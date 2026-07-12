// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

/// @file GaussianFit.hh
/// @brief Gaussian fitting routines for charge sharing position reconstruction.
///
/// Provides 1D and 2D Gaussian fitting for extracting hit positions from
/// charge distributions. Adapted from the simulation's GaussianFitter.cc
/// for use in EICrecon reconstruction.

#ifndef CHARGESHARING_FIT_GAUSSIANFIT_HH
#define CHARGESHARING_FIT_GAUSSIANFIT_HH

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace chargesharing::fit {

// ============================================================================
// Constants
// ============================================================================

namespace constants {
constexpr double kInvalidValue = std::numeric_limits<double>::quiet_NaN();
constexpr double kMinAmplitude = 1e-18;
constexpr double kDefaultErrorPercent = 5.0;
} // namespace constants

// ============================================================================
// Distance-Weighted Error Configuration
// ============================================================================

/// @brief Configuration for the distance-weighted error model used in chi-square fits.
///
/// Models the per-pixel fit uncertainty as a function of distance from the
/// cluster center, using a power-law dependence.  In "inverse" mode the
/// central pixels (closest to the hit) receive the smallest errors and
/// therefore the highest weight in the fit, which improves position resolution.
/// All percentage values are fractions of the maximum charge Q_max.
struct DistanceWeightedErrorConfig {
    bool enabled{false};          ///< Master switch for distance-weighted errors
    double scalePixels{1.0};      ///< Distance scale in pixel units
    double exponent{1.0};         ///< Power law exponent
    double floorPercent{1.0};     ///< Minimum error as % of Q_max
    double capPercent{50.0};      ///< Maximum error as % of Q_max
    bool powerInverse{true};      ///< Use inverse power model (high weight at center)
    double pixelSpacing{0.5};     ///< Pixel pitch for distance calculation (mm)
    double truthCenterX{0.0};     ///< True hit position X (for distance calc)
    double truthCenterY{0.0};     ///< True hit position Y (for distance calc)
    bool preferTruthCenter{true}; ///< Use truth position instead of fit estimate
};

// ============================================================================
// Distance-Weighted Error Functions (small, inline)
// ============================================================================

/// Apply sigma bounds (floor and cap) as percentage of max charge
inline double applySigmaBounds(double sigma, double maxCharge, double floorPercent, double capPercent) {
    if (!std::isfinite(sigma) || sigma <= 0.0) {
        return constants::kInvalidValue;
    }
    double floorAbs = 0.0;
    if (std::isfinite(maxCharge) && maxCharge > 0.0 && std::isfinite(floorPercent) && floorPercent > 0.0) {
        floorAbs = floorPercent * 0.01 * maxCharge;
    }
    double capAbs = std::numeric_limits<double>::infinity();
    if (std::isfinite(maxCharge) && maxCharge > 0.0 && std::isfinite(capPercent) && capPercent > 0.0) {
        capAbs = capPercent * 0.01 * maxCharge;
    }
    if (floorAbs > 0.0 && sigma < floorAbs) {
        sigma = floorAbs;
    }
    if (sigma > capAbs) {
        sigma = capAbs;
    }
    return std::isfinite(sigma) && sigma > 0.0 ? sigma : constants::kInvalidValue;
}

/// Compute distance-weighted sigma using power law model
/// sigma = floor * (1 + r)^exponent, clamped to [floor, cap]
inline double distancePowerSigma(double distance, double maxCharge, double pixelSpacing, double scalePixels,
                                 double exponent, double floorPercent, double capPercent) {
    if (!std::isfinite(distance) || distance < 0.0 || !std::isfinite(maxCharge) || maxCharge <= 0.0 ||
        !std::isfinite(pixelSpacing) || pixelSpacing <= 0.0) {
        return constants::kInvalidValue;
    }
    const double distanceScale = std::max(scalePixels * pixelSpacing, 1e-12);
    const double ratio = distance / distanceScale;
    const double baseSigma =
        std::max(0.0, floorPercent) * 0.01 * maxCharge * std::pow(1.0 + ratio, std::max(0.0, exponent));
    return applySigmaBounds(baseSigma, maxCharge, floorPercent, capPercent);
}

/// Compute distance-weighted sigma using inverse power model
/// sigma = cap / (1 + r)^exponent, clamped to [floor, cap]
/// This gives smaller errors (higher weight) at the center
inline double distancePowerSigmaInverse(double distance, double maxCharge, double pixelSpacing, double scalePixels,
                                        double exponent, double floorPercent, double capPercent) {
    if (!std::isfinite(distance) || distance < 0.0 || !std::isfinite(maxCharge) || maxCharge <= 0.0 ||
        !std::isfinite(pixelSpacing) || pixelSpacing <= 0.0) {
        return constants::kInvalidValue;
    }
    const double distanceScale = std::max(scalePixels * pixelSpacing, 1e-12);
    const double ratio = distance / distanceScale;
    double sigmaMax = (capPercent / 100.0) * maxCharge;
    const double denom = std::pow(1.0 + std::max(0.0, ratio), std::max(0.0, exponent));
    const double sigmaRaw = (denom > 0.0) ? (sigmaMax / denom) : sigmaMax;
    return applySigmaBounds(sigmaRaw, maxCharge, floorPercent, capPercent);
}

/// Compute distance-weighted error for a pixel
/// @param pixelX X position of pixel center
/// @param pixelY Y position of pixel center
/// @param config Distance-weighted error configuration
/// @param maxCharge Maximum charge in the distribution
/// @return Sigma (error) value for this pixel
inline double computeDistanceWeightedError(double pixelX, double pixelY, const DistanceWeightedErrorConfig& config,
                                           double maxCharge) {
    if (!config.enabled) {
        return config.floorPercent * 0.01 * maxCharge; // Uniform error
    }

    const double dx = pixelX - config.truthCenterX;
    const double dy = pixelY - config.truthCenterY;
    const double distance = std::sqrt(dx * dx + dy * dy);

    if (config.powerInverse) {
        return distancePowerSigmaInverse(distance, maxCharge, config.pixelSpacing, config.scalePixels, config.exponent,
                                         config.floorPercent, config.capPercent);
    } else {
        return distancePowerSigma(distance, maxCharge, config.pixelSpacing, config.scalePixels, config.exponent,
                                  config.floorPercent, config.capPercent);
    }
}

/// Compute distance-weighted error for 1D fit (distance from center position)
inline double computeDistanceWeightedError1D(double pixelPos, double centerPos,
                                             const DistanceWeightedErrorConfig& config, double maxCharge) {
    if (!config.enabled) {
        return config.floorPercent * 0.01 * maxCharge;
    }

    const double distance = std::abs(pixelPos - centerPos);

    if (config.powerInverse) {
        return distancePowerSigmaInverse(distance, maxCharge, config.pixelSpacing, config.scalePixels, config.exponent,
                                         config.floorPercent, config.capPercent);
    } else {
        return distancePowerSigma(distance, maxCharge, config.pixelSpacing, config.scalePixels, config.exponent,
                                  config.floorPercent, config.capPercent);
    }
}

// ============================================================================
// Fit Result Structures
// ============================================================================

/// @brief Result of a 1D Gaussian + baseline fit to a row or column charge profile.
///
/// The fitted model is A * exp(-0.5*((x - mu)/sigma)^2) + B.
struct GaussFit1DResult {
    bool converged{false};                    ///< True if the minimizer converged
    double A{constants::kInvalidValue};       ///< Peak amplitude above baseline (ADC counts or charge units)
    double mu{constants::kInvalidValue};      ///< Reconstructed hit position along the strip/row (mm)
    double sigma{constants::kInvalidValue};   ///< Gaussian width, related to charge spread (mm)
    double B{constants::kInvalidValue};       ///< Constant baseline offset (electronic pedestal or leakage)
    double chi2{constants::kInvalidValue};    ///< Chi-square of the fit
    double ndf{constants::kInvalidValue};     ///< Number of degrees of freedom (N_points - N_params)
    double muError{constants::kInvalidValue}; ///< 1-sigma uncertainty on mu from the fit covariance matrix (mm)
};

/// @brief Result of a 2D Gaussian + baseline fit to a pixel cluster charge map.
///
/// The fitted model is A * exp(-0.5*(u²/sigmaX² + v²/sigmaY²)) + B,
/// where u = (x-muX)cos(theta) + (y-muY)sin(theta),
///       v = -(x-muX)sin(theta) + (y-muY)cos(theta).
struct GaussFit2DResult {
    bool converged{false};                     ///< True if the minimizer converged
    double A{constants::kInvalidValue};        ///< Peak amplitude above baseline (charge units)
    double muX{constants::kInvalidValue};      ///< Reconstructed hit position X (mm)
    double muY{constants::kInvalidValue};      ///< Reconstructed hit position Y (mm)
    double sigmaX{constants::kInvalidValue};   ///< Gaussian width along major axis (mm)
    double sigmaY{constants::kInvalidValue};   ///< Gaussian width along minor axis (mm)
    double B{constants::kInvalidValue};        ///< Constant baseline offset (charge units)
    double theta{0.0};                         ///< Rotation angle of Gaussian axes (radians, [-pi/2, pi/2])
    double chi2{constants::kInvalidValue};     ///< Chi-square of the fit
    double ndf{constants::kInvalidValue};      ///< Number of degrees of freedom
    double muXError{constants::kInvalidValue}; ///< 1-sigma uncertainty on muX from fit covariance (mm)
    double muYError{constants::kInvalidValue}; ///< 1-sigma uncertainty on muY from fit covariance (mm)
};

/// @brief Combined output of the charge-sharing position reconstruction.
///
/// Holds the final reconstructed (X, Y, Z) position and the underlying
/// fit results from both the 1D row/column projections and the 2D fit.
struct ReconstructionResult {
    double reconX{constants::kInvalidValue}; ///< Reconstructed X position
    double reconY{constants::kInvalidValue}; ///< Reconstructed Y position
    double reconZ{constants::kInvalidValue}; ///< Z position (unchanged from hit)

    GaussFit1DResult fitRowX; ///< Row fit result (for X position)
    GaussFit1DResult fitColY; ///< Column fit result (for Y position)
    GaussFit2DResult fit2D;   ///< 2D fit result (alternative)

    bool valid() const { return std::isfinite(reconX) && std::isfinite(reconY); }
};

// ============================================================================
// Gaussian Functions (exposed for TF1/TF2)
// ============================================================================

/// 1D Gaussian with constant offset: A * exp(-0.5*((x-mu)/sigma)^2) + B
inline double gauss1DPlusB(double* x, double* p) {
    const double A = p[0];
    const double mu = p[1];
    const double sigma = p[2];
    const double B = p[3];
    const double dx = (x[0] - mu) / sigma;
    return A * std::exp(-0.5 * dx * dx) + B;
}

/// 2D Gaussian with rotation and constant offset:
/// A * exp(-0.5 * (u^2/sigx^2 + v^2/sigy^2)) + B
/// where u = (x-mux)*cos(theta) + (y-muy)*sin(theta)
///       v = -(x-mux)*sin(theta) + (y-muy)*cos(theta)
inline double gauss2DPlusB(double* xy, double* p) {
    const double A = p[0];
    const double muX = p[1];
    const double muY = p[2];
    const double sigX = p[3];
    const double sigY = p[4];
    const double B = p[5];
    const double theta = p[6];
    const double cosT = std::cos(theta);
    const double sinT = std::sin(theta);
    const double dxRaw = xy[0] - muX;
    const double dyRaw = xy[1] - muY;
    const double u = ( dxRaw * cosT + dyRaw * sinT) / sigX;
    const double v = (-dxRaw * sinT + dyRaw * cosT) / sigY;
    return A * std::exp(-0.5 * (u * u + v * v)) + B;
}

// ============================================================================
// Utility Functions (small, inline)
// ============================================================================

/// @brief Compute the charge-weighted centroid of a 1D distribution (fallback when fit fails).
/// @param positions Pixel center coordinates (mm).
/// @param charges   Measured charge at each pixel.
/// @param baseline  Pedestal to subtract before weighting (default 0).
/// @return Pair of (centroid position, success flag).
inline std::pair<double, bool> weightedCentroid(const std::vector<double>& positions,
                                                const std::vector<double>& charges, double baseline = 0.0) {

    if (positions.size() != charges.size() || positions.empty()) {
        return {constants::kInvalidValue, false};
    }

    double weightedSum = 0.0;
    double weightTotal = 0.0;

    for (std::size_t i = 0; i < positions.size(); ++i) {
        const double weight = std::max(0.0, charges[i] - baseline);
        if (weight > 0.0) {
            weightTotal += weight;
            weightedSum += weight * positions[i];
        }
    }

    if (weightTotal <= 0.0) {
        return {constants::kInvalidValue, false};
    }

    return {weightedSum / weightTotal, true};
}

/// Estimate initial sigma from weighted variance
inline double estimateSigma(const std::vector<double>& positions, const std::vector<double>& charges, double baseline,
                            double pixelSpacing, double sigmaLo, double sigmaHi) {

    double weightedSum = 0.0;
    double weightTotal = 0.0;

    for (std::size_t i = 0; i < positions.size(); ++i) {
        const double weight = std::max(0.0, charges[i] - baseline);
        if (weight > 0.0) {
            weightTotal += weight;
            weightedSum += weight * positions[i];
        }
    }

    double sigma = constants::kInvalidValue;
    if (weightTotal > 0.0) {
        const double mean = weightedSum / weightTotal;
        double variance = 0.0;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const double weight = std::max(0.0, charges[i] - baseline);
            if (weight > 0.0) {
                const double dx = positions[i] - mean;
                variance += weight * dx * dx;
            }
        }
        variance /= weightTotal;
        sigma = std::sqrt(std::max(variance, 1e-12));
    }

    if (!std::isfinite(sigma) || sigma <= 0.0) {
        sigma = std::max(0.25 * pixelSpacing, 1e-6);
    }

    return std::clamp(sigma, sigmaLo, sigmaHi);
}

// ============================================================================
// 1D Gaussian Fitting
// ============================================================================

/// @brief Configuration controlling the 1D Gaussian fit (bounds, error model, pixel geometry).
struct GaussFit1DConfig {
    double muLo;                                          ///< Lower bound for mu
    double muHi;                                          ///< Upper bound for mu
    double sigmaLo;                                       ///< Lower bound for sigma
    double sigmaHi;                                       ///< Upper bound for sigma
    double qMax;                                          ///< Maximum charge value (for uncertainty estimate)
    double pixelSpacing;                                  ///< Pixel pitch (mm)
    double errorPercent{constants::kDefaultErrorPercent}; ///< Error as % of max (fallback when noise model inactive)
    DistanceWeightedErrorConfig distanceErrorConfig{};    ///< Distance-weighted error config
    double centerPosition{0.0};                           ///< Center position for distance calculation (1D)
    double gainSigma{0.0};          ///< Per-pixel multiplicative gain sigma (0 = use errorPercent instead)
    double noiseElectronSigma{0.0}; ///< Additive electronic noise floor in charge units (Coulombs)
};

/// @brief Fit a 1D Gaussian + baseline to a row or column charge profile to extract hit position.
///
/// Uses ROOT's Minuit2 (Fumili2, falling back to Migrad) to minimize chi-square.
/// Requires at least 3 data points.
///
/// @param positions Pixel center coordinates along the projection axis (mm).
/// @param charges   Integrated charge measured at each pixel (ADC or fC).
/// @param config    Fit bounds, pixel pitch, and error model parameters.
/// @return GaussFit1DResult with fitted parameters; check `converged` before using.
GaussFit1DResult fitGaussian1D(const std::vector<double>& positions, const std::vector<double>& charges,
                               const GaussFit1DConfig& config);

// ============================================================================
// 2D Gaussian Fitting
// ============================================================================

/// @brief Configuration controlling the 2D Gaussian fit (bounds, error model, pixel geometry).
struct GaussFit2DConfig {
    double muXLo, muXHi;     ///< Bounds for mu_x
    double muYLo, muYHi;     ///< Bounds for mu_y
    double sigmaLo, sigmaHi; ///< Bounds for sigma
    double qMax;             ///< Maximum charge value
    double pixelSpacing;     ///< Pixel pitch (mm)
    double errorPercent{constants::kDefaultErrorPercent}; ///< Error as % of max (fallback when noise model inactive)
    DistanceWeightedErrorConfig distanceErrorConfig{}; ///< Distance-weighted error config
    double gainSigma{0.0};          ///< Per-pixel multiplicative gain sigma (0 = use errorPercent instead)
    double noiseElectronSigma{0.0}; ///< Additive electronic noise floor in charge units (Coulombs)
};

/// @brief Fit a rotated 2D Gaussian + baseline to a pixel cluster charge map to extract (X, Y) hit position.
///
/// Uses ROOT's Minuit2 (Fumili2, falling back to Migrad) to minimize chi-square.
/// The Gaussian includes a rotation angle theta. Requires at least 6 pixels (7 parameters).
///
/// @param xPositions X coordinates of pixel centers (mm).
/// @param yPositions Y coordinates of pixel centers (mm).
/// @param charges    Integrated charge measured at each pixel (ADC or fC).
/// @param config     Fit bounds, pixel pitch, and error model parameters.
/// @return GaussFit2DResult with fitted parameters; check `converged` before using.
GaussFit2DResult fitGaussian2D(const std::vector<double>& xPositions, const std::vector<double>& yPositions,
                               const std::vector<double>& charges, const GaussFit2DConfig& config);

} // namespace chargesharing::fit

#endif // CHARGESHARING_FIT_GAUSSIANFIT_HH
