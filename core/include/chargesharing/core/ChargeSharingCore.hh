// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

/// @file ChargeSharingCore.hh
/// @brief Core charge sharing algorithms shared between simulation and reconstruction.
///
/// This header provides the LogA charge sharing model calculations
/// that are used both by the Geant4 simulation and the EICrecon plugin.

#ifndef CHARGESHARING_CORE_CHARGESHARINGCORE_HH
#define CHARGESHARING_CORE_CHARGESHARINGCORE_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace chargesharing::core {

/// @brief Active pixel selection mode controlling which pixels form the normalization denominator.
///
/// Different modes trade spatial coverage against signal-to-noise: wider denominators
/// dilute the fractions but reduce sensitivity to edge effects.
enum class ActivePixelMode {
    Neighborhood,   ///< All pixels in neighborhood
    RowCol,         ///< Cross pattern (center row + center column)
    RowCol3x3,      ///< Cross pattern + 3x3 center block
    ChargeBlock2x2,      ///< 4 pixels with highest weight (contiguous)
    ChargeBlock3x3,      ///< 9 pixels with highest weight (contiguous)
    ThresholdAboveNoise  ///< Pads above N×σ_noise; threshold applied post-noise in caller
};

/// @brief Reconstruction method for extracting hit position from charge fractions.
///
/// These algorithms estimate the sub-pixel hit coordinate from the measured
/// (or simulated) charge distribution across the neighborhood.
enum class ReconMethod {
    Centroid,   ///< Charge-weighted centroid (simple)
    Gaussian1D, ///< Fit 1D Gaussians to row and column slices
    Gaussian2D  ///< Fit 2D Gaussian to full neighborhood
};

// ============================================================================
// Constants
// ============================================================================

namespace constants {
constexpr double kMillimeterPerMicron = 1.0e-3;
constexpr double kGuardFactor = 1.0 + 1e-6;     // Prevents log(0) singularity when distance == d0, with minimal bias
constexpr double kMinD0Micron = 1e-6;           // Floor for d0 to avoid division by zero in log(d/d0)
constexpr double kOutOfBoundsFraction = std::numeric_limits<double>::quiet_NaN(); // Sentinel: pixel outside detector bounds
} // namespace constants

// ============================================================================
// Charge Model Calculations
// ============================================================================

/// @brief Calculate distance from hit point to nearest pad metal edge.
///
/// Following Tornago et al. (NIM A 1003, 2021), $d_i$ is the distance from the
/// hit to the nearest point on the pad metal boundary, not the pad center.
/// For a hit inside the pad, returns 0. This gives the correct signal gradient
/// near pad edges where charge sharing contrast is highest.
/// @param dxMM Offset in X from hit to pixel center (mm).
/// @param dyMM Offset in Y from hit to pixel center (mm).
/// @param padHalfWidthMM Half of pad width in X (mm). Pass 0 to fall back to center-to-center distance.
/// @param padHalfHeightMM Half of pad height in Y (mm). Pass 0 to fall back to center-to-center distance.
/// @return Distance to nearest pad edge in mm (>= 0).
inline double calcDistanceToEdge(double dxMM, double dyMM,
                                   double padHalfWidthMM = 0.0, double padHalfHeightMM = 0.0) {
    if (padHalfWidthMM <= 0.0 && padHalfHeightMM <= 0.0) {
        return std::hypot(dxMM, dyMM);
    }
    const double edgeDx = std::max(0.0, std::abs(dxMM) - padHalfWidthMM);
    const double edgeDy = std::max(0.0, std::abs(dyMM) - padHalfHeightMM);
    return std::hypot(edgeDx, edgeDy);
}

/// @brief Calculate pad view angle (alpha) per Cartiglia (TREDI 2020, slide 9).
///
/// alpha = atan( (l/2)*sqrt(2) / ((l/2)*sqrt(2) + d) )
/// where l is the pad side length and d is the distance from the hit to the
/// nearest pad metal edge (Tornago et al., arXiv:2007.09528, Eq. 4).
/// @param edgeDistMM Distance from hit to nearest pad edge (mm).
/// @param padWidthMM Pad width (mm).
/// @param padHeightMM Pad height (mm).
/// @return View angle alpha (radians), in range (0, pi/4].
inline double calcPadViewAngle(double edgeDistMM, double padWidthMM, double padHeightMM) {
    const double l = (padWidthMM + padHeightMM) / 2.0;
    const double numerator = (l / 2.0) * std::sqrt(2.0);
    const double denominator = numerator + edgeDistMM;
    if (edgeDistMM <= 0.0) {
        return std::atan(1.0); // pi/4 -- hit on pad metal
    }
    return std::atan(numerator / denominator);
}

/// Calculate weight for LogA model: w_i = alpha_i / ln(d_i/d0)
/// @param distanceMM Distance from hit to nearest pad edge (mm)
/// @param alphaPadViewAngle Pad view angle (radians)
/// @param d0MM Reference distance d0 (mm)
inline double calcWeightLogA(double distanceMM, double alphaPadViewAngle, double d0MM) {
    const double minSafeDistance = d0MM * constants::kGuardFactor;
    const double safeDistance = std::max(distanceMM, minSafeDistance);
    const double logValue = std::log(safeDistance / d0MM);
    if (!(logValue > 0.0 && std::isfinite(logValue))) {
        static thread_local bool warned = false;
        if (!warned) {
            std::cerr << "[ChargeSharingCore] Warning: calcWeightLogA returned 0 due to invalid logValue=" << logValue
                      << " (distance=" << distanceMM << ", d0=" << d0MM << ")\n";
            warned = true;
        }
        return 0.0;
    }
    return alphaPadViewAngle / logValue;
}

// ============================================================================
// Neighborhood Data Structures
// ============================================================================

/// @brief Data for a single pixel in the charge-sharing neighborhood.
///
/// Stores geometry, raw weight, and normalized charge fractions computed
/// under different denominator conventions (full neighborhood, row, col, block).
struct NeighborPixel {
    int di{0};            ///< Offset from center pixel (row)
    int dj{0};            ///< Offset from center pixel (col)
    double centerX{0.0};  ///< Pixel center X (mm)
    double centerY{0.0};  ///< Pixel center Y (mm)
    double distance{0.0}; ///< Distance from hit to pixel center (mm)
    double alpha{0.0};    ///< Pad view angle
    double weight{0.0};   ///< Raw weight (before normalization)
    double fraction{0.0}; ///< Normalized charge fraction (full neighborhood denominator)
    double charge{0.0};   ///< Charge in Coulombs (optional)
    bool inBounds{false}; ///< Within detector bounds

    // Mode-specific fractions (different denominators for diagnostics)
    double fractionRow{0.0};   ///< Fraction using row-only denominator
    double fractionCol{0.0};   ///< Fraction using col-only denominator
    double fractionBlock{0.0}; ///< Fraction using block denominator

};

/// @brief Complete result of a neighborhood charge-sharing calculation.
///
/// Contains per-pixel weights and fractions for an NxN grid centered on the
/// hit pixel. Provides helper accessors for extracting row/column slices
/// used by 1-D Gaussian fitting.
struct NeighborhoodResult {
    std::vector<NeighborPixel> pixels;
    int centerPixelI{0};
    int centerPixelJ{0};
    double centerPixelX{0.0};
    double centerPixelY{0.0};
    double totalWeight{0.0};

    /// Get center row slice (di varies, dj=0)
    std::vector<const NeighborPixel*> getCenterRow() const;

    /// Get center column slice (di=0, dj varies)
    std::vector<const NeighborPixel*> getCenterCol() const;
};

// ============================================================================
// Neighborhood Calculation
// ============================================================================

/// @brief Configuration parameters for the neighborhood charge-sharing calculation.
///
/// Groups detector geometry (pad size, pitch, pixel count) and model parameters
/// (d0 and radius) needed by calculateNeighborhood().
struct NeighborhoodConfig {
    ActivePixelMode activeMode{ActivePixelMode::Neighborhood};
    int radius{2};               ///< Neighborhood half-width (2 = 5x5)
    double pixelSizeMM{0.15};    ///< Pad size (mm)
    double pixelSizeYMM{0.15};   ///< Pad size Y (mm), 0 = same as X
    double pixelSpacingMM{0.5};  ///< Pixel pitch (mm)
    double pixelSpacingYMM{0.5}; ///< Pixel pitch Y (mm), 0 = same as X
    double d0Micron{1.0};        ///< LogA d0 parameter (um)
    int numPixelsX{0};           ///< Total pixels in X (for bounds), 0 = unbounded
    int numPixelsY{0};           ///< Total pixels in Y (for bounds), 0 = unbounded
    int minIndexX{0};            ///< Minimum valid DD4hep index in X (can be negative)
    int minIndexY{0};            ///< Minimum valid DD4hep index in Y (can be negative)
};

/// @brief Calculate charge fractions for the pixel neighborhood around a hit position.
///
/// This is the main entry point for the charge-sharing model. Given a hit position
/// and the struck pixel coordinates, it computes weights and normalized fractions
/// for all pixels within the configured radius, applying the chosen signal model
/// and active-pixel mode.
/// @param hitX Hit X position in local coordinates (mm).
/// @param hitY Hit Y position in local coordinates (mm).
/// @param centerI Center (struck) pixel row index.
/// @param centerJ Center (struck) pixel column index.
/// @param centerX Center pixel X position (mm).
/// @param centerY Center pixel Y position (mm).
/// @param config Neighborhood configuration (model, geometry, radius).
/// @return NeighborhoodResult containing per-pixel fractions and diagnostic slices.
NeighborhoodResult calculateNeighborhood(double hitX, double hitY, int centerI, int centerJ, double centerX,
                                         double centerY, const NeighborhoodConfig& config);

} // namespace chargesharing::core

#endif // CHARGESHARING_CORE_CHARGESHARINGCORE_HH
