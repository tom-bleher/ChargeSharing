// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "chargesharing/fit/GaussianFit.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>

#include <Fit/BinData.h>
#include <Fit/Fitter.h>
#include <Math/MinimizerOptions.h>
#include <Math/WrappedMultiTF1.h>
#include <TF1.h>
#include <TF2.h>

namespace chargesharing::fit {

GaussFit1DResult fitGaussian1D(const std::vector<double>& positions, const std::vector<double>& charges,
                               const GaussFit1DConfig& config) {

    GaussFit1DResult result;

    if (positions.size() != charges.size() || positions.size() < 3) {
        return result;
    }

    // Find min/max for initial estimates
    auto [minIt, maxIt] = std::minmax_element(charges.begin(), charges.end());
    const double qMin = *minIt;
    const double qMax = *maxIt;
    const double A0 = std::max(constants::kMinAmplitude, qMax - qMin);
    const double B0 = qMin;

    // Find position of maximum charge
    const std::size_t maxIdx = std::distance(charges.begin(), maxIt);
    const double mu0 = positions[maxIdx];

    // Estimate initial sigma
    const double sigma0 = estimateSigma(positions, charges, B0, config.pixelSpacing, config.sigmaLo, config.sigmaHi);

    // Calculate base uncertainty
    const bool useNoiseModel = config.gainSigma > 0.0 || config.noiseElectronSigma > 0.0;
    const double uniformSigma = std::max(1e-12, config.qMax * config.errorPercent / 100.0);

    // Compute per-pixel errors
    std::vector<double> perPixelErrors(positions.size());
    const bool useDistanceErrors = config.distanceErrorConfig.enabled;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (useNoiseModel) {
            // Physics-based: sigma_i = sqrt((sigma_gain * Q_i)^2 + sigma_noise^2)
            const double gainTerm = config.gainSigma * charges[i];
            perPixelErrors[i] = std::sqrt(gainTerm * gainTerm +
                                          config.noiseElectronSigma * config.noiseElectronSigma);
            if (!std::isfinite(perPixelErrors[i]) || perPixelErrors[i] <= 0.0) {
                perPixelErrors[i] = std::max(config.noiseElectronSigma, 1e-20);
            }
        } else if (useDistanceErrors) {
            perPixelErrors[i] = computeDistanceWeightedError1D(positions[i], config.centerPosition,
                                                               config.distanceErrorConfig, config.qMax);
            if (!std::isfinite(perPixelErrors[i]) || perPixelErrors[i] <= 0.0) {
                perPixelErrors[i] = uniformSigma;
            }
        } else {
            perPixelErrors[i] = uniformSigma;
        }
    }

    // Set up TF1 for fitting
    thread_local TF1 fitFunc("fitGauss1D", gauss1DPlusB, -1e9, 1e9, 4);

    auto [minPos, maxPos] = std::minmax_element(positions.begin(), positions.end());
    const double margin = 0.5 * config.pixelSpacing;
    fitFunc.SetRange(*minPos - margin, *maxPos + margin);

    const double amplitudeMax = std::max(constants::kMinAmplitude, 2.0 * qMax);
    const double baselineMax = std::max(constants::kMinAmplitude, qMax);

    fitFunc.SetParameters(A0, mu0, sigma0, B0);
    fitFunc.SetParLimits(0, constants::kMinAmplitude, amplitudeMax);
    fitFunc.SetParLimits(1, config.muLo, config.muHi);
    fitFunc.SetParLimits(2, config.sigmaLo, config.sigmaHi);
    fitFunc.SetParLimits(3, -baselineMax, baselineMax);

    // Set up ROOT Fitter
    ROOT::Math::WrappedMultiTF1 wrapped(fitFunc, 1);
    ROOT::Fit::BinData data(static_cast<unsigned int>(positions.size()), 1);

    for (std::size_t i = 0; i < positions.size(); ++i) {
        data.Add(positions[i], charges[i], perPixelErrors[i]);
    }

    ROOT::Fit::Fitter fitter;
    fitter.Config().SetMinimizer("Minuit2", "Fumili2");
    fitter.Config().MinimizerOptions().SetStrategy(0);
    fitter.Config().MinimizerOptions().SetTolerance(1e-4);
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.SetFunction(wrapped);

    fitter.Config().ParSettings(0).SetLimits(constants::kMinAmplitude, amplitudeMax);
    fitter.Config().ParSettings(1).SetLimits(config.muLo, config.muHi);
    fitter.Config().ParSettings(2).SetLimits(config.sigmaLo, config.sigmaHi);
    fitter.Config().ParSettings(3).SetLimits(-baselineMax, baselineMax);
    fitter.Config().ParSettings(0).SetValue(A0);
    fitter.Config().ParSettings(1).SetValue(mu0);
    fitter.Config().ParSettings(2).SetValue(sigma0);
    fitter.Config().ParSettings(3).SetValue(B0);

    bool ok = fitter.Fit(data);
    if (!ok) {
        fitter.Config().SetMinimizer("Minuit2", "Migrad");
        fitter.Config().MinimizerOptions().SetStrategy(1);
        fitter.Config().MinimizerOptions().SetTolerance(1e-3);
        ok = fitter.Fit(data);
    }

    if (ok) {
        result.converged = true;
        result.A = fitter.Result().Parameter(0);
        result.mu = fitter.Result().Parameter(1);
        result.sigma = fitter.Result().Parameter(2);
        result.B = fitter.Result().Parameter(3);
        result.chi2 = fitter.Result().Chi2();
        result.ndf = fitter.Result().Ndf();
        result.muError = fitter.Result().Error(1);
    }

    return result;
}

GaussFit2DResult fitGaussian2D(const std::vector<double>& xPositions, const std::vector<double>& yPositions,
                               const std::vector<double>& charges, const GaussFit2DConfig& config) {

    GaussFit2DResult result;

    const std::size_t nPts = xPositions.size();
    if (nPts != yPositions.size() || nPts != charges.size() || nPts < 6) {
        return result;
    }

    // Find min/max for initial estimates
    auto [minIt, maxIt] = std::minmax_element(charges.begin(), charges.end());
    const double qMin = *minIt;
    const double qMax = *maxIt;
    const double A0 = std::max(constants::kMinAmplitude, qMax - qMin);
    const double B0 = qMin;

    // Find position of maximum charge
    const std::size_t maxIdx = std::distance(charges.begin(), maxIt);
    const double muX0 = xPositions[maxIdx];
    const double muY0 = yPositions[maxIdx];

    // Estimate initial sigmas from weighted variance
    auto estimateSigma2D = [&](bool forX) -> double {
        double wsum = 0.0, m = 0.0;
        for (std::size_t i = 0; i < nPts; ++i) {
            const double w = std::max(0.0, charges[i] - B0);
            const double c = forX ? xPositions[i] : yPositions[i];
            wsum += w;
            m += w * c;
        }
        if (wsum <= 0.0) {
            return std::clamp(0.25 * config.pixelSpacing, config.sigmaLo, config.sigmaHi);
        }
        m /= wsum;
        double var = 0.0;
        for (std::size_t i = 0; i < nPts; ++i) {
            const double w = std::max(0.0, charges[i] - B0);
            const double c = forX ? xPositions[i] : yPositions[i];
            const double d = c - m;
            var += w * d * d;
        }
        var /= wsum;
        return std::clamp(std::sqrt(std::max(var, 1e-12)), config.sigmaLo, config.sigmaHi);
    };

    const double sigX0 = estimateSigma2D(true);
    const double sigY0 = estimateSigma2D(false);

    // Calculate base uncertainty
    const bool useNoiseModel = config.gainSigma > 0.0 || config.noiseElectronSigma > 0.0;
    const double uniformSigma = std::max(1e-12, config.qMax * config.errorPercent / 100.0);

    // Compute per-pixel errors
    std::vector<double> perPixelErrors(nPts);
    const bool useDistanceErrors = config.distanceErrorConfig.enabled;
    for (std::size_t i = 0; i < nPts; ++i) {
        if (useNoiseModel) {
            // Physics-based: sigma_i = sqrt((sigma_gain * Q_i)^2 + sigma_noise^2)
            const double gainTerm = config.gainSigma * charges[i];
            perPixelErrors[i] = std::sqrt(gainTerm * gainTerm +
                                          config.noiseElectronSigma * config.noiseElectronSigma);
            if (!std::isfinite(perPixelErrors[i]) || perPixelErrors[i] <= 0.0) {
                perPixelErrors[i] = std::max(config.noiseElectronSigma, 1e-20);
            }
        } else if (useDistanceErrors) {
            perPixelErrors[i] =
                computeDistanceWeightedError(xPositions[i], yPositions[i], config.distanceErrorConfig, config.qMax);
            if (!std::isfinite(perPixelErrors[i]) || perPixelErrors[i] <= 0.0) {
                perPixelErrors[i] = uniformSigma;
            }
        } else {
            perPixelErrors[i] = uniformSigma;
        }
    }

    // Find range
    auto [xMinIt, xMaxIt] = std::minmax_element(xPositions.begin(), xPositions.end());
    auto [yMinIt, yMaxIt] = std::minmax_element(yPositions.begin(), yPositions.end());
    const double margin = 0.5 * config.pixelSpacing;
    const double xMin = *xMinIt - margin;
    const double xMax = *xMaxIt + margin;
    const double yMin = *yMinIt - margin;
    const double yMax = *yMaxIt + margin;

    // Set up TF2 for fitting (7 params: A, muX, muY, sigX, sigY, B, theta)
    thread_local TF2 fitFunc("fitGauss2D", gauss2DPlusB, -1e9, 1e9, -1e9, 1e9, 7);
    fitFunc.SetRange(xMin, yMin, xMax, yMax);

    const double amplitudeMax = std::max(constants::kMinAmplitude, 2.0 * qMax);
    const double baselineMax = std::max(constants::kMinAmplitude, qMax);
    constexpr double kHalfPi = M_PI / 2.0;

    // Set up ROOT Fitter
    ROOT::Math::WrappedMultiTF1 wrapped(fitFunc, 2);
    ROOT::Fit::BinData data(static_cast<unsigned int>(nPts), 2);

    for (std::size_t i = 0; i < nPts; ++i) {
        double xy[2] = {xPositions[i], yPositions[i]};
        data.Add(xy, charges[i], perPixelErrors[i]);
    }

    ROOT::Fit::Fitter fitter;
    fitter.Config().SetMinimizer("Minuit2", "Fumili2");
    fitter.Config().MinimizerOptions().SetStrategy(0);
    fitter.Config().MinimizerOptions().SetTolerance(1e-4);
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.SetFunction(wrapped);

    fitter.Config().ParSettings(0).SetLimits(constants::kMinAmplitude, amplitudeMax);
    fitter.Config().ParSettings(1).SetLimits(config.muXLo, config.muXHi);
    fitter.Config().ParSettings(2).SetLimits(config.muYLo, config.muYHi);
    fitter.Config().ParSettings(3).SetLimits(config.sigmaLo, config.sigmaHi);
    fitter.Config().ParSettings(4).SetLimits(config.sigmaLo, config.sigmaHi);
    fitter.Config().ParSettings(5).SetLimits(-baselineMax, baselineMax);
    fitter.Config().ParSettings(6).SetLimits(-kHalfPi, kHalfPi);
    fitter.Config().ParSettings(0).SetValue(A0);
    fitter.Config().ParSettings(1).SetValue(muX0);
    fitter.Config().ParSettings(2).SetValue(muY0);
    fitter.Config().ParSettings(3).SetValue(sigX0);
    fitter.Config().ParSettings(4).SetValue(sigY0);
    fitter.Config().ParSettings(5).SetValue(B0);
    fitter.Config().ParSettings(6).SetValue(0.0);

    bool ok = fitter.Fit(data);
    if (!ok) {
        fitter.Config().SetMinimizer("Minuit2", "Migrad");
        fitter.Config().MinimizerOptions().SetStrategy(1);
        fitter.Config().MinimizerOptions().SetTolerance(1e-3);
        ok = fitter.Fit(data);
    }

    if (ok) {
        result.converged = true;
        result.A = fitter.Result().Parameter(0);
        result.muX = fitter.Result().Parameter(1);
        result.muY = fitter.Result().Parameter(2);
        result.sigmaX = fitter.Result().Parameter(3);
        result.sigmaY = fitter.Result().Parameter(4);
        result.B = fitter.Result().Parameter(5);
        result.theta = fitter.Result().Parameter(6);
        result.chi2 = fitter.Result().Chi2();
        result.ndf = fitter.Result().Ndf();
        result.muXError = fitter.Result().Error(1);
        result.muYError = fitter.Result().Error(2);
    }

    return result;
}

} // namespace chargesharing::fit
