// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "chargesharing/core/NoiseModel.hh"
#include "chargesharing/core/ChargeSharingCore.hh"

#include <algorithm>

namespace chargesharing::core {

NoiseModel::NoiseModel() : m_generator(kDefaultSeed) {}

NoiseModel::NoiseModel(unsigned int seed) : m_generator(seed) {}

double NoiseModel::applyNoise(double chargeC) {
    if (!m_config.enabled || chargeC <= 0.0) {
        return chargeC;
    }

    // 1. Apply gain variation (multiplicative)
    //    gain_effective = gain * (1 + gaussian(0, sigma))
    //    where sigma is uniformly sampled from [gainSigmaMin, gainSigmaMax]
    double gainSigma =
        m_uniformDist(m_generator) * (m_config.gainSigmaMax - m_config.gainSigmaMin) + m_config.gainSigmaMin;
    double gainFactor = 1.0 + m_gaussDist(m_generator) * gainSigma;
    gainFactor = std::max(0.0, gainFactor); // Prevent negative gain

    double noisyCharge = chargeC * gainFactor;

    // 2. Apply electronic noise (additive)
    //    noise = gaussian(0, noiseElectrons) * elementaryCharge
    double noiseC = m_gaussDist(m_generator) * m_config.electronNoiseCount * m_config.elementaryCharge;
    noisyCharge += noiseC;

    return std::max(0.0, noisyCharge);
}

void NoiseModel::applyNoise(NeighborhoodResult& result) {
    if (!m_config.enabled) {
        return;
    }

    for (auto& pixel : result.pixels) {
        if (pixel.inBounds && pixel.charge > 0.0) {
            pixel.charge = applyNoise(pixel.charge);
        }
    }
}

void NoiseModel::applyNoise(std::vector<double>& charges) {
    if (!m_config.enabled) {
        return;
    }

    for (auto& q : charges) {
        if (q > 0.0) {
            q = applyNoise(q);
        }
    }
}

} // namespace chargesharing::core
