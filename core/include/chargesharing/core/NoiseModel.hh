// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

/// @file NoiseModel.hh
/// @brief Realistic noise injection for charge values (per-pixel gain + electronic noise).

#ifndef CHARGESHARING_CORE_NOISEMODEL_HH
#define CHARGESHARING_CORE_NOISEMODEL_HH

#include <random>
#include <vector>

namespace chargesharing::core {

struct NeighborhoodResult; // forward decl

/// @brief Configuration for realistic noise injection into charge values.
///
/// Models two physical effects: per-pixel gain non-uniformity (multiplicative)
/// and front-end electronic noise (additive, in electron-equivalent units).
struct NoiseConfig {
    bool enabled{true};
    double gainSigmaMin{0.01};                ///< Minimum per-pixel gain variation (1%)
    double gainSigmaMax{0.05};                ///< Maximum per-pixel gain variation (5%)
    double electronNoiseCount{500.0};         ///< Electronic noise RMS (electrons)
    double elementaryCharge{1.602176634e-19}; ///< For converting electrons to Coulombs
};

/// @brief Noise model for realistic charge simulation.
///
/// Applies per-pixel gain variations (multiplicative Gaussian) and additive
/// electronic noise to charge values in Coulombs. Use one instance per thread.
///
/// @note Thread Safety: This class contains an std::mt19937 generator and
/// distribution objects whose operator() mutates internal state. Instances
/// are NOT thread-safe; each thread must own its own NoiseModel.
class NoiseModel {
public:
    /// Fixed default seed: reconstruction output must be reproducible across
    /// runs (EICrecon requirement), so the default constructor must never
    /// draw entropy from std::random_device. Use setSeed()/the seeded
    /// constructor for statistically independent replicas.
    static constexpr unsigned int kDefaultSeed = 20260710U;

    NoiseModel();
    explicit NoiseModel(unsigned int seed);

    void setConfig(const NoiseConfig& config) { m_config = config; }
    const NoiseConfig& config() const { return m_config; }

    void setSeed(unsigned int seed) { m_generator.seed(seed); }

    /// Apply noise to a single charge value (in Coulombs)
    /// @param chargeC Input charge in Coulombs
    /// @return Noisy charge in Coulombs
    double applyNoise(double chargeC);

    /// Apply noise to all pixels in a neighborhood result
    /// Updates the charge field of each pixel
    void applyNoise(NeighborhoodResult& result);

    /// Apply noise to a vector of charges (in Coulombs)
    void applyNoise(std::vector<double>& charges);

private:
    NoiseConfig m_config{};
    std::mt19937 m_generator;
    std::normal_distribution<double> m_gaussDist{0.0, 1.0};
    std::uniform_real_distribution<double> m_uniformDist{0.0, 1.0};
};

} // namespace chargesharing::core

#endif // CHARGESHARING_CORE_NOISEMODEL_HH
