// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include <string>

namespace eicrecon {

/// Configuration for LGADGaussianClustering.
///
/// Union-find clustering of TrackerHits followed by Gaussian-fit sub-pad
/// position extraction. Drop-in shape-compatible with upstream
/// LGADHitClustering but with Gaussian reconstruction bolted on.
struct LGADGaussianClusteringConfig {
    /// DD4hep readout name (required; set per detector).
    std::string readout;

    /// Time gate (ns) for union-find neighbour merge.
    double deltaT{1.0};

    /// Fit uncertainty as a percentage of the cluster's max charge.
    double fitErrorPercent{5.0};

    /// Gaussian surrogate width (mm) used as the fixed value, or as the seed when
    /// the width is floated. 0 -> derive geometrically as half the pad pitch.
    double fitSigmaMM{0.0};

    /// Float a single isotropic width in the 2D fit when the cluster has enough
    /// pads to constrain it; otherwise hold it at fitSigmaMM. Disable for
    /// calibration studies that want a strictly fixed width.
    bool fitFloatSigma{true};
};

} // namespace eicrecon
