// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include <string>

namespace eicrecon {

/// Configuration for LGADChargeSharingRecon.
///
/// Only user-settable digitization fields live here. Readout pitch and grid
/// coordinates come from segmentation; electrode metal size and sensor
/// thickness come from explicit detector constants in the compact geometry.
struct LGADChargeSharingReconConfig {
    /// DD4hep readout name (required; set per detector).
    std::string readout;

    /// Energy deposit threshold for accepting a hit (GeV).
    float minEDepGeV{0.0F};

    /// Neighborhood half-width. 2 -> 5x5 grid.
    int neighborhoodRadius{2};

    /// LogA model: transverse hit size d0 (microns).
    double d0Micron{1.0};

    /// Silicon e/h pair ionization energy (eV).
    double ionizationEnergyEV{3.6};

    /// AC-LGAD amplification factor (gain).
    double amplificationFactor{20.0};

    /// Enable per-pad gain variation + electronic noise injection.
    bool noiseEnabled{true};

    /// Electronic noise RMS (electrons).
    double noiseElectronCount{500.0};

    /// Charge-integration window (ns). SimTrackerHit contributions landing on the
    /// same pad within this window are summed into one channel before electronics
    /// are applied. Deposits separated by more than this become distinct pulses.
    double integrationWindowNs{1.0};

    /// Per-channel readout threshold (electrons). A channel is dropped when its
    /// digitized charge (after noise) falls at or below this value. 0 keeps every
    /// channel with positive charge.
    double thresholdElectrons{0.0};
};

} // namespace eicrecon
