// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include "algorithms/reco/LGADChargeSharingRecon.h"
#include "algorithms/reco/LGADChargeSharingReconConfig.h"

#include "extensions/jana/JOmniFactory.h"
#include "services/algorithms_init/AlgorithmsInit_service.h"

#include <edm4eic/EDM4eicVersion.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>

#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
#include <edm4eic/MCRecoTrackerHitLinkCollection.h>
#endif

#include <memory>
#include <string>

namespace eicrecon {

/// Thin JOmniFactory wrapper around LGADChargeSharingRecon.
///
    /// Exposes only the digitization parameters; everything else is derived from
/// DD4hep inside the algorithm's init().
class LGADChargeSharingRecon_factory
    : public JOmniFactory<LGADChargeSharingRecon_factory, LGADChargeSharingReconConfig> {
public:
    using AlgoT = LGADChargeSharingRecon;

private:
    std::unique_ptr<AlgoT> m_algo;

    PodioInput<edm4hep::SimTrackerHit> m_in_simhits{this};
    PodioOutput<edm4eic::RawTrackerHit> m_out_raw_hits{this};
    PodioOutput<edm4eic::TrackerHit> m_out_hits{this};
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
    PodioOutput<edm4eic::MCRecoTrackerHitLink> m_out_links{this};
#endif
    PodioOutput<edm4eic::MCRecoTrackerHitAssociation> m_out_assocs{this};

    Service<AlgorithmsInit_service> m_algorithms_init{this};

    ParameterRef<std::string> m_readout{this, "readout", config().readout,
                                        "DD4hep readout name for segmentation lookup"};
    ParameterRef<float> m_min_edep{this, "minEDepGeV", config().minEDepGeV,
                                   "Minimum energy deposit threshold (GeV)"};
    ParameterRef<int> m_neighborhood_radius{this, "neighborhoodRadius", config().neighborhoodRadius,
                                            "Neighborhood half-width (2 = 5x5 grid)"};

    ParameterRef<double> m_d0_micron{this, "d0Micron", config().d0Micron,
                                     "LogA model: transverse hit size d0 (um)"};
    ParameterRef<double> m_ionization_energy{this, "ionizationEnergyEV", config().ionizationEnergyEV,
                                             "Silicon e/h pair ionization energy (eV)"};
    ParameterRef<double> m_amplification{this, "amplificationFactor", config().amplificationFactor,
                                         "AC-LGAD amplification factor (gain)"};
    ParameterRef<bool> m_noise_enabled{this, "noiseEnabled", config().noiseEnabled,
                                       "Enable electronic noise injection"};
    ParameterRef<double> m_noise_electrons{this, "noiseElectronCount", config().noiseElectronCount,
                                           "Electronic noise RMS (electrons)"};
    ParameterRef<double> m_threshold_electrons{this, "thresholdElectrons", config().thresholdElectrons,
                                               "Per-channel readout threshold (electrons)"};

public:
    void Configure() {
        m_algo = std::make_unique<AlgoT>(GetPrefix());
        m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
        m_algo->applyConfig(config());
        m_algo->init();
    }

    void ChangeRun(int32_t /* run_number */) {}

    void Process(int32_t /* run_number */, uint64_t /* event_number */) {
        m_algo->process({m_in_simhits()},
                        {m_out_raw_hits().get(), m_out_hits().get(),
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
                         m_out_links().get(),
#endif
                         m_out_assocs().get()});
    }
};

} // namespace eicrecon
