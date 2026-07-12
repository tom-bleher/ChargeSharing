// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include "algorithms/tracking/LGADGaussianClustering.h"
#include "algorithms/tracking/LGADGaussianClusteringConfig.h"

#include "extensions/jana/JOmniFactory.h"
#include "services/algorithms_init/AlgorithmsInit_service.h"

#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackerHitCollection.h>

#include <memory>
#include <string>

namespace eicrecon {

class LGADGaussianClustering_factory
    : public JOmniFactory<LGADGaussianClustering_factory, LGADGaussianClusteringConfig> {
public:
    using AlgoT = LGADGaussianClustering;

private:
    std::unique_ptr<AlgoT> m_algo;

    PodioInput<edm4eic::TrackerHit> m_hits_input{this};
    PodioOutput<edm4eic::Measurement2D> m_clusters_output{this};

    Service<AlgorithmsInit_service> m_algorithms_init{this};

    ParameterRef<std::string> m_readout{this, "readout", config().readout,
                                        "DD4hep readout name for segmentation lookup"};
    ParameterRef<double> m_deltaT{this, "deltaT", config().deltaT,
                                  "Time gate (ns) for union-find neighbour merge"};
    ParameterRef<double> m_fitErrorPercent{this, "fitErrorPercent", config().fitErrorPercent,
                                           "Fit uncertainty as percentage of cluster max charge"};

public:
    void Configure() {
        m_algo = std::make_unique<AlgoT>(GetPrefix());
        m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
        m_algo->applyConfig(config());
        m_algo->init();
    }

    void ChangeRun(int32_t /* run_number */) {}

    void Process(int32_t /* run_number */, uint64_t /* event_number */) {
        m_algo->process({m_hits_input()}, {m_clusters_output().get()});
    }
};

} // namespace eicrecon
