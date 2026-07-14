// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include "algorithms/tracking/LGADTrackTruthAssociation.h"
#include "extensions/jana/JOmniFactory.h"

#include <edm4eic/MCRecoTrackParticleAssociationCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/TrackCollection.h>

#include <memory>

namespace eicrecon {

class LGADTrackTruthAssociation_factory
    : public JOmniFactory<LGADTrackTruthAssociation_factory, NoConfig> {
public:
    using AlgoT = LGADTrackTruthAssociation;

private:
    std::unique_ptr<AlgoT> m_algo;
    PodioInput<edm4eic::Track> m_tracks{this};
    PodioInput<edm4eic::MCRecoTrackerHitAssociation> m_raw_hit_associations{this};
    PodioOutput<edm4eic::MCRecoTrackParticleAssociation> m_track_associations{this};

public:
    void Configure() {
        m_algo = std::make_unique<AlgoT>(GetPrefix());
        m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
        m_algo->init();
    }

    void ChangeRun(int32_t /* runNumber */) {}

    void Process(int32_t /* runNumber */, uint64_t /* eventNumber */) {
        m_algo->process({m_tracks(), m_raw_hit_associations()}, {m_track_associations().get()});
    }
};

} // namespace eicrecon
