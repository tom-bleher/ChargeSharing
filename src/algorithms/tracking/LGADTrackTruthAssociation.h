// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/MCRecoTrackParticleAssociationCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/TrackCollection.h>

#include <string_view>

namespace eicrecon {

using LGADTrackTruthAssociationAlgorithm =
    algorithms::Algorithm<
        algorithms::Input<edm4eic::TrackCollection,
                          edm4eic::MCRecoTrackerHitAssociationCollection>,
        algorithms::Output<edm4eic::MCRecoTrackParticleAssociationCollection>>;

/// Produces generated-ancestor track truth associations from charge-weighted
/// raw-hit contributors. Reconstruction remains truth-blind; this algorithm
/// only defines the validation and matching labels after the fit.
class LGADTrackTruthAssociation : public LGADTrackTruthAssociationAlgorithm {
public:
    LGADTrackTruthAssociation(std::string_view name)
        : LGADTrackTruthAssociationAlgorithm{
              name, {"inputTracks", "inputRawTrackerHitAssociations"},
              {"outputTrackAssociations"},
              "Charge-weighted generated-ancestor track truth associations"} {}

    void init() final {}
    void process(const Input&, const Output&) const final;
};

} // namespace eicrecon
