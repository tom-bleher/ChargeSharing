// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "LGADTrackTruthAssociation.h"

#include "LGADTruthUtils.h"

#include <edm4hep/MCParticle.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace eicrecon {

namespace {

struct MCParticleObjectIDLess {
    bool operator()(const edm4hep::MCParticle& lhs, const edm4hep::MCParticle& rhs) const {
        const auto lhsID = lhs.getObjectID();
        const auto rhsID = rhs.getObjectID();
        return lhsID.collectionID != rhsID.collectionID ? lhsID.collectionID < rhsID.collectionID
                                                        : lhsID.index < rhsID.index;
    }
};

} // namespace

void LGADTrackTruthAssociation::process(const Input& input, const Output& output) const {
    const auto [tracks, rawHitAssociations] = input;
    auto [trackAssociations] = output;

    if (!tracks || !rawHitAssociations || !trackAssociations) {
        return;
    }

    for (const auto& track : *tracks) {
        std::map<edm4hep::MCParticle, double, MCParticleObjectIDLess> ancestorWeights;

        for (const auto& measurement : track.getMeasurements()) {
            std::map<edm4hep::MCParticle, double, MCParticleObjectIDLess> measurementWeights;
            const auto& hits = measurement.getHits();
            const auto& hitWeights = measurement.getWeights();

            for (std::size_t index = 0; index < hits.size(); ++index) {
                const auto rawHit = hits[index].getRawHit();
                if (!rawHit.isAvailable()) {
                    continue;
                }
                const double hitWeight = hitWeights.size() == hits.size() ? hitWeights[index] : 1.0;
                if (!(hitWeight > 0.0) || !std::isfinite(hitWeight)) {
                    continue;
                }

                for (const auto& association : *rawHitAssociations) {
                    if (association.getRawHit() != rawHit || !association.getSimHit().isAvailable()) {
                        continue;
                    }
                    const double associationWeight = std::max(0.0f, association.getWeight());
                    const auto ancestor = generatedAncestor(association.getSimHit().getParticle());
                    if (associationWeight > 0.0 && ancestor.isAvailable()) {
                        measurementWeights[ancestor] += hitWeight * associationWeight;
                    }
                }
            }

            const double measurementTotal = std::accumulate(
                measurementWeights.begin(), measurementWeights.end(), 0.0,
                [](double sum, const auto& entry) { return sum + entry.second; });
            if (!(measurementTotal > 0.0)) {
                continue;
            }
            for (const auto& [ancestor, weight] : measurementWeights) {
                ancestorWeights[ancestor] += weight / measurementTotal;
            }
        }

        const double totalWeight = std::accumulate(
            ancestorWeights.begin(), ancestorWeights.end(), 0.0,
            [](double sum, const auto& entry) { return sum + entry.second; });
        if (!(totalWeight > 0.0)) {
            continue;
        }
        for (const auto& [ancestor, weight] : ancestorWeights) {
            auto association = trackAssociations->create();
            association.setRec(track);
            association.setSim(ancestor);
            association.setWeight(static_cast<float>(weight / totalWeight));
        }
    }
}

} // namespace eicrecon
