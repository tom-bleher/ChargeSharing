// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include "LGADGaussianClusteringConfig.h"

#include <DD4hep/Detector.h>
#include <DD4hep/Segmentations.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DDSegmentation/BitFieldCoder.h>
#include <DDSegmentation/Segmentation.h>
#include <algorithms/algorithm.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackerHitCollection.h>

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "algorithms/interfaces/WithPodConfig.h"
#include "algorithms/tracking/ActsGeometryProvider.h"

namespace eicrecon {

using LGADGaussianClusteringAlgorithm =
    algorithms::Algorithm<algorithms::Input<edm4eic::TrackerHitCollection>,
                          algorithms::Output<edm4eic::Measurement2DCollection>>;

/// Union-find hit clustering with Gaussian-fit sub-pad position extraction.
///
/// Drop-in replacement for upstream LGADHitClustering with the max-ADC /
/// centroid recipe swapped for the shared Gaussian fitter under
/// chargesharing::fit. Uses the same DD4hep segmentation neighbour lookup
/// as upstream.
class LGADGaussianClustering : public LGADGaussianClusteringAlgorithm,
                               public WithPodConfig<LGADGaussianClusteringConfig> {
public:
    LGADGaussianClustering(std::string_view name)
        : LGADGaussianClusteringAlgorithm{name, {"inputHits"}, {"outputClusters"},
                                          "AC-LGAD charge-sharing clustering with Gaussian fit."} {}

    void init() final;
    void process(const Input&, const Output&) const final;

    /// Result of a pure-math cluster position reconstruction, used by tests.
    struct ClusterPosition {
        double reconX{0.0};
        double reconY{0.0};
        double sigma2X{0.0};
        double sigma2Y{0.0};
        bool fitConverged{false};
    };

    /// Pure-math cluster position reconstruction. Extracted from the main
    /// process() loop so it can be unit-tested without DD4hep / Acts services.
    ///
    /// Runs the single 2D Gaussian estimator, with centroid fallback when the
    /// fit does not converge.
    /// Coordinates, pitches and charges are all in millimetres / arbitrary
    /// charge units; the caller is responsible for unit conversions.
    static ClusterPosition reconstructClusterPosition(
        const std::vector<double>& xPos, const std::vector<double>& yPos,
        const std::vector<double>& charges, double centerX, double centerY,
        double maxEdep, double pitchX, double pitchY, double fitErrorPercent);

    /// Disjoint-set data structure used to merge neighbour hits into clusters.
    ///
    /// Exposed here (rather than as a private nested class) so unit tests can
    /// exercise the merge/find logic directly without a DD4hep geometry.
    class UnionFind {
        std::vector<int> m_parent, m_rank;

    public:
        explicit UnionFind(int n);
        int find(int id);
        void merge(int id1, int id2);
    };

private:
    void reconstructCluster(const Output& output, const std::vector<edm4eic::TrackerHit>& hits) const;
    const dd4hep::DDSegmentation::Segmentation*
    getLocalSegmentation(const dd4hep::rec::CellID& cellID) const;

    const dd4hep::rec::CellIDPositionConverter* m_converter = nullptr;
    const dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;
    const dd4hep::Detector* m_detector = nullptr;
    dd4hep::Segmentation m_seg;
    std::shared_ptr<const ActsGeometryProvider> m_acts_context;

    mutable std::unordered_map<const dd4hep::DetElement*,
                               const dd4hep::DDSegmentation::Segmentation*>
        m_segmentation_map;
};

} // namespace eicrecon
