// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include <algorithms/algorithm.h>
#include <algorithms/interfaces/WithPodConfig.h>

#include "chargesharing/core/ChargeSharingCore.hh"
#include "chargesharing/core/NoiseModel.hh"
#include "LGADChargeSharingReconConfig.h"

#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace dd4hep {
class Detector;
namespace DDSegmentation {
class BitFieldCoder;
class Segmentation;
} // namespace DDSegmentation
namespace rec {
class CellIDPositionConverter;
} // namespace rec
} // namespace dd4hep

namespace eicrecon {

namespace core = ::chargesharing::core;

using LGADChargeSharingReconAlgorithm = algorithms::Algorithm<
    algorithms::Input<edm4hep::SimTrackerHitCollection>,
    algorithms::Output<edm4eic::RawTrackerHitCollection, edm4eic::TrackerHitCollection,
                       edm4eic::MCRecoTrackerHitAssociationCollection>>;

/// AC-LGAD charge-sharing digitization: SimTrackerHit -> per-pad TrackerHits.
///
/// Applies the LogA sharing model and emits the induced charge on every active
/// pad. Position extraction happens once, downstream, in clustering.
///
/// Thread safety: m_noise_model owns std::mt19937 state and is NOT safe for
/// concurrent use on a single instance. JOmniFactory creates one algorithm
/// per worker thread which makes this safe inside JANA2. Outside JANA2,
/// callers must own one instance per thread.
class LGADChargeSharingRecon
    : public LGADChargeSharingReconAlgorithm,
      public eicrecon::WithPodConfig<LGADChargeSharingReconConfig> {
public:
    LGADChargeSharingRecon(std::string_view name)
        : LGADChargeSharingReconAlgorithm{
              name,
              {"inputSimTrackerHits"},
              {"outputRawTrackerHits", "outputTrackerHits", "outputHitAssociations"},
              "AC-LGAD per-pad charge-sharing digitization."} {}

    void init() final;
    void process(const Input&, const Output&) const final;

    // ----- helper types kept public for tests -----

    struct SingleHitInput {
        std::array<double, 3> hitPositionMM{};
        std::optional<std::array<double, 3>> pixelHintMM{};
        std::optional<std::pair<int, int>> pixelIndexHint{};
        double energyDepositGeV{0.0};
        std::uint64_t cellID{0};
    };

    struct NeighborData {
        double fraction{std::numeric_limits<double>::quiet_NaN()};
        double chargeC{0.0};
        double distanceMM{0.0};
        double alphaRad{0.0};
        double pixelXMM{0.0};
        double pixelYMM{0.0};
        int pixelId{-1};
        int di{0};
        int dj{0};
    };

    struct SingleHitResult {
        std::array<double, 3> nearestPixelCenterMM{};
        int pixelRowIndex{0};
        int pixelColIndex{0};
        double totalCollectedChargeC{0.0};
        std::vector<NeighborData> neighbors;
        int numActiveNeighbors{0};
        int neighborhoodRadius{0};
    };

    /// Process a single SimTrackerHit (used directly by unit tests).
    SingleHitResult processSingleHit(const SingleHitInput& input) const;

    /// Internal geometry state derived in init(). Exposed for tests that
    /// want to bypass DD4hep and drive the algorithm on a synthetic pad grid.
    struct Geometry {
        double pixelSpacingXMM{0.5};
        double pixelSpacingYMM{0.5};
        double pixelSizeXMM{0.15};
        double pixelSizeYMM{0.15};
        double gridOffsetXMM{0.0};
        double gridOffsetYMM{0.0};
        double detectorThicknessMM{0.05};
        double pixelThicknessMM{0.02};
        double detectorZCenterMM{-10.0};
        int pixelsPerSide{0};
        bool useXZCoordinates{false};
        std::string fieldNameX{"x"};
        std::string fieldNameY{"y"};
        int minIndexX{0};
        int maxIndexX{-1};
        int minIndexY{0};
        int maxIndexY{-1};

        bool hasBoundsX() const { return maxIndexX >= minIndexX; }
        bool hasBoundsY() const { return maxIndexY >= minIndexY; }
    };

    /// Test hook: set geometry directly, bypassing DD4hep segmentation lookup.
    /// Must be called AFTER applyConfig() and BEFORE init(). When set, init()
    /// will skip the GeoSvc segmentation path.
    void setGeometryForTesting(const Geometry& geom) {
        m_geom = geom;
        m_skip_dd4hep_init = true;
    }

private:
    struct PixelLocation {
        std::array<double, 3> center{};
        int indexI{0};
        int indexJ{0};
    };

    SingleHitResult processSingleHitImpl(const SingleHitInput& input) const;

    PixelLocation findNearestPixelFallback(const std::array<double, 3>& positionMM) const;
    PixelLocation pixelLocationFromIndices(int indexI, int indexJ) const;

    static double indexToPosition(int index, double pitch, double offset) {
        return static_cast<double>(index) * pitch + offset;
    }
    static int positionToIndex(double position, double pitch, double offset) {
        return static_cast<int>(std::floor((position + 0.5 * pitch - offset) / pitch));
    }

    /// Derived geometry, populated in init() from DD4hep.
    Geometry m_geom{};

    /// Electronic-noise RNG state. Mutable so it can be perturbed from
    /// logically-const methods; not safe for concurrent access.
    mutable core::NoiseModel m_noise_model{};

    /// Test-only: skip DD4hep segmentation lookup when geometry was pre-seeded.
    bool m_skip_dd4hep_init{false};

    /// DD4hep segmentation state (populated in init() when available).
    const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{nullptr};
    const dd4hep::DDSegmentation::Segmentation* m_segmentation{nullptr};
    const dd4hep::rec::CellIDPositionConverter* m_converter{nullptr};

};

} // namespace eicrecon
