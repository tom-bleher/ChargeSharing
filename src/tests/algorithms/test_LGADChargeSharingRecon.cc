// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <edm4hep/MCParticleCollection.h>

#include "algorithms/reco/LGADChargeSharingRecon.h"
#include "algorithms/tracking/LGADTruthUtils.h"

#include <algorithms/logger.h>

#include <cmath>
#include <memory>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using eicrecon::LGADChargeSharingRecon;
using eicrecon::LGADChargeSharingReconConfig;

namespace {

LGADChargeSharingRecon::Geometry makeTestGeometry(double electrodeSizeMM = 0.15) {
    LGADChargeSharingRecon::Geometry g{};
    g.pixelSpacingXMM = 0.5;
    g.pixelSpacingYMM = 0.5;
    g.pixelSizeXMM = electrodeSizeMM;
    g.pixelSizeYMM = electrodeSizeMM;
    g.detectorThicknessMM = 0.05;
    g.pixelThicknessMM = 0.02;
    g.detectorZCenterMM = -10.0;
    g.numPixelsX = 21;
    g.numPixelsY = 21;
    g.minIndexX = -10;
    g.maxIndexX = 10;
    g.minIndexY = -10;
    g.maxIndexY = 10;
    return g;
}

std::unique_ptr<LGADChargeSharingRecon> makeAlgorithm(int radius = 2,
                                                      double electrodeSizeMM = 0.15) {
    auto algo = std::make_unique<LGADChargeSharingRecon>("test_lgad_digitizer");
    algo->level(algorithms::LogLevel::kError);

    LGADChargeSharingReconConfig cfg;
    cfg.readout = "TestReadout";
    cfg.neighborhoodRadius = radius;
    cfg.d0Micron = 1.0;
    cfg.ionizationEnergyEV = 3.6;
    cfg.amplificationFactor = 20.0;
    cfg.noiseEnabled = false;
    cfg.noiseElectronCount = 0.0;
    algo->applyConfig(cfg);
    algo->setGeometryForTesting(makeTestGeometry(electrodeSizeMM));
    algo->init();
    return algo;
}

LGADChargeSharingRecon::SingleHitInput makeHit(double x, double y, int i = 0, int j = 0) {
    LGADChargeSharingRecon::SingleHitInput hit{};
    hit.hitPositionMM = {x, y, -10.0};
    hit.energyDepositGeV = 1.0e-3;
    hit.pixelIndexHint = std::pair<int, int>{i, j};
    return hit;
}

} // namespace

TEST_CASE("LGADChargeSharingRecon digitizer normalizes pad fractions", "[lgad][digitizer]") {
    const auto result = makeAlgorithm()->processSingleHit(makeHit(0.05, 0.03));

    double fractionSum = 0.0;
    double chargeSum = 0.0;
    for (const auto& pad : result.neighbors) {
        fractionSum += pad.fraction;
        chargeSum += pad.chargeC;
    }

    REQUIRE(result.neighbors.size() == 25);
    REQUIRE(result.numActiveNeighbors == 25);
    REQUIRE_THAT(fractionSum, WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(chargeSum, WithinRel(result.totalCollectedChargeC, 1e-12));
}

TEST_CASE("LGADChargeSharingRecon digitizer preserves spatial charge asymmetry", "[lgad][digitizer]") {
    const auto result = makeAlgorithm()->processSingleHit(makeHit(0.20, 0.0));

    const LGADChargeSharingRecon::NeighborData* left = nullptr;
    const LGADChargeSharingRecon::NeighborData* right = nullptr;
    for (const auto& pad : result.neighbors) {
        if (pad.di == -1 && pad.dj == 0)
            left = &pad;
        if (pad.di == 1 && pad.dj == 0)
            right = &pad;
    }

    REQUIRE(left != nullptr);
    REQUIRE(right != nullptr);
    CHECK(right->chargeC > left->chargeC);
}

TEST_CASE("LGADChargeSharingRecon digitizer supports a one-pad neighborhood", "[lgad][digitizer]") {
    const auto result = makeAlgorithm(0)->processSingleHit(makeHit(0.04, -0.02));

    REQUIRE(result.neighbors.size() == 1);
    CHECK_THAT(result.neighbors.front().fraction, WithinAbs(1.0, 1e-12));
    CHECK(result.neighbors.front().di == 0);
    CHECK(result.neighbors.front().dj == 0);
}

TEST_CASE("LGADChargeSharingRecon digitizer clips neighborhoods at sensor bounds", "[lgad][digitizer]") {
    const auto result = makeAlgorithm()->processSingleHit(makeHit(5.0, 5.0, 10, 10));

    REQUIRE(result.neighbors.size() == 9);
    double fractionSum = 0.0;
    for (const auto& pad : result.neighbors)
        fractionSum += pad.fraction;
    CHECK_THAT(fractionSum, WithinAbs(1.0, 1e-12));
}

TEST_CASE("LGADChargeSharingRecon derives physical bounds for a centered 500 um grid",
          "[lgad][digitizer][geometry]") {
    const auto xRange = LGADChargeSharingRecon::indexRangeForBox(8.0, 0.5, 0.25);
    CHECK(xRange.first == -16);
    CHECK(xRange.second == 15);

    auto hit = makeHit(7.75, 7.75, 15, 15);
    hit.pixelHintMM = std::array<double, 3>{7.75, 7.75, 0.0};
    hit.indexBounds = LGADChargeSharingRecon::SingleHitInput::IndexBounds{-16, 15, -16, 15};

    const auto result = makeAlgorithm()->processSingleHit(hit);
    REQUIRE(result.neighbors.size() == 9);
    double fractionSum = 0.0;
    for (const auto& pad : result.neighbors)
        fractionSum += pad.fraction;
    CHECK_THAT(fractionSum, WithinAbs(1.0, 1e-12));
}

TEST_CASE("LGADChargeSharingRecon LogA charges are finite and non-negative", "[lgad][digitizer]") {
    const auto result = makeAlgorithm()->processSingleHit(makeHit(0.0, 0.0));
    for (const auto& pad : result.neighbors) {
        CHECK(std::isfinite(pad.chargeC));
        CHECK(pad.chargeC >= 0.0);
    }
}

TEST_CASE("LGADChargeSharingRecon keeps electrode width distinct from readout pitch",
          "[lgad][digitizer][geometry]") {
    const auto physical = makeAlgorithm(2, 0.15)->processSingleHit(makeHit(0.20, 0.0));
    const auto fullPitch = makeAlgorithm(2, 0.50)->processSingleHit(makeHit(0.20, 0.0));

    const auto centerFraction = [](const auto& result) {
        for (const auto& pad : result.neighbors) {
            if (pad.di == 0 && pad.dj == 0)
                return pad.fraction;
        }
        return 0.0;
    };

    CHECK(centerFraction(physical) < 0.5);
    CHECK(centerFraction(fullPitch) > 0.99);
}

TEST_CASE("LGADChargeSharingRecon separates encoded indices from sensor-local pad center",
          "[lgad][digitizer][geometry]") {
    auto hit = makeHit(0.24, -0.08, 7, -3);
    hit.pixelHintMM = std::array<double, 3>{0.20, -0.10, 0.0};

    const auto result = makeAlgorithm()->processSingleHit(hit);

    CHECK(result.pixelRowIndex == 7);
    CHECK(result.pixelColIndex == -3);
    CHECK_THAT(result.nearestPixelCenterMM[0], WithinAbs(0.20, 1e-12));
    CHECK_THAT(result.nearestPixelCenterMM[1], WithinAbs(-0.10, 1e-12));
}

// ============================================================================
// Channel aggregation (aggregate-before-electronics)
// ============================================================================

namespace {

LGADChargeSharingRecon::PadContribution makeContribution(std::uint64_t cellID, double timeNs,
                                                         double chargeC, int simHitIndex) {
    LGADChargeSharingRecon::PadContribution c{};
    c.cellID = cellID;
    c.position = edm4hep::Vector3f{0.0f, 0.0f, 0.0f};
    c.timeNs = timeNs;
    c.inducedChargeC = chargeC;
    c.inducedEnergyGeV = chargeC; // proportional; exact factor irrelevant here
    c.simHitIndex = simHitIndex;
    return c;
}

} // namespace

TEST_CASE("aggregateChannels sums coincident same-pad contributions into one channel",
          "[lgad][digitizer][aggregation]") {
    std::vector<LGADChargeSharingRecon::PadContribution> contribs{
        makeContribution(0x10, /*t=*/1.0, /*q=*/3.0, /*sim=*/0),
        makeContribution(0x10, /*t=*/1.2, /*q=*/1.0, /*sim=*/1),
    };

    const auto channels = LGADChargeSharingRecon::aggregateChannels(contribs);

    REQUIRE(channels.size() == 1);
    CHECK(channels[0].cellID == 0x10);
    CHECK_THAT(channels[0].chargeC, WithinRel(4.0, 1e-12));
    CHECK_THAT(channels[0].energyGeV, WithinRel(4.0, 1e-12));
    CHECK_THAT(channels[0].timeNs, WithinRel(1.0, 1e-12));

    REQUIRE(channels[0].contributors.size() == 2);
    double weightSum = 0.0;
    for (const auto& [simIdx, weight] : channels[0].contributors) {
        weightSum += weight;
        if (simIdx == 0)
            CHECK_THAT(weight, WithinRel(0.75, 1e-12));
        if (simIdx == 1)
            CHECK_THAT(weight, WithinRel(0.25, 1e-12));
    }
    CHECK_THAT(weightSum, WithinAbs(1.0, 1e-12));
}

TEST_CASE("aggregateChannels sums same-pad contributions regardless of time",
          "[lgad][digitizer][aggregation]") {
    std::vector<LGADChargeSharingRecon::PadContribution> contribs{
        makeContribution(0x10, /*t=*/1.0, /*q=*/2.0, /*sim=*/0),
        makeContribution(0x10, /*t=*/50.0, /*q=*/5.0, /*sim=*/1),
    };

    const auto channels = LGADChargeSharingRecon::aggregateChannels(contribs);

    REQUIRE(channels.size() == 1);
    CHECK_THAT(channels[0].chargeC, WithinRel(7.0, 1e-12));
    CHECK_THAT(channels[0].timeNs, WithinRel(1.0, 1e-12));
}

TEST_CASE("aggregateChannels keeps distinct pads separate", "[lgad][digitizer][aggregation]") {
    std::vector<LGADChargeSharingRecon::PadContribution> contribs{
        makeContribution(0x10, 1.0, 2.0, 0),
        makeContribution(0x20, 1.0, 3.0, 0),
    };

    const auto channels = LGADChargeSharingRecon::aggregateChannels(contribs);

    REQUIRE(channels.size() == 2);
    CHECK(channels[0].cellID == 0x10);
    CHECK(channels[1].cellID == 0x20);
}

TEST_CASE("generatedAncestor maps a Geant4 daughter to its generated particle",
          "[lgad][truth][secondaries]") {
    edm4hep::MCParticleCollection particles;
    auto primary = particles.create();
    primary.setGeneratorStatus(1);
    auto secondary = particles.create();
    secondary.setGeneratorStatus(0);
    secondary.addToParents(primary);

    CHECK(eicrecon::generatedAncestor(primary) == primary);
    CHECK(eicrecon::generatedAncestor(secondary) == primary);
}
