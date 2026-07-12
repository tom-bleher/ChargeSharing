// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "algorithms/reco/LGADChargeSharingRecon.h"

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
    g.pixelsPerSide = 21;
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
