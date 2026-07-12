// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover
//
// Unit tests for LGADGaussianClustering. The full process() path needs a
// DD4hep detector + Acts surface map, so these tests drive the algorithm's
// two pieces of pure math directly:
//
//   * UnionFind (disjoint-set) correctness.
//   * reconstructClusterPosition(...) for the canonical 2D Gaussian estimator
//     and its centroid fallback.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "algorithms/tracking/LGADGaussianClustering.h"

#include <cmath>
#include <vector>

using eicrecon::LGADGaussianClustering;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ============================================================================
// UnionFind
// ============================================================================

TEST_CASE("LGADGaussianClustering::UnionFind: singletons are their own root",
          "[lgad][clustering][unionfind]") {
    LGADGaussianClustering::UnionFind uf(5);
    for (int i = 0; i < 5; ++i) {
        CHECK(uf.find(i) == i);
    }
}

TEST_CASE("LGADGaussianClustering::UnionFind: merge produces a single root",
          "[lgad][clustering][unionfind]") {
    LGADGaussianClustering::UnionFind uf(6);
    uf.merge(0, 1);
    uf.merge(1, 2);
    uf.merge(4, 5);

    CHECK(uf.find(0) == uf.find(1));
    CHECK(uf.find(1) == uf.find(2));
    CHECK(uf.find(0) == uf.find(2));
    CHECK(uf.find(4) == uf.find(5));

    // Clusters {0,1,2} and {3} and {4,5} should be disjoint
    CHECK(uf.find(0) != uf.find(3));
    CHECK(uf.find(0) != uf.find(4));
    CHECK(uf.find(3) != uf.find(4));
}

TEST_CASE("LGADGaussianClustering::UnionFind: transitive merges unify components",
          "[lgad][clustering][unionfind]") {
    LGADGaussianClustering::UnionFind uf(4);
    uf.merge(0, 1);
    uf.merge(2, 3);
    uf.merge(1, 2);
    const int root = uf.find(0);
    for (int i = 0; i < 4; ++i) {
        CHECK(uf.find(i) == root);
    }
}

TEST_CASE("reconstructClusterPosition: centroid fallback when < 3 hits",
          "[lgad][clustering][centroid]") {
    const double pitch = 0.5;
    std::vector<double> xPos{0.1, 0.4};
    std::vector<double> yPos{0.0, 0.0};
    std::vector<double> q{1.0, 1.0};

    auto result = LGADGaussianClustering::reconstructClusterPosition(
        xPos, yPos, q, 0.0, 0.0, 1.0, pitch, pitch, 5.0);

    CHECK_THAT(result.reconX, WithinAbs(0.25, 1e-9));
    CHECK_THAT(result.reconY, WithinAbs(0.0, 1e-9));
    CHECK_FALSE(result.fitConverged);
}

TEST_CASE("reconstructClusterPosition: empty input falls through to center",
          "[lgad][clustering][edge]") {
    std::vector<double> empty;
    auto result = LGADGaussianClustering::reconstructClusterPosition(
        empty, empty, empty, /*centerX=*/1.5, /*centerY=*/-0.5, 0.0, 0.5, 0.5,
        5.0);

    CHECK_THAT(result.reconX, WithinAbs(1.5, 1e-12));
    CHECK_THAT(result.reconY, WithinAbs(-0.5, 1e-12));
    CHECK_FALSE(result.fitConverged);
}

// ============================================================================
// Synthetic pad charge distributions
// ============================================================================

namespace {

/// Build a 5x5 Gaussian charge distribution over pitch=0.5 mm pads centered
/// on (cx, cy) with peak at (muX, muY). Returns xPos, yPos, q in mm / charge.
struct GaussianPadGrid {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> q;
    double maxQ{0.0};
};

GaussianPadGrid makeGaussianGrid(double cx, double cy, double muX, double muY, double sigma,
                                 double pitch = 0.5, int radius = 2, double amplitude = 1000.0) {
    GaussianPadGrid g;
    for (int di = -radius; di <= radius; ++di) {
        for (int dj = -radius; dj <= radius; ++dj) {
            const double px = cx + di * pitch;
            const double py = cy + dj * pitch;
            const double rx = (px - muX) / sigma;
            const double ry = (py - muY) / sigma;
            const double qv = amplitude * std::exp(-0.5 * (rx * rx + ry * ry));
            g.x.push_back(px);
            g.y.push_back(py);
            g.q.push_back(qv);
            g.maxQ = std::max(g.maxQ, qv);
        }
    }
    return g;
}

} // namespace

// KNOWN DEFECT: the rotated-2D Minuit fit reconstructs a peak biased toward
// the grid center on this input (reconX ~ -0.04 vs truth -0.12). Kept visible
// as [!mayfail] until the estimator rework (model-fraction chi2/LUT primary,
// Gaussian as comparator) replaces the current fit.
TEST_CASE("reconstructClusterPosition: Gaussian2D recovers peak of wide grid",
          "[lgad][clustering][gauss2d][!mayfail]") {
    const double pitch = 0.5;
    const double muX = -0.12;
    const double muY = 0.07;
    const double sigma = 0.18;

    auto g = makeGaussianGrid(/*cx=*/0.0, /*cy=*/0.0, muX, muY, sigma, pitch);

    auto result = LGADGaussianClustering::reconstructClusterPosition(
        g.x, g.y, g.q, 0.0, 0.0, g.maxQ, pitch, pitch, 5.0);

    CHECK_THAT(result.reconX, WithinAbs(muX, 0.1 * pitch));
    CHECK_THAT(result.reconY, WithinAbs(muY, 0.1 * pitch));
    CHECK(result.sigma2X > 0.0);
    CHECK(result.sigma2Y > 0.0);
}

// KNOWN DEFECT: on flat (peakless) input the Minuit fit reports success with a
// spurious offset instead of triggering the centroid fallback. Kept visible as
// [!mayfail] until the estimator rework adds an explicit fit-quality gate.
TEST_CASE("reconstructClusterPosition: Gaussian2D falls back to centroid on unfittable input",
          "[lgad][clustering][gauss2d][fallback][!mayfail]") {
    const double pitch = 0.5;
    // All-equal charges with no distinguishable peak -> fit should fail and
    // centroid fallback should give the geometric mean of pad centers.
    std::vector<double> x{-pitch, 0.0, pitch, -pitch, 0.0, pitch};
    std::vector<double> y{-pitch, -pitch, -pitch, pitch, pitch, pitch};
    std::vector<double> q(6, 1.0);

    auto result = LGADGaussianClustering::reconstructClusterPosition(
        x, y, q, 0.0, 0.0, 1.0, pitch, pitch, 5.0);

    // Uniform weights -> centroid is mean of positions which is (0, 0).
    CHECK_THAT(result.reconX, WithinAbs(0.0, 1e-6));
    CHECK_THAT(result.reconY, WithinAbs(0.0, 1e-6));
    CHECK(result.sigma2X > 0.0);
    CHECK(result.sigma2Y > 0.0);
}
