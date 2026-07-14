// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "LGADGaussianClustering.h"

#include "chargesharing/fit/GaussianFit.hh"

#include <DD4hep/Detector.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DDSegmentation/MultiSegmentation.h>

#include <algorithms/geo.h>
#include <algorithms/interfaces/ActsSvc.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <edm4hep/utils/vector_utils.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace eicrecon {

namespace cfit = ::chargesharing::fit;

// ============================================================================
// Union-Find
// ============================================================================

LGADGaussianClustering::UnionFind::UnionFind(int n) : m_parent(n), m_rank(n, 0) {
    std::iota(m_parent.begin(), m_parent.end(), 0);
}

int LGADGaussianClustering::UnionFind::find(int id) {
    if (m_parent[id] == id)
        return id;
    return m_parent[id] = find(m_parent[id]);
}

void LGADGaussianClustering::UnionFind::merge(int id1, int id2) {
    int root1 = find(id1);
    int root2 = find(id2);
    if (root1 != root2) {
        if (m_rank[root1] > m_rank[root2])
            m_parent[root2] = root1;
        else if (m_rank[root1] < m_rank[root2])
            m_parent[root1] = root2;
        else {
            m_parent[root1] = root2;
            m_rank[root2]++;
        }
    }
}

// ============================================================================
// MultiSegmentation resolution
// ============================================================================

const dd4hep::DDSegmentation::Segmentation*
LGADGaussianClustering::getLocalSegmentation(const dd4hep::rec::CellID& cellID) const {
    auto segmentation_type = m_seg.type();
    const dd4hep::DDSegmentation::Segmentation* segmentation = m_seg.segmentation();
    while (segmentation_type == "MultiSegmentation") {
        const auto* multi =
            dynamic_cast<const dd4hep::DDSegmentation::MultiSegmentation*>(segmentation);
        segmentation = &multi->subsegmentation(cellID);
        segmentation_type = segmentation->type();
    }
    if (segmentation_type != "CartesianGridXY" && segmentation_type != "CartesianGridXZ") {
        throw std::runtime_error("LGADGaussianClustering: segmentation is neither "
                                 "CartesianGridXY nor CartesianGridXZ");
    }
    return segmentation;
}

// ============================================================================
// init()
// ============================================================================

void LGADGaussianClustering::init() {
    m_converter = algorithms::GeoSvc::instance().cellIDPositionConverter();
    m_detector = algorithms::GeoSvc::instance().detector();
    m_seg = m_detector->readout(m_cfg.readout).segmentation();
    m_decoder = m_seg.decoder();
    m_acts_context = algorithms::ActsSvc::instance().acts_geometry_provider();

    info("LGADGaussianClustering initialized: readout={}, timeResolutionNs={} ns",
         m_cfg.readout, m_cfg.timeResolutionNs);
}

// ============================================================================
// process()
// ============================================================================

void LGADGaussianClustering::process(const Input& input, const Output& output) const {
    const auto [hits] = input;
    auto [clusters] = output;

    if (!hits || !clusters)
        return;

    const int nHits = static_cast<int>(hits->size());
    if (nHits == 0)
        return;

    std::unordered_map<dd4hep::rec::CellID, std::vector<int>> hitsByCellID;
    for (int i = 0; i < nHits; ++i) {
        hitsByCellID[(*hits)[i].getCellID()].push_back(i);
    }

    UnionFind uf(nHits);

    for (const auto& [cellID, hitIndices] : hitsByCellID) {
        for (std::size_t a = 0; a < hitIndices.size(); ++a) {
            for (std::size_t b = a + 1; b < hitIndices.size(); ++b) {
                const int i = hitIndices[a];
                const int j = hitIndices[b];
                uf.merge(i, j);
            }
        }

        const auto* seg = getLocalSegmentation(cellID);
        std::set<dd4hep::rec::CellID> neighborCells;
        seg->neighbours(cellID, neighborCells);

        for (const auto& neighborID : neighborCells) {
            auto it = hitsByCellID.find(neighborID);
            if (it == hitsByCellID.end())
                continue;
            for (int i : hitIndices) {
                for (int j : it->second) {
                    uf.merge(i, j);
                }
            }
        }
    }

    std::unordered_map<int, std::vector<edm4eic::TrackerHit>> clusterMap;
    for (int i = 0; i < nHits; ++i) {
        clusterMap[uf.find(i)].push_back((*hits)[i]);
    }

    for (const auto& [root, clusterHits] : clusterMap) {
        reconstructCluster(output, clusterHits);
    }
}

// ============================================================================
// Cluster position reconstruction
// ============================================================================

void LGADGaussianClustering::reconstructCluster(const Output& output,
                                                const std::vector<edm4eic::TrackerHit>& hits) const {
    auto [clusters] = output;

    if (hits.empty())
        return;

    const auto seedCellID = hits.front().getCellID();
    const auto* context = m_converter->findContext(seedCellID);
    if (!context) {
        error("No DetElement context for cellID {:#018x}", seedCellID);
        return;
    }
    const auto& surfaceMap = m_acts_context->surfaceMap();
    const auto is = surfaceMap.find(context->identifier);
    if (is == surfaceMap.end()) {
        error("No Acts surface found for volume ID {:#018x}", context->identifier);
        return;
    }
    const Acts::Surface* surface = is->second;
    const auto& gctx = m_acts_context->getActsGeometryContext();

    std::vector<double> xPos, yPos, charges;
    xPos.reserve(hits.size());
    yPos.reserve(hits.size());
    charges.reserve(hits.size());

    double maxEdep = 0.0;
    double totalEdep = 0.0;
    dd4hep::rec::CellID maxCellID = 0;
    double centerX = 0.0;
    double centerY = 0.0;
    double earliestTime = std::numeric_limits<double>::max();

    for (const auto& hit : hits) {
        const auto cellID = hit.getCellID();
        const double edep = hit.getEdep();

        const auto& global = hit.getPosition();
        const auto local = surface->globalToLocal(
            gctx, Acts::Vector3{global.x, global.y, global.z}, Acts::Vector3::Zero(),
            1 * Acts::UnitConstants::um);
        if (!local.ok()) {
            warning("Pad center is not on the Acts surface for cellID {:#018x}; dropping cluster",
                    cellID);
            return;
        }

        xPos.push_back(local.value()[0]);
        yPos.push_back(local.value()[1]);
        charges.push_back(edep);
        totalEdep += edep;

        if (edep > maxEdep) {
            maxEdep = edep;
            maxCellID = cellID;
            centerX = local.value()[0];
            centerY = local.value()[1];
        }
        if (hit.getTime() < earliestTime) {
            earliestTime = hit.getTime();
        }
    }

    if (totalEdep <= 0.0 || charges.empty())
        return;

    const auto* seg = getLocalSegmentation(maxCellID);
    const auto cellDimensions = seg->cellDimensions(maxCellID);
    if (cellDimensions.size() < 2) {
        error("Segmentation for cellID {:#018x} returned fewer than two cell dimensions",
              maxCellID);
        return;
    }
    const double pitchX = cellDimensions[0] / dd4hep::mm;
    const double pitchY = cellDimensions[1] / dd4hep::mm;

    const auto clusterPos = reconstructClusterPosition(xPos, yPos, charges, centerX, centerY,
                                                       maxEdep, pitchX, pitchY,
                                                       m_cfg.fitErrorPercent, m_cfg.fitSigmaMM,
                                                       m_cfg.fitFloatSigma);
    double reconX = clusterPos.reconX;
    double reconY = clusterPos.reconY;
    double sigma2X = clusterPos.sigma2X;
    double sigma2Y = clusterPos.sigma2Y;

    // A failed Minuit fit can return non-finite position or a non-positive
    // covariance; either is undefined behaviour inside the Acts CKF
    // measurement selection (observed as a segfault). Fall back to the
    // max-charge pad center with binary-readout errors.
    if (!std::isfinite(reconX) || !std::isfinite(reconY)) {
        reconX = centerX;
        reconY = centerY;
        sigma2X = pitchX * pitchX / 12.0;
        sigma2Y = pitchY * pitchY / 12.0;
    }
    if (!std::isfinite(sigma2X) || sigma2X <= 0.0) {
        sigma2X = pitchX * pitchX / 12.0;
    }
    if (!std::isfinite(sigma2Y) || sigma2Y <= 0.0) {
        sigma2Y = pitchY * pitchY / 12.0;
    }

    auto cluster = clusters->create();
    cluster.setSurface(surface->geometryId().value());
    // Fit coordinates were projected from global pad centers onto this Acts
    // surface, so they are already valid Measurement2D local coordinates.
    cluster.setLoc({static_cast<float>(reconX), static_cast<float>(reconY)});
    cluster.setTime(static_cast<float>(earliestTime));

    const float timeErr = static_cast<float>(m_cfg.timeResolutionNs);
    cluster.setCovariance({static_cast<float>(sigma2X), static_cast<float>(sigma2Y), timeErr * timeErr, 0.0f});

    for (const auto& hit : hits) {
        cluster.addToHits(hit);
        cluster.addToWeights(static_cast<float>(hit.getEdep() / totalEdep));
    }
}

// ============================================================================
// Pure-math cluster position reconstruction (unit-testable)
// ============================================================================

LGADGaussianClustering::ClusterPosition LGADGaussianClustering::reconstructClusterPosition(
    const std::vector<double>& xPos, const std::vector<double>& yPos,
    const std::vector<double>& charges, double centerX, double centerY, double maxEdep,
    double pitchX, double pitchY, double fitErrorPercent, double fitSigmaMM, bool fitFloatSigma) {
    ClusterPosition out{};

    double totalEdep = 0.0;
    for (double q : charges)
        totalEdep += q;

    if (totalEdep <= 0.0 || charges.empty()) {
        out.reconX = centerX;
        out.reconY = centerY;
        out.sigma2X = pitchX * pitchX / 12.0;
        out.sigma2Y = pitchY * pitchY / 12.0;
        return out;
    }

    auto centroid = [&](double& rx, double& ry) {
        rx = 0.0;
        ry = 0.0;
        for (std::size_t i = 0; i < charges.size(); ++i) {
            rx += charges[i] * xPos[i];
            ry += charges[i] * yPos[i];
        }
        rx /= totalEdep;
        ry /= totalEdep;
    };

    if (charges.size() < 3) {
        centroid(out.reconX, out.reconY);
        out.sigma2X = pitchX * pitchX / 12.0;
        out.sigma2Y = pitchY * pitchY / 12.0;
        return out;
    }

    // Canonical estimator: Gaussian 2D.
    const double spacing = std::min(pitchX, pitchY);
    cfit::GaussFit2DConfig cfg;
    cfg.muXLo = centerX - pitchX;
    cfg.muXHi = centerX + pitchX;
    cfg.muYLo = centerY - pitchY;
    cfg.muYHi = centerY + pitchY;
    // Bounds for a floated width: the charge-cloud sigma can sit well below the
    // pitch, so allow down to 0.1*pitch and up to ~1.5*pitch.
    cfg.sigmaLo = 0.1 * spacing;
    cfg.sigmaHi = 1.5 * std::max(pitchX, pitchY);
    cfg.qMax = maxEdep;
    cfg.pixelSpacing = spacing;
    cfg.errorPercent = fitErrorPercent;
    // Calibration width (seed when floated); 0 lets the fitter fall back to 0.5*pitch.
    cfg.sigmaFixed = fitSigmaMM;
    cfg.floatSigma = fitFloatSigma;

    auto result = cfit::fitGaussian2D(xPos, yPos, charges, cfg);
    if (result.converged) {
        out.reconX = result.muX;
        out.reconY = result.muY;
        out.sigma2X = result.muXError * result.muXError;
        out.sigma2Y = result.muYError * result.muYError;
        out.fitConverged = true;
    } else {
        centroid(out.reconX, out.reconY);
        out.sigma2X = pitchX * pitchX / 12.0;
        out.sigma2Y = pitchY * pitchY / 12.0;
    }
    return out;
}

} // namespace eicrecon
