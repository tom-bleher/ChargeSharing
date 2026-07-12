// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "LGADChargeSharingMonitor.h"

#include <algorithms/interfaces/ActsSvc.h>
#include <algorithms/tracking/ActsGeometryProvider.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <JANA/JApplication.h>
#include <services/log/Log_service.h>
#include <services/rootfile/RootFile_service.h>

#include <TDirectory.h>
#include <TFile.h>

#include <algorithm>
#include <cmath>

namespace eicrecon {

namespace {
constexpr const char* kMonitorDirName = "LGADChargeSharing";
constexpr const char* kMeasurementCollection = "B0TrackerClusterHits";
constexpr const char* kAssociationCollection = "B0TrackerChargeSharingHitAssociations";
}

LGADChargeSharingMonitor::LGADChargeSharingMonitor() {
    SetTypeName(NAME_OF_THIS);
    // ProcessSequential(const JEvent&) is only dispatched in ExpertMode; the
    // JANA default (LegacyMode) calls the no-op Process() and fills nothing.
    SetCallbackStyle(CallbackStyle::ExpertMode);
}

void LGADChargeSharingMonitor::Init() {
    auto app = GetApplication();
    m_acts_context = algorithms::ActsSvc::instance().acts_geometry_provider();

    // Use the shared benchmark ROOT file activated with -Phistsfile=....
    auto rootfile_svc = app->GetService<RootFile_service>();
    auto* rootfile = rootfile_svc->GetHistFile();

    TDirectory* mainDir = rootfile->mkdir(kMonitorDirName);
    if (mainDir == nullptr) {
        mainDir = rootfile->GetDirectory(kMonitorDirName);
    }
    mainDir->cd();
    m_tree = new TTree("hits", "LGAD charge sharing residuals");

    m_tree->Branch("residualX", &m_residual_x, "residualX/D");
    m_tree->Branch("residualY", &m_residual_y, "residualY/D");
}

void LGADChargeSharingMonitor::ProcessSequential(const JEvent& event) {
    const edm4eic::Measurement2DCollection* measurements = nullptr;
    const edm4eic::MCRecoTrackerHitAssociationCollection* associations = nullptr;
    try {
        measurements = event.GetCollection<edm4eic::Measurement2D>(kMeasurementCollection);
        associations =
            event.GetCollection<edm4eic::MCRecoTrackerHitAssociation>(kAssociationCollection);
    } catch (...) {
        return;
    }
    if (measurements == nullptr || associations == nullptr) {
        return;
    }

    for (const auto& measurement : *measurements) {
        fillData(measurement, *associations);
    }
}

void LGADChargeSharingMonitor::fillData(
    const edm4eic::Measurement2D& measurement,
    const edm4eic::MCRecoTrackerHitAssociationCollection& associations) {
    const auto& hits = measurement.getHits();
    const auto& hitWeights = measurement.getWeights();
    if (hits.empty()) {
        return;
    }

    Acts::Vector3 truthGlobal = Acts::Vector3::Zero();
    double totalWeight = 0.0;

    for (std::size_t i = 0; i < hits.size(); ++i) {
        const auto rawHit = hits[i].getRawHit();
        if (!rawHit.isAvailable()) {
            continue;
        }

        const double hitWeight = hitWeights.size() == hits.size() ? hitWeights[i] : 1.0;
        if (!(hitWeight > 0.0) || !std::isfinite(hitWeight)) {
            continue;
        }

        double associationWeight = 0.0;
        for (const auto& association : associations) {
            if (association.getRawHit() == rawHit && association.getSimHit().isAvailable()) {
                associationWeight += std::max(0.0f, association.getWeight());
            }
        }
        if (!(associationWeight > 0.0)) {
            continue;
        }

        for (const auto& association : associations) {
            if (association.getRawHit() != rawHit || !association.getSimHit().isAvailable()) {
                continue;
            }
            const double weight =
                hitWeight * std::max(0.0f, association.getWeight()) / associationWeight;
            const auto& truth = association.getSimHit().getPosition();
            truthGlobal += weight * Acts::Vector3{truth.x, truth.y, truth.z};
            totalWeight += weight;
        }
    }

    if (!(totalWeight > 0.0) || !m_acts_context) {
        return;
    }
    truthGlobal /= totalWeight;

    const Acts::Surface* surface = nullptr;
    for (const auto& [volumeID, candidate] : m_acts_context->surfaceMap()) {
        (void)volumeID;
        if (candidate != nullptr && candidate->geometryId().value() == measurement.getSurface()) {
            surface = candidate;
            break;
        }
    }
    if (surface == nullptr) {
        return;
    }

    const auto truthLocal = surface->globalToLocal(m_acts_context->getActsGeometryContext(),
                                                   truthGlobal, Acts::Vector3::Zero(),
                                                   1 * Acts::UnitConstants::mm);
    if (!truthLocal.ok()) {
        return;
    }

    m_residual_x = measurement.getLoc().a - truthLocal.value()[0];
    m_residual_y = measurement.getLoc().b - truthLocal.value()[1];

    m_tree->Fill();
}

void LGADChargeSharingMonitor::Finish() {
    auto app = GetApplication();
    auto log = app->GetService<Log_service>()->logger("LGADChargeSharingMonitor");

    log->info("TTree 'hits' contains {} entries", m_tree->GetEntries());
}

} // namespace eicrecon
