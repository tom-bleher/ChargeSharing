// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include <JANA/JEventProcessor.h>

#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/Measurement2DCollection.h>

#include <TTree.h>

#include <memory>

class ActsGeometryProvider;

namespace eicrecon {

/// Writes charge-sharing truth residuals to the shared benchmark ROOT file.
class LGADChargeSharingMonitor : public JEventProcessor {
public:
    LGADChargeSharingMonitor();
    ~LGADChargeSharingMonitor() override = default;

    void Init() override;
    void ProcessSequential(const JEvent& event) override;
    void Finish() override;

private:
    void fillData(const edm4eic::Measurement2D& measurement,
                  const edm4eic::MCRecoTrackerHitAssociationCollection& associations);

    std::shared_ptr<const ActsGeometryProvider> m_acts_context;
    TTree* m_tree{nullptr};
    double m_residual_x{0.0};
    double m_residual_y{0.0};
};

} // namespace eicrecon
