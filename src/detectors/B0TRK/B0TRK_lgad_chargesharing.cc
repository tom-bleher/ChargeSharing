// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include <JANA/JApplication.h>

#include <edm4eic/EDM4eicVersion.h>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/reco/LGADChargeSharingRecon_factory.h"
#include "factories/tracking/LGADGaussianClustering_factory.h"

// Upstream EICrecon tracking factories, reused for a parallel B0 CKF chain
// that consumes the charge-sharing Measurement2D instead of the standard
// pitch/sqrt(12) B0TrackerMeasurements.
#include "factories/tracking/ActsToTracks_factory.h"
#include "factories/tracking/AmbiguitySolver_factory.h"
#include "factories/tracking/CKFTracking_factory.h"

extern "C" {

void InitPlugin(JApplication* app) {
    InitJANAPlugin(app);

    // -----------------------------------------------------------------------
    // B0 tracker: SimTrackerHit -> RawTrackerHit + TrackerHit (+ association)
    // -----------------------------------------------------------------------
    eicrecon::LGADChargeSharingReconConfig b0_cfg;
    b0_cfg.readout = "B0TrackerHits";

    app->Add(new JOmniFactoryGeneratorT<eicrecon::LGADChargeSharingRecon_factory>(
        "B0TrackerChargeSharingHitReco",
        {"B0TrackerHits"},
        {"B0TrackerChargeSharingRawHits", "B0TrackerChargeSharingHits",
         "B0TrackerChargeSharingHitAssociations"},
        b0_cfg, app));

    // -----------------------------------------------------------------------
    // B0 tracker: TrackerHit -> Measurement2D (cluster positions)
    // -----------------------------------------------------------------------
    eicrecon::LGADGaussianClusteringConfig b0_cluster_cfg;
    b0_cluster_cfg.readout = "B0TrackerHits";
    app->Add(new JOmniFactoryGeneratorT<eicrecon::LGADGaussianClustering_factory>(
        "B0TrackerChargeSharingClustering",
        {"B0TrackerChargeSharingHits"},
        {"B0TrackerClusterHits"},
        b0_cluster_cfg, app));

    // -----------------------------------------------------------------------
    // Parallel truth-seeded CKF chain on the charge-sharing measurements.
    // Mirrors upstream's B0TrackerCKFTruthSeeded* chain (tracking.cc) with
    // "CS" collection names so both can coexist in one reconstruction pass.
    // -----------------------------------------------------------------------
    app->Add(new JOmniFactoryGeneratorT<eicrecon::CKFTracking_factory>(
        "B0TrackerCSCKFTruthSeededTrajectories",
        {"B0TrackerTruthSeeds", "B0TrackerClusterHits"},
        {
            "B0TrackerCSCKFTruthSeededActsTrackStatesUnfiltered",
            "B0TrackerCSCKFTruthSeededActsTracksUnfiltered",
        },
        {
            .numMeasurementsMin = 3,
        },
        app));

    app->Add(new JOmniFactoryGeneratorT<eicrecon::ActsToTracks_factory>(
        "B0TrackerCSCKFTruthSeededTracksUnfiltered",
        {
            "B0TrackerClusterHits",
            "B0TrackerTruthSeeds",
            "B0TrackerCSCKFTruthSeededActsTrackStatesUnfiltered",
            "B0TrackerCSCKFTruthSeededActsTracksUnfiltered",
            "B0TrackerChargeSharingHitAssociations",
        },
        {
            "B0TrackerCSCKFTruthSeededTrajectoriesUnfiltered",
            "B0TrackerCSCKFTruthSeededTrackParametersUnfiltered",
            "B0TrackerCSCKFTruthSeededTracksUnfiltered",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
            "B0TrackerCSCKFTruthSeededTrackUnfilteredLinks",
#endif
            "B0TrackerCSCKFTruthSeededTrackUnfilteredAssociations",
        },
        app));

    app->Add(new JOmniFactoryGeneratorT<eicrecon::AmbiguitySolver_factory>(
        "B0TrackerCSTruthSeededAmbiguityResolutionSolver",
        {"B0TrackerCSCKFTruthSeededActsTrackStatesUnfiltered",
         "B0TrackerCSCKFTruthSeededActsTracksUnfiltered"},
        {
            "B0TrackerCSCKFTruthSeededActsTrackStates",
            "B0TrackerCSCKFTruthSeededActsTracks",
        },
        {
            .n_measurements_min = 3,
        },
        app));

    app->Add(new JOmniFactoryGeneratorT<eicrecon::ActsToTracks_factory>(
        "B0TrackerCSCKFTruthSeededTracks",
        {
            "B0TrackerClusterHits",
            "B0TrackerTruthSeeds",
            "B0TrackerCSCKFTruthSeededActsTrackStates",
            "B0TrackerCSCKFTruthSeededActsTracks",
            "B0TrackerChargeSharingHitAssociations",
        },
        {
            "B0TrackerCSCKFTruthSeededTrajectories",
            "B0TrackerCSCKFTruthSeededTrackParameters",
            "B0TrackerCSCKFTruthSeededTracks",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
            "B0TrackerCSCKFTruthSeededTrackLinks",
#endif
            "B0TrackerCSCKFTruthSeededTrackAssociations",
        },
        app));
}

} // extern "C"
