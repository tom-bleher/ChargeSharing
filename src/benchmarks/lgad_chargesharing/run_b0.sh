#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# Truth-residual benchmark for the LGAD charge-sharing reconstruction on the
# B0 tracker. Runs the full pipeline:
#   1. Generate forward protons into B0 acceptance
#   2. ddsim Geant4 simulation
#   3. eicrecon with B0TRK_lgad_chargesharing + LGAD_chargesharing_benchmark
#   4. Hand off to validate.py to compute residual statistics and exit with
#      pass/fail against per-detector thresholds
#
# Must be run inside eic-shell with the plugin built + installed.
# Usage: run_b0.sh [NEVENTS] [OUTDIR]
set -e

NEVENTS=${1:-200}
OUTDIR=${2:-/tmp/lgad_chargesharing_bench/b0}
mkdir -p "$OUTDIR"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
TEST_DIR="${REPO_ROOT}/test"
COMPACT="${REPO_ROOT}/test_b0.xml"
PLUGIN_DIR="${LGAD_PLUGIN_DIR:-${REPO_ROOT}/install}"

if [ ! -d "$PLUGIN_DIR/plugins" ]; then
    echo "ERROR: Plugin install directory not found at $PLUGIN_DIR/plugins"
    echo "       Build with: cmake --build build --target install"
    exit 1
fi

echo "=== B0 LGAD charge-sharing benchmark ($NEVENTS events) ==="
echo "  Compact:  $COMPACT"
echo "  Output:   $OUTDIR"
echo ""

echo "=== Step 1: Generate forward protons ==="
python3 "${TEST_DIR}/b0/gen_b0_particles.py" \
    --nevents "$NEVENTS" \
    --output "$OUTDIR/b0_gen.hepmc" \
    --particle proton \
    --Emin 50.0 --Emax 100.0

echo ""
echo "=== Step 2: ddsim simulation ==="
# The minimal test compact defines no tracker_region_* constants, so the
# truth-trimming user particle handler must be disabled (keeps all MC truth).
ddsim --compactFile "$COMPACT" \
      --numberOfEvents "$NEVENTS" \
      --inputFiles "$OUTDIR/b0_gen.hepmc" \
      --part.userParticleHandler='' \
      --outputFile "$OUTDIR/b0_sim.edm4hep.root"

echo ""
echo "=== Step 3: EICrecon (B0TRK plugin + benchmark monitor) ==="
export EICrecon_MY="$PLUGIN_DIR"
OUTPUT_COLLECTIONS="EventHeader,MCParticles,B0TrackerHits,B0TrackerChargeSharingRawHits,B0TrackerChargeSharingHits,B0TrackerChargeSharingHitAssociations,B0TrackerClusterHits,B0TrackerCSCKFTruthSeededTrajectories,B0TrackerCSCKFTruthSeededTrackParameters,B0TrackerCSCKFTruthSeededTracks,B0TrackerCSCKFTruthSeededTrackAssociations,B0TrackerCSCKFTruthSeededTrackRootAssociations"
if [ "${LGAD_ENABLE_MODERN_LINKS:-1}" = "1" ]; then
    OUTPUT_COLLECTIONS="${OUTPUT_COLLECTIONS},B0TrackerChargeSharingRawHitLinks,B0TrackerCSCKFTruthSeededTrackLinks"
fi
eicrecon \
    -Pplugins=B0TRK_lgad_chargesharing,LGAD_chargesharing_benchmark \
    -Pjana:plugin_path="${PLUGIN_DIR}/plugins" \
    -Pnthreads=1 \
    -Pjana:nevents="$NEVENTS" \
    -Ppodio:output_file="$OUTDIR/b0_reco.edm4eic.root" \
    -Ppodio:output_collections="$OUTPUT_COLLECTIONS" \
    -Pdd4hep:xml_files="$COMPACT" \
    -Phistsfile="$OUTDIR/b0_bench.root" \
    "$OUTDIR/b0_sim.edm4hep.root"

echo ""
echo "=== Step 4: Validate residuals ==="
python3 "${SCRIPT_DIR}/validate.py" \
    --detector B0Tracker \
    --histfile "$OUTDIR/b0_bench.root" \
    --max-rms-x-um 150 \
    --max-rms-y-um 150 \
    --min-entries 50
