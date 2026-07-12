#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# Top-level test runner for the LGAD charge-sharing EICrecon plugin suite.
#
# Phases:
#   1. Catch2 unit tests  (ctest on the build tree)
#   2. B0 tracker truth-residual benchmark
#
# Must be run inside eic-shell with the plugin already built & installed.
#
# Usage:
#   bash run_all_tests.sh [--nevents N] [--skip-unit] [--skip-bench]
set -e

NEVENTS=200
OUTDIR="/tmp/lgad_chargesharing_bench"
RUN_B0=true
RUN_UNIT=true
RUN_BENCH=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --nevents)   NEVENTS="$2"; shift 2 ;;
        --outdir)    OUTDIR="$2";  shift 2 ;;
        --skip-unit) RUN_UNIT=false; shift ;;
        --skip-bench) RUN_BENCH=false; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PLUGIN_DIR="${REPO_ROOT}/install"
BUILD_DIR="${REPO_ROOT}/build"
BENCH_DIR="${REPO_ROOT}/src/benchmarks/lgad_chargesharing"

echo "============================================"
echo "  LGAD Charge-Sharing Test Suite"
echo "============================================"
echo "  Events:     $NEVENTS"
echo "  Output:     $OUTDIR"
echo "  Unit tests: $RUN_UNIT"
echo "  Benchmarks: $RUN_BENCH"
echo ""

# ----------------------------------------------------------------------------
# Environment checks
# ----------------------------------------------------------------------------
echo "=== Checking environment ==="
ERRORS=0

if $RUN_BENCH; then
    if ! command -v ddsim &>/dev/null; then
        echo "  ERROR: ddsim not found. Are you inside eic-shell?"
        ERRORS=$((ERRORS + 1))
    fi

    if ! command -v eicrecon &>/dev/null; then
        echo "  ERROR: eicrecon not found. Are you inside eic-shell?"
        ERRORS=$((ERRORS + 1))
    fi

    if [ ! -d "$PLUGIN_DIR/plugins" ]; then
        echo "  ERROR: Plugins not installed at $PLUGIN_DIR/plugins"
        echo "         Build with: cmake --build build --target install"
        ERRORS=$((ERRORS + 1))
    fi

    if ! python3 -c "import uproot" 2>/dev/null; then
        echo "  WARNING: uproot not installed. Benchmark validation will fail."
        echo "           Install with: pip3 install uproot awkward numpy"
    fi
fi

if [ $ERRORS -gt 0 ]; then
    echo ""
    echo "  $ERRORS error(s) found. Fix before running tests."
    exit 1
fi
echo "  Environment OK"
echo ""

# Track phase exit codes
UNIT_EXIT=0
B0_EXIT=0

# ----------------------------------------------------------------------------
# Phase 1: Catch2 unit tests
# ----------------------------------------------------------------------------
if $RUN_UNIT; then
    echo "============================================"
    echo "  Phase 1: Catch2 unit tests"
    echo "============================================"
    if [ ! -d "$BUILD_DIR" ]; then
        echo "  BUILD_DIR $BUILD_DIR not found; skipping unit tests."
        echo "  Build first with: cmake -S . -B build -DBUILD_TESTING=ON"
        UNIT_EXIT=1
    else
        if ctest --test-dir "$BUILD_DIR" --output-on-failure; then
            UNIT_EXIT=0
        else
            UNIT_EXIT=$?
        fi
    fi
    echo ""
fi

# ----------------------------------------------------------------------------
# Phase 2: B0 benchmark
# ----------------------------------------------------------------------------
if $RUN_BENCH && $RUN_B0; then
    echo "============================================"
    echo "  Phase 2: B0 TRACKER BENCHMARK"
    echo "============================================"
    if bash "${BENCH_DIR}/run_b0.sh" "$NEVENTS" "$OUTDIR/b0"; then
        echo ""
        echo "  B0 benchmark: OK"
    else
        B0_EXIT=$?
        echo ""
        echo "  B0 benchmark: FAILED (exit $B0_EXIT)"
    fi
    echo ""
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
echo "============================================"
echo "  SUMMARY"
echo "============================================"
TOTAL_EXIT=0
status() { if [ "$1" -eq 0 ]; then echo "PASS"; else echo "FAIL"; fi; }
if $RUN_UNIT; then
    echo "  Unit tests:          $(status $UNIT_EXIT)"
    [ $UNIT_EXIT -ne 0 ] && TOTAL_EXIT=1
fi
if $RUN_BENCH && $RUN_B0; then
    echo "  B0 benchmark:        $(status $B0_EXIT)"
    [ $B0_EXIT -ne 0 ] && TOTAL_EXIT=1
fi
echo "============================================"
exit $TOTAL_EXIT
