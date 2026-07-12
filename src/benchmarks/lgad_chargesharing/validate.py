#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
"""Truth-residual validator for the LGAD charge-sharing benchmark.

Reads the hist file produced by -Phistsfile=... (populated by the
LGAD_chargesharing_benchmark plugin) and enforces residual thresholds on the
(residualX, residualY) distributions.

The plugin writes one TTree, ``LGADChargeSharing/hits``, containing the scalar
branches ``residualX`` and ``residualY``.

Exit codes:
  0 - all checks pass
  1 - at least one check failed
  2 - file missing or unreadable
"""

import argparse
import math
import sys
from pathlib import Path


def _import_uproot():
    try:
        import uproot  # noqa: F401

        return uproot
    except ImportError as exc:  # pragma: no cover
        print(f"ERROR: uproot is required for benchmark validation ({exc})")
        print("       Install with: pip install uproot awkward numpy")
        sys.exit(2)


def _load_residuals(uproot, path: Path):
    """Return the residualX and residualY numpy arrays."""
    f = uproot.open(str(path))
    try:
        tree = f["LGADChargeSharing/hits"]
    except KeyError:
        return None

    branches = set(tree.keys())
    required = {"residualX", "residualY"}
    if not required.issubset(branches):
        return None

    arrays = tree.arrays(["residualX", "residualY"], library="np")
    return arrays["residualX"], arrays["residualY"]


# Preserve the benchmark's established core resolution metric. Hits outside
# this window are reported as tails and excluded from the RMS.
CORE_WINDOW_MM = 0.5


def _robust_rms_mm(values):
    """Return (mean, rms, n_core, n_tail) in mm within +-CORE_WINDOW_MM."""
    import numpy as np

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan"), float("nan"), 0, 0
    core = arr[np.abs(arr) <= CORE_WINDOW_MM]
    n_tail = len(arr) - len(core)
    if len(core) == 0:
        return float("nan"), float("nan"), 0, n_tail
    return float(core.mean()), float(core.std(ddof=0)), len(core), n_tail


def validate(path: Path, detector: str, max_rms_x_um: float, max_rms_y_um: float,
             min_entries: int) -> int:
    uproot = _import_uproot()

    if not path.exists():
        print(f"  [FAIL] hist file missing: {path}")
        return 2

    header = f"[{detector}] {path.name}"
    print("=" * 72)
    print(header)
    print("=" * 72)

    tree_data = _load_residuals(uproot, path)
    failures = 0

    if tree_data is None:
        print("  [FAIL] LGADChargeSharing/hits with residualX/residualY not found")
        return 1

    rx, ry = tree_data
    n = min(len(rx), len(ry))
    meanX, rmsX_mm, _, n_tailX = _robust_rms_mm(rx)
    meanY, rmsY_mm, _, n_tailY = _robust_rms_mm(ry)
    rmsX_um = rmsX_mm * 1000.0
    rmsY_um = rmsY_mm * 1000.0
    n_tail = max(n_tailX, n_tailY)

    print(f"  entries:                 {n}")
    print(f"  core window (mm):        +-{CORE_WINDOW_MM} "
          f"({n_tail} tail hit(s) outside, excluded from RMS)")
    print(f"  residualX (um) mean/rms: {meanX*1000:+.2f} / {rmsX_um:.2f}")
    print(f"  residualY (um) mean/rms: {meanY*1000:+.2f} / {rmsY_um:.2f}")
    print(f"  thresholds (um):         <= {max_rms_x_um} (X), <= {max_rms_y_um} (Y)")
    print(f"  minimum entries:         >= {min_entries}")

    def _check(label: str, condition: bool, detail: str) -> None:
        nonlocal failures
        status = " OK " if condition else "FAIL"
        print(f"  [{status}] {label}: {detail}")
        if not condition:
            failures += 1

    _check("entries", n >= min_entries, f"n={n}")
    _check("rms(residualX)", math.isfinite(rmsX_um) and rmsX_um <= max_rms_x_um,
           f"rmsX = {rmsX_um:.2f} um")
    _check("rms(residualY)", math.isfinite(rmsY_um) and rmsY_um <= max_rms_y_um,
           f"rmsY = {rmsY_um:.2f} um")

    print()
    if failures:
        print(f"  SUMMARY: {failures} failure(s) for {detector}")
        return 1
    print(f"  SUMMARY: all checks passed for {detector}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--histfile", required=True, help="Path to -Phistsfile output ROOT file")
    parser.add_argument("--detector", required=True,
                        choices=["B0Tracker"],
                        help="Detector key written by LGADChargeSharingMonitor")
    parser.add_argument("--max-rms-x-um", type=float, default=150.0,
                        help="Maximum allowed RMS of residualX (microns)")
    parser.add_argument("--max-rms-y-um", type=float, default=150.0,
                        help="Maximum allowed RMS of residualY (microns)")
    parser.add_argument("--min-entries", type=int, default=50,
                        help="Minimum number of hits needed for the check")
    args = parser.parse_args()

    return validate(Path(args.histfile), args.detector, args.max_rms_x_um,
                    args.max_rms_y_um, args.min_entries)


if __name__ == "__main__":
    sys.exit(main())
