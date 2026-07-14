#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
"""Per-cluster charge-neighborhood event display.

For each selected cluster, draws the digitized pad charges around the cluster
in the sensor's local grid frame, outlines the pads used by the Gaussian fit,
and marks where the generated primary and its Geant4 secondaries deposited
energy.

Usage (inside eic-shell):
    plot_neighborhood.py --reco b0_reco.edm4eic.root --event 3 --output nbh.pdf
    plot_neighborhood.py --reco b0_reco.edm4eic.root --find-secondaries --output nbh.pdf
"""

import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

E_CHARGE = 1.602176634e-19  # C

HITS_NAME = "B0TrackerChargeSharingHits"
CLUSTERS_NAME = "B0TrackerClusterHits"
ASSOC_NAME = "B0TrackerChargeSharingHitAssociations"

# Readout id spec "system:8,layer:4,module:12,sensor:2,x:32:-16,z:-16":
# bits 0-31 identify the sensor, x/z are signed 16-bit pad indices.
SENSOR_MASK = (1 << 32) - 1
IX_START, IZ_START, INDEX_WIDTH = 32, 48, 16


def signed_field(cell_id, start, width):
    value = (cell_id >> start) & ((1 << width) - 1)
    return value - (1 << width) if value >= 1 << (width - 1) else value


def pad_indices(cell_id):
    return (signed_field(cell_id, IX_START, INDEX_WIDTH),
            signed_field(cell_id, IZ_START, INDEX_WIDTH))


def object_key(obj):
    oid = obj.getObjectID()
    return (oid.collectionID, oid.index)


def generated_ancestor(particle):
    for _ in range(1024):
        if particle.getGeneratorStatus() != 0 or particle.parents_size() == 0:
            break
        particle = particle.getParents(0)
    return particle


def sim_particle(sim_hit):
    getter = getattr(sim_hit, "getParticle", None) or getattr(sim_hit, "getMCParticle")
    return getter()


def collect_truth(assocs):
    """rawHit key -> list of (sim position, weight, is_secondary, ancestor key)."""
    truth = {}
    for assoc in assocs:
        sim = assoc.getSimHit()
        if not sim.isAvailable():
            continue
        particle = sim_particle(sim)
        pos = sim.getPosition()
        entry = (
            (pos.x, pos.y, pos.z),
            max(0.0, assoc.getWeight()),
            particle.isAvailable() and particle.getGeneratorStatus() == 0,
            object_key(generated_ancestor(particle)) if particle.isAvailable() else None,
        )
        truth.setdefault(object_key(assoc.getRawHit()), []).append(entry)
    return truth


class SensorFrame:
    """Affine map between (ix, iz) pad indices and global positions, fitted
    from the sensor's own digitized pads (exact for a planar regular grid)."""

    def __init__(self, pads):
        indices = np.array([[p["ix"], p["iz"], 1.0] for p in pads])
        positions = np.array([p["pos"] for p in pads])
        coeff, _, self.rank, _ = np.linalg.lstsq(indices, positions, rcond=None)
        self.axis_u, self.axis_v, self.origin = coeff[0], coeff[1], coeff[2]
        self.pitch_u = float(np.linalg.norm(self.axis_u))
        self.pitch_v = float(np.linalg.norm(self.axis_v))

    def usable(self):
        return self.rank == 3 and self.pitch_u > 0.0 and self.pitch_v > 0.0

    def local_mm(self, position):
        basis = np.column_stack([self.axis_u, self.axis_v])
        frac, *_ = np.linalg.lstsq(basis, np.asarray(position) - self.origin, rcond=None)
        return frac[0] * self.pitch_u, frac[1] * self.pitch_v


def draw_cluster(pdf, event_index, cluster_index, cluster, pads_by_sensor, truth, margin):
    member_keys = set()
    for hit in cluster.getHits():
        raw = hit.getRawHit()
        if raw.isAvailable():
            member_keys.add(object_key(raw))
    sensors = {p["sensor"] for pads in pads_by_sensor.values() for p in pads
               if p["key"] in member_keys}
    if not sensors:
        return False
    sensor = sensors.pop()  # clusters never span sensors (DD4hep neighbours)
    pads = pads_by_sensor[sensor]
    if len(pads) < 3:
        return False
    frame = SensorFrame(pads)
    if not frame.usable():
        return False

    members = [p for p in pads if p["key"] in member_keys]
    seed = max(members, key=lambda p: p["charge"])
    span_u = max(abs(p["ix"] - seed["ix"]) for p in members)
    span_v = max(abs(p["iz"] - seed["iz"]) for p in members)
    nx, ny = span_u + margin, span_v + margin
    window = [p for p in pads
              if abs(p["ix"] - seed["ix"]) <= nx and abs(p["iz"] - seed["iz"]) <= ny]

    pitch_u, pitch_v = frame.pitch_u, frame.pitch_v
    to_u = lambda ix: ix * pitch_u
    to_v = lambda iz: iz * pitch_v

    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=max(p["charge"] for p in window))
    cmap = matplotlib.colormaps["viridis"]

    fig, ax = plt.subplots(figsize=(7.6, 6.6), layout="constrained")
    for ix in range(seed["ix"] - nx, seed["ix"] + nx + 1):
        for iz in range(seed["iz"] - ny, seed["iz"] + ny + 1):
            ax.add_patch(Rectangle((to_u(ix) - 0.5 * pitch_u, to_v(iz) - 0.5 * pitch_v),
                                   pitch_u, pitch_v, facecolor="none",
                                   edgecolor="0.85", linewidth=0.6, zorder=1))
    for pad in window:
        color = cmap(norm(pad["charge"]))
        ax.add_patch(Rectangle(
            (to_u(pad["ix"]) - 0.5 * pitch_u, to_v(pad["iz"]) - 0.5 * pitch_v),
            pitch_u, pitch_v, facecolor=color, edgecolor="white",
            linewidth=1.0, zorder=2))
        luminance = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]
        ax.text(to_u(pad["ix"]), to_v(pad["iz"]) - 0.3 * pitch_v,
                f"{pad['charge'] * 1e15:.2f}", ha="center", va="center", fontsize=7,
                color="black" if luminance > 0.55 else "white", zorder=4)
    for pad in members:
        ax.add_patch(Rectangle(
            (to_u(pad["ix"]) - 0.5 * pitch_u, to_v(pad["iz"]) - 0.5 * pitch_v),
            pitch_u, pitch_v, facecolor="none", edgecolor="black",
            linewidth=2.0, zorder=3))

    # Truth deposits linked to the window's pads, deduplicated by position.
    seen, primaries, secondaries = set(), [], []
    ancestor_weights = {}
    for pad in window:
        for pos, weight, is_secondary, ancestor in truth.get(pad["key"], []):
            if pad["key"] in member_keys and ancestor is not None:
                ancestor_weights[ancestor] = ancestor_weights.get(ancestor, 0.0) + weight
            if pos in seen:
                continue
            seen.add(pos)
            (secondaries if is_secondary else primaries).append(frame.local_mm(pos))
    for points, marker, color, size in ((primaries, "*", "#d62728", 230),
                                        (secondaries, "X", "#ff7f0e", 95)):
        if points:
            ax.scatter([p[0] for p in points], [p[1] for p in points], marker=marker,
                       s=size, c=color, edgecolors="white", linewidths=1.0, zorder=5)

    ax.set_xlim(to_u(seed["ix"]) - (nx + 0.5) * pitch_u,
                to_u(seed["ix"]) + (nx + 0.5) * pitch_u)
    ax.set_ylim(to_v(seed["iz"]) - (ny + 0.5) * pitch_v,
                to_v(seed["iz"]) + (ny + 0.5) * pitch_v)
    ax.set_aspect("equal")
    ax.set_xlabel("sensor u [mm]")
    ax.set_ylabel("sensor v [mm]")
    purity = (max(ancestor_weights.values()) / sum(ancestor_weights.values())
              if ancestor_weights else float("nan"))
    z_mm = np.mean([p["pos"][2] for p in members])
    ax.set_title(f"Event {event_index}, cluster {cluster_index} "
                 f"(z = {z_mm:.0f} mm, purity = {purity:.2f})")

    fig.colorbar(plt.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(vmin=0.0, vmax=norm.vmax * 1e15), cmap=cmap),
        ax=ax, label="Charge [fC]")
    handles = [
        Line2D([], [], marker="s", ls="none", markerfacecolor="none",
               markeredgecolor="black", markeredgewidth=2, markersize=11,
               label="Pads in Gaussian fit"),
        Line2D([], [], marker="*", ls="none", color="#d62728", markeredgecolor="white",
               markersize=13, label="Primary deposit"),
        Line2D([], [], marker="X", ls="none", color="#ff7f0e", markeredgecolor="white",
               markersize=9, label="Secondary deposit"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.09),
              ncol=3, frameon=False, fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)
    return True


def event_has_secondary(frame):
    for assoc in frame.get(ASSOC_NAME):
        sim = assoc.getSimHit()
        if sim.isAvailable():
            particle = sim_particle(sim)
            if particle.isAvailable() and particle.getGeneratorStatus() == 0:
                return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reco", required=True, help="podio output of eicrecon")
    parser.add_argument("--event", type=int, default=0, help="event index to draw")
    parser.add_argument("--cluster", type=int, default=None,
                        help="only this cluster index (default: all in the event)")
    parser.add_argument("--find-secondaries", action="store_true",
                        help="draw the first event containing a secondary deposit")
    parser.add_argument("--margin", type=int, default=2,
                        help="pads of margin around the cluster (default 2)")
    parser.add_argument("--output", required=True, help="output PDF")
    args = parser.parse_args()

    from podio.reading import get_reader

    events = get_reader(args.reco).get("events")

    event_index = args.event
    if args.find_secondaries:
        event_index = next(
            (i for i, frame in enumerate(events) if event_has_secondary(frame)), None)
        if event_index is None:
            sys.exit("No event with a secondary deposit found")
        events = get_reader(args.reco).get("events")

    frame = events[event_index]
    clusters = frame.get(CLUSTERS_NAME)
    truth = collect_truth(frame.get(ASSOC_NAME))

    pads_by_sensor = {}
    for hit in frame.get(HITS_NAME):
        raw = hit.getRawHit()
        if not raw.isAvailable():
            continue
        cell_id = hit.getCellID()
        ix, iz = pad_indices(cell_id)
        pos = hit.getPosition()
        pads_by_sensor.setdefault(cell_id & SENSOR_MASK, []).append({
            "key": object_key(raw),
            "sensor": cell_id & SENSOR_MASK,
            "ix": ix, "iz": iz,
            "pos": (pos.x, pos.y, pos.z),
            "charge": raw.getCharge() * E_CHARGE,
        })

    drawn = 0
    with PdfPages(args.output) as pdf:
        for index, cluster in enumerate(clusters):
            if args.cluster is not None and index != args.cluster:
                continue
            if draw_cluster(pdf, event_index, index, cluster, pads_by_sensor,
                            truth, args.margin):
                drawn += 1
    if drawn == 0:
        sys.exit(f"No clusters drawn for event {event_index}")
    print(f"Wrote {drawn} cluster page(s) for event {event_index} to {args.output}")


if __name__ == "__main__":
    main()
