#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
"""Generate forward particles (protons/pions) for B0 tracker testing.

Produces a HepMC3 ASCII file with particles in the B0 tracker acceptance.
The B0 tracker sits at z ~ 6.3 m along the ion beam axis, which is tilted
by the crossing angle (~25 mrad) from the z-axis in the x-z plane.

Usage (inside eic-shell):
    python3 gen_b0_particles.py --nevents 200 --output b0_gen.hepmc

Then feed to ddsim:
    ddsim --compactFile $COMPACT --inputFiles b0_gen.hepmc --outputFile b0_sim.edm4hep.root
"""

import argparse
import math
import numpy as np
import pyhepmc

# Physics constants
PROTON_MASS = 0.938272    # GeV
PION_MASS = 0.13957039    # GeV

# EIC crossing angle: ion beam tilted by ~25 mrad in x-z plane
ION_CROSSING_ANGLE = 0.025  # rad

# B0 tracker geometry
B0_Z = 6300.0       # mm, center of B0 tracker
B0_INNER_R = 35.0   # mm
B0_OUTER_R = 150.0  # mm

# Angular acceptance relative to ion beam axis
# At z=6300mm: inner_r=35mm -> theta_min ~ 5.6 mrad, outer_r=150mm -> theta_max ~ 23.8 mrad
THETA_MIN = 0.006   # rad (~6 mrad from ion beam axis)
THETA_MAX = 0.022   # rad (~22 mrad from ion beam axis)

PARTICLES = {
    "proton": (2212, PROTON_MASS),
    "piplus": (211, PION_MASS),
    "piminus": (-211, PION_MASS),
}


def generate_event(particle_pdg, particle_mass, Emin, Emax, rng):
    """Generate a single forward particle aimed at the B0 tracker.

    The particle direction is computed in the ion beam frame (tilted by
    crossing angle), then rotated to the lab frame.
    """
    # Energy: uniform in [Emin, Emax]
    E = rng.uniform(Emin, Emax)
    p = math.sqrt(max(0, E**2 - particle_mass**2))
    if p < 0.001:
        return None

    # Angle from ion beam axis: uniform in theta (flat in acceptance)
    theta_beam = rng.uniform(THETA_MIN, THETA_MAX)
    phi_beam = rng.uniform(0, 2 * math.pi)

    # Direction in ion beam frame
    px_beam = p * math.sin(theta_beam) * math.cos(phi_beam)
    py_beam = p * math.sin(theta_beam) * math.sin(phi_beam)
    pz_beam = p * math.cos(theta_beam)

    # Rotate from ion beam frame to lab frame:
    # Ion beam is tilted by +crossing_angle around y-axis
    # R_y(alpha): px_lab = px_beam*cos(a) + pz_beam*sin(a)
    #             pz_lab = -px_beam*sin(a) + pz_beam*cos(a)
    ca = math.cos(ION_CROSSING_ANGLE)
    sa = math.sin(ION_CROSSING_ANGLE)
    px_lab = px_beam * ca + pz_beam * sa
    py_lab = py_beam
    pz_lab = -px_beam * sa + pz_beam * ca

    return (px_lab, py_lab, pz_lab, E, particle_pdg)


def main():
    parser = argparse.ArgumentParser(description="Generate forward particles for B0 tracker testing")
    parser.add_argument("--nevents", type=int, default=200)
    parser.add_argument("--output", default="b0_gen.hepmc")
    parser.add_argument("--particle", default="proton", choices=list(PARTICLES.keys()),
                        help="Particle species")
    parser.add_argument("--Emin", type=float, default=50.0, help="Min energy (GeV)")
    parser.add_argument("--Emax", type=float, default=100.0, help="Max energy (GeV)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    pdg, mass = PARTICLES[args.particle]
    rng = np.random.default_rng(args.seed)

    with pyhepmc.open(args.output, "w") as f:
        n_written = 0
        for iev in range(args.nevents):
            result = generate_event(pdg, mass, args.Emin, args.Emax, rng)
            if result is None:
                continue

            px, py, pz, E, pid = result

            evt = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
            evt.event_number = n_written
            v = pyhepmc.GenVertex(pyhepmc.FourVector(0, 0, 0, 0))
            evt.add_vertex(v)
            p = pyhepmc.GenParticle(pyhepmc.FourVector(px, py, pz, E), pid, 1)
            p.generated_mass = mass
            v.add_particle_out(p)
            f.write(evt)

            n_written += 1

    print(f"Generated {n_written} {args.particle} events -> {args.output}")
    print(f"  Energy range: [{args.Emin}, {args.Emax}] GeV")
    print(f"  Theta range: [{THETA_MIN*1000:.1f}, {THETA_MAX*1000:.1f}] mrad from ion beam")
    print(f"  Crossing angle: {ION_CROSSING_ANGLE*1000:.1f} mrad")


if __name__ == "__main__":
    main()
