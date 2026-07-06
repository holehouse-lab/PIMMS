#!/usr/bin/env python
"""
Build the starting configuration for the multiphase (core-shell) condensate demo
and write it to a restart file.

Four chain species A, B, C, D (each a sticky 5-bead homopolymer) are packed into a
single spherical droplet as four concentric bands - A in the core, then B, C and D
as successive shells. The bands are laid down by (1) placing many non-overlapping
straight 5-bead rods uniformly inside a sphere, (2) sorting them by how far their
centre sits from the droplet centre, and (3) handing the innermost rods to A, the
next to B, and so on. Fewer chains go to the inner bands and more to the outer bands
(a shell's volume grows with radius), so the four layers come out at comparable
thickness.

The seed is deliberately ROUGH (the bands overlap and the droplet is not yet
spherical). The accompanying KEYFILE.kf runs the Monte Carlo that rounds the droplet
and sharpens the four interfaces - the interactions in params.prm (each species
attracts itself and its immediate neighbours but is REPELLED by non-neighbours) make
the concentric A|B|C|D nesting a free-energy minimum, so the layering is maintained
and cleaned up rather than scrambled.

Run from this directory:  ``python build_onion_seed.py``
Output:                   ``onion_seed.restart``
"""

import pickle

import numpy as np

# ---------------------------------------------------------------------------
# geometry / composition (everything fits inside a 40 x 40 x 40 hardwall box)
# ---------------------------------------------------------------------------
DIMS = [40, 40, 40]
CENTRE = np.array([20, 20, 20])
RADIUS = 13                       # radius of the seed droplet (lattice units)

# chains per species, innermost -> outermost. Inner bands get fewer chains because
# an inner shell has less volume; this gives four layers of similar thickness.
COUNTS = {"A": 90, "B": 150, "C": 180, "D": 220}
TYPE_ID = {"A": 0, "B": 1, "C": 2, "D": 3}     # chainType must be 1:1 with sequence

# unit steps for a straight 5-bead rod
AXES = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

rng = np.random.default_rng(1)


def _in_box(p):
    return all(0 <= p[d] < DIMS[d] for d in range(3))


def build():
    """Place non-overlapping rods in the sphere, sorted by radius for banding."""
    occupied = set()
    rods = []                       # (centre_radius, [positions]) per rod
    n_target = sum(COUNTS.values())

    while len(rods) < n_target:
        for _attempt in range(600):
            start = CENTRE + rng.integers(-RADIUS, RADIUS + 1, size=3)
            axis = np.array(AXES[rng.integers(6)])
            pts = [tuple(int(v) for v in (start + axis * k)) for k in range(5)]

            if len(set(pts)) < 5:
                continue
            if any(np.linalg.norm(np.array(p) - CENTRE) > RADIUS for p in pts):
                continue
            if not all(_in_box(p) for p in pts):
                continue
            if any(p in occupied for p in pts):
                continue

            for p in pts:
                occupied.add(p)
            centre_r = np.linalg.norm(np.mean(pts, axis=0) - CENTRE)
            rods.append((centre_r, [list(p) for p in pts]))
            break
        else:
            # could not place another rod (droplet is jammed) - stop early
            break

    # innermost rods first, so the concentric assignment below is radial
    rods.sort(key=lambda r: r[0])

    chains = {}
    cid = 1
    idx = 0
    for letter, count in COUNTS.items():
        for _ in range(min(count, len(rods) - idx)):
            chains[cid] = [rods[idx][1], letter * 5, TYPE_ID[letter]]
            cid += 1
            idx += 1

    return chains


if __name__ == "__main__":
    chains = build()
    restart = {"DIMENSIONS": DIMS, "ENERGY": 0, "HARDWALL": True, "CHAINS": chains}
    with open("onion_seed.restart", "wb") as fh:
        pickle.dump(restart, fh)

    n_beads = sum(len(c[0]) for c in chains.values())
    print(f"multiphase seed: {len(chains)} chains, {n_beads} beads "
          f"in 4 concentric bands (A core -> D shell)")
    print("wrote onion_seed.restart")
