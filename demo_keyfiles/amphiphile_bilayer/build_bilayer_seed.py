#!/usr/bin/env python
"""
Build the starting configuration for the amphiphile bilayer demo and write it to a
restart file.

Each molecule is a 5-bead amphiphile ``HHTTT`` - two hydrophilic HEAD beads (``H``)
and three hydrophobic TAIL beads (``T``). They are laid down as a flat lipid-style
BILAYER lying in the x-y plane in the middle of a z-elongated, periodic box:

    z (up)
     |      H H . . . . H H          <- upper leaflet heads (outer surface)
     |      T T . . . . T T
     |      T T ...core... T T       <- hydrophobic tail core (both leaflets meet)
     |      T T . . . . T T
     |      H H . . . . H H          <- lower leaflet heads (outer surface)
     +---------------------- x,y (membrane plane)

The two leaflets point their tails inward (so the tails form a hydrophobic core) and
their heads outward (toward the empty "solvent" above and below along z). One
amphiphile is placed per (x, y) column in each leaflet, so the membrane tiles the
periodic x-y plane.

The seed is a perfectly ordered membrane; the accompanying KEYFILE.kf runs the Monte
Carlo that lets it relax to a natural, rough-edged bilayer. The interactions in
params.prm (tail-tail attraction, hydrophilic heads, head/tail demixing) make the
bilayer a stable arrangement, so it stays a membrane rather than dissolving - a clean
demonstration of amphiphile self-organisation.

Run from this directory:  ``python build_bilayer_seed.py``
Output:                   ``bilayer_seed.restart``

Bead types (defined in params.prm): H = hydrophilic head, T = hydrophobic tail.
"""

import pickle

# ---------------------------------------------------------------------------
# geometry (everything fits inside a 16 x 16 x 44 periodic box)
# ---------------------------------------------------------------------------
DIMS = [16, 16, 44]

# z-columns occupied by an amphiphile in each leaflet, listed in bead order
# (H, H, T, T, T from head to tail-tip). The two tail tips (23 and 22) sit one
# lattice site apart so the leaflets' tails touch and form a single core.
UPPER_LEAFLET_Z = [27, 26, 25, 24, 23]   # head at z=27 (outer, top), tails down to 23
LOWER_LEAFLET_Z = [18, 19, 20, 21, 22]   # head at z=18 (outer, bottom), tails up to 22

SEQUENCE = "HHTTT"
CHAIN_TYPE = 0                            # single amphiphile species


def build():
    """Place one amphiphile per (x, y) column in each of the two leaflets."""
    chains = {}
    cid = 1
    for x in range(DIMS[0]):
        for y in range(DIMS[1]):
            for leaflet_z in (UPPER_LEAFLET_Z, LOWER_LEAFLET_Z):
                positions = [[x, y, z] for z in leaflet_z]
                chains[cid] = [positions, SEQUENCE, CHAIN_TYPE]
                cid += 1
    return chains


if __name__ == "__main__":
    chains = build()
    restart = {"DIMENSIONS": DIMS, "ENERGY": 0, "HARDWALL": False, "CHAINS": chains}
    with open("bilayer_seed.restart", "wb") as fh:
        pickle.dump(restart, fh)

    n_beads = sum(len(c[0]) for c in chains.values())
    print(f"amphiphile bilayer seed: {len(chains)} molecules ({n_beads} beads) "
          f"in a bilayer spanning the {DIMS[0]}x{DIMS[1]} periodic plane")
    print("wrote bilayer_seed.restart")
