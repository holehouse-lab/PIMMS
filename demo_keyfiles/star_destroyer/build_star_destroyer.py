#!/usr/bin/env python
"""
Build a 3D Imperial Star Destroyer out of PIMMS chains and write it to a restart
file (plus a matching freeze file).

The Star Destroyer is assembled from horizontal "rows" of beads - each row is a
single connected chain (Chebyshev-1 adjacent beads), stacked in x (length) and z
(height) to fill a solid wedge. The pieces:

    * a flat triangular DECK (the iconic top-down arrowhead silhouette),
    * a WEDGE underside that thickens toward the stern (the side-on dagger shape),
    * a raised command SUPERSTRUCTURE near the stern,
    * a BRIDGE block, a central command TOWER spike and two sensor domes.

Every chain is given a sequential integer chainID starting at 1, and ALL of them
are written into ``freeze.txt`` so the whole ship is held fixed during the
simulation while the EXTRA_CHAIN "spaceships" fly around it.

Run from this directory:  ``python build_star_destroyer.py``
Outputs:                  ``star_destroyer.restart`` and ``freeze.txt``

Bead types (defined in params.prm): H = hull, T = tower/superstructure.
"""

import pickle

# ---------------------------------------------------------------------------
# geometry parameters (everything fits inside a 60 x 60 x 60 hardwall box)
# ---------------------------------------------------------------------------
DIMS = [60, 60, 60]
X_BOW, X_STERN = 12, 46          # length axis: pointed bow .. wide stern
YC = 30                          # centre line (width axis)
Z_DECK = 30                      # height of the flat top deck
MAX_HALF_W = 16                  # half-width of the deck at the stern
KEEL_DEPTH = 6                   # how deep the wedge gets under the stern

occupied = set()                 # (x, y, z) sites already used (no overlaps!)
chains = {}                      # chainID -> [positions, sequence, chainType]
_seq_types = {}                  # sequence -> chainType (PIMMS requires a 1:1 map)
_next_id = 1


def _chain_type(seq):
    """A unique chainType per distinct sequence (PIMMS enforces type<->seq 1:1)."""
    return _seq_types.setdefault(seq, len(_seq_types))


def _add(positions, bead):
    """Register one chain from a list of (already connected) positions."""
    global _next_id
    if not positions:
        return
    seq = bead * len(positions)
    chains[_next_id] = [positions, seq, _chain_type(seq)]
    _next_id += 1


def _row(x, z, y_lo, y_hi, bead):
    """Add one chain: a connected run of beads along y at fixed (x, z)."""
    positions = []
    for y in range(y_lo, y_hi + 1):
        site = (x, y, z)
        if site in occupied:
            raise RuntimeError(f"overlap at {site} - geometry bug")
        occupied.add(site)
        positions.append([x, y, z])
    _add(positions, bead)


def _column(x, y, z_lo, z_hi, bead):
    """Add one chain: a connected vertical run of beads along z at fixed (x, y)."""
    positions = []
    for z in range(z_lo, z_hi + 1):
        site = (x, y, z)
        if site in occupied:
            raise RuntimeError(f"overlap at {site} - geometry bug")
        occupied.add(site)
        positions.append([x, y, z])
    _add(positions, bead)


def half_width(x):
    """Deck half-width at length-position x (0 at the bow, MAX at the stern)."""
    frac = (x - X_BOW) / (X_STERN - X_BOW)
    return round(frac * MAX_HALF_W)


# ---- 1. flat triangular deck (the top-down silhouette), bead type H ---------
for x in range(X_BOW, X_STERN + 1):
    hw = half_width(x)
    _row(x, Z_DECK, YC - hw, YC + hw, "H")

# ---- 2. wedge underside: thickens toward the stern, narrows to a keel --------
for x in range(X_BOW, X_STERN + 1):
    frac = (x - X_BOW) / (X_STERN - X_BOW)
    depth = round(frac * KEEL_DEPTH)
    hw = half_width(x)
    for d in range(1, depth + 1):
        hw_d = max(0, hw - 2 * d)          # each layer down is narrower -> V keel
        _row(x, Z_DECK - d, YC - hw_d, YC + hw_d, "H")

# ---- 3. command superstructure near the stern, bead type T ------------------
SX0, SX1 = 36, 44
for layer in range(1, 4):                  # z = 31, 32, 33
    z = Z_DECK + layer
    shw = max(1, 6 - layer)                # shrinks with height
    for x in range(SX0 + (layer - 1), SX1 - (layer - 1) + 1):
        _row(x, z, YC - shw, YC + shw, "T")

# ---- 4. bridge block (z = 34, 35) -------------------------------------------
for z in (Z_DECK + 4, Z_DECK + 5):
    for x in range(39, 42):
        _row(x, z, YC - 2, YC + 2, "T")

# ---- 5. central command tower spike (z = 36..38) ----------------------------
_column(40, YC, Z_DECK + 6, Z_DECK + 8, "T")

# ---- 6. the twin sensor domes flanking the tower ----------------------------
for dy in (-2, 2):
    _row(40, Z_DECK + 6, YC + dy, YC + dy, "T")

# ---------------------------------------------------------------------------
# write the restart file + the freeze file
# ---------------------------------------------------------------------------
restart = {
    "DIMENSIONS": DIMS,
    "ENERGY": 0,
    "HARDWALL": True,
    "CHAINS": chains,
}
with open("star_destroyer.restart", "wb") as fh:
    pickle.dump(restart, fh)

with open("freeze.txt", "w") as fh:
    fh.write("# freeze every chain that makes up the Star Destroyer hull\n")
    ids = list(chains.keys())
    for i in range(0, len(ids), 20):
        fh.write("C " + " ".join(str(c) for c in ids[i:i + 20]) + "\n")

n_beads = sum(len(c[0]) for c in chains.values())
print(f"Star Destroyer: {len(chains)} chains, {n_beads} beads "
      f"(chainIDs 1..{len(chains)} all frozen)")
print("wrote star_destroyer.restart and freeze.txt")
