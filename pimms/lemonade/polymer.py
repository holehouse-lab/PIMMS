"""
Polymer - a lightweight view of one chain within one frame.

A Polymer holds only ``(store, frame_index, chain_index)``; every position array
and scalar it exposes is a view onto - or a cached lookup into - the trajectory
store's batched arrays, so creating one allocates nothing per bead.
"""

import numpy as np

from . import _analysis


class Polymer:

    __slots__ = ("_store", "_f", "_c", "_a0", "_a1")

    def __init__(self, store, frame_index, chain_index):
        self._store = store
        self._f = frame_index
        self._c = chain_index
        self._a0 = int(store.topology.offsets[chain_index])
        self._a1 = int(store.topology.offsets[chain_index + 1])

    # -- identity ----------------------------------------------------------
    @property
    def chain_index(self):
        return self._c

    @property
    def frame_index(self):
        return self._f

    @property
    def chain_type(self):
        return int(self._store.topology.chain_types[self._c])

    @property
    def sequence(self):
        return self._store.topology.sequences[self._c]

    def __len__(self):
        return self._a1 - self._a0

    # -- positions ---------------------------------------------------------
    @property
    def positions(self):
        """Raw (wrapped) integer lattice positions, ``(L, 3)`` view."""
        return self._store.positions[self._f, self._a0:self._a1]

    @property
    def whole_positions(self):
        """Positions made contiguous across periodic boundaries, ``(L, 3)``."""
        return self._store.whole_positions()[self._f, self._a0:self._a1]

    # -- cached single-chain scalars (indexed from the batched arrays) -----
    @property
    def center_of_mass(self):
        return self._store.centers_of_mass()[self._f, self._c]

    @property
    def radius_of_gyration(self):
        return float(self._store.radius_of_gyration()[self._f, self._c])

    @property
    def asphericity(self):
        return float(self._store.asphericity()[self._f, self._c])

    @property
    def end_to_end_distance(self):
        return float(self._store.end_to_end()[self._f, self._c])

    @property
    def straddles_boundary(self):
        """True if any bond crosses a periodic boundary in the raw positions."""
        nd = self._store.n_dim
        p = self.positions[:, :nd]
        if len(p) < 2:
            return False
        return bool(np.any(np.abs(np.diff(p, axis=0)) > 1))

    # -- per-chain matrices ------------------------------------------------
    def distance_map(self):
        """Full ``(L, L)`` inter-bead Euclidean distance matrix (whole coords)."""
        return _analysis.distance_map(self.whole_positions[:, :self._store.n_dim].astype(np.float64))

    def internal_scaling(self):
        """``(separations, mean_distance)`` internal-scaling profile."""
        return _analysis.internal_scaling(self.whole_positions[:, :self._store.n_dim].astype(np.float64))

    def __repr__(self):
        return (f"<Polymer chain={self._c} frame={self._f} "
                f"seq={self.sequence!r} type={self.chain_type}>")
