"""
Frame - a view of one trajectory snapshot.

A Frame holds only ``(store, frame_index)``. Polymers are created on demand as
thin views; the grid and clusters are built lazily (and only for this frame) the
first time they are asked for - never eagerly at load time.
"""

import numpy as np

from pimms import lattice_analysis_utils as _lau
from .polymer import Polymer
from .cluster import Cluster


class _PosProvider:
    """Minimal chain adapter for PIMMS's ``get_cluster_distribution`` (which reads
    ``.chainID`` and ``.get_ordered_positions()``)."""
    __slots__ = ("chainID", "_p")

    def __init__(self, chainID, positions):
        self.chainID = chainID
        self._p = positions

    def get_ordered_positions(self):
        return self._p


class Frame:

    __slots__ = ("_store", "_f", "_clusters")

    def __init__(self, store, frame_index):
        self._store = store
        self._f = frame_index
        self._clusters = None

    # -- identity ----------------------------------------------------------
    @property
    def index(self):
        return self._f

    @property
    def time(self):
        return float(self._store.times[self._f])

    @property
    def dimensions(self):
        return self._store.dimensions

    @property
    def hardwall(self):
        return self._store.hardwall

    @property
    def n_chains(self):
        return self._store.n_chains

    @property
    def n_atoms(self):
        return self._store.n_atoms

    # -- positions / polymers ---------------------------------------------
    @property
    def positions(self):
        """Raw integer positions of every bead this frame, ``(n_atoms, 3)`` view."""
        return self._store.positions[self._f]

    @property
    def all_bead_positions(self):
        return self.positions

    def polymer(self, chain_index):
        return Polymer(self._store, self._f, chain_index)

    @property
    def polymers(self):
        return [Polymer(self._store, self._f, c) for c in range(self._store.n_chains)]

    def __len__(self):
        return self._store.n_chains

    def __getitem__(self, chain_index):
        n = self._store.n_chains
        i = int(chain_index)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError(f"chain index {chain_index} out of range (0..{n - 1})")
        return Polymer(self._store, self._f, i)

    def __iter__(self):
        for c in range(self._store.n_chains):
            yield Polymer(self._store, self._f, c)

    # -- grid / clusters ---------------------------------------------------
    @property
    def grid(self):
        """A ``dimensions``-shaped int grid for this frame (site = chain index + 1)."""
        return self._store.frame_grid(self._f)

    @property
    def clusters(self):
        """Connected-component clusters (list of :class:`Cluster`), largest first."""
        if self._clusters is None:
            store = self._store
            grid = store.frame_grid(self._f)
            off = store.topology.offsets
            frame = store.positions[self._f]
            chain_dict = {c + 1: _PosProvider(c + 1, frame[off[c]:off[c + 1]].tolist())
                          for c in range(store.n_chains)}
            cluster_lists = _lau.get_cluster_distribution(grid, chain_dict)
            self._clusters = [Cluster(store, self._f, [cid - 1 for cid in cl])
                              for cl in cluster_lists]
        return self._clusters

    @property
    def droplet(self):
        """The largest cluster in this frame (the condensate), or None if empty."""
        clusters = self.clusters
        return clusters[0] if clusters else None

    def __repr__(self):
        return f"<Frame index={self._f} n_chains={self._store.n_chains}>"
