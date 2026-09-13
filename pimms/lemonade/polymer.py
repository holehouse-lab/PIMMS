## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Polymer - a lightweight view of one chain within one frame.

A Polymer holds only ``(store, frame_index, chain_index)``; every position array
and scalar it exposes is a view onto - or a cached lookup into - the trajectory
store's batched arrays, so creating one allocates nothing per bead.
"""

import numpy as np

from . import _analysis


class Polymer:
    """One chain in one frame: a view, not a copy.

    Obtained from a :class:`~pimms.lemonade.Frame` (``frame[chain_index]``) or by
    iterating a :class:`~pimms.lemonade.Cluster`. Its scalars are looked up in
    the store's batched arrays, so they cost nothing per chain.
    """

    __slots__ = ("_store", "_f", "_c", "_a0", "_a1")

    def __init__(self, store, frame_index, chain_index):
        """Bind a view to one chain of one frame and cache its atom range.

        Parameters
        ----------
        store : pimms.lemonade._store.TrajectoryStore
            The trajectory's backing store.
        frame_index : int
            Index of the frame this chain is viewed in.
        chain_index : int
            0-based chain index. Used to look the chain's atom block up in the
            topology's CSR offsets.
        """
        self._store = store
        self._f = frame_index
        self._c = chain_index
        self._a0 = int(store.topology.offsets[chain_index])
        self._a1 = int(store.topology.offsets[chain_index + 1])

    # -- identity ----------------------------------------------------------
    @property
    def chain_index(self):
        """int : 0-based index of this chain in the system."""
        return self._c

    @property
    def frame_index(self):
        """int : Index of the frame this chain is viewed in."""
        return self._f

    @property
    def chain_type(self):
        """int : Integer type label of this chain (the keyfile CHAIN order)."""
        return int(self._store.topology.chain_types[self._c])

    @property
    def sequence(self):
        """str : The chain's 1-letter bead sequence."""
        return self._store.topology.sequences[self._c]

    def __len__(self):
        """Number of beads in the chain.

        Returns
        -------
        int
            The chain length ``L``.
        """
        return self._a1 - self._a0

    # -- positions ---------------------------------------------------------
    @property
    def positions(self):
        """numpy.ndarray : ``(L, 3)`` int32 raw (wrapped) lattice positions.

        A read-only view into the store, in bonded order.
        """
        return self._store.positions[self._f, self._a0:self._a1]

    @property
    def whole_positions(self):
        """numpy.ndarray : ``(L, 3)`` int32 positions made contiguous across PBC.

        A read-only view into the store's unwrapped positions, so distances
        along the chain are plain Euclidean ones.
        """
        return self._store.whole_positions()[self._f, self._a0:self._a1]

    # -- cached single-chain scalars (indexed from the batched arrays) -----
    @property
    def center_of_mass(self):
        """numpy.ndarray : ``(n_dim,)`` float64 centre of mass, in whole coordinates."""
        return self._store.centers_of_mass()[self._f, self._c]

    @property
    def radius_of_gyration(self):
        """float : Radius of gyration of this chain, in lattice units."""
        return float(self._store.radius_of_gyration()[self._f, self._c])

    @property
    def asphericity(self):
        """float : Asphericity of this chain; zero for a spherically symmetric one."""
        return float(self._store.asphericity()[self._f, self._c])

    @property
    def end_to_end_distance(self):
        """float : Distance between the first and last bead, in lattice units."""
        return float(self._store.end_to_end()[self._f, self._c])

    @property
    def straddles_boundary(self):
        """bool : True if any bond crosses a periodic boundary in the raw positions."""
        nd = self._store.n_dim
        p = self.positions[:, :nd]
        if len(p) < 2:
            return False
        return bool(np.any(np.abs(np.diff(p, axis=0)) > 1))

    # -- per-chain matrices ------------------------------------------------
    def distance_map(self):
        """Inter-bead distance matrix of this chain, from its whole coordinates.

        Returns
        -------
        numpy.ndarray
            ``(L, L)`` float64 matrix of pairwise Euclidean distances, in
            lattice units.
        """
        return _analysis.distance_map(self.whole_positions[:, :self._store.n_dim].astype(np.float64))

    def internal_scaling(self):
        """Mean inter-bead distance as a function of sequence separation.

        Returns
        -------
        separations : numpy.ndarray
            ``(L - 1,)`` int64 sequence separations, running ``1 .. L-1``.
        mean_distance : numpy.ndarray
            ``(L - 1,)`` float64 mean distance at each separation, in lattice
            units.
        """
        return _analysis.internal_scaling(self.whole_positions[:, :self._store.n_dim].astype(np.float64))

    def __repr__(self):
        """One-line summary of the chain.

        Returns
        -------
        str
            The chain and frame indices, the sequence and the chain type.
        """
        return (f"<Polymer chain={self._c} frame={self._f} "
                f"seq={self.sequence!r} type={self.chain_type}>")
