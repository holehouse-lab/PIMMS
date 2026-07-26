"""
TrajectoryStore - the columnar backing store for a loaded trajectory.

All bead positions for the whole trajectory live in a single contiguous
``(n_frames, n_atoms, 3)`` int32 array; per-chain structure is described by the
CSR ``offsets`` in the :class:`~pimms.lemonade._topology.Topology`. Frame / Polymer
/ Cluster objects are thin *views* onto this store, so navigating the hierarchy
allocates no per-bead or per-frame Python objects.

Expensive quantities (unwrapped "whole" positions, and the batched Rg / COM /
gyration-tensor / end-to-end arrays) are computed once, for the whole trajectory
at once, and memoised.
"""

import numpy as np

from . import _analysis
from .kernels import _pbc


class _ChainPositions:
    """Minimal chain adapter for PIMMS's ``get_cluster_distribution``.

    That routine only ever reads ``.chainID`` and ``.get_ordered_positions()``.
    """

    __slots__ = ("chainID", "_positions")

    def __init__(self, chainID, positions):
        self.chainID = chainID
        self._positions = positions

    def get_ordered_positions(self):
        return self._positions


class TrajectoryStore:

    def __init__(self, positions, dimensions, spacing, hardwall, topology, times=None,
                 temperature=None):
        self.positions = np.ascontiguousarray(positions, dtype=np.int32)   # (nf, na, 3), wrapped
        self.dimensions = tuple(int(d) for d in dimensions)
        self.n_dim = len(self.dimensions)
        self.spacing = float(spacing)
        self.hardwall = bool(hardwall)
        self.temperature = None if temperature is None else float(temperature)
        self.topology = topology
        self.n_frames = int(self.positions.shape[0])
        self.times = (np.asarray(times, dtype=np.float64)
                      if times is not None else np.arange(self.n_frames, dtype=np.float64))

        # box vector padded to length 3 (z period is 1 for 2D so z never wraps)
        self._dimarr = np.array(list(self.dimensions) + [1] * (3 - self.n_dim), dtype=np.int64)

        # memoised batched results (computed on first access)
        self._whole = None
        self._com = None
        self._rg = None
        self._eig = None
        self._ete = None

        # per-frame connected-component membership, memoised (see cluster_membership)
        self._cluster_members = {}

    # -- sizes -------------------------------------------------------------
    @property
    def n_chains(self):
        return self.topology.n_chains

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    # -- positions ---------------------------------------------------------
    def whole_positions(self):
        """Unwrapped ("whole") positions, ``(nf, na, 3)`` int32; computed once."""
        if self._whole is None:
            self._whole = _pbc.unwrap_chains(self.positions, self.topology.offsets,
                                             self._dimarr, self.n_dim)
        return self._whole

    def _whole_k(self):
        return self.whole_positions()[..., :self.n_dim].astype(np.float64)

    # -- batched single-chain analyses ------------------------------------
    def centers_of_mass(self):
        if self._com is None:
            self._com = _analysis.centers_of_mass(self._whole_k(),
                                                  self.topology.offsets, self.topology.lengths)
        return self._com

    def radius_of_gyration(self):
        if self._rg is None:
            self._rg = _analysis.radius_of_gyration(self._whole_k(), self.topology.offsets,
                                                    self.topology.lengths, com=self.centers_of_mass())
        return self._rg

    def gyration_eigenvalues(self):
        if self._eig is None:
            self._eig = _analysis.gyration_eigenvalues(self._whole_k(), self.topology.offsets,
                                                       self.topology.lengths, com=self.centers_of_mass())
        return self._eig

    def asphericity(self):
        return _analysis.asphericity(self.gyration_eigenvalues())

    def end_to_end(self):
        if self._ete is None:
            self._ete = _analysis.end_to_end(self._whole_k(), self.topology.offsets, self.topology.lengths)
        return self._ete

    # -- grids -------------------------------------------------------------
    def frame_grid(self, f):
        """A freshly painted ``dimensions``-shaped int32 grid for frame ``f``
        (site value = chain index + 1, 0 = empty). Built on demand, never cached."""
        grid = np.zeros(self.dimensions, dtype=np.int32)
        fp = np.ascontiguousarray(self.positions[f], dtype=np.int32)
        ids = (self.topology.atom_chainid + 1).astype(np.int32)
        if self.n_dim == 3:
            _pbc.paint_frame_grid_3d(fp, ids, grid)
        else:
            _pbc.paint_frame_grid_2d(fp, ids, grid)
        return grid

    # -- clusters ----------------------------------------------------------
    def cluster_membership(self, f):
        """Connected components of frame ``f`` as chain-index lists, largest first.

        Memoised on the store rather than on a :class:`~pimms.lemonade.Frame`, because
        ``traj[f]`` mints a *new* Frame on every access - so a per-Frame cache is thrown
        away between passes. ``phase_separation.analyze`` walks the trajectory five
        times (condensed fraction, cluster count, largest cluster, the density profile
        and the droplet shape), and every one of those passes was re-running the whole
        connected-component search: measured at 5 decompositions per frame, about
        three-quarters of the total runtime of ``analyze()``.

        Only the cheap part - the membership lists - is cached. The per-cluster
        geometry (single-image positions, convex hulls) stays on the transient
        :class:`~pimms.lemonade.Cluster` objects so it can still be garbage collected;
        holding those would cost hundreds of MB over a long trajectory.

        Ordering is by bead count, descending (see :attr:`pimms.lemonade.Frame.clusters`).
        """
        members = self._cluster_members.get(f)
        if members is None:
            from pimms import lattice_analysis_utils as _lau

            offsets = self.topology.offsets
            frame = self.positions[f]
            chain_dict = {c + 1: _ChainPositions(c + 1, frame[offsets[c]:offsets[c + 1]].tolist())
                          for c in range(self.n_chains)}

            cluster_lists = _lau.get_cluster_distribution(self.frame_grid(f), chain_dict)

            # chainIDs are 1-based on the grid; lemonade chain indices are 0-based
            members = [[cid - 1 for cid in cluster] for cluster in cluster_lists]

            # bead count of a cluster straight from the CSR offsets - no need to build
            # Cluster objects just to sort them
            def _n_beads(cluster):
                return int(sum(offsets[c + 1] - offsets[c] for c in cluster))

            members.sort(key=_n_beads, reverse=True)
            self._cluster_members[f] = members

        return members

    # -- slicing -----------------------------------------------------------
    def subset(self, key):
        """A sub-store over a slice/array of frames; shares the topology."""
        return TrajectoryStore(self.positions[key], self.dimensions, self.spacing,
                               self.hardwall, self.topology, times=self.times[key],
                               temperature=self.temperature)
