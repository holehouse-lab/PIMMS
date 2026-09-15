## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
TrajectoryStore - the columnar backing store for a loaded trajectory.

All bead positions for the whole trajectory live in a single contiguous
``(n_frames, n_beads, 3)`` int32 array; per-chain structure is described by the
CSR ``offsets`` in the :class:`~pimms.lemonade._topology.Topology`. Frame / Polymer
/ Cluster objects are thin *views* onto this store, so navigating the hierarchy
allocates no per-bead or per-frame Python objects.

Expensive quantities (unwrapped "whole" positions, and the batched Rg / COM /
gyration-tensor / end-to-end arrays) are computed once, for the whole trajectory
at once, and memoised.
"""

import math
import numbers

import warnings

import numpy as np

from . import _analysis
from .kernels import _pbc


class _ChainPositions:
    """Minimal chain adapter for PIMMS's ``get_cluster_distribution``.

    That routine only ever reads ``.chainID`` and ``.get_ordered_positions()``.
    """

    __slots__ = ("chainID", "_positions")

    def __init__(self, chainID, positions):
        """Wrap one chain's positions in the interface PIMMS's clustering expects.

        Parameters
        ----------
        chainID : int
            1-based chain identifier, matching the value painted on the frame
            grid (lemonade's 0-based chain index plus one).
        positions : list of list of int
            The chain's bead positions in bonded order, as ``[x, y, z]`` lists.
        """
        self.chainID = chainID
        self._positions = positions

    def get_ordered_positions(self):
        """Return the chain's bead positions in bonded order.

        Returns
        -------
        list of list of int
            The ``[x, y, z]`` positions this adapter was built with.
        """
        return self._positions


class TrajectoryStore:
    """Columnar backing store holding every frame of one loaded trajectory.

    Frame, Polymer and Cluster objects are views onto this store; the batched
    per-chain analyses and the unwrapped positions are computed once here and
    memoised for the whole trajectory.
    """

    def __init__(self, positions, dimensions, spacing, hardwall, topology, times=None,
                 temperature=None):
        """Validate the trajectory arrays and set up the memoisation slots.

        The positions are copied into a contiguous int32 array and made
        read-only: geometry derived from them is cached, so a caller mutating
        the array afterwards would leave different properties describing
        different trajectories.

        Parameters
        ----------
        positions : array_like
            ``(n_frames, n_beads, 3)`` array of integer lattice coordinates,
            already wrapped into the box. The third column is zero for a 2D
            system.
        dimensions : sequence of int
            Box extent in lattice units: 2 or 3 positive integers.
        spacing : float
            Lattice spacing in angstroms (PIMMS ``LATTICE_TO_ANGSTROMS``).
        hardwall : bool
            ``True`` if the run used hard walls rather than periodic boundaries.
        topology : pimms.lemonade._topology.Topology
            Chain/bead topology; must describe exactly ``n_beads`` beads.
        times : array_like, optional
            One time value per frame (default ``None``, which numbers the frames
            ``0 .. n_frames-1``).
        temperature : float, optional
            Simulation temperature in PIMMS reduced units, used as ``k_B T``
            (default ``None``, meaning unknown).

        Raises
        ------
        ValueError
            If the positions are not an integer ``(frames, beads, 3)`` array, if
            they fall outside the box (or carry a non-zero z in 2D), if the
            topology's bead count disagrees with the positions, if ``spacing``
            or ``temperature`` is not finite and positive, if ``hardwall`` is not
            a bool, or if ``times`` does not hold one finite value per frame.
        """
        raw_positions = np.asarray(positions)
        if raw_positions.ndim != 3 or raw_positions.shape[2] != 3:
            raise ValueError("TrajectoryStore positions must have shape (frames, beads, 3)")
        if not np.issubdtype(raw_positions.dtype, np.integer):
            raise ValueError("TrajectoryStore positions must contain integer lattice coordinates")

        try:
            raw_dimensions = tuple(dimensions)
        except TypeError:
            raise ValueError(
                "TrajectoryStore dimensions must contain 2 or 3 positive integers")
        int32_max = np.iinfo(np.int32).max
        if (len(raw_dimensions) not in (2, 3) or
                any(isinstance(d, (bool, np.bool_)) or
                    not isinstance(d, numbers.Integral) or d <= 0 or d > int32_max
                    for d in raw_dimensions)):
            raise ValueError("TrajectoryStore dimensions must contain 2 or 3 positive integers")
        self.dimensions = tuple(int(d) for d in raw_dimensions)
        if raw_positions.shape[1] != topology.n_beads:
            raise ValueError(
                f"topology describes {topology.n_beads} beads but positions contain "
                f"{raw_positions.shape[1]}")

        for axis, extent in enumerate(self.dimensions):
            coord = raw_positions[..., axis]
            if coord.size and (coord.min() < 0 or coord.max() >= extent):
                raise ValueError(
                    f"TrajectoryStore positions are outside dimension {axis} (size {extent})")
        if len(self.dimensions) == 2 and raw_positions[..., 2].size and np.any(
                raw_positions[..., 2] != 0):
            raise ValueError("2D TrajectoryStore positions must have z == 0")

        try:
            spacing = float(spacing)
        except (TypeError, ValueError, OverflowError):
            raise ValueError("TrajectoryStore spacing must be a finite positive number")
        if not math.isfinite(spacing) or spacing <= 0:
            raise ValueError("TrajectoryStore spacing must be a finite positive number")
        if not isinstance(hardwall, (bool, np.bool_)):
            raise ValueError("TrajectoryStore hardwall must be True or False")
        if temperature is not None:
            try:
                temperature = float(temperature)
            except (TypeError, ValueError, OverflowError):
                raise ValueError(
                    "TrajectoryStore temperature must be a finite positive number")
            if not math.isfinite(temperature) or temperature <= 0:
                raise ValueError(
                    "TrajectoryStore temperature must be a finite positive number")

        positions = np.ascontiguousarray(raw_positions, dtype=np.int32)
        if positions is raw_positions:
            positions = positions.copy()          # never freeze the caller's own array
        self.positions = positions
        # TrajectoryStore memoises geometry derived from this array. Allowing a
        # caller to mutate positions after (say) whole_positions or Rg had been
        # cached would leave different properties describing different
        # trajectories, so loaded trajectory coordinates are immutable.
        self.positions.flags.writeable = False
        self.n_dim = len(self.dimensions)
        self.spacing = spacing
        self.hardwall = bool(hardwall)
        self.temperature = temperature
        self.topology = topology
        self.n_frames = int(self.positions.shape[0])
        self.times = (np.array(times, dtype=np.float64)
                      if times is not None else np.arange(self.n_frames, dtype=np.float64))
        if self.times.shape != (self.n_frames,):
            raise ValueError("TrajectoryStore times must contain one value per frame")
        if not np.all(np.isfinite(self.times)):
            raise ValueError("TrajectoryStore times must contain only finite values")
        self.times.flags.writeable = False

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
        """int : Number of chains in the trajectory."""
        return self.topology.n_chains

    @property
    def n_beads(self):
        """int : Total number of beads per frame."""
        return self.topology.n_beads

    @property
    def n_atoms(self):
        """int : Deprecated alias of :attr:`n_beads`.

        PIMMS is a coarse-grained model and its particles are beads; the name was
        inherited from the PDB/XTC vocabulary and is kept only so scripts written
        against earlier builds keep running.
        """
        warnings.warn("n_atoms is deprecated, use n_beads (PIMMS has beads, not atoms)",
                      DeprecationWarning, stacklevel=2)
        return self.n_beads

    # -- positions ---------------------------------------------------------
    def whole_positions(self):
        """Unwrapped ("whole") positions, computed once and memoised.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_beads, 3)`` int32 read-only array in which every
            chain is contiguous across the periodic boundaries, so intra-chain
            distances are ordinary Euclidean distances.
        """
        if self._whole is None:
            self._whole = _pbc.unwrap_chains(self.positions, self.topology.offsets,
                                             self._dimarr, self.n_dim)
            self._whole.flags.writeable = False
        return self._whole

    def _whole_k(self):
        """Whole positions as floats, truncated to the dimensions in use.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_beads, n_dim)`` float64 copy of the whole positions,
            which is what the batched routines in ``_analysis`` consume.
        """
        return self.whole_positions()[..., :self.n_dim].astype(np.float64)

    # -- batched single-chain analyses ------------------------------------
    def centers_of_mass(self):
        """Per-chain centres of mass for every frame (memoised).

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains, n_dim)`` float64 read-only array.
        """
        if self._com is None:
            self._com = _analysis.centers_of_mass(self._whole_k(),
                                                  self.topology.offsets, self.topology.lengths)
            self._com.flags.writeable = False
        return self._com

    def radius_of_gyration(self):
        """Per-chain radius of gyration for every frame (memoised).

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains)`` float64 read-only array, in lattice units.
        """
        if self._rg is None:
            self._rg = _analysis.radius_of_gyration(self._whole_k(), self.topology.offsets,
                                                    self.topology.lengths, com=self.centers_of_mass())
            self._rg.flags.writeable = False
        return self._rg

    def gyration_eigenvalues(self):
        """Ascending gyration-tensor eigenvalues per chain and frame (memoised).

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains, n_dim)`` float64 read-only array, eigenvalues
            ascending along the last axis.
        """
        if self._eig is None:
            self._eig = _analysis.gyration_eigenvalues(self._whole_k(), self.topology.offsets,
                                                       self.topology.lengths, com=self.centers_of_mass())
            self._eig.flags.writeable = False
        return self._eig

    def asphericity(self):
        """Per-chain asphericity for every frame, from the cached eigenvalues.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains)`` float64 array; zero for a spherically
            symmetric chain.
        """
        return _analysis.asphericity(self.gyration_eigenvalues())

    def end_to_end(self):
        """Per-chain end-to-end distance for every frame (memoised).

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains)`` float64 read-only array, in lattice units.
        """
        if self._ete is None:
            self._ete = _analysis.end_to_end(self._whole_k(), self.topology.offsets)
            self._ete.flags.writeable = False
        return self._ete

    # -- grids -------------------------------------------------------------
    def frame_grid(self, f):
        """A freshly painted occupancy grid for one frame.

        Built on demand, never cached.

        NOTE: the painting kernels index with bounds checking off and require
        in-box positions; ``lemonade.load`` guarantees this (``np.mod`` at load
        time), but a hand-built TrajectoryStore must wrap its positions itself.

        Parameters
        ----------
        f : int
            Frame index to paint.

        Returns
        -------
        numpy.ndarray
            ``dimensions``-shaped int32 grid; each occupied site holds the
            1-based chain index of the bead on it, empty sites hold ``0``.
        """
        grid = np.zeros(self.dimensions, dtype=np.int32)
        fp = np.ascontiguousarray(self.positions[f], dtype=np.int32)
        ids = (self.topology.bead_chainid + 1).astype(np.int32)
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

        Parameters
        ----------
        f : int
            Frame index to decompose.

        Returns
        -------
        list of list of int
            One list of 0-based chain indices per connected component, largest
            (most beads) first.
        """
        members = self._cluster_members.get(f)
        if members is None:
            from pimms import lattice_analysis_utils as _lau

            offsets = self.topology.offsets
            frame = self.positions[f]
            chain_dict = {c + 1: _ChainPositions(c + 1, frame[offsets[c]:offsets[c + 1]].tolist())
                          for c in range(self.n_chains)}

            # honour the box boundary: under HARDWALL two chains against opposite
            # walls are NOT neighbours, so the connected-component search must not
            # treat the box as periodic (it otherwise merged them into one cluster
            # and the single-image gather then dragged one across the wall).
            cluster_lists = _lau.get_cluster_distribution(self.frame_grid(f), chain_dict,
                                                          hardwall=self.hardwall)

            # chainIDs are 1-based on the grid; lemonade chain indices are 0-based
            members = [[cid - 1 for cid in cluster] for cluster in cluster_lists]

            # bead count of a cluster straight from the CSR offsets - no need to build
            # Cluster objects just to sort them
            def _n_beads(cluster):
                """Total bead count of one cluster, read from the CSR offsets.

                Parameters
                ----------
                cluster : list of int
                    0-based chain indices making up the cluster.

                Returns
                -------
                int
                    Number of beads across those chains.
                """
                return int(sum(offsets[c + 1] - offsets[c] for c in cluster))

            members.sort(key=_n_beads, reverse=True)
            self._cluster_members[f] = members

        return members

    # -- slicing -----------------------------------------------------------
    def subset(self, key):
        """A sub-store over a subset of frames; shares the topology.

        Parameters
        ----------
        key : slice or numpy.ndarray
            Anything that indexes the frame axis of the position array: a slice
            or an integer/boolean index array.

        Returns
        -------
        TrajectoryStore
            A new store holding the selected frames (and their times), with the
            same box, spacing, hardwall flag, topology and temperature.
        """
        return TrajectoryStore(self.positions[key], self.dimensions, self.spacing,
                               self.hardwall, self.topology, times=self.times[key],
                               temperature=self.temperature)
