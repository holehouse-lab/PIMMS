## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
LatticeTrajectory - the top of the lemonade object hierarchy.

    LatticeTrajectory  ->  Frame  ->  Polymer            (one chain in one frame)
                                  ->  Cluster -> Polymer  (connected group)

Indexing a trajectory by an integer yields a :class:`Frame` view; slicing yields a
new LatticeTrajectory over that frame range (sharing the topology, no re-parse).
Trajectory-level analysis methods return whole ``(n_frames, n_chains)`` arrays in a
handful of vectorised numpy operations.
"""

import numpy as np

from .frame import Frame


class LatticeTrajectory:
    """A whole PIMMS lattice trajectory: the top of the lemonade hierarchy.

    Index by an integer for a :class:`~pimms.lemonade.Frame`, or slice for a new
    LatticeTrajectory over that frame range. Built by
    :func:`pimms.lemonade.load`.
    """

    __slots__ = ("_store",)

    def __init__(self, store):
        """Wrap a backing store.

        Parameters
        ----------
        store : pimms.lemonade._store.TrajectoryStore
            The columnar store holding the positions, box and topology. Sliced
            trajectories share the topology of the store they came from.
        """
        self._store = store

    # -- metadata ----------------------------------------------------------
    @property
    def store(self):
        """pimms.lemonade._store.TrajectoryStore : The backing store."""
        return self._store

    @property
    def dimensions(self):
        """tuple of int : Box extent in lattice units, one entry per dimension."""
        return self._store.dimensions

    @property
    def n_dim(self):
        """int : Dimensionality of the system, 2 or 3."""
        return self._store.n_dim

    @property
    def hardwall(self):
        """bool : True if the run used hard walls rather than periodic boundaries."""
        return self._store.hardwall

    @property
    def spacing(self):
        """float : Lattice spacing in angstroms."""
        return self._store.spacing

    @property
    def temperature(self):
        """float or None : Simulation temperature (== k_B T in PIMMS reduced units)."""
        return self._store.temperature

    @property
    def n_frames(self):
        """int : Number of frames in the trajectory."""
        return self._store.n_frames

    @property
    def n_chains(self):
        """int : Number of chains in the system."""
        return self._store.n_chains

    @property
    def n_atoms(self):
        """int : Total number of beads per frame."""
        return self._store.n_atoms

    @property
    def topology(self):
        """pimms.lemonade._topology.Topology : The chain/bead topology."""
        return self._store.topology

    @property
    def times(self):
        """numpy.ndarray : ``(n_frames,)`` float64 simulation time of each frame."""
        return self._store.times

    @property
    def sequences(self):
        """list of str : The 1-letter bead sequence of each chain."""
        return list(self._store.topology.sequences)

    @property
    def chain_types(self):
        """numpy.ndarray : ``(n_chains,)`` int32 type label of each chain."""
        return self._store.topology.chain_types

    # -- navigation --------------------------------------------------------
    def __len__(self):
        """Number of frames.

        Returns
        -------
        int
            The frame count, so ``len(traj) == traj.n_frames``.
        """
        return self._store.n_frames

    def __getitem__(self, key):
        """Select a frame, or a sub-trajectory.

        Parameters
        ----------
        key : int or slice or numpy.ndarray
            An integer (negative counts from the end) picks one frame; a slice
            or index array picks a range of frames.

        Returns
        -------
        Frame or LatticeTrajectory
            A :class:`~pimms.lemonade.Frame` view for an integer key, otherwise
            a new LatticeTrajectory over the selected frames, sharing the
            topology.

        Raises
        ------
        IndexError
            If an integer key is out of range.
        """
        if isinstance(key, (int, np.integer)):
            i = int(key)
            if i < 0:
                i += self._store.n_frames
            if not 0 <= i < self._store.n_frames:
                raise IndexError(f"frame index {key} out of range (0..{self._store.n_frames - 1})")
            return Frame(self._store, i)
        # slice or fancy index -> a real sub-trajectory sharing the topology
        return LatticeTrajectory(self._store.subset(key))

    def __iter__(self):
        """Iterate over the frames in order.

        Yields
        ------
        Frame
            A view of each frame, oldest first.
        """
        for i in range(self._store.n_frames):
            yield Frame(self._store, i)

    # -- raw arrays --------------------------------------------------------
    @property
    def positions(self):
        """numpy.ndarray : ``(n_frames, n_atoms, 3)`` int32 raw lattice positions.

        Read-only, wrapped into the box. The z column is zero for a 2D system.
        """
        return self._store.positions

    def whole_positions(self):
        """Positions with every chain made whole across PBC.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_atoms, 3)`` int32 read-only array in which each chain
            is contiguous, so intra-chain distances are plain Euclidean ones.
        """
        return self._store.whole_positions()

    # -- batched analyses (whole trajectory at once) -----------------------
    def radius_of_gyration(self):
        """Per-chain radius of gyration for every frame.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains)`` float64 Rg in lattice units.
        """
        return self._store.radius_of_gyration()

    def center_of_mass(self):
        """Per-chain centre of mass for every frame.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains, n_dim)`` float64 centres of mass, in lattice
            units and in whole (unwrapped) coordinates.
        """
        return self._store.centers_of_mass()

    def asphericity(self):
        """Per-chain asphericity for every frame.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains)`` float64 asphericity; zero for a spherically
            symmetric chain.
        """
        return self._store.asphericity()

    def end_to_end_distance(self):
        """Per-chain end-to-end distance for every frame.

        Returns
        -------
        numpy.ndarray
            ``(n_frames, n_chains)`` float64 distance between the first and last
            bead of each chain, in lattice units.
        """
        return self._store.end_to_end()

    def __repr__(self):
        """One-line summary of the trajectory's size and box.

        Returns
        -------
        str
            The frame, chain and bead counts plus the box dimensions.
        """
        return (f"<LatticeTrajectory {self._store.n_frames} frames, "
                f"{self._store.n_chains} chains, {self._store.n_atoms} beads, "
                f"box {self._store.dimensions}>")
