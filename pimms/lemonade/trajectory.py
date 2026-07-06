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

    __slots__ = ("_store",)

    def __init__(self, store):
        self._store = store

    # -- metadata ----------------------------------------------------------
    @property
    def store(self):
        return self._store

    @property
    def dimensions(self):
        return self._store.dimensions

    @property
    def n_dim(self):
        return self._store.n_dim

    @property
    def hardwall(self):
        return self._store.hardwall

    @property
    def spacing(self):
        return self._store.spacing

    @property
    def temperature(self):
        """Simulation temperature (== k_B T in PIMMS reduced units), or None."""
        return self._store.temperature

    @property
    def n_frames(self):
        return self._store.n_frames

    @property
    def n_chains(self):
        return self._store.n_chains

    @property
    def n_atoms(self):
        return self._store.n_atoms

    @property
    def topology(self):
        return self._store.topology

    @property
    def times(self):
        return self._store.times

    @property
    def sequences(self):
        return list(self._store.topology.sequences)

    @property
    def chain_types(self):
        return self._store.topology.chain_types

    # -- navigation --------------------------------------------------------
    def __len__(self):
        return self._store.n_frames

    def __getitem__(self, key):
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
        for i in range(self._store.n_frames):
            yield Frame(self._store, i)

    # -- raw arrays --------------------------------------------------------
    @property
    def positions(self):
        """Raw integer positions, ``(n_frames, n_atoms, 3)``."""
        return self._store.positions

    def whole_positions(self):
        """Positions with every chain made whole across PBC, ``(n_frames, n_atoms, 3)``."""
        return self._store.whole_positions()

    # -- batched analyses (whole trajectory at once) -----------------------
    def radius_of_gyration(self):
        """``(n_frames, n_chains)`` per-chain Rg for every frame."""
        return self._store.radius_of_gyration()

    def center_of_mass(self):
        """``(n_frames, n_chains, n_dim)`` per-chain COM for every frame."""
        return self._store.centers_of_mass()

    def asphericity(self):
        """``(n_frames, n_chains)`` per-chain asphericity for every frame."""
        return self._store.asphericity()

    def end_to_end_distance(self):
        """``(n_frames, n_chains)`` per-chain end-to-end distance for every frame."""
        return self._store.end_to_end()

    def __repr__(self):
        return (f"<LatticeTrajectory {self._store.n_frames} frames, "
                f"{self._store.n_chains} chains, {self._store.n_atoms} beads, "
                f"box {self._store.dimensions}>")
