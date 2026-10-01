## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Frame - a view of one trajectory snapshot.

A Frame holds only ``(store, frame_index)``. Polymers are created on demand as
thin views; the grid and clusters are built lazily (and only for this frame) the
first time they are asked for - never eagerly at load time.
"""

import warnings
from .polymer import Polymer
from .cluster import Cluster


class Frame:
    """One snapshot of a trajectory: a view, not a copy.

    Obtained by indexing a :class:`~pimms.lemonade.LatticeTrajectory` with an
    integer. Index a Frame with a chain index for a
    :class:`~pimms.lemonade.Polymer`.
    """

    __slots__ = ("_store", "_f", "_clusters")

    def __init__(self, store, frame_index):
        """Bind a view to one frame of a store.

        Parameters
        ----------
        store : pimms.lemonade._store.TrajectoryStore
            The trajectory's backing store.
        frame_index : int
            Index of the frame this view refers to. Already range-checked by
            the caller (``LatticeTrajectory.__getitem__``).
        """
        self._store = store
        self._f = frame_index
        self._clusters = None

    # -- identity ----------------------------------------------------------
    @property
    def index(self):
        """int : Index of this frame within its trajectory."""
        return self._f

    @property
    def time(self):
        """float : Time stamp of this frame from the XTC.

        PIMMS stamps each frame with its index in the file, not a Monte Carlo
        step (see :attr:`LatticeTrajectory.times
        <pimms.lemonade.LatticeTrajectory.times>`).
        """
        return float(self._store.times[self._f])

    @property
    def dimensions(self):
        """tuple of int : Box extent in lattice units, one entry per dimension."""
        return self._store.dimensions

    @property
    def hardwall(self):
        """bool : True if the run used hard walls rather than periodic boundaries."""
        return self._store.hardwall

    @property
    def n_chains(self):
        """int : Number of chains in the system."""
        return self._store.n_chains

    @property
    def n_beads(self):
        """int : Total number of beads in the frame."""
        return self._store.n_beads

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

    # -- positions / polymers ---------------------------------------------
    @property
    def positions(self):
        """numpy.ndarray : ``(n_beads, 3)`` int32 raw positions of every bead.

        A read-only view into the store, wrapped into the box.
        """
        return self._store.positions[self._f]

    @property
    def all_bead_positions(self):
        """numpy.ndarray : Alias of :attr:`positions`, ``(n_beads, 3)`` int32."""
        return self.positions

    def polymer(self, chain_index):
        """Build a view of one chain in this frame.

        Equivalent to ``frame[chain_index]``: the index is range-checked and
        negative values count from the end. (It used to hand a raw negative
        index to the chain table, which produced a Polymer with an empty bead
        range.)

        Parameters
        ----------
        chain_index : int
            0-based chain index; negative values count from the end.

        Returns
        -------
        Polymer
            A view of that chain in this frame.

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        return self[chain_index]

    @property
    def polymers(self):
        """list of Polymer : A view of every chain in this frame, in chain order."""
        return [Polymer(self._store, self._f, c) for c in range(self._store.n_chains)]

    def __len__(self):
        """Number of chains in the frame.

        Returns
        -------
        int
            The chain count, so ``len(frame) == frame.n_chains``.
        """
        return self._store.n_chains

    def __getitem__(self, chain_index):
        """Select one chain of this frame.

        Parameters
        ----------
        chain_index : int
            0-based chain index; negative values count from the end.

        Returns
        -------
        Polymer
            A view of that chain in this frame.

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        n = self._store.n_chains
        i = int(chain_index)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError(f"chain index {chain_index} out of range (0..{n - 1})")
        return Polymer(self._store, self._f, i)

    def __iter__(self):
        """Iterate over the chains of this frame.

        Yields
        ------
        Polymer
            A view of each chain, in chain order.
        """
        for c in range(self._store.n_chains):
            yield Polymer(self._store, self._f, c)

    # -- grid / clusters ---------------------------------------------------
    @property
    def grid(self):
        """numpy.ndarray : ``dimensions``-shaped int32 occupancy grid for this frame.

        Each occupied site holds the 1-based chain index of the bead on it,
        empty sites hold ``0``. Painted afresh on every access.
        """
        return self._store.frame_grid(self._f)

    @property
    def clusters(self):
        """list of Cluster : Connected-component clusters, largest first.

        Two chains are connected when any bead of one is within Chebyshev
        distance 1 of any bead of the other (the short-range contact shell),
        across periodic boundaries but never through a hardwall; long-range
        pairs do not connect chains here.

        "Largest" means **most beads**. PIMMS's own
        :func:`~pimms.lattice_analysis_utils.get_cluster_distribution` orders clusters by
        the number of *chains* they contain, which is the same thing only when every chain
        is the same length. In a multi-component system with chains of different lengths the
        two orderings differ, and everything downstream that treats ``clusters[0]`` as the
        condensate (``droplet``, ``condensed_fraction``, ``largest_cluster_size``, the
        density profiles, the droplet shape and the surface-tension estimators) would then
        be measuring the wrong cluster. So the bead-count ordering is imposed once, in
        :meth:`TrajectoryStore.cluster_membership <pimms.lemonade._store.TrajectoryStore>`,
        for every consumer.

        Clusters with the same bead count are ordered by their lowest chain index,
        lowest first. In a frame where two clusters tie for the largest,
        ``clusters[0]`` is therefore the one holding the lowest-numbered chain: a
        fixed convention that does not move under translation or reflection of the
        system, though it is not a physical distinction.

        The expensive connected-component search is memoised on the *store*, not here:
        ``traj[f]`` builds a fresh Frame each time, so a per-Frame cache would be
        discarded between analysis passes. The Cluster objects themselves stay
        per-Frame, so their geometry caches remain collectable.
        """
        if self._clusters is None:
            store = self._store
            self._clusters = [Cluster(store, self._f, members)
                              for members in store.cluster_membership(self._f)]
        return self._clusters

    @property
    def droplet(self):
        """Cluster or None : The largest cluster in this frame (the condensate).

        ``None`` when the frame holds no cluster at all.
        """
        clusters = self.clusters
        return clusters[0] if clusters else None

    def __repr__(self):
        """One-line summary of the frame.

        Returns
        -------
        str
            The frame index and the number of chains.
        """
        return f"<Frame index={self._f} n_chains={self._store.n_chains}>"
