## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Cluster - a connected group of polymers within one frame.

Built lazily by :attr:`Frame.clusters`. Geometric properties that need the cluster
gathered into a single periodic image (COM, Rg, volume, radial density) route
through PIMMS's single-image ("snakesearch", Cython-accelerated) and gross-property
machinery, computed once and cached.
"""

import numpy as np

from pimms import lattice_analysis_utils as _lau
from .polymer import Polymer


class Cluster:
    """A connected group of chains within one frame.

    Built lazily by :attr:`Frame.clusters <pimms.lemonade.Frame.clusters>`, which
    orders clusters by bead count, so ``frame.clusters[0]`` is the condensate.
    Geometry is computed on the cluster gathered into a single periodic image and
    cached on the object.
    """

    __slots__ = ("_store", "_f", "_chains", "_si", "_gross")

    def __init__(self, store, frame_index, chain_indices):
        """Bind a cluster to a set of chains in one frame.

        Parameters
        ----------
        store : pimms.lemonade._store.TrajectoryStore
            The trajectory's backing store.
        frame_index : int
            Index of the frame this cluster was found in.
        chain_indices : iterable of int
            0-based indices of the chains making up the connected component.
        """
        self._store = store
        self._f = frame_index
        self._chains = list(chain_indices)
        self._si = None
        self._gross = None

    # -- membership --------------------------------------------------------
    @property
    def chain_indices(self):
        """tuple of int : 0-based indices of the chains in this cluster."""
        return tuple(self._chains)

    @property
    def n_chains(self):
        """int : Number of chains in the cluster."""
        return len(self._chains)

    @property
    def n_beads(self):
        """int : Total number of beads across the cluster's chains."""
        off = self._store.topology.offsets
        return int(sum(off[c + 1] - off[c] for c in self._chains))

    @property
    def polymers(self):
        """list of Polymer : A view of each chain in the cluster."""
        return [Polymer(self._store, self._f, c) for c in self._chains]

    def __len__(self):
        """Number of chains in the cluster.

        Returns
        -------
        int
            The chain count, so ``len(cluster) == cluster.n_chains``.
        """
        return len(self._chains)

    def __iter__(self):
        """Iterate over the chains of the cluster.

        Yields
        ------
        Polymer
            A view of each chain, in the order the clustering reported them.
        """
        for c in self._chains:
            yield Polymer(self._store, self._f, c)

    # -- positions ---------------------------------------------------------
    def _raw(self):
        """Gather the raw positions of the cluster's beads, chain by chain.

        Returns
        -------
        numpy.ndarray
            ``(n_beads, 3)`` int32 wrapped positions, concatenated in cluster
            chain order.
        """
        off = self._store.topology.offsets
        frame = self._store.positions[self._f]
        return np.concatenate([frame[off[c]:off[c + 1]] for c in self._chains], axis=0)

    @property
    def positions(self):
        """numpy.ndarray : ``(n_beads, 3)`` int32 raw (wrapped) bead positions."""
        return self._raw()

    def single_image_positions(self):
        """Cluster gathered into one periodic image (cached).

        Under HARDWALL the box is not periodic, clusters cannot straddle a wall,
        and the raw positions already ARE the single image - the periodic gather
        must not run (it used to drag chains across the wall when the periodic
        connected-component search wrongly merged them).

        Returns
        -------
        numpy.ndarray
            ``(n_beads, n_dim)`` float64 positions in a single image, so the
            cluster is spatially contiguous and its COM and hull are meaningful.
        """
        if self._si is None:
            nd = self._store.n_dim
            raw = self._raw()[:, :nd]
            if self._store.hardwall:
                self._si = np.asarray(raw, dtype=np.float64)
            else:
                self._si = np.asarray(_lau.correct_cluster_positions_to_single_image(
                    [raw], list(self._store.dimensions))[0], dtype=np.float64)
        return self._si

    # -- geometry ----------------------------------------------------------
    @property
    def center_of_mass(self):
        """numpy.ndarray : ``(n_dim,)`` float64 COM of the single-image cluster."""
        return self.single_image_positions().mean(axis=0)

    @property
    def radius_of_gyration(self):
        """float : Radius of gyration of the whole cluster, in lattice units."""
        si = self.single_image_positions()
        d = si - si.mean(axis=0)
        return float(np.sqrt(np.mean(np.einsum("ij,ij->i", d, d))))

    @property
    def asphericity(self):
        """float : Asphericity of the cluster; zero for a spherically symmetric one."""
        si = self.single_image_positions()
        d = si - si.mean(axis=0)
        tensor = (d.T @ d) / len(d)
        ev = np.linalg.eigvalsh(tensor)
        if ev.shape[0] == 3:
            return float(ev[2] - 0.5 * (ev[0] + ev[1]))
        return float(ev[-1] - ev[0])

    def _gross_props(self):
        """Convex-hull properties of the single-image cluster (cached).

        Returns
        -------
        list of float
            ``[volume, surface_area, density]`` of the convex hull, all three
            ``-1`` when the hull is degenerate. In 2D the "volume" is the
            polygon area and the "surface area" its perimeter.
        """
        if self._gross is None:
            self._gross = _lau.compute_cluster_gross_properties([self.single_image_positions()])[0]
        return self._gross

    @property
    def volume(self):
        """float : Convex-hull volume in lattice units (the hull area in 2D).

        ``-1`` if the cluster is too small or degenerate to build a hull from.
        """
        return float(self._gross_props()[0])

    @property
    def surface_area(self):
        """float : Convex-hull surface area (the hull perimeter in 2D).

        ``-1`` if the cluster is too small or degenerate to build a hull from.
        """
        return float(self._gross_props()[1])

    @property
    def density(self):
        """float : Beads per unit convex-hull volume (``-1`` if degenerate)."""
        return float(self._gross_props()[2])

    @property
    def sphericity(self):
        """float : Isoperimetric sphericity of the convex hull, in ``(0, 1]``.

        1 is a perfect sphere/circle; ``nan`` if the hull volume/area is
        degenerate.

        3D: ``pi**(1/3) (6 V)**(2/3) / A``.  2D: ``4 pi A / P**2`` (``V`` is area,
        ``A`` is perimeter).
        """
        vol = self.volume
        area = self.surface_area
        if vol <= 0 or area <= 0:
            return float("nan")
        if self._store.n_dim == 3:
            return float(np.pi ** (1.0 / 3.0) * (6.0 * vol) ** (2.0 / 3.0) / area)
        return float(4.0 * np.pi * vol / (area * area))

    def radial_density_profile(self, minimum_cluster_size_in_beads=None):
        """Radial occupancy profile about the cluster COM (see PIMMS).

        Each shell reports the fraction of its lattice sites that are occupied,
        where a shell is the set of sites at a given Chebyshev distance from the
        (rounded) centre of mass.

        **The profile starts at shell 1**, so entry ``k`` is the shell at
        Chebyshev distance ``k + 1`` - plot it against ``1, 2, 3, ...``, not
        against ``0, 1, 2, ...``. Shell 0 is the single site the centre of mass
        rounds onto; its occupancy is one site out of one and carries no density
        information, so it is never emitted. This matches the
        ``CLUSTER_RADIAL_DENSITY_PROFILE.dat`` column convention that PIMMS
        itself writes (see the output-files documentation).

        Parameters
        ----------
        minimum_cluster_size_in_beads : int, optional
            Clusters with fewer beads than this are skipped and ``None`` is
            returned for them (default ``None``, i.e. no size filter).

        Returns
        -------
        list of float or None
            Occupied-site fraction of shells ``1, 2, 3, ...`` about the cluster
            centre of mass, or ``None`` if the cluster was skipped (the
            underlying routine returns nothing for it, so indexing ``[0]`` would
            raise).
        """
        profiles = _lau.compute_cluster_radial_density_profile(
            [self.single_image_positions()], list(self._store.dimensions),
            minimum_cluster_size_in_beads=minimum_cluster_size_in_beads,
            hardwall=bool(self._store.hardwall))
        return profiles[0] if profiles else None

    # -- composition -------------------------------------------------------
    @property
    def chain_type_composition(self):
        """dict : ``chain_type -> number of chains of that type`` in the cluster."""
        types = self._store.topology.chain_types
        out = {}
        for c in self._chains:
            t = int(types[c])
            out[t] = out.get(t, 0) + 1
        return out

    @property
    def bead_type_composition(self):
        """dict : ``bead type character -> bead count`` across the whole cluster."""
        seqs = self._store.topology.sequences
        out = {}
        for c in self._chains:
            for ch in seqs[c]:
                out[ch] = out.get(ch, 0) + 1
        return out

    def __repr__(self):
        """One-line summary of the cluster.

        Returns
        -------
        str
            The frame index and the cluster's chain and bead counts.
        """
        return f"<Cluster frame={self._f} n_chains={self.n_chains} n_beads={self.n_beads}>"
