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

    __slots__ = ("_store", "_f", "_chains", "_si", "_gross")

    def __init__(self, store, frame_index, chain_indices):
        self._store = store
        self._f = frame_index
        self._chains = list(chain_indices)
        self._si = None
        self._gross = None

    # -- membership --------------------------------------------------------
    @property
    def chain_indices(self):
        return tuple(self._chains)

    @property
    def n_chains(self):
        return len(self._chains)

    @property
    def n_beads(self):
        off = self._store.topology.offsets
        return int(sum(off[c + 1] - off[c] for c in self._chains))

    @property
    def polymers(self):
        return [Polymer(self._store, self._f, c) for c in self._chains]

    def __len__(self):
        return len(self._chains)

    def __iter__(self):
        for c in self._chains:
            yield Polymer(self._store, self._f, c)

    # -- positions ---------------------------------------------------------
    def _raw(self):
        off = self._store.topology.offsets
        frame = self._store.positions[self._f]
        return np.concatenate([frame[off[c]:off[c + 1]] for c in self._chains], axis=0)

    @property
    def positions(self):
        """Raw (wrapped) positions of every bead in the cluster, ``(n_beads, 3)``."""
        return self._raw()

    def single_image_positions(self):
        """Cluster gathered into one periodic image, ``(n_beads, n_dim)`` (cached)."""
        if self._si is None:
            nd = self._store.n_dim
            raw = self._raw()[:, :nd]
            self._si = np.asarray(_lau.correct_cluster_positions_to_single_image(
                [raw], list(self._store.dimensions))[0], dtype=np.float64)
        return self._si

    # -- geometry ----------------------------------------------------------
    @property
    def center_of_mass(self):
        return self.single_image_positions().mean(axis=0)

    @property
    def radius_of_gyration(self):
        si = self.single_image_positions()
        d = si - si.mean(axis=0)
        return float(np.sqrt(np.mean(np.einsum("ij,ij->i", d, d))))

    @property
    def asphericity(self):
        si = self.single_image_positions()
        d = si - si.mean(axis=0)
        tensor = (d.T @ d) / len(d)
        ev = np.linalg.eigvalsh(tensor)
        if ev.shape[0] == 3:
            return float(ev[2] - 0.5 * (ev[0] + ev[1]))
        return float(ev[-1] - ev[0])

    def _gross_props(self):
        if self._gross is None:
            self._gross = _lau.compute_cluster_gross_properties([self.single_image_positions()])[0]
        return self._gross

    @property
    def volume(self):
        """Convex-hull volume (``-1`` if the cluster is too small / degenerate)."""
        return float(self._gross_props()[0])

    @property
    def surface_area(self):
        return float(self._gross_props()[1])

    @property
    def density(self):
        """Beads per convex-hull volume (``-1`` if degenerate)."""
        return float(self._gross_props()[2])

    @property
    def sphericity(self):
        """Isoperimetric sphericity of the convex hull, in ``(0, 1]`` (1 = a perfect
        sphere/circle). ``nan`` if the hull volume/area is degenerate.

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
        """Radial occupancy profile about the cluster COM (see PIMMS)."""
        return _lau.compute_cluster_radial_density_profile(
            [self.single_image_positions()], list(self._store.dimensions),
            minimum_cluster_size_in_beads=minimum_cluster_size_in_beads)[0]

    # -- composition -------------------------------------------------------
    @property
    def chain_type_composition(self):
        """dict ``chain_type -> count``."""
        types = self._store.topology.chain_types
        out = {}
        for c in self._chains:
            t = int(types[c])
            out[t] = out.get(t, 0) + 1
        return out

    @property
    def bead_type_composition(self):
        """dict ``bead_type_char -> count`` across all beads in the cluster."""
        seqs = self._store.topology.sequences
        out = {}
        for c in self._chains:
            for ch in seqs[c]:
                out[ch] = out.get(ch, 0) + 1
        return out

    def __repr__(self):
        return f"<Cluster frame={self._f} n_chains={self.n_chains} n_beads={self.n_beads}>"
