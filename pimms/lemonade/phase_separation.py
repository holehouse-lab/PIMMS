"""
Phase-separation and droplet-physics analysis for lemonade trajectories.

This module quantifies liquid-liquid phase separation of a PIMMS lattice system:
the coexistence (binodal) densities of the dense and dilute phases, the condensed
fraction and cluster-size order parameters, interfacial width, and droplet shape.

Two complementary geometries are supported:

* **Droplet (spherical)** - a radial density profile about the largest cluster's
  centre of mass, fit to ``rho(r) = 1/2(rho_d + rho_v) - 1/2(rho_d - rho_v)
  tanh((r - R)/w)`` to extract the dense density ``rho_d``, the dilute (vapour)
  density ``rho_v``, the droplet radius ``R`` and the interface width ``w``.
* **Slab** - a 1D density profile along the box's long axis (the geometry of the
  ``slab_phase_separation`` demo), slabs re-centred per frame and fit to a
  two-interface ``tanh`` to extract the same coexistence quantities.

Densities are volume fractions (occupied lattice sites per available lattice site),
so ``rho`` runs 0..1 and is directly comparable across box sizes.

Typical use::

    from pimms.lemonade import phase_separation as ps
    result = ps.analyze(traj)          # everything, auto-detecting the geometry
    r, rho = ps.radial_density_profile(traj)
    fit = ps.fit_radial_profile(r, rho)
"""

from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# cluster / order-parameter helpers
# ---------------------------------------------------------------------------

def _largest_clusters(traj, min_beads=1):
    """Yield ``(frame_index, clusters_sorted_desc)`` skipping empty frames."""
    for f in range(traj.n_frames):
        clusters = [c for c in traj[f].clusters if c.n_beads >= min_beads]
        yield f, clusters


def condensed_fraction(traj, min_beads=1):
    """Per-frame fraction of all beads that sit in the single largest cluster.

    Returns ``(n_frames,)``. This is the basic phase-separation order parameter:
    ~0 in a well-mixed phase, ->1 when most material is in one condensate.
    """
    total = traj.n_atoms
    out = np.zeros(traj.n_frames)
    for f, clusters in _largest_clusters(traj, min_beads):
        if clusters:
            out[f] = clusters[0].n_beads / total
    return out


def largest_cluster_size(traj, by="beads", min_beads=1):
    """Per-frame size of the largest cluster (``by='beads'`` or ``'chains'``)."""
    out = np.zeros(traj.n_frames, dtype=np.int64)
    for f, clusters in _largest_clusters(traj, min_beads):
        if clusters:
            out[f] = clusters[0].n_beads if by == "beads" else clusters[0].n_chains
    return out


def number_of_clusters(traj, min_beads=2):
    """Per-frame count of clusters with at least ``min_beads`` beads."""
    out = np.zeros(traj.n_frames, dtype=np.int64)
    for f, clusters in _largest_clusters(traj, min_beads):
        out[f] = len(clusters)
    return out


def cluster_size_distribution(traj, by="beads", min_beads=1):
    """All cluster sizes pooled across frames, as one flat array (for histograms)."""
    sizes = []
    for _f, clusters in _largest_clusters(traj, min_beads):
        sizes.extend(c.n_beads if by == "beads" else c.n_chains for c in clusters)
    return np.asarray(sizes, dtype=np.int64)


# ---------------------------------------------------------------------------
# density profiles
# ---------------------------------------------------------------------------

def _min_image(delta, dims):
    """Wrap displacement(s) into ``[-L/2, L/2)`` per axis (broadcasts over dims)."""
    return delta - dims * np.round(delta / dims)


def _shell_site_counts(dimensions, edges):
    """Number of lattice sites whose PBC (min-image) distance from a point falls in
    each radial shell. Translationally invariant under PBC, so computed once."""
    dims = np.asarray(dimensions, dtype=np.float64)
    grids = np.indices(dimensions).reshape(len(dimensions), -1).T.astype(np.float64)
    r = np.sqrt((_min_image(grids, dims) ** 2).sum(axis=1))
    counts, _ = np.histogram(r, bins=edges)
    return counts.astype(np.float64)


def radial_density_profile(traj, bin_width=1.0, r_max=None, min_beads=2):
    """Spherically averaged density profile about the largest cluster's COM.

    For every frame the minimum-image distance of *every* bead from the condensate
    centre of mass is binned into radial shells and normalised by the number of
    lattice sites in each shell, giving a volume-fraction profile that falls from
    the dense core to the dilute background. Profiles are averaged over frames.

    Returns
    -------
    (radii, density) : tuple of 1D arrays
        Shell-centre radii and the frame-averaged occupied fraction ``rho(r)``.
    """
    dims = np.asarray(traj.dimensions, dtype=np.float64)
    nd = traj.n_dim
    if r_max is None:
        r_max = float(min(traj.dimensions)) / 2.0
    edges = np.arange(0.0, r_max + bin_width, bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    site_counts = _shell_site_counts(traj.dimensions, edges)
    safe = site_counts.copy()
    safe[safe == 0] = np.nan

    acc = np.zeros(len(centers))
    n_used = 0
    positions = traj.positions
    for f, clusters in _largest_clusters(traj, min_beads):
        if not clusters:
            continue
        # bin bead distances from the INTEGER COM using the same metric as the
        # (integer-origin, PBC-invariant) shell site counts, so occupied <= available
        # in every shell and the density is a true occupied fraction in [0, 1].
        com = np.mod(np.round(np.asarray(clusters[0].center_of_mass, dtype=np.float64)), dims[:nd])
        d = _min_image(positions[f][:, :nd].astype(np.float64) - com, dims[:nd])
        r = np.sqrt((d * d).sum(axis=1))
        counts, _ = np.histogram(r, bins=edges)
        acc += counts / safe
        n_used += 1

    density = acc / n_used if n_used else acc
    return centers, np.nan_to_num(density)


def slab_density_profile(traj, axis=None, min_beads=2):
    """1D density profile (volume fraction) along ``axis`` (default: the longest
    box axis), with the dense slab re-centred each frame so it does not smear out
    as the slab diffuses.

    Returns
    -------
    (coordinate, density) : tuple of 1D arrays
        Lattice coordinate along the axis and the frame-averaged occupied fraction.
    """
    dims = traj.dimensions
    if axis is None:
        axis = int(np.argmax(dims))
    length = dims[axis]
    cross_section = int(np.prod([d for i, d in enumerate(dims) if i != axis]))
    coord = np.arange(length)
    angle = 2.0 * np.pi * coord / length

    acc = np.zeros(length)
    positions = traj.positions
    for f in range(traj.n_frames):
        col = positions[f][:, axis]
        counts = np.bincount(col, minlength=length).astype(np.float64)
        # circular centre of mass of the 1D density -> shift dense region to L/2
        cx = (counts * np.cos(angle)).sum()
        cy = (counts * np.sin(angle)).sum()
        com = (np.arctan2(-cy, -cx) + np.pi) / (2.0 * np.pi) * length
        shift = int(round(length / 2.0 - com))
        acc += np.roll(counts, shift)
    return coord.astype(np.float64), acc / (traj.n_frames * cross_section)


# ---------------------------------------------------------------------------
# tanh binodal fits
# ---------------------------------------------------------------------------

def _tanh_droplet(r, rho_d, rho_v, radius, width):
    return 0.5 * (rho_d + rho_v) - 0.5 * (rho_d - rho_v) * np.tanh((r - radius) / width)


def _tanh_slab(z, rho_d, rho_v, half_width, width, center):
    return rho_v + 0.5 * (rho_d - rho_v) * (
        np.tanh((z - (center - half_width)) / width) - np.tanh((z - (center + half_width)) / width))


@dataclass
class BinodalFit:
    """Result of a ``tanh`` fit to a density profile."""
    rho_dense: float
    rho_dilute: float
    interface_width: float
    radius: float = float("nan")        # droplet radius (spherical geometry)
    half_width: float = float("nan")    # slab half-width (slab geometry)
    success: bool = True


def fit_radial_profile(radii, density):
    """Fit a spherical droplet profile; returns a :class:`BinodalFit`."""
    from scipy.optimize import curve_fit
    radii = np.asarray(radii, float)
    density = np.asarray(density, float)
    rho_d0 = float(density[:max(1, len(density) // 5)].mean())
    rho_v0 = float(density[-max(1, len(density) // 5):].mean())
    r0 = float(radii[np.argmin(np.abs(density - 0.5 * (rho_d0 + rho_v0)))])
    try:
        p, _ = curve_fit(_tanh_droplet, radii, density,
                         p0=[rho_d0, rho_v0, r0, 1.0],
                         bounds=([0, 0, 0, 0.1], [1, 1, radii.max(), radii.max()]),
                         maxfev=10000)
        return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                          radius=float(p[2]), interface_width=float(p[3]))
    except Exception:
        return BinodalFit(rho_dense=min(max(rho_d0, 0.0), 1.0),
                          rho_dilute=min(max(rho_v0, 0.0), 1.0),
                          interface_width=float("nan"), radius=r0, success=False)


def fit_slab_profile(coord, density):
    """Fit a slab (two-interface) profile; returns a :class:`BinodalFit`."""
    from scipy.optimize import curve_fit
    coord = np.asarray(coord, float)
    density = np.asarray(density, float)
    length = coord[-1] - coord[0] + 1
    center = float(length / 2.0)
    rho_d0 = float(density.max())
    rho_v0 = float(np.median(density[density < 0.5 * density.max()])) if np.any(density < 0.5 * density.max()) else 0.0
    hw0 = float(np.sum(density > 0.5 * (rho_d0 + rho_v0)) / 2.0)
    try:
        p, _ = curve_fit(lambda z, rd, rv, hw, w: _tanh_slab(z, rd, rv, hw, w, center),
                         coord, density, p0=[rho_d0, rho_v0, max(hw0, 1.0), 1.0],
                         bounds=([0, 0, 0, 0.1], [1, 1, length, length]), maxfev=10000)
        return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                          half_width=float(p[2]), interface_width=float(p[3]))
    except Exception:
        return BinodalFit(rho_dense=rho_d0, rho_dilute=rho_v0,
                          interface_width=float("nan"), half_width=hw0, success=False)


# ---------------------------------------------------------------------------
# droplet shape (frame-averaged over the largest cluster)
# ---------------------------------------------------------------------------

def droplet_shape(traj, min_beads=2):
    """Frame-averaged largest-cluster geometry.

    Returns a dict of mean radius of gyration, asphericity, sphericity, convex-hull
    volume and density (each averaged over the frames that contain a cluster).
    """
    rg, asph, sph, vol, dens = [], [], [], [], []
    for _f, clusters in _largest_clusters(traj, min_beads):
        if not clusters:
            continue
        c = clusters[0]
        rg.append(c.radius_of_gyration)
        asph.append(c.asphericity)
        sph.append(c.sphericity)
        vol.append(c.volume)
        dens.append(c.density)

    def _mean(x):
        x = np.asarray(x, float)
        x = x[np.isfinite(x) & (x > -1)]        # drop degenerate (-1) hull values
        return float(x.mean()) if x.size else float("nan")

    return {"radius_of_gyration": _mean(rg), "asphericity": _mean(asph),
            "sphericity": _mean(sph), "volume": _mean(vol), "density": _mean(dens)}


# ---------------------------------------------------------------------------
# top-level summary
# ---------------------------------------------------------------------------

@dataclass
class PhaseSeparationResult:
    geometry: str
    condensed_fraction: float
    condensed_fraction_series: np.ndarray = field(repr=False)
    n_clusters: float
    largest_cluster_beads: float
    binodal: BinodalFit
    shape: dict
    profile: tuple = field(repr=False, default=None)

    @property
    def rho_dense(self):
        return self.binodal.rho_dense

    @property
    def rho_dilute(self):
        return self.binodal.rho_dilute

    @property
    def is_phase_separated(self):
        """Heuristic: a clear density gap and most material condensed."""
        b = self.binodal
        return bool(np.isfinite(b.rho_dense) and np.isfinite(b.rho_dilute)
                    and b.rho_dense > 2.0 * max(b.rho_dilute, 1e-6)
                    and self.condensed_fraction > 0.3)


def analyze(traj, geometry="auto", min_beads=2):
    """Run the full phase-separation analysis and return a
    :class:`PhaseSeparationResult`.

    ``geometry`` is ``'sphere'``, ``'slab'`` or ``'auto'`` (slab if one box axis is
    noticeably longer than the others, else spherical).
    """
    dims = traj.dimensions
    if geometry == "auto":
        geometry = "slab" if max(dims) >= 1.5 * min(dims) else "sphere"

    cf = condensed_fraction(traj, min_beads=min_beads)
    n_clusters = number_of_clusters(traj, min_beads=min_beads)
    largest = largest_cluster_size(traj, by="beads", min_beads=min_beads)

    if geometry == "slab":
        coord, dens = slab_density_profile(traj, min_beads=min_beads)
        fit = fit_slab_profile(coord, dens)
        profile = (coord, dens)
    else:
        coord, dens = radial_density_profile(traj, min_beads=min_beads)
        fit = fit_radial_profile(coord, dens)
        profile = (coord, dens)

    return PhaseSeparationResult(
        geometry=geometry,
        condensed_fraction=float(cf.mean()),
        condensed_fraction_series=cf,
        n_clusters=float(n_clusters.mean()),
        largest_cluster_beads=float(largest.mean()),
        binodal=fit,
        shape=droplet_shape(traj, min_beads=min_beads),
        profile=profile,
    )
