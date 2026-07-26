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
    """Yield ``(frame_index, clusters)`` per frame, largest (most beads) first.

    :attr:`Frame.clusters <pimms.lemonade.Frame>` is ordered by bead count, so ``clusters[0]``
    here is the condensate. Frames with no qualifying cluster yield an empty list.
    """
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
    each radial shell. Translationally invariant under PBC, so computed once.

    Built from separable per-axis offsets and accumulated a slab at a time. Listing
    every site explicitly (``np.indices(dims).reshape(nd, -1).T`` as float64) allocates
    ``n_dim`` float64 values per lattice site up front - about 200 MB for a 200^3 box,
    plus the distance array - which scales with box volume and has nothing to do with
    how many beads are actually being analysed. Histogram counts are additive, so
    accumulating over slabs gives exactly the same result with memory bounded by one
    slab.
    """
    dims = np.asarray(dimensions, dtype=np.float64)

    # per-axis minimum-image offset of every coordinate from the origin
    offsets = [_min_image(np.arange(int(n), dtype=np.float64), d)
               for n, d in zip(dimensions, dims)]
    squared = [off ** 2 for off in offsets]

    counts = np.zeros(len(edges) - 1, dtype=np.int64)

    # sum the per-axis squared offsets by broadcasting, one x-slab at a time. The
    # addition order matches the axis-order reduction of the original.
    if len(dimensions) == 2:
        for x2 in squared[0]:
            r = np.sqrt(x2 + squared[1])
            counts += np.histogram(r, bins=edges)[0]
    else:
        plane = squared[1][:, np.newaxis] + squared[2][np.newaxis, :]
        for x2 in squared[0]:
            r = np.sqrt(x2 + plane)
            counts += np.histogram(r, bins=edges)[0]

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


def slab_density_profile(traj, axis=None):
    """1D density profile (volume fraction) along ``axis`` (default: the longest
    box axis), with the dense slab re-centred each frame so it does not smear out
    as the slab diffuses.

    Every bead in the box is binned, dense and dilute alike - that is what makes the
    result a *density* profile that a coexistence fit can be run against. (It therefore
    takes no ``min_beads``: unlike the radial profile, this one does not need clusters
    at all. It previously accepted a ``min_beads`` argument that was never read.)

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
    """Result of a ``tanh`` fit to a density profile.

    ``success`` means **the fit is well posed** - the data actually constrains the
    parameters - not merely that the optimiser converged. A converged-but-degenerate
    fit (see :func:`_fit_is_usable`) is reported with ``success=False``, because
    trusting it silently is the more dangerous failure. When ``success`` is ``False``,
    ``rho_dense`` and ``rho_dilute`` fall back to robust percentiles of the *observed*
    profile, so they stay bounded and physically meaningful (for a homogeneous system
    they simply coincide) instead of being unconstrained extrapolations.

    ``reason`` is empty on success and otherwise names the check that failed.

    .. important::

       ``success`` is a **numerical** guarantee, not a physical one. It says the fit is
       meaningful *as a fit*; it does **not** decide whether the system is phase
       separated. A well-posed fit can still describe a shallow density modulation in a
       one-phase system - and near the critical point the coexistence gap closes
       continuously, so no numerical check can draw that line for you.

       For the physical question, use :attr:`PhaseSeparationResult.is_phase_separated`,
       or apply your own density-contrast threshold to ``rho_dense`` / ``rho_dilute``.
       Be aware too that :func:`slab_density_profile` re-centres the slab every frame,
       which aligns the fluctuations of even a *homogeneous* system into a shallow
       central hump - so a small, well-fitted density gap is not by itself evidence of
       coexistence.
    """
    rho_dense: float
    rho_dilute: float
    interface_width: float
    radius: float = float("nan")        # droplet radius (spherical geometry)
    half_width: float = float("nan")    # slab half-width (slab geometry)
    success: bool = True
    reason: str = ""


#: Fraction of its own asymptotic gap that a fitted profile must actually traverse
#: inside the data for the fit to count as a genuine two-phase description.
_MIN_REALISED_FRACTION = 0.5

#: Absolute floor on the coexistence gap, in occupied-site fraction.
_MIN_DENSITY_GAP = 1e-3

#: The gap must also stand clear of the scatter in the profile: a tanh will happily
#: fit a hump in the noise of a flat profile and report the hump as a coexistence gap.
_MIN_GAP_OVER_NOISE = 3.0


def _observed_binodal(density):
    """Robust dense/dilute estimates read straight off the profile.

    Used as the fallback when the tanh fit is unusable. For a homogeneous profile the
    two values coincide, which is the correct answer for a one-phase system.
    """
    lo, hi = np.percentile(density, [5, 95])
    return float(np.clip(hi, 0.0, 1.0)), float(np.clip(lo, 0.0, 1.0))


def _fit_is_usable(density, model_values, rho_dense, rho_dilute):
    """Is this converged fit actually a two-phase description of the profile?

    Two independent ways for a ``tanh`` fit to converge on nonsense, and one check for
    each:

    **1. The asymptotes are never reached.** A ``tanh`` with a very wide interface is,
    over a finite window, almost a straight line. The optimiser can then place
    ``rho_dense`` and ``rho_dilute`` almost anywhere - typically pinned at the bounds,
    1 and 0 - because nothing in the data constrains them. The fit converges, the
    residual is small, and the reported coexistence densities are pure extrapolation,
    attained nowhere in the box. The signature is that the model, evaluated over the
    data, traverses only a small part of the gap between its own asymptotes. A genuine
    profile spans essentially all of it: it really does reach the dense plateau in the
    middle and the dilute plateau at the edges.

    **2. The gap is noise.** Even on a flat profile a ``tanh`` can find a small hump in
    the scatter and fit it as a very thin slab. The asymptotes *are* then realised, so
    check 1 passes - but the "coexistence gap" is the size of the noise. So the gap
    must also stand clear of the residual scatter.

    Returns ``(ok, reason)``.
    """
    gap = abs(float(rho_dense) - float(rho_dilute))
    if gap < _MIN_DENSITY_GAP:
        return False, "no density gap: the profile is homogeneous"

    span = float(np.max(model_values) - np.min(model_values))
    if span < _MIN_REALISED_FRACTION * gap:
        return False, (
            f"degenerate fit: the model spans only {span:.3g} of its own {gap:.3g} "
            f"asymptotic gap inside the box, so rho_dense/rho_dilute are "
            f"extrapolations that the data does not constrain"
        )

    noise = float(np.std(np.asarray(density, float) - np.asarray(model_values, float)))
    if gap < _MIN_GAP_OVER_NOISE * noise:
        return False, (
            f"density gap ({gap:.3g}) is not significant against the scatter in the "
            f"profile ({noise:.3g}): this is noise, not coexistence"
        )
    return True, ""


def fit_radial_profile(radii, density):
    """Fit a spherical droplet profile; returns a :class:`BinodalFit`.

    A fit that converges but does not describe a real two-phase profile is returned
    with ``success=False`` - see :class:`BinodalFit`.
    """
    from scipy.optimize import curve_fit
    radii = np.asarray(radii, float)
    density = np.asarray(density, float)
    rho_d0 = float(density[:max(1, len(density) // 5)].mean())
    rho_v0 = float(density[-max(1, len(density) // 5):].mean())
    r0 = float(radii[np.argmin(np.abs(density - 0.5 * (rho_d0 + rho_v0)))])
    obs_d, obs_v = _observed_binodal(density)
    try:
        p, _ = curve_fit(_tanh_droplet, radii, density,
                         p0=[rho_d0, rho_v0, r0, 1.0],
                         bounds=([0, 0, 0, 0.1], [1, 1, radii.max(), radii.max()]),
                         maxfev=10000)
    except Exception:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float("nan"), radius=r0, success=False,
                          reason="curve_fit did not converge")

    ok, reason = _fit_is_usable(density, _tanh_droplet(radii, *p), p[0], p[1])
    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float(p[3]), radius=float(p[2]),
                          success=False, reason=reason)

    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      radius=float(p[2]), interface_width=float(p[3]))


def fit_slab_profile(coord, density):
    """Fit a slab (two-interface) profile; returns a :class:`BinodalFit`.

    A fit that converges but does not describe a real two-phase profile is returned
    with ``success=False`` - see :class:`BinodalFit`. Above the critical temperature
    the profile is flat, and an unguarded ``tanh`` fit will happily report a large,
    entirely fictitious density gap for it.
    """
    from scipy.optimize import curve_fit
    coord = np.asarray(coord, float)
    density = np.asarray(density, float)
    length = coord[-1] - coord[0] + 1
    # the slab is re-centred in the middle of the window; anchor on coord[0] rather than
    # assuming the coordinate axis starts at zero (it does for slab_density_profile, but
    # not necessarily for a profile the caller built themselves)
    center = float(coord[0] + length / 2.0)
    rho_d0 = float(density.max())
    rho_v0 = float(np.median(density[density < 0.5 * density.max()])) if np.any(density < 0.5 * density.max()) else 0.0
    hw0 = float(np.sum(density > 0.5 * (rho_d0 + rho_v0)) / 2.0)
    obs_d, obs_v = _observed_binodal(density)

    def _model(z, rd, rv, hw, w):
        return _tanh_slab(z, rd, rv, hw, w, center)

    try:
        p, _ = curve_fit(_model, coord, density,
                         p0=[rho_d0, rho_v0, max(hw0, 1.0), 1.0],
                         bounds=([0, 0, 0, 0.1], [1, 1, length, length]), maxfev=10000)
    except Exception:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float("nan"), half_width=hw0, success=False,
                          reason="curve_fit did not converge")

    ok, reason = _fit_is_usable(density, _model(coord, *p), p[0], p[1])

    # A "slab" that fills the box has no dilute phase to coexist with: the fitted
    # dilute density is then unconstrained, whatever the profile happens to look like.
    if ok and 2.0 * float(p[2]) >= length:
        ok, reason = False, "degenerate fit: the slab fills the box, leaving no dilute phase"

    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float(p[3]), half_width=float(p[2]),
                          success=False, reason=reason)

    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      half_width=float(p[2]), interface_width=float(p[3]))


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
        """Heuristic: a *usable* binodal fit, a clear density gap, and most material
        condensed.

        The ``binodal.success`` requirement is load-bearing. Without it, a homogeneous
        (supercritical) system passes: the unguarded tanh fit invents a large density
        gap, and ``condensed_fraction`` is no help either, because at moderate density
        the contact-based clustering *percolates* and reports nearly all the material
        as "condensed" even when there is no condensate at all. A spanning network is
        not a droplet.
        """
        b = self.binodal
        return bool(b.success
                    and np.isfinite(b.rho_dense) and np.isfinite(b.rho_dilute)
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
        coord, dens = slab_density_profile(traj)
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
