## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Surface-tension estimation from interfacial undulations (capillary-wave theory).

Two geometries, both driven by the fluctuation spectrum of the condensate's
interface:

* **Slab** (`slab_surface_tension`) - the condensate spans the in-plane
  directions and is bounded along one axis, giving two nearly-flat interfaces whose
  height field ``h(x, y)`` obeys ``<|h(q)|^2> = kT / (gamma A q^2)``. Fitting the
  low-``q`` capillary spectrum gives ``gamma``. This is the robust method. The
  modes are plane waves under periodic boundaries and cosine modes between hard
  walls.

* **Droplet** (`droplet_surface_tension`) - a compact cluster whose radius
  ``R(theta, phi)`` fluctuates in spherical-harmonic modes with
  ``<|u_lm|^2> = kT / (gamma R0^2 (l-1)(l+2))`` for ``l >= 2``. Best-effort: it needs
  a single, well-formed, reasonably large droplet and many frames, and even then
  it reads low (see :func:`droplet_surface_tension`).

Because PIMMS uses ``exp(-dE / T)`` (k_B = 1, energies in interaction units), the
temperature is ``k_B T`` directly and ``gamma`` comes out in **reduced units**
(interaction energy per lattice area). Temperature is taken from the trajectory
(``traj.temperature``: the keyfile ``TEMPERATURE``, or ``QUENCH_END`` for a quench)
unless passed explicitly. Both estimators need a 3D system.
"""

from dataclasses import dataclass, field

import numpy as np

# the spanning detector is shared with the phase-separation module: a cluster
# that winds the box is a network or a slab, not a droplet, and the droplet
# estimator has to know that before it gathers the cluster into one image
from .phase_separation import (_cluster_spans_box, _largest_cluster_spanning_axes,
                               _normalise_axis, _vote_slab_normal)


@dataclass
class SurfaceTension:
    """Result of a capillary-wave surface-tension estimate (reduced units).

    Attributes
    ----------
    gamma : float
        Surface tension in reduced units (interaction energy per lattice area).
        ``nan`` when no frame could be used, ``inf`` for a perfectly flat
        interface.
    method : str
        Which estimator produced it, ``'slab'`` or ``'droplet'``.
    temperature : float
        The ``k_B T`` used, in PIMMS reduced units.
    n_modes : int
        Number of interfacial modes the estimate averages over: independent
        Fourier wavevectors for the slab (``+q`` and ``-q`` count once),
        spherical-harmonic degrees ``l`` for the droplet. ``0`` when the
        estimate failed (``gamma`` is ``nan``) or the slab interface was
        perfectly flat (``gamma`` is ``inf``).
    gamma_std : float
        Spread of the per-mode gammas, an uncertainty proxy rather than a
        standard error.
    spectrum : tuple or None
        The fitted spectrum for plotting: ``(q, P(q))`` for the slab method,
        where ``q`` is the lattice wavenumber
        ``sqrt((2 - 2 cos qx) + (2 - 2 cos qy))`` and ``P`` the frame- and
        face-averaged ``|FFT(delta h)|^2`` (under a hardwall the cosine-mode
        power, scaled so that ``P q^2 = N kT / gamma`` holds there too);
        ``(l, <|u_l|^2>)`` for the droplet method. ``None`` when ``gamma`` is
        ``nan``.
    n_polar : int
        Droplet method only: number of polar bins in the angular grid used.
    n_azim : int
        Droplet method only: number of azimuthal bins in the angular grid used.
    """
    gamma: float
    method: str
    temperature: float
    n_modes: int
    gamma_std: float = float("nan")           # spread across modes (uncertainty proxy)
    spectrum: tuple = field(repr=False, default=None)   # (q_or_l, power) for plotting/inspection
    n_polar: int = None                       # droplet method: angular grid actually used
    n_azim: int = None


def _resolve_kT(traj, temperature):
    """Settle which temperature to use, and check it is usable.

    PIMMS works in reduced units with ``k_B = 1``, so the temperature is
    ``k_B T`` directly.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory, whose ``temperature`` (from the keyfile) is the fallback.
    temperature : float or None
        An explicit override. ``None`` falls back to the trajectory.

    Returns
    -------
    float
        The temperature to use as ``k_B T``.

    Raises
    ------
    ValueError
        If neither source supplies a temperature, or if it is not a finite
        positive number.
    """
    kT = temperature if temperature is not None else traj.temperature
    if kT is None:
        raise ValueError("temperature is unknown - pass temperature=... or load with a keyfile")
    try:
        kT = float(kT)
    except (TypeError, ValueError):
        raise ValueError("temperature must be a finite positive number")
    if not np.isfinite(kT) or kT <= 0:
        raise ValueError("temperature must be a finite positive number")
    return kT


# ---------------------------------------------------------------------------
# slab capillary waves
# ---------------------------------------------------------------------------

def _interface_heights(cluster_positions, axis, in_plane, dims, hardwall=False):
    """Instantaneous upper/lower interface height fields h(x, y) for a slab.

    The slab is re-centred along ``axis`` (circular mean) so it does not wrap the
    boundary, then for each in-plane column the top-most / bottom-most occupied
    site gives the two interfaces. Empty columns are filled with the field mean.

    Under ``hardwall`` the box is not periodic and the slab cannot wrap, so the
    heights are taken as they are (no re-centring): a face lying on a wall then
    keeps its true height ``0`` or ``L-1`` and the caller can recognise it.

    Parameters
    ----------
    cluster_positions : numpy.ndarray
        ``(n_beads, 3)`` float64 positions of the condensate's beads.
    axis : int
        Index of the slab normal: the axis the interfaces are measured along.
    in_plane : tuple of int
        The two remaining axis indices, in the order the height field is
        indexed by.
    dims : tuple of int
        Box extent in lattice units, one entry per axis.
    hardwall : bool, optional
        ``True`` skips the circular re-centring, because a hardwall slab cannot
        wrap the boundary (default ``False``).

    Returns
    -------
    h_upper : numpy.ndarray or None
        ``(Lx, Ly)`` float64 height of the upper interface, or ``None`` if fewer
        than half the columns are occupied this frame.
    h_lower : numpy.ndarray or None
        ``(Lx, Ly)`` float64 height of the lower interface, or ``None`` under
        the same condition.
    """
    ax_i, (px, py) = axis, in_plane
    L = dims[ax_i]
    Lx, Ly = dims[px], dims[py]

    z = cluster_positions[:, ax_i].astype(np.float64)
    if hardwall:
        zc = z
    else:
        # circular centre along the axis -> shift slab to the middle
        ang = 2.0 * np.pi * z / L
        com = (np.arctan2(-(np.sin(ang).sum()), -(np.cos(ang).sum())) + np.pi) / (2.0 * np.pi) * L
        zc = np.mod(z - com + L / 2.0, L)

    ix = np.mod(np.round(cluster_positions[:, px]).astype(int), Lx)
    iy = np.mod(np.round(cluster_positions[:, py]).astype(int), Ly)

    upper = np.full((Lx, Ly), np.nan)
    lower = np.full((Lx, Ly), np.nan)
    # per-column max / min occupied height
    for x, y, zz in zip(ix, iy, zc):
        if np.isnan(upper[x, y]) or zz > upper[x, y]:
            upper[x, y] = zz
        if np.isnan(lower[x, y]) or zz < lower[x, y]:
            lower[x, y] = zz

    filled = ~np.isnan(upper)
    if filled.sum() < 0.5 * Lx * Ly:            # too holey to be a slab this frame
        return None, None
    upper[~filled] = np.nanmean(upper)
    lower[~filled] = np.nanmean(lower)
    return upper, lower


def slab_surface_tension(traj, axis=None, min_beads=2, n_modes=8, temperature=None):
    """Estimate surface tension from the slab interface capillary spectrum.

    ``gamma = N kT / <P(q) q^2>`` with ``N = Lx*Ly``, ``P`` the
    frame/interface-averaged ``|FFT(delta h)|^2`` and ``q^2`` the lattice
    dispersion ``(2 - 2 cos qx) + (2 - 2 cos qy)``, averaged over the
    ``n_modes`` lowest independent modes.

    **Hardwall.** Between hard walls the height field is not periodic in-plane,
    and its capillary modes are not plane waves: they are the modes of the
    lattice Laplacian with free ends, ``cos(pi m (x + 1/2) / Lx) cos(pi n (y +
    1/2) / Ly)``, with ``q^2 = (2 - 2 cos(pi m / Lx)) + (2 - 2 cos(pi n /
    Ly))``. The height field is projected on those, and equipartition gives the
    same ``gamma = N kT / <P q^2>`` once each mode's power is divided by its
    norm (``P = N a^2`` with ``a`` the amplitude on the normalised mode). The
    periodic transform used to be applied here as well. A plane wave is not a
    mode of the walled interface: its lowest wavevectors pick up power from
    several cosine modes and from the step between one wall and the other, and
    the estimate read ``0.55 gamma`` on synthetic slabs built from a known
    capillary spectrum, where the cosine projection reads the same as the
    periodic estimator does on periodic ones (within a few per cent of
    ``gamma``).

    .. warning::

       That is the only validation the hardwall estimator has: synthetic
       height fields. On real PIMMS hardwall slabs it is NOT reliable. Between
       walls the condensate's faces are not flat on average - they carry a
       static dome, which the projection reads as capillary power - and one
       real run read ``2.7 +- 56`` where its periodic twin reads ``41.9``.
       Subtracting the time-averaged face does not repair it. Measure the
       surface tension in a periodic box; a warning is raised when this
       estimator is used under HARDWALL.

    **Skipped frames.** The estimator is only meaningful for a slab: a largest
    cluster that spans both in-plane axes (connected to its own periodic image
    through those faces, or touching both walls under a hardwall) and does not
    span the slab normal. Any other frame is left out and a warning says how
    many there were. A cluster that also spans the normal is a network - the
    contact clustering of a homogeneous solution at moderate volume fraction
    looks like this - and its per-column top and bottom beads are not two
    interfaces; one that misses an in-plane axis is a droplet or a strip, whose
    empty columns would be filled in with the mean height. The per-column
    coverage test (fewer than half the columns occupied) still applies to the
    frames that pass.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse. Must be 3D.
    axis : int, optional
        Index of the slab normal; negative values count from the last axis.
        Default ``None``: the axis the largest cluster does not span when it
        spans the other two in more than half of the frames, and the longest
        box axis otherwise (the rule of
        :func:`pimms.lemonade.phase_separation.slab_normal`).
    min_beads : int, optional
        Ignore clusters smaller than this many beads when picking the condensate
        (default ``2``).
    n_modes : int, optional
        Number of independent lowest-``q`` modes to average over (the capillary
        regime). Under periodic boundaries the conjugate ``+q/-q`` Fourier
        coefficients of the real height field count once; the cosine modes used
        under a hardwall are all independent (default ``8``).
    temperature : float, optional
        Override the trajectory's temperature, used as ``k_B T`` (default
        ``None``).

    Returns
    -------
    SurfaceTension
        The estimate, with ``method='slab'``. ``gamma`` is ``nan`` and
        ``n_modes`` zero if no frame yielded a usable interface, and ``inf`` for
        a perfectly flat interface.

    Raises
    ------
    ValueError
        If ``n_modes`` is not a positive integer, if the system is not 3D, if
        ``axis`` is not an integer in ``-3 .. 2``, or if no temperature is
        available.
    """
    import warnings

    if isinstance(n_modes, (bool, np.bool_)) or not isinstance(n_modes, (int, np.integer)) or n_modes < 1:
        raise ValueError("n_modes must be a positive integer")

    kT = _resolve_kT(traj, temperature)
    dims = traj.dimensions
    if traj.n_dim != 3:
        raise ValueError("slab capillary-wave analysis requires a 3D system")
    # The axes each frame's largest cluster spans, found once: the same answer
    # picks the slab normal (when the caller leaves it open) and decides below
    # which frames hold a slab. Only the axes are kept, not the gathered
    # clusters, so memory does not grow with the number of frames. (The gather
    # that finds the axes warns when a periodic cluster winds the box, which a
    # slab always does in-plane; the helper mutes that warning.)
    frame_spans = _largest_cluster_spanning_axes(traj, min_beads)
    if axis is None:
        # the longest box axis used to be taken, always; a slab across a short
        # axis, or in a cubic box, was then measured along an in-plane direction
        axis = _vote_slab_normal(frame_spans, tuple(dims)[:3])
    else:
        # a negative axis used to leave in_plane = (0, 1, 2): three in-plane axes,
        # an unpacking error for some and a nan with a baffling warning for others
        axis = _normalise_axis(axis, 3, "slab_surface_tension")
    in_plane = tuple(i for i in range(3) if i != axis)
    Lx, Ly = dims[in_plane[0]], dims[in_plane[1]]
    N = Lx * Ly

    # in-plane wavevectors and the LATTICE dispersion (lattice units). The height
    # field lives on the integer lattice, so the exact Gaussian capillary spectrum
    # involves 2-2cos(q) per axis, not the continuum q^2: using q^2 over-weights
    # each mode by q^2/(2-2cos q) and systematically under-estimates gamma by
    # 2-10% at the box sizes PIMMS typically uses (5% for the lowest mode of an
    # 8-wide cross-section). 2-2cos(q) -> q^2 in the continuum limit, so nothing
    # changes for large boxes.
    hardwall = bool(getattr(traj, "hardwall", False))
    L_axis = dims[axis]
    if hardwall:
        warnings.warn(
            "slab_surface_tension: under HARDWALL this estimator is validated only on "
            "synthetic capillary spectra and is not reliable on real hardwall slabs (their "
            "faces carry a static dome that is read as capillary power); measure the surface "
            "tension in a periodic box.", stacklevel=2)
        # Between hard walls the interface is not periodic in-plane: its modes
        # are those of the lattice Laplacian with free (Neumann) ends, the
        # cosines cos(pi m (x + 1/2) / Lx), with eigenvalue 2 - 2 cos(pi m / Lx)
        # per axis. Each mode is one degree of freedom, so equipartition reads
        # gamma * lambda_mn * a_mn^2 * |phi_mn|^2 = kT for the amplitude a_mn on
        # the mode phi_mn. Storing P_mn = N (sum dh phi_mn)^2 / |phi_mn|^2 keeps
        # the periodic formula gamma = N kT / <P q^2> valid as it stands.
        cos_x = np.cos(np.pi * np.outer(np.arange(Lx), np.arange(Lx) + 0.5) / Lx)
        cos_y = np.cos(np.pi * np.outer(np.arange(Ly), np.arange(Ly) + 0.5) / Ly)
        norm_x = (cos_x ** 2).sum(axis=1)
        norm_y = (cos_y ** 2).sum(axis=1)
        mode_norm = np.outer(norm_x, norm_y)
        q2 = ((2.0 - 2.0 * np.cos(np.pi * np.arange(Lx) / Lx))[:, np.newaxis]
              + (2.0 - 2.0 * np.cos(np.pi * np.arange(Ly) / Ly))[np.newaxis, :])
    else:
        qx = 2.0 * np.pi * np.fft.fftfreq(Lx)
        qy = 2.0 * np.pi * np.fft.fftfreq(Ly)
        QX, QY = np.meshgrid(qx, qy, indexing="ij")
        q2 = (2.0 - 2.0 * np.cos(QX)) + (2.0 - 2.0 * np.cos(QY))

    power = np.zeros((Lx, Ly))
    n_used = 0
    n_network = 0                          # frames whose cluster also spans the normal
    n_not_slab = 0                         # frames whose cluster misses an in-plane axis
    n_with_cluster = 0
    for f, spans in enumerate(frame_spans):
        if spans is None:
            continue
        n_with_cluster += 1
        # Any cluster covering half the columns used to be accepted, so a
        # percolating network (30 % random occupancy of an elongated box) came
        # back with a finite gamma and no warning. The spanning axes say what the
        # cluster is.
        if axis in spans:
            n_network += 1
            continue
        if not all(a in spans for a in in_plane):
            n_not_slab += 1
            continue
        cluster = [c for c in traj[f].clusters if c.n_beads >= min_beads][0]
        pos = cluster.positions.astype(np.float64)
        for h in _interface_heights(pos, axis, in_plane, dims, hardwall=hardwall):
            if h is None:
                continue
            # Under HARDWALL a face pressed against a wall is not an interface:
            # it is flat because the wall is flat. Averaging its zero capillary
            # power in with the real face halved <P(q)> and doubled gamma for a
            # condensate wetting a wall.
            if hardwall and (np.median(h) <= 0.5 or np.median(h) >= L_axis - 1.5):
                continue
            dh = h - h.mean()
            if hardwall:
                power += N * (cos_x @ dh @ cos_y.T) ** 2 / mode_norm
            else:
                power += np.abs(np.fft.fft2(dh)) ** 2
            n_used += 1
    if n_network:
        warnings.warn(
            "slab_surface_tension: %d of %d frames with a cluster were skipped because the "
            "largest cluster also spans the slab normal (axis %d) - a network, not a slab, "
            "with no pair of interfaces to measure; see phase_separation.spanning_fraction."
            % (n_network, n_with_cluster, axis), stacklevel=2)
    if n_not_slab:
        warnings.warn(
            "slab_surface_tension: %d of %d frames with a cluster were skipped because the "
            "largest cluster does not span both in-plane axes %s - a droplet or a strip, not "
            "a slab; see phase_separation.spanning_fraction."
            % (n_not_slab, n_with_cluster, in_plane), stacklevel=2)
    if n_used == 0:
        return SurfaceTension(gamma=float("nan"), method="slab", temperature=kT, n_modes=0)
    power /= n_used

    # select the lowest-|q| non-zero modes (the capillary regime)
    flat_q2 = q2.ravel()
    flat_P = power.ravel()
    candidates = np.argsort(flat_q2)
    independent = []
    seen_pairs = set()
    for flat_idx in candidates[flat_q2[candidates] > 1e-12]:
        ix, iy = np.unravel_index(int(flat_idx), (Lx, Ly))
        # a cosine mode has no conjugate partner: every (m, n) is its own mode
        conjugate = (ix, iy) if hardwall else ((-ix) % Lx, (-iy) % Ly)
        pair = min((ix, iy), conjugate)
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        independent.append(int(flat_idx))
        if len(independent) == n_modes:
            break
    order = np.asarray(independent, dtype=np.int64)
    # In the capillary regime P(q) q^2 = N kT / gamma is constant, so average it
    # over the low-q modes and invert (a less noise-sensitive estimator than
    # averaging per-mode gammas). The per-mode spread is reported as uncertainty.
    pq2 = flat_P[order] * flat_q2[order]
    if np.mean(pq2) == 0.0:
        # a fluctuation-free (perfectly flat) interface has zero capillary
        # power: the surface-tension limit is +infinity. Return the sentinel
        # explicitly rather than letting numpy emit divide-by-zero warnings.
        return SurfaceTension(gamma=float("inf"), method="slab", temperature=kT,
                              n_modes=0, gamma_std=float("nan"),
                              spectrum=(np.sqrt(flat_q2[order]), flat_P[order]))
    gamma = float(N * kT / np.mean(pq2))
    with np.errstate(divide="ignore"):
        gamma_modes = N * kT / pq2
    return SurfaceTension(gamma=gamma, method="slab", temperature=kT,
                          n_modes=len(order),
                          gamma_std=float(np.std(gamma_modes[np.isfinite(gamma_modes)]))
                          if np.isfinite(gamma_modes).any() else float("nan"),
                          spectrum=(np.sqrt(flat_q2[order]), flat_P[order]))


# ---------------------------------------------------------------------------
# droplet shape (spherical-harmonic) fluctuations
# ---------------------------------------------------------------------------

def droplet_surface_tension(traj, l_max=5, n_polar=None, n_azim=None, min_beads=30,
                            temperature=None):
    """Estimate surface tension from droplet shape (spherical-harmonic) fluctuations.

    For each frame the largest cluster is centred on its COM and its interface
    radius ``R(theta, phi)`` sampled on an angular grid; the dimensionless
    fluctuation ``R/R0 - 1`` is projected onto real solid-angle-weighted spherical
    harmonics. ``<|u_lm|^2>`` (averaged over ``m`` and frames) is fit to
    ``kT / (gamma R0^2 (l-1)(l+2))`` for ``l = 2..l_max``.

    **Angular grid.** The interface radius in each angular bin is the radius of the
    outermost bead in it. By default the grid is chosen from the droplet size so
    that each bin holds of the order of one surface bead: ``n_polar = 2 R0``
    (rounded and clamped to ``8..64``) and ``n_azim = 2 n_polar``, with
    ``R0 = (3 n / 4 pi)^(1/3)`` and ``n`` the median bead count of the largest
    cluster over every frame that holds one of at least ``min_beads`` beads. Pass
    ``n_polar``/``n_azim`` to override; the grid actually used is reported on the
    result.

    **Accuracy.** The estimate reads low. Near the poles of the grid the
    azimuthal bins are narrower than a lattice site, so many of them hold no
    surface bead and report an inner one; the radius there reads short in every
    frame, and that deficit, being symmetric about the grid axis, lands in the
    even modes (``l = 2, 4``) as apparent fluctuation. On lattice droplets filled
    from a known capillary spectrum (modes ``l = 2 .. 1.5 R0``, the droplet
    centre at random positions relative to the lattice, 200 frames, 10-15 seeds)
    the automatic grid read ``0.89``, ``0.87`` and ``0.83 gamma`` for
    ``R0 = 8, 12, 18`` at ``gamma = kT``, and ``0.90`` and ``0.84 gamma`` at
    ``gamma = 0.5`` and ``1.5 kT`` (``R0 = 12``), with 2 - 3 % scatter between
    seeds; a droplet centred exactly on a lattice site (the test-suite oracle) is
    a special case that reads ``0.95 - 0.98 gamma``. No fixed grid does reliably
    better: ``8 x 16`` read ``0.97 - 0.99 gamma`` on those droplets but
    ``1.14 - 1.37 gamma`` on smooth droplets carrying only ``l <= 5``, and grids
    twice as fine as the automatic one ``0.5 - 1.7 gamma``. On real PIMMS
    droplets of the ``slab_phase_separation`` demo's chains (``R0`` about 9 and
    11) the automatic grid gave about three quarters of the slab estimate for
    the same chains and temperature, and other grids a half to one and a half
    times it. Treat a droplet estimate as a rough, low number and prefer the slab
    estimator whenever the geometry allows it.

    NOTE: meaningful only for a single, compact, reasonably large droplet sampled
    over many frames; small/rough/multi-droplet systems give noisy estimates.

    **Skipped frames.** Two kinds of frame are left out of the average, and a
    warning says how many of each there were. A frame whose largest cluster spans
    the box (winds through a periodic face, or touches both walls under a
    hardwall) is a network or a slab rather than a droplet: it has no closed
    interface to expand, and under periodic boundaries the single image the
    gather would return for it depends on the search order, so fitting a spectrum
    to it would be fitting noise. The usual case is frame 0 of a run that saved its
    equilibration, which is the random starting placement. A frame whose droplet
    fills fewer than half of the angular bins (too small or too rough for the
    grid) is skipped for the second reason; pass ``n_polar``/``n_azim`` to size
    the grid yourself if that happens on frames you expected to count.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse. Must be 3D.
    l_max : int, optional
        Highest spherical-harmonic degree fitted; modes run ``l = 2..l_max``
        (default ``5``). The fit needs at least two modes, so ``l_max`` below
        3 always returns ``nan``.
    n_polar : int, optional
        Number of polar bins in the angular grid, at least 2. Default ``None``,
        which sizes the grid from the median droplet size as described above.
    n_azim : int, optional
        Number of azimuthal bins, at least 4. Default ``None``, i.e.
        ``2 * n_polar``.
    min_beads : int, optional
        Ignore clusters smaller than this many beads when picking the droplet
        (default ``30``).
    temperature : float, optional
        Override the trajectory's temperature, used as ``k_B T`` (default
        ``None``).

    Returns
    -------
    SurfaceTension
        The estimate, with ``method='droplet'``, ``spectrum`` set to
        ``(l, <|u_l|^2>)`` and the angular grid actually used reported on
        ``n_polar`` / ``n_azim``. ``gamma`` is ``nan`` and ``n_modes`` zero when
        no frame was usable or fewer than two modes survived.

    Raises
    ------
    ValueError
        If the system is not 3D, if ``n_polar < 2`` or ``n_azim < 4``, or if no
        temperature is available.
    """
    import warnings

    kT = _resolve_kT(traj, temperature)
    if traj.n_dim != 3:
        raise ValueError("droplet spherical-harmonic analysis requires a 3D system")

    if n_polar is None or n_azim is None:
        # size the grid from the MEDIAN largest-cluster size over the frames that
        # will be analysed: a compact cluster of n beads has R0 ~ (3n / 4pi)^(1/3)
        # lattice units. Sizing from the first qualifying frame alone let a
        # not-yet-condensed leading frame (frame 0 is the start configuration,
        # and SAVE_EQ keeps the equilibration frames) pick a grid far too fine
        # for the real droplets, which the coverage test below then rejected
        # frame by frame, silently.
        sizes = []
        for f in range(traj.n_frames):
            clusters = [c for c in traj[f].clusters if c.n_beads >= min_beads]
            if clusters:
                sizes.append(clusters[0].n_beads)
        n_largest = float(np.median(sizes)) if sizes else 0.0
        r0_est = (3.0 * n_largest / (4.0 * np.pi)) ** (1.0 / 3.0) if n_largest else 0.0
        auto_polar = int(np.clip(round(2.0 * r0_est), 8, 64))
        if n_polar is None:
            n_polar = auto_polar
        if n_azim is None:
            n_azim = 2 * n_polar
    n_polar, n_azim = int(n_polar), int(n_azim)
    if n_polar < 2 or n_azim < 4:
        raise ValueError("n_polar must be >= 2 and n_azim >= 4")

    # scipy.special.sph_harm(m, l, azimuth, polar); the newer sph_harm_y has a
    # swapped/reordered signature, so keep the stable one and mute its deprecation.
    def _ylm(m, l, azimuth, polar):
        """Complex spherical harmonic, across both scipy signatures.

        Parameters
        ----------
        m : int
            Order, ``-l <= m <= l``.
        l : int
            Degree.
        azimuth : numpy.ndarray
            Float64 azimuthal angles in ``[0, 2 pi)``.
        polar : numpy.ndarray
            Float64 polar angles in ``[0, pi]``.

        Returns
        -------
        numpy.ndarray
            Complex128 ``Y_lm`` on the angular grid.
        """
        from scipy import special
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if hasattr(special, "sph_harm_y"):
                return special.sph_harm_y(l, m, polar, azimuth)
            return special.sph_harm(m, l, azimuth, polar)

    ls = np.arange(2, l_max + 1)
    # angular grid + solid-angle weights
    polar = (np.arange(n_polar) + 0.5) * np.pi / n_polar
    azim = (np.arange(n_azim) + 0.5) * 2.0 * np.pi / n_azim
    P, Az = np.meshgrid(polar, azim, indexing="ij")
    weight = np.sin(P) * (np.pi / n_polar) * (2.0 * np.pi / n_azim)

    # precompute conj(Y_lm) on the grid for every (l, m)
    Yconj = {}
    for l in ls:
        for m in range(-l, l + 1):
            Yconj[(l, m)] = np.conj(_ylm(m, l, Az, P))

    u2_sum = {l: 0.0 for l in ls}          # sum over frames of mean_m |u_lm|^2
    r0_sq_sum = 0.0
    n_used = 0
    n_skipped = 0                          # frames whose droplet did not cover the grid
    n_spanning = 0                         # frames whose largest cluster spans the box
    for f in range(traj.n_frames):
        clusters = [c for c in traj[f].clusters if c.n_beads >= min_beads]
        if not clusters:
            continue
        # A cluster that spans the box - winding through a periodic face, or
        # touching both walls under a hardwall - is a network or a slab, not a
        # droplet. It has no closed interface, so there is no radius to expand in
        # spherical harmonics, and under periodic boundaries the single image the
        # gather would hand back is search-order dependent anyway. Skip the frame
        # and report it below rather than fit a spectrum to noise. The usual
        # culprit is frame 0 of a run that saved its equilibration: the random
        # starting placement percolates as a contact network.
        if _cluster_spans_box(clusters[0], False):
            n_spanning += 1
            continue
        pos = clusters[0].single_image_positions()
        rel = pos - pos.mean(axis=0)
        r = np.sqrt((rel ** 2).sum(axis=1))
        good = r > 1e-6
        rel, r = rel[good], r[good]
        theta = np.arccos(np.clip(rel[:, 2] / r, -1.0, 1.0))
        phi = np.mod(np.arctan2(rel[:, 1], rel[:, 0]), 2.0 * np.pi)

        # interface radius per angular bin = max r in the bin
        pi_idx = np.clip((theta / np.pi * n_polar).astype(int), 0, n_polar - 1)
        ai_idx = np.clip((phi / (2.0 * np.pi) * n_azim).astype(int), 0, n_azim - 1)
        R = np.zeros((n_polar, n_azim))
        np.maximum.at(R, (pi_idx, ai_idx), r)
        filled = R > 0
        if filled.sum() < 0.5 * n_polar * n_azim:
            n_skipped += 1
            continue
        R[~filled] = R[filled].mean()

        w = weight / weight.sum() * (4.0 * np.pi)          # normalise weights to 4 pi
        R0 = (R * w).sum() / (4.0 * np.pi)
        delta = R / R0 - 1.0
        r0_sq_sum += R0 * R0
        for l in ls:
            acc = 0.0
            for m in range(-l, l + 1):
                u = (delta * Yconj[(l, m)] * w).sum()
                acc += (u.real ** 2 + u.imag ** 2)
            u2_sum[l] += acc / (2 * l + 1)                 # average over m
        n_used += 1

    if n_spanning:
        warnings.warn(
            "droplet_surface_tension: %d of %d frames with a cluster were skipped because "
            "the largest cluster spans the box - a network or a slab, not a droplet, with "
            "no closed interface to expand; see phase_separation.spanning_fraction."
            % (n_spanning, n_spanning + n_skipped + n_used), stacklevel=2)
    if n_skipped:
        warnings.warn(
            "droplet_surface_tension: %d of %d frames with a cluster were skipped because the "
            "droplet filled fewer than half of the %dx%d angular bins (too small, too rough, or "
            "not a single droplet); pass n_polar/n_azim to size the grid yourself."
            % (n_skipped, n_skipped + n_used, n_polar, n_azim), stacklevel=2)
    if n_used == 0:
        return SurfaceTension(gamma=float("nan"), method="droplet", temperature=kT, n_modes=0,
                              n_polar=n_polar, n_azim=n_azim)

    u2 = np.array([u2_sum[l] / n_used for l in ls])
    R0_sq = r0_sq_sum / n_used
    # 1/<|u_l|^2> = (gamma R0^2 / kT) (l-1)(l+2)  -> slope through origin
    x = (ls - 1) * (ls + 2)
    with np.errstate(divide="ignore"):
        y = 1.0 / u2
    valid = np.isfinite(y) & (u2 > 0)
    if valid.sum() < 2:
        return SurfaceTension(gamma=float("nan"), method="droplet", temperature=kT, n_modes=0,
                              n_polar=n_polar, n_azim=n_azim)
    slope = float(np.sum(x[valid] * y[valid]) / np.sum(x[valid] ** 2))   # least-squares through origin
    gamma = slope * kT / R0_sq
    # per-mode gammas for an uncertainty proxy
    gamma_l = kT / (R0_sq * x[valid] * u2[valid])
    return SurfaceTension(gamma=float(gamma), method="droplet", temperature=kT,
                          n_modes=int(valid.sum()), gamma_std=float(np.std(gamma_l)),
                          spectrum=(ls, u2), n_polar=n_polar, n_azim=n_azim)


# ---------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------

def surface_tension(traj, geometry="auto", temperature=None, **kwargs):
    """Estimate surface tension, dispatching on geometry.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    geometry : str, optional
        ``'slab'``, ``'droplet'`` (or its synonym ``'sphere'``), or ``'auto'``
        (default), which picks slab when one box axis is at least 1.5 times the
        shortest and droplet otherwise.
    temperature : float, optional
        Override the trajectory's temperature, used as ``k_B T`` (default
        ``None``).
    kwargs : dict
        Keyword arguments (``**kwargs``) passed straight through to the chosen
        estimator
        (:func:`slab_surface_tension` or :func:`droplet_surface_tension`), so
        ``axis``, ``n_modes``, ``l_max``, ``n_polar``, ``n_azim`` and
        ``min_beads`` are set here.

    Returns
    -------
    SurfaceTension
        Whatever the chosen estimator returned.

    Raises
    ------
    ValueError
        If ``geometry`` is not one of the accepted values.
    """
    dims = traj.dimensions
    if geometry == "auto":
        geometry = "slab" if max(dims) >= 1.5 * min(dims) else "droplet"
    if geometry == "slab":
        return slab_surface_tension(traj, temperature=temperature, **kwargs)
    if geometry in ("droplet", "sphere"):
        # "sphere" accepted as a synonym (phase_separation.analyze's vocabulary)
        return droplet_surface_tension(traj, temperature=temperature, **kwargs)
    raise ValueError(
        "surface_tension: unknown geometry %r (use 'auto', 'slab', 'droplet' or "
        "the synonym 'sphere')" % (geometry,))
