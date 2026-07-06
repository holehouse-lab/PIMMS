"""
Surface-tension estimation from interfacial undulations (capillary-wave theory).

Two geometries, both driven by the fluctuation spectrum of the condensate's
interface:

* **Slab** (`slab_surface_tension`) - the condensate spans the periodic in-plane
  directions and is bounded along one axis, giving two nearly-flat interfaces whose
  height field ``h(x, y)`` obeys ``<|h(q)|^2> = kT / (gamma A q^2)``. Fitting the
  low-``q`` capillary spectrum gives ``gamma``. This is the robust method.

* **Droplet** (`droplet_surface_tension`) - a compact cluster whose radius
  ``R(theta, phi)`` fluctuates in spherical-harmonic modes with
  ``<|u_lm|^2> = kT / (gamma R0^2 (l-1)(l+2))`` for ``l >= 2``. Best-effort: it needs
  a single, well-formed, reasonably large droplet and many frames to be reliable.

Because PIMMS uses ``exp(-dE / T)`` (k_B = 1, energies in interaction units), the
temperature is ``k_B T`` directly and ``gamma`` comes out in **reduced units**
(interaction energy per lattice area). Temperature is taken from the trajectory
(the keyfile ``TEMPERATURE``) unless passed explicitly.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SurfaceTension:
    """Result of a capillary-wave surface-tension estimate (reduced units)."""
    gamma: float
    method: str
    temperature: float
    n_modes: int
    gamma_std: float = float("nan")           # spread across modes (uncertainty proxy)
    spectrum: tuple = field(repr=False, default=None)   # (q_or_l, power) for plotting/inspection


def _resolve_kT(traj, temperature):
    kT = temperature if temperature is not None else traj.temperature
    if kT is None:
        raise ValueError("temperature is unknown - pass temperature=... or load with a keyfile")
    return float(kT)


# ---------------------------------------------------------------------------
# slab capillary waves
# ---------------------------------------------------------------------------

def _interface_heights(cluster_positions, axis, in_plane, dims):
    """Instantaneous upper/lower interface height fields h(x, y) for a slab.

    The slab is re-centred along ``axis`` (circular mean) so it does not wrap the
    boundary, then for each in-plane column the top-most / bottom-most occupied
    site gives the two interfaces. Empty columns are filled with the field mean.
    Returns ``(h_upper, h_lower)`` as ``(Lx, Ly)`` arrays, or ``(None, None)`` if
    too few columns are occupied.
    """
    ax_i, (px, py) = axis, in_plane
    L = dims[ax_i]
    Lx, Ly = dims[px], dims[py]

    z = cluster_positions[:, ax_i].astype(np.float64)
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

    Returns a :class:`SurfaceTension`. ``n_modes`` is the number of lowest-``q``
    Fourier modes used (the capillary regime); ``gamma = N kT / <P(q) q^2>`` with
    ``N = Lx*Ly`` and ``P`` the frame/interface-averaged ``|FFT(delta h)|^2``.
    """
    kT = _resolve_kT(traj, temperature)
    dims = traj.dimensions
    if traj.n_dim != 3:
        raise ValueError("slab capillary-wave analysis requires a 3D system")
    if axis is None:
        axis = int(np.argmax(dims))
    in_plane = tuple(i for i in range(3) if i != axis)
    Lx, Ly = dims[in_plane[0]], dims[in_plane[1]]
    N = Lx * Ly

    # in-plane wavevectors and |q|^2 (lattice units)
    qx = 2.0 * np.pi * np.fft.fftfreq(Lx)
    qy = 2.0 * np.pi * np.fft.fftfreq(Ly)
    QX, QY = np.meshgrid(qx, qy, indexing="ij")
    q2 = QX ** 2 + QY ** 2

    power = np.zeros((Lx, Ly))
    n_used = 0
    for f in range(traj.n_frames):
        clusters = [c for c in traj[f].clusters if c.n_beads >= min_beads]
        if not clusters:
            continue
        pos = clusters[0].positions.astype(np.float64)
        for h in _interface_heights(pos, axis, in_plane, dims):
            if h is None:
                continue
            dh = h - h.mean()
            power += np.abs(np.fft.fft2(dh)) ** 2
            n_used += 1
    if n_used == 0:
        return SurfaceTension(gamma=float("nan"), method="slab", temperature=kT, n_modes=0)
    power /= n_used

    # select the lowest-|q| non-zero modes (the capillary regime)
    flat_q2 = q2.ravel()
    flat_P = power.ravel()
    order = np.argsort(flat_q2)
    order = order[flat_q2[order] > 1e-12][:n_modes]
    # In the capillary regime P(q) q^2 = N kT / gamma is constant, so average it
    # over the low-q modes and invert (a less noise-sensitive estimator than
    # averaging per-mode gammas). The per-mode spread is reported as uncertainty.
    pq2 = flat_P[order] * flat_q2[order]
    gamma = float(N * kT / np.mean(pq2))
    gamma_modes = N * kT / pq2
    return SurfaceTension(gamma=gamma, method="slab", temperature=kT,
                          n_modes=len(order), gamma_std=float(np.std(gamma_modes)),
                          spectrum=(np.sqrt(flat_q2[order]), flat_P[order]))


# ---------------------------------------------------------------------------
# droplet shape (spherical-harmonic) fluctuations
# ---------------------------------------------------------------------------

def droplet_surface_tension(traj, l_max=5, n_polar=8, n_azim=16, min_beads=30,
                            temperature=None):
    """Estimate surface tension from droplet shape (spherical-harmonic) fluctuations.

    For each frame the largest cluster is centred on its COM and its interface
    radius ``R(theta, phi)`` sampled on an angular grid; the dimensionless
    fluctuation ``R/R0 - 1`` is projected onto real solid-angle-weighted spherical
    harmonics. ``<|u_lm|^2>`` (averaged over ``m`` and frames) is fit to
    ``kT / (gamma R0^2 (l-1)(l+2))`` for ``l = 2..l_max``.

    NOTE: reliable only for a single, compact, reasonably large droplet sampled over
    many frames; small/rough/multi-droplet systems give noisy estimates. Returns a
    :class:`SurfaceTension` whose ``spectrum`` is ``(l, <|u_l|^2>)``.
    """
    import warnings

    kT = _resolve_kT(traj, temperature)
    if traj.n_dim != 3:
        raise ValueError("droplet spherical-harmonic analysis requires a 3D system")

    # scipy.special.sph_harm(m, l, azimuth, polar); the newer sph_harm_y has a
    # swapped/reordered signature, so keep the stable one and mute its deprecation.
    def _ylm(m, l, azimuth, polar):
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
    for f in range(traj.n_frames):
        clusters = [c for c in traj[f].clusters if c.n_beads >= min_beads]
        if not clusters:
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

    if n_used == 0:
        return SurfaceTension(gamma=float("nan"), method="droplet", temperature=kT, n_modes=0)

    u2 = np.array([u2_sum[l] / n_used for l in ls])
    R0_sq = r0_sq_sum / n_used
    # 1/<|u_l|^2> = (gamma R0^2 / kT) (l-1)(l+2)  -> slope through origin
    x = (ls - 1) * (ls + 2)
    y = 1.0 / u2
    valid = np.isfinite(y) & (u2 > 0)
    if valid.sum() < 2:
        return SurfaceTension(gamma=float("nan"), method="droplet", temperature=kT, n_modes=0)
    slope = float(np.sum(x[valid] * y[valid]) / np.sum(x[valid] ** 2))   # least-squares through origin
    gamma = slope * kT / R0_sq
    # per-mode gammas for an uncertainty proxy
    gamma_l = kT / (R0_sq * x[valid] * u2[valid])
    return SurfaceTension(gamma=float(gamma), method="droplet", temperature=kT,
                          n_modes=int(valid.sum()), gamma_std=float(np.std(gamma_l)),
                          spectrum=(ls, u2))


# ---------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------

def surface_tension(traj, geometry="auto", temperature=None, **kwargs):
    """Estimate surface tension, dispatching on geometry (``'slab'`` / ``'droplet'``
    / ``'auto'``; auto uses the box shape - slab if one axis is >= 1.5x the others).
    """
    dims = traj.dimensions
    if geometry == "auto":
        geometry = "slab" if max(dims) >= 1.5 * min(dims) else "droplet"
    if geometry == "slab":
        return slab_surface_tension(traj, temperature=temperature, **kwargs)
    return droplet_surface_tension(traj, temperature=temperature, **kwargs)
