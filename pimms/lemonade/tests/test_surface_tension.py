"""
Tests for capillary-wave surface-tension estimation.

Surface tension is inherently noisy for small lattice condensates, so the tests on
real PIMMS output assert the machinery, the units plumbing (temperature = kT), and
a *positive, finite* estimate on a genuinely phase-separated slab / droplet. The
value itself is checked against synthetic interfaces drawn from the capillary-wave
spectrum of a known ``gamma`` (the oracle), for both the slab and the droplet
estimator, and the slab estimator must refuse a structure that is not a slab.
"""

import warnings

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import surface_tension as st
from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology
from pimms.lemonade.trajectory import LatticeTrajectory


def test_temperature_is_loaded(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert traj.temperature == pytest.approx(22.0)


def test_slab_surface_tension_positive(traj_slab_files):
    xtc, pdb, keyfile = traj_slab_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    # frame 0 is the compact configuration equilibrated under hard walls in the
    # 12^3 box before the resize: it does not wind the in-plane faces, so it is
    # not a slab yet and is skipped (it used to pass on column coverage alone)
    with pytest.warns(UserWarning, match="does not span both in-plane axes"):
        result = st.slab_surface_tension(traj)
    assert isinstance(result, st.SurfaceTension)
    assert result.method == "slab"
    assert result.temperature == pytest.approx(26.0)
    assert np.isfinite(result.gamma) and result.gamma > 0
    assert result.n_modes > 0


def test_slab_counts_equal_power_wavevectors_as_distinct_modes(monkeypatch):
    class Cluster:
        n_beads = 100
        positions = np.zeros((100, 3), dtype=float)

        def spanning_axes(self):
            # a slab: spans both in-plane axes, not the normal (axis 2)
            return [0, 1]

    class Frame:
        clusters = [Cluster()]

    class Trajectory:
        dimensions = np.array([8, 8, 20])
        n_dim = 3
        n_frames = 1
        temperature = 2.0

        def __getitem__(self, index):
            return Frame()

    x = np.arange(8)[:, None]
    y = np.arange(8)[None, :]
    heights = np.cos(2 * np.pi * x / 8) + np.cos(2 * np.pi * y / 8)
    monkeypatch.setattr(st, "_interface_heights", lambda *args, **kwargs: (heights, None))

    result = st.slab_surface_tension(Trajectory(), n_modes=2)

    # The x and y modes happen to have identical power, but they are independent
    # wavevectors. Value-based np.unique used to collapse them to one mode.
    assert result.n_modes == 2
    assert len(result.spectrum[0]) == 2
    assert np.isfinite(result.gamma)


@pytest.mark.parametrize("n_modes", [0, -1, 1.5, True])
def test_slab_rejects_invalid_mode_count(traj_slab_files, n_modes):
    xtc, pdb, keyfile = traj_slab_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    with pytest.raises(ValueError, match="positive integer"):
        st.slab_surface_tension(traj, n_modes=n_modes)


def test_droplet_surface_tension_runs(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    # frame 0 of this fixture is the percolating start configuration, which the
    # estimator skips and reports; the test below pins that behaviour
    with pytest.warns(UserWarning, match="spans the box"):
        result = st.droplet_surface_tension(traj, l_max=4, min_beads=20)
    assert result.method == "droplet"
    assert result.temperature == pytest.approx(22.0)
    assert np.isfinite(result.gamma)
    # spectrum is (l, <|u_l|^2>) for l = 2..l_max
    ls, u2 = result.spectrum
    assert list(ls) == [2, 3, 4]
    assert np.all(u2 > 0)


def test_dispatch_by_geometry(traj_condensed_files, traj_slab_files):
    cxtc, cpdb, ckey = traj_condensed_files
    sxtc, spdb, skey = traj_slab_files
    cubic = lemonade.load(xtc=cxtc, pdb=cpdb, keyfile=ckey)
    slab = lemonade.load(xtc=sxtc, pdb=spdb, keyfile=skey)
    with pytest.warns(UserWarning, match="spans the box"):    # percolating frame 0
        assert st.surface_tension(cubic).method == "droplet"  # cubic box
    with pytest.warns(UserWarning, match="does not span both in-plane axes"):  # frame 0
        assert st.surface_tension(slab).method == "slab"      # elongated box


def test_requires_temperature(traj_condensed_files):
    xtc, pdb, _keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb)                    # no keyfile -> no temperature
    assert traj.temperature is None
    with pytest.raises(ValueError, match="temperature"):
        st.droplet_surface_tension(traj)
    # explicit temperature override works
    with pytest.warns(UserWarning, match="spans the box"):    # percolating frame 0
        result = st.droplet_surface_tension(traj, temperature=22.0, l_max=4, min_beads=20)
    assert result.temperature == pytest.approx(22.0)


@pytest.mark.parametrize("temperature", [0, -1, float("nan"), float("inf"), "hot"])
def test_rejects_nonphysical_temperature(traj_condensed_files, temperature):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    with pytest.raises(ValueError, match="finite positive"):
        st.droplet_surface_tension(
            traj, temperature=temperature, l_max=4, min_beads=20)


def test_droplet_surface_tension_skips_frames_whose_cluster_spans_the_box(traj_condensed_files):
    """A spanning cluster is not a droplet, so the estimator must skip it and say so.

    A cluster that winds the periodic box has no closed interface, so there is no
    radius to expand in spherical harmonics, and the single image the gather
    returns for it is search-order dependent. The estimator used to gather it
    anyway (and the gather would warn that the shape was meaningless). Frame 0 of
    this fixture is the random starting placement, which percolates as a contact
    network, so it is exactly the case: the frame is skipped, reported, and the
    estimate from the remaining frames is still finite.
    """
    from pimms.lemonade import phase_separation as ps

    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert ps.spanning_fraction(traj) > 0, "fixture no longer spans; pick another case"

    with pytest.warns(UserWarning, match="spans the box"):
        result = st.droplet_surface_tension(traj, l_max=4, min_beads=20)
    assert result.method == "droplet"
    assert np.isfinite(result.gamma)


# ---------------------------------------------------------------------------
# oracle interfaces of known surface tension
# ---------------------------------------------------------------------------

def _capillary_slab(seed: int, gamma: float, kT: float, Lx: int, Ly: int, Lz: int,
                    n_frames: int, half: int) -> LatticeTrajectory:
    """Lattice slabs whose two faces are drawn from the capillary spectrum of ``gamma``.

    Each face is an independent Gaussian height field on the ``Lx x Ly`` lattice
    with ``<|h(q)|^2> = N kT / (gamma lambda(q))``, ``lambda(q) = sum_i (2 - 2 cos
    q_i)`` (the lattice Laplacian, the exact spectrum of a discrete Gaussian
    interface), rounded to whole sites. Every column is then filled between its two
    faces. This is the definition of ``gamma`` the estimator claims to invert,
    written out independently of it.

    Parameters
    ----------
    seed : int
        Seed for the height fields.
    gamma : float
        Surface tension the faces are drawn from, in reduced units.
    kT : float
        Temperature, also stored on the trajectory.
    Lx, Ly, Lz : int
        Box extent; ``z`` is the slab normal.
    n_frames : int
        Number of frames.
    half : int
        Half the mean slab thickness, in sites.

    Returns
    -------
    LatticeTrajectory
        Periodic trajectory of single-bead chains, the same bead count in every
        frame (any surplus is removed from the slab's centre plane, which leaves
        both faces intact).
    """
    rng = np.random.default_rng(seed)
    qx = 2.0 * np.pi * np.fft.fftfreq(Lx)
    qy = 2.0 * np.pi * np.fft.fftfreq(Ly)
    QX, QY = np.meshgrid(qx, qy, indexing="ij")
    lam = (2.0 - 2.0 * np.cos(QX)) + (2.0 - 2.0 * np.cos(QY))
    amp = np.zeros_like(lam)
    amp[lam > 0] = np.sqrt(kT / (gamma * lam[lam > 0]))
    X, Y = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing="ij")
    centre = Lz // 2
    frames = []
    for _ in range(n_frames):
        faces = [np.rint(np.fft.ifft2(np.fft.fft2(rng.normal(size=(Lx, Ly))) * amp).real)
                 .astype(int) for _face in range(2)]
        lo = (centre - half - faces[0]).ravel()
        hi = (centre + half + faces[1]).ravel()
        n_col = hi - lo
        first = np.repeat(np.cumsum(n_col) - n_col, n_col)
        z = np.repeat(lo, n_col) + (np.arange(n_col.sum()) - first)
        frames.append(np.stack([np.repeat(X.ravel(), n_col), np.repeat(Y.ravel(), n_col), z],
                               axis=1))
    n = min(len(f) for f in frames)
    frames = [np.delete(f, np.nonzero(f[:, 2] == centre)[0][:len(f) - n], axis=0)
              for f in frames]
    store = TrajectoryStore(np.stack(frames).astype(np.int32), (Lx, Ly, Lz), 3.65, False,
                            Topology(["A"] * n), temperature=kT)
    return LatticeTrajectory(store)


def _ylm(l: int, m: int, polar: np.ndarray, azimuth: np.ndarray) -> np.ndarray:
    """Complex spherical harmonic ``Y_lm``, across both scipy signatures.

    Parameters
    ----------
    l : int
        Degree.
    m : int
        Order, ``0 <= m <= l`` here.
    polar : numpy.ndarray
        Polar angles in ``[0, pi]``.
    azimuth : numpy.ndarray
        Azimuthal angles in ``[0, 2 pi)``.

    Returns
    -------
    numpy.ndarray
        Complex ``Y_lm`` at each angle.
    """
    from scipy import special
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if hasattr(special, "sph_harm_y"):
            return special.sph_harm_y(l, m, polar, azimuth)
        return special.sph_harm(m, l, azimuth, polar)


def _capillary_droplet(seed: int, gamma: float, kT: float, R0: float, L: int,
                       n_frames: int) -> LatticeTrajectory:
    """Lattice droplets whose radius fluctuates with the capillary spectrum of ``gamma``.

    The relative radius ``R(theta, phi) / R0 - 1`` is a sum of real spherical
    harmonics with ``<|u_lm|^2> = kT / (gamma R0^2 (l - 1)(l + 2))`` for ``l = 2 ..
    1.5 R0``, and a site is filled when it lies inside that surface. Every droplet
    is centred exactly on a lattice site, which is the geometry where the
    estimator's angular-grid bias is smallest (about 0.95-0.98 gamma, against about
    0.85 for a droplet centred at a random sub-site position), so this oracle pins
    the estimator's formula rather than its accuracy on real droplets.

    Parameters
    ----------
    seed : int
        Seed for the mode amplitudes.
    gamma : float
        Surface tension the modes are drawn from, in reduced units.
    kT : float
        Temperature, also stored on the trajectory.
    R0 : float
        Mean droplet radius in sites.
    L : int
        Edge of the cubic periodic box.
    n_frames : int
        Number of frames.

    Returns
    -------
    LatticeTrajectory
        Trajectory of single-bead chains, the same bead count in every frame (any
        surplus is removed from the droplet's innermost sites).
    """
    rng = np.random.default_rng(seed)
    g = np.arange(L) - L // 2
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    band = (r > 0.4 * R0) & (r < 1.7 * R0)
    rb = r[band]
    polar = np.arccos(np.clip(Z[band] / rb, -1.0, 1.0))
    azimuth = np.mod(np.arctan2(Y[band], X[band]), 2.0 * np.pi)
    basis = []
    for l in range(2, int(round(1.5 * R0)) + 1):
        var = kT / (gamma * R0 ** 2 * (l - 1) * (l + 2))
        basis.append((var, [(m, _ylm(l, m, polar, azimuth)) for m in range(l + 1)]))
    frames = []
    for _ in range(n_frames):
        delta = np.zeros_like(rb)
        for var, harmonics in basis:
            for m, y in harmonics:
                if m == 0:
                    delta += rng.normal(0.0, np.sqrt(var)) * y.real
                else:
                    a, b = rng.normal(0.0, np.sqrt(var / 2.0), size=2)
                    delta += 2.0 * ((a + 1j * b) * y).real
        inside = r <= 0.4 * R0
        inside[band] = rb < R0 * (1.0 + delta)
        frames.append(np.stack([X[inside], Y[inside], Z[inside]], axis=1) + L // 2)
    n = min(len(f) for f in frames)
    trimmed = []
    for f in frames:
        innermost = np.argsort(((f - L // 2) ** 2).sum(axis=1))[:len(f) - n]
        trimmed.append(np.delete(f, innermost, axis=0))
    store = TrajectoryStore(np.stack(trimmed).astype(np.int32), (L, L, L), 3.65, False,
                            Topology(["A"] * n), temperature=kT)
    return LatticeTrajectory(store)


def test_slab_surface_tension_recovers_a_known_gamma():
    """The slab estimator must return the gamma its capillary spectrum was drawn from.

    An 8 x 8 cross-section is where the lattice dispersion matters most: the
    continuum ``q^2`` in place of ``2 - 2 cos q`` reads the same slabs as 0.88
    gamma, while the estimator reads them as 0.98 - 1.01 gamma across seeds (the
    residual is the rounding of the faces to whole sites). A soft slab (gamma =
    0.1 kT) keeps that rounding small. The tolerance of 7 % is about four times
    the scatter of the estimate from 200 frames.
    """
    gamma = 0.1
    traj = _capillary_slab(seed=1, gamma=gamma, kT=1.0, Lx=8, Ly=8, Lz=48,
                           n_frames=200, half=8)
    with warnings.catch_warnings():
        warnings.simplefilter("error")                     # a clean slab: nothing skipped
        result = st.slab_surface_tension(traj)
    assert result.n_modes == 8
    assert result.gamma / gamma == pytest.approx(1.0, abs=0.07)


def test_droplet_surface_tension_recovers_a_known_gamma():
    """The droplet estimator must return the gamma its harmonic spectrum was drawn from.

    Two checks, because the fit can go wrong two ways. The estimate itself must be
    within 10 % of gamma on these site-centred droplets (it reads 0.94 - 0.98 gamma
    across seeds; real, randomly centred droplets read lower, see the
    ``droplet_surface_tension`` docstring). And every mode must give the same gamma: the spread
    of the per-mode estimates is below 5 % of gamma from the correct ``(l - 1)(l +
    2)`` stiffness, while ``l(l + 1)`` in its place spreads them by 9 - 14 % and
    reads the whole estimate as 0.86 - 0.89 gamma.
    """
    gamma = 1.0
    traj = _capillary_droplet(seed=1, gamma=gamma, kT=1.0, R0=6.0, L=24, n_frames=300)
    with warnings.catch_warnings():
        warnings.simplefilter("error")                     # one compact droplet per frame
        result = st.droplet_surface_tension(traj)
    assert result.n_modes == 4                             # l = 2..5
    assert result.gamma / gamma == pytest.approx(1.0, abs=0.10)
    assert result.gamma_std / result.gamma < 0.07


# ---------------------------------------------------------------------------
# the slab estimator refuses what is not a slab
# ---------------------------------------------------------------------------

def _occupied_traj(frames: list, dims: tuple) -> LatticeTrajectory:
    """Wrap per-frame site lists (same length) in a periodic trajectory at kT = 1.

    Parameters
    ----------
    frames : list of numpy.ndarray
        ``(n_beads, 3)`` integer sites per frame.
    dims : tuple of int
        Box extent.

    Returns
    -------
    LatticeTrajectory
        One single-bead chain per bead.
    """
    arr = np.stack(frames).astype(np.int32)
    return LatticeTrajectory(TrajectoryStore(arr, dims, 3.65, False,
                                             Topology(["A"] * arr.shape[1]), temperature=1.0))


def test_slab_surface_tension_refuses_a_network():
    """A cluster that also spans the slab normal is a network, not a slab.

    At 30 % random occupancy of an elongated box the contact clustering is one
    network that winds all three axes. Any cluster covering half the columns used
    to be accepted, so this returned gamma = 0.14 +/- 0.14 with no warning.
    """
    rng = np.random.default_rng(0)
    dims = (12, 12, 36)
    n_sites = int(np.prod(dims))
    frames = [np.stack(np.unravel_index(rng.choice(n_sites, int(0.3 * n_sites), replace=False),
                                        dims), axis=1) for _ in range(6)]
    traj = _occupied_traj(frames, dims)
    with pytest.warns(UserWarning, match="6 of 6 frames .* also spans the slab normal"):
        result = st.slab_surface_tension(traj)
    assert np.isnan(result.gamma) and result.n_modes == 0


def test_slab_surface_tension_refuses_a_strip():
    """A strip that spans only one in-plane axis is not a slab.

    It covers two thirds of the columns, which passed the old coverage test, and the
    empty columns were filled in with the mean height as if they were interface.
    """
    rng = np.random.default_rng(1)
    dims = (12, 12, 36)
    frames = []
    for _ in range(4):
        top = 20 + rng.integers(0, 2, size=(12, 8))
        frames.append(np.array([(x, y, z) for x in range(12) for y in range(8)
                                for z in range(14, int(top[x, y]))]))
    n = min(len(f) for f in frames)
    traj = _occupied_traj([f[:n] for f in frames], dims)
    with pytest.warns(UserWarning, match="4 of 4 frames .* does not span both in-plane axes"):
        result = st.slab_surface_tension(traj)
    assert np.isnan(result.gamma) and result.n_modes == 0
