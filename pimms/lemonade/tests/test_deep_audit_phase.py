"""
Regression tests for the deep-audit findings in the lemonade phase-separation and
surface-tension code (auditor A2: A2-1, A2-2, A2-3, A2-4, A2-9, A2-10, A2-11, and
the four surviving mutants of auditor T1's finding T7).

Every trajectory here is built by hand with numpy and has an answer that follows
from how it was built: an athermal solution (independent self-avoiding walks,
no interaction of any kind) is one phase, a slab of stated densities has those
densities, a height field drawn from a capillary spectrum of stated tension has
that tension. Nothing is compared against the routine under test.
"""

import warnings

import numpy as np
import pytest

from pimms.lemonade import phase_separation as ps
from pimms.lemonade import surface_tension as st
from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology
from pimms.lemonade.trajectory import LatticeTrajectory

# the 26 Chebyshev-neighbour steps a PIMMS bond may take in 3D
_STEPS = np.array(
    [
        (a, b, c)
        for a in (-1, 0, 1)
        for b in (-1, 0, 1)
        for c in (-1, 0, 1)
        if (a, b, c) != (0, 0, 0)
    ]
)


def _make(frames, sequences, dims, hardwall=False, temperature=1.0):
    """Wrap hand-built integer positions in a LatticeTrajectory.

    Parameters
    ----------
    frames : array_like
        ``(n_frames, n_beads, n_dim)`` integer lattice positions.
    sequences : list of str
        One sequence per chain, in bead order.
    dims : tuple of int
        Box extent in lattice units.
    hardwall : bool, optional
        Whether the box has hard walls (default ``False``).
    temperature : float, optional
        Temperature stored on the trajectory (default ``1.0``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory.
    """
    frames = np.asarray(frames, dtype=np.int32)
    if frames.shape[-1] == 2:
        frames = np.concatenate(
            [frames, np.zeros(frames.shape[:-1] + (1,), np.int32)], axis=-1
        )
    store = TrajectoryStore(
        frames,
        tuple(dims),
        3.65,
        hardwall,
        Topology(list(sequences)),
        temperature=temperature,
    )
    return LatticeTrajectory(store)


def _athermal_frame(rng, dims, n_chains, length, hardwall):
    """One configuration of an athermal polymer solution.

    Each chain is a self-avoiding walk of unit Chebyshev steps from a uniformly
    random start, grown independently of every other chain apart from excluded
    volume. There is no interaction, so the solution is one phase by
    construction, whatever a density profile of it looks like.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    dims : tuple of int
        Box extent in lattice units (3D).
    n_chains : int
        Number of chains.
    length : int
        Beads per chain.
    hardwall : bool
        ``True`` keeps every bead inside the box; ``False`` wraps periodically.

    Returns
    -------
    numpy.ndarray
        ``(n_chains * length, 3)`` integer positions, chain by chain.
    """
    dims = np.asarray(dims)
    occupied = set()
    beads = []
    for _ in range(n_chains):
        while True:
            walk = [tuple(int(v) for v in rng.integers(0, dims))]
            if walk[0] in occupied:
                continue
            taken = {walk[0]}
            stuck = False
            for _ in range(length - 1):
                for step in _STEPS[rng.permutation(len(_STEPS))]:
                    site = np.asarray(walk[-1]) + step
                    if hardwall:
                        if np.any(site < 0) or np.any(site >= dims):
                            continue
                    else:
                        site = site % dims
                    site = tuple(int(v) for v in site)
                    if site not in occupied and site not in taken:
                        walk.append(site)
                        taken.add(site)
                        break
                else:
                    stuck = True
                    break
            if not stuck:
                break
        occupied |= taken
        beads.extend(walk)
    return np.asarray(beads)


def _athermal_traj(seed, dims, n_chains, length, hardwall, n_frames=30):
    """An athermal (one-phase) trajectory of independent configurations.

    Parameters
    ----------
    seed : int
        Seed of the generator.
    dims : tuple of int
        Box extent in lattice units (3D).
    n_chains : int
        Number of chains.
    length : int
        Beads per chain.
    hardwall : bool
        Whether the box has hard walls.
    n_frames : int, optional
        Number of frames (default ``30``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory.
    """
    rng = np.random.default_rng(seed)
    frames = [
        _athermal_frame(rng, dims, n_chains, length, hardwall) for _ in range(n_frames)
    ]
    return _make(frames, ["A" * length] * n_chains, dims, hardwall=hardwall)


def _slab_frame(rng, dims, axis, planes, rho_dense, rho_dilute):
    """One frame of a slab: random occupancy at two stated densities.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    dims : tuple of int
        Box extent in lattice units.
    axis : int
        The slab normal.
    planes : iterable of int
        Coordinates along ``axis`` of the dense planes.
    rho_dense : float
        Occupied fraction of every dense plane (exact, per plane).
    rho_dilute : float
        Occupied fraction of every other plane (exact, per plane).

    Returns
    -------
    numpy.ndarray
        ``(n_beads, n_dim)`` integer positions.
    """
    dims = tuple(dims)
    others = [d for i, d in enumerate(dims) if i != axis]
    cross = int(np.prod(others))
    planes = set(planes)
    out = []
    for z in range(dims[axis]):
        n = int(round((rho_dense if z in planes else rho_dilute) * cross))
        idx = rng.choice(cross, n, replace=False)
        sub = np.stack(np.unravel_index(idx, others), axis=-1)
        out.append(np.insert(sub, axis, z, axis=1))
    return np.concatenate(out)


def _slab_traj(
    seed, dims, axis, plane_sets, rho_dense=0.75, rho_dilute=0.0625, hardwall=False
):
    """A slab trajectory with the dense planes given frame by frame.

    Parameters
    ----------
    seed : int
        Seed of the generator.
    dims : tuple of int
        Box extent in lattice units.
    axis : int
        The slab normal.
    plane_sets : list of iterable of int
        Per frame, the coordinates of the dense planes (the same number in
        every frame, so the bead count is constant).
    rho_dense : float, optional
        Occupied fraction of the dense planes (default ``0.75``).
    rho_dilute : float, optional
        Occupied fraction of the other planes (default ``0.0625``).
    hardwall : bool, optional
        Whether the box has hard walls (default ``False``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory, every bead a one-bead chain.
    """
    rng = np.random.default_rng(seed)
    frames = [
        _slab_frame(rng, dims, axis, planes, rho_dense, rho_dilute)
        for planes in plane_sets
    ]
    return _make(frames, ["A"] * len(frames[0]), dims, hardwall=hardwall)


def _quiet_analyze(traj, **kwargs):
    """``ps.analyze`` with warnings muted.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    kwargs : dict
        Keyword arguments passed through to :func:`ps.analyze`.

    Returns
    -------
    pimms.lemonade.phase_separation.PhaseSeparationResult
        The result.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ps.analyze(traj, **kwargs)


# ---------------------------------------------------------------------------
# A2-1: a one-phase solution in slab geometry is not phase separated
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "hardwall, dims, n_chains, length",
    [
        (False, (12, 12, 40), 4, 40),  # periodic: re-centring builds a hump
        (False, (12, 12, 40), 6, 40),
        (False, (10, 10, 60), 5, 40),
        (True, (12, 12, 40), 12, 40),  # hardwall: depletion layers at the walls
    ],
)
def test_athermal_solution_is_never_phase_separated_in_slab_geometry(
    hardwall, dims, n_chains, length
):
    """An athermal solution is one phase, so ``analyze`` must say so, seed after seed.

    The baseline reported ``is_phase_separated=True`` for 3-4 of these five
    seeds in each periodic case and for one in the hardwall case.
    """
    for seed in range(5):
        res = _quiet_analyze(_athermal_traj(seed, dims, n_chains, length, hardwall))
        assert res.geometry == "slab"
        assert res.binodal.success is False, (seed, res.binodal)
        assert res.binodal.reason
        assert res.is_phase_separated is False, (seed, res.binodal)
        # the fallback densities are percentiles of the observed profile
        profile = res.profile[1][np.isfinite(res.profile[1])]
        assert profile.min() <= res.rho_dilute <= res.rho_dense <= profile.max()


def test_wall_depletion_layer_is_not_a_dilute_phase():
    """The unaligned profile of a real athermal PIMMS run between hard walls.

    These 48 numbers are the plain per-plane mean of an athermal run (61 chains
    of 20 beads, 16 x 16 x 48, HARDWALL, 151 frames; auditor A2's run
    ``onephase_confirm/c4``): a flat bulk at 0.105 with a ramp of three or four
    planes at each wall. The two-interface tanh used to fit the ramps as the
    edges of a slab filling the box and report ``success`` with a dilute density
    far below anything in the bulk.
    """
    rho = np.array(
        [
            0.022,
            0.056,
            0.083,
            0.099,
            0.103,
            0.109,
            0.109,
            0.109,
            0.103,
            0.106,
            0.103,
            0.101,
            0.104,
            0.105,
            0.107,
            0.105,
            0.110,
            0.107,
            0.103,
            0.105,
            0.108,
            0.109,
            0.111,
            0.110,
            0.111,
            0.112,
            0.111,
            0.110,
            0.107,
            0.106,
            0.105,
            0.103,
            0.102,
            0.103,
            0.103,
            0.102,
            0.102,
            0.101,
            0.105,
            0.104,
            0.107,
            0.107,
            0.108,
            0.107,
            0.101,
            0.082,
            0.056,
            0.023,
        ]
    )
    z = np.arange(rho.size, dtype=float)
    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success is False
    assert "no dilute plateau" in fit.reason
    # percentile fallback, bounded by what was observed
    assert rho.min() <= fit.rho_dilute < fit.rho_dense <= rho.max()


def test_gap_must_be_resolved_in_single_frames():
    """The frame-to-frame scatter decides between a slab and an aligned hump.

    The profile is an exact two-interface tanh with a gap of 0.04, which the
    averaged-profile checks cannot fault. If one plane's density scatters by
    0.03 from frame to frame the gap is 0.94 standard deviations of the
    difference between a dense and a dilute plane (0.04 / (0.03 sqrt 2)), far
    below three: a hump. If the scatter is 0.002 the gap is 14 of them: a slab.
    """
    z = np.arange(40, dtype=float)
    rho = 0.05 + 0.5 * 0.04 * (np.tanh((z - 12.0) / 1.5) - np.tanh((z - 28.0) / 1.5))
    noisy = ps.fit_slab_profile(z, rho, hardwall=False, frame_scatter=np.full(40, 0.03))
    assert noisy.success is False
    assert "not resolved in single frames" in noisy.reason

    quiet = ps.fit_slab_profile(
        z, rho, hardwall=False, frame_scatter=np.full(40, 0.002)
    )
    assert quiet.success is True
    assert quiet.rho_dense == pytest.approx(0.09, abs=1e-4)
    assert quiet.rho_dilute == pytest.approx(0.05, abs=1e-4)
    assert quiet.half_width == pytest.approx(8.0, abs=1e-3)

    # the threshold itself: three standard deviations of the plane difference
    edge = 0.04 / (3.0 * np.sqrt(2.0))
    assert (
        ps.fit_slab_profile(
            z, rho, hardwall=False, frame_scatter=np.full(40, 0.98 * edge)
        ).success
        is True
    )
    assert (
        ps.fit_slab_profile(
            z, rho, hardwall=False, frame_scatter=np.full(40, 1.02 * edge)
        ).success
        is False
    )

    with pytest.raises(ValueError, match="frame_scatter"):
        ps.fit_slab_profile(z, rho, hardwall=False, frame_scatter=np.zeros(7))


def test_frame_scatter_is_the_standard_deviation_over_frames():
    """``slab_density_profile_with_scatter`` against its definition.

    A one-phase hardwall trajectory is not aligned, so plane ``z`` of the
    profile is the mean over frames of (beads in plane z) / (sites per plane),
    and the scatter is the standard deviation of the same numbers.
    """
    dims = (12, 12, 40)
    traj = _athermal_traj(11, dims, 12, 40, hardwall=True)
    pos = traj.positions
    per_frame = np.array(
        [np.bincount(pos[f][:, 2], minlength=40) for f in range(traj.n_frames)]
    )
    per_frame = per_frame / 144.0
    z, rho, scatter = ps.slab_density_profile_with_scatter(traj, axis=2)
    assert np.array_equal(z, np.arange(40.0))
    assert rho == pytest.approx(per_frame.mean(axis=0), abs=1e-12)
    assert scatter == pytest.approx(per_frame.std(axis=0, ddof=1), abs=1e-12)
    z2, rho2 = ps.slab_density_profile(traj, axis=2)
    assert np.array_equal(rho, rho2) and np.array_equal(z, z2)


@pytest.mark.parametrize(
    "hardwall, plane_sets",
    [
        (False, [range(15, 25)] * 12),  # periodic slab
        (True, [range(15, 25)] * 12),  # free hardwall slab
        (
            True,
            [range(lo, lo + 10) for lo in (3, 9, 14, 20, 26, 28, 6, 17, 22, 11, 25, 4)],
        ),
        (True, [range(0, 10)] * 12),  # wets the low wall
        (True, [range(30, 40)] * 12),  # wets the high wall
        (False, [range(18, 21)] * 12),  # three planes
    ],
)
def test_genuine_slabs_are_still_phase_separated(hardwall, plane_sets):
    """The new checks must not cost a single real slab: densities 0.75 / 0.0625."""
    dims = (8, 8, 40)
    traj = _slab_traj(5, dims, 2, plane_sets, hardwall=hardwall)
    res = _quiet_analyze(traj)
    thickness = len(list(plane_sets[0]))
    assert res.binodal.success is True, res.binodal
    assert res.slab_axis == 2
    assert res.rho_dense == pytest.approx(0.75, abs=0.02)
    assert res.rho_dilute == pytest.approx(0.0625, abs=0.01)
    assert 2.0 * res.binodal.half_width == pytest.approx(thickness, abs=0.6)
    assert res.is_phase_separated is True


@pytest.mark.parametrize(
    "hardwall, dims, n_chains, length, geometry",
    [
        (False, (16, 16, 16), 3, 40, "auto"),  # cubic box: auto picks droplet geometry
        (False, (20, 20, 20), 6, 40, "auto"),
        (False, (12, 12, 40), 4, 40, "sphere"),
        (True, (12, 12, 40), 4, 40, "sphere"),
    ],
)
def test_athermal_solution_is_never_phase_separated_in_droplet_geometry(
    hardwall, dims, n_chains, length, geometry
):
    """The same one-phase solutions, profiled about their largest cluster.

    The largest cluster of a dilute solution is a chain or a few touching
    chains, and a chain is denser than its surroundings: the baseline reported a
    droplet of radius 3 and density 0.15 - 0.2, with ``is_phase_separated=True``,
    for two to four of these six seeds in every case.
    """
    for seed in range(6):
        res = _quiet_analyze(
            _athermal_traj(seed, dims, n_chains, length, hardwall), geometry=geometry
        )
        assert res.geometry == "sphere"
        assert res.binodal.success is False, (seed, res.binodal)
        assert res.is_phase_separated is False, (seed, res.binodal)


def test_radial_fit_gap_must_be_resolved_in_single_frames():
    """An exact tanh droplet profile with a gap of 0.15, the size of a coil's.

    With shells that scatter by 0.08 from frame to frame the gap is 1.3 standard
    deviations of the difference between a core shell and an outer shell: not a
    droplet. With 0.004 it is 27 of them.
    """
    radii = np.arange(0.5, 14.0, 1.0)
    rho = 0.5 * (0.17 + 0.02) - 0.5 * (0.17 - 0.02) * np.tanh((radii - 4.0) / 1.2)
    noisy = ps.fit_radial_profile(radii, rho, frame_scatter=np.full(radii.size, 0.08))
    assert noisy.success is False
    assert "not resolved in single frames" in noisy.reason

    quiet = ps.fit_radial_profile(radii, rho, frame_scatter=np.full(radii.size, 0.004))
    assert quiet.success is True
    assert quiet.rho_dense == pytest.approx(0.17, abs=1e-4)
    assert quiet.rho_dilute == pytest.approx(0.02, abs=1e-4)
    assert quiet.radius == pytest.approx(4.0, abs=1e-3)

    with pytest.raises(ValueError, match="frame_scatter"):
        ps.fit_radial_profile(radii, rho, frame_scatter=np.zeros(3))


def test_radial_scatter_is_the_standard_deviation_over_frames():
    """``radial_density_profile_with_scatter`` against its definition.

    A solid 5 x 5 x 5 cube is the largest cluster in every frame and its centre
    of mass is the site (10, 10, 10). The vapour, redrawn every frame and kept
    two sites clear of the cube, sets the scatter. Each shell's density is the
    number of beads at minimum-image distance ``[k, k + 1)`` from that site over
    the number of lattice sites there, both counted here by brute force.
    """
    dims = (20, 20, 20)
    rng = np.random.default_rng(21)
    sites = np.indices(dims).reshape(3, -1).T
    offset = sites - np.array([10, 10, 10])
    offset = offset - 20 * np.round(offset / 20.0)
    distance = np.sqrt((offset**2).sum(axis=1))
    cube = sites[np.abs(offset).max(axis=1) <= 2]
    far = sites[np.abs(offset).max(axis=1) >= 5]
    frames = [
        np.concatenate([cube, far[rng.choice(len(far), 300, replace=False)]])
        for _ in range(8)
    ]
    traj = _make(frames, ["A"] * len(frames[0]), dims)

    edges = np.arange(0.0, 11.0, 1.0)
    shell_sites = np.histogram(distance, bins=edges)[0]
    per_frame = []
    for frame in frames:
        d = frame - np.array([10, 10, 10])
        d = d - 20 * np.round(d / 20.0)
        per_frame.append(
            np.histogram(np.sqrt((d**2).sum(axis=1)), bins=edges)[0] / shell_sites
        )
    per_frame = np.array(per_frame)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        radii, rho, counts, scatter = ps.radial_density_profile_with_scatter(traj)
    assert radii == pytest.approx(0.5 * (edges[:-1] + edges[1:]))
    assert counts == pytest.approx(shell_sites)
    assert rho == pytest.approx(per_frame.mean(axis=0), abs=1e-12)
    assert scatter == pytest.approx(per_frame.std(axis=0, ddof=1), abs=1e-12)
    assert scatter[:2] == pytest.approx(0.0, abs=1e-12)  # the cube's core never changes
    assert scatter[6:].min() > 0.0  # the vapour does


# ---------------------------------------------------------------------------
# A2-2: a one-phase hardwall profile is neither translated nor split
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_chains, length", [(12, 40), (40, 30), (75, 30)])
def test_hardwall_one_phase_profile_is_the_plain_mean(n_chains, length):
    """Volume fractions 0.08, 0.21 and 0.39: the profile is the per-plane mean.

    The baseline translated each frame by the centroid of its noise (so the
    depleted wall planes were smeared into the bulk) and, at the higher
    densities, split the frames into "wall" and "free" and left some out with a
    warning about a condensate touching a wall.
    """
    dims = (12, 12, 40)
    traj = _athermal_traj(3, dims, n_chains, length, hardwall=True, n_frames=20)
    pos = traj.positions
    plain = (
        np.mean(
            [np.bincount(pos[f][:, 2], minlength=40) for f in range(traj.n_frames)],
            axis=0,
        )
        / 144.0
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        z, rho = ps.slab_density_profile(traj, axis=2)
    assert not [w for w in caught if "slab_density_profile" in str(w.message)]
    assert not np.isnan(rho).any()
    assert rho == pytest.approx(plain, abs=1e-12)


def test_hardwall_frame_classifier_needs_a_dilute_region():
    """A block missing only its wall planes is a full box, not a slab."""
    full = np.full(40, 30.0)
    full[[0, 39]] = 9.0  # depleted wall planes, below half the peak
    assert ps._hardwall_slab_frame_kind(full) == ("free", False)
    full[1] = full[38] = 12.0  # two depleted planes per wall
    assert ps._hardwall_slab_frame_kind(full) == ("free", False)

    slab = np.full(40, 2.0)
    slab[15:25] = 30.0  # a free slab, 15 dilute planes each side
    assert ps._hardwall_slab_frame_kind(slab) == ("free", True)
    film = np.full(40, 2.0)
    film[:10] = 30.0  # a film on the low wall, 30 dilute planes
    assert ps._hardwall_slab_frame_kind(film) == ("wall", True)
    near = np.full(40, 2.0)
    near[4:36] = 30.0  # four dilute planes each side: still a slab
    assert ps._hardwall_slab_frame_kind(near) == ("free", True)


# ---------------------------------------------------------------------------
# A2-3: the slab normal is the axis the condensate does not span
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "dims, axis, geometry",
    [
        ((20, 20, 8), 1, "auto"),  # normal along one of two equal long axes
        ((20, 20, 8), 0, "auto"),
        ((8, 20, 20), 2, "auto"),
        ((8, 20, 20), 1, "auto"),
        ((14, 14, 14), 2, "slab"),  # cubic box: no longest axis at all
        ((14, 14, 14), 1, "slab"),
        ((30, 8, 8), 0, "auto"),  # the case that always worked
    ],
)
def test_slab_normal_is_found_from_the_condensate(dims, axis, geometry):
    """A six-plane slab at 0.9 / 0.02, whichever way it lies in the box."""
    lo = dims[axis] // 2 - 3
    traj = _slab_traj(
        2, dims, axis, [range(lo, lo + 6)] * 4, rho_dense=0.9, rho_dilute=0.02
    )
    assert ps.slab_normal(traj) == axis
    res = _quiet_analyze(traj, geometry=geometry)
    assert res.geometry == "slab"
    assert res.slab_axis == axis
    assert res.profile[0].size == dims[axis]
    assert res.binodal.success is True, res.binodal
    assert res.rho_dense == pytest.approx(0.9, abs=0.03)
    assert res.rho_dilute == pytest.approx(0.02, abs=0.01)
    assert res.is_phase_separated is True


def test_analyze_axis_argument_overrides_the_choice():
    """``axis=`` is honoured, negative values included, and range-checked."""
    dims = (8, 20, 20)
    traj = _slab_traj(2, dims, 2, [range(7, 13)] * 4, rho_dense=0.9, rho_dilute=0.02)
    along_y = _quiet_analyze(traj, axis=1)  # an in-plane direction: flat
    assert along_y.slab_axis == 1
    assert along_y.is_phase_separated is False
    assert _quiet_analyze(traj, axis=-1).slab_axis == 2
    assert _quiet_analyze(traj, axis=-1).is_phase_separated is True
    assert _quiet_analyze(traj, geometry="sphere").slab_axis is None
    for bad in (3, -4, 1.0, True):
        with pytest.raises(ValueError, match="axis"):
            ps.analyze(traj, axis=bad)


def test_slab_normal_falls_back_to_the_longest_axis():
    """With no slab the clusters do not say, and the longest box axis is used."""
    traj = _athermal_traj(1, (10, 12, 30), 4, 20, hardwall=False, n_frames=4)
    assert ps.slab_normal(traj) == 2
    assert ps._vote_slab_normal([(0, 1), (0, 1), (), None], (9, 9, 9)) == 2
    assert (
        ps._vote_slab_normal([(0, 1), (), (), None], (9, 12, 9)) == 1
    )  # 1 of 3: no majority
    assert ps._vote_slab_normal([(0, 1, 2)] * 3, (9, 9, 30)) == 2  # a network


# ---------------------------------------------------------------------------
# A2-9: negative axis
# ---------------------------------------------------------------------------


def test_negative_axis_is_normalised_in_the_slab_profile():
    """``axis=-1`` is the last axis of the system, in 3D and in 2D."""
    traj = _slab_traj(4, (8, 8, 40), 2, [range(15, 25)] * 3)
    z, rho = ps.slab_density_profile(traj, axis=2)
    z_neg, rho_neg = ps.slab_density_profile(traj, axis=-1)
    assert np.array_equal(rho, rho_neg)
    # every bead is counted once: the profile sums to beads / sites per plane
    assert rho_neg.sum() * 64 == pytest.approx(traj.n_beads)
    assert np.array_equal(
        ps.slab_density_profile(traj, axis=-3)[1],
        ps.slab_density_profile(traj, axis=0)[1],
    )

    flat = _slab_traj(
        4, (20, 60), 1, [range(25, 35)] * 3, rho_dense=0.9, rho_dilute=0.05
    )
    y, stripe = ps.slab_density_profile(flat, axis=-1)
    assert stripe.size == 60
    assert stripe[25:35] == pytest.approx(0.9)
    assert np.array_equal(stripe, ps.slab_density_profile(flat, axis=1)[1])

    for bad in (3, -4, 2.0, "z"):
        with pytest.raises(ValueError, match="axis"):
            ps.slab_density_profile(traj, axis=bad)
    with pytest.raises(ValueError, match="axis"):
        ps.slab_density_profile(flat, axis=2)  # a 2D system has no axis 2
    with pytest.raises(ValueError, match="axis"):
        ps.slab_density_profile(flat, axis=-3)


# ---------------------------------------------------------------------------
# A2-10: thin films and films on both walls
# ---------------------------------------------------------------------------


def test_one_plane_slabs_are_rejected_with_a_slab_reason():
    """One plane cannot give a density: free, at the low wall or at the high wall."""
    dims = (8, 8, 40)
    for hardwall, planes in ((False, [20]), (True, [20]), (True, [0]), (True, [39])):
        traj = _slab_traj(6, dims, 2, [planes] * 20, hardwall=hardwall)
        z, rho = ps.slab_density_profile(traj, axis=2)
        fit = ps.fit_slab_profile(z, rho, hardwall=hardwall)
        assert fit.success is False, (hardwall, planes, fit)
        assert "too thin" in fit.reason
        assert "droplet" not in fit.reason


def test_two_plane_films_get_the_same_verdict_at_either_wall():
    """A film and its mirror image are the same film."""
    dims = (8, 8, 40)
    fits = []
    for planes in ([0, 1], [38, 39]):
        traj = _slab_traj(6, dims, 2, [planes] * 20, hardwall=True)
        z, rho = ps.slab_density_profile(traj, axis=2)
        fits.append(ps.fit_slab_profile(z, rho, hardwall=True))
    low, high = fits
    assert low.success is True and high.success is True
    for fit in fits:
        assert fit.rho_dense == pytest.approx(0.75, abs=0.02)
        assert fit.rho_dilute == pytest.approx(0.0625, abs=0.01)
        assert 2.0 * fit.half_width == pytest.approx(2.0, abs=0.5)


def test_condensate_wetting_both_walls_is_fitted():
    """Two five-plane films with the dilute phase between them."""
    dims = (8, 8, 40)
    planes = list(range(5)) + list(range(35, 40))
    traj = _slab_traj(7, dims, 2, [planes] * 20, hardwall=True)
    res = _quiet_analyze(traj)
    fit = res.binodal
    assert fit.success is True, fit
    assert fit.rho_dense == pytest.approx(0.75, abs=0.02)
    assert fit.rho_dilute == pytest.approx(0.0625, abs=0.01)
    assert fit.half_width == pytest.approx(2.5, abs=0.3)  # half of one film's 5 planes
    assert res.is_phase_separated is True

    # films of unequal thickness: 4 and 8 planes, mean 6, half of it 3
    uneven = list(range(4)) + list(range(32, 40))
    traj = _slab_traj(7, dims, 2, [uneven] * 20, hardwall=True)
    z, rho = ps.slab_density_profile(traj, axis=2)
    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success is True, fit
    assert fit.half_width == pytest.approx(3.0, abs=0.3)


def test_slab_thinner_than_its_interfaces_is_rejected():
    """An exact tanh hump whose half-width is 0.8 of its interface width.

    The profile at its centre reaches tanh(0.8) = 66% of the way to the dense
    asymptote: more than the half that the generic check asks for, so the fit
    used to succeed with ``rho_dense = 0.65`` where the profile never exceeds
    0.45. A slab with half-width 1.5 interface widths (91%) is accepted.
    """
    z = np.arange(40, dtype=float)
    rho = 0.05 + 0.5 * 0.6 * (np.tanh((z - 16.8) / 4.0) - np.tanh((z - 23.2) / 4.0))
    assert rho.max() < 0.46
    fit = ps.fit_slab_profile(z, rho, hardwall=False)
    assert fit.success is False
    assert "thinner than its interfaces" in fit.reason
    assert rho.min() <= fit.rho_dilute < fit.rho_dense <= rho.max()

    thick = 0.05 + 0.5 * 0.6 * (np.tanh((z - 17.0) / 2.0) - np.tanh((z - 23.0) / 2.0))
    fit = ps.fit_slab_profile(z, thick, hardwall=False)
    assert fit.success is True
    assert fit.rho_dense == pytest.approx(0.65, abs=1e-3)
    assert fit.half_width == pytest.approx(3.0, abs=1e-2)


# ---------------------------------------------------------------------------
# A2-11: the radial centre rounds halves one way
# ---------------------------------------------------------------------------


def test_radial_profile_is_invariant_under_a_unit_translation():
    """A cluster whose centre of mass sits on a half-integer along x.

    A 5 x 5 plate in the plane ``x = c`` and a 25-bead rod in the plane
    ``x = c + 1`` have their centre of mass at ``x = c + 0.5`` exactly. The
    profile's origin is the lattice site nearest to it, and a half has to go one
    way. ``np.round`` sends halves to the even neighbour: the origin was the
    plate for even ``c`` and the rod for odd ``c``, so moving the whole system
    by one site changed the profile (the plate's corners are 2.83 from the
    plate's centre and 3.0 from the rod's, a different shell).
    """
    dims = (32, 32, 32)
    plate = [(0, y, z) for y in range(5) for z in range(5)]
    rod = [(1, 2, z) for z in range(-10, 15)]
    cluster = np.array(plate + rod)
    assert cluster.mean(axis=0).tolist() == [0.5, 2.0, 2.0]
    profiles = []
    for corner in (10, 11, 12, 13):
        traj = _make([cluster + np.array([corner, 10, 12])], ["A"] * len(cluster), dims)
        radii, rho = ps.radial_density_profile(traj)
        profiles.append(rho)
    for other in profiles[1:]:
        assert other == pytest.approx(profiles[0], abs=1e-12)
    assert profiles[0][0] == pytest.approx(1.0)  # the origin site is occupied


# ---------------------------------------------------------------------------
# A2-4: hardwall slab surface tension on cosine modes
# ---------------------------------------------------------------------------


def _neumann_field(rng, n_x, n_y, gamma, kT=1.0):
    """A height field drawn from the capillary spectrum of a walled interface.

    The energy is ``(gamma / 2) sum (h_i - h_j)^2`` over nearest-neighbour
    columns with free ends. Its normal modes are ``cos(pi m (x + 1/2) / n_x)
    cos(pi n (y + 1/2) / n_y)`` with eigenvalue ``lam = (2 - 2 cos(pi m / n_x))
    + (2 - 2 cos(pi n / n_y))``, and equipartition gives each amplitude the
    variance ``kT / (gamma lam |phi|^2)``.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    n_x, n_y : int
        Cross-section in lattice units.
    gamma : float
        Surface tension.
    kT : float, optional
        Temperature (default ``1.0``).

    Returns
    -------
    numpy.ndarray
        ``(n_x, n_y)`` float64 heights, zero mean.
    """
    x = np.arange(n_x) + 0.5
    y = np.arange(n_y) + 0.5
    h = np.zeros((n_x, n_y))
    for m in range(n_x):
        for n in range(n_y):
            if m == 0 and n == 0:
                continue
            phi = np.outer(np.cos(np.pi * m * x / n_x), np.cos(np.pi * n * y / n_y))
            lam = (2 - 2 * np.cos(np.pi * m / n_x)) + (2 - 2 * np.cos(np.pi * n / n_y))
            h += (
                rng.standard_normal()
                * np.sqrt(kT / (gamma * lam * (phi**2).sum()))
                * phi
            )
    return h


def _film_traj(heights, dims, hardwall, base=0):
    """Solid columns from ``base`` up to a given top height, one bead per chain.

    Parameters
    ----------
    heights : list of numpy.ndarray
        Per frame, the ``(n_x, n_y)`` integer height of the top bead of every
        column.
    dims : tuple of int
        Box extent in lattice units.
    hardwall : bool
        Whether the box has hard walls.
    base : int, optional
        Height of the lowest bead of every column (default ``0``, on the wall).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory. Frames are brought to a common bead count by removing
        beads from the plane just above ``base``, which is never a top bead.
    """
    frames = []
    for top in heights:
        frames.append(
            np.array(
                [
                    (x, y, z)
                    for x in range(dims[0])
                    for y in range(dims[1])
                    for z in range(base, int(top[x, y]) + 1)
                ]
            )
        )
    n = min(len(f) for f in frames)
    rng = np.random.default_rng(0)
    out = []
    for f in frames:
        extra = len(f) - n
        if extra:
            interior = np.nonzero(f[:, 2] == base + 1)[0]
            f = np.delete(f, rng.choice(interior, extra, replace=False), axis=0)
        out.append(f)
    return _make(out, ["A"] * n, dims, hardwall=hardwall, temperature=1.0)


@pytest.mark.parametrize("gamma", [0.5, 1.0])
def test_hardwall_slab_surface_tension_recovers_a_known_gamma(gamma):
    """A film on the low wall whose free face carries a known capillary spectrum.

    The periodic transform applied to this walled interface read 0.55 gamma.
    """
    n = 8
    rng = np.random.default_rng(17)
    heights = [
        np.clip(np.rint(12 + _neumann_field(rng, n, n, gamma)), 4, 22).astype(int)
        for _ in range(160)
    ]
    traj = _film_traj(heights, (n, n, 24), hardwall=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = st.slab_surface_tension(traj, axis=2)
    assert res.n_modes == 8
    assert res.gamma == pytest.approx(gamma, rel=0.12)


def test_hardwall_slab_spectrum_is_the_cosine_projection():
    """One frame, checked mode by mode against the definition.

    ``gamma = kT / mean(lam a^2 |phi|^2)`` over the eight softest cosine modes,
    with ``a = sum(dh phi) / |phi|^2`` the amplitude of the height fluctuation
    on the mode.
    """
    n_x, n_y = 6, 7
    rng = np.random.default_rng(8)
    top = rng.integers(8, 13, size=(n_x, n_y))
    traj = _film_traj([top], (n_x, n_y, 20), hardwall=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = st.slab_surface_tension(traj, axis=2, temperature=1.0)

    dh = top - top.mean()
    x = np.arange(n_x) + 0.5
    y = np.arange(n_y) + 0.5
    modes = []
    for m in range(n_x):
        for k in range(n_y):
            if m == 0 and k == 0:
                continue
            phi = np.outer(np.cos(np.pi * m * x / n_x), np.cos(np.pi * k * y / n_y))
            lam = (2 - 2 * np.cos(np.pi * m / n_x)) + (2 - 2 * np.cos(np.pi * k / n_y))
            norm = (phi**2).sum()
            amplitude = (dh * phi).sum() / norm
            modes.append((lam, lam * amplitude**2 * norm))
    modes.sort(key=lambda pair: pair[0])
    expected = 1.0 / np.mean([energy for _lam, energy in modes[:8]])
    assert res.gamma == pytest.approx(expected, rel=1e-9)
    assert res.spectrum[0] == pytest.approx(
        np.sqrt([lam for lam, _e in modes[:8]]), rel=1e-12
    )


def test_slab_surface_tension_axis_is_normalised_and_chosen_from_the_slab():
    """``axis=-1`` is axis 2; out of range is refused; a cubic box needs no hint."""
    rng = np.random.default_rng(4)
    dims = (10, 10, 10)
    heights = [rng.integers(5, 8, size=(10, 10)) for _ in range(3)]
    traj = _film_traj(heights, dims, hardwall=False, base=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        explicit = st.slab_surface_tension(traj, axis=2)
        negative = st.slab_surface_tension(traj, axis=-1)
        default = st.slab_surface_tension(traj)  # argmax(dims) would be axis 0
    assert np.isfinite(explicit.gamma)
    assert negative.gamma == explicit.gamma
    assert default.gamma == explicit.gamma
    for bad in (3, 5, -4, 0.0):
        with pytest.raises(ValueError, match="axis"):
            st.slab_surface_tension(traj, axis=bad)


# ---------------------------------------------------------------------------
# T7: mutants that survived the existing tests
# ---------------------------------------------------------------------------


def test_slab_surface_tension_is_invariant_under_translation_along_the_normal():
    """The circular centre must find the slab wherever it sits (T7, sine sign).

    With the sign of the sine sum flipped the centre of a slab at L/4 comes out
    at 3L/4, the "re-centred" slab straddles the boundary and the heights are
    those of two half slabs.
    """
    rng = np.random.default_rng(6)
    dims = (8, 8, 40)
    low = [rng.integers(15, 18, size=(8, 8)) for _ in range(4)]
    high = [rng.integers(23, 26, size=(8, 8)) for _ in range(4)]
    frames = [
        np.array(
            [
                (x, y, z)
                for x in range(8)
                for y in range(8)
                for z in range(lo[x, y], hi[x, y] + 1)
            ]
        )
        for lo, hi in zip(low, high)
    ]
    n = min(len(f) for f in frames)
    keep = []
    for f in frames:
        middle = np.nonzero(f[:, 2] == 20)[0]
        keep.append(np.delete(f, middle[: len(f) - n], axis=0))
    gammas = []
    for shift in (0, 10, 20, 30, 37):  # slab centre at L/2, 3L/4, 0, L/4, ...
        moved = [
            np.column_stack([f[:, 0], f[:, 1], (f[:, 2] + shift) % 40]) for f in keep
        ]
        traj = _make(moved, ["A"] * n, dims)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gammas.append(st.slab_surface_tension(traj, axis=2).gamma)
    assert np.isfinite(gammas[0])
    assert gammas[1:] == pytest.approx([gammas[0]] * 4, rel=1e-9)


def test_slab_surface_tension_counts_the_frames_it_skips():
    """Three droplet frames are three skipped frames (T7, ``n_not_slab``)."""
    dims = (12, 12, 30)
    ball = np.array(
        [(x, y, z) for x in range(4, 8) for y in range(4, 8) for z in range(10, 14)]
    )
    traj = _make([ball] * 3, ["A"] * len(ball), dims)
    with pytest.warns(
        UserWarning,
        match="slab_surface_tension: 3 of 3 frames with a cluster "
        "were skipped because the largest cluster does not "
        "span both in-plane",
    ):
        res = st.slab_surface_tension(traj, axis=2)
    assert np.isnan(res.gamma) and res.n_modes == 0


def test_flat_slab_returns_infinite_gamma_and_the_spectrum_it_measured():
    """A perfectly flat slab: ``gamma = inf`` and zero power on the softest modes (T7)."""
    dims = (6, 6, 20)
    slab = np.array(
        [(x, y, z) for x in range(6) for y in range(6) for z in range(8, 12)]
    )
    traj = _make([slab] * 2, ["A"] * len(slab), dims)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = st.slab_surface_tension(traj, axis=2, n_modes=4)
    assert res.gamma == float("inf") and res.n_modes == 0
    q, power = res.spectrum
    # the four softest independent wavevectors of a 6 x 6 lattice: (1,0), (0,1)
    # with 2 - 2 cos(2 pi / 6) = 1, then (1,1) and (1,-1) with 2
    assert np.sort(q) == pytest.approx(np.sqrt([1.0, 1.0, 2.0, 2.0]), abs=1e-12)
    assert len(power) == 4 and np.all(power == 0.0)


def test_periodic_slab_fit_pins_the_centre_and_hardwall_fit_does_not():
    """An off-centre slab (T7, pinned versus fitted centre).

    ``hardwall=False`` describes a profile that ``slab_density_profile`` has
    re-centred, so the model's centre is the middle of the window and a slab
    sitting elsewhere cannot be fitted. ``hardwall=True`` fits the centre.
    """
    z = np.arange(40, dtype=float)
    centred = np.where((z >= 15) & (z <= 24), 0.8, 0.05)
    shifted = np.where((z >= 5) & (z <= 14), 0.8, 0.05)

    pinned = ps.fit_slab_profile(z, centred, hardwall=False)
    assert pinned.success is True
    assert pinned.half_width == pytest.approx(5.0, abs=0.3)

    free = ps.fit_slab_profile(z, shifted, hardwall=True)
    assert free.success is True
    assert free.half_width == pytest.approx(5.0, abs=0.3)
    assert free.rho_dense == pytest.approx(0.8, abs=0.01)

    wrong = ps.fit_slab_profile(z, shifted, hardwall=False)
    assert not (
        wrong.success
        and wrong.half_width == pytest.approx(5.0, abs=0.3)
        and wrong.rho_dense == pytest.approx(0.8, abs=0.01)
    )


# ---------------------------------------------------------------------------
# second round (review REV_F8a): chains as the unit, thresholds pinned
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "dense_chain_length, separated", [(8, True), (96, True), (160, False)]
)
def test_slab_must_hold_significantly_more_chains_than_a_uniform_solution(
    dense_chain_length, separated
):
    """The same 600 beads at the same sites, cut into fewer and longer chains.

    Ten planes of 40 are 75% full (six rows of eight sites each), and the 480
    beads in them form one snake path cut into chains of 8, 96 or 160 beads; 15
    more chains of 8 are the vapour. With ``Lw`` the bead-weighted mean chain
    length (``0.8 l + 1.6``) the system holds ``600 / Lw`` chains, the slab
    ``480 / Lw``, and a uniform solution would put a quarter of them there:
    ``sigma = (330 / Lw) / sqrt(600 / Lw * 0.25 * 0.75) = 31.11 / sqrt(Lw)``,
    i.e. 11.0, 3.51 and 2.73. Three are needed: a lump of three chains is not a
    phase, however clean its averaged profile.
    """
    dims = (8, 8, 40)
    snake = []
    for y in range(6):
        xs = range(8) if y % 2 == 0 else range(7, -1, -1)
        snake.extend((x, y) for x in xs)
    dense = []
    for k, z in enumerate(range(15, 25)):
        path = snake if k % 2 == 0 else snake[::-1]
        dense.extend((x, y, z) for x, y in path)
    vapour_planes = [z for z in range(40) if not 15 <= z < 25]
    vapour = [(x, 3, z) for z in vapour_planes[::2] for x in range(8)]
    frames = [np.array(dense + vapour)] * 4
    sequences = ["A" * dense_chain_length] * (480 // dense_chain_length) + [
        "A" * 8
    ] * 15
    res = _quiet_analyze(_make(frames, sequences, dims))
    mean_length = 0.8 * dense_chain_length + 1.6
    assert res.chain_excess_sigma == pytest.approx(
        31.11 / np.sqrt(mean_length), rel=1e-3
    )
    assert res.is_phase_separated is separated
    assert res.binodal.success is separated
    if separated:
        assert res.rho_dense == pytest.approx(0.75, abs=0.01)
        assert res.rho_dilute == pytest.approx(0.0625, abs=0.01)
    else:
        assert (
            "chains" in res.binodal.reason
            and "cannot be told apart" in res.binodal.reason
        )


def test_chain_excess_sigma_against_the_binomial_formula():
    """``_chain_excess_sigma`` on a hand-made profile: 6 dense planes of 30."""
    density = np.full(30, 0.02)
    density[12:18] = 0.80
    lengths = np.full(40, 25)  # 40 chains of 25 beads = 1000 beads
    cross = 1000.0 / density.sum()  # so that the profile holds exactly 1000 beads
    sigma, n_dense, n_expected = ps._chain_excess_sigma(
        density, cross, 0.80, 0.02, lengths
    )
    in_dense = 6 * 0.80 * cross / 25.0
    assert n_dense == pytest.approx(in_dense)
    assert n_expected == pytest.approx(40 * 0.2)
    assert sigma == pytest.approx((in_dense - 8.0) / np.sqrt(40 * 0.2 * 0.8))
    # no dense region, or nothing but dense region: no excess to speak of
    assert ps._chain_excess_sigma(np.full(30, 0.5), cross, 0.8, 0.02, lengths)[0] == 0.0
    assert ps._chain_excess_sigma(np.full(30, 0.9), cross, 0.8, 0.02, lengths)[0] == 0.0


def test_dilute_plateau_width_and_tolerance_are_pinned():
    """Exact tanh slabs that fill most of a 24-plane box.

    Interfaces at 2.5 and 21.5 with width 0.3: planes 0, 1, 2, 22, 23 are within
    10% of the gap of ``rho_dilute``, five planes. With width 1.0 only planes 0,
    1 and 23 are (plane 2 and plane 22 sit half a site from an interface, 27% of
    the gap up). Interfaces at 1.5 and 22.5 with width 0.6: plane 0 alone.
    """
    z = np.arange(24, dtype=float)

    def slab(half_width, width):
        return 0.02 + 0.5 * 0.8 * (
            np.tanh((z - (12.0 - half_width)) / width)
            - np.tanh((z - (12.0 + half_width)) / width)
        )

    sharp = slab(9.5, 0.3)
    assert ps.fit_slab_profile(z, sharp, hardwall=False).success is True
    assert (
        ps.fit_slab_profile(z, sharp, hardwall=False, min_dilute_planes=5).success
        is True
    )
    six = ps.fit_slab_profile(z, sharp, hardwall=False, min_dilute_planes=6)
    assert six.success is False and "no dilute plateau" in six.reason

    broad = slab(9.5, 1.0)
    assert (
        ps.fit_slab_profile(z, broad, hardwall=False, min_dilute_planes=3).success
        is True
    )
    assert (
        ps.fit_slab_profile(z, broad, hardwall=False, min_dilute_planes=4).success
        is False
    )

    one = slab(10.5, 0.6)  # plane 0 only; the default asks for three
    assert ps.fit_slab_profile(z, one, hardwall=False).success is False
    assert (
        ps.fit_slab_profile(z, one, hardwall=False, min_dilute_planes=2).success
        is False
    )
    assert (
        ps.fit_slab_profile(z, one, hardwall=False, min_dilute_planes=1).success is True
    )


def test_frame_scatter_reaches_the_wetting_and_both_wall_fits():
    """The optional single-plane test is applied on every slab model."""
    z = np.arange(40, dtype=float)
    film = np.where(z < 10, 0.75, 0.0625)
    both = np.where((z < 6) | (z > 33), 0.75, 0.0625)
    for profile in (film, both):
        quiet = ps.fit_slab_profile(
            z, profile, hardwall=True, frame_scatter=np.full(40, 0.001)
        )
        assert quiet.success is True
        noisy = ps.fit_slab_profile(
            z, profile, hardwall=True, frame_scatter=np.full(40, 0.5)
        )
        assert noisy.success is False and "cannot be told apart" in noisy.reason


def test_radial_fit_must_reach_its_dilute_asymptote():
    """A profile that ends on the interface has no dilute density to report.

    An exact tanh droplet of radius 5 and width 2, seen out to r = 6.5 only: the
    last shell is still 18% of the gap above the dilute asymptote.
    """
    radii = np.arange(0.5, 7.0, 1.0)
    rho = 0.5 * (0.8 + 0.02) - 0.5 * (0.8 - 0.02) * np.tanh((radii - 5.0) / 2.0)
    short = ps.fit_radial_profile(radii, rho)
    assert short.success is False
    assert "extrapolation" in short.reason
    assert rho.min() <= short.rho_dilute <= short.rho_dense <= rho.max()

    radii = np.arange(0.5, 14.0, 1.0)
    rho = 0.5 * (0.8 + 0.02) - 0.5 * (0.8 - 0.02) * np.tanh((radii - 5.0) / 2.0)
    full = ps.fit_radial_profile(radii, rho)
    assert full.success is True
    assert full.rho_dilute == pytest.approx(0.02, abs=1e-4)


def test_hardwall_slab_surface_tension_warns_that_it_is_not_validated_on_real_runs():
    """Under HARDWALL the estimator says what it is worth."""
    rng = np.random.default_rng(2)
    heights = [rng.integers(8, 11, size=(6, 6)) for _ in range(2)]
    traj = _film_traj(heights, (6, 6, 20), hardwall=True)
    with pytest.warns(UserWarning, match="not reliable on real hardwall slabs"):
        st.slab_surface_tension(traj, axis=2)
