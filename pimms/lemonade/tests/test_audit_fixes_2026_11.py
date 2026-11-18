"""
Regression tests for the November 2026 audit of lemonade.

F1: ``geometry='auto'`` chose slab or droplet from the box shape alone, so a
droplet in an elongated box came back as a "successful" slab with a dense
density about four times too low, and a slab in a cubic box as not phase
separated (surface tension ``nan``). F2: a second slab on the dilute plateau
was averaged into ``rho_dilute`` with ``success=True``. F3: a periodic
trajectory loaded under hard walls was not detected. F4: non-integer chain
indices, a bool temperature and bool frame selections were coerced silently.
F6: the both-walls fit returned ``nan`` rather than its starting guess when the
optimiser failed. F7: upper-casing a keyfile sequence changed its length.

A review of those fixes found three defects they had added, and this file
pins their corrections too. D1: the cluster vote sent a droplet pressed
against the four side walls of a hardwall box (which "spans" both axes by
touching their walls), or a layer full of holes, in a box that is not
elongated, to the slab fit, which reported about half their density. D2:
``traj[True:]`` still dropped frame 0. D3: between hard walls one dense plane
on each wall counted as two neighbouring planes. And frames already reported
as wrapped through the walls were reported a second time by the F3 check.

Every trajectory is built by hand with numpy (or run with PIMMS), with an
answer that follows from how it was built. Each test failed on the tree before
its fix, except the ones that guard a fix against its own side effect, which
pass on the tree before the first fixes and fail without the correction:
``test_auto_does_not_call_a_lone_coil_in_an_elongated_box_a_droplet`` (the
false positive a simpler form of the F1 fix would have added),
``test_auto_does_not_profile_a_partly_filled_layer_as_a_slab_in_a_cubic_box``
(D1), ``test_one_dense_plane_on_each_wall_is_not_a_second_slab`` (D3),
``test_frames_wrapped_through_the_walls_are_reported_once`` and
``test_hardwall_trajectory_does_not_trip_the_face_check`` (false alarms of the
F3 check).
"""

import warnings

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import _load
from pimms.lemonade import phase_separation as ps
from pimms.lemonade import surface_tension as st
from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology
from pimms.lemonade.tests.conftest import _make_trajectory
from pimms.lemonade.tests.test_deep_audit_phase import (_athermal_traj, _make,
                                                        _quiet_analyze, _slab_traj)
from pimms.lemonade.tests.test_review5_fixes import resized_equilibration_run  # noqa: F401 (fixture)

# Densities of the synthetic condensates: occupied-site fractions chosen well
# inside the two phases of a lattice fluid, so that the truth is unambiguous.
_RHO_DENSE = 0.9
_RHO_DILUTE_DROPLET = 0.01
_RHO_DILUTE_SLAB = 0.02

# The bead types PIMMS makes of the keyfile chain 'aßıſéε' under
# CASE_INSENSITIVE_CHAINS: 'a', 'é' and 'ε' upper-case to one character
# each; 'ß' (upper case 'SS'), dotless 'ı' (upper case 'I') and long
# 'ſ' (upper case 'S') are kept as written.
_UNICODE_BEADS = ("A", "ß", "ı", "ſ", "É", "Ε")

# A parameter file in which every one of those bead types interacts weakly
# with every other one, so that the run is an ordinary one-phase solution.
_UNICODE_PARAMS = (
    "".join("ANGLE_PENALTY\t%s\t0\t0\t0\n" % b for b in _UNICODE_BEADS)
    + "".join("%s\t%s\t-1\n" % (a, b) for i, a in enumerate(_UNICODE_BEADS)
              for b in _UNICODE_BEADS[i:])
    + "".join("%s\t0\t0\n" % b for b in _UNICODE_BEADS))


def _droplet_traj(seed, dims, radius, n_frames=12, hardwall=False):
    """A droplet of monomers at a random centre each frame, in a dilute background.

    Sites within ``radius`` of the centre are occupied at ``_RHO_DENSE``, the
    rest at ``_RHO_DILUTE_DROPLET``, the same number of beads in every frame.
    Under a hardwall the centre is kept a site clear of the walls.

    Parameters
    ----------
    seed : int
        Seed of the generator.
    dims : tuple of int
        Box extent in lattice units (2D or 3D).
    radius : int
        Droplet radius in lattice units.
    n_frames : int, optional
        Number of frames (default ``12``).
    hardwall : bool, optional
        Whether the box has hard walls (default ``False``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory, every bead a one-bead chain.
    """
    rng = np.random.default_rng(seed)
    box = np.asarray(dims)
    grid = np.stack(np.meshgrid(*[np.arange(d) for d in dims], indexing="ij"),
                    axis=-1).reshape(-1, len(dims))
    frames = []
    n_inside = n_outside = None
    for _ in range(n_frames):
        if hardwall:
            centre = rng.integers(radius + 1, box - radius - 1)
        else:
            centre = rng.integers(0, box)
        delta = grid - centre
        if not hardwall:
            delta = (delta + box // 2) % box - box // 2
        inside = (delta ** 2).sum(axis=1) <= radius * radius
        if n_inside is None:
            n_inside = int(round(_RHO_DENSE * inside.sum()))
            n_outside = int(round(_RHO_DILUTE_DROPLET * (~inside).sum()))
        pick = np.concatenate([rng.choice(np.flatnonzero(inside), n_inside, replace=False),
                               rng.choice(np.flatnonzero(~inside), n_outside, replace=False)])
        frames.append(grid[pick])
    return _make(frames, ["A"] * (n_inside + n_outside), dims, hardwall=hardwall)


def _wandering_planes(seed, length, thickness, n_frames=12):
    """Per frame, the planes of a slab at a random position along a periodic axis.

    Parameters
    ----------
    seed : int
        Seed of the generator.
    length : int
        Box length along the slab normal.
    thickness : int
        Planes in the slab.
    n_frames : int, optional
        Number of frames (default ``12``).

    Returns
    -------
    list of list of int
        The dense planes of each frame.
    """
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, length, size=n_frames)
    return [[(int(z0) + k) % length for k in range(thickness)] for z0 in starts]


# ---------------------------------------------------------------------------
# F1: geometry='auto' follows the condensate
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("dims, radius, hardwall", [
    ((16, 16, 40), 5, False),
    ((16, 16, 40), 5, True),
    ((30, 80), 8, False),
])
def test_auto_profiles_a_droplet_in_an_elongated_box_as_a_droplet(dims, radius, hardwall):
    """The box used to decide: slab, and a dense density of 0.27 (3D) or 0.44 (2D)."""
    traj = _droplet_traj(3, dims, radius, hardwall=hardwall)
    res = _quiet_analyze(traj)
    assert res.rho_dense == pytest.approx(_RHO_DENSE, abs=0.03), res.binodal
    assert res.rho_dilute == pytest.approx(_RHO_DILUTE_DROPLET, abs=0.005)
    assert res.geometry == "sphere"
    assert res.binodal.success is True
    assert res.is_phase_separated is True


@pytest.mark.parametrize("dims, axis", [((16, 16, 16), 1), ((16, 16, 16), 2), ((32, 32), 0)])
def test_auto_profiles_a_slab_in_a_cubic_box_as_a_slab(dims, axis):
    """The box used to decide: droplet, a spanning cluster, so not phase separated."""
    planes = _wandering_planes(4, dims[axis], 6)
    traj = _slab_traj(4, dims, axis, planes, rho_dense=_RHO_DENSE, rho_dilute=_RHO_DILUTE_SLAB)
    # _slab_traj fills each plane with a whole number of beads
    cross_section = int(np.prod([d for i, d in enumerate(dims) if i != axis]))
    dense = round(_RHO_DENSE * cross_section) / cross_section
    dilute = round(_RHO_DILUTE_SLAB * cross_section) / cross_section
    res = _quiet_analyze(traj)
    assert res.is_phase_separated is True, res.binodal
    assert res.rho_dense == pytest.approx(dense, abs=0.01)
    assert res.rho_dilute == pytest.approx(dilute, abs=0.005)
    assert res.geometry == "slab"
    assert res.slab_axis == axis


def test_surface_tension_auto_follows_the_condensate():
    """A cubic-box slab used to go to the droplet estimator (``nan``), and a droplet
    in an elongated box to the slab estimator (``nan``, every frame skipped)."""
    slab = _slab_traj(5, (16, 16, 16), 1, _wandering_planes(5, 16, 6),
                      rho_dense=_RHO_DENSE, rho_dilute=_RHO_DILUTE_SLAB)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        on_slab = st.surface_tension(slab)
    assert np.isfinite(on_slab.gamma)
    assert on_slab.method == "slab"

    droplet = _droplet_traj(6, (16, 16, 40), 5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        on_droplet = st.surface_tension(droplet)
    assert np.isfinite(on_droplet.gamma)
    assert on_droplet.method == "droplet"


def test_explicit_slab_geometry_refuses_a_compact_droplet():
    """A droplet profiled as a slab is diluted plane by plane; it used to pass
    every check with a dense density of 0.27."""
    traj = _droplet_traj(3, (16, 16, 40), 5)
    res = _quiet_analyze(traj, geometry="slab")
    assert res.binodal.success is False, res.binodal
    assert "not a slab" in res.binodal.reason
    assert res.is_phase_separated is False
    assert res.geometry == "slab"


def test_auto_does_not_call_a_lone_coil_in_an_elongated_box_a_droplet():
    """Two athermal 20-bead chains in a 16 x 16 x 48 hardwall box: one phase.

    The largest cluster is a lone coil, compact in every frame, and the droplet
    fit (which has no chain-number test) takes it for a droplet of density
    0.37. Sending every compact cluster to the droplet analysis would have made
    this a false positive; the slab analysis is asked first in an elongated box.
    This one passed before the fix too.
    """
    traj = _athermal_traj(1, (16, 16, 48), 2, 20, True)
    spans = ps._largest_cluster_spanning_axes(traj)
    assert all(axes == () for axes in spans if axes is not None)
    res = _quiet_analyze(traj)
    assert res.is_phase_separated is False, res.binodal
    assert res.geometry == "slab"


def test_auto_geometry_rule():
    """The shared rule on hand-written spanning axes, and on frames where it needs them.

    A box stand-in suffices wherever the rule does not look at the frames: an
    elongated box, or no decisive slab vote.
    """
    from types import SimpleNamespace

    cube = SimpleNamespace(dimensions=(20, 20, 20), n_dim=3)
    long_box = SimpleNamespace(dimensions=(12, 12, 40), n_dim=3)
    assert ps._auto_geometry(long_box, [(0, 1)] * 3 + [(), None], 2) == ("slab", True)
    assert ps._auto_geometry(long_box, [()] * 3 + [(0,)], 2) == ("sphere", True)
    assert ps._auto_geometry(cube, [()] * 3 + [(0,)], 2) == ("sphere", True)
    # a cylinder (one axis of three), a network, or no majority: the box decides
    for frames in ([(0,)] * 4, [(0, 1, 2)] * 4, [(), (0,), (0, 1), (2,)], [None] * 3):
        assert ps._auto_geometry(cube, frames, 2) == ("sphere", False)
        assert ps._auto_geometry(long_box, frames, 2) == ("slab", False)
    # the same normal must win: slabs in alternating orientations do not
    assert ps._auto_geometry(cube, [(0, 1), (1, 2)] * 2, 2) == ("sphere", False)

    # in a box that is not elongated the vote must be backed by filled planes
    slab = _slab_traj(4, (16, 16, 16), 1, _wandering_planes(4, 16, 6),
                      rho_dense=_RHO_DENSE, rho_dilute=_RHO_DILUTE_SLAB)
    assert ps._auto_geometry(slab, ps._largest_cluster_spanning_axes(slab), 2) == ("slab", True)
    stripe = _slab_traj(4, (32, 32), 0, _wandering_planes(4, 32, 6),
                        rho_dense=_RHO_DENSE, rho_dilute=_RHO_DILUTE_SLAB)
    assert ps._auto_geometry(stripe, ps._largest_cluster_spanning_axes(stripe), 2) == ("slab", True)
    pressed = _wall_droplet_traj(5)
    spans = ps._largest_cluster_spanning_axes(pressed)
    assert ps._decisive_slab_normal(spans, 3) == 2
    assert ps._auto_geometry(pressed, spans, 2) == ("sphere", False)


def _wall_droplet_traj(seed, n_frames=12):
    """A droplet pressed against the four side walls of a 20 x 20 x 26 hardwall box.

    A sphere of radius 10 (occupancy ``_RHO_DENSE``) centred on the box's axis,
    at a random height, in a ``_RHO_DILUTE_DROPLET`` background: it touches both
    walls of x and of y, so by the hardwall definition it spans them, and the
    slab vote picks z. The reviewer's reproducer of D1.

    Parameters
    ----------
    seed : int
        Seed of the generator.
    n_frames : int, optional
        Number of frames (default ``12``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory, every bead a one-bead chain.
    """
    dims = (20, 20, 26)
    radius = 10.0
    rng = np.random.default_rng(seed)
    grid = np.stack(np.meshgrid(*[np.arange(d) for d in dims], indexing="ij"),
                    axis=-1).reshape(-1, 3)
    frames = []
    n_inside = n_outside = None
    for _ in range(n_frames):
        centre = np.array([10, 10, rng.integers(11, 15)])
        inside = np.sqrt(((grid - centre) ** 2).sum(axis=1)) <= radius
        if n_inside is None:
            n_inside = int(round(_RHO_DENSE * inside.sum()))
            n_outside = int(round(_RHO_DILUTE_DROPLET * (~inside).sum()))
        pick = np.concatenate([rng.choice(np.flatnonzero(inside), n_inside, replace=False),
                               rng.choice(np.flatnonzero(~inside), n_outside, replace=False)])
        frames.append(grid[pick])
    return _make(frames, ["A"] * (n_inside + n_outside), dims, hardwall=True)


def _perforated_layer_traj(seed, n_frames=12):
    """A six-plane layer at ``_RHO_DENSE`` with two 8 x 8 holes, in a periodic 20^3 box.

    The layer is connected to its own image across x and y, so it wins the
    slab vote (normal z), but a third of each of its planes is empty: its slab
    profile reads about 0.61, not 0.9. A stand-in for the near-critical
    perforated layer of the review (real run ``J2_box18x18x24_T175``).

    Parameters
    ----------
    seed : int
        Seed of the generator.
    n_frames : int, optional
        Number of frames (default ``12``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        The trajectory, every bead a one-bead chain.
    """
    dims = (20, 20, 20)
    rng = np.random.default_rng(seed)
    frames = []
    n_layer = n_rest = None
    for _ in range(n_frames):
        layer = np.zeros(dims, dtype=bool)
        z0 = int(rng.integers(0, 20))
        layer[:, :, [(z0 + k) % 20 for k in range(6)]] = True
        shift = rng.integers(0, 20, size=2)
        for corner in ((0, 0), (10, 10)):
            xs = [(shift[0] + corner[0] + k) % 20 for k in range(8)]
            ys = [(shift[1] + corner[1] + k) % 20 for k in range(8)]
            layer[np.ix_(xs, ys, np.arange(20))] &= False
        sites, rest = np.argwhere(layer), np.argwhere(~layer)
        if n_layer is None:
            n_layer = int(round(_RHO_DENSE * len(sites)))
            n_rest = int(round(_RHO_DILUTE_SLAB * len(rest)))
        pick = np.concatenate([sites[rng.choice(len(sites), n_layer, replace=False)],
                               rest[rng.choice(len(rest), n_rest, replace=False)]])
        frames.append(pick)
    return _make(frames, ["A"] * (n_layer + n_rest), dims)


@pytest.mark.parametrize("make_traj", [_wall_droplet_traj, _perforated_layer_traj])
def test_auto_does_not_profile_a_partly_filled_layer_as_a_slab_in_a_cubic_box(make_traj):
    """Review D1: the slab vote sent these to the slab fit, which reported success
    with rho_dense 0.69 (droplet pressed against the walls) and 0.61 (layer with
    holes) for a dense phase of 0.9. Before the vote they were refused as spanning
    clusters, and they are again; the surface tension no longer uses the slab
    estimator on them either."""
    traj = make_traj(5)
    res = _quiet_analyze(traj)
    assert not (res.binodal.success and res.rho_dense < 0.8), res.binodal
    assert res.is_phase_separated is False
    assert res.geometry == "sphere"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert st.surface_tension(traj).method == "droplet"


def test_plane_coverage_on_hand_made_planes():
    """Full, holed and disc-shaped planes, periodic and between walls."""
    full = np.ones((10, 10), dtype=bool)
    full[::3, ::3] = False                          # scattered vacancies are not gaps
    assert ps._plane_coverage(full, hardwall=False) == (1.0, 1.0)
    holed = np.ones((10, 10), dtype=bool)
    holed[2:7, 2:7] = False                         # a 5 x 5 hole: its 3 x 3 core stays open
    assert ps._plane_coverage(holed, hardwall=False)[0] == pytest.approx(1.0 - 9 / 100)
    x, y = np.meshgrid(np.arange(20), np.arange(20), indexing="ij")
    disc = (x - 9.5) ** 2 + (y - 9.5) ** 2 <= 10.0 ** 2
    # a disc touching the four walls covers most of the plane, but the edge only
    # near the four points of contact
    plane, wall = ps._plane_coverage(disc, hardwall=True)
    assert plane > ps._MIN_PLANE_COVERAGE
    assert wall < ps._MIN_WALL_COVERAGE


# ---------------------------------------------------------------------------
# F2: a second slab on the dilute plateau
# ---------------------------------------------------------------------------


def test_a_second_slab_is_not_averaged_into_the_dilute_density():
    """Slabs of 8 and 4 planes at 0.9 in 0.02: rho_dilute used to come back 0.13."""
    rng = np.random.default_rng(0)
    planes = []
    for z0 in rng.integers(0, 40, size=12):
        first = [(int(z0) + k) % 40 for k in range(8)]
        second = [(int(z0) + 22 + k) % 40 for k in range(4)]
        planes.append(first + second)
    traj = _slab_traj(0, (12, 12, 40), 2, planes, rho_dense=_RHO_DENSE,
                      rho_dilute=_RHO_DILUTE_SLAB)
    res = _quiet_analyze(traj)
    assert res.binodal.success is False, res.binodal
    assert "second condensate" in res.binodal.reason
    assert res.is_phase_separated is False


def test_second_slab_check_on_a_profile():
    """A clean profile with a second 3-plane slab is refused; scattered dense
    planes in the dilute region (one chain filling half a narrow plane) are not."""
    z = np.arange(40, dtype=float)
    rho = np.full(40, _RHO_DILUTE_SLAB)
    rho[16:24] = _RHO_DENSE
    second = rho.copy()
    second[2:5] = _RHO_DENSE
    fit = ps.fit_slab_profile(z, second, hardwall=False)
    assert fit.success is False, fit
    assert "second condensate" in fit.reason

    scattered = rho.copy()
    scattered[[3, 9]] = 0.5
    fit = ps.fit_slab_profile(z, scattered, hardwall=False)
    assert fit.success is True, fit
    assert fit.rho_dense == pytest.approx(_RHO_DENSE, abs=0.02)


@pytest.mark.parametrize("walls", ["low", "both"])
def test_second_slab_check_on_the_hardwall_fits(walls):
    """The wetting-film and both-walls fits apply the check too: without it a
    second slab of 3 or 4 planes came back with success and rho_dilute 0.09 - 0.12."""
    z = np.arange(48, dtype=float)
    rho = np.full(48, _RHO_DILUTE_SLAB)
    if walls == "low":
        rho[:8] = _RHO_DENSE
        rho[28:31] = _RHO_DENSE
    else:
        rho[:6] = rho[-6:] = _RHO_DENSE
        rho[22:26] = _RHO_DENSE
    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success is False, fit
    assert "second condensate" in fit.reason


def test_one_dense_plane_on_each_wall_is_not_a_second_slab():
    """Review D3: between hard walls the two end planes are not neighbours. A free
    slab with one adsorbed plane on each wall (above the fitted midpoint, below
    half the peak, so it is not fitted as two wall films) used to be refused as
    "2 neighbouring planes"; the same layer on one wall passed."""
    z = np.arange(40, dtype=float)
    rho = np.full(40, _RHO_DILUTE_SLAB)
    rho[15:25] = 0.85
    rho[20] = 0.97
    rho[0] = rho[-1] = 0.47
    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success is True, fit


def test_second_slab_must_stand_clear_of_the_plateau_scatter():
    """Two neighbouring planes above the midpoint, but within three robust standard
    deviations of a plateau that scatters that much, are not a second slab."""
    model = np.full(40, 0.1)
    model[15:25] = 0.9
    noisy = model.copy()
    noisy[np.r_[0:15, 25:40]] = np.tile([0.0, 0.25], 15)     # robust sigma about 0.37
    noisy[[3, 4]] = 0.55
    ok, reason = ps._dense_planes_on_dilute_plateau(noisy, model, 0.9, 0.1)
    assert ok is True, reason
    quiet = model.copy()
    quiet[[3, 4]] = 0.55
    ok, reason = ps._dense_planes_on_dilute_plateau(quiet, model, 0.9, 0.1)
    assert ok is False
    assert "2 neighbouring planes" in reason


# ---------------------------------------------------------------------------
# F3: a periodic trajectory loaded under hard walls
# ---------------------------------------------------------------------------


def test_periodic_trajectory_loaded_under_hard_walls_warns(traj3d_files):
    """Bonds through a face of the box cannot exist under hard walls."""
    xtc, pdb, keyfile = traj3d_files
    with pytest.warns(UserWarning, match="go through a face of the box"):
        traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile, hardwall=True)
    assert traj.hardwall is True
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert not [w for w in caught if "go through a face" in str(w.message)]


def test_frames_wrapped_through_the_walls_are_reported_once(tmp_path):
    """A hardwall frame wider than the box is wrapped, which tears its bonds; the
    out-of-box warning says so, and the face check used to say it again."""
    from pimms.lemonade.tests.test_deep_audit_load import _write

    chain = np.zeros((1, 12, 3), dtype=np.int32)
    chain[0, :, 0] = np.arange(12) - 1
    xtc, pdb = _write(tmp_path, chain, ["A" * 12], (10, 6, 5))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        lemonade.load(xtc, pdb, hardwall=True)
    messages = [str(w.message) for w in caught]
    assert any("outside the hardwall box" in m for m in messages)
    assert not any("go through a face" in m for m in messages)


def test_hardwall_trajectory_does_not_trip_the_face_check(resized_equilibration_run):  # noqa: F811
    """The compact phase of a RESIZED_EQUILIBRATION run is a real hardwall run."""
    import os
    _sim, rundir, keyfile = resized_equilibration_run
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(xtc=os.path.join(rundir, "eq_traj.xtc"),
                             pdb=os.path.join(rundir, "eq_START.pdb"), keyfile=keyfile)
    assert traj.hardwall is True
    assert not [w for w in caught if "go through a face" in str(w.message)]


def test_face_check_on_hand_made_frames():
    """One chain of four beads across the x face of a 10-box, and the same chain whole."""
    torn = np.array([[[8, 5, 5], [9, 5, 5], [0, 5, 5], [1, 5, 5]]], dtype=np.int32)
    whole = torn.copy()
    whole[0, :, 0] = [3, 4, 5, 6]
    offsets = np.array([0, 4])
    with pytest.warns(UserWarning, match=r"1 of 1 frames .* \(1 in frame 0\)"):
        assert _load._check_bonds_inside_walls(torn, offsets, (10, 10, 10)) is True
    # every frame is checked, not only the first
    with pytest.warns(UserWarning, match=r"1 of 2 frames .* \(1 in frame 1\)"):
        assert _load._check_bonds_inside_walls(np.concatenate([whole, torn]), offsets,
                                               (10, 10, 10)) is True
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _load._check_bonds_inside_walls(whole, offsets, (10, 10, 10)) is False
        # a frame the out-of-box check has already reported is not reported again
        assert _load._check_bonds_inside_walls(torn, offsets, (10, 10, 10),
                                               already_reported=np.array([True])) is False
        # a bond that is not a single step even through the minimum image is the
        # other check's business, not this one's
        assert _load._check_bonds_inside_walls(
            np.array([[[1, 5, 5], [4, 5, 5]]], dtype=np.int32), np.array([0, 2]),
            (10, 10, 10)) is False


# ---------------------------------------------------------------------------
# F4: silent type coercion
# ---------------------------------------------------------------------------


def test_chain_index_must_be_an_integer(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    frame = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)[0]
    for bad in (2.9, "1", 1.0, np.float64(2.0)):
        with pytest.raises(TypeError, match="chain index must be an integer"):
            frame[bad]
        with pytest.raises(TypeError, match="chain index must be an integer"):
            frame.polymer(bad)
    assert frame[np.int64(1)].chain_index == 1
    assert frame[-1].chain_index == frame.n_chains - 1


def test_bool_temperature_is_refused(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    with pytest.raises(ValueError, match="temperature"):
        lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile, temperature=True)
    with pytest.raises(ValueError, match="temperature"):
        TrajectoryStore(np.zeros((1, 1, 3), dtype=np.int32), (8, 8, 8), 3.65, False,
                        Topology(["A"]), temperature=np.bool_(True))
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    with pytest.raises(ValueError, match="temperature"):
        st.slab_surface_tension(traj, temperature=True)
    with pytest.raises(ValueError, match="temperature"):
        st.droplet_surface_tension(traj, temperature=True)


def test_bool_slice_of_a_trajectory_is_refused(traj3d_files):
    """Review D2: ``traj[True:]`` dropped frame 0 and ``traj[:True]`` kept one frame."""
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    for key in (slice(True, None), slice(None, True), slice(None, None, True),
                slice(np.bool_(False), None)):
        with pytest.raises(TypeError, match="not a bool"):
            traj[key]
    assert traj[1:].n_frames == traj.n_frames - 1
    assert traj[np.int64(1):].n_frames == traj.n_frames - 1


@pytest.mark.parametrize("selection", [{"start": True}, {"stop": True}, {"step": True},
                                       {"start": np.bool_(False)}])
def test_bool_frame_selection_is_refused(traj3d_files, selection):
    """``start=True`` was a numpy slice from 1, which dropped frame 0 silently."""
    xtc, pdb, keyfile = traj3d_files
    with pytest.raises(TypeError, match="not a bool"):
        lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile, **selection)
    with pytest.raises(TypeError, match="not a bool"):
        lemonade.load(pdb=pdb, keyfile=keyfile, **selection)


# ---------------------------------------------------------------------------
# F6: the both-walls fit reports its starting guess when the optimiser fails
# ---------------------------------------------------------------------------


def test_both_walls_fit_reports_its_starting_guess_when_the_optimiser_fails(monkeypatch):
    """Every other fit returns its starting guess, as the BinodalFit docstring says."""
    z = np.arange(40, dtype=float)
    rho = np.full(40, _RHO_DILUTE_SLAB)
    rho[:6] = rho[-6:] = _RHO_DENSE        # a film on each wall, 28 dilute planes

    def failing_curve_fit(*args, **kwargs):
        """Stand-in for ``scipy.optimize.curve_fit`` that never converges.

        Parameters
        ----------
        args : tuple
            Ignored.
        kwargs : dict
            Ignored.

        Raises
        ------
        RuntimeError
            Always.
        """
        raise RuntimeError("no convergence")

    monkeypatch.setattr("scipy.optimize.curve_fit", failing_curve_fit)
    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.reason == "curve_fit did not converge"
    assert fit.success is False
    # the starting guess: a dilute region half as wide as the 28 planes below
    # the midpoint, so films of mean thickness (40 - 28) / 2 = 6
    assert fit.half_width == pytest.approx(0.25 * (40 - 28))


# ---------------------------------------------------------------------------
# F7: upper-casing keeps a sequence's length
# ---------------------------------------------------------------------------


def test_keyfile_sequences_are_upper_cased_as_pimms_does(tmp_path):
    """The keyfile chains must match the START.pdb the engine wrote for the same keyfile.

    ``'aßıſéε'.upper()`` is the seven beads ``'ASSISÉΕ'``, and the engine keeps
    ``ß``, ``ı`` and ``ſ`` as written: the keyfile stopped matching the PDB and
    its chain types were dropped with a warning that it belongs to another run.
    """
    from pimms.keyfile_parser import upper_case_chain

    chain = "aßıſéε"
    xtc, pdb, keyfile = _make_trajectory(
        tmp_path, [10, 10, 10], [(3, chain), (2, "εa")], seed=3,
        params=_UNICODE_PARAMS, n_steps=60, xtc_freq=20)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert not [w for w in caught if "does not match the PDB topology" in str(w.message)]
    assert traj.topology.sequences == ["".join(_UNICODE_BEADS)] * 3 + ["ΕA"] * 2
    # lemonade's copy of the rule agrees with the engine's, character by character
    sample = ("abcxyzABC019_ßıſéÿεωΐẞ"
              "ﬁŉǅⅰⓐ")
    for bead in sample:
        assert _load._upper_case_chain(bead) == upper_case_chain(bead), bead
    assert _load._upper_case_chain(chain) == "".join(_UNICODE_BEADS)
