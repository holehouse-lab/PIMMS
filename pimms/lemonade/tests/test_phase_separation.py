"""
Tests for the phase-separation / droplet-physics analysis.

Machinery (shapes, bounds, fits) is checked on the general 3D fixture; the physics
(a real dense/dilute split, elevated condensed fraction) is checked on a strongly
self-attracting system that phase separates.
"""

import warnings

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import phase_separation as ps

# What radial_density_profile says when it leaves out a frame whose largest
# cluster spans the box (frame 0 of the condensed fixture does: see
# test_radial_profile_is_occupied_fraction).
SPANNING_FRAME_LEFT_OUT = "left out because the largest cluster spans the box"


# ---------------------------------------------------------------------------
# order parameters
# ---------------------------------------------------------------------------

def test_order_parameter_shapes_and_bounds(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)

    cf = ps.condensed_fraction(traj)
    assert cf.shape == (traj.n_frames,)
    assert np.all((cf >= 0) & (cf <= 1))

    nc = ps.number_of_clusters(traj, min_beads=1)
    assert nc.shape == (traj.n_frames,)
    assert np.all(nc >= 1)                      # at least one cluster per frame

    lc = ps.largest_cluster_size(traj, by="beads")
    assert lc.shape == (traj.n_frames,)
    assert np.all(lc <= traj.n_beads)

    sizes = ps.cluster_size_distribution(traj)
    assert sizes.ndim == 1 and sizes.sum() > 0


# ---------------------------------------------------------------------------
# density profiles are occupied fractions in [0, 1]
# ---------------------------------------------------------------------------

def test_radial_profile_is_occupied_fraction(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    # frame 0 of this fixture is the pre-equilibration start, a random placement
    # that percolates as a contact network. A cluster connected to its own
    # periodic image has no centre to take a profile about, so the standalone
    # profile leaves that frame out and says so - which is asserted here
    with pytest.warns(UserWarning, match=SPANNING_FRAME_LEFT_OUT):
        r, rho = ps.radial_density_profile(traj)
    assert r.shape == rho.shape
    assert np.all(rho >= -1e-9) and np.all(rho <= 1.0 + 1e-9)
    # a condensate: dense near the COM, dilute far away
    assert rho[0] > rho[-1]


def test_slab_profile_bounds_and_length(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    z, rho = ps.slab_density_profile(traj, axis=2)
    assert len(z) == traj.dimensions[2]
    assert np.all(rho >= -1e-9) and np.all(rho <= 1.0 + 1e-9)


# ---------------------------------------------------------------------------
# tanh fits
# ---------------------------------------------------------------------------

def test_fits_return_ordered_binodal(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)

    # see test_radial_profile_is_occupied_fraction for why this warns
    with pytest.warns(UserWarning, match=SPANNING_FRAME_LEFT_OUT):
        r, rho = ps.radial_density_profile(traj)
    fit = ps.fit_radial_profile(r, rho)
    assert isinstance(fit, ps.BinodalFit)
    assert 0 <= fit.rho_dilute <= fit.rho_dense <= 1
    z, rhoz = ps.slab_density_profile(traj, axis=2)
    sfit = ps.fit_slab_profile(z, rhoz)
    assert 0 <= sfit.rho_dilute <= sfit.rho_dense <= 1


# ---------------------------------------------------------------------------
# the fits must REFUSE a profile that is not two-phase
#
# These are regression tests for a real bug. ``curve_fit`` does not raise on a flat
# profile - it converges to a very wide tanh, which over a finite window is nearly a
# straight line, and then places rho_dense/rho_dilute wherever it likes (usually at
# the bounds). The fit "succeeded" and reported rho_dense = 1.0 with an interface
# 17 sites wide for a system that was provably homogeneous. Anything above the
# critical temperature hit this.
# ---------------------------------------------------------------------------

def _slab_profile(length=60, rho_d=0.9, rho_v=0.02, half_width=10.0, width=1.5):
    z = np.arange(length, dtype=float)
    centre = length / 2.0
    rho = ps._tanh_slab(z, rho_d, rho_v, half_width, width, centre)
    return z, rho


def test_slab_fit_accepts_a_real_slab():
    z, rho = _slab_profile()
    fit = ps.fit_slab_profile(z, rho)

    assert fit.success and fit.reason == ""
    assert fit.rho_dense == pytest.approx(0.9, abs=0.02)
    assert fit.rho_dilute == pytest.approx(0.02, abs=0.02)
    assert fit.interface_width == pytest.approx(1.5, rel=0.2)


def test_slab_fit_rejects_a_flat_profile():
    """A homogeneous (supercritical) system must not yield a coexistence gap."""
    z = np.arange(60, dtype=float)
    rho = np.full_like(z, 0.125)

    fit = ps.fit_slab_profile(z, rho)

    assert not fit.success
    assert fit.reason                                   # says why
    # and the fallback values are the observed density, not an extrapolation
    assert fit.rho_dense == pytest.approx(0.125, abs=1e-6)
    assert fit.rho_dilute == pytest.approx(0.125, abs=1e-6)
    assert fit.rho_dense - fit.rho_dilute < 1e-3        # no invented gap


def test_slab_fit_rejects_a_noisy_flat_profile():
    """The realistic version: flat to within noise, as a real supercritical run is."""
    rng = np.random.default_rng(0)
    z = np.arange(60, dtype=float)
    rho = 0.125 + rng.normal(0, 0.004, size=z.size)

    fit = ps.fit_slab_profile(z, rho)

    assert not fit.success
    # the pre-fix bug reported rho_dense ~ 1.0 here
    assert fit.rho_dense < 0.2
    assert fit.rho_dense - fit.rho_dilute < 0.05


def test_slab_fit_rejects_a_slab_that_fills_the_box():
    """With no dilute region in view, the dilute density is unconstrained."""
    z, rho = _slab_profile(length=40, half_width=25.0)   # 2*hw > length
    fit = ps.fit_slab_profile(z, rho)
    assert not fit.success


def test_radial_fit_rejects_a_flat_profile():
    r = np.arange(1, 25, dtype=float)
    rho = np.full_like(r, 0.1)

    fit = ps.fit_radial_profile(r, rho)

    assert not fit.success
    assert fit.rho_dense - fit.rho_dilute < 1e-3


def test_radial_fit_accepts_a_real_droplet():
    r = np.arange(1, 25, dtype=float)
    rho = ps._tanh_droplet(r, 0.85, 0.01, 10.0, 1.2)

    fit = ps.fit_radial_profile(r, rho)

    assert fit.success and fit.reason == ""
    assert fit.rho_dense == pytest.approx(0.85, abs=0.02)
    assert fit.radius == pytest.approx(10.0, rel=0.1)


def test_is_phase_separated_is_false_when_the_fit_is_degenerate():
    """The guard has to hold at the top level too, not just in the fit."""
    z = np.arange(60, dtype=float)
    rho = np.full_like(z, 0.125)
    fit = ps.fit_slab_profile(z, rho)

    result = ps.PhaseSeparationResult(
        geometry="slab",
        condensed_fraction=0.97,        # percolating network: high, but NOT a condensate
        condensed_fraction_series=np.full(10, 0.97),
        n_clusters=1.0,
        largest_cluster_beads=3000.0,
        binodal=fit,
        shape={},
        profile=(z, rho),
    )
    assert not result.is_phase_separated


# ---------------------------------------------------------------------------
# full analysis on a phase-separated system
# ---------------------------------------------------------------------------

def test_analyze_detects_phase_separation(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    result = ps.analyze(traj)

    assert result.geometry == "sphere"                      # cubic box
    assert 0 <= result.rho_dilute <= result.rho_dense <= 1
    assert result.rho_dense > 2 * max(result.rho_dilute, 1e-6)  # a real density gap
    assert result.condensed_fraction > 0.3                  # most material condensed
    assert result.largest_cluster_beads > 0
    # droplet shape is computed
    assert np.isfinite(result.shape["radius_of_gyration"])
    assert result.is_phase_separated


def test_analyze_geometry_override(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    slab_result = ps.analyze(traj, geometry="slab")
    assert slab_result.geometry == "slab"
    assert 0 <= slab_result.rho_dilute <= slab_result.rho_dense <= 1


# ---------------------------------------------------------------------------
# droplet shape
# ---------------------------------------------------------------------------

def test_droplet_and_sphericity(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    droplet = traj[-1].droplet
    assert droplet is not None
    assert droplet.n_beads > 0
    sph = droplet.sphericity
    # sphericity is either a valid (0, 1.3] value or nan for a degenerate hull
    assert np.isnan(sph) or (0 < sph < 1.3)


# ---------------------------------------------------------------------------
# "largest cluster" must mean MOST BEADS
#
# Regression test for a real bug. PIMMS's get_cluster_distribution orders clusters
# by the number of CHAINS they contain. lemonade then took clusters[0] as the
# condensate everywhere - condensed_fraction, largest_cluster_size, the density
# profiles, droplet_shape, frame.droplet and the surface-tension estimators. Those
# two orderings agree only when every chain is the same length; in a
# multi-component system with unequal chain lengths lemonade silently measured the
# wrong cluster. The trajectory below is built so the two orderings disagree
# maximally: one long chain (10 beads, 1 chain) against three touching monomers
# (3 beads, 3 chains).
# ---------------------------------------------------------------------------

def _two_cluster_traj():
    """A one-frame 3D trajectory with a bead-rich cluster and a chain-rich one."""
    from pimms.lemonade._topology import Topology
    from pimms.lemonade._store import TrajectoryStore
    from pimms.lemonade.trajectory import LatticeTrajectory

    positions = [[x, 0, 0] for x in range(10)]          # chain 0: 10 beads, 1 chain
    positions += [[0, 5, 0], [1, 5, 0], [2, 5, 0]]      # chains 1-3: 3 beads, 3 chains

    topology = Topology(["A" * 10, "A", "A", "A"])
    store = TrajectoryStore(np.array([positions], dtype=np.int32),
                            (12, 12, 12), 3.65, False, topology)
    return LatticeTrajectory(store)


def test_clusters_are_ordered_by_bead_count_not_chain_count():
    traj = _two_cluster_traj()
    clusters = traj[0].clusters

    assert len(clusters) == 2
    # the chain-rich cluster has MORE chains but FEWER beads - beads must win
    assert [c.n_beads for c in clusters] == [10, 3]
    assert [c.n_chains for c in clusters] == [1, 3]
    assert traj[0].droplet.n_beads == 10


def test_order_parameters_follow_the_bead_rich_cluster():
    traj = _two_cluster_traj()

    assert ps.largest_cluster_size(traj, by="beads")[0] == 10
    assert ps.largest_cluster_size(traj, by="chains")[0] == 1
    # 10 of the 13 beads sit in the condensate
    assert ps.condensed_fraction(traj)[0] == pytest.approx(10 / 13)


# ---------------------------------------------------------------------------
# slab fitting on a coordinate axis that does not start at zero
# ---------------------------------------------------------------------------

def test_slab_fit_handles_an_offset_coordinate_axis():
    """The tanh centre is anchored on coord[0], not assumed to be at the origin."""
    z, rho = _slab_profile(length=60, rho_d=0.9, rho_v=0.02, half_width=10.0, width=1.5)

    shifted = ps.fit_slab_profile(z + 100.0, rho)

    assert shifted.success and shifted.reason == ""
    assert shifted.rho_dense == pytest.approx(0.9, abs=0.02)
    assert shifted.rho_dilute == pytest.approx(0.02, abs=0.02)
    assert shifted.half_width == pytest.approx(10.0, rel=0.1)


def test_slab_density_profile_takes_no_min_beads():
    """It bins every bead in the box, so the old (never-read) argument is gone."""
    import inspect

    assert "min_beads" not in inspect.signature(ps.slab_density_profile).parameters


def test_analyze_clusters_each_frame_exactly_once(monkeypatch):
    """The connected-component search must not be repeated per analysis pass.

    ``traj[f]`` mints a fresh Frame every time it is indexed, so caching clusters on the
    Frame meant each of analyze()'s passes (condensed fraction, cluster count, largest
    cluster, density profile, droplet shape) re-ran the whole decomposition - five times
    per frame. Membership is memoised on the store instead.
    """
    from pimms import lattice_analysis_utils as lau

    traj = _two_cluster_traj()

    calls = []
    real = lau.get_cluster_distribution
    monkeypatch.setattr(lau, "get_cluster_distribution",
                        lambda *a, **k: (calls.append(1), real(*a, **k))[1])

    ps.analyze(traj)

    assert len(calls) == traj.n_frames


def test_cluster_membership_cache_is_shared_across_frame_objects():
    traj = _two_cluster_traj()

    # two independently created Frame views of the same frame
    first = traj[0].clusters
    second = traj[0].clusters

    assert [c.n_beads for c in first] == [c.n_beads for c in second]
    # the memoised membership is the same object for both
    assert traj.store.cluster_membership(0) is traj.store.cluster_membership(0)


def test_shell_site_counts_match_the_exhaustive_site_enumeration():
    """The slab-wise accumulation must reproduce the all-sites-at-once histogram.

    Listing every lattice site explicitly as float64 made this scale with box VOLUME
    (~124 MB for a 120^3 box) regardless of how many beads were being analysed.
    Histogram counts are additive, so accumulating a slab at a time is exact.
    """
    def reference(dimensions, edges):
        dims = np.asarray(dimensions, dtype=np.float64)
        grids = np.indices(dimensions).reshape(len(dimensions), -1).T.astype(np.float64)
        r = np.sqrt((ps._min_image(grids, dims) ** 2).sum(axis=1))
        return np.histogram(r, bins=edges)[0].astype(np.float64)

    for dims in [(8, 8), (10, 14), (6, 6, 6), (12, 12, 12), (9, 13, 17)]:
        edges = np.arange(0.0, min(dims) / 2 + 1.0, 1.0)
        counts = ps._shell_site_counts(dims, edges)

        assert np.array_equal(counts, reference(dims, edges))
        assert counts.sum() <= np.prod(dims)          # never more sites than the box holds


# ---------------------------------------------------------------------------
# hardwall slab profile: wall-touching and free frames, and vacated planes
# ---------------------------------------------------------------------------

def _hardwall_slab_traj(starts: list, dims: tuple, thickness: int, n_vapour: int,
                        seed: int, hardwall: bool = True) -> "lemonade.LatticeTrajectory":
    """Hand-build a trajectory of one full slab per frame plus a known vapour.

    Frame ``k`` holds a slab filling every site of planes ``starts[k] ..
    starts[k] + thickness - 1`` along the last axis, and exactly ``n_vapour``
    beads scattered over the sites outside it, so the true dense density is 1
    and the true dilute density is ``n_vapour`` over the number of sites outside
    the slab, in every frame.

    Parameters
    ----------
    starts : list of int
        First slab plane of each frame, along the last (longest) axis.
    dims : tuple of int
        Box extent, 2 or 3 entries; the slab normal is the last axis.
    thickness : int
        Number of planes the slab fills.
    n_vapour : int
        Number of vapour beads placed at random outside the slab in each frame.
    seed : int
        Seed for the vapour placement.
    hardwall : bool, optional
        Boundary condition of the returned trajectory (default ``True``).

    Returns
    -------
    pimms.lemonade.LatticeTrajectory
        One single-bead chain per bead.
    """
    from pimms.lemonade._store import TrajectoryStore
    from pimms.lemonade._topology import Topology
    from pimms.lemonade.trajectory import LatticeTrajectory

    rng = np.random.default_rng(seed)
    sites = np.array(list(np.ndindex(*dims)), dtype=np.int32)
    frames = []
    for z0 in starts:
        in_slab = (sites[:, -1] >= z0) & (sites[:, -1] < z0 + thickness)
        outside = sites[~in_slab]
        vapour = outside[rng.choice(len(outside), n_vapour, replace=False)]
        pts = np.concatenate([sites[in_slab], vapour])
        if len(dims) == 2:
            pts = np.concatenate([pts, np.zeros((len(pts), 1), dtype=np.int32)], axis=1)
        frames.append(pts)
    store = TrajectoryStore(np.stack(frames).astype(np.int32), tuple(dims), 3.65, hardwall,
                            Topology(["A"] * len(frames[0])))
    return LatticeTrajectory(store)


def test_hardwall_slab_profile_does_not_mix_wall_and_free_frames():
    """A slab that touches a wall in one frame of six must not fake a wall film.

    Each frame used to be handled on its own: the free slabs were moved to the
    window centre and the wall-touching one was left at the wall, so the average
    was a film of density 1/6 at the wall plus a slab of density 5/6, which the fit
    reported as a successful rho_dense 0.84 / rho_dilute 0.06 against a truth of
    1.0 / 0.016. Only the majority kind (here the five free frames) may be averaged,
    and the frame left out must be reported.
    """
    dims, thickness, n_vapour = (8, 8, 48), 10, 40
    starts = [0, 14, 20, 26, 30, 12]           # frame 0 touches the low wall
    traj = _hardwall_slab_traj(starts, dims, thickness, n_vapour, seed=0)
    true_dilute = n_vapour / (8 * 8 * (48 - thickness))

    with pytest.warns(UserWarning, match="the 1 wall-touching frames were left out"):
        z, rho = ps.slab_density_profile(traj)
    # every kept frame is a full slab aligned onto the same planes, so the plateau
    # is exactly 1 and nothing of the wall-touching frame survives at the wall
    assert rho.max() == pytest.approx(1.0)
    assert np.sum(np.isclose(rho, 1.0)) == thickness
    assert np.nanmax(rho[:thickness]) < 0.1

    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success
    assert fit.rho_dense == pytest.approx(1.0, abs=0.01)
    assert fit.rho_dilute == pytest.approx(true_dilute, rel=0.15)
    # a perfectly sharp step pins each interface only to within one plane
    assert 2.0 * fit.half_width == pytest.approx(thickness, abs=1.0)


def test_hardwall_slab_profile_keeps_the_wall_film_when_it_is_the_majority():
    """The converse split: five wall-wetting frames, one free - the film is kept."""
    dims, thickness = (8, 8, 48), 10
    traj = _hardwall_slab_traj([0, 0, 0, 0, 0, 20], dims, thickness, 40, seed=1)
    with pytest.warns(UserWarning, match="the 1 free-slab frames were left out"):
        z, rho = ps.slab_density_profile(traj)
    assert np.allclose(rho[:thickness], 1.0)
    assert np.nanmax(rho[thickness:]) < 0.1


def test_hardwall_homogeneous_profile_is_not_split():
    """A one-phase hardwall system must not trigger the wall/free split.

    Its dense planes are noise, so which end of the box they happen to reach is
    random; roughly a third of the frames of a 5 % random occupancy used to count
    as free slabs. Splitting those off would warn on every supercritical run.
    """
    from pimms.lemonade._store import TrajectoryStore
    from pimms.lemonade._topology import Topology
    from pimms.lemonade.trajectory import LatticeTrajectory

    rng = np.random.default_rng(3)
    dims = (12, 12, 36)
    n_sites = int(np.prod(dims))
    n = int(0.05 * n_sites)
    sites = np.array(list(np.ndindex(*dims)), dtype=np.int32)
    frames = np.stack([sites[rng.choice(n_sites, n, replace=False)] for _ in range(60)])
    traj = LatticeTrajectory(TrajectoryStore(frames, dims, 3.65, True, Topology(["A"] * n)))

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _z, rho = ps.slab_density_profile(traj)
    assert np.isfinite(rho).all()
    assert np.mean(rho) == pytest.approx(n / n_sites, rel=0.02)


@pytest.mark.parametrize("dims,thickness,n_vapour", [((16, 16, 48), 10, 19),
                                                    ((40, 60), 12, 25)])
def test_hardwall_free_slab_dilute_density_is_not_biased_low(dims, thickness, n_vapour):
    """Planes vacated by the alignment must not be padded with the median.

    At about half a bead per plane the median dilute count is 0, so padding the
    vacated planes with it pulled the hardwall dilute density down to 0.76 of the
    truth in 3D (0.79 for a 2D stripe) while the periodic profile of the very same
    frames read 1.00. The truth is exact here: every frame holds the same number
    of vapour beads outside a slab of known thickness.
    """
    rng = np.random.default_rng(5)
    length = dims[-1]
    starts = [int(rng.integers(3, length - thickness - 3)) for _ in range(40)]
    true_dilute = n_vapour / (int(np.prod(dims[:-1])) * (length - thickness))

    fits = {}
    for hardwall in (True, False):
        traj = _hardwall_slab_traj(starts, dims, thickness, n_vapour, seed=7,
                                   hardwall=hardwall)
        fits[hardwall] = ps.analyze(traj, geometry="slab").binodal
    for fit in fits.values():
        assert fit.success
        assert fit.rho_dense == pytest.approx(1.0, abs=0.01)
        assert fit.rho_dilute == pytest.approx(true_dilute, rel=0.05)


def test_hardwall_profile_marks_uncovered_planes_nan_and_the_fit_skips_them():
    """A plane no frame covers after alignment is nan, and fitting ignores it.

    Here the free slab always sits near the low wall, so every frame is moved up by
    the same amount and the lowest planes are never covered: they have no data,
    which the old padding hid by inventing some.
    """
    dims, thickness, n_vapour = (8, 8, 48), 10, 30
    traj = _hardwall_slab_traj([2] * 8, dims, thickness, n_vapour, seed=2)
    z, rho = ps.slab_density_profile(traj)
    covered = np.isfinite(rho)
    assert not covered.all()
    assert not covered[0]                                  # vacated by every frame
    assert np.sum(np.isclose(rho[covered], 1.0)) == thickness

    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success
    assert fit.rho_dense == pytest.approx(1.0, abs=0.01)
    assert fit.rho_dilute == pytest.approx(n_vapour / (8 * 8 * (48 - thickness)), rel=0.15)
