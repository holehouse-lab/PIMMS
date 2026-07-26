"""
Tests for the phase-separation / droplet-physics analysis.

Machinery (shapes, bounds, fits) is checked on the general 3D fixture; the physics
(a real dense/dilute split, elevated condensed fraction) is checked on a strongly
self-attracting system that phase separates.
"""

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import phase_separation as ps


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
    assert np.all(lc <= traj.n_atoms)

    sizes = ps.cluster_size_distribution(traj)
    assert sizes.ndim == 1 and sizes.sum() > 0


# ---------------------------------------------------------------------------
# density profiles are occupied fractions in [0, 1]
# ---------------------------------------------------------------------------

def test_radial_profile_is_occupied_fraction(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
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
