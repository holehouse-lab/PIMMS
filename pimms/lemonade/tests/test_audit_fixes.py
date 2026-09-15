"""
Regression tests for the phase-separation / surface-tension defects found by the
third 1.0.8 audit. Every trajectory here is synthetic with a hand-computable
answer; nothing is compared against another lemonade or PIMMS routine.
"""

import warnings

import numpy as np
import pytest

from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology
from pimms.lemonade.trajectory import LatticeTrajectory
from pimms.lemonade import phase_separation as ps
from pimms.lemonade import surface_tension as st


def _make(frames, seqs, dims, hardwall=False, temperature=1.0):
    store = TrajectoryStore(np.asarray(frames, dtype=np.int32), tuple(dims), 3.65, hardwall,
                            Topology(seqs), temperature=temperature)
    return LatticeTrajectory(store)


def _random_monomers(rng, dims, density, exclude=None):
    sites = np.indices(dims).reshape(len(dims), -1).T
    if exclude is not None:
        keep = ~np.array([tuple(s) in exclude for s in sites])
        sites = sites[keep]
    n = int(round(density * len(sites)))
    return sites[rng.choice(len(sites), n, replace=False)]


# ---------------------------------------------------------------------------
# HIGH: a homogeneous solution must not be reported as phase separated
# ---------------------------------------------------------------------------

def test_homogeneous_random_occupancy_is_not_phase_separated():
    rng = np.random.default_rng(0)
    dims = (12, 12, 12)
    frames = [_random_monomers(rng, dims, 0.3) for _ in range(3)]
    seqs = ["A"] * len(frames[0])
    traj = _make(frames, seqs, dims)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ps.analyze(traj, geometry="sphere")
    assert res.percolation_fraction == 1.0
    assert res.spanning_fraction == 1.0
    assert res.binodal.success is False
    assert "spans the box" in res.binodal.reason
    assert res.is_phase_separated is False
    # the fallback densities are percentiles of the observed (~0.3) profile, not 1.0
    assert res.rho_dense < 0.8


def test_innermost_single_site_shell_cannot_manufacture_a_droplet():
    """Synthetic profile: a flat 0.3 background with the one-site r<1 shell reading
    0.8 - what the largest cluster's own centre looks like in a homogeneous box."""
    radii = np.arange(0.5, 6.0, 1.0)
    density = np.full(len(radii), 0.3)
    density[0] = 0.8
    sites = np.array([1, 26, 98, 218, 386, 602], dtype=float)[:len(radii)]
    fit = ps.fit_radial_profile(radii, density, site_counts=sites)
    assert fit.success is False
    # and even with equal weights the fitted radius (< 1) is rejected
    fit_unweighted = ps.fit_radial_profile(radii, density)
    assert fit_unweighted.success is False


def test_spanning_fraction_distinguishes_slab_from_network():
    dims = (8, 8, 24)
    slab = [[x, y, z] for x in range(8) for y in range(8) for z in range(8, 16)]
    traj = _make([slab], ["A" * 8] * 64, dims)
    assert ps.spanning_fraction(traj, all_axes=False) == 1.0    # spans x and y
    assert ps.spanning_fraction(traj, all_axes=True) == 0.0     # not z


def test_spanning_needs_a_pair_that_touches_through_the_face():
    """A contact staircase from x = 0 to x = 5 in a 6-box reaches the box length
    on x but its two ends differ by 3 in y, so no pair meets through the face: the
    cluster does not wind and has an unambiguous single image. PIMMS's gather
    already knew that; lemonade's extent-only detector called it spanning. Under a
    hardwall the same staircase touches both x walls, which is what spanning means
    there."""
    stairs = [[0, 0, 3], [1, 0, 3], [1, 1, 3], [2, 1, 3], [2, 2, 3],
              [3, 2, 3], [3, 3, 3], [4, 3, 3], [5, 3, 3]]
    seqs = ["A"] * len(stairs)
    assert ps.spanning_fraction(_make([stairs], seqs, (6, 6, 6))) == 0.0
    assert ps.spanning_fraction(_make([stairs], seqs, (6, 6, 6), hardwall=True)) == 1.0
    rod = [[x, 2, 3] for x in range(6)]
    assert ps.spanning_fraction(_make([rod], ["A"] * 6, (6, 6, 6))) == 1.0


def _droplet_frames_plus_one_spanning_frame(rng):
    """Nine frames of a 5x5x5 cube droplet in a 16-box and one frame whose largest
    cluster is a ring around the x axis (connected to its own image) with an arm,
    every other bead an isolated monomer. Returns (frames, seqs, droplet_frames)."""
    L, n_beads = 16, 125
    droplet = [[x, y, z] for x in range(5, 10) for y in range(5, 10) for z in range(5, 10)]
    ring = [[x, 8, 8] for x in range(L)]
    arm = [[3, 8 + k, 8] for k in range(1, 6)]
    occ = set(map(tuple, ring + arm))
    rest = []
    while len(ring) + len(arm) + len(rest) < n_beads:
        s = tuple(int(v) for v in rng.integers(0, L, size=3))
        if any(max(abs(s[i] - o[i]) for i in range(3)) <= 1 for o in occ):
            continue
        occ.add(s)
        rest.append(list(s))
    return [droplet] * 9 + [ring + arm + rest], ["A"] * n_beads, [droplet] * 9


def test_droplet_shape_leaves_out_frames_whose_largest_cluster_spans_the_box():
    """One spanning frame in ten used to move the mean asphericity of a perfect
    cube from 0.0 to about 1.5; it is now left out, and a warning says so."""
    frames, seqs, droplet_only = _droplet_frames_plus_one_spanning_frame(np.random.default_rng(3))
    traj = _make(frames, seqs, (16, 16, 16))
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        shape = ps.droplet_shape(traj)
    assert any("droplet_shape: 1 of 10 frames" in str(w.message) for w in rec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reference = ps.droplet_shape(_make(droplet_only, seqs, (16, 16, 16)))
    for key in ("radius_of_gyration", "asphericity", "sphericity", "volume", "density"):
        assert shape[key] == pytest.approx(reference[key])
    assert shape["asphericity"] == pytest.approx(0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ps.analyze(traj, geometry="sphere")
    assert res.shape["asphericity"] == pytest.approx(0.0)
    assert res.spanning_fraction == pytest.approx(0.1)


def test_radial_profile_leaves_out_spanning_frames_unless_every_frame_spans():
    """The same spanning frame is left out of the radial profile (the profile of
    the nine droplet frames is reproduced exactly). When every frame spans, the
    profile about the arbitrary centre is still returned, so analyze() keeps its
    percentile fallback for the two densities on a network."""
    frames, seqs, droplet_only = _droplet_frames_plus_one_spanning_frame(np.random.default_rng(3))
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        r, rho = ps.radial_density_profile(_make(frames, seqs, (16, 16, 16)))
    assert any("radial_density_profile: 1 of 10 frames" in str(w.message) for w in rec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r_ref, rho_ref = ps.radial_density_profile(_make(droplet_only, seqs, (16, 16, 16)))
    assert np.array_equal(r, r_ref)
    np.testing.assert_allclose(rho, rho_ref, equal_nan=True)

    spanning_only = [frames[-1]] * 3
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        r_s, rho_s = ps.radial_density_profile(_make(spanning_only, seqs, (16, 16, 16)))
    assert any("every one of the 3 frames" in str(w.message) for w in rec)
    assert np.isfinite(rho_s).any()


def test_hardwall_spanning_frames_are_left_out_as_wall_bounded_condensates():
    """Under a hardwall the ring frame's largest cluster touches both x walls. Its
    centre of mass is exact, so the old 'no single image' rationale does not apply;
    it is left out of the radial profile and the shape averages because it is a
    wall-bounded film or network rather than a droplet, and the warnings say so."""
    frames, seqs, droplet_only = _droplet_frames_plus_one_spanning_frame(np.random.default_rng(3))
    traj = _make(frames, seqs, (16, 16, 16), hardwall=True)
    assert traj[9].clusters[0].spanning_axes() == [0]
    assert traj[0].clusters[0].spanning_axes() == []
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        r, rho = ps.radial_density_profile(traj)
        shape = ps.droplet_shape(traj)
    msgs = [str(w.message) for w in rec]
    assert any("radial_density_profile: 1 of 10 frames" in m and "touches both walls" in m
               for m in msgs)
    assert any("droplet_shape: 1 of 10 frames" in m and "touches both walls" in m for m in msgs)
    assert not any("single-image gather" in m for m in msgs)     # nothing is gathered
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = _make(droplet_only, seqs, (16, 16, 16), hardwall=True)
        r_ref, rho_ref = ps.radial_density_profile(ref)
        shape_ref = ps.droplet_shape(ref)
    assert np.array_equal(r, r_ref)
    np.testing.assert_allclose(rho, rho_ref, equal_nan=True)
    for key in ("radius_of_gyration", "asphericity", "sphericity", "volume", "density"):
        assert shape[key] == pytest.approx(shape_ref[key])
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        r_all, rho_all = ps.radial_density_profile(_make([frames[-1]] * 2, seqs, (16, 16, 16),
                                                         hardwall=True))
    assert any("every one of the 2 frames" in str(w.message)
               and "wall-bounded condensate" in str(w.message) for w in rec)
    assert np.isfinite(rho_all).any()


def test_cluster_spanning_axes_runs_the_percolation_test_once_per_gather(monkeypatch):
    """The gather's own percolation test is switched off for lemonade clusters and
    run once on the cached image: one visit per axis per cluster, the documented
    gather warning still raised on the first gather, and the answer cached so the
    spanning detector and the shape calls never re-test."""
    from pimms import cluster_utils
    frames, seqs, _ = _droplet_frames_plus_one_spanning_frame(np.random.default_rng(3))
    calls = []
    real = cluster_utils._axis_percolates

    def counting(arr, dimensions, d, *args, **kwargs):
        calls.append(d)
        return real(arr, dimensions, d, *args, **kwargs)

    monkeypatch.setattr(cluster_utils, "_axis_percolates", counting)
    ring = _make([frames[-1]], seqs, (16, 16, 16))[0].clusters[0]
    with pytest.warns(UserWarning, match="single-image gather: cluster percolates.*axis 0"):
        ring.single_image_positions()
    assert calls == [0, 1, 2]
    assert ring.spanning_axes() == [0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ring.single_image_positions()
        ring.spanning_axes()
        assert ps._cluster_spans_box(ring, False) is True
        assert ps._cluster_spans_box(ring, True) is False
        ring.radius_of_gyration
    assert calls == [0, 1, 2]
    # the detector on a fresh object gathers once with the warning muted
    calls.clear()
    fresh = _make([frames[-1]], seqs, (16, 16, 16))[0].clusters[0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert ps._cluster_spans_box(fresh, False) is True
    assert calls == [0, 1, 2]
    # a compact droplet: no warning, no spanning axis
    droplet = _make([frames[0]], seqs, (16, 16, 16))[0].clusters[0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert droplet.spanning_axes() == []


# ---------------------------------------------------------------------------
# MEDIUM: hardwall radial profiles
# ---------------------------------------------------------------------------

def _droplet_plus_vapour(rng, dims, n_frames):
    block = {(8 + x, 8 + y, 8 + z) for x in range(5) for y in range(5) for z in range(5)}
    frames = []
    for _ in range(n_frames):
        vapour = _random_monomers(rng, dims, 0.03, exclude=block)
        frames.append(np.vstack([np.asarray(sorted(block)), vapour]))
    # keep the vapour count identical across frames for a clean topology
    n_vap = min(len(f) for f in frames) - 125
    frames = [np.vstack([f[:125], f[125:125 + n_vap]]) for f in frames]
    seqs = ["A"] * (125 + n_vap)
    rho_v = n_vap / (np.prod(dims) - 125)
    return frames, seqs, rho_v


def test_hardwall_radial_profile_marks_empty_shells_as_nan_not_zero():
    rng = np.random.default_rng(1)
    dims = (20, 20, 20)
    frames, seqs, _ = _droplet_plus_vapour(rng, dims, 2)
    traj = _make(frames, seqs, dims, hardwall=True)
    r, rho = ps.radial_density_profile(traj)
    # beyond the farthest in-box site from the centre there is nothing to occupy
    assert np.isnan(rho[-1])
    assert not np.any(rho[np.isfinite(rho)] < 0)
    r2, rho2, sites = ps.radial_density_profile_with_site_counts(traj)
    assert np.array_equal(r, r2)
    assert sites[0] == 1.0                    # the r<1 shell is one site


@pytest.mark.parametrize("hardwall", [False, True])
def test_dilute_density_recovered_with_and_without_hardwall(hardwall):
    rng = np.random.default_rng(2)
    dims = (20, 20, 20)
    frames, seqs, rho_v_hand = _droplet_plus_vapour(rng, dims, 12)
    traj = _make(frames, seqs, dims, hardwall=hardwall)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ps.analyze(traj, geometry="sphere")
    assert res.binodal.success is True
    assert res.rho_dense > 0.9
    assert res.rho_dilute == pytest.approx(rho_v_hand, rel=0.2)
    assert res.is_phase_separated is True


# ---------------------------------------------------------------------------
# MEDIUM: slab tools under hardwall
# ---------------------------------------------------------------------------

def _rough_face(rng, Lx, Ly):
    return rng.integers(0, 3, size=(Lx, Ly))


def test_wetting_slab_surface_tension_counts_only_the_free_face():
    """A slab wetting z=0 under hardwall has ONE interface. Its gamma must equal
    that of a periodic slab whose two faces carry the same height field (mirrored),
    since both faces then have identical capillary power. The wall face used to be
    averaged in with zero power, doubling gamma."""
    rng = np.random.default_rng(3)
    Lx, Ly, Lz = 8, 8, 40
    dims = (Lx, Ly, Lz)
    top = _rough_face(rng, Lx, Ly)

    wet, wet_seqs = [], []
    for x in range(Lx):
        for y in range(Ly):
            col = [[x, y, z] for z in range(0, 10 + top[x, y])]
            wet.extend(col)
            wet_seqs.append("A" * len(col))
    wet_traj = _make([wet], wet_seqs, dims, hardwall=True)

    mirrored, mir_seqs = [], []
    for x in range(Lx):
        for y in range(Ly):
            col = [[x, y, z] for z in range(15 - top[x, y], 25 + top[x, y])]
            mirrored.extend(col)
            mir_seqs.append("A" * len(col))
    ref_traj = _make([mirrored], mir_seqs, dims, hardwall=False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wet_res = st.slab_surface_tension(wet_traj, temperature=1.0)
        ref_res = st.slab_surface_tension(ref_traj, temperature=1.0)
    assert np.isfinite(wet_res.gamma) and wet_res.n_modes == 8
    assert wet_res.gamma == pytest.approx(ref_res.gamma, rel=1e-6)


def test_wetting_slab_profile_is_not_recentred_and_gets_a_single_interface_fit():
    rng = np.random.default_rng(4)
    Lx, Ly, Lz = 8, 8, 40
    dims = (Lx, Ly, Lz)
    top = _rough_face(rng, Lx, Ly)
    wet, seqs = [], []
    for x in range(Lx):
        for y in range(Ly):
            col = [[x, y, z] for z in range(0, 10 + top[x, y])]
            wet.extend(col)
            seqs.append("A" * len(col))
    traj = _make([wet], seqs, dims, hardwall=True)
    z, rho = ps.slab_density_profile(traj, axis=2)
    assert rho[0] == 1.0 and rho[-1] == 0.0          # dense phase left at the wall
    fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert fit.success is True
    assert fit.rho_dense == pytest.approx(1.0, abs=0.02)
    assert fit.rho_dilute == pytest.approx(0.0, abs=0.02)
    assert 2.0 * fit.half_width == pytest.approx(10.0 + top.mean(), abs=1.0)
    # the shape-based default makes the same decision
    assert ps.fit_slab_profile(z, rho).success is True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ps.analyze(traj, geometry="slab")
    assert res.binodal.success is True and res.is_phase_separated is True


# ---------------------------------------------------------------------------
# MEDIUM: droplet surface tension grid
# ---------------------------------------------------------------------------

def test_droplet_grid_defaults_scale_with_droplet_size():
    dims = (40, 40, 40)
    g = np.arange(40)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    inside = (X - 20) ** 2 + (Y - 20) ** 2 + (Z - 20) ** 2 <= 8.0 ** 2
    sphere = np.stack([X[inside], Y[inside], Z[inside]], axis=1)
    traj = _make([sphere], ["A"] * len(sphere), dims)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = st.droplet_surface_tension(traj, temperature=1.0)
    assert res.n_polar == 16 and res.n_azim == 32          # 2 * R0 with R0 ~ 8
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        forced = st.droplet_surface_tension(traj, temperature=1.0, n_polar=8, n_azim=16)
    assert forced.n_polar == 8 and forced.n_azim == 16


# ---------------------------------------------------------------------------
# final-review fixes
# ---------------------------------------------------------------------------

def _slab_frames(dims, centres, half=5):
    """Perfect slabs (density 1) of thickness 2*half along z, one per centre."""
    nx, ny, nz = dims
    frames = []
    for c in centres:
        pts = [(x, y, z) for x in range(nx) for y in range(ny) for z in range(c - half, c + half)]
        frames.append(pts)
    return np.asarray(frames, dtype=np.int32)


def test_hardwall_free_slab_is_aligned_before_averaging():
    """A perfect slab diffusing between the walls used to be averaged without
    alignment: the profile was a broad hump (rho_dense 0.29, thickness 34 for a
    slab of thickness 10) reported with success=True."""
    dims = (8, 8, 48)
    frames = _slab_frames(dims, [10, 15, 22, 28, 33, 38])
    seqs = ["A"] * frames.shape[1]
    traj = _make(frames, seqs, dims, hardwall=True)
    res = ps.analyze(traj, geometry="slab")
    assert res.binodal.success, res.binodal.reason
    assert res.binodal.rho_dense == pytest.approx(1.0, abs=0.05)
    assert res.binodal.rho_dilute == pytest.approx(0.0, abs=0.05)
    assert 2 * res.binodal.half_width == pytest.approx(10.0, abs=0.7)


def test_hardwall_off_centre_slab_fits_with_a_free_centre():
    z = np.arange(40, dtype=float)
    dens = np.zeros(40)
    dens[5:15] = 1.0
    fit = ps.fit_slab_profile(z, dens, hardwall=True)
    assert fit.success, fit.reason
    assert fit.rho_dense == pytest.approx(1.0, abs=0.05)
    assert fit.rho_dilute == pytest.approx(0.0, abs=0.05)
    assert fit.half_width == pytest.approx(5.0, abs=0.6)


def test_inverted_two_interface_fit_is_rejected():
    """A condensate wetting both walls looks, to the two-interface model, like a
    dilute slab in a dense background; the fit used to converge with the dense
    and dilute densities swapped and success=True."""
    z = np.arange(40, dtype=float)
    dens = np.ones(40)
    dens[12:28] = 0.05
    fit = ps.fit_slab_profile(z, dens, hardwall=False)
    assert not fit.success
    assert "inverted" in fit.reason


def test_spanning_fallback_reads_the_dense_density_off_usable_shells():
    """A uniform 30 % occupancy in a 12-box: the fallback percentiles used to
    include the one-site innermost shell and reported rho_dense = 0.6."""
    rng = np.random.default_rng(3)
    dims = (12, 12, 12)
    occ = rng.random(dims) < 0.3
    pts = np.argwhere(occ).astype(np.int32)
    traj = _make(pts[None, :, :], ["A"] * len(pts), dims)
    res = ps.analyze(traj, geometry="droplet")
    assert not res.is_phase_separated
    assert res.binodal.rho_dense < 0.45, res.binodal.rho_dense
    assert res.binodal.rho_dilute == pytest.approx(0.3, abs=0.06)


def test_2d_radial_fit_uses_a_2d_shell_floor():
    """A perfect disc of radius ~3.5 in a 12x12 box: with the 3D floor of 20 sites
    only three shells survived and the fit was refused as 'fewer than four usable
    radial shells'."""
    dims = (12, 12)
    pts = [(x, y, 0) for x in range(12) for y in range(12) if (x - 6) ** 2 + (y - 6) ** 2 <= 3.5 ** 2]
    traj = _make(np.asarray([pts], dtype=np.int32), ["A"] * len(pts), dims)
    res = ps.analyze(traj, geometry="droplet")
    assert res.binodal.success, res.binodal.reason
    assert res.binodal.radius == pytest.approx(3.5, abs=0.6)


def test_droplet_grid_is_sized_from_the_typical_droplet_not_the_first_frame():
    """A leading not-yet-condensed frame (its 'largest cluster' is most of the
    box) used to size the angular grid far too fine, after which every real
    droplet frame failed the coverage test and gamma was nan without warning."""
    import warnings
    rng = np.random.default_rng(5)
    dims = (30, 30, 30)
    soup = np.argwhere(rng.random(dims) < 0.5).astype(np.int32)
    centre = np.array([15, 15, 15])
    ball = np.array([p for p in np.ndindex(dims) if ((np.asarray(p) - centre) ** 2).sum() <= 5.0 ** 2], dtype=np.int32)
    n = min(len(soup), len(ball))
    frames = np.stack([soup[:n], ball[:n], ball[:n], ball[:n], ball[:n]])
    traj = _make(frames, ["A"] * n, dims, temperature=1.0)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = st.droplet_surface_tension(traj)
    assert res.n_polar == int(np.clip(round(2.0 * (3.0 * n / (4.0 * np.pi)) ** (1.0 / 3.0)), 8, 64))
    assert np.isfinite(res.gamma)
    assert any("skipped" in str(w.message) for w in rec)         # the soup frame


def test_blank_pdb_chain_column_is_not_a_type_label():
    from types import SimpleNamespace
    from pimms.lemonade._topology import Topology

    def chain(seq, label):
        atoms = [SimpleNamespace(residue=SimpleNamespace(name=r)) for r in seq]
        return SimpleNamespace(atoms=atoms, chain_id=label)

    top = SimpleNamespace(chains=[chain("AA", " "), chain("GG", " ")])
    types = [int(t) for t in Topology.from_mdtraj(top).chain_types]
    assert types[0] != types[1]


def test_trajectory_store_does_not_freeze_the_callers_arrays():
    arr = np.zeros((2, 3, 3), dtype=np.int32)
    arr[1, :, 0] = 1
    times = np.arange(2, dtype=np.float64)
    store = TrajectoryStore(arr, (5, 5, 5), 3.65, False, Topology(["A", "A", "A"]), times=times)
    assert arr.flags.writeable and times.flags.writeable
    assert not store.positions.flags.writeable


def test_frame_polymer_accepts_negative_indices_and_range_checks():
    # frame.polymer(-1) used to hand the raw negative index to the chain
    # table, giving a Polymer with an empty bead range whose len() raised.
    arr = np.zeros((1, 5, 3), dtype=np.int32)
    arr[0, :, 0] = np.arange(5)
    store = TrajectoryStore(arr, (9, 9, 9), 3.65, False,
                            Topology(["AA", "BBB"]),
                            times=np.zeros(1, dtype=np.float64))
    frame = LatticeTrajectory(store)[0]

    assert frame.polymer(-1).chain_index == frame[-1].chain_index == 1
    assert len(frame.polymer(-1)) == 3
    assert frame.polymer(-2).chain_index == 0
    with pytest.raises(IndexError):
        frame.polymer(2)
    with pytest.raises(IndexError):
        frame.polymer(-3)


def test_load_rejects_a_frame_selection_that_keeps_nothing(traj3d_files):
    # a window past the end of a short run used to surface as an IndexError from
    # the box lookup rather than a message about the selection
    from pimms import lemonade
    xtc, pdb, keyfile = traj3d_files
    with pytest.raises(ValueError, match="keeps no frames"):
        lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile, start=10_000, stop=20_000)
    with pytest.raises(ValueError, match="keeps no frames"):
        lemonade.load(xtc=xtc, pdb=pdb, start=10_000)
