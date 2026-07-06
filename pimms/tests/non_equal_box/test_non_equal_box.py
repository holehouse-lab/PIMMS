"""
Non-cubic (unequal-dimension) box tests.

PIMMS supports boxes where x != y != z (3D) or x != y (2D), but such boxes are
rarely used in practice, so a routine that mishandled a single axis would be a
latent bug. These tests establish that:

1. The Cython kernels keep the tracked energy consistent with a from-scratch
   recompute on non-cubic boxes whose SHORT axis forces cross-boundary
   interactions and chain wrapping (``ENERGY_CHECK`` would raise on a per-axis
   PBC bug).                                              -- fast
2. Permuting the axes of a non-cubic box is an exact physical invariance, so the
   equilibrium ensembles must match; each permutation lands the short axis on a
   different physical axis, stressing every axis's PBC path.   -- slow
3. A cubic box and a volume-matched non-cubic box give the same bulk ensemble in
   the dilute limit (chains small vs every axis).              -- slow
4. ``compute_cluster_radial_density_profile`` sizes its shells from the SMALLEST
   axis (fixed bug: it used the largest), matching an independent reference.  -- fast

The equal-vs-unequal comparisons are statistical (not bit-exact); tolerances follow
the detailed-balance suite's ``assert_same_equilibrium`` convention.
"""

import os
import contextlib

import numpy as np
import pytest

from pimms.tests import kernel_test_utils as U
from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation
from pimms import lattice_utils
from pimms import lattice_analysis_utils as lau


# A 16-bead chain; with a short axis of 8 it cannot avoid wrapping/interacting
# across the periodic boundary once it diffuses.
_WRAPPING_CHAINS = [(6, "AABB" * 4)]
_MOVES = {"MOVE_CRANKSHAFT": 0.7, "MOVE_SLITHER": 0.3}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _run_sim(dirpath, dim, ff, hardwall, moves, *, box, chains, n_steps, seed,
             temperature=45, with_rg=False):
    """Run a short in-process simulation and return (energy_trace, mean_rg, sim).

    ``energy_trace`` is the full ENERGY.dat energy column; ``mean_rg`` is the
    mean radius of gyration over the 2nd half of RG.dat (or None if with_rg is
    False). ``sim`` is the finished Simulation (its LATTICE holds the final
    configuration, e.g. for a boundary-straddle check).
    """
    os.makedirs(dirpath, exist_ok=True)
    extra = {
        "EN_FREQ": 10,
        "PRINT_FREQ": 1000000,
        "XTC_FREQ": 1000000,
        "RESTART_FREQ": 1000000,
        # cluster analysis is expensive and not needed here
        "ANA_CLUSTER": 1000000,
    }
    if with_rg:
        extra["ANALYSIS_FREQ"] = 1000000
        extra["ANA_POL"] = 20
    else:
        extra["ANALYSIS_FREQ"] = 1000000

    U.write_param_file(os.path.join(dirpath, "params.prm"), ff)
    U.write_keyfile(os.path.join(dirpath, "KEYFILE.kf"), dim, hardwall, moves,
                    box=box, chains=chains, temperature=temperature,
                    n_steps=n_steps, equilibration=n_steps // 4, seed=seed, extra=extra)

    cwd = os.getcwd()
    os.chdir(dirpath)
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(keyfile.keyword_lookup)
            sim.run_simulation()
        energy = np.loadtxt("ENERGY.dat", delimiter="\t")[:, 1]
        mean_rg = None
        if with_rg and os.path.isfile("RG.dat"):
            rows = []
            with open("RG.dat") as fh:
                for line in fh:
                    parts = line.split()          # split() ignores the trailing tab
                    if len(parts) >= 2:
                        rows.append([float(x) for x in parts[1:]])  # drop the step column
            if rows:
                arr = np.array(rows[len(rows) // 2:])   # 2nd half = production
                mean_rg = float(arr.mean())
    finally:
        os.chdir(cwd)
    return energy, mean_rg, sim


def _chebyshev_reference_profile(cluster, dims):
    """Independent reference: density at Chebyshev shell k = (# beads at shell k) /
    (# lattice sites in shell k), with the scan depth bounded by the SMALLEST axis
    (matching the fixed compute_cluster_radial_density_profile)."""
    com = np.array(lattice_utils.center_of_mass_from_positions(cluster.tolist(), dims))
    d = len(dims)
    offset_max = int(min(dims) / 2) - 1
    occ = set(map(tuple, cluster.tolist()))
    n_beads = len(cluster)
    max_beads = n_beads - (1 if tuple(com.tolist()) in occ else 0)
    prof, found = [], 0
    for k in range(1, offset_max + 1):
        total = (2 * k + 1) ** d - (2 * k - 1) ** d
        occ_k = sum(1 for p in cluster.tolist() if max(abs(np.array(p) - com)) == k)
        prof.append(occ_k / total)
        found += occ_k
        if found == max_beads:
            break
    if len(prof) < offset_max:
        prof.extend([0] * (offset_max - len(prof)))
    return prof


def _compact_cluster(center, radius):
    """A filled Chebyshev ball of the given radius about an integer center."""
    c = np.array(center)
    d = len(center)
    pts = []
    rng = range(-radius, radius + 1)
    if d == 2:
        for dx in rng:
            for dy in rng:
                pts.append((c + np.array([dx, dy])).tolist())
    else:
        for dx in rng:
            for dy in rng:
                for dz in rng:
                    pts.append((c + np.array([dx, dy, dz])).tolist())
    return np.array(pts)


# ---------------------------------------------------------------------------
# 1. energy consistency on non-cubic PBC boxes (fast) - front-line kernel check
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("ff", ["SR", "LR", "SLR"])
def test_energy_consistency_noncubic_pbc(tmp_path, ff):
    # short axis = 8 forces continual cross-boundary interactions AND chain
    # wrapping; ENERGY_CHECK recomputes the energy from scratch (per-axis-correct
    # Python Hamiltonian) and raises if the Cython-kernel-tracked energy drifts.
    trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, ff, hardwall=False, moves=_MOVES,
        box=[8, 16, 32], chains=[(8, "AABB" * 4)],
        n_steps=600, energy_check=25, seed=17, temperature=45, return_sim=True)
    assert trace.shape[0] > 0
    assert sim.LATTICE.any_chains_straddle_boundary(), (
        "expected at least one chain to cross the short-axis PBC boundary "
        "(test would be vacuous otherwise)")


@pytest.mark.parametrize("ff", ["SR", "SLR"])
def test_energy_consistency_noncubic_2d_pbc(tmp_path, ff):
    trace, sim = U.run_sim_with_energy_check(
        tmp_path, 2, ff, hardwall=False, moves=_MOVES,
        box=[8, 48], chains=[(8, "AABB" * 4)],
        n_steps=600, energy_check=25, seed=19, temperature=45, return_sim=True)
    assert trace.shape[0] > 0
    assert sim.LATTICE.any_chains_straddle_boundary()


def test_energy_consistency_noncubic_hardwall(tmp_path):
    # hardwall confines rather than wraps: energy must still stay consistent, and
    # no chain may straddle a boundary (there is no periodic image to cross into).
    trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SLR", hardwall=True, moves=_MOVES,
        box=[8, 16, 32], chains=[(8, "AABB" * 4)],
        n_steps=600, energy_check=25, seed=17, temperature=45, return_sim=True)
    assert trace.shape[0] > 0
    assert not sim.LATTICE.any_chains_straddle_boundary()


# ---------------------------------------------------------------------------
# 2. axis-permutation equilibrium equivalence (slow) - the rigorous test
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_axis_permutation_equilibrium_pbc(tmp_path):
    # [8,16,32] and its cyclic permutations are the SAME physical system with axes
    # relabelled, so their equilibrium ensembles must be identical. The short axis
    # (8) sits on x, then z, then y - exercising each axis's PBC wrap in turn.
    n_steps = 4000
    boxes = {"xyz": [8, 16, 32], "zxy": [32, 8, 16], "yzx": [16, 32, 8]}
    results = {}
    for name, box in boxes.items():
        e, rg, sim = _run_sim(str(tmp_path / name), 3, "SLR", False, _MOVES,
                              box=box, chains=_WRAPPING_CHAINS, n_steps=n_steps,
                              seed=23, with_rg=True)
        results[name] = (e[len(e) // 2:], rg, sim)

    ref_e, ref_rg, _ = results["xyz"]
    for name in ("zxy", "yzx"):
        test_e, test_rg, sim = results[name]
        U.assert_same_equilibrium(ref_e, test_e, f"axis-perm {name} vs xyz",
                                  rel_floor=0.05, k_sigma=3.0)
        # mean Rg is a scalar invariant too; a broken per-axis unwrap would inflate
        # the Rg of a boundary-crossing chain on the affected axis.
        assert abs(test_rg - ref_rg) <= 0.15 * abs(ref_rg) + 0.5, (
            f"axis-perm {name}: mean Rg {test_rg:.3f} vs xyz {ref_rg:.3f}")

    # non-vacuous: the short axis forces boundary crossings in every permutation
    for name in boxes:
        assert results[name][2].LATTICE.any_chains_straddle_boundary(), (
            f"{name}: expected a chain to cross the short-axis PBC boundary")


# ---------------------------------------------------------------------------
# 3. cubic vs volume-matched non-cubic, dilute limit (slow)
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_cubic_vs_noncubic_dilute_pbc(tmp_path):
    # Volume-matched (12^3 == 8*12*18 == 1728) and dilute: 4-bead chains have an
    # extent well below the short axis (8), so there is no self-image interaction
    # and the bulk ensemble is box-shape-independent. Kept dispersed with a warm
    # temperature so chains do not aggregate into a box-spanning droplet.
    n_steps = 4000
    chains = [(10, "AABB")]
    e_cubic, _, _ = _run_sim(str(tmp_path / "cubic"), 3, "SLR", False, _MOVES,
                             box=[12, 12, 12], chains=chains, n_steps=n_steps,
                             seed=31, temperature=70)
    e_noncubic, _, _ = _run_sim(str(tmp_path / "noncubic"), 3, "SLR", False, _MOVES,
                                box=[8, 12, 18], chains=chains, n_steps=n_steps,
                                seed=31, temperature=70)
    U.assert_same_equilibrium(e_cubic[len(e_cubic) // 2:], e_noncubic[len(e_noncubic) // 2:],
                              "cubic vs non-cubic (dilute)", rel_floor=0.06, k_sigma=3.0)


# ---------------------------------------------------------------------------
# 4. radial density profile: shells sized from the SMALLEST axis (the fix)
# ---------------------------------------------------------------------------
def test_radial_density_profile_noncubic_uses_min_axis():
    # box [10,10,40]: min axis 10 -> offset_max = int(10/2)-1 = 4. Before the fix
    # this used max axis 40 -> offset_max 19, tying the profile to the wrong axis.
    dims = [10, 10, 40]
    cluster = _compact_cluster([5, 5, 20], radius=2)
    profile = lau.compute_cluster_radial_density_profile([cluster], dims)[0]
    assert len(profile) == int(min(dims) / 2) - 1 == 4
    ref = _chebyshev_reference_profile(cluster, dims)
    assert np.allclose(profile, ref), f"{profile} != {ref}"


def test_radial_density_profile_axis_swap_invariant():
    # the same cluster, with box + coordinates permuted, must give the same profile
    base_dims = [10, 10, 40]
    cluster = _compact_cluster([5, 5, 20], radius=2)
    profiles = []
    for perm in [(0, 1, 2), (2, 0, 1), (1, 2, 0)]:
        dims = [base_dims[p] for p in perm]
        permuted = cluster[:, list(perm)]
        profiles.append(lau.compute_cluster_radial_density_profile([permuted], dims)[0])
    for p in profiles[1:]:
        assert np.allclose(profiles[0], p)


def test_radial_density_profile_cubic_unchanged():
    # cubic box: min == max, so behaviour is exactly as before the fix; still must
    # match the independent Chebyshev reference.
    dims = [20, 20, 20]
    cluster = _compact_cluster([10, 10, 10], radius=3)
    profile = lau.compute_cluster_radial_density_profile([cluster], dims)[0]
    assert len(profile) == int(min(dims) / 2) - 1 == 9
    ref = _chebyshev_reference_profile(cluster, dims)
    assert np.allclose(profile, ref)


def test_radial_density_profile_2d_noncubic():
    dims = [8, 40]
    cluster = _compact_cluster([4, 20], radius=1)
    profile = lau.compute_cluster_radial_density_profile([cluster], dims)[0]
    assert len(profile) == int(min(dims) / 2) - 1 == 3
    ref = _chebyshev_reference_profile(cluster, dims)
    assert np.allclose(profile, ref)
