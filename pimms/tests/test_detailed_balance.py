"""
Detailed-balance tests for the optimized moves.

A move that respects detailed balance samples the same Boltzmann equilibrium as
the trusted crankshaft reference (the fast serial crankshaft kernel is bit-exact
to the reference mega_crank kernel - see test_kernel_correctness). Each test
equilibrates a system with crankshaft, then from that SAME configuration runs
both crankshaft and the move-under-test and checks their equilibrium energies
agree. The known detailed-balance bugs (the parallel frozen-halo bug, the
endpoint-only TSMMC acceptance) drove the energy far from equilibrium and would
fail these comfortably.

Everything is seeded -> deterministic, so these do not flake.

These are heavier than the correctness tests; they are marked `slow` so they can
be deselected with `-m "not slow"` when iterating quickly.
"""
import os
import contextlib

import numpy as np
import pytest

from pimms.tests import kernel_test_utils as U
from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation

pytestmark = pytest.mark.slow

HARDWALLS = (True, False)
HW_IDS = {True: "HW", False: "PBC"}


# ---------------------------------------------------------------------------
# slither (reptation) megamove - 2D and 3D, hardwall + PBC, full SLR forcefield
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_slither_detailed_balance(tmp_path, dim, hardwall):
    st = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_SLITHER": 1.0})
    ref, test = U.db_compare(st, U.slither_megastep, equilibrate=120, sample=120)
    U.assert_same_equilibrium(ref, test, f"slither {dim}D {HW_IDS[hardwall]} SLR")


# ---------------------------------------------------------------------------
# pull (cooperative reptation) megamove - 2D and 3D, hardwall + PBC, full SLR.
#
# The pull is not ergodic on its own (it cannot translate whole chains or move
# single-bead chains), so detailed balance is tested the way the move is actually
# used: the pull as the DOMINANT move with a little crankshaft for ergodicity must
# reach the SAME Boltzmann equilibrium as crankshaft alone. A pull that violated
# detailed balance would continuously bias the chain and shift that equilibrium
# away from the crankshaft reference.
# ---------------------------------------------------------------------------
def _pull_dominant_step(state, g, t, idx, energy, seed):
    e, _ = U.pull_megastep(state, g, t, idx, energy, seed, substeps=30)
    return U.crank_megastep(state, g, t, idx, e, seed + 777, substeps=400)


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_pull_detailed_balance(tmp_path, dim, hardwall):
    st = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0})
    ref, test = U.db_compare(st, _pull_dominant_step, equilibrate=150, sample=150)
    U.assert_same_equilibrium(ref, test, f"pull {dim}D {HW_IDS[hardwall]} SLR")


# ---------------------------------------------------------------------------
# parallel checkerboard kernel - 3D, hardwall + PBC, full SLR forcefield.
# A dispersed box so the domain decomposition forms multiple blocks (the regime
# where the historical frozen-halo detailed-balance bug appeared).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_parallel_detailed_balance(tmp_path, hardwall):
    st = U.build_state(tmp_path, 3, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                       box=[30, 30, 30],
                       chains=[(22, "AABB"), (22, "AAAA"), (18, "A")])

    def parallel_step(state, g, t, i, e, seed):
        return U.parallel_megastep(state, g, t, i, e, seed, nthreads=4)

    ref, test = U.db_compare(st, parallel_step, equilibrate=110, sample=110,
                             crank_substeps=4500)
    U.assert_same_equilibrium(ref, test, f"parallel 3D {HW_IDS[hardwall]} SLR")


# ---------------------------------------------------------------------------
# parallel checkerboard kernel - 2D (mega_crank_parallel_2D). A dispersed 2D box
# that decomposes into multiple blocks, so the frozen-halo handling is exercised;
# the 2D parallel kernel must reach the same equilibrium as the serial crankshaft.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_parallel_2D_detailed_balance(tmp_path, hardwall):
    st = U.build_state(tmp_path, 2, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                       box=[40, 40],
                       chains=[(22, "AABB"), (22, "AAAA"), (18, "A")])

    def parallel_step(state, g, t, i, e, seed):
        return U.parallel_megastep_2D(state, g, t, i, e, seed, nthreads=4)

    ref, test = U.db_compare(st, parallel_step, equilibrate=110, sample=110,
                             crank_substeps=4500)
    U.assert_same_equilibrium(ref, test, f"parallel 2D {HW_IDS[hardwall]} SLR")


# ---------------------------------------------------------------------------
# TSMMC moves - coordinated by the Simulation, so compared via two full runs
# (crankshaft-only reference vs TSMMC+crankshaft) reaching equilibrium.
# ---------------------------------------------------------------------------
def _sim_equilibrium(tmp_path, sub, dim, ff, hardwall, moves, *, n_steps, seed):
    """Run a short simulation and return the 2nd-half (mean, std) of ENERGY.dat."""
    d = tmp_path / sub
    d.mkdir()
    extra = {
        "TSMMC_JUMP_TEMP": 120,
        "TSMMC_STEP_MULTIPLIER": 20,
        "TSMMC_NUMBER_OF_POINTS": 10,
        "EN_FREQ": 10,
        "PRINT_FREQ": 1000000,
        "XTC_FREQ": 1000000,
        "ANALYSIS_FREQ": 1000000,
        "RESTART_FREQ": 1000000,
    }
    U.write_param_file(os.path.join(str(d), "params.prm"), ff)
    U.write_keyfile(os.path.join(str(d), "KEYFILE.kf"), dim, hardwall, moves,
                    temperature=40, n_steps=n_steps, equilibration=n_steps // 4,
                    seed=seed, extra=extra)
    cwd = os.getcwd()
    os.chdir(str(d))
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            Simulation(keyfile.keyword_lookup).run_simulation()
        e = np.loadtxt("ENERGY.dat", delimiter="\t")[:, 1]
    finally:
        os.chdir(cwd)
    e = e[len(e) // 2:]
    return e.mean(), e.std()


@pytest.mark.parametrize("tsmmc_move", ["MOVE_CTSMMC", "MOVE_SYSTEM_TSMMC", "MOVE_MULTICHAIN_TSMMC"])
def test_tsmmc_detailed_balance(tmp_path, tsmmc_move):
    # 3D, short-range, PBC. Reference is crankshaft only; the test run mixes the
    # TSMMC move with crankshaft. TSMMC is a temperature-excursion move, so a
    # broken acceptance biases the sampled energy well away from the reference.
    n_steps = 4000
    ref_mean, ref_std = _sim_equilibrium(
        tmp_path, "ref", 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
        n_steps=n_steps, seed=21)
    test_mean, test_std = _sim_equilibrium(
        tmp_path, "test", 3, "SR", False,
        {"MOVE_CRANKSHAFT": 0.7, tsmmc_move: 0.3}, n_steps=n_steps, seed=21)

    tol = 2.5 * max(ref_std, test_std) + 0.05 * abs(ref_mean)
    assert abs(test_mean - ref_mean) <= tol, (
        f"{tsmmc_move}: detailed balance violated - crank E={ref_mean:.1f}+/-{ref_std:.1f}, "
        f"tsmmc E={test_mean:.1f}+/-{test_std:.1f}, |diff|={abs(test_mean - ref_mean):.1f} > tol={tol:.1f}")


# ---------------------------------------------------------------------------
# VMMC (virtual-move Monte Carlo collective move, code 14). VMMC recruits and
# rigidly translates a cluster of chains, so it only does real work in a dense,
# strongly-attractive system - exactly the regime where it matters and where a
# wrong forward/reverse link ratio or cutoff would bias the sampled energy. A
# dense 3D SLR box (small VMMC displacement so collective moves clear hard-core
# clashes) lets clusters of 2-4 chains move; crank+VMMC must reach the SAME
# Boltzmann equilibrium as crankshaft alone.
# ---------------------------------------------------------------------------
def _vmmc_equilibrium(tmp_path, sub, moves, *, seed, n_steps):
    d = tmp_path / sub
    d.mkdir()
    extra = {
        "EN_FREQ": 10,
        "PRINT_FREQ": 1000000,
        "XTC_FREQ": 1000000,
        "ANALYSIS_FREQ": 1000000,
        "RESTART_FREQ": 1000000,
        "VMMC_MAX_DISPLACEMENT": 2,
    }
    U.write_param_file(os.path.join(str(d), "params.prm"), "SLR")
    U.write_keyfile(os.path.join(str(d), "KEYFILE.kf"), 3, False, moves,
                    box=[16, 16, 16], chains=[(8, "AABB"), (8, "AAAA"), (8, "AABBA")],
                    temperature=45, n_steps=n_steps, equilibration=n_steps // 4,
                    seed=seed, extra=extra)
    cwd = os.getcwd()
    os.chdir(str(d))
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(keyfile.keyword_lookup)
            sim.run_simulation()
        e = np.loadtxt("ENERGY.dat", delimiter="\t")[:, 1]
    finally:
        os.chdir(cwd)
    e = e[len(e) // 2:]
    return e.mean(), e.std(), sim


def test_vmmc_detailed_balance(tmp_path):
    n_steps = 4000
    ref_mean, ref_std, _ = _vmmc_equilibrium(
        tmp_path, "ref", {"MOVE_CRANKSHAFT": 1.0}, seed=21, n_steps=n_steps)
    test_mean, test_std, sim = _vmmc_equilibrium(
        tmp_path, "test", {"MOVE_CRANKSHAFT": 0.5, "MOVE_VMMC": 0.5}, seed=21, n_steps=n_steps)

    # the test is only meaningful if VMMC actually performed collective moves
    assert sim.vmmc_accepted_multichain > 5, (
        f"VMMC accepted too few multi-chain moves ({sim.vmmc_accepted_multichain}) - "
        f"the system is not exercising the collective move, so this is not a real test")
    assert sim.vmmc_max_accepted_cluster >= 2

    tol = 2.5 * max(ref_std, test_std) + 0.05 * abs(ref_mean)
    assert abs(test_mean - ref_mean) <= tol, (
        f"VMMC detailed balance violated - crank E={ref_mean:.1f}+/-{ref_std:.1f}, "
        f"vmmc E={test_mean:.1f}+/-{test_std:.1f}, |diff|={abs(test_mean - ref_mean):.1f} > tol={tol:.1f} "
        f"(accepted {sim.vmmc_accepted_multichain} multi-chain moves, max cluster {sim.vmmc_max_accepted_cluster})")


# ---------------------------------------------------------------------------
# Jump-and-relax (code 13). The move is composed of three sub-steps that each
# individually preserve the Boltzmann distribution (relax -> Metropolis-accepted
# jump -> relax), so crank+jump-and-relax must reach the SAME equilibrium as
# crankshaft alone. The earlier deferred-acceptance formulation (one accept/reject
# on the post-relaxation energy) broke detailed balance and would bias this energy.
# Run in a dilute-ish box so the jump (step 2) actually lands sometimes.
# ---------------------------------------------------------------------------
def _jr_equilibrium(tmp_path, sub, moves, *, seed, n_steps):
    d = tmp_path / sub
    d.mkdir()
    extra = {
        "EN_FREQ": 10,
        "PRINT_FREQ": 1000000,
        "XTC_FREQ": 1000000,
        "ANALYSIS_FREQ": 1000000,
        "RESTART_FREQ": 1000000,
        "CRANKSHAFT_SUBSTEPS": 400,
    }
    U.write_param_file(os.path.join(str(d), "params.prm"), "SLR")
    U.write_keyfile(os.path.join(str(d), "KEYFILE.kf"), 3, False, moves,
                    box=[20, 20, 20], chains=[(12, "AABB"), (12, "AAAA")],
                    temperature=50, n_steps=n_steps, equilibration=n_steps // 4,
                    seed=seed, extra=extra)
    cwd = os.getcwd()
    os.chdir(str(d))
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(keyfile.keyword_lookup)
            sim.run_simulation()
        e = np.loadtxt("ENERGY.dat", delimiter="\t")[:, 1]
    finally:
        os.chdir(cwd)
    e = e[len(e) // 2:]
    return e.mean(), e.std(), sim


def test_jump_and_relax_detailed_balance(tmp_path):
    n_steps = 6000
    ref_mean, ref_std, _ = _jr_equilibrium(
        tmp_path, "ref", {"MOVE_CRANKSHAFT": 1.0}, seed=21, n_steps=n_steps)
    test_mean, test_std, sim = _jr_equilibrium(
        tmp_path, "test", {"MOVE_CRANKSHAFT": 0.5, "MOVE_JUMP_AND_RELAX": 0.5}, seed=21, n_steps=n_steps)

    # only meaningful if the jump actually fires
    assert sim.ACC.accepted_count[13] > 0, "jump-and-relax accepted no jumps - not a real test"

    tol = 2.5 * max(ref_std, test_std) + 0.05 * abs(ref_mean)
    assert abs(test_mean - ref_mean) <= tol, (
        f"jump-and-relax detailed balance violated - crank E={ref_mean:.1f}+/-{ref_std:.1f}, "
        f"J&R E={test_mean:.1f}+/-{test_std:.1f}, |diff|={abs(test_mean - ref_mean):.1f} > tol={tol:.1f}")
