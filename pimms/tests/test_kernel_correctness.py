"""
Kernel correctness tests - bit-exactness and energy-consistency for the
optimized Cython kernels across the full matrix of:

    dimensionality (2D/3D) x forcefield (SR / SR+LR / SR+LR+SLR) x HARDWALL (T/F)

Kernels covered: fast serial crankshaft (mega_crank / mega_crank_2D), the
parallel checkerboard kernel (mega_crank_parallel, 3D), the slither/reptation
megamove (mega_slither / mega_slither_2D) and the TSMMC moves (via a short
simulation guarded by ENERGY_CHECK).

See kernel_test_utils for the forcefield/system construction and the meaning of
the two properties (energy-consistency and detailed balance).
"""
import time

import numpy as np
import pytest

from pimms.tests import kernel_test_utils as U

DIMS = (2, 3)
FORCEFIELDS = ("SR", "LR", "SLR")
HARDWALLS = (True, False)

# every (dim, ff, hardwall) combination
ALL_CASES = [(d, ff, hw) for d in DIMS for ff in FORCEFIELDS for hw in HARDWALLS]
CASE_IDS = [f"{d}D-{ff}-{'HW' if hw else 'PBC'}" for (d, ff, hw) in ALL_CASES]


def _crank_moves():
    return {"MOVE_CRANKSHAFT": 1.0}


# ---------------------------------------------------------------------------
# The serial kernels must use the platform-independent splitmix64 generator,
# NOT the macOS libc rand() (Park-Miller MINSTD LCG: x -> 16807*x mod 2^31-1),
# which is weak (short period ~2.1e9, lattice structure) and platform-dependent.
# (mega_crank's crand_test/seed_C_rand expose the shared generator; the fast
# kernel is proven identical to it by test_fast_crank_bit_exact_vs_reference.)
# ---------------------------------------------------------------------------
def test_serial_prng_is_splitmix64_not_minstd():
    import pimms.mega_crank as mc

    # fixed, platform-independent range (libc RAND_MAX is 2^31-1 on macOS but
    # 32767 on Windows; ours is pinned)
    assert mc.RAND_MAX_test() == 2147483647

    # the MINSTD signature is that seeds 1,2,3 yield first outputs 16807,33614,
    # 50421 (= 16807 * seed). splitmix64 must not do that.
    firsts = []
    for s in (1, 2, 3):
        mc.seed_C_rand(s)
        firsts.append(mc.crand_test())
    assert firsts != [16807, 33614, 50421], "serial PRNG reverted to MINSTD!"
    assert firsts[1] != 2 * firsts[0]

    # uniformity, full range, and ~zero lag-1 autocorrelation over a sample
    mc.seed_C_rand(987654321)
    x = np.array([mc.crand_test() for _ in range(100000)], dtype=np.float64) / 2147483647.0
    assert 0.49 < x.mean() < 0.51
    assert x.min() < 0.02 and x.max() > 0.98
    assert abs(np.corrcoef(x[:-1], x[1:])[0, 1]) < 0.01

    # reproducible: same seed -> same stream
    mc.seed_C_rand(42); a = [mc.crand_test() for _ in range(5)]
    mc.seed_C_rand(42); b = [mc.crand_test() for _ in range(5)]
    assert a == b


# ---------------------------------------------------------------------------
# meta-test: the SR / LR / SLR forcefields must actually populate the
# corresponding energy terms, otherwise the range-specific tests are vacuous.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("ff", FORCEFIELDS)
def test_forcefield_exercises_expected_ranges(tmp_path, dim, ff):
    st = U.build_state(tmp_path, dim, ff, True, _crank_moves())
    _tot, sr, lr, slr, _ang = st.energy_terms
    assert sr != 0, "short-range term should always be populated"
    if ff in ("LR", "SLR"):
        assert lr != 0, "LR forcefield must populate the long-range term"
        assert st.has_LR()
    else:
        assert lr == 0 and slr == 0 and not st.has_LR()
    if ff == "SLR":
        assert slr != 0, "SLR forcefield must populate the super-long-range term"
    else:
        assert slr == 0


# ---------------------------------------------------------------------------
# bit-exactness: the fast serial crankshaft kernel must be byte-for-byte
# identical to the reference mega_crank kernel (same RNG stream).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_fast_crank_bit_exact_vs_reference(tmp_path, dim, ff, hardwall):
    st = U.build_state(tmp_path, dim, ff, hardwall, _crank_moves())
    substeps = 4000
    # identical bead selector + seed for both kernels
    bsel = U.make_bead_selector(st, substeps, seed=123)

    g1, t1, i1 = st.fresh()
    e_ref = U.crank_megastep(st, g1, t1, i1, st.energy, 4242,
                             substeps=substeps, fast=False, bsel=bsel)
    g2, t2, i2 = st.fresh()
    e_fast = U.crank_megastep(st, g2, t2, i2, st.energy, 4242,
                              substeps=substeps, fast=True, bsel=bsel)

    assert e_fast == e_ref, "incremental energy differs between fast and reference"
    assert np.array_equal(np.asarray(g1), np.asarray(g2)), "grid differs"
    assert np.array_equal(np.asarray(t1), np.asarray(t2)), "type_grid differs"
    assert np.array_equal(np.asarray(i1), np.asarray(i2)), "idx_to_bead differs"
    # and the bit-exact result is itself energy-consistent
    assert e_fast == U.recompute_energy(st, g2, t2, i2)


# ---------------------------------------------------------------------------
# energy-consistency: incremental energy == from-scratch recompute, over
# several megamoves (so state accumulates), for crankshaft and slither.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_crank_energy_consistency(tmp_path, dim, ff, hardwall):
    st = U.build_state(tmp_path, dim, ff, hardwall, _crank_moves())
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    for m in range(6):
        e = U.crank_megastep(st, g, t, i, e, 100 + m)
        assert e == U.recompute_energy(st, g, t, i), f"crank energy drift at megastep {m}"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"


@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_slither_energy_consistency(tmp_path, dim, ff, hardwall):
    st = U.build_state(tmp_path, dim, ff, hardwall, {"MOVE_SLITHER": 1.0})
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    for m in range(6):
        e = U.slither_megastep(st, g, t, i, e, 100 + m)
        assert e == U.recompute_energy(st, g, t, i), f"slither energy drift at megastep {m}"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"


@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_pull_energy_consistency(tmp_path, dim, ff, hardwall):
    # the cooperative-pull cascade + revert must keep the incrementally-tracked
    # energy exactly equal to a from-scratch recompute, and conserve beads.
    st = U.build_state(tmp_path, dim, ff, hardwall, {"MOVE_PULL": 1.0})
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    total_acc = 0
    for m in range(6):
        e, acc = U.pull_megastep(st, g, t, i, e, 100 + m)
        total_acc += acc
        assert e == U.recompute_energy(st, g, t, i), f"pull energy drift at megastep {m}"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"
    assert total_acc > 0, "pull accepted no moves (kernel may be all-reject)"


# ---------------------------------------------------------------------------
# VMMC (virtual-move Monte Carlo collective move, code 14) - a Simulation-level
# move, so energy-consistency is checked via a short ENERGY_CHECK run rather than
# a kernel driver. The move computes dE as a from-scratch total-energy recompute
# and reverts on reject, so a clean run proves the SR+LR+SLR dE stays exact
# through accepts and reverts across every dim/forcefield/boundary combination.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_vmmc_energy_consistency(tmp_path, dim, ff, hardwall):
    U.run_sim_with_energy_check(tmp_path, dim, ff, hardwall,
                                {"MOVE_CRANKSHAFT": 0.6, "MOVE_VMMC": 0.4},
                                n_steps=150, energy_check=3, seed=7,
                                extra={"VMMC_MAX_DISPLACEMENT": 2})


def test_vmmc_performs_collective_moves(tmp_path):
    # The whole point of VMMC is collective motion: in a dense, attractive system
    # it must recruit and accept multi-chain cluster translations, not just
    # degenerate single-chain moves. (Detailed balance for those collective moves
    # is checked separately in test_detailed_balance.py.)
    _, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SLR", False, {"MOVE_CRANKSHAFT": 0.5, "MOVE_VMMC": 0.5},
        n_steps=1500, energy_check=200, seed=7, temperature=45,
        box=[16, 16, 16], chains=[(8, "AABB"), (8, "AAAA"), (8, "AABBA")],
        extra={"VMMC_MAX_DISPLACEMENT": 2}, return_sim=True)
    assert sim.vmmc_accepted_multichain > 0, "VMMC accepted no multi-chain cluster moves"
    assert sim.vmmc_max_accepted_cluster >= 2


# ---------------------------------------------------------------------------
# Jump-and-relax (code 13) - a Simulation-level composite move (relax -> jump ->
# relax). Each sub-step is committed/reverted in place, so a clean ENERGY_CHECK
# run proves the tracked energy stays exact through the relaxations, the accepted
# jumps and the reverted jumps across every dim/forcefield/boundary combination.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_jump_and_relax_energy_consistency(tmp_path, dim, ff, hardwall):
    U.run_sim_with_energy_check(tmp_path, dim, ff, hardwall,
                                {"MOVE_CRANKSHAFT": 0.5, "MOVE_JUMP_AND_RELAX": 0.5},
                                n_steps=150, energy_check=3, seed=7,
                                extra={"CRANKSHAFT_SUBSTEPS": 300})


def test_jump_and_relax_relocates_chains(tmp_path):
    # The move must actually relocate chains: in a dilute-ish system its jump
    # (step 2) should be accepted at least sometimes. (Detailed balance is checked
    # in test_detailed_balance.py.)
    _, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SLR", False, {"MOVE_CRANKSHAFT": 0.5, "MOVE_JUMP_AND_RELAX": 0.5},
        n_steps=1500, energy_check=200, seed=7, temperature=55,
        box=[20, 20, 20], chains=[(10, "AABB"), (10, "AAAA")],
        extra={"CRANKSHAFT_SUBSTEPS": 300}, return_sim=True)
    # accepted_count[13] counts accepted jumps within the jump-and-relax move
    assert sim.ACC.accepted_count[13] > 0, "jump-and-relax never accepted a jump"


def test_pull_throughput(tmp_path):
    # the pull is a fast Cython megamove: a large batch must complete quickly and
    # accept a non-trivial fraction (guards against an all-reject regression, e.g.
    # an nR inversion or a broken empty/occupied predicate).
    st = U.build_state(tmp_path, 3, "SLR", False, {"MOVE_PULL": 1.0},
                       box=[20, 20, 20],
                       chains=[(18, "AABB"), (18, "AAAA"), (12, "AABBA")])
    g, t, i = st.fresh()
    e = st.energy
    accepted = 0
    pulls = 0
    npull = 30
    nchains = len(U.chain_meta(i)[0])
    t0 = time.perf_counter()
    for m in range(20):
        e, acc = U.pull_megastep(st, g, t, i, e, m, substeps=npull)
        accepted += acc
        pulls += npull * nchains
    dt = time.perf_counter() - t0
    assert accepted > 0, "no pull moves accepted"
    # very generous wall-clock budget (thousands of pulls should be well under a second)
    assert dt < 30.0, f"pull throughput too slow: {pulls} pulls in {dt:.2f}s"


# the parallel checkerboard kernel is 3D-only; check several thread counts so
# multi-block decomposition (where cross-block races would show up) is exercised.
@pytest.mark.parametrize("ff", FORCEFIELDS)
@pytest.mark.parametrize("hardwall", HARDWALLS)
@pytest.mark.parametrize("nthreads", (1, 2, 4))
def test_parallel_energy_consistency(tmp_path, ff, hardwall, nthreads):
    # a larger, dispersed box so the domain decomposition forms multiple blocks
    st = U.build_state(tmp_path, 3, ff, hardwall, _crank_moves(),
                       box=[30, 30, 30],
                       chains=[(20, "AABB"), (20, "AAAA"), (15, "A")])
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    for m in range(4):
        e = U.parallel_megastep(st, g, t, i, e, 555 + m, nthreads=nthreads)
        assert e == U.recompute_energy(st, g, t, i), \
            f"parallel energy drift at megastep {m} (nthreads={nthreads})"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"


# ---------------------------------------------------------------------------
# TSMMC energy-consistency: TSMMC moves are coordinated by the Simulation, so we
# run a short in-process simulation with ENERGY_CHECK enabled. A from-scratch
# energy recomputation that disagreed with the tracked energy raises
# SimulationEnergyException, so simply completing the run proves consistency.
# ---------------------------------------------------------------------------
TSMMC_MOVES = {
    "ctsmmc": {"MOVE_CRANKSHAFT": 0.5, "MOVE_CTSMMC": 0.5},
    "multichain": {"MOVE_CRANKSHAFT": 0.5, "MOVE_MULTICHAIN_TSMMC": 0.5},
    "system": {"MOVE_CRANKSHAFT": 0.5, "MOVE_SYSTEM_TSMMC": 0.5},
}


@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
@pytest.mark.parametrize("tsmmc", list(TSMMC_MOVES), ids=list(TSMMC_MOVES))
def test_tsmmc_energy_consistency(tmp_path, dim, ff, hardwall, tsmmc):
    # Completes iff every ENERGY_CHECK recomputation matched the tracked energy,
    # otherwise SimulationEnergyException is raised.
    U.run_sim_with_energy_check(tmp_path, dim, ff, hardwall, TSMMC_MOVES[tsmmc],
                                n_steps=60, energy_check=5)


# ---------------------------------------------------------------------------
# PARALLELIZE keyword end-to-end. The parallel checkerboard kernel is 3D-only;
# in 2D the PARALLELIZE flag must gracefully fall back to the serial 2D kernel.
# Either way a short run with ENERGY_CHECK must stay energy-consistent.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_parallelize_simulation_energy_consistency(tmp_path, dim, ff, hardwall):
    U.run_sim_with_energy_check(tmp_path, dim, ff, hardwall,
                                {"MOVE_CRANKSHAFT": 1.0}, n_steps=40, energy_check=5,
                                extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 4})
