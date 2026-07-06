"""
Kernel correctness tests - bit-exactness and energy-consistency for the
optimized Cython kernels across the full matrix of:

    dimensionality (2D/3D) x forcefield (SR / SR+LR / SR+LR+SLR) x HARDWALL (T/F)

Kernels covered: fast serial crankshaft (mega_crank / mega_crank_2D), the
parallel checkerboard kernel (mega_crank_parallel / mega_crank_parallel_2D), the slither/reptation
megamove (mega_slither / mega_slither_2D) and the TSMMC moves (via a short
simulation guarded by ENERGY_CHECK).

See kernel_test_utils for the forcefield/system construction and the meaning of
the two properties (energy-consistency and detailed balance).
"""
import time
import os
import contextlib

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


# 3D parallel checkerboard kernel; check several thread counts so the multi-block
# decomposition (where cross-block races would show up) is exercised.
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


# 2D parallel checkerboard kernel (mega_crank_parallel_2D). Same multi-block /
# multi-thread coverage as the 3D test, over a 2D box.
@pytest.mark.parametrize("ff", FORCEFIELDS)
@pytest.mark.parametrize("hardwall", HARDWALLS)
@pytest.mark.parametrize("nthreads", (1, 2, 4))
def test_parallel_2D_energy_consistency(tmp_path, ff, hardwall, nthreads):
    # a larger, dispersed 2D box so the domain decomposition forms multiple blocks
    st = U.build_state(tmp_path, 2, ff, hardwall, _crank_moves(),
                       box=[40, 40],
                       chains=[(20, "AABB"), (20, "AAAA"), (15, "A")])
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    for m in range(4):
        e = U.parallel_megastep_2D(st, g, t, i, e, 555 + m, nthreads=nthreads)
        assert e == U.recompute_energy(st, g, t, i), \
            f"2D parallel energy drift at megastep {m} (nthreads={nthreads})"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"


def test_parallel_2D_thread_count_independent(tmp_path):
    # The block decomposition is independent of the thread count, and each block
    # is a disjoint, deterministically-seeded sub-problem, so the SAME state + seed
    # must produce the IDENTICAL result with 1 vs 4 threads (the per-block work is
    # deterministic and the integer delta sum is order-free).
    st = U.build_state(tmp_path, 2, "SLR", False, _crank_moves(),
                       box=[40, 40], chains=[(20, "AABB"), (20, "AAAA"), (15, "A")])
    g1, t1, i1 = st.fresh()
    g4, t4, i4 = st.fresh()
    e1 = U.parallel_megastep_2D(st, g1, t1, i1, st.energy, 99, nthreads=1)
    e4 = U.parallel_megastep_2D(st, g4, t4, i4, st.energy, 99, nthreads=4)
    assert e1 == e4, f"energy differs across thread counts: {e1} vs {e4}"
    assert np.array_equal(np.asarray(g1), np.asarray(g4)), "grid differs across thread counts"
    assert np.array_equal(np.asarray(i1), np.asarray(i4)), "idx_to_bead differs across thread counts"


# parallel SLITHER (chain-level frozen-halo: mega_slither_parallel / _2D). The box
# is large enough to decompose into multiple blocks AND for the (compact) chains to
# fit inside a block interior, so chains are genuinely distributed across blocks.
@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("ff", FORCEFIELDS)
@pytest.mark.parametrize("hardwall", HARDWALLS)
@pytest.mark.parametrize("nthreads", (1, 2, 4))
def test_parallel_slither_energy_consistency(tmp_path, dim, ff, hardwall, nthreads):
    box = [40, 40, 40] if dim == 3 else [40, 40]
    st = U.build_state(tmp_path, dim, ff, hardwall, {"MOVE_SLITHER": 1.0},
                       box=box, chains=[(20, "AABB"), (20, "AAAA"), (15, "A")])
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    for m in range(4):
        e = U.slither_parallel_megastep(st, g, t, i, e, 321 + m, nthreads=nthreads)
        assert e == U.recompute_energy(st, g, t, i), \
            f"parallel slither energy drift at megastep {m} (dim={dim}, nthreads={nthreads})"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"


@pytest.mark.parametrize("dim", DIMS)
def test_parallel_slither_thread_count_independent(tmp_path, dim):
    # As for the crankshaft kernel, the chain-level decomposition is independent of
    # the thread count, so 1 vs 4 threads must give a bit-identical result.
    box = [40, 40, 40] if dim == 3 else [40, 40]
    st = U.build_state(tmp_path, dim, "SLR", False, {"MOVE_SLITHER": 1.0},
                       box=box, chains=[(20, "AABB"), (20, "AAAA"), (15, "A")])
    g1, t1, i1 = st.fresh()
    g4, t4, i4 = st.fresh()
    e1 = U.slither_parallel_megastep(st, g1, t1, i1, st.energy, 88, nthreads=1)
    e4 = U.slither_parallel_megastep(st, g4, t4, i4, st.energy, 88, nthreads=4)
    assert e1 == e4, f"energy differs across thread counts: {e1} vs {e4}"
    assert np.array_equal(np.asarray(g1), np.asarray(g4)), "grid differs across thread counts"
    assert np.array_equal(np.asarray(i1), np.asarray(i4)), "idx_to_bead differs across thread counts"


# parallel PULL (chain-level frozen-halo: mega_pull_parallel / _2D). Chains are
# length >= 3 (pull-eligible) and compact, in a multi-block box.
@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("ff", FORCEFIELDS)
@pytest.mark.parametrize("hardwall", HARDWALLS)
@pytest.mark.parametrize("nthreads", (1, 2, 4))
def test_parallel_pull_energy_consistency(tmp_path, dim, ff, hardwall, nthreads):
    box = [40, 40, 40] if dim == 3 else [40, 40]
    st = U.build_state(tmp_path, dim, ff, hardwall, {"MOVE_PULL": 1.0},
                       box=box, chains=[(20, "AABBAB"), (20, "AAAAAA"), (15, "AB")])
    g, t, i = st.fresh()
    e = st.energy
    beads0 = int(np.count_nonzero(np.asarray(g)))
    for m in range(4):
        e = U.pull_parallel_megastep(st, g, t, i, e, 654 + m, nthreads=nthreads)
        assert e == U.recompute_energy(st, g, t, i), \
            f"parallel pull energy drift at megastep {m} (dim={dim}, nthreads={nthreads})"
    assert int(np.count_nonzero(np.asarray(g))) == beads0, "bead count changed"


@pytest.mark.parametrize("dim", DIMS)
def test_parallel_pull_thread_count_independent(tmp_path, dim):
    box = [40, 40, 40] if dim == 3 else [40, 40]
    st = U.build_state(tmp_path, dim, "SLR", False, {"MOVE_PULL": 1.0},
                       box=box, chains=[(20, "AABBAB"), (20, "AAAAAA"), (15, "AB")])
    g1, t1, i1 = st.fresh()
    g4, t4, i4 = st.fresh()
    e1 = U.pull_parallel_megastep(st, g1, t1, i1, st.energy, 91, nthreads=1)
    e4 = U.pull_parallel_megastep(st, g4, t4, i4, st.energy, 91, nthreads=4)
    assert e1 == e4, f"energy differs across thread counts: {e1} vs {e4}"
    assert np.array_equal(np.asarray(g1), np.asarray(g4)), "grid differs across thread counts"
    assert np.array_equal(np.asarray(i1), np.asarray(i4)), "idx_to_bead differs across thread counts"


# FROZEN CHAINS under the parallel kernels: a frozen chain must never move (its
# beads stay exactly put and act as fixed obstacles), while the rest of the system
# still evolves and the tracked energy stays consistent.
_PARALLEL_DRIVERS = {
    "crank":   lambda st, *a, **k: (U.parallel_megastep if st.dim == 3 else U.parallel_megastep_2D)(st, *a, **k),
    "slither": U.slither_parallel_megastep,
    "pull":    U.pull_parallel_megastep,
}


@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("move", ["crank", "slither", "pull"])
def test_parallel_frozen_chains_do_not_move(tmp_path, dim, move):
    box = [40, 40, 40] if dim == 3 else [40, 40]
    # all chains length >= 3 so the system is valid for pull as well
    st = U.build_state(tmp_path, dim, "SLR", False, {"MOVE_CRANKSHAFT": 1.0},
                       box=box, chains=[(15, "AABBAB"), (15, "AAAAAA"), (12, "ABA")])
    g, t, i = st.fresh()
    idx0 = np.asarray(i)
    chain_ids = np.unique(idx0[:, 4])
    frozen = set(int(c) for c in chain_ids[::2])          # freeze every other chain

    frozen_rows = np.isin(idx0[:, 4], list(frozen))
    pos_before = idx0[frozen_rows][:, 5:5 + dim].copy()
    grid_before = np.asarray(g).copy()

    driver = _PARALLEL_DRIVERS[move]
    substeps = {"crank": 8000, "slither": 20, "pull": 20}[move]
    e = st.energy
    for m in range(4):
        e = driver(st, g, t, i, e, 123 + m, nthreads=4, frozen=frozen, substeps=substeps)
        assert e == U.recompute_energy(st, g, t, i), \
            f"energy drift with frozen chains ({move} {dim}D) at megastep {m}"

    idx_after = np.asarray(i)
    pos_after = idx_after[np.isin(idx_after[:, 4], list(frozen))][:, 5:5 + dim]
    assert np.array_equal(pos_before, pos_after), \
        f"a frozen chain moved under the parallel {move} kernel ({dim}D)"
    # the run must not be vacuous: the non-frozen part of the system moved
    assert np.any(grid_before != np.asarray(g)), "nothing moved at all (vacuous test)"


def test_parallelize_with_freeze_file_e2e(tmp_path):
    """End-to-end: PARALLELIZE together with a FREEZE_FILE, through the real
    Simulation, exercising all three parallel moves (crankshaft/slither/pull). The
    run must stay energy-consistent (ENERGY_CHECK raises otherwise) and the frozen
    chains must end exactly where they started."""
    box = [40, 40]
    chains = [(8, "AABBAB"), (8, "AAAAAA")]
    moves = {"MOVE_CRANKSHAFT": 0.34, "MOVE_SLITHER": 0.33, "MOVE_PULL": 0.33}
    base = {"ENERGY_CHECK": 4, "PRINT_FREQ": 10**9, "XTC_FREQ": 10**9,
            "ANALYSIS_FREQ": 10**9, "RESTART_FREQ": 10**9, "EN_FREQ": 10,
            "PARALLELIZE": True, "PARALLEL_THREADS": 4}
    U.write_param_file(os.path.join(str(tmp_path), "params.prm"), "SLR")

    def write(extra):
        U.write_keyfile(os.path.join(str(tmp_path), "KEYFILE.kf"), 2, False, moves,
                        box=box, chains=chains, n_steps=48, equilibration=8, seed=11,
                        extra=extra)

    # construct once (no freeze file) just to discover the chainIDs
    write(base)
    with U._chdir(str(tmp_path)), contextlib.redirect_stdout(open(os.devnull, "w")):
        sim0 = U.Simulation(U.KeyFileParser("KEYFILE.kf").keyword_lookup)
        chain_ids = sorted(sim0.LATTICE.chains.keys())
    frozen_ids = chain_ids[:5]
    with open(os.path.join(str(tmp_path), "freeze.txt"), "w") as fh:
        fh.write("C " + " ".join(str(c) for c in frozen_ids) + "\n")

    # real run with the freeze file + PARALLELIZE
    base_fz = dict(base, FREEZE_FILE="freeze.txt")
    write(base_fz)

    def snapshot(sim):
        return {c: [tuple(int(x) for x in p)
                    for p in sim.LATTICE.chains[c].get_ordered_positions()]
                for c in frozen_ids}

    with U._chdir(str(tmp_path)), contextlib.redirect_stdout(open(os.devnull, "w")):
        sim = U.Simulation(U.KeyFileParser("KEYFILE.kf").keyword_lookup)
        assert set(sim.frozen_chains) == set(frozen_ids), "freeze file not applied"
        before = snapshot(sim)
        sim.run_simulation()                 # raises if the tracked energy drifts
        after = snapshot(sim)

    assert before == after, "a frozen chain moved during a PARALLELIZE run"


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
# PARALLELIZE keyword end-to-end (2D and 3D now both have a parallel checkerboard
# kernel). A short run with ENERGY_CHECK must stay energy-consistent whichever
# kernel (parallel or serial fallback) is selected.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dim,ff,hardwall", ALL_CASES, ids=CASE_IDS)
def test_parallelize_simulation_energy_consistency(tmp_path, dim, ff, hardwall):
    U.run_sim_with_energy_check(tmp_path, dim, ff, hardwall,
                                {"MOVE_CRANKSHAFT": 1.0}, n_steps=40, energy_check=5,
                                extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 4})
