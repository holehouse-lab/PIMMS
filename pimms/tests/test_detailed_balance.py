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

Everything is seeded -> deterministic, so these do not flake. Deterministic is
not the same as robust, though: a fixture that passes on the shipped seed and
fails on the next one breaks the first time an unrelated change shifts a random
stream. The pull fixtures were in that state (each failed about one seed in ten,
see their comments) and are now tuned so that they pass, and their positive
controls bite, on eight other seed offsets as well.

Every kernel-level fixture carries a POSITIVE CONTROL: the same move re-run with
its inverse temperature scaled by a fixture-specific factor, which the fixture's
own assertion must FAIL. Without that, an equilibrium test degrades silently into
a vacuous green - which is where these fixtures were. Measured with the shipped
trace lengths, no kernel-level fixture could resolve a 20 % error in the
acceptance of the move it was testing, and the pull and parallel-slither fixtures
could not resolve a 50 % error either, because the tolerance was built partly
from the SEM of the trace it was judging (a broken move drifts, its
autocorrelation time balloons, and it buys itself a pass with the very pathology
that shows it is broken). All of these now pass ``ref_only=True``, run longer
traces, assert the test trace's integrated autocorrelation time is small enough
for the comparison to mean anything, and compare the mean squared radius of
gyration as well as the energy - the mean energy is blind by construction to any
bias that leaves it unchanged.

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

# Relative floor for the <Rg^2> comparison. Rg^2 is a much narrower distribution
# than the energy, so 0.5 % of it is below the level at which two independent
# equilibrium runs of the same system routinely differ; 2 % gives these fixtures a
# resolution of roughly 2.5-3.5 % of <Rg^2>, which is what they can honestly claim.
RG_FLOOR = 0.02

# ---------------------------------------------------------------------------
# slither (reptation) megamove - 2D and 3D, hardwall + PBC, full SLR forcefield.
# Control: invtemp x 1.2 on the slither only. Measured control margins (the
# factor by which the control run exceeds the tolerance) 2.1-3.2, so this fixture
# genuinely resolves a 20 % acceptance error; at the old sample of 120 with the
# two-SEM tolerance it resolved none of the five seeds the review tried.
# ---------------------------------------------------------------------------
SLITHER_CONTROL = 1.2

# ---------------------------------------------------------------------------
# pull (cooperative reptation) megamove - 2D and 3D, hardwall + PBC, full SLR.
#
# The pull is not ergodic on its own (it cannot translate whole chains or move
# single-bead chains), so detailed balance is tested the way the move is actually
# used: the pull as the DOMINANT move with a little crankshaft for ergodicity must
# reach the SAME Boltzmann equilibrium as crankshaft alone. A pull that violated
# detailed balance would continuously bias the chain and shift that equilibrium
# away from the crankshaft reference.
#
# The composite step was 30 pull substeps per chain followed by 400 crankshaft
# substeps, i.e. four crankshaft sweeps of relaxation after every pull burst, and
# the trusted move simply re-equilibrated whatever the pull injected: an invtemp
# x 1.5 pull shifted the mean energy by 1.6 % against a 0.5 % relative floor. 60
# pull substeps against 100 crankshaft substeps keeps the composite chain mixing
# (tau of the test trace stays below 8) while leaving the injected bias visible.
#
# Even so, this fixture cannot resolve a 20 % error: at invtemp x 1.2 the measured
# control margins are 0.15 (2D PBC) to 2.1 (3D HW), i.e. below 1 in three of the
# four cases, because most accepted pulls are close to energy-neutral and the bias
# is comparable to the relative floor. The control is therefore set at 1.5, and a
# 20 % Metropolis error in the pull remains outside what an equilibrium comparison
# on this system can see - the forward/reverse transition counting of single pull
# and slither sub-moves in test_whole_chain_transition_counts.py is the pin for
# that (its positive control rejects beta x 1.2), not this test.
#
# The tolerance is built from the crankshaft reference alone (tau ~0.7-0.9
# megamoves), but one composite step is a far weaker update than a 2500-substep
# crankshaft megamove: the pull-dominated trace has tau ~3.5 (2D) to ~7 (3D)
# composite steps, so with one step per sample its mean carried 2.2-2.8 times the
# standard error the tolerance assumes (measured over 264 and 80 disjoint
# fixture-sized windows: 9.7 against 4.4 in 2D hardwall, 19.7 against 7.0 in 3D
# periodic). A correct kernel then failed about one seed in ten, on the mean or
# on the stationarity guard - the shipped seed passed, three of the next eight
# did not. It is neither a bias nor slow equilibration: over 8 seeds the 2D
# hardwall test and reference means differ by -0.11 +- 0.54 at 40000 samples
# (0.01 % of |E|) and the 3D periodic ones by 0.05 % at 12000, with no transient
# in either trace, and the 2D hardwall pull kernel reproduces exactly enumerated
# Boltzmann distributions in a 7x7 box to within 0.05 % of the energy. Each
# sample is therefore taken every 8 composite steps (the parallel-slither
# remedy), which brings the test trace's tau to ~0.7 (2D) and ~1.1 (3D), level
# with the reference, and leaves it a pure pull-dominated trace. The 2D systems
# also run twice the samples: their relative floor is barely larger than the
# reference SEM, so the headroom of the control has to come from the trace
# length. Replayed on the long traces, the worst test/tolerance ratio is 0.57
# (2D hardwall, 16 windows) and 0.21 (3D periodic, 8 windows).
# ---------------------------------------------------------------------------
PULL_CONTROL = 1.5
PULL_SWEEPS_PER_SAMPLE = 8
PULL_SAMPLES = {2: 2400, 3: 1200}

# ---------------------------------------------------------------------------
# parallel checkerboard kernel - 3D, hardwall + PBC, full SLR forcefield.
# A dispersed box so the domain decomposition forms multiple blocks (the regime
# where the historical frozen-halo detailed-balance bug appeared). NB: the box
# must actually split for an LR system - the old 30^3 box was a SINGLE block
# under the LR halo, so this test never exercised the halo logic it describes.
# The layout is asserted so the test can never silently become vacuous again.
#
# This fixture is limited by the 0.5 % relative floor rather than by its trace
# length: an invtemp x 1.2 parallel crank shifts the mean energy by 0.92 % of |E|,
# so the control margin saturates near 1.9 however long the traces are (measured
# 1.32 at sample 1000). Lengthening it further buys almost nothing.
# ---------------------------------------------------------------------------
PARALLEL_CONTROL = 1.2

# ---------------------------------------------------------------------------
# parallel SLITHER kernel (mega_slither_parallel / _2D), 2D + 3D. A box large
# enough to decompose into multiple blocks with the (compact) chains fitting
# inside block interiors, so chains are distributed across blocks and the
# chain-level frozen-halo handling is exercised; must reach the same equilibrium
# as the serial crankshaft.
#
# The 3D fixture used to be 210 beads in a 40^3 box. With the chain-level halo
# (W = 5) that box splits into 20-site blocks with a 10-site interior, so only
# about 4 of its 66 chains were eligible in any sweep and they sat in two of the
# eight blocks: a valid but thin exercise of the block-parallel path. It is now
# a 48^3 box (24-site blocks, 14-site interior) at 2 % occupancy, where roughly a
# fifth of the 592 chains move per sweep, spread over all eight blocks. The
# control margin here is set by the relative-energy floor, not the trace length
# (1500 samples at the old density gained nothing), and it rises with density:
# 1.0-1.2 at 0.3 %, 1.4 at 0.7 %, 1.8-1.9 at 1.3 %, 2.1 at 2 %. Measured control
# margins at invtemp x 1.2: 2.1 (3D), 1.5-3.2 (2D, unchanged).
#
# The tolerance is built from the crankshaft reference alone (tau ~0.6 megamoves),
# but a parallel slither sweep only moves the chains inside the block interiors of
# that sweep's shift, so the parallel trace decorrelates per SWEEP, not per
# sub-move: its tau is ~5 sweeps whatever the substep count. In 2D, where the
# relative floor is small, that left the test mean with two to three times the
# error the tolerance assumes, and a correct kernel failed about one seed in four
# (long runs of 3000 sweeps put the parallel mean within 1 SEM of the reference,
# with the sign of the offset varying between seeds, in both geometries). Each 2D
# sample is therefore taken every 8 sweeps, which brings its tau to ~0.8 and keeps
# the trace a pure parallel-slither trace: over 8 seeds per geometry the worst
# test/tolerance ratio is 0.5 and the control margins are 1.7-3.0. The 3D fixture
# (tau ~6, but a much larger relative floor) stays at one sweep per sample, with a
# worst test/tolerance ratio of 0.23 over 8 seeds.
# ---------------------------------------------------------------------------
PARALLEL_SLITHER_CONTROL = 1.2
PARALLEL_SLITHER_SWEEPS_PER_SAMPLE = {2: 8, 3: 1}

# ---------------------------------------------------------------------------
# parallel PULL kernel (mega_pull_parallel / _2D), 2D + 3D, multi-block box. Pull
# rearranges sub-segments but does not translate chains freely, so (as for the
# serial pull DB test) the step mixes parallel pull with serial crankshaft for
# ergodicity; the parallel pull must not bias the crankshaft equilibrium.
#
# The 3D system was 342 beads in a 40^3 box - 0.5 % occupancy, so nearly every
# pull was energy-neutral and an invtemp x 1.5 pull shifted the mean energy by
# less than half the tolerance (measured control margins 0.28 and 0.45: the
# fixture had no power at all). It then became a 32^3 box at 2 % occupancy, run
# warm enough (T = 85) that the dense system does not condense - at T = 55 the
# same box phase-separates and the crankshaft REFERENCE picks up an integrated
# autocorrelation time of 70+ megamoves, which destroys the comparison from the
# other side.
#
# That 32^3 box, however, only decomposed for the CRANKSHAFT layout (halo 2,
# eight blocks), which is what the guard asserted. The pull and slither kernels
# use the wider chain-level halo (W = R_int + 2 = 5 with SLR interactions) and
# blocks of at least 4W sites, so for them 32 // 20 = 1 block per axis: the 3D
# parallel pull ran single-block, with no halo, no interior-restricted targets
# and one random stream, and its multi-block detailed balance was never tested.
# The 3D box is now 48^3 (2 x 2 x 2 blocks of 24, interior 14 sites, so every
# 6-mer fits) at the same 2 % occupancy and temperature, and the guard asserts
# the PULL layout. Only the chains inside a block interior move in a sweep,
# about a fifth of them, so the pull carries less weight against the crankshaft
# relaxation than it did single-block; scaling the relaxation with the bead
# count left the control at 0.8-1.1. The pull therefore gets 120 substeps per
# chain against a 300-substep relaxation, which measured control margins of
# 1.8 and 2.4 at invtemp x 1.5 (2D, unchanged: 1.7-2.2) with an integrated
# autocorrelation time of 15-20 megamoves, hence the 900-sample trace.
#
# That set the fixture's power, not its robustness to the seed. On eight other
# seeds the two 2D fixtures failed three runs of sixteen - twice on the positive
# control (margins 0.83 and 0.88), once on the mean - and the 3D periodic one
# tripped the autocorrelation guard once (tau 49 against a limit of 45). None of
# it is a bias: over 8 seeds and 12000 samples the 2D hardwall parallel-pull and
# reference means differ by +0.8 +- 4.3 (0.02 % of |E|). The causes:
#   * one composite step per sample left the test and control traces with tau
#     ~6-7 samples in 2D against a tolerance built from the reference alone, so
#     their means carried 1.7 times the reference's standard error;
#   * the 2D reference itself mixes slowly at T = 55: its 4500-substep megamove
#     has an effective tau of ~3 megamoves, which the autocorrelation estimator
#     (it truncates at the first negative lag) reads as 1.1-2.7, so the
#     reference SEM, and with it the tolerance, came out 1.5 times too small -
#     11.1 estimated against 16.5 measured over 80 fixture-sized windows;
#   * with both of those, 4 of 32 fixture-sized windows of the control fell
#     below a margin of 1.
# The 2D fixture therefore records a sample every 6 composite steps, uses an
# 18000-substep reference megamove (the equilibration lengthens with it) and
# runs 2400 samples. The 3D fixture's tolerance is dominated by the relative
# floor (0.5 % of |E| ~ 90000), so its mean was never at risk; its weak point is
# the autocorrelation guard, because the dense 3D system has a slow mode in both
# traces (the reference's own tau is 3-9 megamoves). A sample every 2nd
# composite step took the test trace's tau to 4-28 samples over eight other
# seeds at 900 samples, every case passing with control margins of 1.4-3.2 -
# but to 43.7, against a limit of 45, on the shipped seed once the kernels'
# random streams changed. The 3D traces are therefore 1800 samples as well,
# which doubles that limit and shortens the reference SEM the control is
# measured against. The 3D control is still the thinnest headroom in this file
# and these two cases are its longest-running.
# ---------------------------------------------------------------------------
PARALLEL_PULL_CONTROL = 1.5
PARALLEL_PULL_SWEEPS_PER_SAMPLE = {2: 6, 3: 2}
PARALLEL_PULL_SAMPLES = {2: 2400, 3: 1800}
PARALLEL_PULL_REFERENCE_SUBSTEPS = {2: 18000, 3: 4500}

# ---------------------------------------------------------------------------
# The four Python single-chain moves (codes 2-5: chain translate, chain rotate,
# chain pivot, head pivot). None of them had an equilibrium fixture at all, in
# either dimensionality, although all four are enabled in the shipped demo
# keyfiles. crank + the four of them must reach the same Boltzmann equilibrium as
# crank alone.
#
# What this can and cannot see. It compares the mean energy of two full
# simulations, and its positive control (the same move set run at TEMPERATURE / k)
# is resolved at k = 1.2. That is a real but blunt instrument: a rigid move with a
# one-way proposal shifts no energy at all, and an unbounded proposal asymmetry in
# the pivot was measured to move the mean Rg of a lone 14-mer by well under 1 % once
# crankshaft runs alongside it. The primary evidence for these four moves is
# therefore the free-draw forward/reverse transition counting in
# test_proposal_symmetry.py, which kills that class of bug deterministically; this
# fixture guards the acceptance criterion and the energy bookkeeping around them.
# ---------------------------------------------------------------------------
PYTHON_MOVES = {"MOVE_CHAIN_TRANSLATE": 0.2, "MOVE_CHAIN_ROTATE": 0.2,
                "MOVE_CHAIN_PIVOT": 0.2, "MOVE_HEAD_PIVOT": 0.2}
PYTHON_MOVE_CONTROL = 1.2


def _assert_detailed_balance(res, label):
    """Assert one kernel-level detailed-balance comparison, controls included.

    Four things are checked, in the order that localises a failure best: the test
    trace carries enough independent samples for any of this to mean something;
    the mean energy matches the crankshaft reference; the mean squared radius of
    gyration matches it too (a conformational observable, blind to the energy);
    and the deliberately mis-tempered control run does NOT match, which is what
    establishes that the first three assertions had power.

    Parameters
    ----------
    res : kernel_test_utils.DBResult
        The three traces from :func:`kernel_test_utils.db_compare_with_control`.

    label : str
        Fixture name used in every failure message.
    """
    U.assert_trace_is_well_sampled(res.test_energy, label)
    U.assert_same_equilibrium(res.ref_energy, res.test_energy, label, ref_only=True)
    U.assert_same_equilibrium(res.ref_rg2, res.test_rg2, f"{label} <Rg^2>",
                              ref_only=True, rel_floor=RG_FLOOR)
    U.assert_control_is_resolved(res.ref_energy, res.control_energy, label,
                                 res.control_scale)


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_slither_detailed_balance(tmp_path, dim, hardwall):
    st = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_SLITHER": 1.0})

    def slither_step(state, g, t, i, e, seed):
        return U.slither_megastep(state, g, t, i, e, seed)

    def control_step(state, g, t, i, e, seed):
        with U.scaled_invtemp(state, SLITHER_CONTROL):
            return U.slither_megastep(state, g, t, i, e, seed)

    res = U.db_compare_with_control(st, slither_step, control_step, equilibrate=120,
                                    sample=900, control_scale=SLITHER_CONTROL)
    _assert_detailed_balance(res, f"slither {dim}D {HW_IDS[hardwall]} SLR")


def _pull_dominant_step(state, g, t, idx, energy, seed, *, scale=1.0):
    """One sample of the pull-dominated composite move.

    Runs ``PULL_SWEEPS_PER_SAMPLE`` composite steps, each a burst of 60 pull
    substeps per chain followed by 100 crankshaft substeps for ergodicity, and
    returns the energy after the last one.

    Parameters
    ----------
    state : kernel_test_utils.State
        The system, as returned by :func:`kernel_test_utils.build_state`.

    g, t, idx : numpy.ndarray
        Occupancy grid, type grid and bead table, mutated in place.

    energy : int
        Total energy on entry.

    seed : int
        Seed of this sample; every composite step inside it gets its own
        derived seed, distinct from those of every other sample.

    scale : float, optional
        Multiplier applied to the inverse temperature the PULL is handed (the
        crankshaft always runs at the true one). ``1.0`` is the move under
        test, anything else is the positive control.

    Returns
    -------
    int
        Total energy after the composite steps.
    """
    e = energy
    for sweep in range(PULL_SWEEPS_PER_SAMPLE):
        sweep_seed = seed * PULL_SWEEPS_PER_SAMPLE + sweep
        with U.scaled_invtemp(state, scale):
            e, _ = U.pull_megastep(state, g, t, idx, e, sweep_seed, substeps=60)
        e = U.crank_megastep(state, g, t, idx, e, sweep_seed + 777, substeps=100)
    return e


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_pull_detailed_balance(tmp_path, dim, hardwall):
    st = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0})

    def control_step(state, g, t, i, e, seed):
        return _pull_dominant_step(state, g, t, i, e, seed, scale=PULL_CONTROL)

    res = U.db_compare_with_control(st, _pull_dominant_step, control_step,
                                    equilibrate=150, sample=PULL_SAMPLES[dim],
                                    control_scale=PULL_CONTROL)
    _assert_detailed_balance(res, f"pull {dim}D {HW_IDS[hardwall]} SLR")


@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_parallel_detailed_balance(tmp_path, hardwall):
    from pimms import mega_crank_fast
    box = [40, 40, 40]
    assert mega_crank_fast.parallel_crank_layout_info(*box, True)["num_blocks"] > 1
    st = U.build_state(tmp_path, 3, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                       box=box,
                       chains=[(40, "AABB"), (40, "AAAA"), (30, "A")])

    def parallel_step(state, g, t, i, e, seed, *, scale=1.0):
        with U.scaled_invtemp(state, scale):
            return U.parallel_megastep(state, g, t, i, e, seed, nthreads=4)

    res = U.db_compare_with_control(
        st, parallel_step,
        lambda s, g, t, i, e, sd: parallel_step(s, g, t, i, e, sd, scale=PARALLEL_CONTROL),
        equilibrate=110, sample=1000, crank_substeps=4500,
        control_scale=PARALLEL_CONTROL)
    _assert_detailed_balance(res, f"parallel 3D {HW_IDS[hardwall]} SLR")


# ---------------------------------------------------------------------------
# parallel checkerboard kernel - 2D (mega_crank_parallel_2D). A dispersed 2D box
# that decomposes into multiple blocks, so the frozen-halo handling is exercised;
# the 2D parallel kernel must reach the same equilibrium as the serial crankshaft.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_parallel_2D_detailed_balance(tmp_path, hardwall):
    from pimms import mega_crank_fast
    assert mega_crank_fast.parallel_crank_layout_info(40, 40, 1, True)["num_blocks"] > 1
    st = U.build_state(tmp_path, 2, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                       box=[40, 40],
                       chains=[(22, "AABB"), (22, "AAAA"), (18, "A")])

    def parallel_step(state, g, t, i, e, seed, *, scale=1.0):
        with U.scaled_invtemp(state, scale):
            return U.parallel_megastep_2D(state, g, t, i, e, seed, nthreads=4)

    res = U.db_compare_with_control(
        st, parallel_step,
        lambda s, g, t, i, e, sd: parallel_step(s, g, t, i, e, sd, scale=PARALLEL_CONTROL),
        equilibrate=110, sample=600, crank_substeps=4500,
        control_scale=PARALLEL_CONTROL)
    _assert_detailed_balance(res, f"parallel 2D {HW_IDS[hardwall]} SLR")


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_parallel_slither_detailed_balance(tmp_path, dim, hardwall):
    from pimms import mega_crank_fast
    if dim == 3:
        box, chains, sample = [48, 48, 48], [(216, "AABB"), (216, "AAAA"), (160, "A")], 1000
    else:
        box, chains, sample = [40, 40], [(24, "AABB"), (24, "AAAA"), (18, "A")], 600
    # the chain-level layout the slither kernel actually uses, not the crankshaft's
    layout = mega_crank_fast.parallel_layout_info(
        box[0], box[1], box[2] if dim == 3 else 1, True)
    assert layout["num_blocks"] > 1, layout
    assert max(len(seq) for _, seq in chains) <= min(
        layout["block_size"][d] - 2 * layout["W"]
        for d in range(dim) if layout["blocks"][d] > 1), layout
    st = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                       box=box, chains=chains)

    sweeps = PARALLEL_SLITHER_SWEEPS_PER_SAMPLE[dim]

    def slither_step(state, g, t, i, e, seed, *, scale=1.0):
        with U.scaled_invtemp(state, scale):
            for sweep in range(sweeps):
                e = U.slither_parallel_megastep(state, g, t, i, e, seed * sweeps + sweep,
                                                substeps=60, nthreads=4)
        return e

    res = U.db_compare_with_control(
        st, slither_step,
        lambda s, g, t, i, e, sd: slither_step(s, g, t, i, e, sd,
                                               scale=PARALLEL_SLITHER_CONTROL),
        equilibrate=140, sample=sample, crank_substeps=4500,
        control_scale=PARALLEL_SLITHER_CONTROL)
    _assert_detailed_balance(res, f"parallel slither {dim}D {HW_IDS[hardwall]} SLR")


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_parallel_pull_detailed_balance(tmp_path, dim, hardwall):
    from pimms import mega_crank_fast
    if dim == 3:
        box, chains, temperature = [48, 48, 48], \
            [(152, "AABBAB"), (152, "AAAAAA"), (152, "ABA")], 85
        relax, substeps = 300, 120
    else:
        box, chains, temperature = [40, 40], \
            [(24, "AABBAB"), (24, "AAAAAA"), (18, "ABA")], 55
        relax, substeps = 342, 60
    sweeps = PARALLEL_PULL_SWEEPS_PER_SAMPLE[dim]
    # the chain-level layout the pull kernel actually uses, not the crankshaft's
    layout = mega_crank_fast.parallel_layout_info(
        box[0], box[1], box[2] if dim == 3 else 1, True)
    assert layout["num_blocks"] > 1, layout
    assert max(len(seq) for _, seq in chains) <= min(
        layout["block_size"][d] - 2 * layout["W"]
        for d in range(dim) if layout["blocks"][d] > 1), layout
    st = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                       box=box, chains=chains, temperature=temperature)

    def pull_step(state, g, t, i, e, seed, *, scale=1.0):
        for sweep in range(sweeps):
            sweep_seed = seed * sweeps + sweep
            with U.scaled_invtemp(state, scale):
                e = U.pull_parallel_megastep(state, g, t, i, e, sweep_seed, substeps=substeps,
                                             nthreads=4)
            e = U.crank_megastep(state, g, t, i, e, sweep_seed + 777, substeps=relax)
        return e

    res = U.db_compare_with_control(
        st, pull_step,
        lambda s, g, t, i, e, sd: pull_step(s, g, t, i, e, sd, scale=PARALLEL_PULL_CONTROL),
        equilibrate=150, sample=PARALLEL_PULL_SAMPLES[dim],
        crank_substeps=PARALLEL_PULL_REFERENCE_SUBSTEPS[dim],
        control_scale=PARALLEL_PULL_CONTROL)
    _assert_detailed_balance(res, f"parallel pull {dim}D {HW_IDS[hardwall]} SLR")


def _python_move_equilibrium(tmp_path, sub, moves, *, dim, hardwall, seed, n_steps,
                             temperature):
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
    box = [14, 14, 14] if dim == 3 else [20, 20]
    chains = [(10, "AABB"), (10, "AAAA"), (6, "A")]
    U.write_param_file(os.path.join(str(d), "params.prm"), "SLR")
    U.write_keyfile(os.path.join(str(d), "KEYFILE.kf"), dim, hardwall, moves,
                    box=box, chains=chains, temperature=temperature,
                    n_steps=n_steps, equilibration=n_steps // 4, seed=seed, extra=extra)
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
    return e[len(e) // 2:], sim


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", HARDWALLS, ids=[HW_IDS[h] for h in HARDWALLS])
def test_python_single_chain_moves_detailed_balance(tmp_path, dim, hardwall):
    n_steps = 8000
    temperature = 45
    mixed = dict(PYTHON_MOVES, MOVE_CRANKSHAFT=0.2)
    ref, _ = _python_move_equilibrium(
        tmp_path, "ref", {"MOVE_CRANKSHAFT": 1.0}, dim=dim, hardwall=hardwall,
        seed=21, n_steps=n_steps, temperature=temperature)
    test, sim = _python_move_equilibrium(
        tmp_path, "test", mixed, dim=dim, hardwall=hardwall, seed=21,
        n_steps=n_steps, temperature=temperature)

    # only meaningful if all four moves actually fire
    for code, name in ((2, "chain translate"), (3, "chain rotate"),
                       (4, "chain pivot"), (5, "head pivot")):
        assert sim.ACC.accepted_count[code] > 20, (
            f"{name} was accepted only {sim.ACC.accepted_count[code]} times - "
            f"this fixture is not exercising it")

    label = f"python single-chain moves {dim}D {HW_IDS[hardwall]}"
    U.assert_same_equilibrium(ref, test, label, ref_only=True)

    # positive control: the SAME move mix at a 20 % colder temperature must fail.
    # A Simulation-level fixture cannot mis-temper one move in isolation, so this
    # measures the fixture's resolving power on the mean energy rather than
    # isolating the Python moves - which is the honest claim it can make.
    control, _ = _python_move_equilibrium(
        tmp_path, "control", mixed, dim=dim, hardwall=hardwall, seed=21,
        n_steps=n_steps, temperature=temperature / PYTHON_MOVE_CONTROL)
    U.assert_control_is_resolved(ref, control, label, PYTHON_MOVE_CONTROL)


# ---------------------------------------------------------------------------
# TSMMC moves - coordinated by the Simulation, so compared via two full runs
# (crankshaft-only reference vs TSMMC+crankshaft) reaching equilibrium.
# ---------------------------------------------------------------------------
def _sim_equilibrium(tmp_path, sub, dim, ff, hardwall, moves, *, n_steps, seed, return_sim=False):
    """Run a short simulation and return the 2nd half of its ENERGY.dat trace.

    With ``return_sim`` the finished Simulation is returned as well, so a caller
    can check that the move under test actually did something.
    """
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
        sim = Simulation(keyfile.keyword_lookup)
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim.run_simulation()
        e = np.loadtxt("ENERGY.dat", delimiter="\t")[:, 1]
    finally:
        os.chdir(cwd)
    if return_sim:
        return e[len(e) // 2:], sim
    return e[len(e) // 2:]


@pytest.mark.parametrize("tsmmc_move", ["MOVE_CTSMMC", "MOVE_SYSTEM_TSMMC", "MOVE_MULTICHAIN_TSMMC"])
def test_tsmmc_detailed_balance(tmp_path, tsmmc_move):
    # 3D, short-range, PBC. Reference is crankshaft only; the test run mixes the
    # TSMMC move with crankshaft. TSMMC is a temperature-excursion move, so a
    # broken acceptance biases the sampled energy well away from the reference.
    n_steps = 4000
    ref = _sim_equilibrium(
        tmp_path, "ref", 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
        n_steps=n_steps, seed=21)
    test, sim = _sim_equilibrium(
        tmp_path, "test", 3, "SR", False,
        {"MOVE_CRANKSHAFT": 0.7, tsmmc_move: 0.3}, n_steps=n_steps, seed=21,
        return_sim=True)

    # a TSMMC move that never accepts is a no-op and matches the crankshaft
    # reference trivially, so the excursions must really have been taken (the
    # main-chain counters hold one entry per excursion, whatever the variant)
    code = {"MOVE_CTSMMC": 9, "MOVE_MULTICHAIN_TSMMC": 10, "MOVE_SYSTEM_TSMMC": 12}[tsmmc_move]
    assert sim.ACC.accepted_count[code] >= 20, (
        "%s accepted only %d of %d excursions - too few to test anything"
        % (tsmmc_move, sim.ACC.accepted_count[code], sim.ACC.move_count[code]))
    assert sim.ACC.accepted_count[code] < sim.ACC.move_count[code], (
        "%s accepted every excursion, so the acceptance test was never exercised" % tsmmc_move)

    # reference-only autocorrelation SEM + 0.5 % floor, with a stationarity
    # guard on the test trace: the old 2.5 * std + 5 % criterion let a beta x1.5
    # acceptance error pass every one of these cases
    U.assert_same_equilibrium(ref, test, tsmmc_move, ref_only=True)


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
    return e[len(e) // 2:], sim


def test_vmmc_detailed_balance(tmp_path):
    n_steps = 4000
    ref, _ = _vmmc_equilibrium(
        tmp_path, "ref", {"MOVE_CRANKSHAFT": 1.0}, seed=21, n_steps=n_steps)
    test, sim = _vmmc_equilibrium(
        tmp_path, "test", {"MOVE_CRANKSHAFT": 0.5, "MOVE_VMMC": 0.5}, seed=21, n_steps=n_steps)

    # the test is only meaningful if VMMC actually performed collective moves
    assert sim.vmmc_accepted_multichain > 5, (
        f"VMMC accepted too few multi-chain moves ({sim.vmmc_accepted_multichain}) - "
        f"the system is not exercising the collective move, so this is not a real test")
    assert sim.vmmc_max_accepted_cluster >= 2

    U.assert_same_equilibrium(
        ref, test,
        f"VMMC (accepted {sim.vmmc_accepted_multichain} multi-chain moves, max cluster "
        f"{sim.vmmc_max_accepted_cluster})", ref_only=True)


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
    return e[len(e) // 2:], sim


def test_jump_and_relax_detailed_balance(tmp_path):
    n_steps = 6000
    ref, _ = _jr_equilibrium(
        tmp_path, "ref", {"MOVE_CRANKSHAFT": 1.0}, seed=21, n_steps=n_steps)
    test, sim = _jr_equilibrium(
        tmp_path, "test", {"MOVE_CRANKSHAFT": 0.5, "MOVE_JUMP_AND_RELAX": 0.5}, seed=21, n_steps=n_steps)

    # only meaningful if the jump actually fires
    assert sim.ACC.accepted_count[13] > 0, "jump-and-relax accepted no jumps - not a real test"

    U.assert_same_equilibrium(ref, test, "jump-and-relax", ref_only=True)


# ---------------------------------------------------------------------------
# Cluster translation / rotation (codes 7 and 8). Both are energy-neutral rigid
# moves whose correctness rests entirely on proposal symmetry; neither had an
# equilibrium test. The complete review found a directional spurious clash in
# cluster_translate (the cluster was deleted chain by chain) and a one-way
# refusal in hardwall cluster_rotate on non-cubic boxes (the periodic winding
# guard). crank + cluster moves must reach the same equilibrium as crank alone,
# under periodic boundaries (cubic box) and under a hardwall (non-cubic box).
# ---------------------------------------------------------------------------
def _cluster_equilibrium(tmp_path, sub, moves, *, hardwall, box, seed, n_steps):
    d = tmp_path / sub
    d.mkdir()
    extra = {
        "EN_FREQ": 10,
        "PRINT_FREQ": 1000000,
        "XTC_FREQ": 1000000,
        "ANALYSIS_FREQ": 1000000,
        "RESTART_FREQ": 1000000,
    }
    U.write_param_file(os.path.join(str(d), "params.prm"), "SR")
    U.write_keyfile(os.path.join(str(d), "KEYFILE.kf"), 3, hardwall, moves,
                    box=box, chains=[(6, "AABB"), (6, "AAAA"), (8, "A")],
                    temperature=30, n_steps=n_steps, equilibration=n_steps // 4,
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
    return e[len(e) // 2:], sim


@pytest.mark.parametrize("hardwall,box", [(False, [12, 12, 12]), (True, [10, 12, 14])],
                         ids=["pbc-cubic", "hardwall-noncubic"])
def test_cluster_moves_detailed_balance(tmp_path, hardwall, box):
    """An equilibrium check with modest power (it resolves a ~1 % shift of the
    mean energy); a purely directional proposal bias in an energy-neutral move
    shifts no energy at all, so the direct proposal-symmetry pins in
    test_full_review_fixes.py are the primary evidence for the two cluster-move
    fixes. This test guards the energy bookkeeping around them."""
    n_steps = 4000
    ref, _ = _cluster_equilibrium(
        tmp_path, "ref", {"MOVE_CRANKSHAFT": 1.0}, hardwall=hardwall, box=box, seed=21, n_steps=n_steps)
    test, sim = _cluster_equilibrium(
        tmp_path, "test",
        {"MOVE_CRANKSHAFT": 0.6, "MOVE_CLUSTER_TRANSLATE": 0.2, "MOVE_CLUSTER_ROTATE": 0.2},
        hardwall=hardwall, box=box, seed=21, n_steps=n_steps)

    # only meaningful if both cluster moves actually fire
    assert sim.ACC.accepted_count[7] > 5, sim.ACC.accepted_count[7]
    assert sim.ACC.accepted_count[8] > 5, sim.ACC.accepted_count[8]

    U.assert_same_equilibrium(ref, test, f"cluster moves ({'hardwall' if hardwall else 'pbc'})", ref_only=True)
