"""Regression tests for the compiled-kernel findings of the deep audit.

Five things are pinned here, each against a value worked out from its
definition in this file rather than against the kernel's own output.

**The serial Metropolis test (M3-1 / E2-1).** The serial kernels used to accept
an uphill move iff ``float32(r / (2^31 - 1)) < float32(exp(-beta dE))`` with r a
31-bit draw. r == 0 passes that for any positive Boltzmann factor, so every
uphill move was accepted at least 2^-31 = 4.66e-10 of the time: 33 times too
often at beta dE = 25 and 10^8 times too often at 40. The test now has a second
stage. When the first stage accepts and the factor lies inside the draw's own
bin ``[r, r + 1) / 2^31``, further draws place the draw inside that bin, so the
accepted probability is the float32 factor itself. Nothing else changes: the
first draw, the first comparison and every decision outside that one bin are
what they were, which is why the regression fixtures did not move.

**Entry guards (M3-2).** Bounds checks are off in the kernels, so a selector
that is too short or names a row that does not exist, a 2D table handed to a 3D
kernel, or chain arrays that have drifted apart, used to read and write the
wrong memory without a word. Every public kernel now raises ``ValueError``
before it touches anything. Each bad input below is built as a view of a larger
array, so the same call is memory-safe on a kernel without the guard as well.

**Thread clamp (M3-3).** The thread budget handed to OpenMP is limited to the
number of blocks.

**Block shift (M2-1, M2-2).** The parallel kernels draw their block-origin shift
from one generator that is reseeded per megamove, not from a generator
constructed per megamove, and the two crankshaft kernels draw it over the whole
box, as the whole-chain kernels do, so that on an axis with a remainder a site's
chance of being movable no longer depends on where it sits.
"""

import math
from fractions import Fraction
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import pytest

import pimms.mega_crank as ref_kernel
import pimms.mega_crank_2D as ref_kernel_2D
import pimms.mega_crank_fast as fk
from pimms.tests import kernel_test_utils as ktu

# ---------------------------------------------------------------------------
# splitmix64, written out from its definition so a test can choose a state
# ---------------------------------------------------------------------------

_MASK: int = (1 << 64) - 1
_GAMMA: int = 0x9E3779B97F4A7C15
_C1: int = 0xBF58476D1CE4E5B9
_C2: int = 0x94D049BB133111EB
_C1_INV: int = pow(_C1, -1, 1 << 64)
_C2_INV: int = pow(_C2, -1, 1 << 64)
_GAMMA_INV: int = pow(_GAMMA, -1, 1 << 64)
_PRNG_MAX: int = 2147483647
_TWO31: int = 1 << 31

_SERIAL_CRANK = {3: "mega_crank", 2: "mega_crank_2D"}
_PARALLEL_CRANK = {3: "mega_crank_parallel", 2: "mega_crank_parallel_2D"}
_SERIAL_CHAIN = {
    3: ["mega_slither", "mega_pull"],
    2: ["mega_slither_2D", "mega_pull_2D"],
}
_PARALLEL_CHAIN = {
    3: ["mega_slither_parallel", "mega_pull_parallel"],
    2: ["mega_slither_parallel_2D", "mega_pull_parallel_2D"],
}

_REMAINDER_BOX: int = 17
_REMAINDER_BETA: float = 0.05
_CONTROL_SCALE: float = 1.1


def _mix(z: int) -> int:
    """The splitmix64 output function.

    Parameters
    ----------
    z : int
        A generator state, already advanced by the increment.

    Returns
    -------
    int
        The 64-bit output for that state.
    """
    z = ((z ^ (z >> 30)) * _C1) & _MASK
    z = ((z ^ (z >> 27)) * _C2) & _MASK
    return z ^ (z >> 31)


def _undo_xorshift(z: int, shift: int) -> int:
    """Invert ``x ^ (x >> shift)`` on 64 bits.

    Parameters
    ----------
    z : int
        The value after the xor-shift.

    shift : int
        The shift distance.

    Returns
    -------
    int
        The x for which ``x ^ (x >> shift) == z``.
    """
    out = z
    for _ in range(64 // shift + 1):
        out = z ^ (out >> shift)
    return out


def _unmix(z: int) -> int:
    """Invert the splitmix64 output function.

    Parameters
    ----------
    z : int
        A 64-bit output.

    Returns
    -------
    int
        The (advanced) state that produces it.
    """
    z = _undo_xorshift(z, 31)
    z = (z * _C2_INV) & _MASK
    z = _undo_xorshift(z, 27)
    z = (z * _C1_INV) & _MASK
    return _undo_xorshift(z, 30)


def _state_for_draw(r: int, low: int, n: int = 1) -> int:
    """A generator state from which the n-th ``mc_rand()`` returns exactly r.

    ``mc_rand`` returns the top 31 bits of one splitmix64 step, so the other 33
    bits are free: `low` picks one of the 2^33 states that give the same r.

    Parameters
    ----------
    r : int
        The 31-bit draw wanted, in [0, 2^31 - 1].

    low : int
        The 33 bits of the output that ``mc_rand`` throws away.

    n : int, optional
        Which draw (1-based) should return r.

    Returns
    -------
    int
        The state (equivalently, the kernel seed) to start from.
    """
    output = (r << 33) | (low & ((1 << 33) - 1))
    return (_unmix(output) - n * _GAMMA) & _MASK


def _draws_between(before: int, after: int) -> int:
    """Number of generator steps that take the state from `before` to `after`.

    Parameters
    ----------
    before, after : int
        The two states.

    Returns
    -------
    int
        The step count; every step adds the same odd increment, so it is the
        difference times that increment's inverse modulo 2^64.
    """
    return ((after - before) * _GAMMA_INV) & _MASK


def _second_stage(state: int, q: Fraction) -> Tuple[int, int]:
    """Decision and draw count of the second stage, from its definition.

    The first stage has consumed one step. The second stage accepts iff a
    uniform ``U = 0.d1 d2 d3 ...`` written in base 2^32 is below q, where each
    digit is the top 32 bits of the next generator step; it stops at the first
    digit that differs from the matching digit of q.

    Parameters
    ----------
    state : int
        The state before the FIRST draw of the acceptance test.

    q : Fraction
        The probability the second stage has to realise, in (0, 1).

    Returns
    -------
    tuple of (int, int)
        (decision, total draws consumed including the first-stage draw).
    """
    st = (state + _GAMMA) & _MASK
    draws = 1
    while q > 0:
        st = (st + _GAMMA) & _MASK
        draws += 1
        digit = _mix(st) >> 32
        q_digit = math.floor(q * (1 << 32))
        if digit < q_digit:
            return 1, draws
        if digit > q_digit:
            return 0, draws
        q = q * (1 << 32) - q_digit
    return 0, draws


def _float32_factor(invtemp: float, d_energy: int) -> np.float32:
    """The float32 Boltzmann factor exactly as ``accept_or_reject`` forms it.

    Parameters
    ----------
    invtemp : float
        Inverse temperature; the kernels take it as a C float.

    d_energy : int
        The (positive) energy change.

    Returns
    -------
    numpy.float32
        ``float32(exp(float32(-dE) * float32(invtemp)))``.
    """
    exponent = np.float32(-d_energy) * np.float32(invtemp)
    return np.float32(math.exp(float(exponent)))


def _first_stage(r: int, factor: np.float32) -> int:
    """The historical single-stage test for a 31-bit draw r.

    Parameters
    ----------
    r : int
        The draw.

    factor : numpy.float32
        The float32 Boltzmann factor.

    Returns
    -------
    int
        1 iff ``float32(r / (2^31 - 1)) < factor``.
    """
    return 1 if np.float32(np.float64(r) / np.float64(_PRNG_MAX)) < factor else 0


def _expected_crank(state: int, invtemp: float, d_energy: int) -> Tuple[int, int]:
    """Decision and draw count ``accept_or_reject`` must give from `state`.

    Parameters
    ----------
    state : int
        Generator state before the call.

    invtemp : float
        Inverse temperature.

    d_energy : int
        The (positive) energy change.

    Returns
    -------
    tuple of (int, int)
        (decision, draws): the first stage as it always was, then the second
        stage only when the first accepted and the factor lies inside the
        draw's own bin.
    """
    factor = _float32_factor(invtemp, d_energy)
    r = _mix((state + _GAMMA) & _MASK) >> 33
    if not _first_stage(r, factor):
        return 0, 1
    scaled = Fraction(float(factor)) * _TWO31
    if r + 1 <= scaled:
        return 1, 1
    return _second_stage(state, scaled - r)


# ---------------------------------------------------------------------------
# M3-1 / E2-1: the serial acceptance test
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("beta_de", [17, 20, 21, 25, 30, 36, 60])
def test_serial_acceptance_is_exact_in_the_threshold_bin(beta_de: int) -> None:
    """The bin that holds the Boltzmann factor is split in the right proportion.

    For draws in the bin below the factor's bin the move is accepted on one
    draw, for draws in the bin above it is rejected on one draw, and inside the
    factor's bin the decision and the number of draws are those of the
    digit-by-digit comparison written out in ``_second_stage``. Both the fast
    and the reference kernel must agree with that, state by state. Summed over
    the bins this makes the accepted probability the float32 factor exactly,
    which is within float32 rounding of exp(-beta dE).

    Parameters
    ----------
    beta_de : int
        beta * dE, realised as invtemp = 1 and an integer dE.
    """
    factor = _float32_factor(1.0, beta_de)
    scaled = Fraction(float(factor)) * _TWO31
    r_star = math.floor(scaled)
    assert scaled != r_star, "pick a case whose factor is strictly inside a bin"
    # the float32 factor is the target; it sits within float32 rounding of the truth
    assert abs(float(factor) / math.exp(-beta_de) - 1.0) < 2.0**-23

    second_stage_accepts = 0
    second_stage_total = 0
    for r in (r_star - 1, r_star, r_star + 1):
        if r < 0:
            continue
        for low in range(1, 301):
            state = _state_for_draw(r, low * 7919 + r)
            want = _expected_crank(state, 1.0, beta_de)
            if r < r_star:
                assert want == (1, 1)
            elif r > r_star:
                assert want == (0, 1)
            else:
                assert want[1] >= 2
                second_stage_total += 1
                second_stage_accepts += want[0]
            for name, probe in (
                ("fast", fk.accept_or_reject_probe),
                ("reference", ref_kernel.accept_or_reject_probe),
            ):
                decision, after = probe(state, 1.0, 0, beta_de)
                assert (decision, _draws_between(state, after)) == want, (
                    f"{name} kernel, beta dE = {beta_de}, first draw r = {r} (factor's bin is "
                    f"{r_star}): got decision {decision} after "
                    f"{_draws_between(state, after)} draws, definition gives {want}"
                )
    # the old test accepted all 300 of these; the definition accepts a fraction q
    q = float(scaled - r_star)
    sigma = math.sqrt(second_stage_total * q * (1.0 - q))
    assert abs(second_stage_accepts - second_stage_total * q) <= 5.0 * sigma + 1.0


@pytest.mark.parametrize("beta_de, n_states", [(25, 20000), (30, 200000)])
def test_uphill_acceptance_below_the_old_floor_follows_boltzmann(
    beta_de: int, n_states: int
) -> None:
    """Below 2^-31 the accepted fraction is exp(-beta dE), not 4.66e-10.

    A move this far uphill can only be accepted when the first draw is 0, which
    happens 2^-31 of the time. We force that draw and count acceptances over
    many states that differ in the 33 discarded bits: the old test accepted
    every one of them, the exact test accepts a fraction q = factor * 2^31, so
    that 2^-31 * q is the Boltzmann factor.

    Parameters
    ----------
    beta_de : int
        beta * dE, realised as invtemp = 1 and an integer dE.

    n_states : int
        Number of forced states; chosen so the expected count is tens to hundreds.
    """
    factor = float(_float32_factor(1.0, beta_de))
    q = factor * _TWO31
    assert 0.0 < q < 1.0
    accepted = 0
    for low in range(1, n_states + 1):
        accepted += fk.accept_or_reject_probe(
            _state_for_draw(0, low * 104729), 1.0, 0, beta_de
        )[0]
    sigma = math.sqrt(n_states * q * (1.0 - q))
    assert abs(accepted - n_states * q) <= 5.0 * sigma, (
        f"beta dE = {beta_de}: {accepted} of {n_states} forced zero draws accepted, expected "
        f"{n_states * q:.1f} +- {sigma:.1f}; the old 2^-31 floor would accept all of them"
    )
    # and therefore the unconditional probability is the Boltzmann factor
    implied = (accepted / n_states) * 2.0**-31
    assert implied < 0.5 * 2.0**-31
    assert (
        abs(implied / math.exp(-beta_de) - 1.0) <= 5.0 * sigma / (n_states * q) + 1e-6
    )


def test_serial_acceptance_first_stage_is_untouched() -> None:
    """Away from the factor's bin every decision is the historical one, on one draw.

    This is what keeps existing serial trajectories bit-identical: 40000
    arbitrary states over a spread of temperatures and energy changes, each
    compared with the single-stage float32 test written out here, for the fast
    and the reference kernel.
    """
    rng = np.random.default_rng(20261001)
    states = rng.integers(0, 1 << 63, size=40000, dtype=np.uint64)
    betas = (0.0, 0.01, 0.37, 1.0, 5.0)
    in_threshold_bin = 0
    for i, raw in enumerate(states):
        state = int(raw)
        beta = betas[i % len(betas)]
        d_energy = 1 + (i % 23)
        want = _expected_crank(state, beta, d_energy)
        if want[1] > 1:
            in_threshold_bin += 1
        else:
            r = _mix((state + _GAMMA) & _MASK) >> 33
            assert want[0] == _first_stage(r, _float32_factor(beta, d_energy))
        got_fast = fk.accept_or_reject_probe(state, beta, 0, d_energy)
        got_ref = ref_kernel.accept_or_reject_probe(state, beta, 0, d_energy)
        assert (got_fast[0], _draws_between(state, got_fast[1])) == want
        assert got_ref == got_fast
    # a draw lands in the factor's bin 2^-31 of the time; none of these should
    assert in_threshold_bin == 0


def test_serial_acceptance_edge_cases() -> None:
    """Underflow never accepts, downhill never draws, and the ends of the range hold.

    The first comparison is deliberately left as it was, float32 rounding
    included, so a draw of 2^31 - 1 at beta = 0 is still rejected (it rounds to
    1.0, which is not below a factor of 1.0).
    """
    zero_draw = _state_for_draw(0, 12345)
    max_draw = _state_for_draw(_PRNG_MAX, 12345)
    for probe in (fk.accept_or_reject_probe, ref_kernel.accept_or_reject_probe):
        # exp(-200) underflows float32 to 0: rejected even on a zero draw, one draw used
        decision, after = probe(zero_draw, 1.0, 0, 200)
        assert (decision, _draws_between(zero_draw, after)) == (0, 1)
        # beta = 0: the factor is exactly 1
        decision, after = probe(zero_draw, 0.0, 0, 5)
        assert (decision, _draws_between(zero_draw, after)) == (1, 1)
        decision, after = probe(max_draw, 0.0, 0, 5)
        assert (decision, _draws_between(max_draw, after)) == (0, 1)
        # flat and downhill moves are accepted without a draw
        for new_energy in (0, -3):
            assert probe(zero_draw, 1.0, 0, new_energy) == (1, zero_draw)


@pytest.mark.parametrize(
    "beta_de, n_forward, n_reverse",
    [(20, 3, 7), (25, 3, 7), (30, 5, 2), (60, 1, 9), (700, 8, 1)],
)
def test_pull_acceptance_is_exact_in_the_threshold_bin(
    beta_de: int, n_forward: int, n_reverse: int
) -> None:
    """The pull's Metropolis-Hastings test gets the same second stage.

    Its threshold is ``(nF / nR) exp(-beta dE)`` in double precision and its
    first comparison is ``r / (2^31 - 1) < threshold`` in double; both are kept,
    and the bin that holds the threshold is resolved exactly.

    Parameters
    ----------
    beta_de : int
        beta * dE, realised as invtemp = 1 and an integer dE.

    n_forward, n_reverse : int
        The forward and reverse first-target counts.
    """
    threshold = (float(n_forward) / float(n_reverse)) * math.exp(-float(beta_de))
    assert 0.0 < threshold < 1.0
    scaled = Fraction(threshold) * _TWO31
    r_star = math.floor(scaled)
    refined = 0
    for r in (r_star - 1, r_star, r_star + 1):
        if r < 0:
            continue
        for low in range(1, 301):
            state = _state_for_draw(r, low * 7919 + r)
            if not (float(r) / float(_PRNG_MAX)) < threshold:
                want = (0, 1)
            elif r + 1 <= scaled:
                want = (1, 1)
            else:
                want = _second_stage(state, scaled - r)
                refined += 1
            decision, after = fk.accept_or_reject_ratio_probe(
                state, 1.0, 0, beta_de, n_forward, n_reverse
            )
            assert (decision, _draws_between(state, after)) == want, (
                f"pull, beta dE = {beta_de}, nF/nR = {n_forward}/{n_reverse}, first draw {r}"
            )
    assert refined == 300
    state = _state_for_draw(0, 999)
    # an impossible reverse move and a ratio of 1 or more are settled without a draw
    assert fk.accept_or_reject_ratio_probe(state, 1.0, 0, 5, 3, 0) == (0, state)
    assert fk.accept_or_reject_ratio_probe(state, 1.0, 0, 0, 7, 3) == (1, state)


@pytest.mark.parametrize("dim", [2, 3])
def test_forced_zero_draw_no_longer_accepts_a_forbidden_move(
    tmp_path, dim: int
) -> None:
    """End to end through the production crankshaft kernel.

    Two A monomers sit in contact, and seeds are constructed so that the
    acceptance draw of the first substep is exactly 0. For every such seed that
    proposes an uphill move we run the same substep at beta dE = 60, where the
    Boltzmann factor is 8.8e-27. The old test accepted all of them.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory for the keyfile.

    dim : int
        Dimensionality; a monomer proposal takes `dim` draws before the
        acceptance draw.
    """
    box = [9] * dim
    state = ktu.build_state(
        tmp_path, dim, "SR", False, {"MOVE_CRANKSHAFT": 1.0}, box=box, chains=[(2, "A")]
    )
    idx = np.asarray(state.idx0).copy()
    idx[0, 5 : 5 + dim] = [4] * dim
    idx[1, 5 : 5 + dim] = [5] + [4] * (dim - 1)
    grid = np.zeros(box, dtype=state.lattice.grid.dtype)
    type_grid = np.zeros(box, dtype=state.lattice.type_grid.dtype)
    for row in idx:
        grid[tuple(row[5 : 5 + dim])] = row[4]
        type_grid[tuple(row[5 : 5 + dim])] = row[2]
    e0 = ktu.recompute_energy(state, grid.copy(), type_grid.copy(), idx.copy())
    kernel = fk.mega_crank if dim == 3 else fk.mega_crank_2D
    selector = np.zeros(1, dtype=np.int64)

    tested = 0
    for low in range(1, 3000):
        seed = _state_for_draw(0, low, n=dim + 1)
        g, t, i = grid.copy(), type_grid.copy(), idx.copy()
        e_hot, accepted_hot = kernel(
            g, t, i, *state.tables, e0, 0.0, 1, selector, seed, 0
        )
        d_energy = e_hot - e0
        if accepted_hot != 1 or d_energy <= 0:
            continue
        assert ktu.recompute_energy(state, g, t, i) == e_hot
        g, t, i = grid.copy(), type_grid.copy(), idx.copy()
        e_cold, accepted_cold = kernel(
            g, t, i, *state.tables, e0, 60.0 / d_energy, 1, selector, seed, 0
        )
        assert (e_cold, accepted_cold) == (e0, 0), (
            f"{dim}D crankshaft accepted a move with beta dE = 60 (dE = +{d_energy}) because "
            f"its acceptance draw was 0"
        )
        assert np.array_equal(i, idx) and np.array_equal(g, grid)
        tested += 1
        if tested >= 5:
            break
    assert tested >= 5, "the construction did not produce enough uphill proposals"


@pytest.mark.parametrize("dim", [2, 3])
def test_fast_and_reference_megamoves_agree_through_the_second_stage(
    tmp_path, dim: int
) -> None:
    """Whole megamoves whose first acceptance test enters the second stage.

    The reference kernel takes a C-int seed, so its generator state cannot be
    chosen freely. Instead 31-bit seeds are scanned for ones whose acceptance
    draw r is small, and the temperature is then tuned so that the float32
    Boltzmann factor falls strictly inside that draw's own bin, alternately a
    fifth and four fifths of the way up it so that the second stage both
    accepts and rejects. For each such seed the decision must be the one the
    definition gives, in both kernels, and a 60-substep megamove from the same
    seed must leave the fast and the reference kernel with the same energy,
    accept count, grids and bead table - which they only do if both consumed the
    same number of draws in the second stage.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory for the keyfile.

    dim : int
        Dimensionality; a monomer proposal takes `dim` draws before the
        acceptance draw.
    """
    box = [9] * dim
    state = ktu.build_state(
        tmp_path, dim, "SR", False, {"MOVE_CRANKSHAFT": 1.0}, box=box, chains=[(2, "A")]
    )
    idx = np.asarray(state.idx0).copy()
    idx[0, 5 : 5 + dim] = [4] * dim
    idx[1, 5 : 5 + dim] = [5] + [4] * (dim - 1)
    grid = np.zeros(box, dtype=state.lattice.grid.dtype)
    type_grid = np.zeros(box, dtype=state.lattice.type_grid.dtype)
    for row in idx:
        grid[tuple(row[5 : 5 + dim])] = row[4]
        type_grid[tuple(row[5 : 5 + dim])] = row[2]
    e0 = ktu.recompute_energy(state, grid.copy(), type_grid.copy(), idx.copy())
    fast = fk.mega_crank if dim == 3 else fk.mega_crank_2D
    reference = ref_kernel.mega_crank if dim == 3 else ref_kernel_2D.mega_crank_2D

    # the acceptance draw (draw number dim + 1) of every seed, vectorised
    seeds = np.arange(1, 4000001, dtype=np.uint64)
    z = seeds + np.uint64(((dim + 1) * _GAMMA) & _MASK)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(_C1)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(_C2)
    z = z ^ (z >> np.uint64(31))
    draws = (z >> np.uint64(33)).astype(np.int64)
    candidates = np.nonzero(draws < (1 << 16))[0]

    one = np.zeros(1, dtype=np.int64)
    selector_rng = np.random.RandomState(5)
    outcomes = {0: 0, 1: 0}
    for c in candidates:
        seed, r = int(seeds[c]), int(draws[c])
        g, t, i = grid.copy(), type_grid.copy(), idx.copy()
        e_hot, accepted_hot = fast(g, t, i, *state.tables, e0, 0.0, 1, one, seed, 0)
        d_energy = e_hot - e0
        if accepted_hot != 1 or d_energy <= 0:
            continue
        # a float32 invtemp whose factor lies strictly inside the bin (r, r + 1) / 2^31
        fraction = 0.2 if sum(outcomes.values()) % 2 == 0 else 0.8
        beta = np.float32(-math.log((r + fraction) / _TWO31) / d_energy)
        for _ in range(40):
            beta = np.nextafter(beta, np.float32(0), dtype=np.float32)
        found = None
        for _ in range(81):
            factor = _float32_factor(float(beta), d_energy)
            scaled = Fraction(float(factor)) * _TWO31
            if r < scaled < r + 1 and _first_stage(r, factor):
                found = float(beta)
                break
            beta = np.nextafter(beta, np.float32(np.inf), dtype=np.float32)
        if found is None:
            continue
        before_test = (seed + dim * _GAMMA) & _MASK
        assert _mix((before_test + _GAMMA) & _MASK) >> 33 == r
        want, n_draws = _expected_crank(before_test, found, d_energy)
        assert n_draws >= 2, "the construction must reach the second stage"
        for name, kernel in (("fast", fast), ("reference", reference)):
            g, t, i = grid.copy(), type_grid.copy(), idx.copy()
            _, accepted = kernel(g, t, i, *state.tables, e0, found, 1, one, seed, 0)
            assert accepted == want, (name, seed, r, found)
        outcomes[want] += 1
        # a longer megamove: the draws after the second stage must line up in both kernels
        selector = selector_rng.randint(0, 2, size=60).astype(np.int64)
        selector[0] = 0
        g_f, t_f, i_f = grid.copy(), type_grid.copy(), idx.copy()
        g_r, t_r, i_r = grid.copy(), type_grid.copy(), idx.copy()
        out_f = fast(g_f, t_f, i_f, *state.tables, e0, found, 60, selector, seed, 0)
        out_r = reference(
            g_r, t_r, i_r, *state.tables, e0, found, 60, selector, seed, 0
        )
        assert tuple(out_f) == tuple(out_r), (seed, out_f, out_r)
        assert np.array_equal(g_f, g_r) and np.array_equal(t_f, t_r)
        assert np.array_equal(i_f, i_r)
        assert out_f[0] == ktu.recompute_energy(state, g_f, t_f, i_f)
    assert outcomes[0] >= 3 and outcomes[1] >= 3, (
        f"both outcomes of the second stage must be exercised, got {outcomes}"
    )


# ---------------------------------------------------------------------------
# M3-2: entry guards
# ---------------------------------------------------------------------------


def _padded(array: np.ndarray) -> np.ndarray:
    """A view of `array` that has one valid element on either side of it.

    Used to build out-of-range inputs that stay memory-safe on a kernel with no
    guard: index -1 and index len(array) both land on a copy of element 0.

    Parameters
    ----------
    array : numpy.ndarray
        The array to wrap; padded along its first axis.

    Returns
    -------
    numpy.ndarray
        A contiguous view equal to `array`, inside a larger owned buffer.
    """
    return np.ascontiguousarray(np.concatenate([array[:1], array, array[:1]]))[1:-1]


class _KernelInputs:
    """Valid arguments for every megamove kernel of one dimensionality.

    Attributes
    ----------
    state : ktu.State
        The built system.

    dim : int
        2 or 3.

    grid, type_grid, idx : numpy.ndarray
        Pristine copies of the two grids and the bead table.

    offsets, lengths, homo : numpy.ndarray
        The per-chain layout arrays.

    selector : numpy.ndarray
        A chain selector naming every chain once.

    mask : numpy.ndarray
        An all-zero frozen mask.
    """

    def __init__(self, state: ktu.State) -> None:
        self.state = state
        self.dim = state.dim
        self.grid, self.type_grid, self.idx = state.fresh()
        self.offsets, self.lengths, self.homo = ktu.chain_meta(self.idx)
        self.selector = np.arange(len(self.offsets), dtype=np.int32)
        self.mask = np.zeros(self.idx.shape[0], dtype=np.int32)

    def arrays(self) -> Dict[str, np.ndarray]:
        """Fresh copies of every array argument, keyed by argument name.

        Returns
        -------
        dict
            Argument name to array; the caller swaps in a corrupted entry.
        """
        return {
            "grid": self.grid.copy(),
            "type_grid": self.type_grid.copy(),
            "idx": _padded(self.idx.copy()),
            "bead_selector": np.zeros(4, dtype=np.int64),
            "offsets": self.offsets.copy(),
            "lengths": self.lengths.copy(),
            "homo": self.homo.copy(),
            "selector": self.selector.copy(),
            "mask": self.mask.copy(),
        }

    def call(
        self, name: str, a: Dict[str, np.ndarray], nsteps: int = 4, threads: int = 1
    ) -> tuple:
        """Call one kernel by name with the arrays in `a`.

        Parameters
        ----------
        name : str
            Name of the kernel in ``pimms.mega_crank_fast``.

        a : dict
            The array arguments, as returned by :meth:`arrays`.

        nsteps : int, optional
            Substep count for the crankshaft kernels.

        threads : int, optional
            Thread budget for the parallel kernels.

        Returns
        -------
        tuple
            Whatever the kernel returns.
        """
        st = self.state
        kernel = getattr(fk, name)
        beta, hw, seed = st.acc.invtemp, st.hardwall_int, 12345
        if name in _SERIAL_CRANK.values():
            return kernel(
                a["grid"],
                a["type_grid"],
                a["idx"],
                *st.tables,
                st.energy,
                beta,
                nsteps,
                a["bead_selector"],
                seed,
                hw,
            )
        if name in _PARALLEL_CRANK.values():
            return kernel(
                a["grid"],
                a["type_grid"],
                a["idx"],
                *st.tables,
                st.energy,
                beta,
                nsteps,
                seed,
                hw,
                threads,
                a["mask"],
            )
        chain_args = (
            a["grid"],
            a["type_grid"],
            a["idx"],
            a["offsets"],
            a["lengths"],
            a["homo"],
            a["selector"],
            *st.tables,
            st.energy,
            beta,
            seed,
            hw,
            int(self.lengths.max()),
        )
        if "parallel" in name:
            return kernel(*chain_args, threads, a["mask"])
        return kernel(*chain_args)

    def all_kernels(self) -> List[str]:
        """Names of the six kernels of this dimensionality.

        Returns
        -------
        list of str
            Serial and parallel crankshaft, slither and pull.
        """
        return (
            [_SERIAL_CRANK[self.dim], _PARALLEL_CRANK[self.dim]]
            + _SERIAL_CHAIN[self.dim]
            + _PARALLEL_CHAIN[self.dim]
        )


@pytest.fixture(scope="module", params=[2, 3])
def kernel_inputs(request, tmp_path_factory) -> _KernelInputs:
    """A small mixed system and valid kernel arguments, in 2D and in 3D.

    Parameters
    ----------
    request : pytest.FixtureRequest
        Carries the dimensionality.

    tmp_path_factory : pytest.TempPathFactory
        Where the keyfile is written.

    Returns
    -------
    _KernelInputs
        The bundle the guard tests corrupt one argument of.
    """
    tmp = tmp_path_factory.mktemp(f"guards{request.param}d")
    return _KernelInputs(
        ktu.build_state(tmp, request.param, "SR", False, {"MOVE_CRANKSHAFT": 1.0})
    )


def _assert_refused(
    inputs: _KernelInputs, name: str, bad: Dict[str, np.ndarray], **kwargs: int
) -> None:
    """Assert a kernel raises ValueError naming itself and leaves its inputs alone.

    Parameters
    ----------
    inputs : _KernelInputs
        The valid bundle.

    name : str
        The kernel to call.

    bad : dict
        Corrupted replacements for one or more array arguments.

    **kwargs : int
        Passed to :meth:`_KernelInputs.call` (``nsteps``).
    """
    args = inputs.arrays()
    args.update(bad)
    before = {key: np.array(value, copy=True) for key, value in args.items()}
    with pytest.raises(ValueError, match=name + ":"):
        inputs.call(name, args, **kwargs)
    for key, value in args.items():
        assert np.array_equal(value, before[key]), (
            f"{name} modified {key} before refusing"
        )


def test_guards_let_the_production_inputs_through(kernel_inputs: _KernelInputs) -> None:
    """Every kernel runs on consistent inputs and keeps the energy consistent.

    Parameters
    ----------
    kernel_inputs : _KernelInputs
        The valid bundle.
    """
    for name in kernel_inputs.all_kernels():
        args = kernel_inputs.arrays()
        out = kernel_inputs.call(name, args)
        assert out[0] == ktu.recompute_energy(
            kernel_inputs.state, args["grid"], args["type_grid"], np.array(args["idx"])
        )


def test_guard_bead_table_columns_and_grid_shapes(kernel_inputs: _KernelInputs) -> None:
    """A table with too few columns, or mismatched grids, is refused by every kernel.

    The narrow table is a column slice of the real one, so a kernel without the
    guard would still find its coordinates in the underlying buffer.

    Parameters
    ----------
    kernel_inputs : _KernelInputs
        The valid bundle.
    """
    dim = kernel_inputs.dim
    for name in kernel_inputs.all_kernels():
        narrow = kernel_inputs.arrays()["idx"][:, : 4 + dim]
        _assert_refused(kernel_inputs, name, {"idx": narrow})
        short_type_grid = kernel_inputs.type_grid.copy()[..., :-1]
        _assert_refused(kernel_inputs, name, {"type_grid": short_type_grid})


def test_guard_bead_selector(kernel_inputs: _KernelInputs) -> None:
    """The serial crankshaft refuses a short selector and an out-of-range one.

    Parameters
    ----------
    kernel_inputs : _KernelInputs
        The valid bundle.
    """
    name = _SERIAL_CRANK[kernel_inputs.dim]
    n_beads = kernel_inputs.idx.shape[0]
    # four entries exist in memory, the view exposes one, and two are asked for
    _assert_refused(
        kernel_inputs, name, {"bead_selector": np.zeros(4, np.int64)[:1]}, nsteps=2
    )
    for bad_index in (n_beads, -1):
        selector = np.array([0, bad_index, 0, 0], dtype=np.int64)
        _assert_refused(kernel_inputs, name, {"bead_selector": selector}, nsteps=4)
    # an out-of-range entry beyond nsteps is never read, so it is not an error
    args = kernel_inputs.arrays()
    args["bead_selector"] = np.array([0, 0, n_beads, -1], dtype=np.int64)
    kernel_inputs.call(name, args, nsteps=2)


def test_guard_chain_arrays_and_selector(kernel_inputs: _KernelInputs) -> None:
    """Slither and pull refuse chain arrays that do not describe the table.

    The layout cases use an empty chain selector, so a kernel without the guard
    would attempt nothing; the guard does not depend on the selector.

    Parameters
    ----------
    kernel_inputs : _KernelInputs
        The valid bundle.
    """
    dim = kernel_inputs.dim
    n_beads = kernel_inputs.idx.shape[0]
    n_chains = len(kernel_inputs.offsets)
    empty = np.zeros(0, dtype=np.int32)
    for name in _SERIAL_CHAIN[dim] + _PARALLEL_CHAIN[dim]:
        _assert_refused(
            kernel_inputs,
            name,
            {"homo": kernel_inputs.homo.copy()[:-1], "selector": empty},
        )
        too_long = kernel_inputs.lengths.copy()
        too_long[-1] += 1
        assert kernel_inputs.offsets[-1] + too_long[-1] == n_beads + 1
        _assert_refused(kernel_inputs, name, {"lengths": too_long, "selector": empty})
        negative = kernel_inputs.offsets.copy()
        negative[0] = -1
        _assert_refused(kernel_inputs, name, {"offsets": negative, "selector": empty})
    for name in _SERIAL_CHAIN[dim]:
        for bad_chain in (n_chains, -1):
            bad = {
                "offsets": _padded(kernel_inputs.offsets.copy()),
                "lengths": _padded(kernel_inputs.lengths.copy()),
                "homo": _padded(kernel_inputs.homo.copy()),
                "selector": np.array([0, bad_chain], dtype=np.int32),
            }
            _assert_refused(kernel_inputs, name, bad)


def test_guard_frozen_mask_length(kernel_inputs: _KernelInputs) -> None:
    """The parallel kernels refuse a frozen mask that is not one entry per bead.

    A length-1 mask is the dangerous case: numpy would broadcast it over every
    bead without complaint.

    Parameters
    ----------
    kernel_inputs : _KernelInputs
        The valid bundle.
    """
    dim = kernel_inputs.dim
    for name in [_PARALLEL_CRANK[dim]] + _PARALLEL_CHAIN[dim]:
        _assert_refused(kernel_inputs, name, {"mask": np.zeros(1, dtype=np.int32)})
        _assert_refused(
            kernel_inputs,
            name,
            {"mask": np.zeros(kernel_inputs.idx.shape[0] + 1, dtype=np.int32)},
        )


# ---------------------------------------------------------------------------
# M3-3 and M2-1: thread clamp and the shared shift generator
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def multiblock_inputs(tmp_path_factory) -> Dict[str, _KernelInputs]:
    """Systems on boxes that every parallel kernel splits into several blocks.

    Returns
    -------
    dict
        Kernel name to its bundle: a 26 x 26 box in 2D (9 crankshaft blocks, 4
        whole-chain blocks), 16^3 for the 3D crankshaft (8 blocks) and 24^3 for
        the 3D whole-chain kernels (8 blocks).
    """
    moves = {"MOVE_CRANKSHAFT": 1.0}

    def build(tag: str, dim: int, box: Sequence[int]) -> _KernelInputs:
        """Build one bundle.

        Parameters
        ----------
        tag : str
            Directory name.

        dim : int
            Dimensionality.

        box : sequence of int
            Box dimensions.

        Returns
        -------
        _KernelInputs
            The bundle for that box.
        """
        return _KernelInputs(
            ktu.build_state(
                tmp_path_factory.mktemp(tag), dim, "SR", False, moves, box=list(box)
            )
        )

    two_d = build("mb2d", 2, [26, 26])
    crank_3d = build("mb3dcrank", 3, [16, 16, 16])
    chain_3d = build("mb3dchain", 3, [24, 24, 24])
    assert fk.parallel_crank_layout_info(26, 26, 1, False)["num_blocks"] == 9
    assert fk.parallel_layout_info(26, 26, 1, False)["num_blocks"] == 4
    assert fk.parallel_crank_layout_info(16, 16, 16, False)["num_blocks"] == 8
    assert fk.parallel_layout_info(24, 24, 24, False)["num_blocks"] == 8
    out = {_PARALLEL_CRANK[2]: two_d, _PARALLEL_CRANK[3]: crank_3d}
    out.update({name: two_d for name in _PARALLEL_CHAIN[2]})
    out.update({name: chain_3d for name in _PARALLEL_CHAIN[3]})
    return out


def test_thread_budget_is_clamped_to_the_block_count(
    multiblock_inputs: Dict[str, _KernelInputs],
) -> None:
    """An absurd thread budget is cut to the block count and changes nothing.

    3,000,000,000 does not fit a C int, so before the clamp the kernel raised
    OverflowError at the first parallel megamove of a run whose keyfile had
    passed the parser. The result must be the one a single thread gives.

    Parameters
    ----------
    multiblock_inputs : dict
        Kernel name to a bundle on a multi-block box.
    """
    for name, inputs in multiblock_inputs.items():
        one = inputs.arrays()
        out_one = inputs.call(name, one, nsteps=400, threads=1)
        many = inputs.arrays()
        out_many = inputs.call(name, many, nsteps=400, threads=3000000000)
        assert out_many == out_one, name
        assert out_one[2] > 0, f"{name}: nothing was movable, the comparison is empty"
        for key in ("grid", "type_grid", "idx"):
            assert np.array_equal(one[key], many[key]), (name, key)
    # the rule itself: [1, number of blocks]
    assert fk._clamp_threads(0, 4) == 1
    assert fk._clamp_threads(-3, 4) == 1
    assert fk._clamp_threads(2, 4) == 2
    assert fk._clamp_threads(10000, 4) == 4
    assert fk._clamp_threads(3000000000, 64) == 64


def test_parallel_kernels_do_not_construct_a_generator_per_megamove(
    multiblock_inputs: Dict[str, _KernelInputs], monkeypatch
) -> None:
    """The block shift comes from a reseeded module-level generator.

    Constructing ``numpy.random.RandomState`` costs several hundred
    microseconds, which was most of a small parallel megamove.

    Parameters
    ----------
    multiblock_inputs : dict
        Kernel name to a bundle on a multi-block box.

    monkeypatch : pytest.MonkeyPatch
        Used to count constructions.
    """
    constructed: List[int] = []
    real = np.random.RandomState

    def counting(*args, **kwargs):
        """Stand-in constructor that records each call.

        Parameters
        ----------
        *args, **kwargs
            Forwarded to ``numpy.random.RandomState``.

        Returns
        -------
        numpy.random.RandomState
            A real generator, so a kernel that does construct one still runs.
        """
        constructed.append(1)
        return real(*args, **kwargs)

    prepared = [
        (name, inputs, inputs.arrays()) for name, inputs in multiblock_inputs.items()
    ]
    monkeypatch.setattr(np.random, "RandomState", counting)
    for name, inputs, args in prepared:
        inputs.call(name, args, nsteps=50, threads=1)
    monkeypatch.undo()
    assert len(constructed) == 0, (
        f"{len(constructed)} RandomState constructions in {len(prepared)} parallel megamoves"
    )


# ---------------------------------------------------------------------------
# M2-1 / M2-2: the block shift is a whole-box draw fixed by the seed
# ---------------------------------------------------------------------------


def _axis_block(
    coordinate: int, shift: int, extent: int, n_blocks: int, halo: int
) -> int:
    """Block index of a site along one axis, or -1 if it cannot move this sweep.

    Written from the definition of the decomposition: the axis is cut into
    `n_blocks` blocks of ``extent // n_blocks`` sites starting at `shift`, the
    outermost `halo` sites of every block are frozen, and so is the remainder.

    Parameters
    ----------
    coordinate : int
        The site's coordinate along the axis.

    shift : int
        This sweep's origin shift along the axis.

    extent : int
        Box length along the axis.

    n_blocks : int
        Number of blocks along the axis.

    halo : int
        Frozen-halo width.

    Returns
    -------
    int
        The block index, or -1.
    """
    length = extent // n_blocks
    s = (coordinate - shift) % extent
    if s >= n_blocks * length:
        return -1
    within = s % length
    if within < halo or within >= length - halo:
        return -1
    return s // length


def _single_chain_state(
    tmp_path, dim: int, box: Sequence[int], sequence: str
) -> Tuple[ktu.State, Callable[[Sequence[Sequence[int]]], tuple]]:
    """A one-chain system whose chain can be placed by hand.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Where the keyfile is written.

    dim : int
        Dimensionality.

    box : sequence of int
        Box dimensions.

    sequence : str
        The chain's sequence.

    Returns
    -------
    tuple
        The state, and a function mapping a list of bead positions to the
        ``(grid, type_grid, idx)`` arrays with the chain at those positions.
    """
    state = ktu.build_state(
        tmp_path,
        dim,
        "SR",
        False,
        {"MOVE_CRANKSHAFT": 1.0},
        box=list(box),
        chains=[(1, sequence)],
    )

    def place(positions: Sequence[Sequence[int]]) -> tuple:
        """Build the kernel arrays with the chain at `positions`.

        Parameters
        ----------
        positions : sequence of sequence of int
            One in-box position per bead.

        Returns
        -------
        tuple
            ``(grid, type_grid, idx)``.
        """
        idx = np.asarray(state.idx0).copy()
        grid = np.zeros(list(box), dtype=state.lattice.grid.dtype)
        type_grid = np.zeros(list(box), dtype=state.lattice.type_grid.dtype)
        for row, position in zip(idx, positions):
            row[5 : 5 + dim] = position
            grid[tuple(position)] = row[4]
            type_grid[tuple(position)] = row[2]
        return grid, type_grid, idx

    return state, place


@pytest.mark.parametrize("dim", [2, 3])
def test_crank_block_shift_is_drawn_over_the_whole_box(tmp_path, dim: int) -> None:
    """Which sweeps can move a bead follows from the seed and a whole-box shift.

    A 26-site axis without long-range beads has halo 1 and three blocks of 8,
    leaving a remainder of 2. One monomer is placed at several sites and the
    parallel crankshaft is asked for one attempt per seed; it reports an
    attempt iff the bead is inside a block interior. The oracle draws the shift
    as ``RandomState(seed & 0x7FFFFFFF).randint(0, 26)`` per axis, in x, y, z
    order. With the old one-block-length draw the kernel disagrees with it, and
    the movable fraction of a site depends on its coordinate (0.50 to 0.75 per
    axis); with the whole-box draw it is 18/26 per axis everywhere.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory for the keyfile.

    dim : int
        Dimensionality.
    """
    extent, halo, n_blocks = 26, 1, 3
    assert n_blocks == min(4, extent // (8 * halo)) and extent % n_blocks != 0
    box = [extent] * dim
    state, place = _single_chain_state(tmp_path, dim, box, "A")
    kernel = fk.mega_crank_parallel if dim == 3 else fk.mega_crank_parallel_2D
    mask = np.zeros(1, dtype=np.int32)
    n_seeds = 500
    uniform = ((n_blocks * (extent // n_blocks - 2 * halo)) / extent) ** dim
    sigma = math.sqrt(uniform * (1.0 - uniform) / n_seeds)
    for site in ([0] * dim, [7] * dim, [13] * dim, [25] * dim, [8, 24, 16][:dim]):
        movable_count = 0
        for seed in range(1, n_seeds + 1):
            seed = seed * 2654435761 + 97
            generator = np.random.RandomState(seed & 0x7FFFFFFF)
            shifts = [int(generator.randint(0, extent)) for _ in range(dim)]
            movable = all(
                _axis_block(c, s, extent, n_blocks, halo) >= 0
                for c, s in zip(site, shifts)
            )
            grid, type_grid, idx = place([site])
            _, _, attempted = kernel(
                grid, type_grid, idx, *state.tables, 0, 0.0, 1, seed, 0, 1, mask
            )
            assert attempted == (1 if movable else 0), (
                f"{dim}D parallel crankshaft, bead at {site}, seed {seed}: kernel attempted "
                f"{attempted}, a whole-box shift of {shifts} makes it movable={movable}"
            )
            movable_count += attempted
        assert abs(movable_count / n_seeds - uniform) <= 5.0 * sigma, (
            f"bead at {site}: movable in {movable_count / n_seeds:.3f} of sweeps, "
            f"uniform coverage is {uniform:.3f}"
        )


def test_whole_chain_block_shift_is_unchanged(tmp_path) -> None:
    """The slither kernel's shift is still the same function of the seed.

    The whole-chain kernels already drew over the whole box; this pins that
    moving them onto the shared generator left their shifts alone. A 27-site
    axis without long-range beads has halo 3 and two blocks of 13 (remainder 1);
    a straight 3-mer is movable iff all three beads share one block interior on
    both axes, and the kernel then attempts every entry of the chain selector.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory for the keyfile.
    """
    extent, halo, n_blocks = 27, 3, 2
    assert n_blocks == min(4, extent // (4 * halo))
    state, place = _single_chain_state(tmp_path, 2, [extent, extent], "AAA")
    offsets = np.array([0], dtype=np.int32)
    lengths = np.array([3], dtype=np.int32)
    homo = np.array([1], dtype=np.int32)
    selector = np.zeros(2, dtype=np.int32)
    mask = np.zeros(3, dtype=np.int32)
    seen = set()
    for start in ([0, 0], [11, 20], [25, 9]):
        positions = [[(start[0] + k) % extent, start[1]] for k in range(3)]
        for seed in range(1, 301):
            seed = seed * 2654435761 + 11
            generator = np.random.RandomState(seed & 0x7FFFFFFF)
            shifts = [int(generator.randint(0, extent)) for _ in range(2)]
            blocks = {
                tuple(
                    _axis_block(p[a], shifts[a], extent, n_blocks, halo)
                    for a in range(2)
                )
                for p in positions
            }
            movable = len(blocks) == 1 and min(next(iter(blocks))) >= 0
            grid, type_grid, idx = place(positions)
            _, _, attempted = fk.mega_slither_parallel_2D(
                grid,
                type_grid,
                idx,
                offsets,
                lengths,
                homo,
                selector,
                *state.tables,
                0,
                0.0,
                seed,
                0,
                3,
                1,
                mask,
            )
            assert attempted == (2 if movable else 0), (start, seed, shifts)
            seen.add(movable)
    assert seen == {True, False}


# ---------------------------------------------------------------------------
# M2-2: exact equilibrium of the parallel crankshaft on a box with a remainder
# ---------------------------------------------------------------------------


def _chi2_z(observed: np.ndarray, probabilities: np.ndarray) -> float:
    """Chi-square of counts against exact probabilities, as a standard normal score.

    Bins whose expected count is below 10 are pooled into one. The statistic is
    turned into a z-score with the Wilson-Hilferty cube-root approximation, so
    the caller can threshold it without a p-value table.

    Parameters
    ----------
    observed : numpy.ndarray
        Counts per bin.

    probabilities : numpy.ndarray
        Exact probability of each bin; sums to 1.

    Returns
    -------
    float
        Approximately N(0, 1) when the counts follow the probabilities; large
        and positive when they do not.
    """
    expected = probabilities * observed.sum()
    small = expected < 10.0
    if small.any():
        observed = np.append(observed[~small], observed[small].sum())
        expected = np.append(expected[~small], expected[small].sum())
    dof = len(observed) - 1
    chi2 = float(((observed - expected) ** 2 / expected).sum())
    return ((chi2 / dof) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * dof))) / math.sqrt(
        2.0 / (9.0 * dof)
    )


def _remainder_box_scores(
    tmp_path, seed: int, n_megamoves: int = 400000, stride: int = 10
) -> Tuple[float, float]:
    """Sample one ABA trimer with the parallel crankshaft on a 17 x 17 box.

    Without long-range beads a 17-site axis has halo 1 and two blocks of 8, so
    one site is left over: the remainder case that the whole-box shift changed.
    A single chain in a periodic box has an energy that depends on its
    conformation only, so the exact distribution over energy levels follows from
    the 56 self-avoiding conformations, each weighted by exp(-beta E) with E from
    the from-scratch Hamiltonian.

    The chain's position is deliberately not tested. It is uniform in
    equilibrium, but the chain diffuses across the box far more slowly than its
    conformation relaxes, so position samples a few megamoves apart are strongly
    correlated and a chi-square on them fails on a correct kernel (z up to 10
    over seven seeds at a stride of 10, against at most 1.8 for the energy
    levels).

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory for the keyfile.

    seed : int
        Base seed; megamove m uses ``seed + m``.

    n_megamoves : int, optional
        Parallel crankshaft megamoves of 12 substeps each, on 2 threads.

    stride : int, optional
        A sample is taken every `stride` megamoves, which is what makes the
        samples independent enough for a chi-square test (see the test).

    Returns
    -------
    tuple of (float, float)
        z-scores of the sampled energy-level counts against the exact Boltzmann
        distribution at beta, and against the exact distribution at 1.1 beta
        (the positive control).
    """
    extent = _REMAINDER_BOX
    layout = fk.parallel_crank_layout_info(extent, extent, 1, False)
    assert layout["num_blocks"] > 1 and extent % layout["blocks"][0] != 0, layout
    state, place = _single_chain_state(tmp_path, 2, [extent, extent], "ABA")
    steps = [(dx, dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1) if (dx, dy) != (0, 0)]
    level_weight: Dict[int, List[float]] = {}
    n_conformations = 0
    for d1 in steps:
        for d2 in steps:
            chain = [
                (8, 8),
                (8 + d1[0], 8 + d1[1]),
                (8 + d1[0] + d2[0], 8 + d1[1] + d2[1]),
            ]
            if chain[2] == chain[0]:
                continue
            n_conformations += 1
            energy = ktu.recompute_energy(state, *place([list(p) for p in chain]))
            level_weight.setdefault(energy, [0.0, 0.0])
            level_weight[energy][0] += math.exp(-_REMAINDER_BETA * energy)
            level_weight[energy][1] += math.exp(
                -_CONTROL_SCALE * _REMAINDER_BETA * energy
            )
    assert n_conformations == 56
    levels = sorted(level_weight)
    assert len(levels) >= 3, "the forcefield must separate several energy levels"
    exact = np.array([level_weight[e][0] for e in levels])
    exact /= exact.sum()
    control = np.array([level_weight[e][1] for e in levels])
    control /= control.sum()

    grid, type_grid, idx = place([[8, 8], [9, 8], [10, 8]])
    energy = ktu.recompute_energy(state, grid.copy(), type_grid.copy(), idx.copy())
    mask = np.zeros(3, dtype=np.int32)
    counts = dict.fromkeys(levels, 0)
    attempted_total = 0
    for m in range(n_megamoves):
        energy, _, attempted = fk.mega_crank_parallel_2D(
            grid,
            type_grid,
            idx,
            *state.tables,
            energy,
            _REMAINDER_BETA,
            12,
            seed + m,
            0,
            2,
            mask,
        )
        attempted_total += attempted
        if m % stride == stride - 1:
            counts[energy] += 1
    assert energy == ktu.recompute_energy(state, grid, type_grid, idx)
    assert 0 < attempted_total < 12 * n_megamoves, (
        "some sweeps must find nothing movable"
    )
    observed = np.array([counts[e] for e in levels])
    return _chi2_z(observed, exact), _chi2_z(observed, control)


@pytest.mark.slow
@pytest.mark.parametrize("seed", [1000003, 77000011])
def test_parallel_crankshaft_samples_boltzmann_on_a_remainder_box(
    tmp_path, seed: int
) -> None:
    """The parallel crankshaft reaches the exact distribution on a non-divisible box.

    Every shipped detailed-balance fixture of the parallel crankshaft uses
    40-site axes, which divide evenly into blocks; the whole-box shift changed
    what happens when they do not. Here the energy-level counts of one ABA
    trimer on a 17 x 17 box (two blocks of 8 and a remainder of 1 per axis) are
    compared with the exact Boltzmann distribution by chi-square, and the same
    counts must reject the exact distribution at 1.1 beta. Samples are taken
    every 10 megamoves: over eight seeds the score against the exact
    distribution stayed within 2 and the control score above 10.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory for the keyfile.

    seed : int
        Base seed of the run.
    """
    z_exact, z_control = _remainder_box_scores(tmp_path, seed)
    assert z_exact < 4.0, (
        f"energy levels disagree with the exact distribution (z = {z_exact:.2f})"
    )
    assert z_control > 6.0, (
        f"POSITIVE CONTROL FAILED: the counts are also compatible with beta x {_CONTROL_SCALE} "
        f"(z = {z_control:.2f}), so this test could not see an acceptance error of that size"
    )
