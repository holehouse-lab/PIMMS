"""Proposal-symmetry tests: the half of detailed balance no acceptance rule can repair.

Every move in PIMMS that is accepted with the bare Metropolis criterion
``min(1, exp(-beta dE))`` - the four Python single-chain moves and the crankshaft
kernels - relies on its PROPOSAL being symmetric,

    T_prop(A -> B) == T_prop(B -> A)   for every pair of microstates A, B.

If that fails, nothing downstream can fix it. The strongest violation is
``T(A -> B) > 0`` with ``T(B -> A) == 0``: detailed balance is then broken for
any acceptance function whatsoever, because ``pi(A) T(A->B) a(A->B) > 0`` while
``pi(B) T(B->A) a(B->A) == 0`` identically.

Two families of test live here.

**The four Python single-chain moves** (chain translate, rotate, pivot, head
pivot) had no equilibrium or symmetry test at all. Restricting the pivot's
rotation set to {90, 180} degrees - so the inverse of a 90 degree pivot is never
drawn - passed all 767 tests in the repository while giving 2 of 2 reachable
states in 2D and 4 of 4 in 3D an exactly zero reverse rate. The reversibility
pins that did exist monkeypatch ``random.randint`` to inject the forward angle
and then its inverse, so they test that a CHOSEN rotation is invertible and can
never see a restricted DRAWN set. Everything here therefore uses free draws.

**The crankshaft proposal itself**, which is the universal reference: every
kernel-level detailed-balance test equilibrates with crankshaft and compares
against a crankshaft trace, and the only crank-specific correctness test compares
the fast kernel against the reference kernel it was transcribed from. Breaking
both the same way (dropping the top y-site of the proposal box) left every
bit-exactness and energy-consistency test green and was caught, by accident, in
4 of 20 equilibrium comparisons. The oracle here is written from the definition
of the move in numpy and shares no code with any kernel.

Nothing in this file is marked `slow` (the whole file runs in about two minutes,
and the crank oracle part in about twenty seconds). That is deliberate: the
accidental detection of the crank mutation lived entirely in the one slow-marked
file in the repository, so `-m "not slow"` - which that file's own docstring
recommends for iteration - reported a broken proposal as fully green.
"""
import contextlib
import io
import math
import random
from collections import Counter
from itertools import product
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pytest

import pimms.mega_crank as ref_kernel_3D
import pimms.mega_crank_2D as ref_kernel_2D
import pimms.mega_crank_fast as fk
from pimms import lattice_utils
from pimms.tests import kernel_test_utils as U


Site = Tuple[int, ...]

# ===========================================================================
# Part 1 - the four Python single-chain moves
# ===========================================================================

# Per move: free draws taken from each state, the cap on DISCOVERED destinations,
# and which geometric seed states to add. chain_translate spreads its proposals
# over the whole box, so it needs many draws and few states; the pivots and the
# rotation concentrate on a handful of destinations and need the seeds instead,
# because their characteristic failure is a destination that is never proposed at
# all - which a state set discovered from A alone can never contain.
_MOVE_SETUP: Dict[str, Dict] = {
    # chain_rotate is rigid, so the endpoint seeds are configurations it can
    # neither reach nor leave: they would cost draws and assert nothing. Its own
    # destination set is small enough (3 cardinal rotations in 2D, 23 in 3D) that
    # the discovery pass finds all of it.
    "chain_translate": dict(n_draws=16000, max_discovered=12,
                            seed_kinds=("translate",)),
    "chain_rotate": dict(n_draws=6000, max_discovered=12,
                         seed_kinds=("translate",)),
    "chain_pivot": dict(n_draws=4000, max_discovered=12,
                        seed_kinds=("endpoint",)),
    "head_pivot": dict(n_draws=4000, max_discovered=12,
                       seed_kinds=("endpoint",)),
}

# Discovery pass: enough draws to find every state reachable from A whose
# probability is not vanishingly small.
_N_DISCOVER: int = 4000

# Symmetry threshold in binomial sigma. Under H0 the difference of the two
# directed counts has standard deviation sqrt(nF + nR); with a few hundred
# ordered pairs per case, 4.5 sigma is ~2e-3 expected false alarms per case.
_K_SIGMA: float = 4.5

# A pair with fewer than this many counts in total carries no information about
# the soft symmetry test (the hard one-way test still applies to it).
_MIN_PAIR_COUNTS: int = 40


def _place(state, mobile: Sequence[Sequence[int]],
           obstacles: Sequence[Sequence[int]]) -> None:
    """Reset both grids and put the system in a named configuration.

    Every draw starts from a clean, identical lattice: the moves mutate the grid
    in place when they succeed, so the state has to be rebuilt rather than
    reverted, or a rejected-then-accepted sequence would leak.

    Parameters
    ----------
    state : kernel_test_utils.State
        The built system. Chain 1 is the mobile chain; chains 2, 3, ... are the
        single-bead obstacles.

    mobile : sequence of sequence of int
        Ordered positions for chain 1.

    obstacles : sequence of sequence of int
        One position per obstacle chain, in chain-ID order.
    """
    lat = state.lattice
    lat.grid[:] = 0
    lat.type_grid[:] = 0
    pos = [[int(c) for c in p] for p in mobile]
    lattice_utils.place_chain_by_position(pos, lat.grid, 1, safe=True)
    lat.insert_chain_into_type_grid(1, pos, list(range(len(pos))), safe=True)
    lat.chains[1].set_ordered_positions(pos)
    for k, p in enumerate(obstacles):
        cid = k + 2
        q = [[int(c) for c in p]]
        lattice_utils.place_chain_by_position(q, lat.grid, cid, safe=True)
        lat.insert_chain_into_type_grid(cid, q, [0], safe=True)
        lat.chains[cid].set_ordered_positions(q)


def _key(positions: Iterable[Sequence[int]]) -> Tuple[Site, ...]:
    """Hashable ordered-position key for a chain configuration."""
    return tuple(tuple(int(c) for c in p) for p in positions)


def _draw_outcomes(state, mover: Callable, start: Sequence[Sequence[int]],
                   obstacles: Sequence[Sequence[int]], hardwall: bool,
                   n_draws: int, seed: int) -> Counter:
    """Count the outcome states of `n_draws` FREE proposals from `start`.

    A refused proposal (hard-sphere clash, hardwall violation, degenerate chain)
    is part of the proposal kernel - it maps the draw back to "stay" - so it is
    counted as a self-transition rather than discarded.

    Parameters
    ----------
    state : kernel_test_utils.State
        The built system.

    mover : callable
        One of the bound move methods, called as ``mover(chain, grid, hardwall=)``.

    start : sequence of sequence of int
        The configuration to propose from.

    obstacles : sequence of sequence of int
        Static single-bead chains, re-placed before every draw.

    hardwall : bool
        Passed through to the move.

    n_draws : int
        Number of proposals.

    seed : int
        Seed for the ``random`` module - the moves draw from it directly, and
        nothing is monkeypatched, which is the whole point.

    Returns
    -------
    collections.Counter
        Outcome configuration key -> count.
    """
    random.seed(seed)
    counts: Counter = Counter()
    start_key = _key(start)
    for _ in range(n_draws):
        _place(state, start, obstacles)
        move_event, moved = mover(state.lattice.chains[1], state.lattice.grid,
                                  hardwall=hardwall)
        if not moved:
            counts[start_key] += 1
        else:
            counts[_key(move_event.moved_chain_positions)] += 1
    return counts


def _constructed_states(start: Sequence[Sequence[int]],
                        obstacles: Sequence[Sequence[int]], dims: Sequence[int],
                        hardwall: bool, kinds: Sequence[str]
                        ) -> List[Tuple[Site, ...]]:
    """Legal configurations of the same chain, built geometrically.

    These are seeds for the state set, and they matter for one specific failure
    mode: a move that STOPS proposing a destination entirely. Such a state is
    never discovered by sampling from A - that is precisely what has gone wrong -
    yet A is still reachable from it, so the flow is one-way. Only a state set
    that does not come from the move under test can see that.

    Nothing here calls a move; the configurations are constructed from the
    geometry of the model (a bond is a Chebyshev-1 step, two beads may not share
    a site) and are legal whether or not any move can reach them.

    Parameters
    ----------
    start : sequence of sequence of int
        The reference configuration.

    obstacles : sequence of sequence of int
        Occupied sites belonging to other chains.

    dims : sequence of int
        Box dimensions.

    hardwall : bool
        Under a hard wall, configurations that straddle a face are illegal and
        are not generated.

    kinds : sequence of str
        Any of "endpoint" (relocate each terminal bead to a free neighbour of its
        anchor) and "translate" (rigid shifts of the whole chain).

    Returns
    -------
    list
        Configuration keys, excluding `start` itself.
    """
    dim = len(dims)
    chain = [tuple(int(c) for c in p) for p in start]
    blocked = {tuple(int(c) for c in p) for p in obstacles}
    out = set()

    def legal(positions):
        seen = set()
        for p in positions:
            if p in blocked or p in seen:
                return False
            seen.add(p)
        for a, b in zip(positions, positions[1:]):
            if any(min(abs(a[k] - b[k]), dims[k] - abs(a[k] - b[k])) > 1
                   for k in range(dim)):
                return False
            if hardwall and any(abs(a[k] - b[k]) > 1 for k in range(dim)):
                return False
        return True

    if "endpoint" in kinds:
        for bead, anchor in ((0, 1), (len(chain) - 1, len(chain) - 2)):
            if anchor < 0 or anchor >= len(chain):
                continue
            base = chain[anchor]
            for offset in product((-1, 0, 1), repeat=dim):
                raw = tuple(base[k] + offset[k] for k in range(dim))
                if hardwall and any(not (0 <= raw[k] < dims[k]) for k in range(dim)):
                    continue
                site = tuple(raw[k] % dims[k] for k in range(dim))
                trial = list(chain)
                trial[bead] = site
                if legal(trial):
                    out.add(tuple(trial))

    if "translate" in kinds:
        shifts = [tuple(1 if k == axis else 0 for k in range(dim))
                  for axis in range(dim)]
        shifts += [tuple(-1 for _ in range(dim)), tuple(2 for _ in range(dim))]
        for shift in shifts:
            trial = [tuple((p[k] + shift[k]) % dims[k] for k in range(dim))
                     for p in chain]
            if hardwall and any(not (0 <= p[k] + shift[k] < dims[k])
                                for p in chain for k in range(dim)):
                continue
            if legal(trial):
                out.add(tuple(trial))

    out.discard(tuple(chain))
    return sorted(out)


def _transition_matrix(state, mover: Callable, start: Sequence[Sequence[int]],
                       obstacles: Sequence[Sequence[int]], hardwall: bool,
                       *, n_draws: int, seed: int, max_discovered: int,
                       seed_kinds: Sequence[str]
                       ) -> Tuple[List[Tuple[Site, ...]], np.ndarray]:
    """Directed transition counts among a set of states.

    The set is the union of (a) the states a discovery pass reaches from `start`
    and (b) configurations built geometrically by :func:`_constructed_states`.
    The full pass then draws `n_draws` free proposals from EVERY state in the
    set, including `start`, so the counts in the two directions are two binomials
    with the same number of trials and can be compared without any modelling.

    Parameters
    ----------
    state : kernel_test_utils.State
        The built system.

    mover : callable
        The move under test.

    start : sequence of sequence of int
        The configuration the neighbourhood is grown from.

    obstacles : sequence of sequence of int
        Static obstacles.

    hardwall : bool
        Passed through to the move.

    n_draws : int
        Proposals drawn from each state.

    seed : int
        Base seed; each state gets its own derived stream.

    max_discovered : int
        Cap on the number of DISCOVERED destinations carried into the matrix.
        Only chain_translate needs it (its reachable set is the whole box).

    seed_kinds : sequence of str
        Passed to :func:`_constructed_states`.

    Returns
    -------
    states : list
        The configuration keys, `start` first.

    counts : numpy.ndarray
        Integer matrix, ``counts[i, j]`` = proposals from ``states[i]`` that
        landed on ``states[j]``.
    """
    start_key = _key(start)
    dims = list(state.lattice.dimensions)
    found = _draw_outcomes(state, mover, start, obstacles, hardwall,
                           _N_DISCOVER, seed)
    reachable = sorted(k for k in found if k != start_key)
    if len(reachable) > max_discovered:
        # deterministic, frequency-blind thinning: a broken proposal shows up as
        # a rare or one-way state, so selecting the most frequent would hide it
        step = len(reachable) / max_discovered
        reachable = [reachable[int(k * step)] for k in range(max_discovered)]
    constructed = _constructed_states(start, obstacles, dims, hardwall, seed_kinds)
    extra = sorted(set(constructed) - set(reachable) - {start_key})
    states = [start_key] + reachable + extra
    index = {s: i for i, s in enumerate(states)}
    counts = np.zeros((len(states), len(states)), dtype=np.int64)
    for i, s in enumerate(states):
        outcomes = _draw_outcomes(state, mover, [list(p) for p in s], obstacles,
                                  hardwall, n_draws, seed + 7919 * (i + 1))
        for dest, n in outcomes.items():
            j = index.get(dest)
            if j is not None:
                counts[i, j] = n
    return states, counts


def _assert_proposal_matrix_is_symmetric(states, counts, label: str,
                                         *, min_reachable: int = 2) -> None:
    """Assert a directed proposal-count matrix is symmetric off the diagonal.

    Two assertions, in order of strength. The first is exact and needs no
    statistics: no ordered pair may have a positive forward count and a zero
    reverse count. That single check kills any rotation set, axis set, offset
    range or rejection rule that has stopped being closed under inversion. The
    second is the binomial comparison, which catches a reverse rate that is
    non-zero but wrong.

    Parameters
    ----------
    states : list
        Configuration keys, in the row/column order of `counts`.

    counts : numpy.ndarray
        Directed proposal counts, equal numbers of draws per row.

    label : str
        Case name for the failure messages.

    min_reachable : int, optional
        Fewest destinations the START state must actually reach, so a move that
        silently stops proposing anything cannot pass vacuously. Counted from the
        transitions rather than from the size of the state set, which also holds
        geometrically constructed states the move may never reach.
    """
    reached = int(np.count_nonzero(counts[0])) - int(counts[0, 0] > 0)
    assert reached >= min_reachable, (
        f"{label}: the move reached only {reached} distinct states from the start "
        f"configuration - this case has become vacuous")

    one_way = [(states[i], states[j], int(counts[i, j]))
               for i in range(len(states)) for j in range(len(states))
               if i != j and counts[i, j] > 0 and counts[j, i] == 0]
    assert not one_way, (
        f"{label}: PROPOSAL IS NOT REVERSIBLE - {len(one_way)} ordered pair(s) have "
        f"T(A->B) > 0 and T(B->A) == 0, which breaks detailed balance for ANY "
        f"acceptance rule. First: A={one_way[0][0]} B={one_way[0][1]} "
        f"nF={one_way[0][2]} nR=0")

    worst = (0.0, None)
    for i in range(len(states)):
        for j in range(i + 1, len(states)):
            n_f, n_r = int(counts[i, j]), int(counts[j, i])
            total = n_f + n_r
            if total < _MIN_PAIR_COUNTS:
                continue
            z = abs(n_f - n_r) / math.sqrt(total)
            if z > worst[0]:
                worst = (z, (states[i], states[j], n_f, n_r))
    assert worst[0] <= _K_SIGMA, (
        f"{label}: proposal counts are asymmetric at {worst[0]:.1f} sigma - "
        f"A={worst[1][0]} B={worst[1][1]} nF={worst[1][2]} nR={worst[1][3]}")


# The systems. Boxes are the smallest the keyfile parser allows (7 per axis) or a
# little larger; the obstacles are there to exercise the hard-sphere refusal,
# which is part of the proposal kernel and is where a directional rejection would
# hide. The straddling case starts the chain across a periodic face, the regime
# the 1.0.8 anchor rework touched.
_SYSTEMS_2D: Dict[str, Dict] = {
    "bulk": dict(box=[9, 9], mobile=[[3, 4], [4, 4], [5, 4], [6, 4]],
                 obstacles=[[1, 1], [7, 7], [2, 7], [7, 2], [4, 8]]),
    "straddle": dict(box=[9, 9], mobile=[[7, 4], [8, 4], [0, 4], [1, 4]],
                     obstacles=[[4, 4], [4, 5], [2, 7]]),
}
_SYSTEMS_3D: Dict[str, Dict] = {
    "bulk": dict(box=[7, 7, 7], mobile=[[2, 3, 3], [3, 3, 3], [4, 3, 3], [5, 3, 3]],
                 obstacles=[[1, 1, 1], [5, 5, 5], [1, 5, 1], [5, 1, 5]]),
    "straddle": dict(box=[7, 7, 7], mobile=[[5, 3, 3], [6, 3, 3], [0, 3, 3], [1, 3, 3]],
                     obstacles=[[3, 3, 3], [3, 4, 3], [1, 1, 5]]),
}

_PY_MOVES = ("chain_translate", "chain_rotate", "chain_pivot", "head_pivot")


def _build_move_system(tmp_path, dim: int, hardwall: bool, spec: Dict):
    """Build a Simulation whose lattice we drive by hand.

    Chain 1 is the mobile chain and every other chain is a single-bead obstacle;
    :func:`_place` overwrites all of their positions before each draw, so the
    random initial packing the Simulation produces does not matter.
    """
    n_obstacles = len(spec["obstacles"])
    chains = [(1, "A" * len(spec["mobile"]))]
    if n_obstacles:
        chains.append((n_obstacles, "A"))
    with contextlib.redirect_stdout(io.StringIO()):
        return U.build_state(tmp_path, dim, "SR", hardwall,
                             {"MOVE_CRANKSHAFT": 1.0}, box=spec["box"],
                             chains=chains, temperature=10, seed=7)


@pytest.mark.parametrize("move", _PY_MOVES)
@pytest.mark.parametrize("dim", (2, 3))
# boundary condition and placement are parametrised TOGETHER, over the three
# combinations that exist: a chain cannot straddle a periodic face under hard
# walls, so generating that fourth combination and skipping it at run time would
# only pad the count with cases that were never cases
@pytest.mark.parametrize("hardwall, system", [
    (False, "bulk"), (False, "straddle"), (True, "bulk")],
    ids=["PBC-bulk", "PBC-straddle", "HW-bulk"])
def test_python_single_chain_move_proposal_is_symmetric(tmp_path, move, dim,
                                                        hardwall, system):
    """Free-draw forward/reverse proposal counts must balance, for all four moves.

    This is the test that was missing: nothing in the suite ran chain translate,
    rotate, pivot or head pivot against an equilibrium reference or counted their
    transitions, so a proposal set that is not closed under inversion shipped
    silently. The counts are taken with the module RNG seeded and NOTHING
    monkeypatched, over the whole one-step neighbourhood rather than one
    hand-picked pair, in 2D and 3D, under periodic and hardwall boundaries, and
    with the chain both in the bulk and straddling a periodic face.
    """
    spec = (_SYSTEMS_2D if dim == 2 else _SYSTEMS_3D)[system]
    state = _build_move_system(tmp_path, dim, hardwall, spec)
    mover = getattr(state.sim.MOVER, move)
    states, counts = _transition_matrix(
        state, mover, spec["mobile"], spec["obstacles"], hardwall, seed=1234,
        **_MOVE_SETUP[move])
    _assert_proposal_matrix_is_symmetric(
        states, counts, f"{move} {dim}D {'HW' if hardwall else 'PBC'} {system}")


# ===========================================================================
# Part 2 - the crankshaft proposal, against an independent numpy oracle
# ===========================================================================

_N_TYPES: int = 3


def _zero_tables(dim: int):
    """Interaction and angle tables that are identically zero.

    With every table zeroed, a crank proposal has delta_energy == 0 and is
    accepted without consuming randomness, so the bead's committed position after
    a single substep IS the proposed site (or its old site, if the proposal
    clashed or was refused by the wall). That is what turns the production kernel
    into something whose proposal distribution can be read off directly, and it
    needs no new exports - ``crank_it`` is a cdef with no wrapper, and testing the
    reference kernel alone would prove nothing about the fast one anyway.
    """
    pair = np.zeros((_N_TYPES, _N_TYPES), dtype=np.int32)
    if dim == 3:
        angle = np.zeros((_N_TYPES, 3, 3, 3, 3, 3, 3), dtype=np.int32)
    else:
        angle = np.zeros((_N_TYPES, 3, 3, 3, 3), dtype=np.int32)
    return (pair, pair.copy(), pair.copy(), angle)


def _build_chain(box: Sequence[int], positions: Sequence[Site]):
    """(grid, type_grid, idx) for a single chain occupying `positions`.

    Column 0 of idx is the terminal/internal flag the kernels branch on: 0 for a
    lone monomer, 1 for the N-terminal bead, 3 for the C-terminal bead and 4 for
    an internal (OXO) bead.
    """
    dim = len(box)
    n = len(positions)
    grid = np.zeros(box, dtype=np.int32)
    type_grid = np.zeros(box, dtype=np.int32)
    idx = np.zeros((n, 8), dtype=np.int64)
    if n == 1:
        flags = [0]
    else:
        flags = [1] + [4] * (n - 2) + [3]
    for bead, p in enumerate(positions):
        grid[tuple(p)] = 1
        type_grid[tuple(p)] = 1
        idx[bead, 0] = flags[bead]
        idx[bead, 1] = 0          # not long-range
        idx[bead, 2] = 1          # residue type
        idx[bead, 3] = 1          # chain length flag / bead type
        idx[bead, 4] = 1          # chainID
        idx[bead, 5] = p[0]
        idx[bead, 6] = p[1]
        idx[bead, 7] = p[2] if dim == 3 else 0
    return grid, type_grid, idx


def _candidate_box(positions: Sequence[Site], bead: int,
                   dims: Sequence[int]) -> List[Tuple[int, int]]:
    """Per-axis (min, max) of the sites `bead` may legally occupy, de-periodised.

    A bond is a Chebyshev-1 step, so a bead bonded to neighbours n1 (and n2) may
    sit exactly on the axis-wise intersection of their unit cubes; a lone monomer
    may sit anywhere in its own unit cube. The intersection is taken in a frame
    de-periodised about the bead's CURRENT position - the minimum image of each
    bond vector is the bond vector - which makes the boundary-straddling case fall
    out of the same arithmetic instead of needing its own branch.

    Note that the box is a function of the NEIGHBOURS alone. The bead's own unit
    cube must not be intersected in: that is exactly the dependence on the current
    position that would break inversion closure, and the kernels correctly do not
    have it - a terminal bead routinely lands two sites from where it started.

    Parameters
    ----------
    positions : sequence of tuple of int
        The chain's ordered positions.

    bead : int
        Index of the bead being moved.

    dims : sequence of int
        Box dimensions.

    Returns
    -------
    list of tuple of int
        One (low, high) pair per axis, in the de-periodised frame.
    """
    dim = len(dims)
    own = positions[bead]
    neighbours = [positions[b] for b in (bead - 1, bead + 1)
                  if 0 <= b < len(positions)]
    axes = []
    for k in range(dim):
        length = dims[k]
        lo, hi = None, None
        for neighbour in neighbours:
            delta = (neighbour[k] - own[k] + length // 2) % length - length // 2
            centre = own[k] + delta
            lo = centre - 1 if lo is None else max(lo, centre - 1)
            hi = centre + 1 if hi is None else min(hi, centre + 1)
        if lo is None:
            lo, hi = own[k] - 1, own[k] + 1
        axes.append((lo, hi))
    return axes


def _expected_distribution(box: Sequence[int], positions: Sequence[Site],
                           bead: int, hardwall: bool) -> Dict[Site, float]:
    """Expected distribution of the committed site, from the move's definition.

    Uniform over the candidate box. Sites that are occupied fold onto "no move",
    and under a hard wall so do sites that would only be reachable by wrapping -
    both foldings are symmetric because the box is fixed by the neighbours, not
    by where the bead currently is.

    Returns
    -------
    dict
        Site -> probability. Includes the bead's own site, which collects the
        "propose to stay" and every folded proposal.
    """
    dims = list(box)
    axes = _candidate_box(positions, bead, dims)
    raw = list(product(*[range(lo, hi + 1) for lo, hi in axes]))
    occupied = {tuple(p) for i, p in enumerate(positions) if i != bead}
    own = tuple(positions[bead])
    probs: Dict[Site, float] = {}
    for site in raw:
        wrapped = tuple(v % L for v, L in zip(site, dims))
        outside = any(not (0 <= v < L) for v, L in zip(site, dims))
        if wrapped in occupied or (hardwall and outside):
            destination = own
        else:
            destination = wrapped
        probs[destination] = probs.get(destination, 0.0) + 1.0 / len(raw)
    return probs


def _sample_serial(kernel, box, positions, bead, hardwall, n_trials, seed0):
    """Committed-site histogram over `n_trials` one-substep crank megamoves."""
    dim = len(box)
    tables = _zero_tables(dim)
    selector = np.array([bead], dtype=np.int64)
    counts: Counter = Counter()
    for trial in range(n_trials):
        grid, type_grid, idx = _build_chain(box, positions)
        kernel(grid, type_grid, idx, *tables, 0, np.float32(0.0), 1, selector,
               seed0 + trial, 1 if hardwall else 0)
        counts[tuple(int(idx[bead, 5 + k]) for k in range(dim))] += 1
    return counts


def _sample_parallel(kernel, box, positions, bead, hardwall, n_trials, seed0,
                     n_threads=2):
    """The same histogram for the parallel checkerboard kernels.

    The parallel kernel has no bead selector - it picks beads itself and freezes
    whatever falls in the current halo - so most substeps leave the bead alone.
    Those extra self-transitions are an artefact of the driver, not of the
    proposal, so callers compare the distribution CONDITIONAL on the bead moving.
    """
    dim = len(box)
    tables = _zero_tables(dim)
    counts: Counter = Counter()
    for trial in range(n_trials):
        grid, type_grid, idx = _build_chain(box, positions)
        frozen = np.zeros(idx.shape[0], dtype=np.int32)
        kernel(grid, type_grid, idx, *tables, 0, np.float32(0.0), 1,
               seed0 + trial, 1 if hardwall else 0, n_threads, frozen)
        counts[tuple(int(idx[bead, 5 + k]) for k in range(dim))] += 1
    return counts


def _chi_square(counts: Counter, probs: Dict[Site, float], label: str,
                *, conditional: bool, own: Site) -> Tuple[float, int]:
    """Compare an empirical site histogram with the oracle's distribution.

    Parameters
    ----------
    counts : collections.Counter
        Observed committed sites.

    probs : dict
        Oracle distribution.

    label : str
        Case name for the failure message.

    conditional : bool
        Drop the "stay" bucket and renormalise. Needed for the parallel kernels,
        whose driver adds self-transitions of its own.

    own : tuple of int
        The bead's starting site, i.e. the "stay" bucket.

    Returns
    -------
    chi2 : float
    dof : int
    """
    observed = dict(counts)
    expected = dict(probs)
    if conditional:
        observed.pop(own, None)
        expected.pop(own, None)
        scale = sum(expected.values())
        assert scale > 0, f"{label}: the oracle predicts no move is ever possible"
        expected = {k: v / scale for k, v in expected.items()}
    off_box = {k: v for k, v in observed.items() if k not in expected}
    assert not off_box, (
        f"{label}: the kernel proposed {len(off_box)} site(s) OUTSIDE the "
        f"anchor-intersection box: {sorted(off_box)[:5]}")
    total = sum(observed.values())
    missing = [k for k, v in expected.items() if v > 0 and observed.get(k, 0) == 0]
    assert not missing, (
        f"{label}: {len(missing)} legal site(s) were NEVER proposed in {total} "
        f"draws: {sorted(missing)[:5]} - the proposal box is not closed under "
        f"inversion, which no acceptance rule can repair")
    chi2 = 0.0
    for site, p in expected.items():
        e = p * total
        chi2 += (observed.get(site, 0) - e) ** 2 / e
    return chi2, max(len(expected) - 1, 1)


# Geometries. Straight along each axis, bent in two planes, a two-away diagonal
# whose box degenerates to a line of three sites, a triptic straddling a periodic
# face, and (hardwall) one against the wall and one along it - the places where
# the de-periodisation or the wall refusal could break the symmetry independently
# of the draw.
_CRANK_3D: List[Tuple[str, List[Site], bool]] = [
    ("straight-x", [(3, 5, 5), (4, 5, 5), (5, 5, 5)], False),
    ("straight-y", [(5, 3, 5), (5, 4, 5), (5, 5, 5)], False),
    ("straight-z", [(5, 5, 3), (5, 5, 4), (5, 5, 5)], False),
    ("bent-xy", [(3, 5, 5), (4, 5, 5), (4, 6, 5)], False),
    ("bent-yz", [(5, 3, 5), (5, 4, 5), (5, 4, 6)], False),
    ("diagonal", [(3, 4, 5), (4, 5, 5), (5, 6, 5)], False),
    # the moved bead sits off the line of its anchors, so the anchor box extends
    # two sites away from it on y: a proposal box that (wrongly) also intersected
    # the bead's own unit cube would never offer y = 5 + 1 here
    ("offset-bead", [(3, 5, 5), (4, 4, 5), (5, 5, 5)], False),
    ("straddle-x", [(10, 5, 5), (0, 5, 5), (1, 5, 5)], False),
    ("straddle-z", [(5, 5, 10), (5, 5, 0), (5, 5, 1)], False),
    ("wall-along", [(0, 4, 5), (0, 5, 5), (0, 6, 5)], True),
    ("wall-adjacent", [(0, 4, 5), (1, 4, 5), (2, 4, 5)], True),
]
_CRANK_2D: List[Tuple[str, List[Site], bool]] = [
    ("straight-x", [(3, 5), (4, 5), (5, 5)], False),
    ("straight-y", [(5, 3), (5, 4), (5, 5)], False),
    ("bent-xy", [(3, 5), (4, 5), (4, 6)], False),
    ("diagonal", [(3, 4), (4, 5), (5, 6)], False),
    ("straddle-y", [(5, 10), (5, 0), (5, 1)], False),
    ("wall-along", [(0, 4), (0, 5), (0, 6)], True),
    ("wall-adjacent", [(0, 4), (1, 4), (2, 4)], True),
]

_BOX_3D = (11, 11, 11)
_BOX_2D = (11, 11)

# 20000 draws over a 9- or 12-site box gives >1500 counts per site, which puts a
# single dropped site at chi2 in the thousands (measured 15003 on the mutation
# that the rest of the suite nearly missed) against a threshold of 60.
_N_CRANK_TRIALS: int = 20000
_CHI2_LIMIT: float = 60.0


@pytest.mark.parametrize("kernel_name", ("fast", "reference"))
@pytest.mark.parametrize("name,positions,hardwall", _CRANK_3D,
                         ids=[c[0] for c in _CRANK_3D])
def test_crank_3D_internal_proposal_is_uniform_over_the_anchor_box(
        kernel_name, name, positions, hardwall):
    """The 3D crank proposal must be uniform over exactly the bonded in-box sites.

    The oracle is the axis-wise intersection of the two anchors' unit cubes,
    computed in numpy from the definition of a bond; it shares no code with any
    kernel, which is what the suite lacked. ``test_fast_crank_bit_exact_vs_reference``
    pins the fast kernel against the reference kernel it was transcribed from, so
    a proposal bug applied to both is invisible to it - both kernels are checked
    here, separately, against something neither of them produced.
    """
    kernel = fk.mega_crank if kernel_name == "fast" else ref_kernel_3D.mega_crank
    label = f"crank 3D {kernel_name} {name}"
    counts = _sample_serial(kernel, _BOX_3D, positions, 1, hardwall,
                            _N_CRANK_TRIALS, seed0=1)
    probs = _expected_distribution(_BOX_3D, positions, 1, hardwall)
    chi2, dof = _chi_square(counts, probs, label, conditional=False,
                            own=tuple(positions[1]))
    assert chi2 < _CHI2_LIMIT, (
        f"{label}: proposal distribution is not uniform over the "
        f"anchor-intersection box - chi2 = {chi2:.1f} on {dof} dof")


@pytest.mark.parametrize("kernel_name", ("fast", "reference"))
@pytest.mark.parametrize("name,positions,hardwall", _CRANK_2D,
                         ids=[c[0] for c in _CRANK_2D])
def test_crank_2D_internal_proposal_is_uniform_over_the_anchor_box(
        kernel_name, name, positions, hardwall):
    """The 2D twin of the 3D uniformity oracle (mega_crank_2D)."""
    kernel = fk.mega_crank_2D if kernel_name == "fast" else ref_kernel_2D.mega_crank_2D
    label = f"crank 2D {kernel_name} {name}"
    counts = _sample_serial(kernel, _BOX_2D, positions, 1, hardwall,
                            _N_CRANK_TRIALS, seed0=1)
    probs = _expected_distribution(_BOX_2D, positions, 1, hardwall)
    chi2, dof = _chi_square(counts, probs, label, conditional=False,
                            own=tuple(positions[1]))
    assert chi2 < _CHI2_LIMIT, (
        f"{label}: proposal distribution is not uniform over the "
        f"anchor-intersection box - chi2 = {chi2:.1f} on {dof} dof")


@pytest.mark.parametrize("name,positions,hardwall",
                         [c for c in _CRANK_3D if c[0] in
                          ("straight-x", "bent-xy", "straddle-x", "wall-along")],
                         ids=["straight-x", "bent-xy", "straddle-x", "wall-along"])
def test_crank_3D_terminal_and_monomer_proposals_are_uniform(name, positions, hardwall):
    """Terminal beads and lone monomers get the same treatment as internal beads.

    A terminal bead has one anchor, so its candidate box is that anchor's whole
    unit cube; a monomer has none, so its box is its own. Both go through separate
    branches of the kernel (idx column 0 is 1/3 for termini, 0 for a monomer) and
    neither had any independent check.
    """
    label = f"crank 3D terminal {name}"
    counts = _sample_serial(fk.mega_crank, _BOX_3D, positions, 0, hardwall,
                            _N_CRANK_TRIALS, seed0=31)
    probs = _expected_distribution(_BOX_3D, positions, 0, hardwall)
    chi2, dof = _chi_square(counts, probs, label, conditional=False,
                            own=tuple(positions[0]))
    assert chi2 < _CHI2_LIMIT, f"{label}: chi2 = {chi2:.1f} on {dof} dof"


@pytest.mark.parametrize("dim", (2, 3))
def test_crank_monomer_proposal_is_uniform_over_its_own_unit_cube(dim):
    """A lone monomer's crank proposal covers its whole unit cube, uniformly."""
    box = _BOX_3D if dim == 3 else _BOX_2D
    kernel = fk.mega_crank if dim == 3 else fk.mega_crank_2D
    positions = [(5, 5, 5)] if dim == 3 else [(5, 5)]
    label = f"crank {dim}D monomer"
    counts = _sample_serial(kernel, box, positions, 0, False, _N_CRANK_TRIALS,
                            seed0=57)
    probs = _expected_distribution(box, positions, 0, False)
    assert len(probs) == 3 ** dim, (
        f"{label}: the oracle box has {len(probs)} sites, expected {3 ** dim}")
    chi2, dof = _chi_square(counts, probs, label, conditional=False,
                            own=tuple(positions[0]))
    assert chi2 < _CHI2_LIMIT, f"{label}: chi2 = {chi2:.1f} on {dof} dof"


@pytest.mark.parametrize("kernel_name", ("fast", "reference"))
@pytest.mark.parametrize("name,positions,hardwall",
                         [c for c in _CRANK_3D if c[0] in
                          ("straight-x", "bent-xy", "straddle-x", "wall-along")],
                         ids=["straight-x", "bent-xy", "straddle-x", "wall-along"])
def test_crank_3D_proposal_is_closed_under_inversion(kernel_name, name, positions,
                                                     hardwall):
    """T(s -> s') must equal T(s' -> s) for every pair of free sites in the box.

    The candidate box depends only on the anchors, never on where the bead
    currently sits, so every free site in it must see every other with the same
    rate. This is the property the reviewers' mutation destroyed - a bead on the
    dropped face could leave and nothing could enter - and it is checked here by
    placing the bead at each free site in turn and counting, rather than by
    reasoning about the draw.
    """
    kernel = fk.mega_crank if kernel_name == "fast" else ref_kernel_3D.mega_crank
    label = f"crank 3D {kernel_name} {name} inversion"
    dims = list(_BOX_3D)
    axes = _candidate_box(positions, 1, dims)
    anchors = {tuple(positions[0]), tuple(positions[2])}
    free = sorted({tuple(v % L for v, L in zip(site, dims))
                   for site in product(*[range(lo, hi + 1) for lo, hi in axes])
                   if not (hardwall and any(not (0 <= v < L)
                                            for v, L in zip(site, dims)))}
                  - anchors)
    assert len(free) >= 3, f"{label}: only {len(free)} free sites - vacuous"

    flow: Dict[Tuple[Site, Site], int] = {}
    for start in free:
        placed = [positions[0], start, positions[2]]
        counts = _sample_serial(kernel, _BOX_3D, placed, 1, hardwall, 8000,
                                seed0=100000)
        for site, n in counts.items():
            if site != start:
                flow[(start, site)] = n

    one_way = [(a, b) for (a, b), n in flow.items() if n > 0 and flow.get((b, a), 0) == 0]
    assert not one_way, (
        f"{label}: one-way proposal flow between free sites of the same box: "
        f"{one_way[:4]}")
    worst = (0.0, None)
    for (a, b), n_f in flow.items():
        n_r = flow.get((b, a), 0)
        total = n_f + n_r
        if total < _MIN_PAIR_COUNTS:
            continue
        z = abs(n_f - n_r) / math.sqrt(total)
        if z > worst[0]:
            worst = (z, (a, b, n_f, n_r))
    assert worst[0] <= 5.0, (
        f"{label}: proposal flow is asymmetric at {worst[0]:.1f} sigma - {worst[1]}")


# The parallel checkerboard kernels carry their own copies of the crank proposal
# (crank_it_cp / crank_it_cp_2D), which nothing independent has ever checked. They
# pick their own beads and freeze whatever lands in the current halo, so the
# comparison is conditional on the bead having moved; the box is the smallest that
# still decomposes into more than one block, with the chain inside a block
# interior. Fewer trials because each call spins up a thread team.
_N_PARALLEL_TRIALS: int = 20000


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("hardwall", (False, True), ids=["PBC", "HW"])
def test_parallel_crank_proposal_is_uniform_over_the_anchor_box(dim, hardwall):
    """The parallel kernels' crank proposal obeys the same oracle as the serial one."""
    if dim == 3:
        box, kernel = (32, 32, 32), fk.mega_crank_parallel
        positions = [(8, 8, 8), (9, 8, 8), (10, 8, 8)]
    else:
        box, kernel = (40, 40), fk.mega_crank_parallel_2D
        positions = [(8, 8), (9, 8), (10, 8)]
    assert fk.parallel_crank_layout_info(
        box[0], box[1], box[2] if dim == 3 else 1, True)["num_blocks"] > 1, (
        "the parallel layout must actually split, or this tests the serial path")
    label = f"parallel crank {dim}D {'HW' if hardwall else 'PBC'}"
    counts = _sample_parallel(kernel, box, positions, 1, hardwall,
                              _N_PARALLEL_TRIALS, seed0=11)
    probs = _expected_distribution(box, positions, 1, hardwall)
    chi2, dof = _chi_square(counts, probs, label, conditional=True,
                            own=tuple(positions[1]))
    assert chi2 < _CHI2_LIMIT, (
        f"{label}: proposal distribution is not uniform over the "
        f"anchor-intersection box - chi2 = {chi2:.1f} on {dof} dof")
