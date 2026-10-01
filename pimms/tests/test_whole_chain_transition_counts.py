"""
Forward/reverse transition counting for single slither and pull sub-moves.

The equilibrium comparisons in ``test_detailed_balance.py`` cannot resolve a
20 % error in the pull's Metropolis-Hastings ratio (most accepted pulls are
nearly energy-neutral, so the bias drowns in the relative floor), and the
bit-exact kernel tests cannot see a bug the fast and reference kernels share.
This file pins the ratio directly. From a fixed state x, one sub-move of one
chain is run many times with fresh seeds; for the states y it reaches, the same
sub-move is run many times from y and the returns to x are counted. Detailed
balance requires

    P(x -> y) / P(y -> x) = exp(-beta (E_y - E_x))

for every pair, with E recomputed from scratch by the Hamiltonian rather than
taken from the kernel. A chi-square over the pairs must accept the true beta and
reject beta x 1.2, so the test is known to resolve an error of that size.
"""

from __future__ import annotations

import collections

import numpy as np
import pytest
from scipy import stats

from pimms import mega_crank_fast as fk
from pimms.tests import kernel_test_utils as U

N_FORWARD: int = 60000
N_REVERSE: int = 40000
MIN_COUNT: int = 300
N_PAIRS: int = 6

Key = tuple[tuple[int, int], ...]


def _counts(move: str, hardwall: bool, tmp_path) -> tuple[list[tuple[int, int, int]], float]:
    """Forward and reverse counts for the highest-|dE| pairs reached from the start state.

    Parameters
    ----------
    move : str
        ``'slither'`` or ``'pull'``.
    hardwall : bool
        Whether the 7x7 box has hard walls.
    tmp_path : pathlib.Path
        Directory to build the simulation in.

    Returns
    -------
    tuple
        ``(pairs, invtemp)`` where each pair is ``(n_forward, n_reverse, dE)``.
    """
    state = U.build_state(tmp_path, 2, "SR", hardwall, {"MOVE_PULL": 1.0}, box=[7, 7],
                          chains=[(1, "AABBA")], temperature=30)
    grid0, type_grid0, idx0 = state.fresh()
    offsets, lengths, homo = U.chain_meta(idx0)
    kernel = fk.mega_pull_2D if move == "pull" else fk.mega_slither_2D
    selector = np.zeros(1, dtype=np.int32)
    invtemp = float(state.acc.invtemp)

    def build(key: Key) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        grid = np.zeros_like(grid0)
        type_grid = np.zeros_like(type_grid0)
        idx = idx0.copy()
        for row, site in enumerate(key):
            idx[row, 5:7] = site
            grid[site] = idx[row, 4]
            type_grid[site] = idx[row, 2]
        return grid, type_grid, idx

    def step(key: Key, seed: int) -> tuple[Key, int]:
        grid, type_grid, idx = build(key)
        energy, _accepted = kernel(grid, type_grid, idx, offsets, lengths, homo, selector,
                                   *state.tables, 0, invtemp, seed, state.hardwall_int,
                                   int(lengths.max()))
        return tuple(map(tuple, idx[:, 5:7].tolist())), int(energy)

    x = tuple(map(tuple, idx0[:, 5:7].tolist()))
    rng = np.random.default_rng(1)
    forward: collections.Counter = collections.Counter()
    tracked_dE: dict[Key, int] = {}
    for seed in rng.integers(1, 2**62, N_FORWARD):
        y, dE = step(x, int(seed))
        if y != x:
            forward[y] += 1
            tracked_dE[y] = dE

    # the pairs that carry the most information about beta are the ones with the
    # largest energy change
    candidates = [y for y, n in forward.items() if n >= MIN_COUNT]
    candidates.sort(key=lambda y: (-abs(tracked_dE[y]), -forward[y]))
    E_x = U.recompute_energy(state, *build(x))
    pairs = []
    for y in candidates[:N_PAIRS]:
        back = sum(1 for seed in rng.integers(1, 2**62, N_REVERSE) if step(y, int(seed))[0] == x)
        dE = U.recompute_energy(state, *build(y)) - E_x
        # the kernel's own delta must agree with the from-scratch energies
        assert tracked_dE[y] == dE
        pairs.append((forward[y], back, dE))
    return pairs, invtemp


def _chi_square(pairs: list[tuple[int, int, int]], beta: float) -> float:
    """Chi-square of the log forward/reverse ratios against exp(-beta dE).

    Parameters
    ----------
    pairs : list of tuple
        ``(n_forward, n_reverse, dE)`` per pair.
    beta : float
        The inverse temperature the ratios are tested against.

    Returns
    -------
    float
        Sum over pairs of the squared z-score of the log ratio.
    """
    chi2 = 0.0
    for n_fwd, n_rev, dE in pairs:
        log_ratio = np.log(n_fwd / N_FORWARD) - np.log(n_rev / N_REVERSE)
        chi2 += (log_ratio + beta * dE) ** 2 / (1.0 / n_fwd + 1.0 / n_rev)
    return chi2


@pytest.mark.parametrize("hardwall", (False, True), ids=["PBC", "HW"])
@pytest.mark.parametrize("move", ("slither", "pull"))
def test_single_sub_move_transition_counts_obey_detailed_balance(move: str, hardwall: bool,
                                                                 tmp_path) -> None:
    pairs, beta = _counts(move, hardwall, tmp_path)
    assert len(pairs) >= 2, "too few reachable states to test anything"
    assert all(n_rev > 0 for _n, n_rev, _dE in pairs)
    assert any(dE != 0 for _n, _r, dE in pairs), "every pair was energy-neutral"
    dof = len(pairs)

    chi2 = _chi_square(pairs, beta)
    assert chi2 < stats.chi2.ppf(1.0 - 1e-4, dof), (
        f"{move}: forward/reverse ratios disagree with exp(-beta dE): chi2 = {chi2:.1f} "
        f"on {dof} pairs {pairs}")

    # positive control: the same counts must reject a 20 % error in the exponent
    chi2_wrong = _chi_square(pairs, 1.2 * beta)
    assert chi2_wrong > stats.chi2.ppf(1.0 - 1e-3, dof), (
        f"{move}: the test cannot resolve a 20 % Metropolis error (chi2 = {chi2_wrong:.1f})")
