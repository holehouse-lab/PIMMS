"""
Regression tests for the Chain-level internal-scaling accumulators.

``INTSCAL.dat`` holds <r_ij> and ``INTSCAL_SQUARED.dat`` holds <r_ij^2>, each
averaged over every pair at a sequence gap and over every snapshot. Before 1.0.8
the squared profile was fed the per-snapshot MEAN distance and squared it, so it
held the average of <r>^2 and silently dropped the within-snapshot spread of pair
distances. The pre-existing tests could not see that: they used straight rods
(every pair at a gap has the same distance, so <r^2> == <r>^2) or pinned only the
last row (gap L-1, one pair per chain). These tests drive the real
``Chain.analysis_update_internal_scaling`` over several random self-avoiding
conformations and compare BOTH accumulators, at every gap, against a plain-numpy
mean over pairs and snapshots written from the definition.
"""

import numpy as np
import pytest

from pimms import chain as chain_module


def _self_avoiding_walk(rng: np.random.Generator, n_beads: int, n_dim: int) -> np.ndarray:
    """
    Grow one random self-avoiding walk on the simple cubic (or square) lattice.

    Each step picks uniformly among the unvisited nearest neighbours of the last
    bead; a walk that traps itself is thrown away and regrown.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.

    n_beads : int
        Number of beads in the walk.

    n_dim : int
        Lattice dimensionality, 2 or 3.

    Returns
    -------
    numpy.ndarray
        ``(n_beads, n_dim)`` integer array of unwrapped bead coordinates, first
        bead at the origin.
    """
    steps = np.concatenate([np.eye(n_dim, dtype=int), -np.eye(n_dim, dtype=int)])
    while True:
        walk = [np.zeros(n_dim, dtype=int)]
        visited = {tuple(walk[0])}
        while len(walk) < n_beads:
            options = [walk[-1] + s for s in steps if tuple(walk[-1] + s) not in visited]
            if not options:
                break
            nxt = options[rng.integers(len(options))]
            walk.append(nxt)
            visited.add(tuple(nxt))
        if len(walk) == n_beads:
            return np.array(walk)


def _pair_distances_by_gap(walk: np.ndarray) -> dict[int, np.ndarray]:
    """
    Every pair distance of one conformation, grouped by sequence gap.

    Parameters
    ----------
    walk : numpy.ndarray
        ``(n_beads, n_dim)`` unwrapped (whole-chain) bead coordinates.

    Returns
    -------
    dict of int to numpy.ndarray
        Maps each gap ``g`` in ``1 .. n_beads - 1`` to the ``n_beads - g``
        Euclidean distances ``|x_{i+g} - x_i|``.
    """
    n_beads = len(walk)
    out: dict[int, np.ndarray] = {}
    for gap in range(1, n_beads):
        delta = (walk[gap:] - walk[:-gap]).astype(float)
        out[gap] = np.sqrt((delta ** 2).sum(axis=1))
    return out


@pytest.mark.parametrize("n_dim, hardwall", [(3, False), (2, False), (3, True)])
def test_chain_accumulates_mean_r_and_mean_r_squared_at_every_gap(n_dim: int, hardwall: bool) -> None:
    rng = np.random.default_rng(20260926 + n_dim + int(hardwall))
    n_beads = 12
    n_snapshots = 7
    box = 17
    dims = [box] * n_dim

    # conformations: random self-avoiding walks. Under periodic boundaries each
    # one is dropped at a random offset and wrapped into the box, so some cross a
    # face and the Chain has to make them whole again; under a hardwall they are
    # placed inside the box untouched
    walks = [_self_avoiding_walk(rng, n_beads, n_dim) for _ in range(n_snapshots)]
    stored: list[list[list[int]]] = []
    for walk in walks:
        if hardwall:
            placed = walk - walk.min(axis=0) + 1
            assert placed.max() < box
        else:
            placed = (walk + rng.integers(0, box, size=n_dim)) % box
        stored.append([[int(v) for v in bead] for bead in placed])
    if not hardwall:
        # at least one stored conformation must actually be torn by the box face,
        # or the make-whole path is not exercised
        assert any(np.abs(np.diff(np.array(s), axis=0)).max() > 1 for s in stored)

    chain = chain_module.Chain(
        lattice_grid=np.zeros(dims, dtype=np.int32), dimensions=dims,
        sequence="A" * n_beads, int_seq=[1] * n_beads, LR_int_seq=[1] * n_beads,
        LR_IDX=[], chainID=1, chainType=0, chain_positions=stored[0],
        hardwall=hardwall)
    for positions in stored:
        chain.set_ordered_positions(positions)
        chain.analysis_update_internal_scaling()

    # oracle: pool every pair distance at a gap across all snapshots, then take
    # the plain mean of r and of r**2 (each snapshot contributes the same number
    # of pairs per gap, so this is also the mean of the per-snapshot means)
    per_snapshot = [_pair_distances_by_gap(walk) for walk in walks]
    gaps = range(1, n_beads)
    pooled = {g: np.concatenate([snap[g] for snap in per_snapshot]) for g in gaps}
    expected_r = np.array([pooled[g].mean() for g in gaps])
    expected_r2 = np.array([(pooled[g] ** 2).mean() for g in gaps])

    np.testing.assert_allclose(chain.analysis_get_cumulative_internal_scaling(),
                               expected_r, rtol=1e-12, atol=0)
    np.testing.assert_allclose(chain.analysis_get_internal_scaling_squared(),
                               expected_r2, rtol=1e-12, atol=0)

    # non-vacuity: the conformations must make <r^2> differ from the average of the
    # squared per-snapshot means at every interior gap (gap 1 is a bond, always
    # length 1, and gap L-1 has a single pair), or the pre-1.0.8 bug could pass
    squared_means = np.array([np.mean([snap[g].mean() ** 2 for snap in per_snapshot])
                              for g in gaps])
    interior = slice(1, n_beads - 2)
    assert np.all(expected_r2[interior] - squared_means[interior] > 1e-3)


def test_squared_accumulator_is_not_the_square_of_the_mean_on_an_l_shape() -> None:
    """A hand-checkable case: one L-shaped 4-bead chain, sampled once.

    Three beads run along x and the fourth turns up in y. Gap 1 is three unit
    bonds. Gap 2 has the straight pair (0, 2) at distance 2 and the pair across
    the corner (1, 3) at sqrt(2), so <r> = (2 + sqrt 2) / 2 and
    <r^2> = (4 + 2) / 2 = 3, while the square of the mean is
    (2 + sqrt 2)^2 / 4 = 2.914... Gap 3 is the single end-to-end pair at sqrt(5).
    """
    positions = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [2, 1, 0]]
    chain = chain_module.Chain(
        lattice_grid=np.zeros((9, 9, 9), dtype=np.int32), dimensions=[9, 9, 9],
        sequence="AAAA", int_seq=[1] * 4, LR_int_seq=[1] * 4, LR_IDX=[],
        chainID=1, chainType=0, chain_positions=positions)
    chain.analysis_update_internal_scaling()

    np.testing.assert_allclose(chain.analysis_get_cumulative_internal_scaling(),
                               np.array([1.0, (2.0 + np.sqrt(2.0)) / 2.0, np.sqrt(5.0)]),
                               rtol=1e-12)
    np.testing.assert_allclose(chain.analysis_get_internal_scaling_squared(),
                               np.array([1.0, 3.0, 5.0]), rtol=1e-12)
