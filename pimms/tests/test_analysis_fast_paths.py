"""
Equivalence tests for the optimized analysis paths.

Each of these replaced a Python loop that dominated the runtime of a PIMMS analysis
step. They are only worth having if they compute *the same thing*, so every test here
pins the fast path against a naive reference implementation of what it replaced. Where
the replacement is exact the tests assert bit-level equality; the one place it is not
(the gyration tensor, where the summation order necessarily changes) is pinned far
tighter than the 4-decimal precision the results are ever written out at.
"""

import numpy as np
import pytest

from pimms import lattice_analysis_utils as lau
from pimms import lattice_utils


_BOXES = [[9, 9], [7, 11], [8, 8, 8], [6, 9, 12]]


def _random_positions(rng, dims, n):
    return [[int(rng.integers(0, d)) for d in dims] for _ in range(n)]


# ---------------------------------------------------------------------------
# the scalar distance helper (rewritten from numpy scalars to plain Python)
# ---------------------------------------------------------------------------

def _reference_distance(p1, p2, dims, pbc_correction=True):
    """The original numpy-scalar implementation."""
    x_dif = np.array(p1)[0] - np.array(p2)[0]
    y_dif = np.array(p1)[1] - np.array(p2)[1]
    if pbc_correction:
        if np.abs(x_dif) > dims[0] * 0.5:
            x_dif = dims[0] - np.abs(x_dif)
        if np.abs(y_dif) > dims[1] * 0.5:
            y_dif = dims[1] - np.abs(y_dif)
    if len(dims) == 3:
        z_dif = np.array(p1)[2] - np.array(p2)[2]
        if pbc_correction and np.abs(z_dif) > dims[2] * 0.5:
            z_dif = dims[2] - np.abs(z_dif)
        return np.sqrt(np.power(x_dif, 2) + np.power(y_dif, 2) + np.power(z_dif, 2))
    return np.sqrt(np.power(x_dif, 2) + np.power(y_dif, 2))


@pytest.mark.parametrize("dims", _BOXES)
@pytest.mark.parametrize("pbc", [True, False])
def test_scalar_distance_is_bit_identical_to_the_numpy_version(dims, pbc):
    rng = np.random.default_rng(0)
    for _ in range(400):
        p1 = _random_positions(rng, dims, 1)[0]
        p2 = _random_positions(rng, dims, 1)[0]
        assert float(lau.get_inter_position_distance(p1, p2, dims, pbc)) == \
            float(_reference_distance(p1, p2, dims, pbc))


def test_scalar_distance_handles_float_positions():
    rng = np.random.default_rng(1)
    dims = [12, 12, 12]
    for _ in range(200):
        p1 = list(rng.uniform(0, 12, 3))
        p2 = list(rng.uniform(0, 12, 3))
        assert float(lau.get_inter_position_distance(p1, p2, dims)) == \
            float(_reference_distance(p1, p2, dims))


def test_scalar_distance_takes_the_short_way_round_the_box():
    # 0 and 9 are one lattice step apart across the boundary of a 10-wide box
    assert lau.get_inter_position_distance([0, 0], [9, 0], [10, 10]) == pytest.approx(1.0)
    # ... and 9 apart if the correction is switched off
    assert lau.get_inter_position_distance([0, 0], [9, 0], [10, 10], False) == pytest.approx(9.0)


# ---------------------------------------------------------------------------
# vectorized distance matrix / internal scaling (replacing O(L^2) Python loops)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dims", _BOXES)
@pytest.mark.parametrize("n_pos", [2, 3, 11, 40])
def test_distance_matrix_matches_the_per_pair_helper(dims, n_pos):
    rng = np.random.default_rng(hash((tuple(dims), n_pos)) % 2**31)
    pos = _random_positions(rng, dims, n_pos)

    reference = np.zeros((n_pos, n_pos))
    for i in range(n_pos):
        for j in range(n_pos):
            reference[i][j] = lau.get_inter_position_distance(pos[i], pos[j], dims)

    assert np.array_equal(lau.get_distance_matrix(pos, dims), reference)


def test_distance_matrix_is_blocked_but_still_exact_for_a_long_chain():
    """The row blocking must not change any value (block size < chain length here)."""
    rng = np.random.default_rng(4)
    dims = [30, 30, 30]
    pos = _random_positions(rng, dims, 900)

    matrix = lau.get_distance_matrix(pos, dims)

    assert matrix.shape == (900, 900)
    assert np.array_equal(matrix, matrix.T)                 # symmetric
    for i, j in [(0, 899), (5, 500), (899, 3), (400, 401)]:
        assert matrix[i][j] == lau.get_inter_position_distance(pos[i], pos[j], dims)


@pytest.mark.parametrize("dims", _BOXES)
@pytest.mark.parametrize("n_pos", [3, 12, 45])
def test_internal_scaling_profile_matches_the_per_pair_loop(dims, n_pos):
    rng = np.random.default_rng(hash((tuple(dims), n_pos, 7)) % 2**31)
    pos = _random_positions(rng, dims, n_pos)

    reference = {}
    for gap in range(1, n_pos - 1):
        reference[gap] = np.mean([lau.get_inter_position_distance(pos[i], pos[i + gap], dims)
                                  for i in range(0, n_pos - gap)])

    gaps, means = lau.get_internal_scaling_profile(pos, dims)

    assert gaps == list(reference.keys())
    for gap, mean in zip(gaps, means):
        assert mean == reference[gap]


# ---------------------------------------------------------------------------
# gyration tensor (per-bead Python loop -> single vectorized pass, eig -> eigh)
# ---------------------------------------------------------------------------

def _reference_polymeric_properties(positions, dimensions):
    """The original per-bead accumulation, using the general eigensolver.

    NB: ``np.linalg.eig`` on the (symmetric) gyration tensor can return a *complex*
    array, because it makes no symmetry assumption and the matrix is only symmetric to
    within rounding. That is precisely the hazard the production code avoids by using
    ``eigh``, which is guaranteed real - so the imaginary part is discarded explicitly
    here rather than leaking a ComplexWarning out of the reference implementation.
    """
    com = lattice_utils.center_of_mass_from_positions(positions, dimensions, on_lattice=False)
    tensor = 0
    for pos in positions:
        _, corrected = lattice_utils.pbc_correct(com, pos, dimensions)
        delta = np.array(corrected) - np.array(com)
        tensor = tensor + np.outer(delta, delta)
    eig = np.real(np.linalg.eig(tensor / len(positions))[0])

    if len(dimensions) == 2:
        rg2 = max(0.0, eig[0] + eig[1])
        return [np.sqrt(rg2), 0.0 if rg2 <= 1e-12 else abs(eig[0] - eig[1]) / rg2]

    total = max(0.0, eig[0] + eig[1] + eig[2])
    denom = np.power(total, 2)
    asph = 0.0 if denom <= 1e-12 else 1 - 3 * (
        (eig[0] * eig[1] + eig[1] * eig[2] + eig[2] * eig[0]) / denom)
    return [np.sqrt(total), asph]


@pytest.mark.parametrize("dims", _BOXES)
@pytest.mark.parametrize("n_pos", [3, 25, 300])
def test_polymeric_properties_agree_far_beyond_output_precision(dims, n_pos):
    """Rg and asphericity are written at 3-4 decimal places; agreement here is ~1e-13.

    The vectorized tensor sums in a different order from the per-bead accumulation, so
    exact equality is not on offer - but the difference is many orders of magnitude
    below anything that reaches an output file.
    """
    rng = np.random.default_rng(hash((tuple(dims), n_pos, 13)) % 2**31)
    for _ in range(10):
        pos = _random_positions(rng, dims, n_pos)
        reference = _reference_polymeric_properties(pos, dims)
        fast = lau.get_polymeric_properties(pos, dims)

        assert float(fast[0]) == pytest.approx(float(reference[0]), abs=1e-9)
        assert float(fast[1]) == pytest.approx(float(reference[1]), abs=1e-9)


def test_polymeric_properties_are_real_valued():
    """eigh is used rather than eig, which can go complex on a near-symmetric matrix."""
    rng = np.random.default_rng(21)
    rg, asph = lau.get_polymeric_properties(_random_positions(rng, [10, 10, 10], 50), [10, 10, 10])

    assert isinstance(float(rg), float) and np.isfinite(rg)
    assert isinstance(float(asph), float) and np.isfinite(asph)


# ---------------------------------------------------------------------------
# envelope-pair deduplication
# ---------------------------------------------------------------------------

def test_unique_rows_deduplicates_via_the_packed_key_path():
    rows = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 3, 5]], dtype=np.int32)

    out = lattice_utils._unique_rows(rows)

    assert len(out) == 3
    assert sorted(map(tuple, out.tolist())) == [(1, 2, 3, 4), (1, 2, 3, 5), (5, 6, 7, 8)]


def test_unique_rows_falls_back_when_values_cannot_be_packed():
    # a span this large cannot be packed into an int64 across 6 columns, so the
    # np.void view fallback has to handle it
    rows = np.array([[0, 0, 0, 0, 0, 0],
                     [2 ** 40, 1, 2, 3, 4, 5],
                     [2 ** 40, 1, 2, 3, 4, 5]], dtype=np.int64)

    out = lattice_utils._unique_rows(rows)

    assert len(out) == 2
    assert sorted(map(tuple, out.tolist())) == [(0, 0, 0, 0, 0, 0), (2 ** 40, 1, 2, 3, 4, 5)]


def test_unique_rows_handles_an_empty_array():
    rows = np.empty((0, 4), dtype=np.int32)
    assert len(lattice_utils._unique_rows(rows)) == 0


@pytest.mark.parametrize("dims", [[8, 8], [7, 9, 11]])
def test_skipping_dedupe_preserves_the_set_of_sites_touched(dims):
    """The connected-component searches pass deduplicate=False.

    They only ask which sites are touched, so duplicates are harmless there - but the
    set of sites has to be identical to the deduplicated version, or the search would
    find a different cluster.
    """
    rng = np.random.default_rng(31)
    positions = _random_positions(rng, dims, 12)

    deduped = lattice_utils.build_envelope_pairs(positions, dims)
    raw = lattice_utils.build_envelope_pairs(positions, dims, deduplicate=False)

    assert len(raw) >= len(deduped)
    as_tuples = lambda arr: {tuple(pair.ravel().tolist()) for pair in arr}
    assert as_tuples(raw) == as_tuples(deduped)


# ---------------------------------------------------------------------------
# nearest-position search
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dims", _BOXES)
def test_find_nearest_position_matches_the_scalar_loop(dims):
    rng = np.random.default_rng(hash(tuple(dims)) % 2**31)
    for n in (1, 2, 30, 200):
        pts = _random_positions(rng, dims, n)
        target = _random_positions(rng, dims, 1)[0]

        best_idx, best_dist = -1, float("inf")
        for idx, pos in enumerate(pts):
            d = lau.get_inter_position_distance(target, pos, dims)
            if d < best_dist:
                best_idx, best_dist = idx, d

        idx, dist = lattice_utils.find_nearest_position(target, pts, dims)
        assert idx == best_idx
        assert float(dist) == float(best_dist)


def test_find_nearest_position_keeps_the_first_of_a_tie():
    """np.argmin returns the first minimum, matching the old strict-`<` comparison."""
    pts = [[5, 5], [1, 1], [1, 1], [9, 9]]
    idx, dist = lattice_utils.find_nearest_position([1, 1], pts, [10, 10])

    assert idx == 1
    assert dist == pytest.approx(0.0)
