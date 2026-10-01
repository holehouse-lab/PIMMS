"""
Regression tests for the small follow-up changes made while merging the
October 2026 deep-audit fixes (the pieces that fell between the areas the
individual fixes covered).
"""

from __future__ import annotations

import tracemalloc

import numpy as np
import pytest

from pimms import simulation


@pytest.mark.parametrize("n_chains", [1, 2, 7, 8, 9, 16, 17, 129, 1000])
@pytest.mark.parametrize("seqlen", [1, 2, 5, 33])
def test_mean_over_chains_is_bit_identical_to_the_stacked_mean(n_chains: int, seqlen: int) -> None:
    """The end-of-run distance map must not change by a single bit.

    DISTANCE_MAP.dat is written with four decimals, so a last-bit change in the
    mean could flip a printed digit at a rounding tie. The reference here is the
    expression the engine used before (stack every chain's map, then take the
    mean over the first axis), written out independently of the helper.
    """
    rng = np.random.default_rng(1000 * n_chains + seqlen)
    maps = [rng.random((seqlen, seqlen)) * rng.integers(1, 1000) for _ in range(n_chains)]

    stacked = np.asarray(maps).mean(axis=0)
    running = simulation._mean_over_chains(iter(maps))

    assert running.dtype == np.float64
    assert np.array_equal(stacked, running)
    # the inputs are the chains' own accumulators and must not be modified
    assert running is not maps[0]


def test_mean_over_chains_does_not_stack_the_maps() -> None:
    """Averaging must cost one map of extra memory, not one per chain.

    Forty 200 x 200 maps are 12.8 MB; the stacked mean allocated all of that
    again. The running sum needs one 320 kB array (plus the returned one).
    """
    seqlen = 200
    maps = [np.full((seqlen, seqlen), float(i)) for i in range(40)]
    one_map = maps[0].nbytes

    tracemalloc.start()
    try:
        result = simulation._mean_over_chains(iter(maps))
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()

    assert np.array_equal(result, np.full((seqlen, seqlen), 19.5))
    assert peak < 4 * one_map
