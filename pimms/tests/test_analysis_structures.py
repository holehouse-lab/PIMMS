import numpy as np
import pytest

from pimms import analysis_structures
from pimms.latticeExceptions import AnalysisStructureException


# ---------------------------------------------------------------------------
# the validation paths must raise the documented exception
#
# Regression tests for a real bug: analysis_structures.py raised
# AnalysisStructureException on three validation paths but never imported it, so
# every one of them died with `NameError: name 'AnalysisStructureException' is not
# defined` instead. The failure was invisible because these paths only fire on a
# size mismatch, which is exactly the "this should never happen" case the checks
# exist to report clearly.
# ---------------------------------------------------------------------------

def test_internal_scaling_rejects_mismatched_profile():
    acc = analysis_structures.InternalScaling(seqlen=6)      # gaps 1..4

    with pytest.raises(AnalysisStructureException):
        acc.update_internal_scaling({1: 1.0, 2: 2.0})


def test_internal_scaling_squared_rejects_mismatched_profile():
    acc = analysis_structures.InternalScalingSquared(seqlen=6)

    with pytest.raises(AnalysisStructureException):
        acc.update_internal_scaling({1: 1.0, 2: 2.0})


def test_distance_map_rejects_non_array_and_wrong_shape():
    dmap = analysis_structures.DistanceMap(seqlen=4)

    with pytest.raises(AnalysisStructureException):
        dmap.update_distance_map([[0, 1], [1, 0]])           # not an ndarray

    with pytest.raises(AnalysisStructureException):
        dmap.update_distance_map(np.zeros((3, 3)))           # wrong shape


# ---------------------------------------------------------------------------
# and the happy paths still accumulate a running mean
# ---------------------------------------------------------------------------

def test_internal_scaling_accumulates_a_running_mean():
    acc = analysis_structures.InternalScaling(seqlen=5)       # gaps 1, 2, 3
    acc.update_internal_scaling({1: 1.0, 2: 2.0, 3: 3.0})
    acc.update_internal_scaling({1: 3.0, 2: 4.0, 3: 5.0})

    assert acc.count == 2
    assert acc.get_internal_scaling_array() == pytest.approx([2.0, 3.0, 4.0])


def test_internal_scaling_squared_accumulates_the_mean_of_the_squares():
    acc = analysis_structures.InternalScalingSquared(seqlen=5)
    acc.update_internal_scaling({1: 1.0, 2: 2.0, 3: 3.0})
    acc.update_internal_scaling({1: 3.0, 2: 4.0, 3: 5.0})

    # mean of the SQUARES, not the square of the mean
    assert acc.get_internal_scaling_array() == pytest.approx([5.0, 10.0, 17.0])


def test_distance_map_accumulates_a_running_mean():
    dmap = analysis_structures.DistanceMap(seqlen=3)
    dmap.update_distance_map(np.full((3, 3), 2.0))
    dmap.update_distance_map(np.full((3, 3), 4.0))

    assert dmap.count == 2
    assert dmap.get_distance_map() == pytest.approx(np.full((3, 3), 3.0))
