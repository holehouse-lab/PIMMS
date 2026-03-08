import numpy as np
import pytest

from pimms import numpy_utils


def test_position_in_list_finds_matching_position_with_python_lists():
    assert numpy_utils.position_in_list([1, 2, 3], [[0, 0, 0], [1, 2, 3]]) is True


def test_position_in_list_returns_false_when_absent():
    assert numpy_utils.position_in_list([1, 2, 4], [[0, 0, 0], [1, 2, 3]]) is False


def test_position_in_list_works_with_numpy_arrays():
    pos = np.array([2, 5])
    pos_list = [np.array([0, 0]), np.array([2, 5])]
    assert numpy_utils.position_in_list(pos, pos_list) is True


def test_randneg_returns_only_signed_variants():
    seen = {numpy_utils.randneg(7) for _ in range(200)}
    assert seen.issubset({-7, 7})
    assert seen == {-7, 7}


def test_tetrahedron_volume_batch_and_single_consistency():
    # Unit tetrahedron with known volume 1/6.
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    c = np.array([0.0, 0.0, 1.0])
    d = np.array([0.0, 0.0, 0.0])

    single = numpy_utils.tetrahedron_volume(a, b, c, d)
    assert np.allclose(single, np.array([1.0 / 6.0]))

    batch = numpy_utils.tetrahedron_volume(
        np.array([a, a]),
        np.array([b, b]),
        np.array([c, c]),
        np.array([d, d]),
    )
    assert np.allclose(batch, np.array([1.0 / 6.0, 1.0 / 6.0]))


def test_tetrahedron_volume_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="same shape"):
        numpy_utils.tetrahedron_volume(
            np.array([[1.0, 0.0, 0.0]]),
            np.array([[0.0, 1.0, 0.0]]),
            np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
            np.array([[0.0, 0.0, 0.0]]),
        )


def test_find_nearest_accepts_python_list_input():
    idx, value = numpy_utils.find_nearest([10.0, 1.0, 6.0], 5.4)
    assert idx == 2
    assert value == 6.0


def test_find_nearest_rejects_empty_input():
    with pytest.raises(ValueError, match="empty"):
        numpy_utils.find_nearest([], 1.0)
