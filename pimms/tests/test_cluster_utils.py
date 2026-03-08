import numpy as np
import pytest

from pimms import cluster_utils


def test_find_local_basic_contact_no_wrap():
    target = [5, 5]
    positions = [[5, 6], [8, 8], [4, 4]]

    in_contact, in_contact_si = cluster_utils.find_local(target, positions, [10, 10], space_threshold=1)

    assert in_contact == [[5, 6], [4, 4]]
    assert in_contact_si == [[5, 6], [4, 4]]


def test_find_local_detects_contact_across_pbc_boundary():
    target = [0, 0]
    positions = [[9, 9], [2, 2]]

    in_contact, in_contact_si = cluster_utils.find_local(target, positions, [10, 10], space_threshold=1)

    assert in_contact == [[9, 9]]
    assert in_contact_si == [[-1, -1]]


def test_find_local_preserves_original_target_image_offset():
    # Target is in a translated single-image frame; returned SI contacts should match that frame.
    target = [10, 10]
    positions = [[9, 9], [2, 2]]

    in_contact, in_contact_si = cluster_utils.find_local(target, positions, [10, 10], space_threshold=1)

    assert in_contact == [[9, 9]]
    assert in_contact_si == [[9, 9]]


def test_convert_positions_to_single_image_empty_returns_empty_list():
    assert cluster_utils.convert_positions_to_single_image_snakesearch([], [10, 10]) == []


def test_convert_positions_to_single_image_connected_cluster_keeps_local_connectivity():
    positions = [[9, 0], [0, 0], [1, 0], [1, 1]]

    out = cluster_utils.convert_positions_to_single_image_snakesearch(positions, [10, 10], space_threshold=1)

    assert isinstance(out, np.ndarray)
    assert out.shape == (len(positions), 2)

    # Every original point should still have at least one local neighbor in SI coordinates.
    # This avoids assuming a specific bead order from the search algorithm.
    out_list = out.tolist()
    for p in out_list:
        neighbors = 0
        for q in out_list:
            if p is q:
                continue
            if max(abs(p[0] - q[0]), abs(p[1] - q[1])) <= 1:
                neighbors += 1
        assert neighbors >= 1


def test_convert_positions_to_single_image_disconnected_cluster_raises_clear_error():
    with pytest.raises(ValueError, match="connected"):
        cluster_utils.convert_positions_to_single_image_snakesearch([[0, 0], [5, 5]], [10, 10], space_threshold=1)


def test_build_interface_envelope_pairs_2d_aggregates_nonempty_sites(monkeypatch):
    def fake_pairs_2d(x, y, xdim, ydim, grid):
        if (x, y) == (1, 1):
            return [[1, 1, 1, 2]]
        if (x, y) == (3, 3):
            return [[3, 3, 3, 4], [3, 3, 4, 3]]
        return []

    monkeypatch.setattr(cluster_utils.hyperloop, "get_unique_interface_pairs_2D", fake_pairs_2d)

    positions = [[0, 0], [1, 1], [2, 2], [3, 3]]
    grid = np.zeros((5, 5), dtype=np.int32)

    out = cluster_utils.build_interface_envelope_pairs(positions, [5, 5], grid)

    assert out.shape == (3, 4)
    assert out.tolist() == [[1, 1, 1, 2], [3, 3, 3, 4], [3, 3, 4, 3]]


def test_build_interface_envelope_pairs_3d_aggregates_nonempty_sites_when_first_is_nonempty(monkeypatch):
    def fake_pairs_3d(x, y, z, xdim, ydim, zdim, grid):
        if (x, y, z) == (1, 1, 1):
            return [[1, 1, 1, 1, 1, 2]]
        if (x, y, z) == (4, 4, 4):
            return [[4, 4, 4, 4, 5, 4]]
        return []

    monkeypatch.setattr(cluster_utils.hyperloop, "get_unique_interface_pairs_3D", fake_pairs_3d)

    positions = [[1, 1, 1], [2, 2, 2], [4, 4, 4]]
    grid = np.zeros((6, 6, 6), dtype=np.int32)

    out = cluster_utils.build_interface_envelope_pairs(positions, [6, 6, 6], grid)

    assert out.shape == (2, 6)
    assert out.tolist() == [[1, 1, 1, 1, 1, 2], [4, 4, 4, 4, 5, 4]]


def test_build_interface_envelope_pairs_all_empty_returns_empty_array(monkeypatch):
    monkeypatch.setattr(cluster_utils.hyperloop, "get_unique_interface_pairs_2D", lambda x, y, xdim, ydim, grid: [])

    out = cluster_utils.build_interface_envelope_pairs([[0, 0], [1, 1]], [5, 5], np.zeros((5, 5), dtype=np.int32))
    assert out.shape == (0, 4)


def test_build_interface_envelope_pairs_3d_does_not_duplicate_first_nonempty(monkeypatch):
    def fake_pairs_3d(x, y, z, xdim, ydim, zdim, grid):
        if (x, y, z) == (1, 1, 1):
            return []
        if (x, y, z) == (2, 2, 2):
            return [[2, 2, 2, 2, 2, 3]]
        return []

    monkeypatch.setattr(cluster_utils.hyperloop, "get_unique_interface_pairs_3D", fake_pairs_3d)

    positions = [[1, 1, 1], [2, 2, 2]]
    out = cluster_utils.build_interface_envelope_pairs(positions, [6, 6, 6], np.zeros((6, 6, 6), dtype=np.int32))

    # Intended behavior is one pair (non-duplicated).
    assert out.shape == (1, 6)


def test_build_interface_envelope_pairs_safe_and_slow_returns_expected_shape_2d(monkeypatch):
    monkeypatch.setattr(
        cluster_utils.lattice_utils,
        "build_envelope_pairs",
        lambda positions, dimensions: np.array([[[0, 0], [0, 1]]]),
    )

    out = cluster_utils.build_interface_envelope_pairs_safe_and_slow([[0, 0]], [5, 5])
    assert out.shape == (1, 2, 2)


def test_build_interface_envelope_pairs_rejects_unsupported_dimensionality():
    with pytest.raises(ValueError, match="dimensionality"):
        cluster_utils.build_interface_envelope_pairs([[0]], [5], np.zeros((5,), dtype=np.int32))
