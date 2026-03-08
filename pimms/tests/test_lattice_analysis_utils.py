import warnings

import numpy as np
import pytest

from pimms import lattice_analysis_utils
from pimms.latticeExceptions import AnalysisRoutineException


class _DummyChain:
    def __init__(self, chain_id, positions):
        self.chainID = chain_id
        self._positions = positions

    def get_ordered_positions(self):
        return self._positions


class _DummyLatticeObject:
    def __init__(self, chains):
        self.grid = np.zeros((5, 5), dtype=np.int32)
        self.chains = chains


def test_get_inter_position_distance_2d_with_and_without_pbc():
    p1 = [0, 0]
    p2 = [9, 0]
    dims = [10, 10]

    no_pbc = lattice_analysis_utils.get_inter_position_distance(p1, p2, dims, pbc_correction=False)
    with_pbc = lattice_analysis_utils.get_inter_position_distance(p1, p2, dims, pbc_correction=True)

    assert np.isclose(no_pbc, 9.0)
    assert np.isclose(with_pbc, 1.0)


def test_get_inter_position_distance_3d_with_pbc():
    p1 = [0, 0, 0]
    p2 = [9, 9, 9]
    dims = [10, 10, 10]

    d = lattice_analysis_utils.get_inter_position_distance(p1, p2, dims, pbc_correction=True)
    assert np.isclose(d, np.sqrt(3.0))


def test_get_inter_position_distances_length_mismatch_raises():
    with pytest.raises(AnalysisRoutineException):
        lattice_analysis_utils.get_inter_position_distances([[0, 0]], [[0, 0], [1, 1]], [10, 10])


def test_get_inter_position_distances_matches_scalar_without_pbc():
    p1s = [[0, 0], [2, 2], [5, 1]]
    p2s = [[3, 4], [2, 5], [1, 1]]
    dims = [10, 10]

    vector = lattice_analysis_utils.get_inter_position_distances(p1s, p2s, dims, pbc_correction=False)
    scalar = np.array(
        [
            lattice_analysis_utils.get_inter_position_distance(a, b, dims, pbc_correction=False)
            for a, b in zip(p1s, p2s)
        ]
    )

    assert np.allclose(vector, scalar)


def test_get_cluster_distribution_orders_clusters(monkeypatch):
    chains = {
        10: _DummyChain(10, [[0, 0]]),
        20: _DummyChain(20, [[1, 0]]),
        30: _DummyChain(30, [[4, 4]]),
    }
    grid = np.zeros((5, 5), dtype=np.int32)

    mapping = {10: {10, 20}, 20: {10, 20}, 30: {30}}

    def fake_cc(chain_id, lattice_grid, chain_dict, useChains=True):
        return mapping[chain_id]

    monkeypatch.setattr(lattice_analysis_utils.lattice_utils, "get_all_chains_in_connected_component", fake_cc)

    out = lattice_analysis_utils.get_cluster_distribution(grid, chains)
    assert out[0] == {10, 20}
    assert out[1] == {30}


def test_get_lr_cluster_distribution_orders_clusters(monkeypatch):
    chains = {1: _DummyChain(1, [[0, 0]]), 2: _DummyChain(2, [[1, 0]]), 3: _DummyChain(3, [[2, 0]])}
    lattice_obj = _DummyLatticeObject(chains)

    mapping = {1: {1, 2}, 2: {1, 2}, 3: {3}}

    def fake_lr(chain_id, lo):
        return mapping[chain_id]

    monkeypatch.setattr(lattice_analysis_utils.lattice_utils, "get_all_chains_in_long_range_cluster", fake_lr)

    out = lattice_analysis_utils.get_LR_cluster_distribution(lattice_obj)
    assert out[0] == {1, 2}
    assert out[1] == {3}


def test_get_eigenvalues_of_t_matrix_nonnegative():
    positions = [[0, 0], [2, 0], [4, 0]]
    dims = [20, 20]

    eig, norm = lattice_analysis_utils.get_eigenvalues_of_the_T_matrix(positions, dims, pbc_correction=False)
    eig = np.real(eig)
    assert eig.shape[0] == 2
    assert np.all(eig >= -1e-12)
    assert norm.shape == (2, 2)


def test_get_polymeric_properties_non_degenerate_values_2d():
    positions = [[0, 0], [2, 0], [4, 0]]
    dims = [20, 20]
    rg, asph = lattice_analysis_utils.get_polymeric_properties(positions, dims, pbc_correction=False)

    assert rg > 0
    assert 0 <= asph <= 1.0


def test_get_polymeric_properties_non_degenerate_values_3d():
    positions = [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 2]]
    dims = [20, 20, 20]
    rg, asph = lattice_analysis_utils.get_polymeric_properties(positions, dims, pbc_correction=False)

    assert rg > 0
    assert np.isfinite(asph)


def test_get_polymeric_properties_degenerate_2d_no_runtime_warning():
    positions = [[5, 5], [5, 5], [5, 5]]
    dimensions = [20, 20]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rg, asph = lattice_analysis_utils.get_polymeric_properties(positions, dimensions)

    assert np.isfinite(rg)
    assert np.isfinite(asph)
    assert abs(rg) < 1e-10
    assert abs(asph) < 1e-10
    assert not any(issubclass(w.category, RuntimeWarning) for w in caught)


def test_get_polymeric_properties_degenerate_3d_no_runtime_warning():
    positions = [[2, 2, 2], [2, 2, 2], [2, 2, 2]]
    dimensions = [20, 20, 20]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rg, asph = lattice_analysis_utils.get_polymeric_properties(positions, dimensions)

    assert np.isfinite(rg)
    assert np.isfinite(asph)
    assert abs(rg) < 1e-10
    assert abs(asph) < 1e-10
    assert not any(issubclass(w.category, RuntimeWarning) for w in caught)


def test_extract_positions_from_clusters():
    chain_dict = {
        1: _DummyChain(1, [[0, 0], [1, 0]]),
        2: _DummyChain(2, [[2, 2]]),
        3: _DummyChain(3, [[3, 3]]),
    }
    clusters = [[1, 2], [3]]

    out = lattice_analysis_utils.extract_positions_from_clusters(clusters, chain_dict)
    assert out == [[[0, 0], [1, 0], [2, 2]], [[3, 3]]]


def test_extract_cluster_polymeric_properties_empty_returns_empty():
    out = lattice_analysis_utils.extract_cluster_polymeric_properties([], dimensions=[10, 10])
    assert out == []


def test_extract_cluster_polymeric_properties_dynamic_dimensions(monkeypatch):
    called_dims = []

    def fake_props(cluster, dimensions):
        called_dims.append(dimensions)
        return [1.0, 0.5]

    monkeypatch.setattr(lattice_analysis_utils, "get_polymeric_properties", fake_props)

    clusters = [[[1, 1], [3, 2], [4, 5]]]
    out = lattice_analysis_utils.extract_cluster_polymeric_properties(clusters, dimensions=False)

    assert out == [[1.0, 0.5]]
    assert called_dims == [[14, 15]]


def test_correct_cluster_positions_to_single_image_uses_threshold_1(monkeypatch):
    calls = []

    def fake_convert(cluster, dimensions, space_threshold):
        calls.append(space_threshold)
        return cluster

    monkeypatch.setattr(lattice_analysis_utils.cluster_utils, "convert_positions_to_single_image_snakesearch", fake_convert)

    clusters = [[[0, 0], [1, 0]]]
    out = lattice_analysis_utils.correct_cluster_positions_to_single_image(clusters, [10, 10])

    assert out == clusters
    assert calls == [1]


def test_correct_lr_cluster_positions_to_single_image_uses_threshold_2(monkeypatch):
    calls = []

    def fake_convert(cluster, dimensions, space_threshold):
        calls.append(space_threshold)
        return cluster

    monkeypatch.setattr(lattice_analysis_utils.cluster_utils, "convert_positions_to_single_image_snakesearch", fake_convert)

    clusters = [[[0, 0], [1, 0]]]
    out = lattice_analysis_utils.correct_LR_cluster_positions_to_single_image(clusters, [10, 10])

    assert out == clusters
    assert calls == [2]


def test_compute_cluster_gross_properties_handles_qhull_error():
    # Fewer than 4 non-coplanar points in 3D triggers Qhull error.
    clusters = [np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])]
    out = lattice_analysis_utils.compute_cluster_gross_properties(clusters)
    assert out == [[-1, -1, -1]]


def test_compute_cluster_gross_properties_valid_tetrahedron():
    tetra = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    out = lattice_analysis_utils.compute_cluster_gross_properties([tetra])
    vol, sa, den = out[0]

    assert vol > 0
    assert sa > 0
    assert np.isclose(den, 4.0 / vol)


def test_compute_cluster_radial_density_profile_min_cluster_threshold_skips():
    clusters = [np.array([[1, 1], [2, 1]])]
    out = lattice_analysis_utils.compute_cluster_radial_density_profile(
        clusters, [10, 10], minimum_cluster_size_in_beads=10
    )
    assert out == []


def test_compute_cluster_radial_density_profile_single_bead_2d_and_3d():
    out_2d = lattice_analysis_utils.compute_cluster_radial_density_profile(
        [np.array([[2, 2]])], [10, 10], minimum_cluster_size_in_beads=None
    )
    out_3d = lattice_analysis_utils.compute_cluster_radial_density_profile(
        [np.array([[2, 2, 2]])], [8, 8, 8], minimum_cluster_size_in_beads=None
    )

    assert len(out_2d) == 1
    assert len(out_3d) == 1
    assert all(val == 0 for val in out_2d[0])
    assert all(val == 0 for val in out_3d[0])
