import numpy as np
import pytest

from pimms import lattice_utils
from pimms.latticeExceptions import (
    ChainConnectivityError,
    ChainDeletionFailure,
    ChainInsertionFailure,
    ClusterSizeThresholdException,
    LatticeUtilsException,
    ResidueAugmentException,
    RotationException,
)


class _DummyChain:
    def __init__(self, chain_id, positions, lr_binary=None):
        self.chainID = chain_id
        self._positions = positions
        self._lr_binary = lr_binary if lr_binary is not None else np.zeros(len(positions), dtype=np.int32)

    def get_ordered_positions(self, center_positions=False):
        return self._positions

    def get_LR_binary_array(self):
        return self._lr_binary


class _DummyLattice:
    def __init__(self, dimensions, chains, grid=None, type_grid=None):
        self.dimensions = dimensions
        self.chains = chains
        self.grid = np.zeros(dimensions, dtype=np.int32) if grid is None else grid
        self.type_grid = np.zeros(dimensions, dtype=np.int32) if type_grid is None else type_grid


class _FakeTraj:
    def __init__(self):
        self.topology = object()
        self.time = np.array([0.0])
        self.unitcell_lengths = np.array([[1.0, 1.0, 1.0]])
        self.unitcell_angles = np.array([[90.0, 90.0, 90.0]])
        self.saved = []
        self.saved_xtc = []

    def join(self, other):
        return self

    def save(self, filename):
        self.saved.append(filename)

    def save_xtc(self, filename):
        self.saved_xtc.append(filename)


def test_same_sites_2d_and_3d():
    assert lattice_utils.same_sites([1, 2], [1, 2]) is True
    assert lattice_utils.same_sites([1, 2], [2, 1]) is False
    assert lattice_utils.same_sites([1, 2, 3], [1, 2, 3]) is True
    assert lattice_utils.same_sites([1, 2, 3], [1, 2, 4]) is False


def test_same_sites_rejects_unsupported_dimensionality():
    with pytest.raises(LatticeUtilsException, match="Unsupported dimensionality"):
        lattice_utils.same_sites([1], [1])


def test_get_real_distance_delegates_to_analysis(monkeypatch):
    monkeypatch.setattr(lattice_utils.lattice_analysis_utils, "get_inter_position_distance", lambda a, b, d: 12.5)
    assert lattice_utils.get_real_distance([0, 0], [1, 1], [10, 10]) == 12.5


def test_get_dimensions_and_pbc_convert():
    grid = np.zeros((5, 7, 9), dtype=np.int32)
    assert lattice_utils.get_dimensions(grid) == (5, 7, 9)
    assert lattice_utils.pbc_convert([6, -1, 10], [5, 7, 9]) == [1, 6, 1]


def test_pbc_correct_both_directions():
    a, b = lattice_utils.pbc_correct([9, 0], [1, 0], [10, 10])
    assert a == [9, 0]
    assert b == [11, 0]

    a2, b2 = lattice_utils.pbc_correct([1, 0], [9, 0], [10, 10])
    assert b2 == [-1, 0]


def test_straddle_and_center_positions(monkeypatch):
    assert lattice_utils.do_positions_stradle_pbc_boundary([[0, 0], [1, 0], [2, 0]]) is False
    assert lattice_utils.do_positions_stradle_pbc_boundary([[0, 0], [9, 0]]) is True

    monkeypatch.setattr(lattice_utils, "center_of_mass_from_positions", lambda pos, dims: [0, 0])
    centered = lattice_utils.center_positions([[1, 1], [2, 2]], [10, 10])
    assert centered == [[6.0, 6.0], [7.0, 7.0]]


def test_convert_chain_to_single_image_handles_pbc():
    chain = [[9, 0], [0, 0], [1, 0]]
    out = lattice_utils.convert_chain_to_single_image(chain, [10, 10])
    assert out == [[9, 0], [10, 0], [11, 0]]


def test_adjacent_site_wrappers(monkeypatch):
    monkeypatch.setattr(lattice_utils.hyperloop, "get_adjacent_sites_2D", lambda a, b, x, y, e: np.array([[1, 2], [2, 2]]))
    monkeypatch.setattr(
        lattice_utils.hyperloop,
        "get_adjacent_sites_3D",
        lambda a, b, c, x, y, z, e: np.array([[1, 2, 3], [2, 2, 3]]),
    )

    a2 = lattice_utils.get_adjacent_sites_2D(0, 0, [10, 10])
    a3 = lattice_utils.get_adjacent_sites_3D(0, 0, 0, [10, 10, 10])
    assert a2.shape == (2, 2)
    assert a3.shape == (2, 3)


def test_find_nearest_position_and_empty_error(monkeypatch):
    with pytest.raises(LatticeUtilsException):
        lattice_utils.find_nearest_position([0, 0], [], [10, 10])

    monkeypatch.setattr(lattice_utils, "get_real_distance", lambda t, p, d: p[0])
    idx, dist = lattice_utils.find_nearest_position([0, 0], [[5, 0], [2, 0], [3, 0]], [10, 10])
    assert idx == 1
    assert dist == 2


def test_get_empty_site_adjacent_and_global(monkeypatch):
    grid = np.zeros((4, 4), dtype=np.int32)
    grid[0, 1] = 5

    monkeypatch.setattr(lattice_utils, "get_adjacent_sites_2D", lambda x, y, d: [[0, 1], [1, 0], [0, 3]])
    site, ok = lattice_utils.get_empty_site(grid, adjacentTo=[0, 0], hardwall=True)
    assert ok is True
    assert site in [[1, 0]]

    # global empty-site search
    monkeypatch.setattr(lattice_utils.random, "randint", lambda lo, hi: lo)
    g2 = np.zeros((2, 2), dtype=np.int32)
    assert lattice_utils.get_empty_site(g2) == [0, 0]


def test_get_empty_site_adjacent_failure_returns_dimension_matched_sentinel(monkeypatch):
    grid = np.ones((4, 4), dtype=np.int32)
    monkeypatch.setattr(lattice_utils, "get_adjacent_sites_2D", lambda x, y, d: [[0, 1], [1, 0]])

    site, ok = lattice_utils.get_empty_site(grid, adjacentTo=[0, 0], hardwall=False)
    assert ok is False
    assert site == [-1, -1]


def test_get_empty_site_raises_fast_when_grid_is_full():
    grid = np.ones((3, 3), dtype=np.int32)
    with pytest.raises(LatticeUtilsException, match="fully occupied"):
        lattice_utils.get_empty_site(grid)


def test_insert_chain_and_placement_deletion(monkeypatch):
    grid = np.zeros((5, 5), dtype=np.int32)

    seq = iter([
        [0, 0],  # start
        ([1, 0], True),
        ([2, 0], True),
    ])

    def fake_get_empty_site(lg, adjacentTo=None, hardwall=False):
        v = next(seq)
        if adjacentTo is None:
            return v
        return v

    monkeypatch.setattr(lattice_utils, "get_empty_site", fake_get_empty_site)
    out = lattice_utils.insert_chain(7, 3, grid)
    assert out == [[0, 0], [1, 0], [2, 0]]
    assert np.all(grid[[0, 1, 2], 0] == 7)

    # placement safe failure
    with pytest.raises(ChainInsertionFailure):
        lattice_utils.place_chain_by_position([[0, 0]], grid, 8, safe=True)

    lattice_utils.delete_chain_by_ID(7, grid)
    assert np.all(grid == 0)

    lattice_utils.place_chain_by_position([[0, 0], [0, 1]], grid, 3, safe=False)
    with pytest.raises(ChainDeletionFailure):
        lattice_utils.delete_chain_by_position([[0, 0]], grid, chainID=2)

    lattice_utils.delete_chain_by_position([[0, 0], [0, 1]], grid)
    assert np.all(grid == 0)


def test_grid_get_set_helpers_and_3d_wrapper(monkeypatch):
    g2 = np.zeros((3, 3), dtype=np.int32)
    lattice_utils.set_gridvalue([1, 1], 4, g2)
    assert lattice_utils.get_gridvalue([1, 1], g2) == 4
    assert lattice_utils.get_gridvalue_2D([1, 1], g2) == 4

    g3 = np.zeros((3, 3, 3), dtype=np.int32)
    lattice_utils.set_gridvalue([1, 1, 1], 9, g3)
    assert lattice_utils.get_gridvalue([1, 1, 1], g3) == 9
    monkeypatch.setattr(lattice_utils.hyperloop, "get_gridvalue_3D", lambda l, x, y, z: 42)
    assert lattice_utils.get_gridvalue_3D([1, 1, 1], g3) == 42


def test_grid_get_set_reject_unsupported_dimensionality():
    g4 = np.zeros((2, 2, 2, 2), dtype=np.int32)

    with pytest.raises(LatticeUtilsException, match="Unsupported lattice dimensionality"):
        lattice_utils.get_gridvalue([0, 0, 0, 0], g4)

    with pytest.raises(LatticeUtilsException, match="Unsupported position dimensionality"):
        lattice_utils.set_gridvalue([0, 0, 0, 0], 1, g4)


def test_build_envelope_pairs_2d_and_3d(monkeypatch):
    p = [[1, 1]]

    monkeypatch.setattr(
        lattice_utils.inner_loops,
        "extract_SR_pairs_from_position_2D",
        lambda pos, x, y: np.array([[1, 1, 1, 2], [1, 1, 1, 2]], dtype=np.int32),
    )
    out2 = lattice_utils.build_envelope_pairs(p, [5, 5], hardwall=False)
    assert out2.shape == (1, 2, 2)

    monkeypatch.setattr(
        lattice_utils.inner_loops,
        "extract_SR_pairs_from_position_3D",
        lambda pos, x, y, z: np.array([[1, 1, 1, 1, 1, 2], [1, 1, 1, 1, 1, 2]], dtype=np.int32),
    )
    out3 = lattice_utils.build_envelope_pairs([[1, 1, 1]], [5, 5, 5], hardwall=False)
    assert out3.shape == (1, 2, 3)


def test_build_envelope_pairs_handles_empty_positions():
    out2 = lattice_utils.build_envelope_pairs([], [5, 5], hardwall=False)
    out3 = lattice_utils.build_envelope_pairs([], [5, 5, 5], hardwall=False)

    assert out2.shape == (0, 2, 2)
    assert out3.shape == (0, 2, 3)


def test_build_all_envelope_pairs_2d_and_3d(monkeypatch):
    monkeypatch.setattr(
        lattice_utils.inner_loops,
        "extract_SR_and_LR_pairs_from_position_2D",
        lambda pos, lr, tl, x, y: (
            np.array([[0, 0, 0, 1], [0, 0, 0, 1]], dtype=np.int32),
            np.array([[0, 0, 1, 1]], dtype=np.int32),
            np.array([], dtype=np.int32).reshape(0, 4),
        ),
    )
    sr, lr, slr = lattice_utils.build_all_envelope_pairs([[0, 0]], [1], np.zeros((2, 2), dtype=np.int32), [2, 2])
    assert sr.shape == (1, 2, 2)
    assert lr.shape == (1, 2, 2)
    assert slr.shape == (0, 2, 2)

    monkeypatch.setattr(
        lattice_utils.inner_loops,
        "extract_SR_and_LR_pairs_from_position_3D",
        lambda pos, lr, tl, x, y, z: (
            np.array([[0, 0, 0, 0, 0, 1]], dtype=np.int32),
            np.array([], dtype=np.int32).reshape(0, 6),
            np.array([], dtype=np.int32).reshape(0, 6),
        ),
    )
    sr3, lr3, slr3 = lattice_utils.build_all_envelope_pairs(
        [[0, 0, 0]], [1], np.zeros((2, 2, 2), dtype=np.int32), [2, 2, 2]
    )
    assert sr3.shape == (1, 2, 3)
    assert lr3.shape == (0, 2, 3)
    assert slr3.shape == (0, 2, 3)


def test_build_all_envelope_pairs_handles_empty_positions():
    sr2, lr2, slr2 = lattice_utils.build_all_envelope_pairs([], np.array([], dtype=np.int32), np.zeros((2, 2), dtype=np.int32), [2, 2])
    assert sr2.shape == (0, 2, 2)
    assert lr2.shape == (0, 2, 2)
    assert slr2.shape == (0, 2, 2)

    sr3, lr3, slr3 = lattice_utils.build_all_envelope_pairs([], np.array([], dtype=np.int32), np.zeros((2, 2, 2), dtype=np.int32), [2, 2, 2])
    assert sr3.shape == (0, 2, 3)
    assert lr3.shape == (0, 2, 3)
    assert slr3.shape == (0, 2, 3)


def test_connected_component_and_threshold(monkeypatch):
    g = np.zeros((4, 4), dtype=np.int32)
    g[0, 0] = 1
    g[0, 1] = 2
    g[3, 3] = 3

    chain_dict = {
        1: _DummyChain(1, [[0, 0]]),
        2: _DummyChain(2, [[0, 1]]),
        3: _DummyChain(3, [[3, 3]]),
    }

    def fake_envelope(positions, dimensions, hardwall=False):
        if [0, 0] in positions or [0, 1] in positions:
            return np.array([[[0, 0], [0, 1]]], dtype=np.int32)
        return np.array([], dtype=np.int32).reshape(0, 2, 2)

    monkeypatch.setattr(lattice_utils, "build_envelope_pairs", fake_envelope)

    chains = lattice_utils.get_all_chains_in_connected_component(1, g, chain_dict, useChains=True)
    assert set(chains) == {1, 2}

    with pytest.raises(ClusterSizeThresholdException):
        lattice_utils.get_all_chains_in_connected_component(1, g, chain_dict, threshold=1, useChains=True)


def test_long_range_cluster_component(monkeypatch):
    g = np.zeros((4, 4), dtype=np.int32)
    g[0, 0] = 1
    g[0, 1] = 2
    type_grid = np.zeros((4, 4), dtype=np.int32)
    chains = {
        1: _DummyChain(1, [[0, 0]], lr_binary=np.array([1], dtype=np.int32)),
        2: _DummyChain(2, [[0, 1]], lr_binary=np.array([1], dtype=np.int32)),
    }
    lat = _DummyLattice([4, 4], chains, grid=g, type_grid=type_grid)

    def fake_build_all(positions, lr_binary, type_lattice, dimensions, hardwall=False):
        return (
            np.array([[[0, 0], [0, 1]]], dtype=np.int32),
            np.array([], dtype=np.int32).reshape(0, 2, 2),
            np.array([], dtype=np.int32).reshape(0, 2, 2),
        )

    monkeypatch.setattr(lattice_utils, "build_all_envelope_pairs", fake_build_all)
    out = lattice_utils.get_all_chains_in_long_range_cluster(1, lat)
    assert set(out) == {1, 2}


def test_center_of_mass_from_positions_2d_and_3d():
    com2 = lattice_utils.center_of_mass_from_positions([[0, 0], [0, 0], [0, 0]], [10, 10], on_lattice=True)
    com3 = lattice_utils.center_of_mass_from_positions([[0, 0, 0], [0, 0, 0]], [10, 10, 10], on_lattice=False)
    assert com2 == [0, 0]
    assert len(com3) == 3


def test_center_of_mass_rejects_empty_positions():
    with pytest.raises(LatticeUtilsException, match="positions list is empty"):
        lattice_utils.center_of_mass_from_positions([], [10, 10])


def test_insert_and_delete_residue_safe_and_unsafe():
    g = np.zeros((3, 3), dtype=np.int32)

    lattice_utils.insert_residue([1, 1], g, 7, safe=True)
    assert g[1, 1] == 7

    with pytest.raises(ResidueAugmentException):
        lattice_utils.insert_residue([1, 1], g, 8, safe=True)

    lattice_utils.insert_residue([1, 1], g, 8, safe=False)
    assert g[1, 1] == 8

    with pytest.raises(ResidueAugmentException):
        lattice_utils.delete_residue([1, 1], g, chainID=3)

    lattice_utils.delete_residue([1, 1], g, chainID=8)
    assert g[1, 1] == 0


def test_rotations_and_invalid_inputs():
    out = lattice_utils.run_rotation([[1, 0]], np.array([[0, -1], [1, 0]]))
    assert np.allclose(out[0], [0, 1])

    r2 = lattice_utils.rotate_positions_2D([[1, 0]], 90)
    assert np.allclose(r2[0], [0, 1])

    r3 = lattice_utils.rotate_positions_3D([[1, 0, 0]], "z", 90)
    assert np.allclose(r3[0], [0, 1, 0])

    with pytest.raises(RotationException):
        lattice_utils.rotate_positions_2D([[1, 0]], 45)

    with pytest.raises(RotationException):
        lattice_utils.rotate_positions_3D([[1, 0, 0]], "q", 90)


def test_pdb_wrappers(monkeypatch):
    calls = []
    monkeypatch.setattr(lattice_utils.pdb_utils, "initialize_pdb_file", lambda d, s, f: calls.append(("init", f)))
    monkeypatch.setattr(
        lattice_utils.pdb_utils,
        "build_pdb_file",
        lambda lo, sp, fn, write_connect=False, autocenter=False: calls.append(("build", fn, write_connect, autocenter)),
    )
    monkeypatch.setattr(lattice_utils.pdb_utils, "finalize_pdb_file", lambda f: calls.append(("fin", f)))

    lat = _DummyLattice([5, 5], {})
    lattice_utils.open_pdb_file([5, 5], 3.8, filename="a.pdb")
    lattice_utils.write_lattice_to_pdb(lat, 3.8, filename="b.pdb", write_connect=True, autocenter=True)
    lattice_utils.finish_pdb_file("c.pdb")

    assert calls[0] == ("init", "a.pdb")
    assert calls[1][0] == "build"
    assert calls[2] == ("fin", "c.pdb")


def test_xtc_helpers(monkeypatch):
    fake = _FakeTraj()

    monkeypatch.setattr(lattice_utils, "open_pdb_file", lambda *args, **kwargs: None)
    monkeypatch.setattr(lattice_utils, "write_lattice_to_pdb", lambda *args, **kwargs: None)
    monkeypatch.setattr(lattice_utils, "finish_pdb_file", lambda *args, **kwargs: None)
    monkeypatch.setattr(lattice_utils.os, "remove", lambda fn: None)

    def fake_load(*args, **kwargs):
        return fake

    monkeypatch.setattr(lattice_utils.md, "load", fake_load)

    lat = _DummyLattice([5, 5], {})
    lattice_utils.start_xtc_file(lat, 3.8, pdb_filename="START.pdb", xtc_filename="traj.xtc")
    assert "traj.xtc" in fake.saved_xtc

    lattice_utils.append_to_xtc_file(lat, 3.8, xtc_filename="traj.xtc", autocenter=False)
    assert "traj.xtc" in fake.saved


def test_append_non_redundant_update_master_and_save(monkeypatch):
    fake = _FakeTraj()
    monkeypatch.setattr(lattice_utils.md, "load", lambda *a, **k: fake)
    monkeypatch.setattr(
        lattice_utils.md,
        "Trajectory",
        lambda xyz, topology, time, unitcell_lengths, unitcell_angles: _FakeTraj(),
    )

    chains = {1: _DummyChain(1, [[0, 0], [1, 0]])}
    lat = _DummyLattice([5, 5], chains)

    lattice_utils.append_to_xtc_file_non_redundant(lat, 3.8, pdb_filename="START.pdb", xtc_filename="traj.xtc")
    assert "traj.xtc" in fake.saved

    out = lattice_utils.update_master_traj(lat, 3.8, None, pdb_filename="START.pdb")
    assert isinstance(out, _FakeTraj)

    lattice_utils.save_out_sim(fake, "out.xtc")
    assert "out.xtc" in fake.saved


def test_chain_connectivity_checks():
    lattice_utils.check_chain_connectivity(1, [[0, 0], [1, 0], [2, 0]], [5, 5], verbose=False)

    with pytest.raises(ChainConnectivityError):
        lattice_utils.check_chain_connectivity(1, [[0, 0], [3, 0]], [5, 5], verbose=False)


class _ChainObj:
    def __init__(self, positions):
        self._positions = positions

    def get_ordered_positions(self):
        return self._positions


def test_check_all_chain_connectivity():
    chain_dict = {1: _ChainObj([[0, 0], [1, 0]]), 2: _ChainObj([[2, 2], [2, 3]])}
    lattice_utils.check_all_chain_connectivity(chain_dict, [10, 10], verbose=False)
