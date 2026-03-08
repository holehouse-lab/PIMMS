import importlib.util
import pathlib
import sys
import types

import numpy as np
import pytest

import pimms.analysis_structures as analysis_structures
import pimms.latticeExceptions as latticeExceptions


def _load_chain_module_with_stubs(monkeypatch):
    """Load pimms.chain while stubbing heavy runtime dependencies."""
    config_stub = types.ModuleType("pimms.CONFIG")
    config_stub.NP_INT_TYPE = np.int32

    lattice_utils_stub = types.ModuleType("pimms.lattice_utils")
    lattice_analysis_utils_stub = types.ModuleType("pimms.lattice_analysis_utils")

    lattice_utils_stub.insert_chain = lambda chain_id, length, lattice_grid, default_start=None, hardwall=False: [
        [i, 0] for i in range(length)
    ]
    lattice_utils_stub.center_positions = lambda positions, dimensions: positions
    lattice_utils_stub.convert_chain_to_single_image = lambda positions, dimensions: positions
    lattice_utils_stub.do_positions_stradle_pbc_boundary = lambda positions: False
    lattice_utils_stub.center_of_mass_from_positions = lambda positions, dimensions: [
        float(np.mean([p[0] for p in positions])),
        float(np.mean([p[1] for p in positions])),
    ]

    lattice_analysis_utils_stub.get_inter_position_distance = (
        lambda p1, p2, dimensions: float(np.linalg.norm(np.array(p1) - np.array(p2)))
    )
    lattice_analysis_utils_stub.get_polymeric_properties = lambda positions, dimensions: [1.0, 2.0]

    monkeypatch.setitem(sys.modules, "pimms.CONFIG", config_stub)
    monkeypatch.setitem(sys.modules, "pimms.analysis_structures", analysis_structures)
    monkeypatch.setitem(sys.modules, "pimms.latticeExceptions", latticeExceptions)
    monkeypatch.setitem(sys.modules, "pimms.lattice_utils", lattice_utils_stub)
    monkeypatch.setitem(sys.modules, "pimms.lattice_analysis_utils", lattice_analysis_utils_stub)

    chain_path = pathlib.Path(__file__).resolve().parents[1] / "chain.py"
    spec = importlib.util.spec_from_file_location("pimms.chain", chain_path)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, "pimms.chain", module)
    spec.loader.exec_module(module)

    return module


@pytest.fixture
def chain_module(monkeypatch):
    return _load_chain_module_with_stubs(monkeypatch)


@pytest.fixture
def base_chain(chain_module):
    return chain_module.Chain(
        lattice_grid=np.zeros((10, 10), dtype=np.int32),
        dimensions=[10, 10],
        sequence="ABCD",
        int_seq=[1, 2, 3, 4],
        LR_int_seq=[10, 20, 30, 40],
        LR_IDX=[1, 3],
        chainID=7,
        chainType=2,
        chain_positions=[[0, 0], [1, 0], [2, 0], [3, 0]],
    )


def test_init_with_positions_sets_basic_attributes(chain_module):
    chain = chain_module.Chain(
        lattice_grid=np.zeros((8, 8), dtype=np.int32),
        dimensions=[8, 8],
        sequence="AAAA",
        int_seq=[5, 5, 5, 5],
        LR_int_seq=[9, 9, 9, 9],
        LR_IDX=[0, 2],
        chainID=3,
        chainType=1,
        chain_positions=[[0, 0], [0, 1], [0, 2], [0, 3]],
        fixed=True,
        rigid=True,
    )

    assert chain.chainID == 3
    assert chain.chainType == 1
    assert chain.homopolymer is True
    assert chain.fixed is True
    assert chain.rigid is True
    assert chain.positions == [[0, 0], [0, 1], [0, 2], [0, 3]]


def test_init_raises_when_positions_length_mismatch(chain_module):
    with pytest.raises(latticeExceptions.ChainInitializationException):
        chain_module.Chain(
            lattice_grid=np.zeros((8, 8), dtype=np.int32),
            dimensions=[8, 8],
            sequence="ABCD",
            int_seq=[1, 2, 3, 4],
            LR_int_seq=[1, 2, 3, 4],
            LR_IDX=[],
            chainID=1,
            chainType=1,
            chain_positions=[[0, 0], [1, 0]],
        )


def test_init_raises_on_invalid_lr_index(chain_module):
    with pytest.raises(latticeExceptions.ChainInitializationException, match="Long-range index"):
        chain_module.Chain(
            lattice_grid=np.zeros((8, 8), dtype=np.int32),
            dimensions=[8, 8],
            sequence="ABCD",
            int_seq=[1, 2, 3, 4],
            LR_int_seq=[1, 2, 3, 4],
            LR_IDX=[4],
            chainID=1,
            chainType=1,
            chain_positions=[[0, 0], [1, 0], [2, 0], [3, 0]],
        )


def test_init_uses_insert_chain_with_center_flag(chain_module, monkeypatch):
    calls = {}

    def fake_insert(chain_id, length, lattice_grid, default_start=None, hardwall=False):
        calls["chain_id"] = chain_id
        calls["length"] = length
        calls["default_start"] = default_start
        calls["hardwall"] = hardwall
        return [[2, 2], [3, 2], [4, 2]]

    monkeypatch.setattr(chain_module.lattice_utils, "insert_chain", fake_insert)

    chain = chain_module.Chain(
        lattice_grid=np.zeros((6, 6), dtype=np.int32),
        dimensions=[6, 6],
        sequence="ABC",
        int_seq=[1, 2, 3],
        LR_int_seq=[1, 2, 3],
        LR_IDX=[1],
        chainID=9,
        chainType=4,
        center=True,
        hardwall=True,
    )

    assert chain.positions == [[2, 2], [3, 2], [4, 2]]
    assert calls == {
        "chain_id": 9,
        "length": 3,
        "default_start": [3, 3],
        "hardwall": True,
    }


def test_init_wraps_insertion_failure_message(chain_module, monkeypatch):
    def fail_insert(*args, **kwargs):
        raise latticeExceptions.ChainInsertionFailure("boom")

    monkeypatch.setattr(chain_module.lattice_utils, "insert_chain", fail_insert)

    with pytest.raises(latticeExceptions.ChainInsertionFailure, match="Unable to insert chain 5"):
        chain_module.Chain(
            lattice_grid=np.zeros((6, 6), dtype=np.int32),
            dimensions=[6, 6],
            sequence="ABCDE",
            int_seq=[1, 1, 1, 1, 1],
            LR_int_seq=[1, 1, 1, 1, 1],
            LR_IDX=[],
            chainID=5,
            chainType=0,
            center=False,
        )


def test_length_getters_and_position_selectors(base_chain, chain_module, monkeypatch):
    assert len(base_chain) == 4
    assert base_chain.get_intcode_sequence() == [1, 2, 3, 4]
    assert base_chain.get_LR_positions() == [[1, 0], [3, 0]]
    assert np.array_equal(base_chain.get_LR_binary_array(), np.array([0, 1, 0, 1], dtype=np.int32))
    assert base_chain.get_positions_by_chain_index([0, 2]) == [[0, 0], [2, 0]]

    monkeypatch.setattr(chain_module.lattice_utils, "convert_chain_to_single_image", lambda positions, dimensions: [[x + 10, y + 10] for x, y in positions])
    assert base_chain.get_positions_by_chain_index_single_image_position([1, 3]) == [[11, 10], [13, 10]]


def test_ordered_and_single_image_paths(base_chain, chain_module, monkeypatch):
    monkeypatch.setattr(chain_module.lattice_utils, "do_positions_stradle_pbc_boundary", lambda positions: False)
    assert base_chain.does_chain_stradle_pbc_boundary() is False
    assert base_chain.get_single_image_positions() == base_chain.positions

    monkeypatch.setattr(chain_module.lattice_utils, "do_positions_stradle_pbc_boundary", lambda positions: True)
    monkeypatch.setattr(chain_module.lattice_utils, "convert_chain_to_single_image", lambda positions, dimensions: [[99, 99]] * len(positions))
    monkeypatch.setattr(chain_module.lattice_utils, "center_positions", lambda positions, dimensions: [[p[0] - 1, p[1] - 1] for p in positions])

    assert base_chain.does_chain_stradle_pbc_boundary() is True
    assert base_chain.get_single_image_positions() == [[99, 99]] * 4
    assert base_chain.get_ordered_positions(center_positions=True) == [[98, 98]] * 4


def test_set_ordered_positions_validation(base_chain):
    new_positions = [[9, 0], [8, 0], [7, 0], [6, 0]]
    base_chain.set_ordered_positions(new_positions)
    assert base_chain.positions == new_positions

    with pytest.raises(latticeExceptions.ChainAugmentFailure):
        base_chain.set_ordered_positions([[1, 1]])


def test_get_center_of_mass_uses_lattice_utils(base_chain, chain_module, monkeypatch):
    called = []

    def fake_com(positions, dimensions, on_lattice=True):
        called.append(on_lattice)
        if on_lattice:
            return [123, 456]
        return [123.5, 456.5]

    monkeypatch.setattr(chain_module.lattice_utils, "center_of_mass_from_positions", fake_com)
    assert base_chain.get_center_of_mass(on_lattice=True) == [123, 456]
    assert base_chain.get_center_of_mass(on_lattice=False) == [123.5, 456.5]
    assert called == [True, False]


def test_internal_scaling_instantaneous_and_updates(chain_module, monkeypatch):
    chain = chain_module.Chain(
        lattice_grid=np.zeros((20, 20), dtype=np.int32),
        dimensions=[20, 20],
        sequence="ABCDEFG",
        int_seq=list(range(7)),
        LR_int_seq=list(range(7)),
        LR_IDX=[2, 4],
        chainID=11,
        chainType=0,
        chain_positions=[[i, 0] for i in range(7)],
    )

    monkeypatch.setattr(
        chain_module.lattice_analysis_utils,
        "get_inter_position_distance",
        lambda p1, p2, dimensions: float(abs(p1[0] - p2[0])),
    )

    inst_dict = chain.analysis_get_instantaneous_internal_scaling(mode="dict")
    assert inst_dict == {1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0, 5: 5.0}

    inst_arr = chain.analysis_get_instantaneous_internal_scaling(mode="array")
    assert np.array_equal(inst_arr[0], np.array([1, 2, 3, 4, 5]))
    assert np.array_equal(inst_arr[1], np.array([1.0, 2.0, 3.0, 4.0, 5.0]))

    with pytest.raises(Exception, match="Invalid mode"):
        chain.analysis_get_instantaneous_internal_scaling(mode="bad")

    chain.analysis_update_internal_scaling()
    assert chain.analysis_get_cumulative_internal_scaling() == [1.0, 2.0, 3.0, 4.0, 5.0]
    assert chain.analysis_get_internal_scaling_squared() == [1.0, 4.0, 9.0, 16.0, 25.0]


def test_distance_map_update_and_accessors(chain_module, monkeypatch):
    chain = chain_module.Chain(
        lattice_grid=np.zeros((20, 20), dtype=np.int32),
        dimensions=[20, 20],
        sequence="ABCDE",
        int_seq=list(range(5)),
        LR_int_seq=list(range(5)),
        LR_IDX=[],
        chainID=12,
        chainType=0,
        chain_positions=[[i, 0] for i in range(5)],
    )

    monkeypatch.setattr(
        chain_module.lattice_analysis_utils,
        "get_inter_position_distance",
        lambda p1, p2, dimensions: float(abs(p1[0] - p2[0])),
    )

    dmap = chain.analysis_get_instantaneous_distance_map()
    assert dmap.shape == (5, 5)
    assert np.allclose(np.diag(dmap), 0.0)
    assert dmap[0, 4] == 4.0
    assert dmap[4, 0] == 0.0

    chain.analysis_update_distance_map()
    assert np.array_equal(chain.analysis_get_cumulative_distance_map(), dmap)


def test_end_to_end_and_residue_distance(base_chain, chain_module, monkeypatch):
    monkeypatch.setattr(
        chain_module.lattice_analysis_utils,
        "get_inter_position_distance",
        lambda p1, p2, dimensions: float((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2),
    )

    assert base_chain.analysis_get_end_to_end_distance() == 9.0
    assert base_chain.analysis_get_residue_residue_distance(1, 3) == 4.0


def test_polymeric_properties_and_warning(base_chain, chain_module, monkeypatch, capsys):
    # One call is used by analysis_get_radius_of_gyration, then two by
    # analysis_get_polymeric_properties (minimum-image and single-image).
    props = [[1.0, 2.0], [1.0, 2.0], [1.5, 2.5]]

    def fake_props(positions, dimensions):
        return props.pop(0)

    monkeypatch.setattr(chain_module.lattice_utils, "do_positions_stradle_pbc_boundary", lambda positions: True)
    monkeypatch.setattr(chain_module.lattice_utils, "convert_chain_to_single_image", lambda positions, dimensions: [[p[0] + 20, p[1]] for p in positions])
    monkeypatch.setattr(chain_module.lattice_analysis_utils, "get_polymeric_properties", fake_props)

    assert base_chain.analysis_get_radius_of_gyration() == 1.0
    assert base_chain.analysis_get_polymeric_properties() == [1.0, 2.0]

    captured = capsys.readouterr()
    assert "finite size artefacts" in captured.out


def test_fit_scaling_exponent_short_chain_returns_sentinel(base_chain):
    assert base_chain.analysis_fit_scaling_exponent() == (-1, -1)


def test_analysis_print_methods_emit_output(base_chain, capsys):
    base_chain.analysis_update_internal_scaling()
    base_chain.analysis_print_internal_scaling()
    base_chain.analysis_print_internal_scaling_squared()

    captured = capsys.readouterr()
    assert "1\t" in captured.out
