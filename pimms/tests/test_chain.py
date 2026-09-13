import numpy as np
import pytest

import pimms.latticeExceptions as latticeExceptions


@pytest.fixture
def chain_module():
    """The real pimms.chain module.

    Historically this fixture re-exec'd chain.py against stub modules, but the
    stubs were ALWAYS bypassed: ``from . import lattice_utils`` inside chain.py
    resolves via the already-populated ``pimms`` package attributes, never via
    the patched ``sys.modules`` entries, so every test has always run against
    the real modules. The inert ~50-line stub construction has been removed;
    per-test isolation is done (as before) by monkeypatching attributes on the
    real modules, which pytest restores automatically.
    """
    import pimms.chain
    return pimms.chain


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
        rigid=True,
    )

    assert chain.chainID == 3
    assert chain.chainType == 1
    assert chain.homopolymer is True
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


@pytest.mark.parametrize("overrides, message", [
    ({"sequence": ""}, "non-empty"),
    ({"int_seq": [1]}, "one integer code per residue"),
    ({"LR_int_seq": [1]}, "one integer code per residue"),
    ({"LR_IDX": [1.5]}, "Long-range index"),
    ({"chainID": 0}, "chainID"),
    ({"chainType": -1}, "chainType"),
    ({"dimensions": [8.5, 8]}, "dimensions"),
    ({"chain_positions": [[0, 0], [1, 0], [2, 0], [8, 0]]}, "outside"),
])
def test_init_rejects_malformed_programmatic_chain_inputs(
        chain_module, overrides, message):
    kwargs = dict(
        lattice_grid=np.zeros((8, 8), dtype=np.int32),
        dimensions=[8, 8], sequence="ABCD", int_seq=[1, 2, 3, 4],
        LR_int_seq=[1, 2, 3, 4], LR_IDX=[], chainID=1, chainType=0,
        chain_positions=[[0, 0], [1, 0], [2, 0], [3, 0]],
    )
    kwargs.update(overrides)
    with pytest.raises(latticeExceptions.ChainInitializationException, match=message):
        chain_module.Chain(**kwargs)


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

    # beads sit at [0..6, 0] in a 20-wide box, so no separation reaches the half-box
    # wrap point and the minimum-image distance for a gap of g is exactly g
    inst_dict = chain.analysis_get_instantaneous_internal_scaling(mode="dict")
    assert inst_dict == {1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0, 5: 5.0, 6: 6.0}

    inst_arr = chain.analysis_get_instantaneous_internal_scaling(mode="array")
    assert np.array_equal(inst_arr[0], np.array([1, 2, 3, 4, 5, 6]))
    assert np.array_equal(inst_arr[1], np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))

    with pytest.raises(Exception, match="Invalid mode"):
        chain.analysis_get_instantaneous_internal_scaling(mode="bad")

    chain.analysis_update_internal_scaling()
    assert chain.analysis_get_cumulative_internal_scaling() == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    assert chain.analysis_get_internal_scaling_squared() == [1.0, 4.0, 9.0, 16.0, 25.0, 36.0]


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

    # beads at [0..4, 0] in a 20-wide box: minimum-image distance is just |i - j|
    dmap = chain.analysis_get_instantaneous_distance_map()
    assert dmap.shape == (5, 5)
    assert np.allclose(np.diag(dmap), 0.0)
    assert dmap[0, 4] == 4.0
    assert dmap[4, 0] == 4.0          # full symmetric matrix, not just the upper triangle
    assert np.array_equal(dmap, dmap.T)

    chain.analysis_update_distance_map()
    assert np.array_equal(chain.analysis_get_cumulative_distance_map(), dmap)


def test_end_to_end_and_residue_distance(base_chain, chain_module, monkeypatch):
    monkeypatch.setattr(
        chain_module.lattice_analysis_utils,
        "get_inter_position_distance",
        lambda p1, p2, dimensions, pbc_correction=True:
            float((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2),
    )

    assert base_chain.analysis_get_end_to_end_distance() == 9.0
    assert base_chain.analysis_get_residue_residue_distance(1, 3) == 4.0


def test_polymeric_properties_computed_on_the_chain_made_whole(base_chain, chain_module, monkeypatch, capsys):
    """Under PBC the properties are computed on the chain made whole (bond-walked
    into one image) with the PBC correction OFF - never on the raw wrapped
    positions with a COM-relative image selection, which tore chains spanning
    more than half the box."""
    calls = []

    def fake_props(positions, dimensions, pbc_correction=True):
        calls.append((list(positions), pbc_correction))
        return [1.0, 2.0]

    monkeypatch.setattr(chain_module.lattice_utils, "do_positions_stradle_pbc_boundary", lambda positions: True)
    monkeypatch.setattr(chain_module.lattice_utils, "make_chain_whole",
                        lambda positions, dimensions: [[p[0] + 20, p[1]] for p in positions])
    monkeypatch.setattr(chain_module.lattice_analysis_utils, "get_polymeric_properties", fake_props)

    assert base_chain.analysis_get_radius_of_gyration() == 1.0
    assert base_chain.analysis_get_polymeric_properties() == [1.0, 2.0]

    # both calls used the whole-chain positions and no PBC correction
    assert len(calls) == 2
    for positions, pbc_correction in calls:
        assert pbc_correction is False
        assert positions == [[p[0] + 20, p[1]] for p in base_chain.positions]
    assert "spans more than half the box" not in capsys.readouterr().out


def test_finite_size_warning_fires_once_when_a_chain_spans_over_half_the_box(chain_module, capsys):
    """The finite-size warning used to compare two calculations that tore a chain
    identically, so it could never fire. It now fires (once per chain) when the
    whole chain's extent exceeds half the box on any axis."""
    # a 9-bead rod along x in a 12-wide box: extent 8 > 6, no boundary crossing
    positions = [[i, 4] for i in range(9)]
    chain = chain_module.Chain(
        lattice_grid=np.zeros((12, 12), dtype=np.int32),
        dimensions=[12, 12],
        sequence="A" * len(positions),
        int_seq=[1] * len(positions),
        LR_int_seq=[1] * len(positions),
        LR_IDX=[],
        chainID=7,
        chainType=0,
        chain_positions=positions,
        hardwall=False,
    )
    rg, asph = chain.analysis_get_polymeric_properties()
    # exact Rg of 9 collinear unit-spaced points: sqrt(mean((i-4)^2)) = sqrt(60/9)
    assert rg == pytest.approx(np.sqrt(60.0 / 9.0))
    out = capsys.readouterr().out
    assert out.count("spans more than half the box") == 1
    assert "Chain 7" in out and "axis 0: extent 8 of box 12" in out
    # second call: same answer, no second warning
    assert chain.analysis_get_polymeric_properties()[0] == pytest.approx(rg)
    assert "spans more than half the box" not in capsys.readouterr().out

    # the same rod short of half the box is silent
    short = chain_module.Chain(
        lattice_grid=np.zeros((20, 20), dtype=np.int32), dimensions=[20, 20],
        sequence="A" * 9, int_seq=[1] * 9, LR_int_seq=[1] * 9, LR_IDX=[], chainID=8,
        chainType=0, chain_positions=positions, hardwall=False)
    short.analysis_get_polymeric_properties()
    assert "spans more than half the box" not in capsys.readouterr().out


def test_hardwall_observables_use_cartesian_not_minimum_image(chain_module):
    positions = [[i, 4] for i in range(9)]
    chain = chain_module.Chain(
        lattice_grid=np.zeros((10, 10), dtype=np.int32),
        dimensions=[10, 10],
        sequence="A" * len(positions),
        int_seq=[1] * len(positions),
        LR_int_seq=[1] * len(positions),
        LR_IDX=[],
        chainID=1,
        chainType=0,
        chain_positions=positions,
        hardwall=True,
    )

    arr = np.asarray(positions, dtype=float)
    delta = arr - arr.mean(axis=0)
    reference_rg = np.sqrt(np.trace((delta.T @ delta) / len(arr)))

    assert chain.hardwall is True
    assert chain.analysis_get_end_to_end_distance() == pytest.approx(8.0)
    assert chain.analysis_get_residue_residue_distance(0, 8) == pytest.approx(8.0)
    assert chain.analysis_get_instantaneous_distance_map()[0, 8] == pytest.approx(8.0)
    assert chain.analysis_get_radius_of_gyration() == pytest.approx(reference_rg)
    assert chain.analysis_get_polymeric_properties()[0] == pytest.approx(reference_rg)


def test_lr_binary_array_is_cached_and_read_only(base_chain):
    """LR_IDX never changes, so the per-bead flag array is built once in the constructor.

    It used to be rebuilt on every call with `if i in self.LR_IDX` against a list -
    O(L^2) per call, on a function called once per chain in every full energy evaluation
    and on every single-chain move.
    """
    first = base_chain.get_LR_binary_array()
    second = base_chain.get_LR_binary_array()

    assert list(first) == [0, 1, 0, 1]              # LR_IDX = [1, 3]
    assert second is first                          # same object, not rebuilt
    assert not first.flags.writeable                # callers cannot corrupt the shared copy


def test_fit_scaling_exponent_short_chain_returns_sentinel(base_chain):
    assert base_chain.analysis_fit_scaling_exponent() == (-1, -1)


def test_analysis_print_methods_emit_output(base_chain, capsys):
    base_chain.analysis_update_internal_scaling()
    base_chain.analysis_print_internal_scaling()
    base_chain.analysis_print_internal_scaling_squared()

    captured = capsys.readouterr()
    assert "1\t" in captured.out
