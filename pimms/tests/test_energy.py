import numpy as np
import pytest

from pimms import energy
from pimms.latticeExceptions import EnergyException, ParameterFileException


class _DummyChain:
    def __init__(self, positions, lr_binary, intcodes):
        self._positions = positions
        self._lr_binary = np.array(lr_binary, dtype=np.int64)
        self._intcodes = intcodes

    def get_ordered_positions(self):
        return self._positions

    def get_LR_positions(self):
        return []

    def get_LR_binary_array(self):
        return self._lr_binary

    def get_intcode_sequence(self):
        return self._intcodes


class _DummyLattice:
    def __init__(self, dims, chains):
        self.dimensions = dims
        self.chains = chains
        self.type_grid = np.zeros(tuple(dims), dtype=np.int64)


def _hamiltonian_stub():
    h = energy.Hamiltonian.__new__(energy.Hamiltonian)
    h.hardwall = 0
    h.angle_lookup = np.zeros((3, 3, 3, 3, 3), dtype=np.int64)
    h.parameter_to_int_map = {"0": 0, "A": 1, "B": 2}
    h.LR_parameter_to_int_map = {"B": 2}
    h.human_readable_LR_interaction_table = {"B": {"B": -1}}
    return h


def test_empty_hamiltonian_returns_zero_for_all_methods():
    eh = energy.EmptyHamiltonian()

    assert eh.evaluate_total_energy(None) == 0.0
    assert eh.evaluate_local_energy(None, None) == 0.0
    assert eh.evaluate_local_energy_LR(None, None) == 0.0
    assert eh.evaluate_angle_energy(None, None) == 0.0
    assert eh.convert_sequence_to_integer_sequence(["A", "B"]) == [1, 1]
    assert eh.convert_sequence_to_LR_integer_sequence(["A", "B"]) == []
    assert eh.get_indices_of_long_range_residues(["A", "B"]) == []


def test_convert_sequence_to_integer_sequence_and_unknown_residue():
    h = _hamiltonian_stub()

    assert h.convert_sequence_to_integer_sequence(["0", "A", "B"]) == [0, 1, 2]

    with pytest.raises(ParameterFileException):
        h.convert_sequence_to_integer_sequence(["A", "Z"])


def test_convert_sequence_to_lr_integer_sequence_maps_non_lr_to_minus_one():
    h = _hamiltonian_stub()

    assert h.convert_sequence_to_LR_integer_sequence(["A", "B", "A"]) == [-1, 2, -1]


def test_get_indices_of_long_range_residues_filters_by_lr_table_keys():
    h = _hamiltonian_stub()

    assert h.get_indices_of_long_range_residues(["A", "B", "A", "B"]) == [1, 3]


def test_evaluate_local_energy_short_range_dispatch_2d_and_3d(monkeypatch):
    h = _hamiltonian_stub()
    h.hardwall = 1

    lattice_2d = _DummyLattice([5, 5], {})
    lattice_3d = _DummyLattice([5, 5, 5], {})
    pairs = np.array([[0, 0, 1, 1]], dtype=np.int64)
    table = np.zeros((2, 2), dtype=np.int64)

    monkeypatch.setattr(
        energy.hyperloop,
        "evaluate_local_energy_2D_shortrange",
        lambda type_grid, p, t, hardwall: (100 + hardwall + len(p)),
    )
    monkeypatch.setattr(
        energy.hyperloop,
        "evaluate_local_energy_3D_shortrange",
        lambda type_grid, p, t, hardwall: (200 + hardwall + len(p)),
    )

    assert h._Hamiltonian__evaluate_local_energy_shortrange(lattice_2d, pairs, table) == 102
    assert h._Hamiltonian__evaluate_local_energy_shortrange(lattice_3d, pairs, table) == 202


def test_evaluate_local_energy_non_short_range_dispatch_2d_and_3d(monkeypatch):
    h = _hamiltonian_stub()
    h.hardwall = 0

    lattice_2d = _DummyLattice([5, 5], {})
    lattice_3d = _DummyLattice([5, 5, 5], {})
    pairs = np.array([[0, 0, 1, 1]], dtype=np.int64)
    table = np.zeros((2, 2), dtype=np.int64)

    monkeypatch.setattr(
        energy.hyperloop,
        "evaluate_local_energy_2D_non_shortrange",
        lambda type_grid, p, t, hardwall: (300 + hardwall + len(p)),
    )
    monkeypatch.setattr(
        energy.hyperloop,
        "evaluate_local_energy_3D_non_shortrange",
        lambda type_grid, p, t, hardwall: (400 + hardwall + len(p)),
    )

    assert h._Hamiltonian__evaluate_local_energy_non_shortrange(lattice_2d, pairs, table) == 301
    assert h._Hamiltonian__evaluate_local_energy_non_shortrange(lattice_3d, pairs, table) == 401


def test_evaluate_local_energy_returns_zero_when_no_pairs():
    h = _hamiltonian_stub()
    lattice = _DummyLattice([5, 5], {})

    assert h._Hamiltonian__evaluate_local_energy_shortrange(lattice, [], np.zeros((2, 2))) == 0
    assert h._Hamiltonian__evaluate_local_energy_non_shortrange(lattice, [], np.zeros((2, 2))) == 0


def test_evaluate_local_energy_raises_for_unsupported_dimensions():
    h = _hamiltonian_stub()
    lattice = _DummyLattice([5], {})
    pairs = np.array([[0, 1]], dtype=np.int64)

    with pytest.raises(EnergyException, match="Unsupported dimensionality"):
        h._Hamiltonian__evaluate_local_energy_shortrange(lattice, pairs, np.zeros((2, 2)))

    with pytest.raises(EnergyException, match="Unsupported dimensionality"):
        h._Hamiltonian__evaluate_local_energy_non_shortrange(lattice, pairs, np.zeros((2, 2)))


def test_evaluate_angle_energy_returns_zero_for_short_chain():
    h = _hamiltonian_stub()

    assert h.evaluate_angle_energy([[0, 0], [1, 1]], [1, 1], [10, 10]) == 0.0


def test_evaluate_angle_energy_dispatches_2d_and_3d(monkeypatch):
    h = _hamiltonian_stub()

    def fake_2d(chain_pos, intcodes, lookup, n):
        assert chain_pos.dtype == energy.NP_INT_TYPE
        assert intcodes.dtype == energy.NP_INT_TYPE
        assert n == 3
        return 12

    def fake_3d(chain_pos, intcodes, lookup, n):
        assert chain_pos.dtype == energy.NP_INT_TYPE
        assert intcodes.dtype == energy.NP_INT_TYPE
        assert n == 3
        return 34

    monkeypatch.setattr(energy.hyperloop, "evaluate_angle_energy_2D", fake_2d)
    monkeypatch.setattr(energy.hyperloop, "evaluate_angle_energy_3D", fake_3d)

    p2 = [[0, 0], [1, 1], [2, 2]]
    p3 = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]

    assert h.evaluate_angle_energy(p2, [1, 1, 1], [10, 10]) == 12
    assert h.evaluate_angle_energy(p3, [1, 1, 1], [10, 10, 10]) == 34


def test_evaluate_angle_energy_raises_for_unsupported_dimensions():
    h = _hamiltonian_stub()

    with pytest.raises(EnergyException, match="Unsupported dimensionality"):
        h.evaluate_angle_energy([[0], [1], [2]], [1, 1, 1], [10])


def test_build_angle_interactions_raises_when_angle_missing_for_residue():
    h = _hamiltonian_stub()
    h.parameter_to_int_map = {"0": 0, "A": 1}

    with pytest.raises(ParameterFileException, match="no angle energies defined"):
        h.build_angle_interactions(angle_dict={"0": [0, 0, 0]}, num_dimensions=2, angles_off=False)


def test_build_angle_interactions_handles_only_solvent_without_crashing():
    h = _hamiltonian_stub()
    h.parameter_to_int_map = {"0": 0}

    h.build_angle_interactions(angle_dict={"0": [0, 0, 0]}, num_dimensions=2, angles_off=False)
    assert h.angle_lookup.shape == (1, 3, 3, 3, 3)

    h.build_angle_interactions(angle_dict={"0": [0, 0, 0]}, num_dimensions=3, angles_off=False)
    assert h.angle_lookup.shape == (1, 3, 3, 3, 3, 3, 3)


def test_build_angle_interactions_rejects_unsupported_dimensions():
    h = _hamiltonian_stub()

    with pytest.raises(EnergyException, match="Unsupported dimensionality"):
        h.build_angle_interactions(angle_dict={"A": [1, 2, 3]}, num_dimensions=4, angles_off=False)


def test_evaluate_total_energy_aggregates_terms_and_builds_lr_binary_once(monkeypatch):
    h = _hamiltonian_stub()

    chain1 = _DummyChain(positions=[[0, 0], [1, 0]], lr_binary=[1, 0], intcodes=[1, 1])
    chain2 = _DummyChain(positions=[[2, 0]], lr_binary=[1], intcodes=[2])
    lattice = _DummyLattice([10, 10], {1: chain1, 2: chain2})

    captured = {}

    def fake_build_all_envelope_pairs(all_positions, lr_binary_array, type_grid, dimensions):
        captured["all_positions"] = all_positions
        captured["lr_binary_array"] = lr_binary_array.copy()
        return ("pairs", "lr_pairs", "slr_pairs")

    monkeypatch.setattr(energy.lattice_utils, "build_all_envelope_pairs", fake_build_all_envelope_pairs)

    h.evaluate_local_energy = lambda lattice_obj, pairs: 10
    h.evaluate_local_energy_LR = lambda lattice_obj, pairs: 20
    h.evaluate_local_energy_SLR = lambda lattice_obj, pairs: 30
    h.evaluate_angle_energy = lambda chain_positions, intcodes, dims: len(chain_positions)

    total, e_local, e_lr, e_slr, e_angle = h.evaluate_total_energy(lattice)

    assert captured["all_positions"] == [[0, 0], [1, 0], [2, 0]]
    assert np.array_equal(captured["lr_binary_array"], np.array([1, 0, 1]))

    assert (e_local, e_lr, e_slr, e_angle) == (10, 20, 30, 3)
    assert total == 63


def test_build_interaction_table_non_interacting_zeroes_all(monkeypatch):
    h = _hamiltonian_stub()
    h.residue_names = ["0", "A"]
    h.LR_residue_names = ["A"]
    h.human_readable_interaction_table = {"0": {"0": 1, "A": 2}, "A": {"0": 3, "A": 4}}
    h.human_readable_LR_interaction_table = {"A": {"A": 5}}
    h.human_readable_SLR_interaction_table = {"A": {"A": 6}}
    h.reduced_printing = True

    rit, mapping, lrrit, lr_mapping, slrrit = h.build_interaction_table(non_interacting=True)

    assert mapping == {"0": 0, "A": 1}
    assert lr_mapping == {"A": 1}
    assert np.all(rit == 0)
    assert np.all(lrrit == 0)
    assert np.all(slrrit == 0)
