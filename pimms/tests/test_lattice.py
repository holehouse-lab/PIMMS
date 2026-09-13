import numpy as np
import pytest

from pimms import lattice
from pimms.latticeExceptions import LatticeInitializationException, RestartException


class _DummyChain:
    def __init__(self, chain_id, positions):
        self.chainID = chain_id
        self._positions = positions
        self.positions = positions
        self.int_sequence = [1] * len(positions)

    def get_ordered_positions(self):
        return self._positions

    def set_ordered_positions(self, positions):
        self._positions = positions
        self.positions = positions

    def __len__(self):
        return len(self._positions)


class _DummyHamiltonian:
    def convert_sequence_to_integer_sequence(self, sequence):
        return [1] * len(sequence)

    def convert_sequence_to_LR_integer_sequence(self, sequence):
        return []

    def get_indices_of_long_range_residues(self, sequence):
        return []


class _DummyRestart:
    def __init__(self, dimensions, chains=None, extra_chains=None):
        self.dimensions = dimensions
        self.chains = {} if chains is None else chains
        self.extra_chains = {} if extra_chains is None else extra_chains


def _make_uninitialized_lattice(dimensions):
    obj = lattice.Lattice.__new__(lattice.Lattice)
    obj.dimensions = dimensions
    return obj


def test_fully_defined_initialization_raises_for_lattice_grid_dimension_mismatch():
    lat = _make_uninitialized_lattice([4, 4])

    with pytest.raises(LatticeInitializationException):
        lat._Lattice__fully_defined_initialization(
            dimensions=[4, 4],
            chain_list=[],
            Hamiltonian=None,
            chainsDict={},
            lattice_grid=np.zeros((4, 4, 4), dtype=np.int32),
            type_grid=np.zeros((4, 4), dtype=np.int32),
        )


@pytest.mark.parametrize("dimensions, spacing, hardwall, message", [
    ([4], 3.65, False, "dimensions"),
    ([4.5, 4], 3.65, False, "dimensions"),
    ([4, 4], 0, False, "lattice_to_angstroms"),
    ([4, 4], float("nan"), False, "lattice_to_angstroms"),
    ([4, 4], 3.65, "False", "hardwall"),
])
def test_lattice_constructor_rejects_invalid_geometry(
        dimensions, spacing, hardwall, message):
    with pytest.raises(LatticeInitializationException, match=message):
        lattice.Lattice(
            dimensions, [], Hamiltonian=None, lattice_to_angstroms=spacing,
            chainsDict={}, lattice_grid=np.zeros((4, 4), dtype=np.int32),
            type_grid=np.zeros((4, 4), dtype=np.int32), hardwall=hardwall,
        )


def test_fully_defined_initialization_raises_for_type_grid_dimension_mismatch():
    lat = _make_uninitialized_lattice([4, 4])

    with pytest.raises(LatticeInitializationException):
        lat._Lattice__fully_defined_initialization(
            dimensions=[4, 4],
            chain_list=[],
            Hamiltonian=None,
            chainsDict={},
            lattice_grid=np.zeros((4, 4), dtype=np.int32),
            type_grid=np.zeros((4, 4, 4), dtype=np.int32),
        )


def test_fully_defined_initialization_does_not_hit_undefined_debug_symbol():
    lat = _make_uninitialized_lattice([4, 4])

    lat._Lattice__fully_defined_initialization(
        dimensions=[4, 4],
        chain_list=[],
        Hamiltonian=None,
        chainsDict={},
        lattice_grid=np.zeros((4, 4), dtype=np.int32),
        type_grid=np.zeros((4, 4), dtype=np.int32),
    )

    assert lat.grid.shape == (4, 4)
    assert lat.type_grid.shape == (4, 4)


def test_fully_defined_initialization_rejects_ghost_grid_occupancy():
    lat = _make_uninitialized_lattice([4, 4])
    grid = np.zeros((4, 4), dtype=np.int32)
    grid[1, 1] = 99

    with pytest.raises(LatticeInitializationException, match="occupancy"):
        lat._Lattice__fully_defined_initialization(
            dimensions=[4, 4], chain_list=[], Hamiltonian=None,
            chainsDict={}, lattice_grid=grid,
            type_grid=np.zeros((4, 4), dtype=np.int32),
        )


def test_initialization_from_restart_raises_for_dimension_count_mismatch():
    lat = _make_uninitialized_lattice([4, 4])
    restart = _DummyRestart(dimensions=[4, 4, 4])

    with pytest.raises(RestartException):
        lat._Lattice__initialization_from_restart(_DummyHamiltonian(), restart, hardwall=False)


def test_initialization_from_restart_raises_when_restart_dimensions_exceed_target():
    lat = _make_uninitialized_lattice([4, 4])
    restart = _DummyRestart(dimensions=[5, 4])

    with pytest.raises(RestartException):
        lat._Lattice__initialization_from_restart(_DummyHamiltonian(), restart, hardwall=False)


def test_initialization_from_restart_raises_on_duplicate_extra_chain_ids(monkeypatch):
    lat = _make_uninitialized_lattice([5, 5])

    # restart.chains and restart.extra_chains use the same chainID = 1.
    restart = _DummyRestart(
        dimensions=[5, 5],
        chains={1: [[[1, 1]], ["A"], 0]},
        extra_chains={1: [None, ["A"], 0]},
    )

    monkeypatch.setattr(lattice, "Chain", lambda *args, **kwargs: _DummyChain(args[6], [[1, 1]]))
    monkeypatch.setattr(lattice.lattice_utils, "place_chain_by_position", lambda *args, **kwargs: None)

    with pytest.raises(RestartException):
        lat._Lattice__initialization_from_restart(_DummyHamiltonian(), restart, hardwall=False)


def test_restart_chain_inherits_hardwall_boundary_mode(monkeypatch):
    lat = _make_uninitialized_lattice([5, 5])
    restart = _DummyRestart(
        dimensions=[5, 5],
        chains={1: [[[1, 1]], ["A"], 0]},
    )
    seen = []

    def fake_chain(*args, **kwargs):
        seen.append(kwargs["hardwall"])
        return _DummyChain(args[6], kwargs["chain_positions"])

    monkeypatch.setattr(lattice, "Chain", fake_chain)
    lat._Lattice__initialization_from_restart(
        _DummyHamiltonian(), restart, hardwall=True)

    assert seen == [True]


def test_get_random_chain_raises_when_all_chains_frozen():
    lat = _make_uninitialized_lattice([4, 4])
    lat.chains = {1: "chain1", 2: "chain2"}

    with pytest.raises(LatticeInitializationException, match="all chains are frozen"):
        lat.get_random_chain(frozen_chains=[1, 2])


def _restorable_lattice():
    lat = _make_uninitialized_lattice([4, 4])
    chain = _DummyChain(1, [[0, 0], [1, 0]])
    lat.chains = {1: chain}
    lat.grid = np.zeros((4, 4), dtype=np.int32)
    lat.type_grid = np.zeros((4, 4), dtype=np.int32)
    lat.grid[0, 0] = lat.grid[1, 0] = 1
    lat.type_grid[0, 0] = lat.type_grid[1, 0] = 1
    return lat


def test_restore_from_backup_is_atomic_when_backup_is_inconsistent():
    lat = _restorable_lattice()
    original_grid = lat.grid
    original_type_grid = lat.type_grid
    original_positions = [position[:] for position in lat.chains[1].positions]

    bad_grid = np.zeros((4, 4), dtype=np.int32)
    bad_type_grid = np.zeros((4, 4), dtype=np.int32)
    # The chain claims to occupy [2, 0], but the replacement grid is empty.
    with pytest.raises(LatticeInitializationException, match="does not contain chain"):
        lat.lattice_restorefrombackup(
            bad_grid, bad_type_grid, {1: [[2, 0], [3, 0]]})

    assert lat.grid is original_grid
    assert lat.type_grid is original_type_grid
    assert lat.chains[1].positions == original_positions


def test_restore_from_backup_validates_exact_chain_ids_before_mutating():
    lat = _restorable_lattice()
    original_grid = lat.grid

    with pytest.raises(LatticeInitializationException, match="chain IDs"):
        lat.lattice_restorefrombackup(
            np.zeros((4, 4), dtype=np.int32),
            np.zeros((4, 4), dtype=np.int32),
            {2: [[2, 0], [3, 0]]},
        )

    assert lat.grid is original_grid


def test_restore_from_backup_commits_consistent_state():
    lat = _restorable_lattice()
    new_grid = np.zeros((4, 4), dtype=np.int32)
    new_type_grid = np.zeros((4, 4), dtype=np.int32)
    new_grid[2, 0] = new_grid[3, 0] = 1
    new_type_grid[2, 0] = new_type_grid[3, 0] = 1

    lat.lattice_restorefrombackup(
        new_grid, new_type_grid, {1: [[2, 0], [3, 0]]})

    assert lat.grid is new_grid
    assert lat.type_grid is new_type_grid
    assert lat.chains[1].positions == [[2, 0], [3, 0]]
