import numpy as np
import pytest

from pimms import lattice
from pimms.latticeExceptions import LatticeInitializationException, RestartException


class _DummyChain:
    def __init__(self, chain_id, positions):
        self.chainID = chain_id
        self._positions = positions
        self.int_sequence = [1] * len(positions)

    def get_ordered_positions(self):
        return self._positions


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


def test_get_random_chain_raises_when_all_chains_frozen():
    lat = _make_uninitialized_lattice([4, 4])
    lat.chains = {1: "chain1", 2: "chain2"}

    with pytest.raises(LatticeInitializationException, match="all chains are frozen"):
        lat.get_random_chain(frozen_chains=[1, 2])
