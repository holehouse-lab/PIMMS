import numpy as np
import pytest

from pimms import energy, lattice, lattice_utils
from pimms.latticeExceptions import (
    ChainInsertionFailure,
    LatticeInitializationException,
    RestartException,
)


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
    lat = _make_uninitialized_lattice([7, 7])

    with pytest.raises(LatticeInitializationException):
        lat._Lattice__fully_defined_initialization(
            dimensions=[7, 7],
            chain_list=[],
            Hamiltonian=None,
            chainsDict={},
            lattice_grid=np.zeros((7, 7, 7), dtype=np.int32),
            type_grid=np.zeros((7, 7), dtype=np.int32),
        )


@pytest.mark.parametrize("dimensions, spacing, hardwall, message", [
    ([7], 3.65, False, "dimensions"),
    ([7.5, 7], 3.65, False, "dimensions"),
    ([7, 7], 0, False, "lattice_to_angstroms"),
    ([7, 7], float("nan"), False, "lattice_to_angstroms"),
    ([7, 7], 3.65, "False", "hardwall"),
])
def test_lattice_constructor_rejects_invalid_geometry(
        dimensions, spacing, hardwall, message):
    with pytest.raises(LatticeInitializationException, match=message):
        lattice.Lattice(
            dimensions, [], Hamiltonian=None, lattice_to_angstroms=spacing,
            chainsDict={}, lattice_grid=np.zeros((7, 7), dtype=np.int32),
            type_grid=np.zeros((7, 7), dtype=np.int32), hardwall=hardwall,
        )


@pytest.mark.parametrize("dimensions", [[4, 4], [6, 7], [7, 7, 5], [5, 5, 5], [7, 1]])
@pytest.mark.parametrize("hardwall", [False, True])
def test_lattice_constructor_rejects_axes_shorter_than_seven(dimensions, hardwall):
    # Below 7 sites the LR (Chebyshev 2) and SLR (Chebyshev 3) shells alias
    # across the periodic wrap: in a box of 5 two LR monomers two sites apart
    # were also scored as an SLR pair, (-11, 0, -10, -1, 0) instead of
    # (-10, 0, -10, 0, 0). The keyfile parser always refused such boxes; the
    # Lattice itself used to accept them.
    grid = np.zeros(dimensions, dtype=np.int32)
    with pytest.raises(LatticeInitializationException, match="at least 7"):
        lattice.Lattice(
            dimensions, [], Hamiltonian=None, lattice_to_angstroms=3.65,
            chainsDict={}, lattice_grid=grid, type_grid=grid.copy(),
            hardwall=hardwall,
        )


@pytest.mark.parametrize("dimensions", [[7, 7], [7, 9], [7, 7, 7], [8, 7, 12]])
def test_lattice_constructor_accepts_seven_site_axes(dimensions):
    grid = np.zeros(dimensions, dtype=np.int32)
    lat = lattice.Lattice(
        dimensions, [], Hamiltonian=None, lattice_to_angstroms=3.65,
        chainsDict={}, lattice_grid=grid, type_grid=grid.copy(),
    )
    assert lat.dimensions == dimensions


def test_fully_defined_initialization_raises_for_type_grid_dimension_mismatch():
    lat = _make_uninitialized_lattice([7, 7])

    with pytest.raises(LatticeInitializationException):
        lat._Lattice__fully_defined_initialization(
            dimensions=[7, 7],
            chain_list=[],
            Hamiltonian=None,
            chainsDict={},
            lattice_grid=np.zeros((7, 7), dtype=np.int32),
            type_grid=np.zeros((7, 7, 7), dtype=np.int32),
        )


def test_fully_defined_initialization_does_not_hit_undefined_debug_symbol():
    lat = _make_uninitialized_lattice([7, 7])

    lat._Lattice__fully_defined_initialization(
        dimensions=[7, 7],
        chain_list=[],
        Hamiltonian=None,
        chainsDict={},
        lattice_grid=np.zeros((7, 7), dtype=np.int32),
        type_grid=np.zeros((7, 7), dtype=np.int32),
    )

    assert lat.grid.shape == (7, 7)
    assert lat.type_grid.shape == (7, 7)


def test_fully_defined_initialization_rejects_ghost_grid_occupancy():
    lat = _make_uninitialized_lattice([7, 7])
    grid = np.zeros((7, 7), dtype=np.int32)
    grid[1, 1] = 99

    with pytest.raises(LatticeInitializationException, match="occupancy"):
        lat._Lattice__fully_defined_initialization(
            dimensions=[7, 7], chain_list=[], Hamiltonian=None,
            chainsDict={}, lattice_grid=grid,
            type_grid=np.zeros((7, 7), dtype=np.int32),
        )


def test_initialization_from_restart_raises_for_dimension_count_mismatch():
    lat = _make_uninitialized_lattice([7, 7])
    restart = _DummyRestart(dimensions=[7, 7, 7])

    with pytest.raises(RestartException):
        lat._Lattice__initialization_from_restart(_DummyHamiltonian(), restart, hardwall=False)


def test_initialization_from_restart_raises_when_restart_dimensions_exceed_target():
    lat = _make_uninitialized_lattice([7, 7])
    restart = _DummyRestart(dimensions=[8, 7])

    with pytest.raises(RestartException):
        lat._Lattice__initialization_from_restart(_DummyHamiltonian(), restart, hardwall=False)


def test_initialization_from_restart_raises_on_duplicate_extra_chain_ids(monkeypatch):
    lat = _make_uninitialized_lattice([7, 7])

    # restart.chains and restart.extra_chains use the same chainID = 1.
    restart = _DummyRestart(
        dimensions=[7, 7],
        chains={1: [[[1, 1]], ["A"], 0]},
        extra_chains={1: [None, ["A"], 0]},
    )

    monkeypatch.setattr(lattice, "Chain", lambda *args, **kwargs: _DummyChain(args[6], [[1, 1]]))
    monkeypatch.setattr(lattice.lattice_utils, "place_chain_by_position", lambda *args, **kwargs: None)

    with pytest.raises(RestartException):
        lat._Lattice__initialization_from_restart(_DummyHamiltonian(), restart, hardwall=False)


def test_restart_chain_inherits_hardwall_boundary_mode(monkeypatch):
    lat = _make_uninitialized_lattice([7, 7])
    restart = _DummyRestart(
        dimensions=[7, 7],
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
    lat = _make_uninitialized_lattice([7, 7])
    lat.chains = {1: "chain1", 2: "chain2"}

    with pytest.raises(LatticeInitializationException, match="all chains are frozen"):
        lat.get_random_chain(frozen_chains=[1, 2])


def _restorable_lattice():
    lat = _make_uninitialized_lattice([7, 7])
    chain = _DummyChain(1, [[0, 0], [1, 0]])
    lat.chains = {1: chain}
    lat.grid = np.zeros((7, 7), dtype=np.int32)
    lat.type_grid = np.zeros((7, 7), dtype=np.int32)
    lat.grid[0, 0] = lat.grid[1, 0] = 1
    lat.type_grid[0, 0] = lat.type_grid[1, 0] = 1
    return lat


def test_restore_from_backup_is_atomic_when_backup_is_inconsistent():
    lat = _restorable_lattice()
    original_grid = lat.grid
    original_type_grid = lat.type_grid
    original_positions = [position[:] for position in lat.chains[1].positions]

    bad_grid = np.zeros((7, 7), dtype=np.int32)
    bad_type_grid = np.zeros((7, 7), dtype=np.int32)
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
            np.zeros((7, 7), dtype=np.int32),
            np.zeros((7, 7), dtype=np.int32),
            {2: [[2, 0], [3, 0]]},
        )

    assert lat.grid is original_grid


def test_restore_from_backup_commits_consistent_state():
    lat = _restorable_lattice()
    new_grid = np.zeros((7, 7), dtype=np.int32)
    new_type_grid = np.zeros((7, 7), dtype=np.int32)
    new_grid[2, 0] = new_grid[3, 0] = 1
    new_type_grid[2, 0] = new_type_grid[3, 0] = 1

    lat.lattice_restorefrombackup(
        new_grid, new_type_grid, {1: [[2, 0], [3, 0]]})

    assert lat.grid is new_grid
    assert lat.type_grid is new_type_grid
    assert lat.chains[1].positions == [[2, 0], [3, 0]]


@pytest.mark.parametrize("extra", ["ghost_bead", "ghost_type"])
def test_restore_from_backup_rejects_sites_outside_every_chain(extra):
    # Each chain bead in the backup is consistent with both grids, but the grids
    # hold one more site that no chain owns. That used to be committed silently,
    # leaving an ownerless bead (or an ownerless bead *type*, which every energy
    # evaluation reads) on the lattice.
    lat = _restorable_lattice()
    original_grid = lat.grid
    original_type_grid = lat.type_grid
    original_positions = [position[:] for position in lat.chains[1].positions]

    new_grid = np.zeros((7, 7), dtype=np.int32)
    new_type_grid = np.zeros((7, 7), dtype=np.int32)
    new_grid[2, 0] = new_grid[3, 0] = 1
    new_type_grid[2, 0] = new_type_grid[3, 0] = 1
    new_type_grid[5, 5] = 1
    if extra == "ghost_bead":
        new_grid[5, 5] = 1
        match = "occupancy"
    else:
        match = "empty lattice sites"

    with pytest.raises(LatticeInitializationException, match=match):
        lat.lattice_restorefrombackup(
            new_grid, new_type_grid, {1: [[2, 0], [3, 0]]})

    assert lat.grid is original_grid
    assert lat.type_grid is original_type_grid
    assert lat.chains[1].positions == original_positions


# ---------------------------------------------------------------------------
# Refusing a system with more beads than lattice sites, before placing anything
# ---------------------------------------------------------------------------

@pytest.fixture
def homopolymer_hamiltonian(tmp_path, monkeypatch):
    """A real one-residue Hamiltonian, built in a scratch directory.

    The parameter-file parser writes a copy of the file into the working
    directory, so we build it inside tmp_path.

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytest's per-test scratch directory.
    monkeypatch : pytest.MonkeyPatch
        Used to change into tmp_path for the duration of the test.

    Returns
    -------
    callable
        ``build(num_dimensions)``, which returns an energy.Hamiltonian for a
        homopolymer of residue A in 2 or 3 dimensions.
    """
    monkeypatch.chdir(tmp_path)
    (tmp_path / "p.prm").write_text("A 0 0\nA A -1\nANGLE_PENALTY A 0 0 0\n")

    def build(num_dimensions):
        return energy.Hamiltonian("p.prm", num_dimensions, False, False,
                                  reduced_printing=True)
    return build


def _forbid_chain_growth(monkeypatch):
    """Make any attempt to grow a chain on the lattice fail the test.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Used to replace lattice_utils.insert_chain for the duration of the test.

    Returns
    -------
    None
    """
    def fail(*args, **kwargs):
        raise AssertionError("a chain was placed although the system cannot fit")
    monkeypatch.setattr(lattice_utils, "insert_chain", fail)


@pytest.mark.parametrize("dimensions, chain_list, n_beads", [
    ([7, 7], [[1, "A" * 50]], 50),              # one chain longer than the box volume
    ([7, 7], [[50, "A"]], 50),                   # one monomer too many
    ([7, 7], [[3, "AAAA"], [40, "A"]], 52),      # mixed chain types
    ([7, 7, 7], [[1, "A" * 400]], 400),          # the audit's KEYFILE case
])
def test_de_novo_refuses_more_beads_than_sites_before_placing(
        homopolymer_hamiltonian, monkeypatch, dimensions, chain_list, n_beads):
    # This used to place chains until the box ran out; a lone chain that could
    # never fit then failed in the centre-insertion branch with a message
    # asking the user to report a bug.
    ham = homopolymer_hamiltonian(len(dimensions))
    _forbid_chain_growth(monkeypatch)
    n_sites = int(np.prod(dimensions))

    with pytest.raises(ChainInsertionFailure) as excinfo:
        lattice.Lattice(dimensions, chain_list, ham, 3.6)

    message = str(excinfo.value)
    assert "%i beads" % n_beads in message
    assert "%i sites" % n_sites in message
    assert "overcrowded" in message
    assert "report" not in message.lower()


def test_de_novo_fills_a_box_exactly(homopolymer_hamiltonian):
    # one bead per site is allowed: the capacity check is ``>``, not ``>=``
    ham = homopolymer_hamiltonian(2)
    lat = lattice.Lattice([7, 7], [[49, "A"]], ham, 3.6)
    assert int(np.count_nonzero(lat.grid)) == 49
    assert int(np.count_nonzero(lat.type_grid)) == 49


def test_restart_refuses_extra_chains_that_cannot_fit_before_placing(
        homopolymer_hamiltonian, monkeypatch):
    # The restart chains fit by construction; the EXTRA_CHAINs must be checked
    # against the room left before any of them is grown.
    ham = homopolymer_hamiltonian(2)
    restart = _DummyRestart(
        dimensions=[7, 7],
        chains={1: [[[0, 0], [1, 0], [2, 0], [3, 0]], "AAAA", 0]},
        extra_chains={2: [None, "A" * 46, 1]},
    )
    _forbid_chain_growth(monkeypatch)

    with pytest.raises(ChainInsertionFailure, match="50 beads"):
        lattice.Lattice([7, 7], [], ham, 3.6, restart_object=restart)


def test_restart_places_extra_chains_that_fit(homopolymer_hamiltonian):
    ham = homopolymer_hamiltonian(2)
    restart = _DummyRestart(
        dimensions=[7, 7],
        chains={1: [[[0, 0], [1, 0], [2, 0], [3, 0]], "AAAA", 0]},
        extra_chains={2: [None, "AA", 1]},
    )

    lat = lattice.Lattice([7, 7], [], ham, 3.6, restart_object=restart)

    assert sorted(lat.chains) == [1, 2]
    assert lat.check_grid_consistency() == []
    assert int(np.count_nonzero(lat.grid)) == 6
