import numpy as np
import pytest

from pimms.moves import MoveObject


class _DummyChain:
    def __init__(self, positions, chain_id=None):
        self.positions = positions
        self.chainID = chain_id

    def get_ordered_positions(self):
        return self.positions


class _DummyLattice:
    def __init__(self, chain_ids):
        self.dimensions = [10, 10]
        self.grid = np.zeros((10, 10), dtype=np.int32)
        self.type_grid = np.zeros((10, 10), dtype=np.int32)
        self.chains = {cid: _DummyChain([[0, 0], [1, 0]]) for cid in chain_ids}


class _DummyHamiltonian:
    def __init__(self):
        self.residue_interaction_table = None
        self.LR_residue_interaction_table = None
        self.SLR_residue_interaction_table = None
        self.angle_lookup = None


class _DummyCTSMMC:
    def __init__(self):
        self.inv_temperature_schedule = [1.0]
        self.steps_per_quench_multiplier = 1
        self.inv_target_temperature = 1.0

    def accept_TSMMC(self, new_energy, old_energy, inv_target_temperature, inv_temp):
        return False


class _StopAfterSelection(Exception):
    pass


def test_multichain_tsmmc_selection_uses_integer_upper_bound(monkeypatch):
    mover = MoveObject()
    lattice = _DummyLattice(chain_ids=[1, 2, 3, 4])
    ctsmmc = _DummyCTSMMC()

    def fake_randint(low, high):
        # Regression check: `high` must be a Python int (not numpy.float64).
        assert isinstance(high, int)
        return 1

    def fake_update_idx_to_bead_multiple_chains(lattice_obj, list_of_chains):
        raise _StopAfterSelection

    monkeypatch.setattr("pimms.moves.random.randint", fake_randint)
    monkeypatch.setattr(
        "pimms.moves.crankshaft_list_functions.update_idx_to_bead_multiple_chains",
        fake_update_idx_to_bead_multiple_chains,
    )

    with pytest.raises(_StopAfterSelection):
        mover.multichain_based_TSMMC(
            original_chainID=1,
            latticeObject=lattice,
            current_energy=0.0,
            hamiltonianObject=None,
            CTSMMC=ctsmmc,
            hardwall=False,
            frozen_chains=[],
        )


def test_multichain_tsmmc_noops_when_all_chains_frozen():
    mover = MoveObject()
    lattice = _DummyLattice(chain_ids=[10, 11])
    ctsmmc = _DummyCTSMMC()

    out = mover.multichain_based_TSMMC(
        original_chainID=10,
        latticeObject=lattice,
        current_energy=5.0,
        hamiltonianObject=None,
        CTSMMC=ctsmmc,
        hardwall=False,
        frozen_chains=[10, 11],
    )

    assert out == (lattice, 5.0, 0, False)


def test_multichain_tsmmc_reports_nonzero_total_moves(monkeypatch):
    mover = MoveObject()
    lattice = _DummyLattice(chain_ids=[1, 2])
    ctsmmc = _DummyCTSMMC()
    ctsmmc.accept_TSMMC = lambda *args, **kwargs: True
    ham = _DummyHamiltonian()

    def fake_update_idx_to_bead_multiple_chains(lattice_obj, list_of_chains):
        # 2 beads from chain 1, 2 beads from chain 2.
        return np.array(
            [
                [0, 0, 0, 0, 1, 1, 1],
                [0, 0, 0, 0, 1, 2, 1],
                [0, 0, 0, 0, 2, 1, 2],
                [0, 0, 0, 0, 2, 2, 2],
            ],
            dtype=np.int32,
        )

    def fake_mega_crank(
        grid,
        type_grid,
        idx_to_bead,
        residue_table,
        lr_table,
        slr_table,
        angle_lookup,
        new_energy,
        inv_temp,
        steps_per_temperature,
        bead_selector,
        local_seed,
        hardwall_int,
    ):
        return (new_energy, steps_per_temperature)

    monkeypatch.setattr(
        "pimms.moves.crankshaft_list_functions.update_idx_to_bead_multiple_chains",
        fake_update_idx_to_bead_multiple_chains,
    )
    monkeypatch.setattr("pimms.moves.mega_crank_2D.mega_crank_2D", fake_mega_crank)
    monkeypatch.setattr("pimms.moves.np.random.choice", lambda arr, n, replace=False: np.array([1]))
    monkeypatch.setattr("pimms.moves.random.randint", lambda low, high: 1)

    _, _, total_moves, _ = mover.multichain_based_TSMMC(
        original_chainID=1,
        latticeObject=lattice,
        current_energy=0.0,
        hamiltonianObject=ham,
        CTSMMC=ctsmmc,
        hardwall=False,
        frozen_chains=[],
    )

    # chain_length=4, steps_per_quench_multiplier=1, num_temps=1
    assert total_moves == 4


def test_cluster_translate_handles_noncontiguous_chain_ids(monkeypatch):
    mover = MoveObject()
    chain2 = _DummyChain([[1, 1], [1, 2]], chain_id=2)
    chain4 = _DummyChain([[5, 5], [5, 6]], chain_id=4)
    lattice = type("L", (), {})()
    lattice.dimensions = [10, 10]
    lattice.grid = np.zeros((10, 10), dtype=np.int32)
    lattice.chains = {2: chain2, 4: chain4}

    monkeypatch.setattr(
        "pimms.moves.lattice_utils.get_all_chains_in_connected_component",
        lambda *args, **kwargs: [2],
    )
    monkeypatch.setattr("pimms.moves.lattice_utils.delete_chain_by_position", lambda *args, **kwargs: None)
    monkeypatch.setattr("pimms.moves.lattice_utils.pbc_convert", lambda pos, dims: pos)
    monkeypatch.setattr("pimms.moves.lattice_utils.get_gridvalue", lambda *args, **kwargs: 0)
    monkeypatch.setattr("pimms.moves.lattice_utils.set_gridvalue", lambda *args, **kwargs: None)
    monkeypatch.setattr("pimms.moves.lattice_utils.place_chain_by_position", lambda *args, **kwargs: None)
    monkeypatch.setattr("pimms.moves.lattice_utils.do_positions_stradle_pbc_boundary", lambda *args, **kwargs: False)
    monkeypatch.setattr("pimms.moves.numpy_utils.randneg", lambda v: v)
    monkeypatch.setattr("pimms.moves.random.randint", lambda low, high: 1)

    me, accepted = mover.cluster_translate(
        selected_chain=chain2,
        latticeObject=lattice,
        cluster_move_threshold=None,
        cluster_size_threshold=None,
        hardwall=False,
        frozen_chains=[],
    )

    assert accepted is True
    assert me.move_type == 7


def test_cluster_rotate_handles_noncontiguous_chain_ids(monkeypatch):
    mover = MoveObject()
    chain2 = _DummyChain([[1, 1], [1, 2]], chain_id=2)
    chain4 = _DummyChain([[5, 5], [5, 6]], chain_id=4)
    lattice = type("L", (), {})()
    lattice.dimensions = [10, 10]
    lattice.grid = np.zeros((10, 10), dtype=np.int32)
    lattice.chains = {2: chain2, 4: chain4}

    monkeypatch.setattr(
        "pimms.moves.lattice_utils.get_all_chains_in_connected_component",
        lambda *args, **kwargs: [2],
    )
    monkeypatch.setattr("pimms.moves.lattice_utils.delete_chain_by_position", lambda *args, **kwargs: None)
    monkeypatch.setattr("pimms.moves.lattice_utils.place_chain_by_position", lambda *args, **kwargs: None)
    monkeypatch.setattr("pimms.moves.lattice_utils.center_of_mass_from_positions", lambda pos, dims: [1, 1])
    monkeypatch.setattr("pimms.moves.lattice_utils.rotate_positions_2D", lambda pos, deg: pos)
    monkeypatch.setattr("pimms.moves.lattice_utils.pbc_convert", lambda pos, dims: pos)
    monkeypatch.setattr("pimms.moves.lattice_utils.get_gridvalue", lambda *args, **kwargs: 0)
    monkeypatch.setattr("pimms.moves.lattice_utils.set_gridvalue", lambda *args, **kwargs: None)
    monkeypatch.setattr("pimms.moves.lattice_utils.do_positions_stradle_pbc_boundary", lambda *args, **kwargs: False)
    monkeypatch.setattr("pimms.moves.random.randint", lambda low, high: 0)

    me, accepted = mover.cluster_rotate(
        selected_chain=chain2,
        latticeObject=lattice,
        cluster_move_threshold=None,
        cluster_size_threshold=None,
        hardwall=False,
        frozen_chains=[],
    )

    assert accepted is True
    assert me.move_type == 8


