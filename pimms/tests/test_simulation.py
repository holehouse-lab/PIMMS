from datetime import datetime

import pytest

from pimms import simulation
from pimms.latticeExceptions import SimulationException


class _DummyChain:
    def __init__(self, chain_id=1):
        self.chainID = chain_id
        self.fixed = False

    def get_ordered_positions(self):
        return [[0, 0], [1, 0]]


class _DummyACC:
    def __init__(self, move_selection=2):
        self._move_selection = move_selection
        self.temperature = 100.0
        self.updated_move_logs = []

    def move_selector(self, chain_length):
        return self._move_selection

    def update_move_logs(self, selection, accepted):
        self.updated_move_logs.append((selection, accepted))

    def boltzmann_acceptance(self, old_energy, new_energy):
        return True

    def update_temperature(self, new_temperature):
        self.temperature = new_temperature


class _DummyLattice:
    def __init__(self, num_chains=2):
        self.lattice_to_angstroms = 1.0
        self._num_chains = num_chains
        self._chain = _DummyChain(chain_id=1)
        self.grid = [[0]]
        self.type_grid = [[0]]

    def get_number_of_chains(self):
        return self._num_chains

    def get_random_chain(self, frozen_chains=None):
        return self._chain

    def lattice_backupcopy(self):
        return (None, None, None)

    def lattice_restorefrombackup(self, grid, type_grid, chains):
        return None


class _DummyHamiltonian:
    def evaluate_total_energy(self, lattice_obj):
        return (0.0, 0.0, 0.0, 0.0, 0.0)


def _make_minimal_sim(move_selection=2, num_chains=2):
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.global_start_time = datetime.now()
    sim.n_steps = 1
    sim.QUENCH_RUN = False
    sim.resize_eq = False
    sim.auxillary_chain = False
    sim.reduced_printing = True
    sim.hardwall = False
    sim.frozen_chains = []
    sim.SAVE_AT_END = False
    sim.current_xtc_filename = "traj.xtc"
    sim.current_pdb_filename = "START.pdb"
    sim.equilibration = 0

    sim.ACC = _DummyACC(move_selection=move_selection)
    sim.MOVER = object()
    sim.Hamiltonian = _DummyHamiltonian()
    sim.LATTICE = _DummyLattice(num_chains=num_chains)

    sim.startup_analysis = lambda: None
    sim.simulation_IO = lambda i, old_energy: None
    sim.run_all_analysis = lambda i: None
    sim.end_of_simulation_analysis = lambda: None
    sim.ANAFUNCT_save_restart = lambda i: None

    return sim


def test_run_simulation_invalid_move_selection_raises_simulation_exception(monkeypatch):
    sim = _make_minimal_sim(move_selection=999)

    monkeypatch.setattr(simulation.lattice_utils, "start_xtc_file", lambda *args, **kwargs: None)

    with pytest.raises(SimulationException, match="Invalid option passed"):
        sim.run_simulation()


def test_run_simulation_all_chains_frozen_skips_move_selection(monkeypatch):
    sim = _make_minimal_sim(move_selection=2, num_chains=2)
    sim.n_steps = 3
    sim.frozen_chains = [1, 2]

    def fail_if_called(*args, **kwargs):
        raise AssertionError("get_random_chain should not be called when all chains are frozen")

    monkeypatch.setattr(sim.LATTICE, "get_random_chain", fail_if_called)
    monkeypatch.setattr(simulation.lattice_utils, "start_xtc_file", lambda *args, **kwargs: None)

    # Should complete without proposing any moves.
    sim.run_simulation()


def test_quench_update_heating_increases_temperature(monkeypatch):
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.QUENCH_RUN = True
    sim.QUENCH_FREQ = 1
    sim.QUENCH_START = 100.0
    sim.QUENCH_END = 110.0
    sim.QUENCH_STEPSIZE = 5.0
    sim.reduced_printing = True
    sim.TSMMC_USED = False
    sim.TSMMC_JUMP_TEMP = 110
    sim.TSMMC_INTERPOLATION_MODE = "linear"
    sim.TSMMC_STEP_MULTIPLIER = 1
    sim.TSMMC_NUMBER_OF_POINTS = 3
    sim.TSMMC_FIXED_OFFSET = 0
    sim.ACC = _DummyACC(move_selection=2)

    wrote = {"called": False, "temp": None}

    def fake_write_quench_file(step, temperature, energy):
        wrote["called"] = True
        wrote["temp"] = temperature

    monkeypatch.setattr(simulation.analysis_IO, "write_quench_file", fake_write_quench_file)

    sim.quench_update(i=1, old_energy=-42.0)

    assert sim.ACC.temperature == pytest.approx(105.0)
    assert wrote["called"] is True
    assert wrote["temp"] == pytest.approx(105.0)


def test_quench_update_rejects_non_positive_frequency():
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.QUENCH_FREQ = 0
    sim.QUENCH_START = 100.0
    sim.QUENCH_END = 90.0
    sim.QUENCH_STEPSIZE = 1.0
    sim.reduced_printing = True
    sim.TSMMC_USED = False
    sim.ACC = _DummyACC(move_selection=2)

    with pytest.raises(SimulationException, match="QUENCH_FREQ"):
        sim.quench_update(i=1, old_energy=0.0)


class _DummyTSMMCCoordinator:
    def __init__(self, complete=True, accepted=True, updated_acc=None):
        self._complete = complete
        self._accepted = accepted
        self._updated_acc = updated_acc
        self.system_move_original_energy = 12.0
        self.system_move_original_info = ("grid", "type_grid", "chains")

    def system_move_complete(self):
        return self._complete

    def accept_system_TSMMC(self, old_energy):
        return self._accepted

    def check_in_system_TSMMC(self, acc):
        if self._updated_acc is None:
            return acc
        return self._updated_acc


def test_auxillary_chain_update_complete_and_accepted():
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.reduced_printing = True
    sim.auxillary_chain = True
    sim.ACC = _DummyACC()
    sim.TSMMC_coordinator = _DummyTSMMCCoordinator(complete=True, accepted=True)
    sim.LATTICE = _DummyLattice()

    complete, accepted = sim.auxillary_chain_update(old_energy=10.0)

    assert (complete, accepted) == (True, True)
    assert sim.auxillary_chain is False


def test_auxillary_chain_update_complete_and_rejected_restores_lattice(monkeypatch):
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.reduced_printing = True
    sim.auxillary_chain = True
    sim.ACC = _DummyACC()
    sim.TSMMC_coordinator = _DummyTSMMCCoordinator(complete=True, accepted=False)
    sim.LATTICE = _DummyLattice()

    restored = {"called": False, "args": None}

    def fake_restore(*args):
        restored["called"] = True
        restored["args"] = args

    monkeypatch.setattr(sim.LATTICE, "lattice_restorefrombackup", fake_restore)

    complete, accepted = sim.auxillary_chain_update(old_energy=10.0)

    assert (complete, accepted) == (True, False)
    assert sim.auxillary_chain is False
    assert restored["called"] is True
    assert restored["args"] == ("grid", "type_grid", "chains")


def test_auxillary_chain_update_incomplete_updates_acc():
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.reduced_printing = True
    sim.auxillary_chain = True
    old_acc = _DummyACC()
    new_acc = _DummyACC()
    sim.ACC = old_acc
    sim.TSMMC_coordinator = _DummyTSMMCCoordinator(complete=False, accepted=False, updated_acc=new_acc)
    sim.LATTICE = _DummyLattice()

    complete, accepted = sim.auxillary_chain_update(old_energy=10.0)

    assert (complete, accepted) == (False, False)
    assert sim.ACC is new_acc
    assert sim.auxillary_chain is True


def test_quench_update_target_temperature_turns_off_quench(monkeypatch):
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.QUENCH_RUN = True
    sim.QUENCH_FREQ = 1
    sim.QUENCH_START = 100.0
    sim.QUENCH_END = 100.0
    sim.QUENCH_STEPSIZE = 5.0
    sim.reduced_printing = True
    sim.TSMMC_USED = False
    sim.ACC = _DummyACC()
    sim.ACC.temperature = 100.0

    wrote = {"called": False}
    monkeypatch.setattr(simulation.analysis_IO, "write_quench_file", lambda *args, **kwargs: wrote.__setitem__("called", True))

    sim.quench_update(i=1, old_energy=-1.0)

    assert sim.QUENCH_RUN is False
    assert wrote["called"] is False


def test_quench_update_refreshes_tsmmc_coordinator_when_enabled(monkeypatch):
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.QUENCH_RUN = True
    sim.QUENCH_FREQ = 1
    sim.QUENCH_START = 100.0
    sim.QUENCH_END = 90.0
    sim.QUENCH_STEPSIZE = 5.0
    sim.reduced_printing = True
    sim.TSMMC_USED = True
    sim.TSMMC_JUMP_TEMP = 85
    sim.TSMMC_INTERPOLATION_MODE = "linear"
    sim.TSMMC_STEP_MULTIPLIER = 2
    sim.TSMMC_NUMBER_OF_POINTS = 4
    sim.TSMMC_FIXED_OFFSET = 1
    sim.ACC = _DummyACC()
    sim.ACC.temperature = 100.0

    captured = {"args": None}

    class _FakeTSMMC:
        def __init__(self, *args):
            captured["args"] = args

    monkeypatch.setattr(simulation, "TSMMC", _FakeTSMMC)
    monkeypatch.setattr(simulation.analysis_IO, "write_quench_file", lambda *args, **kwargs: None)

    sim.quench_update(i=1, old_energy=-2.0)

    assert isinstance(sim.TSMMC_coordinator, _FakeTSMMC)
    assert captured["args"] == (95.0, 85, "linear", 2, 4, 1)


def _analysis_keyword_lookup(custom_module=None):
    return {
        "ANA_POL": 10,
        "ANA_INTSCAL": 20,
        "ANA_DISTMAP": 10,
        "ANA_ACCEPTANCE": 10,
        "ANA_CLUSTER": 40,
        "ANA_INTER_RESIDUE": 10,
        "ANA_END_TO_END": 10,
        "ANA_CUSTOM": 30,
        "RESTART_FREQ": 10,
        "ANALYSIS_FREQ": 10,
        "ANALYSIS_MODULE": custom_module,
        "ANA_RESIDUE_PAIRS": [(1, 2)],
    }


def test_setup_analysis_splits_default_and_non_default_frequencies():
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.LATTICE = _DummyLattice()
    sim.ANAFUNCT_polymeric_properties = lambda step: None
    sim.ANAFUNCT_internal_scaling = lambda step: None
    sim.ANAFUNCT_distance_map = lambda step: None
    sim.ANAFUNCT_acceptance = lambda step: None
    sim.ANAFUNCT_cluster_analysis = lambda step: None
    sim.ANAFUNCT_end_to_end = lambda step: None
    sim.ANAFUNCT_save_restart = lambda step: None
    sim.ANAFUNCT_custom_stubb = lambda step: None
    sim.build_R2R_distance_distribution_analysis = lambda pairs: (lambda step: (pairs, step))

    non_default, default = sim.setup_analysis(_analysis_keyword_lookup(custom_module=None))

    assert len(default) == 6
    assert len(non_default) == 3
    assert sorted(non_default.values()) == [20, 30, 40]


def test_setup_analysis_custom_module_receives_step_and_lattice():
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.LATTICE = _DummyLattice()
    sim.ANAFUNCT_polymeric_properties = lambda step: None
    sim.ANAFUNCT_internal_scaling = lambda step: None
    sim.ANAFUNCT_distance_map = lambda step: None
    sim.ANAFUNCT_acceptance = lambda step: None
    sim.ANAFUNCT_cluster_analysis = lambda step: None
    sim.ANAFUNCT_end_to_end = lambda step: None
    sim.ANAFUNCT_save_restart = lambda step: None
    sim.ANAFUNCT_custom_stubb = lambda step: "stub"
    sim.build_R2R_distance_distribution_analysis = lambda pairs: (lambda step: None)

    received = {"step": None, "lattice": None}

    def custom_analysis(step, lattice):
        received["step"] = step
        received["lattice"] = lattice
        return "ok"

    non_default, _default = sim.setup_analysis(_analysis_keyword_lookup(custom_module=custom_analysis))

    custom_fxn = [fxn for fxn, freq in non_default.items() if freq == 30][0]
    result = custom_fxn(99)

    assert result == "ok"
    assert received["step"] == 99
    assert received["lattice"] is sim.LATTICE


def test_run_all_analysis_respects_equilibration_and_frequencies():
    sim = simulation.Simulation.__new__(simulation.Simulation)
    sim.equilibration = 5
    sim.anafreq = 4

    called = []

    def non_default_fxn(step):
        called.append(("non_default", step))

    def default_fxn(step):
        called.append(("default", step))

    sim.non_default_freq_analysis = {non_default_fxn: 3}
    sim.default_freq_analysis = {default_fxn: 4}

    sim.run_all_analysis(step=4)
    assert called == []

    sim.run_all_analysis(step=6)
    assert called == [("non_default", 6)]

    sim.run_all_analysis(step=8)
    assert called[-1] == ("default", 8)
