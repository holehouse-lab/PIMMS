"""
Regression tests for the September 2026 final audit of 1.0.8.

Each test pins one finding: the effective keyfile keeps the custom analysis
module and flags a chain numbering a re-run would not reproduce; an empty
``-k`` argument is an error; a zero-length equilibration switches off the
equilibration offset along with the resize; a system-wide TSMMC excursion is
credited only with what its sub-moves can do; both ENERGY_CHECK failure modes
keep the trajectory and write the abort snapshot; a resized run treats a
previous run's production trajectory as stale; and the parallelization report
describes only the moves and chains that are actually in play.
"""

from __future__ import annotations

import contextlib
import os
import pathlib
import shutil
import subprocess
import sys

import mdtraj as md
import pytest

from pimms import CONFIG
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import SimulationEnergyException, SimulationException
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]

TSMMC_EXTRA: dict[str, object] = {"TSMMC_JUMP_TEMP": 100, "TSMMC_STEP_MULTIPLIER": 2,
                                  "TSMMC_NUMBER_OF_POINTS": 2}


def _write(path: pathlib.Path, moves: dict[str, float], **kw: object) -> None:
    """Write a parameter file and keyfile into ``path`` (a 3D short-range system).

    Parameters
    ----------
    path : pathlib.Path
        Run directory; created if needed.
    moves : dict
        ``{MOVE_KEYWORD: fraction}``.
    **kw
        Passed on to ``kernel_test_utils.write_keyfile``.
    """
    path.mkdir(parents=True, exist_ok=True)
    U.write_param_file(str(path / "params.prm"), "SR")
    U.write_keyfile(str(path / "KEYFILE.kf"), 3, kw.pop("hardwall", False), moves, **kw)


def _simulation(path: pathlib.Path, keyfile: str = "KEYFILE.kf") -> Simulation:
    """Parse ``keyfile`` in ``path`` and build (but do not run) its Simulation.

    Parameters
    ----------
    path : pathlib.Path
        Run directory holding the keyfile and its parameter file.
    keyfile : str, optional
        Keyfile name. Default ``KEYFILE.kf``.

    Returns
    -------
    Simulation
        The constructed simulation, built with ``path`` as the working directory.
    """
    cwd = os.getcwd()
    os.chdir(path)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            return Simulation(KeyFileParser(keyfile).keyword_lookup)
    finally:
        os.chdir(cwd)


def _run(path: pathlib.Path, keyfile: str = "KEYFILE.kf") -> Simulation:
    """Build and run the simulation in ``path``.

    Parameters
    ----------
    path : pathlib.Path
        Run directory.
    keyfile : str, optional
        Keyfile name. Default ``KEYFILE.kf``.

    Returns
    -------
    Simulation
        The finished simulation.
    """
    sim = _simulation(path, keyfile)
    cwd = os.getcwd()
    os.chdir(path)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim.run_simulation()
    finally:
        os.chdir(cwd)
    return sim


# ---------------------------------------------------------------------------
# keyfile_used.kf
# ---------------------------------------------------------------------------

def test_effective_keyfile_keeps_the_custom_analysis_module(tmp_path: pathlib.Path) -> None:
    # the parser replaces ANALYSIS_MODULE by the imported function, which cannot be
    # written back, so keyfile_used.kf used to keep ANA_CUSTOM and drop the module:
    # re-running it warned and silently skipped the analysis
    (tmp_path / "myana.py").write_text(
        "def analysis_function(step, lattice):\n"
        "    with open('custom_out.txt', 'a') as fh:\n"
        "        fh.write('%d\\n' % step)\n")
    _write(tmp_path, {"MOVE_CRANKSHAFT": 1.0}, n_steps=6, equilibration=1,
           extra={"ANALYSIS_MODULE": "myana.py", "ANA_CUSTOM": 2})
    _run(tmp_path)
    used = (tmp_path / CONFIG.EFFECTIVE_KEYFILE_NAME).read_text()
    assert any(line.split()[:3] == ["ANALYSIS_MODULE", ":", "myana.py"]
               for line in used.splitlines())

    rerun = tmp_path / "rerun"
    rerun.mkdir()
    for name in ("myana.py", "params.prm", CONFIG.EFFECTIVE_KEYFILE_NAME):
        shutil.copy(tmp_path / name, rerun / name)
    _run(rerun, CONFIG.EFFECTIVE_KEYFILE_NAME)
    assert (rerun / "custom_out.txt").read_text() == (tmp_path / "custom_out.txt").read_text()


def test_effective_keyfile_flags_chain_ids_a_rerun_would_renumber(tmp_path: pathlib.Path) -> None:
    # an EXTRA_CHAIN that joins an existing type is appended after every original
    # chain, so the chain IDs are not grouped by type, which is how the CHAIN lines
    # of keyfile_used.kf (and so a re-run) number them
    original = tmp_path / "original"
    _write(original, {"MOVE_CRANKSHAFT": 1.0}, box=[12, 12, 12], n_steps=4, equilibration=1,
           chains=[(3, "AAB"), (2, "BB")], extra={"RESTART_FREQ": 2})
    _run(original)
    restarted = tmp_path / "restarted"
    _write(restarted, {"MOVE_CRANKSHAFT": 1.0}, box=[12, 12, 12], n_steps=4, equilibration=1,
           chains=[(3, "AAB"), (2, "BB")], seed=3,
           extra={"RESTART_FILE": "restart.pimms", "EXTRA_CHAIN": "1 AAB"})
    shutil.copy(original / "restart.pimms", restarted / "restart.pimms")
    sim = _run(restarted)
    types = [sim.LATTICE.chains[cid].chainType for cid in sorted(sim.LATTICE.chains)]
    assert types != sorted(types)
    header = (restarted / CONFIG.EFFECTIVE_KEYFILE_NAME).read_text()
    assert "not grouped by chain type" in header

    # and a restart that did not reorder anything carries no such note
    plain = tmp_path / "plain"
    _write(plain, {"MOVE_CRANKSHAFT": 1.0}, box=[12, 12, 12], n_steps=4, equilibration=1,
           chains=[(3, "AAB"), (2, "BB")], seed=3, extra={"RESTART_FILE": "restart.pimms"})
    shutil.copy(original / "restart.pimms", plain / "restart.pimms")
    _run(plain)
    assert "not grouped by chain type" not in (plain / CONFIG.EFFECTIVE_KEYFILE_NAME).read_text()


# ---------------------------------------------------------------------------
# command line and parser
# ---------------------------------------------------------------------------

def test_an_empty_keyfile_argument_is_an_error(tmp_path: pathlib.Path) -> None:
    # PIMMS -k "$KF" with KF unset used to print the usage hint and exit 0, so a
    # batch job that never ran anything looked like a success
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONSAFEPATH")}
    result = subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "PIMMS"), "-k", ""],
                            cwd=str(tmp_path), env=env, capture_output=True, text=True,
                            timeout=120, check=False)
    assert result.returncode == 1
    assert "Could not open file" in result.stdout + result.stderr


def test_zero_equilibration_switches_off_the_offset_with_the_resize(tmp_path: pathlib.Path) -> None:
    # the resize was switched off with a warning, and the offset left behind then
    # failed with "RESIZED_EQUILIBRATION MUST be turned on", blaming the user for a
    # keyword they had given
    _write(tmp_path, {"MOVE_CRANKSHAFT": 1.0}, box=[20, 20, 20], hardwall=True,
           n_steps=4, equilibration=0,
           extra={"RESIZED_EQUILIBRATION": "10 10 10", "EQUILIBRATION_OFFSET": "2 2 2"})
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            keyfile = KeyFileParser("KEYFILE.kf")
    finally:
        os.chdir(cwd)
    assert keyfile.keyword_lookup["RESIZED_EQUILIBRATION"] is False
    assert keyfile.keyword_lookup["EQUILIBRATION_OFFSET"] is False


# ---------------------------------------------------------------------------
# a system-wide TSMMC excursion can do only what its sub-moves can
# ---------------------------------------------------------------------------

def test_system_tsmmc_does_not_rescue_a_move_set_that_cannot_act(tmp_path: pathlib.Path) -> None:
    # dimers cannot be pivoted, and an excursion only ever runs the pivot, so
    # nothing in this run could move - it used to run as a sequence of no-ops
    _write(tmp_path, {"MOVE_CHAIN_PIVOT": 0.5, "MOVE_SYSTEM_TSMMC": 0.5},
           chains=[(4, "AB")], n_steps=4, equilibration=1, extra=TSMMC_EXTRA)
    with pytest.raises(SimulationException, match="No enabled move can act"):
        _simulation(tmp_path)


def test_system_tsmmc_with_rigid_moves_still_warns_that_shapes_are_frozen(tmp_path: pathlib.Path) -> None:
    _write(tmp_path, {"MOVE_CHAIN_TRANSLATE": 0.5, "MOVE_SYSTEM_TSMMC": 0.5},
           chains=[(4, "AABBA")], n_steps=4, equilibration=1, extra=TSMMC_EXTRA)
    _simulation(tmp_path)
    assert "contains no move that can change a chain's shape" in (tmp_path / "log.txt").read_text()


def test_system_tsmmc_alone_is_credited_with_its_crankshaft_fallback(tmp_path: pathlib.Path) -> None:
    # with no non-TSMMC move the excursions run crankshaft megamoves, which can
    # reshape chains, so there is nothing to warn about
    _write(tmp_path, {"MOVE_SYSTEM_TSMMC": 1.0},
           chains=[(4, "AABBA")], n_steps=4, equilibration=1, extra=TSMMC_EXTRA)
    _simulation(tmp_path)
    assert "can change a chain's shape" not in (tmp_path / "log.txt").read_text()


# ---------------------------------------------------------------------------
# ENERGY_CHECK aborts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("failure", ("energy", "grid"))
def test_both_energy_check_failures_keep_the_trajectory_and_dump_the_state(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, failure: str) -> None:
    # the grid-consistency failure used to raise before the abort handling, losing
    # every SAVE_AT_END frame and writing no CONFIG_AT_ENERGY_FAIL snapshot
    _write(tmp_path, {"MOVE_CRANKSHAFT": 1.0}, n_steps=30, equilibration=0,
           extra={"SAVE_AT_END": "True", "XTC_FREQ": 1, "ENERGY_CHECK": 12})
    real_io = Simulation.simulation_IO

    def corrupt_at_step_12(self: Simulation, i: int, energy: int) -> object:
        if i == 12:
            if failure == "grid":
                # retype one occupied site: the tracked and recomputed energies both
                # read type_grid, so only the grid cross-check can see this
                site = tuple(self.LATTICE.chains[1].get_ordered_positions()[0])
                self.LATTICE.type_grid[site] = 1 if self.LATTICE.type_grid[site] != 1 else 2
            else:
                energy = energy + 5
        return real_io(self, i, energy)

    monkeypatch.setattr(Simulation, "simulation_IO", corrupt_at_step_12)
    with pytest.raises(SimulationEnergyException):
        _run(tmp_path)
    assert (tmp_path / "CONFIG_AT_ENERGY_FAIL.pdb").exists()
    assert (tmp_path / "CONFIG_AT_ENERGY_FAIL.xtc").exists()
    frames = md.load(str(tmp_path / "traj.xtc"), top=str(tmp_path / "START.pdb")).n_frames
    assert frames == 13          # frame 0 plus the 12 buffered steps


# ---------------------------------------------------------------------------
# stale outputs of a resized run
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("save_eq", (True, False))
def test_a_resized_run_treats_a_previous_production_trajectory_as_stale(
        tmp_path: pathlib.Path, save_eq: bool) -> None:
    # a resized run opens traj.xtc only at the resize, so one found at start-up is a
    # previous run's; it used to be kept, and a run that died during equilibration
    # sat beside another run's production trajectory
    _write(tmp_path, {"MOVE_CRANKSHAFT": 1.0}, box=[20, 20, 20], hardwall=True,
           n_steps=6, equilibration=3,
           extra={"RESIZED_EQUILIBRATION": "12 12 12", "SAVE_EQ": str(save_eq)})
    stale = _simulation(tmp_path).stale_output_files()
    assert "traj.xtc" in stale and "START.pdb" in stale

    # a run that is not resized has already opened its own production pair
    plain = tmp_path / "plain"
    _write(plain, {"MOVE_CRANKSHAFT": 1.0}, n_steps=6, equilibration=3)
    assert "traj.xtc" not in _simulation(plain).stale_output_files()


# ---------------------------------------------------------------------------
# parallelization report
# ---------------------------------------------------------------------------

def test_parallel_report_describes_only_what_is_in_play(tmp_path: pathlib.Path) -> None:
    # the crankshaft line used to be printed with a kernel and block grid when the
    # crankshaft was not in the move set, and the pull line counted monomers and
    # dimers, which neither pull kernel ever touches, as parallel-kernel chains
    state = U.build_state(tmp_path, 3, "SR", False,
                          {"MOVE_PULL": 0.5, "MOVE_CHAIN_TRANSLATE": 0.5},
                          box=[24, 24, 24], chains=[(10, "AAAA"), (30, "A"), (5, "AB")],
                          extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 2})
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        report = "\n".join(state.sim.report_parallelization())
    assert "Crankshaft (MOVE_CRANKSHAFT): not in the move set" in report
    assert "mega_crank_parallel" not in report
    assert "parallel kernel for 10 chain(s)" in report
    assert "(35 chain(s) shorter than 3 beads are never moved by this move)" in report
