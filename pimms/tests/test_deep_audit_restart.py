"""Deep-audit fixes for restart files and the run lifecycle.

Each test here pins one finding of the deep audit (the IDs are in the test
names' comments) and fails on the code as it stood before the fix:

* a stop by SIGTERM or SIGINT is a clean stop that says where it stopped;
* ``RESTART_CONTINUE`` protects the Hamiltonian (quench settings stored in the
  checkpoint and compared; recomputed energy compared with the recorded one);
* the checkpoint is the last thing written for a step;
* a continuation is refused, with nothing touched, in a directory that still
  holds simulation output;
* a restart file is validated in full while it is read (generator states,
  empty ``CHAINS``, a missing ``TEMPERATURE``);
* the checkpoint is numpy-free, synced before it is renamed, and written through
  a temporary file no other run shares;
* a failed start says so in ``log.txt`` and leaves no ``keyfile_used.kf``;
* an input file that is also an output of the run is refused.
"""

import contextlib
import copy
import os
import pickle
import pickletools
import io
import random
import shutil
import signal
import subprocess
import sys
import threading
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pytest

from pimms import CONFIG
from pimms import lattice_utils
from pimms import restart as restart_module
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import (
    KeyFileException,
    RestartException,
    SimulationException,
)
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

_MOVES = {"MOVE_CRANKSHAFT": 0.6, "MOVE_CHAIN_TRANSLATE": 0.2, "MOVE_CHAIN_PIVOT": 0.2}
_CONTINUE = (
    "RESTART_FILE : restart.pimms",
    "RESTART_CONTINUE : True",
    "RESTART_OVERRIDE_DIMENSIONS : True",
    "RESTART_OVERRIDE_HARDWALL : True",
)

_QUENCH = (
    "QUENCH_RUN : True",
    "QUENCH_START : 60",
    "QUENCH_END : 30",
    "QUENCH_STEPSIZE : 10",
    "QUENCH_FREQ : 6",
    "QUENCH_AS_EQUILIBRATION : False",
)


def _write_keyfile(
    dirpath: str,
    n_steps: int,
    *,
    seed: Optional[int] = 7,
    temperature: float = 40,
    equilibration: int = 2,
    extra: Tuple[str, ...] = (),
    parameter_file: str = "params.prm",
) -> None:
    """Write a small keyfile and its parameter file into ``dirpath``.

    Parameters
    ----------
    dirpath : str
        Run directory (created if absent).
    n_steps : int
        ``N_STEPS``.
    seed : int or None, optional
        ``SEED``; omitted when None (a continuation must not give one).
    temperature : float, optional
        ``TEMPERATURE``.
    equilibration : int, optional
        ``EQUILIBRATION``.
    extra : tuple of str, optional
        Further keyfile lines, written verbatim.
    parameter_file : str, optional
        Name the parameter file is written under and referred to by.

    Returns
    -------
    None
        ``KEYFILE.kf`` and the parameter file are written.
    """
    os.makedirs(dirpath, exist_ok=True)
    U.write_param_file(os.path.join(dirpath, parameter_file), "SR")
    lines = [
        "DIMENSIONS : 12 12 12",
        "PARAMETER_FILE : %s" % parameter_file,
        "TEMPERATURE : %s" % temperature,
        "N_STEPS : %d" % n_steps,
        "EQUILIBRATION : %d" % equilibration,
        "HARDWALL : False",
        "CRANKSHAFT_SUBSTEPS : 100",
        "EN_FREQ : 1",
        "XTC_FREQ : 2",
        "PRINT_FREQ : 100000",
        "ENERGY_CHECK : 0",
        "CHAIN : 6 AABBA",
        "CHAIN : 4 BBBB",
        "CHAIN : 5 A",
    ]
    if seed is not None:
        lines.append("SEED : %d" % seed)
    lines.extend("%s : %s" % kv for kv in _MOVES.items())
    lines.extend(extra)
    with open(os.path.join(dirpath, "KEYFILE.kf"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


@contextlib.contextmanager
def _inside(dirpath: str):
    """Run the body with ``dirpath`` as the working directory and stdout discarded.

    Parameters
    ----------
    dirpath : str
        Directory to work in.

    Yields
    ------
    None
    """
    cwd = os.getcwd()
    os.chdir(dirpath)
    try:
        with open(os.devnull, "w") as sink, contextlib.redirect_stdout(sink):
            yield
    finally:
        os.chdir(cwd)


def _build(dirpath: str, keyfile: str = "KEYFILE.kf") -> Simulation:
    """Parse the keyfile in ``dirpath`` and construct its Simulation.

    Parameters
    ----------
    dirpath : str
        Run directory.
    keyfile : str, optional
        Keyfile name.

    Returns
    -------
    Simulation
        The constructed (not yet run) simulation.
    """
    with _inside(dirpath):
        return Simulation(KeyFileParser(keyfile).keyword_lookup)


def _run(dirpath: str) -> Simulation:
    """Build and run the simulation in ``dirpath``.

    Parameters
    ----------
    dirpath : str
        Run directory.

    Returns
    -------
    Simulation
        The finished simulation.
    """
    sim = _build(dirpath)
    with _inside(dirpath):
        sim.run_simulation()
    return sim


def _run_until(
    dirpath: str,
    monkeypatch: pytest.MonkeyPatch,
    step: int,
    action: Callable[[Simulation], None],
) -> Simulation:
    """Run the simulation in ``dirpath`` and call ``action`` once ``step`` is complete.

    Parameters
    ----------
    dirpath : str
        Run directory.
    monkeypatch : pytest.MonkeyPatch
        Used to wrap ``Simulation.run_all_analysis``, which the master loop
        calls at the end of every completed step.
    step : int
        The step after whose output ``action`` is called.
    action : callable
        Called with the Simulation; typically raises or sends a signal.

    Returns
    -------
    Simulation
        The simulation (only reached if ``action`` did not end the run).
    """
    real = Simulation.run_all_analysis

    def wrapped(self: Simulation, current: int) -> None:
        real(self, current)
        if current == step:
            action(self)

    monkeypatch.setattr(Simulation, "run_all_analysis", wrapped)
    return _run(dirpath)


def _log(dirpath: str) -> str:
    """Return the text of ``log.txt`` in ``dirpath``."""
    with open(os.path.join(dirpath, CONFIG.OUTNAME_LOGFILE)) as fh:
        return fh.read()


def _raw(path: str) -> dict:
    """Return the unpickled dictionary a restart file holds."""
    with open(path, "rb") as fh:
        return pickle.load(fh)


def _rewrite(path: str, change: Callable[[dict], None]) -> None:
    """Apply ``change`` to the dictionary in a restart file and write it back.

    Parameters
    ----------
    path : str
        Restart file to edit in place.
    change : callable
        Mutates the unpickled dictionary.

    Returns
    -------
    None
    """
    data = _raw(path)
    change(data)
    with open(path, "wb") as fh:
        pickle.dump(data, fh)


def _snapshot(dirpath: str) -> Dict[str, bytes]:
    """Return ``{name: bytes}`` for every file in ``dirpath``."""
    out = {}
    for name in sorted(os.listdir(dirpath)):
        with open(os.path.join(dirpath, name), "rb") as fh:
            out[name] = fh.read()
    return out


def _steps(dirpath: str, name: str) -> List[int]:
    """Return the leading step column of a per-step output file (empty if absent)."""
    path = os.path.join(dirpath, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [int(float(line.split()[0])) for line in fh if line.split()]


def _continuation(
    root: str, source: str, name: str = "resumed", n_steps: int = 30, **keywords
) -> str:
    """Make a fresh directory holding a continuation keyfile and ``source``'s checkpoint.

    Parameters
    ----------
    root : str
        Parent directory.
    source : str
        Directory of the run whose ``restart.pimms`` is copied.
    name : str, optional
        Name of the new directory.
    n_steps : int, optional
        Total run length for the continuation.
    **keywords
        Passed to ``_write_keyfile`` (``extra`` is extended with the
        continuation keywords).

    Returns
    -------
    str
        The new directory.
    """
    d = os.path.join(root, name)
    keywords["extra"] = tuple(keywords.get("extra", ())) + _CONTINUE
    _write_keyfile(d, n_steps, seed=None, **keywords)
    shutil.copy(os.path.join(source, "restart.pimms"), os.path.join(d, "restart.pimms"))
    return d


def _t_norm_parameters(dirpath: str) -> None:
    """Make the parameter file's angle penalties temperature-normalised.

    With ``ANGLE_PENALTY_T_NORM`` the penalties are scaled by the equilibrium
    temperature (``QUENCH_END`` in a quench), which is what makes a changed
    quench setting a changed Hamiltonian.

    Parameters
    ----------
    dirpath : str
        Directory holding ``params.prm``.

    Returns
    -------
    None
    """
    path = os.path.join(dirpath, "params.prm")
    with open(path) as fh:
        text = fh.read()
    text = text.replace(
        "ANGLE_PENALTY\tA\t30\t10\t0", "ANGLE_PENALTY_T_NORM\tA\t0.6\t0.2\t0"
    )
    text = text.replace(
        "ANGLE_PENALTY\tB\t50\t20\t0", "ANGLE_PENALTY_T_NORM\tB\t1.0\t0.4\t0"
    )
    assert text.count("T_NORM") == 2
    with open(path, "w") as fh:
        fh.write(text)


@pytest.fixture
def quench_checkpoint(tmp_path) -> str:
    """A quench run (60 -> 30, T_NORM angle penalties) stopped at its step-14 checkpoint.

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytest's per-test temporary directory.

    Returns
    -------
    str
        The directory, whose ``restart.pimms`` was written at step 14 (T = 40).
    """
    d = str(tmp_path / "quench_source")
    _write_keyfile(d, 30, temperature=60, extra=_QUENCH + ("RESTART_FREQ : 7",))
    _t_norm_parameters(d)

    class Stop(Exception):
        pass

    real = Simulation.run_all_analysis

    def stop_after_14(self: Simulation, step: int) -> None:
        real(self, step)
        if step == 14:
            raise Stop

    Simulation.run_all_analysis = stop_after_14
    try:
        with pytest.raises(Stop):
            _run(d)
    finally:
        Simulation.run_all_analysis = real
    raw = _raw(os.path.join(d, "restart.pimms"))
    # two rungs down the ramp (60 -> 50 at step 6, -> 40 at step 12)
    assert raw["STEP"] == 14 and raw["TEMPERATURE"] == 40.0
    return d


@pytest.fixture
def plain_checkpoint(tmp_path) -> str:
    """A 10-step run at T = 40 with no quench; its ``restart.pimms`` is from step 10."""
    d = str(tmp_path / "plain_source")
    _write_keyfile(d, 10)
    _run(d)
    return d


# --------------------------------------------------------------------------- #
# R1-1 / S1-4: SIGTERM and SIGINT are clean stops that say where they stopped
# --------------------------------------------------------------------------- #


def test_sigterm_is_a_clean_stop_with_a_line_saying_where(tmp_path, monkeypatch):
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("RESTART_FREQ : 5",))
    before = signal.getsignal(signal.SIGTERM)
    assert before is signal.SIG_DFL, "the test needs SIGTERM at its default action"

    def send_sigterm(sim: Simulation) -> None:
        # with the default action still in place the signal would kill the test
        # process, so fail instead of sending it
        if signal.getsignal(signal.SIGTERM) is signal.SIG_DFL:
            pytest.fail("run_simulation installed no SIGTERM handler")
        os.kill(os.getpid(), signal.SIGTERM)
        # the handler runs at the next bytecode boundary of the main thread
        for _ in range(1000):
            pass

    with pytest.raises(SystemExit) as stopped:
        _run_until(d, monkeypatch, 13, send_sigterm)
    assert stopped.value.code == 128 + signal.SIGTERM
    assert signal.getsignal(signal.SIGTERM) is before, (
        "the previous handler was not restored"
    )

    lines = [ln for ln in _log(d).splitlines() if "stopped by SIGTERM" in ln]
    assert len(lines) == 1
    # stopped after step 13; the checkpoints fall on steps 5 and 10
    assert "after step 13 of 40" in lines[0] and "is from step 10" in lines[0]
    assert _raw(os.path.join(d, "restart.pimms"))["STEP"] == 10

    # the trajectory was closed and is whole: frame 0 plus every XTC_FREQ (2)
    # multiple up to the last completed step, 13
    import mdtraj as md

    traj = md.load(os.path.join(d, "traj.xtc"), top=os.path.join(d, "START.pdb"))
    assert traj.n_frames == 1 + 13 // 2
    assert not os.path.exists(os.path.join(d, restart_module.RUN_MARKER_FILENAME))


def test_sigint_is_reported_with_the_step_and_the_checkpoint(tmp_path, monkeypatch):
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("RESTART_FREQ : 5",))

    def interrupt(sim: Simulation) -> None:
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        _run_until(d, monkeypatch, 8, interrupt)
    lines = [ln for ln in _log(d).splitlines() if "interrupted by SIGINT" in ln]
    assert len(lines) == 1
    assert "after step 8 of 40" in lines[0] and "is from step 5" in lines[0]


def test_an_interrupt_during_equilibration_says_no_checkpoint_was_written(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, equilibration=20, extra=("RESTART_FREQ : 5",))

    def interrupt(sim: Simulation) -> None:
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        _run_until(d, monkeypatch, 6, interrupt)
    assert not os.path.exists(os.path.join(d, "restart.pimms"))
    line = [ln for ln in _log(d).splitlines() if "interrupted by SIGINT" in ln][0]
    assert "after step 6 of 40" in line and "nothing to resume from" in line


def test_the_sigterm_handler_leaves_an_embedding_programs_handler_alone(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)

    def embedder(signum, frame):  # pragma: no cover - never delivered
        raise AssertionError("not expected")

    seen = []
    previous = signal.signal(signal.SIGTERM, embedder)
    try:
        _run_until(
            d, monkeypatch, 4, lambda sim: seen.append(signal.getsignal(signal.SIGTERM))
        )
        assert seen == [embedder]
        assert signal.getsignal(signal.SIGTERM) is embedder
    finally:
        signal.signal(signal.SIGTERM, previous)


def test_a_run_off_the_main_thread_installs_no_handler_and_completes(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    sim = _build(d)
    outcome = []

    def target() -> None:
        try:
            sim.run_simulation()
            outcome.append(signal.getsignal(signal.SIGTERM))
        except BaseException as error:  # noqa: BLE001 - reported through the list
            outcome.append(error)

    with _inside(d):
        worker = threading.Thread(target=target)
        worker.start()
        worker.join()
    assert outcome == [signal.SIG_DFL]
    assert _steps(d, "ENERGY.dat") == list(range(1, 7))


def test_save_at_end_writes_its_buffered_frames_when_the_run_is_interrupted(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("SAVE_AT_END : True",))

    def interrupt(sim: Simulation) -> None:
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        _run_until(d, monkeypatch, 11, interrupt)
    import mdtraj as md

    traj = md.load(os.path.join(d, "traj.xtc"), top=os.path.join(d, "START.pdb"))
    # frame 0 plus every XTC_FREQ (2) multiple up to the last completed step, 11;
    # without the write on the way out the file held frame 0 alone
    assert traj.n_frames == 1 + 11 // 2


# --------------------------------------------------------------------------- #
# R1-2: RESTART_CONTINUE protects the Hamiltonian
# --------------------------------------------------------------------------- #


def test_the_checkpoint_records_the_quench_settings_and_the_equilibrium_temperature(
    quench_checkpoint, plain_checkpoint
):
    raw = _raw(os.path.join(quench_checkpoint, "restart.pimms"))
    assert raw["EQUILIBRIUM_TEMPERATURE"] == 30.0
    assert raw["QUENCH"]["QUENCH_RUN"] is True
    assert (raw["QUENCH"]["QUENCH_START"], raw["QUENCH"]["QUENCH_END"]) == (60.0, 30.0)
    assert (
        abs(raw["QUENCH"]["QUENCH_STEPSIZE"]) == 10.0
        and raw["QUENCH"]["QUENCH_FREQ"] == 6
    )
    raw = _raw(os.path.join(plain_checkpoint, "restart.pimms"))
    assert raw["EQUILIBRIUM_TEMPERATURE"] == 40.0 and raw["QUENCH"] == {
        "QUENCH_RUN": False
    }


@pytest.mark.parametrize(
    "changed, named",
    [
        (("QUENCH_END : 30", "QUENCH_END : 20"), "QUENCH_END is 20 here and was 30"),
        (("QUENCH_FREQ : 6", "QUENCH_FREQ : 3"), "QUENCH_FREQ is 3 here and was 6"),
        (
            ("QUENCH_STEPSIZE : 10", "QUENCH_STEPSIZE : 5"),
            "QUENCH_STEPSIZE is 5 here and was 10",
        ),
        (
            ("QUENCH_START : 60", "QUENCH_START : 50"),
            "QUENCH_START is 50 here and was 60",
        ),
    ],
)
def test_a_continuation_with_different_quench_settings_is_refused(
    tmp_path, quench_checkpoint, changed, named
):
    quench = tuple(changed[1] if line == changed[0] else line for line in _QUENCH)
    d = _continuation(
        str(tmp_path), quench_checkpoint, n_steps=60, temperature=60, extra=quench
    )
    _t_norm_parameters(d)
    with pytest.raises(KeyFileException, match=named):
        _build(d)


def test_a_quench_checkpoint_is_not_resumed_as_a_plain_run_at_its_temperature(
    tmp_path, quench_checkpoint
):
    # the checkpoint was written at T = 40 on a ramp ending at 30. TEMPERATURE : 40
    # with no quench used to be accepted (and TEMPERATURE : 60 was refused with the
    # advice to "set TEMPERATURE : 40"), giving a Hamiltonian scaled at 40, not 30
    for temperature in (40, 60):
        d = _continuation(
            str(tmp_path),
            quench_checkpoint,
            name="plain_%d" % temperature,
            temperature=temperature,
        )
        _t_norm_parameters(d)
        with pytest.raises(KeyFileException) as refused:
            _build(d)
        message = str(refused.value)
        assert "written by a quench run" in message and "QUENCH_END 30" in message
        assert "set TEMPERATURE :" not in message


def test_a_plain_checkpoint_is_not_resumed_as_a_quench(tmp_path, plain_checkpoint):
    quench = (
        "QUENCH_RUN : True",
        "QUENCH_START : 50",
        "QUENCH_END : 30",
        "QUENCH_STEPSIZE : 10",
        "QUENCH_FREQ : 6",
        "QUENCH_AS_EQUILIBRATION : False",
    )
    d = _continuation(str(tmp_path), plain_checkpoint, temperature=50, extra=quench)
    with pytest.raises(KeyFileException, match="was not a quench"):
        _build(d)


def test_a_continuation_under_a_different_energy_function_is_refused(
    tmp_path, plain_checkpoint
):
    # ANGLES_OFF changes the Hamiltonian and nothing the parser can see. The
    # recorded energy includes the angle penalties, so it must differ from the
    # energy recomputed without them.
    recorded = _raw(os.path.join(plain_checkpoint, "restart.pimms"))["ENERGY"]
    d = _continuation(str(tmp_path), plain_checkpoint, extra=("ANGLES_OFF : True",))
    with pytest.raises(
        SimulationException, match="recorded ENERGY %s" % recorded
    ) as refused:
        _build(d)
    assert "PARAMETER_FILE" in str(refused.value) and "ANGLES_OFF" in str(refused.value)
    # refused before the run wrote or removed any simulation output
    assert sorted(os.listdir(d)) == sorted(
        [
            "KEYFILE.kf",
            "params.prm",
            "restart.pimms",
            "log.txt",
            CONFIG.OUTPUT_USED_PARAMETER_FILE,
        ]
    )
    assert "START FAILED" in _log(d)


def test_a_legitimate_quench_continuation_passes_both_checks(
    tmp_path, quench_checkpoint
):
    d = _continuation(
        str(tmp_path), quench_checkpoint, n_steps=30, temperature=60, extra=_QUENCH
    )
    _t_norm_parameters(d)
    sim = _run(d)
    assert sim.continue_from_step == 14 and _steps(d, "ENERGY.dat") == list(
        range(15, 31)
    )
    # the ramp reached QUENCH_END at step 18, after which the run no longer steps
    # the temperature; its later checkpoints are still checkpoints of a quench run
    # and must be resumable with the same keyfile
    assert sim.ACC.temperature == 30.0
    final = _raw(os.path.join(d, "restart.pimms"))
    assert final["STEP"] == 30 and final["QUENCH"]["QUENCH_RUN"] is True
    again = _continuation(
        str(tmp_path),
        d,
        name="resumed_again",
        n_steps=36,
        temperature=60,
        extra=_QUENCH,
    )
    _t_norm_parameters(again)
    assert _run(again).continue_from_step == 30


def test_a_checkpoint_from_before_the_settings_were_stored_keeps_the_old_rule(
    tmp_path, quench_checkpoint
):
    def strip(data: dict) -> None:
        del data["QUENCH"], data["EQUILIBRIUM_TEMPERATURE"]

    # same settings: accepted, as before
    d = _continuation(
        str(tmp_path), quench_checkpoint, name="same", temperature=60, extra=_QUENCH
    )
    _t_norm_parameters(d)
    _rewrite(os.path.join(d, "restart.pimms"), strip)
    assert _build(d).continue_from_step == 14

    # a ramp that does not contain the checkpoint temperature: refused, and the
    # message says why nothing more could be checked
    quench = (
        "QUENCH_RUN : True",
        "QUENCH_START : 35",
        "QUENCH_END : 15",
        "QUENCH_STEPSIZE : 10",
        "QUENCH_FREQ : 6",
        "QUENCH_AS_EQUILIBRATION : False",
    )
    d = _continuation(
        str(tmp_path), quench_checkpoint, name="off_ramp", temperature=35, extra=quench
    )
    _t_norm_parameters(d)
    _rewrite(os.path.join(d, "restart.pimms"), strip)
    with pytest.raises(
        KeyFileException, match="before the quench settings were stored"
    ):
        _build(d)


# --------------------------------------------------------------------------- #
# R1-3: the checkpoint is the last thing written for a step
# --------------------------------------------------------------------------- #


def test_the_checkpoint_is_written_after_every_other_analysis_of_its_step(
    tmp_path, monkeypatch
):
    # ANALYSIS_FREQ 5 puts the polymer analysis (RG.dat) in the default-frequency
    # group, which used to run AFTER the checkpoint: a run killed in between had
    # a checkpoint at step 10 and no step-10 row in either segment
    d = str(tmp_path / "run")
    _write_keyfile(d, 20, extra=("ANALYSIS_FREQ : 5", "RESTART_FREQ : 10"))
    seen = {}
    real = restart_module.RestartObject.write_to_file

    def record(self) -> None:
        seen.setdefault(self.step, (_steps(".", "RG.dat"), _steps(".", "ENERGY.dat")))
        real(self)

    monkeypatch.setattr(restart_module.RestartObject, "write_to_file", record)
    _run(d)
    rg_rows, energy_rows = seen[10]
    assert rg_rows == [5, 10]
    assert energy_rows == list(range(1, 11))


def test_no_analysis_draws_from_the_global_generators(tmp_path, monkeypatch):
    # moving the checkpoint to the end of the step moves the point at which the
    # generator states are captured past the analyses. That is only harmless if
    # the analyses draw nothing, so pin it: the states stored in each checkpoint
    # are the states the generators had before the step's analyses began.
    d = str(tmp_path / "run")
    _write_keyfile(
        d, 20, extra=("ANALYSIS_FREQ : 5", "RESTART_FREQ : 5", "ANA_CLUSTER : 5")
    )
    before = {}
    real_analysis = Simulation.run_all_analysis

    def wrapped(self: Simulation, step: int) -> None:
        before[step] = (random.getstate(), np.random.get_state())
        real_analysis(self, step)

    stored = {}
    real_write = restart_module.RestartObject.write_to_file

    def record(self) -> None:
        stored[self.step] = (self.rng_python, self.rng_numpy)
        real_write(self)

    monkeypatch.setattr(Simulation, "run_all_analysis", wrapped)
    monkeypatch.setattr(restart_module.RestartObject, "write_to_file", record)
    _run(d)
    assert sorted(s for s in stored if s in before) == [5, 10, 15, 20]
    for step in (5, 10, 15, 20):
        assert stored[step][0] == before[step][0]
        assert stored[step][1][0] == before[step][1][0]
        assert np.array_equal(stored[step][1][1], before[step][1][1])
        assert tuple(stored[step][1][2:]) == tuple(before[step][1][2:])


# --------------------------------------------------------------------------- #
# R1-4: a continuation is refused in a directory that still holds output
# --------------------------------------------------------------------------- #


def test_a_continuation_in_the_stopped_segments_directory_is_refused_untouched(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 30, extra=("RESTART_FREQ : 5",))

    class Stop(Exception):
        pass

    def stop(sim: Simulation) -> None:
        raise Stop

    with pytest.raises(Stop):
        _run_until(d, monkeypatch, 10, stop)
    monkeypatch.undo()
    assert _steps(d, "ENERGY.dat") == list(range(1, 11))

    # resume in place: used to delete rows 1..10 of every file, and the trajectory
    _write_keyfile(d, 30, seed=None, extra=("RESTART_FREQ : 5",) + _CONTINUE)
    before = _snapshot(d)
    with pytest.raises(
        KeyFileException, match="already holds simulation output"
    ) as refused:
        _build(d)
    assert "ENERGY.dat" in str(refused.value) and "fresh directory" in str(
        refused.value
    )
    assert _snapshot(d) == before, "the refusal changed the working directory"

    # the same continuation in a fresh directory runs
    fresh = _continuation(str(tmp_path), d, extra=("RESTART_FREQ : 5",))
    assert _run(fresh).continue_from_step == 10


# --------------------------------------------------------------------------- #
# R1-5 / R1-6 / R1-8: the restart file is validated while it is read
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "key, garbage",
    [
        ("RNG_PYTHON", (3, (1, 2, 3), None)),
        ("RNG_PYTHON", ("not", "a", "state")),
        ("RNG_NUMPY", ("MT19937", [1, 2, 3], 0, 0, 0.0)),
        ("RNG_NUMPY", ("MT19937", [-1] * 624, 624, 0, 0.0)),
        ("RNG_NUMPY", ("no such generator", [0] * 624, 624, 0, 0.0)),
        ("RNG_NUMPY", ("MT19937", "garbage")),
    ],
)
def test_a_malformed_generator_state_is_refused_when_the_file_is_read(
    tmp_path, plain_checkpoint, key, garbage
):
    path = str(tmp_path / "bad.pimms")
    shutil.copy(os.path.join(plain_checkpoint, "restart.pimms"), path)

    def corrupt(data: dict) -> None:
        data[key] = garbage

    _rewrite(path, corrupt)
    python_state, numpy_state = random.getstate(), np.random.get_state()
    with pytest.raises(RestartException, match=key):
        restart_module.RestartObject().build_from_file(path)
    # validated on throwaway generators: the global ones were not touched
    assert random.getstate() == python_state
    assert np.array_equal(np.random.get_state()[1], numpy_state[1])


def test_a_failed_read_leaves_the_restart_object_untouched(tmp_path, plain_checkpoint):
    path = str(tmp_path / "bad.pimms")
    shutil.copy(os.path.join(plain_checkpoint, "restart.pimms"), path)

    def corrupt(data: dict) -> None:
        data["RNG_NUMPY"] = ("MT19937", [1, 2, 3], 0, 0, 0.0)

    _rewrite(path, corrupt)
    R = restart_module.RestartObject()
    with pytest.raises(RestartException):
        R.build_from_file(path)
    assert R.chains == {} and R.dimensions == [] and R.step is None


def test_a_continuation_needs_the_checkpoint_temperature(tmp_path, plain_checkpoint):
    d = _continuation(str(tmp_path), plain_checkpoint)

    def strip(data: dict) -> None:
        del data["TEMPERATURE"]

    _rewrite(os.path.join(d, "restart.pimms"), strip)
    with pytest.raises(KeyFileException, match="no TEMPERATURE"):
        _build(d)


def test_an_empty_chains_entry_is_refused_as_such(tmp_path, plain_checkpoint):
    path = str(tmp_path / "empty.pimms")
    shutil.copy(os.path.join(plain_checkpoint, "restart.pimms"), path)

    def empty(data: dict) -> None:
        data["CHAINS"] = {}

    _rewrite(path, empty)
    with pytest.raises(RestartException, match="CHAINS is empty"):
        restart_module.RestartObject().build_from_file(path)


# --------------------------------------------------------------------------- #
# X1-1 / E2-2 / PERF7: the checkpoint does not depend on the numpy that wrote it
# --------------------------------------------------------------------------- #


def _mentions_numpy(path: str) -> bool:
    """Whether the pickle at ``path`` refers to any numpy object."""
    out = io.StringIO()
    with open(path, "rb") as fh:
        pickletools.dis(fh.read(), out=out)
    return "numpy" in out.getvalue()


def test_the_checkpoint_pickle_holds_no_numpy_object(plain_checkpoint):
    path = os.path.join(plain_checkpoint, "restart.pimms")
    assert not _mentions_numpy(path)
    raw = _raw(path)
    assert raw["RNG_NUMPY"][0] == "MT19937" and len(raw["RNG_NUMPY"][1]) == 624
    assert all(type(word) is int for word in raw["RNG_NUMPY"][1])
    assert all(
        type(c) is int
        for entry in raw["CHAINS"].values()
        for bead in entry[0]
        for c in bead
    )


def test_both_forms_of_the_numpy_state_resume_the_identical_stream(
    tmp_path, plain_checkpoint
):
    new = os.path.join(plain_checkpoint, "restart.pimms")
    old = str(tmp_path / "old_form.pimms")
    shutil.copy(new, old)

    def as_array(data: dict) -> None:
        # the form the 1.0.8 development code wrote: numpy.random.get_state() verbatim
        s = data["RNG_NUMPY"]
        data["RNG_NUMPY"] = (s[0], np.asarray(s[1], dtype=np.uint32), s[2], s[3], s[4])

    _rewrite(old, as_array)
    assert _mentions_numpy(old)

    draws = []
    for path in (new, old):
        R = restart_module.RestartObject()
        R.build_from_file(path)
        assert R.rng_numpy[1].dtype == np.uint32
        generator = np.random.RandomState()
        generator.set_state(R.rng_numpy)
        draws.append(
            generator.randint(0, 2**31 - 1, 50).tolist()
            + generator.standard_normal(5).tolist()
        )
    assert draws[0] == draws[1]

    # and the stream is the one the run itself would have drawn next: the
    # checkpoint is the final one, so the run's generator is still in that state
    expected = np.random.RandomState()
    expected.set_state(
        (
            _raw(new)["RNG_NUMPY"][0],
            np.array(_raw(new)["RNG_NUMPY"][1], dtype=np.uint32),
        )
        + tuple(_raw(new)["RNG_NUMPY"][2:])
    )
    assert (
        draws[0]
        == expected.randint(0, 2**31 - 1, 50).tolist()
        + expected.standard_normal(5).tolist()
    )


class _FakeChain:
    """The three attributes of a Chain that ``build_from_lattice`` reads."""

    def __init__(self, positions: list, sequence: str, chain_type: int) -> None:
        self.positions = positions
        self.sequence = sequence
        self.chainType = chain_type


class _FakeLattice:
    """The two attributes of a Lattice that ``build_from_lattice`` reads."""

    def __init__(self, chains: dict, dimensions: list) -> None:
        self.chains = chains
        self.dimensions = dimensions


def test_the_snapshot_holds_python_ints_and_is_independent_of_the_lattice():
    # some moves leave numpy integer scalars in Chain.positions
    numpy_rows = [
        [np.int64(1), np.int64(2), np.int64(3)],
        [np.int32(1), np.int32(2), np.int32(4)],
    ]
    plain_rows = [[5, 5, 5], [5, 5, 6], [5, 6, 6]]
    lattice = _FakeLattice(
        {1: _FakeChain(numpy_rows, "AB", 0), 2: _FakeChain(plain_rows, "AAB", 1)},
        [10, 10, 10],
    )
    R = restart_module.RestartObject()
    R.build_from_lattice(lattice)
    assert R.chains[1][0] == [[1, 2, 3], [1, 2, 4]]
    assert all(
        type(c) is int for chain in R.chains.values() for bead in chain[0] for c in bead
    )
    assert b"numpy" not in pickle.dumps(R.chains)

    # where the positions were already Python ints the pickle is byte for byte
    # what the deep copy used to give
    assert pickle.dumps(R.chains[2]) == pickle.dumps(
        [copy.deepcopy(plain_rows), "AAB", 1]
    )

    # a snapshot: moving a bead on the lattice afterwards does not move it here
    plain_rows[0][0] = 9
    assert R.chains[2][0][0] == [5, 5, 5]


# --------------------------------------------------------------------------- #
# R1-7 / O1-3: how the checkpoint reaches the disk
# --------------------------------------------------------------------------- #


def test_the_checkpoint_is_synced_before_the_rename_and_uses_its_own_temporary(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 4, equilibration=0)
    sim = _run(d)
    events = []
    real_fsync, real_replace = os.fsync, os.replace

    def fsync(fd: int) -> None:
        events.append(("fsync", os.fstat(fd).st_size))
        real_fsync(fd)

    def replace(src: str, dst: str) -> None:
        events.append(("replace", src, dst, os.path.getsize(src)))
        real_replace(src, dst)

    monkeypatch.setattr(restart_module.os, "fsync", fsync)
    monkeypatch.setattr(restart_module.os, "replace", replace)
    with _inside(d):
        sim.ANAFUNCT_save_restart(4)
    assert [e[0] for e in events] == ["fsync", "replace"]
    _, src, dst, size = events[1]
    # the whole file was flushed before the sync, and the temporary is this process's own
    assert events[0][1] == size == os.path.getsize(os.path.join(d, "restart.pimms"))
    assert dst == CONFIG.RESTART_FILENAME
    assert src == "%s.tmp.%d" % (CONFIG.RESTART_FILENAME, os.getpid())
    assert not os.path.exists(os.path.join(d, src))


def _dead_pid() -> int:
    """Return the process number of a process that has already exited."""
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    child.wait()
    return child.pid


@pytest.fixture
def live_pid():
    """The process number of a process that stays alive for the whole test.

    A child that sleeps is spawned for the purpose and stopped afterwards. (The
    parent's process number is not safe to use: when pytest is the first process
    of a container it is 0.)

    Yields
    ------
    int
        The child's process number.
    """
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    try:
        yield child.pid
    finally:
        child.terminate()
        child.wait()


def _temporary(dirpath: str, suffix: str, age_seconds: float = 0.0) -> str:
    """Create a torn temporary checkpoint file and return its path.

    Parameters
    ----------
    dirpath : str
        Run directory.
    suffix : str
        Appended to ``restart.pimms.tmp`` ('' for the fixed name older versions
        used, ``'.<pid>'`` otherwise).
    age_seconds : float, optional
        How long ago the file should appear to have been written.

    Returns
    -------
    str
        The file's path.
    """
    path = os.path.join(dirpath, CONFIG.RESTART_FILENAME + ".tmp" + suffix)
    with open(path, "wb") as fh:
        fh.write(b"torn")
    if age_seconds:
        then = os.path.getmtime(path) - age_seconds
        os.utime(path, (then, then))
    return path


def test_start_up_removes_only_the_temporaries_no_running_process_can_own(
    tmp_path, live_pid
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 4, equilibration=0)
    two_hours = 2 * restart_module.STALE_TEMPORARY_AGE_SECONDS
    dead = _temporary(d, ".%d" % _dead_pid())
    live = _temporary(d, ".%d" % live_pid)
    # (a second, long-lived child: old enough that no write can still be going)
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    live_but_old = _temporary(d, ".%d" % sleeper.pid, two_hours)
    legacy_recent = _temporary(d, "")
    try:
        _run(d)
    finally:
        sleeper.terminate()
        sleeper.wait()
    # known to be gone from this machine, or too old to be a write in progress
    assert not os.path.exists(dead)
    assert not os.path.exists(live_but_old)
    # a running process's, and the fixed-name one of an older PIMMS, which names
    # no process and is recent: left alone, and the log says so
    assert os.path.exists(live), "a running process's temporary was removed"
    assert os.path.exists(legacy_recent)
    left = [ln for ln in _log(d).splitlines() if "temporary checkpoint file" in ln]
    assert len(left) == 1
    assert (
        os.path.basename(live) in left[0] and os.path.basename(legacy_recent) in left[0]
    )

    # an old fixed-name temporary goes at the next start
    then = os.path.getmtime(legacy_recent) - two_hours
    os.utime(legacy_recent, (then, then))
    _run(d)
    assert not os.path.exists(legacy_recent) and os.path.exists(live)


def test_a_temporary_whose_process_cannot_be_checked_is_left_alone(
    tmp_path, monkeypatch
):
    # where the process table cannot be asked (always the case on Windows) a
    # recent temporary used to be treated as stale and removed, which makes a run
    # that is still writing it die on its rename
    d = str(tmp_path / "run")
    os.makedirs(d)
    monkeypatch.chdir(d)
    monkeypatch.setattr(restart_module, "_process_is_running", lambda pid: None)
    recent = _temporary(d, ".%d" % (os.getpid() + 1))
    old = _temporary(
        d, ".%d" % (os.getpid() + 2), 2 * restart_module.STALE_TEMPORARY_AGE_SECONDS
    )
    stale, kept = restart_module.checkpoint_temporaries()
    assert stale == [os.path.basename(old)] and kept == [os.path.basename(recent)]
    assert restart_module.stale_checkpoint_temporaries() == stale


def test_a_process_number_is_not_trusted_when_another_machine_shares_the_directory(
    tmp_path, monkeypatch
):
    # the run marker names a run on another machine: a process number that does
    # not exist HERE may be that run's, so its recent temporary must stay
    d = str(tmp_path / "run")
    os.makedirs(d)
    monkeypatch.chdir(d)
    recent = _temporary(d, ".%d" % _dead_pid())
    assert restart_module.checkpoint_temporaries() == ([os.path.basename(recent)], [])
    with open(restart_module.RUN_MARKER_FILENAME, "w") as fh:
        fh.write("4242 some-other-node 2026-01-01 00:00:00\n")
    assert restart_module.checkpoint_temporaries() == ([], [os.path.basename(recent)])


def test_a_second_run_in_a_directory_is_warned_about_the_first(tmp_path, live_pid):
    d = str(tmp_path / "run")
    _write_keyfile(d, 4, equilibration=0)
    marker = os.path.join(d, restart_module.RUN_MARKER_FILENAME)
    import socket

    # the marker of a run that is still going (a sleeping child stands in for it)
    first = "%d %s 2026-01-01 00:00:00\n" % (live_pid, socket.gethostname())
    with open(marker, "w") as fh:
        fh.write(first)
    _run(d)
    warned = [ln for ln in _log(d).splitlines() if "another PIMMS run" in ln]
    assert len(warned) == 1 and ("process %d" % live_pid) in warned[0]
    assert "WARNING" in warned[0]
    # this run took its own line out again, and left the first run's where it was:
    # the first run is still going, and a third run must be warned about it too
    with open(marker) as fh:
        assert fh.read() == first
    _run(d)
    assert len([ln for ln in _log(d).splitlines() if "another PIMMS run" in ln]) == 1

    # the marker of a run that was killed: dropped without a warning, and the
    # file goes when the last run using the directory ends
    with open(marker, "w") as fh:
        fh.write("%d %s 2026-01-01 00:00:00\n" % (_dead_pid(), socket.gethostname()))
    _run(d)
    assert "another PIMMS run" not in _log(d)
    assert not os.path.exists(marker)


def test_the_marker_keeps_the_line_of_a_run_on_another_machine(tmp_path):
    # nothing here can tell whether a run on another machine is alive, so its
    # line is never dropped: every start is warned until the file is deleted
    d = str(tmp_path / "run")
    _write_keyfile(d, 4, equilibration=0)
    marker = os.path.join(d, restart_module.RUN_MARKER_FILENAME)
    other = "4242 some-other-node 2026-01-01 00:00:00\n"
    with open(marker, "w") as fh:
        fh.write(other)
    _run(d)
    warned = [ln for ln in _log(d).splitlines() if "another PIMMS run" in ln]
    assert len(warned) == 1 and "on some-other-node" in warned[0]
    with open(marker) as fh:
        assert fh.read() == other


def test_the_marker_is_in_place_while_the_run_is_going(tmp_path, monkeypatch):
    d = str(tmp_path / "run")
    _write_keyfile(d, 4, equilibration=0)
    seen = []
    _run_until(
        d,
        monkeypatch,
        2,
        lambda sim: seen.append(
            open(restart_module.RUN_MARKER_FILENAME).read().split()[0]
        ),
    )
    assert seen == [str(os.getpid())]


# --------------------------------------------------------------------------- #
# K1-6: a failed start says so, and leaves no keyfile_used.kf
# --------------------------------------------------------------------------- #


def test_a_start_that_fails_during_construction_is_recorded_in_the_log(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    _run(d)
    earlier = _snapshot(d)

    # a second start in the same directory that cannot be built: its freeze file
    # names a chain that does not exist
    with open(os.path.join(d, "freeze.txt"), "w") as fh:
        fh.write("C 99\n")
    _write_keyfile(d, 8, extra=("FREEZE_FILE : freeze.txt",))
    with pytest.raises(Exception):
        _build(d)
    assert "START FAILED - NO SIMULATION WAS RUN" in _log(d)
    after = _snapshot(d)
    # the earlier run's data and its effective keyfile are exactly as they were
    for name in (
        "ENERGY.dat",
        "traj.xtc",
        "START.pdb",
        "restart.pimms",
        CONFIG.EFFECTIVE_KEYFILE_NAME,
    ):
        assert after[name] == earlier[name]


def test_keyfile_used_is_not_written_by_a_start_that_fails_before_the_trajectory_opens(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)

    def cannot_open(*args, **kwargs):
        raise OSError("no trajectory for you")

    monkeypatch.setattr(lattice_utils, "open_xtc_writer", cannot_open)
    with pytest.raises(OSError, match="no trajectory for you"):
        _run(d)
    assert not os.path.exists(os.path.join(d, CONFIG.EFFECTIVE_KEYFILE_NAME))
    log = _log(d)
    assert (
        "START FAILED - NO SIMULATION WAS RUN" in log and "before the first step" in log
    )
    monkeypatch.undo()

    # a start that succeeds does write it, before the run's first step
    seen = []
    _run_until(
        d,
        monkeypatch,
        1,
        lambda sim: seen.append(os.path.exists(CONFIG.EFFECTIVE_KEYFILE_NAME)),
    )
    assert seen == [True]
    assert "START FAILED" not in _log(d)


# --------------------------------------------------------------------------- #
# R1-9: the log says when an earlier restart.pimms is still in place
# --------------------------------------------------------------------------- #


def test_the_log_says_when_an_earlier_restart_file_stays_in_place(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    _run(d)
    assert "is present in this directory from before this run started" not in _log(d)
    _run(d)
    line = [ln for ln in _log(d).splitlines() if "from before this run started" in ln]
    assert (
        len(line) == 1
        and CONFIG.RESTART_FILENAME in line[0]
        and "first checkpoint" in line[0]
    )


# --------------------------------------------------------------------------- #
# X1-10: an input file that is also an output of the run
# --------------------------------------------------------------------------- #


def test_an_input_file_named_like_an_output_is_refused_and_survives(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6, parameter_file="ENERGY.dat")
    with open(os.path.join(d, "ENERGY.dat"), "rb") as fh:
        parameters = fh.read()
    # refused by the keyfile parser (before log.txt is started) and, for a
    # keyword dictionary that did not come through a full parse, by the Simulation
    with pytest.raises(
        (KeyFileException, SimulationException), match="PARAMETER_FILE"
    ) as refused:
        _build(d)
    assert "ENERGY.dat" in str(refused.value)
    with _inside(d):
        found = Simulation._input_output_collisions({"PARAMETER_FILE": "ENERGY.dat"})
    assert found == [("PARAMETER_FILE", "ENERGY.dat", "ENERGY.dat")]
    with open(os.path.join(d, "ENERGY.dat"), "rb") as fh:
        assert fh.read() == parameters
    assert not os.path.exists(os.path.join(d, CONFIG.OUTPUT_USED_PARAMETER_FILE))


def test_an_output_name_linked_to_an_input_file_is_refused(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    os.symlink("params.prm", os.path.join(d, "START.pdb"))
    with pytest.raises((KeyFileException, SimulationException), match="START.pdb"):
        _build(d)
    with _inside(d):
        found = Simulation._input_output_collisions({"PARAMETER_FILE": "params.prm"})
    assert found == [("PARAMETER_FILE", "params.prm", "START.pdb")]
    assert os.path.islink(os.path.join(d, "START.pdb"))


def test_rerunning_keyfile_used_in_place_is_still_allowed(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    _run(d)
    first = _steps(d, "ENERGY.dat")
    with _inside(d):
        sim = Simulation(KeyFileParser(CONFIG.EFFECTIVE_KEYFILE_NAME).keyword_lookup)
        sim.run_simulation()
    assert _steps(d, "ENERGY.dat") == first


# --------------------------------------------------------------------------- #
# O1-4: start-up records are not written through a symbolic link
# --------------------------------------------------------------------------- #


def test_start_up_records_replace_a_link_rather_than_the_file_it_points_at(tmp_path):
    # a user who links an earlier segment's files into a new run directory used to
    # have the earlier keyfile_used.kf and chain_to_chainid.txt overwritten in place
    earlier = tmp_path / "earlier"
    earlier.mkdir()
    names = (CONFIG.EFFECTIVE_KEYFILE_NAME, CONFIG.OUTPUT_CHAIN_TO_CHAINID)
    for name in names:
        (earlier / name).write_text("the earlier segment's %s\n" % name)

    d = str(tmp_path / "run")
    _write_keyfile(d, 6, extra=("WRITE_CHAIN_TO_CHAINID : True",))
    for name in names:
        os.symlink(str(earlier / name), os.path.join(d, name))
    _run(d)
    for name in names:
        assert (earlier / name).read_text() == "the earlier segment's %s\n" % name
        assert not os.path.islink(os.path.join(d, name))
        with open(os.path.join(d, name)) as fh:
            assert "the earlier segment" not in fh.read()


# --------------------------------------------------------------------------- #
# O1-9: elapsed times come from the monotonic clock
# --------------------------------------------------------------------------- #


def test_elapsed_times_do_not_depend_on_the_wall_clock(tmp_path, monkeypatch, capsys):
    import datetime as datetime_module
    import time

    from pimms import analysis_general
    from pimms import simulation as simulation_module

    d = str(tmp_path / "run")
    _write_keyfile(d, 40)
    sim = _build(d)

    class ClockSetBackAnHour:
        """A wall clock that is set back one hour just after the run starts."""

        readings = 0

        @classmethod
        def now(cls) -> datetime_module.datetime:
            cls.readings += 1
            noon = datetime_module.datetime(2026, 1, 1, 12, 0, 0)
            return (
                noon if cls.readings == 1 else noon - datetime_module.timedelta(hours=1)
            )

    starts = []
    real = analysis_general.evaluate_performance

    def record(step, start_time, *args, **kwargs):
        starts.append(start_time)
        return real(step, start_time, *args, **kwargs)

    monkeypatch.setattr(simulation_module, "datetime", ClockSetBackAnHour)
    monkeypatch.setattr(analysis_general, "evaluate_performance", record)
    monkeypatch.chdir(d)
    before = time.monotonic()
    sim.run_simulation()
    after = time.monotonic()

    # the performance report is handed a monotonic reading taken as the run started
    assert starts and all(isinstance(s, float) for s in starts)
    assert all(s == sim.global_start_monotonic for s in starts)
    assert before <= sim.global_start_monotonic <= after

    # and the closing line measures the run on that clock: the difference of the
    # two wall-clock readings is -1 hour
    line = [
        ln for ln in capsys.readouterr().out.splitlines() if "Simulation time:" in ln
    ]
    assert len(line) == 1
    hours, minutes, seconds = [
        int(x) for x in line[0].split(":", 2)[-1].replace(",", " ").split()[::2]
    ]
    assert (hours, minutes) == (0, 0) and 0 <= seconds <= int(after - before) + 1


# --------------------------------------------------------------------------- #
# RF3-1: a second signal during the clean-up cannot cut it short
# --------------------------------------------------------------------------- #


def _second_signals_during_the_clean_up(monkeypatch: pytest.MonkeyPatch) -> List[str]:
    """Arrange for SIGTERM and SIGINT to be sent again while a stopped run cleans up.

    The stop report and the ``SAVE_AT_END`` rescue are each wrapped so that both
    signals are sent to this process just before the real routine runs. If a
    signal is not being ignored at that point nothing is sent (it would start a
    second stop inside the clean-up, or interrupt the test session) and the
    place is recorded instead.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Used to wrap ``Simulation._report_unfinished_run`` and
        ``lattice_utils.save_out_sim``.

    Returns
    -------
    list of str
        Filled in as the run stops: ``'sent in <place>'`` for every place both
        signals were sent, ``'<place>: <signal> not ignored'`` otherwise.
    """
    events: List[str] = []

    def send_both(place: str) -> None:
        for signum in (signal.SIGTERM, signal.SIGINT):
            if signal.getsignal(signum) is not signal.SIG_IGN:
                events.append(
                    "%s: %s not ignored" % (place, signal.Signals(signum).name)
                )
                return
        for signum in (signal.SIGTERM, signal.SIGINT):
            os.kill(os.getpid(), signum)
        for _ in range(1000):
            pass
        events.append("sent in %s" % place)

    real_report = Simulation._report_unfinished_run
    real_save = lattice_utils.save_out_sim

    def report(self: Simulation, error: BaseException) -> None:
        send_both("report")
        real_report(self, error)

    def save(*args, **kwargs):
        send_both("rescue")
        return real_save(*args, **kwargs)

    monkeypatch.setattr(Simulation, "_report_unfinished_run", report)
    monkeypatch.setattr(lattice_utils, "save_out_sim", save)
    return events


def test_a_second_signal_during_the_clean_up_loses_neither_the_stop_line_nor_the_frames(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("RESTART_FREQ : 5", "SAVE_AT_END : True"))
    before = (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT))
    events = _second_signals_during_the_clean_up(monkeypatch)

    def send_sigterm(sim: Simulation) -> None:
        if signal.getsignal(signal.SIGTERM) is signal.SIG_DFL:
            pytest.fail("run_simulation installed no SIGTERM handler")
        os.kill(os.getpid(), signal.SIGTERM)
        for _ in range(1000):
            pass

    with pytest.raises(SystemExit) as stopped:
        _run_until(d, monkeypatch, 13, send_sigterm)
    assert stopped.value.code == 128 + signal.SIGTERM
    assert events == ["sent in report", "sent in rescue"]

    # the stop line was written although two more signals arrived while it was
    assert (
        len(
            [
                ln
                for ln in _log(d).splitlines()
                if "stopped by SIGTERM after step 13 of 40" in ln
            ]
        )
        == 1
    )
    # and so were the buffered frames: frame 0 plus the XTC_FREQ (2) multiples up to 13
    import mdtraj as md

    traj = md.load(os.path.join(d, "traj.xtc"), top=os.path.join(d, "START.pdb"))
    assert traj.n_frames == 1 + 13 // 2
    # both dispositions are back as they were before the run
    assert (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT)) == before
    assert not os.path.exists(os.path.join(d, restart_module.RUN_MARKER_FILENAME))


def test_a_second_signal_after_ctrl_c_does_not_lose_the_stop_line(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("RESTART_FREQ : 5",))
    before = (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT))
    events = _second_signals_during_the_clean_up(monkeypatch)

    def interrupt(sim: Simulation) -> None:
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        _run_until(d, monkeypatch, 8, interrupt)
    assert events == ["sent in report"]
    assert (
        len([ln for ln in _log(d).splitlines() if "interrupted by SIGINT" in ln]) == 1
    )
    assert (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT)) == before


@pytest.mark.skipif(
    not hasattr(signal, "pthread_sigmask"), reason="needs POSIX signal masks"
)
def test_two_signals_that_arrive_together_give_one_stop_and_its_line(
    tmp_path, monkeypatch
):
    # a kill and a Ctrl-C delivered at the same moment are handled one after the
    # other: the second handler used to raise in the first lines of the clean-up,
    # before the stop had been reported, and the line was lost
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("RESTART_FREQ : 5",))
    before = (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT))

    def send_both_at_once(sim: Simulation) -> None:
        # without handlers of the run's own for BOTH signals nothing is sent: one
        # of them would end the test session
        if (
            signal.getsignal(signal.SIGTERM) is signal.SIG_DFL
            or signal.getsignal(signal.SIGINT) is signal.default_int_handler
        ):
            pytest.fail("run_simulation does not handle both SIGTERM and SIGINT itself")
        both = {signal.SIGTERM, signal.SIGINT}
        signal.pthread_sigmask(signal.SIG_BLOCK, both)
        try:
            os.kill(os.getpid(), signal.SIGTERM)
            os.kill(os.getpid(), signal.SIGINT)
        finally:
            # both are pending now, and are delivered together here
            signal.pthread_sigmask(signal.SIG_UNBLOCK, both)
        for _ in range(1000):
            pass

    with pytest.raises((KeyboardInterrupt, SystemExit)):
        _run_until(d, monkeypatch, 13, send_both_at_once)
    lines = [ln for ln in _log(d).splitlines() if "after step 13 of 40" in ln]
    assert len(lines) == 1
    assert "interrupted by SIGINT" in lines[0] or "stopped by SIGTERM" in lines[0]
    assert "is from step 10" in lines[0]
    assert (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT)) == before


def test_signals_are_back_to_normal_after_a_run_that_finishes(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    before = (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT))
    _run(d)
    assert (signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT)) == before


def test_the_sigterm_handler_does_nothing_once_a_stop_is_under_way(
    tmp_path, monkeypatch
):
    d = str(tmp_path / "run")
    _write_keyfile(d, 6)
    outcome = []

    def call_the_handler(sim: Simulation) -> None:
        handler = signal.getsignal(signal.SIGTERM)
        if handler is signal.SIG_DFL:
            pytest.fail("run_simulation installed no SIGTERM handler")
        sim._stop_in_progress = True
        try:
            outcome.append(handler(signal.SIGTERM, None))
            outcome.append(signal.getsignal(signal.SIGINT)(signal.SIGINT, None))
        finally:
            sim._stop_in_progress = False

    _run_until(d, monkeypatch, 3, call_the_handler)
    # neither SystemExit nor KeyboardInterrupt: the run carried on to the end
    assert outcome == [None, None] and _steps(d, "ENERGY.dat") == list(range(1, 7))


def test_a_stop_inside_a_steps_analyses_is_not_called_after_that_step(
    tmp_path, monkeypatch
):
    # ANALYSIS_FREQ 5: the polymer analysis runs on step 10 and is where the
    # interrupt lands, so step 10's move and energy row exist but its analysis
    # rows do not; "after step 10" would promise more than is on disk
    d = str(tmp_path / "run")
    _write_keyfile(d, 40, extra=("ANALYSIS_FREQ : 5", "RESTART_FREQ : 5"))
    real = Simulation.ANAFUNCT_polymeric_properties

    def interrupted(self: Simulation, step: int) -> None:
        if step == 10:
            raise KeyboardInterrupt
        real(self, step)

    monkeypatch.setattr(Simulation, "ANAFUNCT_polymeric_properties", interrupted)
    with pytest.raises(KeyboardInterrupt):
        _run(d)
    line = [ln for ln in _log(d).splitlines() if "interrupted by SIGINT" in ln][0]
    assert "during the output of step 10 of 40" in line and "after step 10" not in line
    assert "is from step 5" in line
    assert _steps(d, "ENERGY.dat")[-1] == 10 and _steps(d, "RG.dat") == [5]


# --------------------------------------------------------------------------- #
# RF3-3: the documented recipe for joining a stopped segment and its continuation
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("checkpoint_step", [23, 20], ids=["off_cadence", "on_cadence"])
def test_the_documented_joining_recipe_reproduces_the_uninterrupted_trajectory(
    tmp_path, monkeypatch, checkpoint_step
):
    # docs/restart_files.rst: keep the stopped segment's frames up to the
    # checkpoint step (frame 0 plus one per XTC_FREQ multiple at or before it)
    # and ALWAYS drop frame 0 of the continuation, whether or not the checkpoint
    # step is a multiple of XTC_FREQ
    import mdtraj as md

    xtc_freq, total = 2, 60  # XTC_FREQ 2 is what _write_keyfile sets
    extra = ("RESTART_FREQ : %d" % checkpoint_step,)
    reference = str(tmp_path / "reference")
    _write_keyfile(reference, total, extra=extra)
    _run(reference)

    class Stop(Exception):
        pass

    def stop(sim: Simulation) -> None:
        raise Stop

    # stopped three steps past its checkpoint, as a killed run is
    stopped = str(tmp_path / "stopped")
    _write_keyfile(stopped, total, extra=extra)
    with pytest.raises(Stop):
        _run_until(stopped, monkeypatch, checkpoint_step + 3, stop)
    monkeypatch.undo()
    assert _raw(os.path.join(stopped, "restart.pimms"))["STEP"] == checkpoint_step
    resumed = _continuation(str(tmp_path), stopped, n_steps=total, extra=extra)
    _run(resumed)

    def frames(dirpath: str) -> np.ndarray:
        return md.load(
            os.path.join(dirpath, "traj.xtc"), top=os.path.join(dirpath, "START.pdb")
        ).xyz

    keep = 1 + checkpoint_step // xtc_freq
    assert len(frames(stopped)) > keep, (
        "the stopped segment should run past its checkpoint"
    )
    joined = np.concatenate([frames(stopped)[:keep], frames(resumed)[1:]])
    assert (
        joined.shape
        == frames(reference).shape
        == (1 + total // xtc_freq,) + joined.shape[1:]
    )
    assert np.array_equal(joined, frames(reference))

    # the per-step rows join the same way: the stopped segment's rows up to the
    # checkpoint step, then every row of the continuation
    rows = [s for s in _steps(stopped, "ENERGY.dat") if s <= checkpoint_step] + _steps(
        resumed, "ENERGY.dat"
    )
    assert rows == _steps(reference, "ENERGY.dat") == list(range(1, total + 1))
