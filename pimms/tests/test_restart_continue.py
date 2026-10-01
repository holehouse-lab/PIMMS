"""RESTART_CONTINUE resumes a run bit for bit, and keyfile_used.kf records what ran.

A restart file used to be a configuration only. Resuming from one started a new
run: step numbers from 1, a fresh seed, the temperature from the keyfile. The
checkpoint now also carries the step it was written at, both global generator
states at that moment and the temperature in force, and RESTART_CONTINUE restores
all three, so the resumed segment makes exactly the moves the uninterrupted run
would have made. The decisive checks here compare an uninterrupted run against a
run stopped at a checkpoint and resumed, row by row.

keyfile_used.kf is the run's effective configuration: the resolved keywords after
restart overrides, chain merging and seed generation, written with a provenance
header. The checks here re-parse it, prove the recorded seed reproduces a run
that gave no seed, and read the resolved box and boundary condition from it
through lemonade.
"""

import contextlib
import os
import pickle
import shutil

import numpy as np
import pytest

from pimms import CONFIG
from pimms import restart as restart_module
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import KeyFileException, RestartException
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

_BASE_MOVES = {"MOVE_CRANKSHAFT": 0.5, "MOVE_CHAIN_TRANSLATE": 0.15,
               "MOVE_CHAIN_ROTATE": 0.1, "MOVE_CHAIN_PIVOT": 0.1,
               "MOVE_SLITHER": 0.1, "MOVE_CLUSTER_TRANSLATE": 0.05}


def _write_keyfile(dirpath, n_steps, *, seed=7, moves=None, box=(12, 12, 12), hardwall=False,
                   chains=(("6", "AABBA"), ("5", "A"), ("4", "BBBB")), equilibration=4,
                   temperature=40, extra=()):
    """Write a small keyfile by hand, so every keyword (SEED included) is under control."""
    os.makedirs(dirpath, exist_ok=True)
    U.write_param_file(os.path.join(dirpath, "params.prm"), "SR")
    lines = ["DIMENSIONS : %s" % " ".join(str(b) for b in box),
             "PARAMETER_FILE : params.prm",
             "TEMPERATURE : %s" % temperature,
             "N_STEPS : %d" % n_steps,
             "EQUILIBRATION : %d" % equilibration,
             "HARDWALL : %s" % hardwall,
             "CRANKSHAFT_SUBSTEPS : 200",
             "EN_FREQ : 1", "ANA_POL : 1", "ANA_ACCEPTANCE : 1",
             "XTC_FREQ : 5", "PRINT_FREQ : 100000", "ENERGY_CHECK : 0"]
    if seed is not None:
        lines.append("SEED : %d" % seed)
    for n, seq in chains:
        lines.append("CHAIN : %s %s" % (n, seq))
    for k, v in (moves or _BASE_MOVES).items():
        lines.append("%s : %s" % (k, v))
    lines.extend(extra)
    with open(os.path.join(dirpath, "KEYFILE.kf"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _run(dirpath):
    """Run the keyfile in ``dirpath`` and return the Simulation."""
    cwd = os.getcwd()
    os.chdir(dirpath)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(KeyFileParser("KEYFILE.kf").keyword_lookup)
            sim.run_simulation()
    finally:
        os.chdir(cwd)
    return sim


def _rows(dirpath, name, after_step):
    """Rows of a per-step file with a leading step column, for steps > after_step."""
    out = []
    path = os.path.join(dirpath, name)
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if parts and int(float(parts[0])) > after_step:
                out.append(line.rstrip("\n").rstrip("\t").rstrip(","))
    return out


def _final_state(sim):
    return {c: [list(p) for p in sim.LATTICE.chains[c].get_ordered_positions()]
            for c in sim.LATTICE.chains}


def _stop_and_resume(root, n_total, n_stop, monkeypatch, **kw):
    """Uninterrupted run with a checkpoint at n_stop, and a continuation from it.

    The reference run writes its checkpoints every ``n_stop`` steps; the one
    written at step ``n_stop`` is kept aside (the later ones overwrite
    ``restart.pimms``) and a second run resumes from it to ``n_total``.

    Returns (reference_dir, resumed_dir, reference_sim, resumed_sim).
    """
    ref, cont = (os.path.join(root, d) for d in ("reference", "resumed"))
    kw_ref = dict(kw)
    kw_ref["extra"] = tuple(kw.get("extra", ())) + ("RESTART_FREQ : %d" % n_stop,)
    _write_keyfile(ref, n_total, **kw_ref)

    real_write = restart_module.RestartObject.write_to_file

    def keep_the_one_at_n_stop(self):
        real_write(self)
        if self.step == n_stop:
            shutil.copy(CONFIG.RESTART_FILENAME, "restart_%d.pimms" % n_stop)

    monkeypatch.setattr(restart_module.RestartObject, "write_to_file", keep_the_one_at_n_stop)
    try:
        ref_sim = _run(ref)
    finally:
        monkeypatch.setattr(restart_module.RestartObject, "write_to_file", real_write)

    kw_cont = dict(kw)
    kw_cont["seed"] = None
    kw_cont["extra"] = tuple(kw.get("extra", ())) + (
        "RESTART_FILE : restart.pimms", "RESTART_CONTINUE : True",
        "RESTART_OVERRIDE_DIMENSIONS : True", "RESTART_OVERRIDE_HARDWALL : True")
    _write_keyfile(cont, n_total, **kw_cont)
    shutil.copy(os.path.join(ref, "restart_%d.pimms" % n_stop), os.path.join(cont, "restart.pimms"))
    cont_sim = _run(cont)
    return ref, cont, ref_sim, cont_sim


def _assert_same_continuation(ref, cont, ref_sim, cont_sim, n_stop, files=("ENERGY.dat", "RG.dat", "END_TO_END_DIST.dat")):
    for name in files:
        ref_rows = _rows(ref, name, n_stop)
        cont_rows = _rows(cont, name, n_stop)
        assert ref_rows, "the reference wrote no %s rows after step %d" % (name, n_stop)
        assert cont_rows == ref_rows, "%s differs after step %d" % (name, n_stop)
    assert _final_state(cont_sim) == _final_state(ref_sim)
    assert int(cont_sim.Hamiltonian.evaluate_total_energy(cont_sim.LATTICE)[0]) == \
        int(ref_sim.Hamiltonian.evaluate_total_energy(ref_sim.LATTICE)[0])
    # and the resumed segment really did move things (not a frozen replay)
    assert len(_rows(cont, "ENERGY.dat", n_stop)) >= 5


# --------------------------------------------------------------------------- #
# the checkpoint carries what a continuation needs
# --------------------------------------------------------------------------- #

def test_checkpoint_records_step_generators_and_temperature(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 12)
    sim = _run(d)
    with open(os.path.join(d, "restart.pimms"), "rb") as fh:
        raw = pickle.load(fh)
    assert raw["STEP"] == 12
    assert raw["TEMPERATURE"] == pytest.approx(40.0)
    assert isinstance(raw["RNG_PYTHON"], tuple) and isinstance(raw["RNG_NUMPY"], tuple)
    assert raw["RNG_NUMPY"][0] == "MT19937"
    # the four keys every earlier version wrote are still there, unchanged in kind
    assert set(raw) >= {"CHAINS", "DIMENSIONS", "ENERGY", "HARDWALL"}

    R = restart_module.RestartObject()
    R.build_from_file(os.path.join(d, "restart.pimms"))
    assert R.step == 12 and R.temperature == pytest.approx(40.0)
    assert R.rng_python is not None and R.rng_numpy is not None
    assert R.filename.endswith("restart.pimms")


def test_a_restart_object_built_from_a_lattice_still_writes_the_old_four_key_file(tmp_path):
    d = str(tmp_path / "run")
    _write_keyfile(d, 3, equilibration=0)
    sim = _run(d)
    R = restart_module.RestartObject()
    R.build_from_lattice(sim.LATTICE, sim.hardwall)
    R.set_energy(0)
    cwd = os.getcwd()
    os.chdir(d)
    try:
        R.write_to_file()
    finally:
        os.chdir(cwd)
    with open(os.path.join(d, "restart.pimms"), "rb") as fh:
        raw = pickle.load(fh)
    assert set(raw) == {"CHAINS", "DIMENSIONS", "ENERGY", "HARDWALL"}
    # and such a file still loads, as a plain configuration
    R2 = restart_module.RestartObject()
    R2.build_from_file(os.path.join(d, "restart.pimms"))
    assert R2.step is None and R2.rng_python is None


# --------------------------------------------------------------------------- #
# resumed runs reproduce the uninterrupted run
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("hardwall", [False, True], ids=["PBC", "HW"])
def test_resumed_run_reproduces_the_uninterrupted_run(tmp_path, monkeypatch, hardwall):
    ref, cont, ref_sim, cont_sim = _stop_and_resume(str(tmp_path), 40, 10, monkeypatch, hardwall=hardwall)
    _assert_same_continuation(ref, cont, ref_sim, cont_sim, 10)


def test_resumed_quench_run_restores_the_ramp_temperature(tmp_path, monkeypatch):
    quench = ("QUENCH_RUN : True", "QUENCH_START : 60", "QUENCH_END : 30",
              "QUENCH_STEPSIZE : 10", "QUENCH_FREQ : 6", "QUENCH_AS_EQUILIBRATION : False")
    ref, cont, ref_sim, cont_sim = _stop_and_resume(
        str(tmp_path), 42, 14, monkeypatch, temperature=60, extra=quench, equilibration=2)
    _assert_same_continuation(ref, cont, ref_sim, cont_sim, 14,
                              files=("ENERGY.dat", "RG.dat", "QUENCH.dat"))
    # the checkpoint at step 14 sits two rungs into the ramp (60 -> 50 at 6, -> 40 at 12),
    # so the resumed run has to pick the ramp up at 40, not restart it at 60
    # (read the kept step-14 checkpoint: the resumed run has since overwritten
    # its own restart.pimms with its final one)
    R = restart_module.RestartObject()
    R.build_from_file(os.path.join(ref, "restart_14.pimms"))
    assert R.step == 14 and R.temperature == pytest.approx(40.0)
    assert cont_sim.ACC.temperature == pytest.approx(30.0) == ref_sim.ACC.temperature


def test_resumed_tsmmc_run_rebuilds_the_coordinator_on_the_restored_temperature(tmp_path, monkeypatch):
    moves = {"MOVE_CRANKSHAFT": 0.7, "MOVE_CTSMMC": 0.3}
    tsmmc = ("TSMMC_JUMP_TEMP : 80", "TSMMC_STEP_MULTIPLIER : 2", "TSMMC_NUMBER_OF_POINTS : 3")
    ref, cont, ref_sim, cont_sim = _stop_and_resume(str(tmp_path), 30, 10, monkeypatch, moves=moves, extra=tsmmc)
    _assert_same_continuation(ref, cont, ref_sim, cont_sim, 10)


def test_resumed_parallel_run_reproduces_the_uninterrupted_run(tmp_path, monkeypatch):
    moves = {"MOVE_CRANKSHAFT": 0.6, "MOVE_SLITHER": 0.4}
    par = ("PARALLELIZE : True", "PARALLEL_THREADS : 2")
    ref, cont, ref_sim, cont_sim = _stop_and_resume(
        str(tmp_path), 30, 10, monkeypatch, moves=moves, box=(24, 24, 24),
        chains=(("40", "AABB"), ("20", "A")), extra=par)
    _assert_same_continuation(ref, cont, ref_sim, cont_sim, 10)


def test_resumed_run_labels_rows_with_continued_steps_and_starts_its_trajectory_at_the_checkpoint(tmp_path, monkeypatch):
    ref, cont, ref_sim, cont_sim = _stop_and_resume(str(tmp_path), 30, 10, monkeypatch)
    steps = [int(float(r.split()[0])) for r in _rows(cont, "ENERGY.dat", 0)]
    assert steps == list(range(11, 31))
    import mdtraj as md
    ref_traj = md.load(os.path.join(ref, "traj.xtc"), top=os.path.join(ref, "START.pdb"))
    cont_traj = md.load(os.path.join(cont, "traj.xtc"), top=os.path.join(cont, "START.pdb"))
    # XTC_FREQ 5: reference frames are steps 0,5,...,30; the continuation's frame 0 is
    # the checkpoint (step 10), then 15,...,30 - so its frames are the tail of the reference's
    assert cont_traj.n_frames == 5
    np.testing.assert_allclose(cont_traj.xyz, ref_traj.xyz[2:], atol=1e-6)


# --------------------------------------------------------------------------- #
# what a continuation refuses
# --------------------------------------------------------------------------- #

def _continuation_keyfile(dirpath, n_steps, extra=(), seed=None, **kw):
    _write_keyfile(dirpath, n_steps, seed=seed, extra=tuple(extra) + (
        "RESTART_FILE : restart.pimms", "RESTART_CONTINUE : True"), **kw)


def _parse_in(dirpath):
    cwd = os.getcwd()
    os.chdir(dirpath)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            return KeyFileParser("KEYFILE.kf")
    finally:
        os.chdir(cwd)


@pytest.fixture
def checkpoint(tmp_path):
    """A 10-step run's checkpoint (with continuation state) and a bare one (without)."""
    d = str(tmp_path / "source")
    _write_keyfile(d, 10)
    sim = _run(d)
    bare = str(tmp_path / "bare.pimms")
    R = restart_module.RestartObject()
    R.build_from_lattice(sim.LATTICE, sim.hardwall)
    R.set_energy(0)
    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        R.write_to_file()
        shutil.move(os.path.join(str(tmp_path), "restart.pimms"), bare)
    finally:
        os.chdir(cwd)
    return os.path.join(d, "restart.pimms"), bare


def _refuses(tmp_path, checkpoint_file, match, n_steps=30, seed=None, extra=(), **kw):
    d = str(tmp_path / "attempt")
    _continuation_keyfile(d, n_steps, seed=seed, extra=extra, **kw)
    shutil.copy(checkpoint_file, os.path.join(d, "restart.pimms"))
    with pytest.raises((KeyFileException, RestartException), match=match):
        _parse_in(d)


def test_continuation_refuses_a_checkpoint_without_state(tmp_path, checkpoint):
    _, bare = checkpoint
    _refuses(tmp_path, bare, "no continuation state")


def test_continuation_refuses_an_explicit_seed(tmp_path, checkpoint):
    good, _ = checkpoint
    _refuses(tmp_path, good, "SEED must not be given", seed=99)


def test_continuation_refuses_a_total_length_at_or_before_the_checkpoint(tmp_path, checkpoint):
    good, _ = checkpoint
    _refuses(tmp_path, good, "must exceed the step", n_steps=10)


def test_continuation_refuses_resized_equilibration_and_extra_chains(tmp_path, checkpoint):
    good, _ = checkpoint
    _refuses(tmp_path / "a", good, "RESIZED_EQUILIBRATION",
             extra=("RESIZED_EQUILIBRATION : 9 9 9", "RESTART_OVERRIDE_HARDWALL : True"),
             hardwall=True)
    _refuses(tmp_path / "b", good, "EXTRA_CHAIN", extra=("EXTRA_CHAIN : 1 AB",))


def test_continuation_refuses_a_different_box(tmp_path, checkpoint):
    good, _ = checkpoint
    # the source ran periodic in 12^3; a larger periodic box is refused by the
    # existing restart rules, so use a hardwall source to reach the continuation check
    d = str(tmp_path / "hw_source")
    _write_keyfile(d, 6, hardwall=True)
    _run(d)
    _refuses(tmp_path, os.path.join(d, "restart.pimms"), "DIMENSIONS .* differ",
             box=(16, 16, 16), hardwall=True)


def test_continuation_refuses_a_different_temperature(tmp_path, checkpoint):
    # the checkpoint was written at T = 40; the continuation would run at 40 while
    # the keyfile (and so the T_NORM angle scaling and keyfile_used.kf) said 50
    good, _ = checkpoint
    _refuses(tmp_path, good, "TEMPERATURE 50 differs from the restart file's 40", temperature=50)


def test_continuation_refuses_a_checkpoint_temperature_off_the_quench_ramp(tmp_path, checkpoint):
    good, _ = checkpoint
    quench = ("QUENCH_RUN : True", "QUENCH_START : 60", "QUENCH_END : 50",
              "QUENCH_STEPSIZE : 5", "QUENCH_FREQ : 6", "QUENCH_AS_EQUILIBRATION : False")
    _refuses(tmp_path, good, "outside this keyfile's quench range", temperature=60, extra=quench)


def test_continuation_requires_a_restart_file(tmp_path):
    d = str(tmp_path / "attempt")
    _write_keyfile(d, 30, extra=("RESTART_CONTINUE : True",))
    with pytest.raises(KeyFileException, match="requires a RESTART_FILE"):
        _parse_in(d)


def test_a_bare_checkpoint_still_seeds_a_new_run_without_continue(tmp_path, checkpoint):
    _, bare = checkpoint
    d = str(tmp_path / "fresh")
    _write_keyfile(d, 6, seed=3, extra=("RESTART_FILE : restart.pimms",))
    shutil.copy(bare, os.path.join(d, "restart.pimms"))
    sim = _run(d)
    assert sim.continue_from_step == 0
    assert [int(float(r.split()[0])) for r in _rows(d, "ENERGY.dat", 0)] == list(range(1, 7))


# --------------------------------------------------------------------------- #
# keyfile_used.kf
# --------------------------------------------------------------------------- #

def test_effective_keyfile_is_written_reparses_and_records_a_reproducing_seed(tmp_path):
    # a run that gives NO seed: the effective keyfile must record the one that was used
    d = str(tmp_path / "noseed")
    _write_keyfile(d, 12, seed=None)
    _run(d)
    used = os.path.join(d, CONFIG.EFFECTIVE_KEYFILE_NAME)
    assert os.path.exists(used)
    text = open(used).read()
    assert text.startswith("# Effective configuration of this PIMMS run")
    assert "generated at start-up" in text

    cwd = os.getcwd()
    os.chdir(d)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            parsed = KeyFileParser(CONFIG.EFFECTIVE_KEYFILE_NAME)
    finally:
        os.chdir(cwd)
    recorded_seed = parsed.keyword_lookup["SEED"]
    assert isinstance(recorded_seed, int) and recorded_seed > 0
    assert parsed.keyword_lookup["N_STEPS"] == 12
    assert parsed.keyword_lookup["CHAIN"] == [[6, "AABBA"], [5, "A"], [4, "BBBB"]]

    # the recorded seed reproduces the run
    d2 = str(tmp_path / "replay")
    _write_keyfile(d2, 12, seed=recorded_seed)
    _run(d2)
    assert _rows(d2, "ENERGY.dat", 0) == _rows(d, "ENERGY.dat", 0)
    assert _rows(d2, "RG.dat", 0) == _rows(d, "RG.dat", 0)
    assert "from the keyfile" in open(os.path.join(d2, CONFIG.EFFECTIVE_KEYFILE_NAME)).read()


def test_effective_keyfile_holds_the_overridden_box_and_boundary(tmp_path):
    # a hardwall 9-cube source, restarted by a keyfile that literally says 14-cube periodic
    src = str(tmp_path / "source")
    _write_keyfile(src, 6, box=(9, 9, 9), hardwall=True, seed=5)
    _run(src)
    d = str(tmp_path / "restarted")
    _write_keyfile(d, 6, box=(14, 14, 14), hardwall=False, seed=5, extra=(
        "RESTART_FILE : restart.pimms", "RESTART_OVERRIDE_DIMENSIONS : True",
        "RESTART_OVERRIDE_HARDWALL : True", "ANA_CLUSTER : 1"))
    shutil.copy(os.path.join(src, "restart.pimms"), os.path.join(d, "restart.pimms"))
    sim = _run(d)
    assert sim.hardwall is True and list(sim.LATTICE.dimensions) == [9, 9, 9]

    used = os.path.join(d, CONFIG.EFFECTIVE_KEYFILE_NAME)
    text = open(used).read()
    assert "started from       : restart file" in text and "written at step 6" in text
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            parsed = KeyFileParser(CONFIG.EFFECTIVE_KEYFILE_NAME, parse_only=True)
    finally:
        os.chdir(cwd)
    assert parsed.keyword_lookup["DIMENSIONS"] == [9, 9, 9]
    assert parsed.keyword_lookup["HARDWALL"] is True
    assert parsed.keyword_lookup["RESTART_OVERRIDE_DIMENSIONS"] is False
    assert parsed.keyword_lookup["RESTART_OVERRIDE_HARDWALL"] is False

    # lemonade reads the truth from it, with no restart file to consult
    from pimms import lemonade
    from pimms.lemonade import phase_separation as ps
    traj = lemonade.load(xtc=os.path.join(d, "traj.xtc"), pdb=os.path.join(d, "START.pdb"), keyfile=used)
    assert traj.hardwall is True and tuple(traj.dimensions) == (9, 9, 9)
    expected = {}
    for line in open(os.path.join(d, "NUM_CLUSTERS.dat")):
        step, n = line.split()[:2]
        expected[int(step)] = int(n)
    got = ps.number_of_clusters(traj, min_beads=1)
    # frames are steps 0,5 (XTC_FREQ 5, 6 steps); NUM_CLUSTERS rows are every step
    assert int(got[1]) == expected[5]


def test_effective_keyfile_records_the_resized_equilibration_phase(tmp_path):
    d = str(tmp_path / "resized")
    _write_keyfile(d, 10, box=(14, 14, 14), equilibration=4, seed=2,
                   extra=("RESIZED_EQUILIBRATION : 9 9 9", "SAVE_EQ : True"))
    _run(d)
    text = open(os.path.join(d, CONFIG.EFFECTIVE_KEYFILE_NAME)).read()
    assert "steps 1 to 4 ran in box 9 9 9 under HARDWALL : True" in text
    assert "apply from step 5" in text


# --------------------------------------------------------------------------- #
# continuation keyfiles written back out
# --------------------------------------------------------------------------- #

def test_a_continuation_keyfile_written_by_write_keyfile_reparses(tmp_path, checkpoint):
    # after restart processing RESTART_FILE is a loaded object and is not written,
    # so RESTART_CONTINUE : True used to be written without it and the file was
    # refused on re-parse ("RESTART_CONTINUE : True requires a RESTART_FILE")
    good, _ = checkpoint
    d = str(tmp_path / "attempt")
    _continuation_keyfile(d, 30, extra=("RESTART_OVERRIDE_DIMENSIONS : True",))
    shutil.copy(good, os.path.join(d, "restart.pimms"))
    parsed = _parse_in(d)
    out = os.path.join(d, "written.kf")
    parsed.write_keyfile(out)
    text = open(out).read()
    assert "RESTART_CONTINUE :" in text and "RESTART_FILE" not in text
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            reparsed = KeyFileParser("written.kf")
    finally:
        os.chdir(cwd)
    assert reparsed.keyword_lookup["RESTART_CONTINUE"] is False
    assert reparsed.keyword_lookup["RESTART_OVERRIDE_DIMENSIONS"] is False
    assert list(reparsed.keyword_lookup["DIMENSIONS"]) == [12, 12, 12]


def test_the_summary_counts_the_frames_a_continuation_writes(tmp_path, monkeypatch, capsys):
    # the continuation writes the checkpoint as frame 0 plus the XTC_FREQ multiples
    # after its step, not the whole run's frame count
    ref, cont, _ref_sim, cont_sim = _stop_and_resume(str(tmp_path), 30, 10, monkeypatch)
    import mdtraj as md
    n_frames = md.load(os.path.join(cont, "traj.xtc"), top=os.path.join(cont, "START.pdb")).n_frames
    capsys.readouterr()
    # the resumed run has since overwritten restart.pimms with its final checkpoint
    shutil.copy(os.path.join(ref, "restart_10.pimms"), os.path.join(cont, "restart.pimms"))
    parsed = _parse_in(cont)
    parsed.print_summary()
    line = [ln for ln in capsys.readouterr().out.splitlines() if "Expected number of frames" in ln]
    assert line and int(line[0].split(":")[1]) == n_frames == 1 + (30 // 5 - 10 // 5)


def test_the_effective_keyfile_says_a_continuation_did_not_use_its_seed(tmp_path, monkeypatch):
    # a continuation's generators come from the restart file, so the header must not
    # claim that the generated seed reproduces the run
    _ref, cont, _ref_sim, _cont_sim = _stop_and_resume(str(tmp_path), 30, 10, monkeypatch)
    text = open(os.path.join(cont, CONFIG.EFFECTIVE_KEYFILE_NAME)).read()
    seed_line = [ln for ln in text.splitlines() if ln.startswith("# SEED")][0]
    assert "unused" in seed_line and "reproduces this run" not in seed_line
