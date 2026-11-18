"""Fixes from the November 2026 audit of the engine, its parsers and the run life cycle.

Each test fails on the code before its fix for the reason of the bug.
"""

import contextlib
import importlib
import io
import os
import pickle
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

import pimms
from pimms import CONFIG, crankshaft_list_functions, data_structures, keyfile_parser, parameterfile_parser
from pimms import restart as restart_module
from pimms.energy import Hamiltonian
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import AnalysisRoutineException, KeyFileException, MoveException, SimulationException
from pimms.simulation import Simulation

# the short-range interaction and solvation energy of every bead type the tests use
PAIR_ENERGY = -1
SOLVATION_ENERGY = -2
PIMMS_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "PIMMS"
# the conventional exit status of a process stopped by SIGTERM (128 + 15), which
# PIMMS uses after its clean-up
SIGTERM_EXIT_STATUS = 143
SUBPROCESS_TIMEOUT_S = 300


def write_params(path, beads, pair_energies=None):
    """Write a short-range parameter file for ``beads``, every pair at
    ``PAIR_ENERGY`` unless ``pair_energies`` gives it, with an angle penalty
    and a solvation energy for each bead."""
    pair_energies = pair_energies or {}
    lines = ["ANGLE_PENALTY\t%s\t1\t0\t0" % bead for bead in beads]
    lines += ["%s\t0\t%d" % (bead, SOLVATION_ENERGY) for bead in beads]
    for i, first in enumerate(beads):
        for second in beads[i:]:
            energy = pair_energies.get((first, second), pair_energies.get((second, first), PAIR_ENERGY))
            lines.append("%s  %s\t%s" % (first, second, energy))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_keyfile(directory, chains, extra=(), n_steps=10, seed=3, box=(12, 12, 12)):
    """Write a small crankshaft-only keyfile beside params.prm."""
    lines = ["DIMENSIONS : %s" % " ".join(str(b) for b in box), "PARAMETER_FILE : params.prm",
             "TEMPERATURE : 40", "N_STEPS : %d" % n_steps, "EQUILIBRATION : 2", "MOVE_CRANKSHAFT : 1.0",
             "CRANKSHAFT_SUBSTEPS : 50", "XTC_FREQ : 5", "EN_FREQ : 1", "ENERGY_CHECK : 0",
             "PRINT_FREQ : 100000"]
    if seed is not None:
        lines.append("SEED : %d" % seed)
    lines += ["CHAIN : %d %s" % (count, sequence) for count, sequence in chains]
    lines += list(extra)
    Path(directory, "KEYFILE.kf").write_text("\n".join(lines) + "\n", encoding="utf-8")


def quietly(function, *args, **kwargs):
    """Call ``function`` with its standard output captured, returning (result, output)."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        result = function(*args, **kwargs)
    return result, buffer.getvalue()


# --------------------------------------------------------------------------- #
# ANGLE_PENALTY_T_NORM: exact halves are rounded half-to-even, as documented
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("penalty, temperature, applied", [("0.70", 45, 32), ("0.14", 75, 10)])
def test_t_norm_angle_penalty_halves_are_rounded_to_even(tmp_path, monkeypatch, penalty, temperature, applied):
    """0.70 x 45 is 31.5, which the documented rule takes to 32; in binary
    floating point the product was 31.499999999999996 and was applied as 31."""
    monkeypatch.chdir(tmp_path)
    Path("params.prm").write_text("ANGLE_PENALTY_T_NORM\tA\t%s\t0\t0\nA\t0\t-2\nA  A\t-1\n" % penalty)
    with contextlib.redirect_stdout(io.StringIO()):
        angles = parameterfile_parser.parse_angles("params.prm", temperature)
    assert int(round(angles["A"][0])) == applied
    assert angles["A"][0] == float(Fraction(penalty) * temperature)


def test_t_norm_angle_penalty_with_a_huge_exponent_is_read_at_once(tmp_path, monkeypatch):
    """The parser accepts 0e99999999; exact arithmetic must not expand it."""
    monkeypatch.chdir(tmp_path)
    Path("params.prm").write_text("ANGLE_PENALTY_T_NORM\tA\t0e99999999\t1e-9999999\t%s\nA\t0\t-2\nA  A\t-1\n"
                                  % ("1." + "0" * 5000))
    with contextlib.redirect_stdout(io.StringIO()):
        angles = parameterfile_parser.parse_angles("params.prm", 45)
    assert angles["A"] == [0.0, 0.0, 45.0]


# --------------------------------------------------------------------------- #
# CASE_INSENSITIVE_CHAINS never changes a chain's length
# --------------------------------------------------------------------------- #

def test_upper_casing_a_chain_keeps_its_length_and_its_bead_types(tmp_path, monkeypatch):
    """'ß'.upper() is 'SS', so upper-casing the whole sequence turned the
    three-bead AßA into the four-bead ASSA; 'ı'.upper() is 'I', which turned a
    bead of one type into another."""
    monkeypatch.chdir(tmp_path)
    write_params("params.prm", ["A", "I", "S", "ß", "ı"])
    write_keyfile(tmp_path, [(2, "aßa"), (1, "aıa")])
    parser, _ = quietly(KeyFileParser, "KEYFILE.kf")
    assert [sequence for _, sequence in parser.keyword_lookup["CHAIN"]] == ["AßA", "AıA"]
    assert keyfile_parser.upper_case_chain("azßıſéε") == "AZßıſÉΕ"


# --------------------------------------------------------------------------- #
# the freeze-file summary does not expand a B range
# --------------------------------------------------------------------------- #

def test_freeze_summary_prints_a_huge_b_range_without_expanding_it(tmp_path, monkeypatch):
    """A mistyped B range such as 0-99999999999 used to be expanded bead by
    bead in the start-up summary, hanging PIMMS before the range check could
    refuse it."""
    monkeypatch.chdir(tmp_path)
    write_params("params.prm", ["A"])
    Path("freeze.txt").write_text("B 1 0-99999999999\n")
    write_keyfile(tmp_path, [(2, "AAAA")], extra=["FREEZE_FILE : freeze.txt"])
    parser, _ = quietly(KeyFileParser, "KEYFILE.kf")

    real_format = keyfile_parser.format_bead_runs

    def refuse_long_ranges(indices, *args, **kwargs):
        if isinstance(indices, range) and len(indices) > 10**6:
            raise AssertionError("the summary expanded a range of %d beads" % len(indices))
        return real_format(indices, *args, **kwargs)

    monkeypatch.setattr(keyfile_parser, "format_bead_runs", refuse_long_ranges)
    _, summary = quietly(parser.print_summary)
    assert "chain 1, beads 0-99999999999" in summary
    # and the range check then refuses it, once the chain lengths are known
    with pytest.raises(KeyFileException, match="freezes bead 99999999999 of chain 1"):
        quietly(Simulation, parser.keyword_lookup)


# --------------------------------------------------------------------------- #
# frozen_bead_rows: the vectorised mask equals a plain loop
# --------------------------------------------------------------------------- #

class _Chain:
    def __init__(self, length):
        self.positions = [0] * length


class _Lattice:
    def __init__(self, lengths_by_id):
        self.chains = {chainID: _Chain(length) for chainID, length in lengths_by_id.items()}
        self.crankshaft_lists = np.zeros((sum(lengths_by_id.values()), 8), dtype=np.int64)
        self.chain_layout = None


def _oracle_rows(lengths_by_id, frozen_chains, frozen_beads):
    """The frozen-bead mask written out bead by bead: rows in ascending chainID."""
    rows = []
    for chainID in sorted(lengths_by_id):
        for bead in range(lengths_by_id[chainID]):
            if chainID in frozen_beads:
                rows.append(bead in frozen_beads[chainID])
            else:
                rows.append(chainID in frozen_chains)
    return np.array(rows, dtype=bool)


def test_frozen_bead_rows_equals_a_bead_by_bead_oracle():
    rng = np.random.default_rng(2026)
    for _ in range(500):
        ids = sorted(rng.choice(np.arange(1, 60), size=int(rng.integers(1, 15)), replace=False).tolist())
        lengths = {chainID: int(rng.integers(1, 9)) for chainID in ids}
        frozen = [chainID for chainID in ids if rng.random() < 0.6]
        partly = {}
        for chainID in frozen:
            if rng.random() < 0.6:
                size = int(rng.integers(0, lengths[chainID] + 1))
                partly[chainID] = tuple(sorted(rng.choice(lengths[chainID], size=size, replace=False).tolist()))
        rows = crankshaft_list_functions.frozen_bead_rows(_Lattice(lengths), frozen, partly)
        assert np.array_equal(rows, _oracle_rows(lengths, frozen, partly))


def test_frozen_bead_rows_still_refuses_inconsistent_frozen_sets():
    lattice = _Lattice({1: 4, 2: 3})
    with pytest.raises(MoveException, match="Chain 2 has frozen beads but is not among the frozen chains"):
        crankshaft_list_functions.frozen_bead_rows(lattice, [1], {1: (0,), 2: (0,)})
    with pytest.raises(MoveException, match="A frozen bead of chain 2 is outside the chain"):
        crankshaft_list_functions.frozen_bead_rows(lattice, [1, 2], {1: (0,), 2: (3,)})


# --------------------------------------------------------------------------- #
# a free stretch held at both ends by frozen beads is announced
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("frozen, stretches", [((0, 4), [(1, 3)]), ((0, 2, 4), []), ((2,), []),
                                               ((0, 3, 6, 7, 10), [(1, 2), (4, 5), (8, 9)])])
def test_enclosed_free_stretches(frozen, stretches):
    assert data_structures.enclosed_free_stretches(frozen) == stretches


@pytest.mark.parametrize("freeze_line, warned", [("B 1 0 4", True), ("B 1 0 2 4", False), ("B 1 2", False)])
def test_start_up_warns_about_free_stretches_single_bead_moves_can_trap(tmp_path, monkeypatch, freeze_line, warned):
    """Single-bead crankshaft moves cannot reach every arrangement of two or
    more free beads held at both ends (two of them closing a 2x2 square with
    their frozen neighbours never move), which nothing used to say."""
    monkeypatch.chdir(tmp_path)
    write_params("params.prm", ["A", "B"])
    Path("freeze.txt").write_text(freeze_line + "\n")
    write_keyfile(tmp_path, [(1, "AABBA"), (2, "BBA")], extra=["FREEZE_FILE : freeze.txt"])
    parser, _ = quietly(KeyFileParser, "KEYFILE.kf")
    _, output = quietly(Simulation, parser.keyword_lookup)
    assert ("held at both ends by frozen beads" in output) is warned
    if warned:
        assert "chain 1 beads 1-3" in output


# --------------------------------------------------------------------------- #
# keyfile_used.kf re-runs a lower-case restart as the same system
# --------------------------------------------------------------------------- #

def test_effective_keyfile_of_a_lower_case_restart_keeps_its_beads(tmp_path, monkeypatch):
    """pdb2pimms writes nucleotides as lower-case beads. The echo of a run
    started from such a restart file carried CASE_INSENSITIVE_CHAINS : True,
    so re-running it upper-cased the chains into a different system."""
    monkeypatch.chdir(tmp_path)
    write_params("params.prm", ["a", "B"])
    restart = {"CHAINS": {1: [[[1, 1, 1], [2, 1, 1], [3, 1, 1], [4, 1, 1]], "aaaa", 0],
                          2: [[[1, 5, 5], [2, 5, 5], [3, 5, 5]], "BBB", 1]},
               "DIMENSIONS": [12, 12, 12], "ENERGY": 0, "HARDWALL": True}
    with open("start.pimms", "wb") as handle:
        pickle.dump(restart, handle)
    write_keyfile(tmp_path, [], extra=["RESTART_FILE : start.pimms", "HARDWALL : True"])
    parser, _ = quietly(KeyFileParser, "KEYFILE.kf")
    simulation, _ = quietly(Simulation, parser.keyword_lookup)
    quietly(simulation.write_effective_keyfile, parser.keyword_lookup, filename="echo.kf")
    echo, _ = quietly(KeyFileParser, "echo.kf")
    assert sorted(sequence for _, sequence in echo.keyword_lookup["CHAIN"]) == ["BBB", "aaaa"]


# --------------------------------------------------------------------------- #
# RESTART_CONTINUE refuses a changed energy function that leaves the
# checkpoint's energy alone
# --------------------------------------------------------------------------- #

# a box in which a crankshaft-only run cannot bring the two chains of
# _checkpoint_run within interaction range of each other
FAR_APART_BOX = (30, 30, 30)
# the largest Chebyshev distance at which two beads interact (the super-long-range shell)
INTERACTION_RANGE = 3


def _build_in(directory, extra, n_steps, seed, run):
    """Parse the keyfile in ``directory`` and build (and, if ``run``, run) the Simulation."""
    write_keyfile(directory, [], extra=["HARDWALL : True"] + extra, n_steps=n_steps, seed=seed,
                  box=FAR_APART_BOX)
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        parser, _ = quietly(KeyFileParser, "KEYFILE.kf")
        simulation, _ = quietly(Simulation, parser.keyword_lookup)
        if run:
            quietly(simulation.run_simulation)
    finally:
        os.chdir(cwd)


def _checkpoint_run(directory):
    """Run 10 crankshaft-only steps of an A chain and a C chain that start far
    apart, and return the checkpoint path."""
    directory.mkdir()
    write_params(directory / "params.prm", ["A", "C"])
    start = {"CHAINS": {1: [[[x, 2, 2] for x in range(2, 6)], "AAAA", 0],
                        2: [[[x, 20, 20] for x in range(2, 6)], "CCCC", 1]},
             "DIMENSIONS": list(FAR_APART_BOX), "ENERGY": 0, "HARDWALL": True}
    with open(directory / "start.pimms", "wb") as handle:
        pickle.dump(start, handle)
    _build_in(directory, ["RESTART_FILE : start.pimms"], n_steps=10, seed=3, run=True)
    return directory / "restart.pimms"


def _continue(directory, checkpoint, params_pairs):
    directory.mkdir()
    write_params(directory / "params.prm", ["A", "C"], params_pairs)
    shutil.copy(checkpoint, directory / "restart.pimms")
    _build_in(directory, ["RESTART_FILE : restart.pimms", "RESTART_CONTINUE : True"], n_steps=20,
              seed=None, run=False)


def test_continuation_refuses_a_parameter_file_the_energy_check_cannot_see(tmp_path):
    """The A and C chains are out of each other's reach in the checkpoint, so
    changing the A-C interaction leaves its energy unchanged; the continuation
    used to be accepted, and then ran under the new interaction."""
    checkpoint = _checkpoint_run(tmp_path / "source")
    with open(checkpoint, "rb") as handle:
        chains = pickle.load(handle)["CHAINS"]
    a_beads, c_beads = np.array(chains[1][0]), np.array(chains[2][0])
    closest = np.abs(a_beads[:, None, :] - c_beads[None, :, :]).max(axis=2).min()
    assert closest > INTERACTION_RANGE
    _continue(tmp_path / "same", checkpoint, {})
    with pytest.raises(SimulationException, match="energy function this keyfile builds is not the one"):
        _continue(tmp_path / "changed", checkpoint, {("A", "C"): 300})


def test_energy_function_digest_ignores_line_order_and_sees_every_table(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    write_params("one.prm", ["A", "B"])
    lines = Path("one.prm").read_text().splitlines()
    Path("reordered.prm").write_text("\n".join(reversed(lines)) + "\n")
    write_params("changed.prm", ["A", "B"], {("A", "B"): -2})

    def digest(name, **kwargs):
        hamiltonian, _ = quietly(Hamiltonian, name, 3, kwargs.get("non_interacting", False),
                                 kwargs.get("angles_off", False), temperature=40)
        return hamiltonian.energy_function_digest()

    reference = digest("one.prm")
    assert digest("reordered.prm") == reference
    assert digest("changed.prm") != reference
    assert digest("one.prm", angles_off=True) != reference
    assert digest("one.prm", non_interacting=True) != reference


def test_restart_file_refuses_a_malformed_energy_function_digest(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    restart = {"CHAINS": {1: [[[1, 1, 1], [2, 1, 1]], "AA", 0]}, "DIMENSIONS": [12, 12, 12],
               "ENERGY": 0, "HARDWALL": True, "ENERGY_FUNCTION_DIGEST": "not a digest"}
    with open("bad.pimms", "wb") as handle:
        pickle.dump(restart, handle)
    with pytest.raises(Exception, match="ENERGY_FUNCTION_DIGEST must be a SHA-256 digest"):
        restart_module.RestartObject().build_from_file("bad.pimms")


# --------------------------------------------------------------------------- #
# a custom analysis that calls sys.exit() stops the run as an error
# --------------------------------------------------------------------------- #

def test_custom_analysis_calling_sys_exit_is_an_error(tmp_path, monkeypatch):
    """SystemExit is not an Exception, so it went straight through the wrapper
    and ended the run part-way through with exit status 0."""
    monkeypatch.chdir(tmp_path)
    write_params("params.prm", ["A"])
    Path("leaving.py").write_text("import sys\n\ndef analysis_function(step, lattice):\n    sys.exit()\n")
    write_keyfile(tmp_path, [(2, "AAAA")], n_steps=10,
                  extra=["ANALYSIS_MODULE : leaving.py", "ANA_CUSTOM : 5"])
    parser, _ = quietly(KeyFileParser, "KEYFILE.kf")
    simulation, _ = quietly(Simulation, parser.keyword_lookup)
    with pytest.raises(AnalysisRoutineException, match=r"called sys.exit\(\)"):
        quietly(simulation.run_simulation)


def test_sigterm_during_a_custom_analysis_is_still_a_signal_stop(tmp_path):
    """PIMMS's SIGTERM handler raises SystemExit wherever Python is, which can
    be inside the user's analysis function; that is a scheduler's stop, not an
    error in the user's code."""
    write_params(tmp_path / "params.prm", ["A"])
    (tmp_path / "stopping.py").write_text(
        "import os\nimport signal\n\n"
        "def analysis_function(step, lattice):\n"
        "    if step == 10:\n"
        "        os.kill(os.getpid(), signal.SIGTERM)\n")
    write_keyfile(tmp_path, [(2, "AAAA")], n_steps=20,
                  extra=["ANALYSIS_MODULE : stopping.py", "ANA_CUSTOM : 5"])
    result = subprocess.run([sys.executable, str(PIMMS_SCRIPT), "-k", "KEYFILE.kf"], cwd=tmp_path,
                            capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT_S)
    assert result.returncode == SIGTERM_EXIT_STATUS, result.stderr[-2000:]
    assert "stopped by SIGTERM" in (tmp_path / "log.txt").read_text(encoding="utf-8")
    assert "AnalysisRoutineException" not in result.stderr


# --------------------------------------------------------------------------- #
# outputs linked in from elsewhere are replaced, never written through
# --------------------------------------------------------------------------- #

def test_angle_summary_and_run_marker_replace_a_symbolic_link(tmp_path, monkeypatch):
    run = tmp_path / "run"
    other = tmp_path / "other"
    run.mkdir()
    other.mkdir()
    monkeypatch.chdir(run)
    for name in (CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL, restart_module.RUN_MARKER_FILENAME):
        (other / name).write_text("another run's file\n")
        os.symlink(other / name, name)

    with contextlib.redirect_stdout(io.StringIO()):
        parameterfile_parser.write_angle_parameter_summary({"A": [1, 0, 0]}, "params.prm")
    restart_module._write_run_markers([(12345, "host", "now")])

    for name in (CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL, restart_module.RUN_MARKER_FILENAME):
        assert (other / name).read_text() == "another run's file\n"
        assert not os.path.islink(name)
    assert "ANGLE PARAMETER SUMMARY" in Path(CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL).read_text()


# --------------------------------------------------------------------------- #
# the cluster-threshold warning says what is written
# --------------------------------------------------------------------------- #

def test_cluster_threshold_warning_says_the_files_hold_step_only_rows(tmp_path, monkeypatch):
    """The per-cluster files are still written, one step-only row per cluster
    step; the warning said they would never be written."""
    monkeypatch.chdir(tmp_path)
    write_params("params.prm", ["A"])
    write_keyfile(tmp_path, [(3, "AAAA")], extra=["ANA_CLUSTER_THRESHOLD : 3"])
    _, output = quietly(KeyFileParser, "KEYFILE.kf")
    assert "will hold only the step number" in output
    assert "will never be written" not in output


# --------------------------------------------------------------------------- #
# a cluster kernel built against another numpy falls back with a warning
# --------------------------------------------------------------------------- #

class _NumpyMismatchFinder:
    """Make importing pimms.cluster_kernels raise what an extension built
    against another numpy raises."""

    def find_spec(self, name, path=None, target=None):
        if name == "pimms.cluster_kernels":
            raise ValueError("numpy.dtype size changed, may indicate binary incompatibility")
        return None


def test_a_cluster_kernel_built_against_another_numpy_falls_back(monkeypatch):
    pytest.importorskip("pimms.cluster_kernels")
    from pimms import cluster_utils

    monkeypatch.delattr(pimms, "cluster_kernels")
    monkeypatch.delitem(sys.modules, "pimms.cluster_kernels")
    monkeypatch.setattr(sys, "meta_path", [_NumpyMismatchFinder()] + sys.meta_path)
    try:
        with pytest.warns(UserWarning, match="numpy.dtype size changed"):
            importlib.reload(cluster_utils)
        assert cluster_utils._HAVE_CLUSTER_KERNELS is False
    finally:
        monkeypatch.undo()
        importlib.reload(cluster_utils)
    assert cluster_utils._HAVE_CLUSTER_KERNELS is True
