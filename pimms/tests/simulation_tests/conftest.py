from __future__ import annotations

import os
import shutil
import subprocess
import sys
import sysconfig
from collections.abc import Callable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pytest


@dataclass(frozen=True)
class ExpectedOutputFile:
    source_filename: str
    expected_by_test: dict[int, str]


def _import_root() -> Path:
    """Return the directory that holds the ``pimms`` package these tests belong to.

    This file is ``<root>/pimms/tests/simulation_tests/conftest.py``, so the root is
    three levels above its directory. In a source checkout that is the repository
    root; in an installed copy it is ``site-packages``. It is the directory that has
    to come first on the subprocess ``PYTHONPATH`` for the run to import the same
    ``pimms`` the tests were collected from.

    Returns
    -------
    Path
        Absolute path to the parent directory of the ``pimms`` package.
    """
    return Path(__file__).resolve().parents[3]


def _source_tree_root() -> Path | None:
    """Return the repository root if these tests are running from a source checkout.

    A checkout is recognised by ``setup.py`` and ``scripts/PIMMS`` sitting directly
    in :func:`_import_root`. We look in that one directory only. We used to walk up
    every parent until some directory matched, which goes wrong for an installed
    copy whose environment lives underneath a checkout (a virtual environment
    created inside the repository, say): the walk left ``site-packages``, found the
    enclosing checkout, and the regression suite then ran that checkout's
    ``scripts/PIMMS`` against that checkout's ``pimms`` - so it reported on a
    different build from the one installed, without saying so.

    Returns
    -------
    Path or None
        The repository root, or None when the tests are in an installed package.
    """
    root = _import_root()
    if (root / "setup.py").is_file() and (root / "scripts" / "PIMMS").is_file():
        return root
    return None


def _installed_pimms_script() -> Path:
    """Locate the ``PIMMS`` executable that was installed with this interpreter.

    Used only when the tests run from an installed package, where there is no
    ``scripts/PIMMS`` to point at. We look in the scripts directory of the running
    interpreter's environment, then next to the interpreter itself, and only then
    on ``PATH``. Whichever script is found is run with the current interpreter and
    with :func:`_import_root` first on ``PYTHONPATH``, and ``_run_single_testset``
    checks which ``pimms`` that imports before it trusts the run.

    Returns
    -------
    Path
        Path to the ``PIMMS`` script.

    Raises
    ------
    RuntimeError
        If no ``PIMMS`` executable can be found. The message says where we looked.
    """
    candidates = [
        Path(sysconfig.get_path("scripts")) / "PIMMS",
        Path(sys.executable).parent / "PIMMS",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate

    on_path = shutil.which("PIMMS")
    if on_path is not None:
        return Path(on_path)

    raise RuntimeError(
        "The simulation regression tests are running from an installed copy of PIMMS "
        f"({_import_root()}) and need the PIMMS executable, but could not find it. "
        f"Looked for {candidates[0]} and {candidates[1]}, and for PIMMS on PATH. "
        "Reinstall PIMMS into this environment, or run the tests from a source checkout."
    )


def _testsuite_root() -> Path:
    """Return the root directory for simulation regression fixtures.

    Returns
    -------
    Path
        Absolute path to the directory containing this ``conftest.py`` file
        and all ``test_*`` simulation fixture directories.
    """
    # Keep all simulation fixtures local to this test package.
    return Path(__file__).resolve().parent


def _expected_output_root() -> Path:
    """Return the directory containing generated expected-output artifacts.

    Returns
    -------
    Path
        Path to the ``expected_output`` folder inside the simulation test
        suite directory.
    """
    return _testsuite_root() / "expected_output"


def _resolve_pimms_command() -> list[str]:
    """Build the command that runs the PIMMS these tests belong to.

    From a source checkout this deliberately always runs ``scripts/PIMMS`` from the
    repo root with the current interpreter - NOT a ``PIMMS`` found on ``PATH``. A
    ``PIMMS`` on PATH is whatever is pip-installed, which in a non-editable install
    is a DIFFERENT (older) version than the tree under test; the regression suite
    then silently validated the installed package instead of the code being changed.

    From an installed copy (the tests ship in the wheel so that an install can be
    checked) there is no ``scripts/PIMMS``, so we run the executable installed with
    the current interpreter (see :func:`_installed_pimms_script`).

    In both cases ``_run_single_testset`` puts :func:`_import_root` first on the
    subprocess ``PYTHONPATH``, so the script imports the ``pimms`` these tests sit
    in, and asserts as much before running.

    Returns
    -------
    list[str]
        Command token list suitable for ``subprocess.run``.

    Raises
    ------
    RuntimeError
        From an installed copy, if no ``PIMMS`` executable can be found.
    """
    source_root = _source_tree_root()
    if source_root is not None:
        return [sys.executable, str(source_root / "scripts" / "PIMMS")]
    return [sys.executable, str(_installed_pimms_script())]


def _is_generated_output(path: Path) -> bool:
    """Say whether a file in a fixture directory is a run output rather than an input.

    A fixture directory in a working tree can hold both: the tracked inputs
    (keyfile, parameter file, and for some scenarios a restart or freeze file) and
    whatever a run made by hand, or by an older version of this harness, left
    behind. Parameter files (``.prm``), keyfiles (``.kf``) and anything whose name
    starts with ``KEYFILE`` are always inputs. Tables, trajectories, structures,
    logs (``.dat``, ``.xtc``, ``.pdb``, ``.txt``, ``.lat``) and ``restart.pimms``
    are outputs. Anything else (``in.pimms``, ``frz.in``) is an input.

    Parameters
    ----------
    path : Path
        A file inside a ``test_<n>`` fixture directory.

    Returns
    -------
    bool
        True if PIMMS (or the harness) wrote the file, False if a run reads it.
    """
    if path.suffix in {".prm", ".kf"} or path.name.startswith("KEYFILE"):
        return False
    return path.suffix in {".dat", ".xtc", ".pdb", ".txt", ".lat"} or path.name == "restart.pimms"


def _copy_fixture_inputs(test_dir: Path, run_dir: Path) -> None:
    """Copy the input files of one fixture into the directory the run will use.

    The harness used to run PIMMS inside the fixture directory itself, after
    deleting the previous outputs there. That wrote into the source tree on every
    test run (and deleted the one tracked output, ``test_11``'s
    ``absolute_energies_of_angles.txt``), and two test runs on one checkout
    deleted each other's files part-way through. We now leave the fixture
    directory alone and give every run a directory of its own.

    Parameters
    ----------
    test_dir : Path
        The ``test_<n>`` fixture directory. Nothing in it is changed.
    run_dir : Path
        An existing, empty directory to copy the inputs into.

    Returns
    -------
    None
    """
    for path in sorted(test_dir.iterdir()):
        if path.is_file() and not _is_generated_output(path):
            shutil.copy2(path, run_dir / path.name)


def _read_final_nonempty_line(path: Path) -> str:
    """Read and return the last non-empty line from a text file.

    Parameters
    ----------
    path : Path
        Output file to inspect.

    Returns
    -------
    str
        Final non-blank line after stripping whitespace.

    Raises
    ------
    AssertionError
        Raised when the file contains no non-empty lines.
    """
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    assert lines, f"{path.name} is empty for {path.parent.name}"
    return lines[-1]


def _parse_expected_output_file(path: Path) -> ExpectedOutputFile:
    """Parse one consolidated expected-output file into structured data.

    The expected file format is tab-delimited lines of the form:
    ``test_<n>\t<expected_final_line>``. Empty lines are ignored.

    Parameters
    ----------
    path : Path
        Path to a ``*.final_lines.txt`` expected-output file.

    Returns
    -------
    ExpectedOutputFile
        Parsed source filename and expected final-line values indexed by test
        number.

    Raises
    ------
    AssertionError
        Raised when filename or line format is invalid.
    """
    suffix = ".final_lines.txt"
    assert path.name.endswith(suffix), f"Unexpected expected-output filename: {path.name}"

    source_filename = path.name[: -len(suffix)].replace("__", "/")
    expected_by_test: dict[int, str] = {}

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue

        parts = line.split("\t", maxsplit=1)
        assert len(parts) == 2, f"Malformed line in {path.name}: {raw_line!r}"

        test_label, expected_value = parts
        assert test_label.startswith("test_"), f"Malformed test label in {path.name}: {test_label!r}"
        test_num = int(test_label.removeprefix("test_"))
        expected_by_test[test_num] = expected_value

    return ExpectedOutputFile(source_filename=source_filename, expected_by_test=expected_by_test)


@lru_cache(maxsize=1)
def _load_expected_outputs() -> dict[str, dict[int, str]]:
    """Load all expected-output files into an in-memory lookup table.

    Results are cached for the life of the Python process to avoid repeated
    disk reads across multiple tests in the same session.

    Returns
    -------
    dict[str, dict[int, str]]
        Mapping from source filename (for example ``ENERGY.dat``) to a mapping
        of ``test_number -> expected_final_line``.

    Raises
    ------
    AssertionError
        Raised if expected-output directory or files are missing.
    """
    expected_root = _expected_output_root()
    assert expected_root.exists(), f"Expected output directory missing: {expected_root}"

    expected_files = sorted(expected_root.glob("*.final_lines.txt"))
    assert expected_files, f"No expected-output files found in {expected_root}"

    loaded: dict[str, dict[int, str]] = {}
    for path in expected_files:
        parsed = _parse_expected_output_file(path)
        loaded[parsed.source_filename] = parsed.expected_by_test

    return loaded


def _run_single_testset(test_num: int, run_dir: Path) -> tuple[Path, dict[str, str], dict[str, int]]:
    """Execute one simulation fixture and collect observed final output lines.

    This helper copies the fixture's input files into ``run_dir``, runs PIMMS
    there, writes a per-test log file there, asserts successful execution, and
    then captures observed final lines for files that have expected data for
    that test number. The fixture directory in the tree is only read.

    Parameters
    ----------
    test_num : int
        Numeric fixture identifier corresponding to ``test_<test_num>``.
    run_dir : Path
        An existing, empty directory to run the simulation in. Every output,
        and the log ``pytest_test_<test_num>_log.txt``, is written here.

    Returns
    -------
    tuple[Path, dict[str, str], dict[str, int]]
        Tuple containing:
        1. Path to the directory the simulation ran in (``run_dir``).
        2. Mapping of source filename to observed final non-empty line.
        3. Mapping of source filename to its number of non-empty lines.

    Raises
    ------
    AssertionError
        Raised when fixture directory is missing, simulation exits non-zero,
        or an expected output file for that test is not found.
    """
    testsuite = _testsuite_root()

    test_dir = testsuite / f"test_{test_num}"
    assert test_dir.exists(), f"Simulation fixture directory missing: {test_dir}"

    _copy_fixture_inputs(test_dir, run_dir)

    cmd = _resolve_pimms_command() + ["-k", "KEYFILE.kf"]

    # Force the subprocess to import the pimms these tests belong to (the working
    # tree in a checkout, the installed package otherwise), not another copy that may
    # shadow it, and prove it before we trust the run's output.
    import_root = str(_import_root())
    env = dict(os.environ)
    env["PYTHONPATH"] = import_root + os.pathsep + env.get("PYTHONPATH", "")
    probe = subprocess.run(
        [sys.executable, "-c", "import pimms, sys; sys.stdout.write(pimms.__file__)"],
        cwd=str(run_dir), capture_output=True, text=True, env=env, check=False,
    )
    expected_package = os.path.join(import_root, "pimms")
    imported_package = os.path.dirname(os.path.realpath(probe.stdout)) if probe.stdout else ""
    assert probe.returncode == 0 and imported_package == expected_package, (
        "regression subprocess would import the wrong pimms: "
        f"{probe.stdout!r} (expected the package in {expected_package!r}); "
        f"stderr={probe.stderr!r}"
    )

    result = subprocess.run(
        cmd,
        cwd=str(run_dir),
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
        env=env,
    )

    log_path = run_dir / f"pytest_test_{test_num}_log.txt"
    log_path.write_text(result.stdout + "\n" + result.stderr)

    assert result.returncode == 0, (
        f"Simulation test_{test_num} failed with return code {result.returncode}. "
        f"See log: {log_path}"
    )

    observed_final_lines: dict[str, str] = {}
    observed_line_counts: dict[str, int] = {}
    expected_output_data = _load_expected_outputs()
    for source_filename, expected_by_test in expected_output_data.items():
        # Only validate/capture files that have an expected value for this test.
        if test_num not in expected_by_test:
            continue

        source_path = run_dir / source_filename
        assert source_path.exists(), (
            f"Expected {source_filename} not found for test_{test_num} (run directory: {run_dir}, "
            f"log: {log_path})"
        )
        observed_final_lines[source_filename] = _read_final_nonempty_line(source_path)
        observed_line_counts[source_filename] = sum(
            1 for ln in source_path.read_text().splitlines() if ln.strip())

    return run_dir, observed_final_lines, observed_line_counts


@lru_cache(maxsize=1)
def _load_expected_line_counts() -> dict[str, dict[int, int]]:
    """Load the companion ``*.line_counts.txt`` files (may be absent for old
    baselines). Mapping: source filename -> {test_num: non-empty line count}."""
    expected_root = _expected_output_root()
    loaded: dict[str, dict[int, int]] = {}
    for path in sorted(expected_root.glob("*.line_counts.txt")):
        source_filename = path.name[: -len(".line_counts.txt")].replace("__", "/")
        by_test: dict[int, int] = {}
        for raw_line in path.read_text().splitlines():
            line = raw_line.strip()
            if not line:
                continue
            label, value = line.split("\t", maxsplit=1)
            by_test[int(label.removeprefix("test_"))] = int(value)
        loaded[source_filename] = by_test
    return loaded


@pytest.fixture(scope="session")
def expected_line_counts() -> dict[str, dict[int, int]]:
    """Expected per-file non-empty line counts (write-cadence fingerprint)."""
    return _load_expected_line_counts()


@pytest.fixture(scope="session")
def expected_output_data() -> dict[str, dict[int, str]]:
    """Provide cached expected output values for the full test session.

    Returns
    -------
    dict[str, dict[int, str]]
        Mapping of source output filenames to per-test expected final lines.
    """
    return _load_expected_outputs()


@pytest.fixture(scope="session")
def run_simulation_testset(tmp_path_factory: pytest.TempPathFactory) -> Callable[[int], tuple[Path, dict[str, str], dict[str, int]]]:
    """Expose the single-testset simulation runner as a session fixture.

    Each call runs its scenario in a fresh directory made by pytest's
    ``tmp_path_factory``, so nothing is written into the fixture directories and
    two test runs on the same checkout cannot disturb one another. pytest keeps
    the directories of the last few sessions, so the outputs and the log of a
    failed scenario can still be inspected afterwards; the failure message gives
    the path.

    Parameters
    ----------
    tmp_path_factory : pytest.TempPathFactory
        pytest's session-scoped temporary-directory factory.

    Returns
    -------
    Callable[[int], tuple[Path, dict[str, str], dict[str, int]]]
        Callable that accepts a test number, runs that simulation fixture, and
        returns the run directory, the observed final-line values and the
        observed line counts.
    """
    def _run(test_num: int) -> tuple[Path, dict[str, str], dict[str, int]]:
        return _run_single_testset(test_num, tmp_path_factory.mktemp(f"pimms_regression_test_{test_num}_"))

    return _run
