from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pytest


@dataclass(frozen=True)
class ExpectedOutputFile:
    source_filename: str
    expected_by_test: dict[int, str]


def _repo_root() -> Path:
    """Locate the PIMMS repository root from this test module location.

    The root is identified by the presence of both ``setup.py`` and the
    executable script at ``scripts/PIMMS`` while walking upward from this file.

    Returns
    -------
    Path
        Absolute path to the repository root directory.

    Raises
    ------
    RuntimeError
        Raised when no parent directory matches the expected repository
        structure.
    """
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "setup.py").exists() and (parent / "scripts" / "PIMMS").exists():
            return parent
    raise RuntimeError(f"Unable to locate repository root from {current}")


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
    """Build the command used to invoke the PIMMS executable from THIS working tree.

    Deliberately always runs ``scripts/PIMMS`` from the repo root with the current
    interpreter - NOT a ``PIMMS`` found on ``PATH``. A ``PIMMS`` on PATH is whatever
    is pip-installed, which in a non-editable install is a DIFFERENT (older) version
    than the tree under test; the regression suite then silently validated the
    installed package instead of the code being changed. ``_run_single_testset``
    forces the repo onto ``PYTHONPATH`` of the subprocess so the script imports the
    tree's ``pimms`` too, and asserts as much before running.

    Returns
    -------
    list[str]
        Command token list suitable for ``subprocess.run``.
    """
    script_path = _repo_root() / "scripts" / "PIMMS"
    return [sys.executable, str(script_path)]


def _cleanup_generated_outputs(test_dir: Path) -> None:
    """Remove generated simulation outputs from a fixture test directory.

    This keeps the fixture directory deterministic before each run by deleting
    previous output artifacts while preserving input/configuration files
    (``.prm``, ``.kf``, and files starting with ``KEYFILE``).

    Parameters
    ----------
    test_dir : Path
        Path to an individual ``test_<n>`` fixture directory.

    Returns
    -------
    None
        Performs in-place filesystem cleanup.
    """
    for path in test_dir.iterdir():
        if not path.is_file():
            continue

        if path.suffix in {".prm", ".kf"}:
            continue

        if path.name.startswith("KEYFILE"):
            continue

        if path.suffix in {".dat", ".xtc", ".pdb", ".txt", ".lat"} or path.name == "restart.pimms":
            path.unlink(missing_ok=True)


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


def _run_single_testset(test_num: int) -> tuple[Path, dict[str, str]]:
    """Execute one simulation fixture and collect observed final output lines.

    This helper performs fixture cleanup, runs PIMMS in the selected test
    directory, writes a per-test log file, asserts successful execution, and
    then captures observed final lines for files that have expected data for
    that test number.

    Parameters
    ----------
    test_num : int
        Numeric fixture identifier corresponding to ``test_<test_num>``.

    Returns
    -------
    tuple[Path, dict[str, str]]
        Tuple containing:
        1. Path to the executed fixture directory.
        2. Mapping of source filename to observed final non-empty line.

    Raises
    ------
    AssertionError
        Raised when fixture directory is missing, simulation exits non-zero,
        or an expected output file for that test is not found.
    """
    testsuite = _testsuite_root()

    test_dir = testsuite / f"test_{test_num}"
    assert test_dir.exists(), f"Simulation fixture directory missing: {test_dir}"

    _cleanup_generated_outputs(test_dir)

    cmd = _resolve_pimms_command() + ["-k", "KEYFILE.kf"]

    # Force the subprocess to import THIS working tree's pimms, not a pip-installed
    # copy that may shadow it, and prove it before we trust the run's output.
    repo_root = str(_repo_root())
    env = dict(os.environ)
    env["PYTHONPATH"] = repo_root + os.pathsep + env.get("PYTHONPATH", "")
    probe = subprocess.run(
        [sys.executable, "-c", "import pimms, sys; sys.stdout.write(pimms.__file__)"],
        cwd=str(test_dir), capture_output=True, text=True, env=env, check=False,
    )
    assert probe.returncode == 0 and probe.stdout.startswith(repo_root), (
        "regression subprocess would import the wrong pimms: "
        f"{probe.stdout!r} (expected under {repo_root!r}); stderr={probe.stderr!r}"
    )

    result = subprocess.run(
        cmd,
        cwd=str(test_dir),
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
        env=env,
    )

    log_path = testsuite / f"test_{test_num}/pytest_test_{test_num}_log.txt"
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

        source_path = test_dir / source_filename
        assert source_path.exists(), f"Expected {source_filename} not found for test_{test_num}"
        observed_final_lines[source_filename] = _read_final_nonempty_line(source_path)
        observed_line_counts[source_filename] = sum(
            1 for ln in source_path.read_text().splitlines() if ln.strip())

    return test_dir, observed_final_lines, observed_line_counts


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
def run_simulation_testset():
    """Expose the single-testset simulation runner as a session fixture.

    Returns
    -------
    Callable[[int], tuple[Path, dict[str, str]]]
        Callable that accepts a test number, runs that simulation fixture,
        and returns the fixture directory path plus observed final-line values.
    """
    return _run_single_testset
