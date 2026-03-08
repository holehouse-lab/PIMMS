from __future__ import annotations

import shutil
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
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "setup.py").exists() and (parent / "scripts" / "PIMMS").exists():
            return parent
    raise RuntimeError(f"Unable to locate repository root from {current}")


def _testsuite_root() -> Path:
    # Keep all simulation fixtures local to this test package.
    return Path(__file__).resolve().parent


def _expected_output_root() -> Path:
    return _testsuite_root() / "expected_output"


def _resolve_pimms_command() -> list[str]:
    pimms_exe = shutil.which("PIMMS")
    if pimms_exe:
        return [pimms_exe]

    script_path = _repo_root() / "scripts" / "PIMMS"
    return [sys.executable, str(script_path)]


def _cleanup_generated_outputs(test_dir: Path) -> None:
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
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    assert lines, f"{path.name} is empty for {path.parent.name}"
    return lines[-1]


def _parse_expected_output_file(path: Path) -> ExpectedOutputFile:
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
    testsuite = _testsuite_root()

    test_dir = testsuite / f"test_{test_num}"
    assert test_dir.exists(), f"Simulation fixture directory missing: {test_dir}"

    _cleanup_generated_outputs(test_dir)

    cmd = _resolve_pimms_command() + ["-k", "KEYFILE.kf"]
    result = subprocess.run(
        cmd,
        cwd=str(test_dir),
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )

    log_path = testsuite / f"test_{test_num}/pytest_test_{test_num}_log.txt"
    log_path.write_text(result.stdout + "\n" + result.stderr)

    assert result.returncode == 0, (
        f"Simulation test_{test_num} failed with return code {result.returncode}. "
        f"See log: {log_path}"
    )

    observed_final_lines: dict[str, str] = {}
    expected_output_data = _load_expected_outputs()
    for source_filename, expected_by_test in expected_output_data.items():
        # Only validate/capture files that have an expected value for this test.
        if test_num not in expected_by_test:
            continue

        source_path = test_dir / source_filename
        assert source_path.exists(), f"Expected {source_filename} not found for test_{test_num}"
        observed_final_lines[source_filename] = _read_final_nonempty_line(source_path)

    return test_dir, observed_final_lines


@pytest.fixture(scope="session")
def expected_output_data() -> dict[str, dict[int, str]]:
    return _load_expected_outputs()


@pytest.fixture(scope="session")
def run_simulation_testset():
    return _run_single_testset
