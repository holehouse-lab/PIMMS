"""Regression tests for the lightweight command-line information paths."""

from pathlib import Path
import os
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
PIMMS_SCRIPT = REPO_ROOT / "scripts" / "PIMMS"


def _source_tree_environment():
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    return env


def test_version_command_does_not_import_simulation_stack():
    """Version reporting must not import MDTraj or the native MC kernels."""
    probe = (
        "import runpy, sys; "
        f"sys.argv = [{str(PIMMS_SCRIPT)!r}, '--version']; "
        f"runpy.run_path({str(PIMMS_SCRIPT)!r}, run_name='__main__'); "
        "assert 'pimms.simulation' not in sys.modules; "
        "assert 'pimms.keyfile_parser' not in sys.modules"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe],
        env=_source_tree_environment(),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("version ")
