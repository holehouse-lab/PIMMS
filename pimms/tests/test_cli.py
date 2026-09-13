"""
Command-line tool tests: every ``PIMMS`` command-line mode must work when the
package is used the way a user uses it - from an arbitrary working directory,
through the installed package, with no repo-specific ``PYTHONPATH`` help.

These pin a real failure: a stale ``pimms/`` directory left in site-packages by
an old non-editable install (one orphaned backup file, owned by no distribution)
made Python import ``pimms`` as an ``__init__``-less *namespace* package. The
editable install still resolved every submodule, so simulations ran, but
``pimms.__version__`` did not exist and ``pimms --version`` died with an
``AttributeError``. Tests that only ran ``scripts/PIMMS`` with ``PYTHONPATH``
pointed at the repo could not see it.
"""

import os
import pathlib
import shutil
import subprocess
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
PIMMS_SCRIPT = REPO_ROOT / "scripts" / "PIMMS"


def _clean_env():
    """An environment with no repo-specific import help."""
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONSAFEPATH", None)
    return env


def _run_script(args, cwd, timeout=120):
    """Run ``scripts/PIMMS`` with the test interpreter from ``cwd``."""
    return subprocess.run([sys.executable, str(PIMMS_SCRIPT), *args], cwd=str(cwd),
                          env=_clean_env(), capture_output=True, text=True,
                          timeout=timeout, check=False)


# ---------------------------------------------------------------------------
# the import itself, from outside the repo
# ---------------------------------------------------------------------------
def test_pimms_imports_as_a_regular_package_from_outside_the_repo(tmp_path):
    """``import pimms`` from a foreign cwd must give the real package (with
    ``__init__.py`` and ``__version__``), not a namespace stub.

    If this fails, look for a stray ``pimms/`` directory on ``sys.path`` (e.g.
    ``site-packages/pimms``) that is not owned by any installed distribution.
    """
    probe = ("import pimms, sys\n"
             "print(pimms.__file__)\n"
             "print(list(getattr(pimms, '__path__', [])))\n"
             "print(getattr(pimms, '__version__', 'MISSING'))\n")
    result = subprocess.run([sys.executable, "-c", probe], cwd=str(tmp_path),
                            env=_clean_env(), capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    file_line, path_line, version_line = result.stdout.strip().splitlines()
    assert file_line != "None", (
        "pimms resolved to a NAMESPACE package (no __init__.py) via %s - a stale "
        "directory is shadowing the installed package" % path_line)
    assert file_line.endswith("__init__.py")
    assert version_line != "MISSING"
    assert version_line and version_line != "None"


# ---------------------------------------------------------------------------
# scripts/PIMMS, every mode, from a foreign cwd
# ---------------------------------------------------------------------------
def test_version_flag(tmp_path):
    for flag in ("--version", "-v"):
        result = _run_script([flag], tmp_path)
        assert result.returncode == 0, result.stderr
        assert result.stdout.startswith("version "), result.stdout
        assert result.stdout.strip() != "version None"
        assert result.stdout.strip() != "version unknown"


def test_help_flag(tmp_path):
    result = _run_script(["--help"], tmp_path)
    assert result.returncode == 0, result.stderr
    for flag in ("-keyfile", "--info", "--version"):
        assert flag in result.stdout


def test_no_arguments_prints_usage_hint(tmp_path):
    result = _run_script([], tmp_path)
    assert result.returncode == 0, result.stderr
    assert "try --help" in result.stdout


def test_info_list(tmp_path):
    from pimms import CONFIG

    result = _run_script(["--info"], tmp_path)
    assert result.returncode == 0, result.stderr
    # every documented keyword appears somewhere in the grouped listing
    for keyword in CONFIG.KEYWORDS_DESCRIPTION:
        assert keyword in result.stdout, keyword


def test_info_single_keyword(tmp_path):
    from pimms import CONFIG

    result = _run_script(["--info", "TEMPERATURE"], tmp_path)
    assert result.returncode == 0, result.stderr
    assert "NAME    : TEMPERATURE" in result.stdout
    assert "TYPE    :" in result.stdout
    assert "DESCR   :" in result.stdout
    assert CONFIG.KEYWORDS_DESCRIPTION["TEMPERATURE"][1][:30] in result.stdout


def test_info_all(tmp_path):
    from pimms import CONFIG

    result = _run_script(["--info", "ALL"], tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count("NAME:") == len(CONFIG.KEYWORDS_DESCRIPTION)


def test_info_keyword_lookup_is_case_insensitive(tmp_path):
    """Keyfile keywords are case-insensitive, so `--info hardwall` must work too."""
    result = _run_script(["--info", "hardwall"], tmp_path)
    assert result.returncode == 0, result.stderr
    assert "NAME    : HARDWALL" in result.stdout
    result = _run_script(["--info", "all"], tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count("NAME:") > 10


def test_info_unknown_keyword_exits_nonzero(tmp_path):
    """A typo in `--info <keyword>` is an error, not a successful lookup."""
    result = _run_script(["--info", "NOT_A_KEYWORD"], tmp_path)
    assert result.returncode == 1
    assert "Unknown keyword 'NOT_A_KEYWORD'" in result.stdout


def test_missing_keyfile_exits_nonzero(tmp_path):
    result = _run_script(["-k", "does_not_exist.kf"], tmp_path)
    assert result.returncode == 1
    assert "Could not open file" in result.stdout + result.stderr


def test_keyfile_runs_a_simulation_from_a_foreign_cwd(tmp_path):
    """The full path: -k <keyfile> runs to completion outside the repo."""
    from pimms.tests import kernel_test_utils as U

    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    n_steps=3, equilibration=0, temperature=40, seed=3,
                    # EN_FREQ has to divide into a 3-step run or nothing is ever
                    # written to ENERGY.dat and, since output files are created on
                    # their first row, the file asserted on below never appears.
                    # The default EN_FREQ of 1000 would never fire here.
                    extra={"ENERGY_CHECK": 3, "EN_FREQ": 1})
    result = _run_script(["-k", "KEYFILE.kf"], tmp_path, timeout=300)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    assert "Simulation complete" in result.stdout
    assert (tmp_path / "ENERGY.dat").exists()
    assert (tmp_path / "traj.xtc").exists()


# ---------------------------------------------------------------------------
# the INSTALLED console script (what the user actually types)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("command", ["PIMMS", "pimms"])
def test_installed_console_script_version(tmp_path, command):
    exe = shutil.which(command)
    if exe is None:
        pytest.skip("%s is not on PATH (package not installed as a script)" % command)
    result = subprocess.run([exe, "--version"], cwd=str(tmp_path), env=_clean_env(),
                            capture_output=True, text=True, timeout=120, check=False)
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("version "), result.stdout
    assert result.stdout.strip() not in ("version None", "version unknown")
