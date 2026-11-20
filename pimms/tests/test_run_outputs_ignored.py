"""
Every file a PIMMS run writes must be ignored by git in demo_keyfiles.

The demos are run where they sit, so a ``git add -A`` after running one would
otherwise sweep its trajectory, analysis tables and logs into a commit. Two
checks keep the repository's ignore rules complete:

* every output name PIMMS declares (CONFIG's OUTNAME_* and OUTPUT_* names, the
  per-chain-type variants, the restart file and its checkpoint temporaries, the
  run marker and the trajectory files) is ignored in every demo directory, in a
  new one and in a run directory inside a demo; and
* so is every file a real run leaves behind, which catches an output that a
  later PIMMS writes under a name the first check doesn't know about.

Both ask git itself (``git check-ignore --no-index``), so they test the rules
git actually applies, every .gitignore included. They're skipped outside a git
checkout, since an installed PIMMS has no .gitignore.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from pimms import CONFIG
from pimms import restart
from pimms.tests import kernel_test_utils as U


REPO_ROOT = Path(__file__).resolve().parents[2]
PIMMS_SCRIPT = REPO_ROOT / "scripts" / "PIMMS"
DEMO_DIRECTORY = "demo_keyfiles"

# a demo that doesn't exist yet, and a run directory inside one, so that the
# rules are checked for future demos as well as the current ones
NEW_DEMO_DIRECTORIES = ("a_new_demo", "a_new_demo/run_1")

# the trajectory and structure files a run writes, which simulation.py names
# directly rather than through CONFIG (eq_ files come from RESIZED_EQUILIBRATION,
# CONFIG_AT_ENERGY_FAIL ones from a failed energy check)
TRAJECTORY_NAMES = ("START.pdb", "traj.xtc", "eq_START.pdb", "eq_traj.xtc",
                    "CONFIG_AT_ENERGY_FAIL.pdb", "CONFIG_AT_ENERGY_FAIL.xtc")

# CONFIG attributes that name an output file
OUTPUT_NAME_PREFIXES = ("OUTNAME_", "OUTPUT_")
OTHER_OUTPUT_ATTRIBUTES = ("RESTART_FILENAME", "EFFECTIVE_KEYFILE_NAME", "QUENCHFILE_NAME")

# chain-type indices used for the CHAIN_<type>_ variants (one and two digits)
CHAIN_TYPE_INDICES = (0, 12)

# a checkpoint temporary is restart.pimms.tmp.<pid>; an older PIMMS wrote restart.pimms.tmp
CHECKPOINT_SUFFIXES = (".tmp", ".tmp.12345")

# The real run: two chain types (so the CHAIN_<type>_ files are written), every
# analysis on, a quench, a resized equilibration box and the chain-to-chainID
# table, so that it writes as many kinds of output as one short run can
RUN_BOX = [14, 14, 14]
RUN_EQUILIBRATION_BOX = "12 12 12"
RUN_CHAINS = [(4, "AABBA"), (4, "BBABB")]
RUN_N_STEPS = 40
RUN_EQUILIBRATION = 20
RUN_KEYWORDS = {
    "ANALYSIS_FREQ": 5, "ANA_RESIDUE_PAIRS": "1 3", "ANA_CLUSTER_THRESHOLD": 0,
    "WRITE_CHAIN_TO_CHAINID": "True", "RESIZED_EQUILIBRATION": RUN_EQUILIBRATION_BOX,
    "QUENCH_RUN": "True", "QUENCH_START": 60, "QUENCH_END": 55, "QUENCH_STEPSIZE": 5, "QUENCH_FREQ": 10,
    "QUENCH_AS_EQUILIBRATION": "False",
    "XTC_FREQ": 5, "EN_FREQ": 5, "PRINT_FREQ": 20, "ENERGY_CHECK": 10, "RESTART_FREQ": 10}
RUN_INPUTS = ("KEYFILE.kf", "params.prm")
RUN_TIMEOUT_SECONDS = 600


def _require_git_checkout() -> Path:
    """
    The repository root, or skip the calling test if this isn't a git checkout.

    Returns
    -------
    Path
        The root of the git work tree that holds this test.
    """
    if shutil.which("git") is None:
        pytest.skip("needs git")
    result = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=str(REPO_ROOT), capture_output=True,
                            text=True, check=False)
    if result.returncode != 0 or Path(result.stdout.strip()).resolve() != REPO_ROOT:
        pytest.skip("needs the git checkout of PIMMS; %s is an installed copy" % REPO_ROOT)
    return REPO_ROOT


def _demo_directories(root: Path) -> list[str]:
    """
    Every directory a demo run could write into: each current demo, a new one and a run directory inside one.

    Parameters
    ----------
    root : Path
        The repository root.

    Returns
    -------
    list of str
        Paths relative to the root.
    """
    current = sorted(path.name for path in (root / DEMO_DIRECTORY).iterdir() if path.is_dir())
    assert current, "found no demo directories in %s" % DEMO_DIRECTORY
    return ["%s/%s" % (DEMO_DIRECTORY, name) for name in current + list(NEW_DEMO_DIRECTORIES)]


def _not_ignored(root: Path, names: set[str]) -> list[str]:
    """
    The paths, among every name in every demo directory, that git would not ignore.

    Parameters
    ----------
    root : Path
        The repository root.
    names : set of str
        File names a run writes.

    Returns
    -------
    list of str
        ``<demo directory>/<name>`` for each name git would let through, sorted.
    """
    paths = sorted("%s/%s" % (directory, name) for directory in _demo_directories(root) for name in names)
    # --no-index tests the rules alone, so the paths needn't exist; git prints the
    # ignored ones and exits 1 if there are none
    result = subprocess.run(["git", "check-ignore", "--no-index", "--stdin"], cwd=str(root),
                            input="\n".join(paths) + "\n", capture_output=True, text=True, check=False)
    assert result.returncode in (0, 1), result.stderr
    ignored = set(result.stdout.splitlines())
    return [path for path in paths if path not in ignored]


def _declared_output_names() -> set[str]:
    """
    Every output file name PIMMS declares.

    Returns
    -------
    set of str
        The analysis outputs and their CHAIN_<type>_ variants, every CONFIG
        OUTNAME_* and OUTPUT_* name, the restart file and its checkpoint
        temporaries, the effective keyfile, the quench log, the run marker and
        the trajectory files.
    """
    names = set(CONFIG.analysis_output_files(chain_types=list(CHAIN_TYPE_INDICES)))
    for attribute in dir(CONFIG):
        value = getattr(CONFIG, attribute)
        if attribute.startswith(OUTPUT_NAME_PREFIXES) and isinstance(value, str):
            names.add(value)
    names.update(getattr(CONFIG, attribute) for attribute in OTHER_OUTPUT_ATTRIBUTES)
    names.update(CONFIG.RESTART_FILENAME + suffix for suffix in CHECKPOINT_SUFFIXES)
    names.add(restart.RUN_MARKER_FILENAME)
    names.update(TRAJECTORY_NAMES)
    return names


def test_every_output_name_pimms_declares_is_ignored_in_demo_keyfiles():
    root = _require_git_checkout()
    names = _declared_output_names()
    # the scan found the outputs at all
    assert {"ENERGY.dat", "log.txt", "chain_to_chainid.txt", "pimms_running.pid", "traj.xtc"} <= names
    assert _not_ignored(root, names) == []


def test_every_file_a_run_leaves_is_ignored_in_demo_keyfiles(tmp_path):
    """The run's own record of what it wrote, so an output under a new name can't slip through."""
    root = _require_git_checkout()
    U.write_param_file(tmp_path / "params.prm", "SR")
    U.write_keyfile(tmp_path / "KEYFILE.kf", 3, False, {"MOVE_CRANKSHAFT": 1.0}, box=RUN_BOX, chains=RUN_CHAINS,
                    n_steps=RUN_N_STEPS, equilibration=RUN_EQUILIBRATION, extra=RUN_KEYWORDS)
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(root)
    result = subprocess.run([sys.executable, str(PIMMS_SCRIPT), "-k", "KEYFILE.kf"], cwd=str(tmp_path),
                            env=environment, capture_output=True, text=True, timeout=RUN_TIMEOUT_SECONDS,
                            check=False)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]

    written = {path.name for path in tmp_path.iterdir()} - set(RUN_INPUTS)
    # the run wrote what it was asked to, so the check covers them
    expected = {"traj.xtc", "eq_traj.xtc", "START.pdb", "eq_START.pdb", "restart.pimms", "QUENCH.dat",
                "chain_to_chainid.txt", "RES_TO_RES_DIST.dat", "CLUSTER_RG.dat", "CHAIN_1_DISTANCE_MAP.dat"}
    assert expected <= written, sorted(expected - written)
    assert all(path.is_file() for path in tmp_path.iterdir()), "the run made a directory"
    assert _not_ignored(root, written) == []
