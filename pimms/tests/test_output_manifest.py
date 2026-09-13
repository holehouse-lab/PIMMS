"""
Tests for the output manifest and for lazy output-file creation.

PIMMS output files are created lazily: a writer opens its file at the moment it
has a row to put in it, so a file exists if and only if the run wrote to it.
Up to 1.0.8 start-up instead pre-created about 25 files, so every run left a
pile of zero-length files behind - including files for analyses that were
switched off, analyses whose frequency never fired, and analyses that fired but
had nothing to report.

The flip side of lazy creation is that start-up must DELETE any copy of an
output this run could write, otherwise a re-run in a used directory silently
inherits the previous run's data for the files it does not itself write. That
guarantee (test_rerun_inherits_nothing_from_a_previous_run) is the important one
here.
"""

import contextlib
import os
import re

import pytest

from pimms import CONFIG
from pimms.tests import kernel_test_utils as U


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# a small, fast system: four short chains in a 12^3 box
_BOX = [12, 12, 12]
_CHAINS = [(4, "AABBA")]

# frequencies that keep the runs quick and the output small
_QUIET = {"PRINT_FREQ": 1000, "XTC_FREQ": 1000, "ENERGY_CHECK": 0}


@pytest.fixture(autouse=True)
def _restore_cwd():
    """The runs chdir into their tmp_path; put the working directory back so
    later tests do not depend on this module having run."""
    cwd = os.getcwd()
    yield
    os.chdir(cwd)


def _run(tmp_path, extra, chains=None, n_steps=30, equilibration=2):
    """Run one short simulation in tmp_path and return the set of files present.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to build and run the simulation in.

    extra : dict
        Extra keyfile keywords, merged over the quiet defaults.

    chains : list of tuple or None, optional
        ``(count, sequence)`` chain specification. Defaults to a single chain
        type.

    n_steps : int, optional
        Number of simulation steps.

    equilibration : int, optional
        Number of equilibration steps.

    Returns
    -------
    set
        Names of the files present in ``tmp_path`` once the run has finished.
    """
    keywords = dict(_QUIET)
    keywords.update(extra)

    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=_BOX, chains=chains if chains is not None else _CHAINS,
                          n_steps=n_steps, equilibration=equilibration, extra=keywords)
    os.chdir(tmp_path)
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state.sim.run_simulation()

    return set(os.listdir(tmp_path))


# the analysis output files written by ANA_CLUSTER
_CLUSTER_FILES = ("CLUSTERS.dat", "NUM_CLUSTERS.dat", "CLUSTER_RG.dat", "CLUSTER_ASPH.dat",
                  "CLUSTER_VOL.dat", "CLUSTER_AREA.dat", "CLUSTER_DEN.dat",
                  "LR_CLUSTERS.dat", "NUM_LR_CLUSTERS.dat", "LR_CLUSTER_RG.dat",
                  "LR_CLUSTER_ASPH.dat", "LR_CLUSTER_VOL.dat", "LR_CLUSTER_AREA.dat",
                  "LR_CLUSTER_DEN.dat")

# ...and by ANA_POL
_POL_FILES = ("RG.dat", "ASPH.dat", "END_TO_END_DIST.dat")


# ---------------------------------------------------------------------------
# 1. the manifest covers every file the writers can open
# ---------------------------------------------------------------------------

def test_manifest_covers_every_output_file_named_in_analysis_IO():
    """The manifest is only as good as its coverage. Every CONFIG filename that
    the writers in analysis_IO touch must be in it, or start-up would leave that
    file behind from a previous run."""
    source = open(os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "analysis_IO.py")).read()

    referenced = set()
    for name in re.findall(r"CONFIG\.([A-Z_0-9]+)", source):
        value = getattr(CONFIG, name, None)
        if isinstance(value, str) and value.endswith(".dat"):
            referenced.add(name)

    # sanity: the scan found the writers at all
    assert "OUTNAME_ENERGY" in referenced
    assert "QUENCHFILE_NAME" in referenced

    assert referenced <= set(CONFIG.ANALYSIS_OUTPUT_NAMES), \
        sorted(referenced - set(CONFIG.ANALYSIS_OUTPUT_NAMES))


def test_analysis_output_files_expands_the_per_chain_type_variants():
    """The manifest has to know about the CHAIN_<type>_ files, since those are
    what a multi-component run writes instead of (or beside) the plain names."""
    plain = CONFIG.analysis_output_files()
    expanded = CONFIG.analysis_output_files(chain_types=[0, 1])

    assert set(plain) < set(expanded)
    assert "CHAIN_0_" + CONFIG.OUTNAME_CLUSTERS in expanded
    assert "CHAIN_1_" + CONFIG.OUTNAME_DMAP in expanded
    # the per-step files that have no per-type variant must not gain one
    assert "CHAIN_0_" + CONFIG.OUTNAME_ENERGY not in expanded


# ---------------------------------------------------------------------------
# 2. a disabled analysis produces none of its files
# ---------------------------------------------------------------------------

def test_disabled_cluster_analysis_writes_none_of_its_files(tmp_path):
    """ANA_CLUSTER : 0 used to leave 14 zero-length cluster files on disk."""
    produced = _run(tmp_path, {"ANA_CLUSTER": 0, "ANA_POL": 0, "ANA_ACCEPTANCE": 0,
                               "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0})

    for name in _CLUSTER_FILES:
        assert name not in produced, name
    assert "CLUSTER_RADIAL_DENSITY_PROFILE.dat" not in produced


def test_enabled_cluster_analysis_writes_all_of_its_files(tmp_path):
    """The mirror image of the test above: with the analysis on, every one of
    those files is there (so the test above is not passing because the analysis
    is broken)."""
    produced = _run(tmp_path, {"ANA_CLUSTER": 5, "ANA_POL": 0, "ANA_ACCEPTANCE": 0,
                               "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0})

    for name in _CLUSTER_FILES:
        assert name in produced, name
        assert (tmp_path / name).stat().st_size > 0, name


# ---------------------------------------------------------------------------
# 3. an analysis that never fires produces no file
# ---------------------------------------------------------------------------

def test_analysis_that_never_fires_writes_no_file(tmp_path):
    """ANA_POL is enabled but its frequency never lands inside the 30-step
    production window, so nothing is ever written and no file should appear.
    This is the case that cannot be predicted from the keyfile alone."""
    produced = _run(tmp_path, {"ANA_POL": 1000, "ANA_CLUSTER": 0, "ANA_ACCEPTANCE": 0,
                               "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0,
                               "EN_FREQ": 1000})

    for name in _POL_FILES:
        assert name not in produced, name

    # EN_FREQ must be positive, but a 30-step run never reaches step 1000 either
    assert "ENERGY.dat" not in produced


def test_analysis_that_fires_writes_its_file(tmp_path):
    """Same run with a frequency that does land in the production window."""
    produced = _run(tmp_path, {"ANA_POL": 5, "ANA_CLUSTER": 0, "ANA_ACCEPTANCE": 0,
                               "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0,
                               "EN_FREQ": 10})

    for name in _POL_FILES:
        assert name in produced, name
        assert (tmp_path / name).stat().st_size > 0, name
    assert "ENERGY.dat" in produced


# ---------------------------------------------------------------------------
# 4. an analysis with nothing to measure produces no file
# ---------------------------------------------------------------------------

def test_inter_residue_analysis_without_pairs_writes_no_file(tmp_path):
    """ANA_INTER_RESIDUE with no ANA_RESIDUE_PAIRS has nothing to measure, so
    RES_TO_RES_DIST.dat used to be created and then left empty for the whole
    run."""
    produced = _run(tmp_path, {"ANA_INTER_RESIDUE": 5, "ANA_CLUSTER": 0, "ANA_POL": 0,
                               "ANA_ACCEPTANCE": 0, "ANA_INTSCAL": 0, "ANA_DISTMAP": 0})

    assert "RES_TO_RES_DIST.dat" not in produced


def test_inter_residue_analysis_with_a_pair_writes_its_file(tmp_path):
    produced = _run(tmp_path, {"ANA_INTER_RESIDUE": 5, "ANA_RESIDUE_PAIRS": "0 4",
                               "ANA_CLUSTER": 0, "ANA_POL": 0, "ANA_ACCEPTANCE": 0,
                               "ANA_INTSCAL": 0, "ANA_DISTMAP": 0})

    assert "RES_TO_RES_DIST.dat" in produced
    assert (tmp_path / "RES_TO_RES_DIST.dat").stat().st_size > 0


# ---------------------------------------------------------------------------
# 5. a single-chain-type run writes no CHAIN_<type>_ files
# ---------------------------------------------------------------------------

def test_single_chain_type_run_writes_no_per_chain_type_files(tmp_path):
    """The per-type cluster composition files are meaningless with one type, and
    the per-type internal-scaling files are only used by multi-type runs, so a
    single-type run must produce none of them."""
    produced = _run(tmp_path, {"ANA_CLUSTER": 5, "ANA_INTSCAL": 5, "ANA_DISTMAP": 5,
                               "ANA_POL": 5, "ANA_ACCEPTANCE": 5})

    assert not [name for name in produced if name.startswith("CHAIN_")]
    # the unprefixed versions are what a single-type run writes
    assert "CLUSTERS.dat" in produced
    assert "INTSCAL.dat" in produced


def test_multi_chain_type_run_writes_the_per_chain_type_files(tmp_path):
    produced = _run(tmp_path, {"ANA_CLUSTER": 5, "ANA_INTSCAL": 5, "ANA_DISTMAP": 5},
                    chains=[(2, "AABBA"), (2, "AAAA")])

    for name in ("CHAIN_0_CLUSTERS.dat", "CHAIN_1_CLUSTERS.dat",
                 "CHAIN_0_LR_CLUSTERS.dat", "CHAIN_1_LR_CLUSTERS.dat",
                 "CHAIN_0_INTSCAL.dat", "CHAIN_1_INTSCAL.dat",
                 "CHAIN_0_DISTANCE_MAP.dat", "CHAIN_1_DISTANCE_MAP.dat"):
        assert name in produced, name

    # per-type internal scaling REPLACES the unprefixed files
    assert "INTSCAL.dat" not in produced
    assert "DISTANCE_MAP.dat" not in produced


# ---------------------------------------------------------------------------
# 6. the stale-data guarantee
# ---------------------------------------------------------------------------

def test_rerun_inherits_nothing_from_a_previous_run(tmp_path):
    """The most important property in this file.

    Lazy creation means a re-run no longer overwrites (or truncates) the files
    it does not write, so start-up has to delete them instead. Here the first
    run turns every analysis on and the second turns them off; afterwards not
    one byte of the first run may survive, including in the manufactured files
    that NEITHER run writes.
    """
    first = _run(tmp_path, {"ANA_CLUSTER": 5, "ANA_POL": 5, "ANA_ACCEPTANCE": 5,
                            "ANA_INTSCAL": 5, "ANA_DISTMAP": 5, "EN_FREQ": 10},
                 chains=[(2, "AABBA"), (2, "AAAA")])

    # everything the first run actually produced, plus outputs that neither run
    # writes: a per-type file for a chain type the second run does not have, the
    # conditional and never-triggered outputs, and the quench file
    inherited = sorted(name for name in first if name.endswith(".dat"))
    assert len(inherited) > 10          # the first run really did write a lot

    manufactured = ["CHAIN_7_CLUSTERS.dat", "CHAIN_7_INTSCAL.dat",
                    "CLUSTER_RADIAL_DENSITY_PROFILE.dat",
                    "LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat",
                    "RES_TO_RES_DIST.dat", "QUENCH.dat",
                    "CONFIG_AT_ENERGY_FAIL.pdb"]
    for name in manufactured:
        (tmp_path / name).write_text("data from the previous run\n")

    second = _run(tmp_path, {"ANA_CLUSTER": 0, "ANA_POL": 0, "ANA_ACCEPTANCE": 0,
                             "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0,
                             "EN_FREQ": 10})

    for name in manufactured:
        assert name not in second, name

    # nothing the second run left behind may contain the first run's data: every
    # .dat file present was written by the second run itself
    for name in second:
        if not name.endswith(".dat"):
            continue
        assert "previous run" not in (tmp_path / name).read_text(), name

    # and the second run's analyses really were off, so those files are simply
    # gone rather than re-truncated
    for name in _CLUSTER_FILES + _POL_FILES + ("CHAIN_0_CLUSTERS.dat", "INTSCAL.dat"):
        assert name not in second, name


def test_rerun_keeps_the_restart_file_from_the_previous_run(tmp_path):
    """restart.pimms is the documented exception: it is overwritten at the first
    checkpoint, not deleted at start-up, so a run that dies early leaves the
    previous segment's restart file usable."""
    _run(tmp_path, {"ANA_CLUSTER": 0, "ANA_POL": 0, "ANA_ACCEPTANCE": 0,
                    "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0})
    assert (tmp_path / "restart.pimms").exists()

    marker = (tmp_path / "restart.pimms").read_bytes()
    (tmp_path / "restart.pimms").write_bytes(marker)

    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=_BOX, chains=_CHAINS, n_steps=30, equilibration=2,
                          extra=dict(_QUIET))
    os.chdir(tmp_path)
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state.sim.startup_analysis()

    assert (tmp_path / "restart.pimms").exists()


# ---------------------------------------------------------------------------
# 7. PERFORMANCE.dat carries exactly one header, at the top
# ---------------------------------------------------------------------------

def test_performance_file_has_exactly_one_header_line_at_the_top(tmp_path):
    """The header used to be written by the start-up wipe. It now travels with
    whichever row first creates the file, so it must appear once and first."""
    produced = _run(tmp_path, {"ANA_CLUSTER": 0, "ANA_POL": 0, "ANA_ACCEPTANCE": 0,
                               "ANA_INTSCAL": 0, "ANA_DISTMAP": 0, "ANA_INTER_RESIDUE": 0},
                    n_steps=60)

    assert "PERFORMANCE.dat" in produced
    lines = (tmp_path / "PERFORMANCE.dat").read_text().splitlines()

    assert lines[0].startswith("Step\tE or P\t")
    assert len([line for line in lines if line.startswith("Step\t")]) == 1
    assert len(lines) > 1                       # there are data rows under it
    for line in lines[1:]:
        assert line.split("\t")[0].isdigit()
