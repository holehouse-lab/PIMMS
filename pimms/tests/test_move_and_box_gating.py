"""
Tests for keyfile-level gating of moves and box shapes:

* non-cubic boxes are no longer behind EXPERIMENTAL_FEATURES,
* the pull / TSMMC / jump-and-relax moves are no longer experimental (only VMMC is),
* cluster ROTATION is rejected on a non-cubic box *under periodic boundaries* only
  (it is fine on a cube/square, and fine on a non-cubic box under hardwall).
"""

import os
import io
import contextlib

import pytest

from pimms.tests import kernel_test_utils as U
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import KeyFileException


def _parse(tmp_path, dims, hardwall, moves, *, experimental=False, chains="CHAIN : 4 AABB", extra=None):
    """Run the FULL keyfile parse+validation (not parse_only) and return the parser."""
    d = str(tmp_path)
    U.write_param_file(os.path.join(d, "params.prm"), "SR")
    lines = [
        f"DIMENSIONS : {dims}",
        "PARAMETER_FILE : params.prm",
        "TEMPERATURE : 45",
        "N_STEPS : 100",
        "EQUILIBRATION : 10",
        chains,
        f"HARDWALL : {'True' if hardwall else 'False'}",
    ]
    if experimental:
        lines.append("EXPERIMENTAL_FEATURES : True")
    for kw, val in (extra or {}).items():
        lines.append(f"{kw} : {val}")
    for kw, frac in moves.items():
        lines.append(f"{kw} : {frac}")
    with open(os.path.join(d, "KEYFILE.kf"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return KeyFileParser("KEYFILE.kf")
    finally:
        os.chdir(cwd)


# --------------------------------------------------------------------------- #
# non-cubic boxes no longer require EXPERIMENTAL_FEATURES
# --------------------------------------------------------------------------- #
def test_noncubic_pbc_parses_without_experimental(tmp_path):
    p = _parse(tmp_path, "8 16 32", hardwall=False, moves={"MOVE_CRANKSHAFT": 1.0})
    assert p.keyword_lookup["DIMENSIONS"] == [8, 16, 32]


def test_noncubic_hardwall_parses_without_experimental(tmp_path):
    p = _parse(tmp_path, "8 16 32", hardwall=True, moves={"MOVE_CRANKSHAFT": 1.0})
    assert p.keyword_lookup["DIMENSIONS"] == [8, 16, 32]


# --------------------------------------------------------------------------- #
# cluster rotation is PBC-specific on non-cubic boxes
# --------------------------------------------------------------------------- #
def test_noncubic_pbc_cluster_rotate_raises(tmp_path):
    with pytest.raises(KeyFileException, match="MOVE_CLUSTER_ROTATE cannot be used with a non-cubic"):
        _parse(tmp_path, "8 16 32", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_CLUSTER_ROTATE": 0.4})


def test_noncubic_hardwall_cluster_rotate_allowed(tmp_path):
    # a 90 degree rotation is a valid isometry when there is no periodic wrapping
    p = _parse(tmp_path, "8 16 32", hardwall=True,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_CLUSTER_ROTATE": 0.4})
    assert p.keyword_lookup["MOVE_CLUSTER_ROTATE"] == 0.4


def test_cubic_pbc_cluster_rotate_allowed(tmp_path):
    p = _parse(tmp_path, "20 20 20", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_CLUSTER_ROTATE": 0.4})
    assert p.keyword_lookup["MOVE_CLUSTER_ROTATE"] == 0.4


def test_noncubic_resized_equilibration_pbc_cluster_rotate_raises(tmp_path):
    # production box is CUBIC but the resized-equilibration box is non-cubic; cluster
    # rotation runs during equilibration too, so this must be caught (naming the
    # RESIZED_EQUILIBRATION box).
    with pytest.raises(KeyFileException, match="resized-equilibration"):
        _parse(tmp_path, "40 40 40", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_CLUSTER_ROTATE": 0.4},
               extra={"RESIZED_EQUILIBRATION": "20 20 40"})


def test_cubic_resized_equilibration_pbc_cluster_rotate_allowed(tmp_path):
    # both production and resized-equilibration boxes are cubic -> fine under PBC
    p = _parse(tmp_path, "40 40 40", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_CLUSTER_ROTATE": 0.4},
               extra={"RESIZED_EQUILIBRATION": "20 20 20"})
    assert p.keyword_lookup["RESIZED_EQUILIBRATION"] == [20, 20, 20]


# --------------------------------------------------------------------------- #
# graduated moves no longer need EXPERIMENTAL_FEATURES; VMMC still does
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("move", [
    "MOVE_PULL", "MOVE_CTSMMC", "MOVE_MULTICHAIN_TSMMC",
    "MOVE_SYSTEM_TSMMC", "MOVE_JUMP_AND_RELAX",
])
def test_graduated_moves_parse_without_experimental(tmp_path, move):
    p = _parse(tmp_path, "20 20 20", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, move: 0.4})
    assert p.keyword_lookup[move] == 0.4


def test_vmmc_still_requires_experimental(tmp_path):
    with pytest.raises(KeyFileException, match="EXPERIMENTAL_FEATURES"):
        _parse(tmp_path, "20 20 20", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_VMMC": 0.4})


def test_vmmc_parses_with_experimental(tmp_path):
    p = _parse(tmp_path, "20 20 20", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_VMMC": 0.4}, experimental=True)
    assert p.keyword_lookup["MOVE_VMMC"] == 0.4


# --------------------------------------------------------------------------- #
# EQUILIBRATION_OFFSET / FREEZE_FILE / EXTRA_CHAIN are no longer gated
# --------------------------------------------------------------------------- #
def test_equilibration_offset_parses_without_experimental(tmp_path):
    # EQUILIBRATION_OFFSET requires a RESIZED_EQUILIBRATION box; neither is gated.
    p = _parse(tmp_path, "40 40 40", hardwall=True,
               moves={"MOVE_CRANKSHAFT": 1.0},
               extra={"RESIZED_EQUILIBRATION": "20 20 20",
                      "EQUILIBRATION_OFFSET": "5 5 5"})
    assert p.keyword_lookup["EQUILIBRATION_OFFSET"] == [5, 5, 5]


def test_freeze_file_parses_without_experimental(tmp_path):
    # a freeze file that exists and is well-formed parses fine with no gate
    with open(os.path.join(str(tmp_path), "frz.in"), "w") as fh:
        fh.write("C 1 2\n")
    p = _parse(tmp_path, "20 20 20", hardwall=False,
               moves={"MOVE_CRANKSHAFT": 1.0},
               extra={"FREEZE_FILE": "frz.in"})
    # FREEZE_FILE is replaced by a parsed FreezeFile object during validation
    assert p.keyword_lookup["FREEZE_FILE"] is not False
