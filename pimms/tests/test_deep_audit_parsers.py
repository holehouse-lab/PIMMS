"""
Regression tests for the input-parsing findings of the deep audit.

They cover the keyfile parser, the parameter-file parser, the freeze-file
parser, the keyword limits in ``CONFIG`` and the ``PIMMS`` command line:

* K1-1 / K1-2 / X1-11 - a file that is not UTF-8 is refused by name, and a UTF-8
  byte-order mark is dropped;
* K1-4 / S1-7 - the sub-step and TSMMC counts that size an array have an upper
  limit, applied when the keyfile is read;
* K1-5 - a box whose lattice could not be allocated is refused;
* K1-7 / O1-2 / O1-7 - ``LATTICE_TO_ANGSTROMS`` has a validated range;
* K1-8 - keywords that would be silently ignored are warned about, and are still
  validated;
* K1-9 - command-line combinations that used to skip the run silently;
* K1-10 / E1-2 - only plain ASCII number syntax is a number;
* E1-1 - a temperature-normalised angle penalty that overflows is a clean error;
* M3-3 / X1-9 - ``PARALLEL_THREADS`` has a ceiling and ``0`` respects
  ``OMP_NUM_THREADS`` and the CPU affinity.

Every test works at parse level (nothing is simulated), writes only under
``tmp_path`` and derives its expected values from the definitions rather than
from the routine under test.
"""

import builtins
import contextlib
import importlib.machinery
import importlib.util
import io
import json
import os
import pathlib
import stat
import subprocess
import sys
import types
from typing import Dict, Iterable, List, Optional, Tuple

import pytest

from pimms import CONFIG
from pimms import data_structures
from pimms import keyfile_parser
from pimms import parameterfile_parser
from pimms.latticeExceptions import KeyFileException, ParameterFileException

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
PIMMS_SCRIPT = REPO_ROOT / "scripts" / "PIMMS"

PARAMETER_TEXT = (
    "A 0 0\nB 0 1\nA A -3 -1 0\nA B -1\nB B -2 -1 -1\n"
    "ANGLE_PENALTY A 0 1 2\nANGLE_PENALTY B 0 0 0\n"
)

# 4 chains of 5 beads in a 12 x 12 x 12 box, 40 steps: small, valid and silent
BASE: Dict[str, object] = {
    "DIMENSIONS": "12 12 12",
    "PARAMETER_FILE": "params.prm",
    "TEMPERATURE": "10",
    "N_STEPS": "40",
    "EQUILIBRATION": "0",
    "CHAIN": ["4 AABBA"],
    "SEED": "11",
    "MOVE_CRANKSHAFT": "0.5",
    "MOVE_CHAIN_TRANSLATE": "0.5",
    "XTC_FREQ": "5",
    "EN_FREQ": "5",
}

NOT_A_NUMBER = [
    ("N_STEPS", "1_000"),  # digit-group underscore: int() reads 1000
    ("N_STEPS", "４０"),  # full-width 40
    ("N_STEPS", "٤٠"),  # Arabic-Indic 40
    ("N_STEPS", " 40"),  # no-break space before
    ("N_STEPS", "40 "),  # no-break space after
    ("TEMPERATURE", "1_0.5"),
    ("TEMPERATURE", "１０"),
    ("TEMPERATURE", "10 "),
    ("DIMENSIONS", "1_2 12 12"),
    ("DIMENSIONS", "12 12 12"),
    ("CHAIN", ["４ AABBA"]),
    ("ANA_RESIDUE_PAIRS", ["0 ２"]),
]

INT32_MAX = 2**31 - 1

IGNORED = [
    (
        {
            "QUENCH_START": "10",
            "QUENCH_END": "5",
            "QUENCH_STEPSIZE": "1",
            "QUENCH_FREQ": "2",
            "QUENCH_AS_EQUILIBRATION": "True",
        },
        ["QUENCH_START", "QUENCH_FREQ", "QUENCH_RUN"],
    ),
    (
        {"TSMMC_JUMP_TEMP": "60", "TSMMC_NUMBER_OF_POINTS": "7"},
        ["TSMMC_JUMP_TEMP", "TSMMC_NUMBER_OF_POINTS"],
    ),
    (
        {"RESTART_OVERRIDE_DIMENSIONS": "True"},
        ["RESTART_OVERRIDE_DIMENSIONS", "RESTART_FILE"],
    ),
    (
        {"RESTART_OVERRIDE_HARDWALL": "True"},
        ["RESTART_OVERRIDE_HARDWALL", "RESTART_FILE"],
    ),
    ({"PARALLEL_THREADS": "3"}, ["PARALLEL_THREADS", "PARALLELIZE"]),
    ({"SLITHER_SUBSTEPS": "50"}, ["SLITHER_SUBSTEPS", "MOVE_SLITHER"]),
    ({"PULL_SUBSTEPS": "50"}, ["PULL_SUBSTEPS", "MOVE_PULL"]),
    (
        {
            "CRANKSHAFT_SUBSTEPS": "50",
            "MOVE_CRANKSHAFT": None,
            "MOVE_CHAIN_TRANSLATE": "1.0",
        },
        ["CRANKSHAFT_SUBSTEPS", "MOVE_CRANKSHAFT"],
    ),
    ({"XTC_FREQ": "1000"}, ["XTC_FREQ", "N_STEPS", "traj.xtc"]),
    ({"EN_FREQ": "41"}, ["EN_FREQ", "N_STEPS", "ENERGY.dat"]),
    ({"ANA_RESIDUE_PAIRS": ["2 2"]}, ["ANA_RESIDUE_PAIRS", "2 2"]),
    ({"ANA_RESIDUE_PAIRS": ["1 3", "3 1"]}, ["ANA_RESIDUE_PAIRS", "1 3"]),
    ({"ANA_CLUSTER_THRESHOLD": "4"}, ["ANA_CLUSTER_THRESHOLD", "4"]),
]

NOT_IGNORED = [
    # the feature is switched off in so many words: a template, not a trap
    {
        "QUENCH_RUN": "False",
        "QUENCH_START": "10",
        "QUENCH_END": "5",
        "QUENCH_STEPSIZE": "1",
        "QUENCH_FREQ": "2",
        "QUENCH_AS_EQUILIBRATION": "True",
    },
    {"MOVE_CTSMMC": "0", "TSMMC_JUMP_TEMP": "60"},
    {"PARALLELIZE": "False", "PARALLEL_THREADS": "3"},
    {"MOVE_SLITHER": "0.0", "SLITHER_SUBSTEPS": "50"},
    # written, but with the default value: what keyfile_used.kf does for every keyword
    {
        "TSMMC_JUMP_TEMP": "50.0",
        "SLITHER_SUBSTEPS": "10",
        "PULL_SUBSTEPS": "10",
        "PARALLEL_THREADS": "0",
        "RESTART_OVERRIDE_DIMENSIONS": "False",
        "ANA_CLUSTER_THRESHOLD": "1",
    },
    # in use
    {"MOVE_SLITHER": "0.5", "MOVE_CHAIN_TRANSLATE": "0", "SLITHER_SUBSTEPS": "50"},
    {"PARALLELIZE": "True", "PARALLEL_THREADS": "3"},
    {"MOVE_CTSMMC": "0.5", "MOVE_CHAIN_TRANSLATE": "0", "TSMMC_JUMP_TEMP": "60"},
    # jump-and-relax uses CRANKSHAFT_SUBSTEPS for its relaxations
    {
        "CRANKSHAFT_SUBSTEPS": "50",
        "MOVE_CRANKSHAFT": None,
        "MOVE_CHAIN_TRANSLATE": "0.5",
        "MOVE_JUMP_AND_RELAX": "0.5",
    },
    # a threshold that a cluster of all four chains exceeds, frequencies that fire, distinct pairs
    {
        "ANA_CLUSTER_THRESHOLD": "3",
        "XTC_FREQ": "40",
        "EN_FREQ": "40",
        "ANA_RESIDUE_PAIRS": ["1 3", "0 4"],
    },
]

# bead type and comments that neither ASCII nor (for the epsilon) Latin-1 can hold
NON_ASCII_PARAMETERS = (
    "# \u00c5 and \u03b5\n"
    "A 0 0\n\u00c5 0 1\nA A -3\nA \u00c5 -1\n\u00c5 \u00c5 -2\n"
    "ANGLE_PENALTY A 0 1 2\nANGLE_PENALTY \u00c5 0 0 0\n"
)

LOCALE_CHILD = r"""
import contextlib, io, json, locale, shutil
from pimms import keyfile_parser, parameterfile_parser, pimmslogger

with contextlib.redirect_stdout(io.StringIO()):
    parser = keyfile_parser.KeyFileParser("KEYFILE.kf")
    parser.write_keyfile("written.kf")
    again = keyfile_parser.KeyFileParser("written.kf")
    pimmslogger.log_status("freeze file \u03b5.txt")
    first = parameterfile_parser.parse_energy("params.prm")
    angles = parameterfile_parser.parse_angles("params.prm")
    parameterfile_parser.write_angle_parameter_summary(angles, "params.prm")
    shutil.copy("parameters_used.prm", "echo.prm")
    second = parameterfile_parser.parse_energy("echo.prm")
print(json.dumps({
    "encoding": locale.getpreferredencoding(False),
    "chains": [parser.keyword_lookup["CHAIN"], again.keyword_lookup["CHAIN"]],
    "tables_equal": first == second,
    "residues": first[1],
}))
"""

NBSP = "\u00a0"


def _keyfile_text(overrides: Optional[Dict[str, object]] = None) -> str:
    """
    Build the text of a keyfile from ``BASE`` plus overrides.

    Parameters
    ----------
    overrides : dict, optional
        Keyword to value. A list gives one line per entry (``CHAIN``,
        ``ANA_RESIDUE_PAIRS``) and None removes the keyword.

    Returns
    -------
    str
        The keyfile text.
    """
    values = dict(BASE)
    values.update(overrides or {})
    lines: List[str] = []
    for keyword, value in values.items():
        if value is None:
            continue
        for entry in value if isinstance(value, list) else [value]:
            lines.append("%s : %s" % (keyword, entry))
    return "\n".join(lines) + "\n"


def _parse(
    tmp_path: pathlib.Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture,
    overrides: Optional[Dict[str, object]] = None,
) -> Tuple[keyfile_parser.KeyFileParser, str]:
    """
    Write a keyfile and its parameter file into ``tmp_path`` and parse it there.

    The parser starts ``log.txt`` in the working directory, so the working
    directory is moved to ``tmp_path`` for the duration of the test.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The test's private directory.
    monkeypatch : pytest.MonkeyPatch
        Used to change the working directory.
    capsys : pytest.CaptureFixture
        Used to collect what the parser printed.
    overrides : dict, optional
        Changes to ``BASE`` (see :func:`_keyfile_text`).

    Returns
    -------
    tuple
        The parser and everything it printed.
    """
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    (tmp_path / "KEYFILE.kf").write_text(_keyfile_text(overrides), encoding="utf-8")
    capsys.readouterr()
    parser = keyfile_parser.KeyFileParser("KEYFILE.kf")
    return parser, capsys.readouterr().out


def _warnings(output: str) -> List[str]:
    """
    Pick the warning lines out of the parser's output.

    Parameters
    ----------
    output : str
        What the parser printed.

    Returns
    -------
    list of str
        The lines that carry a warning tag.
    """
    return [line for line in output.splitlines() if "WARNING" in line.upper()]


def _plenty_of_memory(
    monkeypatch: pytest.MonkeyPatch, n_bytes: Optional[int] = 2**50
) -> None:
    """
    Make the parser believe the machine has ``n_bytes`` of physical memory.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Used to replace ``keyfile_parser.physical_memory_bytes``.
    n_bytes : int or None, optional
        The memory to report (None means "cannot be read"). Default 1 PiB, so
        that only the limits that do not depend on the machine can fire.
    """
    monkeypatch.setattr(
        keyfile_parser, "physical_memory_bytes", lambda: n_bytes, raising=False
    )


# ---------------------------------------------------------------------------
# the base keyfile itself
# ---------------------------------------------------------------------------


def test_the_base_keyfile_is_valid_and_silent(tmp_path, monkeypatch, capsys):
    """Every other test changes one thing about this keyfile, so it must be clean."""
    parser, out = _parse(tmp_path, monkeypatch, capsys)
    assert parser.keyword_lookup["DIMENSIONS"] == [12, 12, 12]
    assert parser.keyword_lookup["N_STEPS"] == 40
    assert _warnings(out) == []


# ---------------------------------------------------------------------------
# K1-1 / K1-2 / X1-11: encoding
# ---------------------------------------------------------------------------


def test_a_latin1_keyfile_is_refused_by_name_and_offset(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    raw = ("LATTICE_TO_ANGSTROMS : 3.65  # Å per site\n" + _keyfile_text()).encode(
        "latin-1"
    )
    (tmp_path / "latin1.kf").write_bytes(raw)
    offset = raw.index(b"\xc5")
    with pytest.raises(KeyFileException) as error:
        keyfile_parser.KeyFileParser("latin1.kf")
    message = str(error.value)
    assert "latin1.kf" in message
    assert "byte offset %d" % offset in message
    assert "UTF-8" in message


def test_a_utf16_keyfile_is_refused_and_says_so(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    (tmp_path / "utf16.kf").write_bytes(_keyfile_text().encode("utf-16"))
    with pytest.raises(KeyFileException, match="UTF-16"):
        keyfile_parser.KeyFileParser("utf16.kf")


def test_a_byte_order_mark_does_not_change_how_a_keyfile_parses(
    tmp_path, monkeypatch, capsys
):
    plain, _ = _parse(tmp_path, monkeypatch, capsys)
    (tmp_path / "bom.kf").write_bytes(b"\xef\xbb\xbf" + _keyfile_text().encode("utf-8"))
    with_bom = keyfile_parser.KeyFileParser("bom.kf")
    expected = {k: v for k, v in plain.keyword_lookup.items() if k != "__KEYFILE"}
    assert {
        k: v for k, v in with_bom.keyword_lookup.items() if k != "__KEYFILE"
    } == expected


def test_the_parameter_file_is_read_as_utf8_with_or_without_a_byte_order_mark(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "bom.prm").write_bytes(b"\xef\xbb\xbf" + PARAMETER_TEXT.encode("utf-8"))
    pairs, residues, _, _, _ = parameterfile_parser.parse_energy("bom.prm")
    # the first line is "A 0 0": with the mark glued on, the residue was '﻿A'
    assert residues == ["0", "A", "B"]
    assert pairs["A"]["0"] == 0
    assert set(parameterfile_parser.parse_angles("bom.prm")) == {"A", "B"}

    raw = ("# energies in kT (µ = 0)\n" + PARAMETER_TEXT).encode("latin-1")
    (tmp_path / "latin1.prm").write_bytes(raw)
    for reader in (
        parameterfile_parser.parse_energy,
        parameterfile_parser.parse_angles,
    ):
        with pytest.raises(ParameterFileException) as error:
            reader("latin1.prm")
        assert "latin1.prm" in str(error.value)
        assert "byte offset %d" % raw.index(b"\xb5") in str(error.value)


def test_the_freeze_file_is_read_as_utf8_with_or_without_a_byte_order_mark(tmp_path):
    bom = tmp_path / "bom.frz"
    bom.write_bytes(b"\xef\xbb\xbfC 1 3\n")
    assert data_structures.FreezeFile(str(bom)).chains == [1, 3]

    raw = "# chaîne 1\nC 1\n".encode("latin-1")
    latin1 = tmp_path / "latin1.frz"
    latin1.write_bytes(raw)
    with pytest.raises(KeyFileException) as error:
        data_structures.FreezeFile(str(latin1))
    assert "latin1.frz" in str(error.value)
    assert "byte offset %d" % raw.index(b"\xee") in str(error.value)


# ---------------------------------------------------------------------------
# K1-10 / E1-2: plain ASCII numbers only
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("keyword,value", NOT_A_NUMBER)
def test_a_number_that_is_not_plain_ascii_is_refused(
    tmp_path, monkeypatch, capsys, keyword, value
):
    with pytest.raises(KeyFileException) as error:
        _parse(tmp_path, monkeypatch, capsys, {keyword: value})
    assert keyword in str(error.value)


def test_plain_ascii_number_syntax_is_still_accepted(tmp_path, monkeypatch, capsys):
    parser, _ = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {
            "N_STEPS": "+40",
            "TEMPERATURE": "1.5e1",
            "DIMENSIONS": "12\t12  12",
            "LATTICE_TO_ANGSTROMS": "4.",
            "ANA_RESIDUE_PAIRS": ["3 1"],
        },
    )
    lookup = parser.keyword_lookup
    assert lookup["N_STEPS"] == 40
    assert lookup["TEMPERATURE"] == 15.0
    assert lookup["DIMENSIONS"] == [12, 12, 12]
    assert lookup["LATTICE_TO_ANGSTROMS"] == 4.0
    assert lookup["ANA_RESIDUE_PAIRS"] == [[1, 3]]


@pytest.mark.parametrize(
    "token,expected",
    [
        ("1_000", False),
        ("１２", False),
        ("12 ", False),
        ("1.0", False),
        ("", False),
        ("+7", True),
        ("-3", True),
        ("0012", True),
    ],
)
def test_is_ascii_integer(token, expected):
    assert data_structures.is_ascii_integer(token) is expected


@pytest.mark.parametrize(
    "token,expected",
    [
        ("1_0.5", False),
        ("nan", False),
        ("inf", False),
        ("1e", False),
        (".", False),
        ("0x10", False),
        ("10", True),
        ("10.", True),
        (".5", True),
        ("-2.5E-3", True),
    ],
)
def test_is_ascii_float(token, expected):
    assert data_structures.is_ascii_float(token) is expected


@pytest.mark.parametrize("token", ["1_000", "１２"])
def test_the_parameter_file_refuses_numbers_that_are_not_plain_ascii(
    tmp_path, monkeypatch, token
):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "p.prm").write_text(
        "A 0 0\nA A %s\nANGLE_PENALTY A 0 0 0\n" % token, encoding="utf-8"
    )
    with pytest.raises(ParameterFileException, match="Unable to parse line"):
        parameterfile_parser.parse_energy("p.prm")
    (tmp_path / "a.prm").write_text(
        "A 0 0\nA A -1\nANGLE_PENALTY A 0 %s 0\n" % token, encoding="utf-8"
    )
    with pytest.raises(ParameterFileException):
        parameterfile_parser.parse_angles("a.prm")
    (tmp_path / "t.prm").write_text(
        "A 0 0\nA A -1\nANGLE_PENALTY_T_NORM A %s.5 0 0\n" % token, encoding="utf-8"
    )
    with pytest.raises(ParameterFileException):
        parameterfile_parser.parse_angles("t.prm", temperature=1.0)


def test_the_parameter_file_still_reads_signed_integers(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "p.prm").write_text("A 0 +7\nA A -3\nANGLE_PENALTY A 0 0 0\n")
    pairs = parameterfile_parser.parse_energy("p.prm")[0]
    assert pairs["A"]["0"] == 7
    assert pairs["A"]["A"] == -3


@pytest.mark.parametrize("line", ["C 1_0", "C １"])
def test_the_freeze_file_refuses_chain_ids_that_are_not_plain_ascii(tmp_path, line):
    path = tmp_path / "f.frz"
    path.write_text(line + "\n", encoding="utf-8")
    with pytest.raises(KeyFileException, match="Error parsing chains"):
        data_structures.FreezeFile(str(path))


# ---------------------------------------------------------------------------
# E1-1: a T-normalised penalty whose product with the temperature is infinite
# ---------------------------------------------------------------------------


def test_an_overflowing_t_norm_penalty_is_a_parameter_file_exception(tmp_path):
    path = tmp_path / "t.prm"
    path.write_text("ANGLE_PENALTY_T_NORM A 3 0 0\nA 0 0\nA A -1\n")
    # 3 x 1e308 is inf although both factors are finite
    assert 3 * 1e308 == float("inf")
    with pytest.raises(ParameterFileException, match="outside the supported"):
        parameterfile_parser.parse_angles(str(path), temperature=1e308)
    # a product that fits is still scaled as documented
    assert parameterfile_parser.parse_angles(str(path), temperature=2.5)["A"] == [
        7.5,
        0.0,
        0.0,
    ]


# ---------------------------------------------------------------------------
# K1-4 / S1-7: counts that size an array
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "keyword",
    [
        "CRANKSHAFT_SUBSTEPS",
        "SLITHER_SUBSTEPS",
        "PULL_SUBSTEPS",
        "TSMMC_STEP_MULTIPLIER",
    ],
)
def test_a_substep_count_beyond_int32_is_refused_whether_or_not_the_move_is_on(
    tmp_path, monkeypatch, capsys, keyword
):
    _plenty_of_memory(monkeypatch)
    for value in (INT32_MAX + 1, 99999999999, 2**64):
        with pytest.raises(KeyFileException) as error:
            _parse(tmp_path, monkeypatch, capsys, {keyword: str(value)})
        assert keyword in str(error.value)
        assert str(INT32_MAX) in str(error.value)
    assert CONFIG.MAX_SUBMOVE_SELECTOR_LENGTH == INT32_MAX


def test_a_substep_count_at_the_limit_is_accepted_when_memory_allows(
    tmp_path, monkeypatch, capsys
):
    _plenty_of_memory(monkeypatch)
    # slither and pull are off, so nothing is multiplied by the chain count
    parser, out = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {
            "SLITHER_SUBSTEPS": str(INT32_MAX),
            "PULL_SUBSTEPS": str(INT32_MAX),
            "MOVE_SLITHER": "0",
            "MOVE_PULL": "0",
        },
    )
    assert parser.keyword_lookup["SLITHER_SUBSTEPS"] == INT32_MAX
    assert _warnings(out) == []


@pytest.mark.parametrize(
    "keyword,moves,multiplier",
    [
        # every chain
        ("SLITHER_SUBSTEPS", {"MOVE_SLITHER": "0.5"}, 8),
        # the chains of three or more beads: 2 x AABBA and 1 x AAAAAAA
        ("PULL_SUBSTEPS", {"MOVE_PULL": "0.5"}, 3),
        # the longest chain
        ("TSMMC_STEP_MULTIPLIER", {"MOVE_CTSMMC": "0.5", "TSMMC_JUMP_TEMP": "50"}, 7),
        # at most floor(8 / 4) + 1 = 3 chains, at most 7 + 5 + 5 beads
        (
            "TSMMC_STEP_MULTIPLIER",
            {"MOVE_MULTICHAIN_TSMMC": "0.5", "TSMMC_JUMP_TEMP": "50"},
            17,
        ),
    ],
)
def test_the_selector_limit_is_applied_to_the_array_the_move_builds(
    tmp_path, monkeypatch, capsys, keyword, moves, multiplier
):
    """
    2 x AABBA, 5 x AA, 1 x AAAAAAA: 8 chains, 3 of them pullable, longest 7.

    The largest value whose product with the move's multiplier still fits a C int
    is accepted and the next one is refused, which pins the multiplier exactly.
    """
    _plenty_of_memory(monkeypatch)
    system = {
        "CHAIN": ["2 AABBA", "5 AA", "1 AAAAAAA"],
        "MOVE_CRANKSHAFT": "0.5",
        "MOVE_CHAIN_TRANSLATE": "0",
    }
    system.update(moves)
    largest = INT32_MAX // multiplier
    assert largest * multiplier <= INT32_MAX < (largest + 1) * multiplier

    parser, _ = _parse(tmp_path, monkeypatch, capsys, {**system, keyword: str(largest)})
    assert parser.keyword_lookup[keyword] == largest
    with pytest.raises(KeyFileException) as error:
        _parse(tmp_path, monkeypatch, capsys, {**system, keyword: str(largest + 1)})
    assert keyword in str(error.value)
    assert str((largest + 1) * multiplier) in str(error.value)

    # and the helper the check is built on reports the same product
    lengths = [n for k, _, n in parser._submove_selector_lengths() if k == keyword]
    assert lengths == [largest * multiplier]


@pytest.mark.parametrize(
    "overrides,keyword,length",
    [
        (
            {
                "SLITHER_SUBSTEPS": str(10**9),
                "MOVE_SLITHER": "0.5",
                "MOVE_CHAIN_TRANSLATE": "0",
            },
            "SLITHER_SUBSTEPS",
            4 * 10**9,
        ),
        (
            {
                "PULL_SUBSTEPS": str(10**9),
                "MOVE_PULL": "0.5",
                "MOVE_CHAIN_TRANSLATE": "0",
            },
            "PULL_SUBSTEPS",
            4 * 10**9,
        ),
        (
            {
                "TSMMC_STEP_MULTIPLIER": str(10**9),
                "MOVE_CTSMMC": "0.5",
                "MOVE_CHAIN_TRANSLATE": "0",
                "TSMMC_JUMP_TEMP": "50",
            },
            "TSMMC_STEP_MULTIPLIER",
            5 * 10**9,
        ),
    ],
)
def test_a_selector_beyond_int32_is_refused(
    tmp_path, monkeypatch, capsys, overrides, keyword, length
):
    """Each value is below 2^31 on its own; times 4 chains (or 5 beads) it is not."""
    _plenty_of_memory(monkeypatch)
    with pytest.raises(KeyFileException) as error:
        _parse(tmp_path, monkeypatch, capsys, overrides)
    assert keyword in str(error.value)
    assert str(length) in str(error.value)


def test_a_selector_larger_than_physical_memory_is_refused(
    tmp_path, monkeypatch, capsys
):
    # 2 x 10^8 sub-moves at 8 bytes each is 1.6 GB: more than a 1 GB machine has
    _plenty_of_memory(monkeypatch, 10**9)
    with pytest.raises(KeyFileException, match="physical memory") as error:
        _parse(tmp_path, monkeypatch, capsys, {"CRANKSHAFT_SUBSTEPS": str(2 * 10**8)})
    assert "CRANKSHAFT_SUBSTEPS" in str(error.value)
    # 10^8 sub-moves (0.8 GB) fit the same machine: 8 bytes per entry, not more
    parser, _ = _parse(
        tmp_path, monkeypatch, capsys, {"CRANKSHAFT_SUBSTEPS": str(10**8)}
    )
    assert parser.keyword_lookup["CRANKSHAFT_SUBSTEPS"] == 10**8
    assert CONFIG.SUBMOVE_SELECTOR_BYTES_PER_ENTRY == 8


def test_a_large_selector_that_fits_gets_a_warning_and_a_normal_one_does_not(
    tmp_path, monkeypatch, capsys
):
    _plenty_of_memory(monkeypatch)
    parser, out = _parse(
        tmp_path, monkeypatch, capsys, {"CRANKSHAFT_SUBSTEPS": str(10**8 + 1)}
    )
    assert parser.keyword_lookup["CRANKSHAFT_SUBSTEPS"] == 10**8 + 1
    assert len(_warnings(out)) == 1 and "CRANKSHAFT_SUBSTEPS" in _warnings(out)[0]
    # where physical memory cannot be read only the int32 limit refuses
    _plenty_of_memory(monkeypatch, None)
    _, out = _parse(
        tmp_path, monkeypatch, capsys, {"CRANKSHAFT_SUBSTEPS": str(10**8 + 1)}
    )
    assert len(_warnings(out)) == 1
    _, out = _parse(tmp_path, monkeypatch, capsys, {"CRANKSHAFT_SUBSTEPS": str(10**8)})
    assert _warnings(out) == []
    assert CONFIG.SUBMOVE_SELECTOR_WARN_LENGTH == 10**8


def test_crankshaft_substeps_count_when_the_crankshaft_is_only_a_tsmmc_fallback(
    tmp_path, monkeypatch, capsys
):
    """
    With MOVE_SYSTEM_TSMMC as the only move the excursion falls back to crankshaft
    megamoves, so CRANKSHAFT_SUBSTEPS is in use although no crankshaft move is
    given: it must not be reported as ignored, and its selector is still checked.
    """
    fallback = {
        "MOVE_CRANKSHAFT": None,
        "MOVE_CHAIN_TRANSLATE": None,
        "MOVE_SYSTEM_TSMMC": "1.0",
        "TSMMC_JUMP_TEMP": "50",
    }
    _plenty_of_memory(monkeypatch, 10**9)
    with pytest.raises(KeyFileException, match="CRANKSHAFT_SUBSTEPS"):
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {**fallback, "CRANKSHAFT_SUBSTEPS": str(2 * 10**8)},
        )
    _, out = _parse(
        tmp_path, monkeypatch, capsys, {**fallback, "CRANKSHAFT_SUBSTEPS": "50"}
    )
    assert _warnings(out) == []


def test_tsmmc_number_of_points_has_a_ceiling(tmp_path, monkeypatch, capsys):
    with pytest.raises(KeyFileException, match="TSMMC_NUMBER_OF_POINTS"):
        _parse(
            tmp_path, monkeypatch, capsys, {"TSMMC_NUMBER_OF_POINTS": str(10**6 + 1)}
        )
    parser, _ = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {"TSMMC_NUMBER_OF_POINTS": str(10**6), "MOVE_CTSMMC": "0"},
    )
    assert parser.keyword_lookup["TSMMC_NUMBER_OF_POINTS"] == 10**6
    assert CONFIG.MAX_TSMMC_NUMBER_OF_POINTS == 10**6


# ---------------------------------------------------------------------------
# K1-5: lattice memory
# ---------------------------------------------------------------------------


def test_a_box_larger_than_physical_memory_is_refused(tmp_path, monkeypatch, capsys):
    # 3000^3 sites x 8 bytes = 216 GB against a 16 GB machine
    _plenty_of_memory(monkeypatch, 16 * 10**9)
    with pytest.raises(KeyFileException) as error:
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {"DIMENSIONS": "3000 3000 3000", "LATTICE_TO_ANGSTROMS": "1"},
        )
    assert "DIMENSIONS" in str(error.value)
    assert "216 GB" in str(error.value)
    assert str(3000**3) in str(error.value)
    assert CONFIG.LATTICE_BYTES_PER_SITE == 8


@pytest.mark.parametrize("dimensions", ["400 400 400", "2000 2000"])
def test_ordinary_large_boxes_are_accepted_without_a_warning(
    tmp_path, monkeypatch, capsys, dimensions
):
    """512 MB and 32 MB of grids on a machine with 4 GB: nothing to say."""
    _plenty_of_memory(monkeypatch, 4 * 10**9)
    parser, out = _parse(tmp_path, monkeypatch, capsys, {"DIMENSIONS": dimensions})
    assert parser.keyword_lookup["DIMENSIONS"] == [int(d) for d in dimensions.split()]
    assert _warnings(out) == []


def test_a_box_above_half_of_physical_memory_gets_a_warning(
    tmp_path, monkeypatch, capsys
):
    # 1100^3 x 8 bytes = 10.6 GB: between half of 16 GB and 16 GB
    assert 0.5 * 16e9 < 1100**3 * 8 < 16e9
    _plenty_of_memory(monkeypatch, 16 * 10**9)
    parser, out = _parse(
        tmp_path, monkeypatch, capsys, {"DIMENSIONS": "1100 1100 1100"}
    )
    assert parser.keyword_lookup["DIMENSIONS"] == [1100, 1100, 1100]
    assert len(_warnings(out)) == 1 and "DIMENSIONS" in _warnings(out)[0]


def test_the_fixed_ceiling_applies_where_memory_cannot_be_read(
    tmp_path, monkeypatch, capsys
):
    _plenty_of_memory(monkeypatch, None)
    # 3000^3 x 8 = 2.16e11 bytes is below half a TiB; 6000^3 x 8 = 1.7e12 is above one
    assert 3000**3 * 8 < 0.5 * 2**40 and 6000**3 * 8 > 2**40
    _, out = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {"DIMENSIONS": "3000 3000 3000", "LATTICE_TO_ANGSTROMS": "1"},
    )
    assert _warnings(out) == []
    with pytest.raises(KeyFileException, match="DIMENSIONS"):
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {"DIMENSIONS": "6000 6000 6000", "LATTICE_TO_ANGSTROMS": "1"},
        )
    assert CONFIG.LATTICE_MEMORY_FALLBACK_CEILING == 2**40


def test_the_real_memory_probe_refuses_an_impossible_box(tmp_path, monkeypatch, capsys):
    """
    Nothing is monkeypatched here, so this goes through the real probe.

    100000^3 sites x 8 bytes is 8 PB: more than any machine has and more than the
    1 TiB fallback, so it is refused whatever the probe returns. On Linux and
    macOS the probe must return a number: a probe that quietly returned None
    would leave every refusal test above green (they replace it) while real runs
    fell back to the 1 TiB ceiling.
    """
    with pytest.raises(KeyFileException, match="DIMENSIONS"):
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {"DIMENSIONS": "100000 100000 100000", "LATTICE_TO_ANGSTROMS": "0.05"},
        )
    memory = keyfile_parser.physical_memory_bytes()
    if sys.platform.startswith(("linux", "darwin")):
        assert isinstance(memory, int) and memory > 2**20
    else:
        assert memory is None or (isinstance(memory, int) and memory > 2**20)


# ---------------------------------------------------------------------------
# K1-7 / O1-2 / O1-7: LATTICE_TO_ANGSTROMS
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spacing", ["1e-300", "1e-9", "0.004", "0.0099"])
def test_a_spacing_below_the_xtc_resolution_is_refused(
    tmp_path, monkeypatch, capsys, spacing
):
    """XTC stores 0.001 nm = 0.01 Angstroms; 1e-300 used to end in a ZeroDivisionError."""
    with pytest.raises(KeyFileException, match="LATTICE_TO_ANGSTROMS"):
        _parse(tmp_path, monkeypatch, capsys, {"LATTICE_TO_ANGSTROMS": spacing})


def test_the_smallest_spacing_parses_summarises_and_warns(
    tmp_path, monkeypatch, capsys
):
    parser, out = _parse(
        tmp_path, monkeypatch, capsys, {"LATTICE_TO_ANGSTROMS": "0.01"}
    )
    assert parser.keyword_lookup["LATTICE_TO_ANGSTROMS"] == 0.01
    assert len(_warnings(out)) == 1 and "LATTICE_TO_ANGSTROMS" in _warnings(out)[0]
    parser.print_summary()  # divided by a volume that underflowed to zero before
    _, out = _parse(tmp_path, monkeypatch, capsys, {"LATTICE_TO_ANGSTROMS": "0.1"})
    assert _warnings(out) == []


def test_a_box_of_10000_angstroms_or_more_is_refused(tmp_path, monkeypatch, capsys):
    """The %8.3f PDB field ends at 9999.999: 2739 sites of 3.65 fit, 2740 do not."""
    _plenty_of_memory(monkeypatch)
    assert 2739 * 3.65 < 10000 <= 2740 * 3.65
    assert len("%8.3f" % 9999.999) == 8 and len("%8.3f" % 10000.0) == 9
    parser, _ = _parse(tmp_path, monkeypatch, capsys, {"DIMENSIONS": "20 20 2739"})
    assert parser.keyword_lookup["DIMENSIONS"] == [20, 20, 2739]
    for overrides in (
        {"DIMENSIONS": "20 20 2740"},
        {"DIMENSIONS": "3000 3000"},
        {"DIMENSIONS": "10 10 1100", "LATTICE_TO_ANGSTROMS": "10"},
        {"LATTICE_TO_ANGSTROMS": "1000"},
    ):
        with pytest.raises(KeyFileException) as error:
            _parse(tmp_path, monkeypatch, capsys, overrides)
        assert "DIMENSIONS" in str(error.value) and "LATTICE_TO_ANGSTROMS" in str(
            error.value
        )


def test_the_box_limit_is_applied_to_the_production_box_of_a_resized_run(
    tmp_path, monkeypatch, capsys
):
    """This used to crash at the resize, after the whole equilibration had run."""
    _plenty_of_memory(monkeypatch)
    with pytest.raises(KeyFileException, match="LATTICE_TO_ANGSTROMS"):
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {
                "DIMENSIONS": "20 20 3000",
                "RESIZED_EQUILIBRATION": "20 20 20",
                "EQUILIBRATION_OFFSET": "0 0 2900",
                "EQUILIBRATION": "20",
                "HARDWALL": "True",
            },
        )


# ---------------------------------------------------------------------------
# K1-8: keywords that are given and would be ignored
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("overrides,expected", IGNORED)
def test_an_ignored_keyword_gets_one_warning_that_names_it(
    tmp_path, monkeypatch, capsys, overrides, expected
):
    parser, out = _parse(tmp_path, monkeypatch, capsys, overrides)
    warnings = _warnings(out)
    assert len(warnings) == 1, warnings
    for text in expected:
        assert text in warnings[0]
    assert parser.keyword_lookup["N_STEPS"] == 40


def test_a_warning_leaves_the_values_as_written(tmp_path, monkeypatch, capsys):
    parser, out = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {
            "PARALLEL_THREADS": "3",
            "SLITHER_SUBSTEPS": "50",
            "TSMMC_JUMP_TEMP": "60",
            "ANA_RESIDUE_PAIRS": ["1 3", "3 1"],
        },
    )
    assert len(_warnings(out)) == 4
    lookup = parser.keyword_lookup
    assert (
        lookup["PARALLEL_THREADS"],
        lookup["SLITHER_SUBSTEPS"],
        lookup["TSMMC_JUMP_TEMP"],
    ) == (3, 50, 60.0)
    assert lookup["ANA_RESIDUE_PAIRS"] == [[1, 3], [1, 3]]


@pytest.mark.parametrize("overrides", NOT_IGNORED)
def test_no_warning_when_the_keyword_is_used_or_deliberately_off(
    tmp_path, monkeypatch, capsys, overrides
):
    _, out = _parse(tmp_path, monkeypatch, capsys, overrides)
    assert _warnings(out) == []


@pytest.mark.parametrize(
    "overrides,keyword",
    [
        ({"QUENCH_FREQ": "-1"}, "QUENCH_FREQ"),
        ({"QUENCH_FREQ": "0"}, "QUENCH_FREQ"),
        ({"TSMMC_JUMP_TEMP": "-5"}, "TSMMC_JUMP_TEMP"),
        ({"TSMMC_FIXED_OFFSET": "-3"}, "TSMMC_FIXED_OFFSET"),
    ],
)
def test_a_keyword_that_will_be_ignored_is_still_validated(
    tmp_path, monkeypatch, capsys, overrides, keyword
):
    with pytest.raises(KeyFileException, match=keyword):
        _parse(tmp_path, monkeypatch, capsys, overrides)


def test_the_written_keyfile_reparses_identically_and_without_new_warnings(
    tmp_path, monkeypatch, capsys
):
    """keyfile_used.kf spells out every keyword; none of them may now be 'ignored'."""
    first, out = _parse(tmp_path, monkeypatch, capsys, {"ANALYSIS_FREQ": "5"})
    assert _warnings(out) == []
    first.write_keyfile("written.kf")
    second = keyfile_parser.KeyFileParser("written.kf")
    assert _warnings(capsys.readouterr().out) == []
    skip = ("__KEYFILE",)
    assert {k: v for k, v in second.keyword_lookup.items() if k not in skip} == {
        k: v for k, v in first.keyword_lookup.items() if k not in skip
    }


# ---------------------------------------------------------------------------
# K1-10: messages and announcements
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "overrides,names",
    [
        ({"DIMENSIONS": "6 12 12"}, ["DIMENSIONS"]),
        ({"DIMENSIONS": "6 12"}, ["DIMENSIONS"]),
        (
            {"RESIZED_EQUILIBRATION": "8 8", "EQUILIBRATION": "5"},
            ["RESIZED_EQUILIBRATION", "DIMENSIONS"],
        ),
        (
            {"RESIZED_EQUILIBRATION": "8 8 20", "EQUILIBRATION": "5"},
            ["RESIZED_EQUILIBRATION", "DIMENSIONS"],
        ),
        (
            {
                "RESIZED_EQUILIBRATION": "8 8 8",
                "EQUILIBRATION_OFFSET": "1 1",
                "EQUILIBRATION": "5",
            },
            ["RESIZED_EQUILIBRATION", "EQUILIBRATION_OFFSET"],
        ),
        (
            {
                "RESIZED_EQUILIBRATION": "8 8 8",
                "EQUILIBRATION_OFFSET": "1 1 5",
                "EQUILIBRATION": "5",
            },
            ["RESIZED_EQUILIBRATION", "EQUILIBRATION_OFFSET", "DIMENSIONS"],
        ),
        ({"ANA_RESIDUE_PAIRS": ["0 5"]}, ["ANA_RESIDUE_PAIRS"]),
        ({"ANA_RESIDUE_PAIRS": ["-1 3"]}, ["ANA_RESIDUE_PAIRS"]),
    ],
)
def test_refusals_name_the_keyword(tmp_path, monkeypatch, capsys, overrides, names):
    with pytest.raises(KeyFileException) as error:
        _parse(tmp_path, monkeypatch, capsys, overrides)
    for name in names:
        assert name in str(error.value)


def test_default_announcements(tmp_path, monkeypatch, capsys):
    _, out = _parse(tmp_path, monkeypatch, capsys)
    # the obsolete keyword is not announced, and RESTART_FREQ is announced as the
    # number the run will use: N_STEPS / 10 = 40 / 10
    assert "CRANKSHAFT_MODE" not in out
    assert "No RESTART_FREQ set - using default [4]" in out
    assert "10th-percentile" not in out
    parser, out = _parse(
        tmp_path, monkeypatch, capsys, {"N_STEPS": "7", "XTC_FREQ": "1", "EN_FREQ": "1"}
    )
    assert "No RESTART_FREQ set - using default [1]" in out
    assert parser.keyword_lookup["RESTART_FREQ"] == 1


# ---------------------------------------------------------------------------
# M3-3 / X1-9: PARALLEL_THREADS
# ---------------------------------------------------------------------------


def test_parallel_threads_has_a_ceiling(tmp_path, monkeypatch, capsys):
    for value in ("3000000000", "1025"):
        with pytest.raises(KeyFileException, match="PARALLEL_THREADS") as error:
            _parse(
                tmp_path,
                monkeypatch,
                capsys,
                {"PARALLELIZE": "True", "PARALLEL_THREADS": value},
            )
        assert "1024" in str(error.value)
    parser, _ = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {"PARALLELIZE": "True", "PARALLEL_THREADS": "1024"},
    )
    assert parser.keyword_lookup["PARALLEL_THREADS"] == 1024
    assert CONFIG.MAX_PARALLEL_THREADS == 1024


def test_zero_threads_in_a_simulation_follows_omp_num_threads(
    tmp_path, monkeypatch, capsys
):
    """
    The behaviour itself: a Simulation built with PARALLEL_THREADS : 0 uses
    OMP_NUM_THREADS, where it used to use the machine's core count. The value is
    chosen so that it cannot be the core count by accident.
    """
    from pimms.simulation import Simulation

    wanted = (os.cpu_count() or 1) + 3
    monkeypatch.setenv("OMP_NUM_THREADS", str(wanted))
    parser, _ = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {"PARALLELIZE": "True", "PARALLEL_THREADS": "0"},
    )
    simulation = Simulation(parser.keyword_lookup)
    assert simulation.parallel_threads == wanted


def test_zero_threads_resolves_from_omp_num_threads_then_the_available_cpus(
    monkeypatch,
):
    """Unit test of the new helper (the behaviour is pinned by the test above)."""
    monkeypatch.setattr(keyfile_parser, "available_cpu_count", lambda: 5, raising=False)
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    assert keyfile_parser.resolve_parallel_threads(0) == 5
    assert keyfile_parser.resolve_parallel_threads(None) == 5
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    assert keyfile_parser.resolve_parallel_threads(0) == 3
    monkeypatch.setenv("OMP_NUM_THREADS", "2,4")  # nested-region list: the first level
    assert keyfile_parser.resolve_parallel_threads(0) == 2
    for junk in ("", "0", "-2", "many", "1_0"):
        monkeypatch.setenv("OMP_NUM_THREADS", junk)
        assert keyfile_parser.resolve_parallel_threads(0) == 5
    # an explicit request is used as given, whatever the environment says
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    assert keyfile_parser.resolve_parallel_threads(7) == 7
    # and nothing ever exceeds the ceiling the parser enforces
    monkeypatch.setenv("OMP_NUM_THREADS", "99999")
    assert keyfile_parser.resolve_parallel_threads(0) == 1024


def test_the_available_cpu_count_is_the_process_allowance_not_the_machine(monkeypatch):
    """A 4-CPU allocation on a 128-core node must give 4."""
    monkeypatch.setattr(os, "cpu_count", lambda: 128)
    monkeypatch.setattr(os, "process_cpu_count", lambda: 4, raising=False)
    assert keyfile_parser.available_cpu_count() == 4
    monkeypatch.delattr(os, "process_cpu_count")
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1, 2}, raising=False)
    assert keyfile_parser.available_cpu_count() == 3
    monkeypatch.delattr(os, "sched_getaffinity")
    assert keyfile_parser.available_cpu_count() == 128


# ---------------------------------------------------------------------------
# K1-9: command line
# ---------------------------------------------------------------------------


def _run_cli(args: Iterable[str], cwd: pathlib.Path) -> subprocess.CompletedProcess:
    """
    Run ``scripts/PIMMS`` with the test interpreter.

    Parameters
    ----------
    args : iterable of str
        Command-line arguments.
    cwd : pathlib.Path
        Working directory for the run.

    Returns
    -------
    subprocess.CompletedProcess
        With ``returncode``, ``stdout`` and ``stderr``.
    """
    env = {
        k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONSAFEPATH")
    }
    env["PYTHONPATH"] = str(REPO_ROOT)
    env["OMP_NUM_THREADS"] = "1"
    return subprocess.run(
        [sys.executable, str(PIMMS_SCRIPT), *args],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )


def _write_cli_inputs(directory: pathlib.Path) -> pathlib.Path:
    """
    Write a runnable keyfile and its parameter file into ``directory``.

    Parameters
    ----------
    directory : pathlib.Path
        Where to write.

    Returns
    -------
    pathlib.Path
        The keyfile.
    """
    (directory / "params.prm").write_text(PARAMETER_TEXT)
    keyfile = directory / "KEYFILE.kf"
    keyfile.write_text(
        _keyfile_text(
            {
                "N_STEPS": "2",
                "XTC_FREQ": "1",
                "EN_FREQ": "1",
                "PARAMETER_FILE": str(directory / "params.prm"),
            }
        )
    )
    return keyfile


@pytest.mark.parametrize(
    "arguments",
    [
        ["-k", "KEYFILE.kf", "--info"],
        ["--info", "-k", "KEYFILE.kf"],
        ["-k", "KEYFILE.kf", "--info", "SEED"],
        ["-v", "-k", "KEYFILE.kf"],
    ],
)
def test_a_keyfile_with_info_or_version_is_a_usage_error(tmp_path, arguments):
    """These exited 0 without running the keyfile."""
    _write_cli_inputs(tmp_path)
    result = _run_cli(arguments, tmp_path)
    assert result.returncode == 2
    assert "cannot be combined" in result.stderr
    assert not (tmp_path / "log.txt").exists()


def test_a_repeated_keyfile_option_is_a_usage_error(tmp_path):
    """The last one used to win without a word."""
    _write_cli_inputs(tmp_path)
    result = _run_cli(["-k", "missing.kf", "-k", "KEYFILE.kf"], tmp_path)
    assert result.returncode == 2
    assert "more than once" in result.stderr
    assert not (tmp_path / "log.txt").exists()


@pytest.mark.parametrize("option", ["--keyfile", "-keyfile", "-k"])
def test_every_spelling_of_the_keyfile_option_reaches_the_keyfile(tmp_path, option):
    """--keyfile was refused as an unknown argument (exit 2); a missing file is exit 1."""
    result = _run_cli([option, "does_not_exist.kf"], tmp_path)
    assert result.returncode == 1
    assert 'Could not open file "does_not_exist.kf"' in result.stderr


def test_version_with_info_is_a_usage_error(tmp_path):
    """--version used to win and --info was dropped without a word."""
    for arguments in (["-v", "--info"], ["--info", "SEED", "--version"]):
        result = _run_cli(arguments, tmp_path)
        assert result.returncode == 2
        assert "cannot be combined" in result.stderr
        assert "version " not in result.stdout


def test_a_working_directory_that_cannot_take_output_is_refused_before_the_parse(
    tmp_path,
):
    """
    The wiring of the working-directory check, with no dependence on chmod.

    The child process removes its own (empty) working directory and then runs
    ``PIMMS -k``. Output cannot be written there, so the run must stop with one
    line on stderr and exit status 1 before the keyfile is parsed (the parse ends
    by starting log.txt, which used to be where this failed, with a traceback).
    What this does not cover is the other branch of the check, a directory that
    exists but is not writable: its logic is unit-tested below, and the CLI is
    only run against one when the filesystem honours chmod (next test but one).
    """
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    keyfile = _write_cli_inputs(inputs)
    doomed = tmp_path / "doomed"
    doomed.mkdir()
    wrapper = (
        "import os, runpy, sys; os.rmdir(os.getcwd()); "
        "sys.argv = ['PIMMS', '-k', %r]; runpy.run_path(%r, run_name='__main__')"
        % (str(keyfile), str(PIMMS_SCRIPT))
    )
    env = {
        k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONSAFEPATH")
    }
    env["PYTHONPATH"] = str(REPO_ROOT)
    result = subprocess.run(
        [sys.executable, "-c", wrapper],
        cwd=str(doomed),
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    assert not doomed.exists()
    assert result.returncode == 1
    assert "no longer exists" in result.stderr
    assert "Traceback" not in result.stderr
    # the keyfile was never parsed: the parser's first announcement is absent
    assert "Parsing keyfile" not in result.stdout
    assert sorted(f.name for f in inputs.iterdir()) == ["KEYFILE.kf", "params.prm"]


def _load_cli_module() -> types.ModuleType:
    """
    Import ``scripts/PIMMS`` as a module, without running it.

    The script has no ``.py`` extension, so it is loaded by path; everything it
    does on being run sits under its ``__main__`` guard.

    Returns
    -------
    types.ModuleType
        The loaded script.
    """
    loader = importlib.machinery.SourceFileLoader(
        "pimms_command_line", str(PIMMS_SCRIPT)
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_the_working_directory_check_names_the_directory(tmp_path, monkeypatch):
    """Unit test of the check itself (a new function; the CLI test above pins the wiring)."""
    cli = _load_cli_module()
    monkeypatch.chdir(tmp_path)
    assert cli.working_directory_problem() is None

    # a directory that cannot be written to, whatever this filesystem makes of chmod
    monkeypatch.setattr(os, "access", lambda path, mode: not (mode & os.W_OK))
    problem = cli.working_directory_problem()
    assert "not writable" in problem
    assert str(pathlib.Path.cwd()) in problem

    def gone() -> str:
        raise FileNotFoundError("the working directory was deleted")

    monkeypatch.setattr(os, "getcwd", gone)
    assert "no longer exists" in cli.working_directory_problem()


def test_a_read_only_working_directory_where_the_filesystem_honours_chmod(tmp_path):
    """
    The CLI against a directory that exists but is not writable.

    This one is opportunistic, and skips (visibly) in two situations it cannot
    do anything about: a user who can write regardless of the mode bits (root),
    and a filesystem that does not keep the directory read-only for the length
    of the run (a Dropbox folder restores the write bit within about a second).
    The check is covered without it by the two tests above.
    """
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    keyfile = _write_cli_inputs(inputs)
    read_only = tmp_path / "read_only"
    read_only.mkdir()
    read_only.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        if os.access(str(read_only), os.W_OK):
            pytest.skip(
                "this user or filesystem can write to a directory without write permission"
            )
        result = _run_cli(["-k", str(keyfile)], read_only)
        # a sync client (Dropbox, for one) puts the write bit back within a second
        # of the chmod; if that happened under the run there was nothing to refuse
        if os.access(str(read_only), os.W_OK):
            pytest.skip(
                "the directory did not stay read-only for the length of the run"
            )
    finally:
        read_only.chmod(stat.S_IRWXU)
    assert result.returncode == 1
    assert "not writable" in result.stderr
    assert str(read_only) in result.stderr
    assert "Traceback" not in result.stderr
    assert list(read_only.iterdir()) == []


# ---------------------------------------------------------------------------
# RF4-1: what is read as UTF-8 is written as UTF-8, whatever the locale
# ---------------------------------------------------------------------------


def test_echo_files_are_utf8_and_reusable_under_a_non_utf8_locale(tmp_path):
    """
    Parse and echo in a child process whose locale encoding is not UTF-8.

    The inputs are read as UTF-8 whatever the locale, so the files that echo
    them (``parameters_used.prm``, a written keyfile, the angle summary,
    ``log.txt``) have to be written as UTF-8 too: in the locale encoding the
    write failed on a character the locale lacks, or produced a copy that was
    then refused as "not UTF-8" when used as an input.
    """
    (tmp_path / "params.prm").write_text(NON_ASCII_PARAMETERS, encoding="utf-8")
    (tmp_path / "KEYFILE.kf").write_text(
        _keyfile_text({"CHAIN": ["2 A\u00c5\u03b5A  # \u00c5 per site"]}),
        encoding="utf-8",
    )
    (tmp_path / "child.py").write_text(LOCALE_CHILD, encoding="ascii")
    env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith(("LC_", "LANG", "PYTHON"))
    }
    env.update(
        LC_ALL="en_US.ISO8859-1",
        PYTHONUTF8="0",
        PYTHONCOERCECLOCALE="0",
        PYTHONPATH=str(REPO_ROOT),
    )
    result = subprocess.run(
        [sys.executable, "child.py"],
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
        encoding="ascii",
        timeout=300,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-2000:]
    report = json.loads(result.stdout.strip().splitlines()[-1])
    # the point of the test: the child really did run in a non-UTF-8 locale
    assert report["encoding"].lower().replace("-", "") not in ("utf8", "utf8sig")

    # every echo decodes as UTF-8 and holds the characters that went in
    echo = (tmp_path / "parameters_used.prm").read_bytes().decode("utf-8")
    assert "# \u00c5 and \u03b5\n" in echo
    written = (tmp_path / "written.kf").read_bytes().decode("utf-8")
    # chain sequences are upper-cased on read: epsilon becomes capital epsilon
    assert "A\u00c5\u0395A" in written
    assert "\u03b5.txt" in (tmp_path / "log.txt").read_bytes().decode("utf-8")
    summary = (
        (tmp_path / "absolute_energies_of_angles.txt").read_bytes().decode("utf-8")
    )
    assert "\u00c5 -> " in summary

    # and each echo re-parses to what the original did
    assert report["chains"][0] == report["chains"][1] == [[2, "A\u00c5\u0395A"]]
    assert report["tables_equal"] is True
    assert report["residues"] == ["0", "A", "\u00c5"]


def test_every_text_file_that_echoes_input_is_opened_as_utf8(
    tmp_path, monkeypatch, capsys
):
    """
    The same property without a subprocess or a locale: each echo file must be
    opened for writing with an explicit UTF-8 encoding.
    """
    from pimms import pimmslogger

    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(NON_ASCII_PARAMETERS, encoding="utf-8")
    (tmp_path / "KEYFILE.kf").write_text(_keyfile_text(), encoding="utf-8")
    real_open = builtins.open
    opened = {}

    def spy(file, mode="r", *args, **kwargs):
        if isinstance(file, str) and "b" not in mode and mode[0] in "wa":
            encoding = kwargs.get("encoding", args[1] if len(args) > 1 else None)
            opened.setdefault(os.path.basename(file), set()).add(encoding)
        return real_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", spy)
    parser = keyfile_parser.KeyFileParser("KEYFILE.kf")
    parser.write_keyfile("written.kf")
    pimmslogger.log_status("a status line")
    pimmslogger.log_warning("a warning line")
    pimmslogger.log_error("an error line")
    parameterfile_parser.parse_energy("params.prm")
    parameterfile_parser.write_angle_parameter_summary(
        parameterfile_parser.parse_angles("params.prm"), "params.prm"
    )
    monkeypatch.undo()
    capsys.readouterr()
    for name in (
        "log.txt",
        "written.kf",
        "parameters_used.prm",
        "absolute_energies_of_angles.txt",
    ):
        assert opened.get(name) == {"utf-8"}, (name, opened.get(name))


# ---------------------------------------------------------------------------
# RF4-2 / RF4-3 / RF4-5: non-ASCII spaces, invisible characters, ASCII whitespace
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "keyword,value",
    [
        ("CHAIN", ["4" + NBSP + "AABBA"]),  # used to act as the column separator
        ("ANA_RESIDUE_PAIRS", ["1" + NBSP + "3"]),
        ("HARDWALL", "True" + NBSP),
        ("N_STEPS", "40" + NBSP),
        ("DIMENSIONS", "12 12" + NBSP + "12"),
        (
            "N_STEPS",
            "\u200b40",
        ),  # zero-width space: not whitespace to Python, invisible all the same
        ("TEMPERATURE", "10\u3000"),  # ideographic space
    ],
)
def test_a_non_ascii_space_in_a_keyfile_is_refused_by_name(
    tmp_path, monkeypatch, capsys, keyword, value
):
    with pytest.raises(KeyFileException) as error:
        _parse(tmp_path, monkeypatch, capsys, {keyword: value})
    culprit = [
        ch for ch in (value[0] if isinstance(value, list) else value) if ord(ch) > 127
    ][0]
    assert keyword in str(error.value)
    assert "U+%04X" % ord(culprit) in str(error.value)


def test_a_non_ascii_space_before_a_keyword_is_refused_by_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    (tmp_path / "k.kf").write_text(NBSP + _keyfile_text(), encoding="utf-8")
    with pytest.raises(KeyFileException, match="U\\+00A0 NO-BREAK SPACE") as error:
        keyfile_parser.KeyFileParser("k.kf")
    assert "DIMENSIONS" in str(error.value)


def test_a_file_name_may_hold_a_non_ascii_space(tmp_path, monkeypatch, capsys):
    """The one place a no-break space is data: inside the value of a path keyword."""
    name = "my" + NBSP + "params.prm"
    (tmp_path / name).write_text(PARAMETER_TEXT)
    parser, _ = _parse(tmp_path, monkeypatch, capsys, {"PARAMETER_FILE": name})
    assert parser.keyword_lookup["PARAMETER_FILE"] == name


@pytest.mark.parametrize(
    "line",
    ["A A -3" + NBSP, "A" + NBSP + "A" + NBSP + "-3", NBSP + "A A -3", "A A\u2009-3"],
)
def test_a_non_ascii_space_in_the_parameter_file_is_refused_by_name(
    tmp_path, monkeypatch, line
):
    """`A<NBSP>A<NBSP>-3` used to be read as three columns."""
    monkeypatch.chdir(tmp_path)
    culprit = [ch for ch in line if ord(ch) > 127][0]
    (tmp_path / "p.prm").write_text(
        "A 0 0  # a no-break space in a comment is fine:" + NBSP + "\n" + line + "\n"
        "ANGLE_PENALTY A 0 0 0\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterFileException) as error:
        parameterfile_parser.parse_energy("p.prm")
    assert "U+%04X" % ord(culprit) in str(error.value)
    assert "Line 1 " in str(error.value)

    (tmp_path / "a.prm").write_text(
        "A 0 0\nA A -3\nANGLE_PENALTY A 0" + culprit + "0 0\n", encoding="utf-8"
    )
    with pytest.raises(ParameterFileException) as error:
        parameterfile_parser.parse_angles("a.prm")
    assert "U+%04X" % ord(culprit) in str(error.value)


@pytest.mark.parametrize(
    "line", ["C 1" + NBSP + "3", "C" + NBSP + "1 3", "C 1 3" + NBSP]
)
def test_a_non_ascii_space_in_the_freeze_file_is_refused_by_name(tmp_path, line):
    """`C 1<NBSP>3` used to freeze chains 1 and 3."""
    path = tmp_path / "f.frz"
    path.write_text(
        "# a comment may hold one:" + NBSP + "\n" + line + "\n", encoding="utf-8"
    )
    with pytest.raises(KeyFileException) as error:
        data_structures.FreezeFile(str(path))
    assert "U+00A0 NO-BREAK SPACE" in str(error.value)
    assert "Line 2 " in str(error.value) and "f.frz" in str(error.value)


def test_columns_still_split_on_every_ascii_whitespace_character(
    tmp_path, monkeypatch, capsys
):
    """Vertical tab and form feed separated values before, and still do."""
    parser, _ = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {"DIMENSIONS": "12\x0b13\x0c14", "CHAIN": ["4\tAABBA"]},
    )
    assert parser.keyword_lookup["DIMENSIONS"] == [12, 13, 14]
    assert parser.keyword_lookup["CHAIN"] == [[4, "AABBA"]]
    assert data_structures.split_ascii_whitespace(" a\tb\x0bc\x0cd \r\n") == [
        "a",
        "b",
        "c",
        "d",
    ]
    assert data_structures.split_ascii_whitespace("a" + NBSP + "b") == [
        "a" + NBSP + "b"
    ]
    path = tmp_path / "f.frz"
    path.write_text("C 1\t3\x0b5\n")
    assert data_structures.FreezeFile(str(path)).chains == [1, 3, 5]


def test_a_second_byte_order_mark_is_named(tmp_path, monkeypatch):
    """A mark that is not at byte 0 used to give `unsupported keyword - [DIMENSIONS]`."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    (tmp_path / "k.kf").write_bytes(
        b"\xef\xbb\xbf\xef\xbb\xbf" + _keyfile_text().encode("utf-8")
    )
    with pytest.raises(KeyFileException) as error:
        keyfile_parser.KeyFileParser("k.kf")
    assert "U+FEFF" in str(error.value) and "byte-order mark" in str(error.value)
    assert "\ufeff" not in str(error.value)


def test_utf16_without_a_byte_order_mark_is_explained(tmp_path, monkeypatch):
    """It decodes as UTF-8 full of NULs; the keyword was shown as `[D\\x00I\\x00...]`."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    (tmp_path / "k.kf").write_bytes(_keyfile_text().encode("utf-16-le"))
    with pytest.raises(KeyFileException) as error:
        keyfile_parser.KeyFileParser("k.kf")
    assert "NUL" in str(error.value) and "UTF-16" in str(error.value)
    assert "\x00" not in str(error.value)


def test_an_unsupported_keyword_with_non_ascii_letters_says_so(
    tmp_path, monkeypatch, capsys
):
    with pytest.raises(KeyFileException) as error:
        _parse(tmp_path, monkeypatch, capsys, {"N_ST\uff25PS": "40", "N_STEPS": None})
    assert "unsupported keyword" in str(error.value)
    assert "U+FF25" in str(error.value)


def test_a_utf32_file_is_not_called_utf16(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "params.prm").write_text(PARAMETER_TEXT)
    (tmp_path / "k.kf").write_bytes(_keyfile_text().encode("utf-32"))
    with pytest.raises(KeyFileException) as error:
        keyfile_parser.KeyFileParser("k.kf")
    assert "UTF-32" in str(error.value) and "UTF-16" not in str(error.value)


# ---------------------------------------------------------------------------
# the final box and the final chain list, with a real restart file
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def restart_file(tmp_path_factory) -> pathlib.Path:
    """
    A real ``restart.pimms``: 4 x AABBA in a periodic 12 x 12 x 12 box, written
    by a two-step run of the base keyfile.

    Returns
    -------
    pathlib.Path
        The restart file.
    """
    from pimms.simulation import Simulation

    directory = tmp_path_factory.mktemp("restart_source")
    previous = os.getcwd()
    os.chdir(directory)
    try:
        (directory / "params.prm").write_text(PARAMETER_TEXT)
        (directory / "KEYFILE.kf").write_text(
            _keyfile_text({"N_STEPS": "2", "XTC_FREQ": "1", "EN_FREQ": "1"})
        )
        with contextlib.redirect_stdout(io.StringIO()):
            parser = keyfile_parser.KeyFileParser("KEYFILE.kf")
            Simulation(parser.keyword_lookup).run_simulation()
    finally:
        os.chdir(previous)
    path = directory / "restart.pimms"
    assert path.is_file()
    return path


def test_the_box_checks_use_the_restart_box_under_override_dimensions(
    tmp_path, monkeypatch, capsys, restart_file
):
    """RESTART_OVERRIDE_DIMENSIONS replaces the keyfile box with the restart file's 12^3."""
    restart = {
        "RESTART_FILE": str(restart_file),
        "CHAIN": None,
        "RESTART_OVERRIDE_DIMENSIONS": "True",
    }
    # the keyfile box is fine at this spacing (8 x 900 = 7200) and the restart box
    # is not (12 x 900 = 10800): the box the run will use is the one that counts
    with pytest.raises(KeyFileException) as error:
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {**restart, "DIMENSIONS": "8 8 8", "LATTICE_TO_ANGSTROMS": "900"},
        )
    assert "LATTICE_TO_ANGSTROMS" in str(error.value) and "10800" in str(error.value)

    # and the converse: a keyfile box that could never be allocated is not the
    # box of the run, so it is not what the memory check looks at
    _plenty_of_memory(monkeypatch, 16 * 10**9)
    parser, out = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {**restart, "DIMENSIONS": "3000 3000 3000", "LATTICE_TO_ANGSTROMS": "1"},
    )
    assert list(parser.keyword_lookup["DIMENSIONS"]) == [12, 12, 12]
    assert _warnings(out) == []


def test_the_chain_checks_use_the_restart_chains_and_the_extra_chains(
    tmp_path, monkeypatch, capsys, restart_file
):
    """The keyfile's own CHAIN lines are ignored with a restart file: it holds 4 chains."""
    _plenty_of_memory(monkeypatch)
    restart = {
        "RESTART_FILE": str(restart_file),
        "MOVE_SLITHER": "0.5",
        "MOVE_CHAIN_TRANSLATE": "0",
    }
    # 10^9 x the ONE chain of the keyfile would fit a C int; x the 4 restart chains it does not
    with pytest.raises(KeyFileException) as error:
        _parse(
            tmp_path,
            monkeypatch,
            capsys,
            {**restart, "CHAIN": ["1 AABBA"], "SLITHER_SUBSTEPS": str(10**9)},
        )
    assert "SLITHER_SUBSTEPS" in str(error.value) and str(4 * 10**9) in str(error.value)
    # 5 x 10^8 x the NINE chains of the keyfile would not fit; x 4 it does
    parser, _ = _parse(
        tmp_path,
        monkeypatch,
        capsys,
        {**restart, "CHAIN": ["9 AABBA"], "SLITHER_SUBSTEPS": str(5 * 10**8)},
    )
    assert parser.keyword_lookup["CHAIN"] == [[4, "AABBA"]]

    # a threshold of 4 can never be exceeded by 4 chains, and can by 4 + 3 extra ones
    threshold = {
        "RESTART_FILE": str(restart_file),
        "CHAIN": None,
        "ANA_CLUSTER_THRESHOLD": "4",
    }
    _, out = _parse(tmp_path, monkeypatch, capsys, threshold)
    assert len(_warnings(out)) == 1 and "ANA_CLUSTER_THRESHOLD" in _warnings(out)[0]
    parser, out = _parse(
        tmp_path, monkeypatch, capsys, {**threshold, "EXTRA_CHAIN": ["3 AB"]}
    )
    assert sum(count for count, _ in parser.keyword_lookup["CHAIN"]) == 7
    assert _warnings(out) == []
