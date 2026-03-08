import types

import pytest

from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import KeyFileException


def _write_keyfile(tmp_path, text):
    keyfile = tmp_path / "KEYFILE.kf"
    keyfile.write_text(text)
    return keyfile


def test_parse_chain_success_parse_only(tmp_path):
    keyfile = _write_keyfile(
        tmp_path,
        """
CHAIN : 2 ab
N_STEPS : 10
""".strip(),
    )

    parser = KeyFileParser(str(keyfile), parse_only=True)
    assert parser.keyword_lookup["CHAIN"] == [[2, "ab"]]
    assert parser.keyword_lookup["N_STEPS"] == 10


def test_parse_chain_malformed_raises_keyfileexception(tmp_path):
    keyfile = _write_keyfile(tmp_path, "CHAIN : 2\n")

    with pytest.raises(KeyFileException, match="Invalid CHAIN keyword format"):
        KeyFileParser(str(keyfile), parse_only=True)


def test_parse_extra_chain_malformed_raises_keyfileexception(tmp_path):
    keyfile = _write_keyfile(tmp_path, "EXTRA_CHAIN : not_an_int AB\n")

    with pytest.raises(KeyFileException, match="Invalid EXTRA_CHAIN keyword format"):
        KeyFileParser(str(keyfile), parse_only=True)


def test_parse_ana_residue_pairs_malformed_raises_keyfileexception(tmp_path):
    keyfile = _write_keyfile(tmp_path, "ANA_RESIDUE_PAIRS : 2\n")

    with pytest.raises(KeyFileException, match="Invalid ANA_RESIDUE_PAIRS format"):
        KeyFileParser(str(keyfile), parse_only=True)


def test_parse_parameter_file_with_colon_in_value(tmp_path):
    keyfile = _write_keyfile(tmp_path, "PARAMETER_FILE : C:/tmp/params.prm\n")

    parser = KeyFileParser(str(keyfile), parse_only=True)
    assert parser.keyword_lookup["PARAMETER_FILE"] == "C:/tmp/params.prm"


def test_parse_duplicate_non_repeatable_keyword_raises(tmp_path):
    keyfile = _write_keyfile(
        tmp_path,
        """
N_STEPS : 10
N_STEPS : 20
""".strip(),
    )

    with pytest.raises(KeyFileException, match="second occurence"):
        KeyFileParser(str(keyfile), parse_only=True)


def test_set_dynamic_defaults_disables_negative_frequency_and_sets_restart_freq():
    parser = KeyFileParser.__new__(KeyFileParser)
    parser.keyword_lookup = {
        "ANALYSIS_MODULE": False,
        "ANA_CUSTOM": 0,
        "N_STEPS": 100,
        "ANALYSIS_FREQ": 50,
        "ANA_POL": 0,
        "ANA_INTSCAL": -1,
        "ANA_DISTMAP": 0,
        "ANA_ACCEPTANCE": 0,
        "ANA_INTER_RESIDUE": 0,
        "ANA_CLUSTER": 0,
        "ENERGY_CHECK": 0,
        "RESTART_FREQ": "Every 10th-percentile",
    }

    parser.set_dynamic_defaults()

    assert parser.keyword_lookup["ANA_CUSTOM"] == 110
    assert parser.keyword_lookup["ANA_POL"] == 110
    assert parser.keyword_lookup["ANA_INTSCAL"] == 110
    assert parser.keyword_lookup["ENERGY_CHECK"] == 110
    assert parser.keyword_lookup["RESTART_FREQ"] == 10


def test_restart_file_sanity_updates_chain_from_restart_data():
    parser = KeyFileParser.__new__(KeyFileParser)

    class DummyRestart:
        def __init__(self):
            self.chains = {
                1: [None, "AB", 2],
                2: [None, "AB", 2],
                3: [None, "CD", 5],
            }
            self.dimensions = [10, 10, 10]
            self.hardwall = True
            self.extra_chains = []

        def update_lattice_dimensions(self, dims):
            self.dimensions = list(dims)

    restart_obj = DummyRestart()

    parser.keyword_lookup = {
        "RESTART_FILE": restart_obj,
        "RESTART_OVERRIDE_HARDWALL": True,
        "HARDWALL": False,
        "RESIZED_EQUILIBRATION": False,
        "RESTART_OVERRIDE_DIMENSIONS": True,
        "DIMENSIONS": [8, 8, 8],
        "EXTRA_CHAIN": [],
    }

    parser.sanity_check_and_update_with_restart_file()

    assert parser.keyword_lookup["CHAIN"] == [[2, "AB"], [1, "CD"]]
    assert parser.keyword_lookup["DIMENSIONS"] == [10, 10, 10]
