import random

import numpy as np
import pytest

from pimms import CONFIG
from pimms.acceptance import AcceptanceCalculator
from pimms.latticeExceptions import AcceptanceException


_MOVE_KEYS = [
    "MOVE_CRANKSHAFT",
    "MOVE_CHAIN_TRANSLATE",
    "MOVE_CHAIN_ROTATE",
    "MOVE_CHAIN_PIVOT",
    "MOVE_HEAD_PIVOT",
    "MOVE_SLITHER",
    "MOVE_CLUSTER_TRANSLATE",
    "MOVE_CLUSTER_ROTATE",
    "MOVE_CTSMMC",
    "MOVE_MULTICHAIN_TSMMC",
    "MOVE_RATCHET_PIVOT",
    "MOVE_SYSTEM_TSMMC",
    "MOVE_JUMP_AND_RELAX",
]


def _moveset_with_single_active(active_key):
    d = {k: 0.0 for k in _MOVE_KEYS}
    d[active_key] = 1.0
    return d


def _uniform_moveset():
    w = 1.0 / len(_MOVE_KEYS)
    return {k: w for k in _MOVE_KEYS}


def test_init_sets_temperature_and_invtemp():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())
    assert ac.temperature == 300.0
    assert ac.invtemp == pytest.approx(CONFIG.INVTEMP_FACTOR / 300.0)


@pytest.mark.parametrize("key, expected_code", [
    ("MOVE_CRANKSHAFT", 1),
    ("MOVE_CHAIN_TRANSLATE", 2),
    ("MOVE_CHAIN_ROTATE", 3),
    ("MOVE_CHAIN_PIVOT", 4),
    ("MOVE_HEAD_PIVOT", 5),
    ("MOVE_SLITHER", 6),
    ("MOVE_CLUSTER_TRANSLATE", 7),
    ("MOVE_CLUSTER_ROTATE", 8),
    ("MOVE_CTSMMC", 9),
    ("MOVE_MULTICHAIN_TSMMC", 10),
    ("MOVE_RATCHET_PIVOT", 11),
    ("MOVE_SYSTEM_TSMMC", 12),
    ("MOVE_JUMP_AND_RELAX", 13),
])
def test_move_selector_selects_each_move_when_only_move_enabled(monkeypatch, key, expected_code):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_moveset_with_single_active(key))
    monkeypatch.setattr(random, "random", lambda: 0.5)

    assert ac.move_selector(chain_length=10) == expected_code


def test_move_selector_raises_if_probability_mass_is_unassigned(monkeypatch):
    # Sum of move frequencies < 1 leaves dead range and should trigger AcceptanceException.
    lookup = {k: 0.0 for k in _MOVE_KEYS}
    lookup["MOVE_CRANKSHAFT"] = 0.2
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=lookup)

    monkeypatch.setattr(random, "random", lambda: 0.95)
    with pytest.raises(AcceptanceException):
        ac.move_selector(chain_length=10)


@pytest.mark.parametrize("key", [
    "MOVE_CHAIN_ROTATE",
    "MOVE_CHAIN_PIVOT",
    "MOVE_HEAD_PIVOT",
    "MOVE_SLITHER",
    "MOVE_RATCHET_PIVOT",
])
def test_move_selector_remaps_singleton_chain_moves_to_crankshaft(monkeypatch, key):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_moveset_with_single_active(key))
    monkeypatch.setattr(random, "random", lambda: 0.5)

    assert ac.move_selector(chain_length=1) == 1


@pytest.mark.parametrize("key", [
    "MOVE_CTSMMC",
    "MOVE_MULTICHAIN_TSMMC",
    "MOVE_SYSTEM_TSMMC",
])
def test_move_selector_remaps_tsmmc_moves_when_aux_chain(monkeypatch, key):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_moveset_with_single_active(key))
    ac.auxillary_chain = True
    monkeypatch.setattr(random, "random", lambda: 0.5)

    assert ac.move_selector(chain_length=10) == 1


def test_boltzmann_acceptance_accepts_downhill_or_equal():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    assert ac.boltzmann_acceptance(old_energy=10, new_energy=9) is True
    assert ac.boltzmann_acceptance(old_energy=10, new_energy=10) is True


def test_boltzmann_acceptance_uphill_accepts_when_random_below_expterm(monkeypatch):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    monkeypatch.setattr(random, "random", lambda: 0.0)
    assert ac.boltzmann_acceptance(old_energy=0, new_energy=1) is True


def test_boltzmann_acceptance_uphill_rejects_when_random_above_expterm(monkeypatch):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    monkeypatch.setattr(random, "random", lambda: 1.0)
    assert ac.boltzmann_acceptance(old_energy=0, new_energy=1) is False


def test_update_temperature_updates_temperature_and_invtemp():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())
    ac.update_temperature(250.0)

    assert ac.temperature == 250.0
    assert ac.invtemp == pytest.approx(CONFIG.INVTEMP_FACTOR / 250.0)


def test_update_move_logs_updates_main_chain_counts():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    ac.update_move_logs(selection=2, acceptance=False)
    ac.update_move_logs(selection=2, acceptance=True)

    assert ac.move_count[2] == 2
    assert ac.accepted_count[2] == 1


def test_update_move_logs_updates_aux_chain_counts_when_enabled():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())
    ac.auxillary_chain = True

    ac.update_move_logs(selection=2, acceptance=False)
    ac.update_move_logs(selection=2, acceptance=True)

    assert ac.aux_chain_move_count[2] == 2
    assert ac.aux_chain_accepted_count[2] == 1
    assert ac.move_count[2] == 0
    assert ac.accepted_count[2] == 0


def test_megastep_update_move_logs_updates_correct_counters_main_and_aux():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    ac.megastep_update_move_logs(selection=4, number_accepted=3, number_tried=8)
    assert ac.move_count[4] == 8
    assert ac.accepted_count[4] == 3

    ac.auxillary_chain = True
    ac.megastep_update_move_logs(selection=4, number_accepted=2, number_tried=5)
    assert ac.aux_chain_move_count[4] == 5
    assert ac.aux_chain_accepted_count[4] == 2


def test_alt_markov_chain_update_move_logs_tracks_main_and_aux_totals():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    ac.alt_Markov_chain_update_move_logs(7)
    assert ac.alt_Markov_chain_moves == 7
    assert ac.aux_chain_alt_Markov_chain_moves == 0

    ac.auxillary_chain = True
    ac.alt_Markov_chain_update_move_logs(9)
    assert ac.alt_Markov_chain_moves == 7
    assert ac.aux_chain_alt_Markov_chain_moves == 9


def test_get_total_aux_chain_moves_sums_aux_move_and_alt_chain_counts():
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())
    ac.auxillary_chain = True

    ac.update_move_logs(selection=1, acceptance=True)
    ac.megastep_update_move_logs(selection=2, number_accepted=1, number_tried=5)
    ac.alt_Markov_chain_update_move_logs(4)

    assert ac.get_total_aux_chain_moves() == 1 + 5 + 4


@pytest.mark.parametrize("bad_temp", [0.0, -10.0])
def test_init_rejects_non_positive_temperature(bad_temp):
    with pytest.raises(AcceptanceException):
        AcceptanceCalculator(temp=bad_temp, keyword_lookup=_uniform_moveset())


@pytest.mark.parametrize("bad_temp", [0.0, -1.0])
def test_update_temperature_rejects_non_positive_temperature(bad_temp):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    with pytest.raises(AcceptanceException):
        ac.update_temperature(bad_temp)


@pytest.mark.parametrize("bad_selection", [-1, 0, 14, 100])
def test_update_move_logs_rejects_invalid_selection_index(bad_selection):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    with pytest.raises(AcceptanceException):
        ac.update_move_logs(selection=bad_selection, acceptance=True)


@pytest.mark.parametrize("bad_selection", [-1, 0, 14, 100])
def test_megastep_update_move_logs_rejects_invalid_selection_index(bad_selection):
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())

    with pytest.raises(AcceptanceException):
        ac.megastep_update_move_logs(selection=bad_selection, number_accepted=1, number_tried=2)
