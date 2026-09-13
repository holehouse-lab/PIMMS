"""Focused regressions from the fresh 1.0.8 Monte Carlo core re-review."""

import math
import random

import numpy as np
import pytest

from pimms.acceptance import AcceptanceCalculator
from pimms.latticeExceptions import AcceptanceException
from pimms.moves import (
    MoveObject,
    _parallel_can_move_all_chains,
    _vmmc_cutoff_cdf,
    _vmmc_link_probability,
    _vmmc_offset_shell,
)
from pimms.tests import kernel_test_utils as U


_MOVE_KEYS = (
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
    "MOVE_PULL",
    "MOVE_SYSTEM_TSMMC",
    "MOVE_JUMP_AND_RELAX",
    "MOVE_VMMC",
)


def _only_move(key):
    moves = {name: 0.0 for name in _MOVE_KEYS}
    moves[key] = 1.0
    return moves


def test_mixed_system_singleton_does_not_suppress_global_megamoves(tmp_path, monkeypatch):
    """The outer seed is irrelevant to whole-system SLITHER/PULL eligibility."""
    state = U.build_state(
        tmp_path,
        2,
        "SR",
        False,
        {"MOVE_CRANKSHAFT": 1.0},
        box=[12, 12],
        chains=[(1, "A"), (1, "AAA")],
        seed=17,
        temperature=40,
    )
    chains = list(state.lattice.chains.values())
    # the system holds a monomer; the selector takes no chain argument, so the
    # codes below are what the outer loop would draw whichever chain it picked
    assert sorted(len(chain) for chain in chains) == [1, 3]

    monkeypatch.setattr(random, "random", lambda: 0.5)

    # Per-chain moves are returned as drawn even for the monomer: the move itself
    # rejects them as null moves. (They used to be remapped to a whole-system
    # crankshaft megamove, which contradicted MOVE_CRANKSHAFT : 0.)
    for move_name, move_code in (("MOVE_CHAIN_ROTATE", 3), ("MOVE_CHAIN_PIVOT", 4),
                                 ("MOVE_HEAD_PIVOT", 5)):
        calculator = AcceptanceCalculator(40.0, _only_move(move_name))
        assert calculator.move_selector() == move_code

    # SLITHER and PULL operate on the whole system. Their code must survive the
    # arbitrary selection of a monomer so the eligible polymer is still reached.
    for move_name, move_code in (("MOVE_SLITHER", 6), ("MOVE_PULL", 11)):
        calculator = AcceptanceCalculator(40.0, _only_move(move_name))
        assert calculator.move_selector() == move_code

    # Exercise the megamoves too: slither offers both chains (including its
    # valid monomer-translation path), while pull offers only the eligible 3-mer.
    random.seed(7)
    np.random.seed(7)
    lattice, energy, proposed, _ = state.sim.MOVER.system_slither(
        state.lattice, state.energy, state.acc, state.ham, 1
    )
    assert proposed == 2
    assert energy == state.ham.evaluate_total_energy(lattice)[0]

    lattice, energy, proposed, _ = state.sim.MOVER.system_pull(
        lattice, energy, state.acc, state.ham, 1
    )
    assert proposed == 1
    assert energy == state.ham.evaluate_total_energy(lattice)[0]


def test_singleton_direct_rotation_and_pivots_reject_without_mutation(tmp_path):
    state = U.build_state(
        tmp_path,
        2,
        "SR",
        False,
        {"MOVE_CRANKSHAFT": 1.0},
        box=[9, 9],
        chains=[(1, "A")],
        seed=5,
    )
    monomer = next(iter(state.lattice.chains.values()))
    original_grid = state.lattice.grid.copy()
    mover = MoveObject()

    assert mover.chain_rotate(monomer, state.lattice.grid) == (False, False)
    assert mover.chain_pivot(monomer, state.lattice.grid) == (False, False)
    assert mover.head_pivot(monomer, state.lattice.grid) == (False, False)
    np.testing.assert_array_equal(state.lattice.grid, original_grid)


@pytest.mark.parametrize("beta", [1.0e-12, 0.025, 1.0, 1.0e6])
@pytest.mark.parametrize("delta_energy", [1.0e-9, 0.1, 1.0, 50.0])
def test_vmmc_link_probability_matches_original_formula_in_normal_range(
    beta, delta_energy
):
    expected = 1.0 - math.exp(-beta * delta_energy)
    assert _vmmc_link_probability(beta, delta_energy) == pytest.approx(
        expected, rel=2.0e-8, abs=1.0e-16
    )


@pytest.mark.parametrize(
    "delta_energy, expected",
    [
        (-float("inf"), 0.0),
        (-1.0e308, 0.0),
        (-1000.0, 0.0),
        (0.0, 0.0),
        (1000.0, 1.0),
        (1.0e308, 1.0),
        (float("inf"), 1.0),
    ],
)
def test_vmmc_link_probability_has_finite_extreme_energy_limits(
    delta_energy, expected
):
    # The former exp(-beta*dE)-then-clamp implementation raised OverflowError
    # for the large negative values before it could reach the clamp.
    probability = _vmmc_link_probability(1.0, delta_energy)
    assert math.isfinite(probability)
    assert probability == expected
    assert 0.0 <= probability <= 1.0


def _legacy_vmmc_cutoff(cap, uniform_draw):
    total = sum(1.0 / k for k in range(1, cap + 1))
    running = 0.0
    for k in range(1, cap + 1):
        running += (1.0 / k) / total
        if uniform_draw <= running:
            return k
    return cap


def test_cached_vmmc_cutoff_preserves_draw_boundaries(monkeypatch):
    mover = MoveObject()
    _vmmc_cutoff_cdf.cache_clear()
    cap = 37
    first_boundary = _vmmc_cutoff_cdf(cap)[0]
    draws = (
        0.0,
        np.nextafter(first_boundary, 0.0),
        first_boundary,
        np.nextafter(first_boundary, 1.0),
        0.5,
        np.nextafter(1.0, 0.0),
    )

    for draw in draws:
        monkeypatch.setattr(random, "random", lambda value=draw: value)
        assert mover._vmmc_draw_nc(cap, cap) == _legacy_vmmc_cutoff(cap, draw)

    # The harmonic CDF is built once for a repeated cluster cap instead of on
    # every VMMC attempt.
    assert _vmmc_cutoff_cdf.cache_info().hits >= len(draws)


def test_vmmc_neighbour_offsets_are_cached_with_precomputed_radii():
    _vmmc_offset_shell.cache_clear()
    offsets = _vmmc_offset_shell(3, 3)
    assert offsets is _vmmc_offset_shell(3, 3)
    assert len(offsets) == 7**3 - 1
    assert all(max(abs(value) for value in delta) == radius
               for delta, radius in offsets)
    assert _vmmc_offset_shell.cache_info().hits == 1


def test_parallel_gate_ignores_only_permanently_frozen_long_chains():
    """An immobile 513-mer must not force movable chains onto the serial path."""
    long_length = 513
    idx = np.zeros((long_length + 2, 8), dtype=np.int64)
    idx[:long_length, 4] = 1
    idx[long_length:, 4] = 2
    offsets = np.array([0, long_length], dtype=np.int32)
    lengths = np.array([long_length, 2], dtype=np.int32)
    homo = np.array([0, 1], dtype=np.int32)

    assert not _parallel_can_move_all_chains(
        idx, offsets, lengths, [64, 64, 64], False,
        chain_homo=homo, cap_mode="hetero"
    )
    assert _parallel_can_move_all_chains(
        idx, offsets, lengths, [64, 64, 64], False,
        chain_homo=homo, cap_mode="hetero", frozen_chains=(1,)
    )

    # Freezing the eligible short chain instead does not excuse the over-cap
    # long chain.
    assert not _parallel_can_move_all_chains(
        idx, offsets, lengths, [64, 64, 64], False,
        chain_homo=homo, cap_mode="hetero", frozen_chains=(2,)
    )


@pytest.mark.parametrize("temperature", [float("nan"), float("inf"), -float("inf")])
def test_acceptance_rejects_nonfinite_temperatures(temperature):
    with pytest.raises(AcceptanceException, match="finite and > 0"):
        AcceptanceCalculator(temperature, _only_move("MOVE_CRANKSHAFT"))

    calculator = AcceptanceCalculator(300.0, _only_move("MOVE_CRANKSHAFT"))
    with pytest.raises(AcceptanceException, match="finite and > 0"):
        calculator.update_temperature(temperature)

