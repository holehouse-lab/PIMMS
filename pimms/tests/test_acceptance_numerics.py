"""Definitional tests of the Python Metropolis criterion.

``AcceptanceCalculator.boltzmann_acceptance`` is the acceptance stage for every
chain translate, chain rotate, chain pivot, head pivot and cluster move, and for
the jump of jump-and-relax. For a symmetric proposal it must implement, exactly,

    P(accept | dE > 0) = exp(-dE * invtemp),   invtemp = INVTEMP_FACTOR / T

Anything else samples a Boltzmann distribution at a temperature other than the
one in the keyfile - and because the Cython megamove kernels use the correct
beta, the composite Markov chain then has no single stationary distribution at
all. That is exactly what a units regression, a dropped ``INVTEMP_FACTOR`` or a
mis-typed exponent would produce.

The tests that existed before this file drew ``random.random()`` at exactly 0.0
and 1.0 with dE in {-1, 0, +1}. ``0.0 < expterm`` and ``1.0 < expterm`` are true
and false respectively for EVERY expterm in (0, 1), so the exponent was never
compared against a number: a review injected a factor of 1.2 into it and 274
tests stayed green. What follows compares the criterion against
``math.exp(-dE/T)`` computed independently here, three ways - the exact
acceptance threshold, the acceptance frequency over many draws, and the
unconditional acceptance of downhill moves - and then checks that the frequency
test actually has the power to reject a 20 % error in the exponent.
"""
import math
import random
from typing import Dict, List, Tuple

import pytest

from pimms import CONFIG
from pimms.acceptance import AcceptanceCalculator


_MOVE_KEYS: List[str] = [
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
]

# (temperature, dE) pairs. The sensitivity of the criterion to a wrong exponent
# lives entirely where exp(-dE/T) is far from 1, so the small-T / large-dE
# corner has to be in here: at T = 300, dE = 1 even a 20 % exponent error is
# under 4 sigma at 200000 draws and would be invisible at the sample sizes a
# test can afford.
_CASES: List[Tuple[float, float]] = [
    (10.0, 1.0), (10.0, 2.0), (10.0, 5.0), (10.0, 10.0),
    (40.0, 1.0), (40.0, 5.0), (40.0, 10.0), (40.0, 25.0),
    (300.0, 10.0), (300.0, 50.0), (300.0, 200.0),
]

# draws per (T, dE) row. 50000 puts the binomial sigma at ~0.2 % absolute, which
# separates a 20 % exponent error from the truth by 5-31 sigma on every row.
_N_DRAWS: int = 50000

# how far the measured frequency may sit from the exact value, in binomial sigma.
# Eleven rows at a fixed seed: 4 sigma is ~1 in 15000 per row by chance.
_K_SIGMA: float = 4.0

_SEED: int = 20260910


def _uniform_moveset() -> Dict[str, float]:
    """A keyword_lookup with every move enabled at equal frequency.

    ``AcceptanceCalculator`` reads all fourteen MOVE_* entries to build its move
    selection thresholds; the acceptance criterion itself does not care what they
    are, so an equal split keeps the construction valid and irrelevant.

    Returns
    -------
    dict
        Mapping of each MOVE_* keyword to 1/14.
    """
    w = 1.0 / len(_MOVE_KEYS)
    return {k: w for k in _MOVE_KEYS}


def _row_seed(temperature: float, delta_energy: float) -> int:
    """A distinct RNG seed per (T, dE) row.

    Every row must get its own stream. Sharing one seed makes all the rows read
    the same 50000 uniforms, so whatever the empirical CDF of that particular
    stream does at these thresholds is inherited by every row at once - measured,
    that pushed every z-score to the same sign and the worst row to 2.6 sigma
    where independent streams give 1.7.

    Parameters
    ----------
    temperature : float
        Row temperature.

    delta_energy : float
        Row energy change.

    Returns
    -------
    int
        A deterministic seed for this row.
    """
    return _SEED + 7919 * _CASES.index((temperature, delta_energy))


def _exact_acceptance(temperature: float, delta_energy: float) -> float:
    """The Metropolis acceptance probability, computed from the definition.

    Written out here rather than read from the object under test, so this is an
    oracle and not a transcription: ``INVTEMP_FACTOR / T`` is the documented
    meaning of ``invtemp`` (acceptance.py) and the criterion is the textbook
    Metropolis rule for a symmetric proposal.

    Parameters
    ----------
    temperature : float
        The simulation temperature.

    delta_energy : float
        new_energy - old_energy. Only positive values are meaningful here;
        downhill moves are accepted unconditionally.

    Returns
    -------
    float
        exp(-dE * INVTEMP_FACTOR / T).
    """
    return math.exp(-delta_energy * CONFIG.INVTEMP_FACTOR / temperature)


def _acceptance_frequency(ac: AcceptanceCalculator, delta_energy: float,
                          n_draws: int, seed: int) -> float:
    """Measured acceptance frequency of an uphill move over `n_draws` free draws.

    Parameters
    ----------
    ac : AcceptanceCalculator
        The object under test.

    delta_energy : float
        The (positive) energy change offered on every draw.

    n_draws : int
        Number of independent calls.

    seed : int
        Seed for the ``random`` module, so the result is deterministic.

    Returns
    -------
    float
        accepted / n_draws.
    """
    random.seed(seed)
    accepted = 0
    for _ in range(n_draws):
        if ac.boltzmann_acceptance(0.0, delta_energy):
            accepted += 1
    return accepted / n_draws


@pytest.mark.parametrize("temperature,delta_energy", _CASES,
                         ids=[f"T{int(t)}-dE{int(d)}" for t, d in _CASES])
def test_boltzmann_acceptance_frequency_matches_exp_minus_beta_dE(temperature, delta_energy):
    """The acceptance FREQUENCY must equal exp(-dE/T) within binomial error.

    This is the pin that the branch-structure tests could not be: it compares the
    exponent against a number. A 20 % error in the exponent shows up on every row
    of this grid at 5.1 to 31.2 sigma (see the positive control below); on the
    shipped code the worst row sits at 1.65 sigma.
    """
    ac = AcceptanceCalculator(temp=temperature, keyword_lookup=_uniform_moveset())
    expected = _exact_acceptance(temperature, delta_energy)
    measured = _acceptance_frequency(ac, delta_energy, _N_DRAWS,
                                     _row_seed(temperature, delta_energy))
    sigma = math.sqrt(expected * (1.0 - expected) / _N_DRAWS)
    z = (measured - expected) / sigma
    assert abs(z) <= _K_SIGMA, (
        f"boltzmann_acceptance(T={temperature}, dE={delta_energy}): measured acceptance "
        f"{measured:.6f} vs exact exp(-dE/T) = {expected:.6f}, i.e. {z:.1f} sigma over "
        f"{_N_DRAWS} draws - the Metropolis exponent does not match the temperature")


@pytest.mark.parametrize("temperature,delta_energy", _CASES,
                         ids=[f"T{int(t)}-dE{int(d)}" for t, d in _CASES])
def test_boltzmann_acceptance_threshold_is_exactly_exp_minus_beta_dE(monkeypatch,
                                                                    temperature,
                                                                    delta_energy):
    """The uniform draw is compared against exactly exp(-dE/T), not near it.

    Sharper than the frequency test and completely deterministic: a draw one part
    in 10^9 below the exact acceptance probability must be accepted and one the
    same distance above it must be rejected. This pins the threshold VALUE, so it
    also catches a wrong comparison operator or an off-by-a-constant exponent
    that a frequency test would need many draws to see.
    """
    ac = AcceptanceCalculator(temp=temperature, keyword_lookup=_uniform_moveset())
    expected = _exact_acceptance(temperature, delta_energy)

    monkeypatch.setattr(random, "random", lambda: expected * (1.0 - 1e-9))
    assert ac.boltzmann_acceptance(0.0, delta_energy) is True

    monkeypatch.setattr(random, "random", lambda: expected * (1.0 + 1e-9))
    assert ac.boltzmann_acceptance(0.0, delta_energy) is False


def test_boltzmann_acceptance_takes_downhill_moves_with_certainty(monkeypatch):
    """Downhill and energy-neutral moves are accepted for every possible draw.

    Forcing ``random.random`` to 1.0 - the value that rejects every uphill move -
    proves the downhill branch does not consult it at all, rather than the older
    formulation which monkeypatched the draw away and could not tell the two
    cases apart.
    """
    ac = AcceptanceCalculator(temp=40.0, keyword_lookup=_uniform_moveset())
    monkeypatch.setattr(random, "random", lambda: 1.0)

    for old, new in ((10.0, 9.0), (10.0, 10.0), (10.0, -50.0), (0.0, -1.0)):
        assert ac.boltzmann_acceptance(old, new) is True


def test_boltzmann_acceptance_follows_a_temperature_update():
    """update_temperature must move the acceptance frequency with it.

    Quenching and the TSMMC excursions both re-temper the calculator in place, so
    the criterion has to read the CURRENT temperature. Checking the frequency
    after the update catches a cached exponent that the invtemp attribute test
    cannot see.
    """
    ac = AcceptanceCalculator(temp=300.0, keyword_lookup=_uniform_moveset())
    ac.update_temperature(20.0)

    expected = _exact_acceptance(20.0, 10.0)
    measured = _acceptance_frequency(ac, 10.0, _N_DRAWS, _SEED)
    sigma = math.sqrt(expected * (1.0 - expected) / _N_DRAWS)
    assert abs(measured - expected) <= _K_SIGMA * sigma, (
        f"after update_temperature(20): measured {measured:.6f} vs exact {expected:.6f}")


def test_the_frequency_test_would_catch_a_twenty_percent_exponent_error():
    """Positive control: the grid above must REJECT a 1.2x exponent.

    A statistical test that has never been shown to fail is not evidence of
    anything. This runs the same comparison against a deliberately wrong
    criterion - the exact mutation that survived the whole suite - and requires
    it to breach the threshold on the majority of the grid, so the assertion
    above cannot quietly lose its power to a smaller sample or a milder grid.
    """
    class _MistemperedAcceptance:
        """boltzmann_acceptance with the exponent scaled by 1.2."""

        def __init__(self, temperature: float) -> None:
            self.invtemp = CONFIG.INVTEMP_FACTOR / temperature

        def boltzmann_acceptance(self, old_energy: float, new_energy: float) -> bool:
            if new_energy <= old_energy:
                return True
            return random.random() < math.exp(-(new_energy - old_energy)
                                              * self.invtemp * 1.2)

    caught = 0
    worst = 0.0
    for temperature, delta_energy in _CASES:
        broken = _MistemperedAcceptance(temperature)
        expected = _exact_acceptance(temperature, delta_energy)
        measured = _acceptance_frequency(broken, delta_energy, _N_DRAWS,
                                         _row_seed(temperature, delta_energy))
        sigma = math.sqrt(expected * (1.0 - expected) / _N_DRAWS)
        z = abs(measured - expected) / sigma
        worst = max(worst, z)
        if z > _K_SIGMA:
            caught += 1

    # measured on the shipped code: all 11 rows caught, worst row 31.2 sigma
    assert caught >= len(_CASES) - 2, (
        f"the acceptance grid only resolves a 20 % exponent error in {caught} of "
        f"{len(_CASES)} rows (worst |z| = {worst:.1f}) - it has lost its power")
    assert worst > 20.0, f"worst-row separation is only {worst:.1f} sigma"
