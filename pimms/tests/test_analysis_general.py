from datetime import datetime, timedelta

import pytest

from pimms import analysis_general


class _StubAcceptance:
    """Minimal stand-in for AcceptanceCalculator's throughput accessor."""

    def __init__(self, total=0):
        self._total = total

    def total_attempted_moves(self):
        return self._total


@pytest.fixture
def captured_performance(monkeypatch):
    """Capture what evaluate_performance would write, instead of touching disk."""
    written = []
    monkeypatch.setattr(analysis_general.analysis_IO, "write_performance",
                        lambda *args: written.append(args))
    monkeypatch.setattr(analysis_general.pimmslogger, "log_status", lambda *a, **k: None)
    return written


def test_evaluate_performance_survives_step_zero(captured_performance):
    """step=0 gives a zero step rate; the ETA must not divide by it.

    ``steps_per_second`` is ``step / elapsed``, so on step 0 it is exactly 0.0 and
    the time-remaining estimate used to raise ZeroDivisionError.
    """
    start = datetime.now() - timedelta(seconds=5)

    analysis_general.evaluate_performance(0, start, 100, 10, _StubAcceptance(0))

    assert len(captured_performance) == 1
    step, eq_string, _sps, _mps, _elapsed, remaining = captured_performance[0]
    assert step == 0
    assert eq_string == "E"                    # still inside equilibration
    assert remaining == "00:00:00"             # nothing to extrapolate from yet


def test_evaluate_performance_reports_production_phase_and_rates(captured_performance):
    start = datetime.now() - timedelta(seconds=10)

    analysis_general.evaluate_performance(50, start, 100, 10, _StubAcceptance(5000))

    step, eq_string, sps, mps, elapsed, remaining = captured_performance[0]
    assert step == 50
    assert eq_string == "P"                    # past equilibration
    assert sps == pytest.approx(5.0, rel=0.2)  # 50 steps in ~10 s
    assert mps == pytest.approx(500.0, rel=0.2)
    # both timing strings are hh:mm:ss
    assert len(elapsed.split(":")) == 3
    assert len(remaining.split(":")) == 3


def test_evaluate_performance_elapsed_time_does_not_wrap_after_a_day(captured_performance):
    """Elapsed time uses total_seconds(), not timedelta.seconds (which is < 1 day)."""
    start = datetime.now() - timedelta(days=1, hours=2, minutes=3, seconds=4)

    analysis_general.evaluate_performance(10, start, 100, 5, _StubAcceptance(10))

    elapsed = captured_performance[0][4]
    hours = int(elapsed.split(":")[0])
    assert hours >= 26
