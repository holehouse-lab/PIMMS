import pytest

from pimms.chainTSMMC import TSMMC
from pimms.latticeExceptions import MoveException


@pytest.mark.parametrize("args", [
    (0, 120, "LINEAR", 1, 4, False),
    (100, float("nan"), "LINEAR", 1, 4, False),
    (100, 120, "CUBIC", 1, 4, False),
    (100, 120, "LINEAR", 0, 4, False),
    (100, 120, "LINEAR", 1, 0, False),
    (100, 120, "LINEAR", 1, 4, 0),
])
def test_tsmmc_rejects_invalid_schedule_inputs(args):
    with pytest.raises(MoveException):
        TSMMC(*args)


def test_tsmmc_accepts_case_insensitive_linear_mode():
    coordinator = TSMMC(100, 120, "linear", 1, 4, False)
    assert coordinator.mode == "LINEAR"


@pytest.mark.parametrize("work", [float("nan"), float("inf"), "bad", True])
def test_tsmmc_rejects_invalid_accumulated_work(work):
    coordinator = TSMMC(100, 120, "LINEAR", 1, 4, False)
    with pytest.raises(MoveException, match="work"):
        coordinator.accept_tempered_transition(work)
