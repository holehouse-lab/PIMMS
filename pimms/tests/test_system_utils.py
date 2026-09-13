import numpy as np
import pytest

from pimms import system_utils
from pimms.latticeExceptions import SimulationException


def test_grid_integer_capacity_depends_on_chain_count_not_bead_count(monkeypatch):
    monkeypatch.setattr(system_utils, "NUMPY_INT_TYPE_PYTHON", np.int8)

    # One very long chain only writes chainID 1 into the occupancy grid.
    system_utils.check_beads_to_grid_mapping([[1, "A" * 1000]])

    # By contrast, 128 chains require chainID 128, outside signed-int8 range.
    with pytest.raises(SimulationException, match="number of chains.*128"):
        system_utils.check_beads_to_grid_mapping([[128, "A"]])
