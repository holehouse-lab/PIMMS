## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................



import numpy as np

from . import inner_loops
from .CONFIG import NP_INT_TYPE


# Developer sanity-check script, run by hand with `python -m pimms.cython_testing`.
# The body is guarded so that merely importing this module (it ships inside the
# package) does not execute anything at import time. NB: the previous call was also
# missing the type_grid argument, so the script crashed with a TypeError if run.
if __name__ == "__main__":
    type_grid = np.zeros((200, 200, 200), dtype=NP_INT_TYPE)
    print(inner_loops.extract_SR_and_LR_pairs_from_position_3D(
        np.array([2, 3, 2], dtype=NP_INT_TYPE), 0, type_grid, 200, 200, 200))
