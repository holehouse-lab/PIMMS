## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import random

from . import numpy_utils


# Developer sanity-check script, run by hand with `python -m pimms.<name>`.
# The body is guarded so that merely importing this module (it ships inside the
# package) does not execute a long side-effecting loop at import time.
if __name__ == "__main__":
    dim = 80
    cluster_move_threshold = 60

    results_count = {}

    for i in range(0,5000000):
    
        for d in [80,80,80]:
            r = numpy_utils.randneg(random.randint(1, min(d-1, cluster_move_threshold)))
    
            if r not in results_count:
                results_count[r] = 0

        results_count[r] = results_count[r] +1

    
