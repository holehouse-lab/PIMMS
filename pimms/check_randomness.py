## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2020
## ...........................................................................


import time
import random

from . import mega_crank
from . import numpy_utils

dim = 80
cluster_move_threshold = 60

results_count = {}

for i in range(0,5000000):
    
    for d in [80,80,80]:
        r = numpy_utils.randneg(random.randint(1, min(d-1, cluster_move_threshold)))
    
        if r not in results_count:
            results_count[r] = 0

    results_count[r] = results_count[r] +1

    
