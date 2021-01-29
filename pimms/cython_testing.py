## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2021
## ...........................................................................



import numpy as np

from . import inner_loops


# 10 M steps...
inner_loops.extract_SR_and_LR_pairs_from_position_3D(np.array([2,3,2]), 0, 200,200,200)

