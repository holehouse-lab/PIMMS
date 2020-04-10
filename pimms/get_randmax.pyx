## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2020
## ...........................................................................


# Interface into the cstdlib to get the system
# specific maximum random number (RAND_MAX).
from libc.stdlib cimport RAND_MAX
def get_randmax():
    return RAND_MAX
