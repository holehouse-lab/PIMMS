## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2020
## ...........................................................................

# Oct 2019: all of the imports below used to occur, but I cannot see why, suspect
#           this is legacy from an earlier version.

#import numpy as np
#cimport numpy as np
#cimport cython 
#import random
#import mega_crank

#DTYPE = np.int
#ctypedef np.int_t DTYPE_t


# Interface into the cstdlib to get the system
# specific maximum random number (RAND_MAX).
from libc.stdlib cimport RAND_MAX
def get_randmax():
    return RAND_MAX
