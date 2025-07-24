## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2024
## 
## ...........................................................................

cimport numpy as cnp
import numpy as np

def randint_np(int start, int end):
    return np.random.randint(start, end+1)


def seed_randint_np(int seedval):
    np.random.seed(seedval)


def generate_random_numbers(int min_val, int max_val, int size):    
    cdef cnp.ndarray[np.int64_t, ndim=1] random_numbers = np.random.randint(
        min_val, max_val, size, dtype=np.int64
    )
    return random_numbers 
