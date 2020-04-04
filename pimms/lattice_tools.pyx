## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2020
## 
## ...........................................................................


import numpy as np
cimport numpy as np
cimport cython 

DTYPE = np.int
ctypedef np.int_t DTYPE_t

## File that contains discrete, stateless functions
## for high performance lattice operations
##
##


@cython.boundscheck(False)
@cython.wraparound(False) 
def pbc_correct_3D(np.ndarray[DTYPE_t, ndim=1] posA, np.ndarray[DTYPE_t, ndim=1] posB, np.ndarray[DTYPE_t, ndim=1] DIM):
    """
    Function which performs relative PBC correction using positionA as the universal reference and
    repositioning B if necessary

    """
    
    cdef np.ndarray[np.int_t, ndim=1] newB = np.zeros((3), dtype=np.int)

    cdef int i;

    for i in xrange(0,3):
        if posA[i] - posB[i] > DIM[i]/2:
            newB[i] = posB[i] + DIM[i]
        elif posA[i] - posB[i] < -DIM[i]/2:
            newB[i] = posB[i] - DIM[i]
        else:
            newB[i] = posB[i]
            

    """
    print ">>>> IN"
    print posA
    print posB
    print newB
    print ""
    """


    return newB

        

    
                   
                   
    
