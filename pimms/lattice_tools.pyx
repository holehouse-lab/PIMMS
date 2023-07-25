## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2023
## 
## ...........................................................................


import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython 

ctypedef cnp.int_t NUMPY_INT_TYPE


## File that contains discrete, stateless functions
## for high performance lattice operations
##
##


@cython.boundscheck(False)
@cython.wraparound(False) 
def pbc_correct_3D(NUMPY_INT_TYPE[:] posA, NUMPY_INT_TYPE[:] posB, NUMPY_INT_TYPE[:] DIM):
    """
    Function which performs relative PBC correction using positionA as the universal reference and
    repositioning B if necessary

    """
    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] newB = np.zeros((3), dtype=int)

    cdef int i;

    for i in range(0,3):
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

        

    
                   
                   
    
