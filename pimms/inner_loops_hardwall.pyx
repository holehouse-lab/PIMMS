## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2026
## 
## ...........................................................................

import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython 

from pimms.latticeExceptions import InnerLoopException

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b


from pimms.cython_config cimport NUMPY_INT_TYPE
from pimms.CONFIG import NP_INT_TYPE as NUMPY_INT_TYPE_PYTHON

#from numpy cimport int16_t as NUMPY_INT16_TYPE
#ctypedef NUMPY_INT16_TYPE  NUMPY_INT_TYPE

## inner_loops contains functions for geting positions and bead information in 
## the local 2D or 3D environment
##
##
##



##
#################################################################################################
##


@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_SR_and_LR_pairs_from_position_3D_hardwall(NUMPY_INT_TYPE[:] position, 
                                                      NUMPY_INT_TYPE LR_position, 
                                                      NUMPY_INT_TYPE[:,:,:] type_grid,
                                                      NUMPY_INT_TYPE XDIM, 
                                                      NUMPY_INT_TYPE YDIM, 
                                                      NUMPY_INT_TYPE ZDIM):
    """
    Function that takes a single position ($position) and the type_grid and determines
    the set of pairwise interactions between that central position and the positions
    around it. 

    The pairs are inherently numbered (i.e. [A-B] would be A then B). The
    ordering is by the sign of the (pre-PBC) offset from the central
    position: if the first non-zero component of the offset is positive the
    NEIGHBOUR comes first, otherwise the central position comes first. This rule
    is antisymmetric, so the same physical pair seen from either of its two ends
    is ordered identically (which is what makes downstream de-duplication
    correct) - note it is NOT a sort by the numeric coordinate values (across a
    periodic face the two disagree).

    """

    ## Note all cdef have to happen at the start for C scoping reasons

    # declare some variables
    #cdef int SLR_index, SR_index, LR_index, x_off, y_off, z_off;
    #cdef int x_tmp, y_tmp, z_tmp;

    cdef int SLR_index, SR_index, LR_index, x_off, y_off, z_off
    cdef NUMPY_INT_TYPE x_p, y_p, z_p

    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs

    # first set the central x, y and z positions
    cdef NUMPY_INT_TYPE x = position[0]
    cdef NUMPY_INT_TYPE y = position[1]
    cdef NUMPY_INT_TYPE z = position[2]
    
    # short range only
    if LR_position == 0:

        return (extract_SR_pairs_from_position_3D_hardwall(position, XDIM, YDIM, ZDIM), np.array([], dtype=NUMPY_INT_TYPE_PYTHON), np.array([], dtype=NUMPY_INT_TYPE_PYTHON))

        # (a triple-quoted block of superseded legacy code used to sit here, AFTER
        # the return - Cython evaluated it as an unreachable string expression and
        # the compiler flagged it with -Wunreachable-code; removed)

    elif LR_position == 1:
        
    
        SR_pairs  = np.zeros((26,2,3), dtype=NUMPY_INT_TYPE_PYTHON)
        LR_pairs  = np.zeros((98, 2, 3), dtype=NUMPY_INT_TYPE_PYTHON)
        SLR_pairs = np.zeros((218, 2, 3), dtype=NUMPY_INT_TYPE_PYTHON)
                            
        SR_index  = 0
        LR_index  = 0
        SLR_index = 0
        
        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                for z_off in xrange(-3,4):
                    x_p = pbc_hardwall(x + x_off, XDIM)
                    y_p = pbc_hardwall(y + y_off, YDIM)
                    z_p = pbc_hardwall(z + z_off, ZDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2 and abs(z_off) <2:
                        if x_off == 0 and y_off == 0 and z_off == 0:
                            continue
                        if x_p == -1 or y_p == -1 or z_p == -1:
                            continue
                        
                        # if x_off > 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y
                            SR_pairs[SR_index, 1, 2] = z

                            SR_pairs[SR_index, 0, 0] = x_p
                            SR_pairs[SR_index, 0, 1] = y_p
                            SR_pairs[SR_index, 0, 2] = z_p

                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y
                            SR_pairs[SR_index, 1, 2] = z

                            SR_pairs[SR_index, 0, 0] = x_p
                            SR_pairs[SR_index, 0, 1] = y_p
                            SR_pairs[SR_index, 0, 2] = z_p

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y
                            SR_pairs[SR_index, 1, 2] = z
                            
                            SR_pairs[SR_index, 0, 0] = x_p
                            SR_pairs[SR_index, 0, 1] = y_p
                            SR_pairs[SR_index, 0, 2] = z_p

                        else:
                            SR_pairs[SR_index, 0, 0] = x
                            SR_pairs[SR_index, 0, 1] = y
                            SR_pairs[SR_index, 0, 2] = z
                            
                            SR_pairs[SR_index, 1, 0] = x_p
                            SR_pairs[SR_index, 1, 1] = y_p
                            SR_pairs[SR_index, 1, 2] = z_p
                            
                        SR_index = SR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3 and abs(z_off) < 3:
                        if x_p == -1 or y_p == -1 or z_p == -1:
                            continue
                        if type_grid[x_p, y_p, z_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            LR_pairs[LR_index, 0, 2] = z_p                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z

                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            LR_pairs[LR_index, 0, 2] = z_p

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            LR_pairs[LR_index, 0, 2] = z_p

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            LR_pairs[LR_index, 0, 2] = z
                            
                            LR_pairs[LR_index, 1, 0] = x_p
                            LR_pairs[LR_index, 1, 1] = y_p
                            LR_pairs[LR_index, 1, 2] = z_p
                            
                        LR_index = LR_index+1

                    # SUPER LONG RANGE INTERACTIONS...
                    else:
                        if x_p == -1 or y_p == -1 or z_p == -1:
                            continue
                        if type_grid[x_p, y_p, z_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            SLR_pairs[SLR_index, 0, 2] = z_p                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z

                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            SLR_pairs[SLR_index, 0, 2] = z_p

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            SLR_pairs[SLR_index, 0, 2] = z_p

                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            SLR_pairs[SLR_index, 0, 2] = z
                            
                            SLR_pairs[SLR_index, 1, 0] = x_p
                            SLR_pairs[SLR_index, 1, 1] = y_p
                            SLR_pairs[SLR_index, 1, 2] = z_p
                            
                        SLR_index = SLR_index+1


        return (SR_pairs[0:SR_index], LR_pairs[0:LR_index], SLR_pairs[0:SLR_index])
        
                        


    else:
        raise InnerLoopException('Invalid LR option passed')


@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_SR_and_LR_pairs_from_position_2D_hardwall(NUMPY_INT_TYPE[:] position, 
                                                      int LR_position, 
                                                      NUMPY_INT_TYPE[:,:] type_grid,
                                                      int XDIM, 
                                                      int YDIM):
    """
    Extracts the short-range and long-range pairs from a given position in the
    type grid.  This is a 2D version of the function above.

    Parameters
    ----------
    position : NUMPY_INT_TYPE[:]
        The position in the type grid from which to extract the pairs.

    LR_position : int
        Interaction mode. 0 extracts the short-range pairs only; 1 extracts
        short-range, long-range AND super-long-range pairs. Any other value
        raises InnerLoopException (there is no separate "both" mode 2).

    type_grid : NUMPY_INT_TYPE[:,:]
        The type grid.

    XDIM : int
        The x dimension of the type grid.

    YDIM : int
        The y dimension of the type grid.

    Returns
    -------
    SR_pairs : cnp.ndarray[NUMPY_INT_TYPE, ndim=3]
        The short-range pairs.

    """

    # declare some variables
    cdef int SR_index, LR_index, SLR_index, x_off, y_off
    cdef int x_p, y_p

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
                
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs

    
    # short range only
    if LR_position == 0:

        return (extract_SR_pairs_from_position_2D_hardwall(position, XDIM, YDIM), np.array([], dtype=NUMPY_INT_TYPE_PYTHON), np.array([], dtype=NUMPY_INT_TYPE_PYTHON))    
        

    elif LR_position == 1:


        SR_pairs  = np.zeros((8,  2, 2), dtype=NUMPY_INT_TYPE_PYTHON)
        LR_pairs  = np.zeros((16, 2, 2), dtype=NUMPY_INT_TYPE_PYTHON) # 5 x 5 -  3 x 3
        SLR_pairs = np.zeros((24, 2, 2), dtype=NUMPY_INT_TYPE_PYTHON) # 7 x 7 -  5 x 5
        
        SR_index = 0    
        LR_index = 0
        SLR_index = 0

        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                    x_p = pbc_hardwall(x + x_off, XDIM)
                    y_p = pbc_hardwall(y + y_off, YDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2:
                        if x_off == 0 and y_off == 0:
                            continue
                        if x_p == -1 or y_p == -1:
                            continue
                        
                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y

                            SR_pairs[SR_index, 0, 0] = x_p
                            SR_pairs[SR_index, 0, 1] = y_p
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y

                            SR_pairs[SR_index, 0, 0] = x_p
                            SR_pairs[SR_index, 0, 1] = y_p


                        else:
                            SR_pairs[SR_index, 0, 0] = x
                            SR_pairs[SR_index, 0, 1] = y
                            
                            SR_pairs[SR_index, 1, 0] = x_p
                            SR_pairs[SR_index, 1, 1] = y_p
                            
                        SR_index = SR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3:
                        if x_p == -1 or y_p == -1:
                            continue
                        if type_grid[x_p, y_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y                            
                            
                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p

                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y

                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            
                            LR_pairs[LR_index, 1, 0] = x_p
                            LR_pairs[LR_index, 1, 1] = y_p
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # SUPER LONG RANGE INTERACTIONS...
                    else:
                        if x_p == -1 or y_p == -1:
                            continue
                        if type_grid[x_p, y_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            
                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y

                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                
                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            
                            SLR_pairs[SLR_index, 1, 0] = x_p
                            SLR_pairs[SLR_index, 1, 1] = y_p
                            
                        SLR_index = SLR_index+1

        return (SR_pairs[0:SR_index], LR_pairs[0:LR_index], SLR_pairs[0:SLR_index])

                                
    else:
        raise InnerLoopException('Invalid LR option passed')
        

@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_LR_pairs_from_position_3D_hardwall(NUMPY_INT_TYPE[:] position, 
                                      int LR_position, 
                                      NUMPY_INT_TYPE[:,:,:] type_grid,
                                      int XDIM, 
                                      int YDIM, 
                                      int ZDIM):
    """
    Same as extract_all except ONLY returns LR and SLR pairs


    """
    
    # declare some variables
    cdef int LR_index, SLR_index, x_off, y_off, z_off
    cdef int x_p, y_p, z_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    cdef int z = position[2]
    
    # if no long-range interactions required the return empy 
    # arrays
    if LR_position == 0:
        return (np.array([], dtype=NUMPY_INT_TYPE_PYTHON), np.array([], dtype=NUMPY_INT_TYPE_PYTHON))

    elif LR_position == 1:
        LR_pairs = np.zeros((98, 2, 3), dtype=NUMPY_INT_TYPE_PYTHON)                
        SLR_pairs = np.zeros((218, 2, 3), dtype=NUMPY_INT_TYPE_PYTHON)                
        
        LR_index = 0
        SLR_index = 0
        
        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                for z_off in xrange(-3,4):
                    x_p = pbc_hardwall(x + x_off, XDIM)
                    y_p = pbc_hardwall(y + y_off, YDIM)
                    z_p = pbc_hardwall(z + z_off, ZDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2 and abs(z_off) <2:
                        continue

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3 and abs(z_off) < 3:
                        if x_p == -1 or y_p == -1 or z_p == -1:
                            continue
                        if type_grid[x_p, y_p, z_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            LR_pairs[LR_index, 0, 2] = z_p                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z

                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            LR_pairs[LR_index, 0, 2] = z_p

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            LR_pairs[LR_index, 0, 2] = z_p

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            LR_pairs[LR_index, 0, 2] = z
                            
                            LR_pairs[LR_index, 1, 0] = x_p
                            LR_pairs[LR_index, 1, 1] = y_p
                            LR_pairs[LR_index, 1, 2] = z_p
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Super long range interaction
                    else:
                        if x_p == -1 or y_p == -1 or z_p == -1:
                            continue
                        if type_grid[x_p, y_p, z_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            SLR_pairs[SLR_index, 0, 2] = z_p                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z

                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            SLR_pairs[SLR_index, 0, 2] = z_p

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            SLR_pairs[SLR_index, 0, 2] = z_p

                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            SLR_pairs[SLR_index, 0, 2] = z
                            
                            SLR_pairs[SLR_index, 1, 0] = x_p
                            SLR_pairs[SLR_index, 1, 1] = y_p
                            SLR_pairs[SLR_index, 1, 2] = z_p
                            
                        SLR_index = SLR_index+1

        return (LR_pairs[0:LR_index], SLR_pairs[0:SLR_index])

    else:
        raise InnerLoopException('Invalid LR option passed')


##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_SR_pairs_from_position_3D_hardwall(NUMPY_INT_TYPE[:] position,                 
                                               NUMPY_INT_TYPE XDIM, 
                                               NUMPY_INT_TYPE YDIM, 
                                               NUMPY_INT_TYPE ZDIM):
    """


    """

    # NB: these MUST be plain C ints. Declaring them as NUMPY_INT_TYPE (npy_int32) made Cython
    # route every abs(x_off) in the shell loops through Python's number protocol (boxing to a
    # Python int and back, 343 times per call) - the hardwall extractor ran ~5x slower than its
    # PBC twin for no reason. Cython only inlines abs() for C int/long/double.
    cdef int SR_index, x_off, y_off, z_off
    cdef NUMPY_INT_TYPE x_p, y_p, z_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((26,2,3), dtype=NUMPY_INT_TYPE_PYTHON)

    # first set the central x, y and z positions
    cdef NUMPY_INT_TYPE x = position[0]
    cdef NUMPY_INT_TYPE y = position[1]
    cdef NUMPY_INT_TYPE z = position[2]

    
    SR_index = 0
    
    for x_off in xrange(-1,2):
        for y_off in xrange(-1,2):
            for z_off in xrange(-1,2):
                if x_off == 0 and y_off == 0 and z_off == 0:
                    continue

                x_p = pbc_hardwall(x + x_off, XDIM)
                y_p = pbc_hardwall(y + y_off, YDIM)
                z_p = pbc_hardwall(z + z_off, ZDIM)

                if x_p == -1 or y_p == -1 or z_p == -1:
                    continue

                # if x_off < 0 then the non-central position must come first in the pair
                if x_off > 0:
                    SR_pairs[SR_index, 1, 0] = x
                    SR_pairs[SR_index, 1, 1] = y
                    SR_pairs[SR_index, 1, 2] = z

                    SR_pairs[SR_index, 0, 0] = x_p
                    SR_pairs[SR_index, 0, 1] = y_p
                    SR_pairs[SR_index, 0, 2] = z_p

                    
                # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                elif x_off == 0 and y_off > 0:
                    SR_pairs[SR_index, 1, 0] = x
                    SR_pairs[SR_index, 1, 1] = y
                    SR_pairs[SR_index, 1, 2] = z

                    SR_pairs[SR_index, 0, 0] = x_p
                    SR_pairs[SR_index, 0, 1] = y_p
                    SR_pairs[SR_index, 0, 2] = z_p

                # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                elif x_off == 0 and y_off == 0 and z_off > 0:
                    SR_pairs[SR_index, 1, 0] = x
                    SR_pairs[SR_index, 1, 1] = y
                    SR_pairs[SR_index, 1, 2] = z

                    SR_pairs[SR_index, 0, 0] = x_p
                    SR_pairs[SR_index, 0, 1] = y_p
                    SR_pairs[SR_index, 0, 2] = z_p

                else:
                    SR_pairs[SR_index, 0, 0] = x
                    SR_pairs[SR_index, 0, 1] = y
                    SR_pairs[SR_index, 0, 2] = z

                    SR_pairs[SR_index, 1, 0] = x_p
                    SR_pairs[SR_index, 1, 1] = y_p
                    SR_pairs[SR_index, 1, 2] = z_p

                SR_index = SR_index+1

                
    return SR_pairs[0:SR_index]


##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_LR_pairs_from_position_2D_hardwall(NUMPY_INT_TYPE[:] position, 
                                               int LR_position, 
                                               NUMPY_INT_TYPE[:,:] type_grid,
                                               int XDIM, 
                                               int YDIM):
                                             
    """


    """
    
    # declare some variables
    cdef int LR_index, SLR_index, x_off, y_off
    cdef int x_p, y_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    
    # if no long-range interactions required the return empy 
    # array
    if LR_position == 0:
        return (np.array([], dtype=NUMPY_INT_TYPE_PYTHON), np.array([], dtype=NUMPY_INT_TYPE_PYTHON))

    elif LR_position == 1:
        LR_pairs = np.zeros((16, 2, 2), dtype=NUMPY_INT_TYPE_PYTHON)        
        SLR_pairs = np.zeros((24, 2, 2), dtype=NUMPY_INT_TYPE_PYTHON)        

        LR_index = 0
        SLR_index = 0
        
        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                    x_p = pbc_hardwall(x + x_off, XDIM)
                    y_p = pbc_hardwall(y + y_off, YDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2:
                        continue

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3:
                        if x_p == -1 or y_p == -1:
                            continue
                        if type_grid[x_p, y_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y                            
                            
                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y

                            LR_pairs[LR_index, 0, 0] = x_p
                            LR_pairs[LR_index, 0, 1] = y_p


                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            
                            LR_pairs[LR_index, 1, 0] = x_p
                            LR_pairs[LR_index, 1, 1] = y_p
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Super long range interaction
                    else:
                        if x_p == -1 or y_p == -1:
                            continue
                        if type_grid[x_p, y_p] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y                            
                            
                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y

                            SLR_pairs[SLR_index, 0, 0] = x_p
                            SLR_pairs[SLR_index, 0, 1] = y_p


                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            
                            SLR_pairs[SLR_index, 1, 0] = x_p
                            SLR_pairs[SLR_index, 1, 1] = y_p
                            
                        SLR_index = SLR_index+1
                        

        return (LR_pairs[0:LR_index], SLR_pairs[0:SLR_index])

    else:
        raise InnerLoopException('Invalid LR option passed')

##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_SR_pairs_from_position_2D_hardwall(NUMPY_INT_TYPE[:] position,
                                               int XDIM, 
                                               int YDIM):
    """
    Returns the non-redundant set of pairs associated with the 2D position defined
    by the position array and all possible short-range interaction sites. Returned
    ordering is by the sign of the (pre-PBC) offset from the central
    position: if the first non-zero component of the offset is positive the
    NEIGHBOUR comes first, otherwise the central position comes first. This rule
    is antisymmetric, so the same physical pair seen from either of its two ends
    is ordered identically (which is what makes downstream de-duplication
    correct) - note it is NOT a sort by the numeric coordinate values (across a
    periodic face the two disagree).


    """
    
    # declare some variables
    cdef int SR_index, x_off, y_off
    cdef int x_p, y_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((8,2,2), dtype=NUMPY_INT_TYPE_PYTHON)

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    
    SR_index = 0
    for x_off in xrange(-1,2):
        for y_off in xrange(-1,2):
            if x_off == 0 and y_off == 0:
                continue

            x_p = pbc_hardwall(x + x_off, XDIM)
            y_p = pbc_hardwall(y + y_off, YDIM)

            if x_p == -1 or y_p == -1:
                continue

            # if x_off < 0 then the non-central position must come first in the pair
            if x_off > 0:
                SR_pairs[SR_index, 1, 0] = x
                SR_pairs[SR_index, 1, 1] = y

                SR_pairs[SR_index, 0, 0] = x_p
                SR_pairs[SR_index, 0, 1] = y_p


            # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
            elif x_off == 0 and y_off > 0:
                SR_pairs[SR_index, 1, 0] = x
                SR_pairs[SR_index, 1, 1] = y

                SR_pairs[SR_index, 0, 0] = x_p
                SR_pairs[SR_index, 0, 1] = y_p

            else:
                SR_pairs[SR_index, 0, 0] = x
                SR_pairs[SR_index, 0, 1] = y

                SR_pairs[SR_index, 1, 0] = x_p
                SR_pairs[SR_index, 1, 1] = y_p

            SR_index = SR_index+1

    return SR_pairs[0:SR_index]


@cython.cdivision(True)
cdef NUMPY_INT_TYPE pbc_hardwall(NUMPY_INT_TYPE value, NUMPY_INT_TYPE DIM):    
    """
    Takes an offset position (value) and the max dimensions in that axis (DIM) and if
    this new position is outside the lattice returns -1 (i.e. this is one part of a pair
    that will straddle the periodic boundary)
    """
                            
    if (value < 0) or value > (DIM-1):
        return -1
    else:
        return value
