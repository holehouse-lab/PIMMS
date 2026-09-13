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

from cython.view cimport array

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
def extract_SR_and_LR_pairs_from_position_3D(NUMPY_INT_TYPE[:] position, 
                                             int LR_position, 
                                             NUMPY_INT_TYPE[:,:,:] type_grid,
                                             int XDIM, 
                                             int YDIM, 
                                             int ZDIM):
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

    Parameters
    ----------
    position : (3,) NUMPY_INT_TYPE (int32) memoryview
        The central (x, y, z) site the pairs are built around.
    LR_position : int (C int)
        0 to build short-range pairs only, 1 to also build the long-range and
        super-long-range pairs. Any other value raises.
    type_grid : (XDIM, YDIM, ZDIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid, used only to skip empty (value 0) neighbour sites when
        building the LR/SLR pairs, since those contribute nothing.
    XDIM : int (C int)
        Box size along x, used to wrap the neighbour coordinates.
    YDIM : int (C int)
        Box size along y.
    ZDIM : int (C int)
        Box size along z.

    Returns
    -------
    tuple of numpy.ndarray
        (SR_pairs, LR_pairs, SLR_pairs), each NUMPY_INT_TYPE (int32). SR_pairs
        is always the full (26, 2, 3) nearest-neighbour shell. With
        LR_position == 0 the other two are empty 1D arrays. Otherwise LR_pairs
        is (n_LR, 2, 3) drawn from the 5x5x5 shell outside the nearest
        neighbours (at most 98 pairs) and SLR_pairs is (n_SLR, 2, 3) drawn from
        the 7x7x7 shell outside that (at most 218), both trimmed to the
        occupied sites only.

    Raises
    ------
    InnerLoopException
        If LR_position is neither 0 nor 1.
    """

    # declare some variables
    cdef int SLR_index, SR_index, LR_index, x_off, y_off, z_off
    cdef int x_p, y_p, z_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 

    
    # initialize 
    # 3x3x3 neighborhood minus the central self-pair = 26 pairs
    SR_pairs = np.zeros((26,2,3), dtype=NUMPY_INT_TYPE_PYTHON)

    
    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    cdef int z = position[2]
    
    # short range only
    if LR_position == 0:

        SR_index = 0
        for x_off in xrange(-1,2):
            for y_off in xrange(-1,2):
                for z_off in xrange(-1,2):
                    if x_off == 0 and y_off == 0 and z_off == 0:
                        continue

                    x_p = pbc_correction(x + x_off, XDIM)
                    y_p = pbc_correction(y + y_off, YDIM)
                    z_p = pbc_correction(z + z_off, ZDIM)

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


        return (SR_pairs, np.array([], dtype=NUMPY_INT_TYPE_PYTHON), np.array([], dtype=NUMPY_INT_TYPE_PYTHON))

    elif LR_position == 1:
        LR_pairs  = np.zeros((98, 2, 3), dtype=NUMPY_INT_TYPE_PYTHON)
        SLR_pairs = np.zeros((218, 2, 3), dtype=NUMPY_INT_TYPE_PYTHON)
                
        SR_index  = 0
        LR_index  = 0
        SLR_index = 0
        

        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                for z_off in xrange(-3,4):
                    x_p = pbc_correction(x + x_off, XDIM)
                    y_p = pbc_correction(y + y_off, YDIM)
                    z_p = pbc_correction(z + z_off, ZDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2 and abs(z_off) <2:
                        if x_off == 0 and y_off == 0 and z_off == 0:
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
def extract_SR_and_LR_pairs_from_position_2D(NUMPY_INT_TYPE[:] position, 
                                             int LR_position, 
                                             NUMPY_INT_TYPE[:,:] type_grid,
                                             int XDIM, 
                                             int YDIM):
    """
    2D counterpart of extract_SR_and_LR_pairs_from_position_3D: takes a single
    position ($position) and the type_grid and determines the set of pairwise
    interactions between that central position and the positions around it.

    The pairs are inherently numbered (i.e. [A-B] would be A then B). The
    ordering is by the sign of the (pre-PBC) offset from the central
    position: if the first non-zero component of the offset is positive the
    NEIGHBOUR comes first, otherwise the central position comes first. This rule
    is antisymmetric, so the same physical pair seen from either of its two ends
    is ordered identically (which is what makes downstream de-duplication
    correct) - note it is NOT a sort by the numeric coordinate values (across a
    periodic face the two disagree).

    Parameters
    ----------
    position : (2,) NUMPY_INT_TYPE (int32) memoryview
        The central (x, y) site the pairs are built around.
    LR_position : int (C int)
        0 to build short-range pairs only, 1 to also build the long-range and
        super-long-range pairs. Any other value raises.
    type_grid : (XDIM, YDIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid, used only to skip empty (value 0) neighbour sites when
        building the LR/SLR pairs, since those contribute nothing.
    XDIM : int (C int)
        Box size along x, used to wrap the neighbour coordinates.
    YDIM : int (C int)
        Box size along y.

    Returns
    -------
    tuple of numpy.ndarray
        (SR_pairs, LR_pairs, SLR_pairs), each NUMPY_INT_TYPE (int32). SR_pairs
        is the full (8, 2, 2) nearest-neighbour shell. With LR_position == 0 the
        other two are empty 1D arrays. Otherwise LR_pairs is (n_LR, 2, 2) drawn
        from the 5x5 shell outside the nearest neighbours (at most 16 pairs) and
        SLR_pairs is (n_SLR, 2, 2) drawn from the 7x7 shell outside that (at
        most 24), both trimmed to the occupied sites only.

    Raises
    ------
    InnerLoopException
        If LR_position is neither 0 nor 1.
    """
    
    # declare some variables
    cdef int SR_index, LR_index, SLR_index, x_off, y_off
    cdef int x_p, y_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((8,2,2), dtype=NUMPY_INT_TYPE_PYTHON)
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]

    SR_index = 0    
    # short range only
    if LR_position == 0:


        for x_off in xrange(-1,2):
            for y_off in xrange(-1,2):
                    if x_off == 0 and y_off == 0:
                        continue

                    x_p = pbc_correction(x + x_off, XDIM)
                    y_p = pbc_correction(y + y_off, YDIM)

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

        # return is SR, LR, SLR
        return (SR_pairs, np.array([], dtype=NUMPY_INT_TYPE_PYTHON), np.array([], dtype=NUMPY_INT_TYPE_PYTHON))    

    elif LR_position == 1:

        LR_pairs = np.zeros((16, 2, 2), dtype=NUMPY_INT_TYPE_PYTHON)  # 5 x 5 -  3 x 3       
        SLR_pairs = np.zeros((24, 2, 2), dtype=NUMPY_INT_TYPE_PYTHON) # 7 x 7 -  5 x 5
        
        LR_index = 0
        SLR_index = 0

        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                    x_p = pbc_correction(x + x_off, XDIM)
                    y_p = pbc_correction(y + y_off, YDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2:
                        if x_off == 0 and y_off == 0:
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
def extract_LR_pairs_from_position_3D(NUMPY_INT_TYPE[:] position, 
                                      int LR_position, 
                                      NUMPY_INT_TYPE[:,:,:] type_grid,
                                      int XDIM, 
                                      int YDIM, 
                                      int ZDIM):
    """
    Same as extract_all except ONLY returns LR and SLR pairs

    Pair ordering is the same antisymmetric offset-sign rule used by
    extract_SR_and_LR_pairs_from_position_3D, so pairs from this function can be
    de-duplicated against each other.

    Parameters
    ----------
    position : (3,) NUMPY_INT_TYPE (int32) memoryview
        The central (x, y, z) site the pairs are built around.
    LR_position : int (C int)
        0 to return empty arrays (no long-range interactions in play), 1 to
        build the LR and SLR pairs. Any other value raises.
    type_grid : (XDIM, YDIM, ZDIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid, used to skip empty (value 0) neighbour sites.
    XDIM : int (C int)
        Box size along x, used to wrap the neighbour coordinates.
    YDIM : int (C int)
        Box size along y.
    ZDIM : int (C int)
        Box size along z.

    Returns
    -------
    tuple of numpy.ndarray
        (LR_pairs, SLR_pairs), each NUMPY_INT_TYPE (int32). With
        LR_position == 0 both are empty 1D arrays. Otherwise LR_pairs is
        (n_LR, 2, 3) from the 5x5x5 shell outside the nearest neighbours (at
        most 98) and SLR_pairs is (n_SLR, 2, 3) from the 7x7x7 shell outside
        that (at most 218), both trimmed to the occupied sites only.

    Raises
    ------
    InnerLoopException
        If LR_position is neither 0 nor 1.
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
                    x_p = pbc_correction(x + x_off, XDIM)
                    y_p = pbc_correction(y + y_off, YDIM)
                    z_p = pbc_correction(z + z_off, ZDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2 and abs(z_off) <2:
                        continue

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3 and abs(z_off) < 3:
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
def extract_SR_pairs_from_position_3D(NUMPY_INT_TYPE[:] position,                 
                                          int XDIM, 
                                          int YDIM, 
                                          int ZDIM):
    """
    Returns the non-redundant set of pairs associated with the 3D position defined
    by the position array and all possible short-range interaction sites. Returned
    ordering is by the sign of the (pre-PBC) offset from the central
    position: if the first non-zero component of the offset is positive the
    NEIGHBOUR comes first, otherwise the central position comes first. This rule
    is antisymmetric, so the same physical pair seen from either of its two ends
    is ordered identically (which is what makes downstream de-duplication
    correct) - note it is NOT a sort by the numeric coordinate values (across a
    periodic face the two disagree).

    No type_grid is consulted here, so every one of the 26 neighbour sites is
    returned whether or not it is occupied.

    Parameters
    ----------
    position : (3,) NUMPY_INT_TYPE (int32) memoryview
        The central (x, y, z) site the pairs are built around.
    XDIM : int (C int)
        Box size along x, used to wrap the neighbour coordinates.
    YDIM : int (C int)
        Box size along y.
    ZDIM : int (C int)
        Box size along z.

    Returns
    -------
    (26, 2, 3) numpy.ndarray, NUMPY_INT_TYPE (int32)
        The full nearest-neighbour shell, one pair per surrounding site.
    """
    
    # declare some variables
    cdef int SR_index, x_off, y_off, z_off
    cdef int x_p, y_p, z_p
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((26,2,3), dtype=NUMPY_INT_TYPE_PYTHON)

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    cdef int z = position[2]
    
    SR_index = 0
    for x_off in xrange(-1,2):
        for y_off in xrange(-1,2):
            for z_off in xrange(-1,2):
                if x_off == 0 and y_off == 0 and z_off == 0:
                    continue

                x_p = pbc_correction(x + x_off, XDIM)
                y_p = pbc_correction(y + y_off, YDIM)
                z_p = pbc_correction(z + z_off, ZDIM)

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

    return (SR_pairs)


##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_LR_pairs_from_position_2D(NUMPY_INT_TYPE[:] position, 
                                      int LR_position, 
                                      NUMPY_INT_TYPE[:,:] type_grid,
                                      int XDIM, 
                                      int YDIM):
                                             
    """
    2D counterpart of extract_LR_pairs_from_position_3D: same as extract_all
    except ONLY returns LR and SLR pairs.

    Pair ordering is the same antisymmetric offset-sign rule used by
    extract_SR_and_LR_pairs_from_position_2D, so pairs from this function can be
    de-duplicated against each other.

    Parameters
    ----------
    position : (2,) NUMPY_INT_TYPE (int32) memoryview
        The central (x, y) site the pairs are built around.
    LR_position : int (C int)
        0 to return empty arrays (no long-range interactions in play), 1 to
        build the LR and SLR pairs. Any other value raises.
    type_grid : (XDIM, YDIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid, used to skip empty (value 0) neighbour sites.
    XDIM : int (C int)
        Box size along x, used to wrap the neighbour coordinates.
    YDIM : int (C int)
        Box size along y.

    Returns
    -------
    tuple of numpy.ndarray
        (LR_pairs, SLR_pairs), each NUMPY_INT_TYPE (int32). With
        LR_position == 0 both are empty 1D arrays. Otherwise LR_pairs is
        (n_LR, 2, 2) from the 5x5 shell outside the nearest neighbours (at most
        16) and SLR_pairs is (n_SLR, 2, 2) from the 7x7 shell outside that (at
        most 24), both trimmed to the occupied sites only.

    Raises
    ------
    InnerLoopException
        If LR_position is neither 0 nor 1.
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
                    x_p = pbc_correction(x + x_off, XDIM)
                    y_p = pbc_correction(y + y_off, YDIM)
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2:
                        continue

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3:
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
def extract_SR_pairs_from_position_2D(NUMPY_INT_TYPE[:] position,                                              
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

    No type_grid is consulted here, so every one of the 8 neighbour sites is
    returned whether or not it is occupied.

    Parameters
    ----------
    position : (2,) NUMPY_INT_TYPE (int32) memoryview
        The central (x, y) site the pairs are built around.
    XDIM : int (C int)
        Box size along x, used to wrap the neighbour coordinates.
    YDIM : int (C int)
        Box size along y.

    Returns
    -------
    (8, 2, 2) numpy.ndarray, NUMPY_INT_TYPE (int32)
        The full nearest-neighbour shell, one pair per surrounding site.
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

            x_p = pbc_correction(x + x_off, XDIM)
            y_p = pbc_correction(y + y_off, YDIM)

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

    return (SR_pairs)
        
    
    

@cython.cdivision(True)
cdef int pbc_correction(int value, int DIM):    
    """
    Performs intelligent periodic boundary correction
    which FIRST checks to see if we have a negative
    value and IF NOT uses the % operator - this means
    we can use % without checking the sign giving
    a 35% speedup per call

    The negative branch adds a single box length rather than taking a modulo, so
    it is only correct for value >= -DIM. That holds throughout this module,
    where value is always an in-box coordinate offset by at most 3 sites.

    Parameters
    ----------
    value : int (C int)
        Candidate coordinate along one axis, i.e. an in-box coordinate plus a
        small neighbour offset.
    DIM : int (C int)
        Box size along that axis.

    Returns
    -------
    int (C int)
        The coordinate wrapped back into [0, DIM).
    """
                            
    if value < 0:
        return DIM+value
    else:
        return (value % DIM)
