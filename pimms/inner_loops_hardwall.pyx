## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2024
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

    The pairs are inherently numbered (i.e. [A-B] would be A then B). To determine which
    of the two positions is first in the pair we use the numerical value of the x/y/z 
    positions

    """

    ## Note all cdef have to happen at the start for C scoping reasons

    # declare some variables
    #cdef int SLR_index, SR_index, LR_index, x_off, y_off, z_off;
    #cdef int x_tmp, y_tmp, z_tmp;

    cdef NUMPY_INT_TYPE SLR_index, SR_index, LR_index, x_off, y_off, z_off
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


        # old code that seems to just re-implement the extract_SR_pairs_from_position_3D_hardwall function
        """
        SR_index = 0
        for x_off in xrange(-1,2):
            for y_off in xrange(-1,2):
                for z_off in xrange(-1,2):

                    # if x_off > 0 then the non-central position must come first in the pair
                    if x_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y
                        SR_pairs[SR_index, 1, 2] = z

                        SR_pairs[SR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                        SR_pairs[SR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)


                    # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                    elif x_off == 0 and y_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y
                        SR_pairs[SR_index, 1, 2] = z

                        SR_pairs[SR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                        SR_pairs[SR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                    # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                    elif x_off == 0 and y_off == 0 and z_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y
                        SR_pairs[SR_index, 1, 2] = z

                        SR_pairs[SR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                        SR_pairs[SR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                    else:
                        SR_pairs[SR_index, 0, 0] = x
                        SR_pairs[SR_index, 0, 1] = y
                        SR_pairs[SR_index, 0, 2] = z

                        SR_pairs[SR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                        SR_pairs[SR_index, 1, 2] = pbc_hardwall(z + z_off, ZDIM)

                    SR_index = SR_index+1


        # delete the self-pair
        SR_pairs = np.delete(SR_pairs, 13,0)
        
        #for i in SR_pairs:
        #    print "[%i, %i, %i] -- [%i, %i, %i]" %(i[0][0],i[0][1],i[0][2],i[1][0],i[1][1],i[1][2],)
        

        return (SR_pairs, np.array([], dtype=np.int), np.array([], dtype=int))
        """

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
        The long-range position to use.  0 for short-range only, 1 for long-range
        only, 2 for both.

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
        
        # the code below is (best I can tell) just reprodicing the function 
        # extract_SR_pairs_from_position_2D_hardwalll, so I have used this
        """
        for x_off in xrange(-1,2):
            for y_off in xrange(-1,2):

                    # if x_off < 0 then the non-central position must come first in the pair
                    if x_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y

                        SR_pairs[SR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)


                    # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                    elif x_off == 0 and y_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y

                        SR_pairs[SR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)

                    else:
                        SR_pairs[SR_index, 0, 0] = x
                        SR_pairs[SR_index, 0, 1] = y

                        SR_pairs[SR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                        SR_pairs[SR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)

                    SR_index = SR_index+1


        # delete the self-pair
        SR_pairs = np.delete(SR_pairs, 4,0)

        # delete pairs which would extend across a boundary
        SR_pairs = delete_pbc_pairs(SR_pairs, 2)
        
        # return is SR, LR, SLR
        return (SR_pairs, np.array([], dtype=int), np.array([], dtype=int))    
        """

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

    # ~ash 2024-01-22
    # declare some variables; note we declare these here as NUMPY_INT_TYP but to be honest that's
    # probably not necessary, this was done while I was debugging some code and it works and is
    # fast but other functions define these tmp variables as int which I expect is fine. The main
    # reason I tried converting from int to NUMPY_INT_TYPE was to see if it would speed things up
    # because then addition operations don't need to do a potential type conversion; however,
    # I don't think it made much difference in the end.
    cdef NUMPY_INT_TYPE SR_index, x_off, y_off, z_off
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
    positions are sorted with the largest chain location first, where 'largest'
    is defined as comparing x and x and then y and y. 
    
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
        


@cython.boundscheck(False)  # Deactivate bounds checking for performance
@cython.wraparound(False)   # Deactivate negative indexing
def delete_pbc_pairs(cnp.ndarray[NUMPY_INT_TYPE, ndim=3] pairs_list, ndims):
    """
    Postprocessing function that takes a list of residue pairs in format:
    [id, pair_id, dimension] = position

    where 
    id indexes into the list
    pair_id is 1 or 0 (each pair contains two separate poistions)
    dimension is the relevant dimenion (0=X,1=Y,2=Z)

    This function systematically goes through the list in this format and deletes pairs
    where the position was -1 (as this means the pair was crossing a periodic boundary).

    Effectively, this trims a list of positions removing any pairs that straddle the PBC.

    Parameters
    ----------
    pairs_list : cnp.ndarray[NUMPY_INT_TYPE, ndim=3]
        The list of pairs to be trimmed

    ndims : int
        The number of dimensions in the system (2 or 3)

    Returns
    -------
    cnp.ndarray[NUMPY_INT_TYPE, ndim=3]
        The trimmed list of pairs

    """

    cdef NUMPY_INT_TYPE minus_one = -1
    cdef int n_pairs = pairs_list.shape[0]
    cdef int idx = 0

    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] pairs_to_delete = np.empty(n_pairs, dtype=NUMPY_INT_TYPE_PYTHON)
    cdef int delete_count = 0


    # if we are in two dimensions
    if ndims == 2:

        # dynamically resize pairs_list. We are looping over a numpy array and changing its size. This means that rather than
        # simply looping over each position we need to loop until the idx (which is only incremented if no deletion occurs) 
        # matches the array length.
        for idx in range(n_pairs):
            
            # if we find pbc crossing delete (so list gets shorter)
            if (pairs_list[idx,0,0] == minus_one) or \
               (pairs_list[idx,1,0] == minus_one) or \
               (pairs_list[idx,0,1] == minus_one) or \
               (pairs_list[idx,1,1] == minus_one):
                #pairs_list = np.delete(pairs_list, idx, 0)
                pairs_to_delete[delete_count] = idx
                delete_count = delete_count + 1

    # same for three dimensions
    else:

        for idx in range(n_pairs):

            # if we find pbc crossing delete (so list gets shorter)
            if (pairs_list[idx,0,0] == minus_one) or \
               (pairs_list[idx,1,0] == minus_one) or \
               (pairs_list[idx,0,1] == minus_one) or \
               (pairs_list[idx,1,1] == minus_one) or \
               (pairs_list[idx,0,2] == minus_one) or \
               (pairs_list[idx,1,2] == minus_one):

                pairs_to_delete[delete_count] = idx
                delete_count = delete_count + 1

    # resize the pairs_to_delete array to the correct size
    pairs_to_delete = pairs_to_delete[:delete_count]

    # finally run delete on those pairs. NOTE this returns a
    # new array so new memory is allocated. This is not ideal
    # but I can't see a way around it.
    pairs_list = np.delete(pairs_list, pairs_to_delete, axis=0)
            
    return pairs_list



##
## OLD IMPLEMENTATION KEPT FOR NOW IN CASE WE NEED TO REVERT.
## THIS IS MUCH MUCH SLOWER THAN THE NEW IMPLEMENTATION ABOVE.
##
def delete_pbc_pairs_OLD(pairs_list, ndims):
    """
    Postprocessing function that takes a list of residue pairs in format:
    [id, pair_id, dimension] = position

    where 
    id indexes into the list
    pair_id is 1 or 0 (each pair contains two separate poistions)
    dimension is the relevant dimenion (0=X,1=Y,2=Z)

    This function systematically goes through the list in this format and deletes pairs
    where the position was -1 (as this means the pair was crossing a periodic boundary).

    Effectively, this trims a list of positions removing any pairs that straddle the PBC.

    """
    n_pairs = len(pairs_list)
    idx = 0

    pairs_to_delete = []
    if ndims == 2:

        # dynamically resize pairs_list. We are looping over a numpy array and changing its size. This means that rather than
        # simply looping over each position we need to loop until the idx (which is only incremented if no deletion occurs) 
        # matches the array length. 
        while True: 
            if idx == len(pairs_list)-1:
                return pairs_list

            # if we find pbc crossing delete (so list gets shorter)
            if (pairs_list[idx,0,0] == -1) or (pairs_list[idx,1,0] == -1) or (pairs_list[idx,0,1] == -1) or (pairs_list[idx,1,1] == -1):
                pairs_list = np.delete(pairs_list, idx, 0)                

            # else increment gets bigger
            else:
                idx=idx+1

    # same for three dimensions
    else:

        # dynamically resize pairs_list (see 2D implementation for brief discussion of this)
        while True: 
            if idx == len(pairs_list)-1:
                return pairs_list

            # if we find pbc crossing delete (so list gets shorter)
            if (pairs_list[idx,0,0] == -1) or (pairs_list[idx,1,0] == -1) or (pairs_list[idx,0,1] == -1) or (pairs_list[idx,1,1] == -1) or (pairs_list[idx,0,2] == -1) or (pairs_list[idx,1,2] == -1):
                pairs_list = np.delete(pairs_list, idx, 0)
                
            # else increment gets bigger
            else:
                idx=idx+1

    return pairs_list


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
