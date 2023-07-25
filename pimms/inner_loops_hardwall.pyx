## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2020
## 
## ...........................................................................

import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython 

from pimms.latticeExceptions import InnerLoopException

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b


ctypedef cnp.int_t NUMPY_INT_TYPE

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
                                                      int LR_position, 
                                                      NUMPY_INT_TYPE[:,:,:] type_grid,
                                                      int XDIM, 
                                                      int YDIM, 
                                                      int ZDIM):
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
    cdef int SLR_index, SR_index, LR_index, x_off, y_off, z_off;
    cdef int x_tmp, y_tmp, z_tmp;
    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    cdef int z = position[2]
    
    # short range only
    if LR_position == 0:

        return (extract_SR_pairs_from_position_3D_hardwall(position, XDIM, YDIM, ZDIM), np.array([], dtype=int), np.array([], dtype=int))


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
        
    
        SR_pairs  = np.zeros((27,2,3), dtype=int)
        LR_pairs  = np.zeros((98, 2, 3), dtype=int)
        SLR_pairs = np.zeros((218, 2, 3), dtype=int)
                    
        SR_index  = 0
        LR_index  = 0
        SLR_index = 0
        
        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                for z_off in xrange(-3,4):
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2 and abs(z_off) <2:
                        
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

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3 and abs(z_off) < 3:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)
                        z_tmp = pbc_hardwall(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z

                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            LR_pairs[LR_index, 0, 2] = z
                            
                            LR_pairs[LR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 1, 2] = pbc_hardwall(z + z_off, ZDIM)
                            
                        LR_index = LR_index+1

                    # SUPER LONG RANGE INTERACTIONS...
                    else:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)
                        z_tmp = pbc_hardwall(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z

                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            SLR_pairs[SLR_index, 0, 2] = z
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 1, 2] = pbc_hardwall(z + z_off, ZDIM)
                            
                        SLR_index = SLR_index+1


        # delete the self-pair for the short range interaction (no such pair for the
        # long-range interactions)
        SR_pairs = np.delete(SR_pairs, 13,0)
        SR_pairs = delete_pbc_pairs(SR_pairs, 3)

        # do same for LR and SLR (deletion in return call)
        good_LR_pairs  = LR_pairs[0:LR_index]
        good_SLR_pairs = SLR_pairs[0:SLR_index]
        
        return (SR_pairs, delete_pbc_pairs(good_LR_pairs,3), delete_pbc_pairs(good_SLR_pairs,3))
        
                        


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


    """

    # declare some variables
    cdef int SR_index, LR_index, SLR_index, x_off, y_off        

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
                
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs

    
    # short range only
    if LR_position == 0:

        return (extract_SR_pairs_from_position_2D_hardwall(position, XDIM, YDIM), np.array([], dtype=int), np.array([], dtype=int))    
        
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


        SR_pairs  = np.zeros((9,  2, 2), dtype=int)
        LR_pairs  = np.zeros((16, 2, 2), dtype=int) # 5 x 5 -  3 x 3       
        SLR_pairs = np.zeros((24, 2, 2), dtype=int) # 7 x 7 -  5 x 5
        
        SR_index = 0    
        LR_index = 0
        SLR_index = 0

        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2:
                        
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

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)
                        
                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y                            
                            
                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)                            

                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y

                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            
                            LR_pairs[LR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # SUPER LONG RANGE INTERACTIONS...
                    else:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)

                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y

                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                
                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        SLR_index = SLR_index+1

        # delete the self-pair for the short range interaction (no such pair for the
        # long-range interactions) and pairs that cross the PBC
        SR_pairs = np.delete(SR_pairs, 4,0)
        SR_pairs = delete_pbc_pairs(SR_pairs, 2)

        # do same for LR and SLR (deletion in return call)
        good_LR_pairs  = LR_pairs[0:LR_index]
        good_SLR_pairs = SLR_pairs[0:SLR_index]
      
        return (SR_pairs, delete_pbc_pairs(good_LR_pairs,2), delete_pbc_pairs(good_SLR_pairs,2))

                                
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
    cdef int LR_index, SLR_index, x_off, y_off, z_off;
    cdef int x_tmp, y_tmp, z_tmp;
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    cdef int z = position[2]
    
    # if no long-range interactions required the return empy 
    # arrays
    if LR_position == 0:
        return (np.array([], dtype=int),np.array([], dtype=int))

    elif LR_position == 1:
        LR_pairs = np.zeros((98, 2, 3), dtype=int)                
        SLR_pairs = np.zeros((218, 2, 3), dtype=int)                
        
        LR_index = 0
        SLR_index = 0
        
        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                for z_off in xrange(-3,4):
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2 and abs(z_off) <2:
                        continue

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3 and abs(z_off) < 3:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)
                        z_tmp = pbc_hardwall(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z

                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            LR_pairs[LR_index, 0, 2] = z
                            
                            LR_pairs[LR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            LR_pairs[LR_index, 1, 2] = pbc_hardwall(z + z_off, ZDIM)
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Super long range interaction
                    else:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)
                        z_tmp = pbc_hardwall(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:                            
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z

                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_hardwall(z + z_off, ZDIM)

                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            SLR_pairs[SLR_index, 0, 2] = z
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 1, 2] = pbc_hardwall(z + z_off, ZDIM)
                            
                        SLR_index = SLR_index+1

        good_LR_pairs  = LR_pairs[0:LR_index]
        good_SLR_pairs = SLR_pairs[0:SLR_index]
        return (delete_pbc_pairs(good_LR_pairs,3), delete_pbc_pairs(good_SLR_pairs, 3))

    else:
        raise InnerLoopException('Invalid LR option passed')


##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def extract_SR_pairs_from_position_3D_hardwall(NUMPY_INT_TYPE[:] position,                 
                                               int XDIM, 
                                               int YDIM, 
                                               int ZDIM):
    """


    """
    
    # declare some variables
    cdef int SR_index, x_off, y_off, z_off;
    cdef int x_tmp, y_tmp, z_tmp;
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((27,2,3), dtype=int)

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    cdef int z = position[2]
    
    SR_index = 0
    for x_off in xrange(-1,2):
        for y_off in xrange(-1,2):
            for z_off in xrange(-1,2):

                # if x_off < 0 then the non-central position must come first in the pair
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
    
    # delete pairs which would extend across a boundary
    SR_pairs = delete_pbc_pairs(SR_pairs, 3)

    return (SR_pairs)


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
    cdef int LR_index, SLR_index, x_off, y_off;
    cdef int x_tmp, y_tmp;
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    
    # if no long-range interactions required the return empy 
    # array
    if LR_position == 0:
        return (np.array([]),np.array([]))

    elif LR_position == 1:
        LR_pairs = np.zeros((16, 2, 2), dtype=int)        
        SLR_pairs = np.zeros((24, 2, 2), dtype=int)        

        LR_index = 0
        SLR_index = 0
        
        # loop over long range cube around site 
        for x_off in xrange(-3,4):
            for y_off in xrange(-3,4):
                                        
                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # if short range_interaction
                    if abs(x_off) < 2 and abs(y_off) < 2:
                        continue

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3:

                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)

                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y                            
                            
                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y

                            LR_pairs[LR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)


                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            
                            LR_pairs[LR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Super long range interaction
                    else:
                        # is it worth continuing?
                        x_tmp = pbc_hardwall(x + x_off, XDIM)
                        y_tmp = pbc_hardwall(y + y_off, YDIM)

                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y                            
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y

                            SLR_pairs[SLR_index, 0, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_hardwall(y + y_off, YDIM)


                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_hardwall(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_hardwall(y + y_off, YDIM)
                            
                        SLR_index = SLR_index+1
                        

        good_LR_pairs  = LR_pairs[0:LR_index]
        good_SLR_pairs = SLR_pairs[0:SLR_index]
        return (delete_pbc_pairs(good_LR_pairs,2), delete_pbc_pairs(good_SLR_pairs,2))

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
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((9,2,2), dtype=int)

    # first set the central x, y and z positions
    cdef int x = position[0]
    cdef int y = position[1]
    
    SR_index = 0
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

    return (SR_pairs)
        
    


def delete_pbc_pairs(pairs_list, ndims):
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
    idx=0
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
cdef int pbc_hardwall(int value, int DIM):    
    """
    Takes an offset position (value) and the max dimensions in that axis (DIM) and if
    this new position is outside the lattice returns -1 (i.e. this is one part of a pair
    that will straddle the periodic boundary)
    """
                            
    if (value < 0) or value > (DIM-1):
        return -1
    else:
        return value
