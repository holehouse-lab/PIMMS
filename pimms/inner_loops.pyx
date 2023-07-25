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

from cython.view cimport array


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

    The pairs are inherently numbered (i.e. [A-B] would be A then B). To determine which
    of the two positions is first in the pair we use the numerical value of the x/y/z 
    positions

    """
    
    # declare some variables
    cdef int SLR_index, SR_index, LR_index, x_off, y_off, z_off;
    cdef int x_tmp, y_tmp, z_tmp;
    #cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] LR_pairs 
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SLR_pairs 


    #SR_pairs = np.zeros((27,2,3), dtype=int)

    #SR_pairs = np.empty((27,2,3), dtype=int)

    cdef int dim1 = 27
    cdef int dim2 = 2
    cdef int dim3 = 3
    cdef int[:,:,:] SR_pairs = array(shape=(dim1, dim2, dim3), itemsize=sizeof(int), format="i")

    
    cdef int i_1, j_1, k_1
    for i_1 in range(27):
        for j_1 in range(2):
            for k_1 in range(3):
                SR_pairs[i_1, j_1, k_1] = 0
        
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

                    # if x_off > 0 then the non-central position must come first in the pair
                    if x_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y
                        SR_pairs[SR_index, 1, 2] = z

                        SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                        SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)


                    # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                    elif x_off == 0 and y_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y
                        SR_pairs[SR_index, 1, 2] = z

                        SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                        SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                    # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                    elif x_off == 0 and y_off == 0 and z_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y
                        SR_pairs[SR_index, 1, 2] = z

                        SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                        SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                    else:
                        SR_pairs[SR_index, 0, 0] = x
                        SR_pairs[SR_index, 0, 1] = y
                        SR_pairs[SR_index, 0, 2] = z

                        SR_pairs[SR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                        SR_pairs[SR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)

                    SR_index = SR_index+1


        # delete the self-pair
        SR_pairs = np.delete(SR_pairs, 13,0)

        """
        for i in SR_pairs:
            print "[%i, %i, %i] -- [%i, %i, %i]" %(i[0][0],i[0][1],i[0][2],i[1][0],i[1][1],i[1][2],)
        """

        return (SR_pairs, np.array([], dtype=int), np.array([], dtype=int))

    elif LR_position == 1:
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

                            SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y
                            SR_pairs[SR_index, 1, 2] = z

                            SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y
                            SR_pairs[SR_index, 1, 2] = z
                            
                            SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        else:
                            SR_pairs[SR_index, 0, 0] = x
                            SR_pairs[SR_index, 0, 1] = y
                            SR_pairs[SR_index, 0, 2] = z
                            
                            SR_pairs[SR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            SR_pairs[SR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)
                            
                        SR_index = SR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3 and abs(z_off) < 3:

                        # is it worth continuing?
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)
                        z_tmp = pbc_correction(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z

                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            LR_pairs[LR_index, 0, 2] = z
                            
                            LR_pairs[LR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)
                            
                        LR_index = LR_index+1

                    # SUPER LONG RANGE INTERACTIONS...
                    else:

                        # is it worth continuing?
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)
                        z_tmp = pbc_correction(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z

                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            SLR_pairs[SLR_index, 0, 2] = z
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)
                            
                        SLR_index = SLR_index+1


        # delete the self-pair for the short range interaction (no such pair for the
        # long-range interactions)
        SR_pairs = np.delete(SR_pairs, 13,0)


        return (SR_pairs, LR_pairs[0:LR_index], SLR_pairs[0:SLR_index])
        
                        


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


    """
    
    # declare some variables
    cdef int SR_index, LR_index, SLR_index, x_off, y_off
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] SR_pairs = np.zeros((9,2,2), dtype=int)
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

                    # if x_off < 0 then the non-central position must come first in the pair
                    if x_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y

                        SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)


                    # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                    elif x_off == 0 and y_off > 0:
                        SR_pairs[SR_index, 1, 0] = x
                        SR_pairs[SR_index, 1, 1] = y

                        SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)

                    else:
                        SR_pairs[SR_index, 0, 0] = x
                        SR_pairs[SR_index, 0, 1] = y

                        SR_pairs[SR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                        SR_pairs[SR_index, 1, 1] = pbc_correction(y + y_off, YDIM)

                    SR_index = SR_index+1


        # delete the self-pair
        SR_pairs = np.delete(SR_pairs, 4,0)


        # return is SR, LR, SLR
        return (SR_pairs, np.array([], dtype=int), np.array([], dtype=int))    

    elif LR_position == 1:

        LR_pairs = np.zeros((16, 2, 2), dtype=int)  # 5 x 5 -  3 x 3       
        SLR_pairs = np.zeros((24, 2, 2), dtype=int) # 7 x 7 -  5 x 5
        
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

                            SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SR_pairs[SR_index, 1, 0] = x
                            SR_pairs[SR_index, 1, 1] = y

                            SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)


                        else:
                            SR_pairs[SR_index, 0, 0] = x
                            SR_pairs[SR_index, 0, 1] = y
                            
                            SR_pairs[SR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            SR_pairs[SR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            
                        SR_index = SR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Long range interaction
                    elif abs(x_off) < 3 and abs(y_off) < 3:

                        # is it worth continuing?
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)
                        
                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y                            
                            
                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)                            

                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y

                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            
                            LR_pairs[LR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # SUPER LONG RANGE INTERACTIONS...
                    else:

                        # is it worth continuing?
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)

                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y

                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                
                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_correction(y + y_off, YDIM)

                            
                        SLR_index = SLR_index+1


        # delete the self-pair for the short range interaction (no such pair for the
        # long-range interactions)
        SR_pairs = np.delete(SR_pairs, 4,0)

        return (SR_pairs, LR_pairs[0:LR_index], SLR_pairs[0:SLR_index])
                                
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
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)
                        z_tmp = pbc_correction(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z

                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y
                            LR_pairs[LR_index, 1, 2] = z
                            
                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            LR_pairs[LR_index, 0, 2] = z
                            
                            LR_pairs[LR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            LR_pairs[LR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Super long range interaction
                    else:

                        # is it worth continuing?
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)
                        z_tmp = pbc_correction(z + z_off, ZDIM)

                        if type_grid[x_tmp, y_tmp, z_tmp] == 0:                            
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)                                                                                    
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z

                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off == 0 and z_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y
                            SLR_pairs[SLR_index, 1, 2] = z
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            SLR_pairs[SLR_index, 0, 2] = z
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            SLR_pairs[SLR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)
                            
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

                    SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                    SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                    SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                    
                # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                elif x_off == 0 and y_off > 0:
                    SR_pairs[SR_index, 1, 0] = x
                    SR_pairs[SR_index, 1, 1] = y
                    SR_pairs[SR_index, 1, 2] = z

                    SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                    SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                    SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                # if x_off == 0  and y_off is == 0 and z_off < 1 then the non-central position must come first in the pair
                elif x_off == 0 and y_off == 0 and z_off > 0:
                    SR_pairs[SR_index, 1, 0] = x
                    SR_pairs[SR_index, 1, 1] = y
                    SR_pairs[SR_index, 1, 2] = z

                    SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                    SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                    SR_pairs[SR_index, 0, 2] = pbc_correction(z + z_off, ZDIM)

                else:
                    SR_pairs[SR_index, 0, 0] = x
                    SR_pairs[SR_index, 0, 1] = y
                    SR_pairs[SR_index, 0, 2] = z

                    SR_pairs[SR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                    SR_pairs[SR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                    SR_pairs[SR_index, 1, 2] = pbc_correction(z + z_off, ZDIM)

                SR_index = SR_index+1


    # delete the self-pair
    SR_pairs = np.delete(SR_pairs, 13,0)
        

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
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)

                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y                            
                            
                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            LR_pairs[LR_index, 1, 0] = x
                            LR_pairs[LR_index, 1, 1] = y

                            LR_pairs[LR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 0, 1] = pbc_correction(y + y_off, YDIM)


                        else:
                            LR_pairs[LR_index, 0, 0] = x
                            LR_pairs[LR_index, 0, 1] = y
                            
                            LR_pairs[LR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            LR_pairs[LR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            
                        LR_index = LR_index+1

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Super long range interaction
                    else:
                        # is it worth continuing?
                        x_tmp = pbc_correction(x + x_off, XDIM)
                        y_tmp = pbc_correction(y + y_off, YDIM)

                        if type_grid[x_tmp, y_tmp] == 0:
                            continue

                        # if x_off < 0 then the non-central position must come first in the pair
                        if x_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y                            
                            
                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)
                            
                        # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
                        elif x_off == 0 and y_off > 0:
                            SLR_pairs[SLR_index, 1, 0] = x
                            SLR_pairs[SLR_index, 1, 1] = y

                            SLR_pairs[SLR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 0, 1] = pbc_correction(y + y_off, YDIM)


                        else:
                            SLR_pairs[SLR_index, 0, 0] = x
                            SLR_pairs[SLR_index, 0, 1] = y
                            
                            SLR_pairs[SLR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                            SLR_pairs[SLR_index, 1, 1] = pbc_correction(y + y_off, YDIM)
                            
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

                SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)


            # if x_off == 0  and y_off is <0 then the non-central position must come first in the pair
            elif x_off == 0 and y_off > 0:
                SR_pairs[SR_index, 1, 0] = x
                SR_pairs[SR_index, 1, 1] = y

                SR_pairs[SR_index, 0, 0] = pbc_correction(x + x_off, XDIM)
                SR_pairs[SR_index, 0, 1] = pbc_correction(y + y_off, YDIM)

            else:
                SR_pairs[SR_index, 0, 0] = x
                SR_pairs[SR_index, 0, 1] = y

                SR_pairs[SR_index, 1, 0] = pbc_correction(x + x_off, XDIM)
                SR_pairs[SR_index, 1, 1] = pbc_correction(y + y_off, YDIM)

            SR_index = SR_index+1


    # delete the self-pair
    SR_pairs = np.delete(SR_pairs, 4,0)

    return (SR_pairs)
        
    
    

@cython.cdivision(True)
cdef int pbc_correction(int value, int DIM):    
    """
    Performs intelligent periodic boundary correction
    which FIRST checks to see if we have a negative
    value and IF NOT uses the % operator - this means
    we can use % without checking the sign giving
    a 35% speedup per call
    """
                            
    if value < 0:
        return DIM+value
    else:
        return (value % DIM)
