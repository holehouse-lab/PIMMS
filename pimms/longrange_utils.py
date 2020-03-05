## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

import random
import numpy as np

from . import inner_loops
from . import lattice_utils

#-----------------------------------------------------------------
#    
def get_LR_positions(all_positions, relevant_indices, LR_IDX):
    """
    Function which takes a full list of positions, a relevant indices (i.e. the indices of interest), a 
    list of the indicies which experience long-ranger interactions and the dimensions, and returns a LIST
    of positions which correspond to residues that undergo long range interactions.
        
    """
        
    relevant_positions =[]
    for idx in relevant_indices:
        if idx in LR_IDX:
            relevant_positions.append(all_positions[idx])
                
    return relevant_positions


#-----------------------------------------------------------------
#
def build_LR_envelope_pairs(positions, LR_binary_array, type_grid, dimensions):
    """
    Function which builds the non-redundant set of paired interactions between the positions defined in the input list of positions ($positions).
    Long-range interactions are defined as those which extend over TWO lattice sites. A position is defined as lists/tuples of length 2 or 3  
    (x,y or x,y,z coordinates), and the input variable here ($positions) is a LIST of positions.

    Example: if positions was a list with a SINGLE 2D position in it - [ [4,4] ]  - then we'd return 16 pairs corresponding to a pair between 
    [4,4] and one of A,B,C,D,E,F,J,K,O,P,T,U,V,W,Z,Y as shown on the diagram below

                  x---->
          
            2   3   4   5   6
   y      +-------------------+
   |    2 | A | B | C | D | E |
   |    3 | F | G | H | I | J |
   v    4 | K | L | M | N | O | 
        5 | P | Q | R | S | T |
        6 | U | V | W | X | Y | 
          +-------------------+

    The return variable is a list of pairs of the format:
    
    [[A,B],[A,C],[B,A]] 

    where A/B/C/D are tuples of positions (note the A/B/C/D here do not correspond to the 
    letters in the digram above - I'd just run out of letters...

    ALSO note that the positions in the return list are sorted such that in each pair if we
    compare the absolute value of the x/y/z components in order the first position is always
    larger than the second.

    e.g.

    [ [3,1], [1,7] ] # because 3 > 1

    [ [2,2,5], [2,2,3]] # because 2==2 but 5 > 3

    [ [19, 1], [0, 1]] # because 19 > 0

    Obviously when dealing with pairs of positions the order of the positions doesn't matter
    but the fact that we consistently order the pairs in the same way means that if two IDENTICAL
    pairs are found they will appear identical to one another and can easily be removed easily.

    """

    if len(positions) == 0:
        return []
    
    
    LR_list = []         
    SLR_list = []

    # define differences for 2D vs 3D
    dims = len(dimensions)
    if dims == 2:
        reshape_axis = 4
    else:
        reshape_axis = 6
        
    


    # >>>>>>>>>>>>>>>> if 2D
    if len(dimensions) == 2:

        for i in range(0, len(positions)):
            (LR_tmp, SLR_tmp)  = inner_loops.extract_LR_pairs_from_position_2D(np.array(positions[i], dtype=int), LR_binary_array[i], type_grid, dimensions[0], dimensions[1])
            
            if len(LR_tmp) > 0:
                LR_list.append(LR_tmp)

            if len(SLR_tmp) > 0:
                SLR_list.append(SLR_tmp)



        """
        original cide incase something is wrong
        if len (LR_list) > 0:
            long_range_pairs = np.concatenate(LR_list)
        else:
            return np.array([])
                
        num_pairs = len(long_range_pairs)

        reshaped = np.reshape(long_range_pairs, (num_pairs, 4))

        b = np.ascontiguousarray(reshaped).view(np.dtype((np.void, reshaped.dtype.itemsize * reshaped.shape[1])))
        _, idx = np.unique(b, return_index=True)

        duplicate_free = reshaped[idx]

        return np.reshape(duplicate_free, (len(duplicate_free), 2,2))
        """

        ## This section figures out which sets of pairs we're going to return
        if len(LR_list) > 0:
            long_range_pairs = np.concatenate(LR_list)
            return_LR = len(long_range_pairs)
        else:
            return_LR = 0


        if len(SLR_list) > 0:
            super_long_range_pairs = np.concatenate(SLR_list)            
            return_SLR = len(super_long_range_pairs)
        else:
            return_SLR = 0
        

        # for those with pairs we remove duplicates and restructure
                
        if return_LR > 0:
            num_LR_pairs = len(long_range_pairs)        
            reshaped_LR = np.reshape(long_range_pairs, (num_LR_pairs, 4))
            b = np.ascontiguousarray(reshaped_LR).view(np.dtype((np.void, reshaped_LR.dtype.itemsize * reshaped_LR.shape[1])))
            _, idx = np.unique(b, return_index=True)
            LR_duplicate_free = reshaped_LR[idx]

        if return_SLR > 0:
            num_SLR_pairs = len(super_long_range_pairs)        
            reshaped_SLR = np.reshape(super_long_range_pairs, (num_SLR_pairs, 4))
            b = np.ascontiguousarray(reshaped_SLR).view(np.dtype((np.void, reshaped_SLR.dtype.itemsize * reshaped_SLR.shape[1])))
            _, idx = np.unique(b, return_index=True)
            SLR_duplicate_free = reshaped_SLR[idx]


        if return_LR > 0 and return_SLR > 0:            
            return (np.reshape(LR_duplicate_free, (len(LR_duplicate_free), 2,2)), np.reshape(SLR_duplicate_free, (len(SLR_duplicate_free), 2,2)))

        elif return_LR > 0:
            return (np.reshape(LR_duplicate_free, (len(LR_duplicate_free), 2,2)), np.array([]))

        elif return_SLR > 0:
            return (np.array([]), np.reshape(SLR_duplicate_free, (len(SLR_duplicate_free), 2,2)))

        else:
            return (np.array([]), np.array([]))


    # >>>>>>>>>>>>>>> if 3D
    else:

        for i in range(0, len(positions)):
            (LR_tmp, SLR_tmp)  = inner_loops.extract_LR_pairs_from_position_3D(np.array(positions[i], dtype=int), LR_binary_array[i], type_grid, dimensions[0], dimensions[1], dimensions[2])

            if len(LR_tmp) > 0:
                LR_list.append(LR_tmp)

            if len(SLR_tmp) > 0:
                SLR_list.append(SLR_tmp)

        ## This section figures out which sets of pairs we're going to return
        if len(LR_list) > 0:
            long_range_pairs = np.concatenate(LR_list)
            return_LR = len(long_range_pairs)
        else:
            return_LR = 0


        if len(SLR_list) > 0:
            super_long_range_pairs = np.concatenate(SLR_list)            
            return_SLR = len(super_long_range_pairs)
        else:
            return_SLR = 0
        

        # for those with pairs we remove duplicates and restructure
                
        if return_LR > 0:
            num_LR_pairs = len(long_range_pairs)        
            reshaped_LR = np.reshape(long_range_pairs, (num_LR_pairs, 6))
            b = np.ascontiguousarray(reshaped_LR).view(np.dtype((np.void, reshaped_LR.dtype.itemsize * reshaped_LR.shape[1])))
            _, idx = np.unique(b, return_index=True)
            LR_duplicate_free = reshaped_LR[idx]

        if return_SLR > 0:
            num_SLR_pairs = len(super_long_range_pairs)        
            reshaped_SLR = np.reshape(super_long_range_pairs, (num_SLR_pairs, 6))
            b = np.ascontiguousarray(reshaped_SLR).view(np.dtype((np.void, reshaped_SLR.dtype.itemsize * reshaped_SLR.shape[1])))
            _, idx = np.unique(b, return_index=True)
            SLR_duplicate_free = reshaped_SLR[idx]


        if return_LR > 0 and return_SLR > 0:            
            return (np.reshape(LR_duplicate_free, (len(LR_duplicate_free), 2,3)), np.reshape(SLR_duplicate_free, (len(SLR_duplicate_free), 2,3)))

        elif return_LR > 0:
            return (np.reshape(LR_duplicate_free, (len(LR_duplicate_free), 2,3)), np.array([]))

        elif return_SLR > 0:
            return (np.array([]), np.reshape(SLR_duplicate_free, (len(SLR_duplicate_free), 2,3)))

        else:
            return (np.array([]), np.array([]))
            
            
