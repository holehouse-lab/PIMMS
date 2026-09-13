## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import numpy as np

from . import inner_loops
from . import inner_loops_hardwall
from . CONFIG import NP_INT_TYPE


#-----------------------------------------------------------------
#
def build_LR_envelope_pairs(positions, LR_binary_array, type_grid, dimensions, hardwall=False):
    """
    Function which builds the non-redundant set of paired interactions between the positions defined in the input list of positions ($positions).
    Long-range interactions are defined as those which extend over TWO lattice sites. A position is defined as lists/tuples of length 2 or 3  
    (x,y or x,y,z coordinates), and the input variable here ($positions) is a LIST of positions.

    Example: if positions was a list with a SINGLE 2D position in it - ``[[4, 4]]``
    - then we'd return 16 pairs corresponding to a pair between ``[4, 4]`` and one
    of A, B, C, D, E, F, J, K, O, P, T, U, V, W, Z, Y as shown on the diagram
    below::

                      x---->

                2   3   4   5   6
       y      +-------------------+
       |    2 | A | B | C | D | E |
       |    3 | F | G | H | I | J |
       v    4 | K | L | M | N | O |
            5 | P | Q | R | S | T |
            6 | U | V | W | X | Y |
              +-------------------+

    The return variable is a list of pairs of the format::

        [[A, B], [A, C], [B, A]]

    where A/B/C/D are tuples of positions (note the A/B/C/D here do not
    correspond to the letters in the diagram above - I'd just run out of
    letters).

    ALSO note the pair-ordering convention (inherited from the inner_loops
    extractors): each pair is ordered by the SIGN of the pre-PBC offset from the
    scanned central position - if the first non-zero component of the offset is
    positive the neighbour comes first, otherwise the central position comes
    first. The rule is antisymmetric (the same physical pair seen from either
    end is ordered identically), which is exactly what makes the downstream
    de-duplication correct. It is NOT a sort by numeric coordinate values:
    across a periodic face the two orderings disagree (e.g. a pair spanning the
    x wrap can legitimately be returned as [[6,3],[0,3]]).

    Obviously when dealing with pairs of positions the order of the positions doesn't matter
    but the fact that we consistently order the pairs in the same way means that if two IDENTICAL
    pairs are found they will appear identical to one another and can easily be removed easily.

    Internally this dispatches to the 2D or 3D Cython
    ``extract_LR_pairs_from_position`` routine for each position, collects the
    long-range (LR) and super-long-range (SLR) candidate pairs, removes
    duplicates by reshaping/viewing as a void dtype and applying
    ``np.unique``, and finally reshapes the result into arrays of pairs.

    Parameters
    ----------
    positions : list
        The list of bead positions (each a 2- or 3-element sequence of integer
        lattice coordinates) over which long-range envelope pairs are
        constructed.

    LR_binary_array : numpy.ndarray
        1D integer array with one entry per entry in ``positions``, set to 1
        where the bead engages in long-range interactions and 0 otherwise. This
        is the array returned by ``Chain.get_LR_binary_array()``, and beads
        flagged 0 contribute no pairs.

    type_grid : numpy.ndarray
        The lattice type grid (2D or 3D integer array, matching ``dimensions``)
        used to look up occupancy/identity at candidate neighbour sites.

    dimensions : list
        The box dimensions as a list of ints; its length (2 or 3) selects the
        2D or 3D code path and sets the extent of the periodic wrapping.

    hardwall : bool, optional
        If True the hardwall extractors are used, so no pair across a box wall
        is ever emitted (neighbour sites outside the box do not exist). Default
        False (periodic: neighbours wrap).

    Returns
    -------
    tuple of numpy.ndarray
        A 2-tuple ``(LR_pairs, SLR_pairs)`` where each element is a
        duplicate-free numpy array of shape ``(n_pairs, 2, ndim)`` (with
        ``ndim`` equal to 2 or 3); positions that generate no pairs - or an
        empty ``positions`` input - give a ``(0, 2, ndim)`` array, so callers
        can always unpack and concatenate without special-casing.

    """

    if len(positions) == 0:
        # every caller unpacks a (LR, SLR) pair, so return one - with the same
        # (0, 2, ndim) shape build_all_envelope_pairs uses for its empties
        empty = np.empty((0, 2, len(dimensions)), dtype=NP_INT_TYPE)
        return (empty, empty.copy())
    
    
    LR_list = []         
    SLR_list = []

    # Under a hardwall the hardwall extractors are used, which never emit a pair
    # across a wall; the periodic extractors were used unconditionally before,
    # leaving it to the energy kernel to drop straddling pairs (the energies
    # were right, but the hardwall extractors sat unused and every hardwall
    # envelope carried pairs that could never contribute).
    if hardwall:
        extract_2D = inner_loops_hardwall.extract_LR_pairs_from_position_2D_hardwall
        extract_3D = inner_loops_hardwall.extract_LR_pairs_from_position_3D_hardwall
    else:
        extract_2D = inner_loops.extract_LR_pairs_from_position_2D
        extract_3D = inner_loops.extract_LR_pairs_from_position_3D

    # define differences for 2D vs 3D
    dims = len(dimensions)
    if dims == 2:
        reshape_axis = 4
    else:
        reshape_axis = 6
        
    


    # >>>>>>>>>>>>>>>> if 2D
    if len(dimensions) == 2:

        for i in range(0, len(positions)):
            (LR_tmp, SLR_tmp)  = extract_2D(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_grid, dimensions[0], dimensions[1])
            
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
            return (np.reshape(LR_duplicate_free, (len(LR_duplicate_free), 2,2)), np.empty((0, 2, 2), dtype=NP_INT_TYPE))

        elif return_SLR > 0:
            return (np.empty((0, 2, 2), dtype=NP_INT_TYPE), np.reshape(SLR_duplicate_free, (len(SLR_duplicate_free), 2,2)))

        else:
            return (np.empty((0, 2, dims), dtype=NP_INT_TYPE), np.empty((0, 2, dims), dtype=NP_INT_TYPE))


    # >>>>>>>>>>>>>>> if 3D
    else:

        for i in range(0, len(positions)):
            (LR_tmp, SLR_tmp)  = extract_3D(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_grid, dimensions[0], dimensions[1], dimensions[2])

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
            return (np.reshape(LR_duplicate_free, (len(LR_duplicate_free), 2,3)), np.empty((0, 2, 3), dtype=NP_INT_TYPE))

        elif return_SLR > 0:
            return (np.empty((0, 2, 3), dtype=NP_INT_TYPE), np.reshape(SLR_duplicate_free, (len(SLR_duplicate_free), 2,3)))

        else:
            return (np.empty((0, 2, dims), dtype=NP_INT_TYPE), np.empty((0, 2, dims), dtype=NP_INT_TYPE))
            
            
