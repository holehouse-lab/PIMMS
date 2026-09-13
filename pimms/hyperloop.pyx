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
import pimms.inner_loops as inner_loops

from pimms.cython_config cimport NUMPY_INT_TYPE
from pimms.CONFIG import NP_INT_TYPE as NUMPY_INT_TYPE_PYTHON


##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def get_adjacent_sites_2D(int position1, int position2, int X_DIM, int Y_DIM, int extent_range):
    """
    Build the square block of sites centred on (position1, position2) out to
    extent_range in each direction.

    Note this INCLUDES the original site and does perform PBC correction


    Optimization notes

    - Switching out the while loops for straightup list of position assignment did
      not bring any speedup

    Parameters
    ----------
    position1 : int (C int)
        x coordinate of the central site.
    position2 : int (C int)
        y coordinate of the central site.
    X_DIM : int (C int)
        Box size along x, used to wrap the returned x coordinates.
    Y_DIM : int (C int)
        Box size along y, used to wrap the returned y coordinates.
    extent_range : int (C int)
        Half-width of the block, i.e. offsets from -extent_range to
        +extent_range are taken along each axis.

    Returns
    -------
    ((2 * extent_range + 1)**2, 2) numpy.ndarray, NUMPY_INT_TYPE (int32)
        Wrapped site coordinates, x varying slowest and y fastest.
    """

    cdef int range_multiplier;

    # number of sitesin positions
    range_multiplier = ((extent_range*2)+1) * ((extent_range*2)+1) 

    # declare vars
    cdef int array_index;    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=2] positions = np.zeros((range_multiplier,2), dtype=NUMPY_INT_TYPE_PYTHON)
    cdef int x
    cdef int y

    cdef int minval
    cdef int maxval
                
    array_index = 0

    minval = 0 - extent_range;
    maxval = 1 + extent_range;



    array_index = 0
    x = minval

    while x < maxval:
        y = minval
        while y < maxval:
            positions[array_index, 0] = (position1+x) % X_DIM
            positions[array_index, 1] = (position2+y) % Y_DIM
                
            array_index = array_index+1
            y=y+1
        x=x+1


    return positions 


##
#################################################################################################
##
@cython.boundscheck(False)
@cython.wraparound(False) 
def get_adjacent_sites_3D(int position1, int position2, int position3, int X_DIM, int Y_DIM, int Z_DIM, int extent_range):
    """
    Build the cube of sites centred on (position1, position2, position3) out to
    extent_range in each direction.

    Note this INCLUDES the original site and does perform PBC correction


    Optimization notes

    - Switching out the while loops for straightup list of position assignment did
      not bring any speedup

    Parameters
    ----------
    position1 : int (C int)
        x coordinate of the central site.
    position2 : int (C int)
        y coordinate of the central site.
    position3 : int (C int)
        z coordinate of the central site.
    X_DIM : int (C int)
        Box size along x, used to wrap the returned x coordinates.
    Y_DIM : int (C int)
        Box size along y, used to wrap the returned y coordinates.
    Z_DIM : int (C int)
        Box size along z, used to wrap the returned z coordinates.
    extent_range : int (C int)
        Half-width of the cube, i.e. offsets from -extent_range to
        +extent_range are taken along each axis.

    Returns
    -------
    ((2 * extent_range + 1)**3, 3) numpy.ndarray, NUMPY_INT_TYPE (int32)
        Wrapped site coordinates, x varying slowest and z fastest.
    """

    cdef int range_multiplier;

    # number of sitesin positions
    range_multiplier = ((extent_range*2)+1) * ((extent_range*2)+1) * ((extent_range*2)+1) 


    # declare vars
    cdef int array_index;    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=2] positions = np.zeros((range_multiplier,3), dtype=NUMPY_INT_TYPE_PYTHON)
    cdef int x
    cdef int y
    cdef int z
    cdef int minval
    cdef int maxval
                
    array_index = 0

    minval = 0 - extent_range;
    maxval = 1 + extent_range;

    x = minval
    while x < maxval:
        y = minval
        while y < maxval:            
            z = minval
            while z < maxval:
                
                positions[array_index, 0] = (position1 + x) % X_DIM
                positions[array_index, 1] = (position2 + y) % Y_DIM
                positions[array_index, 2] = (position3 + z) % Z_DIM
                
                array_index = array_index+1
                
                z=z+1
            y=y+1
        x=x+1
    
    return positions


###
###
###
@cython.boundscheck(False)
@cython.wraparound(False) 
def get_unique_interface_pairs_3D(int E_x, int E_y, int E_z, int X_DIM, int Y_DIM, int Z_DIM, NUMPY_INT_TYPE[:,:,:] lattice):
    """
    Returns the nearest neighbour set of pairs corresponding to the central position (E_x, E_y, E_z) where at least one of the 
    two positions (central position and somewhere else) is an empty site. This is useful if you KNOW the central position is filled a
    priori, as it'll then return a non-redundant list of ONLY the pairs which correspond to the chain-solvent interface
    site.

    This is especially useful for rigid-body cluster moves because we KNOW that there cannot be components of the cluster in nearest
    neighbour contact with another chain which is not part of the cluster (otherwise it would be part of the cluster!) such that by
    getting the unique_interface_pairs for each position in the cluster we automatically have the full, non-redundant set of pairs
    over which we need to evaluate energy changes on a move - no need for any kind of post-processing.

    Parameters
    ----------
    E_x : int (C int)
        x coordinate of the central (assumed occupied) site.
    E_y : int (C int)
        y coordinate of the central site.
    E_z : int (C int)
        z coordinate of the central site.
    X_DIM : int (C int)
        Box size along x, used to wrap the neighbour shell.
    Y_DIM : int (C int)
        Box size along y.
    Z_DIM : int (C int)
        Box size along z.
    lattice : (X_DIM, Y_DIM, Z_DIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid, read only to decide whether a site is empty (value 0)
        or occupied.

    Returns
    -------
    (n_interface, 2, 3) numpy.ndarray, NUMPY_INT_TYPE (int32)
        The subset of the 26 short-range pairs around the central site in which
        at least one of the two sites is empty. n_interface is at most 26.
    """
    cdef int i
    cdef int good_index = 0
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] center = np.empty((3,), dtype=NUMPY_INT_TYPE_PYTHON)
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] good_pairs = np.zeros((26,2,3), dtype=NUMPY_INT_TYPE_PYTHON)

    center[0] = E_x
    center[1] = E_y
    center[2] = E_z
    pairs = inner_loops.extract_SR_pairs_from_position_3D(center, X_DIM, Y_DIM, Z_DIM)

    for i in range(26):
        if lattice[pairs[i, 0, 0], pairs[i,0,1], pairs[i,0,2]] == 0 or \
           lattice[pairs[i, 1, 0], pairs[i,1,1], pairs[i,1,2]] == 0:
            good_pairs[good_index] = pairs[i]
            good_index = good_index + 1

    return good_pairs[0:good_index]


###
###
###
@cython.boundscheck(False)
@cython.wraparound(False) 
def get_unique_interface_pairs_2D(int E_x, int E_y, int X_DIM, int Y_DIM, NUMPY_INT_TYPE[:,:] lattice):
    """
    Returns the nearest neighbour set of pairs corresponding to the central position (E_x, E_y) where at least one of the
    two positions (central position and somewhere else) is an empty site. This is useful if you KNOW the central position is filled a
    priori, as it'll then return a non-redundant list of ONLY the pairs which correspond to the chain-solvent interface
    site.

    This is especially useful for rigid-body cluster moves because we KNOW that there cannot be components of the cluster in nearest
    neighbour contact with another chain which is not part of the cluster (otherwise it would be part of the cluster!) such that by
    getting the unique_interface_pairs for each position in the cluster we automatically have the full, non-redundant set of pairs
    over which we need to evaluate energy changes on a move - no need for any kind of post-processing.

    Parameters
    ----------
    E_x : int (C int)
        x coordinate of the central (assumed occupied) site.
    E_y : int (C int)
        y coordinate of the central site.
    X_DIM : int (C int)
        Box size along x, used to wrap the neighbour shell.
    Y_DIM : int (C int)
        Box size along y.
    lattice : (X_DIM, Y_DIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid, read only to decide whether a site is empty (value 0)
        or occupied.

    Returns
    -------
    (n_interface, 2, 2) numpy.ndarray, NUMPY_INT_TYPE (int32)
        The subset of the 8 short-range pairs around the central site in which
        at least one of the two sites is empty. n_interface is at most 8.
    """
    cdef int i
    cdef int good_index = 0
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] center = np.empty((2,), dtype=NUMPY_INT_TYPE_PYTHON)
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] pairs
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] good_pairs = np.zeros((8,2,2), dtype=NUMPY_INT_TYPE_PYTHON)

    center[0] = E_x
    center[1] = E_y
    pairs = inner_loops.extract_SR_pairs_from_position_2D(center, X_DIM, Y_DIM)

    for i in range(8):
        if lattice[pairs[i,0,0], pairs[i,0,1]] == 0 or \
           lattice[pairs[i,1,0], pairs[i,1,1]] == 0:
            good_pairs[good_index] = pairs[i]
            good_index = good_index + 1

    return good_pairs[0:good_index]
    



### -------------------------------------------------------------------------------------------------------------------------------------
### ENERGY

###
###
###
@cython.boundscheck(False)
def get_gridvalue_3D(NUMPY_INT_TYPE[:,:,:] lattice, unsigned int pos1, unsigned int pos2, unsigned int pos3):    
    """
    Read a single site out of a 3D grid.

    NOTE: positions are UNSIGNED and indexed with bounds checking off - a caller
    passing a negative/unwrapped coordinate gets a silent far-out-of-bounds read.
    Callers must pass in-box coordinates.

    Parameters
    ----------
    lattice : (X_DIM, Y_DIM, Z_DIM) NUMPY_INT_TYPE (int32) memoryview
        Grid to read from. This is used for both the occupancy grid and the type
        grid, so the meaning of the value returned depends on what is passed.
    pos1 : unsigned int
        x coordinate, already wrapped into the box.
    pos2 : unsigned int
        y coordinate, already wrapped into the box.
    pos3 : unsigned int
        z coordinate, already wrapped into the box.

    Returns
    -------
    NUMPY_INT_TYPE (int32)
        The value stored at that site.
    """
    return lattice[pos1,pos2,pos3]


###
###
###
@cython.boundscheck(False)
def get_gridvalue_2D(NUMPY_INT_TYPE[:,:] lattice, unsigned int pos1, unsigned int pos2):
    """
    Read a single site out of a 2D grid.

    NOTE: positions are UNSIGNED and indexed with bounds checking off - a caller
    passing a negative/unwrapped coordinate gets a silent far-out-of-bounds read.
    No live code path uses this (the hot path indexes numpy directly); callers
    must pass in-box coordinates.

    Parameters
    ----------
    lattice : (X_DIM, Y_DIM) NUMPY_INT_TYPE (int32) memoryview
        Grid to read from. This is used for both the occupancy grid and the type
        grid, so the meaning of the value returned depends on what is passed.
    pos1 : unsigned int
        x coordinate, already wrapped into the box.
    pos2 : unsigned int
        y coordinate, already wrapped into the box.

    Returns
    -------
    NUMPY_INT_TYPE (int32)
        The value stored at that site.
    """
    return lattice[pos1,pos2]


###
###
###
@cython.boundscheck(False)
def evaluate_local_energy_3D_non_shortrange(NUMPY_INT_TYPE[:,:,:] lattice, 
                                            NUMPY_INT_TYPE[:,:,:] pairs_list, 
                                            NUMPY_INT_TYPE[:,:] interaction_table, 
                                            unsigned int hardwall):
    """
    Evaluate 3D energy for LR and SLR interactions. If hardwall and straddling a PBC we *ignore* these interactions.

    A pair is taken to straddle the boundary when any per-axis separation
    exceeds 3 sites, which cannot happen for a genuine long-range pair inside
    the box.

    Parameters
    ----------
    lattice : (X_DIM, Y_DIM, Z_DIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid; the value at a site is the residue integer code, with 0
        meaning solvent (empty).
    pairs_list : (n_pairs, 2, 3) NUMPY_INT_TYPE (int32) memoryview
        Site pairs to sum over, each entry holding the two wrapped (x, y, z)
        coordinates.
    interaction_table : (n_types, n_types) NUMPY_INT_TYPE (int32) memoryview
        Pairwise energies indexed by the two residue integer codes. This is the
        LR or SLR table depending on which set of pairs is being evaluated.
    hardwall : unsigned int
        1 if hardwall boundaries are on, in which case boundary-straddling pairs
        are skipped; 0 for fully periodic, in which case every pair is summed.

    Returns
    -------
    long
        Summed interaction energy over the supplied pairs.
    """

    cdef long ENERGY = 0   # long: the from-scratch oracle must not wrap for very large systems
    cdef unsigned int i, typeA, typeB, num_pairs;
    cdef int dx, dy, dz   # C-int separations: a Python abs() on int32 buffer elements boxed per pair

    num_pairs = len(pairs_list)

    if hardwall == 1:
        for i in range(num_pairs):        

            # if we're looking at a pair that stradles a PBC 
            dx = pairs_list[i,0,0] - pairs_list[i,1,0]; dy = pairs_list[i,0,1] - pairs_list[i,1,1]; dz = pairs_list[i,0,2] - pairs_list[i,1,2]
            if dx > 3 or dx < -3 or dy > 3 or dy < -3 or dz > 3 or dz < -3:
                continue
            else:
                typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2]]
                typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2]]

                ENERGY = ENERGY + interaction_table[typeA,typeB]


    else:
        for i in range(num_pairs):        

            typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2]]
            typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2]]

            ENERGY = ENERGY + interaction_table[typeA,typeB]
        
    return ENERGY


###
###
###

# precision fix
#def evaluate_local_energy_3D_shortrange(cnp.ndarray[NUMPY_INT_TYPE, ndim=3] lattice, cnp.ndarray[NUMPY_INT_TYPE, ndim=3] pairs_list, cnp.ndarray[np.float_t, ndim=2] interaction_table, unsigned int hardwall):
@cython.boundscheck(False)
def evaluate_local_energy_3D_shortrange(NUMPY_INT_TYPE[:,:,:] lattice, 
                                        NUMPY_INT_TYPE[:,:,:] pairs_list, 
                                        NUMPY_INT_TYPE[:,:] interaction_table, 
                                        unsigned int hardwall):
    """
    Evaluate 3D energy for short-range (nearest-neighbour) interactions.

    With hardwall on, a pair that straddles the boundary is not simply dropped:
    each non-solvent bead in that pair instead gets its solvent interaction,
    which is what a bead sitting against a wall actually sees. A pair is taken
    to straddle the boundary when any per-axis separation exceeds 3 sites.

    Parameters
    ----------
    lattice : (X_DIM, Y_DIM, Z_DIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid; the value at a site is the residue integer code, with 0
        meaning solvent (empty).
    pairs_list : (n_pairs, 2, 3) NUMPY_INT_TYPE (int32) memoryview
        Site pairs to sum over, each entry holding the two wrapped (x, y, z)
        coordinates.
    interaction_table : (n_types, n_types) NUMPY_INT_TYPE (int32) memoryview
        Pairwise short-range energies indexed by the two residue integer codes.
        Column 0 holds the solvent interactions used at a hardwall.
    hardwall : unsigned int
        1 if hardwall boundaries are on, in which case boundary-straddling pairs
        are replaced by solvent interactions; 0 for fully periodic, in which
        case every pair is summed as-is.

    Returns
    -------
    long
        Summed interaction energy over the supplied pairs.
    """

    # precision fix
    # cdef float ENERGY = 0
    cdef long ENERGY = 0   # long: the from-scratch oracle must not wrap for very large systems
    cdef unsigned int i, typeA, typeB, num_pairs;
    cdef int dx, dy, dz   # C-int separations: a Python abs() on int32 buffer elements boxed per pair

    num_pairs = len(pairs_list)

    # if hardwall is ON we do not want to compute interactions across the PBD
    if hardwall == 1:
        for i in range(num_pairs):        

            typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2]]
            typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2]]


            # if we're looking at a pair that stradles a PBC 
            dx = pairs_list[i,0,0] - pairs_list[i,1,0]; dy = pairs_list[i,0,1] - pairs_list[i,1,1]; dz = pairs_list[i,0,2] - pairs_list[i,1,2]
            if dx > 3 or dx < -3 or dy > 3 or dy < -3 or dz > 3 or dz < -3:

                # if A is non-solvent then A is interacting with solvent 
                if typeA > 0:
                    ENERGY = ENERGY + interaction_table[typeA, 0]
                # if B is non-solvent then B is interacting with solvent
                if typeB > 0:
                    ENERGY = ENERGY + interaction_table[typeB, 0]

            else:
                ENERGY = ENERGY + interaction_table[typeA,typeB]


    # if hardwall is off then we compute every pair!
    else:
        for i in range(num_pairs):        

            typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2]]
            typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2]]

            #print("%s + %s" %(str(ENERGY), str(interaction_table[typeA,typeB])))
            ENERGY = ENERGY + interaction_table[typeA,typeB]
            #print(ENERGY)
            
    return ENERGY


###
###
###
@cython.boundscheck(False)
def evaluate_local_energy_2D_shortrange(NUMPY_INT_TYPE[:,:] lattice, 
                                        NUMPY_INT_TYPE[:,:,:] pairs_list, 
                                        NUMPY_INT_TYPE[:,:] interaction_table, 
                                        unsigned int hardwall):
    """
    Using a defined set of redundant pairs defines the interaction between the positions on the typegrid lattice as defined by the interaction table. This is a general purpose non-bonded
    interaction function.

    With hardwall on, a pair that straddles the boundary is not simply dropped:
    each non-solvent bead in that pair instead gets its solvent interaction,
    which is what a bead sitting against a wall actually sees. A pair is taken
    to straddle the boundary when either per-axis separation exceeds 3 sites.

    Parameters
    ----------
    lattice : (X_DIM, Y_DIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid; the value at a site is the residue integer code, with 0
        meaning solvent (empty).
    pairs_list : (n_pairs, 2, 2) NUMPY_INT_TYPE (int32) memoryview
        Site pairs to sum over, each entry holding the two wrapped (x, y)
        coordinates.
    interaction_table : (n_types, n_types) NUMPY_INT_TYPE (int32) memoryview
        Pairwise short-range energies indexed by the two residue integer codes.
        Column 0 holds the solvent interactions used at a hardwall.
    hardwall : unsigned int
        1 if hardwall boundaries are on, in which case boundary-straddling pairs
        are replaced by solvent interactions; 0 for fully periodic, in which
        case every pair is summed as-is.

    Returns
    -------
    long
        Summed interaction energy over the supplied pairs.
    """
    cdef long ENERGY = 0   # long: the from-scratch oracle must not wrap for very large systems
    cdef unsigned int i, typeA, typeB, num_pairs;
    cdef int dx, dy   # C-int separations: a Python abs() on int32 buffer elements boxed per pair

    num_pairs = len(pairs_list)

    # if a hardwall boundary is being applied
    if hardwall == 1:
        for i in range(num_pairs):

            typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1]]
            typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1]]

            # if we're looking at a pair that stradles a PBC 
            dx = pairs_list[i,0,0] - pairs_list[i,1,0]; dy = pairs_list[i,0,1] - pairs_list[i,1,1]
            if dx > 3 or dx < -3 or dy > 3 or dy < -3:

                # if A is non-solvent then A is interacting with solvent 
                if typeA > 0:
                    ENERGY = ENERGY + interaction_table[typeA, 0]
                # if B is non-solvent then B is interacting with solvent
                if typeB > 0:
                    ENERGY = ENERGY + interaction_table[typeB, 0]

            else:
                ENERGY = ENERGY + interaction_table[typeA,typeB]

    else:
        for i in range(num_pairs):
            typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1]]
            typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1]]

            ENERGY = ENERGY + interaction_table[typeA,typeB]
        
    return ENERGY



@cython.boundscheck(False)
def evaluate_local_energy_2D_non_shortrange(NUMPY_INT_TYPE[:, :] lattice, 
                                            NUMPY_INT_TYPE[:, :, :] pairs_list, 
                                            NUMPY_INT_TYPE[:, :] interaction_table, 
                                            unsigned int hardwall):
    """
    Using a defined set of redundant pairs defines the interaction between the positions on the typegrid lattice as defined by the interaction table. This is a general purpose non-bonded
    interaction function.

    This is the LR/SLR variant: with hardwall on, a pair that straddles the
    boundary is dropped outright rather than converted into solvent contacts. A
    pair is taken to straddle the boundary when either per-axis separation
    exceeds 3 sites.

    Parameters
    ----------
    lattice : (X_DIM, Y_DIM) NUMPY_INT_TYPE (int32) memoryview
        The type grid; the value at a site is the residue integer code, with 0
        meaning solvent (empty).
    pairs_list : (n_pairs, 2, 2) NUMPY_INT_TYPE (int32) memoryview
        Site pairs to sum over, each entry holding the two wrapped (x, y)
        coordinates.
    interaction_table : (n_types, n_types) NUMPY_INT_TYPE (int32) memoryview
        Pairwise energies indexed by the two residue integer codes. This is the
        LR or SLR table depending on which set of pairs is being evaluated.
    hardwall : unsigned int
        1 if hardwall boundaries are on, in which case boundary-straddling pairs
        are skipped; 0 for fully periodic, in which case every pair is summed.

    Returns
    -------
    long
        Summed interaction energy over the supplied pairs.
    """
    cdef long ENERGY = 0   # long: the from-scratch oracle must not wrap for very large systems
    cdef unsigned int i, typeA, typeB, num_pairs;
    cdef int dx, dy   # C-int separations: a Python abs() on int32 buffer elements boxed per pair

    num_pairs = len(pairs_list)

    # if a hardwall boundary is being applied
    if hardwall == 1:
        for i in range(num_pairs):

            # skip interactions that jump PBC
            dx = pairs_list[i,0,0] - pairs_list[i,1,0]; dy = pairs_list[i,0,1] - pairs_list[i,1,1]
            if dx > 3 or dx < -3 or dy > 3 or dy < -3:
                continue

            else:                
                typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1]]
                typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1]]

                ENERGY = ENERGY + interaction_table[typeA,typeB]

    else:
        for i in range(num_pairs):
            typeA = lattice[pairs_list[i,0,0], pairs_list[i,0,1]]
            typeB = lattice[pairs_list[i,1,0], pairs_list[i,1,1]]

            ENERGY = ENERGY + interaction_table[typeA,typeB]
        
    return ENERGY


###
###
###
@cython.boundscheck(False)
def evaluate_angle_energy_3D(NUMPY_INT_TYPE[:,:] chain_positions, 
                             NUMPY_INT_TYPE[:] intcode_sequence, 
                             NUMPY_INT_TYPE[:,:,:,:,:,:,:] angle_lookup,
                             int chain_length):
    """
    Sum the 3D angle penalty over every i, i+1, i+2 triplet in a chain.

    The two bond vectors are taken relative to the middle bead and clamped back
    to unit steps, which is how a bond that wraps the periodic boundary (an
    apparent step of size DIM - 1) is folded back to its true direction. The
    penalty is looked up per triplet using the integer code of the MIDDLE bead.

    Parameters
    ----------
    chain_positions : (chain_length, 3) NUMPY_INT_TYPE (int32) memoryview
        Bead positions along the chain, in bead order.
    intcode_sequence : (chain_length,) NUMPY_INT_TYPE (int32) memoryview
        Per-bead residue integer codes, used to select the residue-specific
        slice of angle_lookup.
    angle_lookup : (n_types, 3, 3, 3, 3, 3, 3) NUMPY_INT_TYPE (int32) memoryview
        Precomputed angle penalties indexed by
        [intcode, a0+1, a1+1, a2+1, b0+1, b1+1, b2+1], where a and b are the two
        unit bond vectors and the +1 shifts the -1/0/+1 components into 0/1/2.
    chain_length : int (C int)
        Number of beads in the chain. chain_length - 2 triplets are summed, so
        a chain shorter than 3 beads contributes nothing.

    Returns
    -------
    long
        Total angle penalty for the chain.
    """


    cdef long ENERGY = 0   # long: the from-scratch oracle must not wrap for very large systems
    cdef int i
    cdef int a0, a1, a2, b0, b1, b2

    

    for i in range(chain_length-2):

        # get relative i-th position (relative to i+1)
        a0 = chain_positions[i+1, 0] - chain_positions[i, 0]
        a1 = chain_positions[i+1, 1] - chain_positions[i, 1]
        a2 = chain_positions[i+1, 2] - chain_positions[i, 2]
        if a0 < -1:
            a0 = 1
        elif a0 > 1:
            a0 = -1
        if a1 < -1:
            a1 = 1
        elif a1 > 1:
            a1 = -1
        if a2 < -1:
            a2 = 1
        elif a2 > 1:
            a2 = -1

        
        # get relative 
        b0 = chain_positions[i+1, 0] - chain_positions[i+2, 0]
        b1 = chain_positions[i+1, 1] - chain_positions[i+2, 1]
        b2 = chain_positions[i+1, 2] - chain_positions[i+2, 2]
        if b0 < -1:
            b0 = 1
        elif b0 > 1:
            b0 = -1
        if b1 < -1:
            b1 = 1
        elif b1 > 1:
            b1 = -1
        if b2 < -1:
            b2 = 1
        elif b2 > 1:
            b2 = -1


        # the energy looks up the 'i+1'th residue because the angle is over three residues
        # where those three are i, i+1, i+2 (i.e. i+1 is the middle residue)
        #print intcode_sequence[i+1]
        #print np.shape(angle_lookup)

        
        #print "Bead -1 %s" %(a)
        #print "Bead +1 %s" %(b)
        #print "Bead pos: [%i,%i,%i],[%i,%i,%i],[%i,%i,%i] " % (chain_positions[i,0], chain_positions[i,1], chain_positions[i,2], chain_positions[i+1,0], chain_positions[i+1,1], chain_positions[i+1,2], chain_positions[i+2,0], chain_positions[i+2,1], chain_positions[i+2,2])
        #print "Bead %i -> distance [%i, %i, %i]: %i" % (i, abs(a[0]-b[0]), abs(a[1]-b[1]), abs(a[2]-b[2]), angle_lookup[intcode_sequence[i+1], a[0]+1, a[1]+1, a[2]+1, b[0]+1, b[1]+1, b[2]+1])
        
        ENERGY = angle_lookup[intcode_sequence[i+1], a0+1, a1+1, a2+1, b0+1, b1+1, b2+1] + ENERGY
        #print "%i = %i a=[%i,%i,%i], b=[%i,%i,%i]" % (i, energy, a[0]+1, a[1]+1,a[2]+1,b[0]+1, b[1]+1,b[2]+1)

    return ENERGY

@cython.boundscheck(False)
def evaluate_angle_energy_2D(NUMPY_INT_TYPE[:,:] chain_positions, 
                             NUMPY_INT_TYPE[:] intcode_sequence, 
                             NUMPY_INT_TYPE[:,:,:,:,:] angle_lookup,
                             int chain_length):
    """
    Sum the 2D angle penalty over every i, i+1, i+2 triplet in a chain.

    The two bond vectors are taken relative to the middle bead and clamped back
    to unit steps, which is how a bond that wraps the periodic boundary (an
    apparent step of size DIM - 1) is folded back to its true direction. The
    penalty is looked up per triplet using the integer code of the MIDDLE bead.

    Parameters
    ----------
    chain_positions : (chain_length, 2) NUMPY_INT_TYPE (int32) memoryview
        Bead positions along the chain, in bead order.
    intcode_sequence : (chain_length,) NUMPY_INT_TYPE (int32) memoryview
        Per-bead residue integer codes, used to select the residue-specific
        slice of angle_lookup.
    angle_lookup : (n_types, 3, 3, 3, 3) NUMPY_INT_TYPE (int32) memoryview
        Precomputed angle penalties indexed by [intcode, a0+1, a1+1, b0+1, b1+1],
        where a and b are the two unit bond vectors and the +1 shifts the
        -1/0/+1 components into 0/1/2.
    chain_length : int (C int)
        Number of beads in the chain. chain_length - 2 triplets are summed, so
        a chain shorter than 3 beads contributes nothing.

    Returns
    -------
    long
        Total angle penalty for the chain.
    """


    cdef long ENERGY = 0   # long: the from-scratch oracle must not wrap for very large systems
    cdef int i
    cdef int a0, a1, b0, b1

    for i in range(chain_length-2):
        a0 = chain_positions[i+1, 0] - chain_positions[i, 0]
        a1 = chain_positions[i+1, 1] - chain_positions[i, 1]
        if a0 < -1:
            a0 = 1
        elif a0 > 1:
            a0 = -1
        if a1 < -1:
            a1 = 1
        elif a1 > 1:
            a1 = -1

        
        b0 = chain_positions[i+1, 0] - chain_positions[i+2, 0]
        b1 = chain_positions[i+1, 1] - chain_positions[i+2, 1]
        if b0 < -1:
            b0 = 1
        elif b0 > 1:
            b0 = -1
        if b1 < -1:
            b1 = 1
        elif b1 > 1:
            b1 = -1

        ENERGY = angle_lookup[intcode_sequence[i+1], a0+1, a1+1, b0+1, b1+1] + ENERGY

    return ENERGY
