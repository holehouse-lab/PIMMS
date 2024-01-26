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
    Note this INCLUDES the original site and does perform PBC correction
    

    Optimization notes
    
    - Switching out the while loops for straightup list of position assignment did
      not bring any speedup

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
            positions[array_index][0] = (position1+x) % X_DIM
            positions[array_index][1] = (position2+y) % Y_DIM
                
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
    Note this INCLUDES the original site and does perform PBC correction
    

    Optimization notes
    
    - Switching out the while loops for straightup list of position assignment did
      not bring any speedup

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
    

    """
    cdef int i = 0;
    cdef int index;
    #cdef cpplist[int] indices;
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] pairs = np.zeros((26,2,3), dtype=NUMPY_INT_TYPE_PYTHON)
    cdef unsigned int num_pairs;
    
    pairs = inner_loops.extract_SR_pairs_from_position_3D(np.array([E_x, E_y, E_z]), X_DIM, Y_DIM, Z_DIM)
    
    ## get the indices of pairs where one of the pair is an empty lattice site
    good_pairs = []
    while i < 26:

        # if pair position 1 is an empty size
        if lattice[pairs[i, 0, 0], pairs[i,0,1], pairs[i,0,2]] == 0:
            good_pairs.append(pairs[i])
            
        # else if pair position 2 is an empty site
        elif lattice[pairs[i, 1, 0], pairs[i,1,1], pairs[i,1, 2]] == 0:
            good_pairs.append(pairs[i])
            
        i=i+1
        
    return good_pairs


###
###
###
@cython.boundscheck(False)
@cython.wraparound(False) 
def get_unique_interface_pairs_2D(int E_x, int E_y, int X_DIM, int Y_DIM, NUMPY_INT_TYPE[:,:] lattice):
    """
    Returns the nearest neighbour set of pairs corresponding to the central position (E_x, E_y, E_z) where at least one of the 
    two positions (central position and somewhere else) is an empty site. This is useful if you KNOW the central position is filled a
    priori, as it'll then return a non-redundant list of ONLY the pairs which correspond to the chain-solvent interface
    site.

    This is especially useful for rigid-body cluster moves because we KNOW that there cannot be components of the cluster in nearest
    neighbour contact with another chain which is not part of the cluster (otherwise it would be part of the cluster!) such that by
    getting the unique_interface_pairs for each position in the cluster we automatically have the full, non-redundant set of pairs
    over which we need to evaluate energy changes on a move - no need for any kind of post-processing.
    

    """
    cdef int i = 0;
    cdef int index;

    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=3] pairs = np.zeros((8,2,2), dtype=NUMPY_INT_TYPE_PYTHON);
    cdef unsigned int num_pairs;

    pairs = inner_loops.extract_SR_pairs_from_position_3D(np.array([E_x, E_y]), X_DIM, Y_DIM)
    

    # get the indices of pairs where one of the pair is an empty lattice site
    # indicies is adynamically allocated list of the indices from the pairs
    # array where at least one of the two pairs is an empty lattice site
    
    good_pairs = []
    while i < 8:
        # if pair position 1 is an empty size
        if lattice[pairs[i,0,0], pairs[i, 0, 1]] == 0:
            good_pairs.append(pairs[i])
            
        # else if pair position 2 is an empty site
        elif lattice[pairs[i, 1, 0], pairs[i, 1, 1]] == 0:
            good_pairs.append(pairs[i])

        i=i+1
        
    return good_pairs
    



### -------------------------------------------------------------------------------------------------------------------------------------
### ENERGY

###
###
###
@cython.boundscheck(False)
def get_gridvalue_3D(NUMPY_INT_TYPE[:,:,:] lattice, unsigned int pos1, unsigned int pos2, unsigned int pos3):    
    return lattice[pos1,pos2,pos3]


###
###
###
@cython.boundscheck(False)
def get_gridvalue_2D(NUMPY_INT_TYPE[:,:] lattice, unsigned int pos1, unsigned int pos2):
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
    """

    cdef int ENERGY = 0
    cdef unsigned int i, typeA, typeB, num_pairs;

    num_pairs = len(pairs_list)

    if hardwall == 1:
        for i in range(num_pairs):        

            # if we're looking at a pair that stradles a PBC 
            if (abs(pairs_list[i,0,0] - pairs_list[i,1,0]) > 3) or (abs(pairs_list[i,0,1] - pairs_list[i,1,1]) > 3) or (abs(pairs_list[i,0,2] - pairs_list[i,1,2]) > 3):
                continue
            else:
                typeA = get_gridvalue_3D(lattice, pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2])
                typeB = get_gridvalue_3D(lattice, pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2])

                ENERGY = ENERGY + interaction_table[typeA,typeB]


    else:
        for i in range(num_pairs):        

            typeA = get_gridvalue_3D(lattice, pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2])
            typeB = get_gridvalue_3D(lattice, pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2])

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

    # precision fix
    # cdef float ENERGY = 0
    cdef int ENERGY = 0
    cdef unsigned int i, typeA, typeB, num_pairs;

    num_pairs = len(pairs_list)

    # if hardwall is ON we do not want to compute interactions across the PBD
    if hardwall == 1:
        for i in range(num_pairs):        

            typeA = get_gridvalue_3D(lattice, pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2])
            typeB = get_gridvalue_3D(lattice, pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2])


            # if we're looking at a pair that stradles a PBC 
            if (abs(pairs_list[i,0,0] - pairs_list[i,1,0]) > 3) or (abs(pairs_list[i,0,1] - pairs_list[i,1,1]) > 3) or (abs(pairs_list[i,0,2] - pairs_list[i,1,2]) > 3):

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

            typeA = get_gridvalue_3D(lattice, pairs_list[i,0,0], pairs_list[i,0,1], pairs_list[i,0,2])
            typeB = get_gridvalue_3D(lattice, pairs_list[i,1,0], pairs_list[i,1,1], pairs_list[i,1,2])

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

    """
    cdef int ENERGY = 0
    cdef unsigned int i, typeA, typeB, num_pairs;

    num_pairs = len(pairs_list)

    # if a hardwall boundary is being applied
    if hardwall == 1:
        for i in range(num_pairs):

            typeA = get_gridvalue_2D(lattice, pairs_list[i,0,0], pairs_list[i,0,1])
            typeB = get_gridvalue_2D(lattice, pairs_list[i,1,0], pairs_list[i,1,1])

            # if we're looking at a pair that stradles a PBC 
            if (abs(pairs_list[i,0,0] - pairs_list[i,1,0]) > 3) or (abs(pairs_list[i,0,1] - pairs_list[i,1,1]) > 3):

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
            typeA = get_gridvalue_2D(lattice, pairs_list[i,0,0], pairs_list[i,0,1])
            typeB = get_gridvalue_2D(lattice, pairs_list[i,1,0], pairs_list[i,1,1])

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

    """
    cdef int ENERGY = 0
    cdef unsigned int i, typeA, typeB, num_pairs;

    num_pairs = len(pairs_list)

    # if a hardwall boundary is being applied
    if hardwall == 1:
        for i in range(num_pairs):

            # skip interactions that jump PBC
            if (abs(pairs_list[i,0,0] - pairs_list[i,1,0]) > 3) or (abs(pairs_list[i,0,1] - pairs_list[i,1,1]) > 3):
                continue

            else:                
                typeA = get_gridvalue_2D(lattice, pairs_list[i,0,0], pairs_list[i,0,1])
                typeB = get_gridvalue_2D(lattice, pairs_list[i,1,0], pairs_list[i,1,1])

                ENERGY = ENERGY + interaction_table[typeA,typeB]

    else:
        for i in range(num_pairs):
            typeA = get_gridvalue_2D(lattice, pairs_list[i,0,0], pairs_list[i,0,1])
            typeB = get_gridvalue_2D(lattice, pairs_list[i,1,0], pairs_list[i,1,1])

            ENERGY = ENERGY + interaction_table[typeA,typeB]
        
    return ENERGY


###
###
###
#@cython.boundscheck(False)
#@cython.cdivision(True)
def evaluate_angle_energy_3D(NUMPY_INT_TYPE[:,:] chain_positions, 
                             NUMPY_INT_TYPE[:] intcode_sequence, 
                             NUMPY_INT_TYPE[:,:,:,:,:,:,:] angle_lookup,
                             int chain_length):


    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] a = np.zeros([3], dtype=NUMPY_INT_TYPE_PYTHON)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] b = np.zeros([3], dtype=NUMPY_INT_TYPE_PYTHON)    
    cdef int ENERGY = 0;
    cdef int i;

    

    for i in range(chain_length-2):

        # get relative i-th position (relative to i+1)
        a[0] = fix_angle_pbc_issues(chain_positions[i+1, 0] - chain_positions[i, 0])
        a[1] = fix_angle_pbc_issues(chain_positions[i+1, 1] - chain_positions[i, 1])
        a[2] = fix_angle_pbc_issues(chain_positions[i+1, 2] - chain_positions[i, 2])

        
        # get relative 
        b[0] = fix_angle_pbc_issues(chain_positions[i+1, 0] - chain_positions[i+2, 0])
        b[1] = fix_angle_pbc_issues(chain_positions[i+1, 1] - chain_positions[i+2, 1])
        b[2] = fix_angle_pbc_issues(chain_positions[i+1, 2] - chain_positions[i+2, 2])


        # the energy looks up the 'i+1'th residue because the angle is over three residues
        # where those three are i, i+1, i+2 (i.e. i+1 is the middle residue)
        #print intcode_sequence[i+1]
        #print np.shape(angle_lookup)

        
        #print "Bead -1 %s" %(a)
        #print "Bead +1 %s" %(b)
        #print "Bead pos: [%i,%i,%i],[%i,%i,%i],[%i,%i,%i] " % (chain_positions[i,0], chain_positions[i,1], chain_positions[i,2], chain_positions[i+1,0], chain_positions[i+1,1], chain_positions[i+1,2], chain_positions[i+2,0], chain_positions[i+2,1], chain_positions[i+2,2])
        #print "Bead %i -> distance [%i, %i, %i]: %i" % (i, abs(a[0]-b[0]), abs(a[1]-b[1]), abs(a[2]-b[2]), angle_lookup[intcode_sequence[i+1], a[0]+1, a[1]+1, a[2]+1, b[0]+1, b[1]+1, b[2]+1])
        
        ENERGY = angle_lookup[intcode_sequence[i+1], a[0]+1, a[1]+1, a[2]+1, b[0]+1, b[1]+1, b[2]+1] + ENERGY
        #print "%i = %i a=[%i,%i,%i], b=[%i,%i,%i]" % (i, energy, a[0]+1, a[1]+1,a[2]+1,b[0]+1, b[1]+1,b[2]+1)

    return ENERGY

@cython.boundscheck(False)
@cython.cdivision(True)
def evaluate_angle_energy_2D(NUMPY_INT_TYPE[:,:] chain_positions, 
                             NUMPY_INT_TYPE[:] intcode_sequence, 
                             NUMPY_INT_TYPE[:,:,:,:,:] angle_lookup,
                             int chain_length):


    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] a = np.zeros([2], dtype=NUMPY_INT_TYPE_PYTHON)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] b = np.zeros([2], dtype=NUMPY_INT_TYPE_PYTHON)    
    cdef int ENERGY  = 0;
    cdef unsigned int i;

    for i in range(chain_length-2):
        a[0] = fix_angle_pbc_issues(chain_positions[i+1, 0] - chain_positions[i, 0])
        a[1] = fix_angle_pbc_issues(chain_positions[i+1, 1] - chain_positions[i, 1])

        
        b[0] = fix_angle_pbc_issues(chain_positions[i+1, 0] - chain_positions[i+2, 0])
        b[1] = fix_angle_pbc_issues(chain_positions[i+1, 1] - chain_positions[i+2, 1])

        ENERGY = angle_lookup[intcode_sequence[i+1], a[0]+1, a[1]+1, b[0]+1, b[1]+1] + ENERGY

    return ENERGY


    

@cython.cdivision(True)
cdef int fix_angle_pbc_issues(int distance):
    """
    Hack that takes advantage of the fact that the distances all must be -1,0, +1 and 
    we're always computing a X->O vector where O is the bend point

    """
    if distance < -1:
        return 1
    elif distance > 1:
        return -1
    else:
        return distance
        





    





    



