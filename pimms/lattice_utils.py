## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2022
## ...........................................................................


##
## lattice_utils
##
## lattice_utils contains system agnostic, stateless utilities for lattice operations.
## Utilities are relevant for 2D and 3D lattices
##

import random
import copy
import os
import numpy as np

import mdtraj as md

from .latticeExceptions import ChainInsertionFailure, ChainDeletionFailure, ResidueAugmentException, DebuggingException, MoveSetException, ChainConnectivityError, ClusterSizeThresholdException
#from .pdb_utils import build_pdb_file, finalize_pdb_file, initialize_pdb_file
from . import pdb_utils

from . import hyperloop
from . import inner_loops
from . import inner_loops_hardwall

from . import numpy_utils
from . import lattice_tools
from . import lattice_analysis_utils
from . import IO_utils

from . import CONFIG

#from CONFIG import * # note there are things from CONFIG being used...

#-----------------------------------------------------------------
#
def same_sites(site1, site2):
    """
    Explicit function to check two sites are the same. Automatically generalizes
    to 2D or 3D in response to the site dimensions.    


    Parameters
    ----------
    site1 : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied.

    site2: list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied.


    Returns
    ---------
    bool
        Returns True if the two sites are the same and False if not

    """

    # two dimensions case
    if len(site1) == 2:

        if site1[0] == site2[0] and site1[1] == site2[1]:
            return True
        else: 
            return False

    # three dimensions case
    else:

        if site1[0] == site2[0] and site1[1] == site2[1] and site1[2] == site2[2]:
            return True
        else: 
            return False


#-----------------------------------------------------------------
#
def get_real_distance(posA, posB, dimensions):
    """
    Function to calculate the real distance between two positions.


    Parameters
    ----------
    posA : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects a specific position on the lattice.

    posB : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects a specific position on the lattice.

    dimensions : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects the lattice dimensions.


    Returns
    ---------
    float
        Returns a value that reflects the Euclidian distance between two positions
        on the lattice.


    """
    return lattice_analysis_utils.get_inter_position_distance(posA, posB, dimensions)
       


#----------------------------------------------------------------
#
def get_dimensions(lattice_grid):
    """
    Function that returns the dimensions associated with the lattice

    Parameters
    ----------
    lattice_grid : np.array
        A lattice grid numpy array


    Returns
    ---------
    np.array
        Returns a numpy array of length 2 or 3, depending on the dimensionality
        of the system, where the value at each position reflects the size of the
        lattice in that dimension.
        

    """
    return lattice_grid.shape



#-----------------------------------------------------------------
#
def pbc_convert(position, dimensions):
    """
    Returns lattice site positions after carrying out periodic boundary conditions (PBC)    
    conversions.

    Parameters
    -----------
    position : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects a specific position on the lattice.        

    Returns
    ---------
    list
        A list where the length matches the input `position` list, where the new
        positions reflects the periodic-boundary condition corrected positions.

    """

    pbc_pos=[]
    n_dim = len(position)
        
    # cyle through each dimension and correct
    for idx in range(0,n_dim):
        pbc_pos.append(position[idx]%dimensions[idx])

    return pbc_pos



#-----------------------------------------------------------------
#
def pbc_correct(posA, posB, dimensions):
    """
    Returns the two positions after converting them such that, relative
    to one another, they are in the same image (i.e. not cross a border).
    NOTE that this function ONLY allows works for a single periodic image
    but if you had some weird set of positions that are spanning 
    multiple images this won't work. This is a pretty unlikely but
    just FYI.

    For simplicity the method assumes posA is God and posB can be re-
    set. This is abirtrary, but we don't need to futz with them both!!
    
    NOTE: This method was introduced in 0.9 so has not been as throughly
    vetted as a lot of the core code in PIMMS. 

    """

    newB = []
    n_dim = len(posA)

    # So I have a Cython implementation of the algorithm below, but it's about 1.8 * slower - presumably because loading
    # the data into memory is more expensive than the operation
    # newB = lattice_tools.pbc_correct_3D(np.array(posA,dtype=int), np.array(posB,dtype=int), np.array(dimensions,dtype=int))
            
    # for each pair of in each dimension
    for idx in range(0, n_dim):

        # if those positions are over half the boxwidth away then the minimum
        # distance is across the PBC with posA 
        if posA[idx] - posB[idx] > dimensions[idx]/2:
            newB.append(posB[idx] + dimensions[idx])

        elif posA[idx] - posB[idx] < -dimensions[idx]/2:
            newB.append(posB[idx] - dimensions[idx])

        else:
            newB.append(posB[idx])
    
    return (posA, newB)



#-----------------------------------------------------------------
#
def do_positions_stradle_pbc_boundary(chain_positions):
    """
    For a set of positions returns true if the positions straddle a boundary
    else return False. Note this assumes the positions are connected to one
    another (i.e. each adjacent position is only 1 lattice site away from
    the next).

    """

    n_dim = len(chain_positions[0])

    for pidx in range(0, len(chain_positions)-1):
        p1 = chain_positions[pidx]
        p2 = chain_positions[pidx+1]
        for xyz in range(0,n_dim):
            if abs(p1[xyz]-p2[xyz]) > 1:
                return True

    return False



#-----------------------------------------------------------------
#
def convert_chain_to_single_image(chain_of_positions, dimensions):
    """
    If passed a chain of positions converts them so, relative to the first position,
    they're all in the same image (single image convention). This assumes that
    each consecutive position is one lattice site apart from the next lattice site.

    chain_of_positions should be a list of 2D or 3D positions, note that
    the first position in this list will define the PBC frame being used. We may
    want (in the future) to offer the option to use the terminal position, and then
    assess which reference frame is the least perturbative to the chain. 

    Also worth bearing in mind that if this method is applied to a set of chains 
    in a connected cluster it will probably break everything. We have specific
    algorithm (snakesearch) to solve this problem, which can be found in the
    cluster_utils module.

    dimensions should be list giving the dimensions of the box
    ([X,Y], or [X,Y,Z])

    """

    local_positions = copy.deepcopy(chain_of_positions)
    
    n_dim = len(local_positions[0])
    n_pos = len(chain_of_positions)
        
    positions = []
    current = []

    # initially set the first position (i.e. the first residue
    # in the chain is used as a general reference, and then each
    # residue in turn is referenced to the preceding residue)
    positions.append(local_positions[0])

    # then initialize the 'current' position in the correct
    # number of dimensions
    for dim in range(0, n_dim):
        current.append(local_positions[0][dim])
        
    # for each position update the x/y/z positions of each residue
    # using chain connectivity to guide 
    for pidx in range(1, n_pos):
            
        next_pos = [0]*n_dim # this is just initializing 

        for dim in range(0, n_dim):                

            # if the next position and the current are the same
            if current[dim]  == local_positions[pidx][dim]:
                next_pos[dim] = local_positions[pidx][dim] 
                    
            # i.e. if we're at the edge of a boundary and the next postion is the 
            # other the side (e.g 27-28-[29-0]-1-2)
            elif current[dim] - local_positions[pidx][dim] > 1:                                        
                next_pos[dim] = local_positions[pidx][dim] + dimensions[dim]

                # finally this loop lets us account for chains that span multiple periodic images
                while abs(next_pos[dim] - current[dim]) > 1:
                    next_pos[dim] = next_pos[dim] + dimensions[dim]

            # i.e. if we're at the edge of a boundary and the next postion is the 
            # other the side (e.g 2-1-[0-29]-28-27)
            elif current[dim] - local_positions[pidx][dim] < -1:
                next_pos[dim] = local_positions[pidx][dim] - dimensions[dim]

                # finally this loop lets us account for chains that span multiple periodic images
                while abs(next_pos[dim] - current[dim]) > 1:
                    next_pos[dim] = next_pos[dim] - dimensions[dim]

            # just a simple next-position relationship
            else:
                next_pos[dim] = local_positions[pidx][dim] 
                    
        # Having determined position in each dimension we update the current 
        # position and append the next position to the ever growing list of 
        # newly updated positions
        positions.append(next_pos)            
        current = next_pos

    # uncorrected = copy.deepcopy(positions)
    # finally offset so all positions are positive        
    minDimVal = []
    for dim in range(0, n_dim):
        min_position = min(np.transpose(positions)[dim])
        if min_position < 0:                
            minDimVal.append(abs(min_position))
        else:
            minDimVal.append(0)
            
    for pidx in range(0, n_pos):
        for dim in range(0, n_dim):
            positions[pidx][dim] = positions[pidx][dim] + minDimVal[dim]

    return positions



#-----------------------------------------------------------------
#
def get_adjacent_sites_3D(position1, position2, position3, dimensions, extent_range=1):
    """
    Returns a 4D numpy array where each element is a 3D np array position

    """
    return(hyperloop.get_adjacent_sites_3D(position1, position2, position3, dimensions[0], dimensions[1], dimensions[2], extent_range))



#-----------------------------------------------------------------
#
def get_adjacent_sites_2D(position1, position2, dimensions, extent_range=1):
    """
    Returns a 3D numpy array where each element is a 2D np array position

    """
    return(hyperloop.get_adjacent_sites_2D(position1, position2, dimensions[0], dimensions[1], extent_range))



#-----------------------------------------------------------------
#
def find_nearest_position(target, positions_list, dimensions):
    """
    Given some target position (target) what is the index
    of the position in the positions list that is closest
    to that target? If multiple positions are found we 
    simply return the first one in the positions list

    target [list of ints] 
    2D or 3D posiition list

    position_list [list of list of ints]
    A list of 2D or 3D positions to be compared against the 
    target

    Return:
    -------------
    tuple

        (int, float)
        [0] - Integer reporting on the index of the minimum position   
        [1] - Actual distance between target and the closes position    

    """

    if len(positions_list) == 0:
        raise LatticeUtilsException("Error in lattice_utils.find_nearest_position() - possition_list is empty")

    
    # cycle through each position and find the single position
    # closes to the COM!
    minimum_distance = 100000000000
    min_idx = -1
    idx = 0    
    for pos in positions_list:
        dist = get_real_distance(target, pos, dimensions)
        if dist < minimum_distance:
            center_position = pos
            minimum_distance = dist
            min_idx = idx

        idx=idx+1

    return (min_idx, minimum_distance)

 

#-----------------------------------------------------------------
#
def get_empty_site(lattice_grid, adjacentTo=None, hardwall=False):
    """
    Function which returns the position of an empty site on the lattice.
    If adjacentTo is populated then the returned site is adjacent to the 
    the position defined by adjacentTo.


    """
    
    dimensions = get_dimensions(lattice_grid)

    # if we passed a position which we wish to find a site adjancent to
    if not (adjacentTo is None):

        # If we look for a site adjacent to a preperscribed position
        position = adjacentTo
        
        # list of empty sites
        empty_list   = []

        # get all the sites adjacent (note PBC correction is done here)  
        if len(dimensions) == 2:
            initial_adjacent_sites = get_adjacent_sites_2D(position[0],position[1], dimensions)
        else:
            initial_adjacent_sites = get_adjacent_sites_3D(position[0], position[1], position[2], dimensions)


        # if hardwall requested, then we only consider sites that don't
        # stradle the boundary
        adjacent_sites = []
        if hardwall:            
            for site in initial_adjacent_sites:                
                site_is_ok=True
                for d in range(0, len(dimensions)):
                    if abs(site[d] - adjacentTo[d]) > 1:
                        site_is_ok=False

                if site_is_ok:
                    adjacent_sites.append(site)
        else:
            adjacent_sites = initial_adjacent_sites

        # find the empty sites
        for site in adjacent_sites:                            
            # if the associated grid element is 0...
            if get_gridvalue(site, lattice_grid) == 0:
                empty_list.append(site)
                                
        # if we didn't find any empty sites at all..
        if len(empty_list) == 0:
            return ([-1,-1,-1],False)

        # return the list of good sites - note we're casting to a list so as
        # we return lists rather than np.arrays
        else:
            return (list(empty_list[random.randint(0,len(empty_list)-1)]), True)

    # if we're literally just looking for an empty site anywhere on the lattice
    else:
        empty = False
        count = 0
        while not empty:
            count=count+1

            if count % 100 == 0:
                IO_utils.status_message("Tried %i times but unable to insert a single point into an empty space - maybe grid is full?\nWill keep trying though, cos I'm a trooper!",'warning')

                                            
            # 2D
            if len(dimensions) == 2:

                # select a random position
                x = random.randint(0, dimensions[0]-1)
                y = random.randint(0, dimensions[1]-1)
                    
                # if the possition is empty celebrate with a beer!
                if get_gridvalue([x,y], lattice_grid) == 0:
                    position = [x,y]
                    empty=True

            # 3D
            if len(dimensions) == 3:
                x = random.randint(0, dimensions[0]-1)
                y = random.randint(0, dimensions[1]-1)
                z = random.randint(0, dimensions[2]-1)

                if get_gridvalue([x,y,z], lattice_grid) == 0:
                    position = [x,y,z]
                    empty=True

    return position



#-----------------------------------------------------------------
#
def insert_chain(chainID, chain_length, lattice_grid, default_start=None, hardwall=False):
    """
    Function that inserts a chain into the passed lattice


    Parameters
    -----------------
    chainID : int
        Unique ID that identifies a specific chain on the lattice

    chain_length : int
        Number of residues in the chain

    lattice_grid

    """
            
    attempt      = 0
    completed    = False
        
    while attempt < CONFIG.CHAIN_INIT_ATTEMPTS and not completed:

        # randomly select starting position
        position_list = []

        if default_start:
            position = default_start
        else:
            position = get_empty_site(lattice_grid)        
        

        # save the starting position because it gives another
        # end for chain expansion if we get ourselves into a 
        # knot!
        start_pos    = position
        head_to_tail = False
        construction_failure = False

        # 
        set_gridvalue(position, chainID, lattice_grid)
        position_list.append(position)

        # if we're looking at particles instead of chains, then we're done!
        if chain_length == 1:
            return position_list            

            
        for i in range(1, chain_length):
            # not -1 because we asign the first lattice site outside
            # the loop

            (position, site_found) = get_empty_site(lattice_grid, adjacentTo=position, hardwall=hardwall)
                
            # if we couldnt find a single empty site adjacent to
            # position
            if not site_found:

                if not head_to_tail:

                    # if we haven't yet tried extending from the other end of
                    # the chain
                    head_to_tail=True
                    position=start_pos
                        
                    # try now from the other end of the chain
                    (position, site_found) = get_empty_site(lattice_grid, adjacentTo=position, hardwall=hardwall)
                        
                    if not site_found:
                        IO_utils.status_message("Chain (ID=%i) construction failed [TRY %i of %i]" %(chainID, attempt+1, CONFIG.CHAIN_INIT_ATTEMPTS), 'warning')
                        attempt = attempt+1
                        construction_failure = True
                        delete_chain_by_ID(chainID, lattice_grid)
                        break

                    # if site was found we update and stick this new position at the front,
                    # then continue on with the head_to_tail flag set to true so all other positions
                    # are added to the head
                    position_list.insert(0,position)
                    set_gridvalue(position, chainID, lattice_grid)
                else:
                    # if we're here we've got to a dead end and know the other end of the chain was 
                    # also a dead end!!!
                    IO_utils.status_message("Chain (ID=%i) construction failed [TRY %i of %i]" %(chainID, attempt+1, CONFIG.CHAIN_INIT_ATTEMPTS),'warning')
                    attempt = attempt+1
                    construction_failure = True
                    delete_chain_by_ID(chainID, lattice_grid)
                    break

            # a site was found!
            else:              

                # if we're in head-to-tail mode add the site to the front of the growing list of positions
                if head_to_tail:
                    position_list.insert(0,position)                    
                # else add the site to the end
                else:
                    position_list.append(position)
                    
                set_gridvalue(position, chainID, lattice_grid)

        # if we're outside of that FOR loop because the chain completed...                
        if not construction_failure == True:
            completed=True

    # if we're outside 
    if not completed:
        raise ChainInsertionFailure
    return position_list



#-----------------------------------------------------------------
#
def place_chain_by_position(positions, lattice_grid, chainID, safe=False):
    """
    Sets the positions defined in the positions vectors to a chain. Note
    this can be an entire chain or a subset of a chain

    """


    if safe:
        for position in positions:
            if not get_gridvalue(position, lattice_grid) == 0.0:
                raise ChainInsertionFailure('Tried to place chain %i at position '%chainID + str(position) + " but found it was occupied by chain [%i]..."%get_gridvalue(position, lattice_grid))
            else:
                set_gridvalue(position, chainID, lattice_grid)
    else:
        for position in positions:
            set_gridvalue(position, chainID, lattice_grid)


    
#-----------------------------------------------------------------
#
def delete_chain_by_ID(chainID, lattice_grid):
    """
    Deletes a chain from the lattice based on the chain's
    ID

    """
    lattice_grid[lattice_grid == chainID] = 0.0



#-----------------------------------------------------------------
#
def delete_chain_by_position(positions, lattice_grid, chainID=None):
    """
    Deletes a chain based on supplied position. Can be an entire chain
    or simply a portion of a chain

    """

    # safe (checks before deletion
    if chainID is not None:
        for position in positions:
            if not chainID == get_gridvalue(position, lattice_grid):
                raise ChainDeletionFailure('Tried to delete chain %i at position'%chainID + str(position) + " but this position was not occupied by the expected chain")
            set_gridvalue(position, 0, lattice_grid)

    # fast - no checks...
    else:
        for position in positions:
            set_gridvalue(position, 0, lattice_grid)
            


#-----------------------------------------------------------------
#                        
def get_gridvalue(position, lattice_grid):
    """
    Returns the value on the lattice grid
    associated with the position defined by
    the 2/3 place tuple

    """
    
    dimensions = get_dimensions(lattice_grid)
    
    if len(dimensions) == 2:
        return lattice_grid[position[0]][position[1]]

    if len(dimensions) == 3:
        return (lattice_grid[position[0]][position[1]][position[2]])


#-----------------------------------------------------------------
#                        
def get_gridvalue_2D(position, lattice_grid):
    """
    Returns the value on the lattice grid
    associated with the position defined by
    the 2/3 place tuple

    """
    return lattice_grid[position[0]][position[1]]



#-----------------------------------------------------------------
#                        
def get_gridvalue_3D(position, lattice_grid):
    """
    Returns the value on the lattice grid
    associated with the position defined by
    the 2/3 place tuple

    """
    return hyperloop.get_gridvalue_3D(lattice_grid, position[0], position[1], position[2])



#-----------------------------------------------------------------
#
def set_gridvalue(position, value, lattice_grid):
    """
    Sets the position defined at lattice site $position on 
    $lattice_grid to $value. This CHANGES the $lattice_grid
    object (which is assumed to be a numpy 2D or 3D array)
    and returns it.

    """

    if len(position) == 2:            
        lattice_grid[position[0]][position[1]] = value

    if len(position) == 3:
        lattice_grid[position[0]][position[1]][position[2]] = value

    return lattice_grid


#-----------------------------------------------------------------
#
def build_envelope_pairs(positions, dimensions, hardwall=False):
    """
    Expects a LIST of positions. Returns a unique unordered
    list of tuples, where each tuple is a pair of positions.

    The complete set of these positions represents the non-redundant 
    set of positions that make contact with the positions in the past
    $positions variable.

    dimensions is the dimensions of the lattice

    hardwall is a boolean which determines if we allow a pair of 
    positions to straddle the boundary (in a periodic manner) or not
    
    """
    
    # now remove any duplicat pairs in there
    if len(dimensions) == 2:
        # if 2D
        
        short_range_list = []

        if hardwall:
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops_hardwall.extract_SR_pairs_from_position_2D_hardwall(np.array(positions[i]), dimensions[0], dimensions[1]))
                
        else:
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops.extract_SR_pairs_from_position_2D(np.array(positions[i]), dimensions[0], dimensions[1]))

        envelope_pairs = np.concatenate(short_range_list)
        num_pairs = len(envelope_pairs)

        reshaped = np.reshape(envelope_pairs, (num_pairs, 4))

        b = np.ascontiguousarray(reshaped).view(np.dtype((np.void, reshaped.dtype.itemsize * reshaped.shape[1])))
        _, idx = np.unique(b, return_index=True)

        duplicate_free = reshaped[idx]

        return np.reshape(duplicate_free, (len(duplicate_free), 2,2))
    else:

        short_range_list = []

        if hardwall:
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops_hardwall.extract_SR_pairs_from_position_3D_hardwall(np.array(positions[i]), dimensions[0], dimensions[1], dimensions[2]))
        else:
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops.extract_SR_pairs_from_position_3D(np.array(positions[i]), dimensions[0], dimensions[1], dimensions[2]))

        envelope_pairs = np.concatenate(short_range_list)
        num_pairs = len(envelope_pairs)

        reshaped = np.reshape(envelope_pairs, (num_pairs, 6))

        b = np.ascontiguousarray(reshaped).view(np.dtype((np.void, reshaped.dtype.itemsize * reshaped.shape[1])))
        _, idx = np.unique(b, return_index=True)

        duplicate_free = reshaped[idx]

        return np.reshape(duplicate_free, (len(duplicate_free), 2,3))

#-----------------------------------------------------------------
#
#@profile
def build_all_envelope_pairs(positions, LR_binary_array, type_lattice, dimensions, hardwall=False):
    """
    Expects a LIST of positions and a numpy array of positions which engage in
    long-range interactions (or not) - 0 if not and 1 if yes.

    Returns a list of tuples, where each tuple is a pair of positions.

    The complete set of these positions represents the non-redudant set of long
    -range and short-range pairwise interactions associated with the positions
    defined in position
    """


    # I don't know why I created a mode selector, but now I'm too scared to 
    # remove it...
    mode = 1

    if mode == 1:

        super_long_range_list = []
        long_range_list  = []
        short_range_list = []
        
        ####>>>> 2D
        if len(dimensions) == 2:
            for i in range(0, len(positions)):

                # get enveloping pairs
                if hardwall:
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_2D_hardwall(np.array(positions[i], dtype=int), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1])
                else:
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops.extract_SR_and_LR_pairs_from_position_2D(np.array(positions[i], dtype=int), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1])
                
                short_range_list.append(SR_tmp)
                
                if len(LR_tmp) > 0:
                    long_range_list.append(LR_tmp)

                if len(SLR_tmp) > 0:
                    super_long_range_list.append(SLR_tmp)

        ####>>>> 3D
        else:
            for i in range(0, len(positions)):

                if hardwall:
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_3D_hardwall(np.array(positions[i], dtype=int), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1], dimensions[2])
                else:
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops.extract_SR_and_LR_pairs_from_position_3D(np.array(positions[i], dtype=int), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1], dimensions[2])
                                
                short_range_list.append(SR_tmp)
                
                if len(LR_tmp) > 0:
                    long_range_list.append(LR_tmp)

                if len(SLR_tmp) > 0:
                    super_long_range_list.append(SLR_tmp)
                
    short_range_pairs = np.concatenate(short_range_list)

    # note we have to check LR and SLR pairs seperately !
    if len(long_range_list) > 0:
        long_range_pairs = np.concatenate(long_range_list)        
    else:
        long_range_pairs = np.array([], dtype=int)

    if len(super_long_range_list) > 0:
        super_long_range_pairs = np.concatenate(super_long_range_list)        
    else:
        super_long_range_pairs = np.array([], dtype=int)

    # compute the number of each type of pair
    num_pairs_SR = len(short_range_pairs)
    num_pairs_LR = len(long_range_pairs)
    num_pairs_SLR = len(super_long_range_pairs)

    # now remove any duplicat pairs in there
    if len(dimensions) == 2:
        # if 2D
        
        reshaped_SR = np.reshape(short_range_pairs, (num_pairs_SR, 4))
        reshaped_LR = np.reshape(long_range_pairs, (num_pairs_LR, 4))
        reshaped_SLR = np.reshape(super_long_range_pairs, (num_pairs_SLR, 4))

        # short range witchcraft
        b_SR = np.ascontiguousarray(reshaped_SR).view(np.dtype((np.void, reshaped_SR.dtype.itemsize * reshaped_SR.shape[1])))
        _, idx_SR = np.unique(b_SR, return_index=True)
        duplicate_free_SR = reshaped_SR[idx_SR]

        # long range witchcraft
        b_LR = np.ascontiguousarray(reshaped_LR).view(np.dtype((np.void, reshaped_LR.dtype.itemsize * reshaped_LR.shape[1])))
        _, idx_LR = np.unique(b_LR, return_index=True)
        duplicate_free_LR = reshaped_LR[idx_LR]

        # super long range witchcraft
        b_SLR = np.ascontiguousarray(reshaped_SLR).view(np.dtype((np.void, reshaped_SLR.dtype.itemsize * reshaped_SLR.shape[1])))
        _, idx_SLR = np.unique(b_SLR, return_index=True)
        duplicate_free_SLR = reshaped_SLR[idx_SLR]


        return (np.reshape(duplicate_free_SR, (len(duplicate_free_SR), 2,2)), np.reshape(duplicate_free_LR, (len(duplicate_free_LR), 2,2)), np.reshape(duplicate_free_SLR, (len(duplicate_free_SLR), 2,2)))

    else:
        # if 3D
                 
        reshaped_SR  = np.reshape(short_range_pairs, (num_pairs_SR, 6))
        reshaped_LR  = np.reshape(long_range_pairs,  (num_pairs_LR, 6))
        reshaped_SLR = np.reshape(super_long_range_pairs,  (num_pairs_SLR, 6))

        # short range witchcraft
        b_SR = np.ascontiguousarray(reshaped_SR).view(np.dtype((np.void, reshaped_SR.dtype.itemsize * reshaped_SR.shape[1])))
        _, idx_SR = np.unique(b_SR, return_index=True)
        duplicate_free_SR = reshaped_SR[idx_SR]

        # long range witchcraft
        b_LR = np.ascontiguousarray(reshaped_LR).view(np.dtype((np.void, reshaped_LR.dtype.itemsize * reshaped_LR.shape[1])))
        _, idx_LR = np.unique(b_LR, return_index=True)
        duplicate_free_LR = reshaped_LR[idx_LR]

        # super long range witchcraft
        b_SLR = np.ascontiguousarray(reshaped_SLR).view(np.dtype((np.void, reshaped_SLR.dtype.itemsize * reshaped_SLR.shape[1])))
        _, idx_SLR = np.unique(b_SLR, return_index=True)
        duplicate_free_SLR = reshaped_SLR[idx_SLR]

        return (np.reshape(duplicate_free_SR, (len(duplicate_free_SR), 2,3)), np.reshape(duplicate_free_LR, (len(duplicate_free_LR), 2,3)), np.reshape(duplicate_free_SLR, (len(duplicate_free_SLR), 2,3)))



#-----------------------------------------------------------------
#    
def get_all_chains_in_connected_component(chainID, lattice_grid, chainDict, threshold=None, useChains=True, hardwall=False):
    """
    Function which given a chainID, a dictionary of chain-to-position mappings, and a lattice grid
    will return the set of chains in the connected component containing chainID. Note a connected 
    component is a *heterotypic* structure - i.e. we are looking for a connected component made up
    of *any* chains, not a single type of chain.

    Arguments:
    
    chainID [int]
    The chainID of the chain we initially are asking about

    lattice_grid [2D or 3D np.array]
    Standard lattice grid

    chainDict [dictionary mapping chainIDs to a list of positions or to a chain object]
    Dictionary containing a mapping of each chainID to either a list of positions associated
    with that chain, or the Chain object associated with that chainID

    useChains [Bool]
    Boolean flag which defines if the chainDict is a true dictionary mapping chainID
    to a set of positions, or in fact a dictionary of Chain objects (which contain
    positions which must be accessed using the .get_ordered_positions()). This isn't
    so much a feature as the fact that we want this function to be able to accept
    two different types of chain information (dictionary of lists of positions or
    dictionary of chain objects)
    
    Returns:
    A list of chainIDs associated with the chains in the connected component 
    which contans the chain defined by $chainID
    

    """

    chains     = set([])
    new_chains = set([])
    dimensions = get_dimensions(lattice_grid)
    
    chains.add(chainID)
    new_chains.add(chainID)
                    
    if useChains:        
        positions = chainDict[chainID].get_ordered_positions()
    else:
        positions = chainDict[chainID]
        
    # loop until we break with a return statement
    while True:

        # get all the envelope pairs assoiated with the list of positions
        envelope_pairs = build_envelope_pairs(positions, dimensions, hardwall=hardwall)

        # for each position associated with each pair figure out what chain
        # it comes from
        for pair in envelope_pairs:        
            new_chains.add(get_gridvalue(pair[0], lattice_grid))
            new_chains.add(get_gridvalue(pair[1], lattice_grid))

        # having done that for every pair remove the 'solvent' chains        
        try:
            new_chains.remove(0)
        except KeyError:
            # in the case where our grid is at 100% volume fraction of no solvent 
            # don't freak out that we can't remove solvent because no solvent chains
            # were added (e.g if a chain is entirely encapsulated by other chains)
            pass
            
        # if the set of chains hasn't changed then we're done
        if len(new_chains) == len(chains):
            return list(chains)
        # found at least one new chain
        else:

            # note this expression is doing a set operation and generating
            # the set of chains found in $new_chains which was not found in
            # the $chains set
            newly_found_chains = new_chains - chains

            positions = []
            
            # for all the new chains create a new list of positions            
            for chain in newly_found_chains:

                if useChains:
                    positions.extend(chainDict[chain].get_ordered_positions())
                else:
                    positions.extend(chainDict[chain])

                chains.add(chain)
                
            # if we defined a threshold and we're above it...
            if threshold is not None and len(chains) > threshold:
                raise ClusterSizeThresholdException



#-----------------------------------------------------------------
#    
def get_all_chains_in_long_range_cluster(chainID, latticeObject, hardwall=False):

    
    lattice_grid = latticeObject.grid
    type_grid     = latticeObject.type_grid
    chainDict    = latticeObject.chains

    chains     = set([])
    new_chains = set([])
    dimensions = get_dimensions(lattice_grid)
    
    
    chains.add(chainID)
    new_chains.add(chainID)
                    

    positions      = chainDict[chainID].get_ordered_positions()
    LR_binary_array = chainDict[chainID].get_LR_binary_array()
        
    # loop until we break with a return statement

    while True:

        # get all the envelope pairs assoiated with the list of positions
        (SR_pairs, LR_pairs, SLR_pairs) = build_all_envelope_pairs(positions, LR_binary_array, type_grid, dimensions, hardwall)

        envelope_pairs = np.concatenate((SR_pairs, LR_pairs))    
                
        # for each position associated with each pair figure out what chain
        # it comes from
        for pair in envelope_pairs:        
            new_chains.add(get_gridvalue(pair[0], lattice_grid))
            new_chains.add(get_gridvalue(pair[1], lattice_grid))

        # having done that for every pair remove the 'solvent' chains        
        try:
            new_chains.remove(0)
        except KeyError:
            # in the case where our grid is at 100% volume fraction of no solvent 
            # don't freak out that we can't remove solvent because no solvent chains
            # were added (e.g if a chain is entirely encapsulated by other chains)
            pass
            
        # if the set of chains hasn't changed then we're done
        if len(new_chains) == len(chains):
            return list(chains)
        # found at least one new chain
        else:

            # note this expression is doing a set operation and generating
            # the set of chains found in $new_chains which was not found in
            # the $chains set
            newly_found_chains = new_chains - chains

            positions = []
            LR_binary_array = []
            
            # for all the new chains create a new list of positions            
            for chain in newly_found_chains:
                positions.extend(chainDict[chain].get_ordered_positions())
                LR_binary_array.extend(chainDict[chain].get_LR_binary_array())
                chains.add(chain)

                

#-----------------------------------------------------------------
#    
def center_of_mass_from_positions(positions, dimensions, on_lattice=True):
    """
    Return the center of mass from the list of positions.
    Assumes all positions have the same mass!

    on_lattice can be set to True if you want a lattice-based COM
    or set to False if you want the true off-lattice Euclidean COM

    COM is calculated by implementing the agorithm developed by Bai
    and Breen [1] extended to 3D, which means it determines the
    correct center of mass in a periodic box.

    [1] Bai, L., & Breen, D. (2008). Calculating Center of Mass in an 
    Unbounded 2D Environment. Journal of Graphics, GPU, and Game 
    Tools, 13(4), 53 - 60.

    """

    n_dim = len(dimensions)
    xmax = dimensions[0]
    ymax = dimensions[1]

    x_polar_1 = 0
    x_polar_2 = 0

    y_polar_1 = 0
    y_polar_2 = 0

    if n_dim == 2:           

        for position in positions:
            x_polar_1  = np.cos((position[0]/float(xmax))*2*np.pi) + x_polar_1
            x_polar_2  = np.sin((position[0]/float(xmax))*2*np.pi) + x_polar_2

            y_polar_1  = np.cos((position[1]/float(ymax))*2*np.pi) + y_polar_1
            y_polar_2  = np.sin((position[1]/float(ymax))*2*np.pi) + y_polar_2
            
        mean_x_polar_1 = x_polar_1/len(positions)
        mean_x_polar_2 = x_polar_2/len(positions)

        mean_y_polar_1 = y_polar_1/len(positions)
        mean_y_polar_2 = y_polar_2/len(positions)

        x_real = xmax*(np.arctan2(-mean_x_polar_2,-mean_x_polar_1)+np.pi)/(2*np.pi)
        y_real = ymax*(np.arctan2(-mean_y_polar_2,-mean_y_polar_1)+np.pi)/(2*np.pi)
            
        if on_lattice:
            x = int(round(x_real))
            y = int(round(y_real))
        else:
            x = x_real
            y = y_real

        return pbc_convert([x, y], dimensions)

    if n_dim == 3:           
        zmax = dimensions[2]
        z_polar_1 = 0
        z_polar_2 = 0
        for position in positions:
            x_polar_1  = np.cos((position[0]/float(xmax))*2*np.pi) + x_polar_1
            x_polar_2  = np.sin((position[0]/float(xmax))*2*np.pi) + x_polar_2

            y_polar_1  = np.cos((position[1]/float(ymax))*2*np.pi) + y_polar_1
            y_polar_2  = np.sin((position[1]/float(ymax))*2*np.pi) + y_polar_2

            z_polar_1  = np.cos((position[2]/float(zmax))*2*np.pi) + z_polar_1
            z_polar_2  = np.sin((position[2]/float(zmax))*2*np.pi) + z_polar_2
            
        # Deterine mean values
        mean_x_polar_1 = x_polar_1/len(positions)
        mean_x_polar_2 = x_polar_2/len(positions)

        mean_y_polar_1 = y_polar_1/len(positions)
        mean_y_polar_2 = y_polar_2/len(positions)

        mean_z_polar_1 = z_polar_1/len(positions)
        mean_z_polar_2 = z_polar_2/len(positions)

        x_real = xmax*(np.arctan2(-mean_x_polar_2,-mean_x_polar_1)+np.pi)/(2*np.pi)
        y_real = ymax*(np.arctan2(-mean_y_polar_2,-mean_y_polar_1)+np.pi)/(2*np.pi)
        z_real = zmax*(np.arctan2(-mean_z_polar_2,-mean_z_polar_1)+np.pi)/(2*np.pi)
            
        if on_lattice:
            x = int(round(x_real))
            y = int(round(y_real))
            z = int(round(z_real))
        else:
            x = x_real
            y = y_real
            z = z_real

        return pbc_convert([x, y, z], dimensions)

#######################################################################################
##                                                                                   ##
##                            Residue functions are here                             ##
##                                                                                   ##
#######################################################################################
#
# Note the insert and delete residue functions are basically just wrappers around set_gridvalue
# except they offer some sanity checking, which is probably a good idea (especially for moves)
# though less crucial when developing lower level routines.
#


#-----------------------------------------------------------------
#
def delete_residue(position, lattice, chainID=None):

    if chainID is not None:
        ## Safe version

        todel = get_gridvalue(position, lattice)
        if not todel == chainID:
            raise ResidueAugmentException('Trying to delete a residue at position' + str(position) + ' - expected chainID %i, but got chainID...' %(chainID, todel))
        else:
            set_gridvalue(position, 0.0, lattice)

    else:
        ## No checks version...
        set_gridvalue(position, 0.0, lattice)



#-----------------------------------------------------------------
#
def insert_residue(position, lattice, chainID, safe=True):

    if safe:
        insert_location = get_gridvalue(position, lattice)
        if not insert_location == 0.0:
            raise ResidueAugmentException('Trying to insert a residue for chain at position'%chainID + str(position) + ' - site was occupied by residue from chain %i !!'%insert_location)
        else:
            set_gridvalue(position, chainID, lattice)
    else:
        set_gridvalue(position, chainID, lattice)

#######################################################################################
##                                                                                   ##
##                              Rotation operations                                  ##
##                                                                                   ##
#######################################################################################

def run_rotation(positions, rotation_matrix):
    """
    Perform single point rotation 

    """
    rotated_positions = []    
    for position in positions:        
        rotated_positions.append(np.dot(rotation_matrix, position))

    return rotated_positions



#-----------------------------------------------------------------
#
def rotate_positions_3D(positions, dimension, degrees):    
    """
    Functions to carry out cardinal position rotation around the origin.

    The CARDINAL_ROTATION_3D matrix is assigned in CONFIG, affording
    extremely fast rotation. 

    """
    

    if dimension =='x':
        if degrees == 90:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[0][0])
        if degrees == 180:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[1][0])
        if degrees == 270:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[2][0])

    if dimension =='y':
        if degrees == 90:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[0][1])
        if degrees == 180:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[1][1])
        if degrees == 270:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[2][1])

    if dimension =='z':
        if degrees == 90:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[0][2])
        if degrees == 180:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[1][2])
        if degrees == 270:
            return run_rotation(positions, CONFIG.CARDINAL_ROTATION_3D[2][2])

    # If we get here passed a non cardinal dimension or degrees
    raise RotationException('Trying to rotate axis %s around %s degrees - INVALID' % (str(dimension), str(degrees)))



#-----------------------------------------------------------------
#    
def rotate_positions_2D(positions, degrees):    
    """
    Functions to carry out 2D cardinal position rotation around the origin.

    The CARDINAL_ROTATION_2D matrix is assigned in CONFIG, affording
    extremely fast rotation. 

    """

    if degrees == 90:        
        return run_rotation(positions, CONFIG.CARDINAL_ROTATION_2D[0])

    if degrees == 180:
        return run_rotation(positions, CONFIG.CARDINAL_ROTATION_2D[1])

    if degrees == 270:
        return run_rotation(positions, CONFIG.CARDINAL_ROTATION_2D[2])

    # If we get here passed a non cardinal dimension or degrees
    raise RotationException('Trying to positions around %s degrees - INVALID' % (str(degrees)))



#######################################################################################
##                                                                                   ##
##                              I/O functions are here                               ##
##                                                                                   ##
#######################################################################################

#-----------------------------------------------------------------
#
def open_pdb_file(dimensions, spacing, filename="lattice.pdb"):
    """
    Function that initializes a PDB file to be written to.

    Parameters
    -------------
    dimensions : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects the lattice dimensions.

    spacing : float
        Lattice-to-realspace spacing in angstroms. 

    filename : str
        Filename to write to

    """

    pdb_utils.initialize_pdb_file(dimensions, spacing, filename)
    

#-----------------------------------------------------------------
#
def write_lattice_to_pdb(latticeObject, spacing, filename='lattice.pdb', write_connect=False):
    """
    Wrapper function that dumps the current Lattice object to a PDB file

    Parameters
    -------------
    latticeObject : lattice.Lattice
        Current lattice object

    spacing : float
        Lattice-to-realspace spacing in angstroms. 

    filename : str
        Filename to write to

    """
    pdb_utils.build_pdb_file(latticeObject, spacing, filename, write_connect=write_connect)



#-----------------------------------------------------------------
#
def finish_pdb_file(filename):
    """
    Function that finalizes a PDB by adding terminating information.

    Parameters
    -----------------
    filename : str
        Filename to be finalized

    Returns
    ----------
    None
        No return but the file associated with filename is finalized as a
        PDB file.
    
    

    
    """
    pdb_utils.finalize_pdb_file(filename)


#-----------------------------------------------------------------
#
def start_xtc_file(lattice, spacing, pdb_filename='START.pdb', xtc_filename='traj.xtc'):
    """
    Function that initializes a new .xtc file. This deletes an existing XTC file 
    of the same name to avoid any issues.

    Parameters
    ------------
    lattice : lattice.Lattice 
        Current Lattice object

    pdb_filename : str
        New XTC files need a corresponding PDB file. This defines the name of that
        PDB file.

    xtc_filename : str
        New XTC files need a corresponding PDB file. This defines the name of that
        PDB file.


    Returns
    ------------
    None
        No return value, but a newly initialized XTC file is generated

    """
    # delete the xtc file if it exists already
    try:        
        os.remove(xtc_filename)

        # if the file doesn't exit this throws an OSError that we deal with
        # here and so its never an issue! 
        IO_utils.status_message("Deleted existing XTC file [%s]"%xtc_filename,'startup')

    except OSError:
        pass

    # first build the PDB file
    open_pdb_file(lattice.dimensions, spacing, filename=pdb_filename)
    write_lattice_to_pdb(lattice, spacing, filename=pdb_filename, write_connect=True)
    finish_pdb_file(pdb_filename)

    # next read the PDBFILE, and save as an xtcfile    
    traj = md.load(pdb_filename)
    traj.save_xtc(xtc_filename)
    


#-----------------------------------------------------------------
#
def append_to_xtc_file(lattice, spacing, xtc_filename='traj.xtc'):
    """
    Low level function that adds a current lattice to and existing XTC file

    Parameters
    -----------
    lattice : lattice.Lattice
        A lattice object

    spacing : float
        Lattice-to-realspace spacing in angstroms. 

    xtc_filename : str
        Filename to read from and extend

    Returns
    -----------
    None
        No return by the existing XTC file is extended by one frame

    """

    # first build the PDB file    
    open_pdb_file(lattice.dimensions, spacing, filename='frame.pdb')
    write_lattice_to_pdb(lattice, spacing, filename='frame.pdb')
    finish_pdb_file('frame.pdb')

    xtc_traj = md.load(xtc_filename, top='frame.pdb')
    
    pdb_frame = md.load('frame.pdb')

    new = xtc_traj.join(pdb_frame)

    new.save(xtc_filename)





    

#######################################################################################
##                                                                                   ##
##                           SANITY CHECKING FUNCTIONS                               ##
##                                                                                   ##
#######################################################################################

    
def check_chain_connectivity(chainID, chain_positions, dimensions):
    """
    Debugging function which ensures that a set of positions
    correspond to a valid, connected chain. Useful for
    debugging new moves, though not designed for performance
    during real simulations 

    """

    num_positions = len(chain_positions)
    num_dims      = len(chain_positions[0])

    for position in range(0, num_positions-1):
        current_position = chain_positions[position]
        next_position    = chain_positions[position+1]
        

        for i in range(0,num_dims):

            # if diff between two positions is greater than 1 site
            if abs(current_position[i] - next_position[i]) > 1:
                print("(chain %i pos %i) %s---%s" %(chainID, position, current_position, next_position))
                # maybe a PBC issue... correct and try again
                if current_position[i] > next_position[i]:
                    PBC_increased_next = next_position[i] + dimensions[i]

                    if abs(current_position[i] - PBC_increased_next) > 1:
                        raise ChainConnectivityError('Chain %i appears to not be correctly connected at position %i \n %s' % (chainID, position, chain_positions))                        
                        
                else:
                    PBC_increased_current = current_position[i] + dimensions[i]

                    if abs(PBC_increased_current - next_position[i]) > 1:
                        raise ChainConnectivityError('Chain %i appears to not be correctly connected at position %i \n %s' % (chainID, position, chain_positions))                        
                        


                    # test again                                                                    
    print("CONNECTIVITY FINE")


def check_all_chain_connectivity(list_of_chain_objects, dimensions):
    """
    Debugging function that takes a LATTICE.chains and LATTICE.dimensions
    pair from the simulation object to check the chain connectivity over
    all chains in the simulation.


    """
    
    for chainID in list_of_chain_objects:
        check_chain_connectivity(chainID, list_of_chain_objects[chainID].get_ordered_positions(), dimensions)




        
