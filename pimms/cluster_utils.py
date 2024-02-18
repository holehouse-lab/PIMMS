## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................



import numpy as np
import copy

from . import hyperloop
from . import lattice_utils
from . import lattice_analysis_utils
from .latticeExceptions import UnfinishedCodeException
from .pdb_utils import write_positions_to_file

##
## NOTE that as of right now these functions are not used, but we'll keep them around in case of future
## developments...

def build_interface_envelope_pairs(positions, dimensions, grid):
    """
    Function which determines the set of non-redundant pairs that describe the interface of the positions in
    the list $positions with the solvent. The KEY assumption here is that positions FULLY describes all the positions
    of occupied lattice sites in a connected component. This method WILL BREAK if that assumption is wrong.

    This is the reason this method is off in its own set of utilities - if we're not working with a complete (redundant)
    set of positions that describe a single connected component (cluster) everything is gonna go wrong...

    But assuming that's true, this method is faster than a shark on cocaine.

    """

    # 2D case
    if len(dimensions) == 2:
        
        # initialze the first set of values in the envelope pairs array
        
        # the following code is a bit testy, but we want to start with an appropriatly sized vector and then grow it by vstacking
        initial_position = 0
        searching = True
        while searching:
            position = positions[initial_position]
            envelope_pairs = np.array(hyperloop.get_unique_interface_pairs_2D(position[0],position[1], dimensions[0], dimensions[1], grid))
            initial_position=initial_position+1

            # if we found at least one pair
            if envelope_pairs.shape[0] > 0:
                searching = False

                
        
        # then just add to them
        for position in positions[initial_position:]:
            tmp = np.array(hyperloop.get_unique_interface_pairs_2D(position[0],position[1], dimensions[0], dimensions[1],grid))

            if tmp.shape[0] == 0:
                # if no positions founds skip
                continue

            envelope_pairs = np.vstack((envelope_pairs, tmp))

    # 3D case
    if len(dimensions) == 3:

        # the following code is a bit testy, but we want to start with an appropriatly sized vector and then grow it by vstacking
        initial_position = 0
        searching = True

        while searching:
            position = positions[initial_position]
            envelope_pairs = np.array(hyperloop.get_unique_interface_pairs_3D(position[0], position[1], position[2], dimensions[0], dimensions[1], dimensions[2], grid))
            initial_position=initial_position+1

            # if we found at least one pair
            if envelope_pairs.shape[0] > 0:
                searching = False
        
        for position in positions[1:]:
            tmp = np.array(hyperloop.get_unique_interface_pairs_3D(position[0], position[1], position[2], dimensions[0], dimensions[1], dimensions[2], grid))

            if tmp.shape[0] == 0:
                # if no positions founds skip
                continue
            envelope_pairs = np.vstack((envelope_pairs, tmp))

    return envelope_pairs

    
def build_interface_envelope_pairs_safe_and_slow(positions, dimensions):    
    """
    This provides an alternative implementation of the function above, but instead users the softere assumption
    that interface pairs must be those which are pairs and it's not true that both members of the pair come 
    from within the list of positions supplied.
    
    However, again it's crucial that positions be an exhaustative list of the positions making up a connected 
    component, although this time we don't use any lattice information.

    """
    envelope_pairs = build_envelope_pairs(positions, dimensions)
    
    # for each pair ask if both members of the pair come from the list of positions    
    interface_pairs = []

    np_positions = np.array(positions)

    for pair in envelope_pairs:
        if (~(np_positions-pair[0]).any(axis=1)).any() and (~(np_positions-pair[1]).any(axis=1)).any():
            pass
        else:
            interface_pairs.append(pair)

        if numpy_utils.position_in_list(pair[0], positions) and numpy_utils.position_in_list(pair[1], positions):
            # if both pairs of interactions come from inside the positions passed
            pass

    num_pairs = len(interface_pairs)

    #print "Saved vs. looped - RATIO: [%s]" %(str(float(len(envelope_pairs)-num_pairs)/len(envelope_pairs)))
    if num_pairs == 0:
        interface_pairs = np.array([],dtype=int)

    if len(dimensions) == 2:
        reshaped = np.reshape(interface_pairs, (num_pairs, 2,2))
    else:
        reshaped = np.reshape(interface_pairs, (num_pairs, 2,3))

    return reshaped


#-----------------------------------------------------------------
#    
def convert_positions_to_single_image_snakesearch(original_positions, dimensions, space_threshold=1):
    """
    This is a deterministic and somewhat expensive algorithm that allows you to reconstruct a list of positions in
    such a manner that they're all in the same periodic image. This relies on the positions being a 'connected cluster'
    where the distance between those connections can span up to space_threshold distance. Note that this algorithm will
    probably stop working as space threshold gets big. Also where there are clusters that span the entire system there's
    a chance that the snakesearch algorithm is gonna give you some weird extended structures depending on how it searches.
    This could probably be fixed by performing the subsearches in a radially distributed way around the COM but it's also
    probably worth pointing out that when you have system spanning clusters you probably already have finite size
    artefacts, so other than saying 'yes, we've crossed the percolation threshold' computing properties like Rg,
    volume etc is probably not that meaningful.
    
    snakesearch algorithm
    The snakesearch algorithm lets you construct a single image convention structure (i.e a structure which exists in a 
    single periodic image) by taking advantage of the fact that a cluster must be connected. The algorithm is an iterative
    implementation of what feels like a A* recursive search algorithm, which also then corrects as we cross boundary
    walls to reconstruct a single structure. 

    Starting with the residue closest to the center of mass, we find all the residues in 'contact' with that COM residue
    (where contact is defined as being equal to or 1 < than the space threshold in each dimension on the lattice). For
    each of those positions the same operation is performed, meaning we slowly build out the set of components, where 
    we always have a 'reference' point in some periodic image which we then ensure all connected regions are also in.

    The first bead is always the center of mass of the clusters.

    Parameters:
    -----------

    original_positions: list 
        This is a a list of lists, where each sublists contains the position of a bead in the system. 
        
    dimensions: list of integers

        This list defines the dimensions of the lattice, i.e. the number of lattice sites in each dimension.

    space_threshold: integer
        This value defines the maximum distance between two beads in each dimension for them to be 
        considered 'in contact'.

    Returns:
    ----------
    
    single_image_positions: list
        This is a list of lists, where each sublist contains the position of a bead in the system, but all
        beads are in the same periodic image.

    """

    if len(original_positions) == 0:
        return  []

    # set the number of dimensions
    n_dim = len(dimensions)

    # firstly identify the center of mass of all the points - the algorithm implemented for
    # COM is PBC agnostic and will return the correct on-lattice COM position
    COM = lattice_utils.center_of_mass_from_positions(original_positions, dimensions)

    # find the bead closest to the COM
    center_position = original_positions[lattice_utils.find_nearest_position(COM, original_positions, dimensions)[0]]
    
    # initiallize the unsearch list to contain initially only the starting bead 
    unsearched = []
    unsearched.append(copy.copy(center_position))

    unsearched_pbc = []
    unsearched_pbc.append(copy.copy(center_position))

    PS_2_SIS = {}
    PS_2_SIS[tuple(unsearched[0])] = unsearched[0]
    unfound = copy.deepcopy(original_positions)
    unfound.remove(unsearched[0])
    
    while len(unfound) > 0:

        # On each iteration calculate the distance between each position and the COM and then use the position closes to 
        # the COM as the next residue to locally explore. NOTE this is a very bad implementation where we check the same 
        # positions again and again - PLEASE optimize this code!
        num_unsearched = len(unsearched)        
        dist = lattice_analysis_utils.get_inter_position_distances(np.transpose(np.matrix.repeat(np.array(COM), num_unsearched).reshape(n_dim,num_unsearched)) , np.array(unsearched_pbc), dimensions, pbc_correction=False) 

        # extract the next position in PBC space
        target_pbc = unsearched_pbc.pop(np.argmin(dist))
        
        # also get the position in SIC space and delete from the 
        # unsearched list
        target = PS_2_SIS[tuple(target_pbc)]
        unsearched.remove(target)

        # given the target and the list of 'unfound' positions we now identify
        # the set of positions which are immedietly local to the target (i.e. positions
        # which we KNOW are close so any apparent long range interactions are PBC effects
        # and can be corrected because we have the target in the SIC reference frame)
        [in_contact, in_contact_SI] = find_local(target, unfound, dimensions, space_threshold)

        # update the mapping dictionary which maps periodic space positions to single image
        # space positions
        number_found = len(in_contact)
        for idx in range(0, number_found):
            PS_2_SIS[tuple(in_contact[idx])] = in_contact_SI[idx]

        # grow the set of positions we're yet to search around
        unsearched.extend(in_contact_SI)
        unsearched_pbc.extend(in_contact)

        # remove the positions we found from the unfound list
        for pos in in_contact:
            unfound.remove(pos)


    # finally we loop through each of the dimensions of all the positions of the data and
    # if there are any positions < 0 we correct the positions so we're working in a space
    # where every position is positive
    SIS =  np.array(list(PS_2_SIS.values()))    
    for dim in range(0, n_dim):
        dim_positions = np.transpose(SIS)[dim]

        min_dim = min(dim_positions)

        if min_dim < 0:
            np.transpose(SIS)[dim] = dim_positions+abs(min_dim)


    #print "AFTER MINIM CORRECTION"
    #print SIS
    return SIS
        

#-----------------------------------------------------------------
#    
def find_local(original_target, list_of_positions, dimensions, space_threshold):
    """
    Function developed to work with the snakesearch algorithm. Given an original_target
    and a list of positions, find_local identifies the set of positions in the 
    list of positions list that are within space_threshold of the original_target.

    This in itself is not a particularly interesting idea, BUT what's nice (?) is
    that it returns not only the list of positions in PBC space, but additionally
    a list of positions in the same periodic image as the original_target.
    
    Note this could mean positions are negative e.g. if our original target was at
    [0,0,0] in a 10x10x10 box, and there was a residue at [9,9,9] this would be 
    translated to [-1,-1,-1] - i.e. the same image relative to the original target.


    """

    n_dim = len(dimensions)
    target = copy.deepcopy(original_target)

    # first convert target into PBC space - note this will happily deal with
    # a single image convention position in an arbitrarily far image from 
    # the PBC minimum image. NOTE we record the target offset because for
    # each dimensions we want to be able to correct the positions identified 
    # locally back to the single image the original_target came from
    target_offset = [0]*n_dim
    for dim in range(0, n_dim):
        while target[dim] < 0:
            target[dim] = target[dim] + dimensions[dim]
            target_offset[dim] = target_offset[dim] - dimensions[dim]
        while target[dim] > dimensions[dim]-1:
            target[dim] = target[dim] - dimensions[dim]
            target_offset[dim] = target_offset[dim] + dimensions[dim]

    # having now moved the target into minimum image space we find all the positions
    # in the list_of_positiosn variable which are within space_threshold of the
    # target position
    in_contact = []
    for query in list_of_positions:

        closecheck = 0

        for dim in range(0, n_dim):
            d = target[dim] - query[dim]
            
            # if we're close and don't cross
            # a boundary

            if abs(d) <= space_threshold:
                closecheck=closecheck+1
            else:                                            
                # what about if we're close accross
                # a boundary?

                if d < -(0.5*dimensions[dim]):
                    d = d + dimensions[dim]

                if d > 0.5*dimensions[dim]:                
                    d = d - dimensions[dim]
                
                if abs(d) <= space_threshold:
                    closecheck=closecheck+1

        # if the query was close to the
        # target in all n_dim dimensions
        #print "Final closecheck %i" % closecheck
        if closecheck == n_dim:
            #print "ADDING"
            in_contact.append(query)
            

    # we now have a list of positions which are close to the target, finally
    # we have to convert each position into the single period 
    in_contact_SI = []
    for query in in_contact:

        new_pos = [0]*n_dim
        for dim in range(0, n_dim):
            
            d = target[dim] - query[dim]
            
            # if we're close and don't cross
            # a boundary
            if abs(d) < space_threshold:
                new_pos[dim] = query[dim]
                continue
                
                
            ## If we get here have to modify the
            # position to get a single image
            # convention position
                
            # if we're close to the query is
            # the small side of a boundary
            # 28-[29-0]-1 -> so query = 0 target = 29
            # d = 29 - 0 = 29
            tmp = query[dim]

            if d >= space_threshold:
                while d > space_threshold:
                    tmp = tmp + dimensions[dim]
                    d = target[dim] - tmp
                # when we finish this d should be < 0 and > -(space_threshold)

            # if we're close to the query is
            # the small side of a boundary
            # 1-[0-29]-28 -> so query = 29 target = 0
            # d = 0 - 29 = -29
            elif d < -space_threshold:

                while d < -space_threshold:
                    tmp = tmp - dimensions[dim]
                    d = target[dim] - tmp

            new_pos[dim] = tmp

        in_contact_SI.append(new_pos)
        

    for query in in_contact_SI:
        for dim in range(0, n_dim):
            query[dim] = query[dim] + target_offset[dim]
    
    """
    print "<<<<>>>>>>>>>>>>>>>"
    print "ALL POSITIONS: " + str(list_of_positions)
    print ""
    print "PBC target: " + str(target)
    print "Contact positions in PBC: " + str(in_contact)
    print "SIC target: " + str(original_target)
    print "Contact positions in SIC: " +  str(in_contact_SI)
    print "end of section"
    """
                    
    return (in_contact, in_contact_SI)
                

