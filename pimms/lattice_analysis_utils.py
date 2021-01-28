## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2020
## ...........................................................................


###
### 
###

## Set of tools for analysis routines. ALL routines should 
# 1) not change any of the passed data

import numpy as np
from numpy import linalg as LA
from scipy.spatial import ConvexHull # compute volume of clusters
import scipy.spatial.qhull

from . import CONFIG
from . import lattice_utils
from . import cluster_utils
from . import numpy_utils
from .latticeExceptions import AnalysisRoutineException


def get_inter_position_distance(P1, P2, dimensions):
    """
    Returns the distance between two positions on the lattice (in real space)
    accounting for periodic boundary conditions.
   
    Routine optimized for a single distance (i.e. doesn't perform any of the
    setup/teardown used for vectorized implementations which are important
    when a set of positions are being compared)

    Arguments:
    
    P1 [list of ints]
    A position list (e.g. a list of integers specifying the X/Y or X/Y/Z coordinates of a position) 

    P2 [list of ints]
    A position list (e.g. a list of integers specifying the X/Y or X/Y/Z coordinates of a position) 

    dimensions [list of ints, 2 or 3 in length]
    Defines the box size in 2 or 3 dimensions

    
    """
    x_max = dimensions[0]
    y_max = dimensions[1]

    # convert to numpy arrays
    P1 = np.array(P1)
    P2 = np.array(P2)

    # get x/y positions
    P1_x = P1[0]
    P1_y = P1[1]
    P2_x = P2[0]
    P2_y = P2[1]

    # get vector of differences in X and Y dimensions
    x_dif = P1_x - P2_x
    y_dif = P1_y - P2_y

    # perform PBC correction for distances 
    if np.abs(x_dif) > x_max*0.5:
        x_dif = x_max - np.abs(x_dif)

    if np.abs(y_dif) > y_max*0.5:
        y_dif = y_max - np.abs(y_dif)
    
    # if we're in 3D do all the equivalent work for the 3D dimension
    if len(dimensions) == 3:
        z_max = dimensions[2]
        P1_z = P1[2]
        P2_z = P2[2]
        z_dif = P1_z - P2_z

        if np.abs(z_dif) > z_max*0.5:
            z_dif = z_max - np.abs(z_dif)
        
        distance = np.sqrt(np.power(x_dif,2) + np.power(y_dif, 2) + np.power(z_dif, 2) )

    else:
        distance = np.sqrt(np.power(x_dif,2) + np.power(y_dif, 2))

    return distance



def get_inter_position_distances(P1s, P2s, dimensions, pbc_correction=True):
    """
    Returns the list of distances between lists of two positions on the lattice (in real space)
    accounting for periodic boundary conditions.

    Optimized for multiple values - vectorizes calculations.

    Arguments:

    P1s [list of positions]
    A list of positions, where each position is a 2-length or 3-length list/set specificying X/Y/[Z] 
    coordinate positions on the lattice

    P2s [list of positions] 
    A list of positions, where each position is a 2-length or 3-length list/set specificying X/Y/[Z] 
    coordinate positions on the lattice

    dimensions [list of ints, 2 or 3 in length]
    Defines the box size in 2 or 3 dimensions

    """

    # Check lists are the same length!
    if not len(P1s) == len(P2s):
        raise AnalysisRoutineException('Two lists of positions for distance analysis did not match one another in length')


    # extract box size in X/Y dimensions
    x_max = dimensions[0]
    y_max = dimensions[1]

    # convert to numpy arrays
    P1s = np.array(P1s)
    P2s = np.array(P2s)

    # get all the X/Y positions for list 1 and list 2
    P1_x = P1s.transpose()[0]
    P1_y = P1s.transpose()[1]
    P2_x = P2s.transpose()[0]
    P2_y = P2s.transpose()[1]

    # get vector of differences in X and Y dimensions
    x_dif = P1_x - P2_x
    y_dif = P1_y - P2_y

    # perform PBC correction for distances 
    if pbc_correction:
        x_dif[np.where(abs(x_dif)>0.5*x_max)] = x_max - abs(x_dif[np.where(abs(x_dif)>x_max)])
        y_dif[np.where(abs(y_dif)>0.5*y_max)] = y_max - abs(y_dif[np.where(abs(y_dif)>y_max)])
    
    # if we're in 3D do all the equivalent work for the 3D dimension (Z)
    if len(dimensions) == 3:

        z_max = dimensions[2]
        P1_z = P1s.transpose()[2]
        P2_z = P2s.transpose()[2]
        z_dif = P1_z - P2_z
        
        # PBC correction in Z
        if pbc_correction:
            z_dif[np.where(abs(z_dif)>0.5*z_max)] = z_max - abs(z_dif[np.where(abs(z_dif)>z_max)])
        
        distance_vector = np.sqrt(np.power(x_dif,2) + np.power(y_dif, 2) + np.power(z_dif, 2) )

    else:
        distance_vector = np.sqrt(np.power(x_dif,2) + np.power(y_dif, 2))

    return distance_vector


def get_cluster_distribution(lattice_grid, chainDict):
    """
    Returns a list of lists, where each sublist contains the chainIDs associated 
    with a cluster. Cluster sublists are ordered from largest cluster to smallest.

    This is a computationally expensive algorithm that probably could be ported 
    into Cython at some point...

    Parameters
    ---------------

    lattice_grid : np.array (2D or 3D)
        Standard lattice grid

    chainDict : dict
        Standard dicionary mapping chainIDs to chain objects.

    """
    
    allChainIDs=[]    
    for chainID in chainDict:        
        allChainIDs.append(chainDict[chainID].chainID)

    num_chains = len(allChainIDs)

    # will contain lists of chains belonging to each cluster
    cluster_map = []

    # list of chains we've found so we only examine the minimum
    # number of clusters to get full coverage
    unfound_chains = set(allChainIDs)
    
    # until we've found all the chains...
    while len(unfound_chains) > 0:

        # take the first chainID from the set of unfound chains
        chainID = list(unfound_chains)[0]

        # get the set of chains in the connected component associated with chainID 

        cluster_members = lattice_utils.get_all_chains_in_connected_component(chainID, lattice_grid, chainDict, useChains=True)
        cluster_map.append(cluster_members)        
        
        # remove the found chains from the unfound chains set
        unfound_chains = unfound_chains.difference(cluster_members)

    # sort the cluster list - this sorts cluster map by the length of each sublist
    # and then reverses the order to get a list of sublists with the largest cluster
    # first. Finally, cycle through until we find a cluster smaller than the threshold, 
    # which point we're done
    clusters = sorted(cluster_map, key=len)[::-1]

    return clusters

def get_LR_cluster_distribution(latticeObject):
    """
    Returns a list of lists, where each sublist contains the chainIDs associated 
    with a cluster. Cluster sublists are ordered from largest cluster to smallest.
    LR clusters are defined as clusters were interactions are through short-range
    OR long-range interactions

    Arguments:

    lattice_grid [2D or 3D np.array]
    Standard lattice grid

    chainDict [dictionary mapping chainIDs to chain objects]
    Dictionary containing a mapping of chain objects for each chainID. 

    """
    lattice_grid = latticeObject.grid
    chainDict = latticeObject.chains

    allChainIDs=[]    
    for chainID in chainDict:        
        allChainIDs.append(chainID)

    num_chains = len(allChainIDs)

    # will contain lists of chains belonging to each cluster
    cluster_map = []

    # list of chains we've found so we only examine the minimum
    # number of clusters to get full coverage
    unfound_chains = set(allChainIDs)
    
    # until we've found all the chains...
    while len(unfound_chains) > 0:

        # take the first chainID from the set of unfound chains
        chainID = list(unfound_chains)[0]

        # get the set of chains in the connected component associated with chainID 
        #cluster_members = lattice_utils.get_all_chains_in_connected_component(chainID, lattice_grid, chainDict, useChains=True)
        cluster_members = lattice_utils.get_all_chains_in_long_range_cluster(chainID, latticeObject)
        cluster_map.append(cluster_members)        
        
        # remove the found chains from the unfound chains set
        unfound_chains = unfound_chains.difference(cluster_members)

    # sort the cluster list - this sorts cluster map by the length of each sublist
    # and then reverses the order to get a list of sublists with the largest cluster
    # first
    clusters = sorted(cluster_map, key=len)[::-1]
    
    return clusters


def get_polymeric_properties(positions, dimensions, pbc_correction=True):
    """
    Returns a list of polymeric properties calculated over the set of positions

    [0] - radius of gyration 
    [1] - asphericity

    Rg is defined as

    \sqrt(\dfrac{1}{N}\sum_{k=1}^N(r_k-r_{mean})^2)

    Where
    N = number of residues
    r_{mean} = mean residue position (Center of Mass)


    Arguments:

    positions [list of positions]
    A list of positions, where each position is a 2-length or 3-length list/set specificying X/Y/[Z] 
    coordinate positions on the lattice

    dimensions [list of ints, 2 or 3 in length]
    Defines the box size in 2 or 3 dimensions

    pbc_correction [bool, True or False]
    Defines if we should perform PBC correction here or not. For certain types of analysis
    (notably cluster analysis) the PBC correction is dealt with by the algorithms that
    construct the cluster, such that performing it again here is redundant (and generally
    not possible as the snakesearch algorithm re-positions the cluster in terms of non-periodic
    space).
    

    """
        
    COM   = lattice_utils.center_of_mass_from_positions(positions, dimensions, on_lattice=False)
    N_res = len(positions) 
    n_dim = len(dimensions)
    

    #summation=0 # commented out but left for debugging

    T_PRE = 0
    for pos in positions:
        # commented out but left for debugging
        # summation   = summation+np.square(get_inter_position_distance(pos, COM, dimensions)) # commented out but left for debugging

        # note the implementation of pbc_correct will never change A/COM, so we can just stick
        # with the COM - the 'pos' is always corrected
        if pbc_correction:
            (A,newPos) = lattice_utils.pbc_correct(COM, pos, dimensions)
        else:
            newPos = pos
        
        T_PRE = T_PRE + np.outer(np.array(newPos) - np.array(COM), np.array(newPos) - np.array(COM))
    
    T = T_PRE/len(positions)

    # get the eigenvalues of the T matrix
    (EIG, norm) = LA.eig(T)



    # if we're doing a 2D simulation
    if n_dim == 2:
        # radius of gyration from the gyration tensor
        rg =  np.sqrt(EIG[0]+EIG[1])

        #  acylindiricity --> NOTE i'm not 100% sure this is right but I think it is...
        #
        asph = (1/(rg*rg))*(abs(EIG[0]-EIG[1]))


    else:
        # radius of gyration from the gyration tensor
        rg =  np.sqrt(EIG[0]+EIG[1]+EIG[2])

        # asphericity from the gyration tensor
        asph = 1 - 3*((EIG[0]*EIG[1] + EIG[1]*EIG[2] + EIG[2]*EIG[0])/np.power(EIG[0]+EIG[1]+EIG[2],2))

    if CONFIG.DEBUG:
        if (np.sqrt(summation/N_res) - np.sqrt(EIG[0]+EIG[1]+EIG[2]))> 0.0001:
            print('Difference obtained when calculating Rg using tensor based vs. geometry based approaches')
            print("OLD WAY: " + str(np.sqrt(summation/N_res)))
            print("NEW WAY: " + rg)
            raise AnalysisRoutineException("Difference obtained when calculating Rg using tensor based vs. geometry based approaches")
            
        
            
    return [rg, asph]


def extract_positions_from_clusters(cluster_list, chainDict):
    """
    Function which takes a list of clusters (i.e. a list of lists, where 
    each sublist is a list of chainIDs in a specific cluster) and returns
    a list of lists of the same length where each sublist in the return list
    contains the positon of all residues in the cluster

    """
    
    return_list = []
    for cluster in cluster_list:
    
        sublist = []
        for chainID in cluster:
            sublist.extend(chainDict[chainID].get_ordered_positions())

        return_list.append(sublist)

    return return_list



def extract_cluster_polymeric_properties(cluster_position_list, dimensions=False):
    """
    Function which takes a list of cluster positions (i.e. a list of lists, where 
    each sublist is a list of positions associated with the residues in a specific cluster) 
    and returns a list of lists of the same length where each sublist in the return list
    contains the polymeric properties of the actual cluster.
    
    cluster_positions_list      [list of list of ints]

    List where each sublist is a list of positions. Each sub-list is its own 
    cluster. NOTE that each cluster should exist within its own single image convention, such
    that for EACH cluster we can niavely calculate things over those positions without
    needing to do any PBC related stuff.

    dimensions    [list of ints]
    Defines the dimensions of the lattice the positions sit on. However,
    if PBC correction has already been performed this can be omitted and a
    dynamic dimension can be calculated (+10 of largest value in each dimension)
    
    """
    return_list = []

    # if no positions return empty list
    if len(cluster_position_list) == 0:
        return return_list

    local_dimensions = dimensions
    n_dim = len(cluster_position_list[0][0])

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            

        if dimensions == False:
            # if we didn't explicitly define dimensions calculate them using the positions - NOTE if we 
            # don't define dimensions we cannot perform PBC corrected COM calculations so dimensions
            # should ALWAYS be supplied if the cluster_position_list is not already corrected

            # NOTE - the 10 here is kind of arbitrary, but basically we're definig a bounding box that sits 
            # around the cluster, and we're saying this box is +10 bigger than the most extreme value in every 
            # dimension. The box runs between 0 and +10 of the max

            if n_dim == 2:
                local_dimensions = [max(np.transpose(cluster)[0]+10), max(np.transpose(cluster)[1])+10]
            else:
                local_dimensions = [max(np.transpose(cluster)[0])+10, max(np.transpose(cluster)[1]+10), max(np.transpose(cluster)[2])+10]
        
        
        return_list.append(get_polymeric_properties(cluster, local_dimensions))

    return return_list



def correct_cluster_positions_to_single_image(cluster_position_list, dimensions):
    """
    Function which takes a list of cluster positions (i.e. a list of lists, where 
    each sublist is a list of positions associated with the residues in a specific cluster) 
    and for EACH CLUSTER re-configures the cluster position so the cluster is in its own single periodic image

    """

    num_clusters = len(cluster_position_list)
    
    return_list = []

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            

        # then perform single image PBC correction 
        return_list.append(cluster_utils.convert_positions_to_single_image_snakesearch(cluster, dimensions, space_threshold=1))

    return return_list


def correct_LR_cluster_positions_to_single_image(cluster_position_list, dimensions):
    """
    Function which takes a list of cluster positions (i.e. a list of lists, where 
    each sublist is a list of positions associated with the residues in a specific cluster) 
    and for EACH CLUSTER re-configures the cluster position so the cluster is in its own single periodic image

    """
    return_list = []

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            
                           
        # then perform single image PBC correction 
        return_list.append(cluster_utils.convert_positions_to_single_image_snakesearch(cluster, dimensions, space_threshold=2))

    return return_list



def compute_cluster_gross_properties(cluster_position_list):
    """
    Determines the volume of each cluster in a list of clusters. NOTE that positions
    here MUST have been corrected for PBC effects as the ConvexHull algorithm will
    calculate ASSUMES a single non-periodic image. This means that when you have clusters
    that wrap around the PBC the convex hull algorithm is calculating a single instance
    (i.e. using the boundaries as edges) so take care when extrapolating cluster volume 
    for such system spanning clusters.

    Using these positions and the Complex Hull algorithm we compute the volume, area
    and density of the cluster, generating a return_list where each element in the
    return list corresponds to the following info for each cluster

    [0] - volume
    [1] - surface area
    [2] - density

    """
    
    return_list = []

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            

        # run convex hull - if that throws an exception then
        # set everything to -1
        try:
            CH  = ConvexHull(cluster)
        except scipy.spatial.qhull.QhullError:            
            vol = -1
            SA = -1
            den = -1
            return_list.append([vol, SA, den])
            continue
        

        # Things are easy in later versions of scipy where
        # area and volume are directly computetd, but let's facilitate 
        # backwards compatibility
        # cos we're nice... 
        try:
            vol = CH.volume
            SA  = CH.area 
            den = float(len(cluster))/vol # density in residues/VOLUME [whatever unit that is!?]

        except Exception: # should make this more specific

            # earlier versions of scipy make us compute area and volume ourselves. We have implemnted
            # this volume and density in 3D but not in 2D
            if len(cluster[0]) == 3:
                simplices = np.column_stack((np.repeat(CH.vertices[0], CH.nsimplex),CH.simplices))
                tets = CH.points[simplices]
                vol = np.sum(numpy_utils.tetrahedron_volume(tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3]))
                den = float(len(cluster))/vol # density in residues/VOLUME [whatever unit that is!?]
                SA = -1.0  # ugh implemented 3D area of convex hull 
            else:
                # TO DO: Manual 2D/3D polygon area/volume calculations...
                vol = -1.0             
                SA = -1.0 
                den = -1.0

        # update lists
        return_list.append([vol, SA, den])

    return return_list
        



def compute_cluster_radial_density_profile(cluster_position_list, dimensions):
    """

    """


    ## ------------------------------------------------------------------------------------
    ## First local function
    def __position_in_list(position, pos_list):
        """
        Internal function that asks if a position exists in the list of positions.
        Assumes that 'position' is always length 3, BUT if this a 2D request
        than the Z dim will = False

        """

        # if 2D
        if position[2] == False:
            if [position[0], position[1]] in pos_list:                                                    
                return True
        else:
            if [position[0], position[1], position[2]] in pos_list:
                return True

        return False
    ## ------------------------------------------------------------------------------------

            
    ## ------------------------------------------------------------------------------------
    ## second local function
    def __extract_ring_density(COM, offset, cluster_positions, z_pos=False):
        """Internal algorithm that 'walks' around the periphery of a square, where
        that square's boundaries are set such that the center is in the COM and 
        the min/max defined by -/+ the given offset value.

        Cluster positions is a list of sites occupied by lattice beads. 

        z_pos is by default set to false, but if provided explicitly this correctly
        performs the same operation on a plane in 3D space (with Z axis fixed)

        """

        # count of number of sites with beads in
        occupied=0

        # count of TOTAL number of sites 
        total=0

        #print "All cluster positions:%s" % str(cluster_positions)
        # first do bottom row
        # C = COM, start at x, move left
        # o o o o o
        # o o o o o
        # o o C o o 
        # o o o o o 
        # x - - - - 
        #print "Cluster positions: %s"%(str(cluster_positions))
        #print offset
        #print "COM: %s" %COM
        x_pos = (COM[0] - offset)-1
        y_pos = COM[1] - offset
        #print "moving x dim [%s,%s]..." % (x_pos, y_pos)
        while x_pos < COM[0]+offset:
            total=total+1
            x_pos=x_pos+1
            #print "scanning [%i,%i]"%(x_pos,y_pos)
            if __position_in_list([x_pos, y_pos, z_pos], cluster_positions):
                occupied=occupied+1

        # first do bottom row
        # C = COM, start at x, move left
        # o o o o |
        # o o o o |
        # o o C o | 
        # o o o o x 
        # . . . . .
        # next move up left hand side
        #print "moving y dim [%s,%s]..." % (x_pos, y_pos)
        while y_pos < COM[1]+offset:
            total=total+1
            y_pos=y_pos+1
            #print "scanning [%i,%i]"%(x_pos,y_pos)
            if __position_in_list([x_pos, y_pos, z_pos], cluster_positions):
                occupied=occupied+1


        # first do bottom row
        # C = COM, start at x, move left
        # o o o x .
        # o o o o .
        # o o C o . 
        # o o o o . 
        # . . . . .
        # next do remainder of top row (note -1 to scoot one left)                                
        #print "moving x dim [%s,%s]..." % (x_pos, y_pos)
        while x_pos > COM[0]-offset:            
            total=total+1
            x_pos=x_pos-1
            #print "scanning [%i,%i]"%(x_pos,y_pos)
            if __position_in_list([x_pos, y_pos, z_pos], cluster_positions):
                occupied=occupied+1



        # . . . . .
        # x o o o .
        # o o C o . 
        # o o o o . 
        # . . . . .
        # finally do left column (note now we got to > COM[1] - offset as opposed
        # to >=
        #print "moving x dim [%s,%s]..." % (x_pos, y_pos)
        while y_pos > (COM[1]-offset)+1:
            total=total+1
            y_pos=y_pos-1
            #print "scanning [%i,%i]"%(x_pos,y_pos)
            if __position_in_list([x_pos, y_pos, z_pos], cluster_positions):
                occupied=occupied+1
            

        #print "ending at [%s,%s]..." % (x_pos, y_pos)
        
        #print "total scanned: %i" %(total)
        #print "total occupied: %i" %(occupied)
        #exit(1)
        return(occupied, total)
    ## ------------------------------------------------------------------------------------

    # initialize return densities list (will be a list of list, where each sub-list is a list
    # of radial density 
    return_densities =[]


    # for each set of positions associated with each cluster
    for cluster_positions_nd in cluster_position_list:            

        # must convert to list so we can query [x,y] in cluster_positions
        # (this works on lists but not numpy arrays)
        cluster_positions = cluster_positions_nd.tolist()
        
        # set number of beads in cluster
        num_beads = len(cluster_positions)

        # get cluster COM position
        COM = lattice_utils.center_of_mass_from_positions(cluster_positions, dimensions)
        
        # define the max offset we're going to examine
        offset_max = int((max(dimensions)/2)) -1

        # if COM is occupied bead then we'll never find that bead...
        if COM in cluster_positions:
            max_num_beads = num_beads-1
        else:
            max_num_beads = num_beads

        complete=False
        
        if len(dimensions) == 2:
            
            # not we take advantage of the flooring behaviour here for even dimension. Also
            # the -1 is because the COM position takes one space, and we need x_max*2 + 1 to
            # be == or one less than box dimensions
            
            occupied=0
            ring_density=[]
            offset = 0           

            # run until either we find all the beads OR we get bigger than the box
            #print "Num beads: %i"%(num_beads)
            while offset <= offset_max and not complete:

                
                offset = offset + 1
                
                # get density associated with ring at this offset
                (local_occupied, total) = __extract_ring_density(COM, offset, cluster_positions, z_pos=False)

                
                # update
                occupied = occupied + local_occupied
                ring_density.append(float(local_occupied)/total)

                #print "total occupied = %i (of %i) "%( occupied, max_num_beads)
                

                # check if all beads in the cluster have been found
                if occupied == max_num_beads:
                    complete=True

                if occupied > num_beads:
                    raise Exception('This should never happen and must be a bug')

            # and we're done
            if len(ring_density) < offset_max:
                ring_density.extend((offset_max - len(ring_density))*[0])

            return_densities.append(ring_density)

        else:
            
            
            ring_density=[] # variable that will become a list of average densities as a function of distance from COM
            offset = 0      # distance from COM used     
            complete=False  # flag that gets set IF we find all the beads
            occupied = 0    # counter for number of occupied lattice sites found in ALL rings

            # run until either we find all the beads OR we get bigger than the box            
            while offset <= offset_max and not complete:

                offset = offset + 1
                ring_occupied = 0     # counter for number of occupied lattice sites found in this "ring" (shell, really, because we're in 3D) 
                ring_total = 0        # TOTAL number of sites in the shell (used to normalize the volume element for de
                
                ## ------------------------------------------------------------------------
                #
                #               CENTRAL
                #  Z_offset   -2    -1    0     +1     +2
                #            ##### ##### ##### ##### #####
                #            ##### #   # #   # #   # #####
                #            ##### #   # #   # #   # #####
                #            ##### #   # #   # #   # #####
                #            ##### ##### ##### ##### #####
                #                
                # To calculate cluster density we scan through the complete planes in stage 1 and 3 and the rings in stage 2 below
                #
                
                # first fully scan the Z-plane in the -offset plane
                z_plane = COM[2] - offset                

                # stage 1
                for x_pos in range(COM[0]-offset, COM[0]+offset+1):
                    for y_pos in range(COM[1]-offset, COM[1]+offset+1):

                        # is this position found in the list of cluster positions? 
                        if __position_in_list([x_pos, y_pos, z_plane], cluster_positions):

                            # whenever we find a bead incremement the bead-occupied counter
                            ring_occupied = ring_occupied + 1

                        # regardless increment the ring total
                        ring_total = ring_total + 1

                # stage 2
                # next scan each ring 
                for z_plane in range( ( (COM[2] - offset)+1), ( (COM[2] + offset)-1) +1):
                    (local_occupied, local_total) = __extract_ring_density(COM, offset, cluster_positions, z_pos=z_plane)
                    ring_total = ring_total + local_total
                    ring_occupied = ring_occupied + local_occupied
                       
                # stage 3
                # finally fully scan the terminal Z-plane in the +offset plane
                z_plane = COM[2] + offset                
                for x_pos in range(COM[0]-offset, COM[0]+offset+1):
                    for y_pos in range(COM[1]-offset, COM[1]+offset+1):
                        if __position_in_list([x_pos, y_pos, z_plane], cluster_positions):
                            ring_occupied = ring_occupied + 1
                        ring_total=ring_total+1
                ## ------------------------------------------------------------------------

                # update
                ring_density.append(float(ring_occupied)/ring_total)
                occupied = occupied+ring_occupied

                # check if all beads in the cluster have been found
                if occupied == max_num_beads:
                    complete=True


                if occupied > num_beads:
                    raise Exception('This should never happen and must be a bug')

            # and we're done
            if len(ring_density) < offset_max:
                ring_density.extend( (offset_max - len(ring_density))*[0]  )

            return_densities.append(ring_density)
    

    return return_densities

            


                    

                    
                
                

                    
            
            
            
        
        
