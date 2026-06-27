## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


###
### 
###

## Set of tools for analysis routines. ALL routines should 
# 1) not change any of the passed data

import numpy as np
from numpy import linalg as LA
from scipy.spatial import ConvexHull # compute volume of clusters
from scipy.spatial import QhullError

from . import CONFIG
from . import lattice_utils
from . import cluster_utils
from . import numpy_utils
from .latticeExceptions import AnalysisRoutineException


def get_inter_position_distance(P1, P2, dimensions, pbc_correction=True):
    """
    Returns the distance between two positions on the lattice (in real space)
    accounting for periodic boundary conditions.
   
    Routine optimized for a single distance (i.e. doesn't perform any of the
    setup/teardown used for vectorized implementations which are important
    when a set of positions are being compared)

    Parameters
    -----------------
    
    P1 : [list of ints]
        A position list (e.g. a list of integers specifying the X/Y or X/Y/Z 
        coordinates of a position) 

    P2 : [list of ints]
        A position list (e.g. a list of integers specifying the X/Y or X/Y/Z 
        coordinates of a position) 

    dimensions : [list of ints, 2 or 3 in length]
        Defines the box size in 2 or 3 dimensions

    pbc_correction : bool
        Flag which if set to true means a PBC correction is applied. Default
        is True.

    Returns
    -------
    float
        The (optionally PBC-corrected) Euclidean distance between ``P1`` and
        ``P2`` in real space.

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
    if pbc_correction:
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

        if pbc_correction:
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

    Parameters
    ----------
    P1s : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    P2s : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    pbc_correction : bool, optional
        If True (default), the minimum-image PBC correction is applied to each
        per-dimension separation before computing distances.

    Returns
    -------
    numpy.ndarray
        1D array of (optionally PBC-corrected) Euclidean distances, one per
        position pair.

    Raises
    ------
    AnalysisRoutineException
        If ``P1s`` and ``P2s`` do not have the same length.

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
        # minimum image convention: where |d| > L/2, replace with L - |d|. The
        # selection mask must be identical on both sides (previously the RHS used
        # the always-empty mask abs(d) > L, which crashed on any over-half-box
        # separation due to a shape mismatch).
        x_mask = np.abs(x_dif) > 0.5*x_max
        x_dif[x_mask] = x_max - np.abs(x_dif[x_mask])
        y_mask = np.abs(y_dif) > 0.5*y_max
        y_dif[y_mask] = y_max - np.abs(y_dif[y_mask])
    
    # if we're in 3D do all the equivalent work for the 3D dimension (Z)
    if len(dimensions) == 3:

        z_max = dimensions[2]
        P1_z = P1s.transpose()[2]
        P2_z = P2s.transpose()[2]
        z_dif = P1_z - P2_z
        
        # PBC correction in Z
        if pbc_correction:
            z_mask = np.abs(z_dif) > 0.5*z_max
            z_dif[z_mask] = z_max - np.abs(z_dif[z_mask])
        
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

    Returns
    -------
    list of list of int
        List of clusters, where each sublist contains the chainIDs of the
        chains in that connected component. Clusters are ordered from largest
        to smallest.

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

    Parameters
    ---------------

    latticeObject : Lattice
        The lattice object. Its ``grid`` (the lattice grid) and ``chains``
        (mapping of chainIDs to chain objects) attributes are used, along with
        long-range interaction information, to build the long-range clusters.

    Returns
    -------
    list of list of int
        List of clusters, where each sublist contains the chainIDs of the
        chains in that long-range connected component. Clusters are ordered
        from largest to smallest.

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



def get_eigenvalues_of_the_T_matrix(positions, dimensions, pbc_correction=True):
    """
    Compute the eigenvalues and eigenvectors of the gyration (T) tensor.

    Builds the gyration tensor from the supplied positions relative to their
    (optionally PBC-corrected) center of mass, then diagonalizes it. The
    eigenvalues are the principal components used downstream to compute the
    radius of gyration and asphericity.

    Parameters
    ----------
    positions : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    pbc_correction : bool, optional
        If True (default), each position is PBC-corrected relative to the
        center of mass before contributing to the gyration tensor.

    Returns
    -------
    tuple
        ``(EIG, norm)`` where ``EIG`` is the array of eigenvalues of the
        gyration tensor and ``norm`` is the matrix of corresponding
        eigenvectors (as returned by ``numpy.linalg.eig``).
    """

    # NB: we have verified that even though the center_of_mass_from_positions algorithm
    # seems to have some issues with PBC, the gyration tensor is unaffected and so
    # this code continues to return value values for the gyration tensor that are
    # PBC-correct
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

    return (EIG, norm)

    


def get_polymeric_properties(positions, dimensions, pbc_correction=True):
    r"""
    Returns a list of polymeric properties calculated over the set of positions

    [0] - radius of gyration 
    [1] - asphericity

    Rg is defined as

    \sqrt(\dfrac{1}{N}\sum_{k=1}^N(r_k-r_{mean})^2)

    Where
    N = number of residues
    r_{mean} = mean residue position (Center of Mass)


    Parameters
    ----------
    positions : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    pbc_correction : bool, optional
        Defines whether to perform PBC correction here (default True). For
        certain types of analysis (notably cluster analysis) the PBC correction
        is dealt with by the algorithms that construct the cluster, such that
        performing it again here is redundant (and generally not possible, as
        the snakesearch algorithm re-positions the cluster in terms of
        non-periodic space).

    Returns
    -------
    list of float
        A two-element list ``[rg, asph]`` where ``rg`` is the radius of
        gyration and ``asph`` is the asphericity (acylindricity in 2D), both
        derived from the gyration-tensor eigenvalues. Degenerate cases where
        ``rg ~ 0`` return an asphericity of 0.0.

    """

    n_dim = len(dimensions)


    # compute the eigenvalues and normal of the T matrix. NOTE - the function below USED
    # to be part of this function but we extracted it out
    (EIG, norm) = get_eigenvalues_of_the_T_matrix(positions, dimensions, pbc_correction)
        
    # Numerical tolerance for degenerate chains/clusters where Rg ~ 0.
    eps = 1e-12

    # if we're doing a 2D simulation
    if n_dim == 2:
        # radius of gyration from the gyration tensor
        rg2 = max(0.0, EIG[0] + EIG[1])
        rg = np.sqrt(rg2)

        # acylindiricity. For degenerate cases (rg == 0), define asphericity as 0.
        if rg2 <= eps:
            asph = 0.0
        else:
            asph = abs(EIG[0] - EIG[1]) / rg2


    else:
        # radius of gyration from the gyration tensor
        rg_sum = max(0.0, EIG[0] + EIG[1] + EIG[2])
        rg = np.sqrt(rg_sum)

        # asphericity from the gyration tensor
        denom = np.power(rg_sum, 2)
        if denom <= eps:
            asph = 0.0
        else:
            asph = 1 - 3 * ((EIG[0] * EIG[1] + EIG[1] * EIG[2] + EIG[2] * EIG[0]) / denom)

    return [rg, asph]


def extract_positions_from_clusters(cluster_list, chainDict):
    """
    Function which takes a list of clusters (i.e. a list of lists, where 
    each sublist is a list of chainIDs in a specific cluster) and returns
    a list of lists of the same length where each sublist in the return list
    contains the positon of all residues in the cluster

    Parameters
    ----------
    cluster_list : list of list of int
        List of clusters, where each sublist is a list of chainIDs in that
        cluster.

    chainDict : dict
        Dictionary mapping each chainID to its chain object (each chain object
        must expose ``get_ordered_positions()``).

    Returns
    -------
    list of list
        List of the same length as ``cluster_list`` where each sublist contains
        the ordered positions of all residues belonging to that cluster.

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

    Parameters
    ----------
    cluster_position_list : list of list of positions
        List where each sublist is a list of positions; each sublist is its own
        cluster. NOTE that each cluster should exist within its own single-image
        convention, so that for each cluster the properties can be computed
        naively over those positions without any further PBC handling.

    dimensions : list of int or bool, optional
        Defines the dimensions of the lattice the positions sit on. If PBC
        correction has already been performed this can be left as ``False``
        (the default), in which case a per-cluster bounding box is computed
        dynamically (+10 beyond the largest value in each dimension).

    Returns
    -------
    list of list of float
        List of the same length as ``cluster_position_list`` where each entry
        is the ``[rg, asph]`` polymeric properties of the corresponding
        cluster (empty list if no clusters are supplied).

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
    and for EACH CLUSTER re-configures the cluster position so the cluster is in its own single 
    periodic image.

    Parameters
    -------------
    cluster_postion_list : list
        A list of lists, each sublist is a list of cluster positions (where, in fact, each position
        is also itself a list of 2 or 3 positions

    dimensions : list
        A list of 2- or 3- elements that defines the X/Y or X/Y/Z positions

    Returns
    ----------
    list
        List of the same length as ``cluster_position_list`` where each entry
        is the cluster's positions re-expressed in a single (non-periodic)
        image, as returned by
        ``cluster_utils.convert_positions_to_single_image_snakesearch`` with a
        ``space_threshold`` of 1.

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

    Identical to :func:`correct_cluster_positions_to_single_image` but uses a
    ``space_threshold`` of 2, appropriate for long-range (LR) clusters whose
    members may be separated by more than one lattice site.

    Parameters
    ----------
    cluster_position_list : list
        A list of lists; each sublist is a list of cluster positions (where
        each position is itself a 2- or 3-element list).

    dimensions : list
        A list of 2 or 3 elements defining the X/Y or X/Y/Z box dimensions.

    Returns
    -------
    list
        List of the same length as ``cluster_position_list`` where each entry
        is the cluster's positions re-expressed in a single (non-periodic)
        image, using a ``space_threshold`` of 2.

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
    and density of the cluster.
    
    Parameters
    -----------------
    cluster_position_list : list
        list of np.ndarrays where dimensions of the np.ndarray reflect dimensions of the
        lattice. 


    Returns
    -------------
    list of lists
        Returns a list of lists, where each sublist contains three elements that reflect
        the cluster gross properties. Sublist indices match indices for cluster_position_list
        indices
    
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
        except QhullError:
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
        



def compute_cluster_radial_density_profile(cluster_position_list, dimensions, minimum_cluster_size_in_beads=None):
    """
    Compute the radial density profile of each cluster about its center of mass.

    For each cluster the lattice is scanned outward from the cluster center of
    mass in concentric square rings (2D) or cubic shells (3D). For every shell
    at a given offset the fraction of occupied lattice sites is recorded,
    yielding a density-versus-distance profile. Scanning stops once all beads in
    the cluster have been accounted for (or the box half-extent is reached), and
    short profiles are zero-padded to a common length.

    Parameters
    ----------
    cluster_position_list : list of numpy.ndarray
        List of clusters; each entry is an array of lattice positions (each
        position being a 2- or 3-element coordinate) for that cluster. Positions
        are expected to be single-image (PBC-corrected) per cluster.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    minimum_cluster_size_in_beads : int or None, optional
        If supplied, clusters with fewer beads than this threshold are skipped
        (no profile is emitted for them). Default is None (no filtering).

    Returns
    -------
    list of list of float
        One radial density profile per (non-skipped) cluster; each profile is a
        list of occupied-site fractions as a function of distance from the
        cluster center of mass, zero-padded to a uniform length.

    Raises
    ------
    Exception
        If the number of occupied sites found ever exceeds the number of beads
        in the cluster, which indicates an internal bug.
    """


    ## ------------------------------------------------------------------------------------
    ## First local function
    def __position_in_list(position, pos_list):
        """
        Internal function that asks if a position exists in the list of positions.
        Assumes that 'position' is always length 3, BUT if this a 2D request
        than the Z dim will = False

        """

        # if 2D (note: must use 'is False', not '== False', because an integer
        # z-coordinate of 0 compares equal to False and would wrongly trigger the
        # 2D path, dropping the entire z=0 plane from 3D density profiles)
        if position[2] is False:
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

        # IF we've defined a smallest cluster worth computing for skip
        if minimum_cluster_size_in_beads is not None:
            if num_beads < minimum_cluster_size_in_beads:
                continue

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

            


                    

                    
                
                

                    
            
            
            
        
        
