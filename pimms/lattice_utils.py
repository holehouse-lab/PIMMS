## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


##
## lattice_utils
##
## lattice_utils contains system agnostic, stateless utilities for lattice operations.
## Utilities are relevant for 2D and 3D lattices
##

import random
import copy
import math
import os
import numpy as np

import mdtraj as md

from .latticeExceptions import ChainInsertionFailure, ChainDeletionFailure, ResidueAugmentException, MoveSetException, ChainConnectivityError, ClusterSizeThresholdException, LatticeUtilsException, RotationException

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
from . CONFIG import NP_INT_TYPE

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

    if len(site1) != len(site2):
        raise LatticeUtilsException("Position dimensionality mismatch in same_sites")

    if len(site1) not in (2, 3):
        raise LatticeUtilsException(f"Unsupported dimensionality in same_sites: {len(site1)}")

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

    dimensions : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects the lattice dimensions used as the modulus
        for each dimension.

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

    Parameters
    ----------
    posA : list
        A list of length 2 or 3, depending on the dimensionality of the system,
        giving a position on the lattice. This position is treated as fixed.

    posB : list
        A list of length 2 or 3 (matching `posA`) giving a second position on
        the lattice. This position may be shifted by a single box length in any
        dimension so that it lies in the same periodic image as `posA`.

    dimensions : list
        A list of length 2 or 3 giving the lattice dimensions, used as the box
        width in each dimension.

    Returns
    -------
    tuple
        A 2-tuple ``(posA, newB)`` where ``posA`` is the unchanged input
        position and ``newB`` is `posB` shifted (where required) so that, in
        each dimension, the separation between the two positions is the minimum
        image separation.

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

    The check is performed by walking through consecutive positions and asking,
    in each dimension, whether the absolute difference between adjacent
    positions exceeds 1. A difference greater than 1 implies the bond wraps
    across a periodic boundary.

    Parameters
    ----------
    chain_positions : list
        A list of positions (each a length-2 or length-3 list) that are assumed
        to be consecutively connected along a chain.

    Returns
    -------
    bool
        True if any consecutive pair of positions straddles a periodic
        boundary, otherwise False.

    """

    n_dim = len(chain_positions[0])

    for pidx in range(0, len(chain_positions)-1):
        p1 = chain_positions[pidx]
        p2 = chain_positions[pidx+1]
        for xyz in range(0,n_dim):
            if abs(p1[xyz]-p2[xyz]) > 1:
                return True

    return False


def center_positions(positions, dimensions):
    """
    Returns the positions after centering them in the box. This is useful
    for visualisation and also for calculating the radial distribution
    function (RDF) as it ensures that the RDF is not biased by the
    positions being offset from the centre of the box.

    Added in v0.1.34

    Parameters
    -----------
    positions : list
        A list of lists, where each inner list is a position on the lattice.

    dimensions : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects the lattice dimensions.

    Returns
    ---------
    list
        A list of lists, where each inner list is a position on the lattice,
        that has been centered in the box.

    """

    n_dim = len(dimensions)

    com = center_of_mass_from_positions(positions, dimensions)

    offset = []
    for idx in range(0, n_dim):
        offset.append(dimensions[idx]/2 - com[idx])
        
    new_positions = []

    for pos in positions:
        new_pos = []
        for idx in range(0, n_dim):
            new_pos.append(pos[idx] + offset[idx])
        new_positions.append(new_pos)

    return new_positions


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

    Parameters
    ----------
    chain_of_positions : list
        A list of 2D or 3D positions describing a chain, where each consecutive
        position is assumed to be one lattice site apart from the next. The
        first position defines the periodic image used as the reference frame.

    dimensions : list
        A list of length 2 or 3 giving the dimensions of the box, used as the
        box width in each dimension when unwrapping positions across periodic
        boundaries.

    Returns
    -------
    list
        A list of positions, the same length as `chain_of_positions`, unwrapped
        into a single periodic image and then offset so that all coordinates in
        every dimension are non-negative.

    Raises
    ------
    LatticeUtilsException
        Raised if the chain cannot be unwrapped into a single image within the
        internal escape-counter limit, which indicates the input chain contains
        impossible (non-unit) bonds.

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
                escape_counter = 0
                while abs(next_pos[dim] - current[dim]) > 1:                    
                    next_pos[dim] = next_pos[dim] + dimensions[dim]
                    escape_counter += 1
                    if escape_counter > 100000:
                        raise LatticeUtilsException("Error building single image convention - suggests input chain may have impossible bonds")

            # i.e. if we're at the edge of a boundary and the next postion is the 
            # other the side (e.g 2-1-[0-29]-28-27)
            elif current[dim] - local_positions[pidx][dim] < -1:
                next_pos[dim] = local_positions[pidx][dim] - dimensions[dim]

                # finally this loop lets us account for chains that span multiple periodic images
                escape_counter = 0
                while abs(next_pos[dim] - current[dim]) > 1:
                    next_pos[dim] = next_pos[dim] - dimensions[dim]
                    escape_counter += 1
                    if escape_counter > 100000:
                        raise LatticeUtilsException("Error building single image convention - suggests input chain may have impossible bonds")

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
def make_chain_whole(chain_of_positions, dimensions):
    """
    Unwrap a chain into a single periodic image, anchored at the FIRST bead's
    real (in-box) position.

    Each consecutive bead is placed exactly one lattice step from the previous one,
    so the chain is never torn across a periodic boundary; coordinates may fall
    outside the box on either side. Unlike :func:`convert_chain_to_single_image`,
    this does NOT afterwards shift the whole chain to be non-negative - that shift
    would translate every boundary-crossing chain toward one face of the box, which
    is why an unwrapped trajectory built with the single-image routine only ever
    appeared to bulge out of one side. Keeping the first bead in place makes chains
    spill symmetrically out of whichever face they actually cross.

    This is used purely for trajectory visualisation (``TRAJECTORY_PBC_UNWRAP``) and
    is a lightweight, allocation-cheap loop (no deep copies / numpy transposes),
    since it is called once per chain per written frame.

    Parameters
    ----------
    chain_of_positions : list
        List of 2D or 3D integer positions describing a chain, consecutive beads
        one lattice site apart. The first position defines the reference image.

    dimensions : list
        Box dimensions (length 2 or 3), used as the per-axis period.

    Returns
    -------
    list
        The chain positions unwrapped into a single image, first bead unchanged.

    Raises
    ------
    LatticeUtilsException
        If a bond cannot be resolved within the escape-counter limit (indicates an
        impossible / non-unit bond in the input chain).
    """
    n_dim = len(chain_of_positions[0])
    n_pos = len(chain_of_positions)

    out = [list(chain_of_positions[0])]
    current = list(chain_of_positions[0])

    for pidx in range(1, n_pos):
        nxt = [0] * n_dim
        for dim in range(0, n_dim):
            v = chain_of_positions[pidx][dim]

            if current[dim] == v:
                nxt[dim] = v

            # neighbour sits across the +boundary (e.g. 27-28-[29-0]-1-2)
            elif current[dim] - v > 1:
                v = v + dimensions[dim]
                escape_counter = 0
                while abs(v - current[dim]) > 1:
                    v = v + dimensions[dim]
                    escape_counter += 1
                    if escape_counter > 100000:
                        raise LatticeUtilsException("Error making chain whole - suggests input chain may have impossible bonds")
                nxt[dim] = v

            # neighbour sits across the -boundary (e.g. 2-1-[0-29]-28-27)
            elif current[dim] - v < -1:
                v = v - dimensions[dim]
                escape_counter = 0
                while abs(v - current[dim]) > 1:
                    v = v - dimensions[dim]
                    escape_counter += 1
                    if escape_counter > 100000:
                        raise LatticeUtilsException("Error making chain whole - suggests input chain may have impossible bonds")
                nxt[dim] = v

            else:
                nxt[dim] = v

        out.append(nxt)
        current = nxt

    return out



#-----------------------------------------------------------------
#
def get_adjacent_sites_3D(position1, position2, position3, dimensions, extent_range=1):
    """
    Returns the lattice sites adjacent to a 3D position, generalized over an
    optional neighbourhood extent.

    This is a thin wrapper around the compiled ``hyperloop.get_adjacent_sites_3D``
    routine, which performs the periodic-boundary-corrected enumeration of
    neighbouring sites.

    Parameters
    ----------
    position1 : int
        The x-coordinate of the position whose neighbours are required.

    position2 : int
        The y-coordinate of the position whose neighbours are required.

    position3 : int
        The z-coordinate of the position whose neighbours are required.

    dimensions : list
        A list of length 3 giving the lattice dimensions (X, Y, Z).

    extent_range : int
        Half-width of the neighbourhood to enumerate around the position
        (default 1, i.e. immediately adjacent sites).

    Returns
    -------
    numpy.ndarray
        A 4D numpy array where each element is a 3D numpy array describing an
        adjacent lattice position.

    """
    return(hyperloop.get_adjacent_sites_3D(position1, position2, position3, dimensions[0], dimensions[1], dimensions[2], extent_range))



#-----------------------------------------------------------------
#
def get_adjacent_sites_2D(position1, position2, dimensions, extent_range=1):
    """
    Returns the lattice sites adjacent to a 2D position, generalized over an
    optional neighbourhood extent.

    This is a thin wrapper around the compiled ``hyperloop.get_adjacent_sites_2D``
    routine, which performs the periodic-boundary-corrected enumeration of
    neighbouring sites.

    Parameters
    ----------
    position1 : int
        The x-coordinate of the position whose neighbours are required.

    position2 : int
        The y-coordinate of the position whose neighbours are required.

    dimensions : list
        A list of length 2 giving the lattice dimensions (X, Y).

    extent_range : int
        Half-width of the neighbourhood to enumerate around the position
        (default 1, i.e. immediately adjacent sites).

    Returns
    -------
    numpy.ndarray
        A 3D numpy array where each element is a 2D numpy array describing an
        adjacent lattice position.

    """
    return(hyperloop.get_adjacent_sites_2D(position1, position2, dimensions[0], dimensions[1], extent_range))



#-----------------------------------------------------------------
#
def find_nearest_position(target, positions_list, dimensions):
    """
    Given some target position (target) what is the index
    of the position in the positions list that is closest
    to that target? If multiple positions are found we
    simply return the first one in the positions list.

    Parameters
    ----------
    target : list
        A 2D or 3D position (list of ints) to which all positions in
        `positions_list` are compared.

    positions_list : list
        A list of 2D or 3D positions (each a list of ints) to be compared
        against the target.

    dimensions : list
        A list of length 2 or 3 giving the lattice dimensions, used when
        computing the real (Euclidean, periodic) distance.

    Returns
    -------
    tuple
        A 2-tuple ``(int, float)`` where element [0] is the index of the
        position in `positions_list` closest to `target`, and element [1] is the
        actual distance between `target` and that closest position.

    Raises
    ------
    LatticeUtilsException
        Raised if `positions_list` is empty.

    """

    if len(positions_list) == 0:
        raise LatticeUtilsException("Error in lattice_utils.find_nearest_position() - possition_list is empty")

    # Vectorized minimum-image search. This used to be a Python loop calling the scalar
    # distance helper once per position; for a large condensate that made choosing the
    # snakesearch seed cost roughly ten times as much as the compiled BFS it feeds.
    #
    # np.argmin returns the FIRST minimum, which reproduces the strict `<` comparison of
    # the old loop (ties keep the earliest position), and comparing squared distances is
    # order-equivalent to comparing distances, so the selected index is unchanged.
    positions = np.asarray(positions_list)
    dims = np.asarray(dimensions)

    d = np.abs(positions - np.asarray(target))
    d = np.where(d > 0.5 * dims, dims - d, d)
    squared_distances = (d * d).sum(axis=1)

    min_idx = int(np.argmin(squared_distances))

    return (min_idx, math.sqrt(squared_distances[min_idx]))

 

#-----------------------------------------------------------------
#
def get_empty_site(lattice_grid, adjacentTo=None, hardwall=False):
    """
    Function which returns the position of an empty site on the lattice.
    If adjacentTo is populated then the returned site is adjacent to the
    the position defined by adjacentTo.

    When `adjacentTo` is None the function performs random rejection sampling
    of the whole lattice until an empty site is found. When `adjacentTo` is set,
    only the sites neighbouring that position are considered (optionally
    excluding sites that straddle the boundary when `hardwall` is True), and one
    empty neighbour is selected at random.

    Parameters
    ----------
    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array, where a value of 0 denotes an empty
        (solvent) site.

    adjacentTo : list, optional
        If provided, a 2D or 3D position; the returned empty site is constrained
        to be adjacent to this position. If None (default), an empty site
        anywhere on the lattice is returned.

    hardwall : bool, optional
        If True (and `adjacentTo` is set), only neighbouring sites that do not
        straddle a periodic boundary are considered. Default is False.

    Returns
    -------
    list or tuple
        If `adjacentTo` is None, returns a single position (list) of an empty
        site. If `adjacentTo` is set, returns a 2-tuple ``(position, found)``
        where `position` is the chosen empty neighbour (or a position filled
        with -1 if none was found) and `found` is a bool indicating success.

    Raises
    ------
    LatticeUtilsException
        Raised (when `adjacentTo` is None) if the lattice is fully occupied or
        no empty site is found within the attempt limit.

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
            return ([-1] * len(dimensions), False)

        # return the list of good sites - note we're casting to a list so as
        # we return lists rather than np.arrays
        else:
            return (list(empty_list[random.randint(0,len(empty_list)-1)]), True)

    # if we're literally just looking for an empty site anywhere on the lattice
    else:
        # Fail fast if the lattice is already fully occupied.
        if not np.any(lattice_grid == 0):
            raise LatticeUtilsException("Unable to find empty lattice site: lattice appears fully occupied")

        empty = False
        count = 0
        max_attempts = max(1000, int(np.prod(dimensions)) * 10)
        while not empty:
            count=count+1

            if count % 100 == 0:
                IO_utils.status_message("Tried %i times but unable to insert a single point into an empty space - maybe grid is full?\nWill keep trying though, cos I'm a trooper!" % count, 'warning')

            if count > max_attempts:
                raise LatticeUtilsException(
                    f"Unable to find empty lattice site after {count} attempts in dimensions {dimensions}"
                )

                                            
            # 2D
            if len(dimensions) == 2:

                # select a random position
                x = NP_INT_TYPE(random.randint(0, dimensions[0]-1))
                y = NP_INT_TYPE(random.randint(0, dimensions[1]-1))
                    
                # if the possition is empty celebrate with a beer!
                if get_gridvalue([x,y], lattice_grid) == 0:
                    position = [x,y]
                    empty=True

            # 3D
            if len(dimensions) == 3:
                x = NP_INT_TYPE(random.randint(0, dimensions[0]-1))
                y = NP_INT_TYPE(random.randint(0, dimensions[1]-1))
                z = NP_INT_TYPE(random.randint(0, dimensions[2]-1))

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

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array into which the chain is inserted. This
        array is modified in place as residues are placed.

    default_start : list, optional
        If provided, a position used as the fixed starting site for chain
        growth. If None (default) a random empty site is selected as the start.

    hardwall : bool, optional
        If True, chain growth only considers neighbouring sites that do not
        straddle a periodic boundary. Default is False.

    Returns
    -------
    list
        A list of positions describing the inserted chain, ordered along the
        chain. The lattice grid is also updated in place.

    Raises
    ------
    ChainInsertionFailure
        Raised if a valid chain configuration could not be constructed within
        ``CONFIG.CHAIN_INIT_ATTEMPTS`` attempts.

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
                        IO_utils.status_message(f"Chain (ID={chainID}) construction failed [TRY {attempt+1} of {CONFIG.CHAIN_INIT_ATTEMPTS}]", 'warning')
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

    Parameters
    ----------
    positions : list
        A list of positions (each a 2D or 3D list) that will be assigned to the
        chain.

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array, modified in place.

    chainID : int
        The chain ID value written into the lattice grid at each position.

    safe : bool, optional
        If True, each target position is checked to ensure it is currently empty
        before writing; an occupied site raises an exception. If False (default)
        positions are overwritten without checking.

    Returns
    -------
    None

    Raises
    ------
    ChainInsertionFailure
        Raised (only when `safe` is True) if any target position is already
        occupied by another chain.

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

    All lattice sites whose value equals `chainID` are reset to 0.0 (solvent).

    Parameters
    ----------
    chainID : int
        The chain ID whose residues should be removed from the lattice.

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array, modified in place.

    Returns
    -------
    None

    """
    lattice_grid[lattice_grid == chainID] = 0.0



#-----------------------------------------------------------------
#
def delete_chain_by_position(positions, lattice_grid, chainID=None):
    """
    Deletes a chain based on supplied position. Can be an entire chain
    or simply a portion of a chain

    Parameters
    ----------
    positions : list
        A list of positions (each a 2D or 3D list) to be reset to solvent (0).

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array, modified in place.

    chainID : int, optional
        If provided, each position is checked to confirm it is currently
        occupied by this chain before being cleared; a mismatch raises an
        exception. If None (default) positions are cleared without checking.

    Returns
    -------
    None

    Raises
    ------
    ChainDeletionFailure
        Raised (only when `chainID` is provided) if any position is not occupied
        by the expected chain.

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

    The dimensionality (2D or 3D) is inferred from the shape of the lattice
    grid.

    Parameters
    ----------
    position : list
        A 2D or 3D position (list of ints) to look up.

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array.

    Returns
    -------
    int or float
        The grid value stored at `position` (0 denotes an empty/solvent site,
        otherwise the occupying chain ID).

    Raises
    ------
    LatticeUtilsException
        Raised if the lattice grid has an unsupported dimensionality.

    """

    # use the grid's own rank + a single tuple index. This is a hot path in the
    # cluster connected-component search (called ~once per bead-neighbour), so we
    # avoid the per-call get_dimensions() (grid.shape) lookup and the chained
    # __getitem__ (which builds intermediate array views).
    ndim = lattice_grid.ndim

    if ndim == 2:
        return lattice_grid[position[0], position[1]]

    if ndim == 3:
        return lattice_grid[position[0], position[1], position[2]]

    raise LatticeUtilsException(f"Unsupported lattice dimensionality in get_gridvalue: {ndim}")


#-----------------------------------------------------------------
#                        
def get_gridvalue_2D(position, lattice_grid):
    """
    Returns the value on the lattice grid associated with a 2D position.

    This is a dimensionality-specialized variant of :func:`get_gridvalue` that
    assumes a 2D position and avoids the dimensionality check for speed.

    Parameters
    ----------
    position : list
        A 2D position (list of two ints) to look up.

    lattice_grid : numpy.ndarray
        The 2D lattice grid array.

    Returns
    -------
    int or float
        The grid value stored at `position` (0 denotes an empty/solvent site,
        otherwise the occupying chain ID).

    """
    return lattice_grid[position[0]][position[1]]



#-----------------------------------------------------------------
#                        
def get_gridvalue_3D(position, lattice_grid):
    """
    Returns the value on the lattice grid associated with a 3D position.

    This is a dimensionality-specialized variant of :func:`get_gridvalue` that
    delegates to the compiled ``hyperloop.get_gridvalue_3D`` routine for speed.

    Parameters
    ----------
    position : list
        A 3D position (list of three ints) to look up.

    lattice_grid : numpy.ndarray
        The 3D lattice grid array.

    Returns
    -------
    int or float
        The grid value stored at `position` (0 denotes an empty/solvent site,
        otherwise the occupying chain ID).

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

    Parameters
    ----------
    position : list
        A 2D or 3D position (list of ints) at which to write.

    value : int or float
        The value to write into the lattice grid at `position` (e.g. a chain ID
        or 0 for solvent).

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array, modified in place.

    Returns
    -------
    numpy.ndarray
        The same lattice grid object that was passed in, after modification.

    Raises
    ------
    LatticeUtilsException
        Raised if `position` has an unsupported dimensionality.

    """

    if len(position) == 2:
        lattice_grid[position[0]][position[1]] = value

    if len(position) == 3:
        lattice_grid[position[0]][position[1]][position[2]] = value

    if len(position) not in (2, 3):
        raise LatticeUtilsException(f"Unsupported position dimensionality in set_gridvalue: {len(position)}")

    return lattice_grid


#-----------------------------------------------------------------
#
def _unique_rows(rows):
    """Duplicate-free rows of a 2D integer array (order unspecified).

    The classic numpy idiom for this - viewing each row as a single ``np.void`` scalar
    and calling ``np.unique`` - is a full lexicographic sort, and it showed up as one of
    the largest single costs of a PIMMS analysis step once the surrounding Python loops
    were removed.

    Where the coordinates are small enough to pack losslessly into one int64 (which they
    always are for a lattice: every value is a box coordinate), each row is folded into a
    single integer key and deduplicated on that instead - the same sort, but over one
    column of scalars rather than an n-column structured view. Anything that does not fit
    falls back to the void view, so correctness never depends on the box being small.
    """
    if len(rows) == 0:
        return rows

    lo = int(rows.min())
    hi = int(rows.max())
    span = hi - lo + 1
    n_cols = rows.shape[1]

    # can we pack n_cols digits of base `span` into a signed 64-bit integer?
    packable = span > 0
    if packable:
        limit = 1
        for _ in range(n_cols):
            limit *= span
            if limit > 2 ** 62:
                packable = False
                break

    if packable:
        keys = np.zeros(len(rows), dtype=np.int64)
        for col in range(n_cols):
            keys *= span
            keys += rows[:, col].astype(np.int64) - lo
        _, idx = np.unique(keys, return_index=True)
    else:
        view = np.ascontiguousarray(rows).view(
            np.dtype((np.void, rows.dtype.itemsize * n_cols)))
        _, idx = np.unique(view, return_index=True)

    return rows[idx]


#-----------------------------------------------------------------
#
def build_envelope_pairs(positions, dimensions, hardwall=False, deduplicate=True):
    """
    Expects a LIST of positions. Returns a unique unordered
    list of tuples, where each tuple is a pair of positions.

    The complete set of these positions represents the non-redundant 
    set of positions that make contact with the positions in the past
    $positions variable.

    dimensions is the dimensions of the lattice

    hardwall is a boolean which determines if we allow a pair of
    positions to straddle the boundary (in a periodic manner) or not

    Parameters
    ----------
    positions : list
        A list of positions (each a 2D or 3D list) for which the enveloping
        short-range contact pairs are required.

    dimensions : list
        A list of length 2 or 3 giving the lattice dimensions.

    hardwall : bool, optional
        If True, pairs that straddle the periodic boundary are excluded
        (hardwall variant). Default is False.

    deduplicate : bool, optional
        If True (default) duplicate pairs are removed, which is required whenever the
        pairs are summed over (e.g. an energy evaluation, where a repeated pair would be
        double counted). Deduplication is a full lexicographic sort of the pair array
        and is the dominant cost of this function, so callers that only ask *which sites
        are touched* - the connected-component searches, which feed the pairs straight
        into a set - can and should pass False.

    Returns
    -------
    numpy.ndarray
        A numpy array of shape (N, 2, 2) in 2D or (N, 2, 3) in 3D, where each
        element is an unordered pair of positions making short-range contact
        with the input positions. Duplicate pairs are removed unless
        ``deduplicate`` is False. An empty array of the appropriate shape is
        returned if `positions` is empty.

    """

    if len(positions) == 0:
        if len(dimensions) == 2:
            return np.empty((0, 2, 2), dtype=NP_INT_TYPE)
        return np.empty((0, 2, 3), dtype=NP_INT_TYPE)

    # now remove any duplicat pairs in there
    if len(dimensions) == 2:
        # if 2D
        
        short_range_list = []

        if hardwall:
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops_hardwall.extract_SR_pairs_from_position_2D_hardwall(np.array(positions[i], dtype=NP_INT_TYPE), dimensions[0], dimensions[1]))
                
        else:
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops.extract_SR_pairs_from_position_2D(np.array(positions[i], dtype=NP_INT_TYPE), dimensions[0], dimensions[1]))

        envelope_pairs = np.concatenate(short_range_list)
        num_pairs = len(envelope_pairs)

        reshaped = np.reshape(envelope_pairs, (num_pairs, 4))

        if not deduplicate:
            return np.reshape(reshaped, (num_pairs, 2, 2))

        duplicate_free = _unique_rows(reshaped)

        return np.reshape(duplicate_free, (len(duplicate_free), 2,2))
    else:

        short_range_list = []

        if hardwall:

            for i in range(0, len(positions)):
                short_range_list.append(inner_loops_hardwall.extract_SR_pairs_from_position_3D_hardwall(np.array(positions[i], dtype=NP_INT_TYPE), dimensions[0], dimensions[1], dimensions[2]))
            
        else:
            
            for i in range(0, len(positions)):
                short_range_list.append(inner_loops.extract_SR_pairs_from_position_3D(np.array(positions[i], dtype=NP_INT_TYPE), dimensions[0], dimensions[1], dimensions[2]))


        envelope_pairs = np.concatenate(short_range_list)
        num_pairs = len(envelope_pairs)

        reshaped = np.reshape(envelope_pairs, (num_pairs, 6))

        if not deduplicate:
            return np.reshape(reshaped, (num_pairs, 2, 3))

        duplicate_free = _unique_rows(reshaped)

        return np.reshape(duplicate_free, (len(duplicate_free), 2,3))

#-----------------------------------------------------------------
#
#@profile
def build_all_envelope_pairs(positions, LR_binary_array, type_lattice, dimensions, hardwall=False, deduplicate=True):
    """
    Expects a LIST of positions and a numpy array of positions which engage in
    long-range interactions (or not) - 0 if not and 1 if yes.

    Returns a list of tuples, where each tuple is a pair of positions.

    The complete set of these positions represents the non-redudant set of long
    -range and short-range pairwise interactions associated with the positions
    defined in position

    Parameters
    ----------
    positions : list
        list of positions

    LR_binary_array : numpy array
        numpy array of positions which engage in long-range interactions (or not)
        - 0 if not and 1 if yes.

    type_lattice : numpy.ndarray
        The type grid used by the inner-loop routines to determine which
        long-range / super-long-range pairs are generated from each position.

    dimensions : list
        A list of length 2 or 3 giving the lattice dimensions.

    hardwall : bool, optional
        If True, the hardwall inner-loop variants are used so that pairs do not
        straddle the periodic boundary. Default is False.

    deduplicate : bool, optional
        If True (default) duplicate pairs are removed, which is required whenever the
        pairs are summed over (an energy evaluation would otherwise double count a
        repeated pair). Deduplication is a lexicographic sort and the dominant cost of
        this function, so callers that only need to know *which sites are touched* - the
        long-range connected-component search, which feeds the pairs into a set - should
        pass False.

    Returns
    -------
    tuple
        A 3-tuple ``(SR_pairs, LR_pairs, SLR_pairs)`` of numpy arrays giving,
        respectively, the duplicate-free short-range, long-range and
        super-long-range interaction pairs. Each array has shape (N, 2, 2) in 2D
        or (N, 2, 3) in 3D. If `positions` is empty, three empty arrays of the
        appropriate shape are returned.

    """

    if len(positions) == 0:
        if len(dimensions) == 2:
            empty = np.empty((0, 2, 2), dtype=NP_INT_TYPE)
        else:
            empty = np.empty((0, 2, 3), dtype=NP_INT_TYPE)
        return (empty.copy(), empty.copy(), empty.copy())


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
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_2D_hardwall(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1])
                else:
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops.extract_SR_and_LR_pairs_from_position_2D(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1])
                
                short_range_list.append(SR_tmp)
                
                if len(LR_tmp) > 0:
                    long_range_list.append(LR_tmp)
                    
                if len(SLR_tmp) > 0:
                    super_long_range_list.append(SLR_tmp)

        ####>>>> 3D
        else:
            for i in range(0, len(positions)):

                if hardwall:
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_3D_hardwall(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1], dimensions[2])
                else:
                    
                    (SR_tmp, LR_tmp, SLR_tmp)  = inner_loops.extract_SR_and_LR_pairs_from_position_3D(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_lattice, dimensions[0], dimensions[1], dimensions[2])
                    
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
        long_range_pairs = np.array([], dtype=NP_INT_TYPE)

    if len(super_long_range_list) > 0:
        super_long_range_pairs = np.concatenate(super_long_range_list)        
    else:
        super_long_range_pairs = np.array([], dtype=NP_INT_TYPE)

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

        if not deduplicate:
            return (np.reshape(reshaped_SR, (num_pairs_SR, 2, 2)),
                    np.reshape(reshaped_LR, (num_pairs_LR, 2, 2)),
                    np.reshape(reshaped_SLR, (num_pairs_SLR, 2, 2)))

        # short range witchcraft
        duplicate_free_SR = _unique_rows(reshaped_SR)

        # long range witchcraft
        duplicate_free_LR = _unique_rows(reshaped_LR)

        # super long range witchcraft
        duplicate_free_SLR = _unique_rows(reshaped_SLR)


        return (np.reshape(duplicate_free_SR, (len(duplicate_free_SR), 2,2)), np.reshape(duplicate_free_LR, (len(duplicate_free_LR), 2,2)), np.reshape(duplicate_free_SLR, (len(duplicate_free_SLR), 2,2)))

    else:
        # if 3D
                 
        reshaped_SR  = np.reshape(short_range_pairs, (num_pairs_SR, 6))
        reshaped_LR  = np.reshape(long_range_pairs,  (num_pairs_LR, 6))
        reshaped_SLR = np.reshape(super_long_range_pairs,  (num_pairs_SLR, 6))

        if not deduplicate:
            return (np.reshape(reshaped_SR, (num_pairs_SR, 2, 3)),
                    np.reshape(reshaped_LR, (num_pairs_LR, 2, 3)),
                    np.reshape(reshaped_SLR, (num_pairs_SLR, 2, 3)))

        # short range witchcraft
        duplicate_free_SR = _unique_rows(reshaped_SR)

        # long range witchcraft
        duplicate_free_LR = _unique_rows(reshaped_LR)

        # super long range witchcraft
        duplicate_free_SLR = _unique_rows(reshaped_SLR)

        return (np.reshape(duplicate_free_SR, (len(duplicate_free_SR), 2,3)), np.reshape(duplicate_free_LR, (len(duplicate_free_LR), 2,3)), np.reshape(duplicate_free_SLR, (len(duplicate_free_SLR), 2,3)))



#-----------------------------------------------------------------
#
def _grid_values_at(envelope_pairs, lattice_grid):
    """Grid values at every site of every envelope pair, as a flat list of ints.

    ``envelope_pairs`` is the ``(n_pairs, 2, n_dim)`` array returned by
    :func:`build_envelope_pairs`; this reads the occupying chainID at all
    ``2 * n_pairs`` sites in a single fancy-index instead of a Python loop calling
    :func:`get_gridvalue` twice per pair. The connected-component search ran that loop
    over a million times per cluster-analysis step.

    Parameters
    ----------
    envelope_pairs : numpy.ndarray
        ``(n_pairs, 2, n_dim)`` array of position pairs (may be empty).

    lattice_grid : numpy.ndarray
        The 2D or 3D lattice grid array.

    Returns
    -------
    list of int
        The grid value at every site referenced by ``envelope_pairs`` (0 = solvent).
    """
    if len(envelope_pairs) == 0:
        return []

    sites = np.asarray(envelope_pairs).reshape(-1, lattice_grid.ndim)

    if lattice_grid.ndim == 2:
        values = lattice_grid[sites[:, 0], sites[:, 1]]
    elif lattice_grid.ndim == 3:
        values = lattice_grid[sites[:, 0], sites[:, 1], sites[:, 2]]
    else:
        raise LatticeUtilsException(
            f"Unsupported lattice dimensionality in _grid_values_at: {lattice_grid.ndim}")

    return values.tolist()


#-----------------------------------------------------------------
#
def get_all_chains_in_connected_component(chainID, lattice_grid, chainDict, threshold=None, useChains=True, hardwall=False):
    """
    Function which given a chainID, a dictionary of chain-to-position mappings, and a lattice grid
    will return the set of chains in the connected component containing chainID. Note a connected 
    component is a *heterotypic* structure - i.e. we are looking for a connected component made up
    of *any* chains, not a single type of chain.

    Parameters
    ----------

    chainID : int
        The chainID of the chain we initially are asking about

    lattice_grid : 2D or 3D np.array
        Standard lattice grid

    chainDict : dictionary mapping chainIDs to a list of positions or to a chain object
        Dictionary containing a mapping of each chainID to either a list of positions associated
        with that chain, or the Chain object associated with that chainID

    threshold : int
        The max size of the connected component we are looking for. If None, this is ignored,
        but if set and we generate a connected component larger than this, we will raise a
        ClusterSizeThresholdException. Enables us to avoid situations where we've moving massive
        giant clusters around which may not be efficient if 90% of the chains are in the cluster.


    useChains : Bool
        Boolean flag which defines if the chainDict is a true dictionary mapping chainID
        to a set of positions, or in fact a dictionary of Chain objects (which contain
        positions which must be accessed using the .get_ordered_positions()). This isn't
        so much a feature as the fact that we want this function to be able to accept
        two different types of chain information (dictionary of lists of positions or
        dictionary of chain objects)

    hardwall : Bool
        Boolean flag which defines if we are using a hardwall potential or not. If we are, we
        will use a different method to calculate the connected component. If we are not, we will
        use a more efficient method which uses a union-find algorithm to calculate the connected
    
    Return
    ---------

    list  

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
        # deduplicate=False: the pairs go straight into a set below, so paying for the
        # lexicographic dedupe sort would be wasted work
        envelope_pairs = build_envelope_pairs(positions, dimensions, hardwall=hardwall, deduplicate=False)

        # look up which chain occupies each site of each pair. This is one fancy-index
        # into the grid rather than two Python-level get_gridvalue() calls per pair -
        # that loop ran to well over a million calls per cluster analysis and was the
        # dominant cost of the connected-component search.
        new_chains.update(_grid_values_at(envelope_pairs, lattice_grid))

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

    """
    Function which given a chainID, a dictionary of chain-to-position 
    mappings, and a lattice grid will return the set of chains in the
    connected component where connectivity is defined in terms of 
    long-range interactions.  Note a connected component is a
    *heterotypic* structure - i.e. we are looking for a connected 
    component made up of *any* chains, not a single type of chain.

    Parameters
    ----------

    chainID : int
        The chainID of the chain we initially are asking about

    latticeObject : Lattice object
        Lattice object containing the lattice grid, the type grid, 
        and the chain dictionary

    hardwall : Bool
        Boolean flag which defines if we are using a hardwall potential 
        or not. 

    Return
    ---------
    list
        A list of chainIDs associated with the chains in the connected
        component which contans the chain defined by $chainID



    """

    
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
        # deduplicate=False: as above, these only feed a set of chainIDs
        (SR_pairs, LR_pairs, SLR_pairs) = build_all_envelope_pairs(positions, LR_binary_array, type_grid, dimensions, hardwall, deduplicate=False)

        envelope_pairs = np.concatenate((SR_pairs, LR_pairs))    
                
        # look up which chain occupies each site of each pair (see the note in
        # get_all_chains_in_connected_component - one fancy-index rather than two
        # Python-level grid lookups per pair)
        new_chains.update(_grid_values_at(envelope_pairs, lattice_grid))

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

    Parameters
    ----------
    positions : list
        List of positions to calculate the center of mass from.

    dimensions : list
        List of dimensions of the box.

    on_lattice : bool
        If True, the center of mass is calculated on the lattice.
        If False, the center of mass is calculated in Euclidean space.

    Returns
    -------
    list
        The center of mass of the positions (2D or 3D list depending on
        if the input positions are 2D or 3D).

    """

    if len(positions) == 0:
        raise LatticeUtilsException("Cannot compute center of mass: positions list is empty")

    n_dim = len(dimensions)

    # Circular (periodic-aware) mean of the positions, computed per axis and
    # vectorized over all beads at once (this used to be a Python loop calling
    # np.cos/np.sin scalar-by-scalar per bead per axis - a hot path in the cluster
    # analysis via the single-image seed and the radial-density COM). Each
    # coordinate is mapped to an angle on a circle whose circumference is the box
    # size on that axis; averaging the unit vectors and taking the argument gives
    # the mean position that respects the wrap-around.
    pos = np.asarray(positions, dtype=np.float64)          # (N, n_dim)
    dims = np.asarray(dimensions, dtype=np.float64)         # (n_dim,)

    angles = (pos / dims) * (2.0 * np.pi)
    mean_cos = np.cos(angles).mean(axis=0)
    mean_sin = np.sin(angles).mean(axis=0)

    real = dims * (np.arctan2(-mean_sin, -mean_cos) + np.pi) / (2.0 * np.pi)

    if on_lattice:
        coords = [int(round(float(v))) for v in real]
    else:
        coords = [float(v) for v in real]

    return pbc_convert(coords, dimensions)


    
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
    """

    Delete a residue at a given position in the lattice. This function
    will raise an exception if the position is already occupied by a residue
    from a different chain. This is the safe version of the function. If you
    want to overwrite the residue, set safe=False.

    Parameters
    ----------
    position : list
        Position to delete the residue from.

    lattice : np.array
        Lattice to delete the residue from.

    chainID : int
        Chain ID of the residue to delete. If None, the residue will be deleted
        regardless of the chain ID.

    Returns
    -------
    None

    Raises
    ------
    ResidueAugmentException
        Raised (only when `chainID` is provided) if the residue currently at
        `position` does not belong to the expected chain.

    """

    if chainID is not None:
        ## Safe version

        # get id of residue to delete
        todel = get_gridvalue(position, lattice)

        # if mismatch, raise exception
        if not todel == chainID:
            raise ResidueAugmentException(
                f'Trying to delete a residue at position {str(position)} - expected chainID {chainID}, but got chainID {todel}'
            )
        else:
            set_gridvalue(position, 0.0, lattice)

    else:
        ## No checks version...
        set_gridvalue(position, 0.0, lattice)



#-----------------------------------------------------------------
#
def insert_residue(position, lattice, chainID, safe=True):
    """
    Insert a residue at a given position in the lattice. This function
    will raise an exception if the position is already occupied by a residue
    from a different chain. This is the safe version of the function. If you
    want to overwrite the residue, set safe=False.

    Parameters
    ----------
    position : list
        Position to insert the residue

    lattice : np.array
        Lattice to insert the residue into

    chainID : int
        Chain ID to insert the residue for

    safe : bool
        If True, will raise an exception if the position is 
        already occupied. If False, will overwrite the residue.

    Returns
    -------
    None
    

    """

    if safe:
        insert_location = get_gridvalue(position, lattice)

        # todo - this probably should be an int comparison - check and fix at somepoint...
        if not insert_location == 0.0:
            raise ResidueAugmentException(f'Trying to insert a residue for chain at position {str(position)} in {chainID} - site was occupied by residue from chain {insert_location}! This is a bug - please report.')
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
    low-level function that performs single point rotation. This should generally
    not be called but instead the wrapper functions rotate_positions_3D or 
    rotate_positions_2D should be used.

    Parameters
    ----------
    positions : list
        List of positions to rotate

    rotation_matrix : np.array
        Rotation matrix to apply

    Returns
    -------
    list
        List of rotated positions

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

    Parameters
    ----------
    positions : list
        List of positions to rotate

    dimension : str
        Dimension to rotate around. Must be one of 'x', 'y', or 'z'.

    degrees : int
        Degrees to rotate by. Must be one of 90, 180, or 270.

    Returns
    -------
    list
        List of rotated positions

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

    Parameters
    -------------
    positions : list
        A list of 2D positions to be rotated.

    degrees : int
        The number of degrees to rotate the positions by. Must be one 
        of 90, 180, or 270.

    Returns
    -------------
    rotated_positions : list
        A list of 2D positions that have been rotated by the specified 
        number of degrees.

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

    Returns
    -------
    None
        No return value, but a new PDB file is initialized on disk.

    """

    pdb_utils.initialize_pdb_file(dimensions, spacing, filename)
    

#-----------------------------------------------------------------
#
def write_lattice_to_pdb(latticeObject, spacing, filename='lattice.pdb', write_connect=False, autocenter=False, unwrap=False):
    """
    Wrapper function that dumps the current Lattice object to a PDB file

    Parameters
    -------------
    latticeObject : lattice.Lattice
        Current lattice object

    spacing : float
        Lattice-to-realspace spacing in angstroms. 

    filename : str
        Filename to write to (default is lattice.pdb)

    write_connect : bool
        Flag to write connect information to the PDB file (default is False)

    autocenter : bool
        Flag to center the lattice in the PDB file (default is False). Note 
        that this correctly is dealth with in build_pdb_file - if more than
        one chain this is ignored.

    Returns
    ------------
    None

    """
    pdb_utils.build_pdb_file(latticeObject, spacing, filename, write_connect=write_connect, autocenter=autocenter, unwrap=unwrap)



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
def start_xtc_file(lattice, spacing, pdb_filename='START.pdb', xtc_filename='traj.xtc', unwrap=False):
    """
    Function that initializes a new .xtc file. This deletes an existing XTC file 
    of the same name to avoid any issues.

    Parameters
    ------------
    lattice : lattice.Lattice
        Current Lattice object

    spacing : float
        Lattice-to-realspace spacing in angstroms, used when writing the
        corresponding PDB file.

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
        IO_utils.status_message(f"Deleted existing XTC file [{xtc_filename}]", 'startup')

    except OSError:
        pass

    # first build the PDB file
    open_pdb_file(lattice.dimensions, spacing, filename=pdb_filename)
    write_lattice_to_pdb(lattice, spacing, filename=pdb_filename, write_connect=True, unwrap=unwrap)
    finish_pdb_file(pdb_filename)

    # next read the PDBFILE, and save as an xtcfile
    traj = md.load(pdb_filename)
    traj.save_xtc(xtc_filename)


#-----------------------------------------------------------------
#
def _lattice_frame_xyz_and_box(lattice, spacing, autocenter=False, unwrap=False):
    """
    Build the ``(1, n_atoms, 3)`` coordinate array (in nm) and the orthorhombic box
    (in nm) for one trajectory frame from the current lattice.

    Positions are gathered per chain via ``get_output_positions`` (honouring the
    ``autocenter`` / ``unwrap`` conventions); 2D systems are padded with a zero z.

    Returns
    -------
    tuple
        ``(xyz, box)`` where ``xyz`` is float32 shape ``(1, n_atoms, 3)`` and
        ``box`` is float32 shape ``(1, 3, 3)`` (diagonal box vectors, nm).
    """
    # autocenter is only meaningful for a single chain
    if autocenter and len(lattice.chains) > 1:
        autocenter = False

    is_3d = len(lattice.dimensions) == 3
    cvals = []
    for chainID in lattice.chains:
        positions = lattice.chains[chainID].get_output_positions(autocenter=autocenter, unwrap=unwrap)
        if is_3d:
            cvals.extend(positions)
        else:
            cur = np.array(positions)
            cur = np.hstack((cur, np.zeros((len(cur), 1), dtype=cur.dtype)))
            cvals.extend(list(cur))

    xyz = np.array([cvals], dtype=np.float32) * spacing * 0.1

    dims = lattice.dimensions
    lz = (dims[2] if is_3d else 1)
    box = np.array([[[dims[0] * spacing * 0.1, 0.0, 0.0],
                     [0.0, dims[1] * spacing * 0.1, 0.0],
                     [0.0, 0.0, lz * spacing * 0.1]]], dtype=np.float32)
    return xyz, box


#-----------------------------------------------------------------
#
def open_xtc_writer(lattice, spacing, pdb_filename='START.pdb', xtc_filename='traj.xtc', autocenter=False, unwrap=False):
    """
    Write the topology PDB, open a persistent XTC writer, write the first frame, and
    return the open writer handle.

    This is the efficient replacement for the previous per-frame
    ``append_to_xtc_file_non_redundant`` approach, which re-loaded the entire growing
    trajectory from disk and re-saved it on EVERY frame - O(frames^2) in both wall
    time and disk I/O. Here a single ``mdtraj`` XTC file handle is kept open for the
    whole run and each frame is appended with :func:`write_xtc_frame` in O(1); the
    handle is closed with :func:`close_xtc_writer`.

    Parameters
    ----------
    lattice : lattice.Lattice
        Current lattice.
    spacing : float
        Lattice-to-realspace spacing (angstroms).
    pdb_filename : str
        Topology PDB filename to (re)write.
    xtc_filename : str
        Trajectory filename to create.
    autocenter : bool
        Single-chain autocentring (see build_pdb_file). Default False.
    unwrap : bool
        Make chains whole across PBC before writing (TRAJECTORY_PBC_UNWRAP). Default False.

    Returns
    -------
    mdtraj.formats.XTCTrajectoryFile
        The open writer handle (write more frames with write_xtc_frame, then close
        with close_xtc_writer).
    """
    # (re)write the topology PDB
    open_pdb_file(lattice.dimensions, spacing, filename=pdb_filename)
    write_lattice_to_pdb(lattice, spacing, filename=pdb_filename, write_connect=True, unwrap=unwrap)
    finish_pdb_file(pdb_filename)

    # start a fresh XTC file and write the first frame
    if os.path.exists(xtc_filename):
        os.remove(xtc_filename)
    writer = md.formats.XTCTrajectoryFile(xtc_filename, 'w')
    xyz, box = _lattice_frame_xyz_and_box(lattice, spacing, autocenter=autocenter, unwrap=unwrap)
    writer.write(xyz, box=box)
    return writer


#-----------------------------------------------------------------
#
def write_xtc_frame(writer, lattice, spacing, autocenter=False, unwrap=False):
    """
    Append one frame from the current lattice to an open XTC writer (O(1), no reload).

    Parameters
    ----------
    writer : mdtraj.formats.XTCTrajectoryFile
        Open writer handle from :func:`open_xtc_writer`.
    lattice : lattice.Lattice
        Current lattice.
    spacing : float
        Lattice-to-realspace spacing (angstroms).
    autocenter : bool
        Single-chain autocentring. Default False.
    unwrap : bool
        Make chains whole across PBC before writing. Default False.

    Returns
    -------
    None
    """
    xyz, box = _lattice_frame_xyz_and_box(lattice, spacing, autocenter=autocenter, unwrap=unwrap)
    writer.write(xyz, box=box)


#-----------------------------------------------------------------
#
def close_xtc_writer(writer):
    """
    Close an open XTC writer (flushing the file). Safe to call with ``None``.

    Parameters
    ----------
    writer : mdtraj.formats.XTCTrajectoryFile or None
        The writer handle to close.

    Returns
    -------
    None
    """
    if writer is not None:
        writer.close()



#-----------------------------------------------------------------
#
def append_to_xtc_file(lattice, spacing, xtc_filename='traj.xtc', autocenter=False):

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

    autocenter : bool
        Flag which, if set to True and there's a single chain will center the protein in the box. 
        This is useful for visualization purposes but does mean any translational diffusion will
        be lost. Default = False


    Returns
    -----------
    None
        No return by the existing XTC file is extended by one frame

    """

    # first build the PDB file    
    open_pdb_file(lattice.dimensions, spacing, filename='frame.pdb')
    write_lattice_to_pdb(lattice, spacing, filename='frame.pdb', autocenter=autocenter)
                         
    finish_pdb_file('frame.pdb')

    xtc_traj = md.load(xtc_filename, top='frame.pdb')
    
    pdb_frame = md.load('frame.pdb')

    new = xtc_traj.join(pdb_frame)

    new.save(xtc_filename)


    
#-----------------------------------------------------------------
#
def append_to_xtc_file_non_redundant(lattice,
                                     spacing,
                                     pdb_filename='START.pdb',
                                     xtc_filename='traj.xtc',
                                     autocenter=False,
                                     unwrap=False):

    """
    Low level function that adds a current lattice to an existing XTC file.
    This is different than 'append_to_xtc_file' in that it does not make a frame.pdb
    object each time it needs to save the frame. 
    Uses the START.pdb file as the topology. 

    Parameters
    -----------
    lattice : lattice.Lattice
        A lattice object

    spacing : float
        Lattice-to-realspace spacing in angstroms. 

    pdb_filename : str
        Topology filename to read from and extend

    xtc_filename : str
        Trajectory filename to read from and extend

    autocenter : bool
        Flag which, if set to True and there's a single chain will center the protein in the box. 
        This is useful for visualization purposes but does mean any translational diffusion will
        be lost. Default = False


    Returns
    -----------
    None
        No return, but the existing XTC file is extended by one frame and then saved
        to disk.

    """
    # overide autocenter if more than 1 chain
    if autocenter and len(lattice.chains)>1:
        autocenter = False

        
    # load the xtc trajectory that is already started. NB: raise rather than exit() - a
    # bare `except:` here also swallowed KeyboardInterrupt/SystemExit, and exit() is the
    # `site` builtin (absent under `python -S`) which tears down the caller's process
    # instead of letting them handle the failure.
    try:
        xtc_traj = md.load(xtc_filename, top=pdb_filename)
    except Exception as e:
        raise LatticeUtilsException(
            f"Error loading xtc file {xtc_filename} with topology {pdb_filename}: {e}") from e
            
    # coordinate vals = cvals... now we need to get the positions of the chains in the sim. 
    cvals = []

    # if we're in 3D...
    if len(lattice.dimensions) == 3:

        # iterate over chains.
        for chain in lattice.chains:

            # extend cvals by the coord vals for this chain
            cvals.extend(lattice.chains[chain].get_output_positions(autocenter=autocenter, unwrap=unwrap))


    # if we're in 2D...
    else:

        for chain in lattice.chains:

            # extend cvals by the coord vals for this chain
            curchain = np.array(lattice.chains[chain].get_output_positions(autocenter=autocenter, unwrap=unwrap))
            
            # if we have a 2D array, we need to add a third coordinate.
            # to do this we can just hstack zeros on to the cvals array
            zeros = np.zeros((len(curchain),1),dtype=np.int8)
            curchain = np.hstack((curchain,zeros))
            
            cvals.extend(list(curchain))

    # make the newdims an array times spacing and account for angstoms vs nanometers
    newdims = np.array([cvals])*spacing*0.1    
    
    # make frame trajectory using xyz values times spacing divided by 10 to account for angstroms vs. nm.
    current_frame_traj = md.Trajectory(newdims,
                                       xtc_traj.topology,
                                       time=xtc_traj.time[-1]+1,
                                       unitcell_lengths=xtc_traj.unitcell_lengths[0],
                                       unitcell_angles=xtc_traj.unitcell_angles[0])
    
    # make a new traj by adding the traj for the current frame to the xtc_traj we loaded in and save iteratively over.
    new_traj = xtc_traj.join(current_frame_traj)
    
    # save the new traj as xtc_filename.
    new_traj.save(xtc_filename)


#-----------------------------------------------------------------
#
class TrajectoryAccumulator:
    """In-memory buffer of trajectory frames, materialised into one Trajectory at the end.

    This is what ``SAVE_AT_END : True`` accumulates into. It exists purely for
    performance: the previous approach called ``master_traj.join(frame)`` once per
    saved frame, and because ``join`` returns a NEW Trajectory holding a copy of every
    frame seen so far, writing ``n`` frames copied ``1 + 2 + ... + n`` frames' worth of
    coordinates - quadratic in trajectory length, in both time and memory traffic. (The
    default incremental path was moved onto a persistent XTC writer for the same reason;
    the SAVE_AT_END path kept the quadratic behaviour.)

    Here each frame is appended to a list and the single join happens once, in
    :meth:`to_trajectory`, so the total work is linear. Peak memory is unchanged -
    buffering the whole trajectory is the point of SAVE_AT_END - but the transient
    copies are gone.

    Frames carry the same topology, unit cell and one-per-frame time stamps as before.
    """

    __slots__ = ("_base", "_frames", "_last_time")

    def __init__(self, base):
        """
        Parameters
        ----------
        base : mdtraj.Trajectory
            The topology frame (loaded from the START.pdb), which becomes the first
            frame of the finished trajectory.
        """
        self._base = base
        self._frames = []
        self._last_time = float(base.time[-1])

    @property
    def topology(self):
        return self._base.topology

    def append_frame(self, xyz):
        """Buffer one frame's ``(1, n_atoms, 3)`` coordinates (nm)."""
        self._last_time = self._last_time + 1
        self._frames.append(md.Trajectory(xyz,
                                          self._base.topology,
                                          time=self._last_time,
                                          unitcell_lengths=self._base.unitcell_lengths[0],
                                          unitcell_angles=self._base.unitcell_angles[0]))
        return self

    def to_trajectory(self):
        """Materialise the buffered frames into a single ``mdtraj.Trajectory``."""
        if len(self._frames) == 0:
            return self._base
        return self._base.join(self._frames)

    def __len__(self):
        return 1 + len(self._frames)


#-----------------------------------------------------------------
#
def update_master_traj(lattice, spacing, master_traj, pdb_filename, autocenter=False, unwrap=False):

    """
    Low level function that adds a current lattice to an existing XTC file.
    This is different than 'append_to_xtc_file' in that it does not read in an
    existing XTC file but instead appends a frame to the passed master trajectory
    object.

    If the master_traj object has not yet been initialized, this will also read
    in the pdb_filename and initialize the master_traj object using that as a 
    topology file. This means PDB initialization needs to happen BEFORE this
    function is called! This is a deliberate design choice to allow for the
    master_traj object to be passed between functions and have the topology
    file be read in only once.

    Parameters
    -----------
    lattice : lattice.Lattice
        A lattice object

    spacing : float
        Lattice-to-realspace spacing in angstroms. 

    master_traj : mdtraj.Trajectory 
        master trajectory we will build throught the sim. 

    pdb_file_name : current_pdb_filename
        the current_pdb_filename from simulation.py

    autocenter : bool
        Flag which, if set to True and there's a single chain will center the protein in the box. 
        This is useful for visualization purposes but does mean any translational diffusion will
        be lost. Default = False


    Returns
    -----------
    TrajectoryAccumulator
        The updated accumulator (see :class:`TrajectoryAccumulator`) - pass it back in
        on the next call, and hand it to :func:`save_out_sim` at the end. Frames are
        buffered and joined once rather than re-joined per frame, which is what makes
        this linear rather than quadratic in the number of frames.

    """
    # coordinate vals = cvals... now we need to get the positions of the chains in the sim.
    cvals = []

    # overide autocenter if more than 1 chain
    if autocenter and len(lattice.chains)>1:
        autocenter = False


    # if 3D...
    if len(lattice.dimensions) == 3:
        
        # iterate over chains.
        for chain in lattice.chains:

            # extend cvals by the coord vals for this chain
            cvals.extend(lattice.chains[chain].get_output_positions(autocenter=autocenter, unwrap=unwrap))

    # if 2D....
    else:
        for chain in lattice.chains:

            # extend cvals by the coord vals for this chain
            curchain = np.array(lattice.chains[chain].get_output_positions(autocenter=autocenter, unwrap=unwrap))
            
            # if we have a 2D array, we need to add a third coordinate.
            # to do this we can just hstack zeros on to the cvals array
            zeros = np.zeros((len(curchain),1),dtype=np.int8)
            curchain = np.hstack((curchain,zeros))
            cvals.extend(list(curchain))

    # make the newdims an array times spacing and account for angstoms vs nanometers
    newdims = np.array([cvals])*spacing*0.1  
    
    # if the master_traj is not yet initialized, use the pdb file name as start point.
    # (`is None`, not `== None`: mdtraj Trajectory defines __eq__, so `== None` is not
    # guaranteed to be the identity test that is meant here.) As above, raise rather
    # than swallowing everything with a bare except and calling exit().
    if master_traj is None:

        try:
            base = md.load(pdb_filename, top=pdb_filename)
        except Exception as e:
            raise LatticeUtilsException(
                f'Could not load pdb file: {pdb_filename}: {e}') from e

        master_traj = TrajectoryAccumulator(base)

    # tolerate being handed a plain Trajectory (older callers / tests)
    elif not isinstance(master_traj, TrajectoryAccumulator):
        master_traj = TrajectoryAccumulator(master_traj)

    # buffer this frame. The frames are joined once, in save_out_sim - joining here (as
    # this used to) copies the entire accumulated trajectory on every single frame.
    return master_traj.append_frame(newdims)


#-----------------------------------------------------------------
#
def save_out_sim(master_traj, xtc_filename):
    """
    Save out the master trajectory. 

    Parameters
    ------------
    master_traj : TrajectoryAccumulator or mdtraj.Trajectory
        The trajectory built up through the sim. A
        :class:`TrajectoryAccumulator` (what :func:`update_master_traj` returns) is
        materialised into a single Trajectory here; a plain Trajectory is saved as is.

    xtc_filename : str
        Filename to write to disk.

    Returns
    ------------
    None
        No return, but the existing XTC file is saved to disk

    """
    if isinstance(master_traj, TrajectoryAccumulator):
        master_traj = master_traj.to_trajectory()

    # save the new traj as xtc_filename.
    master_traj.save(xtc_filename)
    

#######################################################################################
##                                                                                   ##
##                           SANITY CHECKING FUNCTIONS                               ##
##                                                                                   ##
#######################################################################################


#-----------------------------------------------------------------
#
def check_chain_connectivity(chainID, chain_positions, dimensions, verbose=True):
    """
    Debugging function which ensures that a set of positions
    correspond to a valid, connected chain. Useful for
    debugging new moves, though not designed for performance
    during real simulations. If verbose is set to True, will
    print out a "CONNECTIVITY FINE" message assuming the 
    function completes without issue.

    Parameters
    ------------

    chainID : int
        Chain ID number

    chain_positions : list of lists
        List of lists of positions for the chain

    dimensions : list
        List of dimensions for the lattice

    verbose : bool
        Flag to print out debug information. Default = True

    Returns
    ------------
    None
        No return, but will raise an error if the chain is not connected


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
    if verbose:
        print("CONNECTIVITY FINE")


#-----------------------------------------------------------------
#
def check_all_chain_connectivity(list_of_chain_objects, dimensions, verbose=True):
    """
    Debugging function that takes a LATTICE.chains and LATTICE.dimensions
    pair from the simulation object to check the chain connectivity over
    all chains in the simulation.

    Parameters
    ------------
    list_of_chain_objects : dict
        Dictionary of chain objects

    dimensions : list
        List of dimensions for the lattice

    verbose : bool
        Flag to print out debug information. Default = True


    Returns
    ------------
    None
        No return, but will raise an error if any chain is not connected


    """
    
    for chainID in list_of_chain_objects:
        check_chain_connectivity(chainID, list_of_chain_objects[chainID].get_ordered_positions(), dimensions, verbose=verbose)




        
