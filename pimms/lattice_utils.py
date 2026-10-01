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

import collections
import glob
import itertools
import random
import copy
import math
import os
import shutil
import numpy as np

import mdtraj as md

from .latticeExceptions import ChainInsertionFailure, ChainDeletionFailure, ResidueAugmentException, ChainConnectivityError, ClusterSizeThresholdException, LatticeUtilsException, RotationException

#from .pdb_utils import build_pdb_file, finalize_pdb_file, initialize_pdb_file
from . import pdb_utils

from . import hyperloop
from . import inner_loops
from . import inner_loops_hardwall

from . import lattice_analysis_utils
from . import IO_utils

from . import CONFIG
from . CONFIG import NP_INT_TYPE

#from CONFIG import * # note there are things from CONFIG being used...

#-----------------------------------------------------------------
#
# Number of beads whose envelope pairs build_all_envelope_pairs holds as
# separate per-bead arrays before folding them into one array. Large enough
# that a single-chain or cluster move is one chunk (so nothing changes for
# them), small enough that the padded buffers pinned by one chunk's LR/SLR
# views stay below ~10 MB.
_ENVELOPE_PAIR_CHUNK = 1024

# The unit cell (float32, shape (1, 3, 3), nm) each buffered trajectory was
# started with, keyed by the absolute path of its topology PDB.
_STARTED_TRAJECTORY_BOX = {}

#-----------------------------------------------------------------
#
# AUTOCENTER is only defined for a single chain. With more than one chain it has
# always been switched off without a word; we now say so, once per process.
_AUTOCENTER_IGNORED_WARNED = False


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

    Raises
    ---------
    LatticeUtilsException
        If the two sites have different dimensionality, or if that
        dimensionality is neither 2 nor 3.

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
    Function to calculate the real distance between two positions. This is a thin
    wrapper around lattice_analysis_utils.get_inter_position_distance(), so the
    minimum-image (PBC) convention is always applied.


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
        being studied, that reflects the lattice dimensions. These are the
        periods used for the minimum-image correction.


    Returns
    ---------
    float
        Returns a value that reflects the minimum-image Euclidian distance
        between two positions on the lattice.


    """
    return lattice_analysis_utils.get_inter_position_distance(posA, posB, dimensions)
       

#----------------------------------------------------------------
#
def get_dimensions(lattice_grid):
    """
    Function that returns the dimensions associated with the lattice

    Parameters
    ----------
    lattice_grid : numpy.ndarray
        A 2D or 3D lattice grid numpy array (either the occupancy grid or the
        type grid - only its shape is used).


    Returns
    ---------
    tuple
        The array shape: a tuple of length 2 or 3, depending on the
        dimensionality of the system, where the value at each position reflects
        the size of the lattice in that dimension.


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
    set. This is arbitrary, but we don't need to futz with them both!!

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


def clamp_positions_to_box(positions, dimensions):
    """
    Rigidly translate a set of positions by the smallest whole-site shift that
    brings every one of them inside ``[0, L)`` on every axis.

    This is what keeps ``AUTOCENTER`` honest in a hardwall box. Centring puts
    the chain's centre of mass on the box centre, and for an asymmetric chain
    (most of the beads near one wall, a tail reaching to the other) that shift
    carries the tail through the wall. Under periodic boundaries that is just
    another image of the bead; under ``HARDWALL`` there is no image, so a
    coordinate outside the box describes a state the simulation can never
    visit. We therefore clamp the centring shift: the chain is centred as far
    as the walls allow and no further. The translation is rigid (every bead
    moves by the same vector), so the conformation is untouched.

    Parameters
    ----------
    positions : list
        A list of lists, where each inner list is a (single-image) position on
        the lattice.

    dimensions : list
        A list of length 2 or 3 giving the lattice dimensions.

    Returns
    -------
    list
        A list of lists holding the translated positions. If the positions
        already lie inside the box they are returned unchanged (as a new
        list). An axis along which the positions span more sites than the box
        holds cannot be made to fit, and is left with its lowest bead on
        site 0.

    """
    n_dim = len(dimensions)
    shift = []
    for idx in range(0, n_dim):
        lowest = min(pos[idx] for pos in positions)
        highest = max(pos[idx] for pos in positions)
        if lowest < 0:
            shift.append(-lowest)
        elif highest > dimensions[idx] - 1:
            # never push the low end through the opposite wall
            shift.append(max(dimensions[idx] - 1 - highest, -lowest))
        else:
            shift.append(0)

    return [[pos[idx] + shift[idx] for idx in range(0, n_dim)] for pos in positions]


def center_positions(positions, dimensions, hardwall=False):
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

    hardwall : bool, optional
        If True the centring shift is clamped so that no bead is carried
        outside ``[0, L)`` (see :func:`clamp_positions_to_box`): in a hardwall
        box a coordinate beyond the wall is not a periodic image of anything.
        Default is False, which centres the centre of mass exactly and lets
        the ends of an asymmetric chain fall outside the box (a valid periodic
        image).

    Returns
    ---------
    list
        A list of lists, where each inner list is a position on the lattice,
        that has been centered in the box.

    """

    n_dim = len(dimensions)

    # The positions handed in here are SINGLE-IMAGE (already unwrapped by
    # convert_chain_to_single_image, so coordinates may exceed the box). Their centre
    # is therefore the plain arithmetic mean - NOT the periodic (circular-mean) COM,
    # which re-wraps the unwrapped coordinates and, for a chain whose true centre sits
    # at or beyond a box face, lands at ~0. The offset was then ~+L/2 and the
    # "centred" chain was placed entirely OUTSIDE the box - i.e. AUTOCENTER broke for
    # exactly the boundary-straddling chains it exists to tidy up. Rounded to the
    # lattice to match the on-lattice convention of the old path for non-straddling
    # chains.
    com = [int(round(float(v))) for v in np.mean(np.asarray(positions, dtype=np.float64), axis=0)]

    # integer box centre: dimensions/2 is a half-integer on an odd axis and would
    # put every bead of the centred chain half a site off the lattice in
    # START.pdb and every trajectory frame
    offset = []
    for idx in range(0, n_dim):
        offset.append(dimensions[idx] // 2 - com[idx])
        
    new_positions = []

    for pos in positions:
        new_pos = []
        for idx in range(0, n_dim):
            new_pos.append(pos[idx] + offset[idx])
        new_positions.append(new_pos)

    if hardwall:
        return clamp_positions_to_box(new_positions, dimensions)

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

    This underpins every intra-chain observable (``Chain.get_analysis_positions``)
    as well as trajectory unwrapping (``TRAJECTORY_PBC_UNWRAP``), and
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

    extent_range : int, optional
        Half-width of the neighbourhood to enumerate around the position.
        Default is 1, i.e. the 3x3x3 block of immediately adjacent sites.

    Returns
    -------
    numpy.ndarray
        A ``((2 * extent_range + 1)**3, 3)`` int32 array of PBC-wrapped site
        coordinates. Note this INCLUDES the central position itself.

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

    extent_range : int, optional
        Half-width of the neighbourhood to enumerate around the position.
        Default is 1, i.e. the 3x3 block of immediately adjacent sites.

    Returns
    -------
    numpy.ndarray
        A ``((2 * extent_range + 1)**2, 2)`` int32 array of PBC-wrapped site
        coordinates. Note this INCLUDES the central position itself.

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

    positions_list : list or numpy.ndarray
        A list of 2D or 3D positions (each a list of ints), or an ``(N, n_dim)``
        integer array, to be compared against the target.

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
    of the whole lattice until an empty site is found. The lattice is only
    scanned for "is there an empty site at all?" after every 100 failed draws,
    so the common case costs a few draws rather than a pass over the whole box.
    When `adjacentTo` is set,
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
        Either way the coordinates are Python ints (these positions become
        ``Chain.positions``, which must not hold numpy scalars).

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
            # (and the coordinates to Python ints: the adjacent sites are int32
            # rows, and this position goes straight into Chain.positions)
            return ([int(coordinate) for coordinate in empty_list[random.randint(0,len(empty_list)-1)]], True)

    # if we're literally just looking for an empty site anywhere on the lattice
    else:
        # A full lattice is reported from inside the loop (after every 100 failed
        # draws) rather than up front. The whole-grid scan is O(volume) and this
        # branch runs once per chain at start-up, so scanning before the first
        # draw made chain placement O(chains x volume) - 10.8 s of an 11 s
        # start-up for 10^4 chains in a 200^3 box. A draw that lands on an empty
        # site returns exactly as before, so no draw is added, removed or
        # reordered for any lattice that has an empty site; a full lattice now
        # raises the same exception after 99 draws instead of none.
        empty = False
        count = 0
        max_attempts = max(1000, int(np.prod(dimensions)) * 10)
        while not empty:
            count=count+1

            if count % 100 == 0:
                if not np.any(lattice_grid == 0):
                    raise LatticeUtilsException("Unable to find empty lattice site: lattice appears fully occupied")
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
                    # Python ints, since this position goes into Chain.positions
                    position = [int(x), int(y)]
                    empty=True

            # 3D
            if len(dimensions) == 3:
                x = NP_INT_TYPE(random.randint(0, dimensions[0]-1))
                y = NP_INT_TYPE(random.randint(0, dimensions[1]-1))
                z = NP_INT_TYPE(random.randint(0, dimensions[2]-1))

                if get_gridvalue([x,y,z], lattice_grid) == 0:
                    position = [int(x), int(y), int(z)]
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
            # get_empty_site raises LatticeUtilsException on a full lattice; callers
            # (Chain.__init__) catch ChainInsertionFailure to produce the "lattice is
            # overcrowded" message, and the full-lattice case is exactly that case.
            try:
                position = get_empty_site(lattice_grid)
            except LatticeUtilsException as e:
                raise ChainInsertionFailure(str(e)) from e


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
            # not -1 because we assign the first lattice site outside
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
    int
        The grid value stored at `position`, as a numpy integer scalar of the
        grid's dtype. 0 denotes an empty/solvent site; otherwise this is the
        occupying chainID, or the bead type code if a type grid was passed.

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
    int
        The grid value stored at `position`, as a numpy integer scalar of the
        grid's dtype. 0 denotes an empty/solvent site; otherwise this is the
        occupying chainID, or the bead type code if a type grid was passed.

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
    int
        The grid value stored at `position`, returned by the compiled routine as
        a plain Python int. 0 denotes an empty/solvent site; otherwise this is
        the occupying chainID, or the bead type code if a type grid was passed.

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

    Parameters
    ----------
    rows : numpy.ndarray
        A 2D integer array whose rows are the items to deduplicate - in practice
        the ``(n_pairs, 4)`` or ``(n_pairs, 6)`` flattened pair array built by
        the envelope routines. An empty input is returned unchanged.

    Returns
    -------
    numpy.ndarray
        The subset of ``rows`` with duplicate rows removed. Row order follows
        the sort used to deduplicate and is NOT the input order, so callers must
        not depend on it.

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


def _fold_envelope_pair_chunk(pending_pairs, finished_pairs, row_width, deduplicate):
    """Fold the per-bead pair arrays gathered so far into one array per pair class.

    Helper for ``build_all_envelope_pairs``. For each of the three pair classes
    (short range, long range, super long range) the per-bead arrays waiting in
    ``pending_pairs`` are concatenated into a single flat ``(n_pairs,
    row_width)`` array, optionally de-duplicated, appended to
    ``finished_pairs``, and the pending list is emptied. Concatenating copies
    the rows out of the per-bead arrays, so the padded buffers those arrays are
    views of can be freed as soon as the pending list is cleared - which is the
    point of folding in chunks rather than once at the end.

    Parameters
    ----------
    pending_pairs : tuple of list of numpy.ndarray
        Three lists (short range, long range, super long range), each holding
        the ``(n, 2, n_dim)`` arrays returned by the per-bead extractor since
        the last fold. Emptied in place.

    finished_pairs : tuple of list of numpy.ndarray
        Three lists that receive one flat ``(n_pairs, row_width)`` array per
        fold for each pair class with anything pending. Appended to in place.

    row_width : int
        Number of integers in one flattened pair: ``2 * n_dim`` (4 in 2D, 6 in
        3D).

    deduplicate : bool
        If True each folded array is passed through ``_unique_rows`` before
        being stored, so duplicates within the chunk are dropped now rather
        than carried to the end. If False the rows are stored in bead order.

    Returns
    -------
    None
        ``pending_pairs`` and ``finished_pairs`` are modified in place.

    """
    for pair_class in range(0, 3):
        if len(pending_pairs[pair_class]) == 0:
            continue

        rows = np.concatenate(pending_pairs[pair_class])
        rows = np.reshape(rows, (len(rows), row_width))
        pending_pairs[pair_class].clear()

        if deduplicate:
            rows = _unique_rows(rows)

        finished_pairs[pair_class].append(rows)


#-----------------------------------------------------------------
#
#@profile
def build_all_envelope_pairs(positions, LR_binary_array, type_lattice, dimensions, hardwall=False, deduplicate=True):
    """
    Expects a LIST of positions and a numpy array of positions which engage in
    long-range interactions (or not) - 0 if not and 1 if yes.

    Returns a list of tuples, where each tuple is a pair of positions.

    The complete set of these positions represents the non-redundant set of long
    -range and short-range pairwise interactions associated with the positions
    defined in position

    Parameters
    ----------
    positions : list
        A list of positions (each a 2D or 3D list) for which the enveloping
        interaction pairs are required.

    LR_binary_array : numpy.ndarray
        1D integer array aligned with ``positions``, holding 1 where the bead at
        that position engages in long-range interactions and 0 where it does
        not (as returned by ``Chain.get_LR_binary_array()``).

    type_lattice : numpy.ndarray
        The type grid (2D or 3D integer array holding the residue integer code
        at each site) used by the inner-loop routines to determine which
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


    num_dims = len(dimensions)

    # the per-bead extractor: the same four compiled routines as before, picked
    # once rather than inside the loop
    if num_dims == 2:
        if hardwall:
            extract_pairs = inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_2D_hardwall
        else:
            extract_pairs = inner_loops.extract_SR_and_LR_pairs_from_position_2D
        box = (dimensions[0], dimensions[1])
    else:
        if hardwall:
            extract_pairs = inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_3D_hardwall
        else:
            extract_pairs = inner_loops.extract_SR_and_LR_pairs_from_position_3D
        box = (dimensions[0], dimensions[1], dimensions[2])

    # The LR and SLR arrays an extractor returns are VIEWS of padded per-bead
    # buffers ((98, 2, 3) and (218, 2, 3) int32 in 3D), so a bead with a handful
    # of long-range neighbours still pins ~8 kB for as long as its views are
    # alive. Holding every bead's views until one final concatenate cost about
    # 14.5 kB per bead at every whole-system energy evaluation (start-up, every
    # ENERGY_CHECK, every restart write) - 15 GB at 10^6 beads. So the per-bead
    # arrays are folded into one array per chunk of beads as we go (which frees
    # the buffers), and each chunk is de-duplicated straight away so what is
    # kept between chunks is already close to its final size.
    #
    # index 0 / 1 / 2 = short range / long range / super long range
    pending_pairs  = ([], [], [])
    finished_pairs = ([], [], [])

    for i in range(0, len(positions)):

        # get enveloping pairs
        (SR_tmp, LR_tmp, SLR_tmp) = extract_pairs(np.array(positions[i], dtype=NP_INT_TYPE), LR_binary_array[i], type_lattice, *box)

        pending_pairs[0].append(SR_tmp)

        # note we have to check LR and SLR pairs seperately !
        if len(LR_tmp) > 0:
            pending_pairs[1].append(LR_tmp)

        if len(SLR_tmp) > 0:
            pending_pairs[2].append(SLR_tmp)

        if (i + 1) % _ENVELOPE_PAIR_CHUNK == 0:
            _fold_envelope_pair_chunk(pending_pairs, finished_pairs, 2 * num_dims, deduplicate)

    _fold_envelope_pair_chunk(pending_pairs, finished_pairs, 2 * num_dims, deduplicate)

    # Stitch the chunks together. The result is the same array, row for row, as
    # de-duplicating one concatenation of every bead's pairs: _unique_rows
    # returns the unique rows in sorted order, and neither the set of unique
    # rows nor the minimum and maximum coordinate (which decide how the rows
    # are sorted) is changed by having removed some duplicates early. Without
    # de-duplication the chunks are simply joined back in bead order.
    all_pairs = []
    for pair_class in range(0, 3):
        chunks = finished_pairs[pair_class]

        if len(chunks) == 0:
            rows = np.empty((0, 2 * num_dims), dtype=NP_INT_TYPE)

        elif len(chunks) == 1:
            # a single chunk (every single-chain and cluster move) was already
            # de-duplicated when it was folded, so there is nothing left to do
            rows = chunks[0]

        else:
            rows = np.concatenate(chunks)
            if deduplicate:
                rows = _unique_rows(rows)

        all_pairs.append(np.reshape(rows, (len(rows), 2, num_dims)))

    return (all_pairs[0], all_pairs[1], all_pairs[2])


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

    lattice_grid : numpy.ndarray
        Standard lattice occupancy grid: a 2D or 3D integer array holding the chainID
        at each site (0 = solvent).

    chainDict : dict
        Dictionary containing a mapping of each chainID to either a list of positions
        associated with that chain, or the Chain object associated with that chainID
        (see useChains).

    threshold : int or None, optional
        The max size of the connected component we are looking for. If None (the
        default), this is ignored, but if set and we generate a connected component
        larger than this, we will raise a ClusterSizeThresholdException. Enables us to
        avoid situations where we've moving massive giant clusters around which may not
        be efficient if 90% of the chains are in the cluster.

    useChains : bool, optional
        Boolean flag which defines if the chainDict is a true dictionary mapping chainID
        to a set of positions, or in fact a dictionary of Chain objects (which contain
        positions which must be accessed using the .get_ordered_positions()). This isn't
        so much a feature as the fact that we want this function to be able to accept
        two different types of chain information (dictionary of lists of positions or
        dictionary of chain objects). Default is True (Chain objects).

    hardwall : bool, optional
        Boolean flag which defines if we are using a hardwall potential or not. Both modes
        run the same breadth-first search over the chain-chain contact (envelope) pairs;
        the flag only selects whether contacts are read with periodic wrapping or with
        wall clipping (no contacts across a wall). Default is False.

    Returns
    ---------
    list
        A list of chainIDs associated with the chains in the connected component
        which contains the chain defined by $chainID

    Raises
    ---------
    ClusterSizeThresholdException
        If threshold is set and the connected component grows beyond it.

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
def get_all_chains_in_long_range_cluster(chainID, latticeObject, hardwall=False,
                                         LR_table=None, SLR_table=None):

    """
    Function which given a chainID, a dictionary of chain-to-position
    mappings, and a lattice grid will return the set of chains in the
    connected component where connectivity is defined in terms of
    long-range interactions.  Note a connected component is a
    *heterotypic* structure - i.e. we are looking for a connected
    component made up of *any* chains, not a single type of chain.

    Two chains are connected by any short-range contact (Chebyshev
    distance 1, regardless of the interaction energy) or by a Chebyshev-2
    (LR) / Chebyshev-3 (SLR) pair with **nonzero interaction energy**. When
    the interaction tables are supplied the energy is looked up directly, so a
    parameter file without an SLR column (all-zero SLR table) or an explicit
    zero entry for a residue pair does not connect chains through that shell.
    Without the tables the weaker structural rule is used: both beads of the
    pair must be LR-capable, which is exactly the set of pairs whose table
    entry *can* be nonzero.

    Parameters
    ----------

    chainID : int
        The chainID of the chain we initially are asking about

    latticeObject : Lattice
        Lattice object containing the lattice grid, the type grid,
        and the chain dictionary

    hardwall : bool, optional
        Boolean flag which defines if we are using a hardwall potential
        or not. Default is False (periodic boundaries).

    LR_table : numpy.ndarray or None, optional
        The long-range residue interaction table: an
        ``(n_residues, n_residues)`` integer array indexed by the integer
        residue codes stored in the lattice's ``type_grid``. When provided,
        Chebyshev-2 pairs connect chains only if their entry is nonzero.
        Default is None.

    SLR_table : numpy.ndarray or None, optional
        The super-long-range residue interaction table, same shape and indexed
        the same way. When provided, Chebyshev-3 pairs connect chains only if
        their entry is nonzero. Default is None.

    Returns
    ---------
    list
        A list of chainIDs associated with the chains in the connected
        component which contains the chain defined by $chainID

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

    # Boolean grid marking every occupied site whose bead is LR-flagged. The
    # LR/SLR envelope extractors emit pairs only FROM LR-flagged beads TO any
    # occupied site, which is a DIRECTED relation: treating it as undirected
    # reachability made the "cluster" depend on the BFS seed - in mixed-flag
    # systems the same configuration decomposed into different, overlapping,
    # chain-double-counting "clusters" (not a partition at all). We therefore
    # symmetrise with the ENERGY-BASED definition: an LR/SLR edge exists only
    # when BOTH endpoint beads are LR-flagged (exactly the pairs with nonzero
    # LR/SLR energy - the tables are zero unless both residues are LR), which
    # is manifestly symmetric. SR edges (any contact) are unchanged.
    #
    # The LR-flag grid is only read by the structural rule (a table of None).
    # Building it is a whole-box allocation plus a Python loop over every bead
    # of every chain, and this function runs once per cluster seed, so building
    # it unconditionally made the long-range cluster analysis quadratic in the
    # number of chains for a dilute system - for a grid nothing read, because
    # the simulation always passes both tables.
    lr_flag_grid = None
    if LR_table is None or SLR_table is None:
        lr_flag_grid = np.zeros(lattice_grid.shape, dtype=bool)
        for _cid, _chain in chainDict.items():
            _flags = _chain.get_LR_binary_array()
            for _p, _f in zip(_chain.get_ordered_positions(), _flags):
                if _f == 1:
                    lr_flag_grid[tuple(_p)] = True

    def _both_endpoints_LR(pairs, table=None):
        """Keep only pairs that carry a long-range interaction.

        With ``table`` given, that means a nonzero table entry for the two
        residue types (the documented "nonzero interaction energy" rule, which
        also excludes shells the parameter file never enabled). Without it,
        both sites must hold LR-flagged beads.

        Parameters
        ----------
        pairs : numpy.ndarray
            ``(n_pairs, 2, n_dim)`` array of LR or SLR candidate pairs, as
            returned by build_all_envelope_pairs(). An empty input is returned
            unchanged.

        table : numpy.ndarray or None, optional
            The LR or SLR residue interaction table, an
            ``(n_residues, n_residues)`` integer array indexed by the residue
            codes held in the type grid. If None (the default) the structural
            rule is used instead: both endpoints must be LR-flagged beads.

        Returns
        -------
        numpy.ndarray
            The subset of ``pairs`` whose two endpoints carry a long-range
            interaction under whichever of the two rules applies.

        """
        if len(pairs) == 0:
            return pairs
        arr = np.asarray(pairs)
        sites = arr.reshape(-1, lattice_grid.ndim)
        if lattice_grid.ndim == 2:
            index = (sites[:, 0], sites[:, 1])
        else:
            index = (sites[:, 0], sites[:, 1], sites[:, 2])
        if table is None:
            keep = lr_flag_grid[index].reshape(-1, 2).all(axis=1)
        else:
            types = type_grid[index].reshape(-1, 2)
            keep = np.asarray(table)[types[:, 0], types[:, 1]] != 0
        return arr[keep]

    # loop until we break with a return statement

    while True:

        # get all the envelope pairs assoiated with the list of positions
        # deduplicate=False: as above, these only feed a set of chainIDs
        (SR_pairs, LR_pairs, SLR_pairs) = build_all_envelope_pairs(positions, LR_binary_array, type_grid, dimensions, hardwall, deduplicate=False)

        # A long-range cluster is connected by any enabled interaction shell.
        # Omitting SLR pairs silently split components joined at Chebyshev
        # distance three, despite the public contract including SLR contacts.
        # LR/SLR pairs are filtered to those with nonzero interaction energy
        # (or, without the tables, to both-endpoints-LR - see above) so the
        # connectivity is symmetric and matches the Hamiltonian. The table
        # check matters: a parameter file with no SLR column leaves the SLR
        # table all-zero, so two LR chains at Chebyshev distance three have
        # zero interaction energy and must NOT be reported as one cluster.
        envelope_pairs = np.concatenate((SR_pairs,
                                         _both_endpoints_LR(LR_pairs, LR_table),
                                         _both_endpoints_LR(SLR_pairs, SLR_table)))
                
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

    COM is calculated by implementing the algorithm developed by Bai
    and Breen [1] extended to 3D, which means it determines the
    correct center of mass in a periodic box.

    [1] Bai, L., & Breen, D. (2008). Calculating Center of Mass in an 
    Unbounded 2D Environment. Journal of Graphics, GPU, and Game 
    Tools, 13(4), 53 - 60.

    Parameters
    ----------
    positions : list
        List of positions (each a 2D or 3D coordinate) to calculate the center
        of mass from. Must not be empty.

    dimensions : list
        List of the box dimensions (2 or 3 ints); the length sets how many axes
        are computed and each value is the period used for the circular mean.

    on_lattice : bool, optional
        If True (the default), the center of mass is rounded to the nearest
        lattice site and returned as ints. If False, the center of mass is
        returned as floats in Euclidean space.

    Returns
    -------
    list
        The center of mass of the positions, wrapped back into the box (2D or
        3D list depending on if the input positions are 2D or 3D; ints if
        on_lattice is True, otherwise floats).

    Raises
    -------
    LatticeUtilsException
        If positions is empty.

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
        Position (2D or 3D coordinate) to delete the residue from.

    lattice : numpy.ndarray
        The 2D or 3D lattice occupancy grid to delete the residue from,
        modified in place.

    chainID : int or None, optional
        Chain ID of the residue to delete. If None (the default), the residue
        will be deleted regardless of the chain ID.

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
        Position (2D or 3D coordinate) to insert the residue

    lattice : numpy.ndarray
        The 2D or 3D lattice occupancy grid to insert the residue into,
        modified in place

    chainID : int
        Chain ID to insert the residue for (this is the value written into the
        grid)

    safe : bool, optional
        If True (the default), will raise an exception if the position is
        already occupied. If False, will overwrite the residue.

    Returns
    -------
    None
        No return value, but the lattice grid is updated in place.

    Raises
    ------
    ResidueAugmentException
        Raised (only when `safe` is True) if the target site is already
        occupied.

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
        List of positions to rotate (each a 2- or 3-element coordinate,
        expressed relative to the rotation origin)

    rotation_matrix : numpy.ndarray
        The rotation matrix to apply: a ``(2, 2)`` or ``(3, 3)`` integer array
        from the CONFIG cardinal rotation tables

    Returns
    -------
    list
        List of rotated positions, each a plain list of the same
        dimensionality as the input position. Integer input gives Python
        ints: the rotated coordinates end up in ``Chain.positions`` (chain
        rotate, chain pivot and cluster rotate all build their new positions
        from them), and a numpy scalar there wraps silently in later scalar
        arithmetic and is pickled into the restart file as a numpy object.

    """
    rotated_positions = []
    for position in positions:
        rotated_positions.append(np.dot(rotation_matrix, position).tolist())

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
        List of 3D positions to rotate (expressed relative to the origin the
        rotation is about)

    dimension : str
        Dimension to rotate around. Must be one of 'x', 'y', or 'z'.

    degrees : int
        Degrees to rotate by. Must be one of 90, 180, or 270.

    Returns
    -------
    list
        List of rotated positions, each a list of length 3 (Python ints for
        integer input; see ``run_rotation``)

    Raises
    -------
    RotationException
        If dimension is not 'x'/'y'/'z' or degrees is not 90/180/270.

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
        A list of 2D positions to be rotated (expressed relative to the origin
        the rotation is about).

    degrees : int
        The number of degrees to rotate the positions by. Must be one
        of 90, 180, or 270.

    Returns
    -------------
    list
        A list of 2D positions (each a list of length 2; Python ints for
        integer input, see ``run_rotation``) that have been rotated by the
        specified number of degrees.

    Raises
    -------------
    RotationException
        If degrees is not 90, 180 or 270.

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

    filename : str, optional
        Filename to write to. Default is lattice.pdb. Any existing file of this
        name is replaced by the new (empty) PDB with its CRYST1 header.

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
    latticeObject : Lattice
        Current Lattice object, whose chains are written out as one frame

    spacing : float
        Lattice-to-realspace spacing in angstroms.

    filename : str, optional
        Filename to write to. This file must already have been initialized by
        open_pdb_file(). Default is lattice.pdb.

    write_connect : bool, optional
        Flag to write CONECT records into the PDB file. Default is False.

    autocenter : bool, optional
        Flag to center the chain in the box in the PDB file. Default is False.
        Autocentring only applies to a single-chain system: with more than one
        chain it is ignored, and we say so (once) with a warning that names the
        AUTOCENTER keyword.

    unwrap : bool, optional
        Flag which, if True, writes each chain as a single whole periodic image
        (bond-walked so it is not torn across a box face), so coordinates may
        fall outside the box. Ignored where autocenter applies, since that
        already unwraps. Default is False.

    Returns
    ------------
    None
        No return value, but a MODEL/ATOM/TER/ENDMDL block is appended to the
        PDB file on disk.

    """
    autocenter = _resolve_autocenter(latticeObject, autocenter)
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
def start_xtc_file(lattice, spacing, pdb_filename='START.pdb', xtc_filename='traj.xtc', autocenter=False, unwrap=False):
    """
    Function that initializes a new .xtc file. This deletes an existing XTC file 
    of the same name to avoid any issues.

    Parameters
    ------------
    lattice : Lattice
        Current Lattice object, written out as the first frame

    spacing : float
        Lattice-to-realspace spacing in angstroms, used when writing the
        corresponding PDB file.

    pdb_filename : str, optional
        New XTC files need a corresponding PDB file (the topology). This defines
        the name of that PDB file. Default is START.pdb.

    xtc_filename : str, optional
        Name of the XTC trajectory file to create. Any existing file of this
        name is deleted first. Default is traj.xtc.

    autocenter : bool, optional
        Flag which, if True and the system holds a single chain, centres that
        chain in the box in the written frame. Default is False.

    unwrap : bool, optional
        Flag which, if True, writes each chain as a single whole periodic image
        (bond-walked, so coordinates may fall outside the box). Ignored where
        autocenter applies. Default is False.

    Returns
    ------------
    None
        No return value, but a newly initialized XTC file (and its topology PDB)
        is generated. The XTC holds one frame - the current lattice - written
        through the same code as every frame of the streamed trajectory, so
        its coordinates and its unit cell (``DIMENSIONS x LATTICE_TO_ANGSTROMS``)
        are exactly what ``SAVE_AT_END : False`` would have written.

    Raises
    ------
    PDBException
        If a coordinate does not fit the PDB columns (see
        :func:`write_topology_pdb`). Neither file is left half-written.

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
    write_topology_pdb(lattice, spacing, pdb_filename, autocenter=autocenter, unwrap=unwrap)

    # ...then frame 0 of the XTC, straight from the lattice. This used to be
    # md.load(pdb).save_xtc(), which took the unit cell from the CRYST1 record:
    # rounded to 0.001 A, absent altogether when mdtraj judges the box too dense
    # to be real (below about 1 A per site), and given spurious ~1e-6 nm
    # off-diagonal terms by mdtraj's lengths/angles round trip in boxes wider
    # than about 23 nm.
    writer = _XTCStreamWriter(md.formats.XTCTrajectoryFile(xtc_filename, 'w'))
    try:
        xyz, box = _lattice_frame_xyz_and_box(lattice, spacing, autocenter=autocenter, unwrap=unwrap)
        writer.write(xyz, box=box)
    finally:
        writer.close()

    # remember the unit cell this trajectory was started with, for the
    # SAVE_AT_END buffer that will be created for it later (see
    # TrajectoryAccumulator: it needs the cell if it ever has to rebuild
    # frame 0, and may be created without a lattice to take it from)
    _STARTED_TRAJECTORY_BOX[os.path.abspath(pdb_filename)] = box


#-----------------------------------------------------------------
#
def write_topology_pdb(lattice, spacing, pdb_filename, autocenter=False, unwrap=False):
    """
    Write the topology PDB (``START.pdb`` and friends) for the current lattice,
    all or nothing.

    The file is built under a temporary name in the same directory and moved
    into place with ``os.replace`` only once it is complete. A write that fails
    part-way - the usual cause being a coordinate too wide for the PDB columns -
    therefore never leaves a truncated ``START.pdb`` behind, and never destroys
    the ``START.pdb`` of an earlier run that the rest of that run's files still
    belong to. Because the finished file is moved into place rather than
    written through, a symbolic link at ``pdb_filename`` is replaced by a
    regular file and the link's target is left alone.

    Parameters
    ----------
    lattice : Lattice
        Current Lattice object, written out as the single model of the PDB.

    spacing : float
        Lattice-to-realspace spacing in angstroms.

    pdb_filename : str
        Name of the PDB file to create or replace.

    autocenter : bool, optional
        Single-chain autocentring (see build_pdb_file). Default is False.

    unwrap : bool, optional
        Make chains whole across PBC before writing. Default is False.

    Returns
    -------
    None
        No return value, but the complete PDB file is on disk.

    Raises
    ------
    PDBException
        If the file cannot be built, e.g. a coordinate at or beyond 10000 A
        (or at or below -1000 A) does not fit the eight PDB coordinate
        columns. The message names the keywords responsible (``DIMENSIONS`` and
        ``LATTICE_TO_ANGSTROMS``, or ``TRAJECTORY_PBC_UNWRAP`` / ``AUTOCENTER``
        when the coordinate lies outside the box).

        The temporary file is removed before the exception propagates.

    OSError
        If the finished file cannot be moved into place (a directory sits
        under ``pdb_filename``, say). The temporary file is removed here too.

    """
    tmp_filename = '%s.tmp.%i' % (pdb_filename, os.getpid())

    try:
        open_pdb_file(lattice.dimensions, spacing, filename=tmp_filename)
        write_lattice_to_pdb(lattice, spacing, filename=tmp_filename, write_connect=True, autocenter=autocenter, unwrap=unwrap)
        finish_pdb_file(tmp_filename)

        # (the three writers above can be swapped out for stubs that write
        # nothing, in which case there is nothing to move into place)
        if os.path.exists(tmp_filename):
            os.replace(tmp_filename, pdb_filename)
    except BaseException:
        # do not leave the partial (or unmovable) file behind, and do not let a
        # failure to tidy up hide the error that matters
        try:
            os.remove(tmp_filename)
        except OSError:
            pass
        raise


def _resolve_autocenter(lattice, autocenter):
    """
    Decide whether autocentring applies to this lattice, and say so (once) if
    it was asked for but cannot be honoured.

    Parameters
    ----------
    lattice : Lattice
        The Lattice object being written out.

    autocenter : bool
        Whether ``AUTOCENTER`` was requested.

    Returns
    -------
    bool
        True if autocentring was requested and the lattice holds exactly one
        chain. False otherwise; if it was requested for a multi-chain system a
        warning naming the keyword is printed (and logged to ``log.txt`` if the
        log exists) the first time this happens.

    """
    global _AUTOCENTER_IGNORED_WARNED

    if not autocenter:
        return False

    if len(lattice.chains) == 1:
        return True

    if not _AUTOCENTER_IGNORED_WARNED:
        _AUTOCENTER_IGNORED_WARNED = True
        msg = ("AUTOCENTER : True is ignored because the system holds %i chains. Autocentring is "
               "only defined for a single chain, so START.pdb and the trajectory are written "
               "uncentred. Remove AUTOCENTER, or use TRAJECTORY_PBC_UNWRAP : True if the aim is to "
               "keep chains whole across the periodic boundary." % len(lattice.chains))
        IO_utils.status_message(msg, 'warning')
        if os.path.exists(CONFIG.OUTNAME_LOGFILE):
            # imported here because pimmslogger is not otherwise needed by this module
            from . import pimmslogger
            pimmslogger.log_warning(msg)

    return False


#-----------------------------------------------------------------
#
def get_chain_output_positions(lattice, chainID, autocenter=False, unwrap=False):
    """
    Return one chain's positions under the output convention in force, as they
    should be written to the topology PDB and to every trajectory frame.

    This is ``Chain.get_output_positions`` plus the one thing the chain cannot
    decide alone: in a hardwall box an autocentred chain must not be carried
    through a wall (see :func:`clamp_positions_to_box`). Every writer (PDB,
    streamed XTC, ``SAVE_AT_END`` buffer) goes through here so that the three
    can never disagree.

    Parameters
    ----------
    lattice : Lattice
        The Lattice object that owns the chain.

    chainID : int
        The ID of the chain whose positions are wanted.

    autocenter : bool, optional
        Whether autocentring applies. The caller is expected to have resolved
        this with :func:`_resolve_autocenter` (single chain only). Default is
        False.

    unwrap : bool, optional
        Write the chain as a single whole periodic image. Ignored where
        autocenter applies. Default is False.

    Returns
    -------
    list
        The chain's bead positions, in N->C order.

    """
    positions = lattice.chains[chainID].get_output_positions(autocenter=autocenter, unwrap=unwrap)
    if autocenter and getattr(lattice, 'hardwall', False):
        positions = clamp_positions_to_box(positions, lattice.dimensions)
    return positions


#-----------------------------------------------------------------
#
def _lattice_frame_xyz_and_box(lattice, spacing, autocenter=False, unwrap=False):
    """
    Build the ``(1, n_beads, 3)`` coordinate array (in nm) and the orthorhombic box
    (in nm) for one trajectory frame from the current lattice.

    Positions are gathered per chain via ``get_output_positions`` (honouring the
    ``autocenter`` / ``unwrap`` conventions); 2D systems are padded with a zero z.

    Parameters
    ----------
    lattice : Lattice
        The Lattice object to snapshot. Its chains are visited in dictionary
        order, which is the same order the topology PDB was written in.

    spacing : float
        Lattice-to-realspace spacing in angstroms. Coordinates are multiplied by
        ``spacing * 0.1`` to convert lattice units to nm.

    autocenter : bool, optional
        If True, centre the chain in the box. Only meaningful for a single-chain
        system; switched off (with a one-off warning) when the lattice holds
        more than one chain. In a hardwall box the centring shift is clamped so
        that no bead is written outside the box. Default is False.

    unwrap : bool, optional
        If True, write each chain as a single whole periodic image (bond-walked,
        so coordinates may fall outside the box). Ignored where autocenter
        applies. Default is False.

    Returns
    -------
    tuple
        ``(xyz, box)`` where ``xyz`` is float32 shape ``(1, n_beads, 3)`` and
        ``box`` is float32 shape ``(1, 3, 3)`` (diagonal box vectors, nm).

    """
    # autocenter is only meaningful for a single chain
    autocenter = _resolve_autocenter(lattice, autocenter)

    is_3d = len(lattice.dimensions) == 3
    cvals = []
    for chainID in lattice.chains:
        positions = get_chain_output_positions(lattice, chainID, autocenter=autocenter, unwrap=unwrap)
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
class _XTCStreamWriter:
    """
    Thin wrapper around an open ``mdtraj.formats.XTCTrajectoryFile`` that stamps
    each frame with a monotonically increasing time/step.

    ``XTCTrajectoryFile.write`` defaults ``time`` and ``step`` to zero for every
    frame, which made all streamed frames indistinguishable in tools that read
    the XTC time/step metadata (the buffered ``SAVE_AT_END`` path, in contrast,
    wrote 0, 1, 2, ...). The wrapper keeps the two output paths consistent.
    Frames are numbered 0, 1, 2, ... in the order they are written (the saved
    frames, not the Monte Carlo step numbers).

    Every frame is pushed to the operating system as soon as it is written
    (see :meth:`write` and :meth:`flush`), so the file on disk always holds
    every completed frame. Without that, frames sat in a 4 KiB stdio buffer
    inside the process: a run that was killed lost up to 32 frames of a small
    system and left a file that ended mid-frame.

    What is NOT guaranteed is that the file ends on a frame boundary at every
    instant. A frame larger than that 4 KiB buffer reaches the file in pieces
    while it is being written, so a process reading the trajectory of a run
    that is still going (or the file left by a kill that landed inside a
    write) can find part of a frame after the last complete one. Readers of a
    live or killed trajectory should read frame by frame and stop at the first
    frame that does not decode.
    """

    def __init__(self, fh, first_frame_index=0):
        """
        Wrap an already-open XTC file handle and start the frame counter.

        Parameters
        ----------
        fh : mdtraj.formats.XTCTrajectoryFile
            An XTC file handle opened for writing. The wrapper takes ownership
            only in the sense that close() closes it; the caller is responsible
            for having opened it in the right mode.

        first_frame_index : int, optional
            The time/step stamp given to the first frame written. Default is 0;
            the ``SAVE_AT_END`` writer passes the number of frames already on
            disk so the frames it adds carry on the count.

        Returns
        -------
        None
            No return value; the handle and frame counter are stored on the
            new object.

        """
        self._fh = fh
        self.frame_index = int(first_frame_index)

    def write(self, xyz, box=None):
        """
        Write one frame, stamping it with the next frame index as time and step,
        and flush it so it is on disk before we return.

        Parameters
        ----------
        xyz : numpy.ndarray
            float32 array of shape ``(1, n_beads, 3)`` holding the frame's
            coordinates in nm.

        box : numpy.ndarray or None, optional
            float32 array of shape ``(1, 3, 3)`` giving the box vectors in nm.
            If None (the default) no box is recorded for the frame.

        Returns
        -------
        None
            No return value, but the frame is written to the underlying file,
            flushed, and the frame counter is advanced.

        """
        self._fh.write(xyz,
                       time=np.array([float(self.frame_index)], dtype=np.float32),
                       step=np.array([self.frame_index], dtype=np.int32),
                       box=box)
        self.frame_index += 1
        self.flush()

    def flush(self):
        """
        Hand every frame written so far to the operating system.

        After this returns, the file on disk holds every frame passed to
        :meth:`write`, whole, and nothing else: at that moment another process
        can read all of it. (Between flushes, while the next frame is being
        written, such a reader may find the start of that frame at the end of
        the file - see the class docstring.) The flushed frames survive the
        PIMMS process being killed outright (``SIGKILL``, an out-of-memory
        kill, a scheduler's wall-time limit).
        This is a flush, not an ``fsync``: the bytes are in the operating
        system's hands, which protects against the death of the process but not
        against the machine itself losing power.

        Safe to call at any time and any number of times: it does nothing if
        the writer has been closed, or if nothing has been written since the
        last call. :meth:`write` already calls it after every frame, so an
        explicit call (e.g. before a restart checkpoint is written) is belt and
        braces.

        Returns
        -------
        None
            No return value.

        """
        fh = self._fh
        if fh is None or not getattr(fh, 'is_open', True):
            return
        flush = getattr(fh, 'flush', None)
        if flush is not None:
            flush()

    def close(self):
        """
        Close the underlying XTC file handle.

        Returns
        -------
        None
            No return value; the file is closed and no further frames can be
            written. Calling close() again is harmless.

        """
        if self._fh is not None:
            self._fh.close()
            self._fh = None


#-----------------------------------------------------------------
#
def open_xtc_writer(lattice, spacing, pdb_filename='START.pdb', xtc_filename='traj.xtc', autocenter=False, unwrap=False):
    """
    Write the topology PDB, open a persistent XTC writer, write the first frame, and
    return the open writer handle.

    This is the efficient replacement for the previous per-frame approach (the
    since-removed ``append_to_xtc_file_non_redundant``), which re-loaded the entire growing
    trajectory from disk and re-saved it on EVERY frame - O(frames^2) in both wall
    time and disk I/O. Here a single ``mdtraj`` XTC file handle is kept open for the
    whole run and each frame is appended with :func:`write_xtc_frame` in O(1); the
    handle is closed with :func:`close_xtc_writer`.

    Parameters
    ----------
    lattice : Lattice
        Current Lattice object, written out as frame 0.
    spacing : float
        Lattice-to-realspace spacing (angstroms).
    pdb_filename : str, optional
        Topology PDB filename to (re)write. Default is START.pdb.
    xtc_filename : str, optional
        Trajectory filename to create; an existing file of that name is deleted
        first. Default is traj.xtc.
    autocenter : bool, optional
        Single-chain autocentring (see build_pdb_file). Default False.
    unwrap : bool, optional
        Make chains whole across PBC before writing (TRAJECTORY_PBC_UNWRAP). Default False.

    Returns
    -------
    _XTCStreamWriter
        The open writer handle (write more frames with write_xtc_frame, then close
        with close_xtc_writer). Frames are stamped with sequential time/step
        metadata (0, 1, 2, ...), and each one is flushed to disk as it is
        written.

    Raises
    ------
    PDBException
        If a coordinate does not fit the PDB columns (see
        :func:`write_topology_pdb`); no truncated PDB is left behind.
    """
    # (re)write the topology PDB (same autocenter/unwrap conventions as the frames)
    write_topology_pdb(lattice, spacing, pdb_filename, autocenter=autocenter, unwrap=unwrap)

    # start a fresh XTC file and write the first frame. lexists, not exists: a
    # symbolic link whose target is missing does not "exist", and opening it
    # for writing would create the trajectory wherever the link points
    if os.path.lexists(xtc_filename):
        os.remove(xtc_filename)
    writer = _XTCStreamWriter(md.formats.XTCTrajectoryFile(xtc_filename, 'w'))
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
    writer : _XTCStreamWriter
        Open writer handle from :func:`open_xtc_writer`.
    lattice : Lattice
        Current Lattice object, snapshotted into the new frame.
    spacing : float
        Lattice-to-realspace spacing (angstroms).
    autocenter : bool, optional
        Single-chain autocentring. Default False.
    unwrap : bool, optional
        Make chains whole across PBC before writing. Default False.

    Returns
    -------
    None
        No return value, but one frame is appended to the open XTC file.

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
    writer : _XTCStreamWriter or None
        The writer handle to close. None is accepted and does nothing, so the
        caller does not have to know whether a writer was ever opened.

    Returns
    -------
    None
        No return value; the underlying XTC file is closed and flushed.

    """
    if writer is not None:
        writer.close()


#-----------------------------------------------------------------
#
def flush_xtc_writer(writer):
    """
    Make sure every trajectory frame written so far is on disk. Safe to call
    with anything the run may be holding as its trajectory.

    This is the call to make before a restart checkpoint is written, so that a
    checkpoint never describes a step whose frame is not yet in ``traj.xtc``.
    The streamed writer already flushes after every frame, so for it this is a
    cheap second line of defence; for the ``SAVE_AT_END`` buffer it does
    nothing, by design (that mode exists to keep the trajectory off the disk
    until the end - use :func:`save_out_sim` to write the buffer out).

    Parameters
    ----------
    writer : _XTCStreamWriter or TrajectoryAccumulator or None
        The streamed-trajectory writer from :func:`open_xtc_writer`, the
        ``SAVE_AT_END`` buffer, or None if no trajectory is open. A writer that
        has already been closed is fine too.

    Returns
    -------
    None
        No return value.

    """
    if writer is not None:
        writer.flush()


#-----------------------------------------------------------------
#
def stale_writer_temporaries(filenames=('START.pdb', 'traj.xtc', 'eq_START.pdb', 'eq_traj.xtc',
                                        'CONFIG_AT_ENERGY_FAIL.pdb', 'CONFIG_AT_ENERGY_FAIL.xtc')):
    """
    List the scratch files that killed runs left beside the trajectory and
    topology files.

    The topology PDB is built as ``<name>.tmp.<pid>`` and moved into place
    (:func:`write_topology_pdb`), and a ``SAVE_AT_END`` trajectory is extended
    through a scratch ``<name>.tmp.<pid>`` (:meth:`TrajectoryAccumulator.write_out`).
    Both remove their scratch file on any error, but a run killed outright in
    the middle of one of those writes cannot, and with one name per process no
    later run would ever overwrite it. Start-up therefore removes them, under
    the rule the checkpoint temporaries follow
    (``restart.stale_checkpoint_temporaries``): a scratch file whose process
    still exists on this machine is left alone, as it may be another run's
    write in progress.

    Parameters
    ----------
    filenames : iterable of str, optional
        The output files whose scratch files are looked for. Default is every
        topology and trajectory file a run can write.

    Returns
    -------
    list of str
        The scratch files that no running process owns.

    """
    # imported here: restart is not otherwise needed by this module
    from . import restart

    stale = []
    for filename in filenames:
        for name in glob.glob(glob.escape(filename) + '.tmp.*'):
            suffix = name.rsplit('.', 1)[-1]
            if not suffix.isdigit():
                continue
            if int(suffix) != os.getpid() and restart._process_is_running(int(suffix)):
                continue
            stale.append(name)
    return stale


#-----------------------------------------------------------------
#
class TrajectoryAccumulator:
    """In-memory buffer of the trajectory frames a ``SAVE_AT_END : True`` run has
    saved but not yet written.

    ``SAVE_AT_END`` keeps the trajectory off the disk while the run is going
    (for filesystems where frequent small writes are slow). The first frame is
    written with the topology PDB when the trajectory is started
    (:func:`start_xtc_file`); every later frame is buffered here as a bare
    float32 coordinate array, and :meth:`write_out` adds the buffered frames to
    the end of that file.

    The frames go through the same writer as the streamed trajectory
    (:class:`_XTCStreamWriter`), carrying the same coordinates, the same
    time/step stamps and the same unit cell built from ``DIMENSIONS x
    LATTICE_TO_ANGSTROMS``, so the finished file is byte-for-byte what
    ``SAVE_AT_END : False`` writes. Earlier versions buffered one
    ``mdtraj.Trajectory`` per frame, joined them and saved the result. That took
    the unit cell from the CRYST1 record of the PDB as mdtraj read it back
    (rounded to 0.001 A; missing, and a ``TypeError`` at the first frame, when
    mdtraj judged the box too dense to be real; with ~1e-6 nm off-diagonal terms
    in boxes wider than about 23 nm), and the final join held two to three
    copies of the whole trajectory in memory at once.

    Writing is incremental: :meth:`write_out` adds the buffered frames to the
    file a batch at a time and drops each batch from the buffer once it is
    wholly in the file, so calling it twice writes nothing twice, and it can
    be called from an error path to rescue what has been buffered. A batch
    that cannot be written in full is taken back out of the file, so a failed
    write leaves a valid trajectory (every frame written before the failure)
    and a buffer that still holds exactly the frames that are not on disk.
    """

    __slots__ = ("_pdb_filename", "_frames", "_box", "_n_appended")

    # Upper bound on the coordinates (in bytes of float32) written to the
    # scratch file in one batch. It bounds the disk space write_out needs on
    # top of the finished trajectory: the scratch file holds the compressed
    # form of one batch, never the whole buffer. A single frame larger than
    # this is still written, as a batch of one.
    BATCH_BYTES = 64 * 1024 * 1024

    def __init__(self, pdb_filename=None):
        """
        Start an empty buffer.

        Parameters
        ----------
        pdb_filename : str or None, optional
            The topology PDB that was written when the trajectory was started.
            It identifies the trajectory: the unit cell :func:`start_xtc_file`
            recorded for it is picked up here, so the buffer knows the cell
            before any frame has been appended. The file itself is only read
            if, at :meth:`write_out` time, the XTC file the frames are to be
            added to has gone missing and frame 0 has to be rebuilt. Default is
            None.

        Returns
        -------
        None
            No return value; the empty frame buffer is stored on the new object.

        """
        self._pdb_filename = pdb_filename
        self._frames = collections.deque()
        self._box = None
        if pdb_filename is not None:
            self._box = _STARTED_TRAJECTORY_BOX.get(os.path.abspath(pdb_filename))
        self._n_appended = 0

    def append_frame(self, xyz, box=None):
        """Buffer one frame's ``(1, n_beads, 3)`` coordinates (nm).

        Parameters
        ----------
        xyz : numpy.ndarray
            The frame's coordinates in nm, shape ``(1, n_beads, 3)``, with beads
            in the same order as the topology. Stored as float32.

        box : numpy.ndarray or None, optional
            The box vectors in nm, float32 shape ``(1, 3, 3)``. The box is the
            same for every frame of one trajectory file, so only the most
            recent one is kept. Default is None (keep whatever was passed
            before).

        Returns
        -------
        TrajectoryAccumulator
            This accumulator, so the caller can rebind the same name on every
            frame (which is how update_master_traj is used).

        """
        self._frames.append(np.asarray(xyz, dtype=np.float32))
        if box is not None:
            self._box = box
        self._n_appended += 1
        return self

    def flush(self):
        """
        Do nothing: ``SAVE_AT_END`` frames stay in memory until the end.

        This exists so that the run can call ``flush()`` on whatever it holds
        as its trajectory (see :func:`flush_xtc_writer`) without having to
        know which saving mode is in use.

        Returns
        -------
        None
            No return value.

        """
        return None

    def write_out(self, xtc_filename):
        """
        Add every buffered frame to the end of ``xtc_filename`` and empty the
        buffer.

        The file is expected to exist already and to hold the frames written so
        far (frame 0 from :func:`start_xtc_file`, plus anything an earlier call
        to this method added); the new frames are stamped to carry on from the
        number of frames found there. If the file is missing, frame 0 is first
        rebuilt from the topology PDB.

        An XTC file is a plain sequence of self-contained frames, so adding
        frames is a matter of writing them to a scratch XTC beside the target
        and appending its bytes (mdtraj cannot open an XTC for appending). Peak
        memory is therefore the buffer itself: nothing is copied or joined.

        The frames go in batches of at most ``BATCH_BYTES`` of coordinates, so
        the scratch file never holds more than one batch. A batch is dropped
        from the buffer only once all of its bytes are in the trajectory; if a
        batch fails part-way (the disk fills up, the run is interrupted) the
        trajectory is cut back to where it was before that batch. Whatever
        happens, then, the file holds whole frames only, the buffer holds
        exactly the frames that are not in the file, and calling this again
        carries on from there.

        Parameters
        ----------
        xtc_filename : str
            The trajectory file to add the buffered frames to.

        Returns
        -------
        int
            The number of frames written by this call (0 if the buffer was
            empty).

        Raises
        ------
        LatticeUtilsException
            If the trajectory file is missing and frame 0 cannot be rebuilt
            because no topology PDB is known or it cannot be read.

        OSError
            If the frames cannot be written (no space left, say). The
            trajectory is left as it was after the last complete batch.

        """
        if not os.path.exists(xtc_filename):
            # (a symbolic link whose target is missing must not be written through)
            if os.path.islink(xtc_filename):
                os.remove(xtc_filename)
            self._rebuild_first_frame(xtc_filename)

        if len(self._frames) == 0:
            return 0

        with md.formats.XTCTrajectoryFile(xtc_filename, 'r') as fh:
            n_on_disk = len(fh)

        tmp_filename = '%s.tmp.%i' % (xtc_filename, os.getpid())
        n_written = 0
        try:
            while len(self._frames) > 0:

                # how many frames fit in one batch (always at least one)
                n_batch, batch_bytes = 0, 0
                for xyz in self._frames:
                    if n_batch > 0 and batch_bytes + xyz.nbytes > self.BATCH_BYTES:
                        break
                    n_batch += 1
                    batch_bytes += xyz.nbytes

                writer = _XTCStreamWriter(md.formats.XTCTrajectoryFile(tmp_filename, 'w'),
                                          first_frame_index=n_on_disk + n_written)
                try:
                    for xyz in itertools.islice(self._frames, n_batch):
                        writer.write(xyz, box=self._box)
                finally:
                    writer.close()

                # all or nothing: _append_file takes its bytes back out of the
                # trajectory if it cannot add every one of them
                self._append_file(tmp_filename, xtc_filename)

                # only now are these frames on disk, so only now do they leave
                # the buffer
                for _ in range(n_batch):
                    self._frames.popleft()
                n_written += n_batch
        finally:
            try:
                os.remove(tmp_filename)
            except OSError:
                pass

        return n_written

    @staticmethod
    def _append_file(source, target):
        """
        Append the bytes of one file to the end of another and flush them, all
        or nothing.

        If the copy does not complete - the disk fills up half-way, or the
        process is interrupted - the target is cut back to the length it had
        before, so it never ends in part of a frame.

        Parameters
        ----------
        source : str
            File whose bytes are read.

        target : str
            File the bytes are appended to. It must already exist.

        Returns
        -------
        None
            No return value.

        Raises
        ------
        OSError
            If the bytes cannot be read or written. The target has been
            restored to its previous length (if that too fails, the original
            error is still the one raised).

        """
        size_before = os.path.getsize(target)
        try:
            with open(source, 'rb') as fin, open(target, 'ab') as fout:
                shutil.copyfileobj(fin, fout)
                fout.flush()
        except BaseException:
            try:
                with open(target, 'r+b') as fh:
                    fh.truncate(size_before)
            except OSError:
                pass
            raise

    def _rebuild_first_frame(self, xtc_filename):
        """
        Recreate a missing trajectory file with frame 0 taken from the topology
        PDB.

        This is a fallback: the file is written by :func:`start_xtc_file` when
        the trajectory is started and normally is still there. If something
        removed it during the run we would rather write the trajectory with
        frame 0 read back from ``START.pdb`` (coordinates to 0.001 A) than
        lose the frames in the buffer.

        Parameters
        ----------
        xtc_filename : str
            The trajectory file to create.

        Returns
        -------
        None
            No return value, but ``xtc_filename`` exists and holds one frame.

        Raises
        ------
        LatticeUtilsException
            If no topology PDB is known, or it cannot be loaded.

        """
        if self._pdb_filename is None:
            raise LatticeUtilsException(
                'Cannot write the SAVE_AT_END trajectory: %s does not exist and no topology PDB '
                'is known to rebuild its first frame from' % xtc_filename)
        try:
            base = md.load(self._pdb_filename, top=self._pdb_filename)
        except Exception as e:
            raise LatticeUtilsException(
                f'Could not load pdb file: {self._pdb_filename}: {e}') from e

        # the cell this trajectory was started with (or that of the frames
        # appended since). Only if neither is known do we fall back on what
        # mdtraj made of the PDB, and then on the CRYST1 record itself: mdtraj
        # discards the cell of a box it judges too dense to be real, and a
        # frame without a unit cell is the one thing we must not write.
        box = self._box
        if box is None:
            box = base.unitcell_vectors
        if box is None:
            box = _cryst1_box(self._pdb_filename)
        writer = _XTCStreamWriter(md.formats.XTCTrajectoryFile(xtc_filename, 'w'))
        try:
            writer.write(np.asarray(base.xyz[:1], dtype=np.float32), box=box)
        finally:
            writer.close()

    def __len__(self):
        """
        Number of frames the finished trajectory will hold.

        Returns
        -------
        int
            The number of frames ever appended (whether still buffered or
            already written out) plus one for the frame written when the
            trajectory was started.

        """
        return 1 + self._n_appended


#-----------------------------------------------------------------
#
def _cryst1_box(pdb_filename):
    """
    Read the orthorhombic unit cell straight from the CRYST1 record of a PDB.

    A last resort for rebuilding frame 0 of a ``SAVE_AT_END`` trajectory whose
    file has gone missing, used only when the cell is not known any other way.
    The record holds three decimals in angstroms, so this is the cell to
    0.0001 nm rather than exactly.

    Parameters
    ----------
    pdb_filename : str
        The PDB file to read.

    Returns
    -------
    numpy.ndarray or None
        float32 array of shape ``(1, 3, 3)`` with the box vectors in nm, or
        None if the file has no readable CRYST1 record.

    """
    try:
        with open(pdb_filename, 'r') as fh:
            for line in fh:
                if line.startswith('CRYST1'):
                    edges = [float(line[6:15]) * 0.1, float(line[15:24]) * 0.1, float(line[24:33]) * 0.1]
                    return np.diag(np.array(edges, dtype=np.float32))[np.newaxis]
    except (OSError, ValueError):
        pass
    return None


#-----------------------------------------------------------------
#
def update_master_traj(lattice, spacing, master_traj, pdb_filename, autocenter=False, unwrap=False):

    """
    Low level function that adds the current lattice as one more frame of the
    trajectory. Rather than reading in and rewriting an XTC file on every call,
    it appends the frame to the passed master trajectory object.

    The frame is built by the same routine the streamed writer uses
    (``_lattice_frame_xyz_and_box``), so its coordinates and its unit cell
    (``DIMENSIONS x LATTICE_TO_ANGSTROMS``, not the CRYST1 record read back
    from the PDB) are exactly what ``SAVE_AT_END : False`` would have written.
    The trajectory file itself must already have been started with
    :func:`start_xtc_file`, which writes frame 0; the frames buffered here are
    added after it by :func:`save_out_sim`.

    Parameters
    -----------
    lattice : Lattice
        A Lattice object, whose current state becomes the new frame

    spacing : float
        Lattice-to-realspace spacing in angstroms.

    master_traj : TrajectoryAccumulator or None
        The master trajectory we build through the sim. Pass None on the first
        call to have a new, empty accumulator created.

    pdb_filename : str
        The current PDB filename (the START.pdb written by simulation.py).
        Nothing is read from it here; the accumulator only remembers the name
        so it can rebuild frame 0 if the trajectory file has gone missing by
        the time the buffer is written out.

    autocenter : bool, optional
        Flag which, if set to True and there's a single chain will center the protein in the box.
        This is useful for visualization purposes but does mean any translational diffusion will
        be lost. Default = False

    unwrap : bool, optional
        Flag which, if True, writes each chain as a single whole periodic image
        (bond-walked, so coordinates may fall outside the box). Ignored where
        autocenter applies. Default is False.


    Returns
    -----------
    TrajectoryAccumulator
        The updated accumulator (see :class:`TrajectoryAccumulator`) - pass it back in
        on the next call, and hand it to :func:`save_out_sim` at the end. Each frame
        is buffered as one float32 coordinate array; nothing is joined or copied,
        so the cost is linear in the number of frames.

    Raises
    -----------
    LatticeUtilsException
        If master_traj is neither None nor a TrajectoryAccumulator.

    """
    # the frame, exactly as the streamed writer would have written it
    xyz, box = _lattice_frame_xyz_and_box(lattice, spacing, autocenter=autocenter, unwrap=unwrap)

    # (`is None`, not `== None`: an mdtraj Trajectory defines __eq__, so `== None`
    # is not guaranteed to be the identity test that is meant here.)
    if master_traj is None:
        master_traj = TrajectoryAccumulator(pdb_filename)

    elif not isinstance(master_traj, TrajectoryAccumulator):
        raise LatticeUtilsException(
            'update_master_traj needs None or a TrajectoryAccumulator as the master trajectory, '
            'but was passed a %s' % type(master_traj).__name__)

    # buffer this frame. The buffered frames are written out once, in save_out_sim.
    return master_traj.append_frame(xyz, box=box)


#-----------------------------------------------------------------
#
def start_master_traj(pdb_filename):
    """
    Create an EMPTY SAVE_AT_END accumulator, for a trajectory that holds only
    the frame written when it was started.

    ``update_master_traj`` both initialises the accumulator from the PDB *and*
    appends the current lattice as a new frame. When a run buffered no frame
    at all (``XTC_FREQ`` larger than the run, or ``SAVE_EQ : False`` with no
    production step divisible by ``XTC_FREQ``) the end-of-run code used that
    call just to have something to save, and so wrote a second frame holding the
    final state, stamped as frame 1. Under the documented ``frame * XTC_FREQ =
    step`` convention that frame described a step that never happened; the
    incremental writer produced one frame for the same run. The resized-
    equilibration path had the same defect and, worse, appended the *production*
    lattice (already swapped in) to ``eq_traj.xtc`` under the equilibration unit
    cell. Use this to obtain the frame-0-only trajectory instead.

    Parameters
    ----------
    pdb_filename : str
        The START.pdb (or eq_START.pdb) written when the trajectory was opened.

    Returns
    -------
    TrajectoryAccumulator
        An accumulator with nothing buffered, ready for :func:`save_out_sim`
        (which then leaves the frame-0-only file written by
        :func:`start_xtc_file` as it is). The PDB is not read here; it is only
        read if the trajectory file has gone missing and frame 0 has to be
        rebuilt.
    """
    return TrajectoryAccumulator(pdb_filename)


#-----------------------------------------------------------------
#
def save_out_sim(master_traj, xtc_filename):
    """
    Save out the master trajectory.

    For a :class:`TrajectoryAccumulator` this writes every frame buffered since
    the last call onto the end of ``xtc_filename`` and empties the buffer. It is
    therefore safe to call more than once, and safe to call from an error path:
    frames already on disk are never written again, and a call with nothing
    buffered leaves the file untouched. This is what lets a ``SAVE_AT_END`` run
    that fails keep the frames it had saved up to the failure.

    Parameters
    ------------
    master_traj : TrajectoryAccumulator or mdtraj.Trajectory
        The trajectory built up through the sim. A
        :class:`TrajectoryAccumulator` (what :func:`update_master_traj` returns)
        has its buffered frames added to the file; a plain Trajectory is saved
        as is, replacing the file.

    xtc_filename : str
        Filename to write to disk.

    Returns
    ------------
    None
        No return, but the XTC file on disk holds every frame saved so far.

    Raises
    ------
    LatticeUtilsException
        If the trajectory file is missing and its first frame cannot be
        rebuilt from the topology PDB.

    """
    if isinstance(master_traj, TrajectoryAccumulator):
        master_traj.write_out(xtc_filename)
        return

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
        Chain ID number, used only in the printed/raised messages

    chain_positions : list of lists
        List of the chain's bead positions in chain order, each a 2- or
        3-element coordinate

    dimensions : list
        List of the box dimensions (2 or 3 ints), used to test whether an
        apparent break is really a periodic wrap

    verbose : bool, optional
        Flag to print out debug information. Default = True

    Returns
    ------------
    None
        No return, but will raise an error if the chain is not connected

    Raises
    ------------
    ChainConnectivityError
        If two consecutive beads are more than one lattice site apart, even
        after allowing for a periodic wrap.

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
        Dictionary mapping chainID to Chain object (i.e. LATTICE.chains),
        despite the name

    dimensions : list
        List of the box dimensions (2 or 3 ints)

    verbose : bool, optional
        Flag to print out debug information. Default = True


    Returns
    ------------
    None
        No return, but will raise an error if any chain is not connected

    Raises
    ------------
    ChainConnectivityError
        If any chain has two consecutive beads more than one lattice site
        apart, even after allowing for a periodic wrap.

    """
    
    for chainID in list_of_chain_objects:
        check_chain_connectivity(chainID, list_of_chain_objects[chainID].get_ordered_positions(), dimensions, verbose=verbose)

