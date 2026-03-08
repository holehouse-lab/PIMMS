## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................



import numpy as np
import copy

from collections import deque

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
    Function which determines the set of non-redundant pairs that describe the 
    interface of the positions in the list $positions with the solvent. The 
    KEY assumption here is that positions FULLY describes all the positions
    of occupied lattice sites in a connected component. This method WILL 
    BREAK if that assumption is wrong.

    This is the reason this method is off in its own set of utilities - 
    if we're not working with a complete (redundant) set of positions that 
    describe a single connected component (cluster) everything is gonna 
    go wrong...

    But assuming that's true, this method is faster than a shark on cocaine.

    """

    n_dim = len(dimensions)

    if n_dim not in (2, 3):
        raise ValueError(f"Unsupported dimensionality for interface envelope calculation: {n_dim}")

    if len(positions) == 0:
        return np.empty((0, 4 if n_dim == 2 else 6), dtype=int)

    pair_chunks = []
    for position in positions:
        if n_dim == 2:
            tmp = np.array(
                hyperloop.get_unique_interface_pairs_2D(
                    position[0],
                    position[1],
                    dimensions[0],
                    dimensions[1],
                    grid,
                )
            )
        else:
            tmp = np.array(
                hyperloop.get_unique_interface_pairs_3D(
                    position[0],
                    position[1],
                    position[2],
                    dimensions[0],
                    dimensions[1],
                    dimensions[2],
                    grid,
                )
            )

        if tmp.shape[0] > 0:
            pair_chunks.append(tmp)

    if len(pair_chunks) == 0:
        return np.empty((0, 4 if n_dim == 2 else 6), dtype=int)

    return np.vstack(pair_chunks)

    
def build_interface_envelope_pairs_safe_and_slow(positions, dimensions):    
    """
    This provides an alternative implementation of the function above, but instead users the softere assumption
    that interface pairs must be those which are pairs and it's not true that both members of the pair come 
    from within the list of positions supplied.
    
    However, again it's crucial that positions be an exhaustative list of the positions making up a connected 
    component, although this time we don't use any lattice information.

    """
    envelope_pairs = lattice_utils.build_envelope_pairs(positions, dimensions)
    
    # for each pair ask if both members of the pair come from the list of positions    
    interface_pairs = []

    np_positions = np.array(positions)

    for pair in envelope_pairs:
        if (~(np_positions-pair[0]).any(axis=1)).any() and (~(np_positions-pair[1]).any(axis=1)).any():
            pass
        else:
            interface_pairs.append(pair)

    num_pairs = len(interface_pairs)

    #print "Saved vs. looped - RATIO: [%s]" %(str(float(len(envelope_pairs)-num_pairs)/len(envelope_pairs)))
    if len(dimensions) == 2:
        if num_pairs == 0:
            return np.empty((0, 2, 2), dtype=int)
        reshaped = np.reshape(np.array(interface_pairs, dtype=int), (num_pairs, 2, 2))
    elif len(dimensions) == 3:
        if num_pairs == 0:
            return np.empty((0, 2, 3), dtype=int)
        reshaped = np.reshape(np.array(interface_pairs, dtype=int), (num_pairs, 2, 3))
    else:
        raise ValueError(f"Unsupported dimensionality for interface envelope calculation: {len(dimensions)}")

    return reshaped


#-----------------------------------------------------------------
#    
def convert_positions_to_single_image_snakesearch(original_positions, dimensions, space_threshold=1):
    """
    Reconstruct a list of positions so that they all lie in a single periodic image.

    Uses a BFS (breadth-first search) walk starting from the bead closest to the
    periodic-aware center of mass. For each visited bead, all unvisited neighbours
    within *space_threshold* (per-dimension, PBC-aware) are corrected into the same
    periodic image as the reference bead and enqueued.

    The algorithm is O(N * K) where K is the typical number of neighbours per bead
    (bounded by (2*space_threshold+1)^n_dim), compared to the original O(N^2)
    implementation.

    Parameters
    ----------
    original_positions : list of lists
        Each sublist is the [x, y] or [x, y, z] position of a bead.

    dimensions : list of int
        Lattice dimensions in each axis.

    space_threshold : int
        Maximum per-dimension distance for two beads to be considered connected.

    Returns
    -------
    numpy.ndarray
        (N, n_dim) array of positions, all in the same periodic image, shifted
        so that every coordinate is >= 0.
    """

    if len(original_positions) == 0:
        return []

    n_dim = len(dimensions)
    dims = np.array(dimensions, dtype=np.int64)
    half_dims = dims / 2.0

    # ---- build a fast lookup from PBC-tuple -> original index ----
    N = len(original_positions)
    pos_arr = np.array(original_positions, dtype=np.int64)  # (N, n_dim)

    # Map each PBC-space tuple to its index for O(1) neighbour discovery
    pbc_to_idx = {}
    for k in range(N):
        key = tuple(pos_arr[k])
        # handle rare duplicate positions (shouldn't happen for valid clusters)
        if key not in pbc_to_idx:
            pbc_to_idx[key] = []
        pbc_to_idx[key].append(k)

    # ---- find the seed: bead closest to the PBC-aware COM ----
    COM = lattice_utils.center_of_mass_from_positions(original_positions, dimensions)
    seed_idx = lattice_utils.find_nearest_position(COM, original_positions, dimensions)[0]

    # ---- BFS to build single-image positions ----
    si_positions = np.empty_like(pos_arr)  # output array
    si_positions[seed_idx] = pos_arr[seed_idx]

    visited = np.zeros(N, dtype=np.bool_)
    visited[seed_idx] = True

    queue = deque()
    queue.append(seed_idx)

    # precompute the offset grid for neighbour search
    # e.g. for space_threshold=1 in 2D: offsets = [-1,0,1] x [-1,0,1]
    offset_ranges = [range(-space_threshold, space_threshold + 1) for _ in range(n_dim)]
    if n_dim == 2:
        offsets = [(dx, dy) for dx in offset_ranges[0] for dy in offset_ranges[1]]
    else:
        offsets = [(dx, dy, dz) for dx in offset_ranges[0] for dy in offset_ranges[1] for dz in offset_ranges[2]]

    found_count = 1
    while queue:
        ref_idx = queue.popleft()
        ref_si = si_positions[ref_idx]  # this bead's single-image position

        # Convert the reference single-image position back to PBC space for
        # neighbour lookup
        ref_pbc = ref_si % dims

        for off in offsets:
            # candidate neighbour in PBC space
            nbr_pbc = tuple((ref_pbc[d] + off[d]) % dims[d] for d in range(n_dim))
            candidates = pbc_to_idx.get(nbr_pbc)
            if candidates is None:
                continue

            for nbr_idx in candidates:
                if visited[nbr_idx]:
                    continue
                visited[nbr_idx] = True
                found_count += 1

                # Place this neighbour in the same periodic image as the reference.
                # The single-image position = ref_si + offset
                # (offset is the PBC-minimal displacement from ref to neighbour)
                nbr_orig = pos_arr[nbr_idx]
                delta = nbr_orig - ref_pbc
                # PBC-correct each dimension so delta is in [-half, +half)
                for d in range(n_dim):
                    if delta[d] > half_dims[d]:
                        delta[d] -= dims[d]
                    elif delta[d] < -half_dims[d]:
                        delta[d] += dims[d]
                si_positions[nbr_idx] = ref_si + delta

                queue.append(nbr_idx)

    if found_count != N:
        raise ValueError(
            "Input positions must form a connected cluster within the provided space_threshold"
        )

    # ---- shift so all coordinates >= 0 ----
    for d in range(n_dim):
        col = si_positions[:, d]
        min_val = col.min()
        if min_val < 0:
            col -= min_val  # in-place shift

    return si_positions
        

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
                

