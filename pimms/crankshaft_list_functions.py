## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................


# Functions for creating and updating the crankshaft_list
#
#
#
import numpy as np
from pimms.latticeExceptions import MoveException


# -----------------------------------------------------------------
#
#
def __get_bead_flag(bead_idx, chain_length):
    
    """
    Internal function that returns the bead-flag associated with a bead position
    in a chain (as in, relative location along the chain) given the chain length.
    This is used for the system_shake moves, in which the bead position is needed
    to evaluate angle changes correctly, as well as call the correct changes in
    position (i.e. need to know if a bead is at the end of the chain or inside a chain).

    See the table below for a mapping of the beadflag to the relevant offset for angles.
    Offset here reflects the max and min offset value that can be applied FROM the selected 
    bead that should be used when calculating 3-residue sub-chains which need their angles 
    evaluated upon bead movement. 
        
    L=3
    0XX 1   0,3
    X0X 4  -1,2
    XX0 3  -2,1
    
    L=4
    0XXX 1  0,3
    X0XX 5  -1,3
    XX0X 6  -2,2
    XXX0 3  -2,1
    
    L=5
    0XXXX 1  0,3
    X0XXX 5 -1,3
    XX0XX 2 -2,3
    XXX0X 6  -2,2
    XXXX0 3  -2,1
    
    L = 6 or greater
    all_others_polymer lenths
    N    --> 1
    N+1  --> 5        
    INTERNAL --> 2
    C-1  --> 6
    C    --> 3

    Parameters
    ----------
    bead_idx : int
        The index of the bead in the chain (i.e. rank position
        along the chain length).

    chain_length : int
        The length of the chain.

    Returns
    -------
    bead_flag : int
        The bead flag associated with the bead position in the chain.
        
    
    """

    # polymer length 1
    if chain_length == 1:
        return 0


    elif chain_length == 2:
        if bead_idx   == 0: 
            return 1                  # start
        else:
            return 3                  # end

    # polymer length 3
    elif chain_length == 3:
        if bead_idx   == 0:           # start 
            return 1
        elif bead_idx == 1:           # start +1
            return 4
        else:                         # end
            return 3

    # polymer length 4
    elif chain_length == 4:
        if bead_idx   == 0:           # start
            return 1
        elif bead_idx == 1:           # start + 1
            return 5
        elif bead_idx == 2:           # end - 1
            return 6
        else:                         # end 
            return 3                  

    # polymer length 5
    elif chain_length == 5:
        if bead_idx   == 0:           # start
            return 1
        elif bead_idx == 1:             # start +1
            return 5
        elif bead_idx == 2:             # start + 2
            return 2
        elif bead_idx == 3:             # end - 1
            return 6                  
        else:
            return 3                    # end

    # polymer length over 6 
    else:
        if bead_idx   == 0:              # start
            return 1
        elif bead_idx == 1:              # start + 1
            return 5
        elif bead_idx == chain_length-2: # end - 1
            return 6
        elif bead_idx == chain_length-1: # end  
            return 3
        else:
            return 2


# -----------------------------------------------------------------
#
#
def __single_chain_idx_to_bead(chainID, latticeObject):
    """
    Function that constructs an idx_to_bead array that can be fed into megacrank functions. This function
    builds the the idx_to_bead information from scratch, and should only be called when the latticeObject
    is initialized. Calling it more often will incurr a totally unnecessary penalty, BUT in case we want to
    add non-equilibrium effects later, this function would let you fully reset and update the idx_to_bead
    information.

    Parameters
    ----------
    chainID : int
        The chainID of the chain that we want to construct the idx_to_bead array for.

    latticeObject : lattice
        The lattice object that we want to construct the idx_to_bead array for.

    Returns
    -------
    idx_to_bead : list
        A list of lists that contains the idx_to_bead information for the chain of interest. The idx_to_bead
        information is a list of lists that contains the following information for each bead in the chain:
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y, z]. Note, none of these numbers should be
        very big - i.e. they scale with box dimensions or number of unique beads, but none scale with absolute
        number of beads in the system.



    """

    idx_to_bead = []

    # extract the chain of interest
    c = latticeObject.chains[chainID]

    # get local positions, chain length, and LR and intecode arrays
    local_pos = c.get_ordered_positions()

    # get chain length
    chain_length  = len(local_pos)

    # get LR binary array (i.e. where each bead engages in LR interactions
    # or not)
    local_LR_binary_array = c.get_LR_binary_array()

    # get the intcode sequence for the chain. intcode is a list of integers
    # that represent the bead identity as encoded by the integer to bead
    # type mapping build by the parameter input file.
    local_intcode_seq = c.get_intcode_sequence()

    ## Construct the single chain idx_to_bead matrix which is 
    if chain_length == 1:
        temp = []
        temp.append(0)
        temp.append(local_LR_binary_array[0])
        temp.append(local_intcode_seq[0])
        temp.append(1)                             # skip angles = True 
        temp.append(chainID) 
        temp.extend(local_pos[0])
        idx_to_bead.append(temp)

    # else on a polymer of length 2 or more beads
    else:

        # set all bead flags 
        for p in range(0, chain_length):
            temp = []

            temp.append(__get_bead_flag(p,chain_length))
            temp.append(local_LR_binary_array[p])
            temp.append(local_intcode_seq[p])

            # skip angles if chain_length is 2
            if chain_length == 2:
                temp.append(1)                             # skip angles = False 
            else:
                temp.append(0)                             # skip angles = True 

            temp.append(chainID)
            temp.extend(local_pos[p])
            idx_to_bead.append(temp)

    return idx_to_bead




# -----------------------------------------------------------------
#
#
def initialize_idx_to_bead(latticeObject):
    """
    Function that constructs a new idx_to_bead matrix using the chain information
    from the passed lattice object. This function DOES NOT edit the latticeObject. The position
    elements are set to whatever the positions are at this moment, but those values are really
    not meant to be used but are basically palceholders that get overwritten.


    # Each bead contains the following information (index position included)

    # 0 - bead_flag
    # 1 - LR binary flag
    # 2 - intcode value
    # 3 - skip angles (1 = true, 0 = false)
    # 4 - chainID
    # 5 - X position
    # 6 - Y position
    # 7 - Z position (optional - depends on if we're in 3D or not)


    # Bead Flags
    # we have six types of flags, which we asign to each bead according to its relative position
    # in a chain. Note the code below has the nice property of working in both two and three dimensions. 
    # The flags used are shown below, and are described in more specific detail in the 
    # __get_bead_flag() function
    # 
    #
    # 0 single bead
    # 1 N-terminal bead
    # 2 central bead with residues +2/-2 around
    # 3 C-terminal bead
    # 4 Central bead in a polymer of L=3  
    # 5 N-termnial +1 bead 
    # 6 C-termina -1 bead
    #
    # With these 6 options you can fully describe all possible bead positions to capture
    # angle effects 

    Parameters
    ----------

    latticeObject : lattice
        The lattice object that we want to construct the idx_to_bead array for.

    Returns
    -------
    idx_to_bead : numpy array
        A numpy array that contains the idx_to_bead information for the entire system. The idx_to_bead
        information is a list of lists that contains the following information for each bead in the chain:
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y, z]. 

    """

    # for each chain, if we have not yet initialized the crankshaft_list,
    # do so now
    
    idx_to_bead = []
    
    # cycle over each chainID in order
    for chainID in sorted(latticeObject.chains.keys()):
        tmp = __single_chain_idx_to_bead(chainID, latticeObject)
        idx_to_bead.extend(tmp)

    # finally convert to numpy array (the format of the crankshaft_list matrix
    return np.array(idx_to_bead, dtype=np.int64)



# -----------------------------------------------------------------
#
#
def initialize_chain_to_firstbead_lookup(latticeObject):
    """
    Function that constructs a dictionary that maps the chainID to the index of the first bead in the 
    chain from the perspective of the idx_to_bead array. This is useful for quickly looking up the 
    first bead in a chain, which is needed for the crankshaft move acceptance criteria.

    Recall that the idx_to_array bead is a 7 x n array where n is the number of beads in the 
    system. The rows are ordered in terms of ordered beads in each chain, ordered by chainID.
    With that in mind, the chain_to_firstbead_lookup is a dictionary that maps the chainID to the
    index associated with the row in the idx_to_bead array that corresponds to the first bead in
    the chain. 

    Parameters
    ----------
    latticeObject : lattice
        The lattice object that we want to construct the chain_to_firstbead_lookup for.

    Returns
    -------
    chain_to_firstbead_lookup : dictionary
        A dictionary that maps the chainID to the first bead in the chain. This is useful for quickly
        looking up the first bead in a chain, which is needed for the crankshaft move acceptance criteria.


    """
    chain_to_firstbead_lookup = {}

    bead = 0
    for chainID in sorted(latticeObject.chains.keys()):
        chain_to_firstbead_lookup[chainID] = bead 
        bead = bead+len(latticeObject.chains[chainID].positions)

    return chain_to_firstbead_lookup




# -----------------------------------------------------------------
#
#
def update_idx_to_bead(latticeObject):
    """
    Function that updates the crankshaft_lists object such that the 
    positions are set to the lattice' current positional state. This 
    function DOES update the latticeObjects.crankshaft_lists, but also 
    returns the full crankshaft_list.

    As a reminder, the latticeObject has an object called crankshaft_lists
    which is the idx_to_bead matrix. This matrix is a 7 x n array where n 
    is the number of beads in the system. The columns are 

    # 0 - bead_flag
    # 1 - LR binary flag
    # 2 - intcode value
    # 3 - skip angles (1 = true, 0 = false)
    # 4 - chainID
    # 5 - X position
    # 6 - Y position
    # 7 - Z position (optional - depends on if we're in 3D or not)

    So the code here updates the positions in the crankshaft_lists matrix
    to the current positions in the lattice object.

    Parameters
    ----------
    latticeObject : lattice
        The lattice object that we want to update the idx_to_bead array for.

    Returns
    -------
    idx_to_bead : numpy array
        A numpy array that contains the idx_to_bead information for the entire system. The idx_to_bead
        information is a list of lists that contains the following information for each bead in the chain:
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y, z].

    """

    # cycle over each chain and update the positions in the crankshaft_lists matrix. This then gets passed into
    # the megacrank function
    local_idx = 0
    for chainID in sorted(latticeObject.chains.keys()):
        
        # extract the positions, convert to a numpy array, and assign to the appropiate positions in the crankshaft_lists matrix
        pos_list = latticeObject.chains[chainID].get_ordered_positions()

        
        latticeObject.crankshaft_lists[local_idx:local_idx+len(pos_list),5:] = np.array(pos_list)
        local_idx = local_idx + len(pos_list)

    idx_to_bead = latticeObject.crankshaft_lists

    # UP-2023-5 updated to np.array cast here so these alwayts return a np.array - note we define
    # the dtype explicitly so this can be passed into cython without issue
    return np.array(idx_to_bead, dtype=np.int64)


#-----------------------------------------------------------------
#
#
def update_idx_to_bead_single_chain(latticeObject, chainID):
    """
    Function that updates the crankshaft_lists object such that the positions are set to the lattice' current
    positional state. This function DOES update the latticeObjects.crankshaft_lists, but also returns
    the full crankshaft_list. 

    This returns a subset of the idx_to_bead matrix for JUST a single chain, which can then be passed to
    a megacrank function.


    Parameters
    ----------
    latticeObject : lattice
        The lattice object that we want to update the idx_to_bead array for.

    chainID : int
        The chainID that we want to update the idx_to_bead array for.

    Returns
    -------
    idx_to_bead : numpy array
        A numpy array that contains the idx_to_bead information for a single chain. The idx_to_bead
        information is a list of lists that contains the following information for each bead in the chain:
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y, z].

    """
    
    # get current positions from the chain object
    pos_list = latticeObject.chains[chainID].get_ordered_positions()

    # get the index of the position in the crankshaft_lists matrix
    local_idx = latticeObject.chain_to_firstbead_lookup[chainID]

    # update the crankshaft_lists matrix using the lattice positions
    latticeObject.crankshaft_lists[local_idx:local_idx+len(pos_list),5:] = np.array(pos_list)

    # extract the idx_to_bead information for the chain
    idx_to_bead = latticeObject.crankshaft_lists[local_idx:local_idx+len(pos_list),:]

    # UP-2023-5 updated to np.array cast here so these alwayts return a np.array - note we define
    # the dtype explicitly so this can be passed into cython without issue. 
    return np.array(idx_to_bead, dtype=np.int64)


# -----------------------------------------------------------------
#
#
def update_idx_to_bead_multiple_chains(latticeObject, chain_list):
    """
    Function that updates the crankshaft_lists object such that the positions are set to the lattice' current
    positional state. This function DOES update the latticeObjects.crankshaft_lists, but also returns
    the full crankshaft_list. 

    This returns a subset of the idx_to_bead matrix for multiple chains, as defined in the chain_list.

    Parameters
    ----------
    latticeObject : lattice
        The lattice object that we want to update the idx_to_bead array for.

    chain_list : list
        A list of chainIDs that we want to update the idx_to_bead array for.

    Returns
    -------
    idx_to_bead : numpy array
        A numpy array that contains the idx_to_bead information for multiple chains. The idx_to_bead
        information is a list of lists that contains the following information for each bead in the chain:
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y, z].

    
    """

    idx_to_bead = []
    for chainID in chain_list:
        

        # get current positions from the chain object
        pos_list = latticeObject.chains[chainID].get_ordered_positions()
    
        # get the index of the position in the crankshaft_lists matrix
        local_idx = latticeObject.chain_to_firstbead_lookup[chainID]

        # update the crankshaft_lists matrix using the lattice positions
        latticeObject.crankshaft_lists[local_idx:local_idx+len(pos_list),5:] = np.array(pos_list)

        # extract the idx_to_bead information for the chain
        idx_to_bead.extend(latticeObject.crankshaft_lists[local_idx:local_idx+len(pos_list),:].tolist())
        


    return np.array(idx_to_bead, dtype=np.int64)



# -----------------------------------------------------------------
#
#

def bead_selector_constructor(num_beads, number_of_steps, latticeObject, chain_override_list=[], safecheck=False):
    """
    Function that returns a list of bead indices that we want to attempt to move. 
    By default this randomly samples all possible beads on the lattice, but we can
    restrain specific chains using the chain_override_list.

    Parameters
    ----------
    num_beads : int
        The total number of beads in the lattice.

    number_of_steps : int
        The number of steps that we want to attempt to move beads.

    latticeObject : lattice
        The lattice object that we want to update the idx_to_bead array for; this
        is used to extract the chain information if an override is passed, or to
        sanity check the number of beads if the safecheck flag is set.

    chain_override_list : list
        A list of chainIDs that we want to sample from. If this is empty then we
        sample from all possible beads.

    safecheck : bool
        A flag that is used to check that the number of beads in the lattice object
        matches the number of beads in the idx_to_bead matrix. This is a safety check
        to ensure that the idx_to_bead matrix is not corrupted.

    Returns
    -------
    bead_indices : numpy array
        A numpy array that contains the indices of the beads that we want to attempt
        to move in a random order. This essentially defines a random order in which
        we want to attempt to move beads.

    """


    if safecheck:
        beadcount = 0
        for chainID in latticeObject.chains:
            beadcount = beadcount + len(latticeObject.chains[chainID])

        if beadcount != num_beads:
            raise MoveException("The number of beads in the lattice object does not match the number of beads in the idx_to_bead matrix. This is a bug")
            

        
    # if override_list is empty then we are randomly sampling from all possible beads 
    if len(chain_override_list) == 0:
        return np.random.randint(0,num_beads,number_of_steps)

