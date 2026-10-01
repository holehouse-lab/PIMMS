## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


# Functions for creating and updating the crankshaft_list
#
#
#
import numpy as np
import collections

# The list-of-lists <-> table copies below are the cost of every megamove at
# scale (they were about 4 ms of a 6 ms crankshaft megamove at 12,000 beads when
# written in Python). The compiled loops in pimms.bookkeeping do the same copies
# without creating a numpy array or a list per chain; the pure-Python versions
# are kept as the fallback when the extension is not built, and as the oracle the
# compiled ones are tested against.
try:
    from . import bookkeeping as _bookkeeping
    _HAVE_BOOKKEEPING = True
except ImportError:  # pragma: no cover - only when the extension is not built
    _bookkeeping = None
    _HAVE_BOOKKEEPING = False
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
        along the chain length), running from 0 to ``chain_length - 1``.

    chain_length : int
        The number of beads in the chain.

    Returns
    -------
    int
        The bead flag associated with the bead position in the chain: 0 for a
        single-bead chain, 1 for the N-terminal bead, 2 for a central bead with
        two neighbours on each side, 3 for the C-terminal bead, 4 for the central
        bead of a 3-mer, 5 for the N-terminal+1 bead and 6 for the C-terminal-1
        bead.


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

    latticeObject : Lattice
        The lattice object that we want to construct the idx_to_bead array for.

    Returns
    -------
    list of list of int
        One row per bead of the chain, in chain order, each row being
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y] in 2D or
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y, z] in 3D.
        Note, none of these numbers should be very big - i.e. they scale with box
        dimensions or number of unique beads, but none scale with absolute number
        of beads in the system.



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

            # skip angles if chain_length is 2 (a 2-bead chain has no angle)
            if chain_length == 2:
                temp.append(1)                             # skip angles = True
            else:
                temp.append(0)                             # skip angles = False

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

    latticeObject : Lattice
        The lattice object that we want to construct the idx_to_bead array for.

    Returns
    -------
    numpy.ndarray
        An ``int64`` array of shape ``(num_beads, 7)`` in 2D or
        ``(num_beads, 8)`` in 3D, holding one row per bead in the system ordered
        by ascending chainID and then chain position. Each row is
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y(, z)].

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



ChainLayout = collections.namedtuple(
    'ChainLayout', ['sorted_ids', 'offset', 'length', 'homo', 'offset32', 'length32', 'n_beads'])
"""Everything about a lattice's chains that never changes during a run.

``sorted_ids`` is the ascending list of chainIDs, which is the order of the rows
in the bead table; ``offset`` and ``length`` (``int64``) give each chain's first
row and row count in that order, with ``offset32`` / ``length32`` the same values
as the ``int32`` arrays the whole-chain kernels take; ``homo`` (``int32``) is 1
for a chain whose beads all share one intcode AND one long-range flag, which is
the precondition of the kernels' O(1) homopolymer energy path; ``n_beads`` is the
total row count. Positions are deliberately not here: they change every move.
"""


def initialize_chain_layout(latticeObject):
    """
    Build the static chain layout of a lattice from its bead table.

    No move ever adds, removes or re-sequences a chain, so the row each chain
    occupies in the bead table, its length and whether it is homopolymeric are
    fixed for the run. Working them out once here removes a per-chain Python loop
    from every megamove; ``update_idx_to_bead`` and the whole-chain megamoves
    used to rebuild all of this every time they were called.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose chains are described. Its ``crankshaft_lists`` table
        must already be initialised, since the homopolymer flags are read from
        its intcode and long-range columns.

    Returns
    -------
    ChainLayout
        The layout, in ascending chainID order.
    """
    sorted_ids = sorted(latticeObject.chains.keys())
    lengths = np.array([len(latticeObject.chains[c].positions) for c in sorted_ids],
                       dtype=np.int64)
    n_chains = len(sorted_ids)
    offsets = np.zeros(n_chains, dtype=np.int64)
    if n_chains > 1:
        offsets[1:] = np.cumsum(lengths)[:-1]
    n_beads = int(lengths.sum()) if n_chains else 0

    table = np.asarray(latticeObject.crankshaft_lists)
    if n_chains and table.ndim == 2 and table.shape[0] == n_beads:
        # a chain is homopolymeric when its intcode (column 2) and its long-range
        # flag (column 1) are constant over its rows: the kernels read bead 0's
        # values for the whole chain on their fast path, so both must be uniform.
        # max == min over each chain's rows says exactly that.
        intcode = table[:, 2]
        lr_flag = table[:, 1]
        homo = ((np.maximum.reduceat(intcode, offsets) == np.minimum.reduceat(intcode, offsets))
                & (np.maximum.reduceat(lr_flag, offsets) == np.minimum.reduceat(lr_flag, offsets)))
        homo = homo.astype(np.int32)
    else:
        homo = np.zeros(n_chains, dtype=np.int32)

    return ChainLayout(sorted_ids=sorted_ids,
                       offset=offsets, length=lengths, homo=homo,
                       offset32=offsets.astype(np.int32), length32=lengths.astype(np.int32),
                       n_beads=n_beads)


def chain_layout(latticeObject):
    """
    The lattice's static chain layout, built on first use and cached on it.

    The cache is checked against the chain count and the table size on every
    call, so a lattice whose chains were replaced wholesale (which nothing in
    PIMMS does after construction, but a test might) gets a fresh layout rather
    than a stale one. Only chainIDs are cached, never Chain objects, because the
    system-wide TSMMC restore writes positions into the existing objects and the
    layout must not care either way.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose layout is wanted.

    Returns
    -------
    ChainLayout
        See :data:`ChainLayout`.
    """
    layout = getattr(latticeObject, 'chain_layout', None)
    table = latticeObject.crankshaft_lists
    n_rows = len(table) if hasattr(table, '__len__') else 0
    if (layout is None or len(layout.sorted_ids) != len(latticeObject.chains)
            or layout.n_beads != n_rows):
        layout = initialize_chain_layout(latticeObject)
        latticeObject.chain_layout = layout
    return layout


def _chains_in_table_order(latticeObject, layout):
    """The Chain objects in the order of their rows in the bead table.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose chains are wanted.

    layout : ChainLayout
        The cached per-chain layout; its ``sorted_ids`` give the row order.

    Returns
    -------
    list of Chain
        One Chain per chain, in bead-table order.
    """
    chains = latticeObject.chains
    return [chains[c] for c in layout.sorted_ids]


# -----------------------------------------------------------------
#
#
def initialize_chain_to_firstbead_lookup(latticeObject):
    """
    Function that constructs a dictionary that maps the chainID to the index of the first bead in the 
    chain from the perspective of the idx_to_bead array. This is useful for quickly looking up the 
    first bead in a chain, which is needed for the crankshaft move acceptance criteria.

    Recall that the idx_to_bead array has one row per bead (n rows) and 7 columns
    in 2D or 8 in 3D. The rows are ordered in terms of ordered beads in each
    chain, ordered by chainID.
    With that in mind, the chain_to_firstbead_lookup is a dictionary that maps the chainID to the
    index associated with the row in the idx_to_bead array that corresponds to the first bead in
    the chain. 

    Parameters
    ----------
    latticeObject : Lattice
        The lattice object that we want to construct the chain_to_firstbead_lookup for.

    Returns
    -------
    dict
        Maps each chainID to the row index in the idx_to_bead array at which that
        chain's first bead sits. This is useful for quickly looking up the first
        bead in a chain, which is needed for the crankshaft move acceptance
        criteria.


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
    Refresh the position columns of the bead table from the chains, and return a copy.

    The lattice keeps one bead table, ``latticeObject.crankshaft_lists``: one row
    per bead in ascending chainID order, with the static columns
    ``[bead_flag, LR_binary, intcode, skip_angles, chainID]`` followed by the
    bead's ``x, y(, z)``. Only the position columns can be out of date, since the
    Python moves change chain positions without touching the table, so this
    copies every chain's current positions into columns 5 on and returns a fresh
    ``int64`` copy of the whole table for a kernel to work on.

    The copy is done by the compiled :mod:`pimms.bookkeeping` loop when it is
    available. In pure Python it was one ``np.array(positions)`` per chain and
    dominated the cost of a megamove at scale; see that module.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose chains supply the positions. Its ``crankshaft_lists``
        is updated in place as a side effect.

    Returns
    -------
    numpy.ndarray
        A fresh ``int64`` array of shape ``(num_beads, 7)`` in 2D or
        ``(num_beads, 8)`` in 3D, one row per bead in ascending chainID order and
        then chain position, each row
        ``[bead_flag, LR_binary, intcode, skip_angles, chainID, x, y(, z)]``.
    """
    layout = chain_layout(latticeObject)
    table = latticeObject.crankshaft_lists
    if _HAVE_BOOKKEEPING and isinstance(table, np.ndarray) and table.dtype == np.int64 \
            and table.flags['C_CONTIGUOUS'] and table.ndim == 2:
        _bookkeeping.gather_positions(table, _chains_in_table_order(latticeObject, layout),
                                      layout.offset, layout.length, table.shape[1] - 5)
    else:
        _update_idx_to_bead_python(latticeObject)
    return np.array(latticeObject.crankshaft_lists, dtype=np.int64)


def _update_idx_to_bead_python(latticeObject):
    """
    The pure-Python position refresh: the fallback, and the oracle for the compiled one.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose ``crankshaft_lists`` position columns are refreshed
        from its chains, in place.

    Returns
    -------
    numpy.ndarray
        A fresh ``int64`` copy of the refreshed table.
    """
    local_idx = 0
    for chainID in sorted(latticeObject.chains.keys()):
        pos_list = latticeObject.chains[chainID].get_ordered_positions()
        latticeObject.crankshaft_lists[local_idx:local_idx + len(pos_list), 5:] = np.array(pos_list)
        local_idx = local_idx + len(pos_list)
    return np.array(latticeObject.crankshaft_lists, dtype=np.int64)


def write_back_positions(latticeObject, idx_to_bead):
    """
    Copy the positions a kernel left in a bead table back into the chains.

    This is the other half of :func:`update_idx_to_bead`: the kernels move beads
    by editing the table (and the grids) in place, and the Chain objects, which
    everything else in PIMMS reads, have to be told afterwards. Every chain gets a
    fresh list of ``[x, y(, z)]`` lists through ``set_ordered_positions``, which
    keeps that method's length check.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose chains are updated.

    idx_to_bead : numpy.ndarray
        The table a kernel returned, one row per bead in ascending chainID order,
        positions from column 5 on.

    Returns
    -------
    None
        The chains are updated in place.
    """
    layout = chain_layout(latticeObject)
    table = np.ascontiguousarray(idx_to_bead, dtype=np.int64)
    if _HAVE_BOOKKEEPING and table.ndim == 2:
        _bookkeeping.scatter_positions(table, _chains_in_table_order(latticeObject, layout),
                                       layout.offset, layout.length, table.shape[1] - 5)
    else:
        _write_back_positions_python(latticeObject, table)


def _write_back_positions_python(latticeObject, idx_to_bead):
    """
    The pure-Python write-back: the fallback, and the oracle for the compiled one.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose chains are updated in place.

    idx_to_bead : numpy.ndarray
        The table a kernel returned, positions from column 5 on.

    Returns
    -------
    None
    """
    local_idx = 0
    for chainID in sorted(latticeObject.chains.keys()):
        n_pos = len(latticeObject.chains[chainID].get_ordered_positions())
        latticeObject.chains[chainID].set_ordered_positions(
            idx_to_bead[local_idx:local_idx + n_pos, 5:].tolist())
        local_idx = local_idx + n_pos


#-----------------------------------------------------------------
#
#
def update_idx_to_bead_single_chain(latticeObject, chainID):
    """
    Function that updates the crankshaft_lists object such that the positions are set to the lattice' current
    positional state. This function DOES update the latticeObjects.crankshaft_lists, and returns a copy
    of the rows it touched.

    This returns a subset of the idx_to_bead matrix for JUST a single chain, which can then be passed to
    a megacrank function.


    Parameters
    ----------
    latticeObject : Lattice
        The lattice object that we want to update the idx_to_bead array for.

    chainID : int
        The chainID that we want to update the idx_to_bead array for.

    Returns
    -------
    numpy.ndarray
        An ``int64`` array with one row per bead of the requested chain (shape
        ``(chain_length, 7)`` in 2D or ``(chain_length, 8)`` in 3D). Each row is
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y(, z)].

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
    positional state. This function DOES update the latticeObjects.crankshaft_lists, and returns a copy
    of the rows it touched.

    This returns a subset of the idx_to_bead matrix for multiple chains, as defined in the chain_list.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice object that we want to update the idx_to_bead array for.

    chain_list : sequence of int
        The chainIDs to build the idx_to_bead array for. The rows of the returned
        array follow this order, not ascending chainID order.

    Returns
    -------
    numpy.ndarray
        An ``int64`` array with one row per bead of the requested chains, blocked
        by chain in the order the chainIDs appear in ``chain_list``. Each row is
        [bead_flag, LR_binary, intcode, skip_angles, chainID, x, y(, z)].


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

def bead_selector_constructor(num_beads, number_of_steps, latticeObject, frozen_chains=(), safecheck=True):
    """
    Function that returns a list of bead indices that we want to attempt to move. 
    By default this randomly samples all possible beads on the lattice, but we can
    restrain specific chains using the frozen_chains list.

    Parameters
    ----------
    num_beads : int
        The total number of beads in the lattice, i.e. the number of rows in the
        idx_to_bead matrix the returned indices point into.

    number_of_steps : int
        The number of bead-move attempts to generate, i.e. the length of the
        returned array.

    latticeObject : Lattice
        The lattice object the beads belong to; used to walk the chains when
        frozen chains have to be excluded, and to sanity check the bead count if
        the safecheck flag is set.

    frozen_chains : sequence of int, optional
        chainIDs from which we do not select beads. Default is ``()``, meaning
        every bead is selectable.

    safecheck : bool, optional
        If True, check that the number of beads held by the lattice object
        matches ``num_beads``. This is a safety check to ensure the idx_to_bead
        matrix is not corrupted. Default is True.

    Returns
    -------
    numpy.ndarray
        An integer array of length ``number_of_steps`` holding the idx_to_bead
        row indices of the beads to attempt to move, drawn uniformly with
        replacement from the selectable beads. This essentially defines a random
        order in which we want to attempt to move beads.

    Raises
    ------
    MoveException
        If ``safecheck`` is True and the number of beads counted from
        ``latticeObject.chains`` does not match ``num_beads``.

    """
    layout = chain_layout(latticeObject)

    # the safety check used to walk every chain calling len(); the layout already
    # knows the total, and a mismatch means the table and the chains disagree
    if safecheck and layout.n_beads != num_beads:
        raise MoveException("The number of beads in the lattice object does not match the number of beads in the idx_to_bead matrix. This is a bug")

    # if frozen_chains is empty then we are randomly sampling from all possible beads
    if len(frozen_chains) == 0:
        return np.random.randint(0, num_beads, number_of_steps)

    # otherwise exclude the beads of the frozen chains. The selectable indices are
    # the rows of every non-frozen chain, in ascending chainID order - the SAME
    # order in which the bead table assigns row indices (a restart whose pickled
    # chains dict was not in ascending order used to freeze the wrong beads). The
    # draw itself is unchanged: np.random.choice with replace=True over the
    # selectable rows, so the random stream is identical to the old per-chain loop.
    frozen_chain = np.isin(np.asarray(layout.sorted_ids), np.asarray(list(frozen_chains)))
    selectable = np.flatnonzero(np.repeat(~frozen_chain, layout.length))
    return np.random.choice(selectable, number_of_steps, replace=True)
