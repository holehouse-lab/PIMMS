## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import math
import numbers
import random
import numpy as np


from .chain import Chain
from . import lattice_utils
from . import crankshaft_list_functions

from . import latticeExceptions
from .latticeExceptions import LatticeInitializationException, TypeGridException, RestartException, ChainInsertionFailure

from . CONFIG import NP_INT_TYPE, OUTPUT_CHAIN_TO_CHAINID


# Smallest box edge (in lattice sites) on which the energy model is well defined.
# The long-range and super-long-range shells reach offsets -3..3 along each axis,
# and those seven offsets only land on seven distinct sites when the axis is at
# least 7 long: in a box of 5 the +3 and -2 offsets are the same site (so one pair
# is scored as both LR and SLR), in a box of 6 the +3 and -3 offsets are (so one
# pair is scored twice). The compiled kernels also recognise a pair that reaches
# its partner only through the periodic wrap by a per-axis separation larger than
# 3, which is unambiguous only from 7 up. The keyfile parser has always refused
# smaller boxes; the Lattice enforces the same floor for programmatic use.
MIN_LATTICE_DIMENSION = 7


class Lattice:
    """The simulation box: the occupancy and bead-type grids, the chains that
    live on them, and the operations that keep the two representations
    consistent.

    A Lattice owns the ``grid`` (site -> chainID, 0 = empty), the ``type_grid``
    (site -> bead intcode), and a ``chains`` dict of :class:`~pimms.chain.Chain`
    objects keyed by 1-based chainID. It handles de-novo construction from a
    keyfile ``CHAIN`` specification, rebuilding from a restart file, box
    resizing for ``RESIZED_EQUILIBRATION``, backup/restore for revertible
    megamoves, and frame output (with optional PBC unwrapping / autocentring).
    """

    def __init__(self, dimensions, 
                 chain_list, 
                 Hamiltonian, 
                 lattice_to_angstroms,
                 chainsDict=None, 
                 lattice_grid=None, 
                 type_grid=None, 
                 restart_object=False, 
                 hardwall=False):

        """

        Lattice objects are the main type of object upon which simulations are run. Each simulation has one (and only one) 
        lattice


        Parameters
        -------------
        dimensions : list of size 2 or 3
            The 2D or 3D dimensions upon which the lattice is defined. Note that all dimensions
            may differ (non-cubic / non-square boxes are supported; the only
            restriction is cluster rotation under periodic boundaries). Every axis must
            be at least MIN_LATTICE_DIMENSION (7) sites long.

        chain_list : list of lists
            Each sublist is a tuple where element 0 is the number of chains and element 1 is the
            sequence of the chain.

        Hamiltonian : energy.Hamiltonian
            Hamiltonian object (as defined in energy.py) for the system. The Lattice
            itself only calls its three sequence-conversion methods
            (convert_sequence_to_integer_sequence,
            convert_sequence_to_LR_integer_sequence and
            get_indices_of_long_range_residues) while building chains, so any object
            that provides them (e.g. energy.EmptyHamiltonian) also works here; it is
            not used at all on the fully specified route.

        lattice_to_angstroms : float or int
            Value that defines the conversion factor of lattice units to angstroms.

        chainsDict : dict, optional
            Dictionary where keys are chainIDs and values are Chain objects. Only used
            if chainsDict, lattice_grid and type_grid are ALL provided, in which case
            the lattice state is taken verbatim from them (after validation). Default
            is None.

        lattice_grid : numpy.ndarray, optional
            2D or 3D integer numpy array with shape equal to dimensions - i.e. this is
            the grid upon which all the beads are defined, where positions are either 0
            (empty) or equal to a chainID. Only used alongside chainsDict and type_grid.
            Default is None.

        type_grid : numpy.ndarray, optional
            2D or 3D integer numpy array with shape equal to dimensions - i.e. this is
            the grid upon which all the beads are defined. Values are either 0 (solvent)
            or equal to a non-solvent bead type. Only used alongside chainsDict and
            lattice_grid. Default is None.

        restart_object : RestartObject, optional
            Object built from a restart file that contains all the information needed to
            reconstruct a lattice. Used only if the fully defined route above is not
            taken. Default is False (no restart).

        hardwall : bool, optional
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC. Default is False.

        Returns
        -------
        None
            No return value; the Lattice's grids, chains and crankshaft lookup
            tables are built in place.

        Raises
        -------------
        LatticeInitializationException
            If the dimensions are not 2 or 3 positive int32 values, if any axis is
            shorter than MIN_LATTICE_DIMENSION (7) sites, if hardwall is not a
            boolean, if lattice_to_angstroms is not a finite positive number, or if a
            fully specified (chainsDict/lattice_grid/type_grid) state is inconsistent.

        RestartException
            If a restart_object is supplied whose dimensionality does not match, or
            whose box is larger than, the requested lattice.

        ChainInsertionFailure
            If the chains to be placed (de novo, or the EXTRA_CHAINs of a restart)
            hold more beads in total than the box has sites, or if a chain could not
            be placed by random growth (see :class:`~pimms.chain.Chain`).

        """
        
        try:
            dimensions = tuple(dimensions)
        except TypeError:
            raise LatticeInitializationException(
                "Lattice dimensions must contain 2 or 3 positive integers")
        int32_max = np.iinfo(NP_INT_TYPE).max
        if (len(dimensions) not in (2, 3) or
                any(isinstance(value, (bool, np.bool_)) or
                    not isinstance(value, numbers.Integral) or
                    value <= 0 or value > int32_max for value in dimensions)):
            raise LatticeInitializationException(
                "Lattice dimensions must contain 2 or 3 positive int32 integers")
        # see MIN_LATTICE_DIMENSION: a shorter axis silently aliases the LR/SLR
        # shells, so the energies would be wrong rather than the run failing
        if min(dimensions) < MIN_LATTICE_DIMENSION:
            raise LatticeInitializationException(
                "Every lattice dimension must be at least %i sites (got %s): in a "
                "shorter box the long-range and super-long-range interaction shells "
                "wrap onto each other and the energies are wrong"
                % (MIN_LATTICE_DIMENSION, list(dimensions)))
        if not isinstance(hardwall, (bool, np.bool_)):
            raise LatticeInitializationException("hardwall must be True or False")
        if (isinstance(lattice_to_angstroms, (bool, np.bool_)) or
                not isinstance(lattice_to_angstroms, numbers.Real) or
                not math.isfinite(float(lattice_to_angstroms)) or
                lattice_to_angstroms <= 0):
            raise LatticeInitializationException(
                "lattice_to_angstroms must be a finite positive number")

        # define box dimensions (in lattice units)
        self.dimensions = [int(value) for value in dimensions]

        # Keep the boundary convention on the container as well as on every
        # Chain. This also normalizes programmatically supplied ``chainsDict``
        # objects, whose constructors may not have received the Lattice flag.
        self.hardwall     = bool(hardwall)

        # define conversion factor
        self.lattice_to_angstroms = float(lattice_to_angstroms)

        # (a long-dead cubic-box-only sanity check used to sit here inside a string
        # literal; non-cubic boxes have been fully supported since 1.0.5, and the
        # stray string was being picked up by autodoc as the attribute docstring of
        # lattice_to_angstroms - removed.)

        self.crankshaft_lists = []

        # if we have provided values for these three objects we are fully defining the lattice structure
        if chainsDict is not None and lattice_grid is not None and type_grid is not None:
            self.__fully_defined_initialization(dimensions, chain_list, Hamiltonian, chainsDict, lattice_grid, type_grid)     

        # if we have provided a restart object
        elif restart_object:
            self.__initialization_from_restart(Hamiltonian, restart_object, hardwall)

        else:            
            self.__de_novo_initialization(dimensions, chain_list, Hamiltonian, hardwall)

        for chain_object in self.chains.values():
            chain_object.hardwall = self.hardwall
            
        # either way we dynamically build the ID-TO-TYPE mapping dictionary at the end..
        self.chainIDtoType = {}
        self.chainTypeList = []
        
        for chainID in self.chains:

            # get each chain's type...
            CT = self.chains[chainID].chainType

            self.chainIDtoType[chainID] = CT
            if not CT in self.chainTypeList:
                self.chainTypeList.append(CT)


        # finally, initialize the crankshaft_list matrix for crankshaft moves, and build
        # the chain_to_firstbead_lookup dictionary, which allows us to look up specific
        # chains in the crankshaft_lists
        self.crankshaft_lists = crankshaft_list_functions.initialize_idx_to_bead(self)
        self.chain_to_firstbead_lookup = crankshaft_list_functions.initialize_chain_to_firstbead_lookup(self)
        # and the static per-chain layout (row offsets, lengths, homopolymer
        # flags) that the megamoves used to rebuild on every call
        self.chain_layout = crankshaft_list_functions.initialize_chain_layout(self)
        

                
    #-----------------------------------------------------------------
    #    ## CHANGEME        
    def __de_novo_initialization(self, dimensions, chain_list, Hamiltonian, hardwall):
        """
        Function that performs random initialization of a lattice based on the passed variables. 
        This is generally going to be the default behaviour for most simulations.
        

        Parameters
        -------------

        dimensions : list of size 2 or 3
            The 2D or 3D dimensions upon which the lattice is defined. Note that all dimensions
            may differ (non-cubic / non-square boxes are supported).

        chain_list : list of lists
            Each sublist is a tuple where element 0 is the number of chains and element 1 is the
            sequence of the chain.

        Hamiltonian : energy.Hamiltonian
            Hamiltonian object (as defined in energy.py); only its sequence-conversion
            methods are used here, to build each chain's integer codes.

        hardwall : bool
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.

        Returns
        -------------

        None
            No return value, but self.grid, self.type_grid and self.chains are built,
            with every chain placed at random (or centred, if the system is a single
            chain).

        Raises
        -------------
        ChainInsertionFailure
            If the chains hold more beads in total than the box has sites (checked
            before anything is placed), or if a chain could not be placed.

        """

        # A system with more beads than sites can never be placed, so say so before
        # placing anything. Otherwise we fill the box and then fail on whichever
        # chain comes last - and for a lone chain that failure came from the
        # centre-insertion path, which reported it as a PIMMS bug.
        self.__check_bead_capacity(sum(int(chain[0]) * len(chain[1]) for chain in chain_list))

        # intialize empty grids
        self.grid         = np.zeros(dimensions, dtype=NP_INT_TYPE)
        self.type_grid    = np.zeros(dimensions, dtype=NP_INT_TYPE)

        # initialize empty chains dictionary
        self.chains       = {}

        # initialize the chainID to 1
        chainID = 1

        # initialize the chainType to 0
        chainType = 0

        # set any/all flags used during initialization
        centerflag = False

        # if we're working with a SINGLE chain place in the center of the box,
        # else totally random
        if len(chain_list) == 1 and chain_list[0][0] == 1:
            centerflag = True
                                
        # for each chain tuple in the chain_list list
        for chain in chain_list:
                
            # extract the number and sequence of the chain
            n_chains   = chain[0]
            chain_seq  = chain[1]

            # for each chain in this chaingroup
            for i in range(0, n_chains):

                # create the integer_sequence associated with the chain's chemical makeup
                int_seq    = Hamiltonian.convert_sequence_to_integer_sequence(chain_seq)
                LR_int_seq = Hamiltonian.convert_sequence_to_LR_integer_sequence(chain_seq)
                LR_IDX     = Hamiltonian.get_indices_of_long_range_residues(chain_seq)
                                        
                # build a new ChainObjet
                ChainObject = Chain(self.grid, dimensions, chain_seq, int_seq, LR_int_seq, LR_IDX, chainID, chainType, center=centerflag, hardwall=hardwall)
                    
                # assign to the dictionary
                self.chains[chainID] = ChainObject

                chainID = chainID + 1

            chainType = chainType + 1
                    
        # initially the type grid is set to a numpy matrix of strings
        self.initialize_type_grid()


    #-----------------------------------------------------------------
    #            
    def __fully_defined_initialization(self, dimensions, chain_list, Hamiltonian, chainsDict, lattice_grid, type_grid):
        """
        Function that performs initialization of a lattice based on the passed variables.

        Parameters
        -------------

        dimensions : list of size 2 or 3
            The 2D or 3D dimensions upon which the lattice is defined. Note that all dimensions
            may differ (non-cubic / non-square boxes are supported).

        chain_list : list of lists
            Each sublist is a tuple where element 0 is the number of chains and element 1 is the
            sequence of the chain.

        Hamiltonian : energy.Hamiltonian
            Hamiltonian object (as defined in energy.py). Not used on this route (the
            chains arrive already built); accepted so the three initialisation routes
            share a signature.

        chainsDict : dict
            Dictionary of Chain objects that have been initialized elsewhere, keyed by
            chainID. Every entry is cross-checked against the two grids before anything
            is committed to the Lattice.

        lattice_grid : numpy.ndarray
            A 2D or 3D integer numpy array that represents the lattice grid, elements are
            either empty (0) or occupied by a bead where they report on the chainID of
            the bead.

        type_grid : numpy.ndarray
            A 2D or 3D integer numpy array of the same shape as lattice_grid, elements are
            either empty (0) or occupied by a bead where they report on the type of the
            bead.

        Returns
        -------------

        None
            No return value, but self.chains, self.grid and self.type_grid are set once
            the supplied state has been validated in full.

        Raises
        -------------
        LatticeInitializationException
            If either grid is not an integer numpy array, if chainsDict is not a
            dictionary, if the grid shapes do not match the lattice dimensions, or if
            the chains, occupancy grid and type grid do not describe the same
            configuration (mismatched keys, duplicate beads, ghost occupancy, or a
            bead type that disagrees with the chain sequence).

        """

        if (not isinstance(lattice_grid, np.ndarray) or
                not np.issubdtype(lattice_grid.dtype, np.integer)):
            raise LatticeInitializationException(
                "Provided lattice grid must be an integer numpy array")
        if (not isinstance(type_grid, np.ndarray) or
                not np.issubdtype(type_grid.dtype, np.integer)):
            raise LatticeInitializationException(
                "Provided type grid must be an integer numpy array")
        if not isinstance(chainsDict, dict):
            raise LatticeInitializationException("chainsDict must be a dictionary")

        # check the dimensions match up
        if tuple(self.dimensions) != tuple(lattice_utils.get_dimensions(lattice_grid)):
            raise LatticeInitializationException('Expected lattice dimensions (%s) did not match provided lattice-grid dimensions (%s)' %(str(dimensions), str(lattice_grid.shape)))
                
        if tuple(self.dimensions) != tuple(lattice_utils.get_dimensions(type_grid)):
            raise LatticeInitializationException('Expected type_grid lattice dimensions (%s) did not match provided lattice-grid dimensions (%s)' %(str(dimensions), str(type_grid.shape)))

        expected_ids = set(chainsDict)
        occupied = set()
        bead_count = 0
        for key, chain in chainsDict.items():
            if key != chain.chainID:
                raise LatticeInitializationException(
                    f"chainsDict key {key} does not match Chain.chainID {chain.chainID}")
            positions = chain.get_ordered_positions()
            if len(positions) != len(chain.int_sequence):
                raise LatticeInitializationException(
                    f"Chain {key} position and integer-sequence lengths differ")
            bead_count += len(positions)
            for bead_index, position in enumerate(positions):
                coordinate = tuple(position)
                if coordinate in occupied:
                    raise LatticeInitializationException(
                        f"More than one chain bead occupies {coordinate}")
                occupied.add(coordinate)
                try:
                    grid_value = int(lattice_grid[coordinate])
                    type_value = int(type_grid[coordinate])
                except (IndexError, TypeError):
                    raise LatticeInitializationException(
                        f"Chain {key} contains invalid position {position!r}")
                if grid_value != key:
                    raise LatticeInitializationException(
                        f"Lattice grid at {coordinate} does not contain chain {key}")
                if type_value != int(chain.int_sequence[bead_index]):
                    raise LatticeInitializationException(
                        f"Type grid at {coordinate} does not match bead "
                        f"{bead_index} of chain {key}")

        # A fully supplied state is built once, so scan the grids here to catch
        # ghost occupancy/type sites that do not belong to any Chain object.
        if int(np.count_nonzero(lattice_grid)) != bead_count:
            raise LatticeInitializationException(
                "Lattice grid occupancy does not match the supplied chains")
        if np.any(type_grid[lattice_grid == 0] != 0):
            raise LatticeInitializationException(
                "Type grid contains bead types at empty lattice sites")
        if expected_ids and not set(np.unique(lattice_grid)).issubset(expected_ids | {0}):
            raise LatticeInitializationException(
                "Lattice grid contains a chainID absent from chainsDict")

        # Commit only after the supplied state has been validated in full.
        self.chains = chainsDict
        self.grid = lattice_grid
        self.type_grid = type_grid


    #-----------------------------------------------------------------
    #            
    def __initialization_from_restart(self, Hamiltonian, restart, hardwall):
        """
        Function that sets up a lattice based on a restart object.

        Parameters
        -------------
        Hamiltonian : energy.Hamiltonian
            Hamiltonian object (as defined in energy.py); only its sequence-conversion
            methods are used here, to build each chain's integer codes.

        restart : RestartObject
            RestartObject that contains all the information needed to restart a
            simulation (dimensions, per-chain positions/sequence/chainType, and any
            extra chains to be placed at random).

        hardwall : bool
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.

        Returns
        -------------

        None
            No return value, but self.grid, self.type_grid and self.chains are built
            from the restart object.

        Raises
        -------------
        RestartException
            If the restart object's dimensionality does not match the lattice, if the
            restart box is larger than the lattice, or if a chainID appears twice.

        ChainInsertionFailure
            If the restart chains plus the EXTRA_CHAINs hold more beads in total than
            the box has sites (checked before anything is placed), or if an extra
            chain could not be placed.

        """

        # check dimensions of restart match passed dimensions
        if len(restart.dimensions) != len(self.dimensions):
            raise RestartException('Number of dimensions in restart file do not match number of dimensions in keyfile')

        for A, B in zip(restart.dimensions, self.dimensions):
            if A > B:
                raise RestartException('Dimensions associated with new lattice are smaller than lattice from the restart object. This is not allowed.')

        # the restart chains fit by construction, but EXTRA_CHAINs can push the
        # system past the number of sites - catch that before placing anything
        n_beads = (sum(len(restart.chains[c][1]) for c in restart.chains) +
                   sum(len(restart.extra_chains[c][1]) for c in restart.extra_chains))
        self.__check_bead_capacity(n_beads)

        # intialize empty grids
        self.grid         = np.zeros(self.dimensions, dtype=NP_INT_TYPE)
        self.type_grid    = np.zeros(self.dimensions, dtype=NP_INT_TYPE)

        # initialize empty chains dictionary
        self.chains       = {}

        # for each chain, extract all info and insert into the lattice grid.
        # Iterate in ascending chainID order regardless of the pickle's key order:
        # the trajectory/PDB writers and chain_to_chainid.txt follow dict order,
        # whereas the per-chain analysis columns (RG, ASPH, END_TO_END, RES_TO_RES)
        # are written in sorted order, so an unsorted restart would silently put
        # the two in different orders.
        for chainID in sorted(restart.chains):
            
            if chainID in self.chains:
                raise RestartException(f'Error when adding chain extracted from Restart file to lattice. ChainID={chainID} was already found in the chains list. This is a major bug')
                
            chain_info = restart.chains[chainID]

            # extract info from restart object
            chain_pos  = chain_info[0]
            chain_seq  = chain_info[1]
            chainType = chain_info[2]


            # build internal representation 
            int_seq    = Hamiltonian.convert_sequence_to_integer_sequence(chain_seq)
            LR_int_seq = Hamiltonian.convert_sequence_to_LR_integer_sequence(chain_seq)
            LR_IDX     = Hamiltonian.get_indices_of_long_range_residues(chain_seq)

            # Restart coordinates retain the boundary convention of the lattice.
            # This matters for every coordinate-derived Chain observable.
            ChainObject = Chain(self.grid, self.dimensions, chain_seq, int_seq, LR_int_seq, LR_IDX,
                                chainID, chainType, center=False, chain_positions=chain_pos,
                                hardwall=hardwall)

            # insert into lattice grid
            lattice_utils.place_chain_by_position(chain_pos, self.grid, chainID, safe=False)

            # update the chains dictionary
            self.chains[chainID] = ChainObject

        ### Add in extra chains
        ### 
        # for each extra chain (ascending order, as above):
        for chainID in sorted(restart.extra_chains):
            if chainID in self.chains:
                raise RestartException(f'Error when adding chain defined as a EXTRA_CHAIN to the lattice. ChainID={chainID} was already found in the chains list. This is a major bug')
                
            chain_info = restart.extra_chains[chainID]

            # extract the chain sequence and chain type
            chain_seq  = chain_info[1]
            chainType = chain_info[2]

            # build internal representation 
            int_seq    = Hamiltonian.convert_sequence_to_integer_sequence(chain_seq)
            LR_int_seq = Hamiltonian.convert_sequence_to_LR_integer_sequence(chain_seq)
            LR_IDX     = Hamiltonian.get_indices_of_long_range_residues(chain_seq)

            # add the new object without chain_pos variable
            ChainObject = Chain(self.grid, self.dimensions, chain_seq, int_seq, LR_int_seq, LR_IDX, chainID, chainType, center=False, hardwall=hardwall)
            # update the chains dictionary
            self.chains[chainID] = ChainObject

        # Finally initially the type grid is set to a numpy matrix of strings
        self.initialize_type_grid()


    #-----------------------------------------------------------------
    #
    def __check_bead_capacity(self, n_beads):
        """
        Refuse a system that holds more beads than the box has sites.

        Such a system can never be placed, however the chains are arranged, so we
        check this before placing anything and report it as the overcrowding it
        is. Having at most one bead per site is necessary but not sufficient: a
        system that passes can still fail later if random chain growth cannot find
        room.

        Parameters
        -------------
        n_beads : int
            Total number of beads (summed over every chain) to be put on the
            lattice.

        Returns
        -------------
        None
            Returns only if n_beads is no larger than the number of lattice sites.

        Raises
        -------------
        ChainInsertionFailure
            If n_beads exceeds the number of sites in self.dimensions.

        """
        n_sites = int(np.prod(self.dimensions))
        if n_beads > n_sites:
            raise ChainInsertionFailure(
                '\nUnable to place the chains: they hold %i beads in total but the %s box has only %i sites.\n'
                'The lattice is overcrowded - use a larger box (DIMENSIONS) or fewer/shorter chains...\n'
                % (n_beads, ' x '.join(str(d) for d in self.dimensions), n_sites))


    #-----------------------------------------------------------------
    #
    def get_number_of_chains(self):
        """
        Function that returns the number of chains in the lattice

        Returns
        ------------
        int
            Number of chains in the lattice
        
        """
        return len(self.chains)


    #-----------------------------------------------------------------
    #
    def check_grid_consistency(self):
        """
        Cross-check the occupancy grid and the type grid against the chains.

        The full-energy recompute reads bead *types* from ``type_grid`` - the same
        array the compiled kernels read - so a type-grid corruption would make the
        tracked and recomputed energies wrong identically and the ENERGY_CHECK
        comparison could not see it. This walks every chain and verifies that each
        bead's site holds its chainID on ``grid`` and its integer residue code on
        ``type_grid``, and that no other site is occupied.

        Returns
        -------
        list of str
            Human-readable descriptions of every inconsistency found (empty when
            the grids are consistent with the chains).
        """
        problems = []
        expected_occupied = 0
        for chainID, chain in self.chains.items():
            positions = chain.get_ordered_positions()
            int_sequence = chain.int_sequence
            if len(int_sequence) != len(positions):
                problems.append("chain %s: %i beads but %i residue codes"
                                % (chainID, len(positions), len(int_sequence)))
            expected_occupied += len(positions)
            for bead_idx, position in enumerate(positions):
                site = tuple(position)
                occupant = int(self.grid[site])
                if occupant != chainID:
                    problems.append("chain %s bead %i at %s: grid holds chain %i"
                                    % (chainID, bead_idx, list(position), occupant))
                if bead_idx < len(int_sequence):
                    stored_type = int(self.type_grid[site])
                    if stored_type != int(int_sequence[bead_idx]):
                        problems.append("chain %s bead %i at %s: type_grid holds %i, sequence says %i"
                                        % (chainID, bead_idx, list(position), stored_type,
                                           int(int_sequence[bead_idx])))

        occupied = int(np.count_nonzero(self.grid))
        if occupied != expected_occupied:
            problems.append("grid has %i occupied sites but the chains hold %i beads"
                            % (occupied, expected_occupied))
        typed = int(np.count_nonzero(self.type_grid))
        if typed != expected_occupied:
            problems.append("type_grid has %i typed sites but the chains hold %i beads"
                            % (typed, expected_occupied))
        return problems


    #-----------------------------------------------------------------
    #
    def get_gridvalue(self, position):
        """
        Function that returns the value on the main lattice grid at a given position (i.e.
        will return 0 or the chainID of the chain that occupies that position).

        Parameters
        ------------
        position : list
            Position in the lattice grid (len=2 or len=3).

        Returns
        ------------
        int
            Value on the lattice grid at the given position (0 if empty, else the
            chainID of the chain that occupies it)

        Raises
        ------------
        LatticeUtilsException
            If position is neither 2D nor 3D.

        """
        return lattice_utils.get_gridvalue(position, self.grid)


    #-----------------------------------------------------------------
    #
    def set_gridvalue(self, position, value):
        """
        Function that sets the value on the main lattice grid at a given position (i.e.
        will set 0 or the chainID of the chain that occupies that position). Note this
        does not have any sanity checking. 

        Parameters
        ------------
        position : list
            Position in the lattice grid (len=2 or len=3).

        value : int
            Value to set at the given position (0 for solvent, else a chainID)

        Returns
        ------------
        numpy.ndarray
            The lattice grid (self.grid), which has been modified in place.

        Raises
        ------------
        LatticeUtilsException
            If position is neither 2D nor 3D.

        """
        return lattice_utils.set_gridvalue(position, value, self.grid)


    #-----------------------------------------------------------------
    #
    def save_as_pdb(self, fname):
        """
        Function that saves the current Lattice as a PDB file

        Parameters
        --------------
        fname : str
            Name of PDB file

        Returns
        ------------
        None

        """
        lattice_utils.open_pdb_file(self.dimensions, self.lattice_to_angstroms, fname)
        lattice_utils.write_lattice_to_pdb(self, self.lattice_to_angstroms, fname, write_connect=True)
        lattice_utils.finish_pdb_file(fname)


    #-----------------------------------------------------------------
    #
    def any_chains_straddle_boundary(self):
        """
        Function that scans each chain associated with the lattice and
        asks if ANY chain straddles the boundary. If any chain does
        returns True, if no chain does returns false.

        Returns
        ------------
        bool
            True if any chain straddles the boundary, False otherwise


        """
        for chainID in self.chains:
            if self.chains[chainID].does_chain_stradle_pbc_boundary():
                return True

        return False

    #-----------------------------------------------------------------
    #
    def get_random_chain(self, frozen_chains=None):
        """
        Randomly select and return a chain object from the lattice.

        If ``frozen_chains`` is provided, those chain IDs are excluded and a
        chain is selected uniformly at random from the remaining selectable
        chains.

        Parameters
        ------------
        frozen_chains : list, optional
            List of chainIDs to exclude from the selection (i.e. frozen
            chains). If None or empty, all chains are eligible. Default
            is None.

        Returns
        ------------
        Chain
            A randomly selected Chain object.

        Raises
        ------
        LatticeInitializationException
            If there are no chains available to select from (either no
            chains exist at all, or every chain is frozen).

        """

        if frozen_chains is None or len(frozen_chains) == 0:
            # Nothing is frozen, so the candidates are every chain, in
            # dictionary order. Building that list afresh on every step was
            # O(chains) per step (0.36 ms at 10^4 chains - more than a chain
            # translation costs), so it is kept and rebuilt only when the chains
            # change. New chains are always added at the end of the dictionary,
            # so a change shows up as a different chain count or a different
            # last chainID. random.choice is handed the same sequence either
            # way, so the draw is unchanged.
            candidate_chain_ids = getattr(self, '_all_chain_ids', None)
            if (candidate_chain_ids is None
                    or len(candidate_chain_ids) != len(self.chains)
                    or (candidate_chain_ids and candidate_chain_ids[-1] != next(reversed(self.chains)))):
                candidate_chain_ids = list(self.chains)
                self._all_chain_ids = candidate_chain_ids
        else:
            frozen_set = set(frozen_chains)
            candidate_chain_ids = [chain_id for chain_id in self.chains
                                   if chain_id not in frozen_set]

        if not candidate_chain_ids:
            if not self.chains:
                raise LatticeInitializationException(
                    "No chains are available for random selection")
            raise LatticeInitializationException(
                "No selectable chains are available (all chains are frozen)")

        return self.chains[random.choice(candidate_chain_ids)]
        


    #-----------------------------------------------------------------
    #
    def initialize_type_grid(self):
        """
        Initualizes the type grid using the chain's int_sequence
        types (i.e. int_sequence contains the chain's sequence where
        bead types (letters) are represented by integers)

        Returns
        ------------
        None
            No return value, but self.type_grid is written for every bead of
            every chain in self.chains.

        """

        for chainID in self.chains:
            chain = self.chains[chainID]
            int_sequence  = chain.int_sequence
            positions = chain.get_ordered_positions()
            
            for i in range(0, len(positions)):
                lattice_utils.set_gridvalue(positions[i], int_sequence[i], self.type_grid)
            
    
    #-----------------------------------------------------------------
    #
    def update_type_grid(self, chainID, old_positions, new_positions, indices, safe=True):
        """
        Update the type grid for the chain associated with the passed
        chainID based on the position vectors passed.

        Might implement chain positions as a selective index vector so you only remove/delete
        the specific parts of a chain being modified, but for that to be worth it it'd have
        to become clear this is a bottleneck performance wise.

        Parameters
        ------------
        chainID : int
            ID of the chain to be updated

        old_positions : list
            List of old positions (each a 2- or 3-element coordinate) to be removed
            from the type grid

        new_positions : list
            List of new positions (each a 2- or 3-element coordinate) to be added to
            the type grid

        indices : list
            List of indices into the chain sequence which correspond to the positions.
            The same indices apply to both old_positions and new_positions.

        safe : bool, optional
            If True (the default), the delete and insert steps sanity check the type
            grid as they go and raise a TypeGridException on any inconsistency. If
            False the sites are simply overwritten, which is faster but assumes the
            caller already knows the update is valid.

        Returns
        ------------
        None
            No return value, but self.type_grid is updated in place.

        Raises
        ------------
        TypeGridException
            If safe is True and either the positions/indices lengths mismatch, an old
            site does not hold the expected bead type, or a new site is not empty.

        """
        
        self.delete_chain_from_type_grid(chainID, old_positions, indices, safe)
        self.insert_chain_into_type_grid(chainID, new_positions, indices, safe)

            
    #-----------------------------------------------------------------
    #                
    def delete_chain_from_type_grid(self, chainID, positions, indices, safe=True):
        """
        Function which deletes the positions/indices associated with chainID from the type grid.         

        Indices should be a vector of indices in the chain which correspond to the positions - i.e. if indies were
        [4,5,6,7] then positions would be a list or array of length 4 where the positions correspond to the positions
        of residues, 4, 5, 6, and 7, respectively.

        Parameters
        ------------
        chainID : int
            ID of the chain to be deleted

        positions : list
            List of positions (each a 2- or 3-element coordinate) to be deleted from
            the type grid

        indices : list
            List of indices into the chain sequence which correspond to the positions

        safe : bool, optional
            If True (the default), check that the positions and indices are the same
            length and that each site currently holds either the bead type the chain's
            sequence expects or 0 (solvent), raising a TypeGridException if not. If
            False the sites are wiped without any checking.

        Returns
        ------------
        None
            No return value, but the relevant sites of self.type_grid are set to 0.

        Raises
        ------
        TypeGridException
            If ``safe`` is True and either the positions/indices lengths
            mismatch or an existing site does not match the chain's
            expected bead type.

        """

        # first check the chain exists and extract the sequence
        chain = self.chains[chainID]

        # get the int sequence 
        sequence  = chain.int_sequence

        # if safe, check the new positions are not already occupied

        if safe:                        
            # check the indices and positions match
            if not len(positions) == len(indices):
                raise TypeGridException(f"Trying to delete positions {indices} of ChainID [{chainID}] from the typegrid but indices and positions do not match")

            # next delete the old positions - again ensuring we're only wiping an old chain
            for i in range(0, len(positions)):     
                current = lattice_utils.get_gridvalue(positions[i], self.type_grid)

                if not (current == sequence[indices[i]] or current == 0):
                    raise TypeGridException('Trying to update the type grid for chain %i, residue %i but the current operation would delete at position ['%(chainID,i) + str(positions[i]) + "], which does not match the expected type based on the chain's sequence - chain's sequence wants [%s] while current positions is [%s]"%(sequence[indices[i]], current))
                
                # delete the type by replacing with a 0 (solvent)
                lattice_utils.set_gridvalue(positions[i], 0,  self.type_grid)
        else:

            # if not safe, just delete the old positions without worrying about what's there
            for i in range(0, len(positions)):  
                lattice_utils.set_gridvalue(positions[i],  0, self.type_grid)


    #-----------------------------------------------------------------
    #
    def insert_chain_into_type_grid(self, chainID, positions, indices, safe=True):
        """
        Function which inserts the positions/indices associated with chainID into the typeGrid. 

        Indices should be a vector of indices in the chain which correspond to the positions - 
        i.e. if indies were [4,5,6,7] then positions would be a list or array of length 4 
        where the positions correspond to the positions of residues, 4, 5, 6, and 7, 
        respectively.

        Parameters
        ------------
        chainID : int
            ID of the chain to be inserted

        positions : list
            List of positions (each a 2- or 3-element coordinate) to be inserted into
            the type grid

        indices : list
            List of indices into the chain sequence which correspond to the positions

        safe : bool, optional
            If True (the default), check that the positions and indices are the same
            length and that every target site is currently empty (0), raising a
            TypeGridException if not. If False the bead types are written without any
            checking.

        Returns
        ------------
        None
            No return value, but the relevant sites of self.type_grid are set to the
            chain's integer bead types.

        Raises
        ------------
        TypeGridException
            If ``safe`` is True and either the positions/indices lengths mismatch or a
            target site is already occupied by a bead type.

        """

        chain = self.chains[chainID]
        sequence  = chain.int_sequence

        if safe:

            if not len(positions) == len(indices):
                raise TypeGridException(f"Trying to insert ChainID [{chainID}] into the type grid but set of chain indices does not match set of positions to futz with")

            for i in range(0, len(positions)):                        
                current = lattice_utils.get_gridvalue(positions[i], self.type_grid)

                if not current == 0:
                    raise TypeGridException('Trying to update the type grid  but the current operation would over-write the site at [' + str(positions[i]) + "] with type [%s] when it's currently set to [%s]"%(sequence[indices[i]], current))
                
                lattice_utils.set_gridvalue(positions[i], sequence[indices[i]], self.type_grid)
        else:
            for i in range(0, len(positions)): 
                lattice_utils.set_gridvalue(positions[i], sequence[indices[i]], self.type_grid)

    
    def lattice_backupcopy(self):
        """
        Function which returns a 3-place tuple with
        1) The lattice main grid (in its current state)
        2) The lattice type grid (in its current state)
        3) A dictionary of the chain positions 
        
        In all cases these are deep copies of the original. This function is only really relevant if we want to make
        a system-wide backup to restore to at a later date. This was originally written for TMMMC moves (restore state
        should the move be rejected) but in principle should.

        Note - this will basically double the memory footprint of the simulation temporarily so, be careful!

        Returns
        ------------
        tuple
            3-place tuple ``(grid, type_grid, chain_positions)`` where the first two
            are copies of the 2D/3D integer lattice and type grids and the third is a
            dictionary keyed by chainID whose values are the chain positions as lists
            of lists


        """
        grid_copy = np.copy(self.grid)
        type_grid_copy = np.copy(self.type_grid)
        chain_pos_copy = {} 

        for i in self.chains:
            chain_pos_copy[i] = np.copy(self.chains[i].positions).tolist()

        return (grid_copy, type_grid_copy, chain_pos_copy)


    def lattice_restorefrombackup(self, grid, type_grid, chain_dict):
        """
        Restore the lattice and its chains to a previously backed-up state.

        The complete replacement state is validated before any live object is
        changed.  A malformed backup therefore raises without leaving the grid,
        type grid, and chain objects describing different configurations.

        Parameters
        ------------
        grid : numpy.ndarray
            The lattice main grid to be restored: a 2D or 3D integer array with the
            same shape as the lattice dimensions.

        type_grid : numpy.ndarray
            The lattice type grid to be restored: a 2D or 3D integer array with the
            same shape as the lattice dimensions.

        chain_dict : dict
            A dictionary of the chain positions to be restored (keys = chainID,
            values = list of positions). The keys must match the lattice's chains
            exactly. These three arguments are what lattice_backupcopy() returns.

        Returns
        ------------
        None
            No return value, but self.grid, self.type_grid and every Chain's positions
            are replaced once the whole backup has passed validation.

        Raises
        ------------
        LatticeInitializationException
            If either grid is not an integer array of the right shape, if the chainIDs
            do not match the lattice, if any chain's positions are the wrong length,
            wrong dimensionality, non-integer, outside the box, duplicated, or
            inconsistent with the supplied grids, or if either grid holds an occupied
            or typed site that belongs to no chain bead.

        """

        expected_shape = tuple(self.dimensions)
        if (not isinstance(grid, np.ndarray) or
                not np.issubdtype(grid.dtype, np.integer) or
                grid.shape != expected_shape):
            raise LatticeInitializationException(
                "Backup lattice grid must be an integer array with shape "
                f"{expected_shape}")
        if (not isinstance(type_grid, np.ndarray) or
                not np.issubdtype(type_grid.dtype, np.integer) or
                type_grid.shape != expected_shape):
            raise LatticeInitializationException(
                "Backup type grid must be an integer array with shape "
                f"{expected_shape}")

        expected_ids = set(self.chains)
        if set(chain_dict) != expected_ids:
            raise LatticeInitializationException(latticeExceptions.message_preprocess(
                "Trying to re-set the lattice using lattice_restorefrombackup but "
                "the supplied chain IDs do not exactly match the lattice chains."))

        restored_positions = {}
        occupied = set()
        for chain_id, chain in self.chains.items():
            positions = chain_dict[chain_id]
            if len(positions) != len(chain):
                raise LatticeInitializationException(
                    f"Backup chain {chain_id} has {len(positions)} positions; "
                    f"expected {len(chain)}")

            copied_positions = []
            for bead_index, position in enumerate(positions):
                if len(position) != len(expected_shape):
                    raise LatticeInitializationException(
                        f"Backup position {position!r} for chain {chain_id} has the "
                        "wrong dimensionality")
                if any(isinstance(value, (bool, np.bool_)) or
                       not isinstance(value, numbers.Integral) for value in position):
                    raise LatticeInitializationException(
                        f"Backup position {position!r} for chain {chain_id} must "
                        "contain integer coordinates")
                copied = [int(value) for value in position]
                if any(value < 0 or value >= expected_shape[axis]
                       for axis, value in enumerate(copied)):
                    raise LatticeInitializationException(
                        f"Backup position {position!r} for chain {chain_id} is "
                        "outside the lattice")
                coordinate = tuple(copied)
                if coordinate in occupied:
                    raise LatticeInitializationException(
                        f"Backup contains more than one bead at {coordinate}")
                occupied.add(coordinate)
                if int(grid[coordinate]) != chain_id:
                    raise LatticeInitializationException(
                        f"Backup grid at {coordinate} does not contain chain "
                        f"{chain_id}")
                if int(type_grid[coordinate]) != int(chain.int_sequence[bead_index]):
                    raise LatticeInitializationException(
                        f"Backup type grid at {coordinate} does not match bead "
                        f"{bead_index} of chain {chain_id}")
                copied_positions.append(copied)
            restored_positions[chain_id] = copied_positions

        # Every chain bead has now been matched against both grids, but that says
        # nothing about the rest of the grids. Mirror the fully defined
        # initialisation and require that they hold nothing else: an extra
        # occupied site would be a ghost bead no Chain owns (it blocks moves), and
        # an extra typed site would interact with its neighbours in every energy
        # evaluation while being invisible to the chains.
        if int(np.count_nonzero(grid)) != len(occupied):
            raise LatticeInitializationException(
                "Backup lattice grid occupancy does not match the chains")
        if np.any(type_grid[grid == 0] != 0):
            raise LatticeInitializationException(
                "Backup type grid contains bead types at empty lattice sites")

        # Commit only after the whole replacement state has passed validation.
        self.grid = grid
        self.type_grid = type_grid
        for chain_id, positions in restored_positions.items():
            self.chains[chain_id].set_ordered_positions(positions)


    def write_chain_to_chainid_file(self):
        """
        Function which writes the chainID-to-sequence mapping file. Writes one
        tab-delimited line per chain in the lattice (chainID, sequence length,
        sequence) to the file named by CONFIG.OUTPUT_CHAIN_TO_CHAINID. This is
        how a chainID in the trajectory or per-chain analysis files is mapped
        back to the sequence it represents.

        Returns
        ------------
        None
            No return value; the file is written (and overwritten if present).

        """

        with open(OUTPUT_CHAIN_TO_CHAINID, 'w', encoding='utf-8') as fh:
            for chainID in self.chains:
                seq = self.chains[chainID].sequence
                fh.write(f'{chainID}\t{len(seq)}\t{seq}\n')

        
