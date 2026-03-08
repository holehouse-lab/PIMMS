## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................


import random
import numpy as np


from .chain import Chain
from . import lattice_utils
from . import crankshaft_list_functions

from . import latticeExceptions
from .latticeExceptions import LatticeInitializationException, TypeGridException, ParameterFileException, RestartException
from . import CONFIG

from . CONFIG import NP_INT_TYPE, OUTPUT_CHAIN_TO_CHAINID


class Lattice:

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
            must be equal (although we plan to update this at somepoint soon).

        chain_list : list of lists
            Each sublist is a tuple where element 0 is the number of chains and element 1 is the 
            sequence of the chain.
        
        Hamiltonian : energy.Hamiltonian (or energy.EmptyHamiltonian)
            Hamltionian object (as defined in energy.py) that provides a way to compute the energy 
            of the system

        lattice_to_angstroms : float or int
            Value that defines the conversion factor of lattice units to angstroms.

        chainsDict : dict {False}
            Dictionary where keys are chain indices and values are chain.Chain objects.

        lattice_grid : np.ndarray {False}
            2D or 3D numpy array defined INITIALLY as np.zeros(dimensions, dtype=int) - i.e. this is 
            the grid upon which all the beads are defined, where positions are either 0 (empty) or 
            equal to a chainID.

        type_grid : np.ndarray {False}
            2D or 3D numpy array defined INITIALLY as np.zeros(dimensions, dtype=int) - i.e. this is 
            the grid upon which all the beads are defined. Values are either 0 (solvent) or equal
            to a non-solvent bead type.
        
        restart_object : restart.RestartObject {False}
            Object built from a restart file that contains all the information needed to reconstruct
            a lattice. 

        hardwall : bool {False}
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.

        hardwall : bool {False}
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.
        
        """
        
        # define box dimensions (in lattice units)
        self.dimensions   = dimensions

        # define conversion factor
        self.lattice_to_angstroms = lattice_to_angstroms

        '''
        # ensure lattice dimensions are consistent. This is a temporary check...
        if len(dimensions) == 2:
            if dimensions[0] != dimensions[1]:                
                raise LatticeInitializationException(latticeExceptions.message_preprocess('In the current version of PIMMS the X/Y dimensions must be equal. In fact this will be updated soon, but, for now avoid passing in non-matching X/Y dimensions'))
        else:  ## CHANGEME
            if dimensions[0] != dimensions[1] or dimensions[1] != dimensions[2] or dimensions[0] != dimensions[2]:
                raise LatticeInitializationException(latticeExceptions.message_preprocess('In the current version of PIMMS the X/Y/Z dimensions must be equal. In fact this will be updated soon, but, for now avoid passing in non-matching X/Y/Z dimensions'))
        ''' 

        self.crankshaft_lists = []

        # if we have provided values for these three objects we are fully defining the lattice structure
        if chainsDict is not None and lattice_grid is not None and type_grid is not None:
            self.__fully_defined_initialization(dimensions, chain_list, Hamiltonian, chainsDict, lattice_grid, type_grid)     

        # if we have provided a restart object
        elif restart_object:
            self.__initialization_from_restart(Hamiltonian, restart_object, hardwall)

        else:            
            self.__de_novo_initialization(dimensions, chain_list, Hamiltonian, hardwall)
            
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
            must be equal.

        chain_list : list of lists
            Each sublist is a tuple where element 0 is the number of chains and element 1 is the
            sequence of the chain.

        Hamiltonian : energy.Hamiltonian (or energy.EmptyHamiltonian)
            Hamltionian object (as defined in energy.py) that provides a way to compute the energy
            of the system

        hardwall : bool
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.
        
        Returns
        -------------

        None
            

        """

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
            must be equal.

        chain_list : list of lists
            Each sublist is a tuple where element 0 is the number of chains and element 1 is the
            sequence of the chain.

        Hamiltonian : energy.Hamiltonian (or energy.EmptyHamiltonian)
            Hamltionian object (as defined in energy.py) that provides a way to compute the energy
            of the system

        chainsDict : dict
            Dictionary of Chain objects that have been initialized elsewhere. This is the primary
            way to initialize a lattice from a restart file.

        lattice_grid : numpy array
            A numpy array that represents the lattice grid, elements are either empty (0) or occupied
            by a bead where they report on the chainID of the bead.

        type_grid : numpy array
            A numpy array that represents the lattice grid, elements are either empty (0) or occupied
            by a bead where they report on the type of the bead.

        Returns
        -------------

        None

        """

        # check the dimenions match up
        if tuple(self.dimensions) != tuple(lattice_utils.get_dimensions(lattice_grid)):
            raise LatticeInitializationException('Expected lattice dimensions (%s) did not match provided lattice-grid dimensions (%s)' %(str(dimensions), str(lattice_grid.shape)))
                
        if tuple(self.dimensions) != tuple(lattice_utils.get_dimensions(type_grid)):
            raise LatticeInitializationException('Expected type_grid lattice dimensions (%s) did not match provided lattice-grid dimensions (%s)' %(str(dimensions), str(type_grid.shape)))

        # set everything
        self.chains    = chainsDict
        self.grid      = lattice_grid
        self.type_grid = type_grid
            
        # sanity check
        if CONFIG.DEBUG:
            for chain in chainsDict.values():
                chainID   = chain.chainID
                positions = chain.get_ordered_positions()

                for position in positions:
                    if not lattice_utils.get_gridvalue(position, self.grid) == chainID:
                        raise LatticeInitializationException('Uh oh! - Check debug portion')                    


    #-----------------------------------------------------------------
    #            
    def __initialization_from_restart(self, Hamiltonian, restart, hardwall):
        """
        Function that sets up a lattice based on a restart object.

        Parameters
        -------------
        Hamiltonian : energy.Hamiltonian (or energy.EmptyHamiltonian)
            Hamltionian object (as defined in energy.py) that provides a way to compute the energy
            of the system

        restart : restart.Restart
            Restart object that contains all the information needed to restart a simulation

        hardwall : bool
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.

        Returns
        -------------
        
        None
        
        """

        # check dimensions of restart match passed dimensions 
        if len(restart.dimensions) != len(self.dimensions):
            raise RestartException('Number of dimensions in restart file do not match number of dimensions in keyfile')
        
        for A, B in zip(restart.dimensions, self.dimensions):
            if A > B:
                raise RestartException('Dimensions associated with new lattice are smaller than lattice from the restart object. This is not allowed.')
                
        # intialize empty grids
        self.grid         = np.zeros(self.dimensions, dtype=NP_INT_TYPE)
        self.type_grid    = np.zeros(self.dimensions, dtype=NP_INT_TYPE)

        # initialize empty chains dictionary
        self.chains       = {}

        # for each chain, extract all info and insert into the lattice grid. 
        for chainID in restart.chains:   
            
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

            # build a new chain object. Note we don't need to worry about hardwall here because it'll be honoring whatever the
            # appropriate schema is
            ChainObject = Chain(self.grid, self.dimensions, chain_seq, int_seq, LR_int_seq, LR_IDX, chainID, chainType, center=False, chain_positions=chain_pos)

            # insert into lattice grid
            lattice_utils.place_chain_by_position(chain_pos, self.grid, chainID, safe=False)

            # update the chains dictionary
            self.chains[chainID] = ChainObject

        ### Add in extra chains
        ### 
        # for each extra chain:
        for chainID in restart.extra_chains:
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
            Value on the lattice grid at the given position
        
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
            Value to set at the given position

        Returns
        ------------
        None
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
        Randomly selects and returns a chain object. If no override is
        provided selects any of the possible chains. If override is provided
        it is assumed each position in override is a valid chainID and one 
        is randomly used to extract a chain

        Parameters
        ------------
        override : list
            List of chainIDs to choose from. If empty, any chain is chosen.

        exclude : list
            List of chainIDs to exclude from the selection. If empty, 
            all chains are considered.

        Returns
        ------------
        Chain
            Chain object

        """

        # we have no override list
        if frozen_chains is None:
            frozen_chains = []

        if len(frozen_chains) == 0:

            if not self.chains:
                raise LatticeInitializationException("No chains are available for random selection")

            # randomly choose from the list of all chainIDs
            return self.chains[random.choice(list(self.chains.keys()))]
            
        else:
            frozen_set = set(frozen_chains)
            candidate_chain_ids = [chain_id for chain_id in self.chains if chain_id not in frozen_set]

            if not candidate_chain_ids:
                raise LatticeInitializationException("No selectable chains are available (all chains are frozen)")

            # randomly choose from the list of all chainIDs that are not in the frozen_chains list
            return self.chains[random.choice(candidate_chain_ids)]
        


    #-----------------------------------------------------------------
    #
    def initialize_type_grid(self):
        """
        Initualizes the type grid using the chain's int_sequence
        types (i.e. int_sequence contains the chain's sequence where
        bead types (letters) are represented by integers)     

        Parameters
        ------------
        None


        Returns
        ------------
        None

        
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
            List of old positions to be removed from the type grid

        new_positions : list
            List of new positions to be added to the type grid

        indices : list
            List of indices in the chain which correspond to the positions

        safe : bool
            If True, will check if the new positions are already occupied by another chain
            and if so will not update the type grid. If False, will update the type grid
            regardless of whether the new positions are occupied by another chain.

        Returns
        ------------
        None

        """
        
        self.delete_chain_from_type_grid(chainID, old_positions, indices, safe)
        self.insert_chain_into_type_grid(chainID, new_positions, indices, safe)

            
    #-----------------------------------------------------------------
    #                
    def delete_chain_from_type_grid(self, chainID, positions, indices, safe=True):
        """
        Function which deletes the positions/indices associated with chainID from the type grid.         

        Indices should be a vector of indices in the chain which correspond to the positions - i.e. if indies were
        [4,5,6,7] then posistions would be a list or array of length 4 where the positions correspond to the positions
        of residues, 4, 5, 6, and 7, respectivly.

        Parameters
        ------------
        chainID : int
            ID of the chain to be deleted

        positions : list
            List of positions to be deleted from the type grid

        indices : list
            List of indices in the chain which correspond to the positions

        safe : bool
            If True, will check if the new positions are already occupied by another chain
            and if so will not update the type grid. If False, will update the type grid
            regardless of whether the new positions are occupied by another chain.

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
        i.e. if indies were [4,5,6,7] then posistions would be a list or array of length 4 
        where the positions correspond to the positions of residues, 4, 5, 6, and 7, 
        respectivly.

        Parameters
        ------------
        chainID : int
            ID of the chain to be inserted

        positions : list
            List of positions to be inserted into the type grid

        indices : list
            List of indices in the chain which correspond to the positions

        safe : bool
            If True, will check if the new positions are already occupied by another chain
            and if so will not update the type grid. If False, will update the type grid
            regardless of whether the new positions are occupied by another chain.
        
        Returns
        ------------
        None
        """

        chain = self.chains[chainID]
        sequence  = chain.int_sequence

        if safe:

            if not len(positions) == len(indices):
                raise TypeGridException(f"Trying to insert ChainID [{chainID}] into the type grid but set of chain indices does not match set of positions to futsz with")

            for i in range(0, len(positions)):                        
                current = lattice_utils.get_gridvalue(positions[i], self.type_grid)

                if not current == 0:
                    raise TypeGridException('Trying to update the type grid  but the current operation would over-write the site at [' + str(positions[i]) + "] with type [%s] when it's currently set to [%s]"%(sequence[i], current))
                
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

        Parameters
        ------------
        None

        Returns
        ------------
        tuple
            3-place tuple with the lattice main grid, the lattice type grid, and a dictionary of the chain positions


        """
        grid_copy = np.copy(self.grid)
        type_grid_copy = np.copy(self.type_grid)
        chain_pos_copy = {} 

        for i in self.chains:
            chain_pos_copy[i] = np.copy(self.chains[i].positions).tolist()

        return (grid_copy, type_grid_copy, chain_pos_copy)


    def lattice_restorefrombackup(self, grid, type_grid, chain_dict):
        """
        Function which restores the lattice and it's encompassed chain objects
        BACK to some former state based on the passed info. NOTE: We do not
        perform any sanity checking at all for this restore function, so if
        you're using it make sure everything is consistent. This is a really
        risky move so PLEASE be carefull...

        Parameters
        ------------
        grid : np.array
            The lattice main grid to be restored

        type_grid : np.array
            The lattice type grid to be restored

        chain_dict : dict
            A dictionary of the chain positions to be restored (keys = chainID,
            values = list of positions).

        Returns
        ------------
        None

        """

        # first, delete the old lattice and type grid
        del self.grid
        self.grid = grid

        del self.type_grid
        self.type_grid = type_grid

        if not len(self.chains) == len(chain_dict):
            raise LatticeInitializationException(latticeExceptions.message_preprocess('Trying to re-set the lattice using lattice_restorefrombackup but there is a mismatch between the number of chains expected and the number provided.'))

        for i in self.chains:
            self.chains[i].positions = chain_dict[i]

        # now, garbage collect to free up memory
        import gc
        gc.collect()


    def write_chain_to_chainid_file(self):
        """
        Function which writes the chain to chainID file. This is a simple
        function which writes the chain to a file with the name chainID.txt
        where chainID is the ID of the chain. This is useful for debugging
        and for checking the state of the chain at any given time.

        Parameters
        ------------
        None

        Returns
        ------------
        None

        """

        with open(OUTPUT_CHAIN_TO_CHAINID, 'w') as fh:
            for chainID in self.chains:
                seq = self.chains[chainID].sequence
                fh.write(f'{chainID}\t{len(seq)}\t{seq}\n')

        
