## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................

##
## restart
##
## The RestartObject implemements a way to read and write restart files. This allows PIMMS
## to restart from previous simulations. Other than chain position, no other state is saved.  
##


import pickle

from . import CONFIG
from .latticeExceptions import RestartException
from . import pimmslogger


class RestartObject:
    """
    Object used to read and write restart files. Restart information ONLY contains information on
    chain position, sequence, and type, and grid dimenisons, but does NOT include any information     

    Note that the self.chains object in a RestartObject has the following structure:

    1. Is a dictionary 
    2. Keys are chainID (i.e. each seperate chain has it's own entry)
    3. values is a list with three elements
       [0] : bead positions (N->C)
       [1] : chain sequence (which will be referenced against the parameter file)
       [2] : chainType : a single value that defines the type of chain

    """


    #-----------------------------------------------------------------
    #       
    def __init__(self):
        """
        Initialization function to create an empty RestartObject
        """
        self.energy = 0
        self.dimensions = 0
        self.extra_chains = {}


    #-----------------------------------------------------------------
    #       
    def __apply_position_offset(self, position_offset):
        """
        Function that allows position of each residue to be offset by some fixed amount. 
        This is not relevant for traditional restart operations, but is useful when using 
        a Restart object to initialize a new (resized) lattice. This requires that the 
        restart object dimensions are big enough to contain the newly offset positions.

        Parameters
        --------------
        position_offset : list
            List of integers with length equal to the number of dimensions in the lattice. Each element
            in the list defines the amount by which the position of each residue should be offset.

        Returns
        -------------
        None
            No return type, but the internal self.chains object will be appropriately
            updated.

        Raises
        -------------
        RestartException
            If the dimensions of the restart object are not big enough to contain the 
            newly offset positions.

        
        """
        
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        # Internal function that tests if a position (pos) in dimension (dim) is valid given the
        # restart lattice' dimensions
        def valid_pos(pos, dim):
            pos = pos+position_offset[dim]
            if (pos < 0) or (pos >= self.dimensions[dim]):
                return False
            else: 
                return True
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
        # check offset dimensions match chain dimensions
        if len(position_offset) != len(self.dimensions):
            raise RestartException('Trying to apply position offset to restart object, but dimensions do not match.')

        n_dim = len(self.dimensions)
        # also requires that all new positions 

        # for each chain
        for chainID in self.chains:
            # for each position in each chain
            for position in self.chains[chainID][0]:

                # for each dimension check that this new position would be in
                # bounds and throw exception if not
                for dim in range(0, n_dim):
                    if valid_pos(position[dim], dim):
                        position[dim] = position[dim] + position_offset[dim]
                    else:
                        raise RestartException(f'Trying to offet a position on chain {chainID} from {position[dim]} to {position[dim] + position_offset[dim]} (dim={dim}) but lattice dimensions are {self.dimensions}')
                        


    #-----------------------------------------------------------------
    #       
    def __update_seq2chainType(self, local_chainType, local_seq, log):
        """
        Internal function called by both build_from_lattice() and build_from_file()
        which ensures an updated and dynamically constructed self.seq2chainType dictionary
        exists which enables mapping of protein sequence to a chainType.

        Note we don't not allow two identical chain sequences to have different chainTypes - there
        are some circumstances where this might be preferable, so, the seq2chainType mapping
        is a one-to-many mapping, although IN GENERAL we probably expect this mostly to be
        a 1-to-1 mapping.

        If a one-to-many mapping is found and log=True then this is written via the pimmslogger as a 
        warning 

        Parameters
        --------------
        local_chainType : int
            The chainType associated with the passed chain

        local_seq : str
            The amino acid sequence of the passed chain.

        log : bool
            Flag which, if set to true, means if this seq already has a chainType defined but the 
            passed chainType is a DIFFERENT value it'll warn the user about this.

        Returns
        -------------
        None
            No return type, but the internal self.seq2chainType dictionary will be appropriately
            updated

        """
        
        # if we've seen this sequence before
        if local_seq in self.seq2chainType:

            # If the chainType assigned here is associated with that previous record
            # move on...
            if local_chainType in self.seq2chainType[local_seq]:
                pass
            else:
                self.seq2chainType[local_seq].append(local_chainType)

                # note this is not strictly a problem, just might be good to know about...
                if log:
                    pimmslogger.log_warning(f'When building RestartObject from Lattice found two identical chains [{local_seq}] with different chainType indices. This is not a bug or problem, but may be undesired...')

        # if we've never seen this sequence before this is easy...
        else:
            self.seq2chainType[local_seq] = [local_chainType]



    #-----------------------------------------------------------------
    #           
    def add_extra_chains(self, extra_chains):
        """
        Function which allows extra chains (as read from a keyfile) to be
        added to a RestartObject so that when a new lattice is initialized
        from this RestartObject those extra chains are randomly placed
        somewhere across the simulation box.

        Note extra_chains ONLY have a sequence and chainType associated
        with them, but do NOT have any positions.

        Parameters
        ----------------
        extra_chains : list
            List with two elements
            [0] = number of chains (int)
            [1] = chain sequence (str)

        Returns
        ----------------
        None
            No return type, but the internal self.extra_chains dictionary 
            will be appropriately updated.
        

        """
        # dynamically calculate what next chainID should be. This is done by determining the max
        # chainID in both the chains and extrachains dictionaries, and then seeting the NEW
        # chainID to 1 + that number
        if len(self.extra_chains) == 0:
            chainID = max(list(self.chains.keys())) + 1
        else:
            chainID = max(max(list(self.chains.keys())), max(list(self.extra_chains.keys())))  + 1

        # extract info and raise exception in a civilized way
        try:
            count = extra_chains[0]
            chain_seq = extra_chains[1]
        except Exception:
            raise RestartException('ERROR parsing EXTRA_CHAINS keyword [{extra_chains}] - could not parse into chain count and chain sequence')
            

        if chain_seq in self.seq2chainType:

            # note - this [0] means we always use the first chain type even if there are multiple
            # chain IDs associated with a specifi chain. 
            local_chainType = self.seq2chainType[chain_seq][0]
        else:

            # if a new chain dynamically calculate what the next chainType should be (next increment
            # after current highest number)
            tmp = []
            for s in self.seq2chainType:
                tmp.extend(self.seq2chainType[s])

            local_chainType = max(tmp) + 1
        
        # finally, after all this set up, add to the extra_chains dict
        for c in range(count):
            self.extra_chains[chainID] = [None, chain_seq, local_chainType]
            chainID = chainID + 1
 

    #-----------------------------------------------------------------
    #       
    def set_energy(self, energy):
        """
        Set the RestartObject's energy value
        """
        self.energy = energy


    #-----------------------------------------------------------------
    #       
    def build_from_lattice(self, LATTICE, hardwall=False, log=False):
        """
        Construct a restart object using a lattice object to set the chain
        positions.

        Parameter
        ------------
        LATTICE : pimms.lattice.Lattice 
            A standard PIMMS latticd object

        hardwall : bool (default = False)
            Flag which sets of the current system defines a hardwall or, if false
            PBC.

        log : bool (default = False)
            Flag which if set to True means warnings are written to the standard PIMMS
            logfile

        Returns
        ----------
        None
            No return type, but the internal self.chains dictionary will be appropriately
            updated.


        """
        self.dimensions = LATTICE.dimensions
        self.hardwall   = hardwall

        # reset chain info...
        self.chains = {}
        self.seq2chainType  = {}

        for chainID in LATTICE.chains:
        
            local_chainType = LATTICE.chains[chainID].chainType
            local_seq = LATTICE.chains[chainID].sequence

            # add the chain to the restart object
            self.chains[chainID] = [LATTICE.chains[chainID].positions, local_seq, local_chainType]

            # udpate the self.seq2chainType dictionary
            self.__update_seq2chainType(local_chainType, local_seq, log)


    #-----------------------------------------------------------------
    #       
    def set_dimensions(self, dimensions):
        """
        Function that allows the dimensions to be overridden. Only needed
        if we're actually changing the lattice size.
        """
        self.dimensions = dimensions


    #-----------------------------------------------------------------
    #       
    def update_lattice_dimensions(self, new_dimensions):
        """
        Function that updates the restart object's dimensions AND moves the chains so 
        they're centered in the new lattice.
        """
        
        ## -----------
        # calculate offset so the chains are placed in the center of the new lattice
        x_off = int((new_dimensions[0] - self.dimensions[0])/2)
        y_off = int((new_dimensions[1] - self.dimensions[1])/2)

        if len(new_dimensions) == 3:
            z_off = int((new_dimensions[2] - self.dimensions[2])/2)
            position_offset=[x_off, y_off, z_off]
        else:
            position_offset=[x_off, y_off]
            ## -----------

        # Next construct and instantiate a new restart object which has the new dimensions
        # including applying the possition offset we calculated above

        # next update the restart object's dimensions
        self.dimensions = new_dimensions
        
        # finaly, apply the offset on this 'new' lattice (order matters, as __apply_position_offset
        # assesses if, given self.dimensions, the offset is valid or not)
        self.__apply_position_offset(position_offset)






    #-----------------------------------------------------------------
    #       
    def build_from_file(self, filename, log=False):
        """
        Function that constructs a restart object from a passed filename. Performs some sanity check
        in reading in the file but doesn't actually check that the chain positions make sense on the 
        lattice. We can and should probably make this better going forwards...

        Parameters
        --------------
        filename : str
            Name of the file to be read

        log : bool (default = False)
            Flag which if set to True means warnings are written to the standard PIMMS
            logfile

        Returns
        -------------
            None but updates the current object to contain self.dimensions, self.energy, self.hardwall 
            and self.chains[] info.

        """


        # if IO issue (not IndexError often thrown if a valid file is found
        # but its not actually a pickle file!
        try:
            input_dict = pickle.load( open(filename, "rb" ) )
        except (IOError, IndexError) as e:
            raise RestartException("Error reading restart file. Error:\n\n%s" %(str(e)))
        
        # extract out key info (throw exception if missing)
        try:
            self.dimensions = input_dict['DIMENSIONS']
            self.energy     = input_dict['ENERGY']
            self.hardwall   = input_dict['HARDWALL']
            local_chains    = input_dict['CHAINS']
        except KeyError as e:
            raise RestartException("Invalid restart file - missing entry for %s" % (e.args[0]))

        # reset chain info...
        self.chains = {}
        self.seq2chainType  = {}

        for chainID in local_chains:

            local_seq = local_chains[chainID][1]
            local_chainType = local_chains[chainID][2]

            self.chains[chainID] = local_chains[chainID]

            # update the self.seq2chainType dictionary
            self.__update_seq2chainType(local_chainType, local_seq, log)


            

        
    #-----------------------------------------------------------------
    #       
    def write_to_file(self):

        output={}
        output['CHAINS'] = {}
        for chainID in self.chains:
            output['CHAINS'][chainID] = self.chains[chainID]

        output['DIMENSIONS'] = self.dimensions
        output['ENERGY']     = self.energy
        output['HARDWALL']   = self.hardwall

        pickle.dump( output, open( CONFIG.RESTART_FILENAME, "wb" ) )



    
