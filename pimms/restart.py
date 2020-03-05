## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
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


class RestartObject:
    """
    Object used to read and write restart files. Restart information ONLY contains information on
    chain position, sequence, and type, and grid dimenisons, but does NOT include any information     

    """


    #-----------------------------------------------------------------
    #       
    def __init__(self):
        """
        Initialization function to create an empty RestartObject
        """
        self.energy = 0
        self.dimensions = 0
 

    #-----------------------------------------------------------------
    #       
    def set_energy(self, energy):
        """
        Set the RestartObject's energy value
        """
        self.energy = energy


    #-----------------------------------------------------------------
    #       
    def build_from_lattice(self, LATTICE, hardwall=False):
        """
        Construct a restart object using a lattice object to set the chain
        positions

        """
        self.dimensions = LATTICE.dimensions
        self.hardwall   = hardwall

        self.chains = {}
        for chainID in LATTICE.chains:
            self.chains[chainID] = [LATTICE.chains[chainID].positions, LATTICE.chains[chainID].sequence, LATTICE.chains[chainID].chainType]



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
        Function that updates the restart object's dimenions AND moves the chains so 
        they're centered in the new lattice.
        """
        
        ## -----------
        # calculate offset so the chains are placed in the center of the new lattice
        x_off = int((new_dimensions[0] - self.dimensions[0])/2)
        y_off = int((new_dimensions[1] - self.dimensions[1])/2)

        if len(new_dimensions) ==3:
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
    def __apply_position_offset(self, position_offset):
        """
        Function that allows position of each residue to be offset by some fixed amount. This is not relevant for traditional restart
        operations, but is useful when using a Restart object to initialize a new (resized) lattice. This requires that the restart
        object dimensions are big enough to contain the newly offset positions.
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
            raise RestartException('Trying to apply position offset to restart object, but dimensions do not match')

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
                        raise RestartException('Trying to offet a position on chain %i from %i to %i (dim=%i) but lattice dimensions are %s' % (chainID, position[dim], position[dim] + position_offset[dim], dim, str(self.dimensions)))
                        


    #-----------------------------------------------------------------
    #       
    def build_from_file(self, filename):
        """
        Function that constructs a restart object from a passed filename. Performs some sanity check
        in reading in the file but doesn't actually check that the chain positions make sense on 

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

        self.chains = {}
        for chainID in local_chains:
            self.chains[chainID] = local_chains[chainID]

        
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



    
