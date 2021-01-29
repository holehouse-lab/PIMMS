## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2021
## ...........................................................................


import numpy as np
import copy

from . import analysis_structures
from . import lattice_utils
from . import lattice_analysis_utils
from .latticeExceptions import ChainInsertionFailure, ChainInitializationException, ChainAugmentFailure
from .CONFIG import *

class Chain:

    def __init__(self, lattice_grid, dimensions, sequence, int_seq, LR_int_seq, LR_IDX, chainID, chainType, chain_positions=None, fixed=False, rigid=False, center=False, hardwall=False):
        """
        Constructor for a Chain object. In PIMMS, every individual polymer is defined as a Chain, so 
        there will be many, many Chains per simulation

        Parameters
        ----------------
        lattice_grid : np.ndarray 
            Grid where chain is going to be inserted into (with other chains in it!) - used
            to find a vacant position for the chain

        sequence : str
            Human readable sequence for the string

        int_seq : list of ints
            The energy-file encoded sequence where different residues are coded as integers as
            used in their short-range energy interactions. For example, if a chain was a homopolymer
            this would be a list of values equal to the chain length where each value is the same
            number

        LR_int_seq : list of ints
            The energy-file encoded sequence where different residues are coded as integers as
            used in their long-range energy interactions

        LR_IDX : list of ints
            List of indices corresponding to the position(s) in the sequence which undergo long
            range interactions

        chainID : int
            Unique identifier for the chain. Note this value will always be greater than 1 - 
            i.e. the first chainID is 1, and chainIDs monotonically increase thereafter

        chainType : int 
            ID for a specific chain type - i.e. many different chains could be the same type.
            chainType are defined by unique indices starting at 0 and monotonically increasing.

        chain_positions : list of positions {None}
            If present providedm a 'new' chain is initialized in the positsion defined by 
            chain_of_positions. Note that "positions" here are 2D or 3D sublists that define
            X/Y or X/Y/Z coordinates. 

        fixed : bool {False}
            If set to True this chain cannot be moved. Relevant for preformed structures which 
            are non-mobile.        

        rigid : bool {False}
            If the chain is limited to rigid body movements only. This is not yet implemented
            but will be soon.
        
        center : bool {False}
            Defines if we're going to try and place the chain in the center of the box/square.
            default is False.

        hardwall : bool {False}
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC.

        """
        
        # set the chain ID
        self.chainID      = chainID

        # set the chain type ID
        self.chainType    = chainType

        # set the lattice dimensions
        self.dimensions   = dimensions

        # set the sequence associated with this chain
        self.sequence     = sequence

        # set the sequence length
        self.seq_len      = len(sequence)

        # set the coded integer sequence (for use in the type grid) Each value
        # corresponds to a distinct type of bead
        self.int_sequence = int_seq

        # set the coded long-range inter sequence (for use in the type grid over
        # long range
        self.LR_int_sequence = LR_int_seq

        # set if the chain is not to be moved AT ALL. 
        self.fixed = fixed

        # set if the chain is limited to rigid-body movements or not
        self.rigid = rigid

        # define the index of residues which undergo LR interactions
        self.LR_IDX = LR_IDX

        # automatically determine if sequence is a homopolymer
        if len(set(sequence)) == 1:
            self.homopolymer = True
        else:
            self.homopolymer = False

        # if we passed chain positions in....
        if chain_positions:
            # if we're intializing a chain when we aready know it's positions on
            # the grid
            
            # check to make sure we're not trying to initialize the wrong number of
            # positions
            if len(chain_positions) == self.seq_len:
                self.positions = chain_positions
                
            else:
                raise ChainInitializationException('Tried to initialize a chain [chainID = %i] with a sequence of length %i but had %i positions' % (chainID, self.seq_len, len(chain_positions)))

            # should probably have a debug sanity check here...
        else:            
            
            # construct a new chain on the lattice $lattice_grid
            # of the length of this sequence with the chainID
            # chainID and return the possitions associated with
            # this new chain

            # if the center flag was passed as true start the chain in the middle of the grid
            if center:                
                if len(self.dimensions) == 2:
                    default_start = [int(self.dimensions[0]/2), int(self.dimensions[1]/2)]
                else:
                    default_start = [int(self.dimensions[0]/2), int(self.dimensions[1]/2), int(self.dimensions[2]/2)]

                
                try:
                    self.positions    = lattice_utils.insert_chain(chainID, len(sequence), lattice_grid, default_start=default_start, hardwall=hardwall)
                except ChainInsertionFailure:
                    raise ChainInsertionFailure('\nUnable to insert chain %i (length %i) into the center.\nThis is not right, as center-insertion should only be used if a single chain is being added. Please report this...\n' % (chainID, self.seq_len))

            
            else:
                try:
                    self.positions    = lattice_utils.insert_chain(chainID, len(sequence), lattice_grid, hardwall=hardwall)
                except ChainInsertionFailure:
                    raise ChainInsertionFailure('\nUnable to insert chain %i (length %i) into the lattice\nThis is generally indicative of the lattice being overcrowded - you probably have too many chains for the lattice size...\n' % (chainID, self.seq_len))


        ## >>> analysis initialization 
        ##
        self.internal_scaling = analysis_structures.InternalScaling(self.seq_len)
        self.internal_scaling_squared = analysis_structures.InternalScalingSquared(self.seq_len)
        self.distance_map     = analysis_structures.DistanceMap(self.seq_len)
        


    #-----------------------------------------------------------------
    #
    def get_ordered_positions(self):
        """
        Returns a list of the chain positions

        """
        return self.positions

    #-----------------------------------------------------------------
    #
    def get_intcode_sequence(self):
        """
        Returns a list where each position corresponds to the integer code
        used by the energy calculations to identify a specific residue type

        """
        
        return self.int_sequence



    #-----------------------------------------------------------------
    #
    def get_single_image_positions(self):
        """
        Returns a list of the chain positions where all positions
        in the chain are corrected to lie in the same periodic image (i.e. 
        "single image convention" as opposed to minimum image convention.)

        This has the nice feature of being able to deal with an arbitrary 
        number of periodic images, so if the chain spans many PBCs
        (e.g. imagine a chain that extends out of its main box through
        another box and INTO another box) this can deal with that.                
        """

        # if the chain does not straddle a periodic boundary
        if not self.does_chain_stradle_pbc_boundary():
            return self.positions
        else:
            return lattice_utils.convert_chain_to_single_image(self.positions, self.dimensions)


    #-----------------------------------------------------------------
    #
    def does_chain_stradle_pbc_boundary(self):
        """
        Determines if the chain straddles a periodic boundary or not
        """

        return lattice_utils.do_positions_stradle_pbc_boundary(self.positions)
        
      
    #-----------------------------------------------------------------
    #
    def get_LR_positions(self):
        """
        Returns a list of chain positions which engage in long range interactions 

        """
        return [self.positions[i] for i in self.LR_IDX]



    #-----------------------------------------------------------------
    #
    def get_LR_binary_array(self):
        """
        Returns a numpy array of chain length, where beads that engage in long
        range interactions are set to 1 and all others are set to 0

        """
        return_list = []
        for i in range(0, self.seq_len):
            if i in self.LR_IDX:
                return_list.append(1)
            else:
                return_list.append(0)
        return np.array(return_list, dtype=int)



    #-----------------------------------------------------------------
    #
    def get_positions_by_chain_index(self, index_list):
        """
        Returns a list of lattice positions based on the index positions
        in the index_list

        """
        position_list = []
        for i in index_list:
            position_list.append(self.positions[i])

        return position_list


    #-----------------------------------------------------------------
    #
    def get_positions_by_chain_index_single_image_position(self, index_list):
        """
        Returns a list of lattice positions based on the index positions
        in the index_list using the single image convention (i.e. all positions
        come from a chain that exists in one periodic dimension without crossing
        a PBC boundary

        """
        SIP = lattice_utils.convert_chain_to_single_image(self.positions, self.dimensions)

        position_list = []
        for i in index_list:
            position_list.append(SIP[i])

        return position_list



    #-----------------------------------------------------------------
    #                    
    def set_ordered_positions(self, positions):
        if len(positions) == self.seq_len:
            self.positions = positions
        else:
            raise ChainAugmentFailure('Tried to set chain %i to a set of positions of length %i, but this chain is actually %i residues long' %(self.chainID, len(positions), len(self.positions)))



    #-----------------------------------------------------------------
    #
    def get_center_of_mass(self, on_lattice=True):
        """
        Returns a chain's center of mass (note for now we assume every bead has the same mass such that
        the center of mass ends up becoming the mean position in all 2 or 3 dimensions (depending on the system) 

        """
        return lattice_utils.center_of_mass_from_positions(self.get_ordered_positions(), self.dimensions)


    #####################################################################################################
    ## INTERNAL SCALING ANALYSIS FUNCTIONS
    ##

    def analysis_update_internal_scaling(self):        
        """
        Function which when called will re-evaluate the chain's current internal scaling information based on its current position and then 
        update the local internal_scaling object to include this most recent analysis.

        Note this does not REPLACE the current internal scaling information, but allows a running average which should become more accurate
        the more frequently the function is called.

        """

        ij_dist = {}

        # for each possible gap size
        for gap in range(1, self.seq_len-1):
            tmp_dis = []

            # cycle over all pairs equal to the gap size and calculate the average distances
            for i in range(0, self.seq_len-gap):
                tmp_dis.append(lattice_analysis_utils.get_inter_position_distance(self.positions[i], self.positions[i+gap], self.dimensions))
            ij_dist[gap] = np.mean(tmp_dis)

        # finally update the internal value (note the squaring is done inside the
        # internal scaling squared function!)
        self.internal_scaling.update_internal_scaling(ij_dist)
        self.internal_scaling_squared.update_internal_scaling(ij_dist)


    #-----------------------------------------------------------------
    #
    def analysis_print_internal_scaling(self):
        """

         """
        self.internal_scaling.print_status()


    #-----------------------------------------------------------------
    #
    def analysis_get_internal_scaling(self):
        """

        """
        return self.internal_scaling.get_internal_scaling_array()        


    #-----------------------------------------------------------------
    #
    def analysis_print_internal_scaling_squared(self):
        """

         """
        self.internal_scaling_squared.print_status()


    #-----------------------------------------------------------------
    #
    def analysis_get_internal_scaling_squared(self):
        """

        """
        return self.internal_scaling_squared.get_internal_scaling_array()        



    #-----------------------------------------------------------------
    #
    def analysis_fit_scaling_exponent(self):
        """

        """
        return self.internal_scaling_squared.fit_scaling_exponent()


    #####################################################################################################
    ## DISTANCE MAP ANALYSIS FUNCTIONS
    ##

    def analysis_update_distance_map(self):        

        # initialize the empty distance map
        distance_map = np.zeros((self.seq_len, self.seq_len),dtype=float)

        # build the upper triangle distance map
        for i in range(0,self.seq_len):
            for j in range(0+i, self.seq_len):
                distance_map[i][j] = lattice_analysis_utils.get_inter_position_distance(self.positions[i], self.positions[j], self.dimensions)

        self.distance_map.update_distance_map(distance_map)

    #-----------------------------------------------------------------
    #
    def analysis_get_distance_map(self):        
        """
        Returns the chain's current distance map

        """
        return self.distance_map.get_distance_map()


    #####################################################################################################
    ## END TO END DISTANCE ANALYSIS FUNCTIONS
    ##

    def analysis_get_end_to_end_distance(self):        
        """
        Returns the chain's current end-to-end
        distance on the lattice

        """

        start  = self.positions[0]
        end    = self.positions[-1]

        return lattice_analysis_utils.get_inter_position_distance(start, end, self.dimensions)


    #####################################################################################################
    ## Positional analysis
    ##

    def analysis_get_residue_residue_distance(self, R1, R2):        
        """
        Returns the inter-residue position as defined by the two
        positions here

        """

        start  = self.positions[R1]
        end    = self.positions[R2]

        return lattice_analysis_utils.get_inter_position_distance(start, end, self.dimensions)


        
    #####################################################################################################
    ## RADIUS OF GYRATION ANALYSIS FUNCTIONS
    ##

    def analysis_get_radius_of_gyration(self):        
        return lattice_analysis_utils.get_polymeric_properties(self.positions, self.dimensions)[0]
        #return lattice_analysis_utils.get_polymeric_properties(self.get_single_image_positions(), self.dimensions)[0]
        


    def analysis_get_polymeric_properties(self):                        
        """
        Function that returns the polymeric properties associated with a chain. These properties are in a list so additional properties can be
        added as we go through.

        [0] - Radius of gyration
        [1] - Asphericity
        
        ## NOTE: Finite size detection!
        I implemented a method of extract self-consistent chain positions such that all positions come from the same periodic image, regardless
        of how many images the chain actually expands over. 

        This actually only ends up making a tangible difference when you have chains that extend over and span multiple boxes, at which 
        point you're probably in some serious trouble. As a result, we use the standard naive PBC correction approach for all analysis
        BUT this warn function lets you explicitly check if you're in a regime where finite size artefacts might be a problem. This is
        called at the same frequency.

        Note that right now this 'single image convention' algorithm is ONLY used here for checking up on finite size artefacts. However,
        it's fully functional and not too expensive so could be used to replaced the normal way of getting positions if needed be. Presumably
        the single image convention algorithm has been implemented by someone else somewhere but this was a naieve implementation I developed
        and it seems to work well.       

        """

        # do both for now, maybe add this as a toggle switch in the future as performance becomes more key
        polymeric_props = lattice_analysis_utils.get_polymeric_properties(self.positions, self.dimensions)
        single_image_PBC_props = lattice_analysis_utils.get_polymeric_properties(self.get_single_image_positions(), self.dimensions)

        if abs(polymeric_props[0] - single_image_PBC_props[0]) > 0.001:
            print("\n[WARNING]: Computing the radius of gyration using minimum image convention vs. single image convention yeilded different results [%3.5f vs %3.5f]. This probably suggests you're experiencing substantial finite size artefacts and should rethink the box-size to chain dimensions. This message will appear everytime this issue is noticed.\n" % (polymeric_props[0], single_image_PBC_props[0]))

        if abs(polymeric_props[1] - single_image_PBC_props[1]) > 0.001:
            print("\n[WARNING]: Computing the asphericity using minimum image convention vs. single image convention yeilded different results [%3.5f vs %3.5f]. This probably suggests you're experiencing substantial finite size artefacts and should rethink the box-size to chain dimensions. This message will appear everytime this issue is noticed.\n" % (polymeric_props[1], single_image_PBC_props[1]))

        return polymeric_props
            

        #return lattice_analysis_utils.get_polymeric_properties(self.get_single_image_positions(), self.dimensions)



            
                
                

            
        
    




        

                                
            


            


        
            
            
