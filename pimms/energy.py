## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

import numpy as np

from . import lattice_utils
from . import lattice_analysis_utils
from .latticeExceptions import EnergyException, ParameterFileException
from . import parameterfile_parser
from . import hyperloop
from . import longrange_utils

class EmptyHamiltonian:
    """
    Dummy class which implements a Hamiltonian used when all interactions
    are turned off. Means we can define whatever stub-functionality here
    without adding code-rot to the true Hamiltonian class

    """

    def __init__(self):

        self.LR_residue_names = []

    def evaluate_total_energy(self, x):
        return 0.0

    def evaluate_local_energy(self, x, y):
        return 0.0

    def evaluate_local_energy_LR(self, x, y):
        return 0.0

    def evaluate_angle_energy(self, x,y):
        return 0.0

    def convert_sequence_to_integer_sequence(self, sequence):
        """
        See the true Hamiltonian class for what this function is
        really doing...

        """
        return [1]*len(sequence)

    def convert_sequence_to_LR_integer_sequence(self, sequence):
        """
        See the true Hamiltonian class for what this function is
        really doing...

        """
        return []

    def get_indices_of_long_range_residues(self, sequence):
        """
        See the true Hamiltonian class for what this function is
        really doing...

        """
        return []


class Hamiltonian:
    """
    Main Hamiltonian class for evaluating system energies

    """

    def __init__(self, parameter_file, num_dimensions, non_interacting, angles_off, hardwall=False, temperature=False):
        """


        human_readable_interaction_table  

        residue_names 

        human_readable_LR_interaction_table

        LR_residue_names                    

        residue_interaction_table

        parameter_to_int_map

        LR_residue_interaction_table

        LR_parameter_to_int_map



        """

        # internally hardwall is an integer
        if hardwall:
            self.hardwall = 1
        else:
            self.hardwall = 0

        # read in and parse the parameter file  - this generates a HUMAN readable table
        # where particles are represent by strings as defined in the interaction table, and the full set
        # of particles defined in that table 
        (self.human_readable_interaction_table, self.residue_names, self.human_readable_LR_interaction_table, self.LR_residue_names, self.human_readable_SLR_interaction_table) = parameterfile_parser.parse_energy(parameter_file) 

        # now build the particle interaction table, where we re-code the particle names as integers - this allows
        # our energy functions to use cython code (optimized cython doesn't support string arrays) which buys
        # some significant performance enhancement
        (self.residue_interaction_table, self.parameter_to_int_map, self.LR_residue_interaction_table, self.LR_parameter_to_int_map, self.SLR_residue_interaction_table) = self.build_interaction_table(non_interacting)

        # Finally extract residue specific angle potentials - note this REQUIRES every residue have a defined angle potential
        if angles_off:
            self.build_angle_interactions(False, num_dimensions, angles_off)
        else:
            angle_dict = parameterfile_parser.parse_angles(parameter_file, temperature)
            self.build_angle_interactions(angle_dict, num_dimensions, angles_off)
            parameterfile_parser.write_angle_parameter_summary(angle_dict, parameter_file)
        

    def set_hardwall(self, value=True):
        if value:
            self.hardwall=1
        else:
            self.hardwall=0



    def evaluate_total_energy(self, latticeObject, id_to_typeMap=None):
        """
        Function which evaluates the total energy of the system.

        This function constructs a redundant list of ALL pairwise interactions between each residue
        and every neighbour, and then uses the evaluate_local_energy to define the energy of the 'local'
        system defined by those pairs (where the local system happens to actually be the entire
        system)

        """
                
        all_positions    = []
        all_LR_positions = []
        LR_binary_array  = np.array([], dtype=int)

        # cycle through each chain extracting the positions to generate a long list
        # associated with the positions of every residue on the lattice
        angle_energy = 0

        # for each chain extract out the positions enaging 
        for chainID in latticeObject.chains:
            all_positions.extend(latticeObject.chains[chainID].get_ordered_positions())
            all_LR_positions.extend(latticeObject.chains[chainID].get_LR_positions())

            # replace this with a list which gets concatonated in a single operation at the end - right now this is not
            # efficient as new memory must be allocated on each loop - BOOO!
            LR_binary_array = np.concatenate((LR_binary_array, latticeObject.chains[chainID].get_LR_binary_array()))
            angle_energy = angle_energy + self.evaluate_angle_energy(latticeObject.chains[chainID].get_ordered_positions(), latticeObject.chains[chainID].get_intcode_sequence(), latticeObject.dimensions)
                
        # build the non-redundant set of pairs for long-range and short range interactions
        (pairs, lr_pairs, slr_pairs) = lattice_utils.build_all_envelope_pairs(all_positions, LR_binary_array, latticeObject.type_grid, latticeObject.dimensions)


        # evaluate the energy associated with all those pairs
        energy_local = self.evaluate_local_energy(latticeObject, pairs) 
        energy_LR    = self.evaluate_local_energy_LR(latticeObject, lr_pairs) 
        energy_SLR   = self.evaluate_local_energy_SLR(latticeObject, slr_pairs) 
        energy_angle = angle_energy
        
        total = energy_local + energy_LR + energy_SLR + energy_angle

        # sum all the energy and return
        return (total, energy_local, energy_LR, energy_SLR, angle_energy)


    #-----------------------------------------------------------------
    #    
    #
    def evaluate_local_energy(self, latticeObject, pairs_list):
        """
        This is really the main energy calculating function - it takes a latticeObject (which
        contains the comple information on what residues are where) and a pairs_list which
        defines the pairs of residues that the energy will be calculated over (i.e. defining
        the 'locality' of this operation - LOCAL does not here mean only short range!).

        """
        return self.__evaluate_local_energy_shortrange(latticeObject, pairs_list, self.residue_interaction_table)
            

      
    #-----------------------------------------------------------------
    #    
    #
    def evaluate_local_energy_LR(self, latticeObject, pairs_list):
        """
        This is really the main energy calculating function - it takes a latticeObject (which
        contains the comple information on what residues are where) and a pairs_list which
        defines the pairs of residues that the energy will be calculated over (i.e. defining
        the 'locality' of this operation - LOCAL does not here mean only short range!).

        """
        return self.__evaluate_local_energy_non_shortrange(latticeObject, pairs_list, self.LR_residue_interaction_table)
        


    #-----------------------------------------------------------------
    #    
    #
    def evaluate_local_energy_SLR(self, latticeObject, pairs_list):
        """
        This is really the main energy calculating function - it takes a latticeObject (which
        contains the comple information on what residues are where) and a pairs_list which
        defines the pairs of residues that the energy will be calculated over (i.e. defining
        the 'locality' of this operation - LOCAL does not here mean only short range!).

        """
        return self.__evaluate_local_energy_non_shortrange(latticeObject, pairs_list, self.SLR_residue_interaction_table)



    #-----------------------------------------------------------------
    #    
    #
    def __evaluate_local_energy_shortrange(self, latticeObject, pairs_list, interaction_table):
        """
        Internal general energy evaluation function for short-range interactions.

        """

        num_dims = len(latticeObject.dimensions)

        # if no pairs included return 0
        if len(pairs_list) == 0:
            return 0
        
        # if we're working with a 2D lattice
        if num_dims == 2:                    
            return  hyperloop.evaluate_local_energy_2D_shortrange(latticeObject.type_grid, pairs_list, interaction_table, self.hardwall)
        
        # if we're working with a 3D lattice
        if num_dims == 3:
            return hyperloop.evaluate_local_energy_3D_shortrange(latticeObject.type_grid, pairs_list, interaction_table, self.hardwall)


    #-----------------------------------------------------------------
    #    
    #
    def __evaluate_local_energy_non_shortrange(self, latticeObject, pairs_list, interaction_table):
        """
        Internal general energy evaluation function for long-range and super-long range interactions (i.e.
        not short range interactions).

        """

        num_dims = len(latticeObject.dimensions)

        # if no pairs included return 0
        if len(pairs_list) == 0:
            return 0
        
        # if we're working with a 2D lattice
        if num_dims == 2:                    
            return  hyperloop.evaluate_local_energy_2D_non_shortrange(latticeObject.type_grid, pairs_list, interaction_table, self.hardwall)
        
        # if we're working with a 3D lattice
        if num_dims == 3:
            return hyperloop.evaluate_local_energy_3D_non_shortrange(latticeObject.type_grid, pairs_list, interaction_table, self.hardwall)

 
    #-----------------------------------------------------------------
    #    
    def evaluate_angle_energy(self, chain_positions, intcode_sequence, dimensions):
        """
        Angle energies are determined by hyperloop functions that basically
        compare the angle vector and use a pre-computed lookup table to convert 
        that angle into some energy penalty. It's lightning fast, and residue
        specific!

        chain_positions is a list of 2D or 3D bead positoins we're going to evaluate over
        intcode_sequence is an equal length list of the bead integercodes being used
        dimensions is a list of the box dimensions

        """

        num_positions = len(chain_positions)

        # cannot compute an angle for a chain with only 1 or 2 beads
        if num_positions < 3:
            return 0.0
        
        # for the 3D case
        if len(dimensions) == 3:
            penalty = hyperloop.evaluate_angle_energy_3D(np.array(chain_positions), np.array(intcode_sequence), self.angle_lookup, num_positions)

        # for the 2D case
        else:
            penalty = hyperloop.evaluate_angle_energy_2D(np.array(chain_positions), np.array(intcode_sequence), self.angle_lookup, num_positions)

        return penalty
                
            

    #-----------------------------------------------------------------
    #    
    def build_interaction_table(self, non_interacting=False):
        """
        Carries out dynamic construction of a two AxA float matrix where indicies along
        X and Y axis correspond to residues defined in the parameter file.

        * The residue interaction table (RIT) defines short-range interactions

        * The LRRIT (long range interaction table) defines long-range interactions

        The return values are as follows:

        RIT - 2D numpy array of floats which describes the short range interactions
              The matrix is indexed using integer codes, where each code maps to
              a specific residue type


        LRRIT - 2D numpy array of floats which describes the long range interactions
                 The matrix is indexed using integer codes, where each code maps to
                 a specific residue type (same mapping as the RIT). The LR_RIT and 
                 the RIT are the same size (which allows the same indexing codes to 
                 be used) but MANY residues will not engage in LR interactions, so
                 those sites are set to np.NaN such that if a bug leads to these being
                 used an error will occur.
            

        MAPPING - a residue-to-code mapping allowing for a residue name to be mapped
                  to its integer code (same code for short range and long range). This
                  is a dictionary, where the key is the residue name (string) and the
                  value is the corresponding integer code

        
        LR_MAPPING - mapping where ONLY LR residues are included. This is a subset of
                     the mapping defined in MAPPING but has the nice feature of only
                     including LR residues (so also offers a simple way to determine
                     if a given residue undergoes LR interactions or not).


        SLRRIT - 2D numpy array of floats which describes the super long range 
                 interactions (SLR). The matrix is indexed using integer codes, where 
                 each code maps to a specific residue type (same mapping as the RIT
                 and the LR_RIT).
                 
        """

        # number of different residue types we're messing with
        n_residues    = len(self.residue_names)

        # initialize the residue interaction table [RIT] as a matrix
        # of zeros
        RIT   = np.zeros(shape=(n_residues,n_residues),dtype=float)
        LRRIT = np.zeros(shape=(n_residues,n_residues),dtype=float)
        SLRRIT = np.zeros(shape=(n_residues,n_residues),dtype=float)

        MAPPING = {}
        LR_MAPPING = {}
        
        # Build int-index interaction matrices for short-range and long-range
        # residue interactions 
        # for each residue (human name)
        ## WARNING:
        ## THIS WORKS BECAUSE WE ASSUME THE FIRST RESIDUE IN
        ## RESIDUE NAMES IS SOLVENT - this means 0 is always
        ## salt
        R1_int = 0       
        for R1 in self.residue_names:
            R2_int = 0            

            # mapping defines how we map the human names to integers
            MAPPING[R1] = R1_int

            # if this residue participates in LR interactions building
            # the LR int-to-residue mapping with the SAME integer code
            if R1 in self.LR_residue_names:
                LR_MAPPING[R1] = R1_int

                
            for R2 in self.residue_names:

                # set short range interaction first
                RIT[R1_int][R2_int] = self.human_readable_interaction_table[R1][R2]

                # if long-range interaction between these two residues, set the values (note the SLR
                # is always defined by PIMMS EVEN if it wasn't provided by the parameter file - in 
                # this case the SLR is set to zero)
                if R2 in self.LR_residue_names and R1 in self.LR_residue_names:
                    LRRIT[R1_int][R2_int]  = self.human_readable_LR_interaction_table[R1][R2]
                    SLRRIT[R1_int][R2_int] = self.human_readable_SLR_interaction_table[R1][R2]

                # else 
                else:
                    LRRIT[R1_int][R2_int] = 0
                    SLRRIT[R1_int][R2_int] = 0

                # If the NON_INTERACTING flag is on, then overwrite and set all interactions to zero, but
                # warn about this (!)
                if non_interacting:
                    print("!!! WARNING !!! : This is a non-interacting run (over-riding parameter file for [%s-%s])" % (R1,R2))
                    RIT[R1_int][R2_int] = 0
                    LRRIT[R1_int][R2_int] = 0
                    SLRRIT[R1_int][R2_int] = 0

                R2_int = R2_int + 1

            R1_int = R1_int + 1

        return (RIT, MAPPING, LRRIT, LR_MAPPING, SLRRIT)



    #-----------------------------------------------------------------
    #    
    def build_angle_interactions(self, angle_dict, num_dimensions, angles_off):
        """
        Function that constructs a lookup table that we use to assign 'angle' withstraints. The angle effects are really related
        to the 1_3 interaction, but can also be used 

        """

        # first verify that every residue in the interaction table has an angle (note if we have extra angles
        # these are just ignored). At the same time build up a 
        resnames_through_space = list(self.parameter_to_int_map.keys())

        # if angles are off then get 'angle' names through the interactions
        # and set all to 0
        if angles_off:
            resnames_angles = list(self.parameter_to_int_map.keys())
        else:
            resnames_angles       = list(angle_dict.keys())

        int_to_penalty = {}
        
        for resname in resnames_through_space:

            # skip solvent!
            if resname is '0':
                pass

            else:
                if resname not in resnames_angles:
                    raise ParameterFileException('Residue %s has interaction energies defined but *no angle energies defined*' % resname)

                # if angles are off set all penalties to 0
                if angles_off:
                    print("!!! WARNING !!! : Angles are turned off (over-riding parameter file for %s angles)" % (resname))
                    int_to_penalty[self.parameter_to_int_map[resname]] = [0,0,0]
                else:
                    int_to_penalty[self.parameter_to_int_map[resname]] = angle_dict[resname]

                
        # int_list is a sorted list of the intergers that map to a residue-specific angle pair
        int_list = list(int_to_penalty.keys())
        int_list.sort()
        
        # for the 3D case
        if num_dimensions == 3:

            ## IDX entries for angle lookup are
            ## 1 : residue type
            ## 2 : dx of -1 res
            ## 3 : dy of -1 res
            ## 4 : dz of -1 res
            ## 5 : dx of +1 res
            ## 6 : dy of +1 res
            ## 7 : dz of +1 res
            ##
            ## while the actulat int associated with the intidx reflects the identity of residue i


            self.angle_lookup = np.zeros((int_list[-1]+1, 3, 3, 3, 3, 3, 3))

            AP1_count=0
            AP2_count=0
            AP3_count=0
        
            for x1 in range(-1,2):
                for y1 in range(-1,2):
                    for z1 in range(-1,2):
                        for x2 in range(-1,2):
                            for y2 in range(-1,2):
                                for z2 in range(-1,2):

                                    # now for EACH residue set the residue specific angle lookup for this angle. Angles in parameter file
                                    # are defined as RESIDUE A1 A2 A3, where A1,A2,A3 correspond to the 0,1,2 indexed positions in
                                    # int_to_penalty[intidx] 
                                    for intidx in int_list:
                                        
                                        # if straight line across central bead (A3 angle) - so i-1 and i+1 are 2 apart.
                                        # first three here define scenario where 2 of 3 dims are in plane
                                        #
                                        
                                        
                                        if ((abs(x2-x1) == 2 and abs(y2-y1) == 0 and abs(z2-z1) == 0) or
                                            (abs(y2-y1) == 2 and abs(x2-x1) == 0 and abs(z2-z1) == 0) or
                                            (abs(z2-z1) == 2 and abs(x2-x1) == 0 and abs(y2-y1) == 0) or
                                            (abs(z2-z1) == 2 and abs(x2-x1) == 2 and abs(y2-y1) == 0) or
                                            (abs(z2-z1) == 2 and abs(y2-y1) == 2 and abs(x2-x1) == 0) or
                                            (abs(y2-y1) == 2 and abs(x2-x1) == 2 and abs(z2-z1) == 0) or
                                            (abs(z2-z1) == 2 and abs(x2-x1) == 2 and abs(y2-y1) == 2)):
                                            penalty = int_to_penalty[intidx][2]
                                            AP3_count = AP3_count+1

                                        
                                        # to be uncommented after Switzerland!
                                        #if ((abs(z2-z1) == 2 and abs(x2-x1) == 2 and abs(y2-y1) == 0) or
                                        #    (abs(z2-z1) == 2 and abs(y2-y1) == 2 and abs(x2-x1) == 0) or
                                        #    (abs(y2-y1) == 2 and abs(x2-x1) == 2 and abs(z2-z1) == 0) or
                                        #    (abs(z2-z1) == 2 and abs(x2-x1) == 2 and abs(y2-y1) == 2)):
                                        #    penalty = int_to_penalty[intidx][2]
                                        #    AP3_count = AP3_count+1


                                        # if not straight line but beads are not adjacent (at least one of the
                                        # distances between i-1 and i+1 is 2
                                        # (A2 angle)
                                        elif (abs(x2-x1) == 2) or (abs(y2-y1) == 2) or (abs(z2-z1) == 2):
                                            penalty = int_to_penalty[intidx][1]
                                            AP2_count=AP2_count+1

                                        # else if adjacent 
                                        # (A1 angle)
                                        else:
                                            penalty = int_to_penalty[intidx][0]
                                            AP1_count = AP1_count+1
                                                                                        
                                        self.angle_lookup[intidx, x1+1,y1+1,z1+1,x2+1,y2+1,z2+1] = penalty
                                    
                                        
            #print "AP1: %i" % AP1_count
            #print "AP2: %i" % AP2_count
            #print "AP3: %i" % AP3_count
        else:
            self.angle_lookup = np.zeros((int_list[-1]+1, 3,3,3,3))
            for x1 in range(-1,2):
                for y1 in range(-1,2):                    
                    for x2 in range(-1,2):
                        for y2 in range(-1,2):


                            # now for EACH residue set the residue specific angle lookup for this angle. Angles in parameter file
                            # are defined as RESIDUE A1 A2 A3, where A1,A2,A3 correspond to the 0,1,2 indexed positions in
                            # int_to_penalty[intidx] 
                            for intidx in int_list:
                                        
                                # if straight line across central bead (A3 angle) - so i-1 and i+2
                                if ((abs(x2-x1) == 2 and abs(y2-y1) == 0) or 
                                (abs(y2-y1) == 2 and abs(x2-x1) == 0) or
                                (abs(x2-x1) == 2 and abs(y2-y1) == 2)):
                                    penalty = int_to_penalty[intidx][2]

                                # if not straight line but beads are not adjacent (at least one of the
                                # distances between i-1 and i+2 is 2
                                # (A2 angle)
                                elif abs(x2-x1) == 2 or abs(y2-y1) == 2:
                                    penalty = int_to_penalty[intidx][1]

                                # else if adjacent 
                                # (A1 angle)
                                else:
                                    penalty = int_to_penalty[intidx][0]
                            

                                self.angle_lookup[intidx, x1+1,y1+1,x2+1,y2+1] = penalty
                                    

    #-----------------------------------------------------------------
    #        
    def convert_sequence_to_integer_sequence(self, sequence):
        """
        This takes an human name residue sequence (e.g. amino acid) and converts
        it into it's local integer-code sequence which is then evaluated by the HYPERLOOP
        code - basically this happens because CYTHON can't do lookups with strings/chars but can with
        INTs, so we convert a list of chars into a list of ints and then go to TOWN on that badboy

        """
        int_seq = []
        for i in sequence:
            try:
                int_seq.append(self.parameter_to_int_map[i])
            except KeyError:
                raise ParameterFileException("Tried to convert residue [%s] into it's integer code, but no value in the parameter file was found!" %i)
                    
        return int_seq

    #-----------------------------------------------------------------
    #    
    def convert_sequence_to_LR_integer_sequence(self, sequence):
        """
        This takes an human name residue sequence (e.g. amino acid) and converts
        it into it's local LONG RANGE integer-code sequence. There is a unique mapping between
        each amino acid and it's integer LR residue code. This is dynamically constructed on
        each simulation run and should not be assumed to be the same between different
        simulations. 

        Also important, the LR_int code is not the same as the short range int_code. The LR_int
        code is defined by the self.LR_parameter_to_int_map, while the int code is defined
        by the self.parameter_to_int_map.

        Also importantly, EVERY residue must have a parameter_to_int_map entry but not every needs
        as LR_parameter_to_int_map entry. 

        Basically, this happens because Cython can't do lookups with strings/chars but can with
        INTs, so we convert a list of chars into a list of ints

        """
        
        int_seq = []
        for i in sequence:
            if i in self.LR_parameter_to_int_map:
                int_seq.append(self.LR_parameter_to_int_map[i])
            else:
                int_seq.append(-1)
                                    
        return int_seq
    
        
    def get_indices_of_long_range_residues(self, sequence):
        """

        Takes an amino acid sequence and returns a list with the index of positions
        which undergo long-range interactions.

        """
        
        LR_IDX = []
        idx=0
        
        for res in sequence:
            if res in self.human_readable_LR_interaction_table:
                LR_IDX.append(idx)

            idx=idx+1

        return LR_IDX
                    
        
                    

        


        

        
        
