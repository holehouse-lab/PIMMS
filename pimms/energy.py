## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................



import numpy as np

from . import lattice_utils
from . import lattice_analysis_utils
from .latticeExceptions import EnergyException, ParameterFileException
from . import parameterfile_parser
from . import hyperloop
from . import IO_utils

from . CONFIG import NP_INT_TYPE


class EmptyHamiltonian:
    """
    Dummy class which implements a Hamiltonian used when all interactions
    are turned off. Means we can define whatever stub-functionality here
    without adding code-rot to the true Hamiltonian class

    """

    def __init__(self):
        """
        Initialize an EmptyHamiltonian.

        Sets up the minimal attributes required for the object to act as a
        drop-in replacement for the true :class:`Hamiltonian` when all
        interactions are switched off. In particular, ``LR_residue_names`` is
        set to an empty list so that no residue is ever treated as
        long-range.

        Returns
        -------
        None

        """

        self.LR_residue_names = []

    def evaluate_total_energy(self, x):
        """
        Stub total-energy evaluation that always returns zero.

        Parameters
        ----------
        x : object
            Placeholder argument (typically a lattice object). Ignored.

        Returns
        -------
        float
            Always ``0.0``.

        """
        return 0.0

    def evaluate_local_energy(self, x, y):
        """
        Stub short-range local-energy evaluation that always returns zero.

        Parameters
        ----------
        x : object
            Placeholder argument (typically a lattice object). Ignored.
        y : object
            Placeholder argument (typically a list of pairs). Ignored.

        Returns
        -------
        float
            Always ``0.0``.

        """
        return 0.0

    def evaluate_local_energy_LR(self, x, y):
        """
        Stub long-range local-energy evaluation that always returns zero.

        Parameters
        ----------
        x : object
            Placeholder argument (typically a lattice object). Ignored.
        y : object
            Placeholder argument (typically a list of pairs). Ignored.

        Returns
        -------
        float
            Always ``0.0``.

        """
        return 0.0

    def evaluate_angle_energy(self, x,y):
        """
        Stub angle-energy evaluation that always returns zero.

        Parameters
        ----------
        x : object
            Placeholder argument (typically chain positions). Ignored.
        y : object
            Placeholder argument (typically an intcode sequence). Ignored.

        Returns
        -------
        float
            Always ``0.0``.

        """
        return 0.0

    def convert_sequence_to_integer_sequence(self, sequence):
        """
        Stub residue-to-integer conversion for the non-interacting case.

        Mirrors :meth:`Hamiltonian.convert_sequence_to_integer_sequence` but,
        because there are no interactions, simply maps every residue to the
        same dummy integer code (``1``).

        Parameters
        ----------
        sequence : list or str
            The (human readable) residue sequence to convert.

        Returns
        -------
        list of int
            A list of ``1`` values with the same length as ``sequence``.

        """
        return [1]*len(sequence)

    def convert_sequence_to_LR_integer_sequence(self, sequence):
        """
        Stub long-range residue-to-integer conversion for the non-interacting case.

        Mirrors :meth:`Hamiltonian.convert_sequence_to_LR_integer_sequence`
        but, because no residue undergoes long-range interactions, always
        returns an empty list.

        Parameters
        ----------
        sequence : list or str
            The (human readable) residue sequence to convert. Ignored.

        Returns
        -------
        list
            Always an empty list.

        """
        return []

    def get_indices_of_long_range_residues(self, sequence):
        """
        Stub lookup of long-range residue indices for the non-interacting case.

        Mirrors :meth:`Hamiltonian.get_indices_of_long_range_residues` but,
        because no residue undergoes long-range interactions, always returns
        an empty list.

        Parameters
        ----------
        sequence : list or str
            The (human readable) residue sequence to inspect. Ignored.

        Returns
        -------
        list
            Always an empty list.

        """
        return []


class Hamiltonian:
    """
    Main Hamiltonian class for evaluating system energies

    """

    def __init__(self, parameter_file, num_dimensions, non_interacting, angles_off, hardwall=False, temperature=False, reduced_printing=False):
        """


        Parameters
        -----------------
        parameter_file : str
            Location of a valid PIMMS parameter file. This is actually parsed
            via the parse_energy() function in the parameterfile_parser
            module, so all sanity checking happens there

        num_dimensions : int
            Number of dimensions the simulation is using (should be 2 or 3)

        non_interacting : bool
            Flag which, if set to true, means we over-ride any parameter
            file information and simply define the system setting all
            bead-bead interactions to zero. Note that bead excluded volume
            is still very much a thing.

        angles_off : bool
            Flag which, if set to true, means we over-ride angle parameters
            and turn all angles off.

        hardwall : bool
            Flag which defines if interactions engage via periodic boundary
            interactions or not. Default is False (i.e. PBC is in effect).

        reduced_printing : bool
            Flag which, if set to true, means we don't print any of the
            over-ride warning messages on startup. Can be useful when we
            want to supress input.

        
        Returns
        -------------
        Hamiltonian object
        
            

        # leave below for now..
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

        # set reduced printing
        if reduced_printing:
            self.reduced_printing = reduced_printing
        else:
            self.reduced_printing = reduced_printing
            

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
        """
        Toggle whether interactions use hardwall or periodic boundaries.

        Internally the hardwall flag is stored as an integer (``1`` for
        hardwall, ``0`` for periodic boundary conditions) because it is passed
        through to the Cython energy kernels.

        Parameters
        ----------
        value : bool
            If True the Hamiltonian uses hardwall boundary conditions
            (``self.hardwall = 1``); if False it uses periodic boundary
            conditions (``self.hardwall = 0``). Default is True.

        Returns
        -------
        None

        """
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

        Concretely, it walks every chain on the lattice, gathers the ordered
        bead positions and the per-bead long-range binary array, accumulates
        the angle energy chain-by-chain, builds the non-redundant short-range,
        long-range and super-long-range pair lists via
        ``lattice_utils.build_all_envelope_pairs`` and then evaluates each
        energy contribution.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object holding all chains, the type grid and the box
            dimensions over which the total energy is computed.

        id_to_typeMap : dict, optional
            Currently unused. Retained for interface compatibility. Default is
            None.

        Returns
        -------
        tuple of float
            A 5-tuple ``(total, energy_local, energy_LR, energy_SLR,
            angle_energy)`` where ``total`` is the sum of the short-range
            (local), long-range, super-long-range and angle energy
            contributions, and the remaining elements are those individual
            contributions.

        """
                
        all_positions = []
        lr_binary_chunks = []

        # cycle through each chain extracting the positions to generate a long list
        # associated with the positions of every residue on the lattice
        angle_energy = 0

        # for each chain extract out the positions enaging 
        for chainID in latticeObject.chains:

            all_positions.extend(latticeObject.chains[chainID].get_ordered_positions())
            lr_binary_chunks.append(latticeObject.chains[chainID].get_LR_binary_array())
            angle_energy = angle_energy + self.evaluate_angle_energy(latticeObject.chains[chainID].get_ordered_positions(), latticeObject.chains[chainID].get_intcode_sequence(), latticeObject.dimensions)

        if len(lr_binary_chunks) > 0:
            LR_binary_array = np.concatenate(lr_binary_chunks)
        else:
            LR_binary_array = np.array([], dtype=int)
                
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

        This particular method evaluates the SHORT-RANGE energy by dispatching
        to ``__evaluate_local_energy_shortrange`` using the short-range
        residue interaction table.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object holding the type grid and box dimensions.

        pairs_list : numpy.ndarray or list
            The (non-redundant) set of position pairs over which the
            short-range energy is evaluated.

        Returns
        -------
        float or int
            The total short-range interaction energy over the supplied pairs.
            Returns 0 when ``pairs_list`` is empty.

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

        This particular method evaluates the LONG-RANGE energy by dispatching
        to ``__evaluate_local_energy_non_shortrange`` using the long-range
        residue interaction table.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object holding the type grid and box dimensions.

        pairs_list : numpy.ndarray or list
            The (non-redundant) set of position pairs over which the
            long-range energy is evaluated.

        Returns
        -------
        float or int
            The total long-range interaction energy over the supplied pairs.
            Returns 0 when ``pairs_list`` is empty.

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

        This particular method evaluates the SUPER-LONG-RANGE (SLR) energy by
        dispatching to ``__evaluate_local_energy_non_shortrange`` using the
        super-long-range residue interaction table.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object holding the type grid and box dimensions.

        pairs_list : numpy.ndarray or list
            The (non-redundant) set of position pairs over which the
            super-long-range energy is evaluated.

        Returns
        -------
        float or int
            The total super-long-range interaction energy over the supplied
            pairs. Returns 0 when ``pairs_list`` is empty.

        """
        return self.__evaluate_local_energy_non_shortrange(latticeObject, pairs_list, self.SLR_residue_interaction_table)



    #-----------------------------------------------------------------
    #    
    #
    def __evaluate_local_energy_shortrange(self, latticeObject, pairs_list, interaction_table):
        """
        Internal general energy evaluation function for short-range interactions.

        Dispatches to the appropriate 2D or 3D Cython hyperloop routine
        (``evaluate_local_energy_2D_shortrange`` /
        ``evaluate_local_energy_3D_shortrange``) based on the lattice
        dimensionality.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object holding the type grid and box dimensions.

        pairs_list : numpy.ndarray or list
            The (non-redundant) set of position pairs over which the energy is
            evaluated.

        interaction_table : numpy.ndarray
            The integer-indexed short-range residue interaction table to use
            when scoring each pair.

        Returns
        -------
        float or int
            The short-range interaction energy over the supplied pairs.
            Returns 0 when ``pairs_list`` is empty.

        Raises
        ------
        EnergyException
            If the lattice dimensionality is neither 2 nor 3.

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

        raise EnergyException(f"Unsupported dimensionality for short-range energy evaluation: {num_dims}")


    #-----------------------------------------------------------------
    #    
    #
    def __evaluate_local_energy_non_shortrange(self, latticeObject, pairs_list, interaction_table):
        """
        Internal general energy evaluation function for long-range and super-long range interactions (i.e.
        not short range interactions).

        Dispatches to the appropriate 2D or 3D Cython hyperloop routine
        (``evaluate_local_energy_2D_non_shortrange`` /
        ``evaluate_local_energy_3D_non_shortrange``) based on the lattice
        dimensionality. The same routine is used for both long-range (LR) and
        super-long-range (SLR) interactions; the distinction is made entirely
        by which ``interaction_table`` is passed in.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object holding the type grid and box dimensions.

        pairs_list : numpy.ndarray or list
            The (non-redundant) set of position pairs over which the energy is
            evaluated.

        interaction_table : numpy.ndarray
            The integer-indexed (long-range or super-long-range) residue
            interaction table to use when scoring each pair.

        Returns
        -------
        float or int
            The non-short-range interaction energy over the supplied pairs.
            Returns 0 when ``pairs_list`` is empty.

        Raises
        ------
        EnergyException
            If the lattice dimensionality is neither 2 nor 3.

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

        raise EnergyException(f"Unsupported dimensionality for non-short-range energy evaluation: {num_dims}")

        
 
    #-----------------------------------------------------------------
    #    
    def evaluate_angle_energy(self, chain_positions, intcode_sequence, dimensions):
        """
        Angle energies are determined by hyperloop functions that basically
        compare the angle vector and use a pre-computed lookup table to convert
        that angle into some energy penalty. It's lightning fast, and residue
        specific!

        Chains with fewer than 3 beads have no defined angle and contribute
        zero. For valid chains the work is dispatched to the 2D or 3D Cython
        angle hyperloop (``evaluate_angle_energy_2D`` /
        ``evaluate_angle_energy_3D``) using the pre-computed
        ``self.angle_lookup`` table.

        Parameters
        ----------
        chain_positions : list
            A list of 2D or 3D bead positions (each a 2- or 3-element list of
            integer lattice coordinates) over which the angle energy is
            evaluated.

        intcode_sequence : list of int
            An equal-length list of the per-bead integer codes used to index
            the residue-specific angle lookup table.

        dimensions : list
            A list of the box dimensions; its length (2 or 3) determines the
            lattice dimensionality.

        Returns
        -------
        float or int
            The total angle penalty for the supplied chain. Returns ``0.0``
            when the chain has fewer than 3 beads.

        Raises
        ------
        EnergyException
            If the lattice dimensionality is neither 2 nor 3.

        """

        num_positions = len(chain_positions)

        # cannot compute an angle for a chain with only 1 or 2 beads
        if num_positions < 3:
            return 0.0
        
        num_dims = len(dimensions)
        if num_dims not in (2, 3):
            raise EnergyException(f"Unsupported dimensionality for angle energy evaluation: {num_dims}")

        # for the 3D case
        if num_dims == 3:
            penalty = hyperloop.evaluate_angle_energy_3D(np.array(chain_positions, dtype=NP_INT_TYPE), np.array(intcode_sequence, dtype=NP_INT_TYPE), self.angle_lookup, num_positions)

        # for the 2D case
        else:
            penalty = hyperloop.evaluate_angle_energy_2D(np.array(chain_positions, dtype=NP_INT_TYPE), np.array(intcode_sequence, dtype=NP_INT_TYPE), self.angle_lookup, num_positions)

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

        Parameters
        ----------
        non_interacting : bool
            If True, every entry in the short-range, long-range and
            super-long-range tables is forced to zero (overriding the
            parameter file), and a warning is emitted per residue pair unless
            reduced printing is enabled. Default is False.

        Returns
        -------
        tuple
            A 5-tuple ``(RIT, MAPPING, LRRIT, LR_MAPPING, SLRRIT)`` as
            described above: the short-range interaction table, the
            residue-name-to-integer mapping, the long-range interaction table,
            the long-range-only residue-name-to-integer mapping, and the
            super-long-range interaction table.

        """

        # number of different residue types we're messing with
        n_residues    = len(self.residue_names)

        # initialize the residue interaction table [RIT] as a matrix
        # of zeros
        RIT    = np.zeros(shape=(n_residues, n_residues), dtype=NP_INT_TYPE)
        LRRIT  = np.zeros(shape=(n_residues, n_residues), dtype=NP_INT_TYPE)
        SLRRIT = np.zeros(shape=(n_residues, n_residues), dtype=NP_INT_TYPE)

        MAPPING = {}
        LR_MAPPING = {}
        
        # Build int-index interaction matrices for short-range and long-range
        # residue interactions 
        # for each residue (human name)
        ## WARNING:
        ## THIS WORKS BECAUSE WE ASSUME THE FIRST RESIDUE IN
        ## RESIDUE NAMES IS SOLVENT - this means 0 is always
        ## the solvent!
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
                    if self.reduced_printing is False:
                        IO_utils.status_message("This is a non-interacting run (over-riding parameter file for [%s-%s])" % (R1,R2), 'warning')
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
        Function that constructs a lookup table that we use to assign 'angle'
        withstraints. The angle effects are really related to the 1_3
        interaction, but can also be used

        The resulting integer-typed lookup table is stored on the instance as
        ``self.angle_lookup``. Its shape depends on dimensionality: in 3D it is
        ``(max_intcode + 1, 3, 3, 3, 3, 3, 3)`` and in 2D it is
        ``(max_intcode + 1, 3, 3, 3, 3)``. The first axis is indexed by the
        residue integer code and the remaining axes encode the relative
        offsets of the i-1 and i+1 neighbours. Each residue must have an
        associated angle definition (three penalties A1/A2/A3); a missing
        definition raises an exception. Non-integer penalties are rounded to
        the nearest integer because the lattice energies are integer valued.
        When angles are switched off (or when only solvent exists) the table is
        populated with zeros.

        Parameters
        ----------
        angle_dict : dict or bool
            Mapping from residue name to its ``[A1, A2, A3]`` angle penalties,
            as parsed from the parameter file. When ``angles_off`` is True this
            argument is ignored (and is typically passed as ``False``).

        num_dimensions : int
            The lattice dimensionality. Must be 2 or 3.

        angles_off : bool
            If True, all angle penalties are forced to zero and the residue
            names are taken from the interaction table rather than from
            ``angle_dict``.

        Returns
        -------
        None
            The function does not return a value; it sets ``self.angle_lookup``
            as a side effect.

        Raises
        ------
        EnergyException
            If ``num_dimensions`` is neither 2 nor 3.

        ParameterFileException
            If a residue has through-space interactions defined but no
            corresponding angle energies.

        """

        if num_dimensions not in (2, 3):
            raise EnergyException(f"Unsupported dimensionality for angle interactions: {num_dimensions}")

        #
        # A key first thing we have to to do is  verify that every residue in the
        # interaction table has an angle. If we have extra angles that are lack
        # pairwise bead interactions these are just ignored

        # first get a list of the bead names where we have intermolecular
        # interactions defined. We SHOULD have angle parameters for each
        # of these!
        resnames_through_space = list(self.parameter_to_int_map.keys())

        # if angles are off, then we get the 'angle' names through the interactions
        # and set all to 0, i.e. we don't need to make sure 
        if angles_off:
            resnames_angles = list(self.parameter_to_int_map.keys())

        # otherwise get the angle names based on the angle information we parsed
        # from the parameter file
        else:
            resnames_angles       = list(angle_dict.keys())


        
        int_to_penalty = {}        
        for resname in resnames_through_space:

            # skip solvent!
            if resname == '0':
                pass

            else:

                # if we lack angle information for this residue throw an exception
                if resname not in resnames_angles:
                    raise ParameterFileException(f'Residue {resname} has interaction energies defined but *no angle energies defined*')

                # if angles are off set all penalties to 0
                if angles_off:
                    if self.reduced_printing == False:
                        IO_utils.status_message(f"Angles are turned off (over-riding parameter file for {resname} angles)", 'warning')
                    int_to_penalty[self.parameter_to_int_map[resname]] = [0,0,0]
                else:
                    int_to_penalty[self.parameter_to_int_map[resname]] = angle_dict[resname]

                
        # The angle-energy pipeline (the Cython hyperloop) uses an integer-typed
        # lookup table (angle_lookup is NP_INT_TYPE), so penalties must be integers.
        # T_NORM angle penalties are floats (penalty * temperature); round them to
        # the nearest integer here rather than letting the later int32 array
        # assignment silently truncate toward zero (which systematically
        # under-weights every fractional penalty).
        for intkey in int_to_penalty:
            rounded = [int(round(p)) for p in int_to_penalty[intkey]]
            if rounded != list(int_to_penalty[intkey]) and self.reduced_printing == False:
                IO_utils.status_message(
                    "Non-integer angle penalties %s rounded to %s (lattice energies are integer-valued)"
                    % (list(int_to_penalty[intkey]), rounded), 'warning')
            int_to_penalty[intkey] = rounded

        # int_list is a sorted list of the intergers that map to a residue-specific angle pair
        int_list = list(int_to_penalty.keys())
        int_list.sort()

        # If all residues are solvent or no angle-enabled residues exist, create an
        # empty lookup that is still shape-consistent for downstream code paths.
        if len(int_list) == 0:
            if num_dimensions == 3:
                self.angle_lookup = np.zeros((1, 3, 3, 3, 3, 3, 3), dtype=NP_INT_TYPE)
            else:
                self.angle_lookup = np.zeros((1, 3, 3, 3, 3), dtype=NP_INT_TYPE)
            return
        
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
            ## while the actual int associated with the intidx reflects the identity of residue i


            self.angle_lookup = np.zeros((int_list[-1]+1, 3, 3, 3, 3, 3, 3), dtype=NP_INT_TYPE)

            AP1_count = 0
            AP2_count = 0
            AP3_count = 0
        
            for x1 in range(-1,2):
                for y1 in range(-1,2):
                    for z1 in range(-1,2):
                        for x2 in range(-1,2):
                            for y2 in range(-1,2):
                                for z2 in range(-1,2):

                                    # now for EACH residue set the residue specific angle lookup for this angle.
                                    # Angles in parameter file are defined as RESIDUE A1 A2 A3, where
                                    # A1, A2, A3 correspond to the 0,1,2 indexed positions in int_to_penalty[intidx] 
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

                                        # define penalty lookup
                                        self.angle_lookup[intidx, x1+1,y1+1,z1+1,x2+1,y2+1,z2+1] = penalty
                                                                            
        else:
            self.angle_lookup = np.zeros((int_list[-1]+1, 3,3,3,3), dtype=NP_INT_TYPE)
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

        Parameters
        ----------
        sequence : list or str
            The (human readable) residue sequence to convert. Each element must
            be a key in ``self.parameter_to_int_map``.

        Returns
        -------
        list of int
            The integer-code sequence, one integer per residue in
            ``sequence``.

        Raises
        ------
        ParameterFileException
            If a residue in ``sequence`` has no corresponding integer code in
            the parameter file mapping.

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

        Parameters
        ----------
        sequence : list or str
            The (human readable) residue sequence to convert.

        Returns
        -------
        list of int
            The long-range integer-code sequence, one entry per residue. A
            value of ``-1`` is used for any residue that does not participate
            in long-range interactions (i.e. is absent from
            ``self.LR_parameter_to_int_map``).

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
        Return the sequence indices of residues that undergo long-range interactions.

        Takes an amino acid sequence and returns a list with the index of positions
        which undergo long-range interactions. A residue is considered
        long-range if it appears as a key in
        ``self.human_readable_LR_interaction_table``.

        Parameters
        ----------
        sequence : list or str
            The (human readable) residue sequence to inspect.

        Returns
        -------
        list of int
            The (zero-based) indices into ``sequence`` of residues that engage
            in long-range interactions.

        """

        LR_IDX = []
        idx=0
        
        for res in sequence:
            if res in self.human_readable_LR_interaction_table:
                LR_IDX.append(idx)

            idx=idx+1

        return LR_IDX


        


        

        
        
