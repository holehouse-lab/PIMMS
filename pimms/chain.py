## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import numbers

import numpy as np

from . import analysis_structures
from . import lattice_utils
from . import lattice_analysis_utils
from .latticeExceptions import ChainInsertionFailure, ChainInitializationException, ChainAugmentFailure
from .CONFIG import *

class Chain:
    """A single polymer chain: its sequence, its ordered bead positions on the
    lattice, and its per-chain analyses.

    A Chain records the one-letter sequence, the integer interaction codes
    (SR and LR), the ordered position list, and its boundary convention
    (``hardwall``). It provides the conformational observables PIMMS reports -
    radius of gyration, asphericity, end-to-end distance, distance maps and
    internal scaling. Under periodic boundaries every one of these is computed
    on the chain *made whole* (bond-walked into a single periodic image, see
    :meth:`get_analysis_positions`); under a hardwall the raw positions are
    already a single image and are used directly.
    """

    def __init__(self,
                 lattice_grid,
                 dimensions,
                 sequence,
                 int_seq,
                 LR_int_seq,
                 LR_IDX,
                 chainID,
                 chainType,
                 chain_positions=None,
                 rigid=False,
                 center=False,
                 hardwall=False):
        """
        Constructor for a Chain object. In PIMMS, every individual polymer is defined as a Chain, so 
        there will be many, many Chains per simulation

        Parameters
        ----------------
        lattice_grid : np.ndarray
            Integer array of shape ``dimensions`` (2D or 3D) holding the occupancy grid
            the chain is going to be inserted into (with other chains in it!) - used
            to find a vacant position for the chain

        dimensions : list of int
            Box size in 2 or 3 dimensions (i.e. a length-2 or length-3 list of positive
            integers). Must match the shape of lattice_grid, and defines both the lattice
            bounds used when validating positions and the box used by every periodic
            calculation this chain performs.

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

        chain_positions : list of positions, optional
            Default None. If provided, the chain is initialized directly at these positions rather than
            being inserted into the lattice by a self-avoiding walk. Note that "positions"
            here are 2D or 3D sublists that define X/Y or X/Y/Z coordinates. There must be
            exactly one position per residue, every position must lie inside the lattice,
            and no two beads may share a coordinate.

        rigid : bool, optional
            If the chain is limited to rigid body movements only. Default False. This is not
            yet implemented but will be soon.

        center : bool, optional
            Defines if we're going to try and place the chain in the center of the box/square.
            Only meaningful when chain_positions is None. Default False.

        hardwall : bool, optional
            Flag which defines if the simulation is using periodic boundary conditions (PBC) or
            hardwall boundary conventions. PIMMS by default uses PBC, so this defaults to False.
            Stored on the chain and used by every coordinate-derived observable to decide
            whether periodic corrections apply.

        Raises
        ------
        ChainInitializationException
            If any of the chain-defining arguments are invalid - a non-positive or
            non-integer chainID, a negative or non-integer chainType, an empty or
            non-string sequence, dimensions that are not 2 or 3 positive integers, a
            lattice_grid whose dtype or shape does not match, interaction code lists of
            the wrong length or type, out-of-range long-range indices, or a
            chain_positions list of the wrong length or containing malformed,
            out-of-bounds or duplicated coordinates.

        ChainInsertionFailure
            If chain_positions is not provided and a self-avoiding walk of the required
            length could not be placed on the lattice (typically an overcrowded lattice).

        """
        
        int_limit = np.iinfo(NP_INT_TYPE)
        if (isinstance(chainID, (bool, np.bool_)) or
                not isinstance(chainID, numbers.Integral) or
                chainID <= 0 or chainID > int_limit.max):
            raise ChainInitializationException(
                "chainID must be a positive integer representable by the lattice grid")
        if (isinstance(chainType, (bool, np.bool_)) or
                not isinstance(chainType, numbers.Integral) or chainType < 0 or
                chainType > int_limit.max):
            raise ChainInitializationException(
                "chainType must be a non-negative integer representable by int32")
        if not isinstance(sequence, str) or len(sequence) == 0:
            raise ChainInitializationException("Chain sequence must be a non-empty string")
        try:
            normalized_dimensions = tuple(dimensions)
        except TypeError:
            raise ChainInitializationException(
                "Chain dimensions must contain 2 or 3 positive integers")
        if (len(normalized_dimensions) not in (2, 3) or
                any(isinstance(value, (bool, np.bool_)) or
                    not isinstance(value, numbers.Integral) or value <= 0
                    for value in normalized_dimensions)):
            raise ChainInitializationException(
                "Chain dimensions must contain 2 or 3 positive integers")
        normalized_dimensions = [int(value) for value in normalized_dimensions]
        if (not isinstance(lattice_grid, np.ndarray) or
                not np.issubdtype(lattice_grid.dtype, np.integer) or
                lattice_grid.shape != tuple(normalized_dimensions)):
            raise ChainInitializationException(
                "Chain lattice_grid must be an integer array matching dimensions")

        def _normalize_codes(values, label):
            """
            Validate a per-residue integer code list and return it as a list of ints.

            Parameters
            ----------
            values : iterable of int
                One integer interaction code per residue, so the length must match
                the chain sequence. Values must be representable by NP_INT_TYPE
                because they are written into the lattice type grid.

            label : str
                Name of the argument being checked ('int_seq' or 'LR_int_seq'),
                used to build the exception message.

            Returns
            -------
            list of int
                The codes converted to plain Python ints.

            Raises
            ------
            ChainInitializationException
                If values is not iterable, does not have one entry per residue, or
                contains a non-integer (booleans included) or out-of-range code.
            """
            try:
                values = list(values)
            except TypeError:
                raise ChainInitializationException(
                    f"{label} must contain one integer code per residue")
            if len(values) != len(sequence):
                raise ChainInitializationException(
                    f"{label} must contain one integer code per residue")
            if any(isinstance(value, (bool, np.bool_)) or
                   not isinstance(value, numbers.Integral) or
                   value < int_limit.min or value > int_limit.max for value in values):
                raise ChainInitializationException(
                    f"{label} codes must be integers representable by int32")
            return [int(value) for value in values]

        int_seq = _normalize_codes(int_seq, "int_seq")
        LR_int_seq = _normalize_codes(LR_int_seq, "LR_int_seq")
        try:
            LR_IDX = list(LR_IDX)
        except TypeError:
            raise ChainInitializationException("LR_IDX must be a sequence of indices")
        if any(isinstance(idx, (bool, np.bool_)) or
               not isinstance(idx, numbers.Integral) or
               idx < 0 or idx >= len(sequence) for idx in LR_IDX):
            raise ChainInitializationException(
                f"Long-range index list is invalid for chainID {chainID} with "
                f"sequence length {len(sequence)}")
        LR_IDX = [int(idx) for idx in LR_IDX]

        # set the chain ID
        self.chainID      = int(chainID)

        # set the chain type ID
        self.chainType    = int(chainType)

        # set the lattice dimensions
        self.dimensions   = normalized_dimensions

        # Boundary convention used by every coordinate-derived observable. A
        # hardwall chain lives in one ordinary Cartesian box, so minimum-image
        # corrections would turn physically distant beads near opposite walls
        # into apparent neighbours.
        self.hardwall     = bool(hardwall)

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

        # set if the chain is limited to rigid-body movements or not
        self.rigid = rigid

        # define the index of residues which undergo LR interactions
        self.LR_IDX = LR_IDX

        # per-bead 0/1 flag for long-range participation. Fixed for the life of the
        # chain, so build it once here (see get_LR_binary_array).
        self._LR_binary_array = np.zeros(len(sequence), dtype=NP_INT_TYPE)
        for idx in LR_IDX:
            self._LR_binary_array[idx] = 1
        self._LR_binary_array.flags.writeable = False

        # automatically determine if sequence is a homopolymer
        if len(set(sequence)) == 1:
            self.homopolymer = True
        else:
            self.homopolymer = False

        # if we passed chain positions in....
        if chain_positions is not None:
            # if we're intializing a chain when we aready know it's positions on
            # the grid
            
            # check to make sure we're not trying to initialize the wrong number of
            # positions
            if len(chain_positions) != self.seq_len:
                raise ChainInitializationException('Tried to initialize a chain [chainID = %i] with a sequence of length %i but had %i positions' % (chainID, self.seq_len, len(chain_positions)))

            normalized_positions = []
            occupied = set()
            for position in chain_positions:
                try:
                    position = list(position)
                except TypeError:
                    raise ChainInitializationException(
                        f"Malformed position in chainID {chainID}")
                if (len(position) != len(self.dimensions) or
                        any(isinstance(value, (bool, np.bool_)) or
                            not isinstance(value, numbers.Integral)
                            for value in position)):
                    raise ChainInitializationException(
                        f"Malformed position {position!r} in chainID {chainID}")
                position = [int(value) for value in position]
                if any(value < 0 or value >= self.dimensions[axis]
                       for axis, value in enumerate(position)):
                    raise ChainInitializationException(
                        f"Position {position!r} in chainID {chainID} is outside the lattice")
                coordinate = tuple(position)
                if coordinate in occupied:
                    raise ChainInitializationException(
                        f"ChainID {chainID} contains two beads at {coordinate}")
                occupied.add(coordinate)
                normalized_positions.append(position)
            self.positions = normalized_positions

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
    def __len__(self):
        """
        Return the number of beads (positions) in the chain.

        Returns
        -------
        int
            The number of lattice positions currently held by the chain, i.e.
            the chain length.

        """
        return len(self.positions)
    

    #-----------------------------------------------------------------
    #
    def get_ordered_positions(self, center_positions=False):
        """
        Returns a list of the chain positions.

        If center_positions is True, the positions are centered in the center 
        of the simulation box and any periodic boundary stuff is fixed, so you 
        get a single image chain. This is useful for visualization purposes.

        Parameters
        ----------

        center_positions : bool, optional
            If True, the positions are made into a single image and then centered in the
            simulation box. Default False.

        Returns
        -------
        list
            A list of the chain positions. If center_positions is True, the positions
            are centered in the center of the simulation box and any periodic boundary
            stuff is fixed, so you get a single image chain. This is useful for
            visualization purposes.

        """

        if center_positions is True:
            return lattice_utils.center_positions(self.get_single_image_positions(), self.dimensions)
        else:            
            return self.positions

    #-----------------------------------------------------------------
    #
    def get_intcode_sequence(self):
        """
        Returns a list where each position corresponds to the integer code
        used by the energy calculations to identify a specific residue type

        Returns
        -------
        list of int
            A list where each position corresponds to the integer code
            used by the energy calculations to identify a specific residue
            type.

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

        Returns
        -------
        list
            A list of the chain positions where all positions in the chain are
            corrected to lie in the same periodic image (i.e. "single image
            convention" as opposed to minimum image convention.)
           
        """

        # if the chain does not straddle a periodic boundary
        if not self.does_chain_stradle_pbc_boundary():
            return self.positions
        else:
            return lattice_utils.convert_chain_to_single_image(self.positions, self.dimensions)


    #-----------------------------------------------------------------
    #
    def get_analysis_positions(self):
        """
        Return the positions every intra-chain observable is computed from.

        Under periodic boundaries a chain that crosses a box face is stored with
        its beads wrapped back into the box, so its raw positions are not a
        contiguous object. Every intra-chain observable (radius of gyration,
        asphericity, end-to-end distance, internal scaling, distance maps and
        residue-residue distances) is therefore computed on the chain **made
        whole**: the beads are bond-walked from the first bead into a single
        periodic image (:func:`~pimms.lattice_utils.make_chain_whole`), after which
        plain Cartesian geometry is exact for any chain that does not percolate the
        box. A chain that does not cross a face is already whole and is returned
        unchanged; so is every hardwall chain.

        This replaces two conventions that were both wrong once a chain spans
        more than half the box along any axis. Selecting each bead's image
        relative to the centre of mass tore such chains into two pieces (a bead
        further than half a box from the COM was shifted by a full box length
        even though it was bonded to its neighbour), which under-reported the
        radius of gyration by up to ~40 % for chains that never crossed a
        boundary at all. Minimum-image pair distances silently picked the
        nearer periodic image of the partner bead for any pair separated by more
        than half a box, so a nearly straight chain of length 0.8 L reported an
        end-to-end distance of 0.2 L. Neither effect required a boundary crossing
        and neither was flagged.

        Returns
        -------
        list
            The bead positions in N->C order, contiguous in a single periodic
            image (coordinates may fall outside the box on either face).
        """
        if self.hardwall or not self.does_chain_stradle_pbc_boundary():
            return self.positions
        return lattice_utils.make_chain_whole(self.positions, self.dimensions)


    #-----------------------------------------------------------------
    #
    def get_output_positions(self, autocenter=False, unwrap=False):
        """
        Return the chain positions to write to a trajectory frame / PDB.

        Selects between three visualisation conventions:

        * ``autocenter`` (single-chain only) - single-image positions centred in
          the box (used by AUTOCENTER),
        * ``unwrap`` - "whole" positions anchored at the first bead, i.e. the chain
          is bond-walked so it is contiguous across periodic boundaries (coordinates
          may fall outside the box on either face; used by TRAJECTORY_PBC_UNWRAP),
        * neither - the raw on-lattice positions (the default).

        ``autocenter`` takes precedence over ``unwrap`` (it already makes the chain
        whole before centring). The bead ordering is identical in all three cases,
        so the trajectory stays consistent with the topology.

        Parameters
        ----------
        autocenter : bool, optional
            If True, return single-image positions centred in the box. Default False.

        unwrap : bool, optional
            If True (and ``autocenter`` is False), return the chain made whole by
            bond-walking it into a single periodic image. Ignored when
            ``autocenter`` is True. Default False.

        Returns
        -------
        list
            The chain's bead positions in N->C order under the selected convention.
        """
        if autocenter:
            return self.get_ordered_positions(center_positions=True)
        if unwrap:
            # make the chain whole in place (anchored at bead 0); only pay the
            # unwrap cost for chains that actually cross a boundary
            if self.does_chain_stradle_pbc_boundary():
                return lattice_utils.make_chain_whole(self.positions, self.dimensions)
            return self.positions
        return self.get_ordered_positions()


    #-----------------------------------------------------------------
    #
    def does_chain_stradle_pbc_boundary(self):
        """
        Determines if the chain straddles a periodic boundary or not.

        Returns
        -------
        bool
            True if the chain straddles a periodic boundary, False otherwise.

        """

        return lattice_utils.do_positions_stradle_pbc_boundary(self.positions)
        
      
    #-----------------------------------------------------------------
    #
    def get_LR_positions(self):
        """
        Returns a list of chain positions which engage in long range interactions

        Returns
        -------
        list
            A list of chain positions which engage in long range interactions, i.e.
            list where each element is a 2 or 3 element list of the x, y [and z]
            coordinates of the chain positions which engage in long range interactions.


        """
        return [self.positions[i] for i in self.LR_IDX]



    #-----------------------------------------------------------------
    #
    def get_LR_binary_array(self):
        """
        Returns a numpy array of chain length, where beads that engage in long
        range interactions are set to 1 and all others are set to 0

        Returns
        -------
        numpy.ndarray
            Read-only array of shape ``(seq_len,)`` and dtype NP_INT_TYPE, where
            beads that engage in long range interactions are set to 1 and all
            others are set to 0.

        Notes
        -----
        ``LR_IDX`` is fixed for the lifetime of a chain, so this array is built once in
        the constructor and handed back on every call. It used to be rebuilt per call
        with ``if i in self.LR_IDX`` against a *list*, i.e. O(L^2) per call - and it is
        called once per chain in every full energy evaluation and on every single-chain
        move. The array is marked read-only so a caller cannot corrupt the shared copy.
        """
        return self._LR_binary_array



    #-----------------------------------------------------------------
    #
    def get_positions_by_chain_index(self, index_list):
        """
        Returns a list of lattice positions based on the index positions
        in the index_list.

        Parameters
        ----------
        index_list : list of int
            Sequence indices (0 to seq_len-1) of the beads to return.

        Returns
        -------
        list
            The raw on-lattice positions of the requested beads, in the order the
            indices were given. Each position is a 2 or 3 element list of x, y [and z]
            coordinates.

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
        a PBC boundary.

        Parameters
        ----------
        index_list : list of int
            Sequence indices (0 to seq_len-1) of the beads to return.

        Returns
        -------
        list
            The single-image positions of the requested beads, in the order the
            indices were given. Each position is a 2 or 3 element list of x, y [and z]
            coordinates, and coordinates may lie outside the box.

        """
        SIP = lattice_utils.convert_chain_to_single_image(self.positions, self.dimensions)

        position_list = []
        for i in index_list:
            position_list.append(SIP[i])

        return position_list



    #-----------------------------------------------------------------
    #                    
    def set_ordered_positions(self, positions):
        """
        Sets the chain positions to the given list of positions. This enables an
        entire set of positions in a chain to be updated.

        Parameters
        ----------
        positions : list
            The full replacement set of lattice positions in N->C order, one per
            residue. Note lattice positions are 2 or 3 element lists of the x, y [and z]
            coordinates of the chain positions.

        Returns
        -------
        None

        Raises
        ------
        ChainAugmentFailure
            If the number of supplied positions does not match the chain's
            sequence length (``self.seq_len``).

        """

        if len(positions) == self.seq_len:
            self.positions = positions
        else:
            raise ChainAugmentFailure(
                f'Tried to set chainID {self.chainID} to a set of positions of length {len(positions)}, '
                f'but this chain requires {self.seq_len} positions'
            )



    #-----------------------------------------------------------------
    #
    def get_center_of_mass(self, on_lattice=True):
        """
        Returns a chain's center of mass (note for now we assume every bead 
        has the same mass such that the center of mass ends up becoming the 
        mean position in all 2 or 3 dimensions (depending on the system).

        Parameters
        ----------
        on_lattice : bool, optional
            If True (the default) then the center of mass is returned
            as a lattice position, if False then the center of mass is returned
            as a continous space position.

        Returns
        -------
        list
            A list of the x, y [and z] coordinates of the center of mass of the 
            chain. If on_lattice is True then the center of mass is returned
            as a lattice position (i.e. integer x/y[/z] positions), if False then 
            the center of mass is return as a continous space position (i.e. float
            x/y[/z] positions).        

        """
        
        return lattice_utils.center_of_mass_from_positions(
            self.get_ordered_positions(),
            self.dimensions,
            on_lattice=on_lattice,
        )


    #####################################################################################################
    ## INTERNAL SCALING ANALYSIS FUNCTIONS
    ##

    def analysis_get_instantaneous_internal_scaling(self, mode='dict'):
        """
        Returns the instantaneous internal scaling profile for the chain.

        If mode is set to dict then a dictionary is returned where keys are the gaps
        and values are the average inter-residue distances.

        If mode is set to array then a 2-seq_len np.ndarray is returned with gaps in 
        column 1 and distances in column 2.

        Parameters
        --------------
        mode : str, optional
            Selector that determines the return type, either 'dict' or 'array'.
            Default 'dict'.


        Returns
        ----------
        dict or np.ndarray
            Sequence separation vs. mean spatial separation for the chain. In
            'dict' mode a dictionary keyed by sequence separation; in 'array'
            mode a ``(2, seq_len-1)`` np.ndarray with separations in row 0 and
            distances in row 1.

        Raises
        ------
        Exception
            If ``mode`` is neither ``'dict'`` nor ``'array'``.
        """

        if mode not in ('dict', 'array'):
            # should make this a better exception
            raise Exception('Invalid mode provided')

        # One vectorized pass per sequence separation, rather than a Python loop over
        # every pair calling the scalar distance helper. Bit-identical, but this used
        # to be one of the two dominant costs of an analysis step (it is O(L^2) pairs
        # per chain per call). Computed on the chain made whole, so no minimum-image
        # step is needed (or wanted - see get_analysis_positions).
        (ij_gaps, ij_vals) = lattice_analysis_utils.get_internal_scaling_profile(
            self.get_analysis_positions(), self.dimensions, pbc_correction=False)

        if mode == 'array':
            return np.array([ij_gaps, ij_vals])

        return dict(zip(ij_gaps, ij_vals))

        

    def analysis_update_internal_scaling(self):        
        """
        Function which when called will re-evaluate the chain's current internal scaling 
        information based on its current position and then update the local 
        internal_scaling object to include this most recent analysis. Note this does 
        not REPLACE the current internal scaling information, but allows a running average 
        which should become more accurate the more frequently the function is called.

        Returns
        ----------
        None
        """

        # Compute the first and second pair-distance moments in the SAME vectorized
        # pass. The squared profile must be mean(r_ij**2), not
        # mean(r_ij)**2; those differ whenever distances at a sequence gap are not all
        # identical within a snapshot.
        ij_gaps, ij_vals, ij_squared = lattice_analysis_utils.get_internal_scaling_profile(
            self.get_analysis_positions(), self.dimensions, pbc_correction=False,
            return_squared=True)

        self.internal_scaling.update_internal_scaling(dict(zip(ij_gaps, ij_vals)))
        self.internal_scaling_squared.update_internal_scaling_squared(
            dict(zip(ij_gaps, ij_squared)))


    #-----------------------------------------------------------------
    #
    def analysis_print_internal_scaling(self):
        """
        Prints the current internal scaling profile for the chain.

        Delegates to ``print_status`` on the chain's ``internal_scaling``
        analysis object.

        Returns
        ----------
        None

        """
        self.internal_scaling.print_status()


    #-----------------------------------------------------------------
    #
    def analysis_get_cumulative_internal_scaling(self):
        """
        Returns the cumulative internal scaling profile for the chain.

        Returns
        ----------
        list of float
            The ensemble-averaged inter-residue distances ordered by increasing
            sequence separation (index 0 is a gap of 1 residue).

        """
        return self.internal_scaling.get_internal_scaling_array()        


    #-----------------------------------------------------------------
    #
    def analysis_print_internal_scaling_squared(self):
        """
        Prints the squared current internal scaling profile for the chain.

        Delegates to ``print_status`` on the chain's
        ``internal_scaling_squared`` analysis object.

        Returns
        -------
        None

        """
        self.internal_scaling_squared.print_status()


    #-----------------------------------------------------------------
    #
    def analysis_get_internal_scaling_squared(self):
        """
        Returns the cumulative squared internal scaling profile for the chain.

        Delegates to ``get_internal_scaling_array`` on the chain's
        ``internal_scaling_squared`` analysis object, giving the
        ensemble-averaged squared inter-residue distances as a function of
        sequence separation.

        Returns
        -------
        list of float
            The ensemble-averaged squared inter-residue distances ordered by
            increasing sequence separation (index 0 is a gap of 1 residue).

        """
        return self.internal_scaling_squared.get_internal_scaling_array()



    #-----------------------------------------------------------------
    #
    def analysis_fit_scaling_exponent(self):
        """
        Fits and returns the polymer scaling exponent for the chain.

        Delegates to ``fit_scaling_exponent`` on the chain's
        ``internal_scaling_squared`` analysis object, which fits the
        (squared) internal scaling profile to extract the apparent scaling
        exponent.

        Returns
        -------
        tuple of float
            ``(nu, R0)``, the fitted scaling exponent and prefactor. Returns
            ``(-1, -1)`` when the chain is too short (fewer than 25 sequence
            separations) for a fit to be meaningful.

        """
        return self.internal_scaling_squared.fit_scaling_exponent()


    #####################################################################################################
    ## DISTANCE MAP ANALYSIS FUNCTIONS
    ##

    def analysis_get_instantaneous_distance_map(self):
        """
        Function that return's the chains instantaneous distance map.

        Returns
        ----------
        np.ndarray (seq_len, seq_len) of floats
            Returns the full symmetric inter-residue distance matrix.

        """
        # A distance map is a square symmetric observable. The historical loop only
        # populated its upper triangle and left the lower triangle at zero, which made
        # half the residue pairs look coincident when the output was loaded or shown
        # directly with imshow.
        return lattice_analysis_utils.get_distance_matrix(
            self.get_analysis_positions(), self.dimensions, pbc_correction=False)


    def analysis_update_distance_map(self):
        """
        Re-evaluate the instantaneous distance map and fold it into the running average.

        Computes the chain's current instantaneous inter-residue distance map
        and passes it to the cumulative ``distance_map`` analysis object, which
        maintains a running average over the ensemble. Like the internal
        scaling update, this accumulates rather than replaces, so the
        cumulative map becomes more accurate the more often it is called.

        Returns
        -------
        None

        """

        local_distance_map = self.analysis_get_instantaneous_distance_map()
        self.distance_map.update_distance_map(local_distance_map)

    #-----------------------------------------------------------------
    #
    def analysis_get_cumulative_distance_map(self):        
        """
        Returns the chain's cumulative distance map obtained over the entire ensemble
        from the update operations

        Returns
        -------
        np.ndarray
            The ensemble-averaged inter-residue distance map of shape
            ``(seq_len, seq_len)`` accumulated via repeated calls to
            ``analysis_update_distance_map``.

        """
        return self.distance_map.get_distance_map()


    #####################################################################################################
    ## END TO END DISTANCE ANALYSIS FUNCTIONS
    ##

    def analysis_get_end_to_end_distance(self):        
        """
        Returns the chain's current end-to-end distance on the lattice

        Computed as the Cartesian distance between the first and last beads of
        the chain made whole (see :meth:`get_analysis_positions`).

        Returns
        -------
        float
            The end-to-end distance between the first and last beads of the
            chain.

        """

        positions = self.get_analysis_positions()
        start  = positions[0]
        end    = positions[-1]

        return lattice_analysis_utils.get_inter_position_distance(
            start, end, self.dimensions, pbc_correction=False)


    #####################################################################################################
    ## Positional analysis
    ##

    def analysis_get_residue_residue_distance(self, R1, R2, positions=None):
        """
        Returns the inter-residue position as defined by the two
        positions here

        Computes the Cartesian distance, on the chain made whole, between the beads at sequence
        indices ``R1`` and ``R2``.

        Parameters
        ----------
        R1 : int
            Sequence index of the first residue.

        R2 : int
            Sequence index of the second residue.

        positions : list, optional
            The chain's whole-chain analysis positions (``get_analysis_positions``)
            if the caller already holds them. Default None, in which case they are
            computed here.

        Returns
        -------
        float
            The inter-residue distance between beads ``R1`` and ``R2``.

        """

        # a caller measuring many pairs on the same configuration passes the
        # whole-chain positions in once rather than re-walking the chain per pair
        if positions is None:
            positions = self.get_analysis_positions()
        start  = positions[R1]
        end    = positions[R2]

        return lattice_analysis_utils.get_inter_position_distance(
            start, end, self.dimensions, pbc_correction=False)


        
    #####################################################################################################
    ## RADIUS OF GYRATION ANALYSIS FUNCTIONS
    ##

    def analysis_get_radius_of_gyration(self):
        """
        Returns the chain's current radius of gyration.

        Computes the chain's polymeric properties from the chain made whole
        (see :meth:`get_analysis_positions`) and returns the first element,
        which is the radius of gyration.

        Returns
        -------
        float
            The radius of gyration of the chain.

        """
        return lattice_analysis_utils.get_polymeric_properties(
            self.get_analysis_positions(), self.dimensions, pbc_correction=False)[0]



    def analysis_get_polymeric_properties(self):                        
        """
        Function that returns the polymeric properties associated with a chain. These properties are in a list so 
        additional properties can be added as we go through.

        [0] - Radius of gyration
        [1] - Asphericity

        Both are computed from the gyration tensor of the chain made whole (see
        :meth:`get_analysis_positions`), which is exact for any chain that does
        not percolate the box.

        ## NOTE: Finite size detection!
        Under periodic boundaries a chain whose extent along some axis exceeds
        half the box is almost certainly interacting with its own periodic
        image, and any minimum-image quantity (inter-chain distances, contact
        and cluster analyses) becomes ambiguous for it. That regime is detected
        here directly from the whole-chain extent and reported once per chain.
        The previous check compared a centre-of-mass image selection against a
        single-image reconstruction: both tore such chains in exactly the same
        way, so the two numbers always agreed and the warning could never fire.

        Returns
        -------
        list of float
            ``[radius_of_gyration, asphericity]`` computed on the whole chain.
            As a side effect, the first time a periodic chain is found to span
            more than half the box along any axis a finite-size warning is
            printed (once per chain).

        """

        positions = self.get_analysis_positions()

        polymeric_props = lattice_analysis_utils.get_polymeric_properties(
            positions, self.dimensions, pbc_correction=False)

        # Hardwall coordinates already form a single, non-periodic image and a
        # chain cannot see its own image, so the finite-size check is only
        # meaningful for periodic systems.
        if self.hardwall or getattr(self, '_finite_size_warned', False):
            return polymeric_props

        pos_array = np.asarray(positions)
        extent = pos_array.max(axis=0) - pos_array.min(axis=0)
        over = [(axis, int(extent[axis]), int(self.dimensions[axis]))
                for axis in range(len(self.dimensions))
                if extent[axis] > 0.5 * self.dimensions[axis]]

        if over:
            self._finite_size_warned = True
            detail = ', '.join('axis %i: extent %i of box %i' % item for item in over)
            print("\n[WARNING]: Chain %s spans more than half the box (%s). Its radius of gyration, asphericity, end-to-end distance, internal scaling and distance map are computed on the chain made whole and remain exact, but a chain this extended is almost certainly interacting with its own periodic image, and minimum-image inter-chain quantities (contacts, clusters) are ambiguous for it. Rethink the box size relative to the chain dimensions. This message is printed once per chain.\n"
                  % (str(self.chainID), detail))

        return polymeric_props



            
                
                

            
        
    




        

                                
            


            


        
            
            
