## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

##
## restart
##
## The RestartObject implemements a way to read and write restart files. This allows PIMMS
## to restart from previous simulations. Other than chain position, no other state is saved.  
##


import copy
import math
import numbers
import os
import pickle

import numpy as np

from . import CONFIG
from .latticeExceptions import RestartException
from . import pimmslogger


def _validated_dimensions(dimensions, label="DIMENSIONS"):
    """
    Return a normalized, positive 2D/3D integer dimension list.

    Used to sanity check any set of lattice dimensions before they are written
    into a RestartObject, either from a restart file on disk or from a caller
    resizing the lattice. Every element must be a genuine positive integer
    (booleans are rejected explicitly) and the returned list is a fresh list of
    Python ints, so the caller never keeps a reference to the input object.

    Parameters
    ----------
    dimensions : list
        The candidate lattice dimensions. Any non-string sequence of 2 or 3
        integers is accepted (list, tuple, numpy array).

    label : str, optional
        The name used for ``dimensions`` in any exception message, so the
        caller can say which set of dimensions was bad. Default is
        ``"DIMENSIONS"`` (the restart-file keyword).

    Returns
    -------
    list
        A new list of 2 or 3 positive Python ints.

    Raises
    ------
    RestartException
        If ``dimensions`` is a string/bytes object, is not iterable, does not
        contain exactly 2 or 3 elements, or contains a value that is not a
        positive integer.

    """
    if isinstance(dimensions, (str, bytes)):
        raise RestartException(
            f"Invalid restart file - {label} must be a 2D or 3D integer sequence")
    try:
        values = list(dimensions)
    except (TypeError, ValueError):
        raise RestartException(
            f"Invalid restart file - {label} must be a 2D or 3D integer sequence")

    if len(values) not in (2, 3):
        raise RestartException(
            f"Invalid restart file - {label} must contain exactly 2 or 3 dimensions")

    normalized = []
    for value in values:
        if (isinstance(value, bool) or
                not isinstance(value, numbers.Integral) or value <= 0):
            raise RestartException(
                f"Invalid restart file - {label} values must be positive integers; got {values}")
        normalized.append(int(value))
    return normalized


class RestartObject:
    """
    Object used to read and write restart files. Restart information ONLY contains information on
    chain position, sequence, and type, plus the grid dimensions, the hardwall flag and the last
    recorded energy. It does NOT include any other simulation state (step count, move statistics,
    random number generator state, or any of the parameters read from the keyfile).

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
        Initialize an empty RestartObject.

        Sets up the internal state with zero energy, empty dimensions, a
        non-hardwall flag, and empty chain / sequence-to-chaintype / extra-chain
        containers. These are subsequently populated via one of the
        ``build_from_*`` methods or by adding extra chains.

        Returns
        -------
        None
            No return value; the new object's attributes are initialised in place.
        """
        self.energy = 0
        self.dimensions = []
        self.hardwall = False
        self.chains = {}
        self.seq2chainType = {}
        self.extra_chains = {}

        # continuation state, present only in restart files written by 1.0.8 or
        # later (older files simply leave these None): the master step the
        # snapshot was taken at, the two global random-number-generator states
        # at that moment, the temperature in force, and where the file came from
        self.step = None
        self.rng_python = None
        self.rng_numpy = None
        self.temperature = None
        self.pimms_version = None
        self.filename = None


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
            """
            Check that a single coordinate remains inside the box once offset.

            Parameters
            ----------
            pos : int
                The current coordinate of a bead along dimension ``dim``.

            dim : int
                The index of the dimension being checked (0=x, 1=y, 2=z), used
                to select both the offset and the box length.

            Returns
            -------
            bool
                True if ``pos + position_offset[dim]`` lies within
                ``[0, self.dimensions[dim])``, otherwise False.

            """
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
        # First pass: validate all proposed positions so this update is atomic.
        for chainID in self.chains:
            for position in self.chains[chainID][0]:
                for dim in range(0, n_dim):
                    if not valid_pos(position[dim], dim):
                        raise RestartException(f'Trying to offet a position on chain {chainID} from {position[dim]} to {position[dim] + position_offset[dim]} (dim={dim}) but lattice dimensions are {self.dimensions}')

        # Second pass: apply position updates once we know every position is valid.
        for chainID in self.chains:
            for position in self.chains[chainID][0]:
                for dim in range(0, n_dim):
                    position[dim] = position[dim] + position_offset[dim]
                        


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
    def add_extra_chains(self, extra_chains, log=False):
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
            [0] = number of chains (int, must be positive)
            [1] = chain sequence (str, must be non-empty)

        log : bool, optional
            Flag which, if set to true, means if this seq already has a
            chainType defined but the passed chainType is a DIFFERENT value
            it'll warn the user about this. Default is False.

        Returns
        ----------------
        None
            No return type, but the internal self.extra_chains dictionary
            will be appropriately updated. One entry is added per requested
            chain, keyed by a newly allocated chainID, with a position entry
            of None because extra chains are placed when the lattice is built.

        Raises
        ----------------
        RestartException
            If extra_chains cannot be unpacked into exactly a count and a
            sequence, if the count is not a positive integer, or if the
            sequence is not a non-empty string.

        """
        # extract info and raise exception in a civilized way
        try:
            if len(extra_chains) != 2:
                raise ValueError
            count = extra_chains[0]
            chain_seq = extra_chains[1]
        except (TypeError, ValueError, IndexError, KeyError):
            raise RestartException(f'ERROR parsing EXTRA_CHAINS keyword [{extra_chains}] - could not parse into chain count and chain sequence')

        if (isinstance(count, (bool, np.bool_)) or
                not isinstance(count, numbers.Integral) or count <= 0):
            raise RestartException(f'ERROR parsing EXTRA_CHAINS keyword [{extra_chains}] - chain count must be a positive integer')
        if not isinstance(chain_seq, str) or len(chain_seq) == 0:
            raise RestartException(
                f'ERROR parsing EXTRA_CHAINS keyword [{extra_chains}] - sequence must be a non-empty string')
        count = int(count)

        # Dynamically calculate next chainID from both base and extra chain maps.
        existing_chain_ids = list(self.chains.keys()) + list(self.extra_chains.keys())
        if len(existing_chain_ids) == 0:
            chainID = 1
        else:
            chainID = max(existing_chain_ids) + 1
            
        if chain_seq in self.seq2chainType:

            # note - this [0] means we always use the first chain type even if there are multiple
            # chain IDs associated with a specific sequence. 
            local_chainType = self.seq2chainType[chain_seq][0]
        else:

            # if a new chain dynamically calculate what the next chainType should be (next increment
            # after current highest number)
            tmp = []
            for s in self.seq2chainType:
                tmp.extend(self.seq2chainType[s])

            if len(tmp) == 0:
                existing_types = [x[2] for x in self.chains.values()] + [x[2] for x in self.extra_chains.values()]
                if len(existing_types) == 0:
                    local_chainType = 0
                else:
                    local_chainType = max(existing_types) + 1
            else:
                local_chainType = max(tmp) + 1

            # update the seq2chainType dictionary
            self.__update_seq2chainType(local_chainType, chain_seq, log)
            
        # finally, after all this set up, add to the extra_chains dict
        for c in range(count):
            self.extra_chains[chainID] = [None, chain_seq, local_chainType]
            chainID = chainID + 1
 

    #-----------------------------------------------------------------
    #       
    def set_continuation_state(self, step, rng_python, rng_numpy, temperature, pimms_version=None):
        """
        Record what a later run needs to resume this snapshot exactly.

        A configuration alone is not enough to continue a run: the resumed run
        also has to pick up the step counter where this one stopped, draw the
        same random numbers the uninterrupted run would have drawn next, and be
        at the same temperature if a quench is in progress. This stores all of
        that so that RESTART_CONTINUE can restore it.

        Parameters
        ----------
        step : int
            The master step the snapshot was taken at (after that step's move
            and its scheduled output).

        rng_python : tuple
            ``random.getstate()`` at that moment.

        rng_numpy : tuple
            ``numpy.random.get_state()`` at that moment.

        temperature : float
            The simulation temperature in force at that step (the ramp value in
            a quench run, ``TEMPERATURE`` otherwise).

        pimms_version : str, optional
            The PIMMS version writing the file, for the record.

        Returns
        -------
        None
            The fields are stored on the object and written by ``write_to_file``.
        """
        self.step = int(step)
        self.rng_python = rng_python
        self.rng_numpy = rng_numpy
        self.temperature = float(temperature)
        self.pimms_version = pimms_version

    #-----------------------------------------------------------------
    #
    def set_energy(self, energy):
        """
        Set the RestartObject's stored energy value.

        Parameters
        ----------
        energy : float
            The total system energy to record in the restart object. This is
            stored verbatim and written out under the ENERGY key when the
            restart file is saved; it is informational only and is not used to
            rebuild the lattice.

        Returns
        -------
        None
            No return value, but ``self.energy`` is updated in place.
        """
        self.energy = energy


    #-----------------------------------------------------------------
    #       
    def build_from_lattice(self, LATTICE, hardwall=False, log=False):
        """
        Construct a restart object using a lattice object to set the chain
        positions.

        Parameters
        ------------
        LATTICE : Lattice
            A standard PIMMS Lattice object. Its ``dimensions`` and its
            ``chains`` dictionary (positions, sequence and chainType of every
            Chain) are copied into this RestartObject; the positions are deep
            copied, so subsequent moves on the lattice do not alter the
            restart snapshot.

        hardwall : bool, optional
            Flag which records whether the current system uses hardwall
            boundaries (True) or periodic boundary conditions (False). Stored
            on the object and written into the restart file. Default is False.

        log : bool, optional
            Flag which, if set to True, means warnings (specifically identical
            sequences mapping to different chainTypes) are written to the
            standard PIMMS logfile. Default is False.

        Returns
        ----------
        None
            No return type, but self.dimensions, self.hardwall, self.chains and
            self.seq2chainType are overwritten and self.extra_chains is reset.


        """
        self.dimensions = list(LATTICE.dimensions)
        self.hardwall   = hardwall

        # reset chain info...
        self.chains = {}
        self.seq2chainType  = {}
        self.extra_chains = {}

        for chainID in LATTICE.chains:
        
            local_chainType = LATTICE.chains[chainID].chainType
            local_seq = LATTICE.chains[chainID].sequence

            # add the chain to the restart object
            self.chains[chainID] = [copy.deepcopy(LATTICE.chains[chainID].positions), local_seq, local_chainType]

            # udpate the self.seq2chainType dictionary
            self.__update_seq2chainType(local_chainType, local_seq, log)


    #-----------------------------------------------------------------
    #       
    def update_lattice_dimensions(self, new_dimensions, manual_offset=None):
        """
        Resize the lattice and reposition the chains within the new lattice.

        Updates the restart object's dimensions and shifts every chain position
        by a per-dimension offset. By default the offset is computed so the
        existing chains end up centred in the larger lattice; alternatively an
        explicit ``manual_offset`` can be supplied. If the offset would move any
        bead outside the new lattice, the original dimensions are restored and
        the underlying :class:`RestartException` is re-raised.

        Parameters
        ----------
        new_dimensions : list of int
            The new lattice dimensions (length 2 or 3). Should be greater than or
            equal to the current dimensions for centring to make sense.
        manual_offset : list of int, optional
            Explicit per-dimension offset to apply to every chain position. If
            ``None`` (the default), a centring offset is computed automatically
            from the difference between ``new_dimensions`` and the current
            dimensions.

        Returns
        -------
        None
            No return value, but ``self.dimensions`` and the stored chain
            positions are updated in place.

        Raises
        ------
        RestartException
            If ``new_dimensions`` is not a valid 2D/3D positive integer
            sequence or does not have the same dimensionality as the current
            box, if ``manual_offset`` is not one integer per dimension, or if
            applying the offset would place a bead outside the new lattice (in
            which case the prior dimensions are restored before re-raising).
        """

        new_dimensions = _validated_dimensions(new_dimensions, "new dimensions")
        if len(new_dimensions) != len(self.dimensions):
            raise RestartException(
                'Trying to resize restart object, but old and new dimensions do not match.')

        ## -----------
        if manual_offset is None:
            # calculate offset so the chains are placed in the center of the new lattice
            x_off = int((new_dimensions[0] - self.dimensions[0])/2)
            y_off = int((new_dimensions[1] - self.dimensions[1])/2)

            if len(new_dimensions) == 3:
                z_off = int((new_dimensions[2] - self.dimensions[2])/2)
                position_offset=[x_off, y_off, z_off]
            else:
                position_offset=[x_off, y_off]
                ## -----------

        # Manually provide the offsets to convert from old dimensions -> new dimensions
        # TODO - add check in keyfile parser that manual_offset is reasonable
        else:
            try:
                position_offset = list(manual_offset)
            except (TypeError, ValueError):
                raise RestartException('Manual restart offset must be an integer sequence.')
            if (len(position_offset) != len(new_dimensions) or
                    any(isinstance(value, bool) or
                        not isinstance(value, numbers.Integral)
                        for value in position_offset)):
                raise RestartException(
                    'Manual restart offset must contain one integer per lattice dimension.')
            position_offset = [int(value) for value in position_offset]

        # Next construct and instantiate a new restart object which has the new dimensions
        # including applying the possition offset we calculated above

        # finaly, apply the offset on this 'new' lattice (order matters, as __apply_position_offset
        # assesses if, given self.dimensions, the offset is valid or not)
        # If this fails, restore prior dimensions.
        old_dimensions = list(self.dimensions)
        self.dimensions = list(new_dimensions)
        try:
            self.__apply_position_offset(position_offset)
        except RestartException:
            self.dimensions = old_dimensions
            raise






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
            Name of the restart file to be read. This must be a pickle written
            by write_to_file().

        log : bool, optional
            Flag which, if set to True, means warnings (specifically identical
            sequences mapping to different chainTypes) are written to the
            standard PIMMS logfile. Default is False.

        Returns
        -------------
        None
            No return type, but on success updates self.dimensions,
            self.energy, self.hardwall, self.chains and self.seq2chainType, and
            resets self.extra_chains. Nothing is written until the whole file
            has validated, so a failed read leaves the object untouched.

        Raises
        -------------
        RestartException
            If the file cannot be read or unpickled, if the top-level object is
            not a dictionary, if any of the DIMENSIONS/ENERGY/HARDWALL/CHAINS
            entries are missing or malformed, or if any chain is invalid (bad
            chainID, sequence/position length mismatch, a bead outside the box,
            two beads on one site, or a disconnected chain).

        """
        # if IO issue (not IndexError often thrown if a valid file is found
        # but its not actually a pickle file!
        try:
            with open(filename, "rb") as fh:
                input_dict = pickle.load(fh)
        except Exception as e:
            raise RestartException("Error reading restart file. Error:\n\n%s" %(str(e)))
        
        # the pickle must hold the documented top-level dictionary
        if not isinstance(input_dict, dict):
            raise RestartException(
                "Invalid restart file - top-level object is %s, expected a dictionary"
                % type(input_dict).__name__)

        # Extract into locals first: a failed read must not leave a previously
        # usable RestartObject half overwritten.
        try:
            dimensions = _validated_dimensions(input_dict['DIMENSIONS'])
            energy = input_dict['ENERGY']
            hardwall = input_dict['HARDWALL']

            # local chains is a dictionary where keys are chainIDs and values are lists with three elements
            # [0] : bead positions (N->C)
            # [1] : chain sequence (which will be referenced against the parameter file)
            # [2] : chainType : a single value that defines the type of chain (many chains can have the same chainType, 
            #       but each chain has a unique chainID)
            local_chains    = input_dict['CHAINS'] 
        except KeyError as e:
            raise RestartException("Invalid restart file - missing entry for %s" % (e.args[0]))

        if not isinstance(local_chains, dict):
            raise RestartException("Invalid restart file - CHAINS entry must be a dictionary")

        if (isinstance(energy, bool) or not isinstance(energy, numbers.Real) or
                not math.isfinite(float(energy))):
            raise RestartException("Invalid restart file - ENERGY must be a finite numeric value")
        if not isinstance(hardwall, (bool, np.bool_)):
            raise RestartException("Invalid restart file - HARDWALL must be True or False")
        hardwall = bool(hardwall)

        new_chains = {}
        new_seq2chain_type = {}

        # one entry PER chain (not per chain type). Track occupancy so an
        # overlapping restart (two beads on one site - which would silently
        # desynchronise the occupancy grid from the chain objects and crash
        # deep in the mover) is rejected here with a clear message.
        _occupied = set()
        max_chain_id = np.iinfo(CONFIG.NP_INT_TYPE).max

        for chainID in local_chains:

            # Zero is the occupancy-grid solvent sentinel.  Accepting chainID=0
            # constructs a Chain object whose beads remain indistinguishable
            # from empty lattice sites; non-integral/overflowing IDs likewise
            # corrupt the fixed-width occupancy grid.
            if (isinstance(chainID, bool) or
                    not isinstance(chainID, numbers.Integral) or
                    chainID <= 0 or chainID > max_chain_id):
                raise RestartException(
                    "Invalid restart file - chainID must be a positive integer "
                    f"representable by the lattice grid; got {chainID!r}")
            chainID = int(chainID)

            # extract info for each chain
            try:
                chain_entry = local_chains[chainID]
                if len(chain_entry) != 3:
                    raise ValueError
                local_pos = chain_entry[0]
                local_seq = chain_entry[1]
                local_chainType = chain_entry[2]
            except (TypeError, IndexError, KeyError, ValueError):
                raise RestartException(f"Invalid restart file - malformed chain entry for chainID={chainID}")

            if (isinstance(local_chainType, bool) or
                    not isinstance(local_chainType, numbers.Integral) or
                    local_chainType < 0):
                raise RestartException(
                    "Invalid restart file - chainType must be a non-negative integer "
                    f"(chainID={chainID})")
            local_chainType = int(local_chainType)

            if not isinstance(local_seq, str) or len(local_seq) == 0:
                raise RestartException(
                    f"Invalid restart file - sequence must be a non-empty string (chainID={chainID})")

            if isinstance(local_pos, (str, bytes)):
                raise RestartException(
                    f"Invalid restart file - positions must be a sequence (chainID={chainID})")
            try:
                local_pos = list(local_pos)
            except (TypeError, ValueError):
                raise RestartException(
                    f"Invalid restart file - positions must be a sequence (chainID={chainID})")

            # check sequence and number of positions match
            if len(local_seq) != len(local_pos):
                raise RestartException("Invalid restart file - sequence length does not match number of positions")

            normalized_positions = []
            for position in local_pos:
                if isinstance(position, (str, bytes)):
                    raise RestartException("Invalid restart file - malformed bead position")
                try:
                    position = list(position)
                except (TypeError, ValueError):
                    raise RestartException("Invalid restart file - malformed bead position")
                if len(position) != len(dimensions):
                    raise RestartException("Invalid restart file - chain position dimensionality does not match DIMENSIONS")

                # every coordinate must be inside the box: a negative value
                # would silently WRAP onto a real cell via numpy indexing (a
                # physically wrong configuration that only crashes much later),
                # and a too-large one would die with a raw IndexError
                normalized_position = []
                for d, c in enumerate(position):
                    if (isinstance(c, bool) or
                            not isinstance(c, numbers.Integral)):
                        raise RestartException(
                            "Invalid restart file - bead coordinates must be integers "
                            f"(chainID={chainID}, position={position})")
                    c = int(c)
                    if c < 0 or c >= dimensions[d]:
                        raise RestartException(
                            "Invalid restart file - bead position %s outside box %s (chainID=%s)"
                            % (list(position), list(dimensions), chainID))
                    normalized_position.append(c)

                _key = tuple(normalized_position)
                if _key in _occupied:
                    raise RestartException(
                        "Invalid restart file - two beads occupy the same site %s (second chainID=%s)"
                        % (list(position), chainID))
                _occupied.add(_key)
                normalized_positions.append(normalized_position)

            # Consecutive beads must be neighbours in the stored boundary mode.
            # The lattice uses the Moore neighbourhood, so a diagonal step is
            # valid provided every per-axis minimum-image displacement is <= 1.
            for previous, current in zip(normalized_positions, normalized_positions[1:]):
                for dim, (a, b) in enumerate(zip(previous, current)):
                    displacement = abs(a - b)
                    if not hardwall:
                        displacement = min(displacement, dimensions[dim] - displacement)
                    if displacement > 1:
                        raise RestartException(
                            "Invalid restart file - chain is disconnected between "
                            f"{previous} and {current} (chainID={chainID})")

            new_chains[chainID] = [normalized_positions, local_seq, local_chainType]
            existing_types = new_seq2chain_type.setdefault(local_seq, [])
            if local_chainType not in existing_types:
                if existing_types and log:
                    pimmslogger.log_warning(
                        f'When reading RestartObject found identical chain [{local_seq}] '
                        'with different chainType indices. This may be undesired...')
                existing_types.append(local_chainType)

        self.dimensions = dimensions
        self.energy = float(energy) if not isinstance(energy, numbers.Integral) else int(energy)
        self.hardwall = hardwall
        self.chains = new_chains
        self.seq2chainType = new_seq2chain_type
        self.extra_chains = {}
        self.filename = filename

        # Continuation state. These keys are optional so that restart files from
        # earlier versions still load; RESTART_CONTINUE refuses a file without
        # them rather than guessing. STEP must be a non-negative integer, the
        # generator states are taken as the opaque tuples random.getstate() and
        # numpy.random.get_state() produced, and TEMPERATURE must be a positive
        # finite number.
        step = input_dict.get('STEP')
        if step is not None:
            if isinstance(step, bool) or not isinstance(step, numbers.Integral) or step < 0:
                raise RestartException("Invalid restart file - STEP must be a non-negative integer")
            step = int(step)
        temperature = input_dict.get('TEMPERATURE')
        if temperature is not None:
            if (isinstance(temperature, bool) or not isinstance(temperature, numbers.Real)
                    or not math.isfinite(temperature) or temperature <= 0):
                raise RestartException("Invalid restart file - TEMPERATURE must be a positive finite number")
            temperature = float(temperature)
        rng_python = input_dict.get('RNG_PYTHON')
        rng_numpy = input_dict.get('RNG_NUMPY')
        if (rng_python is None) != (rng_numpy is None):
            raise RestartException("Invalid restart file - RNG_PYTHON and RNG_NUMPY must be present together")
        if rng_python is not None and not isinstance(rng_python, tuple):
            raise RestartException("Invalid restart file - RNG_PYTHON is not a random.getstate() tuple")
        if rng_numpy is not None and not isinstance(rng_numpy, tuple):
            raise RestartException("Invalid restart file - RNG_NUMPY is not a numpy.random.get_state() tuple")
        self.step = step
        self.rng_python = rng_python
        self.rng_numpy = rng_numpy
        self.temperature = temperature
        self.pimms_version = input_dict.get('PIMMS_VERSION')


    #-----------------------------------------------------------------
    #       
    def write_to_file(self):
        """
        Serialize the restart object to disk as a pickle file.

        Writes a dictionary containing the chain information (``CHAINS``), lattice
        dimensions (``DIMENSIONS``), recorded energy (``ENERGY``) and hardwall
        flag (``HARDWALL``) to ``CONFIG.RESTART_FILENAME`` using :mod:`pickle`.
        Note that ``extra_chains`` are not written; only the materialised
        ``self.chains`` are saved.

        Returns
        -------
        None
            No return value; the restart data is written to
            ``CONFIG.RESTART_FILENAME``.
        """

        output={}
        output['CHAINS'] = {}
        for chainID in self.chains:
            output['CHAINS'][chainID] = self.chains[chainID]

        output['DIMENSIONS'] = self.dimensions
        output['ENERGY']     = self.energy
        output['HARDWALL']   = self.hardwall
        # continuation state (see set_continuation_state); written only when it
        # was recorded, so a restart object built from a bare lattice writes the
        # same four-key file it always did
        if self.step is not None:
            output['STEP'] = int(self.step)
        if self.rng_python is not None and self.rng_numpy is not None:
            output['RNG_PYTHON'] = self.rng_python
            output['RNG_NUMPY'] = self.rng_numpy
        if self.temperature is not None:
            output['TEMPERATURE'] = float(self.temperature)
        if self.pimms_version is not None:
            output['PIMMS_VERSION'] = str(self.pimms_version)

        # ATOMIC write: dump to a temp file in the same directory and rename it
        # over the target. A plain 'wb' open truncated the existing restart
        # BEFORE the new content was complete, so a crash mid-write (exactly the
        # scenario restart files exist for) destroyed the previous good
        # checkpoint AND left the new one unreadable. os.replace is atomic on
        # POSIX, so restart.pimms is always either the old or the new complete
        # snapshot, never a torn one.
        _tmp = CONFIG.RESTART_FILENAME + ".tmp"
        try:
            with open(_tmp, "wb") as fh:
                pickle.dump(output, fh)
            os.replace(_tmp, CONFIG.RESTART_FILENAME)
        finally:
            # A serialization error should preserve the previous checkpoint and
            # must not leave a misleading torn .tmp snapshot behind.
            if os.path.exists(_tmp):
                try:
                    os.remove(_tmp)
                except OSError:
                    pass



    
