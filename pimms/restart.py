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
    """Return a normalized, positive 2D/3D integer dimension list."""
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
            [0] = number of chains (int)
            [1] = chain sequence (str)

        log : bool
            Flag which, if set to true, means if this seq already has a 
            chainType defined but the passed chainType is a DIFFERENT value
            it'll warn the user about this.

        Returns
        ----------------
        None
            No return type, but the internal self.extra_chains dictionary 
            will be appropriately updated.
        

        """
        # extract info and raise exception in a civilized way
        try:
            count = int(extra_chains[0])
            chain_seq = str(extra_chains[1])
        except (TypeError, ValueError, IndexError, KeyError):
            raise RestartException(f'ERROR parsing EXTRA_CHAINS keyword [{extra_chains}] - could not parse into chain count and chain sequence')

        if count <= 0:
            raise RestartException(f'ERROR parsing EXTRA_CHAINS keyword [{extra_chains}] - chain count must be a positive integer')

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
    def set_energy(self, energy):
        """
        Set the RestartObject's stored energy value.

        Parameters
        ----------
        energy : float
            The system energy to record in the restart object (written out when
            the restart file is saved).

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
            If applying the offset would place a bead outside the new lattice
            (in which case the prior dimensions are restored before re-raising).
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



    
