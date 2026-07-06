## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................
# 


import numpy as np
import random

from . import CONFIG

from .latticeExceptions import AcceptanceException

class AcceptanceCalculator:
    """
    The AcceptanceCalculator object is a object associated with each 
    simulation (though if Hamiltonian switch MC is implemented later
    then you'd probably want one AcceptanceCalculator per Hamiltonian)

    The object performs the following functions

    1) Keeps track of the simulation temperature 

    2) Implements move acceptance/rejection (boltzmann_acceptance)

    3) Record move attempts and success for post-hoc analysis (update_move_logs()

    4) Perform the actual randomized move selection (move_selector())

    """

    #-----------------------------------------------------------------
    #
    def _validate_temperature(self, temp):
        """
        Ensure a temperature is physically meaningful and safe for inverse scaling.

        Parameters
        ----------
        temp : float
            The candidate temperature to validate.

        Returns
        -------
        None

        Raises
        ------
        AcceptanceException
            If ``temp`` is not strictly greater than zero (which would make the
            inverse-temperature scaling undefined or non-physical).
        """
        if temp <= 0:
            raise AcceptanceException("Temperature must be > 0")


    #-----------------------------------------------------------------
    #
    def _validate_selection(self, selection):
        """
        Validate a move-log index to avoid accidental negative-indexing semantics.

        The move log lists reserve element 0 as an unused placeholder, so valid
        move selections run from 1 to ``len(self.move_count) - 1`` inclusive.

        Parameters
        ----------
        selection : int
            The move-type index to validate (1 = crankshaft, 2 = chain
            translate, etc.).

        Returns
        -------
        None

        Raises
        ------
        AcceptanceException
            If ``selection`` is less than 1 or greater than or equal to the
            length of the move-count list.
        """
        if selection < 1 or selection >= len(self.move_count):
            raise AcceptanceException(
                f"Invalid move selection index {selection}; expected 1-{len(self.move_count)-1}"
            )


    #-----------------------------------------------------------------
    #
    def __init__(self, temp, keyword_lookup):
        """
        Initialize the acceptance calculator for a simulation.

        Sets the simulation temperature (and the derived inverse temperature),
        builds the cumulative random-number thresholds used to select moves at
        the frequencies requested in the keyfile, and initializes the per-move
        attempt / acceptance bookkeeping (including the parallel set of counters
        used while running inside an auxiliary TSMMC Markov chain).

        Parameters
        ----------
        temp : float
            The simulation temperature (must be > 0). Used to set
            ``self.temperature`` and ``self.invtemp``.

        keyword_lookup : dict
            Mapping of keyfile keywords to values. The ``MOVE_*`` entries
            (relative move frequencies/weights) are read to construct the
            move-selection thresholds.

        Returns
        -------
        None

        Raises
        ------
        AcceptanceException
            If ``temp`` is not strictly greater than zero.
        """

        self._validate_temperature(temp)

        self.temperature = temp
        self.auxillary_chain = False
        

        #self.invtemp = 1.0/(temp*8.3144621) # for kj
        # WARNING should update invtemp here and elsewhere to be defined as
        # k/temp and then set k in a CONFIG file rather than using 1
        self.invtemp = CONFIG.INVTEMP_FACTOR/(temp) # for kT

        self.random_thresholds={}
        rangepos = 0.0

        # construct an upper and lower bounds which ensures a random number between 0 and 1 will fall into each
        # of the MOVE_ types at the appropriate frequency based on the width of the interval defined by
        # the bounds
        for MOVE in ['MOVE_CRANKSHAFT', 'MOVE_CHAIN_TRANSLATE', 'MOVE_CHAIN_ROTATE','MOVE_CHAIN_PIVOT','MOVE_HEAD_PIVOT',
                     'MOVE_SLITHER', 'MOVE_CLUSTER_TRANSLATE','MOVE_CLUSTER_ROTATE', 'MOVE_CTSMMC', 'MOVE_MULTICHAIN_TSMMC',
                     'MOVE_PULL', 'MOVE_SYSTEM_TSMMC', 'MOVE_JUMP_AND_RELAX', 'MOVE_VMMC']:
            self.random_thresholds[MOVE] = [rangepos, rangepos+keyword_lookup[MOVE]]
            rangepos =  rangepos + keyword_lookup[MOVE]

        # NOTE update this (hard-coded value) when a new move is added - delibertly left here to
        # ensure this code is updated appropriately
        NUM_MOVES=14
        
        # we need the +1 because the '0th' move is not counted, but to keep things we just leave the 0th
        # vector position empty and occupy positions from 1 to NUM_MOVES
        self.move_count     = [0]*(NUM_MOVES+1)
        self.accepted_count = [0]*(NUM_MOVES+1)

        self.aux_chain_move_count = [0]*(NUM_MOVES+1)
        self.aux_chain_accepted_count = [0]*(NUM_MOVES+1)

        self.aux_chain_alt_Markov_chain_moves = 0
        self.alt_Markov_chain_moves = 0


    #-----------------------------------------------------------------
    #            
    def move_selector(self, chain_length):
        """Based on the MOVESET defined the keyfile, returns a value between
        1 and 8 which corresponds to a specific move type, as defined below.
        Selection of these numbers is defined based on the frequency specificed
        in the keyfile by the MOVE_* parameters

        If the chain selected is of length 1 then the moveset defaults to either
        a chain translation (rotation makes no sense) OR a cluster rotation or 
        translation. This behaviour is actually hardcoded - specifically because
        if you have a system with various chains then this allows optimal move
        behaviour for single-residue chains vs. multiresidue chains.

        In the future it might be wise to allow each chain-type to have a move
        set defined with it.


        CODE | MOVE ______________________________ 
        1    | CRANKSHAFT 
        2    | CHAIN TRANSLATE
        3    | CHAIN ROTATE
        4    | CHAIN PIVOT
        5    | HEAD PIVOT
        6    | SLITHER
        7    | CLUSTER TRANSLATE
        8    | CLUSTER ROTATE
        9    | TSMMC CHAIN re-arrangement
        10   | TSMMC MULTICHAIN re-arrangement
        11   | PULL
        12   | TSMMC SYSTEM re-arrangement
        13   | JUMP AND RELAX
        14   | VMMC

        Parameters
        ----------
        chain_length : int
            The length of the chain selected to be moved. If this is 1, moves
            that are meaningless for a single bead (chain rotate, chain pivot,
            head pivot, slither and pull; codes 3, 4, 5, 6, 11) are remapped to
            a crankshaft move (code 1).

        Returns
        -------
        int
            The selected MoveType code (an integer between 1 and 14). If the
            calculator is currently operating inside an auxiliary TSMMC chain
            (``self.auxillary_chain`` is True), the nested TSMMC moves (codes 9,
            10, 12) are remapped to a crankshaft move to avoid nesting TSMMC
            excursions.

        Raises
        ------
        AcceptanceException
            If no valid move could be selected, which indicates a bug in how
            the move-selection thresholds were constructed.
        """

        # for all other chain lengths...

        SELECTOR = random.random()  
        rval = -1
            
        # crankshaft move
        if 0.0 <= SELECTOR < self.random_thresholds['MOVE_CRANKSHAFT'][1]:
            rval = 1

        # chain translate move
        if  self.random_thresholds['MOVE_CHAIN_TRANSLATE'][0]<= SELECTOR < self.random_thresholds['MOVE_CHAIN_TRANSLATE'][1]:
            rval = 2

        # chain rotation move
        if  self.random_thresholds['MOVE_CHAIN_ROTATE'][0]<= SELECTOR < self.random_thresholds['MOVE_CHAIN_ROTATE'][1]:
            rval = 3

        # chain pivot
        if  self.random_thresholds['MOVE_CHAIN_PIVOT'][0]<= SELECTOR < self.random_thresholds['MOVE_CHAIN_PIVOT'][1]:
            rval = 4

        # head pivot
        if  self.random_thresholds['MOVE_HEAD_PIVOT'][0]<= SELECTOR < self.random_thresholds['MOVE_HEAD_PIVOT'][1]:
            rval = 5
            
        # chain slither
        if  self.random_thresholds['MOVE_SLITHER'][0]<= SELECTOR < self.random_thresholds['MOVE_SLITHER'][1]:
            rval = 6

        # cluster translate
        if  self.random_thresholds['MOVE_CLUSTER_TRANSLATE'][0]<= SELECTOR < self.random_thresholds['MOVE_CLUSTER_TRANSLATE'][1]:
            rval = 7

        # cluster rotate
        if  self.random_thresholds['MOVE_CLUSTER_ROTATE'][0]<= SELECTOR < self.random_thresholds['MOVE_CLUSTER_ROTATE'][1]:
            rval = 8

        # chain Temperature Sweep Metropolis Monte Carlo
        if  self.random_thresholds['MOVE_CTSMMC'][0]<= SELECTOR < self.random_thresholds['MOVE_CTSMMC'][1]:
            rval = 9

        # Multichain-based Temperature Sweep Metropolis Monte Carlo
        if  self.random_thresholds['MOVE_MULTICHAIN_TSMMC'][0]<= SELECTOR < self.random_thresholds['MOVE_MULTICHAIN_TSMMC'][1]:
            rval = 10

        # Pull
        if  self.random_thresholds['MOVE_PULL'][0]<= SELECTOR < self.random_thresholds['MOVE_PULL'][1]:
            rval= 11

        # Ratchet pivot 
        if  self.random_thresholds['MOVE_SYSTEM_TSMMC'][0]<= SELECTOR < self.random_thresholds['MOVE_SYSTEM_TSMMC'][1]:
            rval= 12

        # Jump and relax
        if  self.random_thresholds['MOVE_JUMP_AND_RELAX'][0]<= SELECTOR < self.random_thresholds['MOVE_JUMP_AND_RELAX'][1]:
            rval= 13

        # VMMC (virtual-move Monte Carlo collective cluster move)
        if  self.random_thresholds['MOVE_VMMC'][0]<= SELECTOR < self.random_thresholds['MOVE_VMMC'][1]:
            rval= 14


        if rval == -1:
            print(SELECTOR)
            raise AcceptanceException('ERROR: Found ourselves without a correct selection - suggests a bug in how the moveset randomization is done!')

        # if single particle there are a few moves which become equivalent to the crankshaft...
        if chain_length == 1:
            if rval == 3 or rval == 4 or rval == 5 or rval == 6 or rval == 11:
                rval = 1 

        # if we're running *inside* a TSMMC system spanning move then we DO NOT perform any of the other TSMMC moves
        # and default to a crankshaft move (avoids nesting TSMMC moves!!)
        if self.auxillary_chain:
            if rval == 9 or rval == 10 or rval == 12:
                rval = 1

        return rval
            

            
           

    #-----------------------------------------------------------------
    #
    def boltzmann_acceptance(self, old_energy, new_energy):
        """
        Accept or reject a move based on dE and the Boltzmann acceptance criterion.

        Downhill (or energy-neutral) moves are always accepted. Uphill moves are
        accepted with probability ``exp(-(new_energy - old_energy) * invtemp)``,
        compared against a uniform random draw.

        Parameters
        ----------
        old_energy : float
            The system energy before the proposed move.

        new_energy : float
            The system energy after the proposed move.

        Returns
        -------
        bool
            True if the move is accepted, False otherwise.
        """

        if new_energy <= old_energy:
            return True
        else:
            expterm = np.exp(-(new_energy-old_energy)*self.invtemp)
            
            if random.random() < expterm:
                return True

            else:
                return False


    #-----------------------------------------------------------------
    #
    def update_temperature(self, temp):
        """
        Update the acceptance object's temperature dynamically.

        Re-validates the new temperature and recomputes the cached inverse
        temperature (``self.invtemp``). Used by quenching and by the TSMMC
        temperature-excursion machinery.

        Parameters
        ----------
        temp : float
            The new temperature (must be > 0).

        Returns
        -------
        None

        Raises
        ------
        AcceptanceException
            If ``temp`` is not strictly greater than zero.
        """
        self._validate_temperature(temp)
        self.temperature = temp
        self.invtemp = CONFIG.INVTEMP_FACTOR/(temp)


    #-----------------------------------------------------------------
    #
    def update_move_logs(self, selection, acceptance):
        """
        Update the log of which moves are attempted and which are accepted.

        ``move_count`` and ``accepted_count`` are object lists where
        ``selection`` defines which move is being updated (1 = crankshaft,
        2 = chain translate, etc). There is a 0th element, but it is just
        ignored. If the calculator is currently running inside an auxiliary
        TSMMC chain (``self.auxillary_chain``), the parallel
        ``aux_chain_*`` counters are updated instead.

        Parameters
        ----------
        selection : int
            The MoveType code of the move being logged.

        acceptance : bool
            True if the move was accepted (increments the accepted counter),
            False otherwise (only the attempt counter is incremented).

        Returns
        -------
        None

        Raises
        ------
        AcceptanceException
            If ``selection`` is out of the valid 1..NUM_MOVES range.
        """

        self._validate_selection(selection)

        if self.auxillary_chain:

            # increment the move count
            self.aux_chain_move_count[selection] = self.aux_chain_move_count[selection] + 1

            if acceptance:
                self.aux_chain_accepted_count[selection] = self.aux_chain_accepted_count[selection] + 1

        else:
        
            # increment the move count
            self.move_count[selection] = self.move_count[selection] + 1

            if acceptance:
                self.accepted_count[selection] = self.accepted_count[selection] + 1


    #-----------------------------------------------------------------
    #
    def megastep_update_move_logs(self, selection, number_accepted, number_tried):
        """
        Bulk-update the move logs for a megamove that performed many sub-moves.

        Like :meth:`update_move_logs`, but increments the attempt and accepted
        counters by arbitrary amounts in a single call - used by the megamoves
        (crankshaft / slither / pull) which internally perform and accept/reject
        many sub-moves per call. ``move_count`` and ``accepted_count`` are object
        lists where ``selection`` defines which move is being updated
        (1 = crankshaft, 2 = chain translate, etc); the 0th element is ignored.
        If running inside an auxiliary TSMMC chain the parallel ``aux_chain_*``
        counters are updated instead.

        Parameters
        ----------
        selection : int
            The MoveType code of the move being logged.

        number_accepted : int
            The number of sub-moves that were accepted.

        number_tried : int
            The number of sub-moves that were attempted/proposed.

        Returns
        -------
        None

        Raises
        ------
        AcceptanceException
            If ``selection`` is out of the valid 1..NUM_MOVES range.
        """

        self._validate_selection(selection)

        # increment the move and accepted count for aux chain counts
        if self.auxillary_chain:
            self.aux_chain_move_count[selection] = self.aux_chain_move_count[selection] + number_tried
            self.aux_chain_accepted_count[selection] = self.aux_chain_accepted_count[selection] + number_accepted
        
        # increment the move and accepted count for main chain counts
        else:
            self.move_count[selection] = self.move_count[selection] + number_tried
            self.accepted_count[selection] = self.accepted_count[selection] + number_accepted


    #-----------------------------------------------------------------
    #
    def alt_Markov_chain_update_move_logs(self, number_tried):
        """
        Update the counter tracking proposed moves made in alternative Markov chains.

        These are the accept/reject events performed inside various sub-moves
        (e.g. the TSMMC temperature-excursion moves) that belong to an
        alternative Markov chain rather than the main one. This bookkeeping is
        for performance reasons only (i.e. correctly comparing the number of
        independent accept/reject events between different simulation
        strategies). If running inside an auxiliary TSMMC chain the parallel
        ``aux_chain_alt_Markov_chain_moves`` counter is updated instead.

        Parameters
        ----------
        number_tried : int
            The number of alternative-Markov-chain moves to add to the counter.

        Returns
        -------
        None
        """

        if self.auxillary_chain:
            self.aux_chain_alt_Markov_chain_moves = self.aux_chain_alt_Markov_chain_moves + number_tried
        else:
            self.alt_Markov_chain_moves = self.alt_Markov_chain_moves + number_tried
     

   
    #-----------------------------------------------------------------
    #            
    def get_total_aux_chain_moves(self):
        """
        Calculate the total number of moves attempted in auxiliary Markov chains.

        Sums every per-move auxiliary-chain attempt counter
        (``aux_chain_move_count``) together with the auxiliary alternative
        Markov-chain move counter (``aux_chain_alt_Markov_chain_moves``). Note
        this considers only the number of moves attempted, not accept/reject
        information.

        Returns
        -------
        int
            The total number of moves attempted across all auxiliary Markov
            chains.
        """
        n_moves = len(self.aux_chain_move_count)
        tmp=0
        for i in range(0, n_moves):
            tmp = tmp + self.aux_chain_move_count[i]
        tmp = tmp + self.aux_chain_alt_Markov_chain_moves

        return tmp


    #-----------------------------------------------------------------
    #
    def total_attempted_moves(self):
        """
        Total number of individual accept/reject events attempted across ALL
        sub-loops so far (i.e. every MC move, not just outer master-loop steps).

        Megamoves (crankshaft / slither / pull) contribute their full per-megamove
        substep counts, single-chain moves contribute one each, and the TSMMC
        temperature-excursion moves contribute their alternative-Markov-chain
        substeps. This is the same quantity written to TOTAL_MOVES.dat.

        Returns
        -------
        int
            The cumulative number of attempted MC moves.
        """
        # moves 9 and 10 (chain / multichain TSMMC) are excluded here because
        # their per-substep attempts are accumulated in alt_Markov_chain_moves;
        # move_count[9]/[10] only count whole excursions.
        total = 0
        for move in [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14]:
            total = total + self.move_count[move]
        total = total + self.alt_Markov_chain_moves
        return total


    
