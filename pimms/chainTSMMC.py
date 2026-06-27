## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................


import random
import numpy as np

from . import CONFIG
from .latticeExceptions import MoveException

##
## Temperature Sweep Metropolis Monte Carlo (TSMMC) is a class of move where the temperature is gradually increased and then decreased. Conceptually,
## this is a fairly simple idea,
##
##    |
##  T |                ooo 
##  E |             ooo   ooo
##  M |          ooo         ooo
##  P |       ooo               ooo
##    |    ooo                     ooo
##    | XXXX                         XXX
##    +-------------------------------------
##            SIMULATION PROGRESS
## 
## The X's here are the main simulation chain, and the o are conformations in the auxillary chain TSMMC move. The TSMMC move is performed as 'normal'
## MC moves as we slowly increase and then decrease the temperature. At the end we accept or reject this new conformation (with an appropriate correction
## to maintain detailed balance), meanin the full TSMMC chain is a kind of 'super' move that leads to a re-arrangement of the chain.
##
## For these moves one can control
## a) The number of steps in each subchain at each temperature
## b) The maximum temperature jumped to
## c) The number of temperatures passed through on the up and down slope (more points - more expensive but smoother transition)
##
## This move can be implemented in three distinct ways (all of which are done in PIMMS)
##
## 1) Chain TSMMC
##    For chain based TSMMC we select a single chain, and for that chain only local moves are performed as the temperature increases and then decreases. 
##    In these chain-based moves ONLY local chain pertubations are performed (i.e. chain wiggling) which is probably not the smarest way to do things
##    but these moves are SO efficient that it allows good re-arrangement of a single chain without MASSIVE decorrelation which would cause these moves
##    to be rejected.
## 
##   The chain moves are implemented as a move in the MOVER object (Chain_based_TSMMC) which calls functionality from a TSMMC object for a few things
##
## 2) Multichain TSMMC
##    Exactly the same as the chain based, except a random selection of chains are chosen, rather than a single chain. In previous versions we selected
##    a specific cluster, but this has the unfortunate consequence of breaking detailed balance as you can't randomly select a cluster of the same size
##    if a cluster merges with another cluster. This move provides a good balance for the re-arrangement of multiple chains simultaneouly - useful if chains 
##    are stuck in a coperative minima - without requiring the ENTIRE system to re-arrange. In various tests we found this move in particular offers a very
##    effective way to systematically move 'down' rough energy landscapes by alleviating on of MCs major drawbacks (the lack of concerted movement) WITHOUT
##    making any assumptions regarding the nature or the players in that concerted movement.
##
## 3) System TSMMC
##    Fundementally different from the chain and multichain based approaches, the system level TSMMC leads to the full lattice being backed up, and then
##    we sequentially alter the main chain temperature while not incrementing the main counter or performing analysis/IO. This basically converst the
##    main Markov chain into a series of auxillary chains, where almost all Monte Carlo moves are available (with the exception of the TSMMC moves - i.e.
##    we do not allow nested TSMMC behaviour). Once the full sweep is done the new system configuration is accepted or rejected 
##
##

class TSMMC:

    def __init__(self, target_temperature, jump_temp, interp_mode, step_multiplier, number_points, fixed_offset):
        """
        Class which allows the pre-computation of temperature-switch MMC moves, basically
        means we pre-calculated the switching schedule here once and then re-call it on
        each move rather than recalculating the same thing millions of times for no good
        reason.

        Parameters
        ----------
        target_temperature : float
            The main temperature the system is running at (the temperature the
            excursion departs from and returns to).

        jump_temp : float
            Default temperature to jump to upon TSMMC moves (the top of the
            excursion). Ignored if ``fixed_offset`` is set.

        interp_mode : str
            Interpolation mode through which the temperature shift happens.
            Currently only ``'LINEAR'`` is supported (the schedule is only built
            when this is ``'LINEAR'``).

        step_multiplier : int or float
            Multiplier for the number of MC steps performed per temperature in
            the schedule.

        number_points : int
            Number of points between the target and the maximum (jump)
            temperature used to define the up/down temperature ramp. More points
            give a smoother (but more expensive) transition.

        fixed_offset : float or False
            If truthy, this offset is added to the target temperature to define
            the jump temperature (overriding ``jump_temp``). This is useful when
            a quench is being run, so that as the target temperature changes the
            jump temperature changes in a linearly proportional manner.

        Returns
        -------
        None
        """

        self.mode = interp_mode
        self.steps_per_quench_multiplier = step_multiplier
        self.target_temperature = target_temperature
        self.inv_target_temperature = CONFIG.INVTEMP_FACTOR/(target_temperature)
        self.fixed_offset = fixed_offset
            
        
        if interp_mode == 'LINEAR':
            
            # get the schedule rounded to 2 decimalplaces going from the jump temperature to the target temperature
        
            # If we're using a fixed value to jump to
            if not fixed_offset:
                dT = jump_temp - target_temperature        
                step = dT/float(number_points)       

            # else if we're using a fixed offset from the target temperature
            else:
                dT = (target_temperature+fixed_offset) - target_temperature        
                step = dT/float(number_points)            
                jump_temp=target_temperature+fixed_offset
            
            self.true_temp_schedule = np.around(np.hstack((np.arange(target_temperature+step, jump_temp+(step), step), np.repeat(jump_temp, CONFIG.TOP_TEMP), np.flipud(np.arange(target_temperature+step, jump_temp+(step), step)))),5)

            self.inv_temperature_schedule = CONFIG.INVTEMP_FACTOR/self.true_temp_schedule
            

                
    #-----------------------------------------------------------------
    #
    def accept_tempered_transition(self, log_work):
        """
        Acceptance for a tempered-transitions / NCMC temperature
        excursion (Neal 1996; Nilmeier et al. 2011).

        Parameters
        ----------
        log_work : float
            The accumulated tempered-transitions work,

                log_work = sum_over_temperature_changes (beta_before - beta_after) * U(x)

            where U(x) is the system energy AT THE INSTANT of each temperature
            change (i.e. before the subsequent propagation at the new
            temperature). The sum runs over every temperature change in the
            excursion, including the initial change off the target temperature
            and the final change back onto it. For a no-op excursion (no moves
            accepted) this telescopes to exactly 0, so the move is always
            accepted, as it must be.

        Returns
        -------
        bool
            True if the whole excursion is accepted.
        """
        if log_work >= 0.0:
            return True
        return random.random() < np.exp(log_work)


    #-----------------------------------------------------------------
    #
    def start_system_TSMMC(self, backup_tuple, original_energy, ACC):
        """
        Initialize the bookkeeping for a system-wide TSMMC temperature excursion.

        Records a backup of the lattice state and the starting energy, resets the
        per-excursion move and temperature-index counters, and initializes the
        tempered-transitions / NCMC work accumulator (starting from the target
        temperature). Also captures the running total of auxiliary-chain moves so
        the number performed during this excursion can be computed at the end.

        Parameters
        ----------
        backup_tuple : tuple
            A backup of the lattice state. The first three elements are stored
            (as ``self.system_move_original_info``) so the system can be reverted
            if the excursion is rejected.

        original_energy : float
            The system energy at the start of the excursion.

        ACC : AcceptanceCalculator
            The acceptance calculator, queried for the current total number of
            auxiliary-chain moves attempted so far.

        Returns
        -------
        None
        """

        self.system_move_count = 0
        self.system_move_temp_idx = 0
        self.system_move_original_info = (backup_tuple[0], backup_tuple[1], backup_tuple[2])
        self.system_move_original_energy = original_energy

        # Tempered-transitions / NCMC bookkeeping for the system-wide excursion.
        # The system starts at the target temperature; work is accumulated at
        # each temperature change in check_in_system_TSMMC.
        self.system_log_work = 0.0
        self.system_prev_inv = self.inv_target_temperature

        # compute the total number of moves that had previously been made during
        # all TSMMC moves
        self.system_move_original_summed_aux_moves = ACC.get_total_aux_chain_moves()
            



    #-----------------------------------------------------------------
    #
    def check_in_system_TSMMC(self, ACC, current_energy):
        """
        This is the TSMMC function that is called EVERY move and updates the
        local counter (i.e. number of steps within the auxillary chain) and
        updates the temperature in the acceptance object appropriately.

        All book-keeping associated with the system TSMMC move is done by
        by the TSMMC_coordinator object.

        current_energy is the system energy at this moment; it is needed to
        accumulate the tempered-transitions work at each temperature change.

        Parameters
        ----------
        ACC : AcceptanceCalculator
            The acceptance calculator whose temperature is updated when a
            temperature change is due (every ``steps_per_quench_multiplier``
            moves).

        current_energy : float
            The system energy at this instant, used to accumulate the
            tempered-transitions work for detailed balance at each temperature
            change.

        Returns
        -------
        AcceptanceCalculator
            The (possibly temperature-updated) acceptance calculator passed in.
        """

        # increment the general counters
        self.system_move_count = self.system_move_count+1

        # if the counter is mod-0 to the number of steps per temperature
        # then we update the temperature
        if self.system_move_count % self.steps_per_quench_multiplier == 0:

            new_inv = CONFIG.INVTEMP_FACTOR / self.true_temp_schedule[self.system_move_temp_idx]

            # tempered-transitions work for the change system_prev_inv -> new_inv,
            # evaluated at the current system energy (before propagating at new_inv).
            # This is required for detailed balance (see accept_tempered_transition).
            self.system_log_work = self.system_log_work + (self.system_prev_inv - new_inv) * current_energy
            self.system_prev_inv = new_inv

            ACC.update_temperature(self.true_temp_schedule[self.system_move_temp_idx])
            self.system_move_temp_idx = self.system_move_temp_idx+1

        # regardless of if updated or not return the ACC object
        return ACC

    #-----------------------------------------------------------------
    #        
    def system_move_complete(self):
        """
        Evaluate whether the full system-wide TSMMC excursion is complete.

        The excursion is complete once the current sub-move is the last sub-move
        at the final temperature in the schedule.

        Returns
        -------
        bool
            True if the full system TSMMC move is complete (i.e. the schedule
            has been fully traversed), False otherwise.
        """

        if ((self.system_move_count + 1) % self.steps_per_quench_multiplier == 0) and (len(self.true_temp_schedule) == (self.system_move_temp_idx+1)):
            return True
        else:
            return False


    #-----------------------------------------------------------------
    #
    def system_move_finalize(self, ACC):
        """
        Reset the TSMMC object to a neutral state after a system-wide excursion.

        Folds the number of moves made during this specific TSMMC auxiliary
        chain into the acceptance object's alternative-Markov-chain counter,
        resets the per-excursion counters, attempts to explicitly free the
        backup object's memory (success of which depends on Python's memory
        management), and restores the acceptance object to the target
        temperature.

        Parameters
        ----------
        ACC : AcceptanceCalculator
            The acceptance calculator to update (move counters and temperature).

        Returns
        -------
        AcceptanceCalculator
            The updated acceptance calculator, restored to the target
            temperature.
        """

        # compute the number of moves made during this specific TSMMC auxillary chain
        # by computing the total NOW and subracting off the total before the move 
        ACC.alt_Markov_chain_moves = ACC.alt_Markov_chain_moves + (ACC.get_total_aux_chain_moves() - self.system_move_original_summed_aux_moves)

        self.system_move_count = 0
        self.system_move_temp_idx = 0
        self.system_move_original_energy = 0

        # the following block is to try and encourage the memory 
        # management to free up the memory used by the backup 
        # object - may or may not work..
        del self.system_move_original_info
        import gc        
        gc.collect()

        ACC.update_temperature(self.target_temperature)
    
        return ACC

    #-----------------------------------------------------------------
    #
    def accept_system_TSMMC(self, current_energy):
        """
        Function which accepts or rejects the GLOBAL TSMMC move.

        Uses the tempered-transitions / NCMC acceptance: the work accumulated at
        every temperature change during the excursion (self.system_log_work),
        plus the final change from the last schedule temperature back to the
        target temperature, evaluated at the final energy.

        Parameters
        ----------
        current_energy : float
            The final system energy at the end of the excursion, used to
            evaluate the work of the final temperature change back to the
            target temperature.

        Returns
        -------
        bool
            True if the whole system-wide TSMMC excursion is accepted, False
            otherwise.
        """
        final_log_work = self.system_log_work + (self.system_prev_inv - self.inv_target_temperature) * current_energy
        return self.accept_tempered_transition(final_log_work)

