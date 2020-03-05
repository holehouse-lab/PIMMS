## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

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
    def __init__(self, temp, keyword_lookup):

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
                     'MOVE_RATCHET_PIVOT', 'MOVE_SYSTEM_TSMMC', 'MOVE_JUMP_AND_RELAX']:
            self.random_thresholds[MOVE] = [rangepos, rangepos+keyword_lookup[MOVE]]            
            rangepos =  rangepos + keyword_lookup[MOVE]

        # NOTE update this (hard-coded value) when a new move is added - delibertly left here to
        # ensure this code is updated appropriately
        NUM_MOVES=13
        
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
        11   | RATCHET PIVOT
        12   | TSMMC SYSTEM re-arrangement
        13   | JUMP AND RELAX

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

        # Ratchet pivot 
        if  self.random_thresholds['MOVE_RATCHET_PIVOT'][0]<= SELECTOR < self.random_thresholds['MOVE_RATCHET_PIVOT'][1]:
            rval= 11

        # Ratchet pivot 
        if  self.random_thresholds['MOVE_SYSTEM_TSMMC'][0]<= SELECTOR < self.random_thresholds['MOVE_SYSTEM_TSMMC'][1]:
            rval= 12

        # Jump and relax 
        if  self.random_thresholds['MOVE_JUMP_AND_RELAX'][0]<= SELECTOR < self.random_thresholds['MOVE_JUMP_AND_RELAX'][1]:
            rval= 13

            
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
        Accept or reject a move based on dE and the Boltzmann acceptance criterion 

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
        Function which allows the Acceptance objects temperature to be 
        updated dynamically

        """
        self.temperature = temp
        self.invtemp = CONFIG.INVTEMP_FACTOR/(temp)


    #-----------------------------------------------------------------
    #
    def update_move_logs(self, selection, acceptance):
        """
        Function which will update a log of which moves are attempted
        and which moves are accepted. move_count and accepted_count are
        object lists where the 'selection' defines which move is 
        being updated (1 = crankshaft, 2 = chain translate etc). Note 
        that there is a 0th element, but it's just ignored.

        """

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
        Function which will update a log of which moves are attempted
        and which moves are accepted. move_count and accepted_count are
        object lists where the 'selection' defines which move is 
        being updated (1 = crankshaft, 2 = chain translate etc). Note 
        that there is a 0th element, but it's just ignored.

        """
        
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
        Function which will update the counter keeping track of all
        proposed moves performed in the alternative Markov chains 
        in various submoves (e.g. TSMMC and ratchet_pivot moves).

        This is for performance reasons only (i.e. correctly comparing
        the number of independent accept/reject events between
        different simulation stratergies).

        """

        if self.auxillary_chain:
            self.aux_chain_alt_Markov_chain_moves = self.aux_chain_alt_Markov_chain_moves + number_tried
        else:
            self.alt_Markov_chain_moves = self.alt_Markov_chain_moves + number_tried
     

   
    #-----------------------------------------------------------------
    #            
    def get_total_aux_chain_moves(self):
        """
        Calculate the total number of moves attempted in auxillary Markov chains
        (note this does not consider accept reject information, just number attempted).


        """
        n_moves = len(self.aux_chain_move_count)
        tmp=0
        for i in range(0, n_moves):
            tmp = tmp + self.aux_chain_move_count[i]
        tmp = tmp + self.aux_chain_alt_Markov_chain_moves

        return tmp


    
