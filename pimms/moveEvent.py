## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................

##
## moveEvent
##
## This file represents a class which defines the type of move being made and all the
## associated changes with that move. By storing move information inside an object we
## encapsulate the information for well defined access later on

from .latticeExceptions import SimulationEnergyException, MoveException

class MoveEvent:

    def __init__(self, original_positions, moved_positions, original_chain_positions, moved_chain_positions, moved_indices, move_type, pivot_point=None, chain_list=[]):
        """
        MoveEvents are basically objects that describe a bunch of information about a single_chain move for easy access. This information could just be passed as 
        a list or dictionary, but we greatly improve the code clarity by making it a specific object. It also lets us implement general functions to manipulate 
        move-based data rather than implementing those manipulations elsewhere.

        As of version 0.16 there are 12 different moves implemented in PIMMS. These are outlined below:

        MoveType Codes:

        1  | crankshaft             # Not dealt with here 
        2  | chain translate   
        3  | chain rotate
        4  | chain pivot
        5  | head pivot
        6  | chain slither
        7  | cluster translate
        8  | cluster rotate
        9  | chain-based TSMMC      # Not dealt with here
        10 | multichain-based TSMMC # Not dealt with here
        11 | ratchet pivot          # Not dealt with here
        12 | system_based TSMMC     # Not dealt with here
        13 | jump and relax

        The moves that are dealt with by a move event are classed as single_chain moves. These are moves where the 'movement' part is performed by a function
        implemented in moves.py, but the energy evaluation is done by the function single_chain_move in simulation.py. This is useful because it means there is
        a generic way to implement new moves, and as long as the move function follows the protocol described in moves.py and creates an appropriate moveEvent
        object then adding a new move is very straight forward.

        """
        
        # the original positions of the chain which have been moved (i.e. *only* those which will move)
        self.original_positions = original_positions

        # the the new positions occupied by a chain which has moved (i.e. *only* those which have moved)
        self.moved_positions = moved_positions

        # the FULL set of positions corresponding to the chain
        # which has been moved (might be the same as original_positions
        # or original_positions might just be a subset of original_chain_positions)
        self.original_chain_positions = original_chain_positions
        
        # the FULL set of positions corresponding to the chain
        # which has been moved (might be the same as translated_position
        # or translated_positions might just be a subset)
        self.moved_chain_positions =  moved_chain_positions

        # the set of indices corresponding to positions along the chain
        # which have been moved - i.e. if only the first three residues were
        # moved this would be [0,1,2]
        self.moved_indices = moved_indices

        # the type of move (determines how the energy is calculated) 
        # each move has its own MoveType code 
        self.move_type = move_type
        
        ### OPTIONAL 
        # pivot point will be set only for moves where a chain pivot occurs
        # and defines the CHAIN INDEX where a pivot occurs
        self.pivot_point = pivot_point


    def get_angle_indice(self, chain_length):
        """
        Returns the set of indices which correspond to the chain positions that have been changed by the move.
        Note that (for future reference) the ASSUMPTION here is that ALL returned indices are contigous in the chain. 
        This means that (for example) moves that cause multiple 'bends' in the chain (such as the slither move) end up
        returning a long list where MOST of the positions don't change. This is somewhat inefficient, BUT if (in the future)
        we implement residue-specific angle potentials this would be necessary, so for now this is being left as is. 
        
        """


        return_indices = []

        # for rigid body moves no angles
        # need evaluating
        if self.move_type in [2,3,7,8]:
            pass
            

        #  head pivot
        elif self.move_type == 5:
            moved_position = self.moved_indices[0]            
            for offset in range(-2,3):
                pos = moved_position + offset
                if pos >= 0 and pos < chain_length:
                    return_indices.append(pos)

        # chain pivot
        elif self.move_type == 4:
            moved_position = self.pivot_point
            
            for offset in range(-2,3):
                pos = moved_position + offset
                if pos >= 0 and pos < chain_length:
                    return_indices.append(pos)

        # chain slither
        elif self.move_type == 6:
            return_indices = list(range(0, chain_length))

        else:
            raise MoveException('Trying to get angle indices for move-code %i but no such angle index fetching implemented! Please report this bug' % self.move_type)



        return return_indices
            
                
                

            
            
            



        
        
        

