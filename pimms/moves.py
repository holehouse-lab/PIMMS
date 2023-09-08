## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2023
## ...........................................................................


import random
import numpy as np
import copy
import sys
import copy 

from . import lattice_utils
from . import numpy_utils
from . import CONFIG
from . import mega_crank
from . import mega_crank_2D
from . import crankshaft_list_functions
from . import IO_utils

from .latticeExceptions import MoveException, ClusterSizeThresholdException
from .moveset import MoveSet
from .moveEvent import MoveEvent

## A note on single chain MC moves (cluster moves are fundementally different...)
## MoveCodes 2 3 4 5 and 6 
##
## All single chain MC moves defined here must fulfill the following properties
##
## 1) Input should be the current 2D or 3D lattice grid and the Chain object which is to
##    be moved. Note that we're passing a Grid [e.g. an np.array]! NOT a Lattice object.
## 
## 2) The Chain object should be treated as READ-ONLY. It is not actually READ ONLY 
##    - it *can* be augmented but it should not be. This object is also stored in 
##    the calling Simulation object's Lattice object, where it is assumed to remain 
##    the same through a move.
## 
## 3) The lattice grid should be altered IF a move can be made, and should be 
##    returned as it came if the move cannot. 
##
## 4) Moves defined here DO NOT evaluate energy but DO evaluate hardsphere clashes
##    - i.e. a move should be rejected if after the move happening we find a site is 
##    occupied by two residues
##
## 5) Moves return two variables: 1) a MoveEvent object which holds all the information needed for 
##    further move acceptance/rejection in a well defined structures 2) True or False to
##    declare if the move should be evaluated or not
##
##
## There are three other types of moves which should be written up
##
## 1) Cluster moves, which are fundementally different but do need to be evaluated by
##    functions in Simulation.py (i.e. energy evaluation is done here)
##
## 2) Moves which take avantage of the mega_crank subchain moves - these moves keep 
##    track of their state and energy and so do not need to be subsequently re-
##    evaluated. These moves are also insanely fast - example being system_shake
##
## 3) TSMMC moves - there are three classes of TSMMC moves, one where a single
##    chain is perturbed, one where a random set of chains are perturbed, and 
##    one where the ENTIRE system is perturbed. For the single chain and set 
##    of chains all moves are crankshaft-like so energy evaluation is done in
##    hyperloop. For the system wide one it basically shifts the entire simulation
##    engine into an auxillary chain so subsequent moves are done as normal but 
##    but at a different temperature, and at the end the complete series of changes
##    are accepted or rejected. 

 
## General 'GOTCHAS'
## Despite the fact we set the Cython randmax initially it gets re-set each time
## the Cython code is called. As a result, any Cython function which uses random
## numbers should PASS a randomly-generated seed value in. If the Python random
## seed was set, this passed_seed will ensure reproducible behaviour. If not it
## doesn't hurt.
##
##
##
##
#


#-----------------------------------------------------------------
#    
class MoveObject:

    def __init__(self):
        """
        MoveObject is a stateless class, who's objects implement chain movement
        functionality but do not have any state associated with themselves

        """
        pass


            
    #-----------------------------------------------------------------
    #    
    #
    def system_shake(self, latticeObject, current_energy, acceptanceObject, hamiltonianObject, number_of_steps, mode, hardwall=False):
        """
        The system_shake move peforms a large number of very local chain perturbations. Each pertubration involves randomly selecting
        any bead, ensuring that complete detailed balance is maintained.


        latticeObject (Lattice object)
        latticeObject (as you might expect) the full lattice Object upon which the simulation is being performed. 

        curren_energy (int)
        Current energy value

        acceptanceObject (AcceptanceCalculator object)
        Contains all necessary details to accept or reject a move

        hamiltonianObject (Hamiltonain objec)
        Self contained object that allows for the evaluation of energy functions, and contains the interaction tables which
        can be passed to external (Cython) code for energy evaluation

        number_of_steps (int)
        Number of Monte Carlo moves to be performed on the single chain.

        mode (string)
        Defines the mode to be used for determining the final number of steps to be used. Currently obselete but kept in case
        we want to change how bead selection is done in the future.
        
        hardwall (bool) {False}
        Sets if a hardwall boundary is to be used. If false, periodic boundary conditions are used, but if true a hard-wall
        that is made of solvent but cannot be penetrated is used.

        crankshaft_prob (float) {0.0}
        If 1.0 means 100% of moves are crankshaft moves and we don't need to accept reject individual system shake
        moves, else we do


        """
        
        # construct the idx_to_bead matrix. which gets passed into megacrank. This matrix contains position and identity information
        # for all beads on the lattice
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead(latticeObject)

        # get the number of beads, build a selection vector, and figure out the number of dimensions
        num_beads = len(idx_to_bead)        
        num_dims  = len(latticeObject.dimensions)
        
        total_accepted = 0
        total_proposed = 0 

        # set hardwall flag
        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0

        local_seed = random.randint(1,sys.maxsize-1) % CONFIG.C_RAND_MAX

        bead_selector = np.random.randint(0,num_beads,number_of_steps)

        ##
        ## Both functions alter alter the grids on the back end and do not explicity
        ## reassign these as they're passed by reference as memoryviews (direct access to
        ## the memory)
        ## 

        
        # 2D
        if num_dims == 2:
            (new_energy, accepted_moves)= mega_crank_2D.mega_crank_2D(latticeObject.grid, 
                                                                      latticeObject.type_grid, 
                                                                      idx_to_bead,
                                                                      hamiltonianObject.residue_interaction_table,
                                                                      hamiltonianObject.LR_residue_interaction_table,
                                                                      hamiltonianObject.SLR_residue_interaction_table,
                                                                      hamiltonianObject.angle_lookup,
                                                                      current_energy,
                                                                      acceptanceObject.invtemp,
                                                                      number_of_steps,
                                                                      bead_selector,
                                                                      local_seed,
                                                                      hardwall_int)
                
        else:
            (new_energy, accepted_moves) = mega_crank.mega_crank(latticeObject.grid, 
                                                                 latticeObject.type_grid, 
                                                                 idx_to_bead,
                                                                 hamiltonianObject.residue_interaction_table,
                                                                 hamiltonianObject.LR_residue_interaction_table,
                                                                 hamiltonianObject.SLR_residue_interaction_table,
                                                                 hamiltonianObject.angle_lookup,
                                                                 current_energy,
                                                                 acceptanceObject.invtemp,
                                                                 number_of_steps,
                                                                 bead_selector,
                                                                 local_seed,
                                                                 hardwall_int)

            





        total_accepted = total_accepted + accepted_moves
        total_proposed = total_proposed + number_of_steps
        
        local_idx=0
        for chainID in sorted(latticeObject.chains.keys()):
            n_pos = len(latticeObject.chains[chainID].get_ordered_positions())

            latticeObject.chains[chainID].set_ordered_positions(idx_to_bead[local_idx:local_idx+n_pos,5:].tolist())
            local_idx=local_idx+n_pos

        current_energy = new_energy 

        return (latticeObject, current_energy, total_proposed, total_accepted)


    
    #-----------------------------------------------------------------
    #    
    def chain_translate(self, ChainToMove, lattice, hardwall=False):
        """
        The chain_translate move allows the full chain to be translated in rigid body
        space around the lattice.

        The cost of this move will scale linearly with chain length (note cost comes 
        from the energy evaluation). However this will often be one of the most expensive
        moves which can be made as we have to evaluate both short range and (if relevant)
        long range interactions for EVERY residue in the chain when re-calculating the 
        new energy.
    
        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        MoveType code: 2
        
        """

        chainID         = ChainToMove.chainID
        chain_positions = ChainToMove.get_ordered_positions()
        dimensions      = lattice_utils.get_dimensions(lattice)
        num_dims        = len(dimensions)

        # delete the chain from the lattice
        lattice_utils.delete_chain_by_position(chain_positions, lattice, chainID)

        # define translation operations through a translation vector which
        # we apply to each chain unit
        offset_vector = []
        for i in range(0, num_dims):
            offset_vector.append(random.randint(0, dimensions[i]-1))

        translated_positions = []

        # for the position of each residue
        for position in chain_positions:
            translated_pos = []

            # for each dimension incremement the position
            for dim in range(0, num_dims):
                translated_pos.append(position[dim] + offset_vector[dim] )

            # carry out periodic boundary correction on the new position
            translated_pos = lattice_utils.pbc_convert(translated_pos, dimensions)

            # check if that position is empty - as soon as we find a position
            # in the lattice which is not empty consider the move rejected.
            # This means we *only* carry out as many transation operations as
            # absolutely necessary
            if not lattice_utils.get_gridvalue(translated_pos, lattice) == 0:

                ## REVERT BACK !!
                
                # Delete the positions we insterted so far
                lattice_utils.delete_chain_by_position(translated_positions, lattice, chainID)

                # re-insert the chain back into its old positions
                lattice_utils.place_chain_by_position(chain_positions, lattice, chainID, safe=True)
                                
                return (False, False)
                #return (False, False, False, False, False)

            # if the position was free update the lattice copy object
            lattice_utils.set_gridvalue(translated_pos, chainID, lattice)
            
            # and add the position to the growing list
            translated_positions.append(translated_pos)

        # If we're outside the for-loop translation operation was a success!

        # Now check for hardwall rules
        if hardwall:
            if lattice_utils.do_positions_stradle_pbc_boundary(translated_positions):
                
                lattice_utils.delete_chain_by_position(translated_positions, lattice, chainID)
                lattice_utils.place_chain_by_position(chain_positions, lattice, chainID, safe=True)                                
                return (False, False)
                

        ## Create the MoveEvent object
        # We moved all the residues so moved_positions and full_moved_chain_positions
        # are the same.        
        ME = MoveEvent(original_positions        = chain_positions,
                       moved_positions           = translated_positions,
                       original_chain_positions  = chain_positions,
                       moved_chain_positions     = translated_positions,                
                       moved_indices             = list(range(0,len(chain_positions))),
                       move_type                 = 2)
                               
        return (ME, True)


    

    #-----------------------------------------------------------------
    #    
    def chain_rotate(self, ChainToMove, lattice, hardwall=False):
        """
        The chain_rotate move allows the full chain to be rotate in rigid body
        space around the chain's center of mass (i.e. minimizing translational
        movement). This is currently achieved by first determining the chain's COM,
        translating the whole chain to the origin, rotating the chain positions 
        about the origin, and translating BACK to its original center of mass. This
        is not the most efficienct implementation and may be re-written in the 
        future but now it works and is in no way a bottle neck to performance.
              
        The cost of this move will scale linearly with chain length (note cost comes 
        from the energy evaluation).

        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        MoveType code: 3

        """
        ## A note on rotations and offset. The offset parameter is calculated here
        ## so the chain can be first converted into a single image and then rotated
        ## (rather than rotating something that exists in periodic image space). This
        ## is actually irrelevant if you're using a box or square and the simulation
        ## vessel, but when vertices are unequal in length PBC conditions break and
        ## so this ensures we can still rotate chains in rectangular boxes


        chainID                  = ChainToMove.chainID
        chain_positions_original = ChainToMove.get_ordered_positions()
        chain_positions          = ChainToMove.get_single_image_positions()
        dimensions               = lattice_utils.get_dimensions(lattice)
        num_dims                 = len(dimensions)


        # delete the chain from the lattice
        lattice_utils.delete_chain_by_position(chain_positions_original, lattice, chainID)

        # 1) Dermine on-lattice center of mass of chain
        COM = ChainToMove.get_center_of_mass()

        # 2) translate each position to the origin
        OC_positions      = []
        rotated_positions = []
        
        # OC_ is 'origin centered' 

        ## 2D rotation 
        if num_dims == 2:

            # see start of file to explain offset
            offset = [chain_positions[0][0] - chain_positions_original[0][0],chain_positions[0][1] - chain_positions_original[0][1]]

            # translate to the origin
            for position in chain_positions:
                OC_positions.append([position[0] - COM[0], position[1] - COM[1]])

            # carry out a random rotation in 2D along a cardinal angle (note we will probably have to update this to a
            # set of discrete intervals to avoid rigid bodies being stuck in one of four rotational states)
            OC_rotated_positions = lattice_utils.rotate_positions_2D(OC_positions, [90,180,270][random.randint(0,2)])
            
            # now translate all the positions back from the origin
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert([position[0] + COM[0]+offset[0], position[1] + COM[1]+offset[1]], dimensions))
                
        ## 3D rotation
        if num_dims == 3:

            # see start of file to explain offset
            offset = [chain_positions[0][0] - chain_positions_original[0][0], chain_positions[0][1] - chain_positions_original[0][1], chain_positions[0][2] - chain_positions_original[0][2]]

            # translate to origin
            for position in chain_positions:
                OC_positions.append([position[0] - COM[0], position[1] - COM[1], position[2] - COM[2]])

            # carry out a random rotation in 3D
            OC_rotated_positions = lattice_utils.rotate_positions_3D(OC_positions, ['x','y','z'][random.randint(0,2)], [90,180,270][random.randint(0,2)])
            #OC_rotated_positions = lattice_utils.rotate_positions_3D(OC_positions, 'x', 90)
            
            # translate back from origin
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert([position[0] + COM[0] + offset[0] , position[1] + COM[1] + offset[1], position[2] + COM[2] + offset[2]], dimensions))
        
        # Now check for hardwall rules
        if hardwall:
            if lattice_utils.do_positions_stradle_pbc_boundary(rotated_positions):                

                # no need to delete anything because nothing inserted yet
                
                lattice_utils.place_chain_by_position(chain_positions_original, lattice, chainID, safe=True)                                
                return (False, False)
            else:
                pass


        # having built a new list of rotated positions let's see if any of them clash. Note the inserted_chain
        # keeps track of what's going on so if we find a clash we only have to cycle over a small number of filled
        # positions to delete the part of the chain we inserted
        inserted_chain = []
        for rotated_pos in rotated_positions:

            # if it turns out a position was already occupied
            if not lattice_utils.get_gridvalue(rotated_pos, lattice) == 0:

                # delete whatever part(s) of the chain we've already added
                lattice_utils.delete_chain_by_position(inserted_chain, lattice, chainID)

                # re-insert the chain back where it was
                lattice_utils.place_chain_by_position(chain_positions_original, lattice, chainID, safe=True)

                # return all the failure
                return (False, False)
            else:
                
                # if the position was free update the lattice copy object
                lattice_utils.set_gridvalue(rotated_pos, chainID, lattice)
                inserted_chain.append(rotated_pos)
            
        # if we get here we have succesfully added all the rotated positions to the lattice.
        # Assume all positions moved (maybe they didn't but determining this ends up being 
        # more computationally expensive
        ME = MoveEvent(original_positions        = chain_positions_original,
                       moved_positions           = rotated_positions,
                       original_chain_positions  = chain_positions_original,
                       moved_chain_positions     = rotated_positions,
                       moved_indices             = range(0,len(chain_positions)),
                       move_type                 = 3)

        return (ME, True)
    


    #-----------------------------------------------------------------
    #        
    def chain_pivot(self, ChainToMove, lattice, pivotPoint_range=None, hardwall=False):
        """
        The chain_pivot move allows part of the chain to perform a rigid
        pivot. Specifically, we randomly select (with uniform probability)
        some position through the chain, and then pivot the shorter half in
        a rigid-body manner using the selected position as an anchor point.

        NOTE: this is only applied to chains where are 3 residues or longer
        - the move is automatically rejected for shorter chains.
        
        The cost of this move will scale linearly with chain length (note cost 
        comes  from the energy evaluation).

        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        !!! WARNING !!!
        The pivotPoint_range allows for the user to define a range of positions
        which can be pivoted. THIS FUNCTIONALITY IS NOT GENERALIZABLE YET - it
        was implemented for a specific use case, and while generalizing it is
        not a major issue this has not yet been done - PLEASE DO NOT USE.

        MoveType code: 4

        """
    
        chainID         = ChainToMove.chainID
        chain_positions = ChainToMove.get_ordered_positions()
        chain_length    = len(chain_positions)
        dimensions      = lattice_utils.get_dimensions(lattice)
        num_dims        = len(dimensions)

        # reject if we have a chain 2 or less in length 
        if chain_length < 3:
            return (False, False)
            
        # select a position along the chain to pivot
        pivot_point = random.randint(1,len(chain_positions)-2)
                
        # if no pivot point range was provided (as default)
        if pivotPoint_range is None:

            # pivot the smaller of the two halves
            if pivot_point > chain_length/2:

                # set if we're going to to rotate positions and then
                # add them to the end of a fixed region
                add_to_end=True

                # positions which will be rotated (pivot point to end)
                positions_to_rotate    = chain_positions[pivot_point:]

                # positions which will be held fixed (0 to pivot point)
                positions_held_fixed   = chain_positions[:pivot_point]

                # sequence indices of positions which will be rotated
                indices = list(range(pivot_point, len(chain_positions)))

            else:
                # variable roles given above
                add_to_end=False
                positions_to_rotate  = chain_positions[:pivot_point]
                positions_held_fixed = chain_positions[pivot_point:]
                indices = list(range(0, pivot_point))

        else:
            # randomly select a position from the pivot point range - right now we're always
            # going to treat the pivot point as defining a region C-terminal which is pivoted while
            # the N-terminal region remains fixed
            pivot_point = pivotPoint_range[random.randint(0,len(pivotPoint_range)-1)]
           
            # set if we're going to to rotate positions and then
            # add them to the end of a fixed region
            add_to_end=True

            # positions which will be rotated (pivot point to end)
            positions_to_rotate    = chain_positions[pivot_point:]

            # positions which will be held fixed (0 to pivot point)
            positions_held_fixed   = chain_positions[:pivot_point]

            # sequence indices of positions which will be rotated
            indices = list(range(pivot_point, len(chain_positions)))
                              
        
        # delete the positions we're going to rotate but keep the rest
        lattice_utils.delete_chain_by_position(positions_to_rotate, lattice, chainID)

        # get the head position as a separate copy - this is the pivot position which should
        # remain fixed as we perform lever arm rotation on the rest of the residues in the
        # positions_to_rotate - note that depending on which side we're rotating we set
        # first or last residue as head - e.g. see diagram below
        #
        # - : residue reminaing fixe
        # p : pivot point
        # x : residue to pivot
        #
        # if add_to_end
        #       -...---------PXXXXX...X 
        # else
        #   XXX...XXXXXXXP------...---



        if add_to_end:
            head_position = positions_to_rotate[0][:]
        else:
            head_position = positions_to_rotate[-1][:]

        # this offset is, again, so we can use PBC boxes with non-equal sides (see intro in the rotate
        # _chain move
        original_positions_to_rotate = positions_to_rotate
        positions_to_rotate = lattice_utils.convert_chain_to_single_image(positions_to_rotate, dimensions)

        # translate the shorter halve's positions to the origin
        if num_dims == 2:
            x_ref = head_position[0]
            y_ref = head_position[1]
            
            # build a list of origin centered positions
            OC_positions = []
            for position in positions_to_rotate:
                OC_positions.append([position[0] - x_ref, position[1] - y_ref])

            # carry out a random rotation in 2D
            OC_rotated_positions = lattice_utils.rotate_positions_2D(OC_positions, [90,180,270][random.randint(0,2)])

            #
            # determine rotation offset - basically when we do the rotationn we want the residue which connects BACK to the chain to be in exactly
            # the same position - i.e.
            #
            #  XXP
            #    XXXX
            #
            #    PXX
            #    XXXX
            #
            # We have to make sure P (the pivot point) remians in the same position - if we're pivoting the first half the pivot point is at the
            # end of the OC_rotated_positions, while if we're pivoting the second half the pivot point is at the front of the OC_rotated_positions        
            #

            if add_to_end:
                x_return = x_ref - OC_rotated_positions[0][0]
                y_return = y_ref - OC_rotated_positions[0][1]
            else:
                x_return = x_ref - OC_rotated_positions[-1][0]
                y_return = y_ref - OC_rotated_positions[-1][1]
            
            # now move all the positions back again from the origin such that they move back to link up with the chain
            rotated_positions = []
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert([position[0] + x_return, position[1] + y_return], dimensions))

        ## 3D rotation
        if num_dims == 3:
            x_ref = head_position[0]
            y_ref = head_position[1]
            z_ref = head_position[2]

            # center the pivot section at the origin
            OC_positions = lattice_utils
            
            OC_positions = []
            for position in positions_to_rotate:
                OC_positions.append([position[0] - x_ref, position[1] - y_ref, position[2] - z_ref])
                
            # carry out a random rotation in 3D
            OC_rotated_positions = lattice_utils.rotate_positions_3D(OC_positions, ['x','y','z'][random.randint(0,2)], [90,180, 270][random.randint(0,2)])
            
            # determine rotation offset (sometimes rotation around the 0 axis will still move the head position - see 2D description
            # for more details!
            if add_to_end:
                return_correction = [x_ref - OC_rotated_positions[0][0], y_ref - OC_rotated_positions[0][1], z_ref - OC_rotated_positions[0][2]]
            else:
                return_correction = [x_ref - OC_rotated_positions[-1][0], y_ref - OC_rotated_positions[-1][1], z_ref - OC_rotated_positions[-1][2]]

            # now move all the positions back again from the origin
            rotated_positions = []
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert([position[0] + return_correction[0], position[1] + return_correction[1], position[2] + return_correction[2]], dimensions))
                


        # Now check for hardwall rules
        if hardwall:
            if lattice_utils.do_positions_stradle_pbc_boundary(rotated_positions):                
                # no need to delete anything because nothing inserted yet
                lattice_utils.place_chain_by_position(original_positions_to_rotate, lattice, chainID, safe=True)                                
                return (False, False)


        # having built a new list of rotated positions let's see if any of them clash. Note the inserted_chain
        # keeps track of what's going on so if we find a clash we only have to cycle over a small number of filled
        # positions to delete the part of the chain we inserted
        inserted_chain = []
        for rotated_pos in rotated_positions:

            # if it turns out a position was already occupied
            if not lattice_utils.get_gridvalue(rotated_pos, lattice) == 0:

                # delete whatever part(s) of the chain we've already added
                lattice_utils.delete_chain_by_position(inserted_chain, lattice, chainID)

                # re-insert the section of chain we deleted previously
                lattice_utils.place_chain_by_position(original_positions_to_rotate, lattice, chainID, safe=True)

                # return all the failure
                return (False, False)
            else:                
                # if the position was free update the lattice copy object
                lattice_utils.set_gridvalue(rotated_pos, chainID, lattice)
                inserted_chain.append(rotated_pos)
        
        # if we get here we succesfully pivoted the chain!
        if add_to_end:
            # add the newly rotated positions to the end of the already
            # known positions
            fully_pivoted_chain = positions_held_fixed + inserted_chain
        else:
            fully_pivoted_chain = inserted_chain + positions_held_fixed

        if not len(fully_pivoted_chain) == len(chain_positions):
            raise Exception('Yeah stop right there...')
            
        ME = MoveEvent(original_positions        = original_positions_to_rotate,
                       moved_positions           = inserted_chain,
                       original_chain_positions  = chain_positions,
                       moved_chain_positions     = fully_pivoted_chain,
                       moved_indices             = indices,
                       move_type                 = 4,
                       pivot_point               = pivot_point)
                   

        return (ME, True)                

    #-----------------------------------------------------------------
    #    
    def head_pivot(self, ChainToMove, lattice, hardwall=False):
        """
        The head_pivot move allows the chain 'head' (either the first or last residue)
        to pivot in some direction.

        This is always going to be very cheap as it 'always' translates to a single
        position change irrespective of chain length. We arbitrarily pick the first 
        or last residue to pivot - i.e. chains have 2 heads. 

        This is - honestly - kind of a stupid move. It was one of the first moves I
        coded up as a very simple and easy to debug move, but is probably not going
        to add much. However, the chain_pivot move won't pivot the ends of chains
        if you have a very short chain so it does actually serve a relevant purpose!

        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        MoveType code: 5
    
        """
        chainID           = ChainToMove.chainID
        chain_positions   = ChainToMove.get_ordered_positions()        
        dimensions        = lattice_utils.get_dimensions(lattice)
        num_dims          = len(dimensions)

        # copy because we want to create a new list of positions
        updated_positions = chain_positions[:] 

        # select one end of the chain as the head
        if random.random() > 0.5:
            
            ## **************************
            ## Working with first residue
        
            # Delete the first residue
            lattice_utils.delete_residue(chain_positions[0], lattice, chainID)

            # get possible new positions by building a list of the sites which are adajcent
            # to the second residue in the chain (having just deleted the first)
            if num_dims == 2:
                possible_positions = lattice_utils.get_adjacent_sites_2D(chain_positions[1][0], chain_positions[1][1], dimensions)
            else:
                possible_positions = lattice_utils.get_adjacent_sites_3D(chain_positions[1][0], chain_positions[1][1],chain_positions[1][2], dimensions)

            # randomly select one of the positions from this list
            possible_position = list(possible_positions[random.randint(0, len(possible_positions)-1)])
                

            # if hardwall boundary
            if hardwall:          

                # do the first and second positions now cross a PBC boundary?
                if lattice_utils.do_positions_stradle_pbc_boundary([possible_position, chain_positions[1]]):                

                    # revert back to original position
                    lattice_utils.insert_residue(chain_positions[0], lattice, chainID)
                    return (False, False)



            # if 'moved' the residue to the same position the original residue came from...
            # no move (same thing) - re insert and return false 
            # NOTE: Philiosphically should this be a rejection or not? I really don't know...
            if lattice_utils.same_sites(possible_position, chain_positions[0]):

                # revert to original position
                lattice_utils.insert_residue(possible_position, lattice, chainID)
                return (False, False)
                                
            # if moved into occupied site then reject
            elif not lattice_utils.get_gridvalue(possible_position, lattice) == 0.0:
            
                # revert to original position
                lattice_utils.insert_residue(chain_positions[0], lattice, chainID)
                return (False, False)

            # else move is A-OK!
            else:               
                
                # update the lattice 
                lattice_utils.insert_residue(possible_position, lattice, chainID)
                
                # set the first residue to the new position in the copy list
                updated_positions[0] = possible_position                

                # indicies represent the first residue
                ME = MoveEvent(original_positions        = [chain_positions[0]],
                               moved_positions           = [possible_position],
                               original_chain_positions  = chain_positions,
                               moved_chain_positions     = updated_positions,
                               moved_indices             = [0],
                               move_type                 = 5)                       
        
                return (ME, True)

        else:
            ## Working with last residue            

            # Delete last residue
            lattice_utils.delete_residue(chain_positions[-1], lattice, chainID)

            # get possible new positions 
            if num_dims == 2:
                possible_positions = lattice_utils.get_adjacent_sites_2D(chain_positions[-2][0], chain_positions[-2][1], dimensions)
            else:
                possible_positions = lattice_utils.get_adjacent_sites_3D(chain_positions[-2][0], chain_positions[-2][1], chain_positions[-2][2], dimensions)

            # randomly select one of the positions from this list
            possible_position = list(possible_positions[random.randint(0, len(possible_positions)-1)])


            # if hardwall boundary
            if hardwall:          

                # do the second to last and last positions now cross a PBC boundary?
                if lattice_utils.do_positions_stradle_pbc_boundary([chain_positions[-2],possible_position]):                

                    # revert back to original position
                    lattice_utils.insert_residue(chain_positions[-1], lattice, chainID)
                    return (False, False)
                
            # if 'moved' the same position head came from...
            # no move (same thing) - re insert and return false
            if lattice_utils.same_sites(possible_position, chain_positions[-1]):

                # revert to original position
                lattice_utils.insert_residue(possible_position, lattice, chainID)
                return (False, False)
                
            # if moved into occupied site then reject
            elif not lattice_utils.get_gridvalue(possible_position, lattice) == 0.0:
            
                # revert to original position
                lattice_utils.insert_residue(chain_positions[-1], lattice, chainID)
                return (False, False)

            else:

                lattice_utils.insert_residue(possible_position, lattice, chainID)

                # set the last residue to the new position in the copy list
                updated_positions[-1] = possible_position

                # indices represent the terminal residue
                ME = MoveEvent(original_positions        = [chain_positions[-1]],
                               moved_positions           = [possible_position],
                               original_chain_positions  = chain_positions,
                               moved_chain_positions     = updated_positions,                               
                               moved_indices             = [len(updated_positions)-1],
                               move_type                 = 5)
                       
        
                return (ME, True)





    #-----------------------------------------------------------------
    #    
    def chain_slither(self, ChainToMove, lattice, hardwall=False):
        """
        Slither the chain along in some (random) direction. Note that unlike other moves we include
        the homopolymer keyword here, although it is not currently implemented. However, in theory
        this makes the energy calculation MUCH cheaper for a homo polymer because you're basically
        just moving the bead at one end to the other end. However, a high performance of this
        is not implemented (it was in an early version of PIMMS) so for right now although
        the slither MOVE is cheap the associated energy calculation requires that EVERY bead
        in the chain is re-evaluated for its new energy.
        
        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        MoveType code: 6 
        
        """        
            
        chainID           = ChainToMove.chainID
        chain_positions   = ChainToMove.get_ordered_positions()[:] # copy of chain positions         
        dimensions        = lattice_utils.get_dimensions(lattice)
        num_dims          = len(dimensions)


        # copy because we want to create a new list of positions


        # select end to define as the 'head'
        if random.random() > 0.5:

            ## Working with first residue
        
            # get possible new positions by building a list of the sites which are adajcent
            # to the second residue in the chain (having just deleted the first)
            if num_dims == 2:
                possible_positions = lattice_utils.get_adjacent_sites_2D(chain_positions[0][0], chain_positions[0][1], dimensions)
            else:
                possible_positions = lattice_utils.get_adjacent_sites_3D(chain_positions[0][0], chain_positions[0][1],chain_positions[0][2], dimensions)
                
            # randomly select one of the positions from this list
            possible_position = list(possible_positions[random.randint(0, len(possible_positions)-1)])

            # if hardwall ask if this new position and the current first residue would now cross
            # a boundary
            if hardwall:
                if lattice_utils.do_positions_stradle_pbc_boundary([possible_position, chain_positions[0]]):                
                    return (False, False)

            # if that position is occupied reject
            if not lattice_utils.get_gridvalue(possible_position, lattice) == 0.0:            
                return (False, False)

            else:                

                # create a new list of positions 
                new_chain_positions = []
                new_chain_positions.append(possible_position)
                new_chain_positions.extend(chain_positions[0:-1])

                # update the lattice 
                # delete chain from the lattice
                lattice_utils.delete_chain_by_position(chain_positions, lattice, chainID)

                # insert chain into new position
                lattice_utils.place_chain_by_position(new_chain_positions, lattice, chainID, safe=True)
                
                
                ME = MoveEvent(original_positions        = chain_positions,
                               moved_positions           = new_chain_positions,
                               original_chain_positions  = chain_positions,
                               moved_chain_positions     = new_chain_positions,
                               moved_indices             = list(range(0, len(chain_positions))),
                               move_type                 = 6)

                return (ME, True)                                    
        else:

            # get possible new positions based on positions adjacent to terminal residue
            if num_dims == 2:
                possible_positions = lattice_utils.get_adjacent_sites_2D(chain_positions[-1][0], chain_positions[-1][1], dimensions)
            else:
                possible_positions = lattice_utils.get_adjacent_sites_3D(chain_positions[-1][0], chain_positions[-1][1], chain_positions[-1][2], dimensions)

            # randomly select one of the positions from this list (note case to list)
            possible_position = list(possible_positions[random.randint(0, len(possible_positions)-1)])

            
            # if hardwall ask if this new position and the current last residue would now cross
            # a boundary
            if hardwall:
                if lattice_utils.do_positions_stradle_pbc_boundary([chain_positions[-1],possible_position]):                
                    return (False, False)

                                
            # if moved into occupied site then reject
            if not lattice_utils.get_gridvalue(possible_position, lattice) == 0.0:            
                return (False, False)

            else:

                # create a new list of positions 
                new_chain_positions = []
                new_chain_positions.extend(chain_positions[1:])
                new_chain_positions.append(possible_position)

                # update the lattice 
                # delete chain from the lattice
                lattice_utils.delete_chain_by_position(chain_positions, lattice, chainID)

                # insert chain into new position
                lattice_utils.place_chain_by_position(new_chain_positions, lattice, chainID, safe=True)
                
                
                ME = MoveEvent(original_positions        = chain_positions,
                               moved_positions           = new_chain_positions,
                               original_chain_positions  = chain_positions,
                               moved_chain_positions     = new_chain_positions,
                               moved_indices             = list(range(0, len(chain_positions))),
                               move_type                 = 6)




                return (ME, True)                



    #-----------------------------------------------------------------
    #    
    def cluster_translate(self, selected_chain, latticeObject, cluster_move_threshold=None, cluster_size_threshold=None, hardwall=False):
        """
        The cluster_translate move allows a connected components (cluster) to be 
        translated in rigid body space around the lattice.

        The cost of this move grows linearly as clusters get big - the 
        cluster_threshold variable facilitates a soft threshold on the max cluster 
        size.

        To ensure detailed balance moves a cluster move MUST NOT lead to the 
        incorporation of new chains into the cluster being moved. Explicitly (thank you 
        Tyler!), if a cluster translation leads to two clusters merging, there is 
        no move in our which would allow that cluster to unmerge again in a single 
        move - i.e. a cluster-merging translation move is irreversible, hence breaking
        detailed balance.

        On the plus side, this also means that cluster moves which we can accept (i.e. 
        which does not lead to a cluster merger or clash) must be energy neutral, 
        so we don't have to run any short-range energy evaluations on it. Note that we
        must still run long-range energy calculations, meaning for charged systems 
        this becomes particularly expensive...

        The move is rejected if there's a hard-sphere clash, *or* we change the cluster
        size after the move (i.e. incorporate new residues in). If not rejected  we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        The cluster_size_threshold is default to None, but in the simulations.py file we
        set this to be such that in the case that a cluster contains ALL the chains it is not
        rotated, otherwise it can be.
        
        MoveType code: 7
        """

        original_chainID  = selected_chain.chainID
        dimensions        = latticeObject.dimensions
        num_dims          = len(dimensions)
        lattice           = latticeObject.grid 
        
        # note that get_all_chains_in_connected_component returns all the chainIDs in the connected
        # component chainID is part of *including* chainID!                                
        try:
            list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                       lattice, 
                                                                                       latticeObject.chains, 
                                                                                       threshold=cluster_size_threshold,
                                                                                       useChains=True, 
                                                                                       hardwall=hardwall)

        # this 'exception' occurs if we're scanning a connecte component and discover it's larger than the cluster_threshold 
        # note this isn't really an exception, but lets us implement a size-threshold in an interrupt-style manner (which is always
        # going to be maximally efficient) - note we MAY move a cluster larger than the threshold IFF we've already found the components - i.e.
        # you cannot assume that the largest cluster moved is equal to the threshold (this should probably be explicitly documented because
        # it's not very intuitive BUT makes everything a lot more efficient)
        except ClusterSizeThresholdException:
            return (False, False)

        # these dictionaries hold chainID indexed list of positions associated with a chain in their original
        # and new position
        old_chain_positions = {}
        new_chain_positions = {}

        # determine what the translation operation is gonna be
        offset_vector = []

        # cluster move threshold allows us to define translational movement as occuring in a maximum stepsize in each
        # dimension
        # if not included than we randomly move some distance
    
        if cluster_move_threshold is None:
            for i in range(0, num_dims):
                offset_vector.append(numpy_utils.randneg(random.randint(1, dimensions[i]-1)))
        else:
            for i in range(0, num_dims):
                offset_vector.append(numpy_utils.randneg(random.randint(1, min(dimensions[i]-1, cluster_move_threshold))))
            
        # now cycle through each chain in the connected commponent
        for chainID in list_of_chains_in_CC:

            # first get the chain's original position and then delete that chain from the lattice
            old_chain_positions[chainID] = latticeObject.chains[chainID].get_ordered_positions()
            lattice_utils.delete_chain_by_position(old_chain_positions[chainID], lattice, chainID)

            # move to its new position
            translated_positions = []
            for position in old_chain_positions[chainID]:

                translated_pos = []
                
                # determine the translated position and apply PBC corretions
                for dim in range(0, num_dims):
                    translated_pos.append(position[dim] + offset_vector[dim] )
                translated_pos = lattice_utils.pbc_convert(translated_pos, dimensions)
            
                # if the proposed position is already occupied back the f*ck up
                if not lattice_utils.get_gridvalue(translated_pos, lattice) == 0:
                
                    # Delete the positions we insterted so far in the *current* chain
                    lattice_utils.delete_chain_by_position(translated_positions, lattice, chainID)

                    # re-insert the full version of this chain
                    lattice_utils.place_chain_by_position(old_chain_positions[chainID], lattice, chainID, safe=True)

                    # For any chains we fully moved... we have to remove ALL the chains
                    # then reinsert them - this is because we *COULD* have carried out an operation where one chain
                    # moved into a space occupied by another chain so need to revert by fully deleting!
                    chains_reinserted = list(new_chain_positions.keys())
                    for chainIDs_inserted in chains_reinserted:
                        lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted)

                    for chainIDs_inserted in chains_reinserted:
                        lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted, safe=True)

                    # reject the move!
                    return (False, False)

                # if the position was free update the lattice grid object and add the position to
                # the growing list of new positions for this chain
                lattice_utils.set_gridvalue(translated_pos, chainID, lattice)
                translated_positions.append(translated_pos)            

            # if we get here we succesfully inserted an entire chain into the grid so save 
            # the translated_positions as the chain's positions, BUT FIRST check for hardwall rules and reject 
            # the move IF we're applying a hardwall boundary and the chain breaks that hardwall
            if hardwall:

                if lattice_utils.do_positions_stradle_pbc_boundary(translated_positions):
                    
                    # this is exactly the same protocol as we use to reject the move in the case of the clash above, just not annotated
                    # as heavily...
                    lattice_utils.delete_chain_by_position(translated_positions, lattice, chainID)                    
                    lattice_utils.place_chain_by_position(old_chain_positions[chainID], lattice, chainID, safe=True)

                    chains_reinserted = list(new_chain_positions.keys())
                    for chainIDs_inserted in chains_reinserted:
                        lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted)

                    for chainIDs_inserted in chains_reinserted:
                        lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted, safe=True)

                    # reject the move!
                    return (False, False)

            new_chain_positions[chainID] = translated_positions
            
        # >>>>
        # if we get here we moved the cluster! However we have to determine if the NEW cluster position is also 
        # the same size connected component - if it's larger this move would break detailed balance, 
        # if it's the same this move is fine but it is (by definition) energy neutral so no need to do 
        # short range energy calculations (need long-range ones though!)

        size_of_original_cluster = len(list_of_chains_in_CC)
        
        # we now build a new, bespoke positions dictionary which contains the new positions the moved chains
        # and the original positions of the chains which haven't moved (i.e. just a lattice up-to-date list
        # of chain positions
        chainPositionDict={}
        for i in range(1, len(latticeObject.chains)+1):            

            # recall chainID 1 is at position 0 in the chains list
            # and so on...
            if i in new_chain_positions:                            
                chainPositionDict[i] = new_chain_positions[i]
            else:
                chainPositionDict[i] = latticeObject.chains[i].get_ordered_positions()
        
        try:
            new_list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                           lattice, 
                                                                                           chainPositionDict, 
                                                                                           threshold=size_of_original_cluster,
                                                                                           useChains=False,
                                                                                           hardwall=hardwall)


            # finally make sure that if we re-selected this chain and got a list of chains it's THE SAME list of chains!
            # This is actually important - same size is too lenient, as you could move across a pbc but keep number of chain
            # fixed - the implementation below is a necessary and sufficient check to ensure that:

            # a) All the new chain IDs were in the list of old IDs
            for new_id in new_list_of_chains_in_CC:
                if new_id not in list_of_chains_in_CC:
                    raise ClusterSizeThresholdException

            # b) All the old chain IDs are in the new of new IDs
            for old_id in list_of_chains_in_CC:
                if old_id not in new_list_of_chains_in_CC:
                    raise ClusterSizeThresholdException
               

        # if we find adding more chains than we had before in the CC then an exception is raised and we know our 
        # cluster move caused cluster merging        
        # ****************************************************************************************************
        except ClusterSizeThresholdException:

            # revert back by deleting the chains we insterted and then re-setting the old chain
            chains_reinserted = list(new_chain_positions.keys())
            for chainIDs_inserted in chains_reinserted:
                lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted)

            for chainIDs_inserted in chains_reinserted:
                lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted, safe=True)

            return (False, False)
        # ****************************************************************************************************

        # if we get here move is a go!
        ME = MoveEvent(original_positions        = old_chain_positions,
                       moved_positions           = new_chain_positions,
                       original_chain_positions  = old_chain_positions,
                       moved_chain_positions     = new_chain_positions,
                       moved_indices             = None,
                       move_type                 = 7)
        
        return (ME, True)                



    #-----------------------------------------------------------------
    #    
    def cluster_rotate(self, selected_chain, latticeObject, cluster_move_threshold=None, cluster_size_threshold=None, hardwall=False):
        """
        The cluster_rotate move allows a connected components (cluster) to be 
        rotated in rigid body space around the lattice. Right now rotation occurs only
        over the cardinal directions (0/90/180/270) because anything other than this is 
        hard on a lattice...

        The cost of this move becomes massive as clusters get big - the cluster_threshold 
        variable facilitates a soft threshold on the max cluster size. 

        The cluster_size_threshold is default to None, but in the simulations.py file we
        set this to be such that in the case that a cluster contains ALL the chains it is not
        rotated, otherwise it can be.
        
        The move is rejected if there's a hard-sphere clash, *or* we change the cluster
        size after the move (i.e. incorporate new residues in). If not rejected  we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.
        
        MoveType code: 8
        """

        original_chainID        = selected_chain.chainID
        dimensions              = latticeObject.dimensions
        num_dims                = len(dimensions)
        lattice                 = latticeObject.grid 

        old_chain_positions            = {}
        new_chain_positions_OC         = {}
        new_chain_positions_OC_rotated = {}
        new_chain_positions            = {}
        
        # note that get_all_chains_in_connected_component returns all the chainIDs in the connected
        # component chainID is part of *including* chainID!                
        try:
            list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                       lattice, 
                                                                                       latticeObject.chains, 
                                                                                       threshold=cluster_size_threshold,
                                                                                       useChains=True,
                                                                                       hardwall=hardwall)

        # this 'exception' occurs if we're scanning a connecte component and discover it's larger than the clutser_threshold 
        # note this isn't really an exception, but lets us implement a size-threshold in an interrupt-style manner (which is always
        # going to be maximally efficient) - note we MAY move a cluster larger than the threshold IFF we've already found the coponents - i.e.
        # you cannot assume that the largest cluster moved is equal to the threshold (this should probably be explicitly documented because
        # it's not very intuitive BUT makes everything a lot more efficient)

            
        except ClusterSizeThresholdException:
            return (False, False)

        # these dictionaries hold chainID indexed list of positions associated with a chain in their original
        # and new position - delete these chains from the lattice!
        all_cluster_positions =[]
        for chainID in list_of_chains_in_CC:

            all_cluster_positions.extend(latticeObject.chains[chainID].get_ordered_positions())
            old_chain_positions[chainID] = latticeObject.chains[chainID].get_ordered_positions()

            lattice_utils.delete_chain_by_position(old_chain_positions[chainID], lattice, chainID)


        
        # so now ALL the chains in the cluster have been deleted from the lattice
        # get the cluster center of mass based on those chains' positions
        COM = lattice_utils.center_of_mass_from_positions(all_cluster_positions, dimensions)

        # this runs the snakesearch algorithm on all the cluster components 
        #single_image_positions = cluster_utils.convert_positions_to_single_image_snakesearch(all_cluster_positions, dimensions)
        
        ## ----------------------------------------------------------------------------------------------------
        ## 2D CASE FIRST
        ##        
        if num_dims == 2:

            rotationFactor = [90,180,270][random.randint(0,2)]

            # now cycle through each chain in the connected commponent moving it such that it's centered on the
            # origin (OC = origin centered)
            for chainID in list_of_chains_in_CC:
                                             
                # move chain to origin
                new_chain_positions_OC[chainID] = []

      
                for position in old_chain_positions[chainID]:
                    new_chain_positions_OC[chainID].append([position[0] - COM[0], position[1] - COM[1]])

                # rotate 2D positions by the rotation operation defined
                new_chain_positions_OC_rotated = lattice_utils.rotate_positions_2D(new_chain_positions_OC[chainID], rotationFactor)
                
                #new_chain_positions_OC_rotated = lattice_utils.rotate_positions_2D(old_chain_positions_SIC, rotationFactor)
                
                # move back to original location
                new_chain_positions[chainID] = []
                
                for position in new_chain_positions_OC_rotated:
                    new_chain_positions[chainID].append(lattice_utils.pbc_convert([position[0] + COM[0], position[1] + COM[1]], dimensions))
                


        ## ----------------------------------------------------------------------------------------------------
        ## 3D CASE SECOND
        ##
        else:
            rotationFactor = [90,180,270][random.randint(0,2)]
            rotationDim    = ['x','y','z'][random.randint(0,2)]

            # now cycle through each chain in the connected commponent moving it such that it's centered on the
            # origin (OC = origin centered)
            for chainID in list_of_chains_in_CC:
                                                
                new_chain_positions_OC[chainID] = []
                
                for position in old_chain_positions[chainID]:
                    new_chain_positions_OC[chainID].append([position[0] - COM[0], position[1] - COM[1], position[2] - COM[2]])
                
                # rotate 3D positions by the rotation operation defined
                new_chain_positions_OC_rotated = lattice_utils.rotate_positions_3D(new_chain_positions_OC[chainID], rotationDim, rotationFactor)
                #new_chain_positions_OC_rotated[chainID] = lattice_utils.rotate_positions_3D(old_chain_positions_SIC, rotationDim, rotationFactor)

                # move back to original location
                new_chain_positions[chainID] = []
                
                for position in new_chain_positions_OC_rotated:
                    new_chain_positions[chainID].append(lattice_utils.pbc_convert([position[0] + COM[0], position[1] + COM[1], position[2] + COM[2]], dimensions))
                    
        
                            
        ## ----------------------------------------------------------------------------------------------------
        # having built a new list of rotated positions for each chain let's see if any of them clash. Note the inserted_chain
        # keeps track of what's going on so if we find a clash we only have to cycle over a small number of filled
        # positions to delete the part of the chain we inserted

        # for each position in each chain
        chains_reinserted = []
        for chainID in new_chain_positions:
            
            rotated_positions = []
            for position in new_chain_positions[chainID]:

                # if the position we're rotating into is CURRENTLY occupied 
                if not lattice_utils.get_gridvalue(position, lattice) == 0:

                    IO_utils.status_message("Rejection because of clash",'info')
                    
                    # Delete the positions we insterted so far in the *current* chain and then
                    # delete all the other chains which were fully rotated
                    lattice_utils.delete_chain_by_position(rotated_positions, lattice, chainID)
                    for chainIDs_rotated in chains_reinserted:
                        lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_rotated], lattice, chainIDs_rotated)

                    # now reinsert ALL the chains back....
                    for chainIDs_org in old_chain_positions:
                         lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_org], lattice, chainIDs_org, safe=True)

                    # reject the move!
                    return (False, False)

                # else we're OK
                rotated_positions.append(position)
                lattice_utils.set_gridvalue(position, chainID, lattice)
                                    

            if hardwall:
                if lattice_utils.do_positions_stradle_pbc_boundary(rotated_positions):
                    
                    # this is exactly the same protocol as we use to reject the move in the case of the clash above, just not annotated
                    # as heavily...
                    lattice_utils.delete_chain_by_position(rotated_positions, lattice, chainID)
                    for chainIDs_rotated in chains_reinserted:
                        lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_rotated], lattice, chainIDs_rotated)

                    # now reinsert ALL the chains back....
                    for chainIDs_org in old_chain_positions:
                         lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_org], lattice, chainIDs_org, safe=True)

                    # reject the move!
                    return (False, False)

            # if we get here inserted the whole chain, so add it to the list of [succesfully] re-inserted chains
            chains_reinserted.append(chainID)

        # if we get here we moved the cluster! However we have to determine if the NEW cluster position is also 
        # the same size connected component - if it's larger this move would break detailed balance, 
        # if it's the same this move is fine but it is (by definition) energy neutral so no need to do 
        # energy calculations

        size_of_original_cluster = len(list_of_chains_in_CC)
        
        # we now build a new, bespoke positions dictionary which contains the new positions of the moved chains
        # and the original positions of the chains which haven't moved (i.e. just a lattice up-to-date list
        # of chain positions
        chainPositionDict={}
        for i in range(1, len(latticeObject.chains)+1):            
            # recall chainID 1 is at position 0 in the chains list
            # and so on...
            if i in new_chain_positions:                            
                chainPositionDict[i] = new_chain_positions[i]
            else:
                chainPositionDict[i] = latticeObject.chains[i].get_ordered_positions()
        
        try:
            new_list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                           lattice, 
                                                                                           chainPositionDict, 
                                                                                           threshold=size_of_original_cluster,
                                                                                           useChains=False,
                                                                                           hardwall=hardwall)

            # finally make sure the new cluster isn't SMALLER (could happen in hardwall mode) 
            for new_id in new_list_of_chains_in_CC:
                if new_id not in list_of_chains_in_CC:
                    raise ClusterSizeThresholdException

            for old_id in list_of_chains_in_CC:
                if old_id not in new_list_of_chains_in_CC:
                    raise ClusterSizeThresholdException


        # if we find adding more chains than we had before in the CC then an exception is raised and we know our cluster move caused cluster
        # merging
        # ****************************************************************************************************
        except ClusterSizeThresholdException:

            IO_utils.status_message("Cluster resize rejection",'info')

            # revert back by deleting the chains we insterted and then re-setting the old chain
            chains_reinserted = list(new_chain_positions.keys())
            for chainIDs_inserted in chains_reinserted:
                lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted)

            for chainIDs_inserted in chains_reinserted:
                lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted, safe=True)

            return (False, False)
        # ****************************************************************************************************

        # if we get here move is a go!
        ME = MoveEvent(original_positions        = old_chain_positions,
                       moved_positions           = new_chain_positions,
                       original_chain_positions  = old_chain_positions,
                       moved_chain_positions     = new_chain_positions,
                       moved_indices             = None,
                       move_type                 = 8)
        
        return (ME, True)                



    #-----------------------------------------------------------------
    #    
    def Chain_based_TSMMC(self, chainID, latticeObject, current_energy, hamiltonianObject, CTSMMC, hardwall=False):
        """
        The chain-based Temperature Sweet Metropolis Monte Carlo move involves (for a single chain)
        slowly increasing and then decreasingthe temperature while  

        In terms of big picture - this move involves creating an alternative Monte Carlo chain. This
        alternative chain experiences an initial jump to a high temperature and then a gradual drop in
        temperature back to the simulation temperature. As the temperature drops from the jump temperature
        up to the high temperature a number of MC moves are performed at the interveneing temperatures.
        
        Once the chain has returned back to the original temperature a FINAL Metropolis accept/reject
        query is performed to ask if the chain in it's new position should be accepted or rejected.

        In this way, although through the move we perform a LOT of MC accept/reject moves they're mostly
        from a different Hamiltonian so we have to treat them all as some (smart) pertubation of the chain
        which gets evaluated at the end. 

        HOWEVER, we evaluate this move here and then don't through the standard single chain energy evaluation
        because throughout the actual move we keep track of the system energy so don't have to re-evaluate 
        after. 

        [1] Mittal, A., Lyle, N., Harmon, T.S., and Pappu, R.V. (2014). Hamiltonian Switch Metropolis Monte Carlo 
            Simulations for Improved Conformational Sampling of Intrinsically Disordered Regions Tethered to Ordered 
            Domains of Proteins. J. Chem. Theory Comput. 10, 3550-3562.

        [2] Gelb, L.D. (2003). Monte Carlo simulations using sampling from an approximate potential. J. Chem. Phys. 118, 7747-7750.

        
        MoveType code: 9


        """

        idx_to_bead = crankshaft_list_functions.update_idx_to_bead_single_chain(latticeObject, chainID)
        chain_length = len(idx_to_bead)
        
        # this is a copy because its a list (if this was a numpy array would be by reference and would
        # not be a copy) 
        original_chain_positions = copy.deepcopy(latticeObject.chains[chainID].get_ordered_positions())

        # save old energy
        old_energy = current_energy
        num_dims = len(latticeObject.dimensions)
        num_temps = len(CTSMMC.inv_temperature_schedule)
        
        steps_per_temperature = chain_length * CTSMMC.steps_per_quench_multiplier

        total_moves = steps_per_temperature * num_temps
            
        # these are passed by reference, but we set to them variables so we can iteratively pass them
        # at different temperatures
        #tmp_grid            = latticeObject.grid
        #tmp_type_grid       = latticeObject.type_grid
        new_energy = current_energy

        # set new energy to current energy - this will be updated sequentially as we proceed
        

        # set hardwall flag
        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0
        
        for temp_idx in range(0, num_temps):

            # set previous and current inverse temperatures
            inv_temp = CTSMMC.inv_temperature_schedule[temp_idx]
            local_seed = random.randint(1,sys.maxsize-1) % CONFIG.C_RAND_MAX

            bead_selector = np.random.randint(0, chain_length, steps_per_temperature)


            ##
            ## Both functions alter alter the grids on the back end and do not explicity
            ## reassign these as they're passed by reference as memoryviews (direct access to
            ## the memory)
            ## 
            
            if num_dims == 2:
                (new_energy, accepted_moves)= mega_crank_2D.mega_crank_2D(latticeObject.grid,
                                                                          latticeObject.type_grid,
                                                                          idx_to_bead,
                                                                          hamiltonianObject.residue_interaction_table,
                                                                          hamiltonianObject.LR_residue_interaction_table,
                                                                          hamiltonianObject.SLR_residue_interaction_table, 
                                                                          hamiltonianObject.angle_lookup,
                                                                          new_energy,
                                                                          inv_temp,
                                                                          steps_per_temperature,
                                                                          bead_selector,
                                                                          local_seed,
                                                                          hardwall_int)
                
            else:


                (new_energy, accepted_moves) = mega_crank.mega_crank(latticeObject.grid,
                                                                     latticeObject.type_grid,
                                                                     idx_to_bead,
                                                                     hamiltonianObject.residue_interaction_table,
                                                                     hamiltonianObject.LR_residue_interaction_table,
                                                                     hamiltonianObject.SLR_residue_interaction_table,
                                                                     hamiltonianObject.angle_lookup,
                                                                     new_energy,
                                                                     inv_temp,
                                                                     steps_per_temperature,
                                                                     bead_selector,
                                                                     local_seed,
                                                                     hardwall_int)

        # if move is accepted update the grids, the energy, and the chain positions
        if CTSMMC.accept_TSMMC(new_energy, old_energy, CTSMMC.inv_target_temperature, inv_temp):

            # udpate the chain positions
            current_energy = new_energy 

            # the beauty is this works for the 2D and 3D case
            latticeObject.chains[chainID].positions = idx_to_bead[:,5:].tolist()
            
            return (latticeObject, current_energy, total_moves, True)
            
        # reject the whole move
        else:
            
            # construct a new list of the chain's new positions based on the tmp_chain_positions matrix
            deletable_positions = idx_to_bead[:,5:].tolist()
            
            # revert the lattice to it's pre-move state 
            lattice_utils.delete_chain_by_position(deletable_positions, latticeObject.grid, chainID)                
            lattice_utils.place_chain_by_position(original_chain_positions, latticeObject.grid, chainID, safe=True)
            
            # set chain positions in the chain-list positions                        
            latticeObject.chains[chainID].set_ordered_positions(original_chain_positions)
                
            # update the type_grid variable BACK 
            latticeObject.update_type_grid(chainID, deletable_positions, original_chain_positions, list(range(0,len(original_chain_positions))), safe=True)
            
            # return everything!
            return (latticeObject, old_energy, total_moves, False)
                                    

                
    #-----------------------------------------------------------------
    #    
    def multichain_based_TSMMC(self, original_chainID, latticeObject, current_energy, hamiltonianObject, CTSMMC, hardwall=False):
        """
        Same idea as Chain_based_TSMMC except here we randomly select some number of chains (currently this is 
        defined by the max_number_selectable function, which is set at 25% of the total number of chains on the
        lattice.

        Then, we sequentially raise and then lower the temperature, and at each different temperature cycle through
        the chains and update their positions. At the end the full move is accepted or rejected. See the 
        chain_based_TSMMC write up for more details on what's actually going on in terms of the TSMMC-ness.
        
        MoveType code: 10
        """
                    
        dimensions      = latticeObject.dimensions
        num_dims        = len(dimensions)
        num_temps       = len(CTSMMC.inv_temperature_schedule)
        lattice         = latticeObject.grid 
        old_energy      = current_energy

        
        ## in the current implementation we randomly select between 1 and 25% of the chains in the system
        # First figure out what 25% of the number of chains is
        num_chains      = latticeObject.get_number_of_chains()
        all_chains      = list(latticeObject.chains.keys())

        # this works with 1 through n chains and give sensible values
        max_number_selectable = np.floor(0.25*num_chains)+1
        number_selectable     = random.randint(1,max_number_selectable)
        list_of_chains = np.random.choice(all_chains, number_selectable,replace=False) # DO NOT REPLACE!

        # list_of_chains is now a list of chain IDs that we're going to perturb

        # this dictionary allows us to map the chainID to the old positions so we can rever if needed
        all_original_positions = {}        

        # for each chain, save the original positions via a deepcopy operation
        for chainID in list_of_chains:            
            # save the original positions in case we have to revert
            all_original_positions[chainID] = copy.deepcopy(latticeObject.chains[chainID].get_ordered_positions())

        # construct a specific idx_to_bead matrix that reflects the beads taken from this list of chains in the order
        # they appear in the list_of_chains
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead_multiple_chains(latticeObject, list_of_chains)

        # total number of beads
        num_beads = idx_to_bead.shape[0]

        # calculate the number of steps per temperature
        steps_per_temperature = len(idx_to_bead)*CTSMMC.steps_per_quench_multiplier
                
        # total proposed moves (keep track for peformance analysis post-factor)
        total_moves = 0

        new_energy          = current_energy
        #tmp_grid            = latticeObject.grid
        #tmp_type_grid       = latticeObject.type_grid

        # set hardwall flag
        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0

            
        for temp_idx in range(0, num_temps):

            # set previous and current inverse temperatures
            inv_temp = CTSMMC.inv_temperature_schedule[temp_idx]
            local_seed = random.randint(1,sys.maxsize-1) % CONFIG.C_RAND_MAX


            bead_selector = np.random.randint(0, num_beads, steps_per_temperature)

            ##
            ## Both functions alter alter the grids on the back end and do not explicity
            ## reassign these as they're passed by reference as memoryviews (direct access to
            ## the memory)
            ## 

            if num_dims == 2:
                
                (new_energy, accepted_moves)= mega_crank_2D.mega_crank_2D(latticeObject.grid,
                                                                          latticeObject.type_grid,
                                                                          idx_to_bead,
                                                                          hamiltonianObject.residue_interaction_table,
                                                                          hamiltonianObject.LR_residue_interaction_table,
                                                                          hamiltonianObject.SLR_residue_interaction_table,
                                                                          hamiltonianObject.angle_lookup,
                                                                          new_energy,
                                                                          inv_temp,
                                                                          steps_per_temperature,
                                                                          bead_selector,
                                                                          local_seed,
                                                                          hardwall_int)
                
            else:

                (new_energy, accepted_moves) = mega_crank.mega_crank(latticeObject.grid,
                                                                     latticeObject.type_grid,
                                                                     idx_to_bead,
                                                                     hamiltonianObject.residue_interaction_table,
                                                                     hamiltonianObject.LR_residue_interaction_table,
                                                                     hamiltonianObject.SLR_residue_interaction_table,
                                                                     hamiltonianObject.angle_lookup,
                                                                     new_energy,
                                                                     inv_temp,
                                                                     steps_per_temperature,
                                                                     bead_selector,
                                                                     local_seed,
                                                                     hardwall_int)

            
        # if move was accepted
        if CTSMMC.accept_TSMMC(new_energy, old_energy, CTSMMC.inv_target_temperature, inv_temp):
            current_energy = new_energy 
            
            
            # cycle over each chain, and for each chain update the positions by exracting the updated
            # positions from the tmp_chain_positions matrix. The tmp_chain_poistions matrix contains ONLY
            # bead positions that were moved in a specific order that corresponds to the the order of beads
            # from the chains in the list_of_chains, so this ensures we update position correctly
            idx=0
            for chainID in list_of_chains:
                
                # chain length
                chain_len = len(latticeObject.chains[chainID].positions)

                latticeObject.chains[chainID].positions = idx_to_bead[idx:idx+chain_len,5:].tolist()
                idx = idx + chain_len

            IO_utils.status_message("Multichain re-arrangement accepted [dE = %i]  (number of chains: %i)" %(new_energy - old_energy, len(list_of_chains)))
            return (latticeObject, current_energy, total_moves, True)
                
        else:
            all_new_positions={}

            # same logic as was used to update the chain positions (see text in the move-success branch
            # for an explanation)
            idx=0
            for chainID in list_of_chains:
                
                # chain length
                chain_len = len(latticeObject.chains[chainID].positions)

                # UP-2023
                all_new_positions[chainID] = idx_to_bead[idx:idx+chain_len,5:].tolist()
                idx = idx + chain_len
                
            ## Having obtained the positions of all the moved beads...

            # delete all the chains off the lattice
            for chainID in list_of_chains:
                
                # delete from main grind
                lattice_utils.delete_chain_by_position(all_new_positions[chainID], latticeObject.grid, chainID)

                # delete from type grid
                latticeObject.delete_chain_from_type_grid(chainID, all_new_positions[chainID], list(range(0,len(all_new_positions[chainID]))), safe=True)
                    
            # re-insert the chains back into their original position
            for chainID in list_of_chains:
                original_chain_positions = all_original_positions[chainID]

                # insert into main grid
                lattice_utils.place_chain_by_position(original_chain_positions, latticeObject.grid, chainID, safe=True)

                # insert into type grid
                latticeObject.insert_chain_into_type_grid(chainID, original_chain_positions, list(range(0,len(original_chain_positions))), safe=True)
                                    
                # set chain positions in the chain-list positions
                latticeObject.chains[chainID].set_ordered_positions(original_chain_positions)

                # reset the energy

            current_energy = old_energy

            return (latticeObject, current_energy, total_moves, False)
           


    #-----------------------------------------------------------------
    #    
    def ratchet_pivot(self, chainID, latticeObject, current_energy, acceptanceObject, hamiltonianObject, hardwall=False):
        """
        Code: 11

        ratchet pivoy was an older move that is explicitly no longer supported. We keep this stub in place, and will likely
        replace this move with something else in the future. 
      
          
        """

        return (latticeObject, current_energy, 0, False)


               
    # System_based_TSMMC
    # CODE 12
    #
    # Not actually a move implemented here, but implemented in the simulation object with functionality in the TSMMC object
    # too. This is here mainly to ensure that the next move added uses CODE 13 (OOO unlucky!!!!! Sucks to be you!)

    #-----------------------------------------------------------------
    #    
    def single_chain_shake(self, chainID, latticeObject, current_energy, acceptanceObject, hamiltonianObject, number_of_steps, mode, hardwall):
        """
        
        latticeObject (Lattice object)
        latticeObject (as you might expect) the full lattice Object upon which the simulation is being performed. 

        curren_energy (int)
        Current energy value

        acceptanceObject (AcceptanceCalculator object)
        Contains all necessary details to accept or reject a move

        hamiltonianObject (Hamiltonain objec)
        Self contained object that allows for the evaluation of energy functions, and contains the interaction tables which
        can be passed to external (Cython) code for energy evaluation

        number_of_steps (int)
        Number of Monte Carlo moves to be performed on the single chain.

        mode (string)
        Defines the mode to be used for determining the final number of steps to be used. Currently obselete but kept in case
        we want to change how bead selection is done in the future.
        
        hardwall (bool) {False}
        Sets if a hardwall boundary is to be used. If false, periodic boundary conditions are used, but if true a hard-wall
        that is made of solvent but cannot be penetrated is used.


        """
        
        # get number of dimenisons and set various initial values
        num_dims = len(latticeObject.dimensions)
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead_single_chain(latticeObject, chainID)
        chain_length = len(idx_to_bead)
           
        # set some initial values, the hardwall flag, and set the randoms seed 
        total_accepted = 0
        total_proposed = 0 

        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0
            
        local_seed = random.randint(1,sys.maxsize-1) % CONFIG.C_RAND_MAX

        bead_selector = np.random.randint(0, chain_length, number_of_steps)

        ##
        ## Both functions alter alter the grids on the back end and do not explicity
        ## reassign these as they're passed by reference as memoryviews (direct access to
        ## the memory)
        ## 

        # 2D
        if num_dims == 2:
            (new_energy, accepted_moves) = mega_crank_2D.mega_crank_2D(latticeObject.grid, 
                                                                       latticeObject.type_grid, 
                                                                       idx_to_bead,
                                                                       hamiltonianObject.residue_interaction_table,
                                                                       hamiltonianObject.LR_residue_interaction_table,
                                                                       hamiltonianObject.SLR_residue_interaction_table,
                                                                       hamiltonianObject.angle_lookup,
                                                                       current_energy,
                                                                       acceptanceObject.invtemp,
                                                                       number_of_steps,
                                                                       bead_selector,
                                                                       local_seed,
                                                                       hardwall_int)
                
        else:
            (new_energy, accepted_moves) = mega_crank.mega_crank(latticeObject.grid, 
                                                                 latticeObject.type_grid, 
                                                                 idx_to_bead,
                                                                 hamiltonianObject.residue_interaction_table,
                                                                 hamiltonianObject.LR_residue_interaction_table,
                                                                 hamiltonianObject.SLR_residue_interaction_table,
                                                                 hamiltonianObject.angle_lookup,
                                                                 current_energy,
                                                                 acceptanceObject.invtemp,
                                                                 number_of_steps,
                                                                 bead_selector,
                                                                 local_seed,
                                                                 hardwall_int)

        total_accepted = total_accepted + accepted_moves
        total_proposed = total_proposed + number_of_steps


        # set energy
        current_energy = new_energy 

        latticeObject.chains[chainID].positions = idx_to_bead[:,5:].tolist()

        return (latticeObject, current_energy, total_proposed, total_accepted)


