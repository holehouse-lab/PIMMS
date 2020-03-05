## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

import numpy as np
import random

from .latticeExceptions import CustomInitializationException

from . import chain
from . import lattice
from . import lattice_utils

class NeurofilamentDemo:
    """


    """

    def __init__(self, Hamiltonian):
        """
        Generates a 500x500x500 box with a 20x20x500 central filemante
        and arbitrary spaced sidearms extending outwards


        """
        
        # every n-th layer in the Z direction a sidearm will be randomly
        # assigned along the filament axis projecting straight out in
        # a cardinal direction
        dimensions = [500,500,500]

        central_filament_type='E'
        sidearm_type = 'E'

        grid         = np.zeros(dimensions, dtype=int)
        type_grid    = np.zeros(dimensions, dtype=int)
        central_filament_positions = []

        type_code = Hamiltonian.convert_sequence_to_integer_sequence(central_filament_type)
        sidearm_type_code = Hamiltonian.convert_sequence_to_integer_sequence(sidearm_type)[0]

        print(sidearm_type_code)
        print(type(sidearm_type_code))

        for z in range(0,500):
            for y in range(0,500):
                for x in range(0,500):

   


                    if x == 245 or x == 254:
                        if y >= 245 and y <=  254:
                            grid[x][y][z] = 1
                            type_grid[x][y][z] = type_code[0]
                            central_filament_positions.append([x,y,z])

                    elif y == 245 or y == 254:
                        if x >= 245 and x <=  254:
                            grid[x][y][z] = 1
                            type_grid[x][y][z] = type_code[0]
                            central_filament_positions.append([x,y,z])


        sequence = central_filament_type*len(central_filament_positions)#38*500
        print(len(sequence))
        int_seq  = Hamiltonian.convert_sequence_to_integer_sequence(sequence)

        #print central_filament_positions
        CENTRAL_FILAMENT_CHAIN = chain.Chain(grid, sequence, int_seq, 1, chain_positions=central_filament_positions,fixed=True)
        chains_dict = {}
        chains_dict[1] = CENTRAL_FILAMENT_CHAIN


        # now build sidearms
        chainID=2
        sidearm_length=200
        for z in range(0, 500,2):
            print("On chain %i" % chainID)
            
            # randomly choose a side of the central filament
            # 0 - north
            # 1 - east
            # 2 - south
            # 3 - west
            side = random.randint(0,3)
            
            # depending on what side we choose define a starting position which is a 1 offset lattice
            # position on the correct side at the Z level defined by $z

            filament_location = random.randint(245,254)

            # north
            if side == 0:
                startpos = [filament_location, 255, z]

            # east
            if side == 1:
                startpos = [255, filament_location, z]

            # south
            if side == 2:
                startpos = [filament_location, 244, z]

            # west
            if side == 3:
                startpos = [244, filament_location, z]

            

            # get a list of positions for the sidearm
            sidearm_positions = self.build_chain(sidearm_length, side, startpos)
            
            # for each position update the main grid and the type_grid
            for pos in sidearm_positions:
                if not grid[pos[0]][pos[1]][pos[2]] == 0:
                    raise CustomInitializationException('Trying to assign an occupied position!')

                grid[pos[0]][pos[1]][pos[2]] = chainID
                type_grid[pos[0]][pos[1]][pos[2]] = sidearm_type_code

            # build the associated chain object
            sequence = sidearm_type*sidearm_length
            int_seq  = Hamiltonian.convert_sequence_to_integer_sequence(sequence)
            newchain = chain.Chain(grid, sequence, int_seq, chainID, chain_positions=sidearm_positions)
            
            # update the chain dictionary and the chainID
            chains_dict[chainID] = newchain
            chainID=chainID+1

                    
        
        latticeObject = lattice.Lattice(dimensions,[], Hamiltonian, chainsDict=chains_dict, lattice_grid=grid, type_grid=type_grid)

        self.LATTICE = latticeObject

        latticeObject.save_as_pdb('NEUROFILAMENT.pdb')


    def build_chain(self, length, orientation, start):
        """
        Build a chain extending outwards from the start position of length $length
        in either orientation 

        0 - nort
        1 - east
        2 - south
        3 - west

        """
        

        # north
        if orientation == 0:
            # x and z are fixed and y increases
            positions = []
            for i in range(0, length):
                positions.append([start[0],start[1]+i, start[2]])


        # east
        elif orientation == 1:
            # y and z are fixed and x increases
            positions = []
            for i in range(0, length):
                positions.append([start[0]+i,start[1], start[2]])

        # south
        elif orientation == 2:
            # x and z are fixed and y decreases
            positions = []
            for i in range(0, length):
                positions.append([start[0],start[1]-i, start[2]])

        # west
        elif orientation == 3:
            # y and z are fixed and x decreases
            positions = []
            for i in range(0, length):
                positions.append([start[0]-i,start[1], start[2]])

        return positions



                
                
            

        
                            


        
                            
                        


        
        
