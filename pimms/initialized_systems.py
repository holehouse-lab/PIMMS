## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import numpy as np
import random

from .latticeExceptions import CustomInitializationException

from . import chain
from . import lattice
from . import lattice_utils

from . CONFIG import NP_INT_TYPE

class NeurofilamentDemo:
    """
    Demonstration builder for a neurofilament-like lattice system.

    Constructs a large 3D lattice containing a central filament bundle running
    along the Z axis with randomly oriented sidearms projecting outwards, then
    builds the corresponding :class:`lattice.Lattice` object and writes it out as
    a PDB file. This is intended as a worked example / demo of programmatic
    system construction rather than as a routine simulation entry point.

    Attributes
    ----------
    LATTICE : lattice.Lattice
        The fully constructed lattice object for the neurofilament system.
    """

    def __init__(self, Hamiltonian, dimensions=None, sidearm_length=200, sidearm_z_spacing=2,
                 lattice_to_angstroms=3.65, write_pdb=True, verbose=True):
        """
        Build the neurofilament demo lattice and (optionally) write it to a PDB file.

        Generates a cubic box (500^3 by default) with a central filament tube running
        the full length of the Z axis and regularly spaced sidearms extending outwards
        in randomly chosen cardinal directions. The resulting lattice is stored on
        ``self.LATTICE`` and, by default, also saved to ``NEUROFILAMENT.pdb``.

        .. note::

           This class had rotted against the evolving :class:`~pimms.chain.Chain` and
           :class:`~pimms.lattice.Lattice` constructors (it predates the long-range /
           chain-type arguments and the ``lattice_to_angstroms`` parameter), so
           instantiating it raised ``TypeError`` before it built anything. The
           constructor calls are now current, the geometry scales with the requested
           box so it can be exercised at test size, and the site scan walks only the
           filament region rather than all ``N^3`` lattice sites.

        Parameters
        ----------
        Hamiltonian : energy.Hamiltonian
            A PIMMS Hamiltonian object, used to convert residue sequences into
            integer-coded sequences for the lattice/type grids. Must define the
            ``E`` bead type.
        dimensions : list of int, optional
            The (cubic) box, default ``[500, 500, 500]``. Must be at least ~40 per
            side so the filament tube and sidearms fit.
        sidearm_length : int, optional
            Beads per sidearm (default 200). Sidearms project straight out from the
            tube, so ``center + 5 + sidearm_length`` must fit inside the box.
        sidearm_z_spacing : int, optional
            A sidearm is placed every this-many Z layers (default 2).
        lattice_to_angstroms : float, optional
            Lattice-to-angstrom conversion used by the Lattice / PDB writer
            (default 3.65, PIMMS' standard spacing).
        write_pdb : bool, optional
            Write ``NEUROFILAMENT.pdb`` on completion (default True).
        verbose : bool, optional
            Print per-chain progress (default True).

        Returns
        -------
        None
            No return value; ``self.LATTICE`` is populated (and a PDB file written
            if requested).

        Raises
        ------
        CustomInitializationException
            If a sidearm is built into an already-occupied lattice position or
            would extend outside the box.
        """

        if dimensions is None:
            dimensions = [500, 500, 500]

        # The filament is a hollow square tube, 10 sites across, centred in X/Y and
        # running the full Z extent. For the historical 500-box this reproduces the
        # original hardcoded walls at 245/254 exactly.
        mid = dimensions[0] // 2
        wall_lo = mid - 5           # 245 for a 500 box
        wall_hi = mid + 4           # 254 for a 500 box

        if wall_lo - 1 - sidearm_length < 0 or wall_hi + 1 + sidearm_length >= min(dimensions[0], dimensions[1]):
            raise CustomInitializationException(
                f'Sidearms of length {sidearm_length} do not fit in a box of {dimensions} - '
                'shrink sidearm_length or grow the box')

        central_filament_type = 'E'
        sidearm_type = 'E'

        grid         = np.zeros(dimensions, dtype=NP_INT_TYPE)
        type_grid    = np.zeros(dimensions, dtype=NP_INT_TYPE)
        central_filament_positions = []

        type_code = Hamiltonian.convert_sequence_to_integer_sequence(central_filament_type)
        sidearm_type_code = Hamiltonian.convert_sequence_to_integer_sequence(sidearm_type)[0]

        # walk only the tube's bounding region (the original scanned every site of the
        # whole box - 125 million iterations for the default 500^3 - to place a tube
        # that occupies a 10x10 column). Iteration order (z, then y, then x) matches
        # the original so the bead ordering of the filament chain is unchanged.
        for z in range(0, dimensions[2]):
            for y in range(wall_lo, wall_hi + 1):
                for x in range(wall_lo, wall_hi + 1):
                    if x == wall_lo or x == wall_hi or y == wall_lo or y == wall_hi:
                        grid[x][y][z] = 1
                        type_grid[x][y][z] = type_code[0]
                        central_filament_positions.append([x, y, z])

        sequence = central_filament_type * len(central_filament_positions)
        int_seq    = Hamiltonian.convert_sequence_to_integer_sequence(sequence)
        LR_int_seq = Hamiltonian.convert_sequence_to_LR_integer_sequence(sequence)
        LR_IDX     = Hamiltonian.get_indices_of_long_range_residues(sequence)

        CENTRAL_FILAMENT_CHAIN = chain.Chain(grid, dimensions, sequence, int_seq, LR_int_seq, LR_IDX,
                                             chainID=1, chainType=0,
                                             chain_positions=central_filament_positions, fixed=True)
        chains_dict = {}
        chains_dict[1] = CENTRAL_FILAMENT_CHAIN

        # now build sidearms
        chainID = 2
        for z in range(0, dimensions[2], sidearm_z_spacing):
            if verbose:
                print("On chain %i" % chainID)

            # randomly choose a side of the central filament
            # 0 - north, 1 - east, 2 - south, 3 - west
            side = random.randint(0, 3)

            # depending on what side we choose define a starting position which is a
            # 1-offset lattice position on the correct side at the Z level defined by z
            filament_location = random.randint(wall_lo, wall_hi)

            if side == 0:                                     # north
                startpos = [filament_location, wall_hi + 1, z]
            if side == 1:                                     # east
                startpos = [wall_hi + 1, filament_location, z]
            if side == 2:                                     # south
                startpos = [filament_location, wall_lo - 1, z]
            if side == 3:                                     # west
                startpos = [wall_lo - 1, filament_location, z]

            # get a list of positions for the sidearm
            sidearm_positions = self.build_chain(sidearm_length, side, startpos)

            # for each position update the main grid and the type_grid
            for pos in sidearm_positions:
                if not grid[pos[0]][pos[1]][pos[2]] == 0:
                    raise CustomInitializationException('Trying to assign an occupied position!')

                grid[pos[0]][pos[1]][pos[2]] = chainID
                type_grid[pos[0]][pos[1]][pos[2]] = sidearm_type_code

            # build the associated chain object
            sequence   = sidearm_type * sidearm_length
            int_seq    = Hamiltonian.convert_sequence_to_integer_sequence(sequence)
            LR_int_seq = Hamiltonian.convert_sequence_to_LR_integer_sequence(sequence)
            LR_IDX     = Hamiltonian.get_indices_of_long_range_residues(sequence)
            newchain = chain.Chain(grid, dimensions, sequence, int_seq, LR_int_seq, LR_IDX,
                                   chainID=chainID, chainType=1,
                                   chain_positions=sidearm_positions)

            # update the chain dictionary and the chainID
            chains_dict[chainID] = newchain
            chainID = chainID + 1

        latticeObject = lattice.Lattice(dimensions, [], Hamiltonian, lattice_to_angstroms,
                                        chainsDict=chains_dict, lattice_grid=grid, type_grid=type_grid)

        self.LATTICE = latticeObject

        if write_pdb:
            latticeObject.save_as_pdb('NEUROFILAMENT.pdb')


    def build_chain(self, length, orientation, start):
        """
        Build a straight chain extending outwards from a starting position.

        Generates a list of lattice positions for a chain of the requested length
        that grows in a single cardinal direction (in the XY plane) from the
        given start position, holding the remaining coordinates fixed.

        Orientation codes:

        - 0 : north (increasing y)
        - 1 : east  (increasing x)
        - 2 : south (decreasing y)
        - 3 : west  (decreasing x)

        Parameters
        ----------
        length : int
            Number of beads (lattice positions) in the chain.
        orientation : int
            Direction in which the chain extends (0=north, 1=east, 2=south,
            3=west).
        start : list of int
            The ``[x, y, z]`` starting lattice position of the chain.

        Returns
        -------
        list of list of int
            A list of ``[x, y, z]`` lattice positions describing the chain.
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



                
                
            

        
                            


        
                            
                        


        
        
