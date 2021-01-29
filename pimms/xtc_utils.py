## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2021
## ...........................................................................

import mdtraj as md

from . import lattice_utils

def append_lattice_conformation_to_xtc(lattice, pdb_file, xtc_file):
    
    
    xtc_traj  = md.load(xtc_file, top=pdb_file)
    
    pdb_frame = md.load(pdb_file)

    new       = xtc_traj.join(pdb_frame)
    new.save(xtc_file)
    

def initialize_xtc_file(lattice, pdb_filename, xtc_filename):

    # first build the PDB file
    lattice_utils.open_pdb_file(lattice.dimensions, filename=pdb_filename)
    lattice_utils.write_lattice_to_pdb(lattice, filename=pdb_filename)
    lattice_utils.finish_pdb_file(pdb_filename)

    # next read the PDBFILE, and save as an xtcfile
    
    traj = md.load(pdb_filename)
    traj.save_xtc(xtc_filename)
    
