## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2024
## ...........................................................................
# 

import numpy as np
from . import get_randmax

# Define the number of attempts that should be made for inserting
# a new chain into the molecule. Default is 20, although perhaps you
# might want to change this for some reason?
CHAIN_INIT_ATTEMPTS = 20

# run code in debug mode. Slower, but runs sanity check for functions. 
# Useful if/when testing new things and when developing code
DEBUG = False

# Inverse temperature (1/KbT) 
INVTEMP_FACTOR = 1.0

# During TSMMC number of steps spent at the top temperature - kind 
# of irrelevant but should be a specific value 
TOP_TEMP = 10 # 

# Dynamically set the maximum possible random number. This will depend
# on the system architecture. This  
C_RAND_MAX = get_randmax.get_randmax()

# Default file name for the quench file generated during quenched simulations
QUENCHFILE_NAME='QUENCH.dat'

# assumed terminal width for STDOUT
TERMINAL_WIDTH=60

# threshold number of beads needed to trigger a cluster radial density profile
# to be calculated
RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD = 27

## NB: THIS VALUE CAN BE CHANGED. To reduce PIMMS' memory footprint you
## you can change this to np.intxxx where xxx could be 16, 32 or 64. In principle
## it could be 8 but this would be quite limiting in terms of number of unique
## beads that could be used (=256, maybe fine?). NOTE that if you change
## this value you must change the corresponding CYTHON config in cython_config.pxd
NP_INT_TYPE = np.int32




## ------------------------------------------------------------------------
##                                KEYWORDS
## ------------------------------------------------------------------------

# list of ALL valid keywords. This list here  
EXPECTED_KEYWORDS = ['DIMENSIONS', 'LATTICE_TO_ANGSTROMS','CHAIN', 'TEMPERATURE', 'N_STEPS', 'PARAMETER_FILE', 'EQUILIBRATION', 
                     'RESIZED_EQUILIBRATION', 'EQUILIBRATION_OFFSET', 'HARDWALL', 'EXPERIMENTAL_FEATURES',
                     'PRINT_FREQ', 'REDUCED_PRINTING', 'XTC_FREQ', 'EN_FREQ', 'SEED', 'ENERGY_CHECK', 'ANALYSIS_FREQ', 
                     'NON_INTERACTING', 'ANGLES_OFF',
                     'CRANKSHAFT_SUBSTEPS', 'CRANKSHAFT_MODE',
                     'MOVE_CRANKSHAFT', 'MOVE_CHAIN_TRANSLATE', 'MOVE_CHAIN_ROTATE','MOVE_CHAIN_PIVOT','MOVE_HEAD_PIVOT',
                     'MOVE_SLITHER', 'MOVE_CLUSTER_TRANSLATE','MOVE_CLUSTER_ROTATE', 'MOVE_CTSMMC','MOVE_MULTICHAIN_TSMMC', 
                     'MOVE_RATCHET_PIVOT', 'MOVE_SYSTEM_TSMMC', 'MOVE_JUMP_AND_RELAX',
                     'QUENCH_RUN', 'QUENCH_FREQ', 'QUENCH_STEPSIZE', 'QUENCH_START', 'QUENCH_END', 'QUENCH_AS_EQUILIBRATION',       
                     'TSMMC_JUMP_TEMP', 'TSMMC_STEP_MULTIPLIER', 'TSMMC_INTERPOLATION_MODE', 'TSMMC_NUMBER_OF_POINTS',
                     'TSMMC_FIXED_OFFSET',
                     'ANA_POL', 'ANA_INTSCAL', 'ANA_DISTMAP', 'ANA_ACCEPTANCE', 'ANA_INTER_RESIDUE', 'ANA_CLUSTER',
                     'ANA_RESIDUE_PAIRS','WRITE_CHAIN_TO_CHAINID',
                     'ANALYSIS_MODULE','ANA_CUSTOM','ANA_CLUSTER_THRESHOLD',
                     'RESTART_FREQ','RESTART_FILE', 'RESTART_OVERRIDE_DIMENSIONS', 'RESTART_OVERRIDE_HARDWALL', 'EXTRA_CHAIN',
                     'CASE_INSENSITIVE_CHAINS', 'AUTOCENTER', 'SAVE_AT_END', 'SAVE_EQ',
                     'FREEZE_FILE']

# These keywords are the keywords that MUST be included if the simulation is going to be run, with
# the one exception of the chain keyword, which we do not make required
REQUIRED_KEYWORDS = ['DIMENSIONS', 'TEMPERATURE', 'N_STEPS', 'PARAMETER_FILE', 'EQUILIBRATION']

# list of experimental keywords (subset of EXPECTED_KEYWORDS)
# These keywords 
EXPERIMENTAL_KEYWORDS = ['TSMMC_JUMP_TEMP', 'TSMMC_STEP_MULTIPLIER', 'TSMMC_INTERPOLATION_MODE', 
                         'TSMMC_NUMBER_OF_POINTS', 'MOVE_CTSMMC','MOVE_MULTICHAIN_TSMMC', 
                         'MOVE_SLITHER', 'MOVE_MULTICHAIN_TSMMC', 'MOVE_RATCHET_PIVOT', 'MOVE_SYSTEM_TSMMC', 'MOVE_JUMP_AND_RELAX',
                         'EXTRA_CHAIN', 'FREEZE_FILE', 'EQUILIBRATION_OFFSET']


DEFAULTS = {}

DEFAULTS['SEED']        = 'random seed'  # this is overwritten in keyfile_parser..assign_defaults()
DEFAULTS['CHAIN']                       = []        # This means we can pass a RESTART_FILE
DEFAULTS['EXTRA_CHAIN']                 = []        # This means we can pass a RESTART_FILE
DEFAULTS['TEMPERATURE']                 = 'N/A'     # This means we can pass a RESTART_FILE


# major setup things
DEFAULTS['RESIZED_EQUILIBRATION']       = False
DEFAULTS['EQUILIBRATION_OFFSET']        = False     # lets 
DEFAULTS['HARDWALL']                    = False     
DEFAULTS['EXPERIMENTAL_FEATURES']       = False     # This must be set to true to use experimental features
DEFAULTS['LATTICE_TO_ANGSTROMS']        = 3.65      # note: in 0.1.34 we update this to 3.65 from 4 as used previously this is a breaking default change  
DEFAULTS['NON_INTERACTING']             = False     # use interactions 
DEFAULTS['ANGLES_OFF']                  = False     # use angles
DEFAULTS['CASE_INSENSITIVE_CHAINS']     = True      # means we cast chains to upper cahse if set to True
DEFAULTS['AUTOCENTER']                  = False     # means we do not be default center single chains in middle of box

# Output stuff
DEFAULTS['PRINT_FREQ']                  = 1000
DEFAULTS['REDUCED_PRINTING']            = False   # if set means output is printed to STDOUT at a reduced rate
DEFAULTS['XTC_FREQ']                    = 1000
DEFAULTS['EN_FREQ']                     = 1000
DEFAULTS['ENERGY_CHECK']                = 20000
DEFAULTS['ANALYSIS_FREQ']               = 1000 

# restart file stuff
DEFAULTS['RESTART_FILE']                = False     # Filename used to initialze
DEFAULTS['RESTART_OVERRIDE_DIMENSIONS'] = False     # 
DEFAULTS['RESTART_OVERRIDE_HARDWALL']   = False     # 

# quench defaults - note that other than QUENCH_RUN setting
# these to UNSET is important and is checked during initialzation
# sanity checks
DEFAULTS['QUENCH_RUN']                  =  False  # don't do a temperature change run                    
DEFAULTS['QUENCH_START']                = 'UNSET' # the name 'UNSET' gets explicitly checked in keyfile_parser() so don't change
DEFAULTS['QUENCH_END']                  = 'UNSET' # the name 'UNSET' gets explicitly checked in keyfile_parser() so don't change
DEFAULTS['QUENCH_STEPSIZE']             = 'UNSET' # the name 'UNSET' gets explicitly checked in keyfile_parser() so don't change
DEFAULTS['QUENCH_FREQ']                 = 'UNSET' # the name 'UNSET' gets explicitly checked in keyfile_parser() so don't change
DEFAULTS['QUENCH_AS_EQUILIBRATION']     = 'UNSET' # the name 'UNSET' gets explicitly checked in keyfile_parser() so don't change

# TSMMC stuff
DEFAULTS['TSMMC_JUMP_TEMP']             = 50
DEFAULTS['TSMMC_STEP_MULTIPLIER']       = 50
DEFAULTS['TSMMC_INTERPOLATION_MODE']    = 'LINEAR'
DEFAULTS['TSMMC_NUMBER_OF_POINTS']      = 20
DEFAULTS['TSMMC_FIXED_OFFSET']          = False   # don't use a fixed offset of TSMMC used

## moveset ketword stuff
DEFAULTS['CRANKSHAFT_MODE']             = 'UNIFORM'
DEFAULTS['CRANKSHAFT_SUBSTEPS']         = 500
DEFAULTS['MOVE_CRANKSHAFT']             = 0.00    # 1
DEFAULTS['MOVE_CHAIN_TRANSLATE']        = 0.00    # 2
DEFAULTS['MOVE_CHAIN_ROTATE']           = 0.00    # 3
DEFAULTS['MOVE_CHAIN_PIVOT']            = 0.00    # 4
DEFAULTS['MOVE_HEAD_PIVOT']             = 0.00    # 5
DEFAULTS['MOVE_SLITHER']                = 0.00    # 6
DEFAULTS['MOVE_CLUSTER_TRANSLATE']      = 0.00    # 7
DEFAULTS['MOVE_CLUSTER_ROTATE']         = 0.00    # 8
DEFAULTS['MOVE_CTSMMC']                 = 0.00    # 9
DEFAULTS['MOVE_MULTICHAIN_TSMMC']       = 0.00    # 10
DEFAULTS['MOVE_RATCHET_PIVOT']          = 0.00    # 11 
DEFAULTS['MOVE_SYSTEM_TSMMC']           = 0.00    # 12
DEFAULTS['MOVE_JUMP_AND_RELAX']         = 0.00    # 12

## Analysis keyword stuff
DEFAULTS['ANALYSIS_MODULE']             = False
DEFAULTS['ANA_CUSTOM']                  = 0 # by default DO NOT use custom analysis code!
DEFAULTS['ANA_RESIDUE_PAIRS']           = []
DEFAULTS['ANA_CLUSTER_THRESHOLD']       = 1 # i.e. don't do 'cluster' analysis on single chains
DEFAULTS['ANA_POL']                     = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['ANA_INTSCAL']                 = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['ANA_DISTMAP']                 = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['ANA_ACCEPTANCE']              = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['ANA_INTER_RESIDUE']           = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['ANA_CLUSTER']                 = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['ANA_CLUSTER']                 = DEFAULTS['ANALYSIS_FREQ']
DEFAULTS['WRITE_CHAIN_TO_CHAINID']      = False
DEFAULTS['FREEZE_FILE']                 = False

# will be updated to a real numerical value by the set_dynamic_defaults function unless othewise stated
DEFAULTS['RESTART_FREQ']                = "Every 10th-percentile"  # this gets explicitly checked so do not change

# saving arguments
DEFAULTS['SAVE_AT_END']         = False # By default do not hold the mdtraj object in memory for the entire simulation.
DEFAULTS['SAVE_EQ']         = True # By default, save the equilibration steps


# FINALLY we do some sanity checking here

for k in EXPECTED_KEYWORDS:
    if k not in DEFAULTS:
        if k not in REQUIRED_KEYWORDS:
            raise Exception(f'No default value set for {k} - this is a bug!')


KEYWORDS_DESCRIPTION = {
    'DIMENSIONS': ['int (2 or 3 values, e.g. A B or A B C)',
                   '[REQUIRED] - Size of the simulation box (in lattice units). 2D or 3D (defines if the simulation is a 2D or 3D simulation)'],
    'LATTICE_TO_ANGSTROMS': ['float', 'Conversion factor for converting lattice units to Angstroms. Used only to define PDB dimensions'],    
    'CHAIN': ['See description', "[REQUIRED] - One of the few multi-component keywords in PIMMS and the only keyword that can appear multiple times, the 'CHAIN' keyword defines a specific polymer chain and the number of that chain that will exist in the simulation. The format should be \n\nCHAIN : N  {CHAIN IDENTIY}\n\nWhere 'N' defines the number of the chain and '{CHAIN IDENTITY}' gives polymer sequence in one-letter alphabet code. As an example\n\nCHAIN : 20 QQQQQQQQQQ\n\nWould give 20 poly-glutamine polymers. In later versions of PIMMS we will be updating this to allow the reading of keyfiles that use three-letter codes."],
    'CASE_INSENSITIVE_CHAINS' : ["bool,", "Boolean flag which, if set to False, means that chain sequence is case sensitive. By default, this is True, which means that upon reading a keyfile, chains are converted to upper case. However, sometimes you may wish for more unique beads, in which case a lower-case chain can be useful."],    
    'TEMPERATURE': ["float (positiv)","[REQUIRED] - Simulation temperature must be a positive number greater than 0. In general a temperature between 10 and 200 is generally appropriate for the energy scales convenient for parameter files."],
    'N_STEPS':["int","[REQUIRED] - Total number of steps the simulation should be run for. Must be a positive integer value."],
    'PARAMETER_FILE': ["string", "[REQUIRED] - Filepath that points to the parameter file for the simulation. This can be a relative path or an absolute path. If the file does not exits the simulation will fail."],
    'EQUILIBRATION': ["int", "[REQUIRED] - Number of steps to be used as equilibration. During equilibration, no analysis is performed and no data is written to the trajectory file."],
    'SAVE_EQ': ["bool", "Boolean (true or false) that determines whether PIMMS saves trajectory frames for the equilibration steps of a simulation. If set to False, PIMMS begins to save your trajectory frames *after* the equilibration steps have completed."],
    'RESIZED_EQUILIBRATION': ['int (2 or 3 values, e.g. A B or A B C)', "Defines alternative simulation dimensions to be used during equilibration. MUST be smaller than the dimensions defined by the DIMENSIONS keyword"],
    'EQUILIBRATION_OFFSET': ['int (2 or 3 values, e.g. A B or A B C)', "Defines the offset of the equilibration box relative to the full simulation box. For each dimension, EQUILIBRATION_OFFSET + RESIZED_EQUILIBRATION MUST be <= DIMENSIONS"],
    'HARDWALL' :["bool", "Boolean flag set to True or False that defines whether a hardwall boundary is used or not. By default, periodic boundary conditions (PBC) are used, but if hardwall is set to true the edges of the simulation box are reflective with an infinitely repulsive potential."],
    'NON_INTERACTING' : ["bool", "Boolean flag set to True or False that defines if a non-interacting simulation should be performed or not. If set to true, all parameterfile-defined interactions are set to zero. This is convenient in that the non-interacting behavior (i.e. excluded volume limit) is a convenient reference state."],
    'ANGLES_OFF' : ["bool", "Boolean flag set to True or False that defines if angle potentials are to be used or not. If set to False (or not set), angles from the parameter file will be used. If set to True, angles are ignored and parameter files do not need to define angles."],
    'EXPERIMENTAL_FEATURES' : ["bool", "Boolean flag set to True or False that defines if experimental/non-supported keywords and features are allowed. STRONGLY recommend leaving this as False, and NONE of the features/behaviors allowed here are guaranteed to work."],
    'SEED' : ["int", 'Random seed. If not set, a random seed is generated, but it provided ensures perfect simulation reproducibility'],
    'PRINT_FREQ' : ["int", 'Frequency with which status information is printed to STDOUT'],
    'XTC_FREQ' : ["int", 'Frequency with which trajectory information is written to the traj.xtc file'],
    'EN_FREQ' : ["int", "Frequency with which the instantaneous potential energy is written to the ENERGY.dat file"],
    'ANALYSIS_FREQ' : ["int", "Master control parameter that sets default frequency for any analysis not specified by more fine-grain frequency information"],
    'ANA_POL' : ["int", "Frequency with which single-chain polymeric analysis is performed"],
    'ANA_INTSCAL' : ["int", "Frequency with which internal-scaling analysis is performed"],
    'ANA_DISTMAP' : ["int", "Frequency with which distance map analysis is performed"],
    'ANA_ACCEPTANCE' : ["int", "Frequency with acceptance ratio information is written out"],
    'ANA_INTER_RESIDUE' : ["int", "Frequency with which inter-residue distance analysis is performed (if requested). This only makes sense if ANA_RESIDUE_PAIRS has pairs of residues defined."],
    'ANA_CLUSTER' : ["int", "Frequency with which cluster analysis is performed"],
    'ANA_RESIDUE_PAIRS' : ['int (2 values)', "Two integers used to define a pair of residues, the distance between which is then calculated every ANA_INTER_RESIDUE steps. Indexing occurs from 0 (i.e., the first residue is 0. Note that at present, inter-residue distances are calculated for EVERY chain, which will trigger an error if there are chains that cannot accommodate a given pair."],
    'AUTOCENTER' : ["bool", "Boolean flag which, if set to True and you're simulating a single chain, means that the chain is automatically centered in the middle of the box. Default = False."],
    'REDUCED_PRINTING' : ["bool", "Boolean flag which, if set to True, means that the printing output is reduced"],
    'SAVE_AT_END' : ["bool", "Boolean flag which, if set to True, holds the Trajectory object in memory and only saves to .xtc at the very end. Faster but potentially more memory intensive. "],
    'WRITE_CHAIN_TO_CHAINID': ["bool", "Boolean flag which, if set to True, means we generate a file which maps each chain to its chainID. This can be useful for freeze chain diagnostics. Default = False."],
    'FREEZE_FILE': ["string", "Filepath that points to the freeze file for the simulation. This can be a relative path or an absolute path. If the file does not exits the simulation will fail. The freeze file is a file that contains a list of chain IDs that are to be frozen in place during the simulation"],
    'ENERGY_CHECK' : ["int", "Frequency with which the energy check is performed. The energy check recomputes the total energy of the system and compares it to the energy calculated by the simulation. If the energies differ an exception is raised."],
    'RESTART_FREQ' : ["int", "Frequency with which the simulation state is saved to a restart file. This allows the simulation to be restarted from the last saved state."],
    'RESTART_FILE' : ["string", "Filepath that points to the restart file for the simulation. This can be a relative path or an absolute path. If the file does not exits the simulation will fail. The restart file is a file that contains the state of the simulation at a given point in time."],
    'RESTART_OVERRIDE_DIMENSIONS' : ["bool", "Boolean flag which, if set to True, means that the dimensions of the simulation are overridden by the dimensions in the restart file. Default = False."],
    'RESTART_OVERRIDE_HARDWALL' : ["bool", "Boolean flag which, if set to True, means that the hardwall setting of the simulation is overridden by the hardwall setting in the restart file. Default = False."],
    'EXTRA_CHAIN' : ['See description', "One of the few multi-component keywords in PIMMS that should only be used if a RESTART_FILE is defined. This keyword allows you to add additional chains into the system that were not originally present in the RESTART_FILE. The format follows the same as the CHAIN keyword (so <number of chains>  <chain sequence>) and multiple EXTRA_CHAIN lines can be included for different types of chains. This means you can setup an initial set of simulations, and then run a simulation from the end-state of the original simulation with new chains added. Moreover, this can be repeated an arbitrary number of times. New chains are randomly inserted to not overlap with existing chains."],
    'QUENCH_RUN' : ["bool", "Boolean flag which, if set to True, means that the simulation is a quench run. This means that the simulation starts at one temperature and then systematically changes to a different temperature. Generally this will be higher to cooler, but could be cooler to higher. Note that the starting temperature is set by QUENCH_START and ending temperature by QUENCH_END, so the TEMPERATURE keyword is ignored. Also, all the QUENCH keywords (QUENCH_START, QUENCH_END, QUENCH_FREQ, QUENCH_STEPSIZE and QUENCH_AS_EQUILIBRATION) must all be set. Default = False."],
    'QUENCH_AS_EQUILIBRATION' : ["bool", "Boolean flag which, if set to True, means that the equilibration period is used for a quench run, and after the equilibration period the simulation temperature is fixed at the QUENCH_END temperature."],
    'QUENCH_START' : ["float", "Starting temperature for the quench run."],
    'QUENCH_END' : ["float", "Ending temperature for the quench run."],
    'QUENCH_FREQ' : ["int", "Frequency with which a change in temperature is performed."],
    'QUENCH_STEPSIZE' : ["float", "The amount by which the temperature is changed at each QUENCH_FREQ. Note this should be a positive value."],
    'MOVE_CRANKSHAFT' : ["float", "Probability of a crankshaft move being attempted. Note all provided MOVE_* keywords must add up to 1.0"],
    'CRANKSHAFT_SUBSTEPS' : ["int", "Number of subtrajectory steps to take for a crankshaft move. Generally we recommend 20-50K but this could be much larger if needed."],
    'MOVE_CHAIN_TRANSLATE' : ["float", "Probability of a molecular translation move being attempted. Note all provided MOVE_* keywords must add up to 1.0"],
    'MOVE_CHAIN_ROTATE' : ["float", "Probability of a molecular rotation move being attempted. Note all provided MOVE_* keywords must add up to 1.0"],
    'MOVE_CHAIN_PIVOT' : ["float", "Probability of a molecular pivot move being attempted. Pivot moves randomly select a point on the chain and then pivot one half of the chain. Note all provided MOVE_* keywords must add up to 1.0"],
    'MOVE_HEAD_PIVOT' : ["int", "Probability of a head pivot move being attempted. Head pivot moves randomly select one of the two ends of a chain in pivot that terminus, but this almost never worth doing so recommended setting this to 0. Note all provided MOVE_* keywords must add up to 1.0"],
    'MOVE_CLUSTER_TRANSLATE' : ["float", "Probability of a cluster translation move being attempted. Cluster translation moves are relatively expensive, so in general wise to keep this at a low number (0.01 to 0.05). Note all provided MOVE_* keywords must add up to 1.0"],
    'MOVE_CLUSTER_ROTATE' : ["float", "Probability of a cluster rotation move being attempted. Cluster rotation moves are relatively expensive, so in general wise to keep this at a low number (0.01 to 0.05). Note all provided MOVE_* keywords must add up to 1.0"]}

    
    

 
    

ONE_TO_THREE = {'A':'ALA', 
                'C':'CYS',
                'D':'ASP',
                'E':'GLU',
                'F':'PHE',
                'G':'GLY',
                'H':'HIS', 
                'I':'ILE',
                'K':'LYS',
                'L':'LEU',
                'M':'MET',
                'N':'ASN',
                'P':'PRO',
                'Q':'GLN',
                'R':'ARG',
                'S':'SER',
                'T':'THR',
                'V':'VAL',
                'W':'TRP',
                'Y':'TYR',
                'X':'XXX'}


## CARDINAL 3D ROTATION MATRICES
## 
## In the interest of speed for rotational operations in
## cardinal lattice axes (90/180/270 degrees) we define
## and set the explicit rotation matrices here. This avoids
## any need to run sin/cos functions and ensures we're exactly
## precise rather than introducing a need to round due to machine
## precision issues
##


# indicies correspond to
# 0 = rotation (90/180/270)
# 1 = axis (x/y/z)
# 2 = rotation matrix row
# 3 = rotation matrix column element 
CARDINAL_ROTATION_3D=np.zeros((3,3,3,3), dtype=int)
        
# 90 degree rotation matrix in X
CARDINAL_ROTATION_3D[0][0][0] = [1, 0,  0] # 1,      0,       0 
CARDINAL_ROTATION_3D[0][0][1] = [0, 0, -1] # 0, cos(90), -sin(90)
CARDINAL_ROTATION_3D[0][0][2] = [0, 1,  0] # 0, sin(90),  cos(90)

# 90 degree rotation matrix in Y
CARDINAL_ROTATION_3D[0][1][0] = [0,  0, 1] #  cos90,  0, sin90
CARDINAL_ROTATION_3D[0][1][1] = [0,  1, 0]  # 0     ,  1,     0 
CARDINAL_ROTATION_3D[0][1][2] = [-1, 0, 0]  # -sin90,  0, cos90 

# 90 degree rotation matrix in Z
CARDINAL_ROTATION_3D[0][2][0] = [0, -1, 0] # cos90, -sin90, 0
CARDINAL_ROTATION_3D[0][2][1] = [1, 0, 0]  # sin90, cos90,  0
CARDINAL_ROTATION_3D[0][2][2] = [0, 0, 1]  # 0    ,     0,  1


# 180 degree rotation matrix in X
CARDINAL_ROTATION_3D[1][0][0] = [1,  0,  0] # 1,      0,       0 
CARDINAL_ROTATION_3D[1][0][1] = [0, -1,  0] # 0, cos(180), -sin(180)
CARDINAL_ROTATION_3D[1][0][2] = [0,  0, -1] # 0, sin(180),  cos(180)

# 180 degree rotation matrix in Y
CARDINAL_ROTATION_3D[1][1][0] = [-1, 0,  0] #  cos180,  0, sin180
CARDINAL_ROTATION_3D[1][1][1] = [0,  1,  0]  # 0     ,  1,     0 
CARDINAL_ROTATION_3D[1][1][2] = [0,  0, -1]  #  -sin180,  0, cos180 

# 180 degree rotation matrix in Z
CARDINAL_ROTATION_3D[1][2][0] = [-1, 0, 0] # cos180, -sin180, 0
CARDINAL_ROTATION_3D[1][2][1] = [0, -1, 0]  # sin180, cos180,  0
CARDINAL_ROTATION_3D[1][2][2] = [0,  0, 1]  # 0    ,     0,  1

## 270
# 270 degree rotation matrix in X
CARDINAL_ROTATION_3D[2][0][0] = [1, 0,  0] # 1,      0,       0 
CARDINAL_ROTATION_3D[2][0][1] = [0, 0,  1] # 0, cos(270), -sin(270)
CARDINAL_ROTATION_3D[2][0][2] = [0, -1, 0] # 0, sin(270),  cos(270)

# 270 degree rotation matrix in Y
CARDINAL_ROTATION_3D[2][1][0] = [0,  0, -1] #  cos(270),  0, sin(270)
CARDINAL_ROTATION_3D[2][1][1] = [0,  1,  0] # 0     ,  1,     0 
CARDINAL_ROTATION_3D[2][1][2] = [1,  0,  0] #  -sin(270),  0, cos(270) 

# 270 degree rotation matrix in Z
CARDINAL_ROTATION_3D[2][2][0] = [0,  1, 0] # cos(270, -sin(270),  0
CARDINAL_ROTATION_3D[2][2][1] = [-1, 0, 0] # sin(270,  cos(270),  0
CARDINAL_ROTATION_3D[2][2][2] = [0,  0, 1] # 0      ,         0,  1

CARDINAL_ROTATION_2D=np.zeros((3,2,2), dtype=int)

# 90 degrees
CARDINAL_ROTATION_2D[0][0] = [ 0, -1] # cos(90), -sin(90)
CARDINAL_ROTATION_2D[0][1] = [ 1,  0] # sin(90), cos(90)

# 180 degrees
CARDINAL_ROTATION_2D[1][0] = [-1,  0] # cos(180), -sin(180)
CARDINAL_ROTATION_2D[1][1] = [ 0, -1] # sin(180), cos(180)

# 270 degrees
CARDINAL_ROTATION_2D[2][0] = [ 0,  1] # cos(270), -sin(270)
CARDINAL_ROTATION_2D[2][1] = [-1,  0] # sin(270), cos(270)

##
## Definition of filenames for default output
##

OUTNAME_NUM_CLUSTERS='NUM_CLUSTERS.dat'
OUTNAME_NUM_LR_CLUSTERS='NUM_LR_CLUSTERS.dat'
OUTNAME_CLUSTERS='CLUSTERS.dat'
OUTNAME_LR_CLUSTERS='LR_CLUSTERS.dat'
OUTNAME_CLUSTER_RG='CLUSTER_RG.dat'
OUTNAME_CLUSTER_ASPH='CLUSTER_ASPH.dat'
OUTNAME_CLUSTER_AREA='CLUSTER_AREA.dat'
OUTNAME_CLUSTER_VOL='CLUSTER_VOL.dat'
OUTNAME_CLUSTER_DENSITY='CLUSTER_DEN.dat'
OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE='CLUSTER_RADIAL_DENSITY_PROFILE.dat'

OUTNAME_LR_CLUSTER_RG='LR_CLUSTER_RG.dat'
OUTNAME_LR_CLUSTER_ASPH='LR_CLUSTER_ASPH.dat'
OUTNAME_LR_CLUSTER_AREA='LR_CLUSTER_AREA.dat'
OUTNAME_LR_CLUSTER_VOL='LR_CLUSTER_VOL.dat'
OUTNAME_LR_CLUSTER_DENSITY='LR_CLUSTER_DEN.dat'
OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE='LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat'

OUTNAME_INTERNAL_SCALING='INTSCAL.dat'
OUTNAME_INTERNAL_SCALING_SQUARED='INTSCAL_SQUARED.dat'
OUTNAME_SCALING_INFORMATION='SCALING_INFORMATION.dat'
OUTNAME_DMAP='DISTANCE_MAP.dat'
OUTNAME_ENERGY='ENERGY.dat'
OUTNAME_RG='RG.dat'
OUTNAME_ASPH='ASPH.dat'
OUTNAME_E2E='END_TO_END_DIST.dat'
OUTNAME_R2R='RES_TO_RES_DIST.dat'
OUTNAME_ACCEPTANCE='ACCEPTANCE.dat'
OUTNAME_MOVES='MOVE_FREQS.dat'
OUTNAME_TOTAL_MOVES='TOTAL_MOVES.dat'
OUTNAME_PERFORMANCE='PERFORMANCE.dat'
OUTNAME_INTER_INTRA='INTERACTIONS.dat'
OUTNAME_MIXING='MIXING.dat'
OUTNAME_LOGFILE='log.txt'

OUTPUT_USED_PARAMETER_FILE='parameters_used.prm'
OUTPUT_FULL_ANGLE_POTENTIAL='absolute_energies_of_angles.txt'
OUTPUT_CHAIN_TO_CHAINID='chain_to_chainid.txt'

RESTART_FILENAME='restart.pimms'

