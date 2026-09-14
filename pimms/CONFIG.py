## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................
# 

import os

import numpy as np

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

# The fixed output range of the kernels' splitmix64 draws (mega_crank.PRNG_MAX), and
# the modulus used to seed numpy and the reference kernel's PRNG from the keyfile SEED.
#
# This used to be the platform's libc RAND_MAX (get_randmax.get_randmax()), a leftover
# from when the kernels drew from libc rand(). They now use splitmix64 with a FIXED
# 31-bit output range on every platform, so the value must be fixed too: with the
# platform value, the same SEED produced different streams on Windows (RAND_MAX = 32767
# there, collapsing the seed space to 15 bits) than on macOS/Linux (2**31 - 1). 2**31 - 1
# is the value macOS/Linux always had.
#
# It is NO LONGER applied to the per-megamove kernel seeds. Reducing those modulo 2**31
# left only about two billion distinct streams, so a long run replayed streams it had
# already used - 14 duplicate seeds in 300,000 consecutive megamoves, and of order 23,000
# colliding pairs expected in a run of 10^7. That biased nothing, because the seed is
# drawn independently of the configuration, but two megamoves sharing a stream are not
# independent samples. The kernels now take the full 63-bit draw.
C_RAND_MAX = 2147483647

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
                     'CRANKSHAFT_SUBSTEPS', 'CRANKSHAFT_MODE', 'SLITHER_SUBSTEPS', 'PULL_SUBSTEPS',
                     'VMMC_MAX_DISPLACEMENT', 'VMMC_MAX_CLUSTER',
                     'MOVE_CRANKSHAFT', 'MOVE_CHAIN_TRANSLATE', 'MOVE_CHAIN_ROTATE','MOVE_CHAIN_PIVOT','MOVE_HEAD_PIVOT',
                     'MOVE_SLITHER', 'MOVE_CLUSTER_TRANSLATE','MOVE_CLUSTER_ROTATE', 'MOVE_CTSMMC','MOVE_MULTICHAIN_TSMMC',
                     'MOVE_PULL', 'MOVE_SYSTEM_TSMMC', 'MOVE_JUMP_AND_RELAX', 'MOVE_VMMC',
                     'QUENCH_RUN', 'QUENCH_FREQ', 'QUENCH_STEPSIZE', 'QUENCH_START', 'QUENCH_END', 'QUENCH_AS_EQUILIBRATION',       
                     'TSMMC_JUMP_TEMP', 'TSMMC_STEP_MULTIPLIER', 'TSMMC_INTERPOLATION_MODE', 'TSMMC_NUMBER_OF_POINTS',
                     'TSMMC_FIXED_OFFSET',
                     'ANA_POL', 'ANA_INTSCAL', 'ANA_DISTMAP', 'ANA_ACCEPTANCE', 'ANA_INTER_RESIDUE', 'ANA_CLUSTER',
                     'ANA_RESIDUE_PAIRS','WRITE_CHAIN_TO_CHAINID',
                     'ANALYSIS_MODULE','ANA_CUSTOM','ANA_CLUSTER_THRESHOLD',
                     'RESTART_FREQ','RESTART_FILE', 'RESTART_OVERRIDE_DIMENSIONS', 'RESTART_OVERRIDE_HARDWALL', 'RESTART_CONTINUE', 'EXTRA_CHAIN',
                     'CASE_INSENSITIVE_CHAINS', 'AUTOCENTER', 'SAVE_AT_END', 'SAVE_EQ',
                     'TRAJECTORY_PBC_UNWRAP',
                     'FREEZE_FILE', 'PARALLELIZE', 'PARALLEL_THREADS']

# These keywords are the keywords that MUST be included if the simulation is going to be run, with
# the one exception of the chain keyword, which we do not make required
REQUIRED_KEYWORDS = ['DIMENSIONS', 'TEMPERATURE', 'N_STEPS', 'PARAMETER_FILE', 'EQUILIBRATION']

# list of experimental keywords (subset of EXPECTED_KEYWORDS)
# These keywords require EXPERIMENTAL_FEATURES : True to be set when used away from
# their default value. Only the VMMC move (and its tuning keywords) remains gated;
# the megamoves (pull, slither, the three TSMMC variants, jump-and-relax), non-cubic
# boxes, and the EXTRA_CHAIN / FREEZE_FILE / EQUILIBRATION_OFFSET keywords have all
# graduated out of this list and can now be used without the gate.
EXPERIMENTAL_KEYWORDS = ['MOVE_VMMC', 'VMMC_MAX_DISPLACEMENT', 'VMMC_MAX_CLUSTER']


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
DEFAULTS['RESTART_CONTINUE']            = False     # resume the run the restart file was written by

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
DEFAULTS['TSMMC_JUMP_TEMP']             = 50.0
DEFAULTS['TSMMC_STEP_MULTIPLIER']       = 50
DEFAULTS['TSMMC_INTERPOLATION_MODE']    = 'LINEAR'
DEFAULTS['TSMMC_NUMBER_OF_POINTS']      = 20
DEFAULTS['TSMMC_FIXED_OFFSET']          = False   # don't use a fixed offset of TSMMC used

## moveset ketword stuff
DEFAULTS['CRANKSHAFT_MODE']             = 'UNIFORM'
DEFAULTS['CRANKSHAFT_SUBSTEPS']         = 500
DEFAULTS['SLITHER_SUBSTEPS']            = 10   # number of slither (reptation) moves applied to EACH chain per slither megamove
DEFAULTS['PULL_SUBSTEPS']               = 10   # number of pull moves applied to EACH chain per pull megamove
DEFAULTS['VMMC_MAX_DISPLACEMENT']       = 3    # max |translation| per dimension for a VMMC collective move
DEFAULTS['VMMC_MAX_CLUSTER']            = 1000 # cap on VMMC cluster size (clamped to the number of chains at runtime)
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
DEFAULTS['MOVE_PULL']                   = 0.00    # 11
DEFAULTS['MOVE_SYSTEM_TSMMC']           = 0.00    # 12
DEFAULTS['MOVE_JUMP_AND_RELAX']         = 0.00    # 13
DEFAULTS['MOVE_VMMC']                   = 0.00    # 14

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
DEFAULTS['TRAJECTORY_PBC_UNWRAP'] = False # By default write raw lattice positions (chains may be split across PBC)

# parallelization of the crankshaft (system_shake) move
DEFAULTS['PARALLELIZE']        = False  # By default use the (serial) optimized kernel
DEFAULTS['PARALLEL_THREADS']   = 0      # 0 => auto (use all available CPU cores)


# FINALLY we do some sanity checking here

for k in EXPECTED_KEYWORDS:
    if k not in DEFAULTS:
        if k not in REQUIRED_KEYWORDS:
            raise Exception(f'No default value set for {k} - this is a bug!')


KEYWORDS_DESCRIPTION = {
    'DIMENSIONS': ['int (2 or 3 values, e.g. A B or A B C)',
                   '[REQUIRED] - Size of the simulation box (in lattice units). Providing 2 values runs a 2D simulation, 3 values a 3D simulation. The axes need NOT be equal: non-cubic/non-square boxes (e.g. 10 20 40) are fully supported with either HARDWALL or periodic boundaries. Every axis must be at least 7 lattice sites (the smallest box that supports the super-long-range interaction shell). The only other restriction is that cluster-rotation moves (MOVE_CLUSTER_ROTATE) cannot be combined with a non-cubic box under periodic boundaries, because a 90-degree rigid rotation is only an energy-preserving symmetry of a cube/square (or of any box under HARDWALL, where there is no periodic wrapping).'],
    'LATTICE_TO_ANGSTROMS': ['float (positive)', 'Conversion factor (default 3.65) for converting lattice units to Angstroms when writing the START.pdb topology and traj.xtc trajectory. Must be greater than 0. This is purely cosmetic - it sets the bead spacing seen in a viewer/analysis (mdtraj reports nm, i.e. lattice_units x LATTICE_TO_ANGSTROMS x 0.1) and has NO effect on the simulation itself or its energetics. It is also the value lemonade needs to reconstruct the lattice from a trajectory, so keep a record of it.'],
    'CHAIN': ['See description', "[REQUIRED] - One of the few keywords that can appear multiple times (the others are EXTRA_CHAIN and ANA_RESIDUE_PAIRS), the 'CHAIN' keyword defines a specific polymer chain and the number of that chain that will exist in the simulation. The format should be \n\nCHAIN : N  {CHAIN IDENTIY}\n\nWhere 'N' defines the number of the chain and '{CHAIN IDENTITY}' gives polymer sequence in one-letter alphabet code. As an example\n\nCHAIN : 20 QQQQQQQQQQ\n\nWould give 20 poly-glutamine polymers. N must be an integer of 1 or more, and the sequence may not contain the character 0, which is reserved for solvent. Each CHAIN line is a distinct chain type, numbered from 0 in the order the lines appear; chainIDs are numbered from 1 across the whole system. This keyword is required UNLESS a RESTART_FILE is provided, in which case the chains are taken from the restart file and CHAIN may be omitted (any CHAIN lines are then ignored). Chains are placed at random (a single chain is placed in the middle of the box); if PIMMS cannot find room for a chain it aborts with an 'overcrowded lattice' error, so keep the occupied volume fraction reported at start-up sensible. In later versions of PIMMS we will be updating this to allow the reading of keyfiles that use three-letter codes."],
    'CASE_INSENSITIVE_CHAINS' : ["bool", "Boolean flag which, if set to False, means that chain sequence is case sensitive. By default, this is True, which means that upon reading a keyfile, CHAIN and EXTRA_CHAIN sequences are converted to upper case, so a chain written with 'a' becomes 'A' and a lower-case bead type can never appear in a chain. Sometimes you may want those extra unique beads, in which case setting this to False is useful. Note the parameter file itself is never case-folded. Every bead type used in a chain must also be defined in the parameter file."],
    'TEMPERATURE': ["float (positive)","[REQUIRED] - Simulation temperature; must be a positive number greater than 0. In general a temperature between 10 and 200 is appropriate for the energy scales typical of PIMMS parameter files. Higher temperatures sample more expanded/disordered states; lower temperatures favour collapse/assembly. Note that PIMMS works in reduced units with k = 1, so the Boltzmann factor is exp(-dE/TEMPERATURE) and the temperature is on the same scale as the parameter-file energies. In a QUENCH_RUN this value is overwritten by QUENCH_START (a warning is printed if the two differ), so use QUENCH_START / QUENCH_END to set the ramp."],
    'N_STEPS':["int","[REQUIRED] - Total number of outer-loop steps to run (including the EQUILIBRATION steps). Must be a positive integer, and must be larger than EQUILIBRATION. Note that one step is typically a great deal of Monte Carlo work: each crankshaft step performs CRANKSHAFT_SUBSTEPS single-bead sub-moves in total (spread at random over all beads), and the slither/pull/TSMMC moves are likewise 'megamoves', so the true number of accept/reject operations is far larger than N_STEPS (see TOTAL_MOVES.dat)."],
    'PARAMETER_FILE': ["string", "[REQUIRED] - Filepath (relative or absolute, with a leading ~ expanded) to the parameter file defining the interaction energies and angle penalties. The simulation fails if it does not exist. The exact parameters used are echoed to parameters_used.prm at startup."],
    'EQUILIBRATION': ["int", "[REQUIRED] - Number of initial steps treated as equilibration. Must be 0 or larger and smaller than N_STEPS. Step number EQUILIBRATION is the LAST equilibration step; production begins at EQUILIBRATION + 1. During equilibration no ANA_* analysis output and no restart snapshot is written (ENERGY.dat and PERFORMANCE.dat are still written throughout, and trajectory frames are saved if SAVE_EQ is True). Choose this large enough that ENERGY.dat has plateaued before production begins."],
    'SAVE_EQ': ["bool", "Boolean (true or false) that determines whether PIMMS saves trajectory frames for the equilibration steps of a simulation. If set to False, PIMMS begins to save your trajectory frames *after* the equilibration steps have completed. Note frame 0 of traj.xtc is always the starting configuration, regardless of this setting (for a RESIZED_EQUILIBRATION run it is the post-resize configuration at step EQUILIBRATION, and eq_START.pdb / eq_traj.xtc exist only when SAVE_EQ is True)."],
    'RESIZED_EQUILIBRATION': ['int (2 or 3 values, e.g. A B or A B C)', "Defines a smaller box to use during equilibration; at the end of equilibration the box is grown to the full DIMENSIONS (with chains re-centred, or placed using EQUILIBRATION_OFFSET). Useful for condensing/assembling a system at high effective concentration before expanding to the production box. Must have the same number of values as DIMENSIONS, must be <= DIMENSIONS in every dimension, and every axis must satisfy the same >= 7 floor as DIMENSIONS (this box IS simulated). The equilibration phase is always run under hardwall boundaries (forced internally, so a system is never resized while chains straddle a periodic face); your production HARDWALL setting takes over once the box has grown. While the resized box is in use the trajectory is written to eq_START.pdb / eq_traj.xtc rather than START.pdb / traj.xtc, and those files are only written if SAVE_EQ is True. Setting EQUILIBRATION to 0 deactivates this keyword with a warning. Incompatible with RESTART_OVERRIDE_DIMENSIONS and with periodic (non-hardwall) restart files; with a hardwall restart file the restart box must be <= the RESIZED_EQUILIBRATION box. See also EQUILIBRATION_OFFSET."],
    'EQUILIBRATION_OFFSET': ['int (2 or 3 values, e.g. A B or A B C)', "Defines the offset of the equilibration box relative to the full simulation box, i.e. where the small box sits inside the production box when it is grown (without it the configuration is re-centred). Requires RESIZED_EQUILIBRATION to be set, and must have the same number of values. Every value must be >= 0, and for each dimension EQUILIBRATION_OFFSET + RESIZED_EQUILIBRATION MUST be <= DIMENSIONS."],
    'HARDWALL' :["bool", "Boolean flag set to True or False that defines whether a hardwall boundary is used or not. By default (False) periodic boundary conditions are used, so the box wraps in every dimension. If set to True the box has hard walls: no bead ever sits outside the box, no bond ever crosses a wall, and no chain is ever left straddling a face. Note that this constrains the states, not the paths between them: a whole-chain rigid translation (MOVE_CHAIN_TRANSLATE) draws its offset uniformly over the box and wraps, so a chain sitting against one wall can be relocated in one move to the opposite wall. Every state it visits is a legal confined state and the proposal is symmetric, so the sampled ensemble is exactly the confined Boltzmann distribution; but the move is a relocation, not a physical passage through the wall, and a trajectory will show the chain jumping across the box. If you need every step to be spatially local, leave MOVE_CHAIN_TRANSLATE at 0. The walls are energetically solvent - a bead beside a wall picks up its bead-solvent energy for each out-of-box neighbour site - so a wall excludes volume but adds no special surface energy of its own. Nothing is read across a wall either: contacts and clusters never connect two beads through a face, and coordinates are plain Cartesian because there is no periodic image to reconstruct."],
    'NON_INTERACTING' : ["bool", "Boolean flag set to True or False that defines if a non-interacting simulation should be performed or not. If set to true, all bead-bead (SR/LR/SLR) interaction and solvation energies are set to zero. Angle penalties are NOT affected - combine with ANGLES_OFF : True for a fully ideal excluded-volume reference state."],
    'ANGLES_OFF' : ["bool", "Boolean flag set to True or False that defines if angle potentials are to be used or not. If set to False (or not set), angles from the parameter file will be used, and EVERY bead type with interaction energies must then carry an ANGLE_PENALTY (or ANGLE_PENALTY_T_NORM) line or the parameter file is rejected. If set to True, angles are ignored (every penalty is forced to zero, which is announced at start-up unless REDUCED_PRINTING is on) and parameter files do not need to define angles."],
    'EXPERIMENTAL_FEATURES' : ["bool", "Boolean flag set to True or False that defines if experimental/non-supported keywords and features are allowed. As of this release the only gated keywords are MOVE_VMMC, VMMC_MAX_DISPLACEMENT and VMMC_MAX_CLUSTER; setting any of them away from its default without EXPERIMENTAL_FEATURES : True is an error. STRONGLY recommend leaving this as False, and NONE of the features/behaviors allowed here are guaranteed to work."],
    'SEED' : ["int (positive)", 'Random seed; must be a positive integer. If not set, a random seed is generated (and announced at start-up), but if provided it ensures perfect simulation reproducibility (an identical keyfile + parameter file + seed reproduces the trajectory bit-for-bit on the same platform and PIMMS version).'],
    'PRINT_FREQ' : ["int (positive)", 'Frequency (in steps) with which a status line (step number, percentage complete and current energy) is printed to STDOUT. Must be greater than 0. Throughput and the estimated time remaining are reported separately at every 5 % of the run (plus once at step 20), and written to PERFORMANCE.dat.'],
    'XTC_FREQ' : ["int (positive)", 'Frequency (in steps, default 1000) with which a trajectory frame is written to traj.xtc (with START.pdb as the topology). Must be greater than 0. Equilibration frames are only written if SAVE_EQ is True; if SAVE_AT_END is True the trajectory is buffered in memory and written once at the end.'],
    'EN_FREQ' : ["int (positive)", "Frequency (in steps) with which the instantaneous potential energy is appended to ENERGY.dat (tab-separated: step, energy). Must be greater than 0. ENERGY.dat is written throughout the run, including during equilibration, though as with every output file it appears only once a row has actually been written to it (so a run shorter than EN_FREQ produces no ENERGY.dat)."],
    'ANALYSIS_FREQ' : ["int", "Master control parameter that sets the default frequency for any analysis whose own ANA_* frequency keyword is not explicitly provided. Set the per-analysis ANA_* keywords to override it. Setting an ANA_* frequency (or ANALYSIS_FREQ itself, or ENERGY_CHECK) to 0 or less disables that analysis for the run. Note this 0-means-disabled convention applies ONLY to the ANA_* keywords and ENERGY_CHECK: PRINT_FREQ, XTC_FREQ, EN_FREQ and RESTART_FREQ must all be positive."],
    'ANA_POL' : ["int", "Frequency with which single-chain polymeric analysis is performed. Writes per-chain radius of gyration to RG.dat, asphericity to ASPH.dat and end-to-end distance to END_TO_END_DIST.dat; each row is the step number followed by one column per chain, in ascending chainID order."],
    'ANA_INTSCAL' : ["int", "Frequency with which internal-scaling analysis is performed. Accumulates the mean internal scaling R(s) as a function of sequence separation s over the run and, at the end, writes INTSCAL.dat and INTSCAL_SQUARED.dat and the fitted SCALING_INFORMATION.dat (scaling exponent nu and prefactor R0); the mean DISTANCE_MAP.dat is controlled separately by ANA_DISTMAP. For multi-component systems these are written per chain type as CHAIN_<TYPE>_* (where <TYPE> is the 0-based integer chain-type index, in keyfile CHAIN-line order) INSTEAD OF the unprefixed files - the unprefixed names are used only when the system has a single chain type."],
    'ANA_DISTMAP' : ["int", "Frequency with which the mean inter-residue distance map is accumulated. The seqlen x seqlen mean distance matrix is written to DISTANCE_MAP.dat at the end of the run (per chain type as CHAIN_<TYPE>_DISTANCE_MAP.dat)."],
    'ANA_ACCEPTANCE' : ["int", "Frequency with which move statistics are written: MOVE_FREQS.dat (attempted moves per move code), ACCEPTANCE.dat (accepted moves per move code) and TOTAL_MOVES.dat (total attempted MC moves, including every megamove and TSMMC sub-move). All three are RUNNING TOTALS since the start of the run, not per-interval counts. Each MOVE_FREQS/ACCEPTANCE row is the step number followed by 14 tab-separated columns, one per move code 1-14 (each MOVE_* keyword description states its code); divide ACCEPTANCE by MOVE_FREQS to get the per-move acceptance ratio."],
    'ANA_INTER_RESIDUE' : ["int", "Frequency with which inter-residue distance analysis is performed and appended to RES_TO_RES_DIST.dat. This only makes sense if ANA_RESIDUE_PAIRS has a pair of residues defined - with no pair defined there is nothing to measure and no RES_TO_RES_DIST.dat is written; the distance is computed for EVERY chain, so all chains must be long enough to contain the pair (this is checked at start-up)."],
    'ANA_CLUSTER' : ["int", "Frequency with which cluster analysis is performed. Identifies short-range clusters (chains in Chebyshev-1 contact) and long-range clusters (chains connected by any pair with a NONZERO interaction energy: a short-range contact, or a Chebyshev-2 / Chebyshev-3 pair whose LR / SLR table entry is nonzero) and writes their size distributions (CLUSTERS.dat and NUM_CLUSTERS.dat, plus LR_CLUSTERS.dat and NUM_LR_CLUSTERS.dat) plus per-cluster radius of gyration, asphericity, surface area, volume and density (CLUSTER_RG/ASPH/AREA/VOL/DEN.dat, and LR_CLUSTER_RG/ASPH/AREA/VOL/DEN.dat for the long-range clusters). Its cost is driven by the number of chains and by the size of the largest cluster (it walks the connected components of the whole system and then computes convex hulls and radial profiles for each), so keep its frequency low for large, many-chain or strongly condensed systems. See ANA_CLUSTER_THRESHOLD."],
    'ANA_RESIDUE_PAIRS' : ['int (2 values)', "Two integers used to define a pair of residues, the distance between which is then calculated every ANA_INTER_RESIDUE steps. May be repeated to monitor several pairs (one RES_TO_RES_DIST.dat line per pair per recorded step). Indexing occurs from 0 (i.e., the first residue is 0) and negative indices are rejected; the two values are stored in ascending order, so the order you write them in does not matter. Note that at present, inter-residue distances are calculated for EVERY chain, so the larger index must be inside every chain in the system (including chains from a RESTART_FILE or EXTRA_CHAIN) - if it is not, the run is rejected at start-up."],
    'AUTOCENTER' : ["bool", "Boolean flag which, if set to True and you are simulating a SINGLE chain, re-centres that chain in the middle of the box in every WRITTEN frame and in the START.pdb topology. This is a visualisation convenience only - the simulation, its energetics and every analysis file are unaffected - and it is silently ignored when the system contains more than one chain. Default = False."],
    'REDUCED_PRINTING' : ["bool", "Boolean flag which, if set to True, silences the chatty per-event STDOUT messages - the 'saving coordinates' remark, the per-checkpoint 'writing restart file' line, the per-move clash/cluster-resize/multichain-rearrangement notes, the quench temperature updates and the TSMMC accept/reject lines. The PRINT_FREQ status line, the energy checks, the start-up summary and everything written to file are unaffected. Useful for long runs whose STDOUT is being captured to a file. Default = False."],
    'SAVE_AT_END' : ["bool", "Boolean flag which, if set to True, holds the Trajectory object in memory and only saves to .xtc at the very end. Faster but potentially more memory intensive - the whole trajectory is buffered, so memory grows with the number of frames (N_STEPS / XTC_FREQ) times the number of beads. If the run aborts on a failed ENERGY_CHECK the buffered frames are still flushed to disk. Default = False."],
    'TRAJECTORY_PBC_UNWRAP' : ["bool", "Boolean flag (default False). Under periodic boundaries a chain that crosses a box face is stored split across the two faces, which looks broken in a viewer. If set to True, PIMMS makes every chain WHOLE before writing each trajectory frame (and the START.pdb topology): each chain is shifted into a single periodic image, effectively extending the lattice beyond the box in x/y/z as needed, so no chain is torn across a boundary. This is purely a visualisation convenience - it does not affect the simulation or its energetics, and coordinates may fall outside the box (the unit cell is unchanged). Has no effect with HARDWALL (chains never cross a boundary). Default = False."],
    'WRITE_CHAIN_TO_CHAINID': ["bool", "Boolean flag which, if set to True, writes chain_to_chainid.txt: one tab-separated line per chain giving its chainID, its length and its sequence. Useful for working out which chainIDs to list in a FREEZE_FILE, and for mapping the per-chain columns of RG.dat / ASPH.dat / END_TO_END_DIST.dat back to sequences. Default = False."],
    'FREEZE_FILE': ["string", "Filepath (relative or absolute, with a leading ~ expanded) to a freeze file; the simulation fails if it does not exist, and an empty value is rejected. The freeze file is a plain-text file listing chainIDs to hold fixed for the whole run, one or more lines of the form 'C <id> <id> ...'. Frozen chains never move but still contribute to the energy (other chains feel them), and the collective moves (cluster translate/rotate, VMMC) reject any move whose cluster would contain a frozen chain. Use WRITE_CHAIN_TO_CHAINID to discover chainIDs (they are numbered from 1). Freezing every chain in the system is refused at start-up, since nothing could then change. Frozen chains are honoured by both the serial and the parallel (PARALLELIZE) move kernels, so freezing and parallelization can be used together."],
    'PARALLELIZE': ["bool", "Boolean flag (True/False) which, if set to True, runs the crankshaft (MOVE_CRANKSHAFT), slither (MOVE_SLITHER) and pull (MOVE_PULL) moves on multi-threaded checkerboard kernels instead of the serial kernels (all other moves stay serial). Works in both 2D and 3D. Beneficial for large, spatially dispersed systems; gives little benefit for small boxes (which decompose into a single block) or collapsed/dense single-droplet systems. For the whole-chain moves (slither and pull) a chain only parallelizes if all its beads fit inside one block's interior (chains spanning a block boundary are frozen that sweep). Since only a short enough chain can ever fit, the chains are split ONCE by length at the start of the run: chains no longer than the smallest block interior - and within the kernels' 512-bead per-chain limit (heteropolymer chains for slither, all chains for pull) - are moved by the parallel kernel, and every longer chain is moved by the serial kernel, with both passes run in every megamove. The split depends only on chain lengths, the box and the frozen set, never on the current configuration, which is what keeps PARALLELIZE from changing the equilibrium being sampled. Deciding per megamove from the chains' current extents (the behaviour up to 1.0.8) does change it: the parallel kernel can never extend a chain past a block interior while the serial kernel can, so a state-dependent choice between them biases the run towards compact chains. The block decomposition is independent of the thread count, so results are identical for any number of threads. NOTE the parallel sampler targets the same equilibrium distribution but follows a DIFFERENT (and per-step slower-relaxing) Markov chain than the serial run: each sweep only the beads inside the block interiors can move (the frozen halos are re-drawn every sweep), so a run that has not reached equilibrium - e.g. a collapsing system - will show a different (less relaxed) energy at a given step count than the serial run. Compare equilibrium averages, not energies at a fixed step. Frozen chains (via FREEZE_FILE) are fully supported: their beads are excluded from moves but kept in place as fixed, energy-contributing obstacles, so PARALLELIZE applies even with a freeze file. If the compiled kernels were built without OpenMP the blocks are executed one after another instead (no speed-up, identical sampling); the start-up parallelization report says which kernel each move will really use and why. Default = False."],
    'PARALLEL_THREADS': ["int", "Number of OpenMP threads used when PARALLELIZE is True. Must be 0 or greater; 0 (the default) means use all available CPU cores. Ignored when PARALLELIZE is False. Results do not depend on the thread count."],
    'ENERGY_CHECK' : ["int", "Frequency (in steps) with which a full from-scratch energy recompute is compared against the incrementally tracked energy, and the occupancy/type grids are cross-checked against the chain objects; a mismatch raises an exception (and dumps the state to CONFIG_AT_ENERGY_FAIL.pdb/.xtc). Each check also prints the full energy decomposition (short-range, long-range, super-long-range, angles) to STDOUT. This is an O(N) safety/debugging check - cheap to run occasionally, expensive every step for large systems. Set to 0 to switch it off entirely (not recommended)."],
    'RESTART_FREQ' : ["int (positive)", "Frequency with which the simulation state is saved to restart.pimms. May be a positive integer step frequency, or the default sentinel 'Every 10th-percentile' which writes every N_STEPS/10 steps once production begins (restart snapshots, like all analysis, are suppressed during equilibration). A restart of the final state is always written when the run completes. Each write is atomic, so an interrupted run never leaves a torn restart.pimms behind. See the Restart files documentation."],
    'RESTART_FILE' : ["string", "Filepath (relative or absolute, with a leading ~ expanded) to a restart.pimms file to start the simulation from; the simulation fails if it does not exist, and an empty value is rejected. When set, it supplies the initial configuration and the CHAIN keyword is not required (the chains come from the restart file, and any keyfile CHAIN lines are ignored). Restart files are Python pickles, so only load files you generated yourself or otherwise trust. See the Restart files documentation for the dimension/hardwall compatibility rules."],
    'RESTART_OVERRIDE_DIMENSIONS' : ["bool", "If True, IGNORE the keyfile DIMENSIONS and adopt the restart file's box exactly as it was saved (a convenience for continuing in the original box without repeating its size in the keyfile). DIMENSIONS is still a required keyword and must still be present (any valid box will do), and the box taken from the restart file is re-validated against the >= 7 per-axis floor. It does NOT grow the box, and is incompatible with RESIZED_EQUILIBRATION. If False (default), the keyfile DIMENSIONS is used and reconciled with the restart: for a HARDWALL restart the keyfile box must be >= the restart box in every axis, and a larger box is grown with the configuration re-centred inside it (growing into a bigger box therefore needs NO override); for a periodic (PBC) restart the keyfile DIMENSIONS must match the restart box exactly. The box can never be made smaller than the restart box. Default = False."],
    'RESTART_CONTINUE' : ["bool", "Boolean flag which, if set to True, resumes the run that wrote RESTART_FILE as though it had never stopped: the step counter continues from the step the restart file was written at, the random-number generators are restored to their state at that step, and in a QUENCH_RUN the temperature is restored to its value at that step, so the resumed run reproduces the uninterrupted run bit for bit (the same ENERGY.dat rows, analysis rows and moves, labelled with the same step numbers). N_STEPS is the TOTAL length of the run, so it must be larger than the restart file's step, and the resumed segment writes rows for the steps after that one; frame 0 of its trajectory is the restart configuration. Requires RESTART_FILE, and a restart file written by PIMMS 1.0.8 or later (earlier files carry no generator state and are refused). SEED must not be given, since the generator state comes from the file; the box, boundary condition and chain list must be exactly those of the original run, so RESIZED_EQUILIBRATION and EXTRA_CHAIN are refused and DIMENSIONS / HARDWALL must equal the restart file's (RESTART_OVERRIDE_DIMENSIONS and RESTART_OVERRIDE_HARDWALL are the easy way to guarantee that). Run each segment in its own directory. Default = False."],
    'RESTART_OVERRIDE_HARDWALL' : ["bool", "Boolean flag which, if set to True, means that the hardwall setting of the simulation is overridden by the hardwall setting in the restart file, rather than using the keyfile HARDWALL value. It does NOT relax the other restart rules: a periodic restart file still cannot be combined with RESIZED_EQUILIBRATION, and the cluster-rotation box rule is applied to whichever HARDWALL/DIMENSIONS combination the run ends up with. Default = False."],
    'EXTRA_CHAIN' : ['See description', "One of the few multi-component keywords in PIMMS, and it can ONLY be used when a RESTART_FILE is defined (using it without one is an error). This keyword allows you to add additional chains into the system that were not originally present in the RESTART_FILE. The format follows the same as the CHAIN keyword (so <number of chains>  <chain sequence>) and multiple EXTRA_CHAIN lines can be included for different types of chains. This means you can setup an initial set of simulations, and then run a simulation from the end-state of the original simulation with new chains added. Moreover, this can be repeated an arbitrary number of times. New chains are randomly inserted so they do not overlap with existing chains, and are given chainIDs that follow on from the restart file's. An extra chain whose sequence already exists in the restart file joins that existing chain type rather than defining a new one."],
    'QUENCH_RUN' : ["bool", "Boolean flag which, if set to True, means that the simulation is a quench run. This means that the simulation starts at one temperature and then systematically changes to a different temperature. Generally this will be higher to cooler, but could be cooler to higher. Note that the starting temperature is set by QUENCH_START and ending temperature by QUENCH_END, so the TEMPERATURE keyword is overwritten with QUENCH_START (with a warning if the two disagree). Also, all the QUENCH keywords (QUENCH_START, QUENCH_END, QUENCH_FREQ, QUENCH_STEPSIZE and QUENCH_AS_EQUILIBRATION) must all be set. The temperature trajectory is recorded in QUENCH.dat, and any ANGLE_PENALTY_T_NORM values in the parameter file are scaled once by QUENCH_END (the production temperature), so T-normalised angle penalties are exact in kT only at the end of the ramp. Default = False."],
    'QUENCH_AS_EQUILIBRATION' : ["bool", "Boolean flag which, if set to True, means that the equilibration period is used for a quench run, and after the equilibration period the simulation temperature is fixed at the QUENCH_END temperature. Setting this OVERWRITES whatever EQUILIBRATION you gave with the length of the quench, (1 + number of temperature changes) x QUENCH_FREQ steps, which is announced at start-up. The extra QUENCH_FREQ steps let the system equilibrate at the final temperature before production begins."],
    'QUENCH_START' : ["float (positive)", "Starting temperature for the quench run; must be greater than 0. Also becomes the simulation TEMPERATURE."],
    'QUENCH_END' : ["float (positive)", "Ending (production) temperature for the quench run; must be greater than 0. May be above or below QUENCH_START (a cooling or a heating ramp), and is the temperature that any ANGLE_PENALTY_T_NORM values are scaled by."],
    'QUENCH_FREQ' : ["int (positive)", "Number of steps spent at each temperature before the next change; must be greater than 0. The whole quench occupies (1 + number of temperature changes) x QUENCH_FREQ steps, which must be smaller than N_STEPS or the run is rejected at start-up."],
    'QUENCH_STEPSIZE' : ["float (positive)", "The amount by which the temperature is changed at each QUENCH_FREQ steps. Give it as a magnitude: the sign is ignored (the absolute value is taken and PIMMS works out the direction from QUENCH_START and QUENCH_END). It must be greater than 0 and no larger than the QUENCH_START to QUENCH_END range. The final step is clamped so the ramp lands exactly on QUENCH_END."],
    'MOVE_CRANKSHAFT' : ["float", "Move code 1. Probability of a crankshaft megamove being attempted. When selected, CRANKSHAFT_SUBSTEPS single-bead sub-moves are performed across the whole system, each with its own Metropolis accept/reject, so a single crankshaft step is a large amount of MC work. This is the workhorse move for relaxing local chain conformation and should usually carry most of the probability. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'CRANKSHAFT_SUBSTEPS' : ["int (positive)", "Number of single-bead sub-moves performed per crankshaft megamove, in TOTAL across the whole system (each sub-move picks a bead at random from every bead on the lattice, so the work is not per chain). Must be greater than 0; the default of 500 is deliberately small, and for production runs 20-50K (or more) is a more usual choice. The same value sets the length of the two relaxations inside a MOVE_JUMP_AND_RELAX move."],
    'SLITHER_SUBSTEPS' : ["int (positive)", "Number of slither (reptation) moves applied to EACH chain, in random order, per slither megamove. Must be greater than 0. A slither advances a chain forwards or backwards like a snake. For homopolymers the energy is evaluated in O(1) (only the moved end matters); for heteropolymers every residue is re-evaluated; single-bead chains become a local translation."],
    'PULL_SUBSTEPS' : ["int (positive)", "Number of pull moves applied to EACH eligible chain, in random order, per pull megamove. Must be greater than 0. A pull move displaces an interior bead and cooperatively 'pulls' the rest of the segment along to restore connectivity, letting chains rearrange in dense systems where rigid moves would clash. Chains shorter than 3 beads have no interior bead and are skipped."],

    'MOVE_VMMC' : ["float", "Move code 14. Probability of a Virtual-Move Monte Carlo (VMMC) collective move being attempted (Whitelam and Geissler, J. Chem. Phys. 127, 154101, 2007). A seed chain is given a trial rigid translation; neighbouring chains are recruited into a moving cluster according to interaction-energy gradients (a neighbour is recruited when moving the seed alone would break their mutual attraction), and the whole cluster translates together. This avoids the kinetic traps that single-chain moves hit in strongly-attractive / condensed phases, while maintaining detailed balance. Seeding on, or recruiting, a frozen chain rejects the move, and each attempt costs a full from-scratch energy evaluation (O(N) in the number of beads), so keep the fraction modest for large systems. EXPERIMENTAL - requires EXPERIMENTAL_FEATURES : True. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'VMMC_MAX_DISPLACEMENT' : ["int (positive)", "Maximum magnitude (per dimension, in lattice units) of the rigid translation proposed by a VMMC move. Must be greater than 0. Small values give local collective moves (recommended); large values rarely succeed in dense phases. Default 3. Changing this from the default requires EXPERIMENTAL_FEATURES : True (VMMC is experimental)."],

    'VMMC_MAX_CLUSTER' : ["int (positive)", "Upper bound on the VMMC cluster size used for the 1/n_c move-frequency correction; the recruited cluster is aborted (move rejected) if it would exceed the drawn cutoff. Must be greater than 0, and is clamped to the number of chains at runtime. Default 1000. Changing this from the default requires EXPERIMENTAL_FEATURES : True (VMMC is experimental)."],
    'MOVE_CHAIN_TRANSLATE' : ["float", "Move code 2. Probability of a whole-chain rigid translation move being attempted. Works on chains of any length, including single beads. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CHAIN_ROTATE' : ["float", "Move code 3. Probability of a whole-chain rigid rotation move being attempted (a cardinal 90/180/270 degree rotation about the bead nearest the chain centroid). Needs chains of at least 2 beads: drawn for a single-bead chain it is a null move that is simply rejected, and PIMMS warns at start-up if monomers are present with this move enabled. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CHAIN_PIVOT' : ["float", "Move code 4. Probability of a molecular pivot move being attempted. Pivot moves randomly select a bead on the chain and rotate the SHORTER arm about it, leaving the rest of the chain fixed. Needs chains of at least 3 beads; shorter chains give a null move that is rejected. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_HEAD_PIVOT' : ["float", "Move code 5. Probability of a head pivot move being attempted. Head pivot moves randomly select one of the two ends of a chain and pivot that terminus, but this is almost never worth doing so we recommend setting it to 0. Needs chains of at least 2 beads. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CLUSTER_TRANSLATE' : ["float", "Move code 7. Probability of a cluster translation move being attempted: the connected cluster containing a randomly chosen chain is translated rigidly as a body. Cluster translation moves are relatively expensive, so in general wise to keep this at a low number (0.01 to 0.05). A cluster containing every chain in the system is rejected, so this move does nothing in a fully condensed single-droplet system. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CLUSTER_ROTATE' : ["float", "Move code 8. Probability of a cluster rotation move being attempted: the connected cluster containing a randomly chosen chain is rotated rigidly by a cardinal 90/180/270 degrees. Cluster rotation moves are relatively expensive, so in general wise to keep this at a low number (0.01 to 0.05). Under periodic boundaries the production box must be cubic/square (see DIMENSIONS), and a cluster that winds around the box is rejected. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'CRANKSHAFT_MODE' : ["str (UNIFORM)", "Obsolete and ignored. It used to define how the number of crankshaft sub-moves scaled with chain length. The keyword is still accepted so that old keyfiles keep working, but PIMMS prints a warning and ignores its value: the crankshaft always performs a fixed CRANKSHAFT_SUBSTEPS sub-moves per megamove, independent of chain length (the old UNIFORM mode). Remove it from new keyfiles."],

    'MOVE_SLITHER' : ["float", "Move code 6. Probability of a slither (reptation) megamove being attempted. When selected, every non-frozen chain is slithered SLITHER_SUBSTEPS times - a chain advances forwards or backwards through the lattice like a snake, which efficiently relaxes chain conformations. Works in 2D and 3D. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_PULL' : ["float", "Move code 11. Probability of a pull (cooperative reptation) megamove being attempted. When selected, every non-frozen chain of length >= 3 is pulled PULL_SUBSTEPS times - an interior bead is displaced and the following beads are cooperatively 'pulled' along to restore connectivity, letting chains rearrange in DENSE systems where rigid moves would clash (the chain termini are not moved by this move, so pair it with crankshaft/slither). Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_CTSMMC' : ["float", "Move code 9. Probability of a single-chain TSMMC (Temperature-Switch Monte Carlo) move being attempted. A randomly selected chain is taken on a temperature EXCURSION - heated along a schedule from TEMPERATURE up to TSMMC_JUMP_TEMP and cooled back - to help it escape local energy minima, with a tempered-transitions acceptance that preserves detailed balance. Controlled by the TSMMC_* keywords. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_MULTICHAIN_TSMMC' : ["float", "Move code 10. Probability of a multi-chain TSMMC move being attempted. As MOVE_CTSMMC, but a randomly selected SUBSET of the non-frozen chains (between 1 and about a quarter of them) undergoes the temperature excursion together. Controlled by the TSMMC_* keywords. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_SYSTEM_TSMMC' : ["float", "Move code 12. Probability of a system-wide TSMMC move being attempted. The ENTIRE system undergoes a temperature excursion (heated along a schedule to TSMMC_JUMP_TEMP and cooled back) to help the whole configuration escape local minima. The sub-moves inside the excursion are drawn from your keyfile move mix, except that a nested TSMMC draw is not allowed and is redrawn from the non-TSMMC moves with their fractions renormalised; if no non-TSMMC move is enabled at all they fall back to crankshaft megamoves, which is announced once. Controlled by the TSMMC_* keywords. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_JUMP_AND_RELAX' : ["float", "Move code 13. Probability of a single-chain jump-and-relax move being attempted. A selected chain is relaxed (a crankshaft sub-trajectory of CRANKSHAFT_SUBSTEPS sub-moves), a rigid translation ('jump') is proposed and accepted or rejected on its own Metropolis criterion, then the chain is relaxed again. Each of the three sub-steps preserves the Boltzmann distribution, so the composite move maintains detailed balance. The jump is evaluated with a full from-scratch energy recompute (O(N) in the number of beads), so keep the fraction modest for large systems. Useful for relocating individual chains and letting them settle; for relocation through dense/condensed phases prefer MOVE_VMMC or MOVE_PULL. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'TSMMC_JUMP_TEMP' : ["float", "The peak ('jump') temperature reached during a TSMMC temperature excursion. MUST be greater than the simulation TEMPERATURE (and, in a QUENCH_RUN, greater than both QUENCH_START and QUENCH_END - a ramp that would reach it is rejected at start-up): a TSMMC move heats the selected chain(s)/system from TEMPERATURE up to TSMMC_JUMP_TEMP and back. Ignored if TSMMC_FIXED_OFFSET is set. Used only by the TSMMC moves (MOVE_CTSMMC / MOVE_MULTICHAIN_TSMMC / MOVE_SYSTEM_TSMMC). Default 50."],

    'TSMMC_STEP_MULTIPLIER' : ["int (positive)", "Sets how much sampling happens at EACH temperature point of a TSMMC excursion; must be greater than 0. For MOVE_CTSMMC the number of sub-moves per temperature is TSMMC_STEP_MULTIPLIER x the chain length, for MOVE_MULTICHAIN_TSMMC it is TSMMC_STEP_MULTIPLIER x the number of beads in the selected chains, and for MOVE_SYSTEM_TSMMC it is exactly TSMMC_STEP_MULTIPLIER whole-system moves. Larger values equilibrate the system more thoroughly at each temperature, at the cost of much slower excursions. Used only by the TSMMC moves. Default 50."],

    'TSMMC_NUMBER_OF_POINTS' : ["int (positive)", "Number of temperature points on the heating ramp of a TSMMC excursion; must be greater than 0. The ramp is TSMMC_NUMBER_OF_POINTS equally spaced temperatures rising from just above the simulation temperature to exactly the jump temperature; the full schedule is that ramp, then a short hold at the jump temperature, then the ramp in reverse. More points give a smoother (more gradual, and more expensive) heating/cooling ramp. Used only by the TSMMC moves. Default 20."],

    'TSMMC_INTERPOLATION_MODE' : ["str (LINEAR)", "How the temperature is interpolated between TEMPERATURE and TSMMC_JUMP_TEMP across the excursion schedule. Currently the only supported value is LINEAR (equal temperature increments); anything else is rejected at parse time. Used only by the TSMMC moves. Default LINEAR."],

    'TSMMC_FIXED_OFFSET' : ["float (positive)", "If set, the TSMMC jump temperature is defined RELATIVE to the current simulation temperature as TEMPERATURE + TSMMC_FIXED_OFFSET rather than using the absolute TSMMC_JUMP_TEMP. Must be greater than 0. This is what you want in a QUENCH_RUN, where the jump temperature then tracks the ramp instead of being overtaken by it. The default is False (use the absolute TSMMC_JUMP_TEMP), but False is NOT keyfile syntax - simply omit the keyword. Used only by the TSMMC moves."],

    'ANALYSIS_MODULE' : ["str (path)", "Filepath (relative or absolute, with a leading ~ expanded) to a user-supplied Python module defining an analysis_function(step, lattice) that PIMMS loads and runs during the simulation (the custom-analysis hook; see ANA_CUSTOM). The module is imported and validated at parse time, so a broken module fails immediately rather than part-way through a run, and any exception it raises at runtime is reported as an error in YOUR code. Giving ANALYSIS_MODULE without a positive ANA_CUSTOM is an error (the module would never run). Omit the keyword (the default) to run no custom analysis; an empty value is rejected."],

    'ANA_CUSTOM' : ["int", "Frequency (in steps) at which the user-defined custom analysis function (from ANALYSIS_MODULE) is run. 0 (default) disables custom analysis. Setting it without an ANALYSIS_MODULE prints a warning and does nothing."],

    'ANA_CLUSTER_THRESHOLD' : ["int", "Connected components containing MORE than this many chains get the per-cluster shape/size analysis (CLUSTER_RG/ASPH/AREA/VOL/DEN and radial profiles). The comparison is strict, so the default of 1 skips single chains. Note the size-distribution files (CLUSTERS.dat / NUM_CLUSTERS.dat and the LR_* variants) always include EVERY component regardless of this threshold."]}


# Logical groupings of keywords used to organise the `PIMMS --info` output under
# subheadings (ordered). Every keyword in EXPECTED_KEYWORDS should appear in
# exactly one group; any that do not are shown under "Other" by the CLI.
KEYWORD_GROUPS = [
    ("Core simulation setup (most are required)",
        ['DIMENSIONS', 'PARAMETER_FILE', 'CHAIN', 'TEMPERATURE', 'N_STEPS',
         'EQUILIBRATION', 'SEED', 'HARDWALL']),

    ("System & chain options",
        ['EXTRA_CHAIN', 'CASE_INSENSITIVE_CHAINS', 'LATTICE_TO_ANGSTROMS',
         'AUTOCENTER', 'NON_INTERACTING', 'ANGLES_OFF', 'FREEZE_FILE']),

    # listed in move-code order (1-14), which is the column order of
    # MOVE_FREQS.dat / ACCEPTANCE.dat - see AcceptanceCalculator.MOVE_CODES
    ("Monte Carlo moves (the MOVE_* probabilities must sum to 1.0)",
        ['MOVE_CRANKSHAFT', 'MOVE_CHAIN_TRANSLATE', 'MOVE_CHAIN_ROTATE',
         'MOVE_CHAIN_PIVOT', 'MOVE_HEAD_PIVOT', 'MOVE_SLITHER',
         'MOVE_CLUSTER_TRANSLATE', 'MOVE_CLUSTER_ROTATE', 'MOVE_CTSMMC',
         'MOVE_MULTICHAIN_TSMMC', 'MOVE_PULL', 'MOVE_SYSTEM_TSMMC',
         'MOVE_JUMP_AND_RELAX', 'MOVE_VMMC']),

    ("Move tuning",
        ['CRANKSHAFT_SUBSTEPS', 'CRANKSHAFT_MODE', 'SLITHER_SUBSTEPS',
         'PULL_SUBSTEPS', 'VMMC_MAX_DISPLACEMENT', 'VMMC_MAX_CLUSTER']),

    ("TSMMC (temperature-switch) excursion settings",
        ['TSMMC_JUMP_TEMP', 'TSMMC_STEP_MULTIPLIER', 'TSMMC_NUMBER_OF_POINTS',
         'TSMMC_INTERPOLATION_MODE', 'TSMMC_FIXED_OFFSET']),

    ("Quench / simulated annealing",
        ['QUENCH_RUN', 'QUENCH_FREQ', 'QUENCH_STEPSIZE', 'QUENCH_START',
         'QUENCH_END', 'QUENCH_AS_EQUILIBRATION']),

    ("Output & I/O",
        ['PRINT_FREQ', 'XTC_FREQ', 'EN_FREQ', 'REDUCED_PRINTING', 'SAVE_EQ',
         'SAVE_AT_END', 'TRAJECTORY_PBC_UNWRAP', 'WRITE_CHAIN_TO_CHAINID', 'ENERGY_CHECK']),

    ("Analysis",
        ['ANALYSIS_FREQ', 'ANA_POL', 'ANA_INTSCAL', 'ANA_DISTMAP',
         'ANA_ACCEPTANCE', 'ANA_INTER_RESIDUE', 'ANA_CLUSTER',
         'ANA_CLUSTER_THRESHOLD', 'ANA_RESIDUE_PAIRS', 'ANALYSIS_MODULE',
         'ANA_CUSTOM']),

    ("Restart",
        ['RESTART_FREQ', 'RESTART_FILE', 'RESTART_OVERRIDE_DIMENSIONS',
         'RESTART_OVERRIDE_HARDWALL', 'RESTART_CONTINUE']),

    ("Equilibration options",
        ['RESIZED_EQUILIBRATION', 'EQUILIBRATION_OFFSET']),

    ("Parallelization",
        ['PARALLELIZE', 'PARALLEL_THREADS']),

    ("Experimental features",
        ['EXPERIMENTAL_FEATURES']),
]

    
    

 
    

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
OUTNAME_LOGFILE='log.txt'

OUTPUT_USED_PARAMETER_FILE='parameters_used.prm'
OUTPUT_FULL_ANGLE_POTENTIAL='absolute_energies_of_angles.txt'
OUTPUT_CHAIN_TO_CHAINID='chain_to_chainid.txt'

RESTART_FILENAME='restart.pimms'

# The resolved configuration a run actually used, written at start-up (see
# Simulation.write_effective_keyfile): the keyfile after restart overrides, seed
# generation and chain merging, with a provenance header. Re-parses as a keyfile.
EFFECTIVE_KEYFILE_NAME='keyfile_used.kf'


##
## Manifest of the analysis output files a simulation can write
##
## Analysis output files are created lazily - every writer opens its file in
## append mode (or, for the end-of-run files, write mode) at the moment it has
## something to put in it - so a file exists if and only if the run wrote at
## least one row to it. Nothing is pre-created, which means nothing can be
## predicted at start-up either: whether CLUSTER_RADIAL_DENSITY_PROFILE.dat
## appears, for instance, depends on whether any cluster ever reached the bead
## threshold.
##
## What start-up DOES need is the reverse of creation. A run that will not write
## a given file must still remove any copy of it left by a previous run in the
## same directory, otherwise the second run silently inherits the first run's
## data. The manifest below is the single canonical list that drives that
## removal (see Simulation.startup_analysis).
##
## The entries are the NAMES of the filename constants above rather than the
## filenames themselves, so that redirecting a constant (a test pointing an
## output at a temporary directory, say) redirects the manifest with it.
##

ANALYSIS_OUTPUT_NAMES = (

    # core state, throughput and move bookkeeping (appended during the run)
    'OUTNAME_ENERGY',
    'OUTNAME_PERFORMANCE',
    'QUENCHFILE_NAME',
    'OUTNAME_MOVES',
    'OUTNAME_ACCEPTANCE',
    'OUTNAME_TOTAL_MOVES',

    # per-chain (polymeric) analysis, appended during the run
    'OUTNAME_RG',
    'OUTNAME_ASPH',
    'OUTNAME_E2E',
    'OUTNAME_R2R',

    # short-range cluster analysis, appended during the run
    'OUTNAME_CLUSTERS',
    'OUTNAME_NUM_CLUSTERS',
    'OUTNAME_CLUSTER_RG',
    'OUTNAME_CLUSTER_ASPH',
    'OUTNAME_CLUSTER_VOL',
    'OUTNAME_CLUSTER_AREA',
    'OUTNAME_CLUSTER_DENSITY',
    'OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE',

    # long-range cluster analysis, appended during the run
    'OUTNAME_LR_CLUSTERS',
    'OUTNAME_NUM_LR_CLUSTERS',
    'OUTNAME_LR_CLUSTER_RG',
    'OUTNAME_LR_CLUSTER_ASPH',
    'OUTNAME_LR_CLUSTER_VOL',
    'OUTNAME_LR_CLUSTER_AREA',
    'OUTNAME_LR_CLUSTER_DENSITY',
    'OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE',

    # accumulated over the run and written once at the end
    'OUTNAME_INTERNAL_SCALING',
    'OUTNAME_INTERNAL_SCALING_SQUARED',
    'OUTNAME_SCALING_INFORMATION',
    'OUTNAME_DMAP',
)


## The subset of ANALYSIS_OUTPUT_NAMES which a multi-chain-type run writes once
## per chain type, with a CHAIN_<type>_ prefix on the basename. The cluster
## composition files are written in addition to the unprefixed files, while the
## internal-scaling/distance-map files replace them.
PER_CHAIN_TYPE_OUTPUT_NAMES = (
    'OUTNAME_CLUSTERS',
    'OUTNAME_LR_CLUSTERS',
    'OUTNAME_INTERNAL_SCALING',
    'OUTNAME_INTERNAL_SCALING_SQUARED',
    'OUTNAME_SCALING_INFORMATION',
    'OUTNAME_DMAP',
)


def analysis_output_files(chain_types=None):
    """
    Return every analysis output file a simulation could write.

    This resolves the manifest (``ANALYSIS_OUTPUT_NAMES``) into actual
    filenames, optionally adding the ``CHAIN_<type>_`` prefixed variants for a
    set of chain types. It is the one place that knows the full set of analysis
    outputs, and is what start-up uses to delete a previous run's files.

    Note that this is the set of files a run *could* write, not the set it will
    write - output files are created lazily, so which of them actually appear
    depends on which analyses fire.

    Parameters
    ----------
    chain_types : iterable or None, optional
        Chain types whose ``CHAIN_<type>_`` prefixed variants should be
        included. If ``None`` (the default) only the unprefixed names are
        returned.

    Returns
    -------
    list
        List of filenames, unprefixed names first and then, if requested, the
        per-chain-type variants in chain-type order.
    """

    filenames = [globals()[name] for name in ANALYSIS_OUTPUT_NAMES]

    if chain_types is not None:
        for chain_type in chain_types:
            for name in PER_CHAIN_TYPE_OUTPUT_NAMES:
                base = globals()[name]
                directory = os.path.dirname(base)
                prefixed = 'CHAIN_%s_%s' % (chain_type, os.path.basename(base))
                if directory:
                    prefixed = os.path.join(directory, prefixed)
                filenames.append(prefixed)

    return filenames

