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
# a new chain into the lattice. Default is 20, although perhaps you
# might want to change this for some reason?
CHAIN_INIT_ATTEMPTS = 20

# Debug-mode flag. Nothing in PIMMS currently reads it; it is kept only so that
# code importing it does not break.
DEBUG = False

# Numerator of the inverse temperature used by every acceptance test:
# beta = INVTEMP_FACTOR / TEMPERATURE. 1.0 means Boltzmann's constant is 1, i.e.
# TEMPERATURE is in the same (reduced) units as the parameter-file energies.
INVTEMP_FACTOR = 1.0

# Number of extra rungs a TSMMC excursion holds at the jump temperature, between
# the heating ramp and the cooling ramp, so the full schedule is
# 2 x TSMMC_NUMBER_OF_POINTS + TOP_TEMP rungs.
TOP_TEMP = 10 # 

# The fixed output range of the kernels' splitmix64 draws (mega_crank.PRNG_MAX), and
# the modulus used to seed the reference kernel's PRNG from the keyfile SEED (numpy's
# global generator is seeded with SEED modulo 2**32).
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

# minimum number of beads a cluster needs (>=) for its radial density profile to
# be calculated
RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD = 27

## NB: THIS VALUE CAN BE CHANGED. To reduce PIMMS' memory footprint you
## you can change this to np.intxxx where xxx could be 16, 32 or 64. In principle
## it could be 8 but this would be quite limiting in terms of number of unique
## beads that could be used (=256, maybe fine?). NOTE that if you change
## this value you must change the corresponding CYTHON config in cython_config.pxd
NP_INT_TYPE = np.int32


## ------------------------------------------------------------------------
##                       INPUT SANITY LIMITS
## ------------------------------------------------------------------------
##
## Limits the keyfile parser applies before anything is allocated or written.
## Each is meant to be out of reach of a legitimate run: a refusal here replaces
## a crash, a multi-gigabyte allocation or a truncated output file later on. The
## keyword descriptions below quote these numbers, so change both together.

# The megamoves (crankshaft, slither, pull, single- and multi-chain TSMMC) build
# one "selector" array per megamove holding one entry per sub-move, and the
# compiled kernels index it with a C int. A selector longer than 2**31 - 1 can
# therefore not be run at all, whatever the machine. The four keywords that set
# its length (CRANKSHAFT_SUBSTEPS, SLITHER_SUBSTEPS, PULL_SUBSTEPS,
# TSMMC_STEP_MULTIPLIER) are refused above this on their own, and so is any
# enabled move whose selector (the keyword times the number of chains or beads it
# multiplies) would exceed it.
MAX_SUBMOVE_SELECTOR_LENGTH = 2**31 - 1

# A selector entry is a 64-bit integer, so a selector costs 8 bytes per sub-move
# for the length of the megamove. A selector that would not fit in physical
# memory is refused; one above SUBMOVE_SELECTOR_WARN_LENGTH entries (0.8 GB,
# rebuilt on every megamove) is allowed with a warning, since it is merely
# expensive. For scale, a production CRANKSHAFT_SUBSTEPS is 20-50 thousand.
SUBMOVE_SELECTOR_BYTES_PER_ENTRY = 8
SUBMOVE_SELECTOR_WARN_LENGTH = 10**8

# A TSMMC excursion visits 2 x TSMMC_NUMBER_OF_POINTS + TOP_TEMP temperatures, and
# those temperatures are rounded to five decimal places, so a ramp of more than
# 10**5 points per unit of temperature repeats rungs. A million points is beyond
# any ramp that could be useful and keeps the schedule arrays at a few tens of MB.
MAX_TSMMC_NUMBER_OF_POINTS = 10**6

# The parallel kernels split the box into at most 64 blocks (4 per axis), so no
# more than 64 threads ever have work; 1024 is far beyond that and beyond the
# core count of any single node, while a larger request makes the OpenMP runtime
# try to create that many threads and abort the process when it cannot.
MAX_PARALLEL_THREADS = 1024

# The lattice is two integer grids (occupancy and bead type) with one entry per
# site: 8 bytes per site at the default NP_INT_TYPE. A box whose grids alone
# exceed the machine's physical memory cannot be allocated and is refused; where
# the physical memory cannot be read the ceiling is 1 TiB, more than any node
# PIMMS has been run on. Above half the limit the box is allowed with a warning.
LATTICE_BYTES_PER_SITE = 2 * np.dtype(NP_INT_TYPE).itemsize
LATTICE_MEMORY_FALLBACK_CEILING = 2**40
LATTICE_MEMORY_WARN_FRACTION = 0.5

# START.pdb holds each coordinate in an 8-column field with three decimals
# (%8.3f), which ends at 9999.999 Angstroms: a box of 10000 Angstroms or more
# along any axis has sites that cannot be written.
PDB_BOX_LIMIT_ANGSTROMS = 10000.0

# traj.xtc is written at mdtraj's default XTC precision of 1000, i.e. coordinates
# are stored in units of 0.001 nm = 0.01 Angstroms. A lattice spacing below one
# such unit puts neighbouring sites on the same stored coordinate (the trajectory
# is destroyed); below ten units the rounding error (up to 0.005 Angstroms) is
# more than 5 % of a lattice step, which is allowed with a warning.
MIN_LATTICE_TO_ANGSTROMS = 0.01
LATTICE_TO_ANGSTROMS_WARN_BELOW = 0.1


## ------------------------------------------------------------------------
##                                KEYWORDS
## ------------------------------------------------------------------------

# list of ALL valid keywords. Every keyword here needs a DEFAULTS entry (unless it
# is in REQUIRED_KEYWORDS), a KEYWORDS_DESCRIPTION entry and a place in
# KEYWORD_GROUPS (the check after DEFAULTS enforces the first).
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
DEFAULTS['TEMPERATURE']                 = 'N/A'     # placeholder only: TEMPERATURE is required, even with a RESTART_FILE


# major setup things
DEFAULTS['RESIZED_EQUILIBRATION']       = False
DEFAULTS['EQUILIBRATION_OFFSET']        = False     # False = centre the equilibration box in the production box
DEFAULTS['HARDWALL']                    = False     
DEFAULTS['EXPERIMENTAL_FEATURES']       = False     # This must be set to true to use experimental features
DEFAULTS['LATTICE_TO_ANGSTROMS']        = 3.65      # note: in 0.1.34 we update this to 3.65 from 4 as used previously this is a breaking default change  
DEFAULTS['NON_INTERACTING']             = False     # use interactions 
DEFAULTS['ANGLES_OFF']                  = False     # use angles
DEFAULTS['CASE_INSENSITIVE_CHAINS']     = True      # means we cast chains to upper case if set to True
DEFAULTS['AUTOCENTER']                  = False     # means we do not by default centre a single chain in written frames

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

# parallelization of the crankshaft, slither and pull megamoves
DEFAULTS['PARALLELIZE']        = False  # By default use the (serial) optimized kernel
DEFAULTS['PARALLEL_THREADS']   = 0      # 0 => auto (use all available CPU cores)


# FINALLY we do some sanity checking here

for k in EXPECTED_KEYWORDS:
    if k not in DEFAULTS:
        if k not in REQUIRED_KEYWORDS:
            raise Exception(f'No default value set for {k} - this is a bug!')


KEYWORDS_DESCRIPTION = {
    'DIMENSIONS': ['int (2 or 3 values, e.g. A B or A B C)',
                   '[REQUIRED] - Size of the simulation box (in lattice units). Providing 2 values runs a 2D simulation, 3 values a 3D simulation. The axes need NOT be equal: non-cubic/non-square boxes (e.g. 10 20 40) are fully supported with either HARDWALL or periodic boundaries. Every axis must be at least 7 lattice sites (the smallest box that supports the super-long-range interaction shell). The lattice costs 8 bytes of memory per site (two integer grids) before anything else is allocated, so a box whose grids would exceed the physical memory of the machine (or 1 TiB where that cannot be read) is refused at start-up and one above half of it gets a warning; for scale a 400 x 400 x 400 box is 0.5 GB and a 1000 x 1000 x 1000 box is 8 GB. The longest axis times LATTICE_TO_ANGSTROMS must also stay below 10000 Angstroms (see LATTICE_TO_ANGSTROMS). The only other restriction is that cluster-rotation moves (MOVE_CLUSTER_ROTATE) cannot be combined with a non-cubic box under periodic boundaries, because a 90-degree rigid rotation is only an energy-preserving symmetry of a cube/square (or of any box under HARDWALL, where there is no periodic wrapping).'],
    'LATTICE_TO_ANGSTROMS' : ["float (positive)", "Conversion factor (default 3.65) for converting lattice units to Angstroms when writing the START.pdb topology and traj.xtc trajectory, and to compute the solute concentration (M) reported at start-up for a 3D box. It has a validated range, set by the two file formats. It must be at least 0.01: traj.xtc stores coordinates in units of 0.001 nm (0.01 Angstroms), so a smaller spacing puts neighbouring lattice sites on the same stored coordinate; below 0.1 the run is allowed with a warning, because the rounding is then more than 5 % of a lattice step. And the longest box axis times LATTICE_TO_ANGSTROMS must be below 10000 Angstroms, because START.pdb holds each coordinate in an 8-column field that ends at 9999.999 (at the default spacing that is a box axis of 2739 sites); the check uses the box the run ends up with, after any RESTART_OVERRIDE_DIMENSIONS. Both are refused when the keyfile is read, before any file is written. Within that range this is purely cosmetic - it sets the bead spacing seen in a viewer/analysis (mdtraj reports nm, i.e. lattice_units x LATTICE_TO_ANGSTROMS x 0.1) and has NO effect on the simulation itself or its energetics; every analysis .dat file stays in lattice units. It is also the value lemonade needs to reconstruct the lattice from a trajectory, so keep a record of it (it is recorded in keyfile_used.kf)."],
    'CHAIN' : ["See description", "[REQUIRED] - One of the few keywords that can appear multiple times (the others are EXTRA_CHAIN and ANA_RESIDUE_PAIRS), the 'CHAIN' keyword defines a specific polymer chain and the number of that chain that will exist in the simulation. The format should be \n\nCHAIN : N  {CHAIN IDENTITY}\n\nWhere 'N' defines the number of copies of the chain and '{CHAIN IDENTITY}' gives the polymer sequence, one character per bead, where every character must be a bead type defined in the parameter file (the run is refused otherwise). As an example\n\nCHAIN : 20 QQQQQQQQQQ\n\nWould give 20 poly-glutamine polymers. N must be an integer of 1 or more, and the sequence may not contain the character 0, which is reserved for solvent. Each CHAIN line is a distinct chain type (even if two lines share a sequence), numbered from 0 in the order the lines appear; chainIDs are numbered from 1 across the whole system, in the same order. Each chain type gets its own chain identifier in START.pdb; there are 62 of them (A-Z, a-z, 0-9), so with more than 62 chain types PIMMS warns at start-up that the types past the 62nd all share the identifier 9. This keyword is required UNLESS a RESTART_FILE is provided, in which case the chains are taken from the restart file and CHAIN may be omitted (any CHAIN lines are then ignored). Chains are placed at random, except that a system of exactly one chain grows that chain from the middle of the box. A system with more beads than the box has sites is refused before anything is placed, and a chain that still cannot be placed after 20 attempts aborts the run with an overcrowded-lattice error (a lone chain that is long relative to the box can trap itself while growing from the centre), so keep the occupied volume fraction reported at start-up sensible. In later versions of PIMMS we will be updating this to allow the reading of keyfiles that use three-letter codes."],
    'CASE_INSENSITIVE_CHAINS' : ["bool", "Boolean flag which, if set to False, means that chain sequence is case sensitive. By default, this is True, which means that upon reading a keyfile, CHAIN and EXTRA_CHAIN sequences are converted to upper case, so a chain written with 'a' becomes 'A' and a lower-case bead type can never appear in a chain. Sometimes you may want those extra unique beads, in which case setting this to False is useful. Note the parameter file itself is never case-folded. Every bead type used in a chain must also be defined in the parameter file."],
    'TEMPERATURE': ["float (positive)","[REQUIRED] - Simulation temperature; must be a positive number greater than 0. In general a temperature between 10 and 200 is appropriate for the energy scales typical of PIMMS parameter files. Higher temperatures sample more expanded/disordered states; lower temperatures favour collapse/assembly. Note that PIMMS works in reduced units with k = 1, so the Boltzmann factor is exp(-dE/TEMPERATURE) and the temperature is on the same scale as the parameter-file energies. In a QUENCH_RUN this value is overwritten by QUENCH_START (a warning is printed if the two differ), so use QUENCH_START / QUENCH_END to set the ramp. It is required even with a RESTART_FILE, and with RESTART_CONTINUE (outside a QUENCH_RUN) it must equal the temperature the restart file was written at, or the run is refused."],
    'N_STEPS' : ["int (positive)", "[REQUIRED] - Total number of outer-loop steps to run (including the EQUILIBRATION steps). Must be a positive integer, and must be larger than EQUILIBRATION. With RESTART_CONTINUE it is still the total length of the whole run, counted from step 0, not the length of the resumed segment. Note that one step is typically a great deal of Monte Carlo work: each crankshaft step performs CRANKSHAFT_SUBSTEPS single-bead sub-moves in total (in a serial run spread at random over all the non-frozen beads; under PARALLELIZE over the beads inside the block interiors), and the slither/pull/TSMMC moves are likewise 'megamoves', so the true number of accept/reject operations is far larger than N_STEPS (see TOTAL_MOVES.dat)."],
    'PARAMETER_FILE': ["string", "[REQUIRED] - Filepath (relative or absolute, with a leading ~ expanded; a relative path is resolved from the directory PIMMS is run in, not the keyfile's directory) to the parameter file defining the interaction energies and angle penalties. The simulation fails if it does not exist. A verbatim copy of the file, under a timestamp header, is written to parameters_used.prm at startup."],
    'EQUILIBRATION': ["int (0 or more)", "[REQUIRED] - Number of initial steps treated as equilibration. Must be 0 or larger and smaller than N_STEPS. In a QUENCH_RUN with QUENCH_AS_EQUILIBRATION : True the value you give is replaced by the length of the quench. Step number EQUILIBRATION is the LAST equilibration step; production begins at EQUILIBRATION + 1. During equilibration no ANA_* analysis output and no restart snapshot is written (ENERGY.dat and PERFORMANCE.dat are still written throughout, and trajectory frames are saved if SAVE_EQ is True). Choose this large enough that ENERGY.dat has plateaued before production begins."],
    'SAVE_EQ' : ["bool", "Boolean (true or false) that determines whether PIMMS saves trajectory frames for the equilibration steps of a simulation. If set to False, PIMMS begins to save your trajectory frames *after* the equilibration steps have completed. Note frame 0 of traj.xtc is always the starting configuration, regardless of this setting (for a RESIZED_EQUILIBRATION run it is the post-resize configuration at step EQUILIBRATION, and eq_START.pdb / eq_traj.xtc exist only when SAVE_EQ is True). A RESIZED_EQUILIBRATION run deletes any traj.xtc / START.pdb left by an earlier run at start-up (they appear only once the box is resized), and with SAVE_EQ False it also deletes stale eq_traj.xtc / eq_START.pdb."],
    'RESIZED_EQUILIBRATION': ['int (2 or 3 values, e.g. A B or A B C)', "Defines a smaller box to use during equilibration; at the end of equilibration (before the move of step EQUILIBRATION + 1) the box is grown to the full DIMENSIONS, with the equilibration box centred in the new box (the configuration keeps its place inside it), or placed using EQUILIBRATION_OFFSET. Useful for condensing/assembling a system at high effective concentration before expanding to the production box. Must have the same number of values as DIMENSIONS, must be <= DIMENSIONS in every dimension, and every axis must satisfy the same >= 7 floor as DIMENSIONS (this box IS simulated). The equilibration phase is always run under hardwall boundaries (forced internally, so a system is never resized while chains straddle a periodic face); your production HARDWALL setting takes over once the box has grown. While the resized box is in use the trajectory is written to eq_START.pdb / eq_traj.xtc rather than START.pdb / traj.xtc, and those files are only written if SAVE_EQ is True; START.pdb / traj.xtc are written from the resize onward, and any copies an earlier run left in the directory are deleted at start-up. Setting EQUILIBRATION to 0 deactivates this keyword (and EQUILIBRATION_OFFSET with it) with a warning. Incompatible with RESTART_OVERRIDE_DIMENSIONS, with RESTART_CONTINUE and with periodic (non-hardwall) restart files; with a hardwall restart file the restart box must be <= the RESIZED_EQUILIBRATION box (a smaller one is centred in it). See also EQUILIBRATION_OFFSET."],
    'EQUILIBRATION_OFFSET': ['int (2 or 3 values, e.g. A B or A B C)', "Defines the offset of the equilibration box relative to the full simulation box, i.e. where the small box sits inside the production box when it is grown (without it the small box is centred). Requires RESIZED_EQUILIBRATION to be set, and must have the same number of values. Every value must be >= 0, and for each dimension EQUILIBRATION_OFFSET + RESIZED_EQUILIBRATION MUST be <= DIMENSIONS. It is deactivated, together with RESIZED_EQUILIBRATION, when EQUILIBRATION is 0."],
    'HARDWALL' :["bool", "Boolean flag set to True or False that defines whether a hardwall boundary is used or not. By default (False) periodic boundary conditions are used, so the box wraps in every dimension. If set to True the box has hard walls: no bead ever sits outside the box, no bond ever crosses a wall, and no chain is ever left straddling a face. Note that this constrains the states, not the paths between them: a whole-chain rigid translation (MOVE_CHAIN_TRANSLATE, and the jump inside MOVE_JUMP_AND_RELAX) draws its offset uniformly over the box and wraps, so a chain sitting against one wall can be relocated in one move to the opposite wall. Every state it visits is a legal confined state and the proposal is symmetric, so the sampled ensemble is exactly the confined Boltzmann distribution; but the move is a relocation, not a physical passage through the wall, and a trajectory will show the chain jumping across the box. If you need a trajectory in which no chain ever jumps through a wall, leave MOVE_CHAIN_TRANSLATE and MOVE_JUMP_AND_RELAX at 0 (cluster translation and VMMC never wrap under a hard wall: a proposal that would carry any bead out of the box is rejected). The walls are energetically solvent - a bead beside a wall picks up its bead-solvent energy for each out-of-box neighbour site - so a wall excludes volume but adds no special surface energy of its own. Nothing is read across a wall either: contacts and clusters never connect two beads through a face, and coordinates are plain Cartesian because there is no periodic image to reconstruct."],
    'NON_INTERACTING' : ["bool", "Boolean flag set to True or False that defines if a non-interacting simulation should be performed or not. If set to true, all bead-bead (SR/LR/SLR) interaction and solvation energies are set to zero. Angle penalties are NOT affected - combine with ANGLES_OFF : True for a fully ideal excluded-volume reference state."],
    'ANGLES_OFF' : ["bool", "Boolean flag set to True or False that defines if angle potentials are to be used or not. If set to False (or not set), angles from the parameter file will be used, and EVERY bead type with interaction energies must then carry an ANGLE_PENALTY (or ANGLE_PENALTY_T_NORM) line or the parameter file is rejected. If set to True, angles are ignored (every penalty is forced to zero, which is announced at start-up unless REDUCED_PRINTING is on) and parameter files do not need to define angles."],
    'EXPERIMENTAL_FEATURES' : ["bool", "Boolean flag set to True or False that defines if experimental/non-supported keywords and features are allowed. As of this release the only gated keywords are MOVE_VMMC, VMMC_MAX_DISPLACEMENT and VMMC_MAX_CLUSTER; setting any of them away from its default without EXPERIMENTAL_FEATURES : True is an error. STRONGLY recommend leaving this as False, and NONE of the features/behaviors allowed here are guaranteed to work."],
    'SEED' : ["int (positive)", 'Random seed; must be a positive integer. If not set, a random seed is generated (and announced at start-up), but if provided it ensures perfect simulation reproducibility (an identical keyfile + parameter file + seed reproduces the trajectory bit-for-bit on the same platform and PIMMS version). The seed, given or generated, is recorded in keyfile_used.kf. A run started from a RESTART_FILE begins a new random stream from this seed; SEED must NOT be given with RESTART_CONTINUE, which restores the generator state saved in the restart file instead (the run is refused if both are set).'],
    'PRINT_FREQ' : ["int (positive)", "Frequency (in steps) with which a status line (step number, percentage complete and current energy) is printed to STDOUT. Must be greater than 0. The same line is also printed on every step that saves a trajectory frame (only on the 10 % marks of the run under REDUCED_PRINTING) and on every ENERGY_CHECK step. Throughput and the estimated time remaining are not printed: they are written to PERFORMANCE.dat and log.txt at every 5 % of the run (plus once at step 20), whatever PRINT_FREQ is."],
    'XTC_FREQ' : ["int (positive)", 'Frequency (in steps, default 1000) with which a trajectory frame is written to traj.xtc (with START.pdb as the topology). Must be greater than 0. Equilibration frames are only written if SAVE_EQ is True; if SAVE_AT_END is True the trajectory is buffered in memory and written once at the end. An XTC_FREQ you set that is larger than N_STEPS never fires, so traj.xtc holds only the starting configuration; PIMMS prints a warning.'],
    'EN_FREQ' : ["int (positive)", "Frequency (in steps) with which the instantaneous potential energy is appended to ENERGY.dat (tab-separated: step, energy). Must be greater than 0. ENERGY.dat is written throughout the run, including during equilibration, though as with every output file it appears only once a row has actually been written to it (so a run shorter than EN_FREQ produces no ENERGY.dat; if you set EN_FREQ larger than N_STEPS yourself PIMMS prints a warning)."],
    'ANALYSIS_FREQ' : ["int", "Master control parameter that sets the default frequency for the six analyses ANA_POL, ANA_INTSCAL, ANA_DISTMAP, ANA_ACCEPTANCE, ANA_INTER_RESIDUE and ANA_CLUSTER when their own keyword is not given (ANA_CUSTOM does not inherit it). Analyses fire on steps that are multiples of their frequency, and only in production (after step EQUILIBRATION). Setting any of these frequencies, or ANA_CUSTOM or ENERGY_CHECK, to 0 or less disables it for the run; setting ANALYSIS_FREQ itself to 0 disables every one of the six not given its own positive frequency. This 0-means-disabled convention applies ONLY to the analysis frequencies and ENERGY_CHECK: PRINT_FREQ, XTC_FREQ, EN_FREQ and RESTART_FREQ must all be positive (0 is refused)."],
    'ANA_POL' : ["int", "Frequency with which single-chain polymeric analysis is performed. Writes per-chain radius of gyration to RG.dat, asphericity to ASPH.dat and end-to-end distance to END_TO_END_DIST.dat; each row is the step number followed by one column per chain, in ascending chainID order. Defaults to ANALYSIS_FREQ. Columns are tab-separated, distances in lattice units; every chain is included (a single-bead chain reads 0) and each is measured on the whole, bond-walked chain, so a chain crossing a periodic face is not torn."],
    'ANA_INTSCAL' : ["int", "Frequency with which internal-scaling analysis is performed (defaults to ANALYSIS_FREQ). Accumulates the mean inter-residue distance as a function of sequence separation s over the production steps and, at the end of the run, writes INTSCAL.dat (s and the mean distance, lattice units), INTSCAL_SQUARED.dat (s and the mean squared distance) and SCALING_INFORMATION.dat (one row per chain: fitted scaling exponent nu and prefactor R0, or -1 -1 for a chain of fewer than 26 beads, too short to fit). A chain type of single beads has no separations and writes no INTSCAL.dat or INTSCAL_SQUARED.dat. If no production step is a multiple of the frequency (possible only when it exceeds the production length), nothing is written and a warning is printed. The mean DISTANCE_MAP.dat is controlled separately by ANA_DISTMAP. For multi-component systems these are written per chain type as CHAIN_<TYPE>_* (where <TYPE> is the 0-based chain-type index, in keyfile CHAIN-line order, or the restart file's types for a restart run) INSTEAD OF the unprefixed files - the unprefixed names are used only when the system has a single chain type."],
    'ANA_DISTMAP' : ["int", "Frequency with which the mean inter-residue distance map is accumulated. The seqlen x seqlen mean distance matrix is written to DISTANCE_MAP.dat at the end of the run (per chain type as CHAIN_<TYPE>_DISTANCE_MAP.dat). Defaults to ANALYSIS_FREQ. The matrix is tab-separated, in lattice units, averaged over every chain of the type and every production sample; the per-type files replace DISTANCE_MAP.dat, and nothing is written (with a warning) if the map was never sampled. Memory and file size grow with the square of the chain length (8 bytes x seqlen x seqlen per chain held in memory), which matters for chains of several thousand beads; a warning is logged when the estimate exceeds 1 GiB of memory or 256 MB of file."],
    'ANA_ACCEPTANCE' : ["int", "Frequency with which move statistics are written (defaults to ANALYSIS_FREQ): MOVE_FREQS.dat (attempted moves per move code), ACCEPTANCE.dat (accepted moves per move code) and TOTAL_MOVES.dat (the step and the total number of individual accept/reject events of every kind). All three are RUNNING TOTALS since the start of the run (equilibration included, though rows are written only in production; a RESTART_CONTINUE segment counts from its checkpoint), not per-interval counts. Each MOVE_FREQS/ACCEPTANCE row is the step number followed by 14 tab-separated columns, one per move code 1-14 (each MOVE_* keyword description states its code). For the crankshaft, slither and pull megamoves (codes 1, 6 and 11) the columns count individual sub-moves - the attempts the kernels actually made, so a parallel sweep that found nothing movable adds 0 - while every other code counts one per move; the TSMMC codes 9, 10 and 12 count whole excursions, whose sub-moves (like the relaxations inside a jump-and-relax move) appear only in TOTAL_MOVES.dat. Divide ACCEPTANCE by MOVE_FREQS to get the per-move acceptance ratio."],
    'ANA_INTER_RESIDUE' : ["int", "Frequency with which inter-residue distance analysis is performed and appended to RES_TO_RES_DIST.dat. This only makes sense if ANA_RESIDUE_PAIRS has a pair of residues defined - with no pair defined there is nothing to measure and no RES_TO_RES_DIST.dat is written; the distance is computed for EVERY chain, so all chains must be long enough to contain the pair (this is checked at start-up). Defaults to ANALYSIS_FREQ. Each line is the step, the two residue indices, then the distance (lattice units, whole chain) in every chain in ascending chainID order, tab-separated."],
    'ANA_CLUSTER' : ["int", "Frequency with which cluster analysis is performed (defaults to ANALYSIS_FREQ). Identifies short-range clusters (chains in Chebyshev-1 contact) and long-range clusters (chains joined by any short-range contact, whatever its energy, or by a Chebyshev-2 / Chebyshev-3 pair whose LR / SLR table entry is nonzero), and writes their size distributions (CLUSTERS.dat and NUM_CLUSTERS.dat, plus LR_CLUSTERS.dat and NUM_LR_CLUSTERS.dat; a multi-component system also writes CHAIN_<TYPE>_CLUSTERS.dat and CHAIN_<TYPE>_LR_CLUSTERS.dat, giving the fraction of each cluster's chains that are of that type). For clusters larger than ANA_CLUSTER_THRESHOLD it also writes radius of gyration, asphericity and convex-hull volume, surface area and density (CLUSTER_RG/ASPH/AREA/VOL/DEN.dat and the LR_CLUSTER equivalents) and, for those of at least 27 beads, a radial density profile (CLUSTER_RADIAL_DENSITY_PROFILE.dat and its LR_CLUSTER twin). A cluster that percolates the periodic box has no defined shape: it is written as nan (Rg, asphericity) and -1 (volume, area, density) with no radial profile, and a warning is logged. Its cost is driven by the number of chains and by the size of the largest cluster (it walks the connected components of the whole system and then computes convex hulls and radial profiles for each), so keep its frequency low for large, many-chain or strongly condensed systems."],
    'ANA_RESIDUE_PAIRS' : ['int (2 values)', "Two integers used to define a pair of residues, the distance between which is then calculated every ANA_INTER_RESIDUE steps. May be repeated to monitor several pairs (one RES_TO_RES_DIST.dat line per pair per recorded step). Indexing occurs from 0 (i.e., the first residue is 0) and negative indices are rejected; the two values are stored in ascending order, so the order you write them in does not matter. A pair naming the same residue twice (always distance 0) and a pair given more than once (duplicate rows) are both accepted with a warning. Note that at present, inter-residue distances are calculated for EVERY chain, so the larger index must be inside every chain in the system (including chains from a RESTART_FILE or EXTRA_CHAIN) - if it is not, the run is rejected at start-up."],
    'AUTOCENTER' : ["bool", "Boolean flag which, if set to True and you are simulating a SINGLE chain, re-centres that chain in the middle of the box in every WRITTEN frame and in the START.pdb topology. This is a visualisation convenience only - the simulation, its energetics and every analysis file are unaffected. Under HARDWALL the chain is centred only as far as the walls allow, so no bead is ever written outside the box. With more than one chain the keyword has no effect, and a warning says so. Default = False."],
    'REDUCED_PRINTING' : ["bool", "Boolean flag which, if set to True, silences the chatty per-event STDOUT messages - the per-checkpoint 'writing restart file' line, the per-move clash/cluster-resize/multichain-rearrangement notes, the quench temperature updates, the system-wide TSMMC start and accept/reject lines, and at start-up the per-pair NON_INTERACTING notes and the ANGLES_OFF and angle-rounding warnings - and cuts the status line and 'Saving coordinates' remark printed with each trajectory frame to the 10 % marks of the run. The PRINT_FREQ status line, the energy checks, the start-up keyfile summary and everything written to file are unaffected. Useful for long runs whose STDOUT is being captured to a file. Default = False."],
    'SAVE_AT_END' : ["bool", "Boolean flag which, if set to True, holds the trajectory in memory and writes traj.xtc only at the very end of the run (START.pdb and a one-frame traj.xtc holding the starting configuration are written at start-up). This avoids per-frame disk writes, but the whole trajectory is buffered, so memory grows with the number of frames (N_STEPS / XTC_FREQ) times the number of beads (12 bytes per bead per frame in 3D). If the run stops early on an error, a failed ENERGY_CHECK, Ctrl-C or a SIGTERM from a scheduler, the buffered frames are still written to disk; only a hard kill (SIGKILL, power loss) loses them and leaves frame 0. Default = False."],
    'TRAJECTORY_PBC_UNWRAP' : ["bool", "Boolean flag (default False). Under periodic boundaries a chain that crosses a box face is stored split across the two faces, which looks broken in a viewer. If set to True, PIMMS makes every chain WHOLE before writing each trajectory frame (and the START.pdb topology): each chain is shifted into a single periodic image, effectively extending the lattice beyond the box in x/y/z as needed, so no chain is torn across a boundary. This is purely a visualisation convenience - it does not affect the simulation or its energetics, and coordinates may fall outside the box (the unit cell is unchanged). Has no effect with HARDWALL (chains never cross a boundary). Default = False."],
    'WRITE_CHAIN_TO_CHAINID': ["bool", "Boolean flag which, if set to True, writes chain_to_chainid.txt: one tab-separated line per chain giving its chainID, its length and its sequence. Useful for working out which chainIDs to list in a FREEZE_FILE, and for mapping the per-chain columns of RG.dat / ASPH.dat / END_TO_END_DIST.dat back to sequences. Default = False."],
    'FREEZE_FILE' : ["string", "Filepath (relative or absolute, with a leading ~ expanded; a relative path is resolved from the directory PIMMS is run in) to a freeze file; the simulation fails if it does not exist, and an empty value is rejected. The freeze file is a plain-text file listing chainIDs to hold fixed for the whole run, one or more lines of the form 'C <id> <id> ...'. Lines starting with # are comments. Naming a chainID that is not in the system, or a line not starting with C, is an error (bead-level B lines are not implemented). Frozen chains never move but still contribute to the energy (other chains feel them), and the collective moves (cluster translate/rotate, VMMC) reject any move whose cluster would contain a frozen chain. Use WRITE_CHAIN_TO_CHAINID to discover chainIDs (they are numbered from 1). Freezing every chain in the system is refused at start-up, since nothing could then change. Frozen chains are honoured by both the serial and the parallel (PARALLELIZE) move kernels, so freezing and parallelization can be used together."],
    'PARALLELIZE' : ["bool", "Boolean flag (True/False) which, if set to True, runs the crankshaft (MOVE_CRANKSHAFT), slither (MOVE_SLITHER) and pull (MOVE_PULL) moves on multi-threaded checkerboard kernels instead of the serial kernels, including when they run as the sub-moves of a system-wide TSMMC excursion (MOVE_SYSTEM_TSMMC); every other move stays serial, including the crankshaft relaxations inside MOVE_CTSMMC, MOVE_MULTICHAIN_TSMMC and MOVE_JUMP_AND_RELAX. Works in both 2D and 3D. Beneficial for large, spatially dispersed systems; gives little benefit for small boxes (which decompose into a single block) or collapsed/dense single-droplet systems. For the whole-chain moves (slither and pull) a chain only parallelizes if all its beads fit inside one block's interior (chains spanning a block boundary are frozen that sweep). Since only a short enough chain can ever fit, the chains are split by length alone, and the split is fixed for a given box (it changes only when a RESIZED_EQUILIBRATION box grows to DIMENSIONS): chains no longer than the smallest block interior - and within the kernels' 512-bead per-chain limit (heteropolymer chains for slither, all chains for pull) - are moved by the parallel kernel, and every longer chain is moved by the serial kernel, with both passes run in every megamove, in an order drawn at random (a fair coin) whenever both have chains, which keeps the megamove reversible as a system-wide TSMMC excursion requires. The split depends only on chain lengths, the box and the frozen set, never on the current configuration, which is what keeps PARALLELIZE from changing the equilibrium being sampled. Deciding per megamove from the chains' current extents (the behaviour before 1.0.8) does change it: the parallel kernel can never extend a chain past a block interior while the serial kernel can, so a state-dependent choice between them biases the run towards compact chains. The block decomposition is independent of the thread count, so results are identical for any number of threads. NOTE the parallel sampler targets the same equilibrium distribution but follows a DIFFERENT (and per-step slower-relaxing) Markov chain than the serial run: each sweep only the beads inside the block interiors can move (the frozen halos are re-drawn every sweep), so a run that has not reached equilibrium - e.g. a collapsing system - will show a different (less relaxed) energy at a given step count than the serial run. Compare equilibrium averages, not energies at a fixed step. The parallel kernels log the sub-moves they actually attempt, so a sweep whose random block shift leaves nothing movable adds 0 attempts to MOVE_FREQS.dat; and the chains on the parallel kernel share the slither/pull budget, picked at random within each block, rather than each getting exactly SLITHER_SUBSTEPS or PULL_SUBSTEPS attempts. Frozen chains (via FREEZE_FILE) are fully supported: their beads are excluded from moves but kept in place as fixed, energy-contributing obstacles, so PARALLELIZE applies even with a freeze file. If the compiled kernels were built without OpenMP the blocks are executed one after another instead (no speed-up, identical sampling); the start-up parallelization report says which kernel each move will really use and why. Default = False."],
    'PARALLEL_THREADS': ["int", "Number of OpenMP threads used when PARALLELIZE is True. Must be between 0 and 1024. 0 (the default) means use every CPU available to the process: the value of the OMP_NUM_THREADS environment variable if it is set to a positive integer, otherwise the number of CPUs the process is allowed to run on (its CPU affinity, which is what a batch scheduler such as SLURM restricts, where the platform exposes it), otherwise the machine's core count. A value above 1024 is refused when the keyfile is read: the box is split into at most 64 blocks, so no more than 64 threads ever have work, and an absurd request makes the OpenMP runtime abort the process when it cannot create the threads. Ignored when PARALLELIZE is False; setting it in a keyfile with no PARALLELIZE line prints a warning saying so. Results do not depend on the thread count."],
    'ENERGY_CHECK' : ["int", "Frequency (in steps) with which a full from-scratch energy recompute is compared against the incrementally tracked energy, and the occupancy/type grids are cross-checked against the chain objects. It runs throughout the run, including equilibration, and each passing check prints the full energy decomposition (short-range, long-range, super-long-range, angles) to STDOUT. Either failure (an energy mismatch, or grids that disagree with the chains) stops the run with a SimulationEnergyException after closing traj.xtc (writing out the buffered trajectory under SAVE_AT_END) and dumping the current configuration to CONFIG_AT_ENERGY_FAIL.pdb/.xtc. This is an O(N) safety/debugging check - cheap to run occasionally, expensive every step for large systems. Set to 0 to switch it off entirely (not recommended)."],
    'RESTART_FREQ' : ["int (positive)", "Frequency with which the simulation state is saved to restart.pimms, as a positive integer step frequency. Snapshots are taken on production steps that are multiples of RESTART_FREQ, each replacing restart.pimms (restart snapshots, like all analysis, are suppressed during equilibration); if omitted the default is N_STEPS/10 (rounded down, at least 1), and the resulting number is announced at start-up. A restart of the final state is always written when the run completes. Each write is atomic, so an interrupted run never leaves a torn restart.pimms behind. See the Restart files documentation."],
    'RESTART_FILE' : ["string", "Filepath (relative or absolute, with a leading ~ expanded; a relative path is resolved from the directory PIMMS is run in) to a restart.pimms file to start the simulation from; the simulation fails if it does not exist, and an empty value is rejected. When set, it supplies the initial configuration and the CHAIN keyword is not required (the chains come from the restart file, and any keyfile CHAIN lines are ignored; EXTRA_CHAIN adds new ones). By default this starts a NEW run from that configuration: its steps are numbered from 1 again (step 0 is the restart configuration) and the random numbers come from this run's SEED. To resume the run that wrote the file exactly, add RESTART_CONTINUE. Restart files are Python pickles, so only load files you generated yourself or otherwise trust. See the Restart files documentation for the dimension/hardwall compatibility rules."],
    'RESTART_OVERRIDE_DIMENSIONS' : ["bool", "If True, IGNORE the keyfile DIMENSIONS and adopt the restart file's box exactly as it was saved (a convenience for continuing in the original box without repeating its size in the keyfile). DIMENSIONS is still a required keyword and must still be present (any valid box will do), and the box taken from the restart file is re-validated against the >= 7 per-axis floor. It does NOT grow the box, and is incompatible with RESIZED_EQUILIBRATION. Setting it to True without a RESTART_FILE has no effect and prints a warning. If False (default), the keyfile DIMENSIONS is used (with the same number of dimensions as the restart file) and reconciled with the restart: for a HARDWALL restart the keyfile box must be >= the restart box in every axis, and a larger box is grown with the restart box centred inside it (growing into a bigger box therefore needs NO override); for a periodic (PBC) restart the keyfile DIMENSIONS must match the restart box exactly. The box can never be made smaller than the restart box. Default = False."],
    'RESTART_CONTINUE' : ["bool", "Boolean flag which, if set to True, resumes the run that wrote RESTART_FILE as though it had never stopped: the step counter continues from the step the restart file was written at, the random-number generators are restored to their state at that step, and the temperature is restored to its value at that step (the keyfile must describe the thermal protocol of the run being resumed: for a quench, QUENCH_RUN with the QUENCH_START, QUENCH_END, QUENCH_STEPSIZE and QUENCH_FREQ recorded in the restart file; otherwise a TEMPERATURE equal to the restart file's temperature. The angle-penalty scaling is built from the keyfile, so anything else would be a different energy function), so the resumed run makes exactly the moves the uninterrupted run would have made and reproduces its per-step output bit for bit (the same ENERGY.dat rows, trajectory frames and per-step analysis rows, labelled with the same step numbers). Two kinds of output are not continued, because the restart file does not store them: the move counters in MOVE_FREQS.dat, ACCEPTANCE.dat and TOTAL_MOVES.dat start again from zero at the checkpoint, and the end-of-run averages (INTSCAL.dat, INTSCAL_SQUARED.dat, SCALING_INFORMATION.dat and DISTANCE_MAP.dat) cover the resumed segment only. N_STEPS is the TOTAL length of the run, so it must be larger than the restart file's step, and the resumed segment writes rows for the steps after that one; frame 0 of its trajectory is the restart configuration. Requires RESTART_FILE, and a restart file written by PIMMS 1.0.8 or later (earlier files carry no generator state and are refused). SEED must not be given, since the generator state comes from the file; the box, boundary condition and chain list must be exactly those of the original run, so RESIZED_EQUILIBRATION and EXTRA_CHAIN are refused and DIMENSIONS / HARDWALL must equal the restart file's (RESTART_OVERRIDE_DIMENSIONS and RESTART_OVERRIDE_HARDWALL are the easy way to guarantee that). Each of these conditions is checked at start-up and every violation is listed in a single refusal, except that a periodic restart file with HARDWALL : True or with RESIZED_EQUILIBRATION is refused first by the general restart rules. As a last guard the energy of the restart configuration is recomputed under this keyfile and must equal the energy stored in the restart file, so a different PARAMETER_FILE, ANGLES_OFF or NON_INTERACTING is refused as well. A continuation must be run in a fresh directory (holding the keyfile, the input files and a copy of the restart file): start-up deletes the output files an earlier run left behind and the resumed segment only writes the steps after the checkpoint, so a continuation in a directory that still holds simulation output is refused before anything is written or deleted. Default = False."],
    'RESTART_OVERRIDE_HARDWALL' : ["bool", "Boolean flag which, if set to True, means that the hardwall setting of the simulation is overridden by the hardwall setting in the restart file, rather than using the keyfile HARDWALL value. Without it the keyfile HARDWALL is used: a hardwall restart file can start either a hardwall or a periodic run, but a periodic restart file cannot start a hardwall run (refused, since its chains may cross a periodic face). It does NOT relax the other restart rules: a periodic restart file still cannot be combined with RESIZED_EQUILIBRATION, and the cluster-rotation box rule is applied to whichever HARDWALL/DIMENSIONS combination the run ends up with. Setting it to True without a RESTART_FILE has no effect and prints a warning. Default = False."],
    'EXTRA_CHAIN' : ['See description', "One of the few multi-component keywords in PIMMS, and it can ONLY be used when a RESTART_FILE is defined (using it without one is an error). This keyword allows you to add additional chains into the system that were not originally present in the RESTART_FILE. The format follows the same as the CHAIN keyword (so <number of chains>  <chain sequence>) and multiple EXTRA_CHAIN lines can be included for different types of chains. This means you can setup an initial set of simulations, and then run a simulation from the end-state of the original simulation with new chains added. Moreover, this can be repeated an arbitrary number of times. New chains are randomly inserted so they do not overlap with existing chains (the same overcrowding refusals as CHAIN apply), and are given chainIDs that follow on from the restart file's. Sequences are upper-cased unless CASE_INSENSITIVE_CHAINS is False. An extra chain whose sequence already exists in the restart file joins that existing chain type rather than defining a new one. Refused with RESTART_CONTINUE, since adding chains makes it a different system."],
    'QUENCH_RUN' : ["bool", "Boolean flag which, if set to True, means that the simulation is a quench run. This means that the simulation starts at one temperature and then systematically changes to a different temperature. Generally this will be higher to cooler, but could be cooler to higher. Note that the starting temperature is set by QUENCH_START and ending temperature by QUENCH_END, so the TEMPERATURE keyword is overwritten with QUENCH_START (with a warning if the two disagree). Also, all the QUENCH keywords (QUENCH_START, QUENCH_END, QUENCH_FREQ, QUENCH_STEPSIZE and QUENCH_AS_EQUILIBRATION) must all be set; conversely they have no effect unless QUENCH_RUN is True, and giving any of them in a keyfile with no QUENCH_RUN line at all prints a warning that names them and says the run is at constant TEMPERATURE (a keyfile that says QUENCH_RUN : False has switched the quench off in so many words and gets no warning; either way their values are still checked: a QUENCH_START or QUENCH_END of 0 or less, a QUENCH_STEPSIZE of 0 or a QUENCH_FREQ of 0 or less is refused even then). The temperature trajectory is recorded in QUENCH.dat, one tab-separated line (step, new temperature, energy at the end of that step) per temperature change, and any ANGLE_PENALTY_T_NORM values in the parameter file are scaled once by QUENCH_END (the production temperature), so T-normalised angle penalties are exact in kT only at the end of the ramp. Default = False."],
    'QUENCH_AS_EQUILIBRATION' : ["bool", "Boolean flag, required (True or False) whenever QUENCH_RUN is True. If set to True, the equilibration period is used for the quench, and after the equilibration period the simulation temperature is fixed at the QUENCH_END temperature. Setting this OVERWRITES whatever EQUILIBRATION you gave with the length of the quench, (1 + number of temperature changes) x QUENCH_FREQ steps, where the number of temperature changes is the QUENCH_START to QUENCH_END range divided by QUENCH_STEPSIZE, rounded up; the new EQUILIBRATION is announced at start-up. The extra QUENCH_FREQ steps let the system equilibrate at the final temperature before production begins. If set to False, EQUILIBRATION is left as you gave it and the ramp simply runs from the start of the run, so it continues into production if EQUILIBRATION is shorter than the quench."],
    'QUENCH_START' : ["float (positive)", "Starting temperature for the quench run; must be greater than 0. Also becomes the simulation TEMPERATURE."],
    'QUENCH_END' : ["float (positive)", "Ending (production) temperature for the quench run; must be greater than 0. May be above or below QUENCH_START (a cooling or a heating ramp), and is the temperature that any ANGLE_PENALTY_T_NORM values are scaled by."],
    'QUENCH_FREQ' : ["int (positive)", "Number of steps spent at each temperature before the next change; must be greater than 0. The temperature changes (before that step's move) on every step that is a multiple of QUENCH_FREQ, until QUENCH_END is reached. The whole quench occupies (1 + number of temperature changes) x QUENCH_FREQ steps, which must be smaller than N_STEPS or the run is rejected at start-up."],
    'QUENCH_STEPSIZE' : ["float (positive)", "The amount by which the temperature is changed at each QUENCH_FREQ steps. Give it as a magnitude: the sign is ignored (the absolute value is taken and PIMMS works out the direction from QUENCH_START and QUENCH_END). It must be greater than 0 and no larger than the QUENCH_START to QUENCH_END range (so QUENCH_START and QUENCH_END must differ). The step need not divide the range: the final step is clamped so the ramp lands exactly on QUENCH_END."],
    'MOVE_CRANKSHAFT' : ["float", "Move code 1. Probability of a crankshaft megamove being attempted. When selected, CRANKSHAFT_SUBSTEPS single-bead sub-moves are performed across the whole system, each with its own Metropolis accept/reject, so a single crankshaft step is a large amount of MC work. This is the workhorse move for relaxing local chain conformation and should usually carry most of the probability. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'CRANKSHAFT_SUBSTEPS' : ["int (positive)", "Number of single-bead sub-moves performed per crankshaft megamove, in TOTAL across the whole system (each sub-move picks a bead uniformly at random from all the non-frozen beads in the system, so the work is not per chain; with PARALLELIZE the sub-moves are shared out over the beads inside the block interiors). Must be greater than 0; the default of 500 is deliberately small, and for production runs 20-50K (or more) is a more usual choice. It may not exceed 2147483647 (2^31 - 1): the kernel indexes the sub-moves of a megamove with a 32-bit integer and PIMMS builds an array of one 8-byte entry per sub-move for every megamove, so a larger value is refused when the keyfile is read, as is one whose array would not fit in physical memory, and a value above 10^8 (an 0.8 GB array per megamove) gets a warning. To do more work, raise N_STEPS instead. The same value sets the length of each of the two single-chain relaxations inside a MOVE_JUMP_AND_RELAX move. Setting it away from the default in a keyfile that gives neither MOVE_CRANKSHAFT nor MOVE_JUMP_AND_RELAX prints a warning that it has no effect (unless it is used by the crankshaft fallback of a MOVE_SYSTEM_TSMMC excursion)."],
    'SLITHER_SUBSTEPS' : ["int (positive)", "Number of slither (reptation) moves applied to EACH chain, in random order, per slither megamove. Must be greater than 0. A slither advances a chain forwards or backwards like a snake. For homopolymers (every bead the same type and long-range flag) the interaction-energy change is a single O(1) end-bead evaluation plus a whole-chain angle term; for heteropolymers every residue is re-evaluated; single-bead chains become a local translation. With PARALLELIZE the chains handled by the parallel kernel share a budget of SLITHER_SUBSTEPS per chain and are picked at random within each block, so the exact per-chain count holds only on the serial kernel. The megamove builds an array of SLITHER_SUBSTEPS x (number of chains) 8-byte entries, indexed with a 32-bit integer, so both SLITHER_SUBSTEPS and that product may not exceed 2147483647 (2^31 - 1) or physical memory (refused when the keyfile is read), and a product above 10^8 (0.8 GB per megamove) gets a warning. Setting it away from the default in a keyfile with no MOVE_SLITHER line prints a warning that it has no effect."],
    'PULL_SUBSTEPS' : ["int (positive)", "Number of pull moves applied to EACH eligible chain, in random order, per pull megamove. Must be greater than 0. A pull move displaces an interior bead and cooperatively 'pulls' the rest of the segment along to restore connectivity, letting chains rearrange in dense systems where rigid moves would clash. Chains shorter than 3 beads have no interior bead and are skipped. With PARALLELIZE the chains handled by the parallel kernel share a budget of PULL_SUBSTEPS per chain and are picked at random within each block, so the exact per-chain count holds only on the serial kernel. The megamove builds an array of PULL_SUBSTEPS x (number of chains of 3 or more beads) 8-byte entries, indexed with a 32-bit integer, so both PULL_SUBSTEPS and that product may not exceed 2147483647 (2^31 - 1) or physical memory (refused when the keyfile is read), and a product above 10^8 (0.8 GB per megamove) gets a warning. Setting it away from the default in a keyfile with no MOVE_PULL line prints a warning that it has no effect."],

    'MOVE_VMMC' : ["float", "Move code 14. Probability of a Virtual-Move Monte Carlo (VMMC) collective move being attempted (Whitelam and Geissler, J. Chem. Phys. 127, 154101, 2007). A seed chain is given a trial rigid translation; neighbouring chains are recruited into a moving cluster according to interaction-energy gradients (a neighbour is recruited when moving the seed alone would break their mutual attraction), and the whole cluster translates together. This avoids the kinetic traps that single-chain moves hit in strongly-attractive / condensed phases, while maintaining detailed balance. The seed is never a frozen chain, and recruiting one rejects the move; each attempt that survives recruitment and the overlap check costs a full from-scratch energy evaluation (O(N) in the number of beads), so keep the fraction modest for large systems. It is a rigid-body move that never changes a chain's shape, so a move set made only of rigid moves (chain translate/rotate, cluster translate/rotate, VMMC) gets a start-up warning (when any unfrozen chain has three or more beads) that every chain keeps the conformation it was built with. EXPERIMENTAL - requires EXPERIMENTAL_FEATURES : True. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'VMMC_MAX_DISPLACEMENT' : ["int (positive)", "Maximum magnitude (per dimension, in lattice units) of the rigid translation proposed by a VMMC move. Must be greater than 0. Every axis is displaced by between 1 and this many sites (never more than the box length minus 1), with a random sign. Small values give local collective moves (recommended); large values rarely succeed in dense phases. Default 3. Changing this from the default requires EXPERIMENTAL_FEATURES : True (VMMC is experimental)."],

    'VMMC_MAX_CLUSTER' : ["int (positive)", "Upper bound on the VMMC cluster size used for the 1/n_c move-frequency correction; the recruited cluster is aborted (move rejected) if it would exceed the drawn cutoff. Must be greater than 0, and is clamped to the number of chains at runtime, so it has no upper limit of its own. Default 1000. Changing this from the default requires EXPERIMENTAL_FEATURES : True (VMMC is experimental)."],
    'MOVE_CHAIN_TRANSLATE' : ["float", "Move code 2. Probability of a whole-chain rigid translation move being attempted. The offset is drawn uniformly over the whole box on each axis and wrapped, so a translation relocates the chain anywhere in the box rather than taking a local step (under HARDWALL a translated chain that would straddle a wall is rejected; see HARDWALL). Works on chains of any length, including single beads. It is a rigid-body move that never changes a chain's shape, so a move set made only of rigid moves (chain translate/rotate, cluster translate/rotate, VMMC) gets a start-up warning (when any unfrozen chain has three or more beads) that every chain keeps the conformation it was built with. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CHAIN_ROTATE' : ["float", "Move code 3. Probability of a whole-chain rigid rotation move being attempted (a cardinal 90/180/270 degree rotation about the bead nearest the chain centroid). Needs chains of at least 2 beads: drawn for a single-bead chain it is a null move that is simply rejected, and PIMMS warns at start-up if monomers are present with this move enabled. It is a rigid-body move that never changes a chain's shape, so a move set made only of rigid moves (chain translate/rotate, cluster translate/rotate, VMMC) gets a start-up warning (when any unfrozen chain has three or more beads) that every chain keeps the conformation it was built with. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CHAIN_PIVOT' : ["float", "Move code 4. Probability of a molecular pivot move being attempted. Pivot moves randomly select an interior bead (any but the two termini) and rotate the SHORTER arm about it, leaving the rest of the chain fixed. Needs chains of at least 3 beads; shorter chains give a null move that is rejected, and PIMMS warns at start-up if single-bead chains are present with this move enabled. Because the shorter arm always moves (the N-terminal arm on a tie), the two beads at a chain's midpoint never move relative to the rest of the chain, so a move set whose only shape-changing moves are pivots gets a start-up warning. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_HEAD_PIVOT' : ["float", "Move code 5. Probability of a head pivot move being attempted. Head pivot moves randomly select one of the two ends of a chain and pivot that terminus, but this is almost never worth doing so we recommend setting it to 0. Needs chains of at least 2 beads; PIMMS warns at start-up if single-bead chains are present with this move enabled, and if it is the only shape-changing move (only the termini would ever move). Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CLUSTER_TRANSLATE' : ["float", "Move code 7. Probability of a cluster translation move being attempted: the cluster containing a randomly chosen chain (every chain linked to it through short-range, Chebyshev-1, contacts) is translated rigidly as a body by a random nonzero offset on each axis. The move is rejected if the moved cluster would touch any other chain (the cluster must remain the same cluster), if it contains a frozen chain, or if it contains every chain in the system, so in a fully condensed single-droplet system it does nothing (a system of one chain is the exception: its lone chain is translated like any other cluster, subject to the same rejections). Under HARDWALL a proposal that would carry any bead out of the box is rejected rather than wrapped. Cluster translation moves are relatively expensive, so in general wise to keep this at a low number (0.01 to 0.05). It is a rigid-body move that never changes a chain's shape, so a move set made only of rigid moves (chain translate/rotate, cluster translate/rotate, VMMC) gets a start-up warning (when any unfrozen chain has three or more beads) that every chain keeps the conformation it was built with. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],
    'MOVE_CLUSTER_ROTATE' : ["float", "Move code 8. Probability of a cluster rotation move being attempted: the cluster containing a randomly chosen chain (as for MOVE_CLUSTER_TRANSLATE) is rotated rigidly by a cardinal 90/180/270 degrees about the bead nearest its centroid. The same rejections apply (contact with another chain, a frozen member, a cluster of every chain in a system of more than one chain), and under periodic boundaries a cluster that winds around the box is also rejected. Under periodic boundaries the production box must be cubic/square (see DIMENSIONS). A draw that maps the cluster exactly onto itself is rejected as a null move; every draw on an isolated single-bead cluster is one, and PIMMS warns at start-up when single-bead chains are present with this move enabled. Cluster rotation moves are relatively expensive, so in general wise to keep this at a low number (0.01 to 0.05). It is a rigid-body move that never changes a chain's shape, so a move set made only of rigid moves (chain translate/rotate, cluster translate/rotate, VMMC) gets a start-up warning (when any unfrozen chain has three or more beads) that every chain keeps the conformation it was built with. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'CRANKSHAFT_MODE' : ["str (UNIFORM)", "Obsolete and ignored. It used to define how the number of crankshaft sub-moves scaled with chain length. The keyword is still accepted so that old keyfiles keep working, but PIMMS prints a warning and ignores its value: the crankshaft always performs a fixed CRANKSHAFT_SUBSTEPS sub-moves per megamove, independent of chain length (the old UNIFORM mode). Remove it from new keyfiles."],

    'MOVE_SLITHER' : ["float", "Move code 6. Probability of a slither (reptation) megamove being attempted. When selected, every non-frozen chain is slithered SLITHER_SUBSTEPS times (with PARALLELIZE, the chains on the parallel kernel share that budget instead; see SLITHER_SUBSTEPS) - a chain advances forwards or backwards through the lattice like a snake, which efficiently relaxes chain conformations. Works in 2D and 3D. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_PULL' : ["float", "Move code 11. Probability of a pull (cooperative reptation) megamove being attempted. When selected, every non-frozen chain of length >= 3 is pulled PULL_SUBSTEPS times (with PARALLELIZE, the chains on the parallel kernel share that budget instead; see PULL_SUBSTEPS) - an interior bead is displaced and the following beads are cooperatively 'pulled' along to restore connectivity, letting chains rearrange in DENSE systems where rigid moves would clash (the chain termini are not moved by this move, so pair it with crankshaft/slither). Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_CTSMMC' : ["float", "Move code 9. Probability of a single-chain TSMMC (Temperature-Switch Monte Carlo) move being attempted. A randomly selected chain is taken on a temperature EXCURSION - heated along a schedule from the current simulation temperature up to the jump temperature (TSMMC_JUMP_TEMP, or the current temperature plus TSMMC_FIXED_OFFSET) and cooled back, with crankshaft sub-moves on that chain alone at every rung - to help it escape local energy minima, with a tempered-transitions acceptance that preserves detailed balance. Controlled by the TSMMC_* keywords. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_MULTICHAIN_TSMMC' : ["float", "Move code 10. Probability of a multi-chain TSMMC move being attempted. As MOVE_CTSMMC, but a randomly selected SUBSET of the N non-frozen chains (a number drawn uniformly between 1 and floor(N/4) + 1, never more than N) undergoes the temperature excursion together. Controlled by the TSMMC_* keywords. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_SYSTEM_TSMMC' : ["float", "Move code 12. Probability of a system-wide TSMMC move being attempted. The ENTIRE system undergoes a temperature excursion (heated along a schedule from the current simulation temperature up to the jump temperature and cooled back) to help the whole configuration escape local minima. The sub-moves inside the excursion are drawn from your keyfile move mix, except that a nested TSMMC draw is not allowed and is redrawn from the non-TSMMC moves with their fractions renormalised; if no non-TSMMC move is enabled at all they fall back to crankshaft megamoves, which is announced once. The whole excursion counts as one step and one MOVE_FREQS.dat attempt: its sub-moves do not advance the step counter and nothing is written until it completes. Because it can only do what its sub-moves can, the start-up move-set checks treat it as those moves: a run whose sub-moves cannot act on any chain is refused, and an excursion over rigid moves alone still gets the fixed-conformation warning. Controlled by the TSMMC_* keywords. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'MOVE_JUMP_AND_RELAX' : ["float", "Move code 13. Probability of a single-chain jump-and-relax move being attempted. A selected chain is relaxed (a crankshaft sub-trajectory of CRANKSHAFT_SUBSTEPS sub-moves), a rigid translation ('jump', drawn like a MOVE_CHAIN_TRANSLATE offset, so it can land anywhere in the box) is proposed and accepted or rejected on its own Metropolis criterion, then the chain is relaxed again. Each of the three sub-steps preserves the Boltzmann distribution, so the composite move maintains detailed balance. The jump is evaluated with a full from-scratch energy recompute (O(N) in the number of beads), so keep the fraction modest for large systems. Useful for relocating individual chains and letting them settle; for relocation through dense/condensed phases prefer MOVE_VMMC or MOVE_PULL. Note all provided MOVE_* keywords must be >= 0 and add up to 1.0"],

    'TSMMC_JUMP_TEMP' : ["float", "The peak ('jump') temperature reached during a TSMMC temperature excursion. When a TSMMC move is enabled and TSMMC_FIXED_OFFSET is not set, it MUST be greater than the simulation TEMPERATURE (and, in a QUENCH_RUN, greater than both QUENCH_START and QUENCH_END - a ramp that would reach it is rejected at start-up): a TSMMC move heats the selected chain(s)/system from the current simulation temperature up to TSMMC_JUMP_TEMP and back. Ignored if TSMMC_FIXED_OFFSET is set. It must be greater than 0 whether or not a TSMMC move is enabled. Used only by the TSMMC moves (MOVE_CTSMMC / MOVE_MULTICHAIN_TSMMC / MOVE_SYSTEM_TSMMC): in a keyfile that gives none of those three, setting any TSMMC keyword away from its default prints one warning naming the keywords that have no effect (a keyfile that writes a TSMMC move as 0 has switched it off in so many words and gets no warning). Default 50."],

    'TSMMC_STEP_MULTIPLIER' : ["int (positive)", "Sets how much sampling happens at EACH temperature point of a TSMMC excursion; must be greater than 0. For MOVE_CTSMMC the number of sub-moves per temperature is TSMMC_STEP_MULTIPLIER x the chain length, for MOVE_MULTICHAIN_TSMMC it is TSMMC_STEP_MULTIPLIER x the number of beads in the selected chains, and for MOVE_SYSTEM_TSMMC it is exactly TSMMC_STEP_MULTIPLIER sub-moves, each an ordinary move drawn from the keyfile move mix (so a crankshaft draw is a whole CRANKSHAFT_SUBSTEPS megamove). An excursion therefore costs 2 x TSMMC_NUMBER_OF_POINTS + 10 times that many sub-moves. Larger values equilibrate the system more thoroughly at each temperature, at the cost of much slower excursions. The single-chain and multi-chain moves build an array of one 8-byte entry per sub-move at every temperature, indexed with a 32-bit integer, so TSMMC_STEP_MULTIPLIER, and for an enabled MOVE_CTSMMC or MOVE_MULTICHAIN_TSMMC its product with the longest chain or with the largest set of chains the move can select, may not exceed 2147483647 (2^31 - 1) or physical memory (refused when the keyfile is read); a product above 10^8 gets a warning. Used only by the TSMMC moves: with no TSMMC move enabled it is still range-checked, and setting it away from the default in a keyfile that gives none of MOVE_CTSMMC, MOVE_MULTICHAIN_TSMMC and MOVE_SYSTEM_TSMMC prints a warning that it has no effect. Default 50."],

    'TSMMC_NUMBER_OF_POINTS' : ["int (positive)", "Number of temperature points on the heating ramp of a TSMMC excursion; must be greater than 0. The ramp is TSMMC_NUMBER_OF_POINTS equally spaced temperatures, rising from one step above the current simulation temperature to exactly the jump temperature (the simulation temperature itself is not a rung). The full schedule is that ramp, ten further rungs at the jump temperature, then the ramp in reverse: 2 x TSMMC_NUMBER_OF_POINTS + 10 rungs in all, each sampled for the number of sub-moves set by TSMMC_STEP_MULTIPLIER. More points give a smoother (more gradual, and more expensive) heating/cooling ramp. It may not exceed 1000000: the schedule temperatures are rounded to five decimal places, so a finer ramp only repeats rungs, and every rung is at least one kernel call. Used only by the TSMMC moves: with no TSMMC move enabled it is still range-checked, and setting it away from the default in a keyfile that gives none of MOVE_CTSMMC, MOVE_MULTICHAIN_TSMMC and MOVE_SYSTEM_TSMMC prints a warning that it has no effect. Default 20."],

    'TSMMC_INTERPOLATION_MODE' : ["str (LINEAR)", "How the temperature is interpolated between the current simulation temperature and the jump temperature across the excursion schedule. Currently the only supported value is LINEAR (equal temperature increments); anything else is rejected at parse time. Used only by the TSMMC moves. Default LINEAR."],

    'TSMMC_FIXED_OFFSET' : ["float (positive)", "If set, the TSMMC jump temperature is defined RELATIVE to the current simulation temperature (the temperature in force, which in a quench is the ramp's current rung) plus TSMMC_FIXED_OFFSET, rather than using the absolute TSMMC_JUMP_TEMP. Must be greater than 0 (checked whether or not a TSMMC move is enabled). This is what you want in a QUENCH_RUN, where the jump temperature then tracks the ramp instead of being overtaken by it. The default is False (use the absolute TSMMC_JUMP_TEMP), but False is NOT keyfile syntax - simply omit the keyword. Used only by the TSMMC moves."],

    'ANALYSIS_MODULE' : ["str (path)", "Filepath (relative or absolute, with a leading ~ expanded; a relative path is resolved from the directory PIMMS is run in) to a user-supplied Python module defining a top-level analysis_function(step, lattice), which PIMMS calls every ANA_CUSTOM steps during production with the step number and the live Lattice object (read it, do not change it; any return value is ignored). The module is imported and validated at parse time (it must import cleanly and define a callable analysis_function accepting those two arguments), so a broken module fails immediately rather than part-way through a run, and any exception it raises at runtime stops the run with an error naming YOUR code. Giving ANALYSIS_MODULE without a positive ANA_CUSTOM is an error (the module would never run). Omit the keyword (the default) to run no custom analysis; an empty value is rejected."],

    'ANA_CUSTOM' : ["int", "Frequency (in steps) at which the user-defined custom analysis function (from ANALYSIS_MODULE) is run, on production steps only. 0 (default) disables it; unlike the other analysis frequencies it does not default to ANALYSIS_FREQ, so it must be set whenever ANALYSIS_MODULE is given. Setting it without an ANALYSIS_MODULE prints a warning and does nothing."],

    'ANA_CLUSTER_THRESHOLD' : ["int", "Connected components containing MORE than this many chains get the per-cluster shape/size analysis (CLUSTER_RG/ASPH/AREA/VOL/DEN and radial profiles). Must be 0 or greater; the comparison is strict, so the default of 1 skips single chains and 0 includes them. A threshold you set at or above the number of chains in the system can never be exceeded, so those files would never be written; PIMMS prints a warning. The radial profiles additionally need at least 27 beads. Note the size-distribution files (CLUSTERS.dat / NUM_CLUSTERS.dat and the LR_* variants) always include EVERY component regardless of this threshold."]}


# Logical groupings of keywords used to organise the `PIMMS --info` output and
# the generated docs/keywords.rst under subheadings (ordered). Every keyword in
# EXPECTED_KEYWORDS should appear in exactly one group; any that do not are shown
# under "Other" by both.
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


# indices correspond to
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

# Defaults the keyfile parser works out at start-up rather than taking verbatim from
# DEFAULTS: every unset per-analysis frequency inherits the keyfile's ANALYSIS_FREQ
# (KeyFileParser.assign_default), RESTART_FREQ becomes max(1, N_STEPS // 10)
# (KeyFileParser.set_dynamic_defaults) and an unset SEED is drawn at random. Showing
# the static DEFAULTS value for these would be wrong as soon as ANALYSIS_FREQ or
# N_STEPS is set, so the keyword reference and `PIMMS --info` show these instead.
DERIVED_DEFAULTS = {
    'ANA_POL': 'ANALYSIS_FREQ',
    'ANA_INTSCAL': 'ANALYSIS_FREQ',
    'ANA_DISTMAP': 'ANALYSIS_FREQ',
    'ANA_ACCEPTANCE': 'ANALYSIS_FREQ',
    'ANA_INTER_RESIDUE': 'ANALYSIS_FREQ',
    'ANA_CLUSTER': 'ANALYSIS_FREQ',
    'RESTART_FREQ': 'N_STEPS / 10',
    'SEED': 'random',
}


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


def display_default(keyword):
    """
    Return the default of a keyfile keyword as a user would write or read it.

    Most defaults are shown as stored in ``DEFAULTS``. The exceptions are the
    defaults the parser derives at start-up (``DERIVED_DEFAULTS``) and the
    placeholders ``DEFAULTS`` uses internally: ``False`` for a keyword whose value
    is not a boolean (a file path, a box, a temperature offset) and ``'UNSET'``
    both mean "not set" and are shown as ``unset``, and an empty list or ``'N/A'``
    is shown as ``none``. Writing ``False`` for a path keyword in a keyfile would
    not mean "unset", which is why the raw value must not be shown.

    Parameters
    ----------
    keyword : str
        An upper-case keyfile keyword.

    Returns
    -------
    str or None
        The display string, or None if the keyword has no default (it is in
        ``REQUIRED_KEYWORDS``, or has no ``DEFAULTS`` entry).
    """
    # a required keyword has no default even when DEFAULTS holds a placeholder
    # for it (TEMPERATURE's 'N/A')
    if keyword not in DEFAULTS or keyword in REQUIRED_KEYWORDS:
        return None
    if keyword in DERIVED_DEFAULTS:
        return DERIVED_DEFAULTS[keyword]
    value = DEFAULTS[keyword]
    type_string = KEYWORDS_DESCRIPTION.get(keyword, [''])[0]
    if value is False and not type_string.startswith('bool'):
        return "unset"
    if isinstance(value, list) and len(value) == 0:
        return "none"
    if value == 'UNSET':
        return "unset"
    if value == 'N/A':
        return "none"
    return str(value)
