pimms
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/pimms.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/pimms)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pimms/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pimms/branch/master)

Lattice simulation package for biomolecule.

### Preamble
If you are reading this I sent you PIMMS because I know and trust you. Congrats! The current version (1.24) is a pre-release beta candidate, that appears to functioning correctly but is being actively worked on.

This is all to say - use this version to test-drive ideas, but I would not recommend publishing work you do with PIMMS until we finalize a distributable Python 3 version. There _shouldn't_ be any bugs in this version... but there could be. So we're working on it!

### Installation

From this directory first make sure all the requiste code is installed. This requires `conda` because some of the packages are (in my experience) only easily installed via conda.

Assuming `conda` is installed and your in the relevant environment:

	# first install the baseline packages
	conda install numpy scipy cython pandas 
	
	# then install mdtraj, which provides the xtc library backend
	conda install mdtraj
	
	# finally install pimms from this directory by running. These flags are not strictly necessary for the first install, but they ensure the cython always gets recompiled on a new installation
	pip install -e . --upgrade --force-reinstall
	
	
### Usage

Once installed, the PIMMS binary should 


### Keyfile keywords
Keyword | Format (type) | Description
:---: | :---: | :---: 
DIMENSIONS | **INT** <br>(2 or 3) <br><br>`A x B`<br> or<br> `A x B x C`  | Size of the simulation box (in lattice units). 2D or 3D (defines if the simulation is a 2D or 3D simulation)
CHAIN | *See description* | One of the few multi component keywords in PIMMS and the only keyword that can appear multiple times, the `CHAIN` keyword defines a specific polymer chain and the number of that chain that will exist in the simulation. The format should be <br> <br> `CHAIN : N  {CHAIN IDENTIY}` <br><br> Where `N` defines the number of the chain and `{CHAIN IDENTITY}` gives polymer sequence in one-letter alphabet code. As an example <br> <br> `CHAIN : 20 QQQQQQQQQQ`<br><br> Would give 20 poly-glutamine polymers. In later versions of PIMMS we will be updating this to allow the reading of keyfiles that use three-letter codes 
TEMPERATURE | **INT** | Simulation temperature to be used (units are arbitrary and depend on the unist of the parameter file)
N_STEPS | **INT**|Number of MAIN CHAIN Monte Carlo steps
PARAMETER_FILE | **STRING** | File location for the parameter file to be read in
EQUILIBRATION | **INT** | Number of steps to be run as equilibration (i.e. before any analysis or trajectory output is generated)
PRINT_FREQ | **INT** | Step frequency at which the system status is printed to STDOUT
XTC_FREQ | **INT** | Step frequency at which the system system configuration is written to XTC file
EN_FREQ | **INT** | Step frequency at which the system energy is written to ENERGY.dat
SEED | **INT** | Random seed to allow reproducible runs
ENERGY_CHECK | **INT** | Step frequency at which global energy is recalculated and compared to the current current as determined locally on each step. All energy calculations are exact, so this is primarily for sanity checking. For large systems this can be an expensive operation.
NON_INTERACTING | **BOOL** | Boolean that defines if the Hamiltonian is non-interacting or not. This is a convenient way to generate 

MOVE_ keywords define the frequency with which different moves are performed during the simulation. These words must sum up to 1.0, and all must be defined.

Keyword | Format (type) | Description
:---: | :---: | :---: 
MOVE\_CHAIN_TRANSLATE |  **FLOAT** | Single chain rigid body translation
MOVE\_CHAIN_ROTATE |  **FLOAT** | Single chain rigid body rotation
MOVE\_CHAIN_PIVOT |  **FLOAT** | Chain pivot at a random potion
MOVE\_HEAD_PIVOT |  **FLOAT** | Pivot the head residue of the chain (this is a somewhat redundant move - recommend setting to 0)
MOVE\_SLITHER |  **FLOAT** |  Slither the chain through the system (BROKEN, do not use)
MOVE\_CLUSTER_TRANSLATE | **FLOAT** | Translate a randomly selected contigous cluster of chains
MOVE_CLUSTER_ROTATE |  **FLOAT** |  Rotate a randomly selected contigous cluster of chains
MOVE\_CTSMMC |  **FLOAT** |  Single chain TSMMC
MOVE\_MULTICHAIN_TSMMC  | **FLOAT** |  Multiple randomly selected chains undergo TSMMC
MOVE\_SYSTEM_TSMMC | **FLOAT** | Entire system undergoes TSMMC
MOVE\_RATCHET_PIVOT |  **FLOAT** | A single chain undergoes a directed pivot move. Haven't yet convinced myself this doesn't break detailed balance


The following block of keywords define various options that controll temperature-sweep Metropolis Monte Carlo moves

Keyword | Format (type) | Description
:---: | :---: | :---: 
MOVE_CRANKSHAFT | **FLOAT** | 
CRANKSHAFT_SUBSTEPS | **INT** |
CRANKSHAFT_MODE | **STRING** |
TSMMC\_JUMP_TEMP | **FLOAT** | 
TSMMC\_STEP_MULTIPLIER | **INT** | 
TSMMC\_INTERPOLATION_MODE | **STRING** | 
TSMMC\_NUMBER_OF_POINTS | **INT** |
TSMMC\_FIXED_OFFSET | **INT** |

The following block of keywords define various options that controll quench based simulations, which allow the temperature to systematically drop throughout the simulation.

Keyword | Format (type) | Description
:---: | :---: | :---: 
QUENCH_RUN | **BOOL** | 
QUENCH_FREQ |  **INT** |
QUENCH_STEPSIZE | **FLOAT** | 
QUENCH_START | **FLOAT** | 
QUENCH_END | **FLOAT** | 
QUENCH\_AS_EQUILIBRATION | **BOOL** | 

The following block of keywords define various options that controll on-the-fly analysis done in CAMPARI

Keyword | Format (type) | Description
:---: | :---: | :---: 
ANALYSIS_FREQ  | **INT** | Step frequeny at which all default analyis is performed. If no 
ANA_POL | **INT** |
ANA_INTSCAL | **INT** |
ANA_DISTMAP | **INT** |
ANA_ACCEPTANCE | **INT** |
ANA\_INTER_RESIDUE | **INT** |
ANA\_CLUSTER | **INT** |
ANA\_RESIDUE_PAIRS | **INT** |


### Copyright

Copyright (c) 2015-2020, Alex Holehouse


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
