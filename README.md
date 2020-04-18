# PIMMS: Polymer Interactions in Multi-component MixtureS
---

[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/pimms.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/pimms)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pimms/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pimms/branch/master)


--- 

###### PIMMS version 0.1.26.4 (April 18th 2020)

## Preamble
If you are reading this, I've sent you PIMMS because I trust you will use it responsibly at this stage. The current version (1.25) is a pre-release beta candidate. This means that while appears to be functioning correctly it is being actively worked on. With that in mind, you are being enlisted to help us find bugs! Please (for now) email any bugs to alex. "Bugs" here include things that 

1. Are wrong/don't make sense
2. Don't behave in a way you expect
3. Errors/exceptions when you wouldn't expect them

More generally - while I'd love for people to use this version to test-drive ideas, I would not recommend publishing work you do with this specific version of PIMMS until we finalize a stable release. There _shouldn't_ be any bugs in this version... but there could be. So we're working on it! This includes the development and deployment of a large unit-test suite, which takes time, but we're working hard! 

## Background
#### What is PIMMS?
PIMMS is a lattice-based simulation engine that allows both 2D and 3D simulations to be performed. Useful features include:

1. Easy to use! Upon installation a command-line executable (`PIMMS`) is available, should be in you `$PATH` variable, and can be used to run simulations. No messing around, it (should) just work! 
2. Easy to define interaction parameters through a simple parameter file (example included in `/demo_keyfiles/demo_1/params.prm`)
3. Easily run fast 2D or 3D lattice based simulations
4. Run simulations with many distinct components
5. Run simulations of a single homo or heteropolymer 
5. Run simulations of many copies of polymers to explore phase behaviour
7. Drive interactions over three distinct length-scales
8. Various other things

#### How is PIMMS written?
PIMMS is written almost fully in Python (=>3.7) with the most computationally intensive parts written in fully optimized `Cython` that compiles down to native C. We're still ironing out kinks, but by having most of the complex behavior in Python maintenance and development is fast and efficient. However, certain functionality (i.e. analysis of large systems) is as a result disproportionately expensive vs. the actual simulation, so you may wish to alter the frequency at which certain analysis routines are performed based on your interests.

PIMMS is a relatively large codebase of ~20K lines of (mostly) Python code. As mentioned, it is under active development, including streamlining and optimization. There are a number of features currently built into PIMMS that are not documented here, either because they are not quite ready or because they are still in development. Again, we're working on finalizing all this up.

#### Who develops PIMMS?
PIMMS was developed by Alex Holehouse during his time in the [Pappu lab](http://pappulab.wustl.edu/). Alex is continuing to develop this package in [his own lab](http://holehouse.wustl.edu/), and there is a longer-term roadmap ahead, but for now one of the main goals is to get the main PIMMS paper out and published.

### Has PIMMS been used in any publications to date?
Why yes it has, thank you for asking! Please check out:

Martin, E. W.\*, Holehouse, A. S.\*, Peran, I.\*, Farag, M., Incicco, J. J., Bremer, A., Grace, C. R., Soranno, A., Pappu, R. V. & Mittag, T. Valence and patterning of aromatic residues determine the phase behavior of prion-like domains. Science 367, 694–699 (2020).

Boeynaems, S., Holehouse, A. S., Weinhardt, V., Kovacs, D., Van Lindt, J., Larabell, C., Van Den Bosch, L., Das, R., Tompa, P. S., Pappu, R. V. & Gitler, A. D. Spontaneous driving forces give rise to protein-RNA condensates with coexisting phases and complex material properties. Proc. Natl. Acad. Sci. U. S. A. 116, 7889–7898 (2019).

Martin, E. W.\*, Holehouse, A. S.\*, Grace, C. R., Hughes, A., Pappu, R. V. & Mittag, T. Sequence Determinants of the Conformational Properties of an Intrinsically Disordered Protein Prior to and upon Multisite Phosphorylation. J. Am. Chem. Soc. 138, 15323–15335 (2016).

#### This documentation seems a little... hastily put together?
Yeah... There's a pandemic going on, so, you know, we're trying... 

 `¯\_(ツ)_/¯ `

## Installation

**NB**: Installation assumes you have set up a correct `conda` environment with Python 3.7 or higher (any 3.7 is fine). Using 3.7 is important as there are some language features in 3.7 that we use that were not in earlier versions. 

If `conda` and `pip` are  new to you, there is a lot of documentation on this online, and I'd suggest taking a look at [this page here](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) as a first step. Assuming `conda` is set up installation _should_ be easy!  

#### Step one: Make sure dependencies are installed

Assuming `conda` is installed and your in the relevant environment, the first thing to do is ensure the channel **conda-forge** is available. 

Specifically run

	conda config --add channels conda-forge

Assuming this works correctly, next install some standard packages:

	conda install numpy scipy cython pandas 

And assuming these work, install `mdtraj`:
	
	# Then install mdtraj, which provides the xtc library backend - this is
	# what lets us write VMD-compatible trajectories
	conda install mdtraj
	
This should all work out of the box without issue. _At this stage_ if anything goes wrong its outside of my hands (although I'm happy to offer advice).

#### Step 2: Installing PIMMS

Assuming the packages above installed correctly, the next step is to actually install PIMMS.

As of version 0.1.23 we've made this even easier. Simply download the release candidate from the (secret) [internal page](https://www.dropbox.com/sh/ozrpqymi5se0xwk/AAC5Ng0BctaF9RrkSLhzVi7sa?dl=0) and then run:

	pip install pimms-<version>.tar.gz
	
For example, for release candidate 0.1.26.4 that would be
	
This _should_ just work!

If all seems to have gone off without a hitch **open a new terminal**, start up the `conda` environment you just installed PIMMS in, and run (from _any_ directory) the command:

	PIMMS --version
	
If it worked, you should see:

	version <current version number>
	
If this part fails, please contact Alex [alex.holehouse@wustl.edu] and we'll try and figure out what's goin' on.

This installation has been tested on both Linux and macOS.

## Usage

### Running a simulation
PIMMS simulations require two files

1. A **keyfile**, which defines the components of the simulation and all aspects of that simulation.
2. A **parameter file** which defines the interactions between distinct components. The **keyfile** also defines the location of the **parameter file**.

Simulations are run as follows:

	PIMMS -k <keyfile.kf>

For convenient some correctly formatted and annotated keyfiles and parameter files are available in the `/demo` directory. We recommend that you use these as a starting point for your own. These files have been annotated, but for completeness a more expansive description of the keywords is provided below. 

### Keyfile keywords

Keyfiles define everything about the system and simulation, including

* What polymers are present, their sequence, and how many of them are there
* How long the simulation should run for and how big the simulation box is
* The frequency of different analysis and output
* The move-set 

Below we outline the keywords you may wish to change. Note that there are additional keywords that control some advanced functionality, but we're still finalizing that behaviour.

#### System setup keywords
Keyword | Format (type) | Description
:---: | :---: | :---: 
DIMENSIONS | **INT** <br>(2 or 3) <br><br>`A x B`<br> or<br> `A x B x C`  | Size of the simulation box (in lattice units). 2D or 3D (defines if the simulation is a 2D or 3D simulation)
CHAIN | *See description* | One of the few multi component keywords in PIMMS and the only keyword that can appear multiple times, the `CHAIN` keyword defines a specific polymer chain and the number of that chain that will exist in the simulation. The format should be <br> <br> `CHAIN : N  {CHAIN IDENTIY}` <br><br> Where `N` defines the number of the chain and `{CHAIN IDENTITY}` gives polymer sequence in one-letter alphabet code. As an example <br> <br> `CHAIN : 20 QQQQQQQQQQ`<br><br> Would give 20 poly-glutamine polymers. In later versions of PIMMS we will be updating this to allow the reading of keyfiles that use three-letter codes 
TEMPERATURE | **FLOAT** | Simulation temperature to be used (units are arbitrary and depend on the units of the parameter file)
N_STEPS | **INT**|Number of main chain Monte Carlo steps
PARAMETER_FILE | **STRING** | Relative or absolute path of the parameter file
EQUILIBRATION | **INT** | Number of steps to be run as equilibration (i.e. before any analysis or trajectory output is generated)
PRINT_FREQ | **INT** | Step frequency at which the system status is printed to STDOUT
XTC_FREQ | **INT** | Step frequency at which the system system configuration is written to XTC file
EN_FREQ | **INT** | Step frequency at which the system energy is written to ENERGY.dat
SEED | **INT** | Random seed to allow reproducible runs
ENERGY_CHECK | **INT** | Step frequency at which global energy is recalculated and compared to the current current as determined locally on each step. All energy calculations are exact, so this is primarily for sanity checking. For large systems this can be an expensive operation.
NON_INTERACTING | **BOOL** | Boolean that defines if the Hamiltonian is non-interacting or not. This is a convenient way to generate "EV" ensembles for the same system configuration.

#### MOVE\_*  keywords
MOVE_ keywords define the frequency with which different moves are performed during the simulation. These values associated with each of these keywords must sum up to 1.0, and all must be defined.

Keyword | Format (type) | Description
:---: | :---: | :---: 
MOVE_CRANKSHAFT | **FLOAT** | Crankshaft moves drive local chain perturbations and are coded in optimized C (and so very fast). In general a large fraction of your simulation moveset should be these moves.
CRANKSHAFT_SUBSTEPS | **INT** | Defines a multiplier for the number of substeps performed. So each time a crankshaft move is selected, the underlying code performs `CRANKSHAFT_SUBSTEPS` multiplied by some scaling factors (defined by `CRANKSHAFT_MODE`) worth of moves for each bead in the system. In this way a single crankshaft move can actually encompass millions of individual MC moves!
CRANKSHAFT_MODE | **KEYWORD** | [`PROPORTIONAL`] defines how chain-length influences the multiplier for the crankshaft moves. Use `PROPORTIONAL`.
MOVE\_CHAIN_TRANSLATE |  **FLOAT** | Single chain rigid body translation
MOVE\_CHAIN_ROTATE |  **FLOAT** | Single chain rigid body rotation
MOVE\_CHAIN_PIVOT |  **FLOAT** | Chain pivot at a random potion
MOVE\_CLUSTER_TRANSLATE | **FLOAT** | Translate a randomly selected contigous cluster of chains
MOVE\_CLUSTER_ROTATE |  **FLOAT** |  Rotate a randomly selected contiguous cluster of chains

The following block of keywords define various options that control quench based simulations. In quench simulations, the simulation starts at a temperature defined by `QUENCH_START` and progressively decreases (or increases) to `QUENCH_END`. This is particularly useful to achieve convergence of complex systems, and is simply an annealing simulation.

Keyword | Format (type) | Description
:---: | :---: | :---: 
QUENCH_RUN | **BOOL** | Boolean (true or false) that defines if a quench run will be used. The 'quench' part of a question run always happens first in the simulation, although quenches can be low to high or high to low. 
QUENCH_FREQ |  **INT** | Frequency at which the temperature is updated.
QUENCH_STEPSIZE | **FLOAT** | Step (in temperature) that is taken each time the temperature is updated
QUENCH_START | **FLOAT** | Starting temperature
QUENCH_END | **FLOAT** | Ending temperature
QUENCH\_AS_EQUILIBRATION | **BOOL** | Boolean (true or false) that defines if the quench is treated as an equilibration period.

The following block of keywords define various options that control on-the-fly analysis done in PIMMS

Keyword | Format (type) | Description
:---: | :---: | :---: 
ANALYSIS_FREQ  | **INT** | Step frequency at which all default analysis is performed. This sets the default for all other types of analysis unless explicitly defined.
ANA_POL | **INT** | Step frequency at which polymeric analysis is performed.
ANA_INTSCAL | **INT** | Step frequency at which internal scaling analysis is performed.
ANA_DISTMAP | **INT** | Step frequency at which distance-map analysis is performed.
ANA_ACCEPTANCE | **INT** | Step frequency at which acceptance information is printed out
ANA\_INTER_RESIDUE | **INT** | Step frequency at which inter-residue interaction analysis is performed
ANA\_CLUSTER | **INT** | Step frequency at which cluster-based analysis is performed 
ANA\_RESIDUE_PAIRS | **INT** **INT** | Defines pairs of residues (i.e. "1 5") which are analyzed for inter-residue distance.



#### Forbidden keywords

The following block of keywords define various options that control temperature-sweep Metropolis Monte Carlo moves. This functionality is not fully ready so I do not recommend using for now until we confirm some key things! To keep things secret, we haven't even included a description of the keywords!

Keyword | Format (type) | Description
:---: | :---: | :---: 
TSMMC\_JUMP_TEMP | **FLOAT** | 
TSMMC\_STEP_MULTIPLIER | **INT** | 
TSMMC\_INTERPOLATION_MODE | **STRING** | 
TSMMC\_NUMBER_OF_POINTS | **INT** |
TSMMC\_FIXED_OFFSET | **INT** |

Similarly, there are some legacy moves which should also not be altered but must be included. Some of these may be removed for the final release, or updated, depending on our ongoing tests

Keyword | Format (type) | Description
:---: | :---: | :---: 
MOVE\_SLITHER |  **FLOAT** |  Slither the chain through the system (BROKEN, do not use). **Must be set to 0.0**
MOVE\_HEAD_PIVOT |  **FLOAT** | Pivot the head residue of the chain (this is a somewhat redundant move) **Must be set to 0.0**
MOVE\_CTSMMC |  **FLOAT** |  Single chain TSMMC **Must be set to 0.0**
MOVE\_MULTICHAIN_TSMMC  | **FLOAT** |  Multiple randomly selected chains undergo TSMMC **Must be set to 0.0**
MOVE\_SYSTEM_TSMMC | **FLOAT** | Entire system undergoes TSMMC **Must be set to 0.0**
MOVE\_RATCHET_PIVOT |  **FLOAT** | A single chain undergoes a directed pivot move. **Must be set to 0.0**

## Parameterfile
The parameter file defines the interactions experienced by the system. Note - EVERY bead defined on a `CHAIN` must be included in the parameter file and fully defined, no exceptions.

A parameterfile has three sections:

#### The angle section
PIMMS has a rudimentary 'backbone' angle term. The start of the parameter file includes a section where those angle strengths are defined. The format is

	ANGLE_PENALTY <residue name> X X X
	
For your purposes these `Xs` should be `0` - i.e. there is no angle restraint applied. We're still optimizing the implementation here, so I wouldn't use this (yet)

#### The bead-bead section
Next, for EVERY bead one must define the bead - bead interaction. PIMMS allows three different distance ranges for interactions (**short range**, **long range**, and **super long range**). These are defined by three distinct values. In this way, of you wanted to define A-B interaction as short, medium and long one might write

	A	B	-30	-10	-5

Which would mean beads A and B are directly adjacent to one another an interaction strength of -30 is realized. When they are 1 site apart -10 and two sites apart -5. These pairwise interactions must be defined for every bead in the system.

#### The bead-solvent section

Finally, we must ALSO define bead-solvent interactions explicitly, The solvent reserves the bead type `0`, such that bead solvent interactions are

	A	0	-5
	
This would say every solvent exposed face of bead A provides -5 energy.

#### Putting it all together
With this, a simple example of a parameterfile might be

	# angle section
	ANGLE_PENALTY <residue name> X X X
	ANGLE_PENALTY <residue name> X X X
	
	# bead-bead interactions
	A	A -5
	A	B -10
	B	B 0
	
	# bead-solvent interaction
	A	0 0
	B	B 0





### Output files
PIMMS generates a ton of output files. Below is  a brief overview of those files. All files are overwritten when the simulation starts, so simulations can be re-run in the same directory. Output files are subdivided below into distinct types. 


#### General 

The following files provide a description of different aspects of the system and simulation.

Filename | Explanation
:---: | :---: 
log.txt | Contains information on system setup. The current code underutilizes this and we are expanding the info that gets written here.
parameters\_used.prm | We've discovered it's useful to explicitly save which parameters were used WITH a simulation for cross-referencing in the future. This means one can 100% reproduce a simulation from the output files.
absolute\_energies\_of\_angles.txt | There are two modes that angle energies can be defined, one of which scales the energies by T, in which case the absolute energies depend on the simulation temperature. This file reports those absolute energies, and is honestly best used for debugging stuff

#### System output

The following files report on the 2D or 3D orientation of the simulated system, and the frequency at which they are written out is determined by `XTC_FREQ`.

Filename | Explanation
:---: | :---:
START.pdb | This file defines the topology of the simulation system and can be viewed in all good molecular viewers (i.e. VMD:  `vmd START.pdb`)
traj.xtc | This trajectory file defines the molecular evolution of the system, and can be viewed in most good molecular viewers in conjunction with the topology file [START.pdb] (i.e. for VMD: `vmd START.pdb traj.xtc`
frame.pdb | The is the instantaneous output of the system generate during simulation (this may disappear in the final release)

The following files report the step-dependent value of different aspects of the simulation. For _all_ the following files the first column reports on the relative step number (where each of the move-types is treated as a single step (i.e. ignoring sub-steps).

Filename | Explanation
:---: | :---:
ENERGY.dat | Reports on the instantaneous potential energy of the system.
MOVE_FREQS.dat | Reports the frequency with which each move type is proposed. Note that for crankshaft the TOTAL number of moves is reported (i.e. including subset MC steps) such that these values can be interpreted as the absolute number of accept/reject moves proposed. 
ACCEPTANCE.dat | Shows the same information as in `MOVE_FREQS.dat` but with accepted moves, allowing the user to back-calculate the acceptance ratio on any move type.
TOTAL\_MOVES.dat | Tracks the total number of moves (again, moves here include all sub-moves in crankshaft steps, so this is a count of the number of accept/reject operations performed.
STEPS\_PER\_SECOND.dat | For evaluating performance, this file writes out the wall-clock time that it takes to complete the `ANALYSIS_FREQ` number of steps. Useful for assessing performance.


#### Instantaneous single-chain analysis files
These files describe instantaneous analysis that will be most relevant for thinking about a single chain. This means output from this analysis is written at some regular interval as defined by ANALYSIS_FREQ or If simulations with many chains are run, the associated analysis is performed for every chain, which - if you don't care about it - can be computationally expensive. 


Filename | Explanation
:---: | :---:
RG.dat | Reports on the instantaneous radius of gyration of the chain(s).
ASPH.dat | Reports on the instantaneous asphericity of the chain(s).
END\_TO\_END\_DIST.dat| Reports on the instantaneous end-to-end distance of the chain(s)
RES\_TO\_RES\_DIST.dat |Reports on the instantaneous residue-to-residue distance, where residue pairs are defined by the keyword `ANA_RESIDUE_PAIRS`

#### Summary single-chain analysis files

These files describe analysis that will be most relevant for thinking about a single chain, but analysis that is reported ONLY at the end of the simulation and represents an ensemble average. 

If simulations with many chains are run, the associated analysis is performed for every chain, which - if you don't care about it - can be computationally expensive. The `ANA_*` keywords control the frequency of analysis moves - in particular the `ANA_POL` will scale all polymeric analysis, which can be useful when you're running simulations in which the behaviour of individual chains is of no interest.


Filename | Explanation
:---: | :---:
INTSCAL.dat| Reports the ensemble-average instantaneous internal scaling profile.
INTSCAL\_SQUARED.dat | Reports the ensemble-average instantaneous root-mean squared (RMS) internal scaling profile.
SCALING\_INFORMATION.dat| Result of analytical fits to the root-mean square scaling profile to extract the apparent scaling exponent and pre-factor that best describes the 1D scaling information. For details on this see the associated discussion in Peran, Holehouse *et al.* PNAS 2019
DISTANCE\_MAP.dat | Reports the ensemble-average inter-residue distance map. A `NON_INTERACTING` simulation can be run to generate the excluded-volume equivalent, allowing for a scaling map to be easily generated.


#### Instantaneous multi-chain output

These files describe analysis that will be most relevant for thinking about a multiple chains interacting together. Two types of multi-chain assemblies are analyzed - clusters and long-range (LR) clusters.

Filename | Explanation
:---: | :---:
CLUSTERS.dat| Lists the number of chains in each possible cluster. Chains not in a cluster are counted as clusters of "1" chain. Clusters are defined as chain in a continuous connected network that are 1 or 2 lattice sites away from one another.
LR\_CLUSTERS.dat| Lists the number of chains in each long-range cluster. Chains not in a long-range cluster are counted as clusters of "1" chain. Long-range clusters are defined as chain in a continuous connected network that are 1 or 2 lattice sites away from one another.
NUM\_[LR]\_CLUSTERS.dat| Total number of (long-range) clusters at a given moment.
[LR]\_CLUSTER\_RG.dat | Radius of gyration of each (long range) cluster. The clusters defined in `CLUSTERS.dat` (or `LR_CLUSTERS.dat`) map to those analyzed here.
[LR]\_CLUSTER\_ASPH.dat | Asphericity of each (long-range) cluster. The clusters defined in `CLUSTERS.dat` (or `LR_CLUSTERS.dat`) map to those analyzed here.
[LR]\_CLUSTER\_VOL.dat | Volume of each (long-range) cluster. The clusters defined in `CLUSTERS.dat` (or `LR_CLUSTERS.dat`) map to those analyzed here.
[LR]\_CLUSTER\_AREA.dat | Volume of each (long-range) cluster. The clusters defined in `CLUSTERS.dat` (or `LR_CLUSTERS.dat`) map to those analyzed here.
[LR]\_CLUSTER\_DEN.dat | Density of each (long-range) cluster. The clusters defined in `CLUSTERS.dat` (or `LR_CLUSTERS.dat`) map to those analyzed here.
[LR]\_CLUSTER\_RADIAL\_DENSITY\_PROFILE.dat | Radial density profile of each (long-range) cluster. The clusters defined in `CLUSTERS.dat` (or `LR_CLUSTERS.dat`) map to those analyzed here.
CHAIN\_\<n>\_CLUSTERS.dat | Fraction of each cluster defined in `CLUSTERS.dat` that consist of CHAIN <n>.
CHAIN\_\<n>\_LR\_CLUSTERS.dat | Fraction of each long-range cluster defined in `LR_CLUSTERS.dat` that consist of CHAIN <n>.

#### Other

Output files that (for now) are not useful/useable

Filename | Explanation
:---: | :---:
restart.pimms | PIMMS allows simulations to be restart, although this functionality is not yet ready. However, the `restart.pimms` is the only file required for restart.

## Example keyfiles
The pimms.tar.gz tarball comes with two examples under

	demo_keyfiles/demo_1  # multi-chain simulation
	demo_keyfiles/demo_2  # single-chain simulation

The keyfiles here (`KEYFILE.kf` are heavily annotated and a separate `readme.md` is found.

### Copyright

Copyright (c) 2015-2020, Alex Holehouse 
