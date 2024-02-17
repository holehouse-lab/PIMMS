# PIMMS: Polymer Interactions in Multi-component MixtureS
---

[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/pimms.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/pimms)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pimms/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pimms/branch/master)


--- 

###### PIMMS version v0.1.35 (February 2024)

# Preamble
PIMMS is still in development, but _in general_ the `master` branch on this repository can be considered (mostly) stable. As of version 0.1.35, we consider this to be the 'gold' release, and no major changes are expected before publication.

As you use this version of PIMMS please report any/all issues where things: 
 
1. Are wrong/don't make sense
2. Don't behave in a way you expect
3. Errors/exceptions when you wouldn't expect them

If you have access to the [github repository](https://github.com/holehouse-lab/PIMMS) please log issues in the [issue tracker](https://github.com/holehouse-lab/PIMMS/issues). If you don't have access,  please contact Alex directly.

There _shouldn't_ be any bugs in this version... but there could be. So we're working on it! This includes the development and deployment of a large unit test suite, which takes time, but we're working hard! 

# Background
### What is PIMMS?
PIMMS is a lattice-based simulation engine that allows both 2D and 3D simulations to be performed. Useful features include:

1. Easy to use! Upon installation a command-line executable (`PIMMS`) is available, should be in you `$PATH` variable, and can be used to run simulations. No messing around, it (should) just work! 
2. Easy to define interaction parameters through a simple parameter file (example included in `/demo_keyfiles/demo_1/params.prm`)
3. Easily run fast 2D or 3D lattice based simulations
4. Run simulations with many distinct components
5. Run simulations of a single homo or heteropolymer 
5. Run simulations of many copies of polymers to explore phase behavior
7. Drive interactions over three distinct length-scales
8. Various other things

### How is PIMMS written?
PIMMS is written almost fully in Python (=>3.7), with the most computationally intensive parts written in fully optimized `Cython` that compiles down to native C. We're still ironing out kinks, but by having most of the complex behavior in Python, maintenance and development is fast and efficient. However, certain functionality (i.e., analysis of large systems) is, as a result, disproportionately expensive vs. the actual simulation, so you may wish to alter the frequency at which certain analysis routines are performed based on your interests.

PIMMS is a relatively large codebase of ~20K lines of (mostly) Python code. As mentioned, it is under active development, including streamlining and optimization. There are several features currently built into PIMMS that are not documented here, either because they are not quite ready or because they are still in development. Again, we're working on finalizing all this up.

### Who develops PIMMS?
An initial version of PIMMS was developed by Alex Holehouse during his time in the [Pappu lab](http://pappulab.wustl.edu/). Since starting [his own lab](http://holehouse.wustl.edu/), much of PIMMS has been re-written, and Dr. Ryan Emenecker has joined as a core developer. Maintenance of PIMMS is now maintained exclusively by the Holehouse lab

### Has PIMMS been used in any publications to date?
Why yes, it has, thank you for asking! Please check out:

Alston, J. J. & Soranno, A. Condensation goes viral: a polymer physics perspective. J. Mol. Biol. 167988 (2023).

Soranno, A., Incicco, J. J., De Bona, P., Tomko, E. J., Galburt, E. A., Holehouse, A. S. & Galletto, R. Shelterin Components Modulate Nucleic Acids Condensation and Phase Separation in the Context of Telomeric DNA. J. Mol. Biol. 434, 167685 (2022).

Sankaranarayanan, M., Emenecker, R. J., Wilby, E. L., Jahnel, M., Trussina, I. R. E. A., Wayland, M., Alberti, S., Holehouse, A. S. & Weil, T. T. Adaptable P body physical states differentially regulate bicoid mRNA storage during early Drosophila development. Dev. Cell 56, 2886–2901.e6 (2021).

Moses, D., Yu, F., Ginell, G. M., Shamoon, N. M., Koenig, P. S., Holehouse, A. S. & Sukenik, S. Revealing the Hidden Sensitivity of Intrinsically Disordered Proteins to their Chemical Environment. J. Phys. Chem. Lett. 11, 10131–10136 (2020).

Holehouse, A. S., Ginell, G. M., Griffith, D. & Böke, E. Clustering of Aromatic Residues in Prion-like Domains Can Tune the Formation, State, and Organization of Biomolecular Condensates. Biochemistry 60, 3566–3581 (2021).

Martin, E. W.\*, Holehouse, A. S.\*, Peran, I.\*, Farag, M., Incicco, J. J., Bremer, A., Grace, C. R., Soranno, A., Pappu, R. V. & Mittag, T. Valence and patterning of aromatic residues determine the phase behavior of prion-like domains. Science 367, 694–699 (2020).

Boeynaems, S., Holehouse, A. S., Weinhardt, V., Kovacs, D., Van Lindt, J., Larabell, C., Van Den Bosch, L., Das, R., Tompa, P. S., Pappu, R. V. & Gitler, A. D. Spontaneous driving forces give rise to protein-RNA condensates with coexisting phases and complex material properties. Proc. Natl. Acad. Sci. U. S. A. 116, 7889–7898 (2019).

Martin, E. W.\*, Holehouse, A. S.\*, Grace, C. R., Hughes, A., Pappu, R. V. & Mittag, T. Sequence Determinants of the Conformational Properties of an Intrinsically Disordered Protein Prior to and upon Multisite Phosphorylation. J. Am. Chem. Soc. 138, 15323–15335 (2016).

### Is this the final documentation?
No. We are actively working on a full Sphinx-based readthedocs documentation suite for PIMMS, but for now, this readme file serves as the core PIMMS documentation. 

# Installation

**NB**: Installation assumes you have set up a correct `conda` environment with Python 3.7 or higher (any 3.7 is fine). Using 3.7 or higher important as there are some language features in 3.7 that we use that were not in earlier versions that PIMMS requires. We recommend using `Python 3.10`.

If `conda` and `pip` are  new to you, there is a lot of documentation on this online, and I'd suggest taking a look at [this page here](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) as a first step. Assuming `conda` is set up, installation _should_ be easy!  

We have also put together some introductory material [which can be found here]().

### Step one: Make sure dependencies are installed

Assuming `conda` is installed and your in the relevant environment, the first thing to do is ensure the channel **conda-forge** is available. 

Specifically run

	conda config --add channels conda-forge
	
We then recommend creating a clean environment with Python 3.10

	conda create -n pimms  python=3.10 -y
	conda activate pimms
	
Note you're welcome to install PIMMS into an existing environment if you want, but we strongly recommend doing this FIRST to ensure things actually install correctly into a vanilla and empty environment. Having a dedicated environment also avoids any dependency clashes and is generally a safer bet.

Assuming this works correctly, next install some standard packages:

	conda install numpy scipy cython pandas 

And assuming these work, install `mdtraj`:
	
	# Then install mdtraj, which provides the xtc library backend - this is
	# what lets us write VMD-compatible trajectories
	conda install mdtraj
	
This should all work out of the box without issue. _At this stage_ if anything goes wrong it's outside of my hands (although I'm happy to offer advice).

### Step two: Installing PIMMS

Assuming the packages above installed correctly, the next step is to actually install PIMMS. 

To do this, take the `pimms-0.1.35.tar.gz` archive you unpacked to access this file and run

	pip install pimms-0.1.35.tar.gz
	
This _should_ just work! 

If all seems to have gone off without a hitch **open a new terminal**, start up the `conda` environment you just installed PIMMS in, and run (from _any_ directory) the command:

	PIMMS --version
	
If it worked, you should see:

	version <current version number>
	
NOTE - the VERY first time you do this may take 5-10 seconds due to the internal Python environment things initializing, but after that should be basically instantaneous.	
	
If this part fails, please contact Alex [alex.holehouse@wustl.edu] and we'll try and figure out what's goin' on.

This installation has been tested and works on both Linux and macOS. If someone has a Windows machine and wants to test this out they are more the welcome, but, PIMMS has (AFAIK) never been run on Windows so I would anticipate things not working well out of the box...

### Installing the current development version
*Relevant if you have access to the PIMMS GitHub repository*
Alternatively you can clone this repository and install from source. I strongly recommend cloning into a sensible location - e.g. not in ~/Downloads but a directory location that makes sense, e.g. for me PIMMS is found in

	/home/alex/Dropbox/tools/pimms/

But you (presumably) have your own organizational approach to directory management, so create a PIMMS directory somewhere that makes sense (!).

Once this directory is created/ready, you can run:

	git clone git@github.com:holehouse-lab/PIMMS.git
	
Then, navigate to the main directory (where `setup.py` is) and run:

	pip install -e . --upgrade --force-reinstall
	
To install directly from source.

Note we include the `--upgrade` and `--force-reinstall` flags to ensure things get built if you pull an update from GitHub and want to reinstall. 
	
To receive updates to code the from GitHub, simply run

	git pull
	
This will download any changes made to the remote repository to your local clone; you can then re-install by re-running the `pip install` command above.

# Usage

### Running a simulation
PIMMS simulations require two files

1. A **keyfile**, which defines the components of the simulation and all aspects of that simulation.
2. A **parameter file** which defines the interactions between distinct components. The **keyfile** also defines the location of the **parameter file**.

Simulations are run as follows:

	PIMMS -k <keyfile.kf>

For convenience some correctly formatted and annotated keyfiles and parameter files are available in the `/demo` directory. We recommend that you use these as a starting point for your own. These files have been annotated, but for completeness a more expansive description of the keywords is provided below. 

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
RESIZED_EQUILIBRATION **[OPTIONAL]**	| **INT** <br>(2 or 3) <br><br>`A x B`<br> or<br> `A x B x C`  | Defines alternative simulation dimensions to be used during equilibration; e.g., if you wish to equilibrate chains at a higher concentration to facilitate single condensate forming. MUST be smaller than the dimensions defined by the DIMENSIONS keyword.
CHAIN | *See description* | One of the few multi component keywords in PIMMS and the only keyword that can appear multiple times, the `CHAIN` keyword defines a specific polymer chain and the number of that chain that will exist in the simulation. The format should be <br> <br> `CHAIN : N  {CHAIN IDENTIY}` <br><br> Where `N` defines the number of the chain and `{CHAIN IDENTITY}` gives polymer sequence in one-letter alphabet code. As an example <br> <br> `CHAIN : 20 QQQQQQQQQQ`<br><br> Would give 20 poly-glutamine polymers. In later versions of PIMMS we will be updating this to allow the reading of keyfiles that use three-letter codes 
TEMPERATURE | **FLOAT** | Simulation temperature to be used (units are arbitrary and depend on the units of the parameter file)
N_STEPS | **INT**|Number of main chain Monte Carlo steps
PARAMETER_FILE | **STRING** | Relative or absolute path of the parameter file
EQUILIBRATION | **INT** | Number of steps to be run as equilibration (i.e. before any analysis or trajectory output is generated)
HARDWALL **[OPTIONAL]**	 | **BOOL** | Boolean flag set to True or False that defines if a hardwall boundary is used or not. By default periodic boundary conditions (PBC) are used, but if hardwall is set to true the edges of the simulation box are reflective with an infinitely repulsive potential. **Default = False**
PRINT_FREQ | **INT** | Step frequency at which the system status is printed to STDOUT
XTC_FREQ | **INT** | Step frequency at which the system system configuration is written to XTC file
EN_FREQ | **INT** | Step frequency at which the system energy is written to ENERGY.dat
SEED **[OPTIONAL]**	 | **INT** | Random seed to allow reproducible runs
ENERGY_CHECK | **INT** | Step frequency at which global energy is recalculated and compared to the current current as determined locally on each step. All energy calculations are exact, so this is primarily for sanity checking. For large systems this can be an expensive operation.
NON_INTERACTING **[OPTIONAL]**	 | **BOOL** | Boolean that defines if the Hamiltonian is non-interacting or not. This is a convenient way to generate "EV" ensembles for the same system configuration. **Default = False**
ANGLES_OFF | **BOOL** |Boolean flag set to True or False that defines if angle potentials are to be used or not. If set to False (or not set), angles from the parameter file will be used. If set to True, angles are ignored and parameter files do not need to define angles. **Default = False**
EXPERIMENTAL_FEATURES **[OPTIONAL]**	 | **BOOL** | Boolean flag set to True or False that defines if experimental features are allowed. Strongly recommend if your name is not Ryan or Alex to keep this at False, and even then, check if your last name is Holehouse or Emenenecker because IF NOT you should still keep it False, probably. **Default = False**
CASE\_INSENSITIVE\_CHAINS **[OPTIONAL]**	 | **BOOL** | Boolean flag which, if set to False, means that chain sequence is case sensitive. By default, this is True, which means upon reading a keyfile chains are converted to upper case. However, sometimes you may wish for more unique beads, in which case a lower-case chain can be useful. **Default = True**
AUTOCENTER **[OPTIONAL]**	| **BOOL** | Only relevant for single-chain simulations, but this flag, if set to True, ensures every frame of the resulting trajectory is centered in the middle of the box. This is especially useful if you want to avoid the need to align a trajectory after for analysis or visualization. **Default = False**
LATTICE_TO_ANGSTROMS **[OPTIONAL]**	| **FLOAT** | Defines the conversion for lattice units to Angstroms for the output trajectory file that's generated. NB This ONLY influences the XTC trajectory, not any of the internal analysis which is always returned in lattice units. **Default = 3.65**

#### MOVE\_*  keywords
MOVE_ keywords define the frequency with which different moves are performed during the simulation. The values associated with each of these keywords must sum up to 1.0, and all must be defined.

Keyword | Format (type) | Description
:---: | :---: | :---: 
MOVE_CRANKSHAFT | **FLOAT** | Crankshaft moves drive local chain perturbations and are coded in optimized C (and so very fast). In general, a large fraction of your simulation moveset should be these moves.
CRANKSHAFT_SUBSTEPS | **INT** | Defines a multiplier for the number of substeps performed. So each time a crankshaft move is selected, the underlying code performs `CRANKSHAFT_SUBSTEPS` multiplied by some scaling factors (defined by `CRANKSHAFT_MODE`) worth of moves for each bead in the system. In this way a single crankshaft move can actually encompass millions of individual MC moves!
CRANKSHAFT_MODE | **KEYWORD** | [`PROPORTIONAL`] defines how chain length influences the multiplier for the crankshaft moves. Use `PROPORTIONAL`.
MOVE\_CHAIN_TRANSLATE |  **FLOAT** | Single chain rigid body translation
MOVE\_CHAIN_ROTATE |  **FLOAT** | Single chain rigid body rotation
MOVE\_CHAIN_PIVOT |  **FLOAT** | Chain pivot at a random potion
MOVE\_CLUSTER_TRANSLATE | **FLOAT** | Translate a randomly selected contiguous cluster of chains
MOVE\_CLUSTER_ROTATE |  **FLOAT** |  Rotate a randomly selected contiguous cluster of chains

The following block of keywords defines various options that control quench-based simulations. In quench simulations, the simulation starts at a temperature defined by `QUENCH_START` and progressively decreases (or increases) to `QUENCH_END`. This is particularly useful to achieve convergence of complex systems and is simply an annealing simulation.

Keyword | Format (type) | Description
:---: | :---: | :---: 
QUENCH_RUN | **BOOL** | Boolean (true or false) that defines if a quench run will be used. The 'quench' part of a question run always happens first in the simulation, although quenches can be low to high or high to low. **Default = False**
QUENCH_FREQ |  **INT** | Frequency at which the temperature is updated.
QUENCH_STEPSIZE | **FLOAT** | Step (in temperature) that is taken each time the temperature is updated
QUENCH_START | **FLOAT** | Starting temperature
QUENCH_END | **FLOAT** | Ending temperature
QUENCH\_AS\_EQUILIBRATION | **BOOL** | Boolean (true or false) that defines if the quench is treated as an equilibration period.

The following block of keywords defines various options that control on-the-fly analysis done in PIMMS

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

#### Save keywords

Keywords for changing how PIMMS saves your trajectory file. 

Keyword | Format (type) | Description
:---: | :---: | :---: 
SAVE\_AT\_END | **BOOL** | Boolean (true or false) that determines whether PIMMS saves your .xtc file at the end of the simulation or saves at each 'save step'. The default is False. If set to true, this will mean the simulation has more sustained RAM usage, but may increase simulation performance substantially depending on hardware configuration and setup.
SAVE_EQ | **BOOL** | Boolean (true or false) that determines whether PIMMS saves trajectory frames for the equilibration steps of a simulation. The default is True. If set to False, PIMMS begins to save your trajectory frames *after* the equilibration steps have completed. 

#### Forbidden keywords

The following block of keywords defines various options that control temperature-sweep Metropolis Monte Carlo moves. This functionality is not fully ready, so I do not recommend using it for now until we confirm some key things! To keep things secret, we haven't even included a description of the keywords!

Keyword | Format (type) | Description
:---: | :---: | :---: 
TSMMC\_JUMP_TEMP | **FLOAT** | 
TSMMC\_STEP_MULTIPLIER | **INT** | 
TSMMC\_INTERPOLATION_MODE | **STRING** | 
TSMMC\_NUMBER_OF_POINTS | **INT** |
TSMMC\_FIXED_OFFSET | **INT** |

Similarly, there are some legacy moves that should also not be altered but must be included. Some of these may be removed for the final release or updated, depending on our ongoing tests

Keyword | Format (type) | Description
:---: | :---: | :---: 
MOVE\_SLITHER |  **FLOAT** |  Slither the chain through the system (BROKEN, do not use). **Must be set to 0.0**
MOVE\_HEAD_PIVOT |  **FLOAT** | Pivot the head residue of the chain (this is a somewhat redundant move) **Must be set to 0.0**
MOVE\_CTSMMC |  **FLOAT** |  Single chain TSMMC **Must be set to 0.0**
MOVE\_MULTICHAIN_TSMMC  | **FLOAT** |  Multiple randomly selected chains undergo TSMMC **Must be set to 0.0**
MOVE\_SYSTEM_TSMMC | **FLOAT** | Entire system undergoes TSMMC **Must be set to 0.0**
MOVE\_RATCHET_PIVOT |  **FLOAT** | A single chain undergoes a directed pivot move. **Must be set to 0.0**

## Parameterfile
The parameter file defines the interactions experienced by the system. Note - EVERY bead defined on a `CHAIN` must be included in the parameter file and fully defined, with no exceptions.

A parameter file has three sections:

#### The angle section
PIMMS has a rudimentary 'backbone' angle term. The start of the parameter file includes a section where those angle strengths are defined. The format is

	ANGLE_PENALTY <residue name> X X X
	
For your purposes, these `Xs` should be `0` - i.e., there is no angle restraint applied. We're still optimizing the implementation here, so I wouldn't use this (yet)

#### The bead-bead section
Next, for EVERY bead, one must define the bead - bead interaction. PIMMS allows three different distance ranges for interactions (**short range**, **long range**, and **super long range**). These are defined by three distinct values. In this way, of you wanted to define A-B interaction as short, medium, and long, one might write

	A	B	-30	-10	-5

This would mean beads A and B are directly adjacent to one another, and an interaction strength of -30 is realized. When they are one site apart -10 and two sites apart -5. These pairwise interactions must be defined for every bead in the system.

#### The bead-solvent section

Finally, we must ALSO define bead-solvent interactions explicitly; The solvent reserves the bead type `0`, such that bead-solvent interactions are

	A	0	-5
	
This would say every solvent-exposed face of bead A provides -5 energy.

#### Putting it all together
With this, a simple example of a parameter file might be

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
absolute\_energies\_of\_angles.txt | There are two modes that angle energies can be defined, one of which scales the energies by T, in which case the absolute energies depend on the simulation temperature. This file reports those absolute energies and is honestly best used for debugging stuff

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
TOTAL\_MOVES.dat | Tracks the total number of moves (again, moves here include all sub-moves in crankshaft steps, so this counts the number of accept/reject operations performed.
STEPS\_PER\_SECOND.dat | For evaluating performance, this file writes out the wall-clock time it takes to complete the `ANALYSIS_FREQ` number of steps. Useful for assessing performance.


#### Instantaneous single-chain analysis files
These files describe the instantaneous analysis that will be most relevant for thinking about a single chain. This means output from this analysis is written at some regular interval as defined by ANALYSIS_FREQ or, if simulations with many chains are run, the associated analysis is performed for every chain, which - if you don't care about it - can be computationally expensive. 

Filename | Explanation
:---: | :---:
RG.dat | Reports on the instantaneous radius of gyration of the chain(s).
ASPH.dat | Reports on the instantaneous asphericity of the chain(s).
END\_TO\_END\_DIST.dat| Reports on the instantaneous end-to-end distance of the chain(s)
RES\_TO\_RES\_DIST.dat |Reports on the instantaneous residue-to-residue distance, where residue pairs are defined by the keyword `ANA_RESIDUE_PAIRS`

#### Summary single-chain analysis files

These files describe the analysis that will be most relevant for thinking about a single chain, but the analysis that is reported ONLY at the end of the simulation represents an ensemble average. 

If simulations with many chains are run, the associated analysis is performed for every chain, which - if you don't care about it - can be computationally expensive. The `ANA_*` keywords control the frequency of analysis moves - in particular, the `ANA_POL` will scale all polymeric analysis, which can be useful when running simulations in which the behavior of individual chains is of no interest.


Filename | Explanation
:---: | :---:
INTSCAL.dat| Reports the ensemble-average instantaneous internal scaling profile.
INTSCAL\_SQUARED.dat | Reports the ensemble-average instantaneous root-mean squared (RMS) internal scaling profile.
SCALING\_INFORMATION.dat| Result of analytical fits to the root-mean square scaling profile to extract the apparent scaling exponent and pre-factor that best describes the 1D scaling information. For details on this, see the associated discussion in Peran, Holehouse *et al.* PNAS 2019
DISTANCE\_MAP.dat | Reports the ensemble-average inter-residue distance map. A `NON_INTERACTING` simulation can be run to generate the excluded-volume equivalent, allowing for a scaling map to be easily generated.


#### Instantaneous multi-chain output

These files describe the analysis that will be most relevant for thinking about multiple chains interacting together. Two types of multi-chain assemblies are analyzed - clusters and long-range (LR) clusters.

Filename | Explanation
:---: | :---:
CLUSTERS.dat| Lists the number of chains in each possible cluster. Chains not in a cluster are counted as clusters of "1" chain. Clusters are defined as chains in a continuous connected network that are 1 or 2 lattice sites away from one another.
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

## Changelog

#### 0.1.35 (February 2024)
* Updated docs
* Updated demo keyfiles
* Fixed mismatch in integer type breaking PIMMS entirely...

#### 0.1.35 (January 2024)
* Major update to Cython backend to facilitated better control over memory usage. In previous versions, PIMMS defined all back-end grids (chain grids and type grids) as [n x n x n] matrices, where each element was a 64-bit number. Because of the cubic term here, as grids become larger the memory footprint associated with PIMMS becomes very large; for a [200 x 200 x 200] grid, the memory footprint can reach 100s of MBs. In version 0.1.35, the backend memory management has been made dynamic; that is, we can compile PIMMS versions that specify the number of bits associated with the elements in the grids. Right now, the default version compiles with 64-bit integers still, but to recompile with a smaller memory footprint you can change just two flags in `CONFIGS.py` and the newly created `cython_config.pxd`; specifically :

		ctypedef cnp.int64_t NUMPY_INT_TYPE 
		
		# can be changed to
		
		ctypedef cnp.int32_t NUMPY_INT_TYPE 
		
		# or 
		
		ctypedef cnp.int16_t NUMPY_INT_TYPE 
		
	While in `CONFIGS.py` 
	
		NP_INT_TYPE = np.int64
		
		# can be changed to
		
		NP_INT_TYPE = np.int32
		
		# or
		
		NP_INT_TYPE = np.int16
		
	Right now, we continue to default to 64-bit numbers as this is test driven, but ultimately the plan is to transition to a 16-bit backend which basically reduces the memory footprint down to 25% of what it would have been with a 64-bit backend. To what extent this improves performance (steps/second) is unclear, but it HUGELY helps running many parallel jobs on many CPU systems where we actually become memory limited!
	
*  In addition to the memory re-write, we re-wrote `delete_pbc_pairs()` in `inner_loops_hardwall.pyx` to substantially improve performance by fully typing the function - for non 64-bit numbers this adds a ~8x improvement in performance, and maybe 1-2x for native (64-bit) memory implementations.
* Despite the big improvement in memory utilization, arguably the most important update in version 0.1.35 is a re-write of how trajectory saving is done. In particular, we previously had a punishingly inefficient approach for writing new trajectory frames that was so stupid it almost makes you wonder if it was a deliberate act of sabotage by someone. In any case, we (Ryan) has re-written this code to [firstly] ensure trajectory writing is done in a single output operation of XTC only data (instead of the 3x I/O operations we had previously, [don't ask...]). 
* Beyond this update, we (Ryan) also added the `SAVE_AT_END` keyword (default = False). If set to True, this means the simulation only writes the entire XTC at the end of the simulation. If you are worried about simulations crashing this is not ideal. However, where this is not a major concern, avoiding many I/O operations offers big gains, especially for larger systems.
* Added `SAVE_EQ` keyword. Default True. If set to False, equilibrations steps are not saved. This works for both `RESIZED_EQUILIBRATION` experiments (an eq.traj file is still made but it only contains a single frame) and when `RESIZED_EQUILIBRATION` is not used (standard sims). When `RESIZED_EQUILIBRATION` is not used, PIMMS will begin saving (or updating the trajobj if `SAVE_AT_END` is set to True) after the equilibration step but does not save before. 
* KEYWORDS ADDED: `SAVE_AT_END`, `SAVE_EQ` (discussed above).
* Version 0.1.35 is the final architectural change prior to the bump to 0.2.0 which will be the first live PIMMS release. Get psyched. 

#### 1.0.34 (September 2023)
* Major update to Cython backend to improve performance. All numpy arrays are now passed as memory views instead of as new arrays, which reduces the overhead on large arrays substantially 
* This big re-write has been tested extensively without any issues identified
* The default lattice-to-realspace value (`LATTICE_TO_ANGSTROMS`) has been updated from 4.0 nm to 3.65 nm
* KEYWORDS ADDED: `CASE_INSENSITIVE_CHAINS`, `AUTOCENTER`

#### 0.1.33-patch
* Update so logfile is always a new file instead of appended to (pimmslogger.py)

#### 0.1.33 (October 2022)
* Restructured to define the DEFAULTS dictionary which sets and explains default parameters for keyfiles. This means default options are encoded directly in
* Updated internal documentation
* Added reduceD_printing mode
* KEYWORDS ADDED: `REDUCED_PRINTING`
* Fixed bug which could lead to an error when a non-essential keyword was unset
* If residues are unknown to single letter-to-three letter conversion ensure that first character in the unknown residue type is not a number, because this causes PDB readers to fail.

#### 0.1.32 (April 2022)
* Added `EXTRA_CHAINS` keyword
* Fixed bug in how pdb chain ID was being written for `TER` lines (always using chain A)
* Added and improved internal `RestartObject` code and functionality (including improved parsing)
* Improved information printed when a RESTART file is used to make it easier to see what is going on.
* Added CONECT records to output PDB files, so bonds between chains are easily visualized
* Added `LATTICE_TO_ANGSTROMS` keyword such that PDB file dimensions are controllable. Default=4 (same as before) so this will not change anything compared to prior simulations.
* Improved code documentation and removed `xtc_utils` due to redundancy.

#### 0.1.31 (March 2022)
* Changed so PDB chains defined by the internal chainType - that is, all chains of same type have same PDB chain ID, which is convenient for visualization 

### Copyright

Copyright (c) 2015-2023, Alex Holehouse 
