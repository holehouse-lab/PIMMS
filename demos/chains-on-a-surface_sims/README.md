## PIMMS Simulations on a Surface ##

This guide walks through the setup and execution of PIMMS coarse-grained simulations in which the chains are affixed to a surface. The specific setup described here does the following:

1. generates a 'surface' of a single bead type based on simulation box dimensions (*d*) at y = 1 
2. positions *n* identical chains on the perimeter of a circle that fits inside the simulation box
3. creates `surface_restart_attached.pimms` file that will be called in the KEYFILE for simualtion setup
	NB: surface beads are frozen, whereas 'protein' chains are permitted to move

### What you'll need:

1. stable PIMMS build (see [installation docs](https://github.com/holehouse-lab/PIMMS))
2. `build_surface_attach.ipynb` - generates `surface_restart_attached.pimms`
3. `KEYFILE.kf` - specifies all simulation details
4. `freezefile.in` - specifies which chains to be immobile (surface chain)
5. parameter file containing appropriate intereaction definitions
	- the `.prm` file used here is a modified version of [Ryan's amino acid parameters](https://github.com/holehouse-lab/PIMMS/blob/master/demo_keyfiles/single_chain_protein_demo/gcf_rje23_v14.prm)
	- the version for surfaces used in this tutorial is `gcf_rje23_v14_ETU_surf_JO.prm`
		- 'JO' refers to the two new bead types introduced for the suface bead type (J) and the chain bead type that attaches to the surface (O)
6. your simulation details: desired chain sequences, # of chains, box size, etc. 
	- this tutorial describes (and has been tested for) the setup of 50, 30-residue chains in a 50x50x50 box

### Toy system: secreted salivary protein, Statherin (P02808)

Statherin sequence excluding signal sequence (residues 1-43) = `DSSEEKFLRRIGRFGYGYGPYQPVPEQPLYPQPYQPQYQQYTF`

From published experiments on this protein, residues 1-14 (`DSSEEKFLRRIGRF`) are responsible for contacting the surface of the tooth; the last residue that stably contacts the tooth/hydroxyapatite surface is F14. We will simulate the C-terminal Y/Q/G-rich sequence, residues 15-43 (`GYGYGPYQPVPEQPLYPQPYQPQYQQYTF`), and use the residue O in position 14 as the site of attachment to the surface. 

The sequence to be simulated is:
	`OGYGYGPYQPVPEQPLYPQPYQPQYQQYTF`
where O is a unique bead type for which we will specify interaction parameters with all other bead types (see below).

### 2. Generating `surface_restart_attached.pimms` 
Use the iPython notebook `build_surface_attach.ipynb` to build the 'surface' chain and each of the protein chains in arranged on the perimeter of a circle for symmetric chain placement. 

Before running, the user must specify: the box side length ('d'), the number of non-surface chains to build ('num_chains'), the padding around the circle within the box side ('shrink'), and the 1-letter amino acid sequence ('seqin').

The resulting file from running this notebook will contain 1 long chain of 'surface' beads (CHAINID 1) and num_chains copies of the sequence to be simulated on the surface (CHAINIDs 2 thru 51).

##### USER INPUTS:
	d = 50      # 50x50x50 size is used in this tutorial
	num_chains = 50      # 50 chains are used in this tutorial
	shrink = 3     # 'shrink' sets the diameter of the generated circle to be (d - shrink)
	seqin = 'OGYGYGPYQPVPEQPLYPQPYQPQYQQYTF'    # string of 1-letter amino acids
*For additional information about how this script builds the surface, please see cell-by-cell comments*

### 3. Editing `KEYFILE.kf`
The simulation options and keywords in the KEYFILE are described in depth [here](https://github.com/holehouse-lab/PIMMS), but we will survey a few particularly important settings below. I usually save the KEYFILE in a subdirectory for a given simualtion attempt, such that the parameter and freeze files are in the directory above. 

##### Input definitions:
	DIMENSIONS : 50 50 50            # Cube dimensions - MUST be the same as that used to generate restart file
	PARAMETER_FILE : ../gcf_rje23_v14_ETU_surf_JO.prm       # path to parameter file relative to KEYFILE.kf location

##### Chain definitions:
	RESTART_FILE : ../surface_restart_attached.pimms		# path to location of restart file made in the iPython notebook
	HARDWALL : True                  # required TRUE for surface simulations
	FREEZE_FILE : ../freezefile.in    	# path to 'freeze_file.in' (we will make this file in the next step)
	WRITE_CHAIN_TO_CHAINID : TRUE
*Note that we do not define any chains in the keyfile itself. The fixed surface and the 50 mobile chains are encoded in the `surface_restart_attached.pimms` file.*

##### Simulation definitions:
	N_STEPS       :  5000            # Number of simulation steps
	EQUILIBRATION :   250            # Number of steps for equilibration
	TEMPERATURE   :  50              # Simulation temperature - units are arbitrary; amino acid parameters are optimized for 50
	EXPERIMENTAL_FEATURES : True

	SAVE_AT_END : True     	# this prevents memory-intensive writeout until the end of the simulation

##### Output information:
	PRINT_FREQ    :	 1              # Frequency at which info is printed to stdout
	XTC_FREQ      :  10              # Frequency at which traj.xtc file is written to
	EN_FREQ       :  5               # Frequency at which energy file is written to
	RESTART_FREQ  : 1000000        # number of steps between creation of 'restart' files (i.e., never for simualtion of 5000 steps)
	ANALYSIS_FREQ : 1000000		   # number of steps between analyses performed (i.e., never for simulation of 5000 steps)

### 4. Creating `freezefile.in`
For a simulation with only one immobile chain, the *freezefile* is very straightforward:

	## freezefile for a surface of beads
	C 1
where C specifies that we are freezing a 'CHAIN' and 1 is the 1-indexed CHAINID. Yes, this is an entire text file with one line it. Any line following an octothorp will be treated as a comment. Because the frozen chains are called based on CHAINID, the same *freezefile* can be used for any simualtions in which CHAINID 1 should be immobile.

### 5. Modifying the `.prm` file for surface simulations
We will go through the parameter file section by section to ensure completeness.

##### Angle penalties:
	ANGLE_PENALTY   A       160     20      0
	ANGLE_PENALTY   M       120     100     0
	ANGLE_PENALTY   I       180     90      0
	ANGLE_PENALTY   V       160     70      0
	ANGLE_PENALTY   L       140     90      0
	ANGLE_PENALTY   W       200     70      0
	ANGLE_PENALTY   F       200     60      0
	ANGLE_PENALTY   Y       180     70      0
	ANGLE_PENALTY   T       180     50      0
	ANGLE_PENALTY   C       140     60      0
	ANGLE_PENALTY   P       300     160     0
	ANGLE_PENALTY   Q       120     100     0
	ANGLE_PENALTY   N       180     50      0
	ANGLE_PENALTY   S       180     20      0
	ANGLE_PENALTY   G       100     10      0
	ANGLE_PENALTY   H       140     80      0
	ANGLE_PENALTY   E       120     110     0
	ANGLE_PENALTY   D       140     60      0
	ANGLE_PENALTY   K       160     70      0
	ANGLE_PENALTY   R       200     70      0
	ANGLE_PENALTY   X       0       0       0
	ANGLE_PENALTY   J       0       0       0          # bead type J (surface) has no angle penalties
	ANGLE_PENALTY   O       0       0       0          # bead type O (surface sticker) has no angle penalties

##### Bead-bead interaction energies:
The full section is not included for conciseness, but the following lines define the bead-bead interactions. At the end of the canonical amino acids, we've added J and O bead interactions with each other and all other residues. 

	J	W	0	0	0
	J	F	0	0	0
	J	H	0	0	0
	J	Y	0	0	0
	J	G	0	0	0
	J	R	0	0	0
	J   K   0   0   0
	J	Q	0	0	0
	J	S	0	0	0
	J	N	0	0	0
	J	D	0	0	0
	J	E	0	0	0
	J	M	0	0	0
	J	T	0	0	0
	J	P	0	0	0
	J	C	0	0	0
	J	A	0	0	0
	J	V	0	0	0
	J	L	0	0	0
	J	I	0	0	0
	J	X	0	0	0
	J   J   0   0   0
	J	O	-1000   0   0         # the letters in the first two columns define the 'interacting pair'; the three numerical columns are short-range, long-range, and super long-range interaction strengths

and similarly for O-beads:

	O	W	0	0	0
	O	F	0	0	0
	O	H	0	0	0
	O	Y	0	0	0
	O	G	0	0	0
	O	R	0	0	0
	O   K   0   0   0
	O	Q	0	0	0
	O	S	0	0	0
	O	N	0	0	0
	O	D	0	0	0
	O	E	0	0	0
	O	M	0	0	0
	O	T	0	0	0
	O	P	0	0	0
	O	C	0	0	0
	O	A	0	0	0
	O	V	0	0	0
	O	L	0	0	0
	O	I	0	0	0
	O	X	0	0	0
	O   O   0   0   0           # all O-bead interactions EXCEPT J-O (defined above) are set to zero

The logic is as follows: 

- J : J interaction should be zero; J bead behavior to other J beads is neither attractive nor repulsive
- O : O interaction should be zero; each chain has an N-terminal O bead to affix it to the surface, but O should not interact with itself (attractive or repulsive)
- J : O short-range interaction is set to -1000; long- and super long-range interactions are set to zero because we want to prevent the mobile chains from dissociating from the surface. In addition, because the surface is comprised of J-beads that all have identical interactions strengths with the O-beads, the mobile chains will be able to translate in x-y plane (on the surface), but will remain bound.

##### Bead-solvent interaction energies:
Similar to above, J and O *only* have interaction abilities with each other, so we will also set solvent interactions to zero for these bead types.

	J   0   0  # for ETU surface
	O   0   0  # for ETU surface

### 6. Running a surface PIMMS simulation!
Per the installation instructions for the PIMMS package, execute your chains-on-a-surface simulation with the following command within the directory that contains the `KEYFILE.kf`:

	> PIMMS -k KEYFILE.kf

This should take about 45 minutes to run locally (M3 Macbook Pro), then visualize using VMD:

	> vmd traj.xtc START.pdb

The first frame of the simulation should look like you're summoning dark beings:

![Chains on surface in VMD!](/assets/vmd_surface_chains_topview.png "chains on a surface in VMD")



































