# pdb2pimms conversion report

`AF-P35637-F1-model_v6.cif` was converted into the PIMMS restart file `AF-P35637-F1-model_v6_restart.pimms` on 2026-11-18 04:32 with PIMMS 1.0.7.post11+g42ead77.d20261030.

```
pdb2pimms AF-P35637-F1-model_v6.cif --protein-mode
```

## Summary

|  |  |
|---|---|
| Input file | `/Users/alex/Dropbox/WASHU/Holehouse/tools/PIMMS/pimms/demo_keyfiles/fus/AF-P35637-F1-model_v6.cif` (mmCIF, entry AF-P35637-F1) |
| Input SHA-256 | `59be9482ef6dd0629b451a93bb1ebc79a657041ecfd9e6cc1cde32cdb3a47e67` |
| Model converted | 1 (of 1 in the file) |
| PIMMS chains | 1 (1 chain type) |
| Beads | 526 |
| Bead types | 20: `A C D E F G H I K L M N P Q R S T V W Y` |
| Box (DIMENSIONS) | 47 42 42 lattice sites (171.5 x 153.3 x 153.3 A) |
| Lattice spacing | 3.65 A per site |
| Bead position | CA (amino acids) and C4' (nucleotides) |
| Fit RMSD | 1.83 A (largest single-bead deviation 3.03 A, at A:TYR25) |
| Occupied sites | 526 of 82908 (0.63 %) |

## Bead types the parameter file must define

The restart file uses these 20 bead types: `A C D E F G H I K L M N P Q R S T V W Y`. A PIMMS parameter file is case-sensitive, so each must be defined exactly as written here. Every bead type needs a solvation line (an interaction with the solvent type `0`), an interaction line with every bead type including itself, and an `ANGLE_PENALTY` or `ANGLE_PENALTY_T_NORM` line unless the keyfile sets `ANGLES_OFF : True`. For 20 types that is 20 solvation lines, 210 pair lines and 20 angle lines.

| Bead | Kind | Residue names (count) | Beads |
|---|---|---|---|
| `A` | amino acid | ALA (13) | 13 |
| `C` | amino acid | CYS (4) | 4 |
| `D` | amino acid | ASP (26) | 26 |
| `E` | amino acid | GLU (11) | 11 |
| `F` | amino acid | PHE (13) | 13 |
| `G` | amino acid | GLY (152) | 152 |
| `H` | amino acid | HIS (3) | 3 |
| `I` | amino acid | ILE (8) | 8 |
| `K` | amino acid | LYS (14) | 14 |
| `L` | amino acid | LEU (3) | 3 |
| `M` | amino acid | MET (8) | 8 |
| `N` | amino acid | ASN (24) | 24 |
| `P` | amino acid | PRO (27) | 27 |
| `Q` | amino acid | GLN (52) | 52 |
| `R` | amino acid | ARG (37) | 37 |
| `S` | amino acid | SER (68) | 68 |
| `T` | amino acid | THR (19) | 19 |
| `V` | amino acid | VAL (5) | 5 |
| `W` | amino acid | TRP (3) | 3 |
| `Y` | amino acid | TYR (36) | 36 |

## Disorder prediction

`--protein-mode` predicted the folded domains and disordered regions of every protein chain with metapredict 3.1.0+10.ge5920ec (Emenecker, Griffith and Holehouse, Biophys. J. 120, 4312, 2021), using its V3 network with a disorder threshold of 0.5, a minimum disordered region of 12 residues, a minimum folded domain of 50 and a gap closure of 10. Each source chain was predicted as one sequence, with the runs either side of any gap joined. The freeze file leaves the disordered residues free and freezes the folded ones; nucleotides and hetero groups are frozen.

| Chain | Amino acids | Folded domains (residues) | Disordered regions (residues) | Disordered |
|---|---|---|---|---|
| A | 526 | A:THR286-ARG372 (87) | A:MET1-ASN285 (285), A:ALA373-TYR526 (154) | 439 |

## Running the keyfile

`AF-P35637-F1-model_v6_keyfile.kf` runs the converted structure as it stands. Its paths are relative to the directory it is in, so run it from there:

```
PIMMS -k AF-P35637-F1-model_v6_keyfile.kf
```

What it sets, and why:

- **Parameter file.** `gcf_rje23_v14.prm`, a copy of the residue-level amino-acid parameter file that ships with PIMMS (the default with `--protein-mode`), at `TEMPERATURE : 50`, the temperature it was parameterised at. pdb2pimms checked that it defines every bead type of the structure.
- **Box.** `HARDWALL : True` with the restart file's box (47 x 42 x 42 sites). With walls rather than periodic boundaries, no chain can meet a periodic image of itself or of another chain; the price is that the walls confine the chains (a wall counts as solvent). For more room, convert again with a larger `--padding` or `--box`, or raise `DIMENSIONS`: PIMMS centres the restart box in a larger one.
- **Freezing.** `FREEZE_FILE : AF-P35637-F1-model_v6_freeze.txt` holds the frozen beads in place (see the freeze file below), so 439 of the 526 beads move.
- **Moves.** `MOVE_CRANKSHAFT : 1`. The free beads of a partly frozen chain are moved only by the crankshaft and by the experimental segment regrowth, so no other move would do anything here. A free stretch held at both ends by frozen beads can get stuck in some arrangements with the crankshaft alone; PIMMS warns at start-up about every such stretch, and `MOVE_SEGMENT_REGROWTH` (with `EXPERIMENTAL_FEATURES : True`) frees the short ones.
- **Length.** `N_STEPS : 1000000` with `EQUILIBRATION : 100000`, `CRANKSHAFT_SUBSTEPS : 500` (about one attempt per movable bead per step, and never fewer than PIMMS's default of 500) and 1000 trajectory frames after equilibration. That is a starting point, not a converged run: check that `ENERGY.dat` has stopped drifting and that whatever you measure has settled, and lengthen the run as needed.

The keyfile in full:

```
# PIMMS keyfile written by pdb2pimms (PIMMS 1.0.7.post11+g42ead77.d20261030) on
# 2026-11-18 04:32. It runs AF-P35637-F1-model_v6_restart.pimms, converted from
# AF-P35637-F1-model_v6.cif, as it stands. Its paths are relative to the directory it is
# in, so run it from there with PIMMS -k AF-P35637-F1-model_v6_keyfile.kf. The run
# length, output frequencies and move set are a starting point;
# AF-P35637-F1-model_v6_report.md explains them.

# the system
RESTART_FILE         : AF-P35637-F1-model_v6_restart.pimms
DIMENSIONS           : 47 42 42         # the restart file's box; any box at least this big works
HARDWALL             : True             # walls, so no chain meets a periodic image
LATTICE_TO_ANGSTROMS : 3.65             # the conversion spacing, so START.pdb and traj.xtc match the input
PARAMETER_FILE       : gcf_rje23_v14.prm
TEMPERATURE          : 50               # the temperature gcf_rje23_v14.prm was parameterised at
FREEZE_FILE          : AF-P35637-F1-model_v6_freeze.txt  # metapredict's folded domains are frozen, its disordered regions move

# the run
N_STEPS              : 1000000
EQUILIBRATION        : 100000
CRANKSHAFT_SUBSTEPS  : 500              # about one attempt per movable bead (439) per step
MOVE_CRANKSHAFT      : 1                # moves the free beads of a partly frozen chain

# output
XTC_FREQ             : 900              # 1000 frames after equilibration
ANALYSIS_FREQ        : 900
EN_FREQ              : 1000
PRINT_FREQ           : 50000
ENERGY_CHECK         : 100000           # recompute the total energy from scratch this often
```

## Using the restart file

Give the restart file to a run with `RESTART_FILE`, as `AF-P35637-F1-model_v6_keyfile.kf` does. The chains come from it, so no `CHAIN` lines are needed (any are ignored).

The file was written with `HARDWALL : True`, because no chain crosses a face of the box. That makes it valid for a run with either boundary condition, and lets `DIMENSIONS` be any box at least 47 x 42 x 42: PIMMS centres this box inside a larger one. `LATTICE_TO_ANGSTROMS` only scales `START.pdb` and `traj.xtc`; setting it to the spacing used here keeps them on the scale of the input structure.

The restart file holds the configuration only. Its `ENERGY` is a placeholder (0), which no run reads: PIMMS computes the energy from the parameter file at start-up. Because it holds no step, random-number state or temperature, it starts a new run and `RESTART_CONTINUE` refuses it. The run box needs 8 bytes per site for the lattice, so this box needs 0.7 MB before anything else.

The freeze file `AF-P35637-F1-model_v6_freeze.txt` freezes the converted structure except the residues metapredict predicts disordered (see the disorder prediction above): 439 of its 526 beads are free, listed below. A chain with no free residue has a `C` line; a chain with some has a `B` line that lists its frozen beads by their index along the chain, counted from 0. Only the crankshaft move (`MOVE_CRANKSHAFT`) and the experimental segment regrowth (`MOVE_SEGMENT_REGROWTH`) move the free beads of a partly frozen chain - every other move leaves such a chain alone - so a run that uses this file needs one of them. Frozen beads still occupy their sites and still contribute to the energy, so the parameter file must still define their bead types. Mobile chains added with `EXTRA_CHAIN` are numbered from 2; they are optional here, because the free residues can move on their own. `EXTRA_CHAIN` sequences are upper-cased unless the keyfile sets `CASE_INSENSITIVE_CHAINS : False`, which an added chain with lower-case bead types (nucleotides) needs.

- chain 1 (A:MET1-TYR526): 439 of 526 beads free - beads 0-284 (A:MET1-ASN285), beads 372-525 (A:ALA373-TYR526)

## Chains

| chainID | chainType | Source chain | Residues | Beads | RMSD (A) | Rg input / lattice (A) | Placed by |
|---|---|---|---|---|---|---|---|
| 1 | 0 | A | A:MET1-TYR526 | 526 | 1.83 | 51.5 / 51.5 | dynamic programming (radius 2) |

Sequences, one character per bead:

```
>1 A:MET1-TYR526
MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNTIFVQGLGENVTIESVADYFKQIGIIKTNKKTGQPMINLYTDRETGKLKGEATVSFDDPPSAKAAIDWFDGKEFSGNPIKVSFATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY
```

## What was not converted

- Nothing: every atom of the model went into a bead.

## How the conversion was done

Each residue became one bead, placed at its CA atom (amino acids) or C4' atom (nucleotides); hetero groups were placed at their heavy-atom centroid.

Consecutive residues of a chain were treated as bonded when the atoms of the linking bond (C-N for a peptide, O3'-P for a phosphodiester) are within 2.0 A, or, where those atoms are missing, when the CA atoms are within 4.2 A or the C4' atoms within 7.5 A. Hetero groups are never bonded to anything.

A PIMMS chain is a self-avoiding walk in which consecutive beads are Moore neighbours (they differ by at most one site along every axis). Each chain was placed by dynamic programming, which finds the walk with the smallest total squared distance between each bead and its scaled input position, among walks that keep every bead within a couple of sites of its position and avoid the chains already placed. Where that walk visited a site twice, the site was forbidden to the bead further from it and the walk found again. A greedy search with backtracking is the fallback when no such walk exists. Chains were placed in the order of the table above.

Deviations are measured after the best rigid-body superposition of the lattice onto the input. Rounding each position to the nearest site of a cubic lattice gives an RMSD of about half a spacing (1.82 A here) on its own, so a fit close to that is as good as the lattice allows.

| Bonds between consecutive beads |  |
|---|---|
| Bonds | 525 |
| Along one axis (1 spacing) | 275 |
| Along a face diagonal (1.41 spacings) | 202 |
| Along a body diagonal (1.73 spacings) | 48 |
| Bond length in the input | 3.72 +/- 0.26 A (range 3.01 to 3.98 A) |
| Bond length on the lattice | 4.48 +/- 0.92 A (range 3.65 to 6.32 A) |
