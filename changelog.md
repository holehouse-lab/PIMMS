## 1.0.8 (November 2026)

Version 1.0.8 represents a major (although largely-non-breaking or feature-adding) update, primarily in terms of improving test coverage and adding better / updated docstring coverage across the base. We also addressed a number of edge-case bugs and shortcomings; these fixes largely address moves and features added over the last few months, and so touch no previously published work, but a few are long-standing analysis-output errors and are worth calling out:

* The periodic radius of gyration calculation carried a small parallel-axis bias
* Every intra-chain observable of a chain spanning more than half the box was wrong
* `INTSCAL_SQUARED.dat` held the square of a mean rather than the mean of squares
* `DISTANCE_MAP.dat` was only ever half populated.

None of these changes reflect any published work but we're noting this just for completeness! In recent months we did add full seed reproducibility (previously, driven by how random numbers were generated in the Cython/Python codebases this did not happen). In 1.0.8 we have changed how we do the random seeding to (1) remove the risk of repeat patterns in long simulations and (2) ensure full seed-dependent reproducibility. The upshot, however, is that running a simulation with 1.0.7 (or earlier) and 1.0.8 with the same seed will not give the exact same trajectories (although should sample from the same system/distribution).

1.0.8 is also a good deal faster: the per-megamove bookkeeping is compiled and long-range energy evaluations are 2.5 to 3 times faster, which roughly halves a default megamove. Every trajectory is bit-identical to before these speed-ups.

New in 1.0.8: `pdb2pimms` (start a simulation from a PDB or mmCIF structure), freezing single beads of a chain, exact run continuation (`RESTART_CONTINUE`), a record of the configuration every run actually used (`keyfile_used.kf`), an experimental segment-regrowth move for partly frozen chains, and three new demos.

Below, we outline the major changes, minor changes, and new features/updates as three separate sections.

### Major changes

These change what a run writes, the numbers in it, or the trajectory a given seed produces. If you compare 1.0.8 output against an earlier version, these are the reasons it differs.

#### Output files

- **Output files are only created when there's something to write to them.** PIMMS used to pre-create an empty stub of every file it might write, so runs left piles of empty files. A missing file now just means that analysis never wrote anything (switched off, never fired in the production window, or nothing to report); it isn't an error.
- **Start-up deletes every output the run could write**, from one manifest in `CONFIG` (`ANALYSIS_OUTPUT_NAMES`, `analysis_output_files()`), so a stale file from an older run can't pass as this run's. This replaces several ad-hoc clean-up blocks that missed some files.

#### Intra-chain observables

- **Every intra-chain observable is computed on the chain made whole.** Rg, asphericity, end-to-end distance, internal scaling, distance maps and residue-residue distances used a single image rebuilt around the circular-mean centre of mass, which tears any chain longer than half the box (a 30-mer in a 13-box had 17 of 56 `RG.dat` rows wrong), and pair distances used the minimum image, which picks the wrong image past half a box. Chains are now bond-walked into one image and measured with plain geometry. Chains shorter than half the box are unchanged.
- **The periodic Rg lost a small one-sided bias.** The gyration tensor was referenced to the circular COM rather than the mean of the rebuilt coordinates.
- **Hardwall observables use hardwall geometry.** `Chain` discarded its hardwall flag, so they used periodic geometry.

#### Internal scaling and distance maps

- **`INTSCAL_SQUARED.dat` is the mean of squared distances, as documented.** It was the square of each snapshot's mean distance, which throws away the spread within a snapshot.
- **`INTSCAL.dat` and `INTSCAL_SQUARED.dat` gain a final row,** the end-to-end separation, which the profiles stopped one gap short of. The scaling-exponent fit can now use it. Scripts that hard-code the row count need a look.
- **`DISTANCE_MAP.dat` is the full symmetric matrix.** Only its upper triangle was ever filled.
- An analysis that never sampled writes nothing and warns, rather than a profile of zeros or a `nan nan` fit.
- Worth noting that in Martin, Holehouse, Peran et al. Science 2020 all analysis was done outside PIMMS (with SOURSOP), so the scaling exponents reported there are unaffected.

#### Cluster analysis

- **Long-range clusters follow the documented rule.** Two chains are in the same cluster when a pair of their beads touches or has a nonzero LR or SLR energy. The code joined beads by distance whatever their energy, left SLR pairs out, and searched in one direction only, so with mixed bead types clusters overlapped and depended on the seed (25 memberships for 20 chains in one run), contaminating every `LR_*` file. The single-image gather uses the same rule now.
- **Percolating clusters are marked, not measured.** A cluster connected to its own periodic image has no Rg, hull or radial profile, so those columns are `nan` or `-1`, the radial row is left out, and it's logged. The test is exact (it also catches clusters that only close after winding a short axis twice) and linear in memory (it was quadratic: 16 GB for a 10^5-bead network).
- **Cluster shapes.** `extract_cluster_polymeric_properties` re-applied a periodic correction to single-image positions, inflating Rg by up to 15% for elongated clusters. The gather now shifts by whole box periods, so a cluster straddling a low face has the right centre of mass.
- **Radial density profiles.**
  - Periodic profiles bin on the minimum image; beads were being dropped.
  - Hardwall shells are normalised by the sites inside the box (a full wall film read 65% full), and empty shells are `nan` rather than 0.
  - The centre is the plain rounded mean, rounding halves up so it moves with the cluster.
  - Profiles kept after the bead-threshold filter keep their original cluster index (they were renumbered from `C1`).
- `CHAIN_<type>_CLUSTERS.dat` rows lead with the step number, like every other per-step file.

#### Sampling and detailed balance

- **`PARALLELIZE` biased `MOVE_SLITHER` and `MOVE_PULL` towards compact chains.** The parallel kernel can only move a chain that fits inside one block's interior, and PIMMS picked the parallel or serial kernel each megamove according to whether every chain fit right then. That's a one-way street: a compact chain got the kernel that can't extend it, and an extended chain got the one that can compact it. With 8-bead chains at infinite temperature, extended conformations were halved and ⟨Rg²⟩ came out 3% low, with every energy check passing. The split is now fixed at start-up by chain length and both kernels run every megamove.
  - *Who was affected:* boxes of at least 24 sites a side (40 with LR interactions) holding chains whose extent hovers around the block interior, roughly 7 to 12 beads in a 24 to 48 site box. If that's you and chain dimensions mattered, rerun.
- **Two more `PARALLELIZE` fixes, which change parallel trajectories.**
  - The block shift is drawn over the whole box rather than one block length. On an axis that doesn't divide evenly, some positions were never movable, and a crankshaft bead's sweep rate depended on its coordinate.
  - The serial and parallel passes run in coin-flip order. A fixed order isn't reversible, which a system-wide TSMMC excursion needs.
- **`PARALLELIZE` crankshaft relaxes more slowly per step than serial, but samples the same equilibrium.** This isn't a bug, but it's why energies at a fixed step differ between the two; don't compare them. Three fixes cut the penalty and made the accounting honest:
  - The frozen halo shrank from 5 to 2 sites with LR interactions (3 to 1 without). On a 44-site box that lifts the movable fraction per sweep from 16% to about 55%.
  - `CRANKSHAFT_SUBSTEPS` is shared out exactly across blocks; rounding down lost attempts.
  - Each block's RNG stream gets a scrambled seed. Seeding from `block + seed` made adjacent sweeps replay each other almost entirely.
  - The parallel kernels also report the attempts they actually make, so `ACCEPTANCE.dat` is no longer understated.
- **VMMC broke detailed balance in three ways** (all fixed):
  - Failed links to chains left behind used the wrong reverse energy, which strongly favoured contact (two monomers touched 86% of the time; the right answer is 43%).
  - Frustrated links (a failed link to a chain that joined through another one) were left out of the acceptance ratio.
  - Under a hardwall, the commit wrapped a chain that recruitment had treated as absent.

  **NB**: I still don't consider VMMC safe for production ;-), so play with it if it makes sense for your system, but probably don't use it for real science yet. It stays experimental.
- **Rigid rotations weren't invertible.** `chain_rotate` and `cluster_rotate` turned about the rounded centre of mass, which usually lands on a different site after the turn, so about a third of accepted rotations couldn't be reversed.
  - They now turn about the bead nearest the centroid.
  - Under periodic boundaries a cluster is single-imaged first, and one that winds round the box is rejected.
  - A rotation that maps a body onto itself is rejected rather than counted as accepted; it was inflating cluster-rotate acceptance about 2.5x with lots of monomers.
- **`cluster_translate` proposals were directional.** Chains were moved one at a time, so a bead landing on a cluster-mate's not-yet-vacated site clashed depending on the direction of the shift. The whole cluster is now lifted before anything is placed. This only bit when clusters were big relative to the box, which is when cluster translation is least useful anyway.
  - Under a hardwall, a box-spanning cluster could be "translated" through the wall with the wrong energy; out-of-box moves are now rejected.
  - `rigid_cluster_move` no longer counts intra-cluster LR pairs twice.
- **`chain_pivot` was off by one.** Pivot point 1 was a null move logged as accepted, 3-mers never moved, and one end of a 4-mer was unreachable.
- **Monomers walked through hard walls** in every crankshaft and slither kernel: single-bead moves were wrapped and never checked against the wall.
- **TSMMC.**
  - The system-wide excursion now runs a palindromic schedule, which tempered transitions require.
  - The linear schedule no longer overshoots `TSMMC_JUMP_TEMP` (a float `np.arange` endpoint).
  - A TSMMC draw inside an excursion is redrawn from the other moves; it had become a crankshaft megamove even with `MOVE_CRANKSHAFT : 0`.
- **The executed move mix matches the keyfile.** A move drawn for a chain it can't act on (say a rotate on a monomer) used to be silently swapped for a crankshaft megamove, so the real mix depended on composition. It's now a rejected null move, and start-up says what fraction of steps that costs. The last selection interval is also snapped to 1.0; float round-off left it a hair short and crashed about one run in 10^8 steps.
- **Move sets that can't change a conformation warn at start-up.** Examples are rigid moves only (dimers included), pull with any rigid move, head-pivot alone, or chain-pivot alone. It's not a detailed-balance bug, since chains still sample within the shapes they started in, but `RG.dat` and friends just repeat. It's a warning rather than a refusal because rigid-body assembly is a legitimate use.
- **The serial kernels accepted any uphill move with probability at least 2^-31.** A zero acceptance draw passed any positive Boltzmann factor, so a "forbidden" penalty wasn't. A second-stage test fixes it; trajectories change only when a draw lands in that one bin.
- **Cluster-only move sets warn.** `MOVE_CLUSTER_TRANSLATE` and `MOVE_CLUSTER_ROTATE` on their own can never make or break a contact. A run where no mobile chain can move at all is refused.

#### Simulation loop and output cadence

- **Output is written after the move.** A row labelled step `N` now describes the state after step `N`. Before, `ENERGY.dat`, frames and analyses were written before each move, so every row was one move behind and the final state was never written.
- **Step `EQUILIBRATION` is the last equilibration step everywhere.** Analyses and checkpoints start strictly after it, and `RESIZED_EQUILIBRATION` resizes after that step's move, so the last equilibration frame stays in `eq_traj.xtc`.
- **System-wide TSMMC excursions.** One no longer uses up two step numbers or repeats a quench or resize update, and one drawn on the last step finishes. `TOTAL_MOVES.dat` no longer counts the excursion on top of its substeps.
- **`SAVE_AT_END`** no longer appends a mislabelled final frame or writes production coordinates into `eq_traj.xtc`.
- **`AUTOCENTER`** no longer writes a straddling chain outside the box, half a site off the lattice on odd axes, or outside a hardwall box. It now also applies to `START.pdb` and frame 0, and it warns with more than one chain, where it does nothing.
- **Quench ramps.** Rung counting was wrong on binary floats (`1.0 -> 0.7` by `0.1` was four rungs), the target was reached a rung late, cooling below T = 1 was refused, `QUENCH_AS_EQUILIBRATION` sized its window from `int(dT)`, and a ramp reaching `TSMMC_JUMP_TEMP` crashed or inverted the excursion. `QUENCH.dat` now carries six significant figures and the post-move energy.
- **`SIGTERM` and Ctrl-C stop cleanly** (exit status 143 and 130). You get one line giving the step reached and the checkpoint on disk, the trajectory is closed, and `SAVE_AT_END` frames are written.
- **A failed start says so in `log.txt`** and doesn't leave a `keyfile_used.kf` describing a run that never happened.
- **Run-directory safety.**
  - An input that is also an output name is refused.
  - `pimms_running.pid` warns when two runs share a directory.
  - A stale output that can't be deleted is reported.
  - Outputs replace a symbolic link rather than writing through it.
  - Scratch files from a killed run are cleaned up.

#### Trajectory and PDB output

- **Systems of 100,000 beads or more.** Atom serials wrap rather than crash, and `CONECT` records are only written when both serials are unique; reused serials bonded the wrong atoms. The `TER` record also crashed on chains of exactly 9999 beads.
- **More than 26 chain types are labelled** `A-Z`, `a-z`, `0-9`. Every type past 26 used to be `Z`, which lemonade merged.
- **`traj.xtc` is flushed after every frame**, so a killed run keeps every completed frame; after a `SIGKILL` it used to be unreadable. Frames carry real time and step metadata.
- **`SAVE_AT_END` writes the same bytes as streaming.** The unit cell comes from `DIMENSIONS` x `LATTICE_TO_ANGSTROMS`; reading it back from `CRYST1` crashed below about 1 Å per site. A failed write leaves whole frames and can be retried.
- **PDB files follow the wwPDB spec.** Columns are exact, files are written to a temporary name and moved into place, and a box of 10,000 Å or more is refused at input. It used to crash after equilibration and leave a truncated `START.pdb`.

#### lemonade

- **`load()` uses the run's effective configuration, not the literal keyfile.** That means restart overrides of `HARDWALL` and `DIMENSIONS`, the forced hard walls of a resized equilibration, `QUENCH_END` as the temperature, and no restart `CHAIN` lines. A wall-confined run analysed as periodic came out as a percolating condensate. Keyfile chain types are applied only when they match the PDB's own partition.
- **Hardwall trajectories were analysed as periodic,** so clusters joined across walls. Clustering, gathers, radial profiles and slab analyses all honour the walls now; a wall film's surface tension was exactly double.
- **`analyze` no longer calls one-phase solutions phase separated.**
  - *Droplet geometry:* tiny inner shells dominated the radial fit, so random occupancy came out as `rho_dense = 1.0`. Small shells and radii below two sites are now ignored. The density gap must show up in single frames, and `percolation_fraction` and `spanning_fraction` are reported.
  - *Slab geometry:* a slab must hold at least 3 standard deviations more chains than a uniform solution would (`chain_excess_sigma`), with a dilute plateau wider than a chain. Before, 96 of 202 one-phase runs were called phase separated; now none are.
- **`geometry='auto'` picks slab or droplet from the shape of the largest cluster, not the box.** A droplet in an elongated box was profiled as a slab (dense density about 4x low), and a slab in a cubic box was missed. An explicit `geometry='slab'` fit of a compact droplet is refused with a reason.
- **A slab fit is refused when its dilute plateau holds a second slab.** That second slab used to be averaged into `rho_dilute`.
- **Spanning frames are skipped, with a warning.** `droplet_shape`, the radial profile and the droplet surface tension now leave out frames whose largest cluster spans the box, since it has no single image. "Spanning" is PIMMS's exact percolation test rather than reaching the box length.
- **Surface tension.**
  - The capillary fit uses the lattice dispersion `2 - 2cos(q)` (continuum `q^2` was up to 12% off), and modes are paired properly.
  - The droplet grid is sized from the median cluster.
  - Wetting films are measured from the wall face.
  - Hardwall slabs are projected on cosine modes, though they're still flagged as unreliable there.
- **`load()` hardening.**
  - It opens running or killed runs, dropping a torn last frame.
  - It warns on a periodic trajectory loaded as hardwall, a box mismatch, a lattice round-off residual, or a bond spacing that is a whole multiple of the expected one.
  - It no longer flattens a 3D trajectory given a 2D `DIMENSIONS`, and honours lower-case chains under `CASE_INSENSITIVE_CHAINS`.
  - `traj.xtc` with `eq_START.pdb` is the production trajectory; it was being forced to hard walls.
- **Smaller changes.**
  - Equal-size clusters are ordered by lowest chain index; it was search order, so ties for the largest cluster can differ from earlier output.
  - Bools are refused as frames, chains, slices and temperatures, and non-integer chain indices raise `TypeError`.
  - `Cluster.radial_density_profile()` entry `k` is the shell at distance `k + 1`, matching the `.dat` file.
  - Slab analyses no longer report shape statistics.
- **New lemonade options and fields:**
  - `PhaseSeparationResult` gains `percolation_fraction`, `spanning_fraction`, `chain_excess_sigma` and `slab_axis`;
  - `analyze(axis=)` and `slab_normal`;
  - `min_dilute_planes=` and `frame_scatter=` on `fit_slab_profile`;
  - `Cluster.spanning_axes()`.

#### Removed, renamed and reseeded

- **Removed: the obsolete neurofilament demo** (`initialized_systems.py`) and its private `Chain(fixed=True)`. `FREEZE_FILE` is the one way to hold chains still.
- **Removed: dead code,** including:
  - the compiled modules `get_randmax`, `lattice_tools` and `random_number`;
  - `cython_backup/`, which shipped stale kernels in the wheel;
  - `append_to_xtc_file`, which re-saved the whole XTC every frame;
  - assorted uncalled methods and constants.
- **Two signatures changed:** `AcceptanceCalculator.move_selector()` no longer takes a chain length, since the move drawn no longer depends on the chain, and `KeyFileParser.parse` drops its unread `verbose` flag.
- **`CRANKSHAFT_MODE`** is accepted with a warning and ignored, so old keyfiles still run.
- **Atoms are beads.** PIMMS particles are beads, so lemonade's `n_atoms` is now `n_beads` (the old name is a deprecated alias) and `atom_chainid` is `bead_chainid`. "Atom" survives only for real PDB records.
- **New seeding, so new trajectories.**
  - Kernel seeds are passed down at their full 63 bits. They were reduced modulo 2^31 - 1, so by the birthday bound a long run repeated seeds (about 23,000 repeats in 10^7 megamoves), which inflated the variance of time averages.
  - `C_RAND_MAX` is a fixed constant rather than the platform `RAND_MAX`, which is 15 bits on Windows.
  - numpy is seeded modulo 2^32.
  - A 1.0.8 run won't reproduce a 1.0.7 trajectory for the same `SEED`, but reproducibility within 1.0.8 is unaffected.
- **Regression baselines were regenerated from this tree.** The old ones came from whatever `PIMMS` was on `PATH`, in practice a pip-installed 1.0.1.

### Minor changes

Edge-case fixes, validation and hardening. None of these changes the numbers a correctly configured run produces.

#### Keyfile, parameter and freeze files

- **Restart overrides go through the keyfile checks.** `RESTART_OVERRIDE_DIMENSIONS` and `RESTART_OVERRIDE_HARDWALL` used to be applied after the checks, so they could slip a bad box past them.
- **More bad input is refused when it's read, instead of misbehaving later:**
  - chain counts of zero or below, and `0` in a sequence;
  - `nan` or `inf`;
  - empty paths (`RESTART_FILE :` silently ran a fresh simulation);
  - out-of-range resize, quench, TSMMC and lattice-spacing values;
  - energies that overflow the kernels' 32-bit storage;
  - malformed freeze-file lines, and a huge `B` range, which used to hang;
  - `ANALYSIS_MODULE` without `ANA_CUSTOM`;
  - values that would crash or exhaust memory mid-run (huge substep counts, thread counts, lattices bigger than RAM).
- **Runs that can't do anything are refused at start-up:** every chain frozen, or every enabled move needing longer chains than the system has. More beads than sites is refused before placement.
- **Ignored keywords warn.** Keywords that are read but have no effect (quench settings without `QUENCH_RUN`, TSMMC settings without a TSMMC move, and so on) get one warning each, as do `CHAIN` lines beside a `RESTART_FILE`.
- **Text encoding.** Input files are read as UTF-8 on every platform with a byte-order mark ignored, and other encodings are refused by name. Numbers and whitespace must be plain ASCII: `1_000`, full-width digits and no-break spaces were being silently misread.
- **`ANGLE_PENALTY_T_NORM` products are computed exactly in decimal,** so exact halves round to even as documented. In floating point, 0.70 x 45 was applied as 31.
- **`CASE_INSENSITIVE_CHAINS` upper-cases bead by bead** and leaves alone any bead whose upper case would change it. `'ß'.upper()` is `'SS'`, which made chains longer.
- **`write_keyfile` output round-trips,** including lists, freeze files and disabled analyses.
- **Smaller fixes.**
  - `~` is expanded in every path, and error messages name their keywords.
  - `EQUILIBRATION : 0` with `RESIZED_EQUILIBRATION` now switches the resize off rather than failing.
  - A `parameters_used.prm` reused as a parameter file no longer stacks a header per run.

#### Restart files

- **Restart files are validated on read and written atomically.** A bad pickle, out-of-box or overlapping beads, a broken chain, or a non-finite energy raises a named error. Because files are written atomically, a crash mid-checkpoint can't destroy the last good one.
- **Chains load in ascending chain-ID order.** A hand-built file in another order mixed up trajectory and analysis columns and froze the wrong beads.
- **A hardwall restart with a chain crossing a periodic face is an error.**
- **`RESTART_CONTINUE` refuses a continuation that isn't the same run.** That covers different quench settings, a different `EQUILIBRATION`, different frozen chains or beads, or a different energy function. The energy function is checked by a SHA-256 digest, because the energy check alone misses a change to bead types that aren't in contact. A directory that still holds output is refused too; the continuation used to delete the stopped segment.
- **Checkpoints are the last thing written for a step,** synced before the rename and numpy-free. One written under numpy 2 couldn't be read under numpy 1.
- **`keyfile_used.kf` writes `CASE_INSENSITIVE_CHAINS : False`.** Re-running the echo of a lower-case restart upper-cased it into a different system.

#### Engine bookkeeping

- **`ENERGY_CHECK` also cross-checks the occupancy and type grids against the chains.** A type-grid corruption threw both energies off identically. Any failure flushes the trajectory and writes `CONFIG_AT_ENERGY_FAIL` before raising.
- **Energies stay integers.** A few paths returned floats, and the hyperloop energy oracle wrapped silently past 2^31.
- **Construction and restore validate everything before changing anything.** This covers lattices and chains, and `Lattice()` refuses an axis shorter than 7 sites, below which the LR and SLR shells alias.
- **Kernel memory safety.** Latent issues in the compiled kernels are fixed (none was reachable from the engine), and every megamove kernel validates its arguments on entry.
- **An unreachable `RESIZED_EQUILIBRATION` retry is gone.** It would have silently lengthened the run; a straddling chain at the resize now raises.
- **The end-of-run time no longer drops whole days,** and `log.txt` uses one time zone.
- **Housekeeping.**
  - Positions are Python ints.
  - Analysis routines are stored by name.
  - The meaningless "reference-kernel seed" log line is gone.
  - A cluster kernel that won't import (say, one built against another numpy) is announced rather than silently dropping every cluster analysis onto the slow path.

#### Command line and packaging

- **A stale install shadowing the package is named.** `PIMMS --version` detects an old `pimms/` directory in site-packages and names it, and `--version` and `--help` no longer import the simulation stack.
- **Usage errors exit non-zero.** That covers a bare `PIMMS`, conflicting or repeated options, `--info <typo>` and `-k ""`; a bare `PIMMS` used to exit 0, so a job script missing its `-k` line looked successful. `--keyfile` is accepted.
- **`PARALLEL_THREADS : 0`** means `OMP_NUM_THREADS` if it's set, otherwise the CPUs this process may use.
- **Dependencies.** mdtraj 1.10.1 or later (1.10.0 can't read wrapped atom serials), numpy 1.23 or later, Cython 3.0 or later to build, and `python-dateutil` dropped.
- **Packaging.**
  - The sdist and wheel no longer carry Cython-generated C (16 MB of a 23 MB wheel) or annotation HTML.
  - `cythonize(force=True)` always regenerates the C, so a tree shared between machines can't carry C from another Cython.
  - Builds from a `git archive` or GitHub tarball report their real version (`.gitattributes`).
- **`build.sh` and `clean.sh`.** `build.sh` takes `uv` or `pip` and also rebuilds the lemonade kernel. `clean.sh` works from any directory and deletes only the PDB files PIMMS writes, since with `pdb2pimms` a PDB may be an input.
- **Trajectory copies use a fixed 64 KiB chunk,** so the memory `SAVE_AT_END` needs doesn't depend on the Python version; Python 3.14 raised the default.
- **Dead config removed:** Travis, LGTM and unused `setup.cfg` sections are gone, and `ruff check` is clean.

### Improvements and enhancements

#### Starting from a structure: `pdb2pimms`

- **`pdb2pimms` turns a PDB or mmCIF structure into a restart file,** so real structures can go into a simulation with mobile chains added around them by `EXTRA_CHAIN`.
  - Each residue becomes one bead. Amino acids keep their one-letter codes, nucleotides get lower-case `a c g t u`, and modified residues take their parent's code.
  - Each chain is fitted onto the lattice by dynamic programming: 1.8 Å RMSD for proteins at the default 3.65 Å spacing.
  - It also writes a freeze file that holds the structure fixed, and a Markdown report of every choice made and every bead type the parameter file needs.
  - It's documented on a new page, *Starting from a structure: pdb2pimms*, and worked through in `demos/nucleosome_from_pdb`.

#### Freezing part of a chain

- **Freeze single beads with `B <chainID> <beads>` lines** (0-based indices or ranges), for example to hold a folded domain still while its tails move.
  - Only the crankshaft moves a partly frozen chain's free beads, which keeps exact detailed balance given the frozen ones. Every whole-chain or segment move treats the chain as frozen.
  - Single-bead moves can trap a free stretch held at both ends, so start-up warns about them.
  - Checkpoints record frozen beads, and `pdb2pimms --mobile` leaves chosen residues free.
- **Experimental: `MOVE_SEGMENT_REGROWTH` (move code 15).** It regrows a window of 1 to `SEGMENT_REGROWTH_MAX_BEADS` free beads of a partly frozen chain (default 3, at most 4), drawing uniformly from every lattice path between the window's fixed neighbours. That frees up stretches the crankshaft can't untangle.
  - Stretches or tails longer than the window still aren't guaranteed to reach every arrangement in highly frozen structures; the docs explain this.
  - It needs `EXPERIMENTAL_FEATURES : True` and is off by default. With it off, every other move is drawn exactly as before.
  - **`MOVE_FREQS.dat` and `ACCEPTANCE.dat` now have 15 count columns.**

#### New demos

- **`demo_keyfiles/nuclear_pore_complex`:** a schematic human nuclear pore complex in a nuclear envelope, sized from the literature. The envelope and scaffold are frozen whole and 88 FG-Nups are anchored by one bead each. Receptors that bind FG motifs get through (23 of 24 reach the nucleus) and inert cargo mostly doesn't (4 of 24).
- **`demo_keyfiles/star_destroyer_fleet`:** a detailed Star Destroyer, 12 TIE fighters and 8 X-wings, each ship a single chain. The Star Destroyer is frozen and the fighters fly around it as rigid bodies.
- **`demos/nucleosome_from_pdb`:** a nucleosome converted with `pdb2pimms`, with its histone tails free and linker histone H1.0 and prothymosin alpha added around it.

#### Reproducibility and provenance

- **`RESTART_CONTINUE : True` resumes a run exactly.** Checkpoints now store the step, both RNG states and the temperature, so the resumed run makes exactly the moves the uninterrupted one would have, with the same rows under the same step numbers. A plain restart (still the default) starts a fresh run from the configuration, which is what you want for staging and replicas.
- **Every run writes `keyfile_used.kf`.** It holds the resolved keywords the run actually used, under a header recording the version, seed, box and restart source. The keyfile you write isn't always what ran (restart overrides, generated seeds, chain-type merging), and three of the lemonade bugs above came from having to re-derive that. It re-parses as a keyfile, and lemonade recommends loading with it.
- **`PERFORMANCE.dat` throughput counts the steps a resumed segment actually took.**

#### Performance

None of these speed-ups changes a trajectory: the regression fixtures pass unchanged.

- **Megamove bookkeeping is compiled** (`pimms.bookkeeping`). Copying bead positions between the Chain objects and the kernels' bead table in Python was most of a megamove's cost, and it capped `PARALLELIZE`. At 12,000 beads a four-thread parallel crankshaft megamove went from 8.1 to 2.3 ms, and a serial one at default substeps from 5.6 to 2.4 ms.
- **LR and SLR energy evaluations are 2.5 to 3 times faster.** The hardwall test was hoisted out of the innermost loop so the compiler can unroll it. This comes on top of about 35% from precomputing wrapped indices.
- **Slither and pull helpers.** The helpers dropped per-access bounds checks (homopolymer slither 2.9x faster), and the parallel fit report is vectorised (49 to 2 ms at 3,000 chains).
- **`MOVE_JUMP_AND_RELAX`** computes the jump's energy locally, 30 to 170 times faster.
- **Cluster analysis.**
  - The LR search is no longer quadratic in the number of chains (47.7 to 0.39 s for 4,000 chains).
  - The single-image gather picks a grid or a sorted search by size, 30 to 250 times faster where the grid was the wrong choice.
- **Whole-system energy evaluations** fold their pair list in chunks: 12 kB per bead of transient memory becomes 3 kB. That transient was what limited the largest runnable system.
- **Smaller wins:**
  - the crankshaft's draw mask with many partly frozen chains (26.5 to 3.4 ms);
  - start-up placement;
  - VMMC neighbour scans;
  - distance-map accumulators;
  - residue-pair analysis.
- **`lemonade.load()` reads only the selected frames,** in blocks: peak memory 1,155 to 164 MB.

#### Reporting

- **With `PARALLELIZE`, start-up says what you actually get:** OpenMP availability, the halo, block grid and movable fraction, and how many chains go to each slither and pull kernel. It warns when a megamove is too small to pay for the parallel kernel.
- **Start-up warnings for moves that can't do anything.** It warns about null moves (naming each move and the fraction of steps it wastes), monomers a move can't move, and move sets that can't change any conformation. Every percolating cluster is logged.
- **`REDUCED_PRINTING` silences the per-event messages.** Quench and summary messages show the values you wrote, and the run logs a completion line and any `ENERGY_CHECK` failure.
- **A `MOVE_VMMC` run reports its accepted collective moves and largest cluster.**
- **A custom analysis that calls `sys.exit()` stops the run with an error.** It used to end the run part-way through with exit status 0.
- **Exact numbers.** `ENERGY.dat` and `QUENCH.dat` write integer energies exactly, timings use the monotonic clock, and a warning flags distance maps over 1 GiB. `absolute_energies_of_angles.txt` records both the requested and the applied (rounded) angle penalties.
- **lemonade warnings.** lemonade warns on a box mismatch, a lattice round-off residual, more than 62 chain types and an empty frame selection.

#### Tooling

- **Benchmarks.** `benchmark_parallel_2d.py` reproduces the speed-up tables on the parallelization page, `end_to_end.py` works again, and the benchmarks build in a temporary copy rather than over a demo's files.
- **`run_up_test_sims.sh`** regenerates every fixture baseline. It was a copy of the clean-up script and never ran anything.

#### Test suite

- **The detailed-balance tests now have real power.** The old tolerance was about thirty times the standard error and passed a 50% error in the Metropolis exponent. Tests now use:
  - autocorrelation-aware tolerances built from the reference trace alone;
  - a second observable, ⟨Rg²⟩;
  - a positive control that runs at a scaled temperature and must fail.
- **Every move is checked against exact enumeration and transition counts.** Each move is tested against exactly enumerated ensembles and free-draw forward and reverse transition counts, and every parallel kernel in boxes that really split into blocks; several parallel fixtures had been running single-block.
- **Energies and analyses are checked against independent oracles.** Energies go against an implementation written from the model definition, and `INTSCAL`, lemonade shapes and surface tensions against plain-numpy oracles.
- **Mutation testing.** Every survivor of a 179-mutant run now has a test that kills it.
- **The regression harness tests this working tree.** It pins each file's line count and runs in temporary directories, so tests no longer write into the repo.
- **Independent audits.** Several audits ran between September and November 2026, checking every move against exact enumeration and every energy against an independent oracle. Their fixes are the ones listed above, each with a regression test.

#### Documentation

- **New `pimms/kernels.md`** covers every compiled kernel, the conventions they all share, and how to edit them safely.
- **Every function has an accurate numpy-style docstring.** Before this release, 102 had none and 34 documented arguments that didn't exist.
- **Every docs page, both READMEs and the keyword reference were checked claim by claim against the code.** The keyword reference is now generated from `CONFIG.py`. Wrong physics descriptions were fixed: slither's non-existent Hastings correction, pull's upside-down ratio, the pre-1.0.8 rotations, and VMMC's links.
- **New pages and sections.**
  - New pages: *Starting from a structure: pdb2pimms* and *Segment regrowth*.
  - New sections: freezing part of a chain, stopping and resuming a run, and irreducibility in the detailed-balance primer.
- **The output-files page** states the lazy-creation rule, the after-the-move step labels, what survives a kill, and the `nan` and `-1` sentinels.
- **The droplet surface-tension estimator's measured accuracy is documented:** 0.83 to 0.90 of the true value on lattice droplets. Use the slab estimator where you can.

## 1.0.7 (August 2026)

Small bug fixes and docs update, including a few additional tests, plus a clean-up of everything the package build log raises.

- **Bug fix (`pytest pimms/` executed a 10-million-iteration loop during collection).** The same bug class fixed in 1.0.5 (`test_megagrank.py`) had two more instances: `pimms/randint_test.py` and `pimms/randneg_test.py` match pytest's default `*_test.py` collection pattern, and their top-level bodies print 10 million / 1 million random numbers. `testpaths` protects a bare `pytest`, but an explicit `pytest pimms/` imported and ran them at collection time (confirmed: collection hangs printing random integers). Worse, all of these dev scripts ship in the wheel as importable modules whose *import* runs the loops. Renamed to `dev_randint_check.py` / `dev_randneg_check.py`, and every developer script in the package (`check_randomness.py`, `cython_testing.py`, `print_interaction_matrix.py`, `dev_megacrank_rng_check.py` included) is now `__main__`-guarded, so importing any shipped module is side-effect free. Two new package-hygiene tests pin this permanently: one imports every shipped module in a fresh interpreter and asserts nothing is printed and nothing raises, the other asserts no file outside the test packages matches a pytest collection pattern.
- **Bug fix (the NeurofilamentDemo worked example could not run at all).** Many moons ago (2015?) we developed a demo with a neurofilament system. This became dead; `initialized_systems.NeurofilamentDemo` - the documented example of programmatic system construction - had rotted against the evolving constructors: it called `Chain` without the `dimensions` / `LR_int_seq` / `LR_IDX` / `chainType` arguments and `Lattice` without the required `lattice_to_angstroms`, so instantiating it raised `TypeError` before building anything. The constructor calls are now current; the geometry is parametrized (box size, sidearm length/spacing) so the demo can be exercised at test size instead of only at the hardcoded 500^3; and the site scan walks only the filament region rather than all 125 million lattice sites. Two tests build the system small and verify grid/chain consistency. (This obsolete demo was subsequently removed in 1.0.8.)
- **Bug fix.** `ANA_RESIDUE_PAIRS` accepted negative residue indices!!! Very silly; this was parsed as valid integers and are used downstream directly as Python list indices, so `ANA_RESIDUE_PAIRS : -1 3` silently measured the distance from the *last* residue - a number the user never asked for, written to `RES_TO_RES_DIST.dat` with no hint anything was wrong. Negative indices are now rejected at keyfile-parse time.
- **Bug fix.** The parameter-file parser's error messages used `\\n` inside f-strings, so a malformed interaction line produced a message with a literal backslash-n in it instead of a line break.
- **Bug fix.** The quench-overshoot error read "suggests an in the logfile" - a garbled fragment. It now names the actual quantities (`QUENCH_STEPSIZE` vs the start-to-end temperature range) and says what to change.
- **Bug fix (`-Wsometimes-uninitialized` in the reference crankshaft kernel).** In `mega_crank.pyx`'s `get_angle_energy_change`, the bead-flag `if/elif` chain (flags 1-6) had no `else`, so a flag-0 (single-bead) entry would read `offset_start` / `offset_end` **uninitialized** and then index `idx_to_bead` with a garbage range in a `boundscheck(False)` function. The path is unreachable for valid input - single beads always carry `skip_angles = 1`, which returns earlier - but the kernel's memory safety was resting on a bookkeeping convention maintained in a different file, and the compiler rightly flagged it. A defensive `else` now yields an empty range (energy change 0), exactly matching the guard the optimized kernel (`mega_crank_fast`) already carried; the kernels remain bit-identical (verified by the equivalence suite). Relatedly, the two `# skip angles = ...` comments in `crankshaft_list_functions.py` were inverted and now match the code.
- **License metadata migrated to an SPDX expression** (`license = "LGPL-3.0-or-later"` plus `license-files`, build requirement bumped to `setuptools>=77`), replacing the deprecated TOML-table form and the deprecated License classifier. `uv build` previously emitted four `SetuptoolsDeprecationWarning`s about these, with support scheduled to end in February 2027; the build is now free of them.
- **Bug fix (the published package shipped 3.5 MB of stale build junk).** `MANIFEST.in`'s `graft pimms` sweeps the *filesystem*, not the git index, so untracked local artifacts went into the sdist AND the wheel: `.o` object files from a 2021-era Python 3.8 in-place build (`pimms/build/temp.macosx-11.0-arm64-cpython-38/`), nine Cython annotation `.html` files, six benchmark logs, `.pytest_cache`, and fifteen `.gitignore` files. The manifest now prunes `pimms/build` / `pimms/.pytest_cache` and globally excludes `*.o`, `*.log`, top-level annotation HTMLs, `.gitignore` and `.DS_Store`; the stale build directory itself is deleted. (Deliberately shipped files - the test fixtures, `pimms/data`, the tracked `_pbc.html` annotation and the checkpoint archive - are untouched.)
- Removed two blocks of superseded legacy code in `inner_loops_hardwall.pyx` that sat *after* `return` statements inside triple-quoted strings - Cython evaluated them as unreachable string expressions and both Cython and the C compiler flagged them (`-Wunreachable-code`). No behavioural change; the compile is now warning-free apart from three benign `-Wsign-compare` in Cython-generated loop counters and the deliberately-retained `crank_it_good` reference routine.
- Changelog bookkeeping: the "Supported Python and dependency versions" section shipped in 1.0.6 but was filed under the 1.0.5 header (and still contained a stale draft bullet claiming the README badge reads "3.8 – 3.14"); it now sits under its own `## 1.0.6` header with the stale bullet corrected, and the note about `docs/_static/custom.css` needing a `git add` is updated - it is committed.

## 1.0.6 (July 2026)

### Supported Python and dependency versions

**PIMMS now requires Python ≥ 3.10 and is verified on 3.10, 3.11, 3.12, 3.13 and 3.14.** Each of those was checked properly rather than assumed: a clean source tree per version, the Cython kernels compiled from scratch, and the **full 737-test suite run on each** - all passing, with no compiler errors and no new warnings.

The floor moves up from 3.8 because 3.8 and 3.9 are both past end-of-life and the scientific stack PIMMS depends on dropped them some time ago, so those installs could only ever resolve to years-old pins. There is deliberately **no upper bound** on `requires-python`: PIMMS uses no version-gated syntax or standard library API, so a new release is expected to work as soon as `numpy` / `scipy` / `mdtraj` publish wheels for it - that, rather than PIMMS, is the practical gate. (3.15 is not claimed: it is at 3.15.0a5 and `numpy` has no wheels for it, so there is nothing to test against.)

- **Bug fix (the declared `scipy` floor was unusable).** `scipy>=1.7` was not merely loose, it was wrong, for two independent reasons. `scipy.spatial.QhullError` - used by the cluster convex-hull analysis - only moved to that location in 1.8.0, so on 1.7.x **importing PIMMS fails outright**; and the 1.8.x macOS/arm64 wheels **segfault inside their own LAPACK `svd`**, which `curve_fit(..., bounds=...)` reaches on every binodal fit. The floor is now `scipy>=1.9`, the oldest release that actually runs.
- The other floors were raised to versions that can actually be selected: `numpy>=1.21` (1.20 has no Python 3.10 support, so `>=1.20` was unreachable) and `mdtraj>=1.10` (likewise). All of these are now *tested* minimums, not nominal ones - the full suite is run against exactly `numpy 1.21.3 / scipy 1.9.0 / mdtraj 1.10.2 / python-dateutil 2.8.0` on Python 3.10, as well as against current releases.

- **Bug fix.** `pimms/setup.py` - a dead, superseded build script sitting *inside* the package - imported `distutils`, which was removed from the standard library in Python 3.12. Nothing referenced it (the real build script is the one at the repo root, which uses setuptools and additionally handles the OpenMP flags and the explicit extension list), but because it lived inside the package directory it was shipped in the wheel as an importable `pimms.setup` module that could not import at all on any modern Python. Removed.
- The source was audited for every other 3.12-3.14 removal and deprecation - the other retired standard-library modules, `pkgutil.find_loader`/`get_loader` (removed in 3.14), `datetime.utcnow`, `typing.ByteString`, `imp.load_*`, `locale.getdefaultlocale`, the removed `unittest` aliases and the removed `ast` node classes. `distutils` was the only hit.
- Classifiers list 3.10 through 3.14, the README badge reads "3.10 – 3.14", and the installation guide states what is actually tested.
- A `ComplexWarning` surfaced by the 3.13/3.14 runs came from the *reference* implementation in the new equivalence tests, which deliberately reproduces the old `np.linalg.eig` behaviour: `eig` makes no symmetry assumption and can return a complex array for a gyration tensor that is symmetric only to within rounding. That is exactly the hazard the switch to `eigh` removed; the reference now discards the imaginary part explicitly so the warning does not leak.

### Documentation and Read the Docs

- **Bug fix (`.gitignore` hid a docs asset).** The byte-compiled block contained a bare `*c`. Git matches such a pattern against *every path component*, not just filenames, so it also ignored any directory ending in "c" - including `docs/_static`. The practical effect was that `docs/_static/custom.css`, which `conf.py` references via `html_css_files` and which carries the PIMMS brand colours, was never committed: the published docs have been building without it. The pattern is now `*.c`, which is what it was meant to be (its only other effect was on the Cython-generated C sources, and those are tracked anyway, so nothing else changes), and `docs/_static/custom.css` is now committed.
- **Bug fix (Read the Docs reported the version as "unknown").** Read the Docs checks out a *shallow* clone, so `git describe` finds no reachable tag and versioningit returns the `0+unknown` default. A `build.jobs.post_checkout` step now fetches the tags and un-shallows the checkout, and `versioningit` is added to `docs/requirements.txt` (the package itself is deliberately not installed for a docs build, so without it there is nothing left for `conf.py` to read the version from). Verified by simulating the whole RTD pipeline - shallow clone, post-checkout job, install from `docs/requirements.txt` only, no compiled extensions, PIMMS not installed - which now renders the correct version and builds with zero warnings.
- **Bug fix.** Two Sphinx warnings from `Simulation.rigid_cluster_move`'s docstring: the `****...****` separator lines were parsed as reStructuredText transitions, which are not legal in that position. Replaced with a bold sub-heading.
- **The docs index page now shows the version and release date**, matching how SOURSOP does it: `docs/conf.py` resolves the version through a fallback chain (versioningit read straight from the git tags, then a versioningit-written `pimms/_version.py`, then the installed `idptools-pimms` metadata, then `"unknown"`) and looks the release month up in this changelog, exposing both to reStructuredText as a `|version_info|` substitution used on `index.rst`. Reading from the git tags means the version is right even in an environment where PIMMS itself is not installed - notably Read the Docs (see the shallow-clone fix above).
- Release headers in this changelog now carry the release month (`## 1.0.0 (July 2026)`), which is what `conf.py` reads the date from. A development build reports the date of the release it sits on top of.
- `START.pdb`'s `CRYST1` record is now documented in the output-files reference as the periodic unit cell (`L * LATTICE_TO_ANGSTROMS`), so the convention is written down rather than inferred.
- The `lemonade` docs now say that clusters are ordered by **bead** count, not chain count, and why the distinction matters for a multi-component system.
- The `ANA_CLUSTER` keyword description no longer claims to be "the heaviest analysis" - after the work above it is comparable to the per-chain analyses. It now describes what its cost actually scales with (the number of chains and the size of the largest cluster).
- The installation guide's testing instructions reflect the new `testpaths` setting: a bare `pytest` picks up both test packages.

### A note for developers

`docs/generate_keywords.py` regenerates the keyword reference from `CONFIG.py`. When Sphinx runs it (via `conf.py`) it uses the checkout, but running it by hand puts `docs/` on `sys.path` first, so a bare `python docs/generate_keywords.py` imports whatever `pimms` is *installed* rather than the working tree. Set `PYTHONPATH` to the repo root if you run it standalone. (`keywords.rst` is gitignored and rebuilt on every docs build, so this only bites when checking a `CONFIG.py` edit by hand.)

Relatedly, `chain.py` reaches its dependencies with `from . import lattice_utils`, which prefers an attribute already set on the `pimms` package object over anything installed into `sys.modules`. The `test_chain.py` stubs are therefore bypassed whenever an earlier test has imported the real module, which made those tests pass in a full run but fail when the file was run on its own. The stubs are now a faithful superset of what `chain.py` uses.

## 1.0.5 (July 2026)

A round of small bug fixes and analysis-layer performance work found while auditing the codebase. Nothing here changes the physics of a simulation - the fixes are in the output files, the analysis layer, the documentation and the tooling, and the sampling itself is untouched. Every fix has a regression test; the suite grows from 640 tests to 737.

### Bug fixes: trajectory output

- **Bug fix (PDB unit cell off by one lattice site).** The `CRYST1` record in `START.pdb` was written as `(L - 1) * LATTICE_TO_ANGSTROMS`, i.e. the extent spanned by the occupied sites, rather than `L * LATTICE_TO_ANGSTROMS`, the period of the lattice (sites `L-1` and `0` are periodic neighbours one lattice unit apart). The XTC frames have always carried the correct `L * spacing` box, so the topology and the trajectory disagreed. Consequences: any PBC-aware calculation done on `START.pdb` in mdtraj/VMD was wrong by one lattice unit, and `lemonade.load(pdb=...)` without a keyfile (where the box has to be inferred from the file) inferred an `L-1` box and then wrapped the coordinates into it, silently corrupting them. For a 2D system the `c` axis is now one lattice unit, again matching the XTC.

### Bug fixes: `lemonade`

- **Bug fix ("largest cluster" meant most chains, not most beads).** `Frame.clusters` inherited its ordering from `get_cluster_distribution`, which sorts by the number of *chains* in a cluster. Everything downstream treats `clusters[0]` as the condensate — `frame.droplet`, `condensed_fraction`, `largest_cluster_size`, the radial density profile, `droplet_shape` and both surface-tension estimators. The two orderings coincide only when every chain is the same length, so in a multi-component system with unequal chain lengths lemonade was measuring the wrong cluster. Clusters are now ordered by bead count, once, in `Frame.clusters`.
- **Bug fix.** `fit_slab_profile` placed the fitted tanh centre at `length / 2`, which is only the middle of the window when the coordinate axis starts at zero. It is now anchored on `coord[0]`, so a profile the caller built on any coordinate range fits correctly. (`slab_density_profile` always starts at zero, so results from `analyze()` are unchanged.)
- `slab_density_profile` no longer accepts `min_beads`. It bins every bead in the box — which is what makes the result a density profile a coexistence fit can be run against — so the argument never did anything.

### Bug fixes: analysis and parsing

- **Bug fix.** `analysis_structures.py` raised `AnalysisStructureException` on three validation paths but never imported it, so a mismatched internal-scaling profile or distance map died with `NameError` instead of the intended, explanatory error.
- **Bug fix.** The freeze-file `B` (bead) directive split on a single space rather than on whitespace, so the empty field left by the space after the `B` made `int('')` fail: a well-formed `B` line was reported as a malformed parse error instead of reaching the not-yet-implemented guard. The bare `except:` clauses around both directive parsers are now `except ValueError:`, so they no longer swallow unrelated failures.
- **Bug fix.** `FreezeFile.validate_freeze_file`'s error message said the offending chain *was* present in the lattice — the opposite of the condition that triggers it.
- **Bug fix.** The default `RESTART_FREQ` (`"Every 10th-percentile"`) floored to `0` for runs of fewer than 10 steps, and the keyfile sanity check then rejected the keyfile with "Expected RESTART_FREQ to be larger than 0" — about a keyword the user had never set. It is now floored at 1, so short smoke-test runs work out of the box.
- **Bug fix.** `evaluate_performance` divided by a step rate of exactly zero when called on step 0.
- **Bug fix.** The "unable to find an empty lattice site" warning in `get_empty_site` printed a literal `%i` — the attempt count was never substituted in.

### Bug fixes: tooling and safety

- **Bug fix.** `pimms/test_megagrank.py` was a hand-run developer script, not a test module, but its name meant pytest imported it during collection and ran its top-level code — which called a `mega_crank.python_randint` entry point that has never existed, so a bare `pytest` at the repo root died during collection before running a single test. Renamed to `pimms/dev_megacrank_rng_check.py`, pointed at the real `randint_python` entry point, and `testpaths` is now set in `pyproject.toml` so pytest only walks the actual test packages.
- `cluster_kernels.snakesearch_single_image` indexes a flat occupancy grid by raw position with bounds checking disabled. It now validates that every input coordinate lies inside the box, turning what was silent out-of-bounds memory access into a clear `ValueError`. (The pure-Python fallback uses a dict and tolerates arbitrary coordinates, which made the divergence easy to miss.)
- **Bug fix.** Two trajectory-writing helpers in `lattice_utils` (`append_to_xtc_file_non_redundant` and `update_master_traj`, the latter on the live `SAVE_AT_END` path) caught every exception with a bare `except:` and then called `exit(1)`. That swallowed `KeyboardInterrupt`/`SystemExit`, and `exit` is the `site` builtin — absent under `python -S` — so the failure path could itself fail, and in any case it tore down the caller's process rather than letting them handle the error. Both now raise `LatticeUtilsException`. `update_master_traj` also tested `master_traj == None` rather than `is None`, which invokes mdtraj's `__eq__` instead of the identity check that is meant.

### Performance: analysis routines

Profiling a representative run (400 steps, 60 x 20-mers in a 16^3 box, `ANALYSIS_FREQ : 25`) showed **95% of the wall time going to the analysis routines** and only 0.17 s to the Monte Carlo engine itself. That run now takes **0.73 s instead of 4.89 s (6.7x)**, with the analysis block down from 4.67 s to 0.51 s (9.2x). The gain grows with chain length, because the analyses that dominated were quadratic in it: a 200-step run of 120 x 40-mers in a 24^3 box goes from **16.7 s to 1.2 s (13.7x)**.

Every change below is either bit-identical to what it replaces or, where noted, agrees far beyond the precision of the output files. Two independent checks back that up: the 15 simulation regression baselines are unchanged, and running the same simulation under the old and new code produces **byte-identical output for 29 of 31 files** - every analysis `.dat`, and bit-identical XTC coordinates. The two that differ are `log.txt` / `PERFORMANCE.dat` (timestamps and timings) and `START.pdb`, whose single differing line is the intentional `CRYST1` bug fix above.

- **`get_inter_position_distance` was 56% of the entire run** (396,000 calls). To measure the distance between two 2- or 3-element lattice positions it built two numpy arrays and dispatched `np.power` / `np.sqrt` per call, so the numpy scalar machinery cost roughly 16x the arithmetic: **6.51 us/call, now 0.42 us** in plain Python. Bit-identical (integer coordinates stay exact through the squaring; `math.sqrt` and `np.sqrt` are both the correctly-rounded IEEE-754 double root), verified over 51,000 random cases including float inputs.
- **The distance-map and internal-scaling analyses were O(L^2) Python double loops** over that function - together 2.7 s of the 4.9 s. New vectorized `lattice_analysis_utils.get_distance_matrix` and `get_internal_scaling_profile` do one numpy pass per row block / per sequence separation instead. Bit-identical. `get_distance_matrix` works in row blocks so the intermediate does not scale as `L^2 * n_dim` for long chains.
- **`get_eigenvalues_of_the_T_matrix` looped in Python over every bead**, calling `pbc_correct` and allocating an `np.outer` per bead. Now a single vectorized pass: **7.17 ms -> 0.74 ms for a 2000-bead cluster (10x)**. It also now uses `eigh` rather than the general `eig` - the gyration tensor is symmetric by construction, and `eig` can return a complex array for a matrix that is only symmetric to within rounding. Every downstream quantity is a symmetric function of the eigenvalues, so the ordering difference is immaterial; results agree to ~1e-14, against output written at 4 decimal places.
- **`ANA_POL` computed every chain's properties twice** - once minimum-image, once single-image - purely to emit a finite-size warning. For a chain that does not straddle a periodic boundary the two inputs are identical by construction and the warning cannot fire, so that case now computes once.
- **The connected-component search made over 1.4 million Python-level `get_gridvalue` calls** per analysis step, two per envelope pair. Those are now a single fancy-index into the grid per round. `get_cluster_distribution` also materialised the whole unfound-chain set (`list(...)[0]`) on every round, which made it O(n_chains^2); it uses `next(iter(...))` now, which selects the same chain.
- **Envelope-pair deduplication is a full lexicographic sort**, and once the surrounding Python loops were gone it became the single largest remaining cost. `build_envelope_pairs` / `build_all_envelope_pairs` take a new `deduplicate` argument, and the two connected-component searches - which feed the pairs straight into a set, where duplicates are harmless - pass `False`. Where deduplication *is* needed (any energy evaluation, which would otherwise double count a repeated pair) rows are now packed into a single int64 key and sorted on that, with the previous `np.void` view kept as a fallback for coordinates too large to pack.
- **`find_nearest_position`** was a Python loop over every bead; for a large condensate, choosing the snakesearch seed cost about ten times as much as the compiled BFS it feeds. Now vectorized, and `cluster_utils` hands it the position array it has already built. `np.argmin` reproduces the old strict-`<` tie-breaking exactly.
- **`Chain.get_LR_binary_array` was O(L^2) per call** (`if i in self.LR_IDX` against a list) and is called once per chain in every full energy evaluation and on every single-chain move. `LR_IDX` is fixed for the life of a chain, so the array is built once in the constructor and returned read-only.

### Performance: `lemonade`

- **`phase_separation.analyze` ran the connected-component decomposition five times per frame.** `traj[f]` mints a fresh `Frame`, so the per-Frame cluster cache was discarded between each of the five passes (condensed fraction, cluster count, largest cluster, density profile, droplet shape). Membership is now memoised on the `TrajectoryStore`: **5 decompositions per frame -> 1, and `analyze()` is 3.5x faster** (0.33 s -> 0.093 s on a 12-frame test). Only the cheap membership lists are cached - the per-cluster geometry stays on the transient `Cluster` objects so it can still be collected.
- **`_analysis.gyration_eigenvalues` materialised an `(n_frames, n_atoms, k, k)` array** before reducing it - 0.72 GB for a 1000-frame, 10k-bead trajectory and 5.8 GB at 2000 frames / 40k beads. It now accumulates the tensor one component at a time (and only the `k(k+1)/2` unique ones), holding a single `(n_frames, n_atoms)` scratch array. Bit-identical.
- **`phase_separation._shell_site_counts` listed every lattice site explicitly** as float64, so its memory scaled with box volume regardless of how many beads were being analysed. Now built from separable per-axis offsets, accumulated a slab at a time: **124 MB -> 0.5 MB on a 120^3 box**, with identical counts.

### Performance: trajectory writing

- **`SAVE_AT_END : True` was O(frames^2).** Each saved frame was added with `master_traj.join(frame)`, and `join` returns a *new* trajectory holding a copy of everything so far, so writing n frames copied `1 + 2 + ... + n` frames' worth of coordinates. Frames are now buffered in a `TrajectoryAccumulator` and joined once when the trajectory is written. (The default incremental path was moved to a persistent XTC writer in 1.0.0 for the same reason; this path had kept the quadratic behaviour.) Frame count, ordering, time stamps and unit cell are unchanged.

## 1.0.0 (July 2026)

The first stable release of PIMMS. This is a large release: the Monte Carlo engine is rewritten around fast compiled kernels and a multi-threaded sampler, several new and corrected moves are added, a new `lemonade` analysis package lands, non-cubic boxes become a first-class configuration, and the whole codebase gains a complete documentation site and an extensive correctness / detailed-balance test suite. The changes are grouped by area below.

### Compiled kernels and the fast engine

- New `pimms/mega_crank_fast.pyx`: allocation-free serial crankshaft kernels (`mega_crank` / `mega_crank_2D`) that are bit-exact drop-in replacements for the reference kernels and are wired into production. The crankshaft is the hot loop of a PIMMS run, so this speeds up essentially every simulation.
- **PRNG:** the serial kernels now draw from **splitmix64** (period 2^64, passes BigCrush, identical on every platform) instead of the platform's libc `rand()`. On macOS that was the Park–Miller MINSTD LCG — a short-period (~2.1e9), lattice-structured generator that is a poor choice for Monte Carlo. The fast and reference kernels are kept bit-identical to each other so their equivalence test still holds; the 15 simulation regression baselines were regenerated under the new generator.
- **Keyfile parser robustness:** duplicate-keyword detection for every keyword, plus typed sanity checks at startup so malformed int/float/bool values fail immediately with a clear message rather than obscurely downstream.

### Parallelization (`PARALLELIZE` / `PARALLEL_THREADS`)

- New OpenMP checkerboard kernels run the **crankshaft, slither and pull** moves multi-threaded, in both **2D and 3D**. The box is split into blocks separated by a frozen halo (width `W = R_int + 2`) so no two blocks' move footprints can touch the same site; blocks run concurrently with private splitmix64 streams and per-block integer energy deltas. The decomposition depends only on box geometry, so the result is **independent of the thread count** and targets the same Boltzmann distribution as the serial sampler.
- Whole-chain moves (slither, pull) use a chain-level decomposition: a chain parallelizes only if all of its beads fit inside one block's interior; chains spanning a boundary are frozen for that sweep.
- **Frozen chains compose with parallelization** — the parallel kernels take a per-bead frozen mask, so `FREEZE_FILE` and `PARALLELIZE` can be used together (frozen beads stay as fixed, energy-contributing obstacles but are never selected to move).
- `PARALLEL_THREADS : 0` (the default) uses all cores. Measured near-linear scaling on large dilute crankshaft/slither/pull-dominated boxes (~8× on 8 threads); little benefit for small boxes (halo-dominated) or a single concentrated droplet in a large box (load imbalance). Benchmark tables are in the docs.

### Cluster analysis performance

#### `pimms/cluster_kernels.pyx` (new compiled kernel)

- Added a Cython kernel, `snakesearch_single_image`, that reimplements the single-image ("snakesearch") reconstruction of a cluster in typed C.
	- Neighbour discovery now uses a flat occupancy-index grid (a `prod(dimensions)` array mapping each PBC position to its bead index) for O(1) lookups, replacing the Python dict of coordinate tuples and the millions of per-neighbour tuple constructions / dict lookups / per-dimension Python loops that dominated the previous BFS.
	- The seed (bead nearest the PBC-aware centre of mass) is still chosen in Python exactly as before, so the kernel is a byte-for-byte drop-in for the pure-Python routine (verified on real condensate clusters and on random self-avoiding walks, 2D/3D, both interaction thresholds).
	- Wired into `cluster_utils.convert_positions_to_single_image_snakesearch` as a fast path that transparently falls back to the pure-Python implementation if the extension has not been built.
	- End-to-end this takes cluster analysis from ~0.80 s to ~0.16 s per call on a 500-chain condensate (~5×), of which the kernel accounts for ~3×.
	- Registered as `pimms.cluster_kernels` in `setup.py`; a rebuild (`build.sh`) is required to pick it up.

#### `pimms/lattice_analysis_utils.py` — radial density profile

- Rewrote `compute_cluster_radial_density_profile` from an O(offset_max^n_dim) concentric-shell site scan to an O(num_beads) histogram: each bead's Chebyshev distance from the cluster COM is binned with a single `np.bincount`, and the density at shell k is (beads at distance k) / (lattice sites in shell k). This removes the per-site membership-test / ring-scan machinery entirely.
- **Bug fix:** the previous ring-scan loop had an off-by-one (`while offset <= offset_max` with a top-of-body increment) that emitted one extra shell at `offset_max + 1`, whose extent (2k+1 = min(box)+1) spills outside the box. Profiles are now correctly capped at `offset_max` shells. The two `*_RADIAL_DENSITY_PROFILE` regression fixtures were regenerated for the box-spanning clusters; every other row is unchanged.

#### `pimms/lattice_utils.py`

- Vectorised `center_of_mass_from_positions`: the circular (PBC-aware) mean was a per-bead Python loop calling scalar `np.cos`/`np.sin`; it is now a single vectorised computation over all beads. The integer COM output is bit-identical (checked over thousands of random cases).
- Micro-optimised `get_gridvalue`: it called `get_dimensions` (an array `.shape` lookup) on every one of its ~1M calls in the connected-component search; it now branches on `grid.ndim` and uses a single tuple index.

### Trajectory writing

- **Performance:** trajectory writing was O(frames²) — every frame reloaded and rewrote the entire XTC. Replaced the reload-append pattern with a persistent XTC writer handle (`open_xtc_writer` / `write_xtc_frame` / `close_xtc_writer`) so each frame is an O(1) append. A ~450-frame run that previously slowed to a crawl now adds ~0.6 s total.
- **New keyword `TRAJECTORY_PBC_UNWRAP`** (default `False`): when enabled, each chain is made whole across periodic boundaries before it is written, so molecules that cross a box face appear contiguous in the trajectory rather than being torn in two.
- **Bug fix (one-sided unwrap):** the single-image routine used for output shifted every boundary-crossing chain to be non-negative, which translated all such chains toward one face of the box (chains only ever appeared to bulge out of one side). Added `make_chain_whole`, which anchors the first bead in place and never applies that shift, so chains now spill symmetrically out of whichever face they actually cross.

### Non-cubic (unequal-dimension) boxes

- PIMMS now supports boxes with `x != y != z` (3D) or `x != y` (2D) under periodic boundaries as a documented, tested configuration (previously blocked behind a hardwall-off guard in the keyfile parser). An exhaustive `ENERGY_CHECK` sweep confirmed the kernels are per-axis correct for every move except cluster rotation.
- **Bug fix:** `compute_cluster_radial_density_profile` sized its shells from `max(dimensions)`, so on a non-cubic box the shells ran out to the longest axis and wrapped the short axes; it now uses `min(dimensions)` — the largest shell that fits inside every axis. Cubic boxes are unchanged.
- Cluster *rotation* is genuinely incompatible with non-cubic PBC (a cardinal 90°/270° rotation is only a symmetry of a cube), so `MOVE_CLUSTER_ROTATE` on a non-cubic periodic box now raises a clear `KeyFileException` explaining the issue rather than silently drifting the tracked energy.

### Monte Carlo moves

- **Slither (`MOVE_SLITHER`)** is now an optimized whole-system reptation megamove (2D + 3D): O(1) interaction energy for homopolymers, O(N) for heteropolymers, and single-bead chains reduce to a local translation. Detailed-balanced. `SLITHER_SUBSTEPS` sets the number of reptations applied to each chain per megamove.
- **Pull (`MOVE_PULL`, code 11)** is a new cooperative-reptation megamove for rearranging dense systems, replacing the defunct ratchet-pivot. An interior bead is displaced and the following beads are pulled along to restore connectivity (the termini are not moved), letting chains rearrange where rigid moves would clash. Detailed balance via a Metropolis–Hastings proposal-multiplicity correction; new `PULL_SUBSTEPS` keyword; requires chains of length ≥ 3.
- **VMMC (`MOVE_VMMC`, code 14)** is a new Virtual-Move Monte Carlo collective move (Whitelam & Geissler, J. Chem. Phys. 127, 154101, 2007) - hat tip to [Eric Deeds](https://deedslab.ibp.ucla.edu/) for this suggestion at BPS in 2016 (!): a seed chain is given a trial translation, neighbouring chains are recruited into a moving cluster by interaction-energy gradients, and the whole cluster translates together — escaping the kinetic traps single-chain moves hit in condensed phases. Metropolis–Hastings detailed balance with a symmetric `1/n_c` cluster-size cutoff and an exact from-scratch ΔE. New `VMMC_MAX_DISPLACEMENT` / `VMMC_MAX_CLUSTER` keywords. Remains experimental.
- **Jump-and-relax (`MOVE_JUMP_AND_RELAX`, code 13).** Previously experimental move that is now decomposed into three π-preserving sub-steps (relax → Metropolis-accepted jump → relax) whose composition is detailed-balanced. The 8 regression baselines using this move were regenerated.
- **TSMMC** temperature-excursion moves (`MOVE_CTSMMC` / `MOVE_MULTICHAIN_TSMMC` / `MOVE_SYSTEM_TSMMC`) use correct **tempered-transitions work accumulation** for detailed balance; the old `systemTSMMC` module was removed.
- The self-contained moves (slither, pull, VMMC, jump-and-relax, TSMMC) all live in `moves.py` (`MoveObject`) behind a uniform parameters-in / mutated-lattice-out interface.
- **The experimental gate is now VMMC-only.** Pull, the TSMMC family, jump-and-relax, the cluster moves, non-cubic boxes, and the `EXTRA_CHAIN` / `FREEZE_FILE` / `EQUILIBRATION_OFFSET` keywords are all first-class and no longer require `EXPERIMENTAL_FEATURES : True`. Only `MOVE_VMMC` (and its `VMMC_*` tuning keywords) remains gated.

### Analysis: the `lemonade` package (new)

- New `pimms.lemonade` package for post-hoc analysis of PIMMS trajectories, loaded from an XTC + PDB + keyfile. It presents a hierarchical, index-only object model — `LatticeTrajectory` → `Frame` → `Polymer` / `Cluster` — over a structure-of-arrays backing store, with compiled kernels for the expensive steps (batched PBC unwrap, grid painting), so loading and analysis stay fast on large trajectories. Lemonade was previously implemented as a separate analysis package but has been brought into PIMMS in 1.0.0..
- Vectorised whole-trajectory conformational analysis — radius of gyration, centre of mass, asphericity, end-to-end distance, distance maps and internal scaling — computed for every chain in every frame in a handful of array operations, with the same quantities available as scalars at the single-chain level.
- Phase-separation and droplet physics: condensed fraction and cluster-size order parameters, coexistence (binodal) densities from `tanh` fits to radial (droplet) or slab density profiles, droplet shape, and interfacial tension estimated from capillary-wave (slab) / spherical-harmonic (droplet) fluctuation spectra, plus a one-call `analyze()` that auto-detects slab vs. droplet geometry. Densities are reported as occupied-site fractions in `[0, 1]`.

### Robustness and correctness

- **Quench bug fix.** Heating quenches (`QUENCH_START < QUENCH_END`) actually cooled: `QUENCH_STEPSIZE` was sign-flipped twice — once by the keyfile parser and again at runtime — so the two negations cancelled and the temperature always decreased. Removed the redundant runtime negation, so heating now heats; cooling is unchanged (all cooling regression sims stay byte-identical). Added heating and cooling regression tests.
- Assorted correctness fixes surfaced during review: angle-penalty rounding, z = 0 cluster-plane handling, the radial-density off-by-one, and the one-sided trajectory unwrap (each detailed in its own section).

### Status reporting

- The startup memory report now shows true data-buffer sizes (via `.nbytes`, not the misleading `sys.getsizeof` for numpy arrays) plus the actual process resident memory. `PERFORMANCE.dat` and the status log report both the outer master-loop steps/second and the overall MC accept/reject throughput across all sub-loops.

### New demos (`demo_keyfiles/`)

- `amphiphile_bilayer` — a lipid-style bilayer self-assembled from 5-bead HHTTT amphiphiles (two hydrophilic heads, three hydrophobic tails) in a z-elongated periodic box. Seeded and relaxed (as in the onion demo); tail–tail cohesion, head solvation and head/tail demixing hold the tails-in / heads-out membrane together.
- `slab_phase_separation` — a slab-geometry phase-separated condensate of a single sticky homopolymer, grown from a small equilibration box via `RESIZED_EQUILIBRATION` into an elongated production box.
- `multiphase_core_shell` — a four-layer core/shell ("onion") condensate seeded from a restart, demonstrating multiphase co-assembly.
- `star_destroyer` — a showcase of restart-loaded frozen chains composing with `EXTRA_CHAIN` and `PARALLELIZE`: a ~170-chain Star Destroyer hull is loaded from a restart and frozen in place while 150 small mobile "spaceships" (added via `EXTRA_CHAIN`) fly around it under pure excluded volume.

### Tests

- Added `pimms/tests/test_pbc.py`: periodic-boundary edge cases for chains large enough to straddle the box across multiple faces / wrap it more than once. It establishes that single-image reconstruction recovers such chains exactly (any number of images, one axis or several) and round-trips bit-for-bit; that minimum-image Rg collapses for chains larger than ~half the box while the single-image value stays exact and the finite-size warning fires; that impossible bonds and disconnected clusters fail with bounded, clear exceptions rather than hanging; and that every move keeps the energy self-consistent under heavy straddling (`ENERGY_CHECK`).
- Added `pimms/tests/non_equal_box/`: per-axis PBC correctness via `ENERGY_CHECK`, axis-permutation ensemble equivalence, and a cubic-vs-non-cubic dilute-bulk check, all with chains that straddle the short axis.
- Added kernel-vs-Python equivalence tests for the snakesearch kernel and a unit test locking in the radial-density off-by-one fix.
- Added a comprehensive kernel correctness + detailed-balance suite (`pimms/tests/test_kernel_correctness.py`, `test_detailed_balance.py`, `kernel_test_utils.py`) covering the fast serial, parallel, slither, pull and TSMMC paths across 2D/3D × SR/LR/SLR × hardwall, plus per-move detailed-balance gates (a move mixed with crankshaft must reach the crankshaft-alone equilibrium).
- Added `pimms/tests/test_move_and_box_gating.py` (which moves and box shapes are / are not gated) and a `pimms/lemonade/tests/` suite for the analysis package (loading, conformational analysis, phase separation, surface tension).

### Documentation

- A complete Sphinx documentation site replaces the old stubs: a User Guide (installation, an overview of the lattice / energy model and the move set, restart files, output files), an auto-generated keyword reference (built from `CONFIG` so it never drifts from `PIMMS --info`), a per-move **Moves** section with detailed-balance derivations, an **Advanced Features** section (quench, freeze files, parallelization, TSMMC, custom analysis, reference controls), the **Analysis (lemonade)** section, and a developer API reference. The custom-analysis path was also hardened to validate user code and fail gracefully.
- `PIMMS --info` now documents **every** keyword (type + description) and groups them under logical subheadings; `--info <keyword>` and `--info ALL` are supported.
- Comprehensive NumPy-style docstrings were added across the package (~225 functions).
- A full accuracy pass reconciled every doc page and keyword description against actual code behaviour — defaults and requirements, the restart-override semantics (`RESTART_OVERRIDE_DIMENSIONS` / `RESTART_OVERRIDE_HARDWALL`), the parameter-file rules (integer energies, the complete / non-redundant short-range matrix, `ANGLE_PENALTY_T_NORM`), the crankshaft sub-move accounting, and the catalogue of output files.

### Housekeeping

- Removed 8 dead functions/classes surfaced by a repo-wide usage audit.
- Git-ignore the files a simulation writes (`*.dat`, `*.xtc`, `*.pdb`, `restart.pimms`, `log.txt`, `parameters_used.prm`, `absolute_energies_of_angles.txt`) so run outputs can't be committed by accident; `clean.sh` sweeps those outputs from the repo root and every demo / validation directory.
- Copyright bumped to 2026.

## 2026-03-14

* Fixed a bug where if `SAVE_AT_END` and `RESIZED_EQUIBRIUM` were set writing a PDB/XTC file failed
* Added tests to address this

## March 2026 0.1.40 update

In preparation for the PIMMS paper, we have conducted a large-scale modernization of the PIMMS codebase. The major 

### Comprehensive unit and simulation tests

* Added 3419 lines of unit tests across many (though not all) modules
* Added large-scale simulation tests across 14 different scenarios that evaluate a large number of outputs (both system state, energies, analysis etc)
* Added performance benchmarking code to assess if changes impact performance.

### Cython Performance Refactors

#### `pimms/inner_loops.pyx`

- Refactored SR/LR pair extractors in 2D and 3D to eliminate per-call `np.delete` cleanup.
	- Self-pair is now skipped inline during the loop instead of allocated and deleted afterward.
	- Corrected PBC coordinates are cached once per neighbor instead of recomputed in each branch.
	- SR arrays are preallocated to their exact maximum size and returned as slices, avoiding post-hoc filtering.

#### `pimms/inner_loops_hardwall.pyx`

- Applied the same extract-inline / skip-invalid pattern to all six hardwall extractor functions (mixed SR/LR, LR-only, SR-only × 2D/3D).
	- Hardwall sentinel (`-1`) checks are now performed inside the neighbor loop, eliminating post-loop `np.delete` and boolean-mask filtering on hot paths.
	- SR output arrays are preallocated to exact capacity and sliced on return.

#### `pimms/hyperloop.pyx`

- Replaced Python-list `append` with typed preallocated arrays and index-sliced returns in `get_unique_interface_pairs_3D` and `get_unique_interface_pairs_2D`.
- Fixed `get_unique_interface_pairs_2D` to call the 2D SR extractor (`extract_SR_pairs_from_position_2D`) instead of the 3D variant.
- Replaced `get_gridvalue_3D`/`get_gridvalue_2D` helper calls with direct `lattice[...]` indexing inside all energy evaluation loops, removing per-iteration Python function-call overhead.
- Refactored `evaluate_angle_energy_3D` and `evaluate_angle_energy_2D` to use scalar `cdef int` locals with inlined PBC-clamp logic instead of temporary NumPy vectors.
- Removed the now-unused `fix_angle_pbc_issues` helper function.

#### `pimms/mega_crank_2D.pyx`

- Eliminated per-move NumPy array allocations in the Monte Carlo loop.
	- `single_bead_crank_2D` and `crank_it_2D` now write into a caller-provided buffer and return `int` (success/fail) instead of allocating a fresh `np.zeros([2])` on every call.
	- A single `new_position` buffer is preallocated once before the main loop and reused across all iterations.
- Cached per-bead fields (`bead_flag`, `old_x`, `old_y`, `lr_vs_sr`, `bead_id`) into C locals at the top of each MC step, reducing repeated `idx_to_bead[bead_index, ...]` indexing.
- Scalarized angle-energy vector math in `get_angle_energy_change_2D`.
	- Removed temporary NumPy arrays `a` and `b`; replaced with `cdef int` scalars `a0`, `a1`, `b0`, `b1`.
	- Cached `bead_flag` and added safe defaults for `offset_start`/`offset_end` to eliminate uninitialized-variable warnings.
	- Replaced legacy `xrange` calls with `range`.

### Cluster Utilities Performance (`pimms/cluster_utils.py`)

- Rewrote `convert_positions_to_single_image_snakesearch` with a BFS-based algorithm.
	- Old algorithm was O(N²): recomputed COM distances over all unsearched beads every iteration, used pure-Python `find_local` with nested loops, and performed O(N) `list.remove()` calls.
	- New algorithm is O(N·K) where K = (2·threshold+1)^n_dim: uses a hash-map (`pbc_to_idx`) for O(1) neighbor lookup, a precomputed offset grid, and a boolean visited array instead of list removals.
	- Eliminates `copy.deepcopy` of the positions list and per-iteration `get_inter_position_distances` allocations.
	- Preserves identical seed selection (PBC-aware COM → nearest bead) and output contract (all coordinates ≥ 0, same periodic image).

### Lattice Core Fixes (`pimms/lattice.py`)

- Fixed silent exception paths in lattice initialization routines.
	- `__fully_defined_initialization(...)` now raises `LatticeInitializationException` when lattice/type-grid dimensions do not match expected dimensions.
	- `__initialization_from_restart(...)` now raises `RestartException` for dimension mismatches, restart-overflow dimensions, and duplicate chain IDs while rebuilding chain tables.
- Fixed fully-defined initialization debug path correctness.
	- Replaced undefined `DEBUG` symbol usage with `CONFIG.DEBUG`.
	- Corrected debug iteration from dictionary keys to chain objects (`chainsDict.values()`), and ensured sanity-check failures actually raise.
- Fixed dimension equality semantics in `__fully_defined_initialization(...)`.
	- Dimension checks now compare normalized tuples, avoiding false mismatches from list-vs-tuple representation differences.
- Hardened `get_random_chain(...)` behavior and reduced per-call overhead.
	- Replaced mutable default argument (`[]`) with `None`.
	- Added explicit failure for empty/all-frozen chain-selection cases with clear `LatticeInitializationException` messages.
	- Reduced repeated membership-check overhead by using a set for frozen-chain filtering.
- Removed unused `scipy.misc` import from `lattice.py`.

### Simulation Core Fixes (`pimms/simulation.py`)

- Fixed invalid move-selection error path in `run_simulation(...)`.
	- Unknown move selector outputs now raise `SimulationException` directly with a clear message.
	- This prevents accidental `NameError` failures on the defensive invalid-selection branch.

- Added explicit all-frozen guard in the main simulation loop.
	- When every chain is frozen, `run_simulation(...)` now skips move proposal for that step instead of attempting random chain selection.
	- This avoids avoidable failures and unnecessary per-step proposal overhead in fully-frozen workflows.

- Hardened and corrected quench update behavior in `quench_update(...)`.
	- Added explicit validation that `QUENCH_FREQ` is positive, raising `SimulationException` for invalid values.
	- Fixed heating-quenches by passing a sign-corrected quench step to the temperature updater so temperature increases toward target during heating runs.

### Packaging and Versioning Migration (`versioneer` -> `versioningit`)

- Migrated package version management from `versioneer` to `versioningit`.
- Added `pyproject.toml` with:
	- PEP 517 build backend (`setuptools.build_meta`)
	- Build requirements including `versioningit>=2`
	- `tool.versioningit` configuration with fallback version `0+unknown`
- Updated `setup.py`:
	- Removed `versioneer` integration (`get_version()`, `get_cmdclass()`)
	- Added `versioningit.get_version()` for package version resolution
- Updated `pimms/__init__.py`:
	- Removed import of versioneer-generated `_version.py`
	- Now uses `importlib.metadata.version("pimms")` with source-tree fallback to `0+unknown`
	- `__git_revision__` is now set to `"unknown"`
- Removed legacy `versioneer` artifacts:
	- deleted `versioneer.py`
	- deleted `pimms/_version.py`
- Updated distribution/config files:
	- removed `versioneer.py` from `MANIFEST.in`
	- removed `[versioneer]` block from `setup.cfg`

### Chain Module Fixes (`pimms/chain.py`)

- Added validation for `LR_IDX` during `Chain` initialization.
	- Invalid long-range residue indices (negative or >= sequence length) now raise `ChainInitializationException` immediately with a clear message.
	- This prevents deferred, less-informative `IndexError` failures in downstream accessors like `get_LR_positions()`.

- Fixed `get_center_of_mass(on_lattice=...)` to correctly forward the `on_lattice` argument.
	- `Chain.get_center_of_mass()` now calls `lattice_utils.center_of_mass_from_positions(..., on_lattice=on_lattice)`.
	- This restores the documented behavior difference between lattice and continuous COM outputs.

- Improved `set_ordered_positions()` error reporting.
	- `ChainAugmentFailure` now reports the user-provided length and required sequence length correctly.
	- This makes debugging incorrect position-array assignments much easier.

### Move Engine Fixes (`pimms/moves.py`)

- Fixed `multichain_based_TSMMC()` chain-selection bound type.
	- `max_number_selectable` is now cast to a Python `int` before `random.randint(...)`.
	- This resolves crashes on Python 3.12 where `randint()` rejects `numpy.float64` bounds.
- Added a guard for the all-frozen/no-selectable-chains case.
	- If there are no selectable chains, the move now exits cleanly with `(latticeObject, current_energy, 0, False)`.

### Cluster Utilities Fixes (`pimms/cluster_utils.py`)

- Fixed `build_interface_envelope_pairs()` to handle edge cases safely and deterministically.
	- Added explicit dimensionality validation (2D/3D only) with a clear `ValueError` for unsupported dimensions.
	- Added empty-input handling to return correctly shaped empty arrays instead of crashing.
	- Reworked pair aggregation to avoid out-of-bounds search loops when all queried positions produce zero interface pairs.
	- Fixed 3D aggregation so the first non-empty site is not duplicated.

- Fixed `build_interface_envelope_pairs_safe_and_slow()` runtime errors.
	- Corrected call site to use `lattice_utils.build_envelope_pairs(...)`.
	- Removed stale/undefined symbol dependency on `numpy_utils` in this routine.
	- Added robust zero-pair behavior returning properly shaped empty arrays.
	- Added unsupported-dimensionality validation with a clear `ValueError`.

- Improved `convert_positions_to_single_image_snakesearch()` failure behavior for disconnected input.
	- Added an explicit connectivity guard: when search frontier is exhausted while positions remain, the function now raises a clear `ValueError` indicating the cluster is not connected under the current `space_threshold`.

### Acceptance Module Fixes (`pimms/acceptance.py`)

- Added explicit temperature validation in `AcceptanceCalculator`.
	- `__init__` and `update_temperature()` now reject non-positive temperatures with `AcceptanceException`.
	- This prevents divide-by-zero errors and non-physical inverse-temperature state.

- Added explicit move-selection index validation in move-log updaters.
	- `update_move_logs()` and `megastep_update_move_logs()` now reject out-of-range move indices.
	- This fixes silent Python negative-index behavior (e.g. `selection=-1` mutating the last move bucket).

### Analysis Output I/O Fixes (`pimms/analysis_IO.py`)

- Fixed output path prefixing for prefixed analysis files.
	- Prefixes are now applied to the basename while preserving any configured output directory.
	- This resolves invalid paths generated by raw string concatenation when config paths are absolute (e.g. `"pre_/abs/path/file"`).

- Improved cluster composition output robustness and performance.
	- Added deterministic chain-type ordering for stable output file generation.
	- Replaced repeated nested scans with precomputed per-cluster type fractions.
	- Added safe handling for empty clusters (writes `0.0000` fractions instead of dividing by zero).

- Added input validation for list-length mismatches.
	- `write_scaling_information(...)` now raises `ValueError` when `all_nu` and `all_R0` lengths differ.
	- `write_residue_residue_distance(...)` now raises `ValueError` when `R2R_info` and `all_data` lengths differ (instead of silent truncation via `zip`).

### Energy Module Fixes (`pimms/energy.py`)

- Added explicit dimensionality validation in energy evaluation paths.
	- Short-range and non-short-range local evaluators now raise `EnergyException` for unsupported dimensions.
	- Angle-energy evaluators and angle-lookup construction now explicitly validate supported dimensions (2D/3D).

- Fixed angle-lookup edge case for solvent-only/no-angle-enabled systems.
	- `build_angle_interactions(...)` now creates a valid minimal zero lookup and returns cleanly when no residue penalties are defined.
	- This prevents crashes from indexing an empty penalty key list.

- Reduced overhead in total-energy evaluation.
	- Replaced repeated per-chain `np.concatenate(...)` calls with chunk collection followed by a single concatenation.
	- This avoids repeated array reallocations in the chain loop.

- Minor cleanup.
	- Removed unused import (`longrange_utils`) from `energy.py`.

### Lattice Utilities Fixes (`pimms/lattice_utils.py`)

- Added explicit dimensionality guards in core lattice helpers.
	- `same_sites(...)` now raises `LatticeException` for mismatched dimensionality and unsupported dimensions.
	- `get_gridvalue(...)` now rejects unsupported lattice dimensions with a clear `LatticeException`.
	- `set_gridvalue(...)` now rejects unsupported position dimensions with a clear `LatticeException`.

- Added robust empty-input handling for envelope and center-of-mass utilities.
	- `build_envelope_pairs(...)` and `build_all_envelope_pairs(...)` now return correctly shaped empty arrays/tuples for empty position lists.
	- `center_of_mass_from_positions(...)` now raises `LatticeException` for empty positions instead of failing indirectly.

### PDB Utilities Fixes (`pimms/pdb_utils.py`)

- Hardened `write_positions_to_file(...)` input validation and boundary logic.
	- Added explicit guards for empty input, invalid shape, and unsupported dimensionality (2D/3D only).
	- Corrected bounds checking to reject coordinates equal to box dimensions (`>=`), preventing out-of-range coordinates from being accepted.
	- Improved inferred-dimension behavior by computing `max + 1` per axis for valid 0-indexed box extents.
	- Reduced overhead by using vectorized axis maxima instead of repeated transpose/index extraction.

- Added dimensionality validation to `build_cryst_line(...)`.
	- The function now raises `PDBException` for non-2D/3D dimension vectors with a clear message.

### Keyfile Parser Fixes (`pimms/keyfile_parser.py`)

- Fixed restart-chain propagation bug in restart sanity processing.
	- `sanity_check_and_update_with_restart_file()` now updates the canonical `CHAIN` keyword instead of writing restart-derived chains to `CHAINS`.
	- This ensures downstream concentration/reporting and chain-based checks use restart-updated chain composition.

- Hardened parser behavior for keyword/value splitting and malformed multi-field keywords.
	- `parse(...)` now splits each input line once on `:` (`split(':', 1)`), allowing valid values that contain colons (for example, path-like values).
	- Added explicit format and conversion validation for `CHAIN`, `EXTRA_CHAIN`, and `ANA_RESIDUE_PAIRS`, raising `KeyFileException` with clear messages instead of leaking raw conversion/index errors.

- Minor parser efficiency improvement.
	- `parse(...)` now iterates line-by-line directly from the file handle rather than materializing the full file content first.

### Move Engine Fixes (`pimms/moves.py`)

- Fixed sparse/non-contiguous chain-ID handling in cluster move validation paths.
	- `cluster_translate(...)` and `cluster_rotate(...)` previously rebuilt chain-position dictionaries using `range(1, len(chains)+1)`, which can fail when chain IDs are not contiguous.
	- Both paths now iterate over actual `latticeObject.chains` keys, preventing `KeyError` and ensuring connected-component checks are correct for arbitrary chain IDs.

- Fixed incorrect move accounting in multichain TSMMC.
	- `multichain_based_TSMMC(...)` previously returned `total_moves = 0` even when proposals were made.
	- `total_moves` now correctly reports `steps_per_temperature * num_temps`.

- Hardened multichain TSMMC chain-selection edge cases.
	- Added an explicit early return when all chains are frozen (or no selectable chains remain): `(latticeObject, current_energy, 0, False)`.
	- Ensured `max_number_selectable` is computed as a bounded Python `int` before `random.randint(...)`.

- Reduced repeated membership-check overhead in move hot paths.
	- Converted repeated `chainID in frozen_chains` list checks to `set` lookups in cluster-translate, cluster-rotate, and multichain-TSMMC preprocessing.

- Minor module cleanup.
	- Removed duplicated `import copy` statement.

### Restart Engine Fixes (`pimms/restart.py`)

- Hardened `RestartObject` initialization to avoid uninitialized-state failures.
	- `__init__` now initializes `hardwall`, `chains`, and `seq2chainType` in addition to existing fields.
	- `dimensions` now starts as an empty list instead of scalar `0`, matching expected list semantics throughout the module.

- Fixed non-atomic position-offset behavior during lattice-dimension updates.
	- `__apply_position_offset(...)` now validates all proposed shifted positions first, then applies updates in a second pass.
	- This prevents partial in-place mutation when one shifted coordinate is invalid.

- Made `update_lattice_dimensions(...)` transactional for failures.
	- If offset application fails, prior dimensions are restored before re-raising `RestartException`.

- Improved robustness of `add_extra_chains(...)` parsing and ID/type assignment.
	- Added strict parsing for chain count and sequence with clear `RestartException` messages.
	- Added positive-count validation for extra-chain requests.
	- Fixed edge-case failures when base chains are empty and when no existing chain types are recorded.
	- Corrected malformed error-string formatting so the offending payload is interpolated.

- Removed aliasing to mutable lattice position data.
	- `build_from_lattice(...)` now deep-copies chain positions instead of storing live references.
	- This prevents accidental restart-state mutation when the source lattice mutates later.

- Prevented stale state leakage when rebuilding restart objects.
	- `build_from_lattice(...)` and `build_from_file(...)` now reset `extra_chains` when reconstructing object state.

- Strengthened restart-file loading and schema validation.
	- `build_from_file(...)` now uses context-managed file I/O.
	- Expanded pickle-read exception handling to include common corruption/parse failures (e.g. unpickling and EOF errors).
	- Added explicit validation that `CHAINS` is a dictionary.
	- Added explicit validation for malformed chain entries and position dimensionality mismatches relative to `DIMENSIONS`.

- Strengthened restart-file writing safety.
	- `write_to_file()` now uses context-managed file I/O for reliable file-handle cleanup.

### Numpy Utilities Fixes (`pimms/numpy_utils.py`)

- Repaired `position_in_list(...)` (previously hard-disabled).
	- Removed unconditional `BrokenException` raise and restored functional element-wise position matching.
	- Added robust handling for empty inputs and both list/NumPy-array call patterns.

- Hardened `tetrahedron_volume(...)` input handling.
	- Added explicit shape-consistency checks across all four point arguments.
	- Added 3D-coordinate validation with clear `ValueError` for invalid dimensionality.
	- Normalized single-tetrahedron and batched inputs via `np.atleast_2d(...)` for consistent behavior.

- Hardened `find_nearest(...)` for broader and safer usage.
	- Added conversion from generic sequence input to NumPy arrays, so Python lists are handled reliably.
	- Added explicit empty-input rejection with a clear `ValueError`.
	- Flattened array input before nearest-index selection to avoid shape-dependent surprises.

### Parameter File Parser Fixes (`pimms/parameterfile_parser.py`)

- Replaced unsafe process exits in `parse_energy(...)` with structured exceptions.
	- Invalid numeric interaction values now raise `ParameterFileException` instead of calling `exit(1)` from library code.
	- This makes parser failures catchable and testable by callers.

- Added robust parsing for interaction integer fields.
	- Introduced centralized integer parsing validation for short-range, long-range, and semi-long-range interaction terms.
	- Non-numeric and float values now report explicit, line-localized parser errors.

- Fixed blank-line and comment-only-line safety in both parsers.
	- `parse_energy(...)` now skips lines that are empty after comment removal, preventing index errors.
	- `parse_angles(...)` now does the same.

- Fixed comment-aware tokenization bug in `parse_angles(...)`.
	- Angle parsing now tokenizes `un_comment` content, ensuring inline comments do not corrupt field parsing.

- Minor parser efficiency/clarity cleanup.
	- Removed repeated `list(dict.keys())` containment checks in hot parsing loops in favor of direct dict membership tests.
