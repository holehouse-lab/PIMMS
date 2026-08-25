## 1.0.8 (August 2026)

For version 1.0.8 we coducted a deep audit of the whole codebase. This identifed a number of edge-case but important bugs and shortcomings that have now been fixed, although we note none of these affect any previously published work (i.e. all in new moves/features introduced over the last few months).

### Detailed balance and sampling correctness

- **Bug fix (VMMC boundary links violated detailed balance).** For links that were tested but did NOT form, with the partner chain left outside the cluster, `vmmc_move` computed the reverse link probability from the seed-displaced pair energy (`Er`) - the quantity that is only correct for links whose partner moves WITH the cluster. The correct boundary reverse probability comes from the forward-displaced energy (`1 - exp(-beta(E0 - Ef))`). The bias systematically over-favoured contact states: exact enumeration of two single-bead chains gave P(contact) = 0.86 observed vs 0.43 exact at T=3. After the one-line correction every enumerated case matches the exact Boltzmann distribution to within statistical error (verified at three temperatures and with three-chain systems that exercise the formed-link path, which was and remains correct).
- **Bug fix (hardwall VMMC wrapped wall-crossing proposals).** The virtual link energies correctly treated a bead translated outside a hardwall box as absent, but the commit path always applied `pbc_convert`, turning that same invalid proposal into an in-box periodic wrap. Single-bead chains evaded the later straddle check entirely. The committed state therefore differed from the state used to construct the Metropolis-Hastings ratio and biased a two-particle hardwall equilibrium. Hardwall VMMC now rejects every raw translated coordinate outside the box before placement; an independent full-state enumeration recovers the exact Boltzmann contact probability and checks every sampled coordinate remains in bounds.
- **Bug fix (`chain_rotate` was not reversible).** Rotations were applied about the *rounded* circular-mean centre of mass; the rounded COM of the rotated chain is generally a different lattice point (guaranteed rounding ties for even-length chains, circular-mean drift for boundary-straddling ones), so the inverse rotation did not return the original configuration - about a third of accepted rotations were not invertible, and an exact-reference Metropolis test showed a consistent equilibrium bias. The move now rotates the chain's displacement vectors about the bead nearest the (single-image) centroid and re-anchors them on that bead's original in-box position; this is exactly invertible for any box shape and for straddling chains. The equilibrium bias vanishes (translate+rotate now matches exact enumeration as tightly as the translate-only control).
- **Bug fix (`cluster_rotate` had the same defect, plus a straddle bug).** Cluster rotations used the rounded circular-mean COM of the *raw* wrapped positions, so a boundary-straddling cluster was not even rotated as a rigid body. The move now single-images the cluster (snakesearch) and rotates displacement vectors about the bead nearest the cluster centroid: 0/1949 accepted rotations are non-invertible (was 636/1949). This matters more than the chain case because energy-neutral cluster rotations are accepted unconditionally, so nothing compensated the bias.
- **Bug fix (`chain_pivot` off-by-one: guaranteed null moves, immobile 3-mers, asymmetric termini).** The N-terminal branch anchored on bead `pivot_point - 1`, so `pivot_point = 1` "rotated" only the anchor - a guaranteed null move that was still fully energy-evaluated and logged as accepted. For 3-mers that was the ONLY pivot point, so `chain_pivot` never moved a 3-mer at all; for L=4 half of all pivots were null and the C-terminal bead could never move. Both arms now anchor on bead `pivot_point` and the shorter-arm comparison uses `(L-1)/2`, so every eligible pivot moves beads and both termini are reachable (null-move fraction measured at exactly 0 for L=3..8, was 1/(L-2)).
- **Bug fix (system-TSMMC protocol was not palindromic).** The temperature-excursion bookkeeping ran the first M-1 sub-moves at the *target* temperature and never visited the last schedule temperature, executing `[target x(M-1), s0 xM, ..., s_(L-2) xM]` - not a palindrome, which the tempered-transitions acceptance requires for detailed balance with Metropolis sub-kernels. `start_system_TSMMC` now switches to the first schedule temperature immediately (accumulating that work term at the starting energy), the check-in advances after each block of M sub-moves, and completion fires at exactly L*M sub-moves: the executed protocol is now literally `[s0 xM, ..., s_(L-1) xM]`, verified move-by-move for M = 1, 5, 20.
- **Bug fix (TSMMC linear schedule overshot the jump temperature).** `np.arange` with a float step and an exclusive float endpoint produced `number_points + 1` ramp elements for ~6% of integer (target, jump, points) combinations, adding a temperature a full step ABOVE `TSMMC_JUMP_TEMP` and making the schedule non-monotonic. The ramp is now `np.linspace(target, jump, points+1)[1:]` - both endpoints exact, always monotone, still palindromic.
- **Bug fix (monomer HARDWALL teleport).** In every crankshaft/slither kernel (serial reference, fast serial, and parallel, 2D and 3D), a single-bead chain's proposed move was PBC-wrapped and never straddle-checked under HARDWALL, so monomers walked straight through the wall (reproduced: a monomer confined to a corner crossed the box in one step). All ten kernel proposal sites now apply the same straddle rejection multi-bead moves already had; fast and reference kernels remain bit-identical, and a new test drives a hardwall monomer 3000 steps per kernel asserting the per-step displacement never exceeds one site.
- **Bug fix (PARALLELIZE silently never moved some chains).** The parallel slither/pull kernels only move a chain when it fits inside one checkerboard block interior, and hard-skip chains longer than their 512-bead stack buffers (hetero chains for slither, all chains for pull). A chain failing either condition was silently frozen on EVERY sweep - reproduced: a 40-mer in a 50-box accepted ~6000 serial slithers and exactly 0 parallel ones. The Python dispatch now measures each chain's periodic extent against the block-interior geometry (via a new `parallel_layout_info` introspection helper) and the stack caps, and falls back to the trusted serial kernel whenever any chain could never move - so PARALLELIZE again changes only speed, never the sampling. Regression test included.

### Simulation loop, quench and acceptance

- **Bug fix (system-TSMMC master-step bookkeeping).** An excursion was started after that step's quench/resize updates, but completion incremented the master counter again and repeated those updates, assigning one proposal two step numbers. If selected on the final step, the loop instead stopped with the auxiliary move unfinished. Auxiliary substeps now hold the master counter fixed, completion reports the resulting state at the original step without duplicating pre-move updates, and the loop remains alive past `N_STEPS` only long enough to finish an in-progress excursion.
- **Bug fix (scheduled output described the wrong state).** `ENERGY.dat`, trajectory frames and analysis were written before each move while labelled with that move's step, so every row/frame was one move behind and the actual final state was absent. Scheduled output now runs after the completed move (including self-contained mega/TSMMC/VMMC moves). The temperature still changes before a quench step, but its `QUENCH.dat` energy is deferred until after that move too. Fully frozen and individually fixed no-op steps still run quench/resize/output/analysis instead of silently skipping all scheduled work.
- **Bug fix (`move_selector` could crash a long run).** Move thresholds are a cumulative float sum; the parser tolerates |sum - 1| up to 1e-7 and float accumulation itself leaves the top of the last interval a hair below 1.0 (ten 0.1s sum to 0.9999999999999999). A random draw in that gap raised `AcceptanceException` and killed the run - expected roughly once per 1e8 steps. The last non-zero interval's upper bound is now snapped to exactly 1.0, but ONLY when it is already within 1e-6 of 1.0, so genuinely unassigned probability mass still raises as intended.
- **Bug fix (quench arithmetic truncated fractional temperature ranges).** `QUENCH_AS_EQUILIBRATION` sized the equilibration window from `int(dT)`, so e.g. a 1.0 -> 0.2 quench in 0.1-steps computed ONE quench period instead of nine and production analysis ran while still quenching. Now `ceil(dT / stepsize)`. The `print_summary` quench-step count also used the sign-flipped stepsize and printed a negative count for heating runs - fixed with the same formula.
- **Bug fix (cooling quenches could not pass below T=1).** `update_temperature_in_quench` raised `TemperatureException` below 1, although TEMPERATURE is a float and the acceptance machinery only requires it positive; a legitimate 1.5 -> 0.5 quench aborted mid-run. The guard is now `<= 0`.
- **Bug fix (heating quench + TSMMC crashed or silently inverted).** When a heating quench reaches `TSMMC_JUMP_TEMP`, rebuilding the TSMMC schedule hit `ZeroDivisionError` (dT = 0), and past it the "excursion" silently became a cooling ramp. `TSMMC.__init__` now raises a clear error whenever the jump temperature is not above the target, and `TSMMC_FIXED_OFFSET` must be positive.
- **Bug fix (run-time message dropped whole days).** The end-of-run "Simulation time" used `relativedelta` fields, which normalise >24 h into `.days` - a 50-hour run reported "2 hours". Now computed from `total_seconds()`.
- The `update_dimensions` "forced move" comment and docstring described a chain-forcing mechanism that does not exist (the returned override list was never consumed); both now describe the actual behaviour (the equilibration window is extended by 100 steps).

- **Bug fix (frozen beads could silently be the wrong beads).** The crankshaft kernels' frozen-bead selector iterated the chains dict in insertion order while the global bead indices are assigned in ascending chainID order; a restart file whose pickled chains dict was not in ascending order therefore froze the *wrong* beads without any warning. The selector now iterates `sorted(chains.keys())`, matching the index assignment.

### Keyfile parser and inputs

- **Bug fix (restart overrides bypassed the box sanity checks).** `RESTART_OVERRIDE_DIMENSIONS` replaces `DIMENSIONS` *after* the keyfile-level checks have already run, so a restart file with a too-small (< 7) box, or a non-cubic box combined with `MOVE_CLUSTER_ROTATE` under PBC (the energy-corruption guard added in 1.0.5), sailed straight through. The overridden box is now re-validated at the point of override.
- **Bug fix (`CHAIN : 0 SEQ` and negative counts accepted).** A zero/negative chain count parsed fine and crashed much later inside mdtraj's XTC writer with an opaque buffer error. Counts below 1 are now rejected at parse time (CHAIN and EXTRA_CHAIN).
- **Missing range checks added**: every `RESIZED_EQUILIBRATION` axis must satisfy the same >= 7 floor as `DIMENSIONS` (that box IS simulated); `EQUILIBRATION_OFFSET` values must be >= 0; `QUENCH_FREQ` must be positive (0 silently produced a zero-length quench-equilibration, negatives failed only at runtime); `LATTICE_TO_ANGSTROMS` must be positive.
- **Bug fix (`ANA_RESIDUE_PAIRS` never validated against restart/EXTRA_CHAIN chains).** The bound check looped over the keyfile's own `CHAIN` lines, which are empty in a restart run, so out-of-range pairs reached the analysis code. The check now re-runs after restart processing over the rebuilt `CHAIN` list plus `EXTRA_CHAIN`.
- **Bug fix (`ANALYSIS_MODULE` without `ANA_CUSTOM` silently never ran).** A custom analysis module was loaded and validated at parse time and then never executed, because `ANA_CUSTOM` defaults to 0. This is now a parse-time error telling the user to set the frequency.
- **Bug fix (`write_keyfile` output could not be re-parsed).** List values were written as Python reprs (`[20, 20, 20]`) and loaded objects (restart, freeze file, the analysis callable) verbatim. Multi-entry keywords now serialise one `KEY : n seq` line per entry, lists as space-separated values, object-valued and derived `__` keywords are skipped - a written keyfile now round-trips through the parser.
- Removed the over-restrictive `RESIZED_EQUILIBRATION` branch of the cluster-rotate cubic-box check: the resized-equilibration phase always runs hardwall (forced in `Simulation.__init__`), where a rigid rotation is a valid isometry on any box shape. Only the production box needs to be cubic/square for `MOVE_CLUSTER_ROTATE` under PBC.
- The parameter-file docs/error message now mention the 5-column `A B X Y Z` (SLR) line form the parser has always accepted.
- A completely full lattice now aborts chain insertion with the intended "lattice is overcrowded" `ChainInsertionFailure` message instead of an unhandled `LatticeUtilsException` crash from the empty-site search.

### Analysis and output correctness

- **Bug fix (AUTOCENTER broke for boundary-straddling chains).** The centring shift was computed from the periodic circular-mean COM of positions that had already been unwrapped to a single image, so for a chain whose true centre sat at or beyond a box face the reference re-wrapped to ~0 and the "centred" chain was written entirely outside the box - exactly the chains AUTOCENTER exists to tidy up. The centre is now the plain (rounded) arithmetic mean of the single-image positions.
- **Bug fix (cluster Rg/asphericity biased for asymmetric clusters).** `extract_cluster_polymeric_properties` receives single-image positions but computed the gyration tensor with the periodic correction ON, in a bounding box barely larger than the cluster - re-wrapping beads more than half that tiny box from the COM and inflating Rg by up to ~15% for elongated clusters. The single-image path now computes plain (non-periodic) gyration properties, and the misleading "the gyration tensor is unaffected" comment is corrected.
- **Bug fix (hardwall observables were reported with periodic distances).** `Chain` discarded its hardwall flag, so end-to-end distances, residue distances, internal scaling, distance maps, Rg and asphericity all used minimum-image geometry. Main-simulation SR/LR cluster discovery also omitted the hardwall argument and could connect chains through opposite faces. The boundary convention is now retained through de-novo, restart and resized lattices; hardwall observables use Cartesian coordinates, hardwall cluster properties skip periodic single-imaging, and both cluster graphs honour the walls.
- **Bug fix (super-long-range contacts missing from LR clusters).** The long-range connected-component builder requested SR, LR and SLR envelope pairs but concatenated only the first two, splitting clusters connected at Chebyshev distance three. SLR pairs are now included and LR-cluster single-imaging uses the matching three-site connectivity threshold.
- **Bug fix (radial density profile silently all-zero for edge clusters).** The profile centred on the *periodic* COM of single-image input; when the rounded COM landed past the box edge it wrapped to ~0 while the beads sat at ~DIM, every Chebyshev distance overflowed the profile range, and the whole profile zeroed out. The centre is now the plain rounded mean of the single-image positions; profiles are verified translation-invariant across the boundary.
- **Bug fix (filtered radial profiles were assigned to the wrong cluster).** Small clusters are intentionally skipped by the bead threshold, but the writers renumbered the surviving profiles from `C1`; if an earlier heterogeneous cluster was below threshold, its later profile was silently attributed to another cluster. The original one-based cluster indices are now carried to both SR and LR radial-profile files.
- **Bug fix (PDB writes crashed at scale).** Atom serials overflowed the 5-column field at 100,000 beads (`PDBException` at topology-write time on large systems) - serials now wrap modulo 100000 (readers rebuild indices sequentially; CONECT records for wrapped serials are skipped as ambiguous). Separately the TER record carried a residue number one PAST the terminal residue, which both violated the PDB spec and crashed outright for chains of exactly 9999 beads - it now equals the last residue's number.
- **Bug fix (scaling-exponent fit).** With no internal-scaling sample ever folded in (analysis frequency beyond the production length) the fit ran `log(0)` and wrote `nan nan` to `SCALING_INFORMATION.dat` - it now returns the `(-1, -1)` sentinel. The log-spaced point selection also never included the largest sequence separation (the most informative point of a scaling fit); the endpoint is now always in the fit set.
- **`CHAIN_<type>_CLUSTERS.dat` rows now lead with the step number**, like every other per-step analysis file (they previously relied on line-index alignment with `CLUSTERS.dat`).
- **Bug fix (log.txt mixed two clocks).** The logger header used local time while every entry used UTC, so entries appeared hours before/after the header. All timestamps are now local time.
- `ANA_CLUSTER_THRESHOLD` docs corrected: the comparison is strict (`>`, so the default of 1 *excludes* single chains) and gates only the per-cluster property analysis - the size-distribution files always include every component. The 2D convex-hull conventions are now documented (`CLUSTER_VOL` = polygon area, `CLUSTER_AREA` = perimeter in 2D), the dead pre-scipy-0.17 manual hull fallback is removed, and the 2D "asphericity" is documented as kappa (the 3D value is kappa squared - the historical 2D form is kept for continuity).
- `SAVE_EQ` docs now state that frame 0 of `traj.xtc` is always the starting configuration regardless of the setting, and `print_summary`'s expected-frame formula now matches the writer's real cadence when EQUILIBRATION is not a multiple of XTC_FREQ.
- The `MOVE_FREQS.dat`/`ACCEPTANCE.dat` loadtxt example in the docs dropped its `delimiter="\t"` (rows end with a trailing tab, which that delimiter parses as an empty column and rejects).

### lemonade

- **Bug fix (HARDWALL trajectories analysed as periodic).** `TrajectoryStore` stored `hardwall` and never consulted it: two chains against opposite walls were merged into one cluster and the single-image gather dragged one of them across the wall, corrupting every condensate statistic for hardwall runs. The connected-component search now honours the wall (via a new `hardwall` argument threaded through `get_cluster_distribution`) and `single_image_positions` is the identity under hardwall. Regression test included.
- **Bug fix (wrong lattice spacing silently corrupted the load).** The lattice round-off residual warning fired only under `verbose=True`; a wrong/omitted `LATTICE_TO_ANGSTROMS` silently produced broken bonds and a mis-sized box. The residual check now always warns above 0.05.
- **Bug fix (`Cluster.radial_density_profile(minimum_cluster_size_in_beads=n)` raised `IndexError`)** for clusters smaller than `n`; it now returns `None`.
- **Bug fix (lower-case keyfile chains silently dropped).** `parse_only=True` skips PIMMS's sanitiser, so lower-case `CHAIN` sequences never matched the (upper-case) PDB residue names and the keyfile chain types were quietly discarded. The loader now honours `CASE_INSENSITIVE_CHAINS` (default True) and upper-cases the specs.
- **`lemonade.load` now cross-checks the box**: a keyfile `DIMENSIONS` that disagrees with the trajectory's own CRYST1/XTC box record (wrong keyfile, or the `eq_` trajectory of a `RESIZED_EQUILIBRATION` run) triggers a loud warning instead of silently breaking every periodic calculation.
- **Surface tension estimator corrections**: the capillary-wave fit now uses the exact lattice dispersion `2 - 2cos(q)` per axis instead of the continuum `q^2` (which under-estimated gamma by 2-10% at typical PIMMS box sizes; verified against synthetic Gaussian height fields). Independent modes are now selected by canonical `+q/-q` wavevector pairs before applying `n_modes`; the old value-based `unique(gamma)` could collapse unrelated modes that merely had equal power and report the wrong sample count/spread. `n_modes` itself is now validated (a bool, non-integer or value below 1 raises `ValueError` instead of silently slicing an empty or garbage mode set).
- `analyze(geometry='slab')` no longer reports `shape` statistics: the convex hull of a box-spanning slab is exactly the quantity the docs warn is meaningless. Density/asphericity documentation now flags the two places lemonade's conventions differ from PIMMS's own output.

### Performance

- **The energy-delta kernel - the hot core of every crankshaft/slither/pull substep - is ~35% faster.** `get_energy_change` recomputed three modulo wraps and three hardwall tests per SITE of the 7x7x7 (and 3x3x3) shells; the wrapped index along each axis depends only on the offset, so 42 per-axis values now get precomputed per call and the inner 343-site loops just index them. Measured: SLR substeps 1206 -> 760 ns, SR 171 -> 129 ns, bit-exact with the reference kernels (the full 250-case equivalence suite passes unchanged).
- The four hyperloop hardwall energy branches and the `inner_loops_hardwall` 3D pair extractor were calling Python-level `abs()` per pair (boxing through PyNumber_Absolute, ~5-10x slower than the PBC twins); both now use C comparisons. `get_adjacent_sites_2D` chained indexing (`positions[i][0]`) defeated its typed buffer and is now a single 2D index.
- `C_RAND_MAX` is now the fixed constant 2^31 - 1 instead of the platform libc `RAND_MAX`. The kernels have used splitmix64 (fixed internal range) since 1.0.5, so the libc value was only a leftover seed-modulus - but on Windows it is 32767, collapsing the per-move seed space to 15 bits and making seed streams platform-dependent. macOS/Linux behaviour (and reproducibility of existing seeds there) is unchanged; `get_randmax` is no longer consulted.

### Test suite

- **The detailed-balance tolerance is now statistically honest.** The old criterion (2.5 x per-sample std + 3% of |E|) was ~30x the standard error of the mean on the actual fixtures, and used the *test* trace's own spread - so a broken move inflated its own tolerance: an injected acceptance bug with beta scaled x1.5 (a gross Metropolis error) passed EVERY detailed-balance case. The tolerance is now 4 x the autocorrelation-aware SEM of the difference plus a 0.5% floor, computed from each trace's own mean uncertainty; a meta-test pins that the beta x1.5 bug now fails (and that a correct kernel still passes against itself).
- **The simulation regression suite now provably tests THIS working tree.** It previously preferred whatever `PIMMS` executable was on PATH - in a typical dev environment a stale pip-installed version, so all 15 end-to-end tests were validating the wrong code (confirmed: the run logs said version 1.0.1). The harness now always runs `scripts/PIMMS` with the repo forced onto the subprocess's `PYTHONPATH`, and asserts the subprocess imports `pimms` from inside the repo before trusting any output.
- **The regression comparisons are no longer blind to write cadence**: alongside each file's final line the suite now pins the file's non-empty line count (an EN_FREQ/ANALYSIS_FREQ off-by-one changes the count but not the final line), and `QUENCH.dat` joined the captured set. All expected outputs were regenerated from the fixed tree - the detailed-balance and TSMMC-protocol fixes above legitimately change every trajectory (each regenerated run passes its internal from-scratch ENERGY_CHECK).
- **New regression tests for every fix in this release** (`test_review_fixes.py` and focused module tests): monomer hardwall confinement across all kernels, chain/cluster rotation reversibility, chain-pivot mobility for L=3..6 and both termini, periodic and hardwall VMMC against exact Boltzmann enumeration, independent Cartesian cluster/chain tensors, hardwall clustering, SLR connectivity, filtered radial-profile labels, independent capillary wavevectors, post-move output/frozen scheduling, move-selector threshold closure, empty LR-envelope shapes, PARALLELIZE serial fallback, TSMMC schedule/protocol/guards, and the detailed-balance-threshold meta-test.
- Test hygiene: the quench unit tests no longer append a `log.txt` wherever pytest is invoked from; `test_chain`'s ~50-line stub-module loader - which was provably inert (the stubs were always bypassed by `from . import` resolution) - is removed in favour of the real modules it was actually testing all along; the PBC-aware centre-of-mass finally has value-level tests (including boundary-straddling cases); and the fixture maintenance scripts are fixed (`run_up_test_sims.sh` was an accidental byte-for-byte copy of the cleanup script with an invalid shebang, and both stopped at test_13 - it now actually runs all fixtures against the working tree and regenerates the baselines).

### Analysis-code audit (oracle-verified)

The same maximum-rigor treatment applied to every analysis function: four parallel audits (per-chain observables, the cluster pipeline, lemonade's core, and the phase-separation/surface-tension estimators), each function verified against independently written oracles - closed-form configurations, union-find clustering, hand-counted shells, synthetic profiles and fluctuation fields with known ground truth - or reported as a concrete bug. Five genuine numerical bugs found and fixed:

- **Bug fix (every PBC radius of gyration was inflated by a parallel-axis error).** The gyration tensor was referenced to the circular (Bai-Breen) COM instead of the arithmetic mean of the reconstructed single-image coordinates; by the parallel-axis theorem every Rg^2 written under periodic boundaries carried a strictly one-sided bias of exactly |mean - circular|^2 (up to ~1e-2 relative at box 7, visible in the second decimal of RG.dat; asphericity perturbed both ways; verified over 1959 configurations to 1.6e-15 as pure parallel-axis). One line - recompute the reference from the reconstructed coordinates - makes the PBC path satisfy the definition exactly and bit-match the Cartesian path for non-straddling chains. Independently confirmed by a second audit re-running an identical-seed simulation: RG.dat now matches the definition to pure print rounding.
- **Bug fix (single-image gather broke box congruence).** snakesearch's shift-everything-to->=0 step was not a whole number of box periods, so the "single image" of any cluster straddling a LOW box face was in an arbitrary translated frame: lemonade's `Cluster.center_of_mass` (reduced mod box) landed up to a cluster radius away from the truth, and `radial_density_profile` - which bins raw wrapped positions about that COM - was smeared to the point of reporting solid-droplet cores as 35-54% empty. Found independently by two audits. The shift is now a whole number of box periods in both the kernel and the Python fallback (all shape-based consumers are translation-invariant and unchanged); straddled and centred droplets now produce identical profiles.
- **Bug fix (the LR-cluster search was directional - its "clusters" were not even a partition).** LR/SLR envelope pairs are only emitted FROM LR-flagged beads, but the component search treated the relation as undirected reachability from arbitrary seeds: in any mixed-flag system (some residue types without LR columns) the same configuration decomposed into seed-dependent, overlapping clusters that double-counted chains (a real run wrote 25 memberships for 20 chains, 7 clusters where the truth is 4), contaminating every LR_* output file. Connectivity is now symmetrised with the energy-based definition - an LR/SLR edge exists only when BOTH beads are LR-capable (exactly the pairs with nonzero energy; SR contacts connect as before) - verified against a union-find oracle on 70/70 randomized systems (all-flagged and mixed-flag, PBC and hardwall) with the partition property holding everywhere. The docs now state this definition precisely.
- **Bug fix (lemonade radial profiles applied the periodic metric to hardwall trajectories).** Material at true distance beyond half the box was folded into inner shells. Hardwall profiles now use plain Cartesian distances with per-COM wall-truncated shell normalisation, and the default range extends to the box diagonal; far material appears at its true radius.
- **Bug fix (`extract_cluster_polymeric_properties` crashed on numpy-array dimensions)** - ambiguous truth test, now `is False`.

Hardening and reporting from the audit: percolating (box-winding) clusters are now DETECTED at the single-image gather and warned about explicitly - their gathered coordinates are BFS-order dependent and any shape quantity computed from them is meaningless (detection is exact: a connected non-winding cluster always fits one image); a fluctuation-free slab returns the explicit gamma = +inf sentinel instead of numpy divide-by-zero warnings; unknown `geometry`/`by` selector strings now raise instead of silently running a different analysis (with droplet/sphere accepted as cross-module synonyms); zero-length chain sequences are rejected at Topology construction (reduceat garbage otherwise); the empty residue-pair analysis returns instead of pointlessly reopening its output file every step; and stale docstrings were corrected. Regression tests pin every fix; the per-chain, cluster, lemonade-core and phase/surface-tension ledgers record positive verification for everything else - including exact recovery of known surface tensions from synthetic capillary fields (slab, within 2.2% at stated ~5% detection power, exact to 1e-12 on deterministic fields) and known gamma from synthetic fluctuating spheres (droplet, within ~1%), tanh-fit parameter recovery with stated power, and every production .dat value matched to oracle recomputation on real simulations across the full configuration matrix.

### Physics-core audit (energy + moves, oracle-verified)

A five-way maximum-rigor audit of every energy-calculating function and every move, with the strictest evidence standard of the release: each function had to earn either a positive correctness argument against an INDEPENDENTLY written oracle (from-scratch numpy Hamiltonians coded straight from the physical definitions - never validated against other in-repo code) or a concrete reproducible bug. Coverage: exact table enumeration; 504 oracle-matched total-energy configurations plus ~75 hand-computed adversarial cases; move-by-move delta verification over thousands of real moves with bitwise revert checks; exact proposal-set enumeration; and full-state-space exact-Boltzmann equilibrium tests with a validated-power control (an injected beta x1.35 error is flagged at |z| ~ 23).

**The energy pipeline is clean.** Zero correctness bugs in table construction, `evaluate_total_energy` (all five decomposition slots bit-equal to the oracle in every configuration), all incremental delta paths (`single_chain_move`/`rigid_cluster_move` bookkeeping, kernel delta algebra, MoveEvent angle windows - integer-exact after every accepted move, bitwise restoration after every reject), and the crank/slither/pull kernels as samplers (exact detailed balance verified transition-by-transition for pull - pi(x)T(x,y) == pi(y)T(y,x) to 1e-10 - and exact-Boltzmann equilibrium for all three, homo and hetero paths).

**Four cluster-move bugs found and fixed** (all verified by reproducer before and after):

- **Bug fix (hardwall cluster translation wrapped through the wall).** `cluster_translate` always applied the periodic wrap, and the per-chain straddle check only inspects consecutive bonds - so a box-spanning cluster under HARDWALL could be "translated" by permuting its chains through the wall: occupancy unchanged, committed as energy-neutral, while the true hardwall SR/solvation energy changed (reproduced: 780 of 780 accepted moves on a spanning-column system carried the wrong dE, permanently corrupting the tracked energy). Hardwall cluster translations now reject any bead whose raw translated coordinate leaves the box - the same convention VMMC uses.
- **Bug fix (rotating a box-winding cluster is not a rigid motion).** A cluster connected to its own periodic image has no consistent single image; rotating its spanning-tree unwrapping maps the winding closure vector onto an axis with a different period, changing intra-cluster minimum-image LR/SLR relations while the move's dE assumed rigidity (found spontaneously in a random sweep: an accepted 8-chain rotation committed dE = -4 while the true dE was +3). Cluster rotations now reject whenever the single-image extent reaches the box length on any axis - a symmetric constraint, so detailed balance is preserved.
- **Bug fix (rigid_cluster_move double-counted intra-cluster LR/SLR pairs).** The per-chain envelope lists were concatenated, so every cross-chain pair inside the cluster was counted once from each end: the computed delta was external + 2x(intra) instead of external + intra. This cancels exactly while intra-cluster geometry is preserved (which is why it survived every equilibrium test) but compounded the winding bug above. The concatenated lists are now exactly de-duplicated (the antisymmetric pair-row ordering makes duplicate rows identical).
- **Bug fix (rotation pivot ties broken by floating-point noise).** The bead-anchored rotations selected the pivot as the argmin of FLOAT distances to the float centroid; two beads exactly equidistant had their tie broken by sub-ulp rounding noise that is not rotation-invariant, so the reverse move could select the other tied bead and fail to invert - a complete detailed-balance violation on tied shapes (measured: 109/10080 forced chain rotations, 9/9 on a tied cluster; a 5-mer counterexample gives pi(x)T(x,y) = 1/3 vs pi(y)T(y,x) = 0 exactly). Pivot selection now uses exact integer arithmetic (argmin of ||n*p_i - sum(p)||^2), which is tie-stable under rotation, and `cluster_rotate` concatenates its chains in sorted-chainID order so the bead ordering (and therefore the tie-break) is canonical in both directions.

Also fixed from the audit: the `chain_homo` fast-path precondition now checks LR-flag uniformity as well as intcode uniformity (implied today by the type-to-LR coupling, now enforced); `EmptyHamiltonian` finally honours the real Hamiltonian's contracts (5-tuple energy decomposition, per-residue -1 LR sentinel); the residue-to-intcode assignment is deterministic across processes (sorted, was set-iteration order - physically unobservable, but a trap for anyone persisting intcodes); `absolute_energies_of_angles.txt` now records both the requested and the APPLIED (rounded-integer) angle penalties; and a series of stale docstrings (NaN claims, pair-ordering claims, the LR-intcode identity, VMMC's boundary-only link convention, the strict cluster-threshold semantics, dead MoveEvent branches) now describe what the code actually does. Regression tests pin all four cluster-move fixes; the six cluster-move regression fixtures were regenerated from the fixed tree (each run passes its internal from-scratch ENERGY_CHECK).

The audit's positive findings are as valuable as the fixes: the angle convention is now stated precisely (displacement-class A1/A2/A3 keyed on the i-1 to i+1 vector, middle-bead lookup, lattice-point-group invariant - which is what makes the rigid moves' empty angle windows legitimate), and every move now carries an explicit verification: proposal symmetry or exact Hastings correction, bitwise reject-restoration, and exact-Boltzmann equilibrium with stated ergodicity closures.

### Cython kernel audit

A dedicated four-way line-by-line audit of all ~11,700 lines of Cython (fast serial kernels, parallel/megamove kernels, reference kernels + hyperloop, and the support modules), against a fixed checklist: bounds safety under `boundscheck(False)`, C-division semantics, integer width/overflow, uninitialized locals, stack-buffer caps, `prange` thread-privacy, RNG state/stream correctness, and proposal-vs-commit consistency. 22 findings, all fixed:

- **Bug fix (per-block parallel PRNG seeds were affine in block + sweep-seed).** Any two (block, sweep) pairs with equal `block + seed` replayed a verbatim-identical random stream: adjacent sweeps shared all-but-one whole-block streams, and by the birthday bound duplicated streams were *expected* within a few thousand parallel megamoves. All six parallel drivers now derive block seeds through a splitmix64-style finalizer (verified collision-free over 32k sweep windows; a regression test pins it). Thread-count independence is unaffected and now pinned in 3D as well as 2D.
- **Bug fix (`randint` was wrong for any start >= 2).** The kernel RNG helper's special-cased formula only honoured its inclusive-range contract for start 0/1 (the only values the kernels use - but the helper is exported). Replaced in BOTH the reference and fast kernels with the general inclusive formula, which is provably bit-identical for start 0/1 - confirmed by the full 250-case bit-exactness suite passing unchanged - and correct for every start.
- **Memory-safety hardening (all latent - none reachable through the engine's real call paths):** a negative substep count wrapped an unsigned loop counter to ~2^32 iterations and read the bead selector far out of bounds (all four crank entry points now no-op on `nsteps <= 0`); the reference kernels' 1.0.7 defensive unknown-flag `else` left an unguarded `[-1]` angle-buffer write that the fast kernels already guarded (now guarded in the references too); the serial slither/pull revert buffers had no in-kernel cap against an understated `max_chain_len` (now skip, mirroring the parallel kernels' 512 guard); `snakesearch_single_image` validated positions but not `seed_idx`/empty input (now raises); lemonade's `unwrap_chains` read one row past the buffer for a zero-length final chain (now skipped).
- **Bug fix (hyperloop's from-scratch energy oracle accumulated in 32-bit int).** The incremental kernels accumulate in `long`; the ENERGY_CHECK oracle wrapped silently past 2^31. All six accumulators widened.
- **Bug fix (packaging: `pimms/cython_backup/` shipped in both the sdist and the wheel).** Three stale pre-fix snapshots of audited kernels - never compiled, pure confusion hazard - were swept in by `graft pimms`. The directory is deleted (git history preserves it), `MANIFEST.in` prunes it defensively, and a hygiene test pins its absence. (A stale `build/` staging tree from an old in-place build was re-injecting the deleted files into fresh wheels via `include_package_data`; it is cleared - delete `build/` after removing any shipped file.) Three fully dead compiled modules (`get_randmax`, `lattice_tools`, `random_number` - zero live callers since the splitmix64 migration) are removed along with their extensions and imports, as is the dead `delete_pbc_pairs`.
- **Structural alignment**: the 2D reference kernel now zeroes the y slot on failed proposals exactly like the fast kernel (the buffers previously carried a stale coordinate into a hardwall branch whose divergence was provably unobservable but pinned by nothing); the PBC 2D LR extractor's empty returns are now typed int32 like all its twins.
- **Doc/dead-code cleanup in the kernels**: the truncated "MOST IMPORTANT COMMENT" and the wrong "x by 6" `idx_to_bead` shape in the reference docstring; the duplicated dead hard-sphere check in `crank_it_good`; unused locals/imports (`num_beads`, `import random`, dead `ok` declarations, unused `num_chains` parameter, 2D `dz`); no-op `cdivision` decorators (and `evaluate_angle_energy_3D` gained the `boundscheck(False)` its 2D twin already had - both its decorators were commented out); the hardwall extractor docstring's nonexistent mode 2; the snakesearch "byte-for-byte" claim scoped to singly-occupied clusters; hetero-only 512-cap documented on the parallel slither; a bounds warning on the unsigned-position `get_gridvalue` helpers.

Six new regression tests pin the fixes (`test_review_fixes.py`): negative-substep no-op across all four crank entry points, the inclusive `randint` contract for every start, collision-freedom of the parallel block seeds, 3D thread-count independence of the checkerboard kernel (previously only 2D was pinned), `snakesearch` seed/empty-input validation, and the absence of `cython_backup` from the tree and manifest.

Everything the audit could not fault is recorded in the audit trail: halo-width sufficiency for every kernel read and write (including across the periodic seam and remainder strips), detailed-balance closure of the interior restriction, the homopolymer O(1) slither shortcut's exactness, RNG consumption-order identity between fast and reference kernels on every path, shell-buffer cardinalities, `cdivision` negative-operand freedom, and the memoryview dtype contracts at every Python call site. Bit-exactness was additionally stress-verified at the minimal box (7) and on non-cubic minimal boxes, neither of which the shipped fixtures previously exercised.

### Documentation audit (docs fully reconciled with the code)

Four parallel doc-vs-code audit passes (moves pages, advanced pages, core pages + keyword reference, lemonade) checked every factual claim in the documentation against the current code, running every code example. 54 findings were fixed. Highlights:

- **Wrong physics descriptions corrected**: the slither page described a Metropolis-Hastings multiplicity correction the kernel does not have (the proposal is symmetric - uniform draw over a fixed-size candidate box - so plain Metropolis is correct, and that is what is implemented); the pull page's Hastings ratio was printed upside-down (`n_fwd/n_rev`, not `n_rev/n_fwd`); the cluster-translate page claimed a fully-connected droplet is diffused bodily (an all-chain cluster is in fact rejected).
- **Stale pre-1.0.8 mechanics replaced**: chain/cluster rotation pages (and the `chain_rotate` docstring) now describe the bead-anchored displacement rotation; the VMMC page distinguishes formed-link vs boundary-link reverse probabilities and documents hardwall rejection; the crankshaft page's proposal description now matches its own detailed-balance argument; the chain-pivot page says the *shorter* arm swings about the (fixed) pivot bead.
- **New 1.0.8 behaviour documented everywhere it was missing**: the PARALLELIZE serial fallback (parallelization page, both megamove pages, and the keyword description), all the new parser validations (QUENCH_FREQ > 0, RESIZED_EQUILIBRATION >= 7, EQUILIBRATION_OFFSET >= 0, TSMMC_FIXED_OFFSET > 0, the heating-quench TSMMC abort), the post-move output cadence and frame-0 convention in the keyword descriptions, the cluster-file leading step columns, the SCALING_INFORMATION per-chain rows and -1/-1 sentinel, the PDB serial wrap, and lemonade's always-on load warnings, slab `shape=None`, lattice-dispersion surface tension and hardwall clustering semantics.
- **Bug fix (the ANALYSIS_MODULE/ANA_CUSTOM parse error was dead code).** The 1.0.8 check ran after `set_dynamic_defaults` had already rewritten a disabled ANA_CUSTOM to "beyond the run length", so it could never fire and a validated module still silently never ran. The check now tests the raw pre-defaults value; regression test added (error raised for 0/unset, clean parse with a real frequency). The converse case also warns: `ANA_CUSTOM` set without an `ANALYSIS_MODULE` prints a parse-time warning that no custom analysis will run, instead of being silently ignored.
- **The API reference now actually documents the core objects.** `Lattice`, `Chain`, `MoveObject` and `TSMMC` had no class-level docstring, so autodoc silently omitted them from the developer reference (their `:class:` cross-references rendered as dead text). All four have proper class docstrings and appear in the built inventory; a stray string literal in `Lattice.__init__` (dead cubic-box-check code inside quotes) that autodoc was picking up as an attribute docstring is removed.
- **Factual corrections across the core pages**: multi-type runs write `CHAIN_<TYPE>_*` internal-scaling files *instead of* (not "also" alongside) the unprefixed ones; `ANA_RESIDUE_PAIRS` is repeatable (not just CHAIN/EXTRA_CHAIN); NON_INTERACTING zeroes only the pairwise/solvation tables (angles need ANGLES_OFF); ANGLE_PENALTY_T_NORM values are floats (the integers-only rule applies to interaction energies and absolute angle penalties); every interacting bead type must carry an angle line unless ANGLES_OFF; restart snapshots are suppressed during equilibration; ENERGY.dat/PERFORMANCE.dat are written throughout equilibration; AUTOCENTER re-centres only *written* frames (single-chain runs only); freeze semantics include the collective-move veto; the docs build instructions now install the real requirements; a stale "required for EXTRA_CHAIN" EXPERIMENTAL_FEATURES line was removed from the star_destroyer demo; the chains-on-a-surface demo keyfile pointed at three files that don't exist (`../`-prefixed parameter, restart and freeze paths) and could not run - it now references the files shipped in the demo directory and the restart file its notebook generates, with the README corrected to match; the FreezeFile docstring example used a chainID (0) that cannot exist.
- `keywords.rst` regenerated from the corrected CONFIG descriptions; the full Sphinx build is warning-free.
- Copyright notices unified to 2015-2026 across the codebase: `scripts/PIMMS` (was 2015-2020), the docs `conf.py` (was 2016-2026) and the archived `initial_docs` README were stale, and 17 shipped modules (the optimised kernel, the whole `lemonade` package and the `fast_kernels` tooling) carried no notice at all.

### File I/O audit (readers + writers, oracle-verified)

A dedicated four-agent audit of every file reader and writer (keyfile/parameter/freeze/restart parsers, restart checkpoints, PDB/XTC trajectory output, and all `.dat` writers), verified against round-trips and the wwPDB spec. 28 findings; the fixes:

- **Restart files are now validated on read and written atomically.** `build_from_file` rejects non-dictionary pickles, beads outside the box (negative coordinates previously wrapped silently onto real lattice sites), and two beads on the same site - each with a `RestartException` naming the offending chain and position. `write_to_file` writes to a temp file and `os.replace`s it into place, so a crash mid-checkpoint can no longer destroy the previous good `restart.pimms`. The pickle security caveat (only load restart files you trust) is now documented.
- **Streamed XTC frames carry real time/step metadata.** The incremental writer stamped every frame with time=0/step=0 (mdtraj's default), while the buffered `SAVE_AT_END` path wrote 0, 1, 2, ...; tools reading the metadata saw all streamed frames as simultaneous. Frames are now stamped sequentially on both paths.
- **AUTOCENTER now applies to `START.pdb` and frame 0.** The topology PDB and the writer's first frame were written without the autocenter shift, so the topology disagreed with every centred frame appended after it.
- **An `ENERGY_CHECK` abort no longer discards a `SAVE_AT_END` trajectory.** The buffered frames are flushed to `traj.xtc` before the exception is raised (the incremental writer was already closed cleanly).
- **PDB output is now column-exact against the wwPDB spec**: ATOM x/y/z are `%8.3f` right-justified in columns 31-54 (values such as `1.5` were previously centre-padded as `  1.5   `, which strict fixed-column readers misparse) and the MODEL serial is right-justified in columns 11-14.
- **Re-runs in a used directory no longer inherit stale outputs.** Startup now also removes all per-chain-type `CHAIN_<T>_*` files (clusters, LR clusters, INTSCAL, INTSCAL_SQUARED, SCALING_INFORMATION, DISTANCE_MAP), `chain_to_chainid.txt` when the keyword is off, and `QUENCH.dat` on non-quench runs - previously a shorter re-run silently mixed two runs' data in anything that globs those files. The one-run-per-directory lifecycle (and the resume-in-a-fresh-directory rule) is now documented in the output and restart pages.
- **One equilibration boundary convention everywhere.** Analysis (and restart snapshots) now treat step == EQUILIBRATION as the last equilibration step (`<=`), matching the XTC saving convention and the PERFORMANCE.dat E/P column, instead of analysing that one step.
- **Parser hardening**: a fully-initialised parser's `write_keyfile` output now round-trips (derived/sentinel keywords - `UNSET` strings, False-valued flags, internal keys - are skipped on write, and the obsolete `CRANKSHAFT_MODE` keyword raises instead of being silently written back); list keywords reject an empty value list; freeze files raise on unrecognised directives with 1-based line numbers instead of silently freezing nothing; `write_acceptance_statistics` and the cluster writers validate their inputs before opening any file, so a failed write can no longer leave a truncated table behind; the angle-summary header now names the parameter file actually used.
- **Console noise control**: per-event rejection messages (clash, cluster resize, multichain rearrangement) and the per-checkpoint restart message respect `REDUCED_PRINTING`; a stray debug `print` of the move count is gone; the logger docstrings say local time because that is what it writes; the run now logs an explicit completion line with the step count.
- Documented: a `SIGKILL`'d run leaves `traj.xtc` valid up to the last complete frame (recoverable frame-wise via `md.formats.XTCTrajectoryFile`; a whole-file `md.load` may refuse the torn tail), and the XTC time/step metadata counts saved frames, not MC steps.

### Independent re-review (second external pass)

A second, from-scratch external review of the 1.0.8 tree (Monte Carlo core, simulation/I/O, analysis) ran after the audits above; its findings were verified and completed here (every shipped restart file re-validated against the new reader checks, the full suite and fixture baselines regenerated). The fixes:

- **Bug fix (`INTSCAL_SQUARED.dat` was the square of a mean, not a mean of squares).** The squared internal-scaling accumulator was fed the per-snapshot *pair-averaged* distance and squared it, giving `<r>^2` per snapshot instead of `<r^2>` and discarding the within-snapshot variance at every sequence gap. The chain analysis now computes both moments in one vectorised pass (`get_internal_scaling_profile(..., return_squared=True)`) and feeds the second moment to a dedicated `update_internal_scaling_squared`. `sqrt(INTSCAL_SQUARED)` is now the RMS internal-scaling profile it was always documented as.
- **Internal-scaling profiles now span every gap 1 .. L-1.** The accumulators and profile stopped one gap short (the end-to-end separation was excluded), so `INTSCAL*.dat` gain a final row and the scaling-exponent fit can use the most informative point.
- **`DISTANCE_MAP.dat` is now the full symmetric matrix.** Only the upper triangle was ever populated, so half the residue pairs looked coincident when the file was loaded or `imshow`n directly.
- **Output-format heads-up**: because of the two items above, `INTSCAL.dat`/`INTSCAL_SQUARED.dat` (and their `CHAIN_<T>_` variants) carry one more row and `DISTANCE_MAP.dat` has a populated lower triangle. Downstream scripts that hard-code the old row count or read only the upper triangle need a look; `INTSCAL_SQUARED` values also change in meaning (see above). All 15 fixture baselines were regenerated, and every changed trajectory (as opposed to output format) was attributed to the SLITHER/PULL monomer change below - test_3/8/9/12 are exactly the monomer-containing fixtures with SLITHER enabled.
- **Disabled analyses stay disabled.** A frequency below 1 was represented only as `N_STEPS + 10`; a resized equilibration that extends the run past that point silently re-activated the "disabled" analysis (and `ENERGY_CHECK`). The parser now records the disabled set explicitly (`__DISABLED_FREQUENCIES`), the simulation skips those routines and no longer writes plausible-looking zero/`-1` final `INTSCAL`/`SCALING_INFORMATION`/`DISTANCE_MAP` files for them, and `write_keyfile` writes a disabled frequency back as `0` so it round-trips as disabled rather than as an enabled analysis at the sentinel cadence.
- **`SLITHER`/`PULL` are no longer remapped to a crankshaft when the outer-loop seed chain is a monomer.** Both are whole-system megamoves whose per-chain eligibility is decided inside the kernel, so the remap made their effective frequency composition-dependent and starved the polymers in mixed monomer/polymer systems. Per-chain rotations and pivots keep the remap (and now also reject singleton chains directly without touching the lattice).
- **VMMC numerics/perf**: link probabilities go through a stable `-expm1` form (the old `exp`-then-clamp overflowed for large negative energy differences), the harmonic cluster-cutoff CDF and Chebyshev neighbour shells are cached (bisection instead of a per-attempt linear scan), and the parallel-kernel gate ignores permanently frozen over-cap chains so a frozen 513-mer no longer forces every movable chain onto the serial path.
- **Analysis API details**: `extract_cluster_polymeric_properties` no longer re-applies minimum-image wrapping when explicit box `dimensions` are passed - its input contract is already-single-image cluster coordinates, and re-wrapping could fold an extended but non-percolating cluster (the production cluster analysis never passed `dimensions`, so `CLUSTER_RG`/`ASPH` outputs are unaffected); `InternalScalingSquared.write_status` defaults to `INTSCAL_SQUARED.dat` instead of overwriting `INTSCAL.dat`; the running-mean accumulators (internal scaling, distance map) use the numerically stable in-place form `mean += (x - mean)/(n+1)`, which also avoids two full temporary matrices per distance-map update.
- **Restart reader hardening**: every field is validated before the object is touched (a failed read cannot leave a previously usable `RestartObject` half-overwritten): dimensions must be a 2D/3D positive-integer list, `ENERGY` finite, `HARDWALL` boolean, chainIDs positive integers representable in the occupancy grid (0 is the solvent sentinel), sequences non-empty strings, coordinates integers, and consecutive beads must be lattice neighbours under the stored boundary mode (a disconnected chain is rejected instead of corrupting the mover). Resize offsets and dimensions are validated the same way; a failed atomic write also removes its `.tmp` file. The PBC-restart-with-`RESIZED_EQUILIBRATION` check now also applies under `RESTART_OVERRIDE_HARDWALL`, which used to bypass it.
- **Parser/inputs**: numeric keywords reject `nan`/`inf` (which defeat every range check); freeze-file directives must carry at least one ID and frozen sets are sorted for deterministic logs; non-finite temperatures are rejected by the acceptance calculator; the `frozen_chains`/`chain_list` mutable-list defaults on the move and event signatures are now immutable tuples.
- **Custom analysis modules**: the module's directory is placed *first* on `sys.path` only while it loads and while `analysis_function` runs (a same-named module earlier on the path used to shadow a local helper, and the directory used to stay on the path permanently); module names include a path digest so `a-b/analysis.py` and `a_b/analysis.py` no longer collide; a `SystemExit` at import is reported as a keyfile error; failed loads unregister their module.
- **Bug fix (`pimms --version` crashed with `AttributeError: module 'pimms' has no attribute '__version__'`).** A stale `pimms/` directory left in site-packages by an old non-editable install (one orphaned backup file, owned by no distribution) was picked up by Python's path finder ahead of the editable-install finder and imported as an `__init__`-less *namespace* package: every submodule still resolved to the repo, so simulations ran, but nothing defined in `pimms/__init__.py` existed. The CLI now detects a namespace-package `pimms` and exits with a message naming the shadowing directory, and reads the version through a helper that falls back to the distribution metadata. New `test_cli.py` exercises every command-line mode (`--version`, `--help`, `--info`/`--info <KW>`/`--info ALL`/unknown keyword, no arguments, missing keyfile, and a full `-k` run) from a *foreign* working directory with no `PYTHONPATH`, checks that `import pimms` from outside the repo is the real package, and runs the installed `PIMMS`/`pimms` console scripts when present - the previous CLI test pointed `PYTHONPATH` at the repo and so could not see this.
- **CLI/packaging**: `PIMMS --version`/`--help` no longer import the simulation stack (NumPy/MDTraj/compiled kernels); `MANIFEST.in` globally excludes generated simulation outputs (`.dat`/`.xtc`/`.pdb`/logs/restarts left in the fixture directories by local test runs used to be swept into the sdist by `graft pimms`) while re-including the one intentional data file; the trajectory writer is always closed via a `run_simulation` guard even when a move or analysis raises; the running-mean accumulators use the in-place stable form.
- Regression tests: `test_rereview_mc_core.py`, `test_rereview_cli.py`, `test_package_hygiene.py`, plus new cases in `test_custom_analysis.py`, `test_data_structures.py` and `test_review_fixes.py` (disabled-frequency round trip, no final files for disabled analyses, mean-of-squares profile).

### PARALLELIZE: same equilibrium, slower relaxation - halo tightened, claims corrected

A user report that `PARALLELIZE : True` and `False` reach "quite different final energies" was run to ground. It is **not** an energy-evaluation or equilibrium bug: the parallel crankshaft's tracked energy is exactly equal to a from-scratch recompute in every multi-block regime (3D SR/LR/SLR, hardwall and periodic, 2D), and high-statistics detailed-balance comparisons in the multi-block regime give the same equilibrium energy as the serial kernel (e.g. 44^3 SLR periodic: −23187 vs −23203, z = −0.3). What differs is the **dynamics**: each sweep only the beads inside the block interiors can move and proposals that would land in a frozen halo are rejected, so the parallel Markov chain relaxes more slowly per step. From the same starting configuration, at a collapse temperature on a 44^3 long-range box, the serial kernel reached −32480 after 150 sweeps and the parallel kernel only −27747 with identical attempt counts; through the full `Simulation` path at default settings the gap after 300 steps was ~1600 energy units and growing. A run judged by its energy at a fixed step count therefore differs; the equilibrium does not. Changes:

- **Crankshaft halo tightened to the minimal race-free width.** The halo `W = R_int + 2` (3 without / 5 with long-range beads) predates the closed-interior rule (moves may never land in the halo); with that rule every write is >= W inside a block face and every read reaches at most `max(R_int, 2)` past the interior edge, while the next interior starts `2W` past it, so `W = 2` with long-range interactions and `W = 1` without is sufficient. Blocks are also kept >= 8W long. On a 44^3 long-range box the movable fraction per sweep rises from 16% to ~55%; on a 30^3 short-range box from 22% to ~51%. The whole-chain slither/pull kernels keep their chain-level `R_int + 2` halo. Results remain independent of the thread count.
- **Documentation corrected.** "Changes only speed, never sampling" overstated it: the equilibrium is unchanged but the relaxation rate is not. The parallelization page, the `PARALLELIZE` keyword description and the kernel comments now say so explicitly, explain why energies at a fixed step count differ for un-equilibrated runs, and point to `mega_crank_fast.parallel_crank_layout_info(X, Y, Z, has_LR)`, which reports the halo, block layout and movable fraction for a box.
- **Startup parallelization report.** With `PARALLELIZE` on, startup now prints and logs exactly which implementation this box/system gets: thread budget and whether the compiled kernels have OpenMP at all (`mega_crank_fast.openmp_info()`; a build without OpenMP is now stated plainly rather than silently running serially), the crankshaft halo/block grid/block size and the fraction of the box movable per sweep (or the single-block warning), and for slither and pull whether every chain fits a block interior - i.e. parallel kernel or serial fallback, with the reason (`moves.parallel_chain_fit_report`). The report is re-issued after a resized equilibration changes the box, and is covered by `test_parallel_report.py`.
- **Test coverage was vacuous for exactly this regime.** The 3D parallel detailed-balance and energy-consistency tests used a 30^3 box, which under the long-range halo was a *single block*, so the halo logic they describe was never exercised for LR systems. They now use a 40^3 box and assert the layout has more than one block; new tests pin bookkeeping exactness in every multi-block regime and the halo width / movable-fraction floor.

### Documented decisions (behaviour deliberately kept)

- The angle-penalty classes (A1/A2/A3) are keyed to the displacement pattern of the i-1 -> i+1 vector, not to the geometric bend angle - some 70-110 degree bends share the A3 class with the straight-through geometry. Changing this would change the physics of every ANGLE simulation, so the code comments and keyfile docs now describe the criterion actually implemented instead of claiming "straight line" / "distinct bend angles".
- The 2D asphericity remains kappa (the 3D value is kappa squared); frame 0 of `traj.xtc` remains the starting configuration regardless of `SAVE_EQ`. Both are now documented.
- `EmptyHamiltonian` gained the missing `evaluate_local_energy_SLR` and the correct three-argument `evaluate_angle_energy` signature, making it the drop-in its docstring always claimed.

## 1.0.7 (August 2026)

Small bug fixes and docs update, including a few additional tests, plus a clean-up of everything the package build log raises.

- **Bug fix (`pytest pimms/` executed a 10-million-iteration loop during collection).** The same bug class fixed in 1.0.5 (`test_megagrank.py`) had two more instances: `pimms/randint_test.py` and `pimms/randneg_test.py` match pytest's default `*_test.py` collection pattern, and their top-level bodies print 10 million / 1 million random numbers. `testpaths` protects a bare `pytest`, but an explicit `pytest pimms/` imported and ran them at collection time (confirmed: collection hangs printing random integers). Worse, all of these dev scripts ship in the wheel as importable modules whose *import* runs the loops. Renamed to `dev_randint_check.py` / `dev_randneg_check.py`, and every developer script in the package (`check_randomness.py`, `cython_testing.py`, `print_interaction_matrix.py`, `dev_megacrank_rng_check.py` included) is now `__main__`-guarded, so importing any shipped module is side-effect free. Two new package-hygiene tests pin this permanently: one imports every shipped module in a fresh interpreter and asserts nothing is printed and nothing raises, the other asserts no file outside the test packages matches a pytest collection pattern.
- **Bug fix (the NeurofilamentDemo worked example could not run at all).** Many moons ago (2015?) we developed a demo with a neurafilament system. This became dead; `initialized_systems.NeurofilamentDemo` - the documented example of programmatic system construction - had rotted against the evolving constructors: it called `Chain` without the `dimensions` / `LR_int_seq` / `LR_IDX` / `chainType` arguments and `Lattice` without the required `lattice_to_angstroms`, so instantiating it raised `TypeError` before building anything. The constructor calls are now current; the geometry is parametrized (box size, sidearm length/spacing) so the demo can be exercised at test size instead of only at the hardcoded 500^3; and the site scan walks only the filament region rather than all 125 million lattice sites. Two tests build the system small and verify grid/chain consistency.
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
- **Jump-and-relax (`MOVE_JUMP_AND_RELAX`, code 13).** Previously expermental move that is now decomposed into three π-preserving sub-steps (relax → Metropolis-accepted jump → relax) whose composition is detailed-balanced. The 8 regression baselines using this move were regenerated.
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
