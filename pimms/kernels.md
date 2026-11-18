# The compiled kernels

A map of PIMMS's Cython: what each module does, the data it works on, where its random numbers come from, and the rules every kernel follows. Per-function detail is in the `.pyx` docstrings. Sphinx mocks these modules, so they're not in the API reference. Paths are relative to `pimms/`.

## The modules

| Module | What it holds | Called from | Threads |
|---|---|---|---|
| `mega_crank_fast.pyx` | The production megamove kernels (serial and parallel crankshaft, slither and pull, 2D and 3D), their energy and angle helpers, and their random-number generators | `moves.py` (`system_shake`, `system_slither`, `system_pull`, the TSMMC and single-chain relaxations); its layout and OpenMP helpers also feed the start-up parallelization report | OpenMP `prange` in the parallel kernels |
| `mega_crank.pyx`, `mega_crank_2D.pyx` | The original crankshaft kernels, kept only as the bit-exactness oracle for the fast kernel | tests (`kernel_test_utils.py`, `test_kernel_correctness.py`, `test_proposal_symmetry.py`, `test_review_fixes.py`, `test_deep_audit_kernels.py`), `fast_kernels/` benchmarks, `dev_*_check.py`. Never a production move | none |
| `bookkeeping.pyx` | The two copies between the Chain objects and the bead table | `crankshaft_list_functions.update_idx_to_bead` / `write_back_positions` | none |
| `inner_loops.pyx`, `inner_loops_hardwall.pyx` | Pair-list extractors: every short-, long- and super-long-range neighbour of a site | `lattice_utils.build_envelope_pairs` / `build_all_envelope_pairs`, `longrange_utils.build_LR_envelope_pairs` | none |
| `hyperloop.pyx` | Sums a pair list against an interaction table; a chain's angle energy; small site helpers | `energy.Hamiltonian`; the site helpers from `lattice_utils` and `cluster_utils` | none |
| `cluster_kernels.pyx` | The single-image gather for cluster analysis | `cluster_utils.convert_positions_to_single_image_snakesearch` | none |
| `system_utils.pyx` | Start-up checks on integer types and chain count | `simulation.py` | none |
| `lemonade/kernels/_pbc.pyx` | lemonade's chain unwrap and grid paint | `lemonade/_store.py` | none |
| `cython_config.pxd` | The two integer typedefs every kernel shares | cimported by all of the above | |

Only `mega_crank_fast` gets OpenMP flags (`openmp_flags` in `setup.py`). On macOS without Homebrew's libomp it's built without OpenMP and `prange` runs serially, with identical results. Elsewhere `-fopenmp` is always passed, so a compiler without OpenMP fails the build. `mega_crank_fast.openmp_info()` says which you've got.

## The data

**Integer types.** `cython_config.pxd` defines `NUMPY_INT_TYPE` (int32, twinned by `CONFIG.NP_INT_TYPE`; `system_utils.check_dtype_consistency` checks they agree at start-up) and `NUMPY_INT_TYPE_long` (int64). Grids and interaction and angle tables are int32; the bead table and bead selector are int64. Energies are integers throughout, accumulated in C `long`, which is exact to 2^63 on macOS and Linux. Where `long` is 32 bits (Windows), the energies in `mega_crank_fast`, `mega_crank`, `mega_crank_2D` and `hyperloop` would silently wrap at 2^31 and the parallel kernels wouldn't run; a port has to move them all to `NUMPY_INT_TYPE_long` first. Nothing checks this at build time.

**The grids.** `latticeObject.grid` holds the chainID at each site, and `latticeObject.type_grid` the bead's integer residue code (its intcode). Both use 0 for empty and are C-contiguous int32 arrays of the box shape. Kernels read and write both in place and never see a Chain object. `system_utils.check_beads_to_grid_mapping` refuses more chains than int32 can label.

**The bead table.** `latticeObject.crankshaft_lists` is int64, one row per bead, ordered by chainID and then position along the chain, so rows `i - 1` and `i + 1` are a bead's bonded neighbours.

| column | content |
|---|---|
| 0 | bead flag (below) |
| 1 | long-range flag: 1 if the bead has long-range interactions |
| 2 | intcode, which indexes the interaction and angle tables |
| 3 | skip-angle flag: 1 for beads of 1- and 2-bead chains |
| 4 | chainID |
| 5, 6, 7 | x, y, z (2D tables have 7 columns and no z) |

The bead flag sets which proposal a kernel makes and which angle triplets a move changes:

| flag | bead | angle window (rows relative to the bead) |
|---|---|---|
| 0 | a single-bead chain | none |
| 1 | first bead | 0 to +2 |
| 3 | last bead | -2 to 0 |
| 4 | middle bead of a 3-mer | -1 to +1 |
| 5 | second bead of a chain of four or more | -1 to +2 |
| 6 | second-to-last bead of a chain of four or more | -2 to +1 |
| 2 | any other interior bead | -2 to +2 |

Columns 0 to 4 never change during a run. `ChainLayout` (`crankshaft_list_functions.py`) holds the matching per-chain layout: sorted chainIDs, each chain's first row and length, the bead count, and a homopolymer flag (1 when a chain's intcode and long-range columns are each constant). It's built in `Lattice.__init__` and cached on the lattice. `chain_layout()` rebuilds it if the chain or row count changes, which only tests do: nothing in PIMMS adds, removes or rewrites a chain after construction (`EXTRA_CHAIN` is applied at construction).

**The round trip.** Every megamove:

1. `update_idx_to_bead` gathers positions from the Chain objects into the lattice's table (`bookkeeping.gather_positions`) and hands the kernel a fresh copy.
2. The kernel updates the copy's positions and both grids.
3. `write_back_positions` scatters the copy back into the Chain objects (`bookkeeping.scatter_positions`, through `Chain.set_ordered_positions`, so the length check still runs).

Between steps 1 and 3 only the kernel's copy is current. Both compiled copies raise `ValueError` if the table and the chains are out of step. The Python loops they replaced (`_update_idx_to_bead_python`, `_write_back_positions_python`) are kept as the fallback and the oracle; `test_bookkeeping.py` runs megamoves both ways. The TSMMC and single-chain relaxations use `update_idx_to_bead_single_chain` and `_multiple_chains` instead, which take just the rows they need (via `chain_to_firstbead_lookup`) and write positions back themselves.

**Interaction and angle tables.** The Hamiltonian provides one `(n_intcodes, n_intcodes)` int32 table per shell (short-, long- and super-long-range), with row and column 0 standing for an empty site. The angle table is indexed by the middle bead's intcode and the two bond vectors of a triplet, each component shifted from -1..1 to 0..2 (seven dimensions in 3D, five in 2D). `ANGLES_OFF` zeroes it; the lookups still happen.

## Random numbers

There are three generators. It matters which one a draw comes from, because the restart file saves the first two and the third is reseeded from the first on every call.

1. **Python's `random`**, seeded from `SEED`. Every megamove draws `local_seed = random.randint(1, sys.maxsize - 1)` from it for the kernel. Move selection and the Python-level moves use it too.
2. **numpy's global generator**, seeded from `SEED` modulo 2^32. It draws the crankshaft's bead selector (`bead_selector_constructor`: `np.random.randint`, or `np.random.choice` over the movable rows when any bead is frozen), the relaxations' bead selectors and the multi-chain TSMMC's chains, and it shuffles the slither and pull chain selectors.
3. **The kernels' own streams**, all from `passed_seed`:
   - **Serial kernels:** one module-global splitmix64 state, reset by `mc_seed(passed_seed)` at the start of every call, so a call is a pure function of its inputs. `mc_rand` returns the top 31 bits of a step. `randint(start, end)` is inclusive and takes one step. `accept_or_reject` takes a step only for an uphill move (plus, rarely, the second stage below).
   - **Parallel kernels:** one splitmix64 state per block, seeded by mixing the block index and `passed_seed` through a splitmix64 finaliser, and a 53-bit uniform. Don't replace the mixing with an affine formula: under one, (block, seed) pairs with the same sum share a stream.
   - **Block shift:** a legacy numpy `RandomState` (`_SHIFT_RNG`), which `_block_shift_rng` reseeds with `passed_seed & 0x7FFFFFFF` each call. It's reseeded rather than rebuilt because building one costs several hundred microseconds.

**Draw counts are part of the contract.** `mega_crank.pyx` has a byte-for-byte copy of the serial generator and of `randint`, evaluated in double precision, which makes the fast serial crankshaft bit-exact with it. Every serial substep takes a fixed number of draws before the acceptance test, even if the proposal then fails on a clash or a wall:

- 3D crankshaft: three (two in 2D);
- 3D slither: one for the direction plus three target offsets (just the offsets for a single bead).

One more follows only for an uphill move. Pull differs: it draws a bead and a direction, then a first target only if one exists. It makes an acceptance draw only if the cascade reconnects the chain, a reverse first target exists and `(nF / nR) exp(-beta dE) < 1`. Don't reorder or add draws without regenerating the regression fixtures.

**The serial acceptance test.** `accept_or_reject` (crankshaft and slither) compares `float32(r / (2^31 - 1))`, with `r` the 31-bit draw, against the float32 Boltzmann factor. On its own that can't resolve probabilities below 2^-31: `r == 0` passes any positive factor, so every uphill move was once accepted at least 4.66e-10 of the time.

A second stage, `mc_rand_below`, fixes this. It runs only when the first comparison accepts and the factor lies inside the draw's own bin `[r, r + 1) / 2^31`. It then compares a uniform against the fraction of that bin below the factor, one base-2^32 digit per generator step. The result:

- the accepted probability is exactly the float32 factor;
- a factor that underflows to zero is never accepted;
- above a factor of about 2^-7 the second stage is never reached, and the accepted probability is within one float32 rounding (a relative 6e-8) below the factor;
- every decision outside the bin is unchanged, so the regression fixtures didn't move.

`accept_or_reject_ratio` (pull) wraps the same second stage around `acc = (nF / nR) exp(-beta dE)`, in double precision. It accepts with probability `acc` exactly, except when `0 < frac(acc x 2^31) <= acc`, where it falls short by at most a relative 2^-31 and is never above.

The reference kernel has byte-for-byte copies of both stages. The probes (`accept_or_reject_probe` in both modules, `accept_or_reject_ratio_probe`) run one test from a chosen generator state for `test_deep_audit_kernels.py`; production never calls them. The parallel tests (`accept_p`, `accept_p_ratio`) use a 53-bit uniform, floored at 2^-53, with no second stage.

## Energy conventions

**Shells.** Interactions are counted over Chebyshev shells around a site:

- short-range (SR): the 3x3x3 block minus the centre;
- long-range (LR): the 5x5x5 minus the 3x3x3;
- super-long-range (SLR): the 7x7x7 minus the 5x5x5.

A bead whose long-range flag is 0 scores only the SR shell. The same radii define long-range clusters in `lattice_utils`, the cluster gather's links and the VMMC neighbour shell.

**Solvent.** Empty sites are real partners in the SR shell: a bead next to one scores `table[type, 0]`, which is how solvation enters. LR and SLR have no solvent term; the extractors skip empty sites there, and the megamove kernels rely on column 0 of the LR and SLR tables being zero.

**Bonds and self-terms.** A bond is always a Chebyshev-1 contact, so bonded pairs are ordinary SR pairs that add a constant to the total and cancel in every delta. The megamove delta visits the centre site of both shells, then subtracts the bead's self-interaction and self-solvent terms.

**Angles.** For each triplet, the two bond vectors are taken from the middle bead, and each component is folded back into -1..1 (a bond across the periodic boundary shows up as `L - 1`). The penalty is `angle_lookup[middle intcode, a + 1, b + 1]`. `hyperloop.evaluate_angle_energy_*`, the reference kernel and `mega_crank_fast` all use this fold.

**Hard walls.** Two methods are in use, and each is only correct where it's used now:

- The `mega_crank_fast` deltas and the `hyperloop` summers treat a pair whose partner is reachable only through the periodic wrap (a per-axis separation above 3) as a wall. In the SR shell the partner reads as empty, so a wall acts like solvent; in LR and SLR the pair is skipped.
- The `inner_loops_hardwall` extractors never emit a pair that leaves the box.

So under a hard wall, SR pair lists for energies must come from the periodic extractor, and they do. `evaluate_total_energy` and `commit_single_chain_move` build their SR lists without the hardwall flag; the hardwall extractors serve only LR and SLR lists and connectivity. The separation-above-3 test needs every box axis to be at least 7 sites, which the keyfile parser and `Lattice()` (`MIN_LATTICE_DIMENSION`) enforce; the kernels don't check. VMMC's neighbour scan skips out-of-box sites in every shell, which is fine because link energies have no solvent term.

**Which move uses which energy path.**

| Move | Energy path |
|---|---|
| crankshaft, slither, pull, the TSMMC and single-chain relaxations | `mega_crank_fast` |
| chain translate, rotate, pivot, head pivot; segment regrowth | `moves.commit_single_chain_move` (via `Simulation.single_chain_move` for the first four): `inner_loops` pair lists summed by `hyperloop`, comparing the old and new regions in each shell, plus the angle change |
| jump-and-relax | two crankshaft-kernel relaxations either side of a chain translation scored by `commit_single_chain_move` |
| cluster translate, rotate | `Simulation.rigid_cluster_move`: LR and SLR only (a rigid cluster's SR neighbours are solvent before and after), and 0 straight away when no residue has long-range interactions |
| VMMC | recruitment link energies: a Python loop or a numpy scan of the offset shell, whichever is cheaper, cross-chain terms only. The move's energy change: `evaluate_total_energy` after the cluster is committed |
| from-scratch total, test oracles | `Hamiltonian.evaluate_total_energy` over whole-system `inner_loops` pair lists summed by `hyperloop`; `tests/test_energy_oracle.py` checks it term by term against an independent numpy implementation |

## The megamove kernels (`mega_crank_fast.pyx`)

### The energy delta

`get_energy_change_c` (and `_2D_c`) is what every substep pays for. Given a bead's old and new sites, its long-range flag and the tables, it:

1. precomputes wrapped coordinates for offsets -3..3 on each axis;
2. sums the old shell;
3. moves the bead in `type_grid` and sums the new shell;
4. subtracts the self-terms and puts the bead back;
5. returns new minus old, solvent changes included.

Callers see `type_grid` unchanged, and the new site must already be empty in `grid`.

It dispatches to two bodies, `_energy_change_periodic_c` and `_energy_change_wall_c`, and the periodic one has no hardwall code. That's for speed: with an `if hardwall` inside the 7x7x7 loops the compiler can't unroll them, which made every LR and SLR evaluation 2.5 to 3 times slower, and a profiler only showed time "in the energy function". Give any new per-site condition its own body.

`get_angle_energy_change_c` returns 0 for a skip-angle bead. Otherwise it scores the flag's angle window before and after the move.

### Serial crankshaft: `mega_crank`, `mega_crank_2D`

Inputs: both grids, the bead table copy, the interaction and angle tables, the current energy, `invtemp`, the substep count, the bead selector, the seed and the hardwall flag. Returns `(energy, accepted)` and updates the grids and table in place.

Each substep takes the next bead from the selector and proposes by flag:

- **single bead (flag 0):** three offsets in -1..1 around its own site;
- **end bead (flags 1, 3):** three offsets around its one bonded neighbour, so it can land anywhere in that neighbour's 3x3x3 cube;
- **any other bead:** a uniform site in the overlap of its two neighbours' cubes (neighbours that straddle the wrap are first placed side by side in an extended box).

An occupied target is rejected, and so, under a hard wall, is a bond (or a single bead's step) across the wrap; either way the draws are already spent. Otherwise the interaction plus angle delta gets a Metropolis test.

`system_shake` draws the bead selector in Python (`bead_selector_constructor`, numpy's generator) so that bead choice and perturbation use different streams. Drawing both from one stream once correlated the bead with the perturbed axis.

### Parallel crankshaft: `mega_crank_parallel`, `mega_crank_parallel_2D`

**Layout.** The box is cut into at most four blocks per axis, each at least `8W` sites long. `W` is the halo width: 2 if any bead has the long-range flag, 1 otherwise. It's half the gap between two moving regions, not the interaction range. Neighbouring interiors are `2W` apart, writes stay `W` inside a block face, and the widest read reaches `max(R_int, 2)` past the interior edge (the SLR shell, the occupancy check, the angle window). So the condition is `2W >= max(R_int, 2)`. `parallel_crank_layout_info` reports the layout and the movable fraction, and `system_shake` uses the parallel kernel only when there's more than one block.

**A sweep.**

1. Draw an origin shift on each split axis, uniform over the whole box.
2. Put each bead in a block's movable list or mark it immovable (in a halo, in the leftover strip, or frozen).
3. Share the `nsteps` attempts among blocks in proportion to their movable beads, so exactly `nsteps` are made (none if nothing is movable, which returns at once).
4. Run `run_block` over the blocks with a dynamic `prange`. Each attempt picks a bead from the block, proposes as the serial kernel does from the block's own stream, rejects a site outside the interior, and accepts on the local delta alone (`accept_p`).
5. Sum the block deltas and counts.

The parallel kernels return `(energy, accepted, attempted)`, and the wrappers log `attempted`, so an empty sweep isn't logged as `nsteps` rejections.

**Frozen beads.** `frozen_mask` has one entry per bead, from `crankshaft_list_functions.frozen_bead_rows`: all beads of a fully frozen chain plus the listed beads of a partly frozen one. Masked beads are never picked but stay in the grids as obstacles.

**Why it's correct.** Within a sweep the movable set is fixed, proposals are symmetric, the interior is closed (a move into the halo couldn't be reversed, so it's rejected), and acceptance is Metropolis on the exact delta. So each sweep preserves the Boltzmann distribution, and the shifting halos let every bead move over a run. Draw the shift over the whole box, not one block length: on an axis with a leftover strip, a one-block draw made how often a site was movable depend on where it sat (0.50 to 0.75 along a 26-site axis). With the whole-box draw it's `nb (L - 2W) / DIM` everywhere, the `movable_fraction` reported.

**Why it relaxes more slowly.** Only interior beads move in a sweep, so compare equilibrium averages between serial and parallel runs, never energies at a fixed step.

**Threads.** The blocks depend only on the box and `W`, so results don't depend on the thread count. `_clamp_threads` (in all six parallel kernels) limits the request to `[1, num_blocks]`, since an unbounded one can make the OpenMP runtime abort.

### Slither: `mega_slither`, `mega_slither_2D` and the parallel twins

In place of the crankshaft's substep count and bead selector, slither (and pull) take each chain's offset, length and homopolymer flag from the layout; a chain selector (each eligible chain repeated `SLITHER_SUBSTEPS` times and shuffled); and `max_chain_len`, which sizes the buffers used to undo a rejection.

**One attempt.** A single-bead chain gets the crankshaft's single-bead proposal. Otherwise draw a direction (head-grow moves the last bead forward, tail-grow the first bead backward) and a target in the 3x3x3 cube around the growing end. Reject the target if it's occupied (even by the other end, which is about to leave) or, under a hard wall, reachable only through the wrap. Then score it:

- **Homopolymer:** the occupied sites change by one out and one in, all of one type, so the interaction change is a single-bead delta from the vacated end to the new site. The angle change comes from `_chain_angle_mode`, without touching the table.
- **Heteropolymer:** every bead's type moves one site along, so the move is done as `L` single-bead moves in sequence, each into the site the last bead left, each scored against the live grids. The deltas sum to the exact total, and a rejection restores the chain from the buffers.

**The parallel kernel.** A slither moves a whole chain, so:

- a chain moves in a block only if all its beads are in that block's interior after the shift;
- the halo is chain-level, `W = R_int + 2` (5 with long-range interactions, 3 without), and blocks are at least `4W` long;
- the shift is uniform over the whole box, so a chain as long as the interior fits from any starting coordinate;
- a new end in the halo is rejected, so chains never leave their blocks and blocks run without locks.

`frozen_mask` marks the beads of every chain the parallel kernel mustn't touch (frozen chains and chains sent to the serial pass). The chain selector is read only for its length, the attempt count. Each attempt picks a chain from its block uniformly with replacement, so a chain isn't tried exactly `SLITHER_SUBSTEPS` times. Heteropolymers longer than the 512-bead per-thread buffer are skipped in the kernel, and the Python side keeps them out anyway.

**The length partition.** `system_slither` and `system_pull` split chains between the two kernels by length (`parallel_chain_partition`). Chains with any frozen bead go to neither. The rest go parallel if they fit in the smallest block interior and are under the buffer cap (for slither, the cap applies only to heteropolymers); a box with one block sends everything serial. The split only changes when a resized equilibration grows the box.

Never choose the kernel from the current configuration. The parallel kernel can't take a chain out of the "fits in a block" region, while the serial kernel can, so choosing by current extents pumped probability into compact conformations and silently biased the mean squared radius of gyration low by a few percent. `parallel_chain_fit_report` is for the start-up summary only.

Both passes run every megamove, chained through the running energy, in an order set by a coin flip. Each pass is reversible but they don't commute, so a fixed order wouldn't be reversible, and a system TSMMC excursion needs every move it runs to be.

### Pull: `mega_pull`, `mega_pull_2D` and the parallel twins

Chains shorter than three beads are skipped. One attempt:

1. Draw an interior bead `i` in `[1, L - 2]` and a direction, which sets the anchor (`i - 1` or `i + 1`).
2. The first targets are the `nF` empty sites adjacent to both `i` and the anchor. If there are none, stop.
3. Move `i` to one of them.
4. Cascade away from the anchor: each next bead takes the previous bead's old site, until the next bead is already adjacent to the last one moved. If that never happens before the chain end, reject.
5. Count the reverse first targets `nR` the same way in the final state, from the last moved bead in the opposite direction.
6. Accept with `(nF / nR) exp(-beta dE)`.

The bead and direction probabilities cancel and the cascade is fixed once the target is chosen, so `nF / nR` is the whole Hastings correction. A rejection clears every current site before rewriting the old ones, so overlapping old and new sites are safe.

The parallel kernel uses the slither's block scheme, but its 512-bead cap applies to every chain. Its one algorithmic difference is `pull_first_targets_interior`, which keeps first targets inside the block interior and counts `nR` the same way. That's required: the reverse move really is drawn from the restricted set, so both counts must be restricted.

### The 2D twins

Every kernel has a 2D version, line for line the same logic with two coordinates, a 7-column table, a five-dimensional angle table and 7x7 and 3x3 rings. It takes one fewer draw per crankshaft or slither proposal (pull takes the same three), and the parallel kernels have no z shift.

## The reference kernels (`mega_crank.pyx`, `mega_crank_2D.pyx`)

The original crankshaft kernels, slow because they allocate numpy arrays inside the substep loop. They're still built because `test_kernel_correctness.py` runs them and the fast kernel from the same selector and seed and requires identical results, and `test_proposal_symmetry.py` uses them to check proposals are uniform. The 2D module borrows the 3D module's generator through Python. Their seed is a C `int`, so bit-exactness can only be tested with seeds below 2^31.

## Pair lists and summers (`inner_loops*.pyx`, `hyperloop.pyx`)

The Python-level moves and the from-scratch total build explicit pair lists. `inner_loops` returns, per site, int32 arrays of shape `(n, 2, ndim)`: the 26 SR pairs (8 in 2D), empty partners included, and, for a bead with the long-range flag, the LR and SLR pairs to occupied sites. A pair's two sites are ordered by the sign of the first non-zero component of their unwrapped offset, so the same pair comes out identically from either end and duplicates can be dropped by exact row match (`lattice_utils._unique_rows`). `inner_loops_hardwall` has the same six functions but drops partners outside the box.

`hyperloop` sums a pair list against a table (2D or 3D, SR or LR/SLR) using the hardwall method above, and computes a chain's angle energy. `evaluate_total_energy` uses it for the whole system, `commit_single_chain_move` for the old and new regions, and `rigid_cluster_move` for LR and SLR only. Its adjacent-site lists feed `lattice_utils.get_empty_site` and the head pivot.

## The cluster gather (`cluster_kernels.pyx`)

`snakesearch_single_image` returns a cluster's beads in one periodic image. Starting from the bead nearest the periodic centre of mass (chosen in Python), it walks breadth first, placing each bead at its parent's image plus their minimum-image separation.

- **Links** follow the rule `lattice_utils.get_all_chains_in_long_range_cluster` builds clusters with: Chebyshev 1, Chebyshev 2 with a non-zero LR entry, or Chebyshev 3 with a non-zero SLR entry. Walking a plain distance instead tore long clusters across a periodic face.
- **Site lookup** is a flat grid or a binary search over sorted site keys, picked by `snakesearch_uses_grid` from a measured crossover. The walk is written out four times (2D/3D, grid/search) because a branch in it cost 8 to 16 percent.
- **Output** is shifted by whole box lengths so no coordinate is negative. A cluster that doesn't connect raises `ValueError`. Whether it wraps the box is decided by `cluster_utils.cluster_percolates`, with the same rule.

The pure-Python gather in `cluster_utils` gives byte-identical results and is the fallback and the oracle. The walk holds the GIL.

## lemonade's kernels (`lemonade/kernels/_pbc.pyx`)

`unwrap_chains` makes every chain whole in every frame by walking its bonds, shifting each bead by whole box lengths until it's within one site of the previous one. Shifts build up along the chain, so chains longer than the box work and coordinates can leave the box. It's the batched `lattice_utils.make_chain_whole`, and because it follows bonds rather than taking a centre-of-mass image, lemonade's intra-chain observables never tear a chain. `paint_frame_grid_3d` and `_2d` write a per-bead value (chain index plus one, so 0 stays empty) into an empty grid for one frame. Bounds checks are off, so positions must already be in the box, which `lemonade.load` guarantees.

## Rules when you edit a kernel

- **Bounds checks are mostly off.** `mega_crank_fast` turns them off for the whole module, the others function by function, so check the decorators. Every index must be in range by construction, and a stray one silently reads the wrong memory. The public megamove kernels call the `_guard_*` helpers on entry, before anything changes, which raise `ValueError` if the grids, table, selectors, chain arrays or frozen mask don't fit together. They don't check the table's contents. New entry points must call them too.
- **No runtime flags inside the shell loops.** Write a separate body, as the periodic and wall energy bodies do.
- **Keep the serial crankshaft bit-exact with the reference kernel; the regression fixtures depend on it.** The generator, the order and number of draws, double-precision `randint` and both acceptance stages are all part of that, so change `accept_or_reject` or `mc_rand_below` in both modules. Anything that changes a random stream must regenerate `pimms/tests/simulation_tests/expected_output/` and say so in the changelog.
- **Keep the hardwall methods paired.** SR energies come from periodic pair lists; the summers and kernels turn SR wall contacts into solvent and drop LR and SLR ones.
- **Fix everything in a parallel sweep before the first draw.** The movable set, block layout and attempt counts must not depend on what happens during the sweep. A new whole-chain kernel takes offsets, lengths and homopolymer flags from `chain_layout()` and is chosen by the length partition, never by the current configuration.
- **Rebuild with `./build.sh uv` (or `./build.sh pip`) after editing a `.pyx`.** `setup.py` always regenerates the C, but a plain reinstall can reuse stale object files, and a stale `.so` can shadow a failed build. Then run `pytest pimms/tests/test_kernel_correctness.py pimms/tests/test_detailed_balance.py pimms/tests/test_proposal_symmetry.py pimms/tests/test_bookkeeping.py pimms/tests/test_parallel_dispatch.py` before the full suite. If timings changed, regenerate the tables in `docs/advanced/parallelization.rst` with `pimms/fast_kernels/benchmark_parallel_2d.py`. The kernel micro-benchmarks are `pimms/fast_kernels/benchmark.py` and `benchmark_parallel.py`.
