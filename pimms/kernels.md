# The compiled kernels

This is the map of everything in PIMMS that is written in Cython: what each module is for, what data it works on, how that data gets from the Python objects into the compiled code and back, where the random numbers come from, and the conventions every kernel has to agree on. The per-function detail lives in the docstrings of the `.pyx` files; this document is about the logic and the information flow, so that you can find your way around and know what you are allowed to change. The Sphinx build mocks these modules, which is why they are not in the API reference.

## The modules at a glance

| Module | What it holds | Called from | Threads |
|---|---|---|---|
| `pimms/mega_crank_fast.pyx` | The production megamove kernels: serial and parallel crankshaft, slither and pull, in 2D and 3D, plus the shared energy-delta and angle helpers and the two random-number generators | `moves.py` (`system_shake`, `system_slither`, `system_pull`, the TSMMC and single-chain relaxations) | OpenMP `prange` in the parallel kernels |
| `pimms/mega_crank.pyx`, `pimms/mega_crank_2D.pyx` | The original, numpy-allocating crankshaft kernels. Kept as the bit-exactness oracle for the fast kernel and never called by a production move | `tests/kernel_test_utils.py`, `tests/test_kernel_correctness.py`, the benchmark scripts, and one harmless seeding call at start-up | none |
| `pimms/bookkeeping.pyx` | The two copies between the Chain objects and the kernels' bead table | `crankshaft_list_functions.update_idx_to_bead` / `write_back_positions` | none |
| `pimms/inner_loops.pyx`, `pimms/inner_loops_hardwall.pyx` | Pair-list extractors: every short-, long- and super-long-range neighbour of one site, periodic and hard-wall variants | `lattice_utils.build_envelope_pairs` / `build_all_envelope_pairs`, `longrange_utils.build_LR_envelope_pairs` | none |
| `pimms/hyperloop.pyx` | Sums a pair list against an interaction table, the chain angle energy, and small site helpers | `energy.Hamiltonian` (the from-scratch total and the per-move deltas of the Python-level moves) | none |
| `pimms/cluster_kernels.pyx` | The breadth-first single-image gather used by the cluster analysis | `cluster_utils.convert_positions_to_single_image_snakesearch` | none |
| `pimms/system_utils.pyx` | Start-up checks that the compiled and Python integer types agree and that the chain count fits the grid dtype | `simulation.py` set-up | none |
| `pimms/lemonade/kernels/_pbc.pyx` | lemonade's batched chain unwrap and per-frame grid paint | `lemonade/_store.py` | none |
| `pimms/cython_config.pxd` | The two integer typedefs every kernel shares | cimported by the modules above | |

Only `mega_crank_fast` is built with the OpenMP flags. If OpenMP is missing at build time the module still compiles and `prange` runs serially with identical results; `mega_crank_fast.openmp_info()` reports which you have.

## The data every kernel works on

**Integer types.** `cython_config.pxd` defines `NUMPY_INT_TYPE` as `int32` and `NUMPY_INT_TYPE_long` as `int64`. The Python twin of the first is `CONFIG.NP_INT_TYPE`; the two must agree, and `system_utils.check_dtype_consistency` proves it at start-up by binding a Python-dtype array to the compiled buffer type. Grids, interaction tables and angle tables are int32; the bead table and the bead selector are int64. Lattice energies are integers end to end: the kernels accumulate into C `long` and return Python ints.

**The two grids.** `latticeObject.grid` holds the occupying chainID at each site, 0 for empty, and `latticeObject.type_grid` holds the occupying bead's integer residue code, again 0 for empty. Both are C-contiguous int32 arrays of the box shape. The kernels read and write both in place; a kernel never sees the Chain objects. Because the grid stores chainIDs, `system_utils.check_beads_to_grid_mapping` refuses a system with more chains than the dtype can label.

**The bead table.** `latticeObject.crankshaft_lists` is an int64 array with one row per bead, rows in ascending chainID and then chain position, so that rows `i - 1` and `i + 1` are a bead's bonded neighbours. The columns are

| column | content |
|---|---|
| 0 | bead flag (below) |
| 1 | long-range flag, 1 if this bead takes part in long-range interactions |
| 2 | intcode, the integer residue type that indexes the interaction and angle tables |
| 3 | skip-angle flag, 1 for beads of 1- and 2-bead chains, which have no angle |
| 4 | chainID |
| 5, 6, 7 | x, y, z (the 2D table has seven columns and no z) |

The bead flag tells a kernel which proposal to make and which angle triplets a move can change:

| flag | bead | angle window (rows relative to the bead) |
|---|---|---|
| 0 | a single-bead chain | none |
| 1 | first bead of a chain | 0 to +2 |
| 3 | last bead | -2 to 0 |
| 4 | middle bead of a 3-mer | -1 to +1 |
| 5 | second bead of a chain of four or more | -1 to +2 |
| 6 | second-to-last bead of a chain of four or more | -2 to +1 |
| 2 | any other interior bead | -2 to +2 |

Columns 0 to 4 never change during a run. The static per-chain layout that goes with the table, `ChainLayout` in `crankshaft_list_functions.py`, is built once in `Lattice.__init__` and cached on the lattice: the sorted chainIDs, each chain's first row and length as int64 and int32 arrays, the number of beads, and a homopolymer flag that is 1 when a chain's intcode and long-range columns are each constant. The cache is rebuilt if the chain count or the row count changes, which is what happens on the scheduled chain insertions and deletions; nothing else in PIMMS rewrites a chain's identity after construction.

**The round trip.** Every megamove does the same three things. `update_idx_to_bead` gathers the current positions out of the Chain objects (whose `positions` are Python lists of lists) into columns 5+ of the lattice's table with `bookkeeping.gather_positions`, then hands the kernel a fresh copy of the whole table. The kernel mutates that copy's position columns along with the two grids. `write_back_positions` scatters the copy back into the Chain objects with `bookkeeping.scatter_positions`, which builds each chain's list of lists in C and installs it through `Chain.set_ordered_positions` so the length check still runs. Between those two calls the lattice-owned table's position columns are stale; only the kernel's copy is current. Both compiled copies raise `ValueError` if the offsets, lengths or a chain's list have gone out of step with the table, and the Python loops they replaced are kept as the fallback and the test oracle (`test_bookkeeping.py` runs megamoves both ways from one seed). The TSMMC and single-chain relaxations use `update_idx_to_bead_single_chain` and `_multiple_chains` instead, which extract just the rows they need through `chain_to_firstbead_lookup` and write positions back themselves.

**Interaction and angle tables.** The Hamiltonian provides three `(n_intcodes, n_intcodes)` int32 tables, one per shell (short-range, long-range, super-long-range), with row and column 0 standing for an empty site, and an angle lookup indexed by the middle bead's intcode and the two bond vectors of a triplet, each component shifted from -1..1 into 0..2: seven dimensions in 3D, five in 2D. `ANGLES_OFF` zeroes that table rather than bypassing the lookups.

## Where the random numbers come from

There are three generators and it matters which one a draw comes from, because the restart file records the state of the first two and the third is reseeded from the first on every call.

1. **Python's global `random`.** Seeded from `SEED` at start-up. Every megamove draws `local_seed = random.randint(1, sys.maxsize - 1)` from it and passes that to the kernel. The move-selection layer and the Python-level moves also draw from it.
2. **numpy's global generator.** Seeded from `SEED` at start-up. `bead_selector_constructor` draws the crankshaft's bead indices from it (`np.random.randint`, or `np.random.choice` over the selectable beads when chains are frozen), and the slither and pull selectors are shuffled with it. Inside the parallel kernels a `np.random.RandomState(passed_seed & 0x7FFFFFFF)` draws the block-origin shift, so that too is a pure function of the seed.
3. **The kernels' own streams.** The serial kernels share one module-global splitmix64 state in `mega_crank_fast`, reset by `mc_seed(passed_seed)` as the first statement of every call, so a call is a pure function of its inputs. `mc_rand` returns the top 31 bits of a splitmix64 step, in `[0, 2^31 - 1]`; `randint(start, end)` is inclusive and takes exactly one step; `accept_or_reject` takes a step only for an uphill move. The parallel kernels give every block its own splitmix64 state on the stack, seeded by a finaliser mix of the block index and `passed_seed` so that no two (block, sweep) pairs share a stream, and use a 53-bit uniform.

The reference kernel in `mega_crank.pyx` carries a byte-for-byte copy of the serial generator and the same `randint` formula, evaluated in double precision, which is what makes the fast serial kernel bit-exact with it: same generator, same number of draws per proposal in the same order, same floating-point path. Every serial substep therefore consumes a fixed number of draws before the acceptance test (three for a 3D crankshaft proposal, two in 2D, one direction plus three target offsets for a slither) whether or not the proposal is rejected for a clash or a wall, and one more only if the move is uphill. Do not reorder or add draws in those kernels without regenerating the regression fixtures.

Start-up also seeds the reference module's generator through `mega_crank.seed_C_rand`. That call has no effect on any production move, since every kernel reseeds itself from `passed_seed`; it is a leftover of the libc `srand` era.

## The energy conventions everything shares

**Shells.** Interactions are counted over Chebyshev shells around a site: the 3x3x3 block minus the centre is the short-range (SR) shell, the 5x5x5 minus the 3x3x3 is the long-range (LR) shell, the 7x7x7 minus the 5x5x5 is the super-long-range (SLR) shell. A bead whose long-range flag is 0 only ever scores the SR shell. The same three radii define long-range cluster membership in `lattice_utils`, the cluster gather's link predicate and the VMMC neighbour shell.

**Solvent.** Empty sites are real partners in the SR shell: a bead next to an empty site scores `table[type, 0]`, which is how solvation enters, and when a bead vacates a site its former SR neighbours gain a solvent contact while the neighbours of the site it fills lose one. The LR and SLR shells have no solvent term at all; the pair extractors skip empty sites there and the megamove kernels rely on column 0 of the LR and SLR tables being zero.

**Bonded pairs and self-terms.** A bond is always a Chebyshev-1 contact, so every bonded pair sits in the SR pair list; it contributes a constant to the total and cancels in every delta, and nothing special-cases it. The megamove energy delta visits the centre site of both shells and subtracts the bead's self-interaction and self-solvent entries afterwards.

**Angles.** For each triplet the two bond vectors are taken relative to the middle bead, each component is folded back into -1..1 (a bond that crosses the periodic boundary shows up as a component of magnitude `L - 1`), and the penalty is `angle_lookup[intcode of the middle bead, a + 1, b + 1]`. `hyperloop.evaluate_angle_energy_*`, the reference kernel and `mega_crank_fast` all implement exactly this fold.

**Hard walls.** Two idioms coexist and they are correct only as they are paired today.

- The energy-delta kernels in `mega_crank_fast` and the pair-list summers in `hyperloop` detect a pair that reaches its partner only through the periodic wrap (any per-axis separation larger than 3) and then read that partner as empty in the SR shell, so a wall behaves like solvent, and skip it entirely in the LR and SLR shells.
- The `inner_loops_hardwall` extractors simply never emit a pair that leaves the box.

Because the summers turn wall contacts into solvent terms, the SR pair list for an energy evaluation under a hard wall must come from the periodic extractor, and it does: `evaluate_total_energy` and `single_chain_move` build their SR lists without the hardwall flag, and the hardwall extractors are used only for LR and SLR lists and for connectivity. The separation-larger-than-3 test is sound only because the keyfile parser refuses any box dimension below 7; the kernels do not check it themselves. VMMC scans its neighbours in pure Python and skips out-of-box sites in every shell, which is consistent because a link energy never contains a solvent term.

**Which Python move uses which kernel.**

| Move | Energy path |
|---|---|
| crankshaft, slither, pull, the TSMMC relaxations, the single-chain relaxation | `mega_crank_fast` |
| jump-and-relax | two crankshaft-kernel relaxations around a jump scored by `evaluate_total_energy` |
| chain translate, chain rotate, chain pivot, head pivot | `simulation.single_chain_move`: pair lists from `inner_loops`, summed by `hyperloop`, the four-way old-region/new-region accounting per shell plus the angle difference |
| cluster translate, cluster rotate | `simulation.rigid_cluster_move`: LR and SLR only (a rigid cluster's SR interface is solvent on both sides), returns 0 outright when no residue has long-range interactions |
| VMMC | pure Python over the offset shell, cross-chain pair terms only |
| the from-scratch total and every test oracle | `Hamiltonian.evaluate_total_energy` over the whole-system pair lists from `inner_loops`, summed by `hyperloop` |

## The megamove kernels (`mega_crank_fast.pyx`)

### The energy delta

`get_energy_change_c` (and `get_energy_change_2D_c`) is what every substep of every megamove pays for. It takes the bead's old and new sites, its long-range flag and the three tables, precomputes the wrapped coordinate along each axis for offsets -3..3 once (42 values instead of three modulo operations per site), sums the old shell, moves the bead in `type_grid`, sums the new shell, subtracts the self-terms, puts the bead back and returns new minus old with the solvent re-scoring included. Callers therefore see `type_grid` unchanged; the new site must already be empty in `grid`.

The function is only a dispatcher over two bodies, `_energy_change_periodic_c` and `_energy_change_wall_c`, with the same signature minus the flag. The periodic body contains no hardwall code at all. That split is a performance decision, not a stylistic one: with a runtime `if hardwall == 1` inside the 7x7x7 loops the C compiler keeps the branch on every one of the 686 site visits and cannot unroll them, which measured as 2.5 to 3 times slower for every LR and SLR evaluation, invisible to a profiler because all of the time was simply "in the energy function". Keep it that way: any new per-site condition belongs in a separate body, never inside the loops.

The angle delta, `get_angle_energy_change_c`, returns 0 for a skip-angle bead, otherwise gathers the window of rows the bead flag dictates, scores the triplets with the current positions, substitutes the proposed position and scores them again.

### Serial crankshaft: `mega_crank`, `mega_crank_2D`

Inputs: the two grids, the bead table copy, the three tables, the angle table, the current total energy, `invtemp`, the number of substeps, the pre-drawn bead selector, the seed and the hardwall flag. Returns `(energy, accepted)` and mutates the grids and the table's position columns in place on every accepted move.

One substep: take the bead index from the selector (no draw); propose by flag. A single bead (flag 0) draws three offsets in -1..1 around its own site; a terminal bead (flags 1 and 3) draws the same three offsets around its one bonded neighbour, so it can land anywhere in that neighbour's 3x3x3 cube; any other bead draws uniformly inside the intersection of its two neighbours' cubes, with the anchors lifted into the extended box first when they straddle the wrap. The proposal is rejected if the target is occupied. Under a hard wall a proposal whose bond would cross the wrap is rejected as well, after the draws have been consumed. A surviving proposal gets the interaction delta plus the angle delta and a Metropolis test on the total; on acceptance the grids, the table and the running energy are updated.

The RNG budget is three draws per proposal in 3D (two in 2D) plus one if the total delta is uphill. `simulation.py` chooses the bead selector in Python precisely so that bead selection and the perturbation come from different streams; deriving both from one stream once correlated the bead index with the perturbed axis.

### Parallel crankshaft: `mega_crank_parallel`, `mega_crank_parallel_2D`

The box is cut into blocks along each axis, at most four per axis and at least `8W` sites long, with `W` the halo width: 2 when any bead carries the long-range flag, 1 otherwise. `W` is not the interaction range. Every block boundary has a halo on both sides, so two concurrently moving interiors are `2W` apart; writes stay at least `W` inside a block face and the widest read (the SLR shell, the occupancy check two sites out, the angle window two bonds along the chain) reaches at most `max(R_int, 2)` past it, so the condition is `2W >= max(R_int, 2)`. `parallel_crank_layout_info` exposes the layout and the movable fraction, and `system_shake` only uses the parallel kernel when it yields more than one block.

Each sweep draws a random origin shift per split axis, buckets every bead into a block or into "immovable" (halo, trailing remainder, or frozen), builds CSR lists per block, splits the `nsteps` attempts across blocks in proportion to their movable-bead counts with a largest-remainder top-up so exactly `nsteps` are made, derives a seed per block, and runs `run_block` over the blocks with a dynamic `prange`. Inside a block each attempt picks a bead uniformly from the block's list, proposes exactly as the serial kernel does with the block's own stream, rejects any proposal whose new site would leave the interior, scores it with the same two delta helpers, and accepts on the delta alone (`accept_p`), which is what lets each block keep a private accumulator. The per-block deltas and acceptance counts are summed serially afterwards.

Why it samples the same equilibrium: within a sweep the movable set is fixed, the proposals are the serial kernel's symmetric proposals, the interior is closed (a halo bead is never selected, so an interior-to-halo move could never be reversed and is rejected), and acceptance is Metropolis on the exact local delta, so each sweep preserves the Boltzmann distribution; the random shift moves the halos every sweep, so over a run every bead moves. Why it relaxes more slowly: only interior beads move in a sweep and the attempts are spent on them alone, so compare equilibrium averages between serial and parallel runs, never energies at a fixed step count. The decomposition depends only on geometry and `W`, never on the thread count, so results are identical for any number of threads. Frozen chains enter only through `frozen_mask`, which removes their beads from the movable set while leaving them in the grids as obstacles.

### Slither: `mega_slither`, `mega_slither_2D` and the parallel twins

Inputs are the crankshaft's plus the per-chain offsets, lengths and homopolymer flags from the layout, a chain selector (each eligible chain repeated `SLITHER_SUBSTEPS` times, shuffled) and `max_chain_len`, the size of the malloc'd revert buffers.

One attempt: a single-bead chain becomes a local translation through the crankshaft's monomer proposal. Otherwise draw a direction (head-grow moves the last bead forward, tail-grow moves the first bead backward), draw a target uniformly in the 3x3x3 cube around the growing end, reject if it is occupied (that includes the end itself and the site the other end is about to vacate) or, under a hard wall, reachable only through the wrap. For a homopolymer the interaction change is exactly the single-bead delta of moving one bead of the chain's type from the vacated end to the new site, because the set of occupied sites afterwards is the old set minus one plus one and every bead has the same type; the angle change is the whole-chain angle energy in the shifted arrangement minus the current one, computed by `_chain_angle_mode` without touching the table. For a heteropolymer every bead's type moves one site along the path, so the reptation is decomposed into `L` sequential single-bead moves, each landing on the site the previous bead just vacated and each scored against the live grids; the deltas telescope to the exact total, and a rejection restores the whole chain from the buffers. On acceptance the table's rows for the chain are shifted by one.

The parallel kernel cannot use a per-bead halo, since a slither moves a whole chain. Instead a chain may move in a block only if every one of its beads lies in that block's interior after the random shift, with the wider chain-level halo `W = R_int + 2` (5 with long-range interactions, 3 without) and blocks of at least `4W` sites. The only new occupied site a slither creates is the growing end, and a proposal whose new end would enter the halo is rejected, so a chain stays inside its block for the whole sweep, footprints are disjoint and blocks run lock-free. `frozen_mask` marks every bead of every chain the parallel kernel must not touch; attempts per block are proportional to the movable chain count and each attempt picks a chain uniformly with replacement, so per-chain attempt counts differ from the serial kernel's exact `SLITHER_SUBSTEPS`. Heteropolymers longer than the 512-bead per-thread buffer are skipped inside the kernel; the Python side keeps them out of the parallel set.

**The length partition.** `system_slither` and `system_pull` split the chains between the parallel and the serial kernel once, from chain lengths alone (`parallel_chain_partition`): a chain goes to the parallel kernel if its length fits inside the smallest block interior and it is under the buffer cap. Both kernels are stationary on their own, but choosing between them from the current configuration is not: the parallel kernel's movable region is closed, so the rate out of "every chain is compact enough" is zero through it while the serial kernel crosses that boundary freely, and gating on current extents pumped probability into compact conformations and biased the mean squared radius of gyration low by a few percent, silently. `parallel_chain_fit_report` still describes the current configuration for the start-up summary and must never gate the dispatch. Both passes run every megamove, serial after parallel, chained through the running energy on the same buffers.

### Pull: `mega_pull`, `mega_pull_2D` and the parallel twins

Chains shorter than three beads are skipped. One attempt: draw an interior bead `i` uniformly in `[1, L - 2]` and a direction, which fixes the anchor (`i - 1` or `i + 1`). Enumerate the first targets: every empty site Chebyshev-adjacent to both the bead's current site and the anchor, `nF` of them; if there are none the attempt ends with nothing moved. Draw one, move bead `i` there and score it, then cascade away from the anchor: each following bead moves into the previous bead's old site until the next bead is already adjacent to the last moved one, which restores connectivity. If the cascade reaches the chain end without restoring, the proposal is rejected and reverted. Otherwise count the reverse first targets `nR` with the same predicate in the final state (the last moved bead, the opposite direction) and accept with the Metropolis-Hastings factor `(nF / nR) exp(-beta dE)`; the bead and direction probabilities cancel and the cascade is deterministic given the first target, so that ratio is the whole Hastings correction. A rejection restores the moved range from the buffers by clearing every current site first and then rewriting the old ones, which is what makes the overlap between old and new sites safe.

The parallel kernel uses the same chain-level block scheme as the slither. The one difference is `pull_first_targets_interior`: first targets are restricted to the block interior, so the single new site a pull creates never enters the halo, and `nR` is counted with the same restricted predicate. That is required rather than convenient: in the parallel chain the reverse proposal really is drawn from the restricted set, so the Hastings ratio must use the restricted counts on both sides.

### The 2D twins

Every kernel has a 2D version with two coordinates, a 7-column bead table, a five-dimensional angle table, 7x7 and 3x3 rings, one fewer draw per proposal and the hardwall tests called with zero for the z arguments. The logic is line for line the same.

## The reference kernels (`mega_crank.pyx`, `mega_crank_2D.pyx`)

These are the original crankshaft kernels. They allocate small numpy arrays inside the substep loop, which is why they are slow and why `mega_crank_fast` exists. They are still built because `test_kernel_correctness.py` runs both from the same bead selector and seed and requires identical energy, grids and bead table, and because `test_proposal_symmetry.py` uses them as the proposal-uniformity oracle. The 2D reference has no generator of its own and calls the 3D module's through Python. Their seed argument is a C `int`, so the bit-exactness test can only be run with seeds below 2^31; production seeds are 63 bits and only ever reach the fast kernel.

## Pair lists and their summers (`inner_loops*.pyx`, `hyperloop.pyx`)

The Python-level moves and the from-scratch total do not use the megamove kernels. They build explicit pair lists. `inner_loops` exposes, per site, extractors that return int32 arrays of shape `(n, 2, ndim)`: the 26 (8 in 2D) SR pairs including empty partners, and the LR and SLR pairs restricted to occupied partners, either all at once for a bead with the long-range flag or SR only. The two sites of a pair are ordered by the sign of the first non-zero component of the pre-wrap offset, an antisymmetric rule that makes the same physical pair come out identically from either end, so the Python side can de-duplicate the concatenated lists by exact row matching (`lattice_utils._unique_rows`). `inner_loops_hardwall` provides the same six functions but drops any partner outside the box.

`hyperloop` sums a pair list against a table into a C `long` in four variants (2D/3D, short-range/non-short-range) with the hardwall idiom described above, and evaluates a chain's angle energy from its positions and intcodes. `Hamiltonian.evaluate_total_energy` is the whole-system version; `simulation.single_chain_move` uses the same pieces for the old-region/new-region accounting of translate, rotate, pivot and head pivot; `rigid_cluster_move` evaluates LR and SLR only. The adjacent-site enumerators feed `lattice_utils.get_empty_site` and the head pivot.

## The cluster gather (`cluster_kernels.pyx`)

`snakesearch_single_image` takes the in-box positions of one cluster and a seed bead (the bead nearest the periodic centre of mass, chosen in Python) and returns a single-image copy: a breadth-first walk that places each newly reached bead at its parent's image plus the minimum-image separation of the two in-box positions. Site lookup is either a flat occupancy grid or a binary search over sorted encoded site keys, chosen by a measured crossover in box volume, bead count and shell size, and the walk is written out four times (2D/3D times grid/search) because a branch inside it cost 8 to 16 percent. With the residue types and the LR and SLR tables supplied, two beads are linked if they are Chebyshev-1 neighbours, or Chebyshev-2 with a non-zero LR entry for their types, or Chebyshev-3 with a non-zero SLR entry; this is the same relation `lattice_utils.get_all_chains_in_long_range_cluster` uses to define the cluster, so the gather walks the relation that built it rather than a plain distance, which used to tear long clusters through a periodic face. The result is shifted by whole box periods so no coordinate is negative, preserving congruence with the box. A cluster that does not connect raises `ValueError`; whether it winds the box is decided in Python by `cluster_utils.cluster_percolates`, which encodes the same link relation, and the pure-Python gather in `cluster_utils` is the byte-identical fallback and test oracle. The helpers are declared `nogil` for inlining, but the walk runs with the GIL held.

## lemonade's kernels (`lemonade/kernels/_pbc.pyx`)

`unwrap_chains` makes every chain whole in every frame by walking its bonds from the first bead: each subsequent bead is shifted by whole box lengths until it is within one site of its unwrapped predecessor. It is cumulative, so a chain longer than the box is handled and coordinates may leave the box; it is the batched twin of `lattice_utils.make_chain_whole`, and it is a bond walk, not a centre-of-mass image, which is why lemonade's intra-chain observables never tear a chain. `paint_frame_grid_3d` and `_2d` write `chainID + 1` into a zero grid for one frame with bounds checking off, so the positions must already be wrapped into the box, which `lemonade.load` guarantees.

## Rules when you edit a kernel

- **Bounds and wraparound checks are off for the whole of `mega_crank_fast`**, `bookkeeping`, `cluster_kernels` and `_pbc`. Every index must be in range by construction: a wrapped lattice coordinate, a table offset from the layout, a bead flag that guarantees the neighbour rows exist. A stray index is a silent read of the wrong memory, not an exception.
- **Never put a runtime flag inside the shell loops.** Specialise the body instead, as the periodic and wall energy bodies do.
- **The serial crankshaft is bit-exact with the reference kernel and the regression fixtures depend on it.** The generator, the order and number of draws per proposal, and double-precision evaluation of `randint` are all part of the contract. A change that alters any random stream must regenerate `pimms/tests/simulation_tests/expected_output/` and say so in the changelog.
- **Keep the hardwall idioms paired.** SR energy from periodic pair lists; the summers and the megamove kernels convert wall contacts to solvent in the SR shell and drop them in the LR and SLR shells.
- **Keep the parallel kernels state-independent within a sweep.** The movable set, the block layout and the attempt allocation are fixed before the first draw and must stay so; the closed-interior rejection is what makes the sweep symmetric. Any new whole-chain kernel takes offsets, lengths and homopolymer flags from `chain_layout()` and is dispatched by the length partition, never by the current configuration.
- **Rebuild with `./build.sh uv` (or `./build.sh pip`)** after editing a `.pyx`; a plain reinstall can reuse stale C. Then run `pytest pimms/tests/test_kernel_correctness.py pimms/tests/test_detailed_balance.py pimms/tests/test_proposal_symmetry.py pimms/tests/test_bookkeeping.py pimms/tests/test_parallel_dispatch.py` before the full suite, and regenerate the tables in `docs/advanced/parallelization.rst` with `pimms/fast_kernels/benchmark_parallel_2d.py` if timings changed. The kernel-only micro-benchmarks are `pimms/fast_kernels/benchmark.py` and `benchmark_parallel.py`.
