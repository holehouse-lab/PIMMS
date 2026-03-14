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
