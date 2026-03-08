# PIMMS Keyfile Keywords

This document summarizes every keyword parsed by `keyfile_parser.py` (i.e., every keyword in `CONFIG.EXPECTED_KEYWORDS`).

## Notes

- Keywords are parsed case-insensitively, but chain sequence handling depends on `CASE_INSENSITIVE_CHAINS`.
- If a keyword is omitted, parser defaults are applied from `CONFIG.DEFAULTS` unless the keyword is required.
- Multi-entry keywords (can appear multiple times): `CHAIN`, `EXTRA_CHAIN`, `ANA_RESIDUE_PAIRS`.
- Experimental keywords require `EXPERIMENTAL_FEATURES : TRUE`.

## Required Keywords

- `DIMENSIONS`
- `TEMPERATURE`
- `N_STEPS`
- `PARAMETER_FILE`
- `EQUILIBRATION`

## Keyword Reference

| Keyword | Type | Description |
| --- | --- | --- |
| `DIMENSIONS` | int list (2 or 3 values) | Simulation box size in lattice units. Also determines 2D vs 3D simulation mode. |
| `LATTICE_TO_ANGSTROMS` | float | Conversion factor used when writing structure/trajectory dimensions in physical units. |
| `CHAIN` | repeated: `<count> <sequence>` | Defines polymer components in the initial system. Can be provided multiple times for different chain types. |
| `TEMPERATURE` | float (>0) | Simulation temperature. If quenching is enabled, startup temperature may be set from `QUENCH_START`. |
| `N_STEPS` | int (>0) | Total number of Monte Carlo steps. |
| `PARAMETER_FILE` | string path | Path to force-field/interaction parameter file. |
| `EQUILIBRATION` | int (>=0) | Number of equilibration steps (analysis/output behavior depends on other settings). |
| `RESIZED_EQUILIBRATION` | int list (2 or 3 values) | Optional smaller box dimensions for equilibration before switching to production dimensions. |
| `EQUILIBRATION_OFFSET` | int list (2 or 3 values) | Positional offset for resized equilibration box placement; must be compatible with production dimensions. |
| `HARDWALL` | bool | Use reflective hardwall boundaries instead of periodic boundary conditions. |
| `EXPERIMENTAL_FEATURES` | bool | Enables use of experimental/non-supported keywords and configurations. |
| `PRINT_FREQ` | int (>0) | Frequency for status/progress messages. |
| `REDUCED_PRINTING` | bool | Reduces runtime verbosity. |
| `XTC_FREQ` | int (>0) | Frequency for trajectory frame output. |
| `EN_FREQ` | int (>0) | Frequency for writing energy values (e.g., `ENERGY.dat`). |
| `SEED` | int (>0) | Random seed used for reproducible initialization and MC behavior. |
| `ENERGY_CHECK` | int (>0) | Frequency for full energy consistency checks. |
| `ANALYSIS_FREQ` | int (>0) | Base analysis frequency used by analysis outputs not explicitly overridden. |
| `NON_INTERACTING` | bool | If true, runs in excluded-volume/non-interacting limit by zeroing interaction terms. |
| `ANGLES_OFF` | bool | If true, disables angle potentials. |
| `CRANKSHAFT_SUBSTEPS` | int (>0) | Number of crankshaft substeps used when a crankshaft move is selected. |
| `CRANKSHAFT_MODE` | string | Parsed keyword, but currently marked obsolete and raises an exception if set. |
| `MOVE_CRANKSHAFT` | float probability | Probability weight for crankshaft moves. |
| `MOVE_CHAIN_TRANSLATE` | float probability | Probability weight for rigid single-chain translation moves. |
| `MOVE_CHAIN_ROTATE` | float probability | Probability weight for rigid single-chain rotation moves. |
| `MOVE_CHAIN_PIVOT` | float probability | Probability weight for single-chain pivot moves. |
| `MOVE_HEAD_PIVOT` | float probability | Probability weight for end-anchored pivot moves. |
| `MOVE_SLITHER` | float probability | Probability weight for slither/reptation-like moves. |
| `MOVE_CLUSTER_TRANSLATE` | float probability | Probability weight for rigid translation of detected clusters. |
| `MOVE_CLUSTER_ROTATE` | float probability | Probability weight for rigid rotation of detected clusters. Requires square/cubic boxes. |
| `MOVE_CTSMMC` | float probability | Probability weight for chain-level TSMMC move attempts. |
| `MOVE_MULTICHAIN_TSMMC` | float probability | Probability weight for multi-chain TSMMC move attempts. |
| `MOVE_RATCHET_PIVOT` | float probability | Probability weight for ratchet-pivot style move attempts. |
| `MOVE_SYSTEM_TSMMC` | float probability | Probability weight for system-level TSMMC move attempts. |
| `MOVE_JUMP_AND_RELAX` | float probability | Probability weight for jump-and-relax style move attempts. |
| `QUENCH_RUN` | bool | Enables temperature-quench mode. Requires quench-specific parameters to be set. |
| `QUENCH_FREQ` | int (>0 when used) | Step interval between quench temperature updates. |
| `QUENCH_STEPSIZE` | float (>0 when used) | Magnitude of temperature change per quench update (direction inferred from start/end). |
| `QUENCH_START` | float | Starting temperature for quench mode. |
| `QUENCH_END` | float | Final target temperature for quench mode. |
| `QUENCH_AS_EQUILIBRATION` | bool | If true, equilibration length is set to quench duration so production starts at final quench temperature. |
| `TSMMC_JUMP_TEMP` | float | High-temperature point used in TSMMC transitions (must exceed simulation temperature unless fixed offset is used). |
| `TSMMC_STEP_MULTIPLIER` | int (>0) | Controls TSMMC schedule granularity/spacing. |
| `TSMMC_INTERPOLATION_MODE` | string | TSMMC interpolation mode. Currently `LINEAR` is accepted. |
| `TSMMC_NUMBER_OF_POINTS` | int (>0) | Number of interpolation points used for TSMMC schedule construction. |
| `TSMMC_FIXED_OFFSET` | float or omitted | Optional TSMMC offset mode parameter. If omitted, defaults to `False` and jump-temperature sanity checks apply. |
| `ANA_POL` | int (>0) | Frequency for per-chain polymer statistics analysis. |
| `ANA_INTSCAL` | int (>0) | Frequency for internal scaling analysis. |
| `ANA_DISTMAP` | int (>0) | Frequency for distance-map analysis. |
| `ANA_ACCEPTANCE` | int (>0) | Frequency for move acceptance diagnostics output. |
| `ANA_INTER_RESIDUE` | int (>0) | Frequency for inter-residue pair-distance analysis. |
| `ANA_CLUSTER` | int (>0) | Frequency for cluster analysis. |
| `ANA_RESIDUE_PAIRS` | repeated: `<res_i> <res_j>` | Residue index pairs used by inter-residue analysis; each entry is normalized so lower index comes first. |
| `WRITE_CHAIN_TO_CHAINID` | bool | Writes mapping file from chain identity to chain IDs. |
| `ANALYSIS_MODULE` | string path | Path to custom Python analysis module to import. |
| `ANA_CUSTOM` | int (>0) | Frequency for invoking custom analysis callbacks from `ANALYSIS_MODULE`. |
| `ANA_CLUSTER_THRESHOLD` | int (>=0) | Minimum cluster size threshold used by cluster analysis routines. |
| `RESTART_FREQ` | int (>0) | Frequency for writing restart snapshots. |
| `RESTART_FILE` | string path | Path to restart file used to initialize system state. |
| `RESTART_OVERRIDE_DIMENSIONS` | bool | If true, use dimensions from restart file instead of keyfile dimensions (subject to other checks). |
| `RESTART_OVERRIDE_HARDWALL` | bool | If true, use hardwall setting from restart file. |
| `EXTRA_CHAIN` | repeated: `<count> <sequence>` | Adds new chains on top of restart state. Valid only when `RESTART_FILE` is provided. |
| `CASE_INSENSITIVE_CHAINS` | bool | If true (default), chain sequences are uppercased during sanity checks. |
| `AUTOCENTER` | bool | If true, enables automatic centering behavior for applicable single-chain setups. |
| `SAVE_AT_END` | bool | If true, holds trajectory in memory and writes at end of run. |
| `SAVE_EQ` | bool | If false, equilibration frames are not written to trajectory output. |
| `FREEZE_FILE` | string path | Path to freeze specification file listing chain IDs constrained during simulation. |

## Parser-Derived (Not Parsed Directly)

These are generated by parser logic and are not keyfile entries:

- `__TSMMC_USED`: internal flag set if any TSMMC move probability is non-zero.
- `EQUILIBRIUM_TEMPERATURE`: derived final/equilibrium temperature (e.g., from quench settings).

## Sanity Rules Worth Remembering

- All `MOVE_*` probabilities must sum to `1.0`.
- `QUENCH_*` settings must be fully specified when `QUENCH_RUN` is enabled.
- Non-square/non-cubic dimensions require experimental features and have additional constraints.
- `EXTRA_CHAIN` is invalid without `RESTART_FILE`.
- Cluster rotation moves require square/cubic dimensions.
