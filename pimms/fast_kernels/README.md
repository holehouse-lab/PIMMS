# fast_kernels — optimized PIMMS crankshaft kernel (drop-in, A/B comparable)

A self-contained, **drop-in optimized replacement** for the hottest part of
PIMMS — the crankshaft Monte Carlo kernel `pimms.mega_crank.mega_crank`, which
profiling showed to be **~85% of total simulation runtime**.

Everything here is kept *separate from the main package* so the new and old
kernels can be compared side by side without touching production code.

## What's the win

| Measurement (two_phase_equilibrium_demo) | Result |
|---|---|
| Kernel micro-benchmark (10k substeps × many megamoves) | **~10× faster** |
| Whole 5000-step simulation, end to end | **~5.4× faster** (≈56s → ≈10s) |
| Correctness | **bit-exact** — identical energy, accepted count, and final grid |

## Why it's faster (and still bit-exact)

The reference kernel is algorithmically fine but allocates numpy arrays
(`np.zeros(...)`) *inside the per-substep hot loop* — in `crank_it`,
`update_position`, `single_bead_crank`, and **four** arrays per call in
`get_angle_energy_change`. Across ~37M substeps per run that is hundreds of
millions of heap allocations routed through Python/numpy.

`mega_crank_fast.pyx` performs the **identical algorithm with the identical
random-number call sequence**, but every transient buffer is a **C stack array**.
Because the RNG stream is preserved exactly, the same inputs + seed yield the
**same** `(energy, accepted_moves)` and the same mutated `grid/type_grid/
idx_to_bead` — so the optimization is verifiable by exact comparison, not just
"looks close".

The kernel helpers are also written `noexcept nogil`, which makes them ready to
be called from `cython.parallel.prange` — the groundwork for the parallel
checkerboard scheme described in `CHECKERBOARD_DESIGN.md`.

## Where the kernel lives

The kernel itself is a normal package module: **`pimms/mega_crank_fast.pyx`**,
built with the package (it is listed in `setup.py`'s `ext_modules`, with graceful
OpenMP for the parallel path). It exposes serial 3D `mega_crank(...)`, serial 2D
`mega_crank_2D(...)`, and parallel 3D `mega_crank_parallel(...)`, each matching the
corresponding reference signature exactly. It is already wired into
`pimms.moves` (every crankshaft call uses it), so there is nothing to "turn on".

This directory holds only the **comparison / validation tooling** for it.

## Files

| File | Purpose |
|---|---|
| `benchmark.py` | A/B **correctness (bit-exact) + speed** on a real system state |
| `benchmark_parallel.py` | parallel kernel **energy-consistency + thread scaling** |
| `end_to_end.py` | runs the **whole** demo simulation both ways and compares final energy + wall time |
| `validation_system/` | a small system exercising paths the two_phase demo doesn't (non-zero angles, PBC/hardwall=0, single beads, OXO-central beads) |
| `validation_large/` | a large, spatially dispersed system for parallel scaling |
| `CHECKERBOARD_DESIGN.md` | design + feasibility for parallelizing local moves |

## Build

The kernel builds with the package - nothing special to do:

```bash
pip install -e .                      # or, to rebuild just the extensions:
python setup.py build_ext --inplace
```

## Run the comparisons

From the repository root:

```bash
# bit-exact correctness + kernel speedup (default: two_phase demo)
python pimms/fast_kernels/benchmark.py

# same, on the angle/PBC/single-bead validation system
python pimms/fast_kernels/benchmark.py pimms/fast_kernels/validation_system 4000 6

# parallel kernel: energy-consistency + thread scaling
python pimms/fast_kernels/benchmark_parallel.py

# whole-simulation end-to-end (fixed seed -> identical trajectory)
python pimms/fast_kernels/end_to_end.py
```

## Validation coverage

Bit-exact agreement has been confirmed on:

- **two_phase_equilibrium_demo** — short-range energy path, `hardwall=1`, chains
  of length 9 (bead flags 1/5/2/6/3), zero angle penalties.
- **validation_system/** — short-range energy path, `hardwall=0` (PBC), non-zero
  angle penalties, single beads (flag 0) and OXO-central beads (flag 4).
- **2D kernel (`mega_crank_2D`)** — bit-exact vs the reference `mega_crank_2D`
  on both PBC and `hardwall=1`, with non-zero angles (5D angle lookup), LR/SLR
  interactions, and single/terminal/internal/OXO beads. ~17x faster than the
  reference 2D kernel (which crossed the Python/C boundary every substep via
  `randint_ext` / `accept_or_reject_ext` / the `def` `get_energy_change_2D`).

The long/super-long-range energy branch (7×7×7) is a faithful line-by-line port
of the reference and is exercised bit-exactly by the regression test suite
(several tests define LR/SLR interactions); pass such a keyfile dir to
`benchmark.py` to check it directly.

## Notes / caveats

- Requires the same toolchain as the main package (Cython + a C compiler +
  numpy headers). `cimport pimms.cython_config` resolves because the kernel
  lives inside the `pimms` package.
- Bit-exactness depends on libc `rand()` behaving identically to the reference
  (it does — both use the same `srand`/`rand` and the same `randint` math).
- The generated `mega_crank_fast.c` and the `.so` are git-ignored.
