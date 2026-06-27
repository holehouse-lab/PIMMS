# Checkerboard (domain-decomposition) parallelization of PIMMS local moves

> Status: **design + feasibility**, not yet implemented. The serial
> allocation-free kernel (`mega_crank_fast.pyx`, ~10× faster, bit-exact) is the
> prerequisite and is already done. This document specifies how to parallelize
> the crankshaft (`system_shake`) substep loop on top of it.

## TL;DR

- Lattice MC with **local** moves is parallelizable by **spatial domain
  decomposition** with a **checkerboard (red/black, 8-colour in 3D) schedule** —
  the standard, well-established approach.
- It is **not** a free lunch: naive static decomposition **violates detailed
  balance at block boundaries**, and PIMMS adds two wrinkles — a wide energy
  halo (±3) and **chain connectivity** (a crankshaft move touches its `i±1`
  neighbours). Both are solvable.
- The serial kernel functions here are already written `noexcept nogil`, so they
  can be called from `cython.parallel.prange` **without the GIL** — the main
  plumbing barrier is already removed.
- Expect **near-linear scaling for large boxes**, but **little benefit for the
  30³ demo** (too few independent blocks). Recommend shipping it as an opt-in
  `parallel=True` path validated by *ensemble observables*, not bit-exactness.

---

## 1. Why local-move MC can be parallelized at all

Two proposed moves are **independent** if neither one's evaluation reads a site
the other writes, and neither writes a site the other writes. If two beads are
far enough apart that their interaction halos don't overlap, their accept/reject
decisions are independent and can be made concurrently — the joint update is
identical to doing them in some serial order.

The **checkerboard** schedule guarantees this distance: partition the box into
cubic blocks, 2-colour them in each dimension (→ 8 colours in 3D), and update
**one colour at a time**. All same-colour blocks are separated by at least one
block width, so threads working different same-colour blocks can never collide.

```
  colour layout (2D slice, • = active this pass):
   ┌───┬───┬───┬───┐
   │ • │   │ • │   │     pass 1: all "•" blocks updated in parallel,
   ├───┼───┼───┼───┤             each by its own thread; then barrier;
   │   │ • │   │ • │             pass 2 shifts to the next colour, etc.
   ├───┼───┼───┼───┤
   │ • │   │ • │   │
   └───┴───┴───┴───┘
```

## 2. The PIMMS-specific interaction range (sets the halo width)

From `get_energy_change` in the kernel:

| Bead class | Neighbourhood scanned | Half-width |
|------------|----------------------|------------|
| short-range only (`LR_vs_SR==0`) | 3×3×3 | 1 |
| long/super-long-range (`LR_vs_SR==1`) | 7×7×7 | **3** |

A crankshaft move also displaces a bead by ±1. So the worst-case **halo width**
a move can read/write from its bead's position is

```
W = R_interaction + move_displacement = 3 + 1 = 4   (LR systems)
W = 1 + 1 = 2                                        (purely short-range systems)
```

**Consequence:** only beads whose position is at least `W` sites *inside* a block
(away from the faces shared with a different-colour neighbour) are safe to move
in a given colour pass. Beads within `W` of an active boundary are **frozen**
that pass. To keep every bead movable across the simulation, the block grid
origin is **randomly shifted each sweep** (see §4).

This also sets a **minimum useful block size**: `L_block ≳ 3·W` (so each block
has a non-trivial movable interior). For LR systems that's ~12–16 sites; for the
30³ demo you'd get only ~2 blocks/dim → 8 colours barely tile the box and there
is almost no parallelism. **Checkerboard pays off for large boxes / many beads**,
which is exactly the regime where runs are slow.

## 3. The chain-connectivity wrinkle (PIMMS-specific)

PIMMS is a **polymer** model: a crankshaft move on internal bead `i` reads its
anchors `i-1` and `i+1`, and terminal-bead moves read `i±1`. So a move's true
footprint is *bead i and its bonded neighbours*, each with its own ±W halo.

A chain can straddle a block boundary, so bead `i` may sit safely inside an
active block while anchor `i+1` lies in a concurrently-updated neighbour. Two
clean mitigations:

- **(A) Footprint test (recommended).** A proposed move is admissible in a pass
  only if `i` *and* its anchors *and* their halos lie within the current block's
  safe interior. Chains (9 beads in the demo) are short relative to a sensible
  block, so only a small boundary fraction is deferred. Deferred moves are
  picked up on later randomly-shifted sweeps.
- **(B) Chain ownership.** Assign each chain to the block holding its
  centre-of-mass bead; only that block may move any of its beads. Simpler
  bookkeeping, but a chain bridging two active blocks serializes onto one of
  them, lowering parallel efficiency near interfaces.

Start with **(A)**; it maximizes parallelism and composes naturally with the
random-shift ergodicity fix.

## 4. Detailed balance / correctness

Static domain decomposition with *fixed* boundaries breaks detailed balance:
boundary beads are systematically under-sampled and the frozen halo biases the
acceptance near interfaces. The standard, rigorous fixes:

1. **Random origin shift per sweep.** Before each full 8-colour sweep, offset the
   block lattice by a uniform random vector in `[0, L_block)^3`. Over many sweeps
   every bead spends time in block interiors, eliminating the boundary bias. This
   is the most common production approach (used in parallel lattice-gas / Ising /
   polymer codes).
2. **Frozen halo within a pass.** Within one colour pass, sites within `W` of an
   active boundary are immutable, so each block's sub-chain of accepted moves is
   exactly a valid serial Markov sub-chain on an isolated region — detailed
   balance holds *locally*, and the random shift restores it *globally*.

3. **Closed movable region — REJECT moves that leave the interior.** This is the
   subtle one (and the bug that was actually shipped first). Only interior beads
   are *selected*, so a halo bead can never make a move. Therefore an interior
   bead must **not** be allowed to move *into* the halo: that forward move would
   have no reverse (the halo bead is frozen), breaking detailed balance and
   biasing the ensemble. Concretely: after proposing a new position, reject the
   move (bead stays put) if the new position is not itself in the movable
   interior. Without this, a condensing system is actively driven apart — e.g.
   the two_phase droplet equilibrated to ~−20k instead of the correct ~−48k.

Net: the parallel sampler targets the **same equilibrium distribution** as the
serial one, but follows a **different trajectory**. It must therefore be
validated by **ensemble observables**, not bit-exact comparison (see §7).
NOTE: energy-consistency (incremental energy == from-scratch recompute) and
thread-independence do **not** detect a detailed-balance violation — both held
while the ensemble was wrong. Always check ensemble agreement with the serial
kernel (e.g. start from a serial-equilibrated config and confirm the parallel
kernel *holds* that energy rather than drifting).

## 5. Energy bookkeeping

The serial kernel threads a single `energy` scalar updated per accepted move.
In parallel:

- Each thread keeps a **private `delta_energy` accumulator** over the moves it
  accepts in its block (regions are non-overlapping → contributions are additive
  and independent).
- After each colour pass (barrier), **reduce** (sum) the per-thread deltas into
  the global energy. `prange(..., reduction)` or an explicit per-thread array +
  serial sum both work.

No locks are needed on `grid`/`type_grid` because the colouring guarantees
disjoint write sets within a pass.

## 6. RNG

`libc rand()`/`srand()` use a **single global state** and are not thread-safe;
calling them from `prange` would both race and serialize. Replace with a
**per-thread counter-based / splittable PRNG** (e.g. xoshiro256\*\* or philox,
seeded per `(sweep, colour, block_id)`), giving independent, reproducible streams.

This is the one place the parallel kernel **cannot** match the serial RNG stream,
so the two will not be bit-identical — by design.

## 7. Validation strategy (since it's not bit-exact)

Run the serial fast kernel and the parallel kernel from the same initial state
for many steps across several seeds, and compare **distributions**, not single
trajectories:

- total energy: mean and histogram (KS test) at equilibrium
- acceptance ratio per move class
- radius of gyration / cluster-size distribution / radial density profile
- detailed-balance spot check: with all interactions off, the bead density must
  stay uniform (no spurious boundary depletion) — a direct test that the
  random-shift fix works.

`benchmark.py` already provides the harness shape; a `--parallel` mode would add
these statistical comparisons.

## 8. Implementation sketch (Cython + OpenMP)

```cython
from cython.parallel cimport prange, threadid

# Precompute, per colour c: a list of (block_id -> bead indices safe to move)
# using the footprint test (§3A) against the current random origin shift.

for sweep in range(n_sweeps):
    shift = random_vector()                      # §4.1 ergodicity
    for colour in range(8):                      # 8 colours in 3D
        blocks = active_blocks[colour]
        # each block is independent -> parallel, no GIL:
        for b in prange(len(blocks), nogil=True, schedule='dynamic'):
            rng = thread_rng[threadid()]
            run_block_substeps(blocks[b], grid, type_grid, idx_to_bead,
                               tables, angle_lookup, rng, &delta[threadid()])
        # barrier (implicit at end of prange)
        energy += sum(delta); reset(delta)
```

`run_block_substeps` is essentially the body already in `mega_crank_fast` —
restricted to a block's safe interior and using `rng` instead of `rand()`. The
energy/move/angle helpers are **already `noexcept nogil`**, so they drop straight
into the `prange` body. Build needs `extra_compile_args=['-fopenmp']` and
`extra_link_args=['-fopenmp']` (on Apple clang, `-Xpreprocessor -fopenmp` +
`libomp`).

## 9. Recommendation

1. **Ship the serial fast kernel now** — 10× kernel / ~5× whole-run, zero
   behaviour change (bit-exact). Low risk, big win.
2. **Treat checkerboard as a separate, opt-in feature** for large boxes and long
   runs. The kernel is already nogil-ready. Budget the work for: per-thread PRNG,
   block/colour bookkeeping with random shifts, footprint admissibility test for
   chains, and the statistical validation suite.
3. **Set expectations on scaling:** great for large systems, marginal for small
   boxes like the 30³ demo. If the science is dominated by many small-box runs,
   **embarrassingly-parallel independent replicas** (one process per core) may
   deliver better aggregate throughput than intra-run parallelism, with none of
   the detailed-balance subtlety.
