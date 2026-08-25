.. _advanced-parallel:

===============
Parallelization
===============

For large systems (2D or 3D) the **crankshaft**, **slither** and **pull** moves can
be run on multi-threaded "checkerboard" kernels instead of the serial ones. Every
other move stays serial.

Quick start
===========

.. code-block:: text

   PARALLELIZE     : True
   PARALLEL_THREADS: 0          # 0 = use all available CPU cores

That is all that is required. With ``PARALLELIZE`` set, the crankshaft always
runs on the parallel kernel; the whole-chain moves (slither/pull) run on it
whenever every chain can fit a block interior, and otherwise automatically fall
back to the serial kernel for that megamove (see below) - so the flag never
changes the equilibrium being sampled. It composes with a :doc:`freeze file <freeze>`
(frozen beads are excluded from the movable set but kept in place as fixed
obstacles). ``PARALLEL_THREADS`` sets the number of
OpenMP threads; ``0`` means "use every core", and the keyword is ignored when
``PARALLELIZE`` is off.

Enabling ``PARALLELIZE`` only changes *which* Markov chain is followed, never the
target distribution, so it can never silently change the physics. (OpenMP must be
available at build time; on macOS this means the Homebrew ``libomp`` package -
without it the kernel simply runs single-threaded.)

.. warning::

   Same equilibrium, **different dynamics**. In every sweep only the beads inside
   the block interiors can move; everything in a frozen halo is scenery for that
   sweep, and a proposal that would land in a halo is rejected. The halos are
   re-drawn at random every sweep, so over a run every bead is moved equally often
   and the equilibrium distribution is exactly the serial one (verified by
   detailed-balance tests in the multi-block regime) - but the system **relaxes
   more slowly per step** than under the serial kernel, because each sweep
   perturbs only part of the system. Consequently a run that has **not** reached
   equilibrium (a collapsing or phase-separating system, or any run judged by
   its energy at a fixed step count) will show a different - less relaxed -
   energy than the serial run at the same step. Compare equilibrium averages,
   not energies at a fixed step, and equilibrate for longer when parallelizing.
   The fraction of the box movable per sweep is reported by
   ``mega_crank_fast.parallel_crank_layout_info(X, Y, Z, has_LR)``.

Startup report
==============

With ``PARALLELIZE`` set, PIMMS prints (and logs) a **parallelization report** at
startup - and again after a resized equilibration swaps in the production box,
since the decomposition depends on the box - describing exactly what this run
gets: the thread budget and whether the compiled kernels actually have OpenMP;
the crankshaft's halo width, block grid, block size and the fraction of the box
movable per sweep (or a warning that the box does not split and the kernel runs
single-threaded); and, for slither and pull, whether every chain fits a block
interior so the move really runs on the parallel kernel, or falls back to the
serial one and why. The same information is available programmatically from
``Simulation.report_parallelization()``,
``mega_crank_fast.parallel_crank_layout_info(X, Y, Z, has_LR)``,
``mega_crank_fast.openmp_info()`` and ``moves.parallel_chain_fit_report(...)``.

How it works
============

The parallel kernels use a **frozen-halo domain decomposition**:

#. **Split the box into blocks.** The simulation box is divided into a grid of
   rectangular blocks (up to 4 per dimension). The block grid depends only on the
   box geometry, **not** on the thread count.

#. **Freeze a halo around each block.** Every block keeps a border of width ``W``
   frozen (no moves there). Because a move is only ever *committed* inside the
   interior (see step 4), every write lies at least ``W`` inside a block face and
   every read (the energy shell of radius ``R_int`` = 1 or 3, the occupancy check
   of a candidate site, the bonded neighbours used by the angle term) reaches at
   most ``max(R_int, 2)`` sites past the interior edge, while the next block's
   interior starts ``2W`` past it. The crankshaft therefore uses the minimal
   race-free width ``W = 2`` with long-range interactions and ``W = 1`` without,
   with blocks kept at least ``8W`` long so at least 75% of each blocked
   dimension is movable. (The whole-chain slither/pull kernels use a wider
   chain-level halo, ``W = R_int + 2``.) Because the halos guarantee that two
   blocks' moves can **never touch the same lattice site**, the blocks are
   completely independent.

   .. note::

      **The halo is not the interaction range - it is half the separation
      between two regions that move at the same time.** Every block boundary
      has a halo on *both* sides, so adjacent interiors are ``2W`` frozen sites
      apart, and it is ``2W`` that must exceed a move's reach. With ``W = 2``
      and blocks of 20 along one axis::

         site:   14 15 16 17 | 18 19 | 20 21 | 22 23 24 25
                 -A interior-  A halo  B halo  -B interior-
                           ^                     ^
                      bead a (thread 1)     bead b (thread 2)

      Bead *a* on A's last interior site reads its energy shell (radius 3, which
      includes the distance-3 SLR layer) out to site 20 - all frozen halo, so an
      SLR partner sitting there is read correctly and cannot move meanwhile.
      Bead *b* reads down to 19. The two shells overlap on 19-20 but both only
      *read* there; *a* writes at or below 17, *b* at or above 22, and the two
      movers are 5 sites apart, beyond SLR range. So ``2W = 4 > 3`` is the whole
      condition; reads reaching *into* a halo are the design, not a violation.

#. **Run the blocks concurrently.** Each block is handed to an OpenMP thread and
   runs a batch of moves with **no locks** and its own independent random-number
   stream. Each block accumulates a private energy change; the global energy is the
   base energy plus the sum of the per-block deltas (an integer sum, so it is
   order-independent). Because the blocks are disjoint and deterministically
   seeded, the result is **bit-identical for any number of threads** - threads only
   change how fast the fixed set of blocks is processed.

#. **Shift the grid each sweep.** Beads/chains sitting in a frozen halo are skipped
   for that sweep. A fresh random origin shift is applied to the block grid on every
   call, so over many sweeps every part of the system spends time in a block
   interior and gets moved - restoring ergodicity. The move set is also kept
   "closed": a move whose result would leave the movable interior is rejected, which
   preserves detailed balance (a halo bead is never selected, so it could never make
   the reverse move).

The parallelized moves differ in the unit they decompose. The **crankshaft** is a
*per-bead* move with a tiny footprint, so the halo applies per bead: any bead at
least ``W`` inside its block is movable. The **slither** and **pull** are
*whole-chain* moves, so their decomposition is at the chain level: a chain is moved
only if **all of its beads** lie inside one block's interior; a chain straddling a
block boundary is frozen for that sweep. (Pull additionally restricts its
cooperative-reptation target search to the block interior, which keeps its
Metropolis-Hastings multiplicity correction self-consistent.)

If any chain's spatial extent could *never* fit inside a block interior - or a
chain exceeds the kernels' 512-bead per-chain buffer (heteropolymer chains for
slither, any chain for pull) - PIMMS detects this at dispatch time and runs that
megamove on the **serial kernel** instead, so no chain is ever silently frozen
on every sweep. You lose the speed-up for that move but never the equilibrium.

Which moves are parallelized
============================

The **crankshaft** (:doc:`/moves/crankshaft`, ``MOVE_CRANKSHAFT``), **slither**
(:doc:`/moves/slither`, ``MOVE_SLITHER``) and **pull** (:doc:`/moves/pull`,
``MOVE_PULL``) moves have parallel kernels, in **both 2D and 3D**:

.. list-table::
   :header-rows: 1
   :widths: 34 16 50

   * - Move
     - Dimensions
     - Parallel kernel
   * - Crankshaft (``MOVE_CRANKSHAFT``)
     - 2D / 3D
     - ``mega_crank_parallel_2D`` / ``mega_crank_parallel``
   * - Slither (``MOVE_SLITHER``)
     - 2D / 3D
     - ``mega_slither_parallel_2D`` / ``mega_slither_parallel``
   * - Pull (``MOVE_PULL``)
     - 2D / 3D
     - ``mega_pull_parallel_2D`` / ``mega_pull_parallel``
   * - all other moves
     - 2D & 3D
     - *(none - always serial)*

Every other move (chain translate/rotate/pivot, the cluster moves, the TSMMC moves,
jump-and-relax and VMMC) runs serially regardless of ``PARALLELIZE``. This is rarely
a limitation, because the crankshaft is the intended workhorse and normally
dominates the move budget.

When it helps
=============

``PARALLELIZE`` speeds a run up only when **all** of the following hold; otherwise it
has little or no effect (but never changes the physics). The relevant variable is
**box geometry, not density** - a dense, uniformly-filled melt in a large box
parallelizes well, whereas a small box does not regardless of how dilute it is.

* **The parallelized moves dominate the move set.** Time is only saved in proportion
  to the fraction of work spent in ``MOVE_CRANKSHAFT``, ``MOVE_SLITHER`` and
  ``MOVE_PULL`` (and their substep counts). A moveset that is mostly
  cluster/TSMMC/VMMC sees little benefit.
* **The box is large relative to the halo.** The block count is capped at 4 per
  dimension, so once the box exceeds ~``32 x W`` sites in a dimension the blocks
  simply grow as ``box / 4`` and the fixed ``2 x W`` frozen halo becomes a small
  fraction of each block - i.e. most of the system is movable each sweep. Boxes
  below ``8 x W`` in a dimension are not split in that dimension (a box that is
  below that in every dimension is a single block and runs serially). For the
  crankshaft ``W`` = 1 for short-range-only systems and 2 with long-range
  interactions (so boxes of 8 / 16 sites split); the whole-chain slither/pull
  kernels use ``W = R_int + 2`` = 3 / 5 and need ``4 x W`` = 12 / 20 sites.
* **Work is spread across the blocks** - and, for slither, **each chain's spatial
  extent is small compared to a block** (so it fits in a block interior; this is
  about chain size, not density - a collapsed chain is compact). A system that fills
  the box evenly (including a dense melt) gives all the threads balanced work. The
  bad case is a single concentrated droplet sitting in a big box: all the beads pile
  into a few blocks, leaving the other threads idle (a load-balance problem, not a
  density one).

In short: parallelization is most useful for **large boxes whose contents are spread
across the box** (dilute *or* dense), dominated by crankshaft and/or slither. It is
least useful for small boxes, a single concentrated droplet in a big box, chains
whose extent rivals the block size, or movesets that lean on the
collective/enhanced-sampling moves.

Measured speed-up
=================

The tables below benchmark the three parallel kernels on a 16-core machine, in 2D
with short-range interactions, on square boxes uniformly filled to ~7.5% with short
chains (4-bead ``AABB`` for crankshaft/slither, 6-bead ``AABBAB`` for pull, so chains
comfortably fit the block interiors). Each megamove performs the same number of move
attempts at every box size; the speed-up is the serial wall-time divided by the
parallel wall-time. (Absolute numbers are hardware-dependent, but the *trends* are
the point. The crankshaft rows were measured with the earlier, wider ``W = R_int
+ 2`` crankshaft halo; the current narrower halo changes the block layout at a
given box size but not the qualitative picture.)

.. list-table:: Slither (``mega_slither_parallel_2D``)
   :header-rows: 1
   :widths: 18 18 14 14

   * - Box
     - serial time
     - 4 threads
     - 8 threads
   * - 64 x 64
     - 10 ms
     - 2.6x
     - 2.9x
   * - 96 x 96
     - 23 ms
     - 4.1x
     - 6.8x
   * - 160 x 160
     - 68 ms
     - 4.0x
     - 6.5x
   * - 256 x 256
     - 174 ms
     - 4.0x
     - 7.1x
   * - 400 x 400
     - 434 ms
     - 4.3x
     - 8.1x

.. list-table:: Crankshaft (``mega_crank_parallel_2D``)
   :header-rows: 1
   :widths: 18 18 14 14

   * - Box
     - serial time
     - 4 threads
     - 8 threads
   * - 64 x 64
     - 2.8 ms
     - 3.3x
     - 3.0x
   * - 96 x 96
     - 6.6 ms
     - 3.5x
     - 6.0x
   * - 160 x 160
     - 20 ms
     - 3.9x
     - 5.9x
   * - 256 x 256
     - 54 ms
     - 3.9x
     - 6.7x
   * - 400 x 400
     - 134 ms
     - 4.2x
     - 7.8x

.. list-table:: Pull (``mega_pull_parallel_2D``)
   :header-rows: 1
   :widths: 18 18 14 14

   * - Box
     - serial time
     - 4 threads
     - 8 threads
   * - 64 x 64
     - 8.1 ms
     - 4.0x
     - 3.7x
   * - 96 x 96
     - 19 ms
     - 3.4x
     - 6.1x
   * - 160 x 160
     - 55 ms
     - 4.0x
     - 7.0x
   * - 256 x 256
     - 141 ms
     - 4.1x
     - 6.9x
   * - 400 x 400
     - 348 ms
     - 3.9x
     - 7.4x

All three moves show the same pattern predicted above: **small boxes scale poorly**
(the fixed halo dominates each block), and the speed-up climbs toward the thread
count as the box grows and the halo becomes a small fraction of each block. At the
largest box the scaling is essentially linear - at 1 / 2 / 4 / 8 threads the speed-up
is 1.1x / 2.2x / 4.2x / 8.1x for slither, 1.1x / 2.2x / 4.2x / 7.7x for crankshaft
and 1.1x / 2.1x / 3.9x / 7.4x for pull (the parallel kernel is even slightly faster
than serial on a single thread, from better cache locality). The crankshaft's serial
cost per megamove is roughly 3x lower than the whole-chain moves' for the same number
of attempts, because it is a cheap per-bead move rather than a reptation/cascade -
but all three parallelize equally well.

Checklist
=========

Reach for ``PARALLELIZE`` when:

* the box is large (comfortably more than ~``16 x W`` sites per dimension),
* the contents are spread across the box rather than balled up in one corner,
* crankshaft/slither/pull make up most of the move budget, and
* PIMMS was built with OpenMP available.

If any of those is not true, leave it off - it will not hurt correctness, but it will
not buy you much either.
