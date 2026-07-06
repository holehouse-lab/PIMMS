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

That is all that is required. ``PARALLELIZE`` is used whenever it is set, including
alongside a :doc:`freeze file <freeze>` (frozen beads are excluded from the movable
set but kept in place as fixed obstacles). ``PARALLEL_THREADS`` sets the number of
OpenMP threads; ``0`` means "use every core", and the keyword is ignored when
``PARALLELIZE`` is off.

Enabling ``PARALLELIZE`` only changes *which* Markov chain is followed, never the
target distribution, so it can never silently change the physics - at worst it has
no effect. (OpenMP must be available at build time; on macOS this means the Homebrew
``libomp`` package - without it the kernel simply runs single-threaded.)

How it works
============

The parallel kernels use a **frozen-halo domain decomposition**:

#. **Split the box into blocks.** The simulation box is divided into a grid of
   rectangular blocks (up to 4 per dimension). The block grid depends only on the
   box geometry, **not** on the thread count.

#. **Freeze a halo around each block.** Every block keeps a border of width ``W``
   frozen (no moves there). ``W = R_int + 2`` is chosen so that the full read+write
   footprint of any single move stays inside the block: ``R_int`` (the interaction
   range, 1 or 3) covers the energy evaluation and ``+2`` covers the bead
   displacement. Because the halos guarantee that two blocks' moves can **never
   touch the same lattice site**, the blocks are completely independent.

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
  dimension, so once the box exceeds ~``16 x W`` sites in a dimension the blocks
  simply grow as ``box / 4`` and the fixed ``2 x W`` frozen halo becomes a small
  fraction of each block - i.e. most of the system is movable each sweep. Small
  boxes (below ~``4 x W``) collapse to a single block and run serially; moderate
  boxes spend a large fraction of each sweep frozen in halos and benefit only
  partially. (``W = R_int + 2`` = 3 for short-range-only, 5 with long-range.)
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
the point.)

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
