.. _advanced:

=================
Advanced features
=================

Beyond the basic move set and analysis, PIMMS offers several features for harder
sampling problems and specialised workflows. Several are **experimental** and
require ``EXPERIMENTAL_FEATURES : True``; these are flagged below.

.. _advanced-quench:

Quench / simulated annealing
============================

A **quench** run ramps the temperature during the simulation rather than holding
it fixed - typically cooling a system to drive assembly (simulated annealing), but
heating is supported too. Enable it with ``QUENCH_RUN : True`` and provide the
full set of quench keywords:

.. code-block:: text

   QUENCH_RUN      : True
   QUENCH_START    : 200      # initial temperature (TEMPERATURE is ignored)
   QUENCH_END      : 40       # final temperature
   QUENCH_FREQ     : 100      # change temperature every 100 steps
   QUENCH_STEPSIZE : 5        # by this much each time (always positive)
   QUENCH_AS_EQUILIBRATION : False

When ``QUENCH_START > QUENCH_END`` PIMMS cools; otherwise it heats. The
temperature steps by ``QUENCH_STEPSIZE`` every ``QUENCH_FREQ`` steps and is clamped
exactly to ``QUENCH_END`` (the run then continues at ``QUENCH_END``). If
``QUENCH_AS_EQUILIBRATION : True``, the ramp *is* the equilibration phase
(equilibration is sized to cover the ramp, after which production runs at
``QUENCH_END``).

Constraints (checked at startup): ``QUENCH_STEPSIZE`` must be positive, the
``START``→``END`` span must be at least one step size, and the ramp cannot exceed
``N_STEPS``. The temperature trajectory and the energy response are written to
``QUENCH.dat``.

.. _advanced-freeze:

Freeze files
============

A **freeze file** holds chosen chains rigidly fixed for the entire simulation -
useful for, e.g., a fixed scaffold or surface that mobile chains interact with.
Point the ``FREEZE_FILE`` keyword at a plain-text file listing chainIDs to freeze:

.. code-block:: text

   # freeze.txt  -  lines beginning with # are comments
   C 0 1 2          # freeze chainIDs 0, 1 and 2
   C 10             # more C lines are allowed

Frozen chains never move, but they **still contribute to the energy** (other
chains feel them normally); they are simply skipped when PIMMS picks beads to move.
To discover which chainID is which, run once with ``WRITE_CHAIN_TO_CHAINID : True``
and read ``chain_to_chainid.txt``.

.. note::

   Freezing forces PIMMS to use the serial crankshaft kernel, so ``PARALLELIZE``
   has no effect while any chain is frozen. (Per-bead freezing - a ``B`` directive
   - is reserved but not yet implemented.)

.. _advanced-parallel:

Parallelization
===============

For large systems (2D or 3D) the crankshaft and slither moves can be run on
multi-threaded "checkerboard" kernels:

.. code-block:: text

   PARALLELIZE     : True
   PARALLEL_THREADS: 0          # 0 = use all available CPU cores

It is used whenever ``PARALLELIZE`` is set and no chains are frozen; with a freeze
file present, PIMMS transparently falls back to the (bit-exact) serial kernel.
Enabling ``PARALLELIZE`` therefore can never silently change the physics - at worst
it has no effect. (OpenMP must be available at build time; on macOS this means
Homebrew ``libomp`` - otherwise the kernel runs single-threaded.)

How it works
------------

The parallel kernels use a **frozen-halo domain decomposition**:

#. **Split the box into blocks.** The simulation box is divided into a grid of
   rectangular blocks (up to 4 per dimension). The block grid depends only on the
   box geometry, **not** on the thread count.

#. **Freeze a halo around each block.** Every block keeps a border of width
   ``W`` frozen (no moves there). ``W = R_int + 2`` is chosen so that the full
   read+write footprint of any single move stays inside the block: ``R_int`` (the
   interaction range, 1 or 3) covers the energy evaluation and ``+2`` covers the
   bead displacement. Because the halos guarantee that two blocks' moves can
   **never touch the same lattice site**, the blocks are completely independent.

#. **Run the blocks concurrently.** Each block is handed to an OpenMP thread and
   runs a batch of moves with **no locks** and its own independent random-number
   stream. Each block accumulates a private energy change; the global energy is
   the base energy plus the sum of the per-block deltas (an integer sum, so it is
   order-independent). Because the blocks are disjoint and deterministically
   seeded, the result is **bit-identical for any number of threads** - threads
   only change how fast the fixed set of blocks is processed.

#. **Shift the grid each sweep.** Beads/chains sitting in a frozen halo are skipped
   for that sweep. A fresh random origin shift is applied to the block grid on
   every call, so over many sweeps every part of the system spends time in a
   block interior and gets moved - restoring ergodicity. The move set is also kept
   "closed": a move whose result would leave the movable interior is rejected,
   which preserves detailed balance (a halo bead is never selected, so it could
   never make the reverse move).

The two parallelized moves differ in the unit they decompose. The **crankshaft**
is a *per-bead* move with a tiny footprint, so the halo applies per bead: any bead
at least ``W`` inside its block is movable. The **slither** is a *whole-chain*
move, so the decomposition is at the chain level: a chain is moved only if **all
of its beads** lie inside one block's interior; a chain straddling a block
boundary is frozen for that sweep.

Which moves are parallelized
----------------------------

The **crankshaft** (:doc:`/moves/crankshaft`, ``MOVE_CRANKSHAFT``) and **slither**
(:doc:`/moves/slither`, ``MOVE_SLITHER``) moves have parallel kernels, in **both 2D
and 3D**:

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
   * - all other moves
     - 2D & 3D
     - *(none - always serial)*

Every other move (chain translate/rotate/pivot, pull, the cluster moves, the TSMMC
moves, jump-and-relax and VMMC) runs serially regardless of ``PARALLELIZE``. This
is rarely a limitation, because the crankshaft is the intended workhorse and
normally dominates the move budget.

The two moves parallelize differently because of what they touch. The crankshaft
is a *per-bead* move with a tiny footprint, so it uses a per-bead frozen-halo
decomposition. The slither is a *whole-chain* move, so it uses a **chain-level**
decomposition: a chain is moved in parallel only if **all of its beads fit inside
one block's interior**. A chain that straddles a block boundary is simply frozen
for that sweep (the next sweep's random origin shift gives it another chance) - so,
as for the crankshaft, enabling ``PARALLELIZE`` never changes the physics.

When it helps
-------------

``PARALLELIZE`` speeds a run up only when **all** of the following hold; otherwise
it has little or no effect (but never changes the physics). Note that the relevant
variable is **box geometry, not density** - a dense, uniformly-filled melt in a
large box parallelizes well, whereas a small box does not regardless of how dilute
it is.

* **The parallelized moves dominate the move set.** Time is only saved in
  proportion to the fraction of work spent in ``MOVE_CRANKSHAFT`` and
  ``MOVE_SLITHER`` (and their substep counts). A moveset that is mostly
  pull/cluster/TSMMC/VMMC sees little benefit.
* **The box is large relative to the halo.** The block count is capped at 4 per
  dimension, so once the box exceeds ~``16 x W`` sites in a dimension the blocks
  simply grow as ``box / 4`` and the fixed ``2 x W`` frozen halo becomes a small
  fraction of each block - i.e. most of the system is movable each sweep. Small
  boxes (below ~``4 x W``) collapse to a single block and run serially; moderate
  boxes spend a large fraction of each sweep frozen in halos and benefit only
  partially. (``W = R_int + 2`` = 3 for short-range-only, 5 with long-range.)
* **Work is spread across the blocks** - and, for slither, **each chain's spatial
  extent is small compared to a block** (so it fits in a block interior; this is
  about chain size, not density - a collapsed chain is compact). A system that
  fills the box evenly (including a dense melt) gives all the threads balanced
  work. The bad case is a single concentrated droplet sitting in a big box: all
  the beads pile into a few blocks, leaving the other threads idle (a
  load-balance problem, not a density one).

In short: parallelization is most useful for **large boxes whose contents are
spread across the box** (dilute *or* dense), dominated by crankshaft and/or
slither. It is least useful for small boxes, a single concentrated droplet in a
big box, chains whose extent rivals the block size, or movesets that lean on the
collective/enhanced-sampling moves.

Measured speed-up
-----------------

The tables below benchmark the two parallel kernels on a 16-core machine, in 2D
with short-range interactions, on square boxes uniformly filled to ~7.5% with
4-bead ``AABB`` chains (so chains comfortably fit the block interiors). Each
megamove performs the same number of move attempts at every box size; the speed-up
is the serial wall-time divided by the parallel wall-time. (Absolute numbers are
hardware-dependent, but the *trends* are the point.)

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

Both moves show the same pattern predicted above: **small boxes scale poorly** (the
fixed halo dominates each block), and the speed-up climbs toward the thread count
as the box grows and the halo becomes a small fraction of each block. At the
largest box the scaling is essentially linear - at 1 / 2 / 4 / 8 threads the
speed-up is 1.1x / 2.2x / 4.2x / 8.1x for slither and 1.1x / 2.2x / 4.2x / 7.7x for
crankshaft (the parallel kernel is even slightly faster than serial on a single
thread, from better cache locality). The crankshaft's serial cost per megamove is
roughly 3x lower than slither's for the same number of attempts, because it is a
cheap per-bead move rather than a whole-chain reptation - but both parallelize
equally well.

.. _advanced-tsmmc:

Temperature-switch Monte Carlo (TSMMC)
======================================

**TSMMC** moves help a simulation hop over energy barriers that ordinary moves
cannot cross. The idea is a temperature *excursion*: part of the system is briefly
heated along a schedule up to a high "jump" temperature and cooled back down,
during which it can rearrange; a **tempered-transitions** acceptance criterion
(accumulating the work done over the temperature schedule) keeps the overall move
detailed-balanced, so the sampled distribution is unchanged.

There are three variants, differing in what is heated:

* ``MOVE_CTSMMC`` (code 9) - a single randomly chosen chain.
* ``MOVE_MULTICHAIN_TSMMC`` (code 10) - a random subset of chains.
* ``MOVE_SYSTEM_TSMMC`` (code 12) - the entire system (the most powerful and the
  most expensive).

The excursion is configured by the ``TSMMC_*`` keywords:

.. code-block:: text

   MOVE_CRANKSHAFT       : 0.9
   MOVE_SYSTEM_TSMMC     : 0.1
   TSMMC_JUMP_TEMP       : 120     # peak temperature (must exceed TEMPERATURE)
   TSMMC_NUMBER_OF_POINTS: 20      # temperature points in the schedule (smoother = more)
   TSMMC_STEP_MULTIPLIER : 50      # MC sub-steps per temperature point
   TSMMC_INTERPOLATION_MODE : LINEAR
   EXPERIMENTAL_FEATURES : True

``TSMMC_JUMP_TEMP`` must be greater than the simulation ``TEMPERATURE`` (unless
``TSMMC_FIXED_OFFSET`` is used, which instead sets the jump temperature relative to
the current temperature - convenient inside quench runs). More schedule points and
a larger step multiplier make each excursion gentler and more thorough, but more
expensive. TSMMC is experimental.

.. _advanced-other:

Other features
==============

**Collective moves for dense systems.** The pull move (``MOVE_PULL``) and
Virtual-Move Monte Carlo (``MOVE_VMMC``) are designed to rearrange dense/condensed
phases where rigid single-chain moves clash or get kinetically trapped; the
slither move (``MOVE_SLITHER``) efficiently relaxes chain conformations. See
:ref:`overview-moves` for the move set, and note that ``MOVE_PULL`` and
``MOVE_VMMC`` are experimental.

**Reference and debugging controls.**

* ``NON_INTERACTING : True`` - zero **all** interaction energies and run a pure
  excluded-volume reference ensemble (the parameter file is still required but its
  energies are ignored).
* ``ANGLES_OFF : True`` - disable backbone-angle penalties (no ``ANGLE_PENALTY``
  lines needed in the parameter file).
* ``ENERGY_CHECK : <freq>`` - periodically recompute the total energy from scratch
  and abort if it has drifted from the incrementally tracked value. An O(N)
  safety/debugging check; cheap occasionally, costly every step.

**Box and trajectory controls.**

* ``RESIZED_EQUILIBRATION`` / ``EQUILIBRATION_OFFSET`` - equilibrate in a smaller
  box and grow to the production ``DIMENSIONS`` (requires ``HARDWALL``); see
  :ref:`overview-setup`.
* ``AUTOCENTER : True`` - keep a single chain centred in the box every frame.
* ``CASE_INSENSITIVE_CHAINS : False`` - treat upper- and lower-case bead letters
  as distinct types.

**Experimental gate.** ``EXPERIMENTAL_FEATURES : True`` unlocks the
not-yet-stable keywords/moves (``EXTRA_CHAIN``, ``MOVE_PULL``, ``MOVE_VMMC``, the
TSMMC moves, ``MOVE_JUMP_AND_RELAX``). None of these are guaranteed to work in
every configuration - leave the gate ``False`` unless you need them.
