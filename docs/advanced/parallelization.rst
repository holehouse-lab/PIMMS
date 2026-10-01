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
   PARALLEL_THREADS: 0          # 0 = use every CPU available to the run

That is all that is required. Both keywords default to off / ``0``. With
``PARALLELIZE`` set, the crankshaft runs on the parallel kernel whenever the box
decomposes into at least two blocks (a one-block layout dispatches directly to the
serial kernel). The whole-chain moves (slither and pull) split their work
between the two kernels instead: the chains are partitioned once, by **chain
length**, and every megamove runs a parallel pass over the chains short enough to be
guaranteed to fit a block interior and a serial pass over the rest. Short
chains keep the speed-up, long chains keep moving. ``PARALLELIZE`` composes with a
:doc:`freeze file <freeze>` (frozen beads are excluded from the movable set but kept
in place as fixed obstacles, and do not constrain the split). ``PARALLEL_THREADS``
sets the number of OpenMP threads; ``0`` means "use every CPU available to the run":
the ``OMP_NUM_THREADS`` environment variable if it is set, otherwise the number of CPUs
the process is allowed to run on (its affinity mask, so a four-core batch allocation on
a large node starts four threads; on macOS, which exposes no affinity, this is the
machine's core count and includes the efficiency cores on Apple silicon). A negative
value, or one above 1024, is rejected at parse time, and the keyword has no effect
when ``PARALLELIZE`` is off.

The split is by chain length rather than by a chain's current shape on purpose, and
this matters more than it looks. Both kernels sample the same distribution on their
own, but *choosing* between them from the current configuration does not: the
parallel kernel's interior is closed (a move that would leave it is rejected), so a
compact chain can never be pulled out past the interior by it, while a serial
fallback would move freely - and a switch that flips on exactly that boundary pumps
probability into compact conformations. Chain length is a property of the *system*,
not of the configuration, so a length-based split is fixed for the whole run and both
passes leave the target distribution alone.

Enabling ``PARALLELIZE`` therefore only changes *which* Markov chain is followed,
never the target distribution. (OpenMP must be available at build time; on macOS this
means the Homebrew ``libomp`` package - without it the kernel simply runs
single-threaded, which the start-up report says plainly.)

.. warning::

   Same equilibrium, **different dynamics**. In every sweep only the beads inside
   the block interiors can move; everything in a frozen halo is scenery for that
   sweep, and a proposal that would land in a halo is rejected. The halos are
   re-drawn at random every sweep, so over a run every part of the box spends time
   in a block interior and the equilibrium distribution is exactly the serial one
   (verified by detailed-balance tests on the kernels in the multi-block regime, and
   for the whole-chain moves by an athermal dispatch-level test that drives
   ``system_slither`` and ``system_pull`` with ``parallelize`` on and off and compares
   the chain-extent and radius-of-gyration distributions) - but the system **relaxes
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
single-threaded); and, for slither and pull, the chain-level block layout, the
interior width, and how the chains sit relative to it. The same information is
available programmatically from ``Simulation.report_parallelization()``,
``mega_crank_fast.parallel_crank_layout_info(X, Y, Z, has_LR)``,
``mega_crank_fast.parallel_layout_info(X, Y, Z, has_LR)`` (the wider chain-level
halo used by slither and pull), ``mega_crank_fast.openmp_info()``,
``moves.parallel_chain_fit_report(...)`` and - the actual dispatch decision, one
boolean per chain - ``moves.parallel_chain_partition(...)``.

A 40 x 40 x 40 periodic box of 60 four-bead ``AABB`` chains with long-range beads,
``PARALLEL_THREADS : 4``, crankshaft + slither + pull in the move set and every
substep count left at its default reports the following (real output, as written to
``log.txt`` without its ``> STATUS: [time]:`` prefix; on stdout the same lines carry
a ``[STARTUP]:`` prefix and are wrapped)::

   PARALLELIZATION REPORT
   ----------------------
   Box: 40x40x40 (3D, periodic); 60 chains, 240 beads; interaction radius 3 (long-range beads present)
   Threads: 4 OpenMP threads per parallel megamove (from PARALLEL_THREADS; 16 CPUs available to this process)
   OpenMP: compiled in (runtime default thread budget 16)
   Crankshaft (MOVE_CRANKSHAFT): kernel mega_crank_parallel - per-bead frozen-halo checkerboard
       halo W=2; block grid 2x2x2 = 8 blocks of 20x20x20 sites; 51% of the box movable per sweep (random block shift every sweep)
   PARALLELIZE: megamove too small for the parallel kernel - at CRANKSHAFT_SUBSTEPS : 500, 4 threads can save at most about 0.038 ms per megamove (75% of 0.05 ms of sampling), against a fixed cost of about 0.12 ms for entering the parallel kernel (0.1 ms + 0.1 us per bead, 240 beads), so this move runs slower than it would without PARALLELIZE. Raise CRANKSHAFT_SUBSTEPS to about 20000 or more (break-even is near 2000), or drop PARALLELIZE. These costs are approximate and machine dependent.
   Slither (MOVE_SLITHER): kernel mega_slither_parallel - chain-level halo W=5, block grid 2x2x2 (20x20x20 sites, interiors 10x10x10); parallel kernel for 60 chain(s) of length <= 10, serial kernel for 0 longer chain(s) - both run every megamove
       (the split is fixed by chain length for the whole run - it is never re-decided from the configuration)
   Pull (MOVE_PULL): kernel mega_pull_parallel - chain-level halo W=5, block grid 2x2x2 (20x20x20 sites, interiors 10x10x10); parallel kernel for 60 chain(s) of length <= 10, serial kernel for 0 longer chain(s) - both run every megamove
       (the split is fixed by chain length for the whole run - it is never re-decided from the configuration)
   Note: the parallel kernels sample the SAME equilibrium as the serial ones but relax more slowly per step (only block interiors move each sweep) - judge equilibration by the observable plateau, not by step count.

The ``PARALLELIZE: megamove too small ...`` line is a warning, so in ``log.txt`` it
carries ``> WARNING: [time]:`` where the others carry ``> STATUS: [time]:``. It
appears under a move whose megamove is too small to pay for the parallel kernel's
fixed cost and names the ``*_SUBSTEPS`` keyword to raise and a value to raise it to
(the rule is given under :ref:`When it helps <advanced-parallel-when>`); here the
slither and pull megamoves, 60 chains at 10 substeps each, are large enough and get
no such line. Its companion is ``PARALLELIZE: only one thread will run the parallel
kernels ... set PARALLEL_THREADS above 1 or drop PARALLELIZE``, issued once when
``PARALLEL_THREADS`` resolves to a single thread, or ``PARALLELIZE: the kernels were
built without OpenMP ... rebuild with OpenMP or drop PARALLELIZE`` on a build that
has no OpenMP: either way every parallel megamove pays the fixed cost for no gain.
The ``Threads:`` line says ``from PARALLEL_THREADS`` for an explicit value and
``PARALLEL_THREADS : 0 -> OMP_NUM_THREADS if set, otherwise every available CPU``
for the default.

A move that is not in the move set is reported as ``not in the move set`` rather
than being described. The other lines you may see are a ``Frozen chains: ...
excluded from every parallel move, kept as fixed obstacles`` count when a
:doc:`freeze file <freeze>` is in use; ``OpenMP: NOT compiled into ...`` (blocks
run one after another, no speed-up) when OpenMP is missing; ``box does not split
... ONE block -> runs single-threaded`` for the crankshaft (the line gives the
``16 x W`` sites needed),
or ``box does not split for the chain-level halo ... every chain runs on the SERIAL
kernel`` for slither and pull (with the ``8 x W`` needed); and, under the slither or
pull line, how many serial-side chains exceed the 512-bead kernel buffer and (for
pull) how many chains are too short to pull at all. The report re-issued after a
resized equilibration is headed ``PARALLELIZATION REPORT - after resized
equilibration (production box)``.

The slither and pull lines tell you how the chains were split between the two
kernels and where the length threshold falls, so a run whose chains are all longer
than the interior will honestly report every chain on the serial side and no
speed-up for that move. The split is a property of the system, not of the
configuration the report happened to be printed from: see
:ref:`How it works <advanced-parallel-how>` below.

.. _advanced-parallel-how:

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
   interior starts ``2W + 1`` past it (so the reach must not exceed ``2W``). The
   crankshaft therefore uses the minimal
   race-free width ``W = 2`` with long-range interactions and ``W = 1`` without,
   with blocks kept at least ``8W`` long so at least 75% of each block's length
   along every split axis is movable (at that minimum that is 56% of a block's
   area in 2D and 42% of its volume in 3D; when the box does not divide evenly
   the few remainder sites are frozen too, so the movable fraction of an axis can
   dip slightly below 75%). (The whole-chain slither/pull kernels use a wider
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
      movers are 5 sites apart, beyond SLR range. So ``2W = 4 >= 3`` is the whole
      condition (without long-range beads ``2W = 2 >= 2``, the occupancy check
      and bonded-neighbour reach); reads reaching *into* a halo are the design,
      not a violation.

#. **Run the blocks concurrently.** Each block is handed to an OpenMP thread and
   runs a batch of moves with **no locks** and its own independent random-number
   stream (a splitmix64 state derived from the megamove's seed and the block index
   through a splitmix64 finalizer). The megamove's total attempt budget is shared out
   in proportion to each block's movable count (beads for the crankshaft, whole
   chains for slither and pull), with the floored shares topped
   up largest-remainder-first so exactly the requested number of attempts is made -
   unless the sweep's random shift leaves nothing movable at all (a small box, or a
   few chains that all straddle halos), in which case the sweep is a no-op and
   makes no attempts. The kernels report the attempts they actually made, and that
   is what ``MOVE_FREQS.dat`` counts, so such a sweep logs 0 attempts rather than
   the requested number. For slither and pull this means the parallel pass does
   **not** move each of its chains exactly ``SLITHER_SUBSTEPS`` (or
   ``PULL_SUBSTEPS``) times, as the serial pass does: it spends a total budget of
   that many attempts per chain on its side, shared among the blocks as above,
   and each block draws its chains uniformly *with replacement* from those that
   sit inside its interior on that sweep. A chain that straddles a halo on a given
   sweep gets no attempts, and the chains that are movable get proportionally
   more.
   Each block accumulates a private energy change; the global energy is the base
   energy plus the sum of the per-block deltas (an integer sum, so it is
   order-independent). Because the blocks are disjoint and deterministically
   seeded, the result is **bit-identical for any number of threads** - threads only
   change how fast the fixed set of blocks is processed. A block is the unit of
   work, so the kernels never ask OpenMP for more threads than there are blocks
   (at most 64): a ``PARALLEL_THREADS`` larger than the block count is cut to it
   inside the kernel, which changes nothing but the number of idle threads.

#. **Shift the grid each sweep.** Beads/chains sitting in a frozen halo are skipped
   for that sweep, as are any sitting in the trailing remainder when a box
   dimension is not an exact multiple of its block size. A fresh random origin
   shift is applied to the block grid on every
   call, so over many sweeps every part of the system spends time in a block
   interior and gets moved - restoring ergodicity. In every parallel kernel the
   shift is drawn uniformly over the whole box along each split axis, not over
   one block length. For the crankshaft that makes every site movable in the
   same fraction of sweeps, ``nb x (L - 2W) / box`` along an axis of ``nb``
   blocks of length ``L``, wherever the site sits; a shift confined to one block
   length left the frozen remainder of a non-divisible axis on the same stretch
   of the box most of the time, so a bead's sweep rate depended on its absolute
   coordinate (from 0.50 to 0.75 along a 26-site axis without long-range
   interactions). The equilibrium was never affected, only how evenly the box
   relaxed. For the whole-chain slither and pull it is what lets a chain as long
   as the interior be placed inside one from every starting position even when
   the box does not divide evenly into blocks. The
   move set is also kept
   "closed": a move whose result would leave the movable interior is rejected, which
   preserves detailed balance (a halo bead is never selected, so it could never make
   the reverse move).

Bit-identical across thread counts does **not** mean bit-identical to the serial
kernel: the parallel kernels draw from a different random stream and restrict every
sweep to the block interiors, so a parallel run and a serial run with the same
``SEED`` follow different trajectories. They target the same distribution, which is
what the detailed-balance tests check.

The parallelized moves differ in the unit they decompose. The **crankshaft** is a
*per-bead* move with a tiny footprint, so the halo applies per bead: any bead at
least ``W`` inside its block is movable. The **slither** and **pull** are
*whole-chain* moves, so their decomposition is at the chain level: a chain is moved
only if **all of its beads** lie inside one block's interior; a chain straddling a
block boundary is frozen for that sweep. (Pull additionally restricts its
cooperative-reptation target search to the block interior, which keeps its
Metropolis-Hastings multiplicity correction self-consistent.)

Because of that chain-level restriction, a chain **longer** than a block interior
could never fit one for any block shift, and would be silently frozen on every sweep
if it were handed to the parallel kernel. So slither and pull partition the chains
once, before any of them move:

* **short chains** - length at most the smallest block interior
  (``block_size - 2W`` on the split axes), and no longer than the kernels' 512-bead
  per-chain buffer (a limit that applies to heteropolymer chains for slither, to
  every chain for pull) - go to the **parallel kernel**;
* **everything else** (frozen chains apart, which neither kernel moves) goes to the
  **serial kernel**.

Every megamove then runs both passes, one after the other, over its own chains: the
chains held out of each pass stay in the grid as fixed obstacles for it, in exactly
the way a frozen chain does. Nothing is skipped, and the long chains keep reptating at
full serial speed. Whenever both passes have chains, which one goes first is decided
by a fair coin every megamove (with only one non-empty pass there is nothing to
order and no coin is drawn).
Each pass is reversible on its own, but the two do not commute (their chains
interact), so a fixed order would give a megamove that preserves the Boltzmann
distribution without being reversible; that is enough for the main loop but not for a
system-wide TSMMC excursion, whose acceptance assumes every sub-move is reversible.
The coin makes the megamove reversible.

The partition is deliberately a function of chain **length**, not of a chain's
current extent. Length is a run constant - as are the frozen set, the box and the
interaction range - so the split is fixed for the whole run (it is recomputed only
when the box itself changes, i.e. after a resized equilibration, which is a scheduled
event and not a property of the configuration). That is what makes the composition of
the two passes leave the target distribution alone. Deciding per megamove from the
current configuration would not: the parallel kernel's interior is closed, so the
rate at which a chain crosses out past the interior is identically zero through it,
while a serial fallback crosses freely - a switch that flips on that same boundary
over-samples compact conformations and biases the radius of gyration low by a few
percent, one-signed, with nothing to show for it in the energy or the move logs.
PIMMS did exactly that before 1.0.8; ``pimms/tests/test_parallel_dispatch.py`` now
pins both halves (the partition never moves, and the parallelized arm reproduces the
serial arm's chain-extent and ``Rg^2`` distributions in an athermal box).

The practical consequence for speed is that ``PARALLELIZE`` buys slither and pull
nothing for chains longer than a block interior. In a short-range box of 24-48 sites
per axis the interior is 6 to 11 sites, so anything past a ~10-mer is on the serial
side; with long-range beads a 40-80 site box gives interiors of 10 to 19. Enlarge the
box if you want the whole-chain moves to parallelize - the interior grows as
``box / 4 - 2W`` once the four-block cap is reached.

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
jump-and-relax and VMMC) runs serially regardless of ``PARALLELIZE``. The one
indirect exception is a system-wide TSMMC excursion: its sub-moves are ordinary
moves, so any crankshaft, slither or pull megamove it draws runs on the parallel
kernels exactly as in the main loop (the chain and multichain TSMMC moves use the
serial crankshaft kernel on their selected chains). This is rarely a limitation,
because the crankshaft is the intended workhorse and normally dominates the move
budget.

.. _advanced-parallel-when:

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
  dimension, so once the box exceeds about ``32 x W`` sites in a dimension the blocks
  simply grow as ``box / 4`` and the fixed ``2 x W`` frozen halo becomes a small
  fraction of each block - i.e. most of the system is movable each sweep. Boxes
  below ``16 x W`` in every dimension cannot form two crankshaft blocks and run
  serially (each crankshaft block is kept at least ``8 x W`` long). For the
  crankshaft ``W`` = 1 for short-range-only systems and 2 with long-range
  interactions (so two blocks require 16 / 32 sites along at least one axis);
  the whole-chain slither/pull kernels use ``W = R_int + 2`` = 3 / 5, keep each
  block at least ``4 x W`` long, and therefore require ``8 x W`` = 24 / 40 sites
  along at least one axis to split into two blocks.
* **The megamoves are big enough.** Launching the threads and bucketing the beads
  into blocks is a fixed cost per megamove, and the bookkeeping around every
  megamove (refreshing the bead table from the chains, drawing the bead selector,
  writing the moved positions back) is compiled but not threaded.
  A megamove has to do enough kernel work to amortise both: keep
  ``CRANKSHAFT_SUBSTEPS`` in the tens of thousands, and raise ``SLITHER_SUBSTEPS``
  / ``PULL_SUBSTEPS`` well above their default of 10 if those moves are to gain
  (see the measurements below).

  The start-up report checks this for you. Entering a parallel kernel costs
  roughly 0.1 ms plus 0.1 microseconds per bead before the first sub-move (about
  0.6 ms at 5,000 beads), which the serial kernel does not pay, and a sub-move
  costs roughly 0.1 microseconds for the crankshaft and 0.4 for a slither or pull
  of a short chain. With ``T`` threads (never more than there are blocks) the most
  a megamove can save is ``sub-moves x cost x (1 - 1/T)``, where the sub-moves are
  ``CRANKSHAFT_SUBSTEPS``, or ``SLITHER_SUBSTEPS`` / ``PULL_SUBSTEPS`` times the
  number of chains on the parallel side. When that is less than the fixed cost the
  parallel kernel is slower than the serial one, and the report adds a warning
  under that move naming the keyword and a value to raise it to (ten times the
  fixed cost in sampling work), for example::

     megamove too small for the parallel kernel - at CRANKSHAFT_SUBSTEPS : 500, 4 threads can save at most about 0.038 ms per megamove (75% of 0.05 ms of sampling), against a fixed cost of about 0.6 ms for entering the parallel kernel (0.1 ms + 0.1 us per bead, 5000 beads), so this move runs slower than it would without PARALLELIZE. Raise CRANKSHAFT_SUBSTEPS to about 60000 or more (break-even is near 8000), or drop PARALLELIZE. These costs are approximate and machine dependent.

  The message gives what the threads can save rather than the sampling time alone,
  because that is what has to beat the fixed cost: at two threads a megamove with
  0.3 ms of sampling can save at most 0.15 ms. At the default
  ``CRANKSHAFT_SUBSTEPS : 500`` the result is a crankshaft several times *slower*
  under ``PARALLELIZE``. The cost figures are approximate and machine dependent, so
  treat the suggested value as an order of magnitude. The rule uses only run
  constants (substeps, bead count, block layout, thread count); it is advice in the
  report and never changes which kernel runs. A single thread, or a build without
  OpenMP, gets its own warning, since it pays the fixed cost for no gain.
* **Work is spread across the blocks** - and, for slither and pull, **the chains are
  shorter than a block interior**, since that is what puts them on the parallel side
  of the length partition (this is about chain size, not density). A system that fills
  the box evenly (including a dense melt) gives all the threads balanced work. The
  bad case is a single concentrated droplet sitting in a big box: all the beads pile
  into a few blocks, leaving the other threads idle (a load-balance problem, not a
  density one).

In short: parallelization is most useful for **large boxes whose contents are spread
across the box** (dilute *or* dense), dominated by crankshaft and/or slither. It is
least useful for small boxes, a single concentrated droplet in a big box, chains
longer than a block interior (which slither and pull always move serially), or
movesets that lean on the collective/enhanced-sampling moves.

Measured speed-up
=================

The tables below were measured with ``pimms/fast_kernels/benchmark_parallel_2d.py``
(see the end of this section) on a 16-core development machine. The systems are 2D
with short-range interactions, on square boxes uniformly filled to ~7.5% with short
chains (4-bead ``AABB`` for crankshaft and slither, 6-bead ``AABBAB`` for pull, so
every chain is short enough for the parallel side of the length partition). They
report two things for each move, because they differ:

* **Wall time per megamove** through the normal ``PARALLELIZE`` path, i.e. what a
  run actually gains. Two megamove sizes are shown: the benchmark's *default* size
  (a crankshaft megamove of 50 000 attempts - ``CRANKSHAFT_SUBSTEPS`` in the 20-50K
  range usual for production, not the keyword's own default of 500; slither and pull
  at their keyword default of 10 substeps per chain), and a *heavy* megamove
  (500 000 crankshaft attempts; 100 substeps per chain).
* **Kernel time only**, for the heavy megamove: the compiled kernel by itself, with
  the per-megamove bookkeeping around it excluded.

The gap between the two is what PIMMS does around every megamove: refreshing the
bead table from the chain objects, drawing the bead selector, and writing the moved
positions back. It costs the same with or without ``PARALLELIZE``, grows with the
number of beads, and is not threaded, but since those copies were compiled
(``pimms.bookkeeping``) it is about a millisecond at 12 000 beads rather than the
five to eleven it cost as per-chain Python loops, which is what used to cap the
wall-time speed-up well below the kernels'. Speed-ups are the serial wall (or
kernel) time divided by the parallel one; each number is the median of seven
megamoves. Absolute times are hardware-dependent and even the kernel-only ratios
move by tens of percent between runs on a busy machine; on a quiet machine the
trends are the point. Under heavy load even the trends blur: a re-run of the script
on a 16-core machine busy with other jobs (load average about 10) gave serial times
10-70% longer and mostly lower speed-ups, and while the default-size slither and pull
megamoves still ran *slower* than serial in the 64 x 64 box, the crankshaft and pull
kernel-only speed-ups no longer climbed steadily with box size.

.. list-table:: Crankshaft (``mega_crank_parallel_2D``)
   :header-rows: 2
   :widths: 14 14 12 14 12 12 16

   * - Box
     - default megamove
     -
     - heavy megamove
     -
     -
     - kernel only
   * -
     - serial
     - 8 threads
     - serial
     - 4 threads
     - 8 threads
     - 8 threads, heavy
   * - 64 x 64
     - 3.5 ms
     - 2.6x
     - 34 ms
     - 2.7x
     - 3.5x
     - 5.3x
   * - 96 x 96
     - 3.4 ms
     - 2.8x
     - 33 ms
     - 2.8x
     - 3.9x
     - 5.7x
   * - 160 x 160
     - 3.6 ms
     - 3.3x
     - 34 ms
     - 3.5x
     - 5.1x
     - 6.1x
   * - 256 x 256
     - 4.2 ms
     - 2.4x
     - 38 ms
     - 2.9x
     - 3.8x
     - 6.0x
   * - 400 x 400
     - 5.0 ms
     - 1.8x
     - 38 ms
     - 2.9x
     - 4.1x
     - 6.3x

.. list-table:: Slither (``mega_slither_parallel_2D``)
   :header-rows: 2
   :widths: 14 14 12 14 12 12 16

   * - Box
     - default megamove
     -
     - heavy megamove
     -
     -
     - kernel only
   * -
     - serial
     - 8 threads
     - serial
     - 4 threads
     - 8 threads
     - 8 threads, heavy
   * - 64 x 64
     - 0.21 ms
     - 0.7x
     - 1.8 ms
     - 2.3x
     - 2.4x
     - 3.1x
   * - 96 x 96
     - 0.46 ms
     - 1.0x
     - 4.1 ms
     - 2.8x
     - 3.4x
     - 4.6x
   * - 160 x 160
     - 1.3 ms
     - 1.7x
     - 11 ms
     - 3.0x
     - 4.2x
     - 6.4x
   * - 256 x 256
     - 3.2 ms
     - 2.2x
     - 29 ms
     - 3.2x
     - 4.6x
     - 7.3x
   * - 400 x 400
     - 8.3 ms
     - 2.4x
     - 72 ms
     - 3.3x
     - 4.7x
     - 7.7x

.. list-table:: Pull (``mega_pull_parallel_2D``)
   :header-rows: 2
   :widths: 14 14 12 14 12 12 16

   * - Box
     - default megamove
     -
     - heavy megamove
     -
     -
     - kernel only
   * -
     - serial
     - 8 threads
     - serial
     - 4 threads
     - 8 threads
     - 8 threads, heavy
   * - 64 x 64
     - 0.18 ms
     - 0.6x
     - 1.5 ms
     - 2.1x
     - 2.1x
     - 2.3x
   * - 96 x 96
     - 0.40 ms
     - 0.9x
     - 3.3 ms
     - 2.5x
     - 3.1x
     - 4.0x
   * - 160 x 160
     - 1.1 ms
     - 1.6x
     - 9.4 ms
     - 2.9x
     - 3.9x
     - 5.0x
   * - 256 x 256
     - 2.7 ms
     - 2.0x
     - 24 ms
     - 2.9x
     - 4.3x
     - 6.0x
   * - 400 x 400
     - 7.0 ms
     - 2.4x
     - 60 ms
     - 3.1x
     - 4.8x
     - 7.0x

Three things to read off the tables.

* **The kernels scale.** Small boxes scale poorly (the fixed halo is a large
  fraction of each block) and the kernel-only speed-up climbs toward the thread
  count as the box grows: at 400 x 400 the heavy megamove's kernel runs 1.0x / 2.0x /
  3.7x / 6.3x faster at 1 / 2 / 4 / 8 threads for the crankshaft, 1.2x / 2.2x / 4.1x /
  7.7x for the slither and 1.1x / 2.0x / 3.9x / 7.0x for the pull (only the 8-thread
  kernel figures appear in the tables; the whole-chain parallel kernels are slightly
  faster than serial even on one thread).
* **What sits outside the kernel is now small, but it is not zero.** Around the
  kernel, each megamove refreshes the bead table, draws the bead selector and
  writes the moved positions back; at 400 x 400 (12 000 beads) that is about 1 ms
  for the crankshaft and 1.5 ms for the whole-chain moves, serial or parallel, now
  that the copies are compiled (they used to be 5 to 11 ms as Python loops). With
  heavy megamoves the kernel dominates and the wall-time speed-up reaches 4.1x to
  4.8x on 8 threads, against 6.3x to 7.7x for the kernels alone; with default-size
  megamoves the kernel does only a few milliseconds of work, the fixed costs are a
  comparable fraction, and the gain is 1.8x to 2.4x.
* **Small megamoves do not pay.** Launching the threads and bucketing the beads
  into blocks is a fixed cost of a few tenths of a millisecond. A slither or pull
  megamove at the default 10 substeps per chain in a 64 x 64 box does only 0.2 ms of
  kernel work, so ``PARALLELIZE`` makes it *slower* (0.6x to 0.7x). Raise
  ``SLITHER_SUBSTEPS`` / ``PULL_SUBSTEPS`` (and keep ``CRANKSHAFT_SUBSTEPS`` in the
  tens of thousands) if you want those moves to benefit.

To re-measure on your own hardware run
``python pimms/fast_kernels/benchmark_parallel_2d.py`` from the repository root
(it reproduces these tables and checks the incrementally tracked energy against a
from-scratch recompute after every configuration); ``python
pimms/fast_kernels/benchmark_parallel.py`` is the 3D correctness and scaling harness
for the crankshaft kernel.

Checklist
=========

Reach for ``PARALLELIZE`` when:

* the box is large (comfortably more than about ``16 x W`` sites per dimension),
* the contents are spread across the box rather than balled up in one corner,
* crankshaft/slither/pull make up most of the move budget, and
* PIMMS was built with OpenMP available.

If any of those is not true, leave it off - it will not hurt correctness, but it will
not buy you much either.
