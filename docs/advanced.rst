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

For large, spatially **dispersed** systems (2D or 3D) the crankshaft move can be
run on a multi-threaded "checkerboard" kernel:

.. code-block:: text

   PARALLELIZE     : True
   PARALLEL_THREADS: 0          # 0 = use all available CPU cores

The parallel kernel partitions space into independent blocks (separated by a
frozen halo wide enough that no two blocks' moves can ever touch the same site)
and updates them concurrently with OpenMP threads. The decomposition depends only
on the box geometry, **not** on the thread count, so a run gives the identical
result for any number of threads. It targets the **same equilibrium distribution**
as the serial kernel but follows a different Markov chain.

It is used whenever ``PARALLELIZE`` is set and no chains are frozen; with a freeze
file present, PIMMS transparently falls back to the (bit-exact) serial kernel.
Enabling ``PARALLELIZE`` therefore can never silently change the physics - at worst
it has no effect (e.g. a box too small to decompose into more than one block). The
benefit is greatest for large dispersed boxes and negligible for small or collapsed
single-droplet systems. (OpenMP must be available at build time; on macOS this
means Homebrew ``libomp`` - otherwise the kernel runs single-threaded.)

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
