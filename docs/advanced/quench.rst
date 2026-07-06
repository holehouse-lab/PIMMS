.. _advanced-quench:

============================
Quench / simulated annealing
============================

A **quench** run changes the temperature *during* the simulation rather than
holding it fixed. The usual purpose is **simulated annealing** - start hot so the
system can explore freely, then cool slowly so it settles into a low-energy,
well-equilibrated state (assembled droplets, folded structures, ordered phases).
The same machinery runs in reverse, so you can also **heat** a configuration to
melt or dissolve it.

Because a quench deliberately walks the temperature along a schedule, it is one of
the few features here that *does* change the ensemble being sampled - that is the
point. Within any single temperature window the moves are the ordinary
detailed-balance moves; the quench simply resets the temperature every so often.

Turning it on
=============

Set ``QUENCH_RUN : True`` and provide the full set of quench keywords. When a
quench is active the plain ``TEMPERATURE`` keyword is **ignored** - the starting
temperature comes from ``QUENCH_START`` instead:

.. code-block:: text

   QUENCH_RUN      : True
   QUENCH_START    : 200      # initial temperature (TEMPERATURE is ignored)
   QUENCH_END      : 40       # final temperature
   QUENCH_FREQ     : 100      # change the temperature every 100 steps
   QUENCH_STEPSIZE : 5        # by this much each change (always a positive number)
   QUENCH_AS_EQUILIBRATION : False

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Keyword
     - Type
     - Meaning
   * - ``QUENCH_RUN``
     - bool
     - Master switch. When ``True``, all of the keywords below must be set.
   * - ``QUENCH_START``
     - float
     - Temperature the run begins at (replaces ``TEMPERATURE``).
   * - ``QUENCH_END``
     - float
     - Target temperature the ramp finishes at.
   * - ``QUENCH_FREQ``
     - int
     - Number of steps between successive temperature changes.
   * - ``QUENCH_STEPSIZE``
     - float
     - Size of each temperature change. **Always positive** - the direction is
       inferred from ``START`` vs ``END`` (see below).
   * - ``QUENCH_AS_EQUILIBRATION``
     - bool
     - If ``True``, the ramp *is* the equilibration phase (see below).

Cooling vs heating
==================

You never specify a direction explicitly - PIMMS infers it from the endpoints:

* ``QUENCH_START > QUENCH_END`` → a **cooling** run (simulated annealing).
* ``QUENCH_START < QUENCH_END`` → a **heating** run (melting/dissolving).

``QUENCH_STEPSIZE`` is given as a positive magnitude in both cases; internally the
step is negated for a heating run so the temperature moves the right way. Every
``QUENCH_FREQ`` steps the temperature is nudged by one ``QUENCH_STEPSIZE`` toward
``QUENCH_END``. When the next step would **overshoot** the target, the temperature
is clamped **exactly** to ``QUENCH_END`` instead, and from that point on the quench
flag is switched off and the remainder of ``N_STEPS`` runs at a constant
``QUENCH_END`` - so a quench always finishes with a stretch of ordinary
fixed-temperature production at the final temperature.

Sizing the ramp
===============

The ramp occupies roughly

.. math::

   \text{ramp steps} \;\approx\; \frac{|\,\text{QUENCH\_START} - \text{QUENCH\_END}\,|}
                                        {\text{QUENCH\_STEPSIZE}} \times \text{QUENCH\_FREQ}

steps, after which the simulation continues at ``QUENCH_END`` for whatever is left
of ``N_STEPS``. For a good anneal you generally want the ramp to be *slow* relative
to how quickly the system relaxes: prefer many small steps (small
``QUENCH_STEPSIZE``, generous ``QUENCH_FREQ``) over a few large jumps, and leave
enough steps after the ramp for the system to equilibrate at the final temperature.

Startup constraints (all checked before the run begins):

* ``QUENCH_STEPSIZE`` must be positive.
* The ``START`` → ``END`` span must be at least one ``QUENCH_STEPSIZE``.
* The ramp must fit inside ``N_STEPS`` (it cannot request more steps than the run
  has).

Using the ramp as equilibration
===============================

With ``QUENCH_AS_EQUILIBRATION : True`` the temperature ramp *is* the equilibration
phase. The equilibration period is sized to cover the ramp, and once the target is
reached the temperature is held at ``QUENCH_END`` for the production phase. This is
the natural choice for "anneal, then measure": all of your analysis then comes from
the constant-temperature production stretch at ``QUENCH_END`` rather than from the
non-equilibrium ramp.

With ``QUENCH_AS_EQUILIBRATION : False`` the ramp runs through the normal
production accounting, so output written during the ramp reflects the changing
temperature.

Output
======

Every temperature change is logged to ``QUENCH.dat``, one tab-separated row per
change (``step``, ``temperature``, ``energy``; no header line is written):

.. code-block:: text

   100       195.00         -1423.0000
   200       190.00         -1502.0000
   ...

This lets you plot the energy against temperature directly - the classic view for
spotting a transition (a sharp drop in energy, or a peak in its fluctuations, as
the system assembles on cooling).

Interaction with TSMMC
======================

If you combine a quench with the :doc:`TSMMC <tsmmc>` moves, the TSMMC coordinator
is rebuilt at the new base temperature each time the quench updates - so the
temperature excursions always heat *relative to the current* simulation
temperature. The ``TSMMC_FIXED_OFFSET`` keyword is especially convenient here,
because it defines the jump temperature as an offset above the current temperature
rather than as a fixed absolute value that might fall below the (falling or rising)
base temperature during the ramp.

Worked example: anneal to assemble
==================================

Cool a mixture from a well-mixed hot state down to an assembly temperature, using
the ramp as equilibration and then measuring at the bottom:

.. code-block:: text

   DIMENSIONS      : 40 40 40
   PARAMETER_FILE  : params.prm
   CHAIN           : 200 AABB

   QUENCH_RUN      : True
   QUENCH_START    : 250
   QUENCH_END      : 30
   QUENCH_FREQ     : 500
   QUENCH_STEPSIZE : 2
   QUENCH_AS_EQUILIBRATION : True

   N_STEPS         : 4000000
   MOVE_CRANKSHAFT : 0.8
   MOVE_SLITHER    : 0.2

Here the ramp spends ``(1 + (250-30)/2) × 500 = 55500`` steps cooling from 250 to 30
as the equilibration phase (the ``1 +`` accounts for the initial temperature point),
then the remaining steps run as production at ``T = 30``, where the analysis is
collected.
