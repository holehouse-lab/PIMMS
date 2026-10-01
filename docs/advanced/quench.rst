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
     - Master switch (default ``False``). When ``True``, all of the keywords
       below must be set - none of them has a default.
   * - ``QUENCH_START``
     - float
     - Temperature the run begins at (replaces ``TEMPERATURE``).
   * - ``QUENCH_END``
     - float
     - Target temperature the ramp finishes at.
   * - ``QUENCH_FREQ``
     - int
     - Number of steps between successive temperature changes. Must be
       positive.
   * - ``QUENCH_STEPSIZE``
     - float
     - Size of each temperature change. **Always positive** - the direction is
       inferred from ``START`` vs ``END`` (see below). A negative value is read
       as its magnitude.
   * - ``QUENCH_AS_EQUILIBRATION``
     - bool
     - If ``True``, the ramp *is* the equilibration phase (see below).

Two side effects happen as the keyfile is parsed. ``TEMPERATURE`` is overwritten
with ``QUENCH_START``, which is announced on stdout as ``[ WARNING ] : Resetting
the starting temperature to the QUENCH_START value ...`` when the two values
disagree (write ``TEMPERATURE : <QUENCH_START>`` and nothing is printed). With
``QUENCH_AS_EQUILIBRATION : True`` the ``EQUILIBRATION`` keyword is overwritten
with the quench window computed below, always announced as ``UPDATING
EQUILIBRATION TO [N] ...``.
Both keywords are still **required** in the keyfile even though the quench
replaces their values. PIMMS also
records ``QUENCH_END`` as the run's equilibrium temperature, which is what the
temperature-normalised angle penalty (``ANGLE_PENALTY_T_NORM``) and the
:doc:`lemonade </lemonade/index>` analysis layer use, so both are exact for the
constant-temperature production stretch and only approximate during the ramp.

Cooling vs heating
==================

You never specify a direction explicitly - PIMMS infers it from the endpoints:

* ``QUENCH_START > QUENCH_END`` → a **cooling** run (simulated annealing).
* ``QUENCH_START < QUENCH_END`` → a **heating** run (melting/dissolving).

``QUENCH_STEPSIZE`` is given as a positive magnitude in both cases; internally the
step is negated for a heating run so the temperature moves the right way, and the
start-up summary's ``QUENCH STEP`` line reports the magnitude you asked for. On every
step that is a multiple of ``QUENCH_FREQ`` the temperature is nudged by one
``QUENCH_STEPSIZE`` toward ``QUENCH_END``, and each change is announced on stdout
(``QUENCH: Updating temperature from 1.000 to 0.900``, with ``(target reached)``
appended when a step lands exactly on the target; suppressed by
``REDUCED_PRINTING``). When the next change would **overshoot** the target, the
temperature is clamped **exactly** to ``QUENCH_END`` instead, and that last change
is announced as ``QUENCH: Trying to update the temperature from ... Setting to
target temperature now...``. On the *next* multiple of ``QUENCH_FREQ`` after that PIMMS
notices the target has been reached, prints ``Reached target temperature of [...] -
no change`` (and logs ``Target temperature reached on step ...``) and switches the
quench off; the remainder of ``N_STEPS`` then runs at a constant ``QUENCH_END``. So a
quench always finishes with a stretch of ordinary fixed-temperature production at the
final temperature. (A fractional ramp is snapped to ``QUENCH_END`` once it is within
a relative ``1e-9`` of it, so accumulated floating-point error never leaves the run
one rung short of its target.)

Sizing the ramp
===============

Write the number of temperature changes ("rungs") the ramp needs as

.. math::

   R \;=\; \left\lceil \frac{|\,\text{QUENCH\_START} - \text{QUENCH\_END}\,|}
                            {\text{QUENCH\_STEPSIZE}} \right\rceil

(an integer ratio up to floating-point round-off counts as that integer rather than
being rounded up). Then:

* the ramp itself takes :math:`R \times \text{QUENCH\_FREQ}` steps - this is the
  ``Number of steps for quenching`` line in the start-up summary, and the number of
  rows written to ``QUENCH.dat``;
* the **quench window**, :math:`(1 + R) \times \text{QUENCH\_FREQ}` steps, is one
  ``QUENCH_FREQ`` longer, because the run has to spend one more full window *at*
  ``QUENCH_END`` before the quench is declared finished. This is the value
  ``QUENCH_AS_EQUILIBRATION`` uses for ``EQUILIBRATION``, and the value the
  fits-in-the-run check below uses. The two numbers differing by one
  ``QUENCH_FREQ`` is deliberate, not an off-by-one.

After the window the simulation continues at ``QUENCH_END`` for whatever is left of
``N_STEPS``. For a good anneal you generally want the ramp to be *slow* relative to
how quickly the system relaxes: prefer many small steps (small ``QUENCH_STEPSIZE``,
generous ``QUENCH_FREQ``) over a few large jumps, and leave enough steps after the
ramp for the system to equilibrate at the final temperature.

Startup constraints (all checked before the run begins, and all fatal):

* Every one of ``QUENCH_FREQ``, ``QUENCH_STEPSIZE``, ``QUENCH_START``,
  ``QUENCH_END`` and ``QUENCH_AS_EQUILIBRATION`` must be present in the keyfile.
* ``QUENCH_STEPSIZE`` must be positive.
* ``QUENCH_FREQ`` must be a positive number of steps.
* The ``START`` → ``END`` span must be at least one ``QUENCH_STEPSIZE``, otherwise
  a single step would overshoot the whole range.
* The quench window must be strictly shorter than ``N_STEPS``
  (:math:`(1 + R) \times \text{QUENCH\_FREQ} < \text{N\_STEPS}`), so the run always
  continues past the end of the ramp. A keyfile that fails this is rejected with
  "This quench will not complete ...".

.. note::

   A plain restart does **not** resume a quench. Started with ``RESTART_FILE``
   alone, a run takes only the configuration from ``restart.pimms``, numbers its
   steps from 1 and begins the whole ramp again at ``QUENCH_START``. To continue
   at the final temperature instead, restart with ``QUENCH_RUN : False`` and
   ``TEMPERATURE : <QUENCH_END>``. To resume the quench *exactly* where it
   stopped, mid-ramp included, use ``RESTART_CONTINUE : True``: the checkpoint
   records the step, the temperature in force and the generator states, and the
   resumed run picks the ramp up at that temperature and reproduces the
   uninterrupted run (see :ref:`restart-continue`). The continuation must keep the
   quench settings of the run it resumes: a checkpoint temperature that does not
   lie between this keyfile's ``QUENCH_START`` and ``QUENCH_END`` is refused (and
   with ``QUENCH_RUN : False`` a ``TEMPERATURE`` that differs from the checkpoint's
   is refused).

Using the ramp as equilibration
===============================

With ``QUENCH_AS_EQUILIBRATION : True`` the temperature ramp *is* the equilibration
phase: whatever you wrote for ``EQUILIBRATION`` is replaced by the quench window
computed above, and once the target is reached the temperature is held at
``QUENCH_END`` for the production phase. This is the natural choice for "anneal, then
measure": since no analysis runs during equilibration, everything in ``RG.dat``,
``CLUSTERS.dat`` and the rest comes from the constant-temperature production stretch
at ``QUENCH_END`` rather than from the non-equilibrium ramp.

With ``QUENCH_AS_EQUILIBRATION : False`` your own ``EQUILIBRATION`` value stands and
the ramp runs through the normal production accounting, so output written during the
ramp reflects the changing temperature. That is what you want if the ramp itself is
the measurement (a melting curve, say) rather than a way of preparing a state.

Output
======

Every temperature change is logged to ``QUENCH.dat``, one tab-separated row per
change (``step``, ``temperature``, ``energy``; no header line is written). The
temperature is printed to six significant digits and the energy right-aligned in a
ten-character field with four decimals, so a 1.0 → 0.7 ramp in steps of 0.1 with
``QUENCH_FREQ : 10`` writes three rows. These are from a real run (ten six-bead
chains; the energy column is whatever the run produced)

.. code-block:: text

   10	0.9	-3480.0000
   20	0.8	-3870.0000
   30	0.7	-4080.0000

Six significant digits means fractional ramps such as ``0.975`` are recorded exactly.
The file holds one row per rung (three here, for :math:`R = 3`); the step at which the
quench is switched off adds no row. ``QUENCH.dat`` is rewritten from scratch at the
start of a quench run, and deleted at the start of a non-quench run, so it never
mixes two runs' data.

The temperature is changed before the labelled Monte Carlo step and the energy
is recorded after that step's move, so each row describes the resulting state at
the displayed temperature.

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
base temperature during the ramp. When ``TSMMC_FIXED_OFFSET`` is set, the absolute
``TSMMC_JUMP_TEMP`` is ignored entirely.

Without ``TSMMC_FIXED_OFFSET``, a keyfile whose ramp would reach or exceed the
absolute ``TSMMC_JUMP_TEMP`` is rejected at parse time (the jump temperature must
exceed both ``QUENCH_START`` and ``QUENCH_END``, not just the starting temperature),
so an excursion can never be silently inverted mid-run - use ``TSMMC_FIXED_OFFSET``
or keep the hot end of the ramp below ``TSMMC_JUMP_TEMP``.

Worked example: anneal to assemble
==================================

Cool a mixture from a well-mixed hot state down to an assembly temperature, using
the ramp as equilibration and then measuring at the bottom:

.. code-block:: text

   DIMENSIONS      : 40 40 40
   PARAMETER_FILE  : params.prm
   CHAIN           : 200 AABB

   # TEMPERATURE and EQUILIBRATION are still REQUIRED keywords even though the
   # quench overwrites both (TEMPERATURE is reset to QUENCH_START, EQUILIBRATION
   # to the quench length with QUENCH_AS_EQUILIBRATION)
   TEMPERATURE     : 250
   EQUILIBRATION   : 0

   QUENCH_RUN      : True
   QUENCH_START    : 250
   QUENCH_END      : 30
   QUENCH_FREQ     : 500
   QUENCH_STEPSIZE : 2
   QUENCH_AS_EQUILIBRATION : True

   N_STEPS         : 4000000
   MOVE_CRANKSHAFT : 0.8
   MOVE_SLITHER    : 0.2

Here :math:`R = (250 - 30)/2 = 110` rungs, so the ramp itself is ``110 × 500 =
55000`` steps and the quench window - the value ``EQUILIBRATION`` is rewritten to -
is ``111 × 500 = 55500`` steps. The remaining ``3944500`` steps run as production at
``T = 30``, where all of the analysis is collected. The start-up summary reports both
numbers::

   Quench running as equilibration: TRUE
   Equilibration length: 55500
   Number of steps for quenching:  55000
   QUENCH FREQ  : 500
   QUENCH START : 250.00
   QUENCH STEP  : 2.00
   QUENCH END   : 30.00
