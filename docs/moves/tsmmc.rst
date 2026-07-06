.. _move-tsmmc:

=======================================
Temperature-switch Monte Carlo (TSMMC)
=======================================

:Keywords: ``MOVE_CTSMMC`` (9), ``MOVE_MULTICHAIN_TSMMC`` (10), ``MOVE_SYSTEM_TSMMC`` (12)
:Status: stable

TSMMC is one algorithm applied at three scopes, so the three move codes share this
page. All three are governed by the same ``TSMMC_*`` keywords (see
:doc:`/advanced/tsmmc`).

How it works
============

A TSMMC move takes part of the system on a **temperature excursion** to hop over
energy barriers that ordinary moves cannot cross. Starting from the simulation
temperature :math:`T`, the temperature is ramped up a schedule to a high "jump"
temperature ``TSMMC_JUMP_TEMP`` and back down again; at each temperature on the
schedule a burst of ordinary Monte Carlo moves is performed. At the high
temperature the system can climb out of a local energy minimum and explore; on the
way back down it re-settles, hopefully into a different basin. The whole excursion
is proposed as a single move and accepted or rejected as a unit.

The three variants differ only in what is heated:

* ``MOVE_CTSMMC`` (code 9) - a single, randomly chosen chain.
* ``MOVE_MULTICHAIN_TSMMC`` (code 10) - a random subset of chains.
* ``MOVE_SYSTEM_TSMMC`` (code 12) - the entire system (most powerful, most
  expensive).

Why detailed balance holds
==========================

A naïve "heat, move, cool, then accept on the final energy" scheme would **not**
be balanced, because the moves made at the elevated temperatures bias the
proposal. TSMMC instead uses **tempered transitions** (Neal, 1996), which restore
exact balance by accumulating the thermodynamic *work* done as the temperature is
switched.

Let the inverse-temperature schedule be
:math:`\beta_0 = \beta,\, \beta_1,\, \dots,\, \beta_{2n} = \beta` (with
:math:`\beta = 1/T`), rising to the jump temperature at the midpoint and returning,
and let :math:`x_k` be the configuration just before the :math:`k`-th temperature
change. Each segment runs ordinary, balanced MC moves at its own temperature. The
excursion is accepted with

.. math::

   A = \min\!\left(1,\; e^{W}\right), \qquad
   W = \sum_{k} \bigl(\beta_k - \beta_{k+1}\bigr)\, E(x_k),

where :math:`W` is the accumulated log-weight ("work") of the temperature
switching. Because the schedule is symmetric (it returns to :math:`\beta`) and the
sub-moves at each temperature are themselves balanced, this acceptance makes the
*whole* excursion satisfy detailed balance with respect to the target Boltzmann
distribution at :math:`T` - the elevated-temperature exploration is corrected for
exactly, so it changes the dynamics but not the sampled distribution. The PIMMS
implementation accumulates :math:`W` as :math:`\sum(\beta_\text{before} -
\beta_\text{after})\,E(x)` over the schedule and is checked against a
crankshaft-only reference by the detailed-balance test suite.

Configuration
=============

The move fractions ``MOVE_CTSMMC`` / ``MOVE_MULTICHAIN_TSMMC`` /
``MOVE_SYSTEM_TSMMC`` (all ``MOVE_*`` must sum to 1.0) select the variants. The
excursion itself is shaped by:

``TSMMC_JUMP_TEMP`` : float
    Peak temperature of the excursion. **Must exceed** ``TEMPERATURE`` (unless
    ``TSMMC_FIXED_OFFSET`` is used).

``TSMMC_NUMBER_OF_POINTS`` : int
    Number of temperature points on the ramp (more = smoother, more expensive).

``TSMMC_STEP_MULTIPLIER`` : int
    MC sub-steps performed at each temperature point.

``TSMMC_INTERPOLATION_MODE`` : str
    How temperatures are spaced; currently only ``LINEAR``.

``TSMMC_FIXED_OFFSET`` : float or False
    If set, the jump temperature is ``TEMPERATURE + TSMMC_FIXED_OFFSET`` rather
    than the absolute ``TSMMC_JUMP_TEMP`` (handy inside quench runs).

TSMMC is most useful for strongly-interacting systems that get stuck; it is
expensive (each move is many sub-moves across the schedule), so it is typically
used at a small fraction alongside the crankshaft.

For the full treatment - the temperature schedule, a step-by-step account of the
tempered-transitions work correction, a separate description of each of the three
variants (what is heated, which sub-moves run, and how the cost scales), and cost /
tuning guidance - see the dedicated :doc:`/advanced/tsmmc` page.
