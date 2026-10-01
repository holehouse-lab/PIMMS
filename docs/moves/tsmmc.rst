.. _move-tsmmc:

======================================
Temperature-switch Monte Carlo (TSMMC)
======================================

:Keywords: ``MOVE_CTSMMC`` (9), ``MOVE_MULTICHAIN_TSMMC`` (10), ``MOVE_SYSTEM_TSMMC`` (12)
:Status: stable
:Scope: one chain (9), a random subset of chains (10), the whole system (12)

TSMMC is one algorithm applied at three scopes, so the three move codes share this
page. All three are governed by the same ``TSMMC_*`` keywords (see
:doc:`/advanced/tsmmc`).

How it works
============

A TSMMC move takes part of the system on a **temperature excursion** to hop over
energy barriers that ordinary moves cannot cross. Starting from the simulation
temperature :math:`T`, the temperature is ramped up a schedule to a high "jump"
temperature ``TSMMC_JUMP_TEMP``, held there briefly, and ramped back down again; at
each rung of the schedule a burst of ordinary Monte Carlo moves is performed. At
the high temperature the system can climb out of a local energy minimum and
explore; on the way back down it re-settles, hopefully into a different basin. The
whole excursion is proposed as a single move and accepted or rejected as a unit,
and a rejected excursion restores the pre-move configuration exactly.

The three variants differ in what is heated and in which sub-moves run:

* ``MOVE_CTSMMC`` (code 9) - the single chain the main loop drew. Only crankshaft
  perturbations of that chain run at each rung.
* ``MOVE_MULTICHAIN_TSMMC`` (code 10) - a fresh random subset of the non-frozen
  chains, of size drawn uniformly between 1 and ``floor(0.25 * n_chains) + 1``,
  where ``n_chains`` is the number of non-frozen chains (clamped to that number),
  and sampled without replacement.
  Only crankshaft perturbations of those chains run at each rung. The chain the
  main loop drew plays no part in this selection.
* ``MOVE_SYSTEM_TSMMC`` (code 12) - the entire system (most powerful, most
  expensive). Here the main loop itself is driven along the schedule, so the
  sub-moves are drawn from the **full move set**. A nested TSMMC draw is not
  allowed: such a draw is redrawn from the non-TSMMC moves with their keyfile
  fractions renormalised, so the excursion samples the same relative move mix as
  the outer loop. A crankshaft fallback is used only if no non-TSMMC move is
  enabled at all, and it is announced once.

Frozen chains are excluded from selection by the chain and multi-chain variants
and are never moved by the system-wide variant's sub-moves.

The system-wide variant can do exactly what its sub-moves can do and nothing more,
and the start-up move-set checks treat it that way: a keyfile whose only other move
needs chains of three beads (say ``MOVE_CHAIN_PIVOT``) is refused for a system of
dimers even with ``MOVE_SYSTEM_TSMMC`` enabled, and ``MOVE_SYSTEM_TSMMC`` plus rigid
moves draws the same "no move can change a chain's shape" warning as the rigid moves
alone (see :ref:`moves-step-anatomy`). The chain and multi-chain variants reshape
chains through their crankshaft perturbations, so they count as conformational
moves.

Why detailed balance holds
==========================

A naive "heat, move, cool, then accept on the final energy" scheme would **not**
be balanced, because the moves made at the elevated temperatures bias the
proposal. TSMMC instead uses **tempered transitions** (Neal, 1996), which restore
exact balance by accumulating the thermodynamic *work* done as the temperature is
switched.

Write :math:`\beta = 1/T` for the simulation's inverse temperature. The excursion
departs from :math:`\beta_0 = \beta`, runs through a schedule of rungs
:math:`\beta_1, \beta_2, \dots, \beta_M`, and returns to :math:`\beta_0` at the
end; :math:`T` itself is **not** a rung. The rungs are an up-ramp of
``TSMMC_NUMBER_OF_POINTS`` steps starting just above :math:`T` and ending exactly
on the jump temperature, a hold of ten further rungs there (a fixed internal
constant, ``CONFIG.TOP_TEMP``), then the up-ramp reversed, so
:math:`M = 2\,\texttt{TSMMC\_NUMBER\_OF\_POINTS} + 10`. Let :math:`x_k` be the
configuration reached after the propagation at rung :math:`\beta_k`, with
:math:`x_0` the starting configuration. Each rung runs ordinary, balanced MC moves
at its own temperature. The excursion is accepted with

.. math::

   A = \min\!\left(1,\; e^{W}\right), \qquad
   W = \sum_{k=1}^{M} \bigl(\beta_{k-1} - \beta_k\bigr)\, E(x_{k-1})
       \;+\; \bigl(\beta_M - \beta_0\bigr)\, E(x_M),

where :math:`W` is the accumulated log-weight ("work") of the temperature
switching. The sum runs over *every* temperature change, including the initial
change off :math:`T` onto the first rung and the final change from the last rung
back onto :math:`T`; for an excursion in which nothing was accepted it telescopes
to exactly 0, so such an excursion is always accepted, as it must be.

Because the schedule is symmetric (it returns to :math:`\beta`), the protocol is
palindromic (the same number of sub-moves at *every* rung -
``TSMMC_STEP_MULTIPLIER`` for the system-wide excursion, that multiple of the number
of heated beads for the chain variants - and none at the target temperature) and
the sub-moves at each temperature are themselves balanced (reversible), this
acceptance makes the *whole* excursion satisfy detailed balance with respect to the
target Boltzmann distribution at :math:`T` - the
elevated-temperature exploration is corrected for exactly, so it changes the
dynamics but not the sampled distribution. The PIMMS implementation accumulates
:math:`W` as :math:`\sum(\beta_\text{before} - \beta_\text{after})\,E(x)` over the
schedule and is checked against a crankshaft-only reference by the
detailed-balance test suite.

Configuration
=============

The move fractions ``MOVE_CTSMMC`` / ``MOVE_MULTICHAIN_TSMMC`` /
``MOVE_SYSTEM_TSMMC`` (all default 0.0; all ``MOVE_*`` must sum to 1.0) select the
variants. The excursion itself is shaped by:

``TSMMC_JUMP_TEMP`` : float
    Peak temperature of the excursion (default 50.0). **Must exceed** the hottest
    temperature the run will visit - ``TEMPERATURE``, or both ``QUENCH_START`` and
    ``QUENCH_END`` in a quench run - unless ``TSMMC_FIXED_OFFSET`` is used. This is
    checked at start-up.

``TSMMC_NUMBER_OF_POINTS`` : int
    Number of temperature points on each ramp; a positive integer (default 20;
    more = smoother, more expensive). The full schedule is this many rungs up, a
    ten-rung hold at the jump temperature, then the mirrored ramp down - 50 rungs
    at the default.

``TSMMC_STEP_MULTIPLIER`` : int
    MC sub-steps performed at each rung of the schedule; a positive integer
    (default 50). For the chain and multi-chain variants this is multiplied by the
    number of beads being heated (the chain length, or the total beads of the
    selected chains); the system-wide variant performs exactly this many full Monte
    Carlo moves at each rung.

``TSMMC_INTERPOLATION_MODE`` : str
    How temperatures are spaced; currently only ``LINEAR`` (the default, read
    case-insensitively). Any other value is rejected when the keyfile is parsed.

``TSMMC_FIXED_OFFSET`` : float
    If set, the jump temperature is ``TEMPERATURE + TSMMC_FIXED_OFFSET`` rather
    than the absolute ``TSMMC_JUMP_TEMP`` (handy inside quench runs, where it
    tracks the moving temperature). It must be positive - a TSMMC excursion always
    heats. To disable it, omit the keyword; there is no ``False`` keyfile value.

TSMMC is most useful for strongly-interacting systems that get stuck; it is
expensive (each move is many sub-moves across the schedule), so it is typically
used at a small fraction alongside the crankshaft. Note that all three variants log
**one** attempt per excursion in ``MOVE_FREQS.dat`` / ``ACCEPTANCE.dat``; the
sub-moves made inside an excursion are counted separately, through the
alternative-Markov-chain counter that feeds ``TOTAL_MOVES.dat``.

Under ``PARALLELIZE`` the chain and multi-chain excursions always run their
crankshaft perturbations on the serial kernel. The system-wide excursion's
sub-moves are ordinary main-loop moves, so a crankshaft, slither or pull drawn
inside it uses the parallel kernel exactly as it would outside (this is why the
parallel slither and pull megamoves run their two passes in a random order: the
excursion's acceptance needs every rung's sub-moves to be reversible).

On STDOUT a system-wide excursion announces itself with
``Performing System TSMMC...`` and ends with ``System TSMMC: ACCEPTED [dE = ...]``
or ``System TSMMC: REJECTED [dE = ...]``, and every accepted multi-chain excursion
prints ``Multichain re-arrangement accepted [dE = ...]`` with the number of chains
it moved. ``REDUCED_PRINTING : True`` silences all of these.

For the full treatment - the temperature schedule, a step-by-step account of the
tempered-transitions work correction, a separate description of each of the three
variants (what is heated, which sub-moves run, and how the cost scales), and cost /
tuning guidance - see the dedicated :doc:`/advanced/tsmmc` page.
