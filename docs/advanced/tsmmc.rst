.. _advanced-tsmmc:

======================================
Temperature-switch Monte Carlo (TSMMC)
======================================

Ordinary local moves relax a system efficiently but struggle to cross large energy
barriers. Once a chain - or the whole system - has fallen into a deep basin, almost
every proposed move is uphill and gets rejected, and the simulation can sit there
essentially stuck, sampling one basin while never discovering a lower one on the
other side of a barrier. **TSMMC** moves are built to get over exactly these
barriers.

The idea is a **temperature excursion**. Instead of proposing a single
perturbation, a TSMMC move takes part of the system (or all of it) on a round trip
in temperature: heat it along a schedule from the simulation temperature :math:`T`
up to a high "jump" temperature, let it rearrange while it is hot, then cool it back
to :math:`T`. While hot, moves that are effectively impossible at :math:`T` become
easy, so the system can climb out of its basin and wander; on the way back down it
settles into whatever basin it has found. The entire excursion - which is internally
thousands of ordinary sub-moves - is then accepted or rejected as a **single**
Monte Carlo move.

This page explains the acceptance correction that makes such a move valid, then
describes the three variants that differ in *what* is taken on the excursion.

Why a naive excursion would be wrong
====================================

The tempting scheme - "heat, move around, cool, then accept on the final energy with
the usual Metropolis rule" - is **not** correct. The moves made at the elevated
temperatures are drawn from hot distributions, not from the target distribution at
:math:`T`, so they bias the proposal in a way the ordinary Metropolis criterion does
not account for. Using it would distort the sampled ensemble: the simulation would
no longer converge on the Boltzmann distribution at :math:`T`.

TSMMC fixes this with a **tempered-transitions** acceptance (Neal, 1996; closely
related to nonequilibrium candidate Monte Carlo, Nilmeier *et al.*, 2011). The trick
is to account for the thermodynamic *work* done as the temperature is switched along
the schedule, and to fold that work into the accept/reject decision. Done correctly,
the elevated-temperature exploration is corrected for **exactly**: TSMMC changes how
the system explores, but not the distribution it converges on.

The temperature schedule
========================

PIMMS works in reduced units where the inverse temperature is simply
:math:`\beta = 1/T`. A TSMMC excursion departs from :math:`\beta_0 = 1/T`, runs
through a fixed schedule of inverse temperatures

.. math::

   \beta_1,\ \beta_2,\ \dots,\ \beta_M,

and returns to :math:`\beta_0` at the end. The schedule is built once when the move
begins (and rebuilt whenever the base temperature changes, for example by a
:doc:`quench <quench>`) out of three pieces:

* an **up-ramp** of ``TSMMC_NUMBER_OF_POINTS`` linearly-spaced rungs starting just
  above :math:`T` and ending exactly on the jump temperature ``TSMMC_JUMP_TEMP``,
* a **hold** of ten further rungs at the jump temperature (a fixed internal
  constant, ``CONFIG.TOP_TEMP``, not a keyword), and
* a **down-ramp** that is the up-ramp reversed, ending on the first rung above
  :math:`T`.

The schedule therefore has exactly

.. math::

   M \;=\; 2\times\texttt{TSMMC\_NUMBER\_OF\_POINTS} \;+\; 10

rungs, and the jump temperature itself is visited twelve times (the two ramp
endpoints plus the ten hold rungs). With the default
``TSMMC_NUMBER_OF_POINTS : 20`` that is 50 rungs. Note that :math:`T` is **not** a
rung: the schedule leaves the simulation temperature immediately and the return to
it is the final term of the work sum below, which is what makes the protocol
palindromic (exactly ``TSMMC_STEP_MULTIPLIER`` sub-moves at every rung, and none at
:math:`T`).

At **each** rung :math:`\beta_k` the system is propagated by a burst of ordinary
Monte Carlo moves that are themselves reversible at :math:`\beta_k` (they
individually satisfy detailed balance for the Boltzmann distribution at that
temperature). The number of sub-moves per rung is set by ``TSMMC_STEP_MULTIPLIER``
(scaled by the number of beads being heated for the chain-based variants - see
below).

The acceptance correction
==========================

Let :math:`x^{(k)}` be the configuration reached after the propagation at rung
:math:`\beta_k`, with :math:`x^{(0)}` the starting configuration, and let
:math:`U(x)` be its energy. The excursion accumulates a single scalar, the
**tempered-transitions work**

.. math::
   :label: tsmmc-work

   W \;=\; \sum_{k=1}^{M}\bigl(\beta_{k-1}-\beta_k\bigr)\,U\!\bigl(x^{(k-1)}\bigr)
       \;+\;\bigl(\beta_{M}-\beta_{0}\bigr)\,U\!\bigl(x^{(M)}\bigr).

Each term is the change in inverse temperature at one switch, multiplied by the
energy of the configuration **at the instant of that switch** (i.e. before the
system is propagated at the new temperature). The sum runs over every temperature
change in the excursion, including the initial departure from :math:`\beta_0` and
the final return to it. The whole excursion is then accepted with

.. math::
   :label: tsmmc-accept

   A \;=\; \min\!\bigl(1,\; e^{W}\bigr).

Two features are worth drawing out:

* **Only the temperature switches contribute.** The propagation *within* a rung uses
  a kernel that already preserves the Boltzmann distribution at :math:`\beta_k`, so
  it needs no correction and does not appear in :eq:`tsmmc-work`. All of the
  bookkeeping lives in the discrete temperature changes, each weighted by the
  instantaneous energy.

* **A no-op is always accepted.** If no sub-move is ever accepted, the energy never
  changes, :math:`U(x^{(k)}) = U_0` for all :math:`k`, and because the schedule
  returns to :math:`\beta_0` the coefficients telescope:

  .. math::

     W \;=\; U_0\!\left[\sum_{k=1}^{M}(\beta_{k-1}-\beta_k) + (\beta_M-\beta_0)\right]
       \;=\; U_0\,(\beta_0-\beta_0) \;=\; 0,

  giving :math:`A = 1`. A move that does nothing is never spuriously rejected - as it
  must not be.

More generally, because the schedule is a closed loop and every intra-rung kernel is
balanced, the correction :eq:`tsmmc-work` exactly cancels the bias introduced by
sampling at the elevated temperatures, so the composite excursion satisfies detailed
balance with respect to the target Boltzmann distribution at :math:`T`. Intuitively,
:math:`W` measures the reversible work of carrying the system around the temperature
loop: excursions that end up in configurations which are *favourable* at :math:`T`
relative to where they started are accepted readily, while those that end
"expensively" are penalised by exactly the amount needed to keep the ensemble
correct. PIMMS accumulates :math:`W` incrementally as
:math:`\sum(\beta_\text{before}-\beta_\text{after})\,U(x)` over the schedule, and the
implementation is checked against a crankshaft-only reference in the
detailed-balance test suite.

The three variants
==================

All three variants share the schedule and the acceptance rule above; they differ
only in **what** is taken on the excursion and, consequently, **which** moves run at
each rung and how expensive the excursion is.

.. list-table::
   :header-rows: 1
   :widths: 8 32 24 36

   * - Code
     - Keyword
     - Scope of excursion
     - Sub-moves used at each rung
   * - 9
     - ``MOVE_CTSMMC``
     - One randomly chosen (non-frozen) chain
     - Local (crankshaft) moves on that chain only
   * - 10
     - ``MOVE_MULTICHAIN_TSMMC``
     - A random subset of the non-frozen chains
     - Local (crankshaft) moves on the selected chains only
   * - 12
     - ``MOVE_SYSTEM_TSMMC``
     - The entire system
     - The full move set (any move except a nested TSMMC)

Chain TSMMC (``MOVE_CTSMMC``)
-----------------------------

A single chain is chosen at random and taken on the temperature excursion by itself.
Only **local crankshaft perturbations of that one chain** are performed at each rung
- the chain "wiggles" more and more freely as it heats, then re-tightens as it cools.
Restricting the excursion to local wiggling is deliberate: it lets the chain
thoroughly rearrange its own conformation without decorrelating so violently that the
final state is almost always rejected.

The amount of work per rung scales with the chain length: the number of sub-moves at
each temperature is

.. math::

   \texttt{chain length}\;\times\;\texttt{TSMMC\_STEP\_MULTIPLIER},

so a single chain-TSMMC move costs exactly :math:`M \times \texttt{chain length}
\times \texttt{TSMMC\_STEP\_MULTIPLIER}` proposed crankshaft sub-moves. Use it to
shake a single stubborn chain out of a bad conformation.

Multichain TSMMC (``MOVE_MULTICHAIN_TSMMC``)
--------------------------------------------

Exactly the chain excursion, but a **random subset of chains** is heated together,
with local moves performed across all of the selected chains. This is the tool for
chains that are stuck in a *cooperative* minimum - tangled or mutually frustrated in
a way no single-chain move can undo - without paying for a full-system excursion.

.. important::

   The subset is chosen as a plain **random selection of chains** - a size drawn
   uniformly between 1 and ``floor(0.25 * N) + 1``, where ``N`` is the number of
   non-frozen chains, then that many of them sampled without replacement - *not* as
   a connected cluster. Selecting
   a cluster would break detailed balance: after the excursion two clusters might
   have merged (or one split), and you could not propose the reverse move with the
   same probability, because the cluster structure has changed. Drawing the subset at
   random - **independently of the configuration** - makes its selection probability
   the same before and after the excursion, so the reverse move is proposable with
   equal probability and the acceptance :eq:`tsmmc-accept` stays valid.

The cost per rung scales with the **total number of beads in the selected subset**
times ``TSMMC_STEP_MULTIPLIER``. In testing this variant is a particularly effective
way to work *down* a rough energy landscape: it supplies the concerted, multi-chain
motion that ordinary single-move MC lacks, without having to guess in advance which
chains need to move together.

System TSMMC (``MOVE_SYSTEM_TSMMC``)
------------------------------------

The system-wide excursion is different in kind. The **entire lattice is backed up**,
and the main simulation is temporarily converted into a chain of auxiliary
simulations: the main-loop temperature is stepped along the schedule, and at each
rung a burst of **ordinary full-system Monte Carlo moves** is run - *any* move is
available except a nested TSMMC (TSMMC excursions are not allowed to recurse). A
nested TSMMC draw is redrawn from the non-TSMMC moves with their keyfile fractions
renormalised, so the excursion samples the same relative move mix as the outer loop;
only if the move set contains nothing *but* TSMMC moves does it fall back to
crankshaft megamoves, and that is announced once.

During the excursion the main step counter is not advanced and no analysis,
trajectory or energy output is produced; it is a self-contained super-move bolted
into the main chain. Its sub-moves go to separate auxiliary-chain counters, so
``MOVE_FREQS.dat`` records one attempt of move 12 per excursion (not the thousands of
sub-moves inside it); the sub-moves are folded into ``TOTAL_MOVES.dat`` when the
excursion completes. A system-wide excursion drawn on the last
step of the run is allowed to finish before the loop exits, so an excursion is never
serialised half-completed.

At the end the whole new system configuration is accepted or rejected with
:eq:`tsmmc-accept` (``System TSMMC: ACCEPTED [dE = ...]`` / ``REJECTED``, suppressed
by ``REDUCED_PRINTING``); on rejection the lattice is restored exactly from the
backup and the step reports the pre-move state. The
number of sub-moves per rung is ``TSMMC_STEP_MULTIPLIER`` (not scaled by size, since
each sub-move is already a full-system move). This is the most powerful variant - it
can rearrange the entire configuration cooperatively - and the most expensive, both
in time and in the memory needed to hold the backup.

Configuring the excursion
=========================

.. code-block:: text

   MOVE_CRANKSHAFT          : 0.9
   MOVE_SYSTEM_TSMMC        : 0.1
   TSMMC_JUMP_TEMP          : 120     # peak temperature (must exceed TEMPERATURE)
   TSMMC_NUMBER_OF_POINTS   : 20      # rungs on each ramp (smoother = more)
   TSMMC_STEP_MULTIPLIER    : 50      # sub-steps per temperature rung
   TSMMC_INTERPOLATION_MODE : LINEAR

With these values the schedule is ``2 x 20 + 10 = 50`` rungs of 50 sub-moves each,
so every draw of ``MOVE_SYSTEM_TSMMC`` costs 2500 full-system Monte Carlo moves.

.. list-table::
   :header-rows: 1
   :widths: 28 12 12 48

   * - Keyword
     - Type
     - Default
     - Meaning
   * - ``TSMMC_JUMP_TEMP``
     - float
     - ``50.0``
     - Peak temperature of the excursion. **Must be greater than the hottest
       temperature the run will visit** - ``TEMPERATURE``, or both
       ``QUENCH_START`` and ``QUENCH_END`` in a quench run - and this is checked
       at start-up. Ignored entirely when ``TSMMC_FIXED_OFFSET`` is set.
   * - ``TSMMC_NUMBER_OF_POINTS``
     - int
     - ``20``
     - Number of temperature rungs on each ramp between :math:`T` and the jump
       temperature. More points = a smoother, more gradual (and more expensive)
       excursion. Must be a positive integer.
   * - ``TSMMC_STEP_MULTIPLIER``
     - int
     - ``50``
     - Monte Carlo sub-steps performed at **each** rung of the schedule (further
       multiplied by the number of beads being heated for the chain/multichain
       variants). More = more thorough equilibration at each temperature, but
       slower. Must be a positive integer.
   * - ``TSMMC_INTERPOLATION_MODE``
     - str
     - ``LINEAR``
     - How the temperature is spaced across the schedule. Currently only
       ``LINEAR`` (equal increments) is supported; anything else is rejected at
       parse time.
   * - ``TSMMC_FIXED_OFFSET``
     - float
     - unset
     - If set, the jump temperature is defined **relative** to the
       current temperature as :math:`T + \texttt{TSMMC\_FIXED\_OFFSET}` instead of
       the absolute ``TSMMC_JUMP_TEMP``. Handy inside a quench (see below). Must
       be positive (a TSMMC excursion always heats); a non-positive offset is
       rejected at parse time. To turn it off, **omit the keyword** - there is no
       ``TSMMC_FIXED_OFFSET : False`` keyfile syntax.

The move fractions ``MOVE_CTSMMC`` / ``MOVE_MULTICHAIN_TSMMC`` /
``MOVE_SYSTEM_TSMMC`` select the variants and, like all ``MOVE_*`` keywords, must
sum to ``1.0`` across the whole move set. All three default to ``0.0``, and the
``TSMMC_*`` keywords only take effect when at least one of them is non-zero. Their
validation is not all conditional, though: ``TSMMC_INTERPOLATION_MODE`` is checked
the moment it is read, and ``TSMMC_STEP_MULTIPLIER`` and ``TSMMC_NUMBER_OF_POINTS``
must be positive whether or not a TSMMC move is in the move set. It is only the
jump-temperature and ``TSMMC_FIXED_OFFSET`` checks that are skipped when no TSMMC
move is enabled.

Cost and tuning
===============

A single TSMMC move is worth exactly

.. math::

   \underbrace{\bigl(2\times\texttt{TSMMC\_NUMBER\_OF\_POINTS} + 10\bigr)}_{\text{rungs, } M}
   \;\times\;\text{(sub-moves per rung)}

ordinary sub-moves, so even a small ``MOVE_*TSMMC`` probability represents a large
amount of work. Guidance:

* **Jump temperature.** Hot enough to melt the barrier you are trying to cross, but
  not so hot that the system fully evaporates and forgets everything useful.
  ``TSMMC_JUMP_TEMP`` must always exceed the base ``TEMPERATURE`` (and, in a quench
  run, both ends of the ramp).
* **Points and step multiplier.** More points and more sub-steps make each excursion
  gentler and more thorough, and therefore more likely to be accepted, at a directly
  proportional cost. If excursions are almost always rejected, make the schedule
  *smoother* (more points) before making it *hotter*.
* **Fraction.** Because each excursion is expensive, TSMMC is normally mixed in at a
  small fraction, with cheap local moves (crankshaft, slither) doing the routine
  sampling in between.
* **Pick the smallest scope that works.** Reach for ``MOVE_CTSMMC`` for a single
  stuck chain, ``MOVE_MULTICHAIN_TSMMC`` for a cooperatively-trapped handful, and
  ``MOVE_SYSTEM_TSMMC`` only when the whole configuration must move at once.

Inside a quench
===============

If you combine TSMMC with a :doc:`quench <quench>`, the excursion coordinator is
rebuilt at the new base temperature every time the quench updates, so the jumps
always heat relative to wherever the ramp currently sits. ``TSMMC_FIXED_OFFSET`` is
the natural choice there: it pins the jump temperature a fixed amount above the
current temperature, so you never risk an absolute ``TSMMC_JUMP_TEMP`` accidentally
falling below the (moving) base temperature during the ramp. A keyfile whose ramp
would reach or exceed the absolute ``TSMMC_JUMP_TEMP`` is rejected at parse time
(the jump temperature must exceed both ``QUENCH_START`` and ``QUENCH_END``), so an
excursion can never be silently inverted mid-run - use ``TSMMC_FIXED_OFFSET``, or
keep the hot end of the ramp below ``TSMMC_JUMP_TEMP``. Should a schedule ever be
asked for with a jump temperature at or below the target, the coordinator refuses to
build it and the run stops with an explicit message rather than quietly running a
flat or cooling "excursion".

References
==========

* R. M. Neal, "Sampling from multimodal distributions using tempered transitions",
  *Statistics and Computing* **6**, 353-366 (1996).
* J. P. Nilmeier, G. E. Crooks, D. D. L. Minh, J. D. Chodera, "Nonequilibrium
  candidate Monte Carlo is an efficient tool for equilibrium simulation",
  *PNAS* **108**, E1009-E1018 (2011).

See also
========

The per-move :doc:`/moves/tsmmc` page places TSMMC in the context of the full move
set and the shared detailed-balance primer.
