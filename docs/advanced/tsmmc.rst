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
:math:`\beta = 1/T`. A TSMMC excursion follows a fixed schedule of inverse
temperatures

.. math::

   \beta_0 \;=\; \tfrac{1}{T},\quad \beta_1,\ \beta_2,\ \dots,\ \beta_M,

built once at the start of the move:

* an **up-ramp** of ``TSMMC_NUMBER_OF_POINTS`` linearly-spaced steps from just above
  :math:`T` to the jump temperature ``TSMMC_JUMP_TEMP``,
* a short **hold** at the jump temperature (a fixed number of points at the top of
  the ramp), and
* a **down-ramp** that mirrors the up-ramp back toward :math:`T`.

So the schedule is a closed loop in temperature - it leaves :math:`T`, reaches the
jump temperature, and returns - of total length :math:`M \approx
2\times\texttt{TSMMC\_NUMBER\_OF\_POINTS}` plus the hold. At **each** rung
:math:`\beta_k` the system is propagated by a burst of ordinary Monte Carlo moves
that are themselves reversible at :math:`\beta_k` (they individually satisfy
detailed balance for the Boltzmann distribution at that temperature). The number of
sub-moves per rung is set by ``TSMMC_STEP_MULTIPLIER`` (scaled by system size for the
chain-based variants - see below).

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
     - One randomly chosen chain
     - Local (crankshaft) moves on that chain only
   * - 10
     - ``MOVE_MULTICHAIN_TSMMC``
     - A random subset of chains
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

so a single chain-TSMMC move costs roughly :math:`M \times \texttt{chain length}
\times \texttt{TSMMC\_STEP\_MULTIPLIER}` crankshaft sub-moves. Use it to shake a
single stubborn chain out of a bad conformation.

Multichain TSMMC (``MOVE_MULTICHAIN_TSMMC``)
--------------------------------------------

Exactly the chain excursion, but a **random subset of chains** is heated together,
with local moves performed across all of the selected chains. This is the tool for
chains that are stuck in a *cooperative* minimum - tangled or mutually frustrated in
a way no single-chain move can undo - without paying for a full-system excursion.

.. important::

   The subset is chosen as a plain **random selection of chains** (a random number
   of them, up to a fraction of the system), *not* as a connected cluster. Selecting
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
available except a nested TSMMC (TSMMC excursions are not allowed to recurse). During
the excursion the main step counter is not advanced and no analysis or trajectory
output is produced; it is a self-contained super-move bolted into the main chain.

At the end the whole new system configuration is accepted or rejected with
:eq:`tsmmc-accept`; on rejection the lattice is restored exactly from the backup. The
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
   TSMMC_NUMBER_OF_POINTS   : 20      # points on each ramp (smoother = more)
   TSMMC_STEP_MULTIPLIER    : 50      # sub-steps per temperature rung
   TSMMC_INTERPOLATION_MODE : LINEAR

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Keyword
     - Type
     - Meaning
   * - ``TSMMC_JUMP_TEMP``
     - float
     - Peak temperature of the excursion. **Must be greater than**
       ``TEMPERATURE``.
   * - ``TSMMC_NUMBER_OF_POINTS``
     - int
     - Number of temperature points on each ramp between :math:`T` and the jump
       temperature. More points = a smoother, more gradual (and more expensive)
       excursion.
   * - ``TSMMC_STEP_MULTIPLIER``
     - int
     - Monte Carlo sub-steps performed at **each** rung of the schedule (further
       multiplied by the number of beads for the chain/multichain variants). More
       = more thorough equilibration at each temperature, but slower.
   * - ``TSMMC_INTERPOLATION_MODE``
     - str
     - How the temperature is spaced across the schedule. Currently only
       ``LINEAR`` (equal increments) is supported.
   * - ``TSMMC_FIXED_OFFSET``
     - float / False
     - If set to a number, the jump temperature is defined **relative** to the
       current temperature as :math:`T + \texttt{TSMMC\_FIXED\_OFFSET}` instead of
       the absolute ``TSMMC_JUMP_TEMP``. Handy inside a quench (see below).

The move fractions ``MOVE_CTSMMC`` / ``MOVE_MULTICHAIN_TSMMC`` /
``MOVE_SYSTEM_TSMMC`` select the variants and, like all ``MOVE_*`` keywords, must
sum to ``1.0`` across the whole move set.

Cost and tuning
===============

A single TSMMC move is worth roughly

.. math::

   \underbrace{\bigl(2\times\texttt{TSMMC\_NUMBER\_OF\_POINTS} + \text{hold}\bigr)}_{\text{rungs, } M}
   \;\times\;\text{(sub-moves per rung)}

ordinary sub-moves, so even a small ``MOVE_*TSMMC`` probability represents a large
amount of work. Guidance:

* **Jump temperature.** Hot enough to melt the barrier you are trying to cross, but
  not so hot that the system fully evaporates and forgets everything useful.
  ``TSMMC_JUMP_TEMP`` must always exceed the base ``TEMPERATURE``.
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
falling below the (moving) base temperature during the ramp.

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
