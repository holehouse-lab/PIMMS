.. _move-segment-regrowth:

================
Segment regrowth
================

:Keyword: ``MOVE_SEGMENT_REGROWTH``
:Move code: 15
:Status: experimental (needs ``EXPERIMENTAL_FEATURES : True``)
:Scope: every partly frozen chain (``FREEZE_FILE`` ``B`` lines); nothing else

Segment regrowth moves the free beads of chains you've partly frozen, a short
window at a time. It exists because the :doc:`crankshaft`, which moves one bead at
a time, can't reach every arrangement of two or more free beads held at both ends
by frozen beads. The simplest case is two free beads that close a 2x2 square with
their frozen neighbours: each one's only candidate sites are the four corners of
the square, all occupied, so neither ever moves. A window that covers both beads
can move them together.

It never touches a chain with no frozen bead, or a chain frozen whole. In a
system with no partly frozen chain it does nothing, and PIMMS says so at start-up
(and refuses the run if it's the only move).

How it works
============

At start-up each partly frozen chain gets a fixed list of **windows**: every run of
1 to ``SEGMENT_REGROWTH_MAX_BEADS`` consecutive free beads. With beads 0 and 4 of
an 8-bead chain frozen and windows of up to 2 beads, for example, the windows are
1, 1-2, 2, 2-3, 3, 5, 5-6, 6, 6-7 and 7.

When the move is drawn it's a **sweep**: every partly frozen chain, taken in a
random order, gets one window move, so a system with 40 partly frozen chains makes
40 window moves per draw. A window move

1. picks one of the chain's windows uniformly at random;
2. lifts the window's beads off the lattice;
3. draws new sites for them uniformly from **every lattice path of that many
   beads** that joins the window's two fixed neighbours - or, for a window at a
   chain end, that leaves its one fixed neighbour - with every bead bonded to the
   next (a bond is a step to any of the 8 neighbouring sites in 2D, 26 in 3D);
4. rejects the proposal if a new site is occupied, or (``HARDWALL``) if the path
   leaves the box; otherwise
5. accepts or rejects the new arrangement on the Metropolis criterion, putting the
   window back if it's rejected.

The paths come from tables built the first time they're needed. They're small for
short windows and grow quickly with length:

.. list-table::
   :header-rows: 1
   :widths: 30 20 25 25

   * - Window
     - Lattice
     - Paths from one end
     - Most paths between two fixed ends
   * - 2 beads
     - 3D
     - 650
     - 192
   * - 3 beads
     - 3D
     - 15,986
     - 3,008
   * - 4 beads
     - 3D
     - 388,994
     - 49,288
   * - 3 beads
     - 2D
     - 368
     - 58

A one-bead window proposes the same new sites as the crankshaft does for that
bead. (The crankshaft also draws, and rejects, the neighbours' own sites, so the
two count their rejected draws differently.) A draw that puts the window back on
the path it already has is rejected as a null move, so ``ACCEPTANCE.dat`` counts
only window moves that changed something.

Why detailed balance holds
==========================

The windows are fixed for the whole run, and the set of paths a window can be given
depends only on where its fixed neighbours are - which the move never changes - not
on where the window's own beads are. So going from arrangement :math:`x` to
:math:`y` is proposed with exactly the same probability as going back:

.. math::

   g(x\to y) = g(y\to x) = \frac{1}{n_\text{windows}} \cdot \frac{1}{n_\text{paths}},

and the plain Metropolis acceptance :math:`A = \min(1, e^{-\Delta E/T})` gives
detailed balance (see :ref:`the primer <moves-db-primer>`). An occupied or
out-of-box path is just a rejected proposal. Taking the chains in a fresh random
order every sweep makes the sweep as a whole reversible too, which a
``MOVE_SYSTEM_TSMMC`` excursion relies on when it runs this move as a sub-move.

The test suite checks the result directly: it enumerates every arrangement of the
free beads in small 2D systems (the trapped square included), computes their exact
Boltzmann weights, and checks the sampled frequencies against them, with a control
that fails when the temperature is 30% off.

What it guarantees, and what it doesn't
=======================================

* A free stretch held at both ends that's **no longer than**
  ``SEGMENT_REGROWTH_MAX_BEADS`` can go from any arrangement to any other in one
  move, so the start-up warning about trapped stretches only names longer ones.
* A **longer** stretch, or a free tail longer than the window, isn't guaranteed.
  It's far better connected than with the crankshaft alone, but no move that only
  touches a few beads at a time can be shown to reach every arrangement of a long
  chain (Madras and Sokal, *J. Stat. Phys.* 47, 573-595, 1987, proved this for
  lattice self-avoiding walks moved by any fixed set of local, length-conserving
  moves). In open space that rarely matters; in a highly frozen structure it can
  (see :ref:`below <segment-regrowth-frozen-limit>`).
* In a **periodic box no wider than** 2 x (stretch + 1) sites, where the stretch
  is the free beads between two frozen ones, the stretch can join its two frozen
  ends going either way round the box. No local move that keeps the frozen beads
  fixed - this one or the crankshaft - can switch a stretch from one way to the
  other, so it's sampled exactly but only among the arrangements that go the way
  it started. That matters only for a loop long enough to reach round the box.

.. _segment-regrowth-frozen-limit:

Highly frozen structures
========================

.. warning::

   If most of the space around a free loop is frozen, a loop longer than
   ``SEGMENT_REGROWTH_MAX_BEADS`` - or a free tail longer than it - can get stuck
   in one group of its arrangements for the whole run. It's still sampled
   correctly within that group, so nothing looks wrong, but what you measure for
   it depends on where it started.

A window move changes a few neighbouring beads at a time while the rest of the
loop stays put. So to get from one shape to another, the loop has to pass through
a series of shapes in between, each one window move from the last. In open space
there are plenty of those. Hemmed in by frozen beads there may be none, because
every in-between shape puts a bead on a frozen site.

Here's the smallest case an exhaustive search turned up: a 4-bead loop between
frozen beads ``A`` and ``B`` in 2D, with ``#`` marking frozen sites:

.. code-block:: none

    . # A # B . #          . # A # B 4 #
    # 1 2 # 4 # #          # . 1 # 3 # #
    # . # 3 # . .          # . # 2 # . .

These are the loop's only two arrangements. Going from one to the other means all
four beads sliding one step along the channel at once, and a window of three
always leaves one of them where it is, so the loop never leaves whichever
arrangement it starts in. A loop wrapped round a frozen part of the structure
gets stuck the same way: reaching the other side means passing through shapes
that overlap frozen beads.

How often this happens depends on how tightly the loop is packed in. We checked
every arrangement of loops surrounded by frozen sites scattered at random (20 to
50% of the nearby sites). With windows of up to 3 beads, 1 to 11% of 4-bead
loops and 7 to 38% of 7-bead loops were split in 2D; with the crankshaft alone it
was 60 to 90%. In 3D it's much rarer: no 4-bead loop was split with windows of 3
(longer 3D loops were too big to check). Random frozen sites aren't a folded
protein, so take these numbers as a sense of scale, not a forecast. With no
frozen beads nearby, every loop and tail we checked could reach every
arrangement.

What you can do:

* Keep enclosed loops no longer than ``SEGMENT_REGROWTH_MAX_BEADS``, raising it to
  4 if you need to. PIMMS warns at start-up about every enclosed stretch longer
  than the window, so read that warning.
* If a loop has to be longer, freezing one bead partway along it splits it into
  two shorter loops. That holds the bead fixed, though, so only do it where that's
  physically reasonable.
* Free tails aren't covered by the start-up warning. A tail longer than the
  window that sits in a tight frozen pocket can get stuck too.
* If you can, start runs from different arrangements of the loop and check that
  what you measure for it agrees. If it doesn't, the loop is stuck.

Cost
====

Each window move is a serial Python move, also under ``PARALLELIZE``. Its energy
change is the same local evaluation a chain pivot uses, so it's cheap - 0.05 to
0.1 ms per window in small 3D test systems, up to about 0.4 ms in a crowded one
with long-range interactions - but it isn't a compiled megamove. It suits tens to
hundreds of partly frozen chains; with a thousand, one sweep takes roughly 0.05
to 0.4 s, so keep ``MOVE_SEGMENT_REGROWTH`` modest. Longer windows are accepted less
often in crowded surroundings, since more of their paths land on an occupied
site.

``MOVE_FREQS.dat`` and ``ACCEPTANCE.dat`` count one attempt per partly frozen chain
per sweep, in column 15.

Configuration
=============

``MOVE_SEGMENT_REGROWTH`` : float
    Probability of drawing a segment-regrowth sweep (all ``MOVE_*`` must sum to
    1.0). Default 0.0. Experimental.

``SEGMENT_REGROWTH_MAX_BEADS`` : int
    The longest window, 1 to 4. Default 3. Experimental.

A run that frees a few loops and tails of a folded structure might use:

.. code-block:: none

   EXPERIMENTAL_FEATURES : True
   FREEZE_FILE           : freeze.txt   # B lines freeze everything but the loops and tails

   MOVE_CRANKSHAFT       : 0.8          # single-bead moves of every free bead
   MOVE_SEGMENT_REGROWTH : 0.2          # one window move per partly frozen chain
   SEGMENT_REGROWTH_MAX_BEADS : 3       # loops of up to 3 free beads are fully sampled
