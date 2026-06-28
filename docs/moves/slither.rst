.. _move-slither:

=====================
Slither (reptation)
=====================

:Keyword: ``MOVE_SLITHER``
:Move code: 6
:Status: core

How it works
============

A slither move is **reptation** - the chain advances forwards or backwards through
the lattice "like a snake". A direction is chosen; the bead at the trailing end is
removed and a new bead is grown at the leading end into a randomly chosen empty,
connectivity-preserving neighbour site. Every interior bead effectively shifts one
place along the contour, so the chain crawls a step without any large rigid
motion. This relaxes chain conformations and lets chains thread through crowded
surroundings very efficiently.

A slither *step* is a **megamove**: when selected, every non-frozen chain is
slithered ``SLITHER_SUBSTEPS`` times in random order (each substep with its own
accept/reject), in an optimised Cython kernel. A single-bead chain has no contour
to reptate along, so its slither degenerates to a local translation.

Why detailed balance holds
==========================

Reptation is *not* a symmetric proposal: the number of empty directions in which
the leading end can grow, :math:`n_\text{fwd}`, generally differs from the number
of ways the reverse crawl could regrow the other end, :math:`n_\text{rev}`. The
move therefore uses the full Metropolis-Hastings acceptance with this
proposal-multiplicity ratio,

.. math::

   A(x\to y) = \min\!\left(1,\; \frac{n_\text{rev}}{n_\text{fwd}}\;
   e^{-\Delta E / T}\right),

which satisfies detailed balance (see :ref:`the primer <moves-db-primer>`). The
energy change :math:`\Delta E` is cheap to evaluate for a **homopolymer** - only
the removed and added end beads change their surroundings, an :math:`O(1)`
calculation - but for a **heteropolymer** the bead identities shift relative to
their positions, so every bead's interactions change and :math:`\Delta E` is an
:math:`O(N)` calculation. The kernel implements both paths and the multiplicity
bookkeeping; the result is checked by the detailed-balance test suite (a
crankshaft-only run and a crankshaft+slither run reach the same equilibrium).

Configuration
=============

``MOVE_SLITHER`` : float
    Probability of selecting a slither megamove (all ``MOVE_*`` must sum to 1.0).

``SLITHER_SUBSTEPS`` : int
    Number of slither moves applied to each chain, in random order, per megamove
    (default 10).

Slither works in both 2D and 3D and is one of the most effective moves for
relaxing chain conformations; a healthy fraction alongside the crankshaft usually
improves mixing markedly. Along with the crankshaft, it is one of the two moves
with a multi-threaded kernel: see :ref:`PARALLELIZE <advanced-parallel>` (it
parallelizes the chains that fit within a block interior).
