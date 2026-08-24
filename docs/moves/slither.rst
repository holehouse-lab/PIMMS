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

The slither proposal is **symmetric**. The direction (head-first or tail-first)
is drawn with probability :math:`1/2`, and the new end site is drawn uniformly
from the *fixed-size* Chebyshev-1 offset box around the growing end
(:math:`3^d` sites, the current site included); a draw landing on an occupied
site is rejected as a hard-sphere clash rather than redrawn. Because the
candidate set has the same fixed size in both directions,
:math:`g(x\to y) = g(y\to x) = \tfrac12 \cdot 3^{-d}`, and the plain
Metropolis criterion

.. math::

   A(x\to y) = \min\!\left(1,\; e^{-\Delta E / T}\right)

satisfies detailed balance (see :ref:`the primer <moves-db-primer>`). The
energy change :math:`\Delta E` is cheap to evaluate for a **homopolymer** - only
the removed and added end beads change their surroundings, an :math:`O(1)`
calculation - but for a **heteropolymer** the bead identities shift relative to
their positions, so every bead's interactions change and :math:`\Delta E` is an
:math:`O(N)` calculation. The kernel implements both paths; the result is
checked by the detailed-balance test suite (a crankshaft-only run and a
crankshaft+slither run reach the same equilibrium).

Configuration
=============

``MOVE_SLITHER`` : float
    Probability of selecting a slither megamove (all ``MOVE_*`` must sum to 1.0).

``SLITHER_SUBSTEPS`` : int
    Number of slither moves applied to each chain, in random order, per megamove
    (default 10).

Slither works in both 2D and 3D and is one of the most effective moves for
relaxing chain conformations; a healthy fraction alongside the crankshaft usually
improves mixing markedly. Along with the crankshaft and pull, it is one of the
three moves with a multi-threaded kernel: see :ref:`PARALLELIZE <advanced-parallel>`
(it parallelizes the chains that fit within a block interior; if any chain is too
long ever to fit a block interior - or exceeds the kernel's 512-bead per-chain
buffer for heteropolymers - the whole megamove automatically falls back to the
serial kernel, so ``PARALLELIZE`` never changes the sampling, only the speed).
