.. _move-head-pivot:

==========
Head pivot
==========

:Keyword: ``MOVE_HEAD_PIVOT``
:Move code: 5
:Status: core (but almost always set to 0)
:Scope: the single chain drawn by the main loop

How it works
============

A head pivot moves a **single terminal bead**. One of the two chain ends is chosen
with probability :math:`1/2`, that bead is lifted off the lattice, and a new site
is drawn uniformly from the full :math:`3^d` block of sites around its bonded
neighbour (:math:`d` = 2 or 3) - so the terminus stays bonded wherever it lands.
Only the end bead moves, so the conformational change is tiny.

Three things reject the move outright: a draw that lands on an occupied site
(including the neighbour's own site, which is always in the candidate block), a
draw that lands back on the terminus' original site (treated as a null move rather
than a trivially accepted no-op), and, under ``HARDWALL``, a draw whose new bond
would cross the wall. A **single-bead chain has no terminus to pivot** and rejects
the move immediately.

This is the same proposal the :doc:`crankshaft` already makes for a terminal bead,
which is why the move adds essentially nothing to a move set that includes the
crankshaft. It is retained because the :doc:`chain_pivot` cannot reach the ends of
very short chains.

Why detailed balance holds
==========================

Both the choice of terminus and the candidate block are fixed by the *neighbour's*
position, which the move does not change, so the reverse move draws the original
site from the same block of :math:`3^d` sites with the same probability:

.. math::

   g(x\to y) = g(y\to x) = \tfrac12 \cdot 3^{-d}.

The plain Metropolis acceptance :math:`A = \min(1, e^{-\Delta E/T})` therefore
preserves detailed balance (see :ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_HEAD_PIVOT`` : float
    Probability of selecting a head-pivot step (all ``MOVE_*`` must sum to 1.0).
    Default 0.0.

In practice this move adds little over the crankshaft and we **recommend leaving
it at 0**. It is retained mainly for completeness.
