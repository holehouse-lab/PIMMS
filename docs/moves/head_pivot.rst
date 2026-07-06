.. _move-head-pivot:

==========
Head pivot
==========

:Keyword: ``MOVE_HEAD_PIVOT``
:Move code: 5
:Status: core (but almost always set to 0)

How it works
============

A head pivot is a pivot (:doc:`chain_pivot`) restricted to a **single terminal
bead**: one of the two chain ends is chosen and that terminus is rotated about its
neighbour. Only the end bead moves, so the conformational change is tiny - much
smaller than a normal pivot, and largely redundant with what the crankshaft
already does at the chain ends.

Why detailed balance holds
==========================

As for the other rigid moves, the choice of terminus and rotation is symmetric
between the forward and reverse proposals,

.. math::

   g(x\to y) = g(y\to x),

so the plain Metropolis acceptance :math:`A = \min(1, e^{-\Delta E/T})` preserves
detailed balance (see :ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_HEAD_PIVOT`` : float
    Probability of selecting a head-pivot step (all ``MOVE_*`` must sum to 1.0).

In practice this move adds little over the crankshaft and we **recommend leaving
it at 0**. It is retained mainly for completeness.
