.. _move-chain-pivot:

===========
Chain pivot
===========

:Keyword: ``MOVE_CHAIN_PIVOT``
:Move code: 4
:Status: core

How it works
============

A pivot move chooses a random interior bead and rigidly **rotates the whole
portion of the chain on one side of it** (a cardinal lattice rotation) about that
pivot point, leaving the other side fixed. Because it moves a large, contiguous
section of the chain in one shot, a pivot makes much larger conformational changes
than a crankshaft - it is an efficient way to decorrelate the global shape of a
chain, especially for swollen/dilute chains. A clash of any moved bead (or a
hardwall straddle) rejects the move.

Why detailed balance holds
==========================

The pivot point is chosen uniformly along the chain and the rotation uniformly
from the cardinal rotations; both choices, and their inverses, are equally likely
for the reverse move, so the proposal is symmetric,

.. math::

   g(x\to y) = g(y\to x),

and the move uses the plain Metropolis acceptance
:math:`A = \min(1, e^{-\Delta E/T})` (see :ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CHAIN_PIVOT`` : float
    Probability of selecting a pivot step (all ``MOVE_*`` must sum to 1.0).

Pivots have a low acceptance rate in dense systems (the rotated tail usually
clashes) but are very effective for single chains and dilute solutions, where a
modest fraction substantially speeds up sampling of the overall chain dimensions.
