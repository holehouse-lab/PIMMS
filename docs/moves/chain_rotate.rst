.. _move-chain-rotate:

============
Chain rotate
============

:Keyword: ``MOVE_CHAIN_ROTATE``
:Move code: 3
:Status: core

How it works
============

The whole chain is rotated as a **rigid body** about its centre of mass. On the
lattice only the cardinal rotations are used - 90°, 180° or 270° (in 3D, about a
randomly chosen axis) - because arbitrary angles do not map lattice sites onto
lattice sites. The chain is translated so its centre of mass is at the origin,
rotated, and translated back. As with translation, a rotated bead landing on an
occupied site (or, under ``HARDWALL``, straddling the boundary) rejects the move.
The internal conformation is preserved; only the chain's orientation changes.

Why detailed balance holds
==========================

The proposed rotation is chosen uniformly from the cardinal rotations, and each
rotation's inverse (e.g. 90° ↔ 270°, 180° ↔ 180°) is equally likely to be
proposed for the reverse move, so the proposal is symmetric,

.. math::

   g(x\to y) = g(y\to x).

The move is accepted with the plain Metropolis criterion
:math:`A = \min(1, e^{-\Delta E/T})`, satisfying detailed balance (see
:ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CHAIN_ROTATE`` : float
    Probability of selecting a chain-rotation step (all ``MOVE_*`` must sum to
    1.0).

No other tuning keywords. Like translation, rotation is most useful in dilute
systems; in dense phases most rotations clash.
