.. _move-chain-rotate:

============
Chain rotate
============

:Keyword: ``MOVE_CHAIN_ROTATE``
:Move code: 3
:Status: core

How it works
============

The whole chain is rotated as a **rigid body** about one of its own beads - the
bead nearest the chain's (single-image) centroid. On the lattice only the
cardinal rotations are used - 90°, 180° or 270° (in 3D, about a randomly chosen
axis) - because arbitrary angles do not map lattice sites onto lattice sites.
The chain's displacement vectors relative to the pivot bead are rotated and
re-anchored on the pivot bead's original lattice position, which makes every
rotation exactly invertible for any box shape, including chains straddling a
periodic boundary. As with translation, a rotated bead landing on an occupied
site (or, under ``HARDWALL``, straddling the boundary) rejects the move. The
internal conformation is preserved; only the chain's orientation changes.

Why detailed balance holds
==========================

The proposed rotation is chosen uniformly from the cardinal rotations, and each
rotation's inverse (e.g. 90° ↔ 270°, 180° ↔ 180°) is equally likely to be
proposed for the reverse move. Anchoring the rotation on a *bead* matters here:
a bead is a lattice point that maps exactly to itself under the rotation, and
the bead nearest the centroid is a rotation-invariant choice for a rigid body,
so the reverse move picks the same physical pivot and the inverse rotation
returns the original configuration bit for bit. (Rotating about a rounded
centre of mass - the pre-1.0.8 behaviour - is *not* invertible: the rounded COM
of the rotated chain is generally a different lattice point, which silently
violated detailed balance.) The proposal is therefore symmetric,

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
