.. _move-chain-translate:

===============
Chain translate
===============

:Keyword: ``MOVE_CHAIN_TRANSLATE``
:Move code: 2
:Status: core

How it works
============

The whole chain is moved as a **rigid body** by a single random displacement
vector: a per-dimension offset is drawn uniformly and every bead of the chain is
shifted by it (with periodic wrapping). If any translated bead would land on an
occupied site the move is rejected as a hard-sphere clash; under ``HARDWALL`` a
chain that would straddle the boundary is also rejected. The chain's internal
conformation is unchanged - only its position in the box moves.

This is one of the more expensive single-chain moves, because the energy of every
bead has to be re-evaluated against its new surroundings, but it is the primary
way an intact chain explores different regions of the box.

Why detailed balance holds
==========================

The displacement is drawn uniformly, and on the lattice a shift by :math:`+v` and
the reverse shift by :math:`-v` are equally likely, so the proposal is symmetric,

.. math::

   g(x\to y) = g(y\to x).

The move is therefore accepted with the plain Metropolis criterion
:math:`A = \min(1, e^{-\Delta E/T})`, which satisfies detailed balance (see
:ref:`the primer <moves-db-primer>`). Hard-sphere clashes correspond to
:math:`\pi = 0` states and are rejected with certainty, consistent with balance.

Configuration
=============

``MOVE_CHAIN_TRANSLATE`` : float
    Probability of selecting a chain-translation step (all ``MOVE_*`` must sum to
    1.0).

There are no other tuning keywords. In dense systems most translations clash and
are rejected; for relocating chains through crowded/condensed phases prefer the
collective moves (:doc:`vmmc`, :doc:`pull`) instead.
