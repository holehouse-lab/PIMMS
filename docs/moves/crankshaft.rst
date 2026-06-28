.. _move-crankshaft:

==========
Crankshaft
==========

:Keyword: ``MOVE_CRANKSHAFT``
:Move code: 1
:Status: core (recommended workhorse)

How it works
============

The crankshaft is PIMMS' fundamental local move and should make up most of a
typical move budget. It perturbs **one bead at a time**: for a chosen bead it
finds the lattice sites that keep the chain connected - for an interior bead, the
empty sites that are Chebyshev-1 adjacent to *both* of its chain neighbours; for a
terminal bead, the empty sites adjacent to its single neighbour - and proposes
moving the bead to one of them chosen uniformly at random. The chain's identity
and bonding are preserved; only the kink at that bead changes.

A single crankshaft *step* is a **megamove**: it sweeps the system performing
``CRANKSHAFT_SUBSTEPS`` such single-bead perturbations per bead, each with its own
accept/reject, all inside an optimised Cython kernel. One step can therefore
encompass millions of elementary Monte Carlo moves.

Why detailed balance holds
==========================

The set of valid destination sites for a bead is determined solely by its
neighbours' (fixed) positions, not by the bead's current position. Hence the
forward proposal "bead at :math:`A \to B`" and the reverse "bead at
:math:`B \to A`" are drawn from the *same* set of size :math:`N`, so

.. math::

   g(x\to y) = g(y\to x) = \tfrac{1}{N}.

The proposal is symmetric, and the move is accepted with the plain Metropolis
criterion

.. math::

   A(x\to y) = \min\!\left(1,\; e^{-\Delta E / T}\right),

which satisfies detailed balance (see :ref:`the primer <moves-db-primer>`). Each
sub-move within the megamove obeys this independently, so the whole sweep leaves
the Boltzmann distribution invariant.

Configuration
=============

``MOVE_CRANKSHAFT`` : float
    Probability of selecting a crankshaft step (all ``MOVE_*`` must sum to 1.0).

``CRANKSHAFT_SUBSTEPS`` : int
    Sub-moves performed per bead each crankshaft step (e.g. 20 000-50 000, larger
    for big systems). This is the main lever on how much work a crankshaft step
    does.

The crankshaft is one of three moves with a multi-threaded kernel (along with the
slither and pull): see :ref:`PARALLELIZE <advanced-parallel>` (2D or 3D, no frozen
chains). It is fast, ergodic for local relaxation, and a good default to dominate
the move set, mixing in small fractions of the other moves for global
rearrangement.
