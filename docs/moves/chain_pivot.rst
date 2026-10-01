.. _move-chain-pivot:

===========
Chain pivot
===========

:Keyword: ``MOVE_CHAIN_PIVOT``
:Move code: 4
:Status: core
:Scope: the single chain drawn by the main loop

How it works
============

A pivot move picks an interior bead uniformly at random - any bead from index 1 to
:math:`L-2` of an :math:`L`-bead chain - and rigidly **rotates the shorter of the
two arms on either side of it** by a cardinal lattice rotation (90, 180 or 270
degrees, drawn uniformly; in 3D about a randomly chosen x, y or z axis). Two arms of
equal length tie in favour of the N-terminal one. The pivot bead itself stays
exactly where it is and the longer arm is untouched, so the two arms are still
bonded through the pivot when the move lands.

The arm is gathered into a single periodic image, its displacement vectors
relative to the pivot bead are rotated, and the result is re-anchored on the pivot
bead's original in-box position and re-wrapped - the same bead-anchored
construction :doc:`chain_rotate` uses, and for the same reason.

The moving arm is removed from the occupancy grid before its rotated sites are
tested, so it may land on sites it currently occupies, but a clash with the fixed
arm or with any other chain rejects the move. Under ``HARDWALL`` a rotation that
would leave the moved arm (or its bond back to the pivot bead) crossing the wall
is rejected too.

**Chains must be at least 3 beads long**; a monomer or a dimer rejects the move
immediately as a null move (logged as a rejected attempt). Because it moves a large,
contiguous section of the chain in one shot, a pivot makes much larger
conformational changes than a crankshaft - it is an efficient way to decorrelate the
global shape of a chain, especially for swollen/dilute chains. The cost is linear in
chain length.

Why detailed balance holds
==========================

The pivot bead is drawn uniformly from the interior of the chain and the rotation
uniformly from the cardinal rotations. The reverse move draws from exactly the
same sets: the chain length is unchanged, so the same pivot index is equally
likely; which arm moves is a function of the pivot index alone (the shorter one),
so the reverse move moves the same arm; and each rotation's inverse is equally
likely to be drawn. Anchoring on the pivot *bead* (a lattice point that maps to
itself under the rotation, and whose position the move does not change) makes the
inverse rotation return the original configuration exactly, even for an arm that
straddles a periodic boundary. So

.. math::

   g(x\to y) = g(y\to x),

and the move uses the plain Metropolis acceptance
:math:`A = \min(1, e^{-\Delta E/T})` (see :ref:`the primer <moves-db-primer>`).

For :math:`L \ge 4` both termini are reachable. In a 3-mer the only pivot point
sits at the midpoint and the tie always swings the N-terminal arm, so bead 2 never
moves under a pivot - that is what :doc:`head_pivot` is for. More generally the two
beads at the middle of the chain - indices :math:`\lfloor (L-1)/2 \rfloor` and
:math:`\lfloor (L-1)/2 \rfloor + 1`, so beads 1 and 2 of a 3-mer - are the pivot or
part of the longer arm for every pivot point and are never moved by this move,
which is why a move set whose only shape-changing move is the pivot draws a
start-up warning (see :ref:`moves-db-primer`).

In 3D some draws map the arm exactly onto itself: an arm whose beads all lie on one
lattice axis through the pivot bead is unchanged by the three rotations about that
axis. Every pivot of a 3- or 4-mer moves a single terminal bead, so any axis-aligned
terminal bond gives such draws. Unlike :doc:`chain_rotate` and
:doc:`cluster_rotate`, the pivot does not screen these out. Such a draw passes the
clash test, is accepted with :math:`\Delta E = 0` and is logged as an *accepted*
pivot in ``ACCEPTANCE.dat``, although nothing moved. The sampled ensemble is
unaffected, but for short or rod-like chains in 3D the code-4 acceptance ratio
overstates how often a pivot actually changes the chain - for a straight,
axis-aligned chain of any length, or a 3- or 4-mer whose terminal bonds lie along
lattice axes, a third of all draws are such identities. In 2D every cardinal
rotation moves a non-empty arm, so the issue does not arise.

Configuration
=============

``MOVE_CHAIN_PIVOT`` : float
    Probability of selecting a pivot step (all ``MOVE_*`` must sum to 1.0).
    Default 0.0.

No other tuning keywords, and no parallel kernel. Pivots have a low acceptance
rate in dense systems (the rotated arm usually clashes) but are very effective for
single chains and dilute solutions, where a modest fraction substantially speeds
up sampling of the overall chain dimensions.
