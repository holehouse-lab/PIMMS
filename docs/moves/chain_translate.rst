.. _move-chain-translate:

===============
Chain translate
===============

:Keyword: ``MOVE_CHAIN_TRANSLATE``
:Move code: 2
:Status: core
:Scope: the single chain drawn by the main loop

How it works
============

The whole chain is moved as a **rigid body** by a single random displacement
vector. The offset in each dimension is drawn uniformly from ``0`` to
``DIMENSIONS[d] - 1``, so this is not a small local step: the chain is relocated
to a uniformly random position anywhere in the box, with periodic wrapping. (A
zero offset in every dimension is one of the possible draws - one in
:math:`L_x L_y` or :math:`L_x L_y L_z` - and simply leaves the chain where it is.
Unlike the identity draws of :doc:`chain_rotate`, it is not singled out: it passes
the clash test, is accepted with :math:`\Delta E = 0` and is logged as an accepted
translation in ``ACCEPTANCE.dat``. In any realistic box the effect on the
acceptance ratio is negligible.)

The chain is removed from the occupancy grid before the translated sites are
tested, so it may land on sites it currently occupies. If any translated bead
would land on a site occupied by *another* chain the move is rejected as a
hard-sphere clash. Under ``HARDWALL`` the move is additionally rejected if the
translated chain would have a bond crossing the wall (consecutive beads more than
one site apart), since bonds cannot pass through a wall. The chain's internal
conformation is unchanged - only its position in the box moves.

Note what the hardwall check does and does not forbid. It rejects any *state* in
which the chain straddles a face, but the offset is still drawn over the whole box
and still wraps, so a chain hugging one wall can be relocated in a single move to
the opposite wall. Both the state it leaves and the state it lands in are legal
confined states, and the wrapped shift is its own inverse in the same sense a
periodic shift is, so this costs nothing in correctness: the sampled ensemble is
the confined Boltzmann distribution. It does mean the move is a relocation rather
than a physical passage through the wall, and a hardwall trajectory will show
chains jumping from one side of the box to the other. Set ``MOVE_CHAIN_TRANSLATE``
to 0 if every step of your run needs to be spatially local. The collective moves
take the stricter line and refuse to wrap at all, because a rigid cluster
translation that wrapped would not be a rigid motion of the confined system.

There is no chain-length restriction: a single-bead chain translates like any
other. Frozen chains are never selected by the main loop, so they are never
translated. The move works the same way in 2D and 3D and in non-cubic boxes (the
offset on each axis is drawn over that axis' own length). Each step logs one
attempt under code 2.

This is one of the more expensive single-chain moves, because the energy of every
bead has to be re-evaluated against its new surroundings (the cost is linear in
chain length), but it is the primary way an intact chain explores different
regions of the box.

Why detailed balance holds
==========================

The displacement is drawn uniformly over the whole box, and on a periodic lattice
a shift by :math:`+v` and the reverse shift by :math:`-v \equiv L - v` are equally
likely draws, so the proposal is symmetric,

.. math::

   g(x\to y) = g(y\to x).

The move is therefore accepted with the plain Metropolis criterion
:math:`A = \min(1, e^{-\Delta E/T})`, which satisfies detailed balance (see
:ref:`the primer <moves-db-primer>`). Hard-sphere clashes correspond to
:math:`\pi = 0` states and are rejected with certainty, consistent with balance;
the hardwall bond check is likewise symmetric, since a configuration with a bond
through the wall is not a legal state in either direction. That symmetry is what
licenses the wrapping described above: a hardwall chain against one face and the
same chain against the opposite face are connected by a shift and its inverse with
equal probability, so balance between them holds and both are weighted by their own
hardwall energy.

Configuration
=============

``MOVE_CHAIN_TRANSLATE`` : float
    Probability of selecting a chain-translation step (all ``MOVE_*`` must sum to
    1.0). Default 0.0.

There are no other tuning keywords, and the move has no parallel kernel. In dense
systems most translations clash and are rejected; for relocating chains through
crowded/condensed phases prefer the collective moves (:doc:`vmmc`, :doc:`pull`)
instead.
