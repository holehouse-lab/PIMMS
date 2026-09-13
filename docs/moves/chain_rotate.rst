.. _move-chain-rotate:

============
Chain rotate
============

:Keyword: ``MOVE_CHAIN_ROTATE``
:Move code: 3
:Status: core
:Scope: the single chain drawn by the main loop

How it works
============

The whole chain is rotated as a **rigid body** about one of its own beads - the
bead nearest the chain's (single-image) centroid. On the lattice only the cardinal
rotations are used - 90, 180 or 270 degrees, in 3D about a randomly chosen x, y or
z axis - because arbitrary angles do not map lattice sites onto lattice sites.

Concretely, the chain is first gathered into a single periodic image, the
displacement vectors of every bead relative to the pivot bead are rotated, and the
result is re-anchored on the pivot bead's original in-box lattice position and
re-wrapped:

.. math::

   \text{rotated}[i] = \mathrm{pbc}\bigl(p_\text{pivot} +
   R\,(s_i - s_\text{pivot})\bigr),

where :math:`s` are the single-image coordinates and :math:`p_\text{pivot}` the
pivot bead's original position. This makes every rotation exactly invertible for
any box shape, including chains straddling a periodic boundary.

As with translation, a rotated bead landing on a site occupied by another chain
rejects the move; under ``HARDWALL`` a rotation that would leave the chain with a
bond crossing the wall is rejected too. The internal conformation is preserved;
only the chain's orientation changes.

A **single-bead chain cannot be rotated**, so the move returns immediately as a
rejected null move. If the system contains monomers and this move is enabled, the
start-up summary says what fraction of steps that will cost (see
:ref:`moves-step-anatomy`).

The same applies to any draw that maps the chain exactly onto itself, which is
what the three rotations about its own axis do to a straight, axis-aligned chain.
Those draws are rejected rather than being reported as accepted rotations. This
does not change what is sampled - accepting an identity commits a bit-identical
state, and detailed balance says nothing about a move from a state to itself - but
it does keep ``ACCEPTANCE.dat`` honest. Before 1.0.8 they were counted as accepted,
which inflated the apparent acceptance of this move for rod-like chains.

Unlike :doc:`cluster_rotate`, this move places no restriction on the box shape: it
evaluates the true energy change of the rotated chain rather than assuming the
rotation is energy neutral, so a non-cubic periodic box is fine.

Why detailed balance holds
==========================

The proposed rotation is chosen uniformly from the cardinal rotations, and each
rotation's inverse (e.g. 90 and 270, 180 and itself) is equally likely to be
proposed for the reverse move. Anchoring the rotation on a *bead* matters here: a
bead is a lattice point that maps exactly to itself under the rotation, and the
bead nearest the centroid is a rotation-invariant choice for a rigid body, so the
reverse move picks the same physical pivot and the inverse rotation returns the
original configuration bit for bit.

The pivot is chosen with exact-integer arithmetic - the bead minimising
:math:`\lVert n\,p_i - \sum_j p_j\rVert^2`, first index wins - rather than by
comparing floating-point distances to the centroid. Two beads exactly equidistant
from the centre would otherwise have their tie broken by rounding noise that is
not invariant under the rotation, so the reverse move could pick the *other* tied
bead and fail to invert.

With a reversible, uniformly-drawn rotation the proposal is symmetric,

.. math::

   g(x\to y) = g(y\to x),

and the move is accepted with the plain Metropolis criterion
:math:`A = \min(1, e^{-\Delta E/T})`, satisfying detailed balance (see
:ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CHAIN_ROTATE`` : float
    Probability of selecting a chain-rotation step (all ``MOVE_*`` must sum to
    1.0). Default 0.0.

No other tuning keywords, and no parallel kernel. Like translation, rotation is
most useful in dilute systems; in dense phases most rotations clash. The cost is
linear in chain length.
