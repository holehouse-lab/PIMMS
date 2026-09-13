.. _move-cluster-rotate:

==============
Cluster rotate
==============

:Keyword: ``MOVE_CLUSTER_ROTATE``
:Move code: 8
:Status: core
:Scope: the connected cluster containing the chain drawn by the main loop

How it works
============

The rotational counterpart of :doc:`cluster_translate`. PIMMS identifies the
connected cluster of chains containing the chain the main loop drew - the same
Chebyshev-1 contact criterion, and the same rejections if the cluster contains a
frozen chain or grows to contain every chain in the system - and rotates the whole
cluster as a rigid body by a cardinal 90/180/270 degree rotation (in 3D about a
randomly chosen x, y or z axis).

Under periodic boundaries the cluster is first mapped into a single periodic image
so a boundary-straddling cluster is rotated as a genuine rigid body. The rotation
is applied to the displacement vectors relative to the cluster *bead* nearest the
cluster's centroid, and the result is re-anchored on that bead's original in-box
lattice position - which makes every rotation exactly invertible. This reorients
an entire aggregate at once, a degree of freedom (whole-cluster tumbling) that
single-chain rotations cannot capture.

Constraints
===========

**Draws that do nothing.** A rotation that maps the cluster exactly onto itself is
rejected rather than reported as an accepted rotation. The case that matters in
practice is an isolated free monomer: its cluster is a single bead, every rotation
is about that bead, and so every draw is an identity. Straight, axis-aligned
clusters lose the rotations about their own axis the same way. Rejecting these does
not change what is sampled, since accepting an identity commits a bit-identical
state, but before 1.0.8 they were logged as accepted, which in a system with free
monomers inflated the apparent acceptance of this move by roughly a factor of two
and a half. If your system has free monomers, note that a cluster rotation drawn on
one can never do anything, and budget the move fraction accordingly.

**Box shape (periodic boundaries only).** Under periodic boundaries this move
requires a **cubic** (3D) / **square** (2D) production box: a 90 degree rotation
swaps two axes, and on unequal periodic axes that swap changes intra-cluster
minimum-image distances, so the rotation stops being energy preserving and the
tracked energy would drift from the true energy. The keyfile parser therefore
rejects ``MOVE_CLUSTER_ROTATE`` with a non-cubic ``DIMENSIONS`` unless ``HARDWALL``
is on (where a rigid rotation is a valid isometry of any box shape). A non-cubic
``RESIZED_EQUILIBRATION`` box is fine - that phase always runs with a hardwall.
The check uses the *final* ``HARDWALL`` and ``DIMENSIONS``, so a restart file that
overrides either of them is validated too.

**Winding clusters (periodic boundaries only).** A cluster that winds around the
box - one that is connected to its own periodic image, so its single-image extent
reaches the box length on some axis - is rejected outright. A cardinal rotation of
such a cluster is not a rigid motion of the periodic system: the winding closure
vector lands on an axis with a different period, so intra-cluster minimum-image
relations change while the move assumes they do not. The rejection is symmetric
(winding is preserved by the move), so detailed balance is unaffected.

**Hardwall.** Under ``HARDWALL`` the cluster's coordinates are used as the plain
Cartesian coordinates they are: no periodic single-image gather and no winding
guard is applied, because a cluster spanning a box axis is a legal, rotatable
configuration there. Instead, a rotation whose raw coordinates would leave the box
is rejected before any wrapping, and a chain left with a bond crossing the wall is
rejected as well.

Why detailed balance holds
==========================

As with cluster translation, two conditions hold:

#. **Symmetric rotation.** The rotation is drawn uniformly from the cardinal
   rotations, whose inverses are equally likely for the reverse move, and it is
   anchored on a cluster bead (a lattice point that maps to itself) so the inverse
   rotation returns the original configuration exactly. The pivot bead is selected
   with exact-integer arithmetic over the cluster's beads concatenated in
   ascending chain-ID order, so both the tie-break and the concatenation order are
   the same in the forward and reverse directions.

#. **Cluster preservation.** The connected component from the same seed chain is
   required to be the same *set* of chain IDs after the move, so the reverse move
   selects and rotates the same cluster.

Together these make the proposal symmetric (:math:`g(x\to y)=g(y\to x)`). As for
cluster translation, no non-cluster bead is Chebyshev-1 adjacent to the cluster
before or after an accepted move, so the **short-range energy change is exactly
zero** and only the long-range Chebyshev-2/3 interfacial terms contribute to
:math:`\Delta E` (with no long-range residues the move is exactly energy neutral
and is committed without an energy evaluation). The plain Metropolis acceptance
:math:`A = \min(1, e^{-\Delta E/T})` then satisfies detailed balance (see
:ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CLUSTER_ROTATE`` : float
    Probability of selecting a cluster-rotation step (all ``MOVE_*`` must sum to
    1.0). Default 0.0. Like cluster translation, keep this small (e.g.
    0.01-0.05) - the cost grows with the size of the cluster.

No other tuning keywords, and no parallel kernel.
