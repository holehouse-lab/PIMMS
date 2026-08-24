.. _move-cluster-rotate:

==============
Cluster rotate
==============

:Keyword: ``MOVE_CLUSTER_ROTATE``
:Move code: 8
:Status: core

How it works
============

The rotational counterpart of :doc:`cluster_translate`. PIMMS identifies the
connected cluster of chains containing the selected chain and rotates the whole
cluster as a rigid body: the cluster is first mapped to a single periodic image
(so a boundary-straddling cluster is rotated as a genuine rigid body), then a
cardinal 90/180/270° rotation is applied about the cluster *bead* nearest the
cluster's centroid, with the rotated displacement vectors re-anchored on that
bead's original lattice position - which makes every rotation exactly
invertible. This reorients an entire aggregate at once - a degree of freedom
(whole-cluster tumbling) that single-chain rotations cannot capture.

Why detailed balance holds
==========================

As with cluster translation, two conditions hold:

#. **Symmetric rotation.** The rotation is drawn uniformly from the cardinal
   rotations, whose inverses are equally likely for the reverse move.

#. **Cluster preservation.** The connected component is required to be unchanged
   by the move, so the reverse move selects and rotates the same cluster.

Together these make the proposal symmetric (:math:`g(x\to y)=g(y\to x)`), and only
the cluster-environment interface contributes to :math:`\Delta E` (intra-cluster
contacts are preserved by the rigid rotation). The plain Metropolis acceptance
:math:`A = \min(1, e^{-\Delta E/T})` then satisfies detailed balance (see
:ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CLUSTER_ROTATE`` : float
    Probability of selecting a cluster-rotation step (all ``MOVE_*`` must sum to
    1.0). Like cluster translation, keep this small (e.g. 0.01-0.05).

    Under periodic boundaries this move requires a **cubic** (3D) / **square**
    (2D) production box: a 90° rigid rotation is only an energy-preserving
    symmetry of a cube/square, so the keyfile parser rejects
    ``MOVE_CLUSTER_ROTATE`` with a non-cubic ``DIMENSIONS`` unless ``HARDWALL``
    is on (where a rigid rotation is a valid isometry of any box shape). A
    non-cubic ``RESIZED_EQUILIBRATION`` box is fine - that phase always runs
    with a hardwall.
