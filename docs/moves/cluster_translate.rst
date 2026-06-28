.. _move-cluster-translate:

=================
Cluster translate
=================

:Keyword: ``MOVE_CLUSTER_TRANSLATE``
:Move code: 7
:Status: core

How it works
============

A cluster move operates on a whole **connected cluster of chains** at once. PIMMS
first identifies the cluster: starting from the selected chain it grows the set of
all chains connected to it through contacts (a connected component). The entire
cluster is then translated as a rigid body by a random displacement, exactly like
a single-chain translation but applied to every chain in the cluster together.

This lets a whole aggregate or droplet diffuse as a unit - motion that no
single-chain move can produce. It is comparatively expensive (identifying and
moving the cluster, and checking for clashes against the rest of the system), so
it is normally used at a small fraction of the move budget.

Why detailed balance holds
==========================

Two ingredients make the move balanced:

#. **Symmetric displacement.** As for :doc:`chain_translate`, a shift by
   :math:`+v` and the reverse :math:`-v` are equally likely.

#. **Cluster preservation.** After the rigid move PIMMS re-identifies the
   connected component and requires it to be the *same* set of chains. This
   guarantees that the reverse move would select and translate the identical
   cluster, keeping the proposal symmetric (:math:`g(x\to y)=g(y\to x)`) and the
   move reversible.

Because the cluster moves rigidly, all *intra*-cluster contacts are unchanged;
only the cluster-environment interface contributes to :math:`\Delta E`. The move
is accepted with the plain Metropolis criterion
:math:`A = \min(1, e^{-\Delta E/T})` on that interfacial energy change, satisfying
detailed balance (see :ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CLUSTER_TRANSLATE`` : float
    Probability of selecting a cluster-translation step (all ``MOVE_*`` must sum
    to 1.0). Keep small (e.g. 0.01-0.05); these moves are expensive.

.. note::

   The cluster here is the **geometric** connected component, so in a single
   fully-connected droplet the "cluster" is the entire condensate and the move
   simply diffuses it bodily. To rearrange chains *within* a dense phase, use the
   energy-gradient collective move :doc:`vmmc`, which recruits and moves
   *sub-clusters*.
