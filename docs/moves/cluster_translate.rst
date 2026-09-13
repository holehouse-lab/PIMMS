.. _move-cluster-translate:

=================
Cluster translate
=================

:Keyword: ``MOVE_CLUSTER_TRANSLATE``
:Move code: 7
:Status: core
:Scope: the connected cluster containing the chain drawn by the main loop

How it works
============

A cluster move operates on a whole **connected cluster of chains** at once. PIMMS
first identifies the cluster: starting from the chain the main loop drew, it grows
the set of all chains reachable through **Chebyshev-1 contact** (the Moore shell -
two chains are linked if any bead of one is a neighbouring site of any bead of the
other), which under ``HARDWALL`` never links chains through a wall. Note this is
purely geometric contact, not an interaction-energy criterion.

The entire cluster is then translated as a rigid body by a random displacement.
The offset in each dimension has a magnitude drawn uniformly from ``1`` to
``DIMENSIONS[d] - 1`` and a random sign, so (unlike :doc:`chain_translate`) the
displacement is never zero on any axis and the cluster is relocated somewhere else
in the box.

The whole cluster is removed from the occupancy grid *before* any translated bead
is placed, so cluster-mates never clash with each other; a translated bead landing
on a non-cluster bead rejects the move. Under ``HARDWALL`` a raw translated
coordinate outside the box is rejected outright rather than wrapped, and a chain
left with a bond crossing the wall is rejected too.

Two further conditions reject the move before it is evaluated:

* the cluster must not contain a **frozen** chain (frozen chains cannot move, so
  the whole cluster is refused);
* the cluster must not exceed the size threshold. The main loop sets that
  threshold to (number of chains) - 1, so a cluster that grows to contain *every*
  chain in the system is never translated.

This lets a whole aggregate or droplet diffuse as a unit - motion that no
single-chain move can produce. It is comparatively expensive (identifying and
moving the cluster, and checking for clashes against the rest of the system), so
it is normally used at a small fraction of the move budget.

Why detailed balance holds
==========================

Two ingredients make the move balanced:

#. **Symmetric displacement.** A shift by :math:`+v` and the reverse :math:`-v`
   are equally likely draws.

#. **Cluster preservation.** After the rigid move PIMMS re-identifies the
   connected component from the same seed chain and requires it to be the *same
   set* of chain IDs (not merely the same size - a cluster could swap members
   across a periodic boundary at fixed size). This guarantees that the reverse
   move would select and translate the identical cluster, keeping the proposal
   symmetric (:math:`g(x\to y)=g(y\to x)`) and the move reversible. A translation
   that merges the cluster with another one is irreversible - no single move
   un-merges them - so it is rejected.

Because the cluster is the *complete* set of chains it touches, no non-cluster
bead is Chebyshev-1 adjacent to it before or after an accepted move, so every
short-range partner is solvent in both configurations and the **short-range energy
change is exactly zero** (walls count as solvent, so this holds under ``HARDWALL``
too). Only the long-range Chebyshev-2/3 terms across the cluster interface can
change, and those are what :math:`\Delta E` is built from. In a system with no
long-range residues the move is exactly energy neutral and is committed without an
energy evaluation at all. Otherwise the plain Metropolis criterion
:math:`A = \min(1, e^{-\Delta E/T})` on the long-range interfacial change satisfies
detailed balance (see :ref:`the primer <moves-db-primer>`).

Configuration
=============

``MOVE_CLUSTER_TRANSLATE`` : float
    Probability of selecting a cluster-translation step (all ``MOVE_*`` must sum
    to 1.0). Default 0.0. Keep small (e.g. 0.01-0.05); these moves are expensive.

No other tuning keywords, and no parallel kernel.

.. note::

   Because a cluster that grows to contain every chain in the system is rejected
   rather than moved (translating the whole system is a no-op for sampling), in a
   single fully condensed system this move never fires. To rearrange chains
   *within* a dense phase, use the energy-gradient collective move :doc:`vmmc`,
   which recruits and moves *sub-clusters*.
