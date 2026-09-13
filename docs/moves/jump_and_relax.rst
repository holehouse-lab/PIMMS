.. _move-jump-and-relax:

==============
Jump and relax
==============

:Keyword: ``MOVE_JUMP_AND_RELAX``
:Move code: 13
:Status: stable
:Scope: the single chain drawn by the main loop

How it works
============

Jump-and-relax concentrates sampling effort on relocating a single chain. It is a
**composite** of three sub-steps applied to the chain the main loop drew:

#. **relax** - a single-chain crankshaft "shake": ``CRANKSHAFT_SUBSTEPS`` local
   perturbations of that chain, each accepted or rejected on its own inside the
   Cython kernel. Always committed.
#. **jump** - a rigid :doc:`chain_translate` of the whole chain, accepted or
   rejected on its own Metropolis criterion (and reverted outright on a
   hard-sphere clash or a hardwall violation).
#. **relax** - a second single-chain shake in the chain's (possibly new) location.
   Always committed.

The relaxations let the chain explore conformations before and after the
relocation attempt, so a chain that lands somewhere viable can settle into its new
surroundings. There is no chain-length restriction, and the move works under
periodic and hardwall boundaries alike.

The accept/reject recorded for move code 13 in ``MOVE_FREQS.dat`` /
``ACCEPTANCE.dat`` is that of the **jump** - the two relaxations are always kept -
and the relaxation sub-moves are counted through the alternative-Markov-chain
counter that feeds ``TOTAL_MOVES.dat``.

Why detailed balance holds
==========================

The argument is **composition**, not a single acceptance test. Each of the three
sub-steps is, on its own, a valid Monte Carlo update that leaves the Boltzmann
distribution :math:`\pi` invariant:

* the relaxation shakes are ordinary single-chain crankshaft sub-chains, each
  obeying detailed balance (:doc:`crankshaft`);
* the jump is a standard single-chain translation accepted with the Metropolis
  criterion :math:`\min(1, e^{-\Delta E/T})` (:doc:`chain_translate`).

A sequence of :math:`\pi`-invariant kernels :math:`K_1, K_2, K_3` has :math:`\pi`
as a stationary distribution of the product :math:`K_1 K_2 K_3`,

.. math::

   (\pi K_1 K_2 K_3)(y) = \pi(y),

so the composite move samples the correct distribution.

Configuration
=============

``MOVE_JUMP_AND_RELAX`` : float
    Probability of selecting a jump-and-relax step (all ``MOVE_*`` must sum to
    1.0). Default 0.0.

``CRANKSHAFT_SUBSTEPS`` : int
    Reused to size the two relaxation shakes (default 500), so one jump-and-relax
    step costs ``2 x CRANKSHAFT_SUBSTEPS`` single-bead sub-moves on the selected
    chain plus one translation.

Performance
===========

There is no parallel kernel for this move. The jump's energy change is obtained
from a **from-scratch total-energy recompute** rather than an incremental local
evaluation, which makes each jump-and-relax step markedly more expensive than a
bare :doc:`chain_translate` in a large system - budget for that when choosing its
fraction.

Because the jump is accepted on its own energy, it does not "rescue" a chain
that lands in a bad spot; for aggressive relocation through dense/condensed phases
prefer :doc:`vmmc` or :doc:`pull`, which are built for that.
