.. _move-jump-and-relax:

==============
Jump and relax
==============

:Keyword: ``MOVE_JUMP_AND_RELAX``
:Move code: 13
:Status: stable

How it works
============

Jump-and-relax concentrates sampling effort on relocating a single chain. It is a
**composite** of three sub-steps applied to the selected chain:

#. **relax** - a single-chain crankshaft "shake" (many local perturbations of the
   chain, each accepted/rejected on its own);
#. **jump** - a rigid translation of the whole chain, accepted or rejected on its
   own merit;
#. **relax** - a second single-chain shake in the chain's (possibly new) location.

The relaxations let the chain explore conformations before and after the
relocation attempt, so a chain that lands somewhere viable can settle into its new
surroundings.

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

.. admonition:: A subtle history

   Earlier versions of this move applied the jump *unconditionally* and deferred a
   **single** accept/reject to the energy *after* both relaxations
   (:math:`E_\text{final}` vs :math:`E_\text{initial}`). Because the relaxations
   bias the proposal, that composite is an *asymmetric* proposal, and a plain
   Metropolis test on :math:`E_\text{final}-E_\text{initial}` is **not** a valid
   Hastings acceptance - it over-accepts downhill moves and breaks detailed
   balance (the bias grows with the relaxation strength). Accepting/rejecting the
   jump on its own merit, between two :math:`\pi`-preserving relaxations, is what
   makes the move correct.

Configuration
=============

``MOVE_JUMP_AND_RELAX`` : float
    Probability of selecting a jump-and-relax step (all ``MOVE_*`` must sum to
    1.0).

``CRANKSHAFT_SUBSTEPS`` : int
    Reused to size the two relaxation shakes.

Because the jump is now accepted on its own energy, it does not "rescue" a chain
that lands in a bad spot; for aggressive relocation through dense/condensed phases
prefer :doc:`vmmc` or :doc:`pull`, which are built for that.
