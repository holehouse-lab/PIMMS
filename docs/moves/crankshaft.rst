.. _move-crankshaft:

==========
Crankshaft
==========

:Keyword: ``MOVE_CRANKSHAFT``
:Move code: 1
:Status: core (recommended workhorse)
:Scope: whole system (megamove)

How it works
============

The crankshaft is PIMMS' fundamental local move and should make up most of a
typical move budget. It perturbs **one bead at a time**, and the set of candidate
sites depends on how that bead is bonded:

* an **interior** bead is proposed a site drawn uniformly from the intersection of
  its two bonded neighbours' Chebyshev-1 neighbourhoods, so the bead stays bonded
  to both;
* a **terminal** bead is proposed a site drawn uniformly from the full
  :math:`3^d` block around its single bonded anchor (:math:`d` = 2 or 3), which
  keeps it bonded to that anchor;
* a **single-bead chain** (a monomer) has no anchor at all, so the block is
  centred on its own site and the move is a local translation.

In every case the bead's current site is inside the candidate set, and a draw that
lands on an occupied site - the bead's own current site, the anchor's site, or any
other bead - is rejected as a hard-sphere clash rather than redrawn. The chain's
identity and bonding are preserved; only the kink at that bead changes.

Under ``HARDWALL`` a proposal is additionally rejected if any bond it would create
crosses the wall (for a monomer, if the step itself would cross it). Without that
check the periodic wrap inside the proposal would let a monomer walk straight
through a wall.

A single crankshaft *step* is a **megamove**: it performs
``CRANKSHAFT_SUBSTEPS`` such single-bead perturbations in total, each one targeting
a bead drawn uniformly at random (with replacement) from all non-frozen beads in
the system, and each with its own accept/reject, all inside an optimised Cython
kernel. With a large ``CRANKSHAFT_SUBSTEPS`` one step therefore encompasses many
thousands of elementary Monte Carlo moves. Frozen chains are never selected but
remain in place as energy-contributing obstacles. (Under ``PARALLELIZE`` the beads
are instead drawn from the block interiors - see `Performance`_ below.)

The same kernel drives the single-chain "shake" used inside :doc:`jump_and_relax`
and the chain and multi-chain :doc:`tsmmc` excursions, restricted to the beads of
the chain(s) being perturbed. Those shakes always use the serial kernel, whatever
``PARALLELIZE`` says.

Every sub-move of a crankshaft step is logged under code 1 in ``MOVE_FREQS.dat`` /
``ACCEPTANCE.dat``, including the ones rejected outright as clashes or hardwall
violations, so the code-1 counters advance by ``CRANKSHAFT_SUBSTEPS`` per crankshaft
step (see `Performance`_ for the one exception under ``PARALLELIZE``). The shakes
inside jump-and-relax and the TSMMC excursions, and crankshaft steps run inside a
system-wide TSMMC excursion, are counted only in ``TOTAL_MOVES.dat``.

Why detailed balance holds
==========================

For a bonded bead the set of valid destination sites is determined solely by its
bonded neighbours' (fixed) positions, not by the bead's current position. Hence
the forward proposal "bead at :math:`A \to B`" and the reverse "bead at
:math:`B \to A`" are drawn from the *same* set of size :math:`N`, so

.. math::

   g(x\to y) = g(y\to x) = \tfrac{1}{N}.

A monomer is the one case where the candidate block travels with the bead, and the
proposal is symmetric there too: both blocks hold :math:`3^d` sites, and :math:`B`
lies in the block centred on :math:`A` exactly when :math:`A` lies in the block
centred on :math:`B`, so :math:`g(x\to y) = g(y\to x) = 3^{-d}`.

The proposal is symmetric, and the move is accepted with the plain Metropolis
criterion

.. math::

   A(x\to y) = \min\!\left(1,\; e^{-\Delta E / T}\right),

which satisfies detailed balance (see :ref:`the primer <moves-db-primer>`). Here
:math:`\Delta E` is the sum of the pair-interaction change (short-range plus, for
long-range beads, the Chebyshev-2 and Chebyshev-3 shells) and the change in the
angle penalty at the moved bead and its bonded neighbours. Each sub-move within
the megamove obeys this independently, so the whole sweep leaves the Boltzmann
distribution invariant.

Configuration
=============

``MOVE_CRANKSHAFT`` : float
    Probability of selecting a crankshaft step (all ``MOVE_*`` must sum to 1.0).
    Default 0.0.

``CRANKSHAFT_SUBSTEPS`` : int
    Total number of single-bead perturbations performed per crankshaft step; must
    be a positive integer (default 500; 20 000-50 000 is a more typical production
    value, larger for big systems). They are spread at random across all non-frozen
    beads in the system, so on average each bead is perturbed
    ``CRANKSHAFT_SUBSTEPS`` / (number of movable beads) times per step. This is the
    main lever on how much work a crankshaft step does, and it also sizes the two
    relaxations of :doc:`jump_and_relax`.

``CRANKSHAFT_MODE`` : str
    **Obsolete and ignored.** It used to scale the substep count with chain
    length. The keyword is still accepted so that old keyfiles keep working, but
    PIMMS prints a warning and ignores its value (any value is accepted and reset
    to ``UNIFORM``): the crankshaft always performs a fixed ``CRANKSHAFT_SUBSTEPS``
    sub-moves per megamove, independent of chain length (the old ``UNIFORM``
    mode). Remove it from new keyfiles.

Performance
===========

The crankshaft is one of three moves with a multi-threaded kernel (along with the
:doc:`slither` and :doc:`pull`): see :ref:`PARALLELIZE <advanced-parallel>`. It
works in 2D and 3D and composes with frozen chains, which the parallel kernel
honours through a per-bead frozen mask. The parallel path is taken only when the
box actually splits into more than one block; a single-block box uses the serial
kernel directly rather than paying the threading overhead.

The parallel kernel decomposes the box into blocks separated by a frozen halo of
width :math:`W` on each blocked face (:math:`W = 2` when any bead carries a
long-range flag, otherwise :math:`W = 1`). A move must both start and land at least
:math:`W` inside its block, so writes never leave the interior while reads may reach
into the neighbouring block's read-only halo; since the next block's interior starts
:math:`2W + 1` sites past the last interior site of this one (a width-:math:`W` halo
on each side of the face), the safety condition is
:math:`2W \ge \max(R_\text{int}, 2)`, with :math:`R_\text{int} = 3` for long-range
systems and 1 otherwise. Blocks are at least :math:`8W` sites long (up to four per
axis), so an axis is only split once it is :math:`16W` sites or longer - 16 without
long-range beads, 32 with them - and the block grid's origin is shifted at random on
every sweep so the halos fall somewhere different each time. Within a block the
beads to move are drawn uniformly, with replacement, from the block's movable beads
(frozen chains excluded). The ``CRANKSHAFT_SUBSTEPS`` attempts are shared between
blocks in proportion to their movable bead counts, with the floored shares topped up
largest-remainder-first so that exactly the requested number of attempts is made -
provided anything is movable at all. On a sweep whose shift leaves no movable bead
inside any block interior the kernel makes no attempts, and that step logs zero
attempts under code 1 rather than ``CRANKSHAFT_SUBSTEPS``.

``PARALLELIZE`` leaves the equilibrium distribution alone, but it is a *different*
Markov chain: only the block interiors move in a given sweep and the kernel draws
from a different random stream, so a parallel run relaxes more slowly per step and
does not reproduce a serial run at the same ``SEED``. Compare equilibrium averages,
not energies at a fixed step.

The crankshaft is fast, ergodic for local relaxation, and a good default to
dominate the move set, mixing in small fractions of the other moves for global
rearrangement.
