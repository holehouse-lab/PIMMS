.. _move-slither:

===================
Slither (reptation)
===================

:Keyword: ``MOVE_SLITHER``
:Move code: 6
:Status: core
:Scope: whole system (megamove)

How it works
============

A slither move is **reptation** - the chain advances forwards or backwards through
the lattice "like a snake". A direction is chosen; a new site is drawn for the
leading end from the :math:`3^d` block around it, and if that site is free the
whole chain shifts one place along its own contour: the bead at the trailing end
vacates its site, every interior bead takes the site of its neighbour towards the
head, and the leading bead occupies the new site. The chain crawls a step without
any large rigid motion, which relaxes chain conformations and lets chains thread
through crowded surroundings very efficiently.

A slither *step* is a **megamove**: when selected, every non-frozen chain is
slithered ``SLITHER_SUBSTEPS`` times, in one globally shuffled order (so a
megamove is ``SLITHER_SUBSTEPS`` x (number of non-frozen chains) sub-moves), each
substep with its own accept/reject, in an optimised Cython kernel. Frozen chains
are never selected but stay in place as energy-contributing obstacles. A
single-bead chain has no contour to reptate along, so its slither degenerates to a
local translation.

Under ``HARDWALL`` a proposed end site that is only reachable through the periodic
wrap is rejected, so a chain can neither grow nor walk through a wall.

Why detailed balance holds
==========================

The slither proposal is **symmetric**. The direction (head-first or tail-first)
is drawn with probability :math:`1/2`, and the new end site is drawn uniformly
from the *fixed-size* Chebyshev-1 offset box around the growing end
(:math:`3^d` sites, the current site included); a draw landing on an occupied
site is rejected as a hard-sphere clash rather than redrawn. Because the
candidate set has the same fixed size in both directions,
:math:`g(x\to y) = g(y\to x) = \tfrac12 \cdot 3^{-d}`, and the plain
Metropolis criterion

.. math::

   A(x\to y) = \min\!\left(1,\; e^{-\Delta E / T}\right)

satisfies detailed balance (see :ref:`the primer <moves-db-primer>`). The
energy change :math:`\Delta E` is cheap to evaluate for a **homopolymer** - only
the removed and added end beads change their surroundings, an :math:`O(1)`
calculation - but for a **heteropolymer** the bead identities shift relative to
their positions, so every bead's interactions change and :math:`\Delta E` is an
:math:`O(N)` calculation, assembled as a telescoping sum of :math:`N` single-bead
moves and reverted bead by bead on a rejection. The kernel implements both paths;
the result is checked by the detailed-balance test suite (a crankshaft-only run
and a crankshaft+slither run reach the same equilibrium).

Configuration
=============

``MOVE_SLITHER`` : float
    Probability of selecting a slither megamove (all ``MOVE_*`` must sum to 1.0).
    Default 0.0.

``SLITHER_SUBSTEPS`` : int
    Number of slither moves applied to each non-frozen chain, in random order, per
    megamove (default 10).

Performance
===========

Slither works in both 2D and 3D and is one of the most effective moves for
relaxing chain conformations; a healthy fraction alongside the crankshaft usually
improves mixing markedly.

Along with the crankshaft and pull, it is one of the three moves with a
multi-threaded kernel: see :ref:`PARALLELIZE <advanced-parallel>`. Because a
slither is a whole-chain move, the parallel kernel decomposes at chain level: a
chain is movable inside a block only if *all* of its beads lie at least
:math:`W = R_\text{int} + 2` sites inside that block (:math:`W = 5` with
long-range beads, 3 without), and a proposal whose new end would enter the frozen
halo is rejected. A chain straddling a block boundary is simply frozen for that
sweep.

A chain can only ever fit a block interior if it is short enough, so PIMMS splits
the chains by *length* once, at the start of the run: chains no longer than the
smallest block interior (and within the kernel's 512-bead per-chain buffer, which
for the slither applies to heteropolymers) go to the parallel kernel, and every
longer chain goes to the serial kernel. Both passes run in every megamove. The
split depends only on the chain lengths, the box and the frozen set, all of which
are fixed for the run, so it never depends on the configuration the system happens
to be in. If the box decomposes into a single block the whole megamove is serial.
The per-block share of the substeps is topped up largest-remainder-first so exactly
the requested number of attempts is made.

The split has to be by length rather than by the chains' current extents, even
though a long chain is usually compact enough to fit. Choosing per megamove from
the current configuration would make the choice of kernel a function of the state,
and because the parallel kernel can never push a chain past a block interior while
the serial kernel can, that choice quietly biases the run towards compact chains.
That was the behaviour up to 1.0.8; see :ref:`PARALLELIZE <advanced-parallel>` for
the measured size of the effect.

``PARALLELIZE`` leaves the equilibrium distribution alone, but it is a *different*
Markov chain: only the block interiors move in a given sweep and the kernel draws
from a different random stream, so a parallel run relaxes more slowly per step and
does not reproduce a serial run at the same ``SEED``. Compare equilibrium averages,
not energies at a fixed step.
