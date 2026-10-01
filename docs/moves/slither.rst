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
vacates its site, each remaining bead takes the site of its neighbour on the leading
side, and the leading bead occupies the new site. The chain crawls a step without
any large rigid motion, which relaxes chain conformations and lets chains thread
through crowded surroundings very efficiently.

A slither *step* is a **megamove**: when selected, every non-frozen chain is
slithered ``SLITHER_SUBSTEPS`` times, in one globally shuffled order (so a
megamove is ``SLITHER_SUBSTEPS`` x (number of non-frozen chains) sub-moves), each
substep with its own accept/reject, in an optimised Cython kernel. The chain the
main loop drew for the step plays no part. Frozen chains are never selected but
stay in place as energy-contributing obstacles. A single-bead chain has no contour
to reptate along, so its slither degenerates to a local translation: a site drawn
uniformly from the :math:`3^d` block around the bead, exactly as the
:doc:`crankshaft` moves a monomer. (Under ``PARALLELIZE`` the chains handed to the
parallel kernel are scheduled differently - see `Performance`_.)

Every substep is logged under code 6 in ``MOVE_FREQS.dat`` / ``ACCEPTANCE.dat``,
including draws rejected outright because the new end site was occupied (which
includes the site the trailing end is about to vacate) or lay across a wall.

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
the result is checked by the detailed-balance test suite, which runs slither-only
megamoves and crankshaft-only megamoves from the same equilibrated configuration
and requires the same equilibrium (with a deliberately mis-tempered slither as a
positive control that must be caught), and by forward/reverse transition counting
of single slither sub-moves.

Configuration
=============

``MOVE_SLITHER`` : float
    Probability of selecting a slither megamove (all ``MOVE_*`` must sum to 1.0).
    Default 0.0.

``SLITHER_SUBSTEPS`` : int
    Number of slither moves applied to each non-frozen chain, in random order, per
    megamove; must be a positive integer (default 10). Under ``PARALLELIZE`` it
    sets the parallel pass's total budget rather than an exact per-chain count (see
    `Performance`_).

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
the chains by *length*: chains no longer than the smallest block interior (and
within the kernel's 512-bead per-chain buffer, which for the slither applies to
heteropolymers) go to the parallel kernel, and every longer chain goes to the
serial kernel. The split depends only on the chain lengths, the box, the
interaction range and the frozen set, so it is a run constant (it changes only if
a resized equilibration changes the box) and never depends on the configuration the
system happens to be in. If the box decomposes into a single block (for the
slither and pull an axis is split only once it is at least 24 sites long without
long-range beads, 40 with them) the whole megamove is serial.

When both sets have chains, each megamove runs both passes, and a fair coin decides
which goes first. Each pass holds the other's chains fixed, and the two do not
commute, so a fixed order would preserve the Boltzmann distribution without being
reversible - enough for the main loop, but not for a
:doc:`system-wide TSMMC <tsmmc>` excursion, whose acceptance assumes reversible
sub-moves. The random order makes the megamove reversible.

The two passes schedule their work differently. The serial pass slithers each of its
chains exactly ``SLITHER_SUBSTEPS`` times in a shuffled order, as above. The
parallel pass has the same total budget (``SLITHER_SUBSTEPS`` x the number of chains
in the parallel set), but it shares it between blocks in proportion to the number of
chains sitting wholly inside each block's interior that sweep - topped up
largest-remainder-first so the whole budget is spent - and each block then picks its
chains uniformly at random, with replacement. So under ``PARALLELIZE`` a given short
chain is slithered a random number of times per megamove rather than exactly
``SLITHER_SUBSTEPS``, and not at all in a sweep that leaves it straddling a halo.
The block grid's origin is shifted at random on every sweep, drawn over the whole
length of each split axis rather than one block length, so that a chain as long as
an interior can be placed inside one from any starting position even when the box
does not divide evenly into blocks. A sweep in which no chain of the parallel set
fits inside any interior makes no attempts, and that pass logs zero attempts rather
than its budget - so the code-6 counts per megamove vary under ``PARALLELIZE``.

The split has to be by length rather than by the chains' current extents, even
though a long chain is usually compact enough to fit. Choosing per megamove from
the current configuration would make the choice of kernel a function of the state,
and because the parallel kernel can never push a chain past a block interior while
the serial kernel can, that choice quietly biases the run towards compact chains.
That was the behaviour before 1.0.8; see :ref:`PARALLELIZE <advanced-parallel>` for
the measured size of the effect.

``PARALLELIZE`` leaves the equilibrium distribution alone, but it is a *different*
Markov chain: only the block interiors move in a given sweep and the kernel draws
from a different random stream, so a parallel run relaxes more slowly per step and
does not reproduce a serial run at the same ``SEED``. Compare equilibrium averages,
not energies at a fixed step.
