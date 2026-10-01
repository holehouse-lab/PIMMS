.. _move-pull:

============================
Pull (cooperative reptation)
============================

:Keyword: ``MOVE_PULL``
:Move code: 11
:Status: stable
:Scope: whole system (megamove)

How it works
============

The pull move is a **cooperative, local reptation of a sub-segment** of a chain.
An interior bead - index :math:`i` drawn uniformly from :math:`1` to :math:`L-2` -
is displaced to a nearby empty site, and the beads on one side of it are "pulled"
along one after another, each taking the site just vacated by its predecessor,
until connectivity is restored. A direction is drawn with probability :math:`1/2`:
C-ward keeps bead :math:`i-1` as the anchor and cascades through beads
:math:`i+1, i+2, \dots`; N-ward keeps bead :math:`i+1` as the anchor and cascades
the other way.

The perturbation propagates bead by bead along the chain rather than moving it
rigidly, so it needs almost no free volume: it can rearrange chains in
**dense/condensed systems** where rigid translate/rotate moves simply clash. Chain
identity and bonding are preserved throughout, and the move works in 2D and 3D.

Two things end a cascade. If the next bead along is already Chebyshev-1 adjacent
to the site its predecessor just moved into, connectivity is restored and the
cascade stops there. If instead the cascade runs all the way to a terminus, the
move is rejected - there is no valid anchor for the reverse move. An accepted pull
therefore **never displaces either chain terminus**, so it is best paired with the
crankshaft or slither, which do move the ends.

A pull requires chains of length :math:`\ge 3` (an interior bead needs neighbours
on both sides). Like the slither, a pull *step* is a megamove: every non-frozen
chain of length :math:`\ge 3` is pulled ``PULL_SUBSTEPS`` times, in one globally
shuffled order, each substep with its own accept/reject inside the Cython kernel
(under ``PARALLELIZE`` the short chains are scheduled per block instead - see
`Performance`_). The chain the main loop drew for the step plays no part. Frozen
chains and chains of one or two beads are never selected but remain as fixed
obstacles. Every substep is logged under code 11 in ``MOVE_FREQS.dat`` /
``ACCEPTANCE.dat``, including the ones that did nothing because the first-target
set was empty or the cascade reached a terminus. A pull megamove in a system with
no eligible chain does nothing and logs no attempts; if ``MOVE_PULL`` is then the
only enabled move (no chain of three or more beads, or every such chain frozen),
the run is refused at start-up instead.

Why detailed balance holds
==========================

The cascade starts by moving the chosen bead to a site drawn uniformly from a
**first-target set** :math:`S`: the empty sites that are Chebyshev-1 neighbours of
*both* the moving bead's current site and the anchor, so the moved bead stays
bonded to the anchor and the next pulled bead can occupy the site it vacates.
(Under ``HARDWALL`` adjacency through the periodic wrap does not count, so no
target on the far side of a wall is ever offered.) If that set is empty the
substep does nothing. Once the first step is fixed, the rest of the cascade is
determined.

The size of this set differs between the forward move and its reverse,
:math:`n_\text{fwd}=|S_\text{fwd}|` and :math:`n_\text{rev}=|S_\text{rev}|`, so
the proposal is asymmetric and the move uses the Metropolis-Hastings acceptance

.. math::

   A(x\to y) = \min\!\left(1,\; \frac{n_\text{fwd}}{n_\text{rev}}\;
   e^{-\Delta E / T}\right).

(The Hastings factor is :math:`g(y\to x)/g(x\to y) = n_\text{fwd}/n_\text{rev}`,
since the forward proposal picks its first target with probability
:math:`1/n_\text{fwd}` and the reverse with :math:`1/n_\text{rev}`; the uniform
draws of the bead index and the direction cancel.) A reverse count of zero makes
the reverse move impossible and rejects the proposal outright.

The cascade is constructed so that the reverse move retraces exactly the forward
path: the reverse of a C-ward pull whose cascade stopped at bead :math:`k` is an
N-ward pull started at bead :math:`k`, and the forward cascade's "did not stop"
conditions force the reverse cascade to stop in the same place. Both counts come
from the *same* predicate, so they can never diverge. The energy change uses the
same telescoping per-bead decomposition as the slither (pair terms plus the angle
penalty, summed over each bead move in the cascade), and the whole construction is
verified by the detailed-balance test suite.

Configuration
=============

``MOVE_PULL`` : float
    Probability of selecting a pull megamove (all ``MOVE_*`` must sum to 1.0).
    Default 0.0.

``PULL_SUBSTEPS`` : int
    Number of pull moves applied to each eligible (non-frozen, length
    :math:`\ge 3`) chain per megamove; must be a positive integer (default 10).
    Under ``PARALLELIZE`` it sets the parallel pass's total budget rather than an
    exact per-chain count (see `Performance`_).

Performance
===========

Pull is aimed squarely at dense-phase rearrangement; for moving whole correlated
groups of chains see :doc:`vmmc`.

Along with the crankshaft and slither, pull has a multi-threaded kernel: see
:ref:`PARALLELIZE <advanced-parallel>`. As a whole-chain move it uses the same
chain-level block decomposition as the slither (a chain parallelizes only if all
its beads sit at least :math:`W = R_\text{int} + 2` sites inside one block); in
addition its first-target search is restricted to the block interior so the
Metropolis-Hastings multiplicity ratio stays self-consistent. As for the slither,
the eligible chains are split by length: those no longer than the smallest block
interior, and within the kernel's 512-bead per-chain buffer (which for pull applies
to *all* chains, not just heteropolymers), are pulled by the parallel kernel, and
every longer chain is pulled by the serial kernel. The split is a run constant and
never depends on the current configuration; a per-megamove choice made from the
chains' current extents biases the sampled conformations towards compact chains,
which is what PIMMS did before 1.0.8. If the box decomposes into a single block the
whole megamove is serial.

Everything else about the two passes is as described for the :doc:`slither`: when
both sets have chains a fair coin decides which pass runs first (which keeps the
megamove reversible for system-wide TSMMC excursions); the parallel pass spends a
total budget of ``PULL_SUBSTEPS`` x (chains in the parallel set), shared between
blocks by how many of those chains sit wholly inside each block's interior that
sweep, with each block picking its chains uniformly at random, so a given short
chain is not pulled exactly ``PULL_SUBSTEPS`` times per megamove; the block origin
is shifted at random over the whole box on every sweep; and a sweep with no movable
chain logs zero attempts rather than its budget.

As for the slither, ``PARALLELIZE`` leaves the equilibrium distribution alone but
follows a different Markov chain from a serial run - only the block interiors move
in a given sweep, from a different random stream - so relaxation per step is slower
and a same-``SEED`` parallel run does not track the serial one. Compare equilibrium
averages, not energies at a fixed step.
