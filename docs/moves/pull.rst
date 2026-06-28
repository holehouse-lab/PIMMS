.. _move-pull:

============================
Pull (cooperative reptation)
============================

:Keyword: ``MOVE_PULL``
:Move code: 11
:Status: experimental (requires ``EXPERIMENTAL_FEATURES : True``)

How it works
============

The pull move is a **cooperative, local reptation of a sub-segment** of a chain.
An interior bead is displaced to a nearby empty site, and the beads following it
are "pulled" along one after another - each taking the site just vacated by its
predecessor - until connectivity is restored. The perturbation propagates
bead-by-bead along the chain rather than moving it rigidly, so it needs almost no
free volume: it can rearrange chains in **dense/condensed systems** where rigid
translate/rotate moves simply clash. Chain identity and bonding are preserved
throughout. Requires chains of length :math:`\ge 3`.

Like the slither, a pull *step* is a megamove: every eligible chain is pulled
``PULL_SUBSTEPS`` times, each with its own accept/reject. The two chain termini are
not displaced by a pull, so it is best paired with the crankshaft/slither, which
do move the ends.

Why detailed balance holds
==========================

The cascade starts by moving the chosen bead to a site drawn uniformly from a
**first-target set** :math:`S` of valid empty sites; once that first step is
fixed, the rest of the cascade is determined. The size of this set differs between
the forward move and its reverse, :math:`n_\text{fwd}=|S_\text{fwd}|` and
:math:`n_\text{rev}=|S_\text{rev}|`, so as with reptation the proposal is
asymmetric and the move uses the Metropolis-Hastings acceptance

.. math::

   A(x\to y) = \min\!\left(1,\; \frac{n_\text{rev}}{n_\text{fwd}}\;
   e^{-\Delta E / T}\right).

The cascade is constructed so that the reverse move retraces exactly the forward
path (the "didn't-stop" conditions of the forward cascade force the reverse to
stop at the same place), which is what guarantees reversibility and makes the
multiplicity ratio the correct Hastings factor. The energy change uses the same
telescoping per-bead decomposition as the slither. The construction is verified by
the detailed-balance test suite.

Configuration
=============

``MOVE_PULL`` : float
    Probability of selecting a pull megamove (all ``MOVE_*`` must sum to 1.0).

``PULL_SUBSTEPS`` : int
    Number of pull moves applied to each (length :math:`\ge 3`) chain per megamove
    (default 10).

Pull is an experimental move aimed squarely at dense-phase rearrangement; for
moving whole correlated groups of chains see :doc:`vmmc`.

Along with the crankshaft and slither, pull has a multi-threaded kernel: see
:ref:`PARALLELIZE <advanced-parallel>`. As a whole-chain move it uses the same
chain-level block decomposition as the slither (a chain parallelizes only if all
its beads fit in a block interior); in addition its cooperative-reptation target
search is restricted to the block interior so the Metropolis-Hastings multiplicity
ratio stays self-consistent.
