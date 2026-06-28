.. _moves:

=====
Moves
=====

PIMMS samples configurations with **Metropolis Monte Carlo (MC)**. At each step it
picks a move at random (according to the ``MOVE_*`` fractions, which must sum to
``1.0``), proposes the corresponding change, and accepts or rejects it so that the
simulation converges on the correct Boltzmann distribution. This section has one
page per move: each describes *how the move works*, *why it preserves detailed
balance*, and the *configuration* relevant to it.

.. toctree::
   :maxdepth: 1

   crankshaft
   chain_translate
   chain_rotate
   chain_pivot
   head_pivot
   slither
   pull
   cluster_translate
   cluster_rotate
   tsmmc
   jump_and_relax
   vmmc

The move set at a glance
========================

.. list-table::
   :header-rows: 1
   :widths: 8 26 10 56

   * - Code
     - Keyword
     - Status
     - Move
   * - 1
     - ``MOVE_CRANKSHAFT``
     - core
     - Local single-bead perturbations (the workhorse)
   * - 2
     - ``MOVE_CHAIN_TRANSLATE``
     - core
     - Rigid translation of a whole chain
   * - 3
     - ``MOVE_CHAIN_ROTATE``
     - core
     - Rigid 90/180/270° rotation of a whole chain
   * - 4
     - ``MOVE_CHAIN_PIVOT``
     - core
     - Pivot one half of a chain about an interior point
   * - 5
     - ``MOVE_HEAD_PIVOT``
     - core
     - Pivot a single terminus (rarely useful)
   * - 6
     - ``MOVE_SLITHER``
     - core
     - Reptation - the chain slides forwards/backwards
   * - 7
     - ``MOVE_CLUSTER_TRANSLATE``
     - core
     - Rigid translation of a connected cluster of chains
   * - 8
     - ``MOVE_CLUSTER_ROTATE``
     - core
     - Rigid rotation of a connected cluster of chains
   * - 9, 10, 12
     - ``MOVE_CTSMMC`` / ``MOVE_MULTICHAIN_TSMMC`` / ``MOVE_SYSTEM_TSMMC``
     - experimental
     - Temperature-switch (tempered-transitions) excursions
   * - 11
     - ``MOVE_PULL``
     - experimental
     - Cooperative reptation of a sub-segment (dense systems)
   * - 13
     - ``MOVE_JUMP_AND_RELAX``
     - experimental
     - Relax → jump → relax a single chain
   * - 14
     - ``MOVE_VMMC``
     - experimental
     - Virtual-Move MC collective cluster move

(The three TSMMC variants share one page, since they are the same algorithm
applied at different scopes.) Experimental moves require
``EXPERIMENTAL_FEATURES : True``.

.. _moves-db-primer:

A detailed-balance primer
=========================

Every move is constructed so that the simulation samples the **Boltzmann
distribution**

.. math::

   \pi(x) \;=\; \frac{1}{Z}\, e^{-E(x)/T},

where :math:`x` is a configuration, :math:`E(x)` its energy, :math:`T` the
``TEMPERATURE`` and :math:`Z` a normalising constant. (PIMMS works in reduced
units where Boltzmann's constant is absorbed into :math:`T`.)

A sufficient condition for converging on :math:`\pi` is **detailed balance**: for
every pair of configurations :math:`x, y`, the move's transition probability
:math:`P` must satisfy

.. math::
   :label: db

   \pi(x)\, P(x \to y) \;=\; \pi(y)\, P(y \to x).

A move is built from a **proposal** :math:`g(x\to y)` (the random change it
suggests) and an **acceptance** probability :math:`A(x\to y)`, so
:math:`P(x\to y) = g(x\to y)\,A(x\to y)` for :math:`y\neq x`. The
**Metropolis-Hastings** acceptance

.. math::
   :label: mh

   A(x\to y) \;=\; \min\!\left(1,\; \frac{g(y\to x)}{g(x\to y)}\,
   e^{-(E(y)-E(x))/T}\right)

satisfies :eq:`db` for *any* proposal. The ratio :math:`g(y\to x)/g(x\to y)`
corrects for an asymmetric proposal.

Two important special cases recur throughout this section:

* **Symmetric proposal.** If forward and reverse proposals are equally likely,
  :math:`g(x\to y) = g(y\to x)`, the ratio is 1 and :eq:`mh` reduces to the plain
  Metropolis criterion :math:`A = \min(1, e^{-\Delta E/T})` with
  :math:`\Delta E = E(y)-E(x)`. The local and rigid-body moves (crankshaft,
  translate, rotate, pivot, cluster moves) all use this.

* **Composition of valid moves.** A sequence of updates that each *individually*
  leave :math:`\pi` invariant also leaves :math:`\pi` invariant. Several PIMMS
  moves are "megamoves" - many sub-moves bundled into one step - or composites
  (e.g. :doc:`jump_and_relax`) that rely on this fact.

Where a move needs more than these (a genuine Hastings ratio, a work-accumulation
factor, or a link-probability product), the relevant page derives it explicitly.

.. note::

   Hard-sphere overlaps are rejected outright: PIMMS never allows two beads on the
   same site, which is equivalent to assigning such configurations infinite
   energy (:math:`\pi = 0`), consistent with :eq:`db`.
