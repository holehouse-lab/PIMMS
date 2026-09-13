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

Every ``MOVE_*`` keyword defaults to ``0.0``, so a keyfile must set the ones it
wants and their total must come to ``1.0`` (the parser checks the sum to within
``1e-7`` and rejects a negative fraction).

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
   :widths: 6 24 9 13 48

   * - Code
     - Keyword
     - Status
     - Scope
     - Move
   * - 1
     - ``MOVE_CRANKSHAFT``
     - core
     - whole system
     - Local single-bead perturbations (the workhorse)
   * - 2
     - ``MOVE_CHAIN_TRANSLATE``
     - core
     - one chain
     - Rigid translation of a whole chain
   * - 3
     - ``MOVE_CHAIN_ROTATE``
     - core
     - one chain
     - Rigid 90/180/270 degree rotation of a whole chain
   * - 4
     - ``MOVE_CHAIN_PIVOT``
     - core
     - one chain
     - Pivot the shorter arm of a chain about an interior bead
   * - 5
     - ``MOVE_HEAD_PIVOT``
     - core
     - one chain
     - Move a single terminus (rarely useful)
   * - 6
     - ``MOVE_SLITHER``
     - core
     - whole system
     - Reptation - the chain slides forwards/backwards
   * - 7
     - ``MOVE_CLUSTER_TRANSLATE``
     - core
     - one cluster
     - Rigid translation of a connected cluster of chains
   * - 8
     - ``MOVE_CLUSTER_ROTATE``
     - core
     - one cluster
     - Rigid rotation of a connected cluster of chains
   * - 9
     - ``MOVE_CTSMMC``
     - stable
     - one chain
     - Temperature-switch excursion of a single chain
   * - 10
     - ``MOVE_MULTICHAIN_TSMMC``
     - stable
     - chain subset
     - Temperature-switch excursion of a random subset of chains
   * - 11
     - ``MOVE_PULL``
     - stable
     - whole system
     - Cooperative reptation of a sub-segment (dense systems)
   * - 12
     - ``MOVE_SYSTEM_TSMMC``
     - stable
     - whole system
     - Temperature-switch excursion of the entire system
   * - 13
     - ``MOVE_JUMP_AND_RELAX``
     - stable
     - one chain
     - Relax, jump, relax a single chain
   * - 14
     - ``MOVE_VMMC``
     - experimental
     - recruited cluster
     - Virtual-Move MC collective cluster move

The three TSMMC variants (codes 9, 10 and 12) share one page, since they are the
same algorithm applied at different scopes. Only the **VMMC** move is still
experimental: setting ``MOVE_VMMC``, ``VMMC_MAX_DISPLACEMENT`` or
``VMMC_MAX_CLUSTER`` away from its default requires ``EXPERIMENTAL_FEATURES :
True``. The temperature-switch, pull and jump-and-relax moves are stable and need
no special flag.

.. _moves-step-anatomy:

What one step actually does
===========================

At the top of every step the main loop draws a **chain** uniformly at random from
the chains that are not frozen (see :ref:`FREEZE_FILE <advanced-freeze>`), and
then independently draws a **move** from the ``MOVE_*`` fractions. The two draws
are unrelated: the move selector never looks at the chain that was picked.

What the move then does with that chain depends on its scope:

* **Per-chain moves** (codes 2, 3, 4, 5, 9, 13) act on the drawn chain. The
  multi-chain excursion (code 10) is the exception: it draws its own random subset
  of the non-frozen chains and ignores the chain the loop picked.
* **Whole-system megamoves** (codes 1, 6 and 11) ignore the drawn chain entirely
  and act on the whole box, performing many accept/reject decisions internally on
  the same Markov chain. The crankshaft draws each of its sub-moves from all the
  non-frozen beads in the system; the slither and the pull sweep every eligible
  non-frozen chain (a pull needs three beads or more). Every one of those sub-moves
  is counted in ``MOVE_FREQS.dat`` / ``ACCEPTANCE.dat``, so for these three codes
  the counters advance in large jumps rather than by one per step.
* **The system-wide TSMMC excursion** (code 12) also acts on the whole box, but it
  logs a single attempt per excursion; the moves made inside the excursion belong
  to an auxiliary Markov chain and appear only in ``TOTAL_MOVES.dat``. The same
  holds for the chain and multi-chain excursions (codes 9 and 10).
* **Cluster moves** (codes 7, 8) act on the connected component containing the
  drawn chain; **VMMC** (code 14) seeds on the drawn chain and recruits outward.
  Both log one attempt per step.

Because the move draw is independent of the chain, a per-chain move can be drawn
for a chain it cannot act on - a rotation, pivot or head pivot of a single-bead
chain, for example. Such a draw is returned as drawn and the move rejects it as a
**null move**, which detailed balance permits.

Two start-up guards go with this:

* If **no enabled move can act on the system** (every chain frozen, or every
  enabled move needs longer chains than the box contains - dimers with only
  ``MOVE_CHAIN_PIVOT``, say) the run is refused with the reason rather than
  producing output for a configuration that never changes.
* If the system contains **monomers** and ``MOVE_CHAIN_ROTATE``,
  ``MOVE_CHAIN_PIVOT`` or ``MOVE_HEAD_PIVOT`` is enabled, the start-up summary
  warns and estimates what fraction of steps will be null moves.

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

Detailed balance is not the whole story
---------------------------------------

Equation :eq:`db` guarantees that :math:`\pi` is *invariant* under the move. It
does not guarantee that the run converges to :math:`\pi`, because it says nothing
about whether the moves can reach the whole configuration space. That second
condition is **irreducibility**, and it is a property of the move *set* you chose
in the keyfile, not of any individual move.

The failure mode to watch for is a move set in which some feature of a chain can
never change. Every rigid move - :doc:`chain_translate`, :doc:`chain_rotate`,
:doc:`cluster_translate`, :doc:`cluster_rotate` and :doc:`vmmc` - relocates or
reorients a chain without ever altering its shape, so a move set built only from
those samples the Boltzmann distribution *restricted to the conformations the
system started in*. The chains will move convincingly around the box while
``RG.dat`` and ``END_TO_END_DIST.dat`` repeat the same numbers on every row, and
the scaling exponent gets fitted to constant data. There are narrower versions of
the same trap: :doc:`pull` moves interior beads but never the termini, so pull
plus a rigid move holds every end-to-end distance fixed; head-pivot alone never
touches an interior bead; chain-pivot alone never moves the beads nearest the
middle of the chain.

None of these is a detailed-balance violation, and none of them will show up in an
``ENERGY_CHECK`` - the energy is tracked perfectly, it is simply the energy of a
restricted ensemble, so the average it converges to depends on the starting
configuration and can be off by any amount in either direction. PIMMS warns at
start-up when it detects one of these move sets, and names the output files whose
numbers should not be trusted. It warns rather than refuses, because a rigid-only
move set is a legitimate way to study rigid-body assembly of chains whose
conformation you deliberately want fixed.

The general rule: to sample conformations you need at least one move that changes
conformation. :doc:`crankshaft` is the usual choice, and :doc:`slither`,
:doc:`pull`, :doc:`chain_pivot`, :doc:`head_pivot` and :doc:`jump_and_relax` also
qualify.
