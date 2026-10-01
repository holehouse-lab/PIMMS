.. _overview:

========
Overview
========

This page explains how PIMMS represents a system and how a simulation works. By
the end you should understand the lattice model, the Monte Carlo move set, the
energy function, and how to set up and judge the convergence of a run. The
exhaustive per-keyword reference lives in :doc:`keywords`.

.. _overview-lattice:

The lattice and polymers
========================

PIMMS is a **lattice** model: space is discretised into a square (2D) or cubic
(3D) grid of sites, set by the ``DIMENSIONS`` keyword (e.g. ``DIMENSIONS : 70 70
70`` for a 70×70×70 cube). Every site is either empty (implicit **solvent**) or
occupied by exactly one **bead** - PIMMS enforces hard-sphere exclusion, so two
beads can never share a site.

A **polymer chain** is a connected string of beads. Consecutive beads along a
chain must be **lattice-adjacent** in the Chebyshev sense: each coordinate differs
by ``-1``, ``0`` or ``+1`` from the previous bead (so a bead has up to 26 valid
successor sites in 3D, 8 in 2D, including diagonals). This flexible connectivity
lets chains fold compactly on the lattice.

Chains are declared with the ``CHAIN`` keyword, which gives a count and a
one-letter sequence and may be repeated to build a **multi-component** mixture::

    CHAIN : 20 QQQQQQQQQQ      # 20 copies of a 10-bead poly-Q homopolymer
    CHAIN : 5  EEEEKKKKEEEE     # 5 copies of a hetero-polymer
    CHAIN : 100 A               # 100 single-bead "particles"

Every bead letter used in a chain must be defined in the parameter file
(:ref:`overview-energy`). A chain of length 1 is a single free particle; a chain
where all beads are identical is a *homopolymer*; otherwise it is a
*heteropolymer*. Sequences are upper-cased on read unless
``CASE_INSENSITIVE_CHAINS : False`` (which lets ``a`` and ``A`` be distinct bead
types).

The simulation box is **periodic** by default (a chain that leaves one face
re-enters the opposite face). Setting ``HARDWALL : True`` instead makes the box
edges hard walls - see :ref:`overview-setup`.

The box need not be cubic/square: unequal axes (e.g. ``DIMENSIONS : 20 20 60``) are
fully supported in 2D and 3D, with periodic or hardwall boundaries. Every axis must
be at least 7 lattice sites (the smallest box that supports the super-long-range
interaction shell). The only other restriction is that cluster-rotation moves
(``MOVE_CLUSTER_ROTATE``) cannot be combined with a non-cubic box under periodic
boundaries, because a 90° rigid rotation is only an energy-preserving symmetry of a
cube/square (or of any box under hardwall).

.. _overview-moves:

Monte Carlo and the move set
============================

PIMMS samples configurations with **Metropolis Monte Carlo (MC)**. Starting from
the current configuration it repeatedly proposes a random change (a *move*),
computes the resulting energy change ``ΔE``, and accepts the move with
probability

.. math::

   P_\text{accept} = \min\!\left(1,\; e^{-\Delta E / T}\right),

where ``T`` is the ``TEMPERATURE`` (PIMMS works in reduced units with Boltzmann's
constant absorbed into ``T``, so energies and temperatures share one arbitrary
scale; for the energy scales typical of PIMMS parameter files a ``TEMPERATURE``
between 10 and 200 is usually the useful range). Downhill and level moves
(``ΔE ≤ 0``) are always accepted; uphill moves are accepted with a
temperature-dependent probability. This plain Metropolis rule is what most moves
use; the few whose proposals are not symmetric add a correction (the pull and VMMC
moves multiply the Boltzmann factor by a proposal-probability ratio, and a TSMMC
excursion is accepted on the work accumulated along its temperature schedule -
each move's page gives its exact criterion). Each move is constructed to satisfy
**detailed balance**, which guarantees that - given enough steps - the simulation
samples the correct Boltzmann distribution of configurations. Hard-sphere overlaps
are rejected outright.

One outer-loop **step** (``N_STEPS`` counts these) is *not* a single move: the
core crankshaft move is a "megamove" that performs ``CRANKSHAFT_SUBSTEPS``
single-bead sub-moves in total (in a serial run, each on a bead drawn at random
from all the non-frozen beads), and several other moves are likewise batched. The
true number of accept/reject operations is therefore far larger than ``N_STEPS``
(it is reported in ``TOTAL_MOVES.dat``).

Which moves are attempted, and how often, is set by the ``MOVE_*`` keywords -
fractions that **must sum to 1.0** (a keyfile whose move fractions do not add up is
rejected at startup). Every move has its own page in :doc:`moves/index`, which
explains the algorithm, why it satisfies detailed balance, and how to configure it.
The available moves (with their internal move codes, the codes used for the columns
of ``MOVE_FREQS.dat`` and ``ACCEPTANCE.dat``) are:

**Local / single-chain moves**

* **Crankshaft** (``MOVE_CRANKSHAFT``, 1) - the workhorse. Fast, local bead
  rotations applied throughout the system in optimised Cython; most of your
  move budget should usually go here. ``CRANKSHAFT_SUBSTEPS`` tunes how much work
  each crankshaft step does.
* **Chain translate / rotate / pivot** (``MOVE_CHAIN_TRANSLATE`` 2,
  ``MOVE_CHAIN_ROTATE`` 3, ``MOVE_CHAIN_PIVOT`` 4) - rigid-body translation or
  rotation of a whole chain, or a pivot of the shorter arm of a chain about a
  randomly chosen bead.
* **Head pivot** (``MOVE_HEAD_PIVOT``, 5) - pivots a single terminus; rarely
  useful (keep at 0).
* **Slither / reptation** (``MOVE_SLITHER``, 6) - advances a chain forwards or
  backwards through the lattice "like a snake", which efficiently relaxes chain
  conformations. Like the crankshaft this is a megamove over the whole system:
  when it is selected *every* non-frozen chain is slithered ``SLITHER_SUBSTEPS``
  times, in random order.
* **Pull** (``MOVE_PULL``, 11) - cooperative reptation of a sub-segment: an
  interior bead is displaced and the following beads are "pulled" along to keep
  the chain connected, letting chains rearrange in **dense** systems where rigid
  moves clash. Also a whole-system megamove: every non-frozen chain of length ≥ 3
  is pulled ``PULL_SUBSTEPS`` times. The two termini are never moved by a pull, so
  pair it with crankshaft and/or slither.
* **Jump-and-relax** (``MOVE_JUMP_AND_RELAX``, 13) - relax a chain, relocate it,
  and relax again.

**Collective / many-chain moves**

* **Cluster translate / rotate** (``MOVE_CLUSTER_TRANSLATE`` 7,
  ``MOVE_CLUSTER_ROTATE`` 8) - rigid-body moves of a whole connected cluster of
  chains. Relatively expensive; keep small (0.01-0.05).
* **Virtual-Move Monte Carlo, VMMC** (``MOVE_VMMC``, 14) - recruits a cluster of
  chains by *interaction-energy gradients* and translates it rigidly, to escape
  the kinetic traps that single-chain moves hit in dense/condensed phases.
* **Temperature-switch MC, TSMMC** (``MOVE_CTSMMC`` 9, ``MOVE_MULTICHAIN_TSMMC``
  10, ``MOVE_SYSTEM_TSMMC`` 12) - take a chain, subset of chains, or the whole
  system on a temperature *excursion* (heated along a schedule up to
  ``TSMMC_JUMP_TEMP``, or the current temperature plus ``TSMMC_FIXED_OFFSET``, and
  cooled back) to hop over energy barriers. See :doc:`advanced/tsmmc`.

The collective and enhanced-sampling moves (TSMMC, pull, jump-and-relax, VMMC) are
powerful for assembly/condensate problems; of these only **VMMC** is still
**experimental** and gated behind ``EXPERIMENTAL_FEATURES : True`` (setting
``MOVE_VMMC``, ``VMMC_MAX_DISPLACEMENT`` or ``VMMC_MAX_CLUSTER`` away from its
default without that flag is an error). A robust default move set for most
problems is mostly crankshaft with a little translate/rotate/pivot and slither.

.. _overview-energy:

Energy: interactions, solvation and angles
==========================================

The total potential energy is a sum of **pairwise bead-bead interactions**,
**bead-solvent (solvation)** terms, and **backbone-angle** penalties. Pairwise
interactions act over **three nested length scales**, defined by how far apart two
beads sit on the lattice (Chebyshev distance):

* **Short range (SR)** - beads in contact (Chebyshev distance 1).
* **Long range (LR)** - Chebyshev distance 2.
* **Super-long range (SLR)** - Chebyshev distance 3.

SR is always present; LR and SLR are optional and only act between bead types that
declare them. This lets you model, e.g., a strong short-ranged "sticker"
attraction plus a weak longer-ranged electrostatic-like tail. Every pair of beads
in range counts, bonded neighbours included: two consecutive beads are always in
SR contact, so each bond contributes its SR energy as a constant that cancels in
every move but is part of the total.

**Solvation** enters through the SR shell only. Every empty site next to a bead
is a bead-solvent contact scored with that bead type's solvation energy, so a
lone bead in bulk carries 26 of them in 3D (8 in 2D), and each bead-bead contact
it forms replaces one. The LR and SLR shells have no solvent term.

All of this is specified in the **parameter file** (the ``PARAMETER_FILE``
keyword), which has a few kinds of line. All interaction energies and *absolute*
``ANGLE_PENALTY`` values must be **integers** (floats are rejected with an error);
the temperature-normalised ``ANGLE_PENALTY_T_NORM`` values are floats:

.. code-block:: text

   ## pairwise interactions:  R1 R2  e_SR [e_LR [e_SLR]]
   A  A   -8                 # A-A short-range contact energy
   A  B   -3  -2             # A-B short-range AND long-range energies
   B  B   -6  -3   3         # B-B short, long and super-long-range energies

   ## solvation (interaction with solvent, denoted 0):  R 0  e_solv
   A  0   -2                 # required for EVERY bead type
   B  0   -1

   ## backbone angle penalties:  ANGLE_PENALTY R  a1 a2 a3
   ANGLE_PENALTY  A   30 10 0
   ANGLE_PENALTY  B   50 20 0

   ## ...or, INSTEAD of the ANGLE_PENALTY line for that bead type, a
   ## temperature-normalised penalty (units of kT, k=1), multiplied by
   ## TEMPERATURE at parse time:  ANGLE_PENALTY_T_NORM R  a1 a2 a3
   ## ANGLE_PENALTY_T_NORM  A   0.5 0.2 0

Notes:

* A bead type takes **one** angle line: either ``ANGLE_PENALTY`` or
  ``ANGLE_PENALTY_T_NORM``, never both. The ``ANGLE_PENALTY_T_NORM`` line above is
  commented out for that reason - as written the file runs, and with that line
  uncommented as well PIMMS stops at start-up with ``Multiple ANGLE_PENALTY
  definitions for residue A``. To use it, delete the ``ANGLE_PENALTY  A`` line.

* ``ANGLE_PENALTY_T_NORM`` is scaled **once**, at parse time, by the temperature the
  production phase runs at: ``TEMPERATURE`` for a fixed-temperature run, and
  ``QUENCH_END`` for a ``QUENCH_RUN``. The Hamiltonian is not rebuilt as a quench
  ramps, so during the ramp the stiffness in units of the *current* kT drifts, and
  the penalties are exact in kT only once the quench reaches ``QUENCH_END``.

* Negative energies are **favourable** (attractive); positive are repulsive.
* The short-range interaction matrix must be **complete and non-redundant**: every
  pair of bead types the file defines (including each type with itself, and
  whether or not a chain uses them) needs exactly one short-range line. For types
  ``{A, B}`` that means ``A A``, ``A B`` and ``B B`` - a missing or duplicated pair
  is an error.
* A **solvation line for every bead type is mandatory** - solvent is the special
  type ``0`` and is part of the short-range matrix, because every empty
  neighbouring site is scored with the bead-solvent energy (see above). The
  solvent-solvent energy is fixed at 0. Long-range (LR/SLR) terms are for
  solute-solute pairs only; a solvent (``0``) entry in an LR/SLR line is an
  error. Unlike the short-range matrix, LR/SLR pairs need not be complete (any
  pair you omit defaults to 0).
* Angle penalties bias the local backbone geometry (three values per residue,
  keyed to displacement classes of the ``i-1`` to ``i+1`` displacement vector -
  each class mixes several geometric bend angles; see :doc:`input_files`). The
  penalty for a bend is taken from the type of the middle bead, so chains shorter
  than three beads have no angle term. Use either ``ANGLE_PENALTY`` (absolute
  integer penalties) or ``ANGLE_PENALTY_T_NORM`` (penalties in units of
  :math:`k_BT` with :math:`k_B=1`, scaled by the production temperature when the
  file is read, which keeps the stiffness fixed relative to temperature); unless
  angles are switched off, every bead type needs one of the two. Set
  ``ANGLES_OFF : True`` to disable angles entirely (then no angle lines are
  needed).
* Set ``NON_INTERACTING : True`` to zero all **pairwise** interaction and
  solvation energies, regardless of the parameter file. Angle penalties are
  unaffected - combine with ``ANGLES_OFF : True`` for a fully ideal
  excluded-volume-only reference simulation.

Internally the energy decomposes into SR, LR, SLR and angle components (see
``Hamiltonian.evaluate_total_energy``); the running total is what
``ENERGY.dat`` records, and ``ENERGY_CHECK`` periodically re-derives it from
scratch as a correctness guard.

.. _overview-setup:

Simulation setup: boundaries and box resizing
==============================================

**Boundary conditions.** With the default periodic boundaries a system behaves as
a bulk phase. ``HARDWALL : True`` gives hard walls instead: any proposal that
would put a bead outside the box or a bond across a face is rejected, and nothing
interacts through a face. This is appropriate for a droplet in a finite container,
or whenever you do not want periodic images. (A whole-chain translation still
draws its offset over the whole box, so it can relocate a chain from one wall to
the opposite one in a single move - see ``HARDWALL`` in :doc:`keywords`.)
Energetically a wall is **solvent**: a bead next to a wall receives its
bead-solvent energy for every out-of-box neighbour site, exactly as it would in
bulk, so walls are neither attractive nor repulsive - they only exclude volume.

**Box resizing for equilibration.** Sometimes you want to *condense* a system at
high effective concentration and then study it in a larger box. ``RESIZED_EQUILIBRATION``
runs the equilibration phase in a smaller box and then grows it to the full
``DIMENSIONS`` for production, placing the small box at the centre of the large one
(the configuration keeps its place inside it); ``EQUILIBRATION_OFFSET`` places the
small box elsewhere. The equilibration phase always runs with a
hardwall regardless of the keyfile; the production ``HARDWALL`` setting may be
True or False. With ``EQUILIBRATION : 0`` there is no phase to resize, so both
keywords are switched off with a warning.

**Centring.** For single-chain runs, ``AUTOCENTER : True`` writes the chain centred
in the middle of the box in every saved frame, so the trajectory needs no post-hoc
alignment. This is a write-time convention only - the sampling itself is unchanged.

(Box resizing also appears when *restarting* from a previous configuration into a
larger box - see :doc:`restart_files`.)

.. _overview-convergence:

Running and converging a simulation
===================================

A simulation is driven by a **keyfile**: a plain-text file of ``KEYWORD : value``
lines (``#`` starts a comment). A minimal but complete keyfile looks like this:

.. code-block:: text

   ## --- system ---
   DIMENSIONS      : 30 30 30
   PARAMETER_FILE  : params.prm
   CHAIN           : 50 AABBAABB     # 50 copies of an 8-bead heteropolymer
   TEMPERATURE     : 60

   ## --- run length ---
   N_STEPS         : 5000
   EQUILIBRATION   : 1000

   ## --- moves (must sum to 1.0) ---
   MOVE_CRANKSHAFT     : 0.8
   CRANKSHAFT_SUBSTEPS : 20000
   MOVE_CHAIN_TRANSLATE: 0.1
   MOVE_SLITHER        : 0.1
   SLITHER_SUBSTEPS    : 200

   ## --- output / analysis ---
   EN_FREQ         : 10
   XTC_FREQ        : 100
   ANA_CLUSTER     : 100

Run it with:

.. code-block:: bash

   PIMMS -k KEYFILE.kf

Every output file is written into the **current working directory**, so run each
simulation in its own directory.

PIMMS first runs ``EQUILIBRATION`` steps and then the remaining
``N_STEPS - EQUILIBRATION`` production steps, writing the requested output files
(see :doc:`output_files`). During equilibration no ``ANA_*`` analysis is performed
and no restart snapshot is written, but ``ENERGY.dat`` and ``PERFORMANCE.dat`` are
written throughout, and trajectory frames are saved unless you set
``SAVE_EQ : False``.

Runs are reproducible: on the same platform and PIMMS version, an identical
keyfile, parameter file and ``SEED`` reproduce the trajectory bit-for-bit. If
``SEED`` is not set a random seed is generated (announced at start-up and recorded
in ``log.txt`` and ``keyfile_used.kf``), so no two runs are alike.

**Judging convergence.** The first thing to check is ``ENERGY.dat``: the potential
energy should fall (or rise) and then **plateau** with stationary fluctuations -
that plateau marks equilibrium. If the energy is still drifting at the end of
equilibration, increase ``EQUILIBRATION`` (and probably ``N_STEPS``). For
assembly/condensate problems, also confirm that structural observables (e.g.
cluster sizes in ``CLUSTERS.dat``, radius of gyration in ``RG.dat``) have
stabilised, and that move **acceptance ratios** (``ACCEPTANCE.dat`` divided by
``MOVE_FREQS.dat``) are reasonable - extremely low acceptance for a move means it
is doing little useful work. ``PERFORMANCE.dat`` reports throughput and an
estimated time-to-completion. As a correctness safeguard, ``ENERGY_CHECK`` (on by
default, every 20000 steps; set it to 0 to disable) makes PIMMS periodically
re-compute the total energy from scratch and cross-check the lattice grids against
the chains, and abort (writing the offending configuration to
``CONFIG_AT_ENERGY_FAIL.pdb``/``.xtc``) if either has drifted.

.. _overview-next:

Beyond a basic run
==================

The keyfile above is the whole workflow for a standard fixed-temperature run.
Everything else PIMMS does is a keyword away, and each feature has its own page:

* **Restarting, resuming and growing a system.** ``RESTART_FREQ`` writes
  ``restart.pimms`` during production (and always at the end of a run);
  ``RESTART_FILE`` starts a new run from it, optionally with new chains added
  (``EXTRA_CHAIN``), and ``RESTART_CONTINUE`` instead resumes the run that wrote
  it exactly, step numbers, random draws and temperature included (so the keyfile
  must describe the same system, box and temperature or quench ramp). A hardwall
  snapshot can also be grown into a larger box by giving the larger
  ``DIMENSIONS``; a periodic snapshot must keep the box it was written with. See
  :doc:`restart_files`.
* **Quench / simulated annealing.** ``QUENCH_RUN`` ramps the temperature from
  ``QUENCH_START`` to ``QUENCH_END`` along a schedule instead of holding
  ``TEMPERATURE`` fixed. See :doc:`advanced/quench`.
* **Frozen chains.** ``FREEZE_FILE`` holds chosen chains rigidly in place as a
  scaffold, surface or template; they never move but still contribute to the
  energy. See :doc:`advanced/freeze`.
* **Parallelisation.** ``PARALLELIZE : True`` runs the crankshaft, slither and pull
  moves on multi-threaded checkerboard kernels (``PARALLEL_THREADS`` sets the thread
  count; chains too long to fit inside a block stay on the serial slither and pull
  kernels). Each sweep only the beads and chains inside the block interiors can
  move, so the attempts are shared among them: the crankshaft sub-moves are not
  spread over every bead, and a chain on the parallel kernel does not get exactly
  ``SLITHER_SUBSTEPS`` or ``PULL_SUBSTEPS`` attempts. The sampled equilibrium is
  the same as a serial run's, but the Markov chain is different and relaxes more
  slowly per step, so compare equilibrium averages rather than energies at a fixed
  step. See :doc:`advanced/parallelization`.
* **Enhanced sampling.** The TSMMC temperature-excursion moves (:doc:`advanced/tsmmc`)
  and the collective moves in :doc:`moves/index`.
* **Reference ensembles and controls.** ``NON_INTERACTING``, ``ANGLES_OFF``,
  ``RESIZED_EQUILIBRATION`` and friends. See :doc:`advanced/reference_controls`.
* **Analysis during the run.** The ``ANA_*`` keywords write the standard
  observables (:doc:`output_files`); ``ANALYSIS_MODULE`` plus ``ANA_CUSTOM`` run
  your own Python analysis against the live lattice (:doc:`advanced/custom_analysis`).
* **Analysis after the run.** The bundled ``lemonade`` package loads a finished
  trajectory and computes conformational, cluster and phase-separation properties.
  See :doc:`lemonade/index`.

The complete list of keywords, with types and defaults, is in :doc:`keywords` (the
same content the ``PIMMS --info`` command prints). The exact file formats are
described in :doc:`input_files` and :doc:`output_files`.
