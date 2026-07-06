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
edges hard, reflecting walls - see :ref:`overview-setup`.

The box need not be cubic/square: unequal axes (e.g. ``DIMENSIONS : 20 20 60``) are
fully supported in 2D and 3D, with periodic or hardwall boundaries. The only
restriction is that cluster-rotation moves (``MOVE_CLUSTER_ROTATE``) cannot be
combined with a non-cubic box under periodic boundaries, because a 90° rigid
rotation is only an energy-preserving symmetry of a cube/square (or of any box under
hardwall).

.. _overview-moves:

Monte Carlo and the move set
============================

PIMMS samples configurations with **Metropolis Monte Carlo (MC)**. Starting from
the current configuration it repeatedly proposes a random change (a *move*),
computes the resulting energy change ``ΔE``, and accepts the move with
probability

.. math::

   P_\text{accept} = \min\!\left(1,\; e^{-\Delta E / T}\right),

where ``T`` is the ``TEMPERATURE``. Downhill moves (``ΔE < 0``) are always
accepted; uphill moves are accepted with a temperature-dependent probability.
Each move is constructed to satisfy **detailed balance**, which guarantees that -
given enough steps - the simulation samples the correct Boltzmann distribution of
configurations. Hard-sphere overlaps are rejected outright.

One outer-loop **step** (``N_STEPS`` counts these) is *not* a single move: the
core crankshaft move is a "megamove" that performs ``CRANKSHAFT_SUBSTEPS``
single-bead sub-moves in total (spread at random over all beads), and several other
moves are likewise batched. The true number of accept/reject operations is
therefore far larger than ``N_STEPS`` (it is reported in ``TOTAL_MOVES.dat``).

Which moves are attempted, and how often, is set by the ``MOVE_*`` keywords -
fractions that **must sum to 1.0**. The available moves (with their internal move
codes) are:

**Local / single-chain moves**

* **Crankshaft** (``MOVE_CRANKSHAFT``, 1) - the workhorse. Fast, local bead
  rotations applied throughout the system in optimised Cython; most of your
  move budget should usually go here. ``CRANKSHAFT_SUBSTEPS`` tunes how much work
  each crankshaft step does.
* **Chain translate / rotate / pivot** (``MOVE_CHAIN_TRANSLATE`` 2,
  ``MOVE_CHAIN_ROTATE`` 3, ``MOVE_CHAIN_PIVOT`` 4) - rigid-body translation or
  rotation of a whole chain, or a pivot of one half of a chain about a randomly
  chosen point.
* **Head pivot** (``MOVE_HEAD_PIVOT``, 5) - pivots a single terminus; rarely
  useful (keep at 0).
* **Slither / reptation** (``MOVE_SLITHER``, 6) - advances a chain forwards or
  backwards through the lattice "like a snake". Efficiently relaxes chain
  conformations; ``SLITHER_SUBSTEPS`` sets how many per chain.
* **Pull** (``MOVE_PULL``, 11) - cooperative reptation of a sub-segment: an
  interior bead is displaced and the following beads are "pulled" along to keep
  the chain connected, letting chains rearrange in **dense** systems where rigid
  moves clash. Requires chains of length ≥ 3.
* **Jump-and-relax** (``MOVE_JUMP_AND_RELAX``, 13) - relax a chain, relocate it,
  and relax again.

**Collective / many-chain moves**

* **Cluster translate / rotate** (``MOVE_CLUSTER_TRANSLATE`` 7,
  ``MOVE_CLUSTER_ROTATE`` 8) - rigid-body moves of a whole connected cluster of
  chains. Relatively expensive; keep small (0.01-0.05).
* **Virtual-Move Monte Carlo, VMMC** (``MOVE_VMMC``, 14) - recruits a cluster of
  chains by *interaction-energy gradients* and moves it rigidly, to escape the
  kinetic traps that single-chain moves hit in dense/condensed phases.
* **Temperature-switch MC, TSMMC** (``MOVE_CTSMMC`` 9, ``MOVE_MULTICHAIN_TSMMC``
  10, ``MOVE_SYSTEM_TSMMC`` 12) - take a chain, subset of chains, or the whole
  system on a temperature *excursion* to hop over energy barriers. See
  :doc:`advanced/index`.

The collective and enhanced-sampling moves (TSMMC, pull, jump-and-relax, VMMC) are
powerful for assembly/condensate problems; of these only **VMMC** is still
**experimental** and gated behind ``EXPERIMENTAL_FEATURES : True``. A robust default
move set for most problems is mostly crankshaft with a little translate/rotate/pivot
and slither.

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
attraction plus a weak longer-ranged electrostatic-like tail.

All of this is specified in the **parameter file** (the ``PARAMETER_FILE``
keyword). All interaction and angle energies must be **integers** (floats are
rejected with an error). It has a few kinds of line:

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

   ## ...or temperature-normalised angle penalties (units of kT, k=1):
   ##   ANGLE_PENALTY_T_NORM R  a1 a2 a3
   ANGLE_PENALTY_T_NORM  A   0.5 0.2 0    # multiplied by TEMPERATURE at parse time

Notes:

* Negative energies are **favourable** (attractive); positive are repulsive.
* The short-range interaction matrix must be **complete and non-redundant**: every
  pair of bead types you use (including each type with itself) needs exactly one
  short-range line. For types ``{A, B}`` that means ``A A``, ``A B`` and ``B B`` -
  a missing or duplicated pair is an error.
* A **solvation line for every bead type is mandatory** - the energy is measured
  relative to a fully solvated reference, so PIMMS needs to know each bead's
  bead-solvent energy. Solvent is the special type ``0``; the solvent-solvent
  energy is fixed at 0. Long-range (LR/SLR) terms are for solute-solute pairs only
  - a solvent (``0``) entry in an LR/SLR line is an error. Unlike the short-range
  matrix, LR/SLR pairs need not be complete (any pair you omit defaults to 0).
* Angle penalties bias the local backbone geometry (three values per residue for
  the distinct lattice bend angles). Use either ``ANGLE_PENALTY`` (absolute integer
  penalties) or ``ANGLE_PENALTY_T_NORM`` (penalties in units of :math:`k_BT` with
  :math:`k_B=1`, scaled by ``TEMPERATURE`` when the file is read - handy for keeping
  the stiffness fixed relative to temperature). Set ``ANGLES_OFF : True`` to disable
  angles entirely (then no angle lines are needed).
* Set ``NON_INTERACTING : True`` to zero **all** interaction energies and run a
  pure excluded-volume reference simulation, regardless of the parameter file.

Internally the energy decomposes into SR, LR, SLR and angle components (see
``Hamiltonian.evaluate_total_energy``); the running total is what
``ENERGY.dat`` records, and ``ENERGY_CHECK`` periodically re-derives it from
scratch as a correctness guard.

.. _overview-setup:

Simulation setup: boundaries and box resizing
==============================================

**Boundary conditions.** With the default periodic boundaries a system behaves as
a bulk phase. ``HARDWALL : True`` gives reflective walls instead - appropriate for
a droplet in a finite container, or whenever you do not want chains wrapping
across the box.

**Box resizing for equilibration.** Sometimes you want to *condense* a system at
high effective concentration and then study it in a larger box. ``RESIZED_EQUILIBRATION``
runs the equilibration phase in a smaller box and then grows it (re-centring the
chains) to the full ``DIMENSIONS`` for production; ``EQUILIBRATION_OFFSET`` places
the small box within the large one. Resizing requires ``HARDWALL : True``.

**Centring.** For single-chain runs, ``AUTOCENTER : True`` keeps the chain in the
middle of the box every frame, so the saved trajectory needs no post-hoc
alignment.

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

PIMMS first runs ``EQUILIBRATION`` steps (no analysis output, and trajectory
frames only if ``SAVE_EQ : True``) and then the remaining production steps,
writing the requested output files (see :doc:`output_files`).

**Judging convergence.** The first thing to check is ``ENERGY.dat``: the potential
energy should fall (or rise) and then **plateau** with stationary fluctuations -
that plateau marks equilibrium. If the energy is still drifting at the end of
equilibration, increase ``EQUILIBRATION`` (and probably ``N_STEPS``). For
assembly/condensate problems, also confirm that structural observables (e.g.
cluster sizes in ``CLUSTERS.dat``, radius of gyration in ``RG.dat``) have
stabilised, and that move **acceptance ratios** (``ACCEPTANCE.dat`` divided by
``MOVE_FREQS.dat``) are reasonable - extremely low acceptance for a move means it
is doing little useful work. ``PERFORMANCE.dat`` reports throughput and an
estimated time-to-completion. As a correctness safeguard you can set
``ENERGY_CHECK`` so PIMMS periodically re-computes the total energy from scratch
and aborts if it has drifted.
