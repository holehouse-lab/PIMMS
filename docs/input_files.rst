.. _input-files:

===========
Input files
===========

A PIMMS run is driven entirely by plain-text input files. There are up to three:

* a **keyfile** - the simulation configuration (box, chains, temperature, moves,
  output frequencies);
* a **parameter file** (``.prm``) - the force field (interaction energies and
  backbone-angle penalties);
* an optional **freeze file** - a list of chains to hold rigidly fixed.

Only the keyfile and the parameter file are required; the freeze file is read only
when you set the ``FREEZE_FILE`` keyword. All three share the same comment
convention: any line whose first non-whitespace character is ``#`` is ignored, as
are blank lines - so ``#`` (or ``##``) can be used freely for headers and
annotations. Inline comments (text after a ``#`` partway along a line) are stripped
too.

.. _input-keyfiles:

Keyfiles
========

The keyfile is the top-level description of a simulation. Each non-comment line
sets one keyword::

    KEYWORD : value

Whitespace around the colon is optional. A few rules govern the file as a whole:

* **Most keywords may appear at most once.** A repeated keyword is an error (PIMMS
  reports the offending keyword rather than silently taking the last value). The two
  exceptions are ``CHAIN`` and ``EXTRA_CHAIN``, which may be repeated to build a
  multi-component system.
* **Required keywords.** ``DIMENSIONS``, ``PARAMETER_FILE``, ``TEMPERATURE``,
  ``N_STEPS`` and ``EQUILIBRATION`` must always be present, plus ``CHAIN`` - unless a
  ``RESTART_FILE`` is provided, in which case the chains come from the restart file
  and ``CHAIN`` may be omitted. Every other keyword falls back to a default.
* **Values are checked on read.** A keyword expecting an integer, float or boolean
  that is given a malformed value fails immediately with a descriptive error, rather
  than deep inside the run.

A few keywords that shape the file are worth calling out here (the full list, with
types and defaults, is the :doc:`keyword reference <keywords>`):

* ``DIMENSIONS`` takes 2 values (a 2D simulation) or 3 (3D). Each value must be at
  least **8** lattice units - the smallest box that can support the super-long-range
  interaction shell - and the axes need not be equal (non-cubic/non-square boxes are
  fully supported).
* ``CHAIN : N SEQUENCE`` declares ``N`` copies of a chain whose one-letter
  ``SEQUENCE`` names the bead types (e.g. ``CHAIN : 20 QQQQQQQQQQ`` for 20 ten-bead
  poly-Q chains). Repeat the keyword for a mixture. Sequences are upper-cased on read
  unless ``CASE_INSENSITIVE_CHAINS : False``; every bead letter used must be defined
  in the parameter file.
* ``PARAMETER_FILE`` points at the ``.prm`` file described :ref:`below
  <input-parameter-files>`; ``FREEZE_FILE`` (optional) points at a
  :ref:`freeze file <input-freeze-files>`.

A minimal but complete keyfile:

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

Every keyword can be inspected from the command line without opening the docs::

    PIMMS --info                 # list every keyword, grouped
    PIMMS --info <KEYWORD>       # full details on one keyword
    PIMMS --info ALL             # print every keyword description

See the :doc:`keyword reference <keywords>` for the exhaustive, auto-generated list.

.. _input-parameter-files:

Parameter files (``.prm``)
==========================

The parameter file, named by the ``PARAMETER_FILE`` keyword, defines the **force
field**: every pairwise interaction energy and every backbone-angle penalty. Bead
types are the one-letter codes used in your ``CHAIN`` sequences; **solvent** is the
special type ``0`` (an empty lattice site).

All interaction and angle energies must be **integers** - a float is rejected with a
clear error. The file has a few kinds of line.

**1. Pairwise interactions** act over three nested length scales, set by the
Chebyshev distance between two beads:

.. code-block:: text

   ## R1 R2  e_SR [e_LR [e_SLR]]
   A  A   -8                 # A-A short-range contact energy (distance 1)
   A  B   -3  -2             # A-B short-range AND long-range (distance 2)
   B  B   -6  -3   3         # B-B short, long AND super-long-range (distance 3)

Short-range (SR) is always present; long-range (LR, distance 2) and super-long-range
(SLR, distance 3) are optional trailing columns.

**2. Solvation** is the bead-solvent energy, written as an interaction with type
``0``:

.. code-block:: text

   ## R 0  e_solv
   A  0   -2
   B  0   -1

**3. Backbone-angle penalties** bias the local chain geometry (three values per
residue, one for each distinct lattice bend angle). Use either absolute
integer penalties or temperature-normalised ones:

.. code-block:: text

   ## absolute integer penalties:
   ANGLE_PENALTY         A   30 10 0

   ## ...or temperature-normalised (units of kT, k=1; multiplied by TEMPERATURE):
   ANGLE_PENALTY_T_NORM  A   0.5 0.2 0

The rules PIMMS enforces:

* **Negative energies are favourable** (attractive); positive energies are repulsive.
* **The short-range matrix must be complete and non-redundant.** Every pair of bead
  types you use needs exactly one short-range line - each type with itself *and*
  every distinct pair. For types ``{A, B}`` that means ``A A``, ``A B`` and ``B B``;
  a missing or duplicated pair is an error.
* **A solvation line for every bead type is mandatory.** Energies are measured
  relative to a fully solvated reference, so PIMMS needs each bead's bead-solvent
  energy. The solvent-solvent energy is fixed at 0.
* **Long-range terms are solute-solute only.** A solvent (``0``) entry in an LR/SLR
  line is an error. Unlike the short-range matrix, LR/SLR pairs need **not** be
  complete: any pair you omit defaults to 0.
* **Angles are optional as a whole.** Use ``ANGLE_PENALTY`` or
  ``ANGLE_PENALTY_T_NORM`` (the T-normalised form keeps stiffness fixed relative to
  temperature). Setting ``ANGLES_OFF : True`` in the keyfile disables angles
  entirely, and no angle lines are then needed.
* ``NON_INTERACTING : True`` in the keyfile zeroes **all** interaction energies for a
  pure excluded-volume reference run, regardless of what the parameter file says.

A complete two-type parameter file:

.. code-block:: text

   ## interactions
   A  A   -8
   A  B   -3  -2
   B  B   -6  -3   3

   ## solvation (required for every bead type)
   A  0   -2
   B  0   -1

   ## angles
   ANGLE_PENALTY  A   30 10 0
   ANGLE_PENALTY  B   50 20 0

The exact parameters PIMMS parsed are echoed to ``parameters_used.prm`` at startup,
so every run is self-documenting. The conceptual role of the energy terms is
discussed in the :ref:`energy model <overview-energy>`.

.. _input-freeze-files:

Freeze files
============

A freeze file holds chosen chains **rigidly fixed** for the whole simulation - a
scaffold, a wall, or a pre-formed template the rest of the system explores around.
Point the ``FREEZE_FILE`` keyword at it; the run aborts at startup if the file does
not exist.

The file is a short list of ``C`` directives, one or more per file:

.. code-block:: text

   # freeze.txt  -  lines beginning with # are comments
   C 1 2 3          # freeze chainIDs 1, 2 and 3 (chainIDs are numbered from 1)
   C 10 11 12 13    # more C lines are allowed; IDs may be split across lines

Each ``C`` line contributes its integer chainIDs to the frozen set; order and the
split across lines do not matter. A frozen chain is **excluded from the pool of
chains PIMMS can move** but is otherwise unchanged: it stays where it was placed,
still excludes volume, and still contributes to the energy, so the mobile chains feel
it exactly as they would any other chain.

To discover which chainID is which, run once with ``WRITE_CHAIN_TO_CHAINID : True``,
which writes ``chain_to_chainid.txt`` mapping every chainID to its length and
sequence.

.. note::

   Freezing is currently at **whole-chain** granularity: a chain is either entirely
   frozen or entirely free. A per-bead freeze directive (a ``B`` line) is reserved in
   the file format but is not yet implemented.

The freeze-file *workflow* - capturing a structure into a restart file, freezing it,
and letting new chains (added with ``EXTRA_CHAIN``) explore around it, including how
freezing composes with :doc:`parallelization <advanced/parallelization>` - is
covered in detail on the :doc:`freeze files <advanced/freeze>` page.
