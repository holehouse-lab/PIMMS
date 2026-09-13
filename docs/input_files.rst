.. _input-files:

===========
Input files
===========

A PIMMS run is driven by plain-text input files. There are up to three:

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

A run can take two further inputs, neither of which follows the conventions above
and both of which are documented elsewhere: a :doc:`restart file <restart_files>`
(``RESTART_FILE``), a binary Python pickle that seeds a run from a previous
configuration, and a user-supplied Python analysis module (``ANALYSIS_MODULE``),
which is imported and run during the simulation.

.. _input-keyfiles:

Keyfiles
========

The keyfile is the top-level description of a simulation. Each non-comment line
sets one keyword::

    KEYWORD : value

Whitespace around the colon is optional, and only the *first* colon separates the
keyword from its value, so a value may itself contain one (a Windows-style path,
say). A few rules govern the file as a whole:

* **Every non-comment line must contain a colon.** A line without one is an error,
  so stray text cannot sit unnoticed in a keyfile.
* **Keyword names are case-insensitive.** ``N_STEPS``, ``n_steps`` and ``N_steps``
  are the same keyword. How the *value* treats case is keyword-specific: booleans
  are read in any case, and chain sequences are upper-cased unless
  ``CASE_INSENSITIVE_CHAINS : False``.
* **An unrecognised keyword is an error.** PIMMS refuses to run rather than
  ignoring a mistyped keyword and quietly using the default in its place.
* **Most keywords may appear at most once.** A repeated keyword is an error (PIMMS
  reports the offending keyword rather than silently taking the last value). The
  exceptions are ``CHAIN``, ``EXTRA_CHAIN`` and ``ANA_RESIDUE_PAIRS``, which may be
  repeated (to build a multi-component system, or to monitor several residue
  pairs).
* **Required keywords.** ``DIMENSIONS``, ``PARAMETER_FILE``, ``TEMPERATURE``,
  ``N_STEPS`` and ``EQUILIBRATION`` must always be present, plus ``CHAIN`` - unless a
  ``RESTART_FILE`` is provided, in which case the chains come from the restart file
  and ``CHAIN`` may be omitted. Every other keyword falls back to a default.
* **Values are checked on read.** A keyword expecting an integer, float or boolean
  that is given a malformed value fails immediately with a descriptive error, rather
  than deep inside the run. Booleans are written ``True`` or ``False`` (any case);
  ``nan`` and ``inf`` are rejected for numeric keywords. The path keywords
  (``PARAMETER_FILE``, ``RESTART_FILE``, ``FREEZE_FILE``, ``ANALYSIS_MODULE``)
  expand a leading ``~`` and reject an empty value - if you do not want the
  feature, remove the keyword rather than leaving it blank.

A few keywords that shape the file are worth calling out here (the full list, with
types and defaults, is the :doc:`keyword reference <keywords>`):

* ``DIMENSIONS`` takes 2 values (a 2D simulation) or 3 (3D). Each value must be at
  least **7** lattice units - the smallest box that can support the super-long-range
  interaction shell - and the axes need not be equal (non-cubic/non-square boxes are
  fully supported).
* ``CHAIN : N SEQUENCE`` declares ``N`` copies (``N`` >= 1) of a chain whose
  one-letter ``SEQUENCE`` names the bead types (e.g. ``CHAIN : 20 QQQQQQQQQQ`` for
  20 ten-bead poly-Q chains). Repeat the keyword for a mixture; each line is a
  separate chain type. Sequences are upper-cased on read unless
  ``CASE_INSENSITIVE_CHAINS : False``; every bead letter used must be defined in the
  parameter file, and ``0`` may not appear in a sequence because it is the solvent
  symbol.
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

All interaction energies and *absolute* ``ANGLE_PENALTY`` values must be
**integers** - a float is rejected with a clear error; the temperature-normalised
``ANGLE_PENALTY_T_NORM`` values are floats. Applied energies are stored as signed
32-bit integers, so values (including temperature-scaled and rounded angle
penalties) must lie between ``-2147483648`` and ``2147483647``; values outside
that range are rejected instead of wrapping to a different energy. The file has
a few kinds of line.

**1. Pairwise interactions** act over three nested length scales, set by the
Chebyshev distance between two beads:

.. code-block:: text

   ## R1 R2  e_SR [e_LR [e_SLR]]
   A  A   -8                 # A-A short-range contact energy (distance 1)
   A  B   -3  -2             # A-B short-range AND long-range (distance 2)
   B  B   -6  -3   3         # B-B short, long AND super-long-range (distance 3)

Short-range (SR) is always present; long-range (LR, distance 2) and super-long-range
(SLR, distance 3) are optional trailing columns, so an interaction line carries 3,
4 or 5 columns and anything else is a parse error. A bead type that appears in any
line with an LR column becomes "LR-capable" and is subsequently tested at
Chebyshev distances 2 and 3 as well as 1.

**2. Solvation** is the bead-solvent energy, written as an interaction with type
``0``:

.. code-block:: text

   ## R 0  e_solv
   A  0   -2
   B  0   -1

**3. Backbone-angle penalties** bias the local chain geometry (three values per
residue, keyed to the *displacement class* of the ``i-1`` to ``i+1`` vector:
``A1`` = the two flanking beads are Chebyshev-adjacent, ``A2`` = mixed
displacement, ``A3`` = every non-zero component of the displacement is 2 -
which includes the straight-through geometry but also some 70-110 degree
bends, so each class mixes several geometric bend angles). Use either absolute
integer penalties or temperature-normalised ones:

.. code-block:: text

   ## absolute integer penalties:
   ANGLE_PENALTY         A   30 10 0

   ## ...or temperature-normalised (units of kT, k=1; multiplied by TEMPERATURE):
   ANGLE_PENALTY_T_NORM  A   0.5 0.2 0

(For a ``QUENCH_RUN`` the multiplier is ``QUENCH_END``, the production
temperature, and it is applied once at parse time - the penalties are exact in kT
only at that temperature, not along the ramp.) Because lattice energies are
integers, a scaled T-normalised penalty is rounded to the nearest integer before
it is used; both the requested and the applied value are written to
``absolute_energies_of_angles.txt``, so you can check what the run actually
applied.

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
  complete: any pair you omit defaults to 0, and PIMMS prints a warning naming each
  pair it filled in, so an accidentally missing LR term is visible at start-up.
* **Angles are optional as a whole - but all-or-nothing.** Use ``ANGLE_PENALTY``
  or ``ANGLE_PENALTY_T_NORM`` (the T-normalised form keeps stiffness fixed
  relative to temperature), one line of exactly five fields per bead type; a bead
  type given two angle lines is a parse error. If angles are enabled (the default),
  **every** bead type with interaction energies must have an angle line - a missing
  one is a parse error. Setting ``ANGLES_OFF : True`` in the keyfile disables angles
  entirely, and no angle lines are then needed.
* ``NON_INTERACTING : True`` in the keyfile zeroes all **pairwise** interaction
  and solvation energies, regardless of what the parameter file says. Angle
  penalties are unaffected - combine with ``ANGLES_OFF : True`` for a fully ideal
  excluded-volume-only reference run.

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

Each ``C`` line contributes its integer chainIDs to the frozen set; order, repeats
and the split across lines do not matter. A ``C`` line with no IDs after it, an ID
that is not an integer, or a line beginning with anything other than ``C`` (or the
reserved ``B``) is an error, reported with the line number. A frozen chain is
**excluded from the pool of chains PIMMS can move** but is otherwise unchanged: it
stays where it was placed, still excludes volume, and still contributes to the
energy, so the mobile chains feel it exactly as they would any other chain. The
collective moves (cluster translate/rotate, VMMC) additionally reject any move
whose cluster would contain a frozen chain. Freezing *every*
chain is refused at start-up, since the configuration could then never change.

To discover which chainID is which, run once with ``WRITE_CHAIN_TO_CHAINID : True``,
which writes ``chain_to_chainid.txt``: one tab-separated line per chain giving its
chainID, its length and its sequence. chainIDs are numbered from 1.

.. note::

   Freezing is currently at **whole-chain** granularity: a chain is either entirely
   frozen or entirely free. A per-bead freeze directive (a ``B`` line) is reserved in
   the file format but is not yet implemented, and using one aborts the run with an
   explicit "not implemented" error.

The freeze-file *workflow* - capturing a structure into a restart file, freezing it,
and letting new chains (added with ``EXTRA_CHAIN``) explore around it, including how
freezing composes with :doc:`parallelization <advanced/parallelization>` - is
covered in detail on the :doc:`freeze files <advanced/freeze>` page.
