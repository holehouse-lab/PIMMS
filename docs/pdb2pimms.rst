.. _pdb2pimms:

==========================================
Starting from a structure: ``pdb2pimms``
==========================================

``pdb2pimms`` turns a PDB or mmCIF structure (an experimental structure, an
AlphaFold model, a frame from another simulation) into a PIMMS
:ref:`restart file <restart-files>`, with one bead per residue, and writes a
Markdown report of everything the conversion did. The report also lists every
bead type the run's parameter file has to define, which is usually the first
thing you need to know. A :ref:`freeze file <advanced-freeze>` listing every
chain is written too, so the structure can be held fixed while other chains move
around it, optionally with chosen residues (a disordered tail, a loop) left free.
Give it a parameter file and it also writes a keyfile that runs the structure.

It is installed alongside ``PIMMS``:

.. code-block:: text

   pdb2pimms 1ubq.cif

writes ``1ubq_restart.pimms``, ``1ubq_report.md`` and ``1ubq_freeze.txt`` into
the current directory and prints a summary, including the bead types:

.. code-block:: text

   Wrote 1ubq_restart.pimms: 1 chain(s), 76 beads, box 14 14 16, fit RMSD 1.83 A
   Wrote 1ubq_report.md
   Wrote 1ubq_freeze.txt: freezes all 1 chain(s); a run using it needs mobile chains too (EXTRA_CHAIN)
   No keyfile written: give --parameter-file and --temperature, or --protein-mode (which uses gcf_rje23_v14.prm), for one that runs the structure
   Bead types the parameter file must define (18): A D E F G H I K L M N P Q R S T V Y

For a protein, ``--protein-mode`` does the rest. It asks
`metapredict <https://github.com/idptools/metapredict>`_ which residues are folded
and which are disordered, freezes the folded ones, and writes a keyfile that runs
the result with the amino-acid forcefield that ships with PIMMS. For the AlphaFold
model of FUS:

.. code-block:: text

   pdb2pimms AF-P35637-F1-model_v6.cif --protein-mode

.. code-block:: text

   Wrote AF-P35637-F1-model_v6_restart.pimms: 1 chain(s), 526 beads, box 47 42 42, fit RMSD 1.83 A
   Wrote AF-P35637-F1-model_v6_report.md
   metapredict 3.1.0: chain A, 526 amino acids: folded A:THR286-ARG372; 439 disordered
   Wrote AF-P35637-F1-model_v6_freeze.txt: freezes 87 of 526 beads (0 chain(s) entirely, 1 partly, 0 left free); ...
   Wrote gcf_rje23_v14.prm: a copy of the gcf_rje23_v14.prm that ships with PIMMS
   Wrote AF-P35637-F1-model_v6_keyfile.kf: run it from <this directory> with: PIMMS -k AF-P35637-F1-model_v6_keyfile.kf
   Bead types the parameter file must define (20): A C D E F G H I K L M N P Q R S T V W Y

The RRM (residues 286-372) is held fixed and the two disordered regions either
side of it move. ``PIMMS -k AF-P35637-F1-model_v6_keyfile.kf`` then runs it (about
5 minutes). See `The keyfile`_ and `Protein mode`_.

Running it
==========

The input may be a PDB file (``.pdb``, ``.ent``) or a PDBx/mmCIF file (``.cif``,
``.mmcif``), either optionally gzipped (``.gz``); a file with some other name is
recognised from its content. No output file is replaced unless you give
``--overwrite``, no output may be the input structure itself, and two names that
differ only by case count as one file (they are one on the default file systems
of macOS and Windows). All of them are written under temporary names and moved
into place only once the restart file, the freeze file and the keyfile have been
read back and accepted by PIMMS's own readers, so a failed conversion leaves
nothing behind.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Option
     - Effect
   * - ``-o``, ``--output FILE``
     - The restart file to write. Default ``<structure name>_restart.pimms`` in
       the current directory.
   * - ``--report FILE``
     - The report to write. Default: next to the restart file, named after it
       (``x_restart.pimms`` or ``x.pimms`` give ``x_report.md``).
   * - ``--freeze FILE``
     - The freeze file to write. Default: next to the restart file, named after
       it (``x_freeze.txt``).
   * - ``--no-freeze``
     - Write no freeze file.
   * - ``--mobile SELECTION``
     - Leave these residues free in the freeze file (repeatable): ``CHAIN`` for a
       whole chain, or ``CHAIN:NUMBER`` / ``CHAIN:FIRST-LAST``, with the
       structure's own chain IDs (``-`` for a blank one) and residue numbers,
       which may be negative (``I:-73--60``). A PDB file that numbers past 9999
       in hybrid-36 (``A000``, ``A001``, ...) is selected by the numbers those
       stand for (10000, 10001, ...). See `Leaving residues free`_.
   * - ``--protein-mode``
     - Freeze the folded domains metapredict predicts and leave the disordered
       regions free, and write a keyfile that runs it with
       ``gcf_rje23_v14.prm``. Needs metapredict (``pip install metapredict``).
       See `Protein mode`_.
   * - ``--parameter-file FILE``
     - The parameter file the keyfile names; needs ``--temperature``. It must
       define every bead type of the structure. Default with
       ``--protein-mode``: a copy of ``gcf_rje23_v14.prm``. See `The keyfile`_.
   * - ``--temperature T``
     - The keyfile's ``TEMPERATURE``. Required with ``--parameter-file``; with
       ``--protein-mode``'s default parameter file it's 50.
   * - ``--keyfile FILE``
     - The keyfile to write. Default: next to the restart file, named after it
       (``x_keyfile.kf``). A keyfile is only written when there's a parameter
       file for it.
   * - ``--no-keyfile``
     - Write no keyfile.
   * - ``--spacing ANGSTROMS``
     - Lattice spacing in Angstroms per site. Default 3.65, the default
       ``LATTICE_TO_ANGSTROMS``.
   * - ``--padding SITES``
     - Empty sites between the structure and each face of the box. Default 3,
       the reach of the super-long-range interaction.
   * - ``--box X Y Z``
     - Use this box instead, with the structure centred in it.
   * - ``--model N``
     - The model to convert. Default: the first model in the file.
   * - ``--chains A,B``
     - Convert only these (author) chains; ``-`` stands for a blank chain ID,
       as in the report.
   * - ``--bead-position``
     - ``backbone`` (default) places amino acids at CA and nucleotides at C4';
       ``centroid`` places every residue at its heavy-atom centroid.
   * - ``--hetero``
     - Also convert every non-water hetero group (ligand, ion, cofactor) into a
       one-bead chain.
   * - ``--map RESNAME=CHARACTER``
     - Use ``CHARACTER`` as the bead type of residue ``RESNAME`` (repeatable).
   * - ``--overwrite``
     - Replace existing output files.

A value that starts with ``-`` (a blank chain ID first) must be joined to its
option with an equals sign, ``--mobile=-:1-5`` or ``--chains=-,A``; otherwise it
is read as an option of its own.

The exit status is 0 on success, 1 if the conversion failed (the reason is
printed on standard error) and 2 for a usage error.

How residues become beads
=========================

Every residue of the chosen model becomes one bead, whose type is a single
character, because a PIMMS sequence has one character per bead:

* **Amino acids** keep their one-letter codes (``ALA`` is ``A``, and so on;
  selenocysteine is ``U``, pyrrolysine ``O`` and the unknown amino acid ``UNK``
  is ``X``).
* **Nucleotides** get lower-case codes, so that they never share a bead type
  with an amino acid: ``a``, ``c``, ``g``, ``u`` for RNA and ``a``, ``c``, ``g``,
  ``t`` for DNA. The parameter file is case-sensitive, so ``a`` and ``A`` are
  different bead types. A structure holding both RNA and DNA therefore uses one
  type for, say, adenosine and deoxyadenosine; the report warns when this
  happens, and ``--map DA=...`` separates them. Residues named ``A``, ``C`` or
  ``G`` are RNA in a current PDB file but DNA in older files, so their sugar
  decides: a deoxyribose makes them DNA. The unknown nucleotides ``N`` and
  ``DN`` are ``n``.
* **Simulation-package names** for standard residues are recognised too:
  AMBER, CHARMM and GROMACS protonation states (``HID``, ``HSE``, ``HISE``,
  ``CYX``, ``LYSN``, ...), AMBER's terminal residues (``NALA``, ``DA5``, ``RA3``,
  and for histidine and bridged cysteine ``NHIE``, ``CHID``, ``CCYX``, ...)
  and R-prefixed RNA (``RA``), CHARMM's base-only nucleotides (``ADE``, ...) and
  the ``T`` of older PDB files. Many of these are also codes of unrelated
  molecules in the PDB Chemical Component Dictionary (``RU`` is a ruthenium ion,
  ``GUA`` glutaric acid), so such a name counts only for a residue that has the
  backbone of what it claims to be (N, CA and C for an amino acid; C4' with a P,
  O3' or O5' for a nucleotide; or that single atom in a CA-only or C4'-only
  model). A calcium ion is also a single atom named CA, but not a carbon (and
  its residue is named CA, which no amino acid is), so it is a hetero group
  like any other ion.
* **Modified residues** take the code of the residue they derive from. The file
  itself says which that is (PDB ``MODRES`` records, mmCIF
  ``_pdbx_struct_mod_residue``), and its word comes before everything but a
  standard name; a built-in table covers common modifications (``MSE`` to
  ``M``, ``SEP`` to ``S``, ``PSU`` to ``u``, ...) when it says nothing. The report lists every residue name behind each bead type and where
  its parent came from.
* A residue whose name is not known but which has an amino-acid backbone (N, CA
  and C) becomes the unknown amino acid ``X``; one with a nucleotide backbone
  becomes the unknown nucleotide ``n``. Both draw a warning.
* **Water** is never converted. ``HOH`` and ``DOD`` are always water; a
  simulation package's water name (``WAT``, ``TIP3``, ``SOL``, ``T3P``, ...)
  counts only for a residue with one oxygen and no carbon, nitrogen, phosphorus
  or sulfur, since ``SOL`` and ``T3P`` are also other molecules' codes.
* A HETATM residue with an amino-acid or nucleotide name (or backbone) that is
  bonded to neither neighbour is a ligand (a free arginine in a binding pocket,
  an AMP, an ATP), not part of a chain, and is treated as a hetero group. A
  modified residue inside a chain is also HETATM, but it is bonded and is kept.
  For this test a residue with more atoms than a CA-only or C4'-only model's
  must have its own atom of the linking bond, so a nucleoside or a nucleoside
  triphosphate (whose alpha phosphorus is named PA, not P) is never joined to
  the end of a chain by the backbone distance alone.
* Other **hetero groups** are left out unless you give ``--hetero``, which makes each one a single bead at its heavy-atom
  centroid, never bonded to anything, with a lower-case letter that no
  nucleotide uses (``b``, ``d``, ``e``, ...) unless ``--map`` names one.

``--map`` overrides any of these. A bead character must be a single printable,
non-space character other than ``0`` (the solvent type) and ``#`` (which starts a
parameter-file comment).

Where alternate locations are given, the location whose atoms have the highest
mean occupancy is kept for each residue (the first listed wins a tie), including
in files that write the alternate locations of several residues as separate
blocks; a ligand that only shares a residue's chain and number stays a residue
of its own, unless that residue has alternate locations too (the two cannot then
be told apart, and the ligand's atoms join the residue). Hydrogens never affect a bead's position, except in a residue that
has no heavy atom, and each belongs to the residue with its chain and number
even when the file lists every hydrogen after the whole chain or names them all
``H`` (Open Babel does both). Where a file gives no elements (a PDB file without
element columns, an mmCIF file without ``_atom_site.type_symbol``), hydrogens
are recognised by their names. A heavy atom whose name repeats within a residue
starts a new residue with the same chain and number, which tells apart
consecutive residues that share a number but also splits a molecule whose atoms
are not uniquely named (Open Babel names a ligand's atoms by element alone);
the report warns whenever that makes a bead, so give such a molecule unique
atom names.

A PDB file whose frames are separated by ``END`` or ``ENDMDL`` records with no
``MODEL`` records is read as one model per frame, so ``--model`` picks a frame.
A ``MODEL`` number used for two blocks of atoms (as in two single-model files
joined together) is refused rather than read as one doubled model.

Chains and bonds
================

Consecutive residues of one chain are bonded when the atoms of the linking
bond are within 2.0 A of each other: C-N for a peptide (1.33 A in a real
structure) and O3'-P for a phosphodiester (1.61 A). Where those atoms are
missing, as in a CA-only model, the CA atoms must be within 4.2 A, or for a
nucleic acid the C4' atoms (or, failing those, the P atoms) within 7.5 A. A PIMMS chain has to be connected, so a chain with missing
residues becomes one PIMMS chain per bonded run, and the report lists each
break and its cause. Chains with identical sequences share a chainType.

Fitting the structure onto the lattice
======================================

Each bead's position is divided by the lattice spacing, and every chain then
has to become a self-avoiding walk in which consecutive beads are Moore
(Chebyshev-1) neighbours. Each chain is placed by dynamic programming: among
the walks that keep every bead within two sites of its scaled position (three
if no such walk exists) and avoid the chains already placed, it finds the one with the smallest total
squared distance between the beads and their positions. Where that walk
visits a site twice, the site is forbidden to the bead further from it and the
walk is found again. A greedy search with backtracking is the fallback when no
walk exists within three sites, or when repairing the walk does not converge,
and a chain that cannot be placed at all (beads packed more densely than the
lattice can hold) stops the conversion with a message suggesting a smaller
spacing.

The report gives the fit as the root-mean-square deviation between the lattice
beads and the input positions, after the best rigid-body superposition, per
chain and overall, together with the largest single deviation, each chain's
radius of gyration before and after, and the bond lengths. Rounding positions
to a cubic lattice on its own costs about half a spacing (1.8 A at the default
3.65 A), and a protein fits at about that. Consecutive C4' atoms of a nucleic
acid are 5 to 7 A apart, more than a lattice step can always span, so nucleic
acids fit less well at the default spacing (3.1 A for a tRNA); a larger
``--spacing`` reduces that at the cost of a coarser lattice. Too large a
spacing fails the other way: a folded protein holds about one residue per
135 cubic Angstroms (a mean residue mass of 110 Da at a partial specific
volume of 0.73 cubic centimetres per gram), so at 6 A spacing (216 cubic Angstroms per site) it has
more residues than sites and has to be squashed. The report warns when the
RMSD is more than 1.5 times the half-spacing floor.

The restart file
================

The restart file is a starting configuration for ``RESTART_FILE``. It holds the
chains, the box and the boundary flag only:

* it is written with ``HARDWALL : True``, because no chain crosses a face of the
  box. That makes it valid for a run with either boundary condition, and lets
  ``DIMENSIONS`` be any box at least as large as the one written: PIMMS centres
  the restart box inside a larger one (see :ref:`restart-complexities`);
* its ``ENERGY`` is a placeholder (0), which no run reads, because PIMMS computes
  the energy from the parameter file at start-up;
* it records no step, random-number state or temperature, so it starts a new run,
  and ``RESTART_CONTINUE`` refuses it.

A run started from it needs, besides the usual keywords:

.. code-block:: text

   RESTART_FILE         : 1ubq_restart.pimms
   DIMENSIONS           : 14 14 16      # or any larger box
   LATTICE_TO_ANGSTROMS : 3.65          # the spacing used for the conversion
   PARAMETER_FILE       : params.prm

No ``CHAIN`` lines are needed (any are ignored). ``LATTICE_TO_ANGSTROMS`` only
scales ``START.pdb`` and ``traj.xtc``; setting it to the conversion's spacing keeps
them on the scale of the input structure. The parameter file must define every
bead type in the report: a solvation line for each, an interaction line for
every pair of types (each type with itself included), and an ``ANGLE_PENALTY``
or ``ANGLE_PENALTY_T_NORM`` line for each unless the keyfile sets
``ANGLES_OFF : True`` (see :ref:`input-parameter-files`).

The freeze file
===============

A converted folded domain keeps its fold only as long as the parameter file
favours it. To hold the structure fixed instead, while other chains move around
it, ``pdb2pimms`` writes a :ref:`freeze file <advanced-freeze>` that lists every
chain of the restart file, one per line, with the residues each came from:

.. code-block:: text

   # PIMMS freeze file written by pdb2pimms. It freezes every chain of 4HHB_restart.pimms,
   # converted from 4HHB.cif. PIMMS will not run with every chain frozen, so a run that
   # uses this file needs mobile chains as well: EXTRA_CHAIN adds them, numbered
   # from 11. Delete a line to let that chain move.
   C 1   # A:VAL1-ARG141, 141 bead(s)
   C 2   # B:VAL1-HIS146, 146 bead(s)
   ...

Frozen chains still occupy their sites and still contribute to the energy, so
the parameter file must define their bead types as before. PIMMS refuses to run
with every chain frozen, so a run that uses the freeze file adds mobile chains
with ``EXTRA_CHAIN``, which are numbered after the restart file's chains and so
leave the freeze file correct:

.. code-block:: text

   RESTART_FILE         : 4HHB_restart.pimms
   DIMENSIONS           : 40 40 40      # room for the added chains
   LATTICE_TO_ANGSTROMS : 3.65
   PARAMETER_FILE       : params.prm
   FREEZE_FILE          : 4HHB_freeze.txt
   EXTRA_CHAIN          : 20 KEKEKEKE

``EXTRA_CHAIN`` sequences are upper-cased unless the keyfile sets
``CASE_INSENSITIVE_CHAINS : False``, which an added chain with lower-case bead
types (nucleotides) needs. To let a chain move, delete its line; to let
everything move, leave out ``FREEZE_FILE``. The file is read back with PIMMS's
own freeze-file reader before it is put in place, and ``--no-freeze`` skips it.

Leaving residues free
---------------------

A folded domain often has disordered tails or loops that should move while the
fold stays put. ``--mobile`` names them, by the structure's own chain IDs and
residue numbers, and the freeze file then freezes everything else: a chain with
no selected residue keeps its ``C`` line, a chain with some gets a ``B`` line
listing its frozen beads by their index along the chain (counted from 0), and a
chain whose residues are all selected is left out. For the nucleosome core
particle of 1KX5, freeing the two H3 N-terminal tails:

.. code-block:: text

   pdb2pimms 1kx5.cif.gz --mobile A:1-26 --mobile E:1-26

gives

.. code-block:: text

   C 1   # I:DA-73-DT73, 147 bead(s)
   C 2   # J:DA-73-DT73, 147 bead(s)
   B 3  26-134  # A:ALA1-ALA135: 109 of 135 bead(s) frozen; free: beads 0-25 (A:ALA1-ARG26)
   ...
   B 7  26-134  # E:ALA1-ALA135: 109 of 135 bead(s) frozen; free: beads 0-25 (E:ALA1-ARG26)
   ...

Only the crankshaft (``MOVE_CRANKSHAFT``) and the experimental segment regrowth
(``MOVE_SEGMENT_REGROWTH``) move the free beads of a partly frozen chain - every
other move leaves such a chain alone (see :ref:`freeze-part-of-a-chain`) - so a
run with this freeze file needs one of them. A chain freed whole is an ordinary mobile chain. The
run needs no ``EXTRA_CHAIN``, since the free residues can move on their own;
added chains are still numbered after the restart file's. A chain ID selects
every converted residue that carries it, so with ``--hetero`` it also frees the
one-bead hetero groups (ions, ligands) of that chain, which then move as chains
of their own; the report lists them. A selection that matches no converted
residue, and a set of selections that leaves nothing frozen, stop the conversion
with an error, and the report lists the free residues of every chain.

The demo in ``demos/nucleosome_from_pdb`` works through all of this: it
converts the nucleosome core particle of PDB entry 1KX5, completes the
parameter file with the DNA bead types the report asks for, freezes the
nucleosome except its histone tails (the N-terminal tail of all eight histones
and the C-terminal tail of both H2As, one ``--mobile`` selection each), adds two
disordered proteins with ``EXTRA_CHAIN`` and maps where they bind.

The keyfile
===========

Given a parameter file (``--parameter-file`` with ``--temperature``, or
``--protein-mode``'s default), ``pdb2pimms`` also writes a keyfile,
``x_keyfile.kf``, that runs the converted structure as it stands. Its paths are
relative to its own directory, so run it from there:

.. code-block:: text

   PIMMS -k 1ubq_keyfile.kf

Every line has a comment saying why it's there, and the report explains the
choices. In short:

* **The system.** ``RESTART_FILE``, ``LATTICE_TO_ANGSTROMS`` (the conversion's
  spacing), ``PARAMETER_FILE`` and ``TEMPERATURE``.
* **The box.** ``HARDWALL : True`` with the restart file's box. With walls rather
  than periodic boundaries, no chain can meet a periodic image of itself, which in
  a box this tight a periodic run would allow. The walls do confine the chains (a
  wall counts as solvent), so for more room convert with a larger ``--padding`` or
  ``--box``, or raise ``DIMENSIONS`` in the keyfile: PIMMS centres the structure in
  a bigger box.
* **Freezing.** ``FREEZE_FILE`` only if the freeze file leaves something free
  (``--mobile`` or ``--protein-mode``). The default freeze file freezes every
  chain, which PIMMS only runs with added chains, so then the keyfile leaves it out
  and the whole structure moves.
* **Moves.** ``MOVE_CRANKSHAFT : 1`` when a chain is partly frozen, since nothing
  else moves its free beads. Otherwise 90% crankshaft and 10% ``MOVE_CHAIN_PIVOT``,
  as in the single-chain protein demo.
* **Length.** 1,000,000 steps, the first 100,000 of them equilibration, with
  ``CRANKSHAFT_SUBSTEPS`` set to the number of movable beads (at least 500), so
  each step tries every movable bead about once. The trajectory gets 1,000 frames.
  For FUS with ``--protein-mode`` that takes 5 minutes and 130 MB, and the energy
  stops drifting after about 200,000 steps.

The run length is a starting point, not a guarantee of convergence. Check that
``ENERGY.dat`` has stopped drifting and that whatever you measure has settled, and
lengthen the run if not. With your own parameter file, ``pdb2pimms`` refuses to
write a keyfile unless that file defines every bead type of the structure (an
interaction with every type, a solvation term and an angle penalty for each).

With ``--protein-mode``'s default, ``gcf_rje23_v14.prm`` (the residue-level
amino-acid forcefield in ``pimms/data``) is copied next to the keyfile, so the run
directory holds everything the run reads. It's run at ``TEMPERATURE : 50``, the
temperature it was parameterised at. A copy that's already there and identical is
kept, so converting several structures into one directory needs no
``--overwrite``; a copy you've edited is never replaced without it.

Protein mode
============

``--protein-mode`` is for proteins with folded domains and disordered regions,
where you want the domains held fixed and the disordered regions to move. It runs
`metapredict <https://github.com/idptools/metapredict>`_ (Emenecker, Griffith and
Holehouse, *Biophys. J.* **120**, 4312, 2021) on each protein chain, then:

* freezes the residues in the folded domains it predicts, with ``B`` lines;
* leaves the residues in the disordered regions free, along with any ``--mobile``
  residues;
* freezes nucleotides and hetero groups, unless ``--mobile`` frees them;
* writes a keyfile using ``gcf_rje23_v14.prm``, unless ``--parameter-file`` names
  another.

metapredict isn't a PIMMS dependency, so install it first (``pip install
metapredict``); without it ``--protein-mode`` stops with an error. It runs with
its ``V3`` network and metapredict 3's default domain settings (a disorder
threshold of 0.5, a minimum disordered region of 12 residues, a minimum folded
domain of 50 and a gap closure of 10). These are given explicitly, so a later
metapredict release can't silently change which residues are frozen. The report
records the version and lists every folded domain and disordered region.

A few details:

* **One sequence per chain.** Each source chain is predicted from its converted
  amino acids in order. Where a chain has missing residues, the runs either side
  of the gap are joined, since metapredict needs one sequence.
* **The chemistry, not the bead type.** metapredict reads each residue's own amino
  acid, so ``--map`` can rename a bead type without changing the prediction.
  Selenocysteine is read as cysteine and pyrrolysine as lysine, with a warning. An
  unknown amino acid (``X``) stops the conversion, unless ``--map`` gives it a
  standard bead type, which metapredict then reads.
* **Nothing folded, or nothing disordered.** If every residue is predicted
  disordered there's nothing to freeze, so no freeze file is written and the
  keyfile moves everything. If none is, the freeze file freezes every chain and the
  keyfile leaves it out, so the whole structure moves; a warning says which.
* **Short chains.** metapredict warns about a sequence no longer than its
  24-residue smoothing window, and the report keeps the warning.

Small folded domains can come out disordered. In FUS, metapredict's defaults call
the small zinc finger (its zinc-binding cysteines are residues 428, 433, 444 and
447) disordered, so it's left free.
``--protein-mode`` can only add free residues, never freeze ones metapredict calls
disordered. To hold such a domain, write the ``B`` lines yourself, or convert
without ``--protein-mode`` and name the free residues with ``--mobile``.

The report
==========

The report (Markdown, readable as plain text) records:

* the input file with its SHA-256 digest, the command line, the PIMMS version and
  the time;
* a summary: model, chains, beads, bead types, box, spacing and the fit;
* every warning (bead types shared by different kinds of residue, unknown
  residues, residues split by a repeated atom name, a poor fit, a box too large
  for the PDB coordinate field);
* the bead types the parameter file must define, with the residue names behind
  each and how many beads have it;
* with ``--protein-mode``, metapredict's version and settings and every folded
  domain and disordered region;
* the keyfile, if one was written, with what it sets and why, and otherwise the
  keyfile lines for using the restart file and the freeze file;
* the residues the freeze file leaves free;
* every PIMMS chain, with its chainID, chainType, source residues, fit and
  sequence;
* every chain break and its cause;
* what was not converted: waters, hetero groups, hydrogens, the alternate
  locations not kept, other models and unselected chains;
* how the conversion was done, with the bond statistics.

Limitations
===========

* The lattice is always 3D; a 2D restart file cannot be made from a structure.
* One bead per residue. A hetero group kept with ``--hetero`` is a single bead
  however large it is.
* Only the bonds along each chain are represented. Disulfides, cross-links and
  bonds to ligands are not (a PIMMS chain is linear).
* For a configuration from a PIMMS run, use that run's own ``restart.pimms``
  rather than converting its ``START.pdb``.
* ``--protein-mode`` predicts disorder from the residues in the file, not from the
  full sequence of the protein, so a structure with missing residues is predicted
  as if its gaps were closed.
