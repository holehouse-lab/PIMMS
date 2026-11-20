.. _advanced-freeze:

============
Freeze files
============

A **freeze file** holds chosen chains, or chosen beads of a chain, rigidly fixed
for the entire simulation. Frozen beads never move, but they are still fully
present on the lattice: they occupy their sites and contribute to the energy, so
the mobile chains feel them exactly as they would any other chain. This is the tool
for building a fixed scaffold, a wall or surface, or a pre-formed template that the
rest of the system explores around - and, by freezing only part of a chain, for a
folded domain whose disordered tail or loop stays free.

Turning it on
=============

Point the ``FREEZE_FILE`` keyword at a plain-text file that lists what to freeze.
This is the only immobilization mechanism. The simulation aborts at startup if the
file does not exist, or if ``FREEZE_FILE`` is given with an empty value.

.. code-block:: text

   FREEZE_FILE : freeze.txt

The freeze file itself is a short list of directives, one per line. Lines beginning
with ``#`` are comments, and a ``#`` anywhere on a line comments out the rest of it:

.. code-block:: text

   # freeze.txt  -  lines beginning with # are comments
   C 1 2 3          # freeze chainIDs 1, 2 and 3 (chainIDs are numbered from 1)
   C 10 11 12 13    # more C lines are allowed; IDs may be split across lines
   B 5 20-134       # freeze beads 20 to 134 of chain 5 (beads are numbered from 0)
   B 6 0-9 40 41    # ranges and single beads may be mixed

A ``C`` line freezes every bead of each chain it names. A ``B`` line freezes some
beads of one chain: the chainID comes first, then the beads, each written as its
index along the chain counted from 0 (as in ``ANA_RESIDUE_PAIRS``) or as an
inclusive range ``first-last``. Both directives may be repeated; the order and the
split across lines do not matter, blank lines are skipped, and a chain or bead named
more than once is frozen once. Every directive must contain at least one ID (and a
``B`` line at least one bead after its chainID); an empty directive, a non-integer
ID, a bead that is neither an index nor an ascending range, or an unknown directive
(including a lower-case ``c`` or ``b``) is a parse-time error rather than a silently
ignored instruction. A chainID that does not exist in the system, or a bead past the
end of its chain, is also an error, raised once the lattice has been built and
naming what was asked for and what is there. A file that contains no directives at
all (only comments or blank lines) is accepted and freezes nothing; the start-up
summary then shows ``Chains to freeze : []``.

The one thing a freeze file cannot do is freeze everything: if every bead in the
system is frozen, PIMMS refuses to start rather than write out ``N_STEPS`` copies
of an unchanging configuration.

The frozen set is echoed in the start-up summary (here for a file freezing chains
1, 2 and 3, and beads 20 to 134 of chain 5) and written to the log::

   --> Freeze File Settings

   Freeze file      : freeze.txt
   Chains to freeze : [1, 2, 3]
   Beads to freeze  : chain 5, beads 20-134

What "frozen" means
===================

A chain with any frozen bead - a chain on a ``C`` line, or a chain with ``B``
lines - is **excluded from the pool of chains PIMMS can pick to move**, at every
layer:

* The outer-loop chain selector never returns such a chain, so none of the
  single-chain moves (translate, rotate, pivot, head pivot, jump-and-relax, chain
  TSMMC) can ever be proposed for one.
* The slither and pull megamoves build their chain selectors from the chains with
  no frozen bead only - in the parallel kernels this is enforced by an explicit
  per-bead frozen mask, in which one frozen bead holds back the whole chain.
* Any *collective* move whose cluster would include such a chain (cluster
  translate, cluster rotate, VMMC seeding and VMMC recruitment) is rejected
  outright, so a frozen chain is never dragged along by its neighbours.
* The multichain TSMMC move draws its random subset from the chains with no frozen
  bead only, and a system-wide TSMMC excursion is made of ordinary sub-moves, each
  of which honours the frozen set as above.

The crankshaft is the exception, and works bead by bead: its sub-moves are drawn
from every bead that is not frozen, so the free beads of a partly frozen chain move
by crankshaft like any other bead, and its frozen beads never do. See
:ref:`freeze-part-of-a-chain` below.

Mobile chains bound to a frozen scaffold still move by every move except the
collective ones that would have to drag the scaffold along. Everything else about a
frozen bead is unchanged:

* It stays exactly where it was placed (from the ``CHAIN`` set-up or, more usually,
  from a ``RESTART_FILE``). The one exception is a ``RESIZED_EQUILIBRATION`` run:
  when the box grows, the whole configuration - frozen beads included - is
  translated rigidly into the production box (by half the difference in box size,
  which places the small box at the centre of the large one, or by
  ``EQUILIBRATION_OFFSET``), so the scaffold keeps its shape and its position
  relative to the mobile chains but not its absolute coordinates.
* It still **excludes volume** - mobile beads cannot overlap it.
* It still **contributes to the energy** - every interaction between a frozen bead
  and a mobile bead is counted normally, so the mobile chains are attracted to or
  repelled by the frozen scaffold just as they would be by a mobile partner.

Frozen means fixed in the box, not just fixed relative to the other frozen beads.
There's no way to freeze a group of beads so that it keeps its shape but can still
translate or rotate as a whole. For one frozen structure in a periodic box this
costs nothing: the mobile chains can go anywhere, so where the structure sits
doesn't change what's sampled. It does matter for two or more frozen structures
(their placement relative to each other is fixed), under a hardwall (the walls are
fixed relative to the structure), and if you care about the structure's own
diffusion.

If you want rigid bodies that move, don't freeze them. Instead enable only rigid
moves (``MOVE_CHAIN_TRANSLATE``, ``MOVE_CHAIN_ROTATE``, the cluster moves,
``MOVE_VMMC``), so that every mobile chain keeps the shape it was built with (see
:ref:`moves-irreducibility`). That applies to every mobile chain in the run:
you can't have rigid chains and flexible chains side by side. The
``star_destroyer_fleet`` demo in ``demo_keyfiles/`` does this, with a frozen Star
Destroyer and fighters that fly around it as rigid bodies.

Freezing adds essentially nothing to the cost of a step: the frozen beads simply
sit in the grid. (Partly frozen chains add a little: each crankshaft megamove
works out which beads it may draw, which takes about a millisecond with 4,000
partly frozen chains.) It does change how the megamoves' work is shared out: the
crankshaft's ``CRANKSHAFT_SUBSTEPS`` attempts are spread over the beads that are not
frozen only, and slither and pull make their per-chain ``SLITHER_SUBSTEPS`` /
``PULL_SUBSTEPS`` attempts for the chains with no frozen bead only.

.. _freeze-part-of-a-chain:

Freezing part of a chain
========================

``B`` lines freeze single beads, so part of a chain can stay free: the disordered
tails of a folded protein, a loop, or a chain tethered by one end. Two moves act on
the free beads of a partly frozen chain: the crankshaft (``MOVE_CRANKSHAFT``), one
bead per sub-move, each constrained by its bonded neighbours (frozen or not), and
the experimental :doc:`segment regrowth <../moves/segment_regrowth>`
(``MOVE_SEGMENT_REGROWTH``), which regrows a short window of free beads between
its fixed neighbours. Every other move leaves such a chain alone, as described
above. This keeps the moves correct: both draw from a set that depends only on the
frozen beads and the window's fixed neighbours, never on where the moving beads
are, which obeys detailed balance for the free beads given the frozen ones. A
pivot, a slither or a rigid move would have to move frozen beads, or be rejected
whenever it would.

Single-bead moves have a limit, though: they cannot always reach every arrangement
of the free beads, so with the crankshaft alone the free beads are sampled
correctly among the arrangements reachable from where they start, which is not
always all of them.

* **A free stretch held at both ends.** Two or more free beads with a frozen bead
  at each end can sit in an arrangement none of them can leave one bead at a time:
  two free beads that close a 2x2 square with their frozen neighbours, for example,
  never move. Such a loop explores only part of its arrangements, and what you
  measure for it can depend on where it started. PIMMS warns at start-up about
  every free stretch of two or more beads held at both ends. ``MOVE_SEGMENT_REGROWTH``
  fixes this for any stretch no longer than ``SEGMENT_REGROWTH_MAX_BEADS`` (3 by
  default, at most 4), which it can regrow whole in one move; with it enabled, the
  warning names only the longer stretches. A single free bead between two frozen
  ones is not affected (its new site is drawn from every site next to both
  neighbours), and neither, in an open space, is a free tail at a chain end.
* **Crowded frozen surroundings.** Frozen beads of any chain that pack tightly
  round a free stretch can divide its arrangements in the same way, a free tail
  included, and segment regrowth doesn't fix this for a stretch or tail longer
  than its window. In 3D it takes densely frozen surroundings; in 2D it happens
  much more easily. See :ref:`segment-regrowth-frozen-limit` for an example and
  what to do about it.

In practice:

* A run with partly frozen chains needs ``MOVE_CRANKSHAFT`` or
  ``MOVE_SEGMENT_REGROWTH`` (directly, or through a ``MOVE_SYSTEM_TSMMC`` whose
  excursions run them). If chains with no frozen bead keep the run going without
  either, PIMMS warns at start-up that the free beads of the partly frozen chains
  will never move.
* A partly frozen chain needs no mobile partner: if every chain holds a frozen bead
  but some beads are free (a protein with its tails released and nothing else), the
  run goes ahead on the crankshaft and segment regrowth alone. The moves that need
  a chain to act on are then rejected null moves, counted as attempts in
  ``MOVE_FREQS.dat``, and PIMMS says at start-up which moves they are and what
  fraction of steps they take; set them to 0. Without either move such a run is
  refused.
* A free segment moves only by local moves - single beads, or windows of up to four
  beads with segment regrowth - which is fine for a tail of a few tens of residues
  but slow for a long free loop. With the crankshaft alone it can trap a loop held
  at both ends, and with segment regrowth too if the loop is longer than the window
  and tightly packed in by frozen beads (above). Pivots on the free segment of a
  partly frozen chain are not implemented.
* ``MOVE_SEGMENT_REGROWTH`` does nothing in a system with no partly frozen chain,
  and PIMMS says so at start-up (and refuses the run if it's the only move).
* A chain whose every bead is frozen by ``B`` lines is treated exactly like a chain
  on a ``C`` line.

The start-up log lists each partly frozen chain and its frozen beads, and a
:ref:`checkpoint <restart-files>` records them, so that ``RESTART_CONTINUE`` can
refuse a continuation whose freeze file freezes different beads.

Finding the chainIDs
====================

ChainIDs are assigned internally, so to know which ID is which, run once with

.. code-block:: text

   WRITE_CHAIN_TO_CHAINID : True

This writes ``chain_to_chainid.txt``, mapping every chainID to its length and
sequence. Read off the IDs you want to pin and list them in the freeze file. When
you add chains to a restart configuration with ``EXTRA_CHAIN``, the new chains are
given fresh IDs after the existing ones, so the original (restart) chains keep the
IDs they had - which is what makes the "freeze the whole restart hull, add mobile
chains around it" pattern below straightforward.

Typical workflow: a frozen template from a restart
==================================================

The most common use is to build or capture a structure, save it to a restart file,
and then re-run with that structure frozen while new chains move around it:

#. Run (or construct) the configuration you want to keep fixed and save a
   ``restart.pimms``.
#. List the chainIDs that make up that structure in a freeze file (all of them, for
   a fully fixed template), or ``B`` lines for the beads to hold where part of a
   chain should stay free.
#. Start a new simulation with ``RESTART_FILE`` + ``FREEZE_FILE``, optionally adding
   mobile chains with ``EXTRA_CHAIN`` (see :doc:`reference_controls`).

For a structure converted from a PDB or mmCIF file with
:ref:`pdb2pimms <pdb2pimms>`, the first two steps are done for you: it writes a
restart file and a freeze file listing every chain in it, and with ``--mobile``
(for example ``--mobile A:1-37``) it leaves the residues you name free and writes
the ``B`` lines for the rest of their chains. With ``--protein-mode`` it asks
metapredict which residues are disordered, freezes the folded domains, leaves the
disordered regions free and writes a keyfile that runs the result (see
:ref:`pdb2pimms`).

.. code-block:: text

   DIMENSIONS            : 60 60 60      # still required; must be compatible with the restart
   RESTART_FILE          : template.restart
   FREEZE_FILE           : freeze.txt
   EXTRA_CHAIN           : 150 AB        # mobile chains added around the frozen template

The ``star_destroyer`` demo in ``demo_keyfiles/`` is exactly this pattern: a 170-chain
Star Destroyer hull loaded from a restart file and frozen in place, with 150 small
mobile chains added by ``EXTRA_CHAIN`` and ``PARALLELIZE : True`` on top.

Two more demos there use freezing on a larger scale. ``nuclear_pore_complex``
freezes a nuclear envelope and a pore scaffold as six whole chains (``C`` lines)
and holds 88 FG-Nups by their first bead alone (``B`` lines), so the rest of each FG
chain moves. ``star_destroyer_fleet`` freezes a 43,725-bead ship with a single
``C 1`` line, with fighters flying round it as rigid bodies (see above).

Works with parallelization
==========================

Freezing composes with :doc:`parallelization`. The parallel checkerboard kernels
exclude frozen beads from the movable set (via an explicit per-bead frozen mask)
but keep them in place as fixed, energy-contributing obstacles, so ``FREEZE_FILE``
and ``PARALLELIZE`` can be used together - the frozen scaffold is respected while
the mobile moves are threaded. The parallel crankshaft reads the mask bead by bead,
so it moves the free beads of a partly frozen chain exactly as the serial one does.

Chains with frozen beads also take no part in the length partition that splits the
whole-chain slither and pull moves between the parallel and serial kernels: neither
kernel is ever asked to move them, so such a chain too long to fit a block interior
is simply an obstacle to both passes. The start-up parallelization report counts
them separately (``Frozen chains: N (M beads) - excluded from every parallel move,
kept as fixed obstacles``, or, when some chains are only partly frozen, ``Frozen
beads: M, in N chain(s) (...)``).
