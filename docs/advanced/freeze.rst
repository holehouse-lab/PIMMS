.. _advanced-freeze:

============
Freeze files
============

A **freeze file** holds chosen chains rigidly fixed for the entire simulation. The
frozen chains never move, but they are still fully present on the lattice: they
occupy their sites and contribute to the energy, so the mobile chains feel them
exactly as they would any other chain. This is the tool for building a fixed
scaffold, a wall or surface, or a pre-formed template that the rest of the system
explores around.

Turning it on
=============

Point the ``FREEZE_FILE`` keyword at a plain-text file that lists the chainIDs to
freeze. This is the only whole-chain immobilization mechanism. The simulation
aborts at startup if the file does not exist, or if ``FREEZE_FILE`` is given with
an empty value.

.. code-block:: text

   FREEZE_FILE : freeze.txt

The freeze file itself is a short list of ``C`` directives, one or more per file.
Lines beginning with ``#`` are comments, and a ``#`` anywhere on a line comments
out the rest of it:

.. code-block:: text

   # freeze.txt  -  lines beginning with # are comments
   C 1 2 3          # freeze chainIDs 1, 2 and 3 (chainIDs are numbered from 1)
   C 10 11 12 13    # more C lines are allowed; IDs may be split across lines

Each ``C`` line contributes its integer chainIDs to the frozen set; the order and
the split across lines do not matter. Blank lines are skipped, duplicate IDs are
removed and the resulting set is processed in numeric order. Every directive must
contain at least one ID; an empty ``C``/``B`` line, a non-integer ID, or an unknown
directive (including a lower-case ``c``) is a parse-time error rather than a
silently ignored instruction. A chainID that does not exist in the system is also
an error, raised once the lattice has been built and naming both the requested and
the available IDs. A file that contains no directives at all (only comments or blank
lines) is accepted and freezes nothing; the start-up summary then shows
``Chains to freeze : []``.

The one thing a freeze file cannot do is freeze everything: if it names *every*
chain in the system, PIMMS refuses to start rather than write out ``N_STEPS`` copies
of an unchanging configuration.

The frozen set is echoed in the start-up summary (here for a file freezing chains
1, 2 and 3) and written to the log::

   --> Freeze File Settings

   Freeze file      : freeze.txt
   Chains to freeze : [1, 2, 3]

What "frozen" means
===================

A frozen chain is **excluded from the pool of chains PIMMS can pick to move**, at
every layer:

* The outer-loop chain selector never returns a frozen chain, so none of the
  single-chain moves (translate, rotate, pivot, head pivot, jump-and-relax, chain
  TSMMC) can ever be proposed for one.
* The whole-system megamoves build their bead/chain selectors from the mobile
  chains only, so the crankshaft, slither and pull never propose a frozen bead -
  in the parallel kernels this is enforced by an explicit per-bead frozen mask.
* Any *collective* move whose cluster would include a frozen chain (cluster
  translate, cluster rotate, VMMC seeding and VMMC recruitment) is rejected
  outright, so a frozen chain is never dragged along by its neighbours.
* The multichain TSMMC move draws its random subset from the mobile chains only,
  and a system-wide TSMMC excursion is made of ordinary sub-moves, each of which
  honours the frozen set as above.

Mobile chains bound to a frozen scaffold still move by every move except the
collective ones that would have to drag the scaffold along. Everything else about a
frozen chain is unchanged:

* It stays exactly where it was placed (from the ``CHAIN`` set-up or, more usually,
  from a ``RESTART_FILE``). The one exception is a ``RESIZED_EQUILIBRATION`` run:
  when the box grows, the whole configuration - frozen chains included - is
  translated rigidly into the production box (by half the difference in box size,
  which places the small box at the centre of the large one, or by
  ``EQUILIBRATION_OFFSET``), so the scaffold keeps its shape and its position
  relative to the mobile chains but not its absolute coordinates.
* It still **excludes volume** - mobile beads cannot overlap it.
* It still **contributes to the energy** - every interaction between a frozen bead
  and a mobile bead is counted normally, so the mobile chains are attracted to or
  repelled by the frozen scaffold just as they would be by a mobile partner.

Freezing adds essentially nothing to the cost of a step: the frozen beads simply
sit in the grid. It does change how the megamoves' work is shared out: the
crankshaft's ``CRANKSHAFT_SUBSTEPS`` attempts are spread over the mobile beads only,
and slither and pull make their per-chain ``SLITHER_SUBSTEPS`` / ``PULL_SUBSTEPS``
attempts for the mobile chains only.

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
   a fully fixed template).
#. Start a new simulation with ``RESTART_FILE`` + ``FREEZE_FILE``, optionally adding
   mobile chains with ``EXTRA_CHAIN`` (see :doc:`reference_controls`).

.. code-block:: text

   DIMENSIONS            : 60 60 60      # still required; must be compatible with the restart
   RESTART_FILE          : template.restart
   FREEZE_FILE           : freeze.txt
   EXTRA_CHAIN           : 150 AB        # mobile chains added around the frozen template

The ``star_destroyer`` demo in ``demo_keyfiles/`` is exactly this pattern: a 170-chain
Star Destroyer hull loaded from a restart file and frozen in place, with 150 small
mobile chains added by ``EXTRA_CHAIN`` and ``PARALLELIZE : True`` on top.

Works with parallelization
==========================

Freezing composes with :doc:`parallelization`. The parallel checkerboard kernels
exclude frozen beads from the movable set (via an explicit per-bead frozen mask)
but keep them in place as fixed, energy-contributing obstacles, so ``FREEZE_FILE``
and ``PARALLELIZE`` can be used together - the frozen scaffold is respected while
the mobile moves are threaded.

Frozen chains also take no part in the length partition that splits the whole-chain
slither and pull moves between the parallel and serial kernels: neither kernel is
ever asked to move them, so a frozen chain too long to fit a block interior is
simply an obstacle to both passes. The start-up parallelization report
counts them separately (``Frozen chains: N (M beads) - excluded from every parallel
move, kept as fixed obstacles``).

.. note::

   Freezing is currently at **whole-chain** granularity: a chain is either entirely
   frozen or entirely free. A per-bead freeze directive (a ``B`` line) is reserved
   in the file format but is **not yet implemented**: a syntactically valid ``B``
   line is parsed and then raises ``UnfinishedCodeException``, so it fails loudly
   rather than being silently ignored.
