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
freeze. The simulation aborts at startup if the file does not exist.

.. code-block:: text

   FREEZE_FILE : freeze.txt

The freeze file itself is a short list of ``C`` directives, one or more per file.
Lines beginning with ``#`` are comments:

.. code-block:: text

   # freeze.txt  -  lines beginning with # are comments
   C 1 2 3          # freeze chainIDs 1, 2 and 3 (chainIDs are numbered from 1)
   C 10 11 12 13    # more C lines are allowed; IDs may be split across lines

Each ``C`` line contributes its integer chainIDs to the frozen set; the order and
the split across lines do not matter.

What "frozen" means
===================

A frozen chain is simply **excluded from the pool of chains PIMMS can pick to
move** - the bead selector skips it. Everything else about it is unchanged:

* It stays exactly where it was placed (from the ``CHAIN`` set-up or, more usually,
  from a ``RESTART_FILE``).
* It still **excludes volume** - mobile beads cannot overlap it.
* It still **contributes to the energy** - every interaction between a frozen bead
  and a mobile bead is counted normally, so the mobile chains are attracted to or
  repelled by the frozen scaffold just as they would be by a mobile partner.

Because frozen chains never move, freezing costs nothing per step beyond keeping
them on the lattice; if anything it speeds the run up slightly, since there are
fewer movable chains to choose from.

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

   RESTART_FILE          : template.restart
   FREEZE_FILE           : freeze.txt
   EXTRA_CHAIN           : 150 AB        # mobile chains added around the frozen template

The ``star_destroyer`` demo in ``demo_keyfiles/`` is exactly this pattern: a large
multi-chain hull loaded from a restart file and frozen in place, with a swarm of
small mobile chains added by ``EXTRA_CHAIN``.

Works with parallelization
==========================

Freezing composes with :doc:`parallelization`. The parallel checkerboard kernels
exclude frozen beads from the movable set but keep them in place as fixed,
energy-contributing obstacles, so ``FREEZE_FILE`` and ``PARALLELIZE`` can be used
together - the frozen scaffold is respected while the mobile moves are threaded.

.. note::

   Freezing is currently at **whole-chain** granularity: a chain is either entirely
   frozen or entirely free. A per-bead freeze directive (a ``B`` line) is reserved
   in the file format but is **not yet implemented**.
