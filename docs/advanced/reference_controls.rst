.. _advanced-other:

==============================================
Reference ensembles & other advanced controls
==============================================

This page collects the remaining advanced keywords: reference/debugging switches,
box and trajectory controls, a couple of chain-handling options, and the
experimental gate that unlocks the not-yet-stable features.

Reference ensembles
===================

Sometimes the most useful comparison is the system with its interactions turned off
- a well-defined reference state that isolates the effect of excluded volume or of
the backbone.

``NON_INTERACTING : True``
   Zero **all** interaction energies and run a pure excluded-volume reference
   ensemble. The parameter file is still required (bead types must be defined) but
   its pairwise and solvation energies are ignored. This is the natural "ideal
   chain in a box" baseline to compare an interacting run against - anything that
   differs from the non-interacting ensemble is a genuine consequence of the
   interactions.

``ANGLES_OFF : True``
   Disable the backbone-angle penalties entirely, so the chains are perfectly
   flexible. With this set you do **not** need ``ANGLE_PENALTY`` lines in the
   parameter file. Useful for isolating the role of chain stiffness, or simply for
   models where stiffness is not wanted.

The two switches are independent and can be combined: ``NON_INTERACTING : True`` with
``ANGLES_OFF : True`` is the fully ideal, freely-jointed, excluded-volume-only chain.

Energy-consistency checking
===========================

PIMMS tracks the total energy **incrementally** - each accepted move adds its energy
change to a running total rather than recomputing the whole Hamiltonian. That is
what makes it fast, but it also means a subtle bookkeeping bug would slowly drift the
tracked energy away from the truth.

``ENERGY_CHECK : <freq>``
   Every ``<freq>`` steps, recompute the total energy from scratch and compare it to
   the incrementally tracked value. If they disagree the run aborts with an exception
   and dumps the offending configuration to ``CONFIG_AT_ENERGY_FAIL.pdb`` /
   ``.xtc`` so you can inspect it. This is an **O(N)** check: cheap to run
   occasionally on a modest system, but expensive if you run it every step on a large
   one. The default frequency is deliberately infrequent; lower it when you are
   developing or debugging and want tight verification, raise it (or leave it) for
   production.

Box and equilibration controls
==============================

``RESIZED_EQUILIBRATION`` / ``EQUILIBRATION_OFFSET``
   Equilibrate in a **smaller** box and then grow to the production ``DIMENSIONS`` at
   the end of equilibration (chains are re-centred as the box expands). This lets you
   condense or assemble a system at high effective concentration and then relax it
   into the full production volume - often much faster than waiting for assembly to
   happen at the production density.

   ``RESIZED_EQUILIBRATION`` gives the equilibration box size (2 or 3 values), and
   ``EQUILIBRATION_OFFSET`` places that smaller box within the full box. Both are
   constrained: ``RESIZED_EQUILIBRATION`` must be ``<= DIMENSIONS`` in every
   dimension, and ``EQUILIBRATION_OFFSET + RESIZED_EQUILIBRATION`` must fit inside
   ``DIMENSIONS``. The equilibration phase is always run under **hardwall**
   boundaries (forced internally, so a system is never resized while chains straddle
   a periodic face); your production ``HARDWALL`` setting takes over once the box has
   grown. The feature is incompatible with ``RESTART_OVERRIDE_DIMENSIONS`` and with
   PBC restart files. See :ref:`overview-setup` for the surrounding set-up keywords.

   .. code-block:: text

      DIMENSIONS            : 60 60 60
      RESIZED_EQUILIBRATION : 30 30 30      # equilibrate at 8x the density...
      EQUILIBRATION_OFFSET  : 15 15 15      # ...centred in the production box

``AUTOCENTER : True``
   For a **single-chain** simulation, keep the chain centred in the box on every
   frame. This stops a lone chain drifting to the edge of a hardwall box or
   wandering across periodic images, which keeps trajectories tidy for
   visualisation and analysis.

Chain-handling options
======================

``CASE_INSENSITIVE_CHAINS : False``
   By default (``True``) chain sequences are upper-cased when the keyfile is read, so
   ``a`` and ``A`` are the same bead type. Set this to ``False`` to treat case as
   significant, which effectively doubles the alphabet of bead types available
   (``A`` and ``a`` become distinct). Every bead letter used in a chain must still be
   defined in the parameter file.

``EXTRA_CHAIN : <count> <sequence>``
   Add chains that were **not** in the original ``RESTART_FILE`` when restarting from
   a saved configuration. The format matches the ``CHAIN`` keyword
   (``<number of chains> <sequence>``), multiple ``EXTRA_CHAIN`` lines are allowed for
   different species, and the new chains are inserted at random so as not to overlap
   anything already present. This is how you take the end-state of one run and
   continue it with additional material - and it can be repeated as many times as you
   like. ``EXTRA_CHAIN`` requires a ``RESTART_FILE`` (there must be an existing
   configuration to add to). It pairs naturally with a :doc:`freeze file
   <freeze>`: freeze the restart configuration and let the extra chains explore
   around it (as in the ``star_destroyer`` demo).

The experimental gate
=====================

``EXPERIMENTAL_FEATURES : True``
   Unlocks the remaining not-yet-stable keywords and moves. As of this writing the
   only gated feature is the ``MOVE_VMMC`` collective move (and its ``VMMC_*`` tuning
   keywords). It is not guaranteed to behave correctly in every configuration, so the
   recommendation is to leave the gate ``False`` unless you specifically need it - and
   to sanity-check the results carefully when you do. (The pull, jump-and-relax and
   temperature-switch (TSMMC) moves, non-cubic boxes, and the ``EXTRA_CHAIN``,
   ``FREEZE_FILE`` and ``EQUILIBRATION_OFFSET`` keywords have all graduated out of the
   experimental gate and need no special flag.)
