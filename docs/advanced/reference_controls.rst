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
   Zero **all** bead-bead (short-range, long-range and super-long-range) and
   solvation energies and run a pure excluded-volume reference ensemble. The
   parameter file is still required (bead types must be defined) but its pairwise
   and solvation energies are ignored. Backbone-angle penalties are **not**
   affected - see ``ANGLES_OFF`` below. The start-up summary says so (``NOTE: This
   is a non-interacting simulation ...``), and unless ``REDUCED_PRINTING`` is on a
   warning is printed for every bead pair whose parameter-file energy is being
   overridden. This is the natural "ideal chain in a box" baseline to compare an
   interacting run against; anything that differs from the non-interacting ensemble
   is a genuine consequence of the interactions.

``ANGLES_OFF : True``
   Disable the backbone-angle penalties entirely, so the chains are perfectly
   flexible (announced at start-up unless ``REDUCED_PRINTING`` is on). With this set
   you do **not** need ``ANGLE_PENALTY`` lines in the parameter file; without it,
   every bead type that has interaction energies needs one, and a missing line is a
   parameter-file error. Useful for isolating the role of chain stiffness, or simply
   for models where stiffness is not wanted.

Both default to ``False``. The two switches are independent and can be combined:
``NON_INTERACTING : True`` with ``ANGLES_OFF : True`` is the fully ideal,
freely-jointed, excluded-volume-only chain.

Energy-consistency checking
===========================

PIMMS tracks the total energy **incrementally** - each accepted move adds its energy
change to a running total rather than recomputing the whole Hamiltonian. That is
what makes it fast, but it also means a subtle bookkeeping bug would slowly drift the
tracked energy away from the truth.

``ENERGY_CHECK : <freq>`` (default ``20000``)
   Every ``<freq>`` steps, recompute the total energy from scratch and compare it to
   the incrementally tracked value. The same pass also cross-checks the occupancy
   and type grids against the chains' own positions and sequences, so a corrupted
   type grid (which would make the tracked and recomputed energies wrong in exactly
   the same way, and therefore invisible to the energy comparison alone) is caught
   too. Either kind of disagreement aborts the run with a
   ``SimulationEnergyException``, and both keep what can be kept: the trajectory
   is closed cleanly up to its last frame (under ``SAVE_AT_END`` the frames
   buffered in memory are written out first rather than lost), and the offending
   configuration is dumped to ``CONFIG_AT_ENERGY_FAIL.pdb`` / ``.xtc`` so you can
   inspect it. A grid inconsistency also reports the problems it found (the first
   ten of them) to stdout and the log. This is an **O(N)** check: cheap to run
   occasionally on a modest system, but expensive if you run it every step on a large
   one. Each check also prints the full energy decomposition (short-range,
   long-range, super-long-range and angle terms) to stdout, which is a convenient
   way to see where the energy of a configuration actually comes from. The default
   frequency is deliberately infrequent; lower it when you are
   developing or debugging and want tight verification, raise it for
   production, or set it to ``0`` to switch it off entirely.

Box and equilibration controls
==============================

``RESIZED_EQUILIBRATION`` / ``EQUILIBRATION_OFFSET``
   Equilibrate in a **smaller** box and then grow to the production ``DIMENSIONS`` at
   the end of equilibration. This lets you
   condense or assemble a system at high effective concentration and then relax it
   into the full production volume - often much faster than waiting for assembly to
   happen at the production density.

   ``RESIZED_EQUILIBRATION`` gives the equilibration box size (2 or 3 values,
   matching the length of ``DIMENSIONS``). Without ``EQUILIBRATION_OFFSET`` the
   small box is placed at the **centre** of the larger one as it grows: every bead
   is shifted by half the difference in box size, so the configuration keeps its
   place inside the small box (it is not re-centred on its own centre of mass).
   With it, every bead is shifted by that explicit per-dimension offset instead,
   which is how you place the small box somewhere other than the centre. Both are
   constrained:

   * ``RESIZED_EQUILIBRATION`` must be ``<= DIMENSIONS`` in every dimension (the
     box can only grow), and every axis must be ``>= 7`` - the equilibration box is
     simulated, so it obeys the same floor as ``DIMENSIONS``.
   * ``EQUILIBRATION_OFFSET`` values must be ``>= 0``, must have the same number of
     entries, and ``EQUILIBRATION_OFFSET + RESIZED_EQUILIBRATION`` must fit inside
     ``DIMENSIONS``. ``EQUILIBRATION_OFFSET`` on its own, without
     ``RESIZED_EQUILIBRATION``, is an error.
   * ``EQUILIBRATION : 0`` makes the feature meaningless, so PIMMS prints a warning
     and switches ``RESIZED_EQUILIBRATION`` off (and ``EQUILIBRATION_OFFSET`` with
     it, if given) rather than resizing at step zero.

   The equilibration phase is always run under **hardwall**
   boundaries (forced internally, so a system is never resized while chains straddle
   a periodic face); your production ``HARDWALL`` setting takes over once the box has
   grown. The box swap happens *after* the step-``EQUILIBRATION`` move, keeping the
   convention that ``EQUILIBRATION`` is the last equilibration step. If
   ``PARALLELIZE`` is on, the parallelization report is re-issued at that point,
   since the block decomposition depends on the box. While the small box is in use
   the trajectory goes to ``eq_START.pdb`` / ``eq_traj.xtc`` (written only when
   ``SAVE_EQ`` is on, the default); the production ``START.pdb`` / ``traj.xtc`` pair is
   opened at the resize, and any copy of it left in the directory by an earlier run
   is deleted at start-up. The feature is incompatible with
   ``RESTART_OVERRIDE_DIMENSIONS``, with ``RESTART_CONTINUE`` and with PBC restart
   files; with a hardwall restart file the restart box must be no larger than the
   ``RESIZED_EQUILIBRATION`` box on any axis.
   See :ref:`overview-setup` for the surrounding set-up keywords.

   .. code-block:: text

      DIMENSIONS            : 60 60 60
      RESIZED_EQUILIBRATION : 30 30 30      # equilibrate at 8x the density...
      EQUILIBRATION_OFFSET  : 15 15 15      # ...centred in the production box

``AUTOCENTER : True`` (default ``False``)
   For a **single-chain** simulation, re-centre the chain in the box in every
   *written* trajectory/PDB frame. The simulation itself is untouched - the
   chain still explores the box and feels any hardwall normally; the centring
   only removes drift from the output, keeping trajectories tidy for
   visualisation and analysis. The chain is also made into a single periodic
   image first, so ``AUTOCENTER`` takes precedence over
   ``TRAJECTORY_PBC_UNWRAP``. Silently ignored when more than one chain is
   present.

Chain-handling options
======================

``CASE_INSENSITIVE_CHAINS : False``
   By default (``True``) chain sequences are upper-cased when the keyfile is read, so
   ``a`` and ``A`` are the same bead type. Set this to ``False`` to treat case as
   significant, which effectively doubles the alphabet of bead types available
   (``A`` and ``a`` become distinct). Only ``CHAIN`` and ``EXTRA_CHAIN`` sequences
   are case-folded - the parameter file never is - and every bead letter used in a
   chain must still be defined there.

``EXTRA_CHAIN : <count> <sequence>``
   Add chains that were **not** in the original ``RESTART_FILE`` when restarting from
   a saved configuration. The format matches the ``CHAIN`` keyword
   (``<number of chains> <sequence>``, with ``<count> >= 1``), it is one of the three
   keywords that may appear on multiple lines (the others being ``CHAIN`` and
   ``ANA_RESIDUE_PAIRS``) so different species can be added together, and the new
   chains are inserted at random so as not to overlap
   anything already present. They are given fresh chainIDs *after* the restart
   chains, so the restart chains keep the IDs they had, and an extra chain whose
   sequence already appears in the restart file joins that existing chain type
   rather than defining a new one. This is how you take the
   end-state of one run and
   continue it with additional material - and it can be repeated as many times as you
   like. ``EXTRA_CHAIN`` requires a ``RESTART_FILE`` (there must be an existing
   configuration to add to) and is rejected without one; it is also refused with
   ``RESTART_CONTINUE``, since adding chains makes it a different system. It pairs
   naturally with a
   :doc:`freeze file <freeze>`: freeze the restart configuration and let the extra
   chains explore around it (as in the ``star_destroyer`` demo).

The experimental gate
=====================

``EXPERIMENTAL_FEATURES : True``
   Unlocks the remaining not-yet-stable keywords and moves. As of this writing the
   gated keywords are exactly ``MOVE_VMMC``, ``VMMC_MAX_DISPLACEMENT`` and
   ``VMMC_MAX_CLUSTER`` - the collective VMMC move and its tuning parameters. The gate
   fires only when one of them is set **away from its default**, so simply writing
   ``MOVE_VMMC : 0.0`` does not require it. VMMC is not guaranteed to behave correctly
   in every configuration, so the
   recommendation is to leave the gate ``False`` unless you specifically need it - and
   to sanity-check the results carefully when you do. (The pull, jump-and-relax and
   temperature-switch (TSMMC) moves, non-cubic boxes, and the ``EXTRA_CHAIN``,
   ``FREEZE_FILE`` and ``EQUILIBRATION_OFFSET`` keywords have all graduated out of the
   experimental gate and need no special flag.)
