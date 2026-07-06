.. _advanced:

=================
Advanced features
=================

Beyond the core move set and the standard analysis, PIMMS has a handful of
features aimed at harder sampling problems and specialised workflows: ramping the
temperature during a run, pinning part of the system in place, spreading the hot
moves across CPU cores, and enhanced-sampling temperature excursions. Each has its
own page, listed under **Advanced Features** in the sidebar and summarised below.

One feature is still **experimental** and requires ``EXPERIMENTAL_FEATURES : True``
in the keyfile (the VMMC move); it is flagged clearly on the page where it appears.
Experimental features are not guaranteed to behave correctly in every configuration,
so leave the gate off unless you specifically need it (see :doc:`reference_controls`).

At a glance
===========

.. list-table::
   :header-rows: 1
   :widths: 24 24 52

   * - Feature
     - Key keyword(s)
     - What it is for
   * - :doc:`quench`
     - ``QUENCH_RUN``
     - Ramp the temperature during a run - simulated annealing (cooling) to
       drive assembly, or heating to melt a structure.
   * - :doc:`freeze`
     - ``FREEZE_FILE``
     - Hold chosen chains rigidly fixed as a scaffold, surface or template that
       the mobile chains still interact with.
   * - :doc:`parallelization`
     - ``PARALLELIZE`` / ``PARALLEL_THREADS``
     - Run the crankshaft, slither and pull moves on multi-threaded
       "checkerboard" kernels for large systems.
   * - :doc:`tsmmc`
     - ``MOVE_CTSMMC`` / ``MOVE_MULTICHAIN_TSMMC`` / ``MOVE_SYSTEM_TSMMC``
     - Temperature-excursion (tempered-transitions) moves that help the system
       hop over energy barriers ordinary moves cannot cross.
   * - :doc:`custom_analysis`
     - ``ANALYSIS_MODULE`` / ``ANA_CUSTOM``
     - Run your own Python analysis code against the live lattice during the run.
   * - :doc:`reference_controls`
     - ``NON_INTERACTING``, ``RESIZED_EQUILIBRATION``, ``EXTRA_CHAIN`` ...
     - Reference ensembles, energy-consistency checks, box/trajectory controls
       and the experimental gate.

None of these features change what PIMMS is sampling unless you ask them to: a
freeze file and ``PARALLELIZE`` leave the target Boltzmann distribution untouched,
a quench deliberately changes the temperature along a schedule you specify, and the
TSMMC moves are constructed to preserve detailed balance. Where a feature *does*
alter the physics (a quench, or ``NON_INTERACTING``) that is the whole point of it,
and the page says so explicitly.
