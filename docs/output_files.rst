.. _output-files:

============
Output files
============

PIMMS writes its results as plain-text ``.dat`` tables plus a molecular trajectory
(``.pdb`` topology + ``.xtc`` frames), all in the directory the simulation is run
from. Which files appear depends on the analysis keywords you enable; the
*frequency* of each is controlled by its ``ANA_*`` keyword (falling back to
``ANALYSIS_FREQ``). This page catalogues every file and then explains how to load
them.

.. note::

   ``.dat`` files are appended each analysis step, so a column-per-quantity,
   row-per-step layout is typical. A handful of analyses (internal scaling,
   distance maps) are **accumulated over the run and written once at the end**;
   these are noted below.

.. _output-catalogue:

The output files
================

Core state & performance
------------------------

``ENERGY.dat``
    Instantaneous potential energy. Two tab-separated columns: ``step``,
    ``energy``. Frequency: ``EN_FREQ``. Always written.

``PERFORMANCE.dat``
    Throughput and timing. Columns: ``step``, ``E``/``P`` (equilibration vs
    production), loop-steps per second, overall MC moves per second (counting all
    sub-loop moves), elapsed time and estimated remaining time (``hh:mm:ss``). A
    header line is written when the file is created.

``QUENCH.dat``
    Only written for quench runs (``QUENCH_RUN : True``). Columns: ``step``,
    ``temperature``, ``energy`` - the temperature ramp and the energy response.
    See :doc:`advanced/quench`.

Move bookkeeping
----------------

These three are written together by ``ANA_ACCEPTANCE``. The per-move columns are
indexed by **move code** (1 = crankshaft, 2 = chain-translate, 3 = chain-rotate,
4 = chain-pivot, 5 = head-pivot, 6 = slither, 7 = cluster-translate,
8 = cluster-rotate, 9 = CTSMMC, 10 = multichain-TSMMC, 11 = pull,
12 = system-TSMMC, 13 = jump-and-relax, 14 = VMMC).

``MOVE_FREQS.dat``
    Number of moves **attempted** of each type. Columns: ``step`` then one count
    per move code.

``ACCEPTANCE.dat``
    Number of moves **accepted** of each type (same layout). Divide by
    ``MOVE_FREQS.dat`` to obtain per-move acceptance ratios.

``TOTAL_MOVES.dat``
    Cumulative total number of accept/reject operations across all sub-loops
    (``step``, ``total``) - the "true" amount of MC work done.

Single-chain (polymeric) analysis
----------------------------------

``RG.dat`` / ``ASPH.dat``
    Per-chain radius of gyration / asphericity. Each row is a ``step`` followed by
    one value per chain. Trigger: ``ANA_POL``.

``END_TO_END_DIST.dat``
    Per-chain end-to-end distance (a ``step`` column followed by one value per
    chain). Written at the same frequency as ``RG.dat``/``ASPH.dat`` (``ANA_POL``).

``RES_TO_RES_DIST.dat``
    Distance between a chosen residue pair, for every chain. Columns: ``step``,
    the two residue indices, then one distance per chain. Trigger:
    ``ANA_INTER_RESIDUE`` with ``ANA_RESIDUE_PAIRS`` set.

``INTSCAL.dat`` / ``INTSCAL_SQUARED.dat``
    Mean internal scaling ``R(|i-j|)`` (and its square) versus sequence
    separation. Columns: ``gap``, ``mean``. The data is accumulated over the run at
    the ``ANA_INTSCAL`` sampling frequency and the file is written once at the end
    of every run.

``SCALING_INFORMATION.dat``
    Fitted polymer-scaling parameters: the apparent scaling exponent ``nu`` and
    prefactor ``R0`` from ``R = R0 · N^nu``. Written once at the end of every run.

``DISTANCE_MAP.dat``
    Mean inter-residue distance map - a ``seqlen × seqlen`` matrix (tab-separated
    rows), accumulated at the ``ANA_DISTMAP`` sampling frequency and written once at
    the end of every run.

For multi-component systems the internal-scaling/distance-map files are also
written per chain type as ``CHAIN_<TYPE>_INTSCAL.dat`` etc.

Cluster analysis
----------------

Enabled by ``ANA_CLUSTER``. PIMMS identifies **short-range clusters** (chains in
direct contact) and **long-range clusters** (chains connected via any
interaction), and reports both their size distributions and per-cluster shape
descriptors. ``ANA_CLUSTER_THRESHOLD`` sets the minimum chain count for a
connected component to be counted as a cluster.

``CLUSTERS.dat`` / ``NUM_CLUSTERS.dat``
    Per-step cluster size distribution (``CLUSTERS.dat``: comma-separated cluster
    sizes) and the number of clusters (``NUM_CLUSTERS.dat``: tab-separated ``step``,
    ``count``).

``CLUSTER_RG.dat`` / ``CLUSTER_ASPH.dat`` / ``CLUSTER_AREA.dat`` / ``CLUSTER_VOL.dat`` / ``CLUSTER_DEN.dat``
    Per-cluster radius of gyration, asphericity, surface area, volume and density
    (one value per cluster per step).

``CLUSTER_RADIAL_DENSITY_PROFILE.dat``
    Radial density profile (density vs distance from the cluster centre of mass),
    computed for sufficiently large clusters.

``LR_CLUSTERS.dat``, ``NUM_LR_CLUSTERS.dat``, ``LR_CLUSTER_RG.dat`` (etc.)
    The same set of files for the **long-range** clusters.

For multi-component systems, ``CHAIN_<TYPE>_CLUSTERS.dat`` (and the long-range
``CHAIN_<TYPE>_LR_CLUSTERS.dat``) records the fraction of each cluster contributed
by that chain type.

Trajectory
----------

``START.pdb``
    Topology file - one ``ATOM`` record per bead, chains labelled by type. Used as
    the topology when loading the trajectory. Its ``CRYST1`` record gives the
    **periodic unit cell**, so an axis of ``L`` lattice sites is written as
    ``L * LATTICE_TO_ANGSTROMS`` angstroms (sites ``L-1`` and ``0`` are periodic
    neighbours one lattice unit apart) - the same box that is written into every XTC
    frame. For a 2D system the ``c`` axis is one lattice unit, since there is no
    periodicity in z.

``traj.xtc``
    The trajectory itself (XTC format, via ``mdtraj``), one frame every
    ``XTC_FREQ`` steps. Equilibration frames are included only if
    ``SAVE_EQ : True``; with ``SAVE_AT_END : True`` the trajectory is buffered in
    memory and written once at the end. Coordinates are scaled by
    ``LATTICE_TO_ANGSTROMS``. (When ``RESIZED_EQUILIBRATION`` is used the
    equilibration phase is written separately as ``eq_START.pdb`` / ``eq_traj.xtc``.)

    By default the raw on-lattice positions are written, so under periodic
    boundaries a chain that crosses a box face appears split across the two faces.
    Set ``TRAJECTORY_PBC_UNWRAP : True`` to make every chain **whole** before each
    frame (and the ``START.pdb`` topology) is written - each chain is shifted into a
    single periodic image, so no chain is torn across a boundary. This is purely a
    visualisation convenience (it does not affect the simulation), and unwrapped
    coordinates may fall outside the box; it has no effect under ``HARDWALL``, where
    chains never cross a boundary.

Echoed inputs & checkpoint
--------------------------

``parameters_used.prm``
    A copy of the parameter file actually used (with a timestamp header), so a run
    is self-documenting.

``absolute_energies_of_angles.txt``
    Human-readable summary of the angle penalties per residue.

``chain_to_chainid.txt``
    Mapping of each chainID to its length and sequence. Written only if
    ``WRITE_CHAIN_TO_CHAINID : True``; handy for choosing chains to put in a
    :ref:`freeze file <advanced-freeze>`.

``log.txt``
    A plain-text run log written on every simulation: the startup banner, the
    resolved configuration and progress/status messages. Handy for reconstructing
    exactly how a run was set up and whether it finished cleanly.

``restart.pimms``
    Configuration checkpoint - see :doc:`restart_files`.

.. _output-analyze:

Analysing the output
====================

Plain-text ``.dat`` files
-------------------------

Most ``.dat`` files are tab-separated; the exceptions are the cluster
size-distribution and per-cluster property files (``CLUSTERS.dat``,
``CLUSTER_RG.dat`` and friends), which are comma-separated. Any tool reads them;
with NumPy:

.. code-block:: python

   import numpy as np

   # energy trace: columns [step, energy]
   step, energy = np.loadtxt("ENERGY.dat", delimiter="\t", unpack=True)
   print("mean production energy:", energy[len(energy)//2:].mean())

   # per-move acceptance ratio
   attempted = np.loadtxt("MOVE_FREQS.dat", delimiter="\t")
   accepted  = np.loadtxt("ACCEPTANCE.dat", delimiter="\t")
   ratio = accepted[:, 1:] / np.clip(attempted[:, 1:], 1, None)   # column k = move code k

The square ``DISTANCE_MAP.dat`` matrix loads directly with
``np.loadtxt("DISTANCE_MAP.dat")`` and can be shown with ``imshow``.

Trajectories (``.pdb`` + ``.xtc``)
----------------------------------

Load the ``.xtc`` frames with the ``START.pdb`` topology using `mdtraj
<https://www.mdtraj.org/>`_:

.. code-block:: python

   import mdtraj as md

   traj = md.load("traj.xtc", top="START.pdb")
   print(traj)                       # n_frames, n_atoms
   # mdtraj coordinates are in nanometres; one lattice unit = LATTICE_TO_ANGSTROMS A,
   # so: lattice_units = traj.xyz * 10.0 / LATTICE_TO_ANGSTROMS
   rg = md.compute_rg(traj)          # radius of gyration per frame, etc.

The same ``START.pdb`` + ``traj.xtc`` pair loads directly in **VMD** (and other
molecular viewers) for visualisation. For 2D simulations the out-of-plane
coordinate is held at zero. If you used ``AUTOCENTER : True`` (single-chain runs)
the chain is already centred each frame, so no alignment is needed before
analysis.
