.. _output-files:

============
Output files
============

PIMMS writes its results as plain-text ``.dat`` tables plus a molecular trajectory
(``.pdb`` topology + ``.xtc`` frames), all in the directory the simulation is run
from. The *frequency* of each analysis is controlled by its ``ANA_*`` keyword
(falling back to ``ANALYSIS_FREQ``); setting any ``ANA_*`` frequency (or
``ENERGY_CHECK``) to 0 disables that analysis for the run, and setting
``ANALYSIS_FREQ`` to 0 disables every analysis whose own ``ANA_*`` keyword was
left unset. This page catalogues every file and then explains how to load them.

**Units.** Every length in a ``.dat`` file (radii of gyration, end-to-end and
residue-residue distances, internal scaling, distance maps) is in **lattice
units**, and areas, volumes and densities are in the matching powers of the
lattice unit; multiply a length by ``LATTICE_TO_ANGSTROMS`` to get angstroms.
Energies are in the parameter file's units (PIMMS works in reduced units with
k = 1, so ``TEMPERATURE`` is on the same scale). The trajectory files
(``START.pdb`` / ``traj.xtc`` and their ``eq_`` and ``CONFIG_AT_ENERGY_FAIL``
counterparts) are the only outputs written in real-space units (angstroms in the
``.pdb`` files, nanometres in the ``.xtc`` files).

.. note::

   **A file appears only if the run wrote at least one row to it.** Output files
   are created lazily, at the moment there is something to put in them, so the
   set of files a run leaves behind is exactly the set of outputs it actually
   produced. An analysis that is switched off produces no file; so does an
   analysis that is switched on but never fires, because its frequency never
   lands in the production window (a run of 500 steps with ``ANA_POL : 1000``,
   or the default ``EN_FREQ : 1000`` in a run shorter than that, writes no
   ``RG.dat`` and no ``ENERGY.dat``). The same holds for an analysis that fires
   but has nothing to say: ``RES_TO_RES_DIST.dat`` needs an
   ``ANA_RESIDUE_PAIRS`` pair, and ``CLUSTER_RADIAL_DENSITY_PROFILE.dat`` needs
   a cluster that reached the bead threshold. The exceptions are deliberate and
   few: the per-cluster property files (``CLUSTER_RG.dat`` and its companions)
   get a row at every cluster-analysis step, holding just the step number when
   no cluster exceeded ``ANA_CLUSTER_THRESHOLD``, so they keep one row per
   ``CLUSTERS.dat`` row; and a chain type of single beads still writes its
   sentinel ``SCALING_INFORMATION.dat`` rows and a 1 x 1 distance map (both
   described below).

   Do not read a missing file as an error, then, and do not read one as evidence
   that an analysis was disabled - check the keyfile (or ``parameters_used.prm``
   and ``log.txt``) for that.

   Per-step files are appended each analysis step, so a column-per-quantity,
   row-per-step layout is typical. A handful of analyses (internal scaling,
   distance maps) are **accumulated over the run and written once at the end**;
   these are noted below.

   A row labelled step ``N`` describes the state **after** Monte Carlo step
   ``N`` has completed. Step 0 exists only as the initial trajectory frame.
   Steps taken inside a TSMMC temperature excursion belong to an auxiliary
   Markov chain and produce no rows of their own; the master step that owns the
   excursion is reported once, when the excursion finishes.

.. note::

   **Equilibration.** All ``ANA_*`` analysis (and the periodic restart snapshot)
   is skipped until production begins. Step ``EQUILIBRATION`` is the *last*
   equilibration step, so the first analysis row is written at the first
   multiple of the analysis frequency strictly greater than ``EQUILIBRATION``.
   ``ENERGY.dat``, ``PERFORMANCE.dat``, ``QUENCH.dat`` and the trajectory are
   written throughout, equilibration included (trajectory frames only if
   ``SAVE_EQ : True``).

.. warning::

   Output files belong to a **single run**: every run start deletes every output
   file it could itself write - the per-step ``.dat`` files, the end-of-run
   files, ``QUENCH.dat``, the trajectory and every stale per-chain-type
   ``CHAIN_<T>_*.dat`` file left by a previous run in the same directory - so a
   re-run in a used directory never inherits an earlier run's data, not even for
   the files it does not itself write. If one of those files cannot be deleted
   (it belongs to another user, or a directory sits under its name) the run
   still starts, but a warning naming the file is printed and written to
   ``log.txt``: rows the new run writes under that name would be appended to
   the old contents, so remove the file by hand or use a clean directory. An
   output file that is a symbolic link is replaced by a regular file; the
   file the link pointed at is never written to. Replacing (like the deletion
   above) is decided by the permissions of the directory, not of the old file:
   a read-only ``START.pdb`` or ``.dat`` file from an earlier run is replaced
   like any other, so write-protecting a file does not protect it from a new
   run in the same directory - keep finished runs in their own directories. The deletion happens as the run starts,
   once the keyfile has been accepted and the system built, so a run refused
   at start-up leaves the previous run's ``.dat`` files and trajectory where
   they were. A refusal while the system is being built (an overcrowded box, or
   a restart snapshot using a bead type the parameter file lacks) comes after
   ``log.txt`` has been restarted and ``parameters_used.prm`` rewritten, so
   those two already belong to the refused run; a keyfile refusal (anything the
   keyfile or restart-file checks reject) does not even touch ``log.txt``, and
   is reported on the terminal only.
   ``restart.pimms`` is the one exception: it is overwritten when the first
   checkpoint is written rather than at start-up, so a run that dies before its
   first checkpoint leaves the previous run's restart file in place. When you
   resume from a restart file, run each segment in its **own directory** -
   resuming in the directory of the previous segment destroys that segment's
   outputs. A run started from a restart file numbers its steps from 1, unless it
   is an exact continuation (``RESTART_CONTINUE : True``), whose rows carry on
   from the checkpoint's step; see :ref:`restart-continue`.

.. _output-catalogue:

The output files
================

Core state & performance
------------------------

``ENERGY.dat``
    Instantaneous potential energy. Two tab-separated columns: ``step``,
    ``energy``. Frequency: ``EN_FREQ``, which must be greater than 0, so every
    run of at least ``EN_FREQ`` steps writes this file, equilibration
    included. The energy is an exact integer written with four decimal zeros
    (``-1234.0000``); it is written from the integer itself, so it stays exact
    however large it is (the same holds for the energy column of
    ``QUENCH.dat``).

``PERFORMANCE.dat``
    Throughput and timing. Columns: ``step``, ``E``/``P`` (equilibration vs
    production), loop-steps per second, overall MC moves per second (counting all
    sub-loop moves), elapsed time and estimated remaining time (``hh:mm:ss``). A
    header line is written when the file is created. Rows are written at every
    multiple of ``round(N_STEPS / 20)`` steps - roughly every 5 % of the run, so
    about twenty rows, but not exactly: ``N_STEPS : 72`` writes every 4 steps
    (18 rows) and ``N_STEPS : 115`` every 6 steps, so its last row is step 114
    and the final step gets no row. When ``N_STEPS`` is below 30 a twentieth of
    the run rounds to a single step and a row is written every step. One extra
    row is written at step 20 so an early throughput estimate is available;
    there is no keyword frequency for this file. Elapsed time and the rates are
    measured on a monotonic clock, so a change to the system clock during the
    run (daylight saving, a time-server correction) does not disturb them.

``QUENCH.dat``
    Only written for quench runs (``QUENCH_RUN : True``); any copy left by a
    previous run is removed at start-up. Three tab-separated columns: ``step``, ``temperature``,
    ``energy`` - the temperature ramp and the post-move energy response at that
    temperature. One row per ``QUENCH_FREQ`` steps at which the temperature
    actually changed (so rows stop once ``QUENCH_END`` is reached), written
    during equilibration as well as production. The temperature column carries
    six significant figures, so a fractional ``QUENCH_STEPSIZE`` round-trips
    exactly.
    See :doc:`advanced/quench`.

Move bookkeeping
----------------

These three are written together by ``ANA_ACCEPTANCE``. The per-move columns are
indexed by **move code** (1 = crankshaft, 2 = chain-translate, 3 = chain-rotate,
4 = chain-pivot, 5 = head-pivot, 6 = slither, 7 = cluster-translate,
8 = cluster-rotate, 9 = CTSMMC, 10 = multichain-TSMMC, 11 = pull,
12 = system-TSMMC, 13 = jump-and-relax, 14 = VMMC).

``MOVE_FREQS.dat``
    Cumulative number of moves **attempted** of each type since the start of
    the run (the counters are never reset, so each row is a running total, like
    ``TOTAL_MOVES.dat``, and the totals include the moves made during
    equilibration even though the first row is written in production).
    Columns: ``step`` then 14 counts, one per move code 1-14. The
    megamoves (crankshaft, slither, pull) contribute every one of their
    per-megamove substeps, so their counts grow much faster than the step
    number (with ``PARALLELIZE : True`` a parallel sweep counts the attempts it
    actually made, so a sweep in which nothing could move adds 0); the TSMMC
    columns (9, 10 and 12) count whole temperature excursions rather than the
    substeps inside them.

``ACCEPTANCE.dat``
    Cumulative number of moves **accepted** of each type (same layout, also a
    running total). Divide by ``MOVE_FREQS.dat`` to obtain per-move acceptance
    ratios.

    Read those ratios as "how often did this move change the configuration",
    with one caveat. The two rotations (codes 3 and 8) reject draws that map the
    body exactly onto itself, so their ratios mean exactly that. The chain pivot
    (code 4) and the chain translation (code 2) do not: a pivot of a straight arm
    and a translation by a zero offset are both counted as accepted although
    nothing moved. Acceptance ratios are therefore comparable between runs of the
    same move, but not across those different moves.

``TOTAL_MOVES.dat``
    Cumulative total number of accept/reject operations across all sub-loops
    (``step``, ``total``) - the "true" amount of MC work done. This includes the
    TSMMC excursion substeps and the relaxation substeps of jump-and-relax,
    which is why it can be orders of magnitude larger than the step number.

Single-chain (polymeric) analysis
----------------------------------

Every single-chain observable below is computed on the chain **made whole**:
under periodic boundaries the beads are bond-walked into one periodic image
before any distance or shape tensor is formed, so a chain that crosses a box
face (or that simply spans more than half the box along an axis) is treated as
the contiguous object it is. Under ``HARDWALL`` the raw positions already form
a single image and are used directly. For chains shorter than half the box this
is identical to minimum-image geometry; beyond that, minimum-image quantities
are wrong for a bonded object (a nearly straight chain of length ``0.8 L`` would
report an end-to-end distance of ``0.2 L``) and the whole-chain values are the
meaningful ones. A chain found to span more than half the box triggers a
finite-size warning (printed once per chain), since it is almost certainly
interacting with its own periodic image.

In every per-chain file the data columns run in ascending ``chainID`` order,
which is the order the chains appear in ``chain_to_chainid.txt``, in the
``START.pdb`` topology and in the trajectory. ``chainID`` numbering starts at 1
and follows keyfile ``CHAIN``-line order. A run started from a restart file
keeps the snapshot's chainIDs, and any ``EXTRA_CHAIN`` chains are numbered
after them, so there the columns of one chain type need not be contiguous
(see :ref:`restart-extra-chains`).

``RG.dat`` / ``ASPH.dat``
    Per-chain radius of gyration / asphericity, both derived from the
    gyration-tensor eigenvalues. Each row is a ``step`` followed by one value per
    chain. Trigger: ``ANA_POL``.

``END_TO_END_DIST.dat``
    Per-chain end-to-end distance (a ``step`` column followed by one value per
    chain). Written at the same frequency as ``RG.dat``/``ASPH.dat`` (``ANA_POL``).

``RES_TO_RES_DIST.dat``
    Distance between a chosen residue pair, for every chain. Columns: ``step``,
    the two residue indices (0-based), then one distance per chain. Trigger:
    ``ANA_INTER_RESIDUE`` with ``ANA_RESIDUE_PAIRS`` set; ``ANA_RESIDUE_PAIRS``
    may be repeated, in which case each recorded step contributes one row per
    pair. With no pair defined there is nothing to measure and the file is not
    written. Because the pair is measured on every chain, every chain must be
    long enough to contain it: a keyfile (or restart file plus ``EXTRA_CHAIN``)
    with a chain too short for a pair is refused at start-up, so alongside, for
    example, a chain type of single beads the only pair accepted is the
    degenerate ``0 0``, whose distance is always zero.

``INTSCAL.dat`` / ``INTSCAL_SQUARED.dat``
    Internal scaling versus sequence separation, one row per gap ``|i-j|`` from
    ``1`` to ``L-1`` (the last row is the end-to-end separation). Columns:
    ``gap``, ``mean``. ``INTSCAL.dat`` holds the mean distance ``<r_ij>``;
    ``INTSCAL_SQUARED.dat`` holds the **mean squared distance** ``<r_ij^2>``
    (averaged over every pair at that gap, every chain the file covers and every
    sample - not the square of the mean), so ``sqrt`` of it is the RMS
    internal-scaling profile. Accumulated over the run at the ``ANA_INTSCAL``
    sampling frequency and written once at the end of the run. Not written when
    ``ANA_INTSCAL`` is disabled, nor when no sample was ever taken because no
    multiple of ``ANA_INTSCAL`` fell in the production window (typically because
    it exceeded the production length; a warning is printed and logged instead of
    a file of zeros), nor for a chain type of single beads, which has no sequence
    separations (its ``SCALING_INFORMATION.dat`` sentinel rows and its 1 x 1
    ``DISTANCE_MAP.dat`` are still written).

``SCALING_INFORMATION.dat``
    Fitted polymer-scaling parameters: one tab-separated row **per chain** giving
    the apparent scaling exponent ``nu`` and prefactor ``R0`` from
    ``R = R0 · N^nu``, fitted to that chain's own ``<r_ij^2>`` profile (not to
    the type-averaged one written to ``INTSCAL_SQUARED.dat``). The fit discards
    the first 15 sequence separations and uses at most 41 log-spaced points, the
    largest separation always among them. A chain too short to fit (fewer than 26
    beads, i.e. fewer than 25 sequence separations) writes the sentinel row
    ``-1.0000 -1.0000``, so in a system of short chains every row is a sentinel.
    Written once at the end of the run, alongside ``INTSCAL.dat``, and on the same
    terms: not written when ``ANA_INTSCAL`` is disabled, nor when nothing was
    ever sampled.

    Before 1.0.8 the never-sampled case was an exception and wrote a file of
    ``-1`` rows. That could not be read as evidence of anything, because it is
    byte-identical to what a fully sampled run of sub-26-bead chains writes. Use
    the warning (on stdout and in ``log.txt``), or the absence of ``INTSCAL.dat``
    alongside the absence of ``SCALING_INFORMATION.dat``, to tell the two apart.

``DISTANCE_MAP.dat``
    Mean inter-residue distance map - a full symmetric ``seqlen × seqlen`` matrix
    (tab-separated rows; zero on the diagonal), accumulated at the ``ANA_DISTMAP``
    sampling frequency and written once at the end of the run. Not written when
    ``ANA_DISTMAP`` is disabled or was never sampled (again with a warning). A
    chain type of single beads gets a 1 x 1 map holding a single ``0.0000``.

    The cost of this analysis grows with the **square of the chain length**.
    Every chain keeps its own ``seqlen × seqlen`` float64 running mean (8 bytes
    per entry, allocated the first time the map is sampled, so a run with
    ``ANA_DISTMAP : 0`` never pays for it) and holds it for the rest of the run.
    On top of those the analysis needs about two more matrices of working
    memory - while a chain is sampled, and again when the maps of one chain type
    are averaged and written - so the peak is about
    ``(chains + 2) × 8 × seqlen²`` bytes per chain type. The file itself is
    about ``8 × seqlen²`` bytes. For chains of a few hundred beads this is
    negligible; a single 10,000-bead chain needs about 2.4 GB of memory and
    writes an 800 MB file. At the first sample PIMMS prints and logs a warning
    naming ``ANA_DISTMAP``, the chain type, its chain count and length and both
    estimates, for every chain type whose estimated peak memory exceeds 1 GiB
    or whose file would exceed 256 MB. Set ``ANA_DISTMAP : 0`` to switch the
    map off for such a system.

For multi-component systems the internal-scaling/distance-map files are written
per chain type as ``CHAIN_<TYPE>_INTSCAL.dat`` etc. **instead of** the unprefixed
files (the unprefixed names are used only for single-chain-type systems).
``<TYPE>`` is the 0-based integer chain-type index, in keyfile ``CHAIN``-line
order (for a run started from a restart file, the snapshot's own chain-type
indices, with any new ``EXTRA_CHAIN`` sequence taking the next free index). In
that case ``CHAIN_<TYPE>_INTSCAL.dat`` and
``CHAIN_<TYPE>_DISTANCE_MAP.dat`` are averaged over the chains of that type
alone, and ``CHAIN_<TYPE>_SCALING_INFORMATION.dat`` holds one row per chain of
that type.

Cluster analysis
----------------

Enabled by ``ANA_CLUSTER``. PIMMS identifies **short-range clusters** (chains in
Chebyshev-1 contact) and **long-range clusters** (chains connected via any pair
with *nonzero interaction energy*: a short-range contact between any beads, or a
Chebyshev-2 / Chebyshev-3 pair whose LR / SLR table entry for the two bead types
is nonzero). The table entry is the whole test in the long-range case, and it is
zero unless both beads are LR-capable, so a pair with only one LR bead never
connects two chains. Nor does a pair of LR-capable beads whose entry happens to be
zero (a parameter file with no SLR column, say), since such a pair carries no
energy. PIMMS reports both the size distributions and per-cluster shape
descriptors.

Clusters are **sorted largest first** in every row of every cluster
file. ``ANA_CLUSTER_THRESHOLD`` gates only the per-cluster shape/size
analysis: components with **more** than this many chains (strict comparison, so
the default of 1 skips single chains) get ``CLUSTER_RG``/``ASPH``/``AREA``/
``VOL``/``DEN`` and radial-profile entries, while the size-distribution files
(``CLUSTERS.dat`` / ``NUM_CLUSTERS.dat`` and the ``LR_*`` variants) always
include every component. Because the clusters are size-sorted, the clusters
above the threshold are a prefix of the list, so column ``k`` of
``CLUSTER_RG.dat`` and friends describes the same cluster as column ``k`` of
``CLUSTERS.dat``.

``CLUSTERS.dat`` / ``NUM_CLUSTERS.dat``
    Per-step cluster size distribution (``CLUSTERS.dat``: the ``step`` followed
    by the comma-separated cluster
    sizes, each the number of **chains** in that cluster) and the number of
    clusters (``NUM_CLUSTERS.dat``: tab-separated ``step``, ``count``).

``CLUSTER_RG.dat`` / ``CLUSTER_ASPH.dat`` / ``CLUSTER_AREA.dat`` / ``CLUSTER_VOL.dat`` / ``CLUSTER_DEN.dat``
    Per-cluster radius of gyration, asphericity, surface area, volume and density
    (one value per cluster per step). A step at which no component exceeded
    ``ANA_CLUSTER_THRESHOLD`` still writes a row, containing just the step number.
    The convex-hull quantities (area, volume,
    density) are ``-1`` for a cluster whose hull is degenerate - all beads
    collinear (3D: coplanar), as for a straight rod - because no hull can be
    built; treat ``-1`` as "undefined", not as a measurement. A cluster that is
    connected to its own periodic image is marked the same way: ``CLUSTER_RG``
    and ``CLUSTER_ASPH`` are written as ``nan``, the three convex-hull columns as
    ``-1``, and no radial profile row is written for it at all (see
    :ref:`percolating-clusters` below). In **2D** simulations the convex-hull
    quantities follow scipy's convention: ``CLUSTER_VOL.dat`` holds the polygon
    *area*, ``CLUSTER_AREA.dat`` holds the *perimeter*, and ``CLUSTER_DEN.dat``
    is beads per unit area. Note also that the 2D "asphericity" is
    :math:`\kappa = |\lambda_1-\lambda_2|/(\lambda_1+\lambda_2)` while the 3D
    value is the relative shape anisotropy :math:`\kappa^2`; square the 2D value
    to compare against 3D.

``CLUSTER_RADIAL_DENSITY_PROFILE.dat``
    Radial density profile (density vs distance from the cluster centre of mass),
    computed for clusters containing at least 27 beads. Each row is the ``step``,
    then a ``C<n>`` label, then the density values; ``<n>`` is the cluster number
    in the corresponding size-sorted ``CLUSTERS.dat`` row, so labels can have gaps
    when a smaller cluster is below the bead threshold, and a step at which no
    cluster reached 27 beads contributes no row at all (a run in which no step
    ever does writes no file). Shell ``k`` holds the
    fraction of the lattice sites at Chebyshev distance ``k`` from the (rounded)
    centre of mass that are occupied by beads of that cluster (a site held by a
    bead of another cluster counts as empty) - the first density column is shell 1 (the
    26 / 8 sites around the centre), and the profile runs out to
    ``min(DIMENSIONS) // 2 - 1`` shells, zero-padded; under ``HARDWALL`` only the sites that lie
    inside the box count, so a cluster wetting a wall is not diluted by sites
    that do not exist. Under periodic boundaries a shell is a set of lattice
    *sites*, so the site at ``COM - k`` is the site at ``COM + L - k`` and beads
    are binned by their minimum-image Chebyshev distance from the centre of
    mass. That only matters for a strongly one-sided cluster with beads more
    than about half a box from its own centre of mass; before 1.0.8 such beads
    were binned by their raw single-image distance and so fell outside every
    shell, which made the outermost shells read too dilute. The shell cap of
    ``min(DIMENSIONS) // 2 - 1`` is what guarantees the minimum-image shells do
    not overlap, so the shell site count is exact.

    The centre is a lattice **site**: the arithmetic mean of the cluster's
    single-image bead positions, rounded to the nearest site, with a mean that
    falls exactly half-way between two sites rounded up on that axis
    (``floor(mean + 0.5)``, the rule lemonade's radial profile uses too). The
    profile is therefore unchanged when the whole configuration is translated
    by whole lattice sites. Before 1.0.8 the half-way case was rounded to the
    nearest *even* coordinate, so the same cluster one site along could get a
    different profile, by up to a few tenths in a shell density; rows of
    ``CLUSTER_RADIAL_DENSITY_PROFILE.dat`` / ``LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat``
    for a cluster whose mean position falls exactly on a half-integer (a few
    percent of small clusters) can therefore differ from those written by
    earlier versions. Exact mirror invariance is not attainable for a
    site-centred profile on a lattice: when the mean lies half-way between two
    sites one of them has to be chosen, and the mirrored cluster ends up centred
    on the mirror image of the other, so a cluster and its mirror image can
    still give different profiles in that half-integer case.

``LR_CLUSTERS.dat``, ``NUM_LR_CLUSTERS.dat``, ``LR_CLUSTER_RG.dat`` (etc.)
    The same set of files for the **long-range** clusters: besides the three
    named above, ``LR_CLUSTER_ASPH.dat``, ``LR_CLUSTER_VOL.dat``,
    ``LR_CLUSTER_AREA.dat``, ``LR_CLUSTER_DEN.dat`` and
    ``LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat``, each with the layout of its
    short-range namesake. Before the shape
    descriptors are computed each cluster is gathered into a single periodic
    image, and the long-range gather walks exactly the relation that defines
    long-range membership - a contact, or a Chebyshev-2 / Chebyshev-3 pair with
    a nonzero LR / SLR table entry. Because of that, whenever the long-range and
    short-range cluster memberships coincide (they always do with a short-range
    parameter file, where both tables are zero) the ``LR_CLUSTER_*`` values are
    identical to the ``CLUSTER_*`` values. Before 1.0.8 the long-range gather
    used a plain "anything within three sites" rule instead, which linked the
    two ends of an elongated cluster through the periodic face even when they
    carried no interaction and so reported a torn, artificially compact object.

For multi-component systems, ``CHAIN_<TYPE>_CLUSTERS.dat`` (and the long-range
``CHAIN_<TYPE>_LR_CLUSTERS.dat``) records, per row, the ``step`` followed by the
fraction of each cluster's chains that are of that type (a fraction of the
chain count, not of the beads), in the same
cluster order as ``CLUSTERS.dat``. (The leading step column
was added in 1.0.8 - older files relied on line alignment with ``CLUSTERS.dat``.)

.. _percolating-clusters:

Clusters that percolate the box
-------------------------------

A cluster that is connected to its own periodic image - a spanning fibril, a
condensate that has grown to touch itself across the boundary - is an unbounded
object in the infinite periodic system the simulation represents. It has no
radius of gyration, no asphericity, no convex hull and no radial density
profile: the gather still returns coordinates, but they are an arbitrary finite
window cut out of an infinite object, and which window you get depends on which
bead the walk started from. PIMMS therefore refuses to report those quantities
as measurements:

* ``CLUSTER_RG.dat`` / ``CLUSTER_ASPH.dat`` (and the ``LR_`` variants) get ``nan``, which is
  what a naive ``mean()`` over the column should choke on rather than quietly
  average.
* ``CLUSTER_VOL.dat`` / ``CLUSTER_AREA.dat`` / ``CLUSTER_DEN.dat`` get ``-1``, the same
  "undefined" convention already used for a degenerate hull.
* ``CLUSTER_RADIAL_DENSITY_PROFILE.dat`` gets no row for that cluster. The ``C<n>``
  labels index the ``CLUSTERS.dat`` column explicitly, so the remaining rows still
  point at the right clusters.
* One line is written to ``log.txt`` (and printed) for **every** analysis step at
  which this happens, naming the step and the clusters involved.

The size distributions themselves (``CLUSTERS.dat``, ``NUM_CLUSTERS.dat`` and the
``LR_`` variants) are untouched: a cluster size is a well-defined integer whatever
the cluster's topology, and it is how you diagnose the percolation in the first
place. A long-range cluster counts as percolating only when a pair of its beads
genuinely *interacts* through the face; merely coming close to its own image is
not enough. Hardwall boxes do not wrap and can never percolate.

The test is exact. A cluster is gathered into one image by walking its contacts
(for a long-range cluster, the contacts and the Chebyshev-2 / Chebyshev-3 pairs
with a nonzero LR / SLR entry that define it), and it percolates if and only if
some pair of beads joined by one of those links is displaced, in the gathered
image, from the lattice offset that links them by a non-zero number of box
lengths: that link closes a loop around the box. It does not matter how many
box lengths the loop covers, so a helical cluster that only meets itself after
winding a short box axis twice is flagged like any other, while a cluster that
merely spans the box without closing is not: a straight rod one bead shorter
than the box, or a staircase that crosses the whole of one axis while its two
ends stay more than a link apart on another (in a box that is longer along
that other axis, say). A diagonal that runs corner to corner of a square or
cubic box does close, through the corner, and is flagged on every axis it
advances along. The test costs a few hundred bytes per bead, so it is
safe on a system-spanning cluster of any size, and each cluster is tested once
per analysis step.

Trajectory
----------

``START.pdb``
    Topology file - one ``ATOM`` record per bead, chains labelled by type (one PDB
    chain identifier per chain type, ``A-Z`` then ``a-z`` then ``0-9``, handed
    out in the order the types first appear in chainID order - so ``A`` for the
    first ``CHAIN`` line, ``B`` for the second, and so on; an ``EXTRA_CHAIN``
    that joins an existing type shares that type's identifier; every type past
    the 62nd shares the identifier ``9``, which is announced at start-up). Used
    as the topology when loading the trajectory.
    Every bead is written as atom name ``CA``; the residue name is the standard
    three-letter code when the bead type is a one-letter amino acid (``A`` ->
    ``ALA``) and otherwise the bead name left-padded with ``X`` (``B`` ->
    ``XXB``). ``CONECT`` records carry the backbone bonds, so viewers draw the
    chains. Residue numbers restart at 1 in every chain (and roll over into a new
    segment past 9999 within a chain longer than that). Atom serials number every
    ``ATOM`` and ``TER`` record and are written modulo 100000 (the 5-column PDB
    limit), so in a file of ``N`` > 99,999 serials (roughly 100,000+ beads) the
    serials 1 to ``N - 100000`` each appear twice. mdtraj/VMD rebuild atom
    indices from the record order, so the beads themselves load correctly, but a
    ``CONECT`` record can only name a serial. A bond is therefore written only
    when both of its serials occur once in the file (``N - 99999`` to 99999),
    and every other bond is omitted: no ``CONECT`` ever joins the wrong
    beads, but the beads at the start and the end of the file carry no bonds, so
    viewers draw those chains unbonded and mdtraj's ``find_molecules()`` splits
    them. Nothing that reads the file back depends on either column:
    ``lemonade`` builds its chains from the ``TER`` blocks and the order of the
    ``ATOM`` records, so a 100,000+ bead trajectory loads and analyses exactly as
    a small one does.
    Its ``CRYST1`` record gives the
    **periodic unit cell**, so an axis of ``L`` lattice sites is written as
    ``L * LATTICE_TO_ANGSTROMS`` angstroms (sites ``L-1`` and ``0`` are periodic
    neighbours one lattice unit apart) - the same box that is written into every XTC
    frame. For a 2D system the ``c`` axis is one lattice unit, since there is no
    periodicity in z. ``CRYST1`` holds three decimals; the XTC frames carry the
    cell at full (single) precision, so prefer the trajectory's box when the
    two differ in the last digit.

    The PDB format gives a coordinate eight columns, so it can only hold values
    from -999.999 to 9999.999 A. The longest box edge (the largest
    ``DIMENSIONS`` value times ``LATTICE_TO_ANGSTROMS``, including the box of a
    ``RESIZED_EQUILIBRATION``) must therefore stay below 10000 A. A bead beyond
    the limit stops the run with an error that names both keywords; the file is
    written under a temporary name and moved into place only when complete, so
    a failed write never leaves a truncated ``START.pdb`` (or replaces an
    earlier run's). ``LATTICE_TO_ANGSTROMS`` only rescales the output
    coordinates, so a smaller spacing is a safe way to fit a very large box.

    With ``AUTOCENTER : True`` the single chain is translated so that its
    centre of mass sits on the box centre, in ``START.pdb`` and in every frame.
    Under periodic boundaries the ends of a lopsided chain may then lie outside
    the box, which is a valid periodic image. Under ``HARDWALL : True`` there is
    no such image, so the shift is clamped: the chain is centred as far as the
    walls allow and every bead is written inside ``[0, L)``. ``AUTOCENTER`` is
    defined for a single chain only; with more than one chain it is ignored and
    a warning says so at start-up.

``traj.xtc``
    The trajectory itself (XTC format, via ``mdtraj``), one frame every
    ``XTC_FREQ`` steps, captured after that step's move. Equilibration frames are included only if
    ``SAVE_EQ : True`` - except **frame 0**, which is always the starting
    configuration (opening the writer records it regardless of ``SAVE_EQ``);
    with ``SAVE_AT_END : True`` frame 0 is written at start-up, the later
    frames are buffered in memory and added to the file once at the end. The
    two modes write the same file, byte for byte: the same coordinates, the
    same frame stamps and the same unit cell (``DIMENSIONS`` times
    ``LATTICE_TO_ANGSTROMS``, exactly orthorhombic). Coordinates are scaled by
    ``LATTICE_TO_ANGSTROMS``. (When ``RESIZED_EQUILIBRATION`` is used with
    ``SAVE_EQ : True`` the equilibration phase is written separately as
    ``eq_START.pdb`` / ``eq_traj.xtc``.)

    By default the raw on-lattice positions are written, so under periodic
    boundaries a chain that crosses a box face appears split across the two faces.
    Set ``TRAJECTORY_PBC_UNWRAP : True`` to make every chain **whole** before each
    frame (and the ``START.pdb`` topology) is written - each chain is shifted into a
    single periodic image, so no chain is torn across a boundary. This is purely a
    visualisation convenience (it does not affect the simulation), and unwrapped
    coordinates may fall outside the box; it has no effect under ``HARDWALL``, where
    chains never cross a boundary.

    Each frame carries sequential ``time``/``step`` metadata (0, 1, 2, ... in
    saved-frame order - these count saved frames, not Monte Carlo steps). With
    ``SAVE_EQ : True`` (the default) frame *k* is step ``k * XTC_FREQ``. With
    ``SAVE_EQ : False`` frame 0 is still the starting configuration (step 0) but
    equilibration frames are skipped, so frame *k* (*k* >= 1) is the *k*-th
    multiple of ``XTC_FREQ`` after ``EQUILIBRATION``, i.e. step
    ``(EQUILIBRATION // XTC_FREQ + k) * XTC_FREQ``. For the production
    ``traj.xtc`` of a ``RESIZED_EQUILIBRATION`` run frame 0 is the post-resize
    configuration at step ``EQUILIBRATION`` and frame *k* is the *k*-th multiple
    of ``XTC_FREQ`` after it. Likewise, for a ``RESTART_CONTINUE`` segment
    frame 0 is the checkpoint configuration (the restart file's step) and frame
    *k* is the *k*-th multiple of ``XTC_FREQ`` after that step, while the
    metadata still count 0, 1, 2, ... from the start of the segment.

    **What survives a run that stops early.** With the default
    ``SAVE_AT_END : False`` every frame is handed to the operating system as
    soon as it has been written, so the ``traj.xtc`` on disk always holds every
    completed frame. If the PIMMS process is killed outright (``SIGKILL``, an
    out-of-memory kill, a scheduler's wall-time limit) all frames written
    before the kill are there, and the file loads with ``md.load`` as usual
    unless the kill landed in the middle of writing a frame. In that case the
    file ends in a partly written (torn) frame after the last complete one.
    The same is true of a trajectory read while the run is still going: every
    completed frame is in the file, but a frame larger than about 4 kB (roughly
    800 beads or more) reaches the file in pieces, so a reader can catch the frame
    that is being written at that moment. ``md.load``, and a plain
    ``md.formats.XTCTrajectoryFile(...).read()``, both refuse a file with a
    torn final frame. ``pimms.lemonade.load`` reads such a file up to its last
    complete frame and warns with the number of frames it recovered; with
    ``mdtraj`` alone, read the frames one at a time and stop at the first that
    does not decode (:ref:`the recipe is given below <reading-a-torn-trajectory>`), or, for a run that
    is still going, simply try again a moment later. The frames are flushed, not
    synced: if the machine itself goes down (power loss, a kernel crash) the
    last few seconds of any file may be lost, as for every other output. With
    ``SAVE_AT_END : True`` the buffered frames exist only in memory. They are
    written if the run stops on an error or is interrupted (``Ctrl-C``,
    ``SIGTERM``, a failed ``ENERGY_CHECK``), but a ``SIGKILL`` or a crash of the
    machine loses them all, and ``traj.xtc`` then holds frame 0 alone. If the
    buffered frames cannot be written (the disk is full, say) ``traj.xtc`` is
    left holding whole frames only - never part of one - and the error is
    reported. What the
    other files of a stopped run hold, and how to continue from it, is described
    in :ref:`restart-stopping`.

``eq_START.pdb`` / ``eq_traj.xtc``
    The equilibration-phase topology and trajectory, written **only** by a
    ``RESIZED_EQUILIBRATION`` run with ``SAVE_EQ : True``. Their ``CRYST1`` record
    and box vectors describe the (smaller) equilibration box, not the production
    box; the production ``START.pdb``/``traj.xtc`` pair opens fresh once the box
    is resized. Any pair left by a previous run is removed at start-up when the
    current run will not write them. Conversely, because a resized run opens its
    production pair only at the resize, it removes any ``START.pdb`` /
    ``traj.xtc`` left by a previous run at start-up, so a run that dies during
    equilibration never sits beside another run's production trajectory.

``CONFIG_AT_ENERGY_FAIL.pdb`` / ``CONFIG_AT_ENERGY_FAIL.xtc``
    A single-frame snapshot of the configuration (as the chain objects hold it)
    at which an ``ENERGY_CHECK`` failed: either the incrementally tracked energy
    disagreed with a from-scratch recompute, or the occupancy/type grids
    disagreed with the chain objects (up to ten of the offending sites of the
    latter are also reported to the terminal and to ``log.txt``). Either failure is recorded in
    ``log.txt`` with the name of this snapshot. The run aborts immediately
    afterwards, but ``traj.xtc`` is closed (or, under ``SAVE_AT_END``, flushed
    from the buffer) first, so the trajectory up to the failure is kept. These
    files should never appear; if they do, keep them and the ``log.txt``. Any
    pair left by a previous run is removed at start-up.

Echoed inputs & checkpoint
--------------------------

``keyfile_used.kf``
    The configuration the run *actually* used, written at start-up once the
    simulation has been built and its trajectory opened, so a run that failed
    to start never leaves one behind (``log.txt`` then carries a
    ``START FAILED - NO SIMULATION WAS RUN`` line instead). The keyfile
    you wrote is not always what ran: ``RESTART_OVERRIDE_HARDWALL`` and
    ``RESTART_OVERRIDE_DIMENSIONS`` replace those keywords with the restart file's
    values, a ``RESIZED_EQUILIBRATION`` forces hard walls for its first phase
    whatever ``HARDWALL`` says, ``EXTRA_CHAIN`` lines merge into existing chain
    types, and a missing ``SEED`` is generated at start-up. This file holds the
    resolved keyword set, one ``KEY : value`` line each, under a comment header
    that records the source keyfile, the PIMMS version, the date, the seed
    and where it came from (yours; generated, and reproducing this run;
    generated but not reproducing it, because the run began from a restart file;
    or generated and unused, because a ``RESTART_CONTINUE`` run restored the
    generators), the resolved box and boundary condition, the restart file the
    run started from and the step it was written at (noting when the overrides
    supplied the box or boundary, and when the run was a continuation), the box
    and boundary condition of a resized equilibration phase, and any
    ``FREEZE_FILE``. It re-parses as a keyfile: the override and continuation
    flags are written as ``False`` because their effect is already in the
    values, and ``ANALYSIS_MODULE`` is written as the module's path. For a
    restarted run it describes the equivalent system started afresh (the chains
    as a plain ``CHAIN`` list, one line per chain type with any ``EXTRA_CHAIN``
    already merged in; a restart snapshot has no keyfile form), and the header
    says which restart file it actually came from. When the snapshot's chainIDs
    are not grouped by chain type (after an ``EXTRA_CHAIN`` joined an existing
    type, say) the header also warns that a re-run from this file would number
    the chains differently. Where it differs from the keyfile you wrote, this is
    the one to trust, and it is the file to hand to ``lemonade.load``.

``parameters_used.prm``
    A copy of the parameter file actually used (with a timestamp header), so a run
    is self-documenting.

``absolute_energies_of_angles.txt``
    Human-readable summary of the angle penalties per residue, listing both the
    requested value and the integer the simulation actually applies. Written only
    when angle penalties are in use; with ``ANGLES_OFF : True`` no file is
    written and any copy left by a previous run is removed at start-up.

``chain_to_chainid.txt``
    Mapping of each chainID to its length and sequence (one tab-separated line
    per chain: ``chainID``, ``length``, ``sequence``), in ascending chainID
    order - the same order as the per-chain analysis columns. Written only if
    ``WRITE_CHAIN_TO_CHAINID : True`` (otherwise any file left by a previous run
    is removed); handy for choosing chains to put in a
    :ref:`freeze file <advanced-freeze>`.

``log.txt``
    A plain-text run log written on every simulation: the run header and start
    time, the resolved run settings (steps, equilibration, ``SAVE_EQ``,
    ``XTC_FREQ``, expected frame count, start and final temperature), the random
    seeds, one progress/throughput line for every row written to
    ``PERFORMANCE.dat``, and a "Simulation complete" line once the last step
    has been taken (just before the end-of-run files and the final
    ``restart.pimms`` are written). Handy for reconstructing exactly how a run was
    set up and whether it finished. It is not a transcript of the terminal,
    though: only warnings and errors that are explicitly logged reach it (the
    never-sampled analysis warnings, the percolating-cluster lines, the
    grid-consistency errors, an ``ENERGY_CHECK`` failure and the move-set
    applicability warnings and refusals, for example), while others are printed
    to the terminal alone, including the per-chain finite-size warning and the
    parameter-file interaction checks. The log is opened only once the
    keyfile has been accepted, so a keyfile or restart-file refusal never
    reaches it (and leaves any ``log.txt`` from a previous run untouched). Keep
    the captured stdout as well as the log. ``REDUCED_PRINTING`` trims what goes
    to the terminal only; it does not change what is written to ``log.txt`` or to
    any other output file.

``restart.pimms``
    Configuration checkpoint, rewritten every ``RESTART_FREQ`` production steps
    and always once more when the run finishes, so the file left behind is the
    final state. See :doc:`restart_files`. Each checkpoint is written to
    ``restart.pimms.tmp.<pid>`` (``<pid>`` being the process number) and moved
    into place only when complete, so ``restart.pimms`` itself is never
    half-written; a temporary left by a killed run is removed by a later
    start-up in that directory, at once if its process is known to be gone from
    that machine and otherwise when it is more than an hour old.

``pimms_running.pid``
    A marker holding one line (process number, host name and start time) for
    each run that is using the directory. A run adds its line when it starts
    and takes it out when it ends, and the file is removed with its last line,
    so it exists only while some run is going. Two runs in one directory delete
    and overwrite each other's output, so a run that finds the line of another
    run that may still be going says so at start-up, on screen and in
    ``log.txt``, and leaves that line where it is (a third run is then warned
    too, even after the second has finished). It is a notice rather than a
    lock: nothing is refused on its account. The line of a run that was killed
    on the same machine is dropped at the next start; the line of a run on
    another machine cannot be checked and stays, with a warning at every
    start, until the file is deleted by hand.

A custom analysis module (``ANALYSIS_MODULE`` / ``ANA_CUSTOM``) can of course
write anything it likes; PIMMS neither creates nor cleans up those files. See
:doc:`advanced/custom_analysis`.

.. _output-analyze:

Analysing the output
====================

Plain-text ``.dat`` files
-------------------------

Most ``.dat`` files are tab-separated; the exceptions are the cluster
size-distribution, composition, per-cluster property and radial-profile files
(``CLUSTERS.dat``, ``CHAIN_<T>_CLUSTERS.dat``, ``CLUSTER_RG.dat`` and friends,
``CLUSTER_RADIAL_DENSITY_PROFILE.dat``, and their ``LR_`` counterparts), which
are comma-separated. ``NUM_CLUSTERS.dat`` and ``NUM_LR_CLUSTERS.dat`` are
tab-separated despite the company they keep.

Most per-step rows end with a **trailing delimiter**, which matters when you load
them: for the tab-separated files simply pass no explicit ``delimiter`` to
``numpy.loadtxt`` (whitespace splitting ignores it), and for the comma-separated
files drop the empty final field. Any tool reads them; with NumPy:

.. code-block:: python

   import numpy as np

   # ndmin=2 keeps every table two-dimensional (rows x columns) even when the
   # file holds a single row, which np.loadtxt would otherwise return as a
   # one-dimensional array (and, with unpack=True, as bare numbers)

   # energy trace: columns [step, energy] (no trailing tab on this one)
   step, energy = np.loadtxt("ENERGY.dat", delimiter="\t", ndmin=2).T
   print("mean production energy:", energy[len(energy)//2:].mean())

   # per-move acceptance ratio
   # note: no explicit delimiter - rows end with a trailing tab, which
   # delimiter="\t" would parse as an empty final column and reject
   attempted = np.loadtxt("MOVE_FREQS.dat", ndmin=2)
   accepted  = np.loadtxt("ACCEPTANCE.dat", ndmin=2)
   ratio = accepted[:, 1:] / np.clip(attempted[:, 1:], 1, None)   # ratio[:, code - 1] = move code

   # per-chain observables: column 0 is the step, then one column per chain
   # in ascending chainID order
   rg = np.loadtxt("RG.dat", ndmin=2)

   # PERFORMANCE.dat has a header line and a text column, so pick columns
   perf = np.loadtxt("PERFORMANCE.dat", skiprows=1, usecols=(0, 2, 3), ndmin=2)

The cluster size distribution is ragged (a different number of clusters each
step), so read it line by line:

.. code-block:: python

   sizes = {}
   with open("CLUSTERS.dat") as fh:
       for line in fh:
           fields = [f for f in line.strip().split(",") if f.strip()]
           if not fields:
               continue
           sizes[int(fields[0])] = [int(v) for v in fields[1:]]   # largest first

The square ``DISTANCE_MAP.dat`` matrix loads directly with
``np.loadtxt("DISTANCE_MAP.dat")`` and can be shown with ``imshow``.

Trajectories (``.pdb`` + ``.xtc``)
----------------------------------

Load the ``.xtc`` frames with the ``START.pdb`` topology using `mdtraj
<https://www.mdtraj.org/>`_:

.. code-block:: python

   import mdtraj as md

   traj = md.load("traj.xtc", top="START.pdb")
   print(traj)                       # frames and beads (mdtraj reports them as n_atoms)
   # mdtraj coordinates are in nanometres; one lattice unit = LATTICE_TO_ANGSTROMS A,
   # so: lattice_units = traj.xyz * 10.0 / LATTICE_TO_ANGSTROMS
   rg = md.compute_rg(traj)          # radius of gyration per frame, etc.

.. _reading-a-torn-trajectory:

The trajectory of a run that was killed, or that is still going, can end in a
partly written frame, which ``md.load`` refuses (see ``traj.xtc`` above).
``pimms.lemonade.load`` copes with that by itself. With ``mdtraj`` alone, read
the frames one at a time and keep the ones that decode:

.. code-block:: python

   import mdtraj as md
   import numpy as np

   frames, boxes = [], []
   with md.formats.XTCTrajectoryFile("traj.xtc") as fh:
       while True:
           try:
               xyz, _time, _step, box = fh.read(n_frames=1)
           except (RuntimeError, OSError):
               break                     # a partly written final frame
           if len(xyz) == 0:
               break                     # the clean end of the file
           frames.append(xyz[0])
           boxes.append(box[0])

   traj = md.Trajectory(np.array(frames), md.load("START.pdb").topology)
   traj.unitcell_vectors = np.array(boxes)

The same ``START.pdb`` + ``traj.xtc`` pair loads directly in **VMD** (and other
molecular viewers) for visualisation. For 2D simulations the out-of-plane
coordinate is held at zero. If you used ``AUTOCENTER : True`` (single-chain runs)
the chain is already centred in ``START.pdb`` and in every frame, so no alignment
is needed before analysis (in a ``HARDWALL`` box it is centred as far as the
walls allow, and never written outside the box).
