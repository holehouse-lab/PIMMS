.. _lemonade-loading:

=====================
Loading a trajectory
=====================

Everything starts with :func:`pimms.lemonade.load`, which reads a finished PIMMS
run and returns a :class:`~pimms.lemonade.LatticeTrajectory`.

.. code-block:: python

   import pimms.lemonade as lemonade
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

What to pass
============

``load`` accepts the files you actually have to hand, in any sensible combination:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Inputs
     - What you get
   * - ``xtc`` + ``pdb``
     - The full trajectory. The PDB provides the topology (it matches the XTC bead
       order exactly); the XTC provides the coordinates over time.
   * - ``xtc`` + ``pdb`` + ``keyfile``
     - As above, **plus** authoritative box dimensions, lattice spacing, hardwall
       flag, temperature and chain *types* taken from the keyfile.
   * - ``pdb`` only
     - A single frame (e.g. ``START.pdb``) - handy for inspecting a starting
       configuration. The box is read from the PDB's ``CRYST1`` record, at any
       lattice spacing.

A PDB is always required alongside an XTC (mdtraj needs a topology to read the
trajectory). Passing the keyfile is optional but recommended - without it, lemonade
*infers* what it can (see below).

Where the numbers come from
===========================

PIMMS writes coordinates in nanometres as ``lattice_index x spacing / 10``.
lemonade inverts that with vectorised arithmetic - ``round(nm / (spacing/10))`` -
recovering the exact integer lattice (the round-off is float32 noise). The XTC is
read a block of frames at a time and each block is converted to ``int32`` before
the next is read, so a load holds the integer lattice of the frames it keeps (12
bytes per bead per frame) plus working memory that does not grow with the length
of the trajectory: about 20 MB (a 4 MB block of decoded coordinates and the arrays
it is converted through), or about five frames' worth of coordinates for a system
so large - above roughly 350,000 beads - that a single frame is bigger than that
block. mdtraj's topology of the PDB, roughly a kilobyte per bead, is held as well
while the file is read. The remaining metadata is resolved in this order:

* **spacing** - from the keyfile ``LATTICE_TO_ANGSTROMS``; otherwise PIMMS's default
  of ``3.65`` angstroms. Override with ``spacing=`` (a finite positive number;
  ``True`` / ``False`` are refused rather than read as ``1`` / ``0``).
* **dimensions** - from the keyfile ``DIMENSIONS``; otherwise inferred from the
  trajectory's box record: the XTC's, or the PDB's ``CRYST1`` record when the XTC
  carries no box or there is no XTC (before 1.0.8 a ``SAVE_AT_END`` trajectory at
  a lattice spacing below about an angstrom was written without one). PIMMS writes a 2D system with a ``z`` period of exactly
  one lattice unit, and that is what marks a trajectory as 2D - a 3D
  configuration that happens to lie in the ``z = 0`` plane stays 3D. Override with
  ``dimensions=(x, y, z)``. A restart override or a resized-equilibration
  trajectory changes where this comes from - see
  :ref:`lemonade-effective-keyfile`.
* **hardwall** - from the keyfile ``HARDWALL``; otherwise ``False``. Override with
  ``hardwall=``. Again, a restart override or a resized-equilibration trajectory
  changes where this comes from - see :ref:`lemonade-effective-keyfile`.
* **temperature** - from the keyfile ``TEMPERATURE`` (needed only for surface
  tension). For a ``QUENCH_RUN`` keyfile PIMMS ignores ``TEMPERATURE``, ramps from
  ``QUENCH_START`` to ``QUENCH_END`` and stays at ``QUENCH_END`` for the rest of
  the run, so ``QUENCH_END`` is what is used (with a warning; frames written
  during the ramp were sampled at intermediate temperatures). Override with
  ``temperature=``, which must be a finite positive number. With neither a
  keyfile nor ``temperature=``, ``traj.temperature`` is ``None``.
* **topology** (chain lengths, sequences, bead types) - always from the PDB, since
  it is written in lockstep with the trajectory. Chains come from the PDB's
  ``TER`` blocks and beads from the order of the ``ATOM`` records; the atom
  serial and residue number columns are never read, so the duplicated serials
  of a very large system (PIMMS writes them modulo 100000 once beads plus
  ``TER`` records pass 99,999) and the per-chain residue numbering make no
  difference. Without a keyfile, PIMMS PDB chain
  identifiers are preserved as chain-type labels (rather than guessing type from
  sequence); a PDB whose chain column is blank carries no type information at
  all, and there chains sharing a sequence are grouped into one type instead.
  PIMMS has 62 identifiers (``A-Z``, ``a-z``, ``0-9``); chain types past the
  62nd share the last one, so for such systems pass the keyfile. If a keyfile is
  given and its expanded ``CHAIN``/``EXTRA_CHAIN`` specification matches the PDB
  - the same sequence on every chain in order and, when the PDB carries chain
  identifiers, the same partition of the chains into types (each identifier
  holding exactly one keyfile type and each keyfile type exactly one
  identifier, except that when the PDB uses all 62 identifiers the keyfile
  types may split the shared last one - and only that one - which is how the
  keyfile recovers merged types) - those authoritative type labels are used,
  numbered in keyfile ``CHAIN`` order. A PDB that has run out of identifiers
  is also matched to keyfile lines that are not in chain order (the
  ``keyfile_used.kf`` of a restart run with ``EXTRA_CHAIN`` chains), by
  sequence, as long as no two lines share a sequence. PDB identifiers are numbered in the order they first appear
  in the file, which is the order PIMMS assigned them. A keyfile that
  describes the same types, i.e. the same ``(count, sequence)`` per type, but
  in a different order keeps the PDB labels without a warning, because they are
  then already right; any other mismatch emits a warning and retains the PDB
  labels. A *restart* keyfile is the exception - see below.

Whatever the trajectory was written with - wrapped positions, or whole molecules via
``TRAJECTORY_PBC_UNWRAP`` - lemonade canonicalises positions back into the box on
load and re-derives whole chains itself, so results do not depend on how the run was
saved.

.. _lemonade-effective-keyfile:

When the keyfile is not the run
===============================

The simplest way to avoid every case in this section is to pass the run's own
``keyfile_used.kf`` rather than the keyfile you wrote. PIMMS writes it at start-up
with every keyword already resolved (the box and boundary condition after any
restart override, the generated seed, the merged chain list) under a header that
records where it all came from, so lemonade can read it literally and be right;
see :doc:`/output_files`. The rest of this section is what lemonade does when it
is given the original keyfile instead.

lemonade reads the keyfile *literally*, which is almost always exactly what the run
did - but PIMMS has a handful of keywords whose whole purpose is to make the run
differ from the literal keyfile, and lemonade resolves those for you rather than
believing the file in front of it:

* **restart overrides.** ``RESTART_OVERRIDE_HARDWALL`` tells PIMMS to ignore the
  keyfile ``HARDWALL`` and adopt the restart file's boundary condition;
  ``RESTART_OVERRIDE_DIMENSIONS`` does the same for the box. When either is set,
  lemonade opens the ``RESTART_FILE`` and takes the value from there. A relative
  ``RESTART_FILE`` is resolved against the *keyfile's* directory first (that is
  where a run directory keeps it), then against the working directory. If the
  restart file cannot be found or read, ``load`` **raises**: the keyfile value is
  known to be the wrong one in that situation, so there is nothing sensible to fall
  back on - pass ``hardwall=`` / ``dimensions=`` yourself.
* **the resized-equilibration trajectory.** The compact phase of a
  ``RESIZED_EQUILIBRATION`` run is always simulated under hard walls, whatever
  ``HARDWALL`` says, and always in the smaller box. Loading ``eq_traj.xtc`` /
  ``eq_START.pdb`` therefore gives you ``hardwall=True`` and the
  ``RESIZED_EQUILIBRATION`` box, and warns saying so. The *filename* is what
  triggers this - PIMMS writes those two names for nothing else - so it works with
  or without a keyfile, and the production ``traj.xtc`` of the same run is
  untouched and keeps its real ``HARDWALL``. Passing ``hardwall=`` or
  ``dimensions=`` still wins.
* **chain types under a restart.** PIMMS discards the keyfile ``CHAIN`` lines when a
  ``RESTART_FILE`` is given (the composition comes from the snapshot instead), and
  an ``EXTRA_CHAIN`` whose sequence already exists in the restart joins that
  existing chain type rather than defining a new one. lemonade therefore does not
  apply the keyfile ``CHAIN``/``EXTRA_CHAIN`` block for a restart keyfile at all,
  and keeps the PDB chain identifiers - which carry PIMMS's real ``chainType``
  order, and are what you need to match its ``CHAIN_<type>_*`` output files.
  The restart run's ``keyfile_used.kf`` has no ``RESTART_FILE`` and lists the
  merged composition as one ``CHAIN`` line per type, while the trajectory keeps
  the snapshot's chains first and appends the ``EXTRA_CHAIN`` chains at the end,
  so its lines need not expand onto the chains in order. The partition check
  above is what catches that: when two types share a sequence the lines expand
  onto the chains with every sequence matching but the wrong types, and they are
  rejected in favour of the PDB labels (silently, since the composition agrees).

  The restart file itself is the one thing a restart keyfile can be checked
  against, so when it can be found (beside the keyfile, then in the working
  directory) lemonade reads its chains, appends the ``EXTRA_CHAIN`` chains, and
  compares the sequences with the PDB's: a mismatch means the keyfile is from
  another run, and warns. For a run with more than 62 chain types the same read
  recovers the types merged under the last PDB identifier, exactly as PIMMS
  assigned them. If the restart file is not there, neither happens and nothing
  is said, except that a 62-identifier PDB gets the collision warning below,
  which then asks for ``keyfile_used.kf``.

Selecting frames
================

You can subsample at load time (cheaper than loading everything and slicing after:
only the frames you keep are decoded from the XTC and held in memory):

.. code-block:: python

   # every 5th frame between 100 and 400
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", start=100, stop=400, step=5)

   # thin an arbitrarily long run down to 200 evenly spaced frames
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", n_frames=200)

``start``/``stop``/``step`` are applied first and ``n_frames`` then thins what
survives, so the two compose. ``n_frames`` is a no-op if the trajectory is
already that short or shorter. They behave exactly as a Python slice of the
frames does (negative values count from the end; a negative ``step`` returns
the frames in reverse).

.. _lemonade-running-trajectory:

Loading a run that is still going
---------------------------------

The ``traj.xtc`` of a run that is still writing, or of one that was killed,
usually ends part-way through a frame. ``load`` reads such a file up to its last
complete frame and warns with the number of frames it recovered; the frame
selection then counts within those frames, so ``start=-10`` is the last ten
*complete* frames. A file that ends cleanly is read directly and gives no such
warning. Finding the last complete frame of a torn file takes one extra pass
over it. This holds for systems of any size, including those of nine beads or
fewer, whose frames the XTC format stores uncompressed. (A frame damaged in the
middle of a file is not a torn tail: if the end of the file is intact, that load
still fails with the reader's own error.)

You can also slice *after* loading - ``traj[100:400:5]`` returns a new
``LatticeTrajectory`` over those frames. It does not re-read the files, and it
shares the parent's topology, but it holds its own copy of the selected frames'
positions and computes its own cached analyses (whole chains, Rg, clusters), so
it costs the memory of those frames; thinning at load time never holds the
frames you drop. An integer index gives a ``Frame`` instead; a slice or an index
array gives a trajectory. Either way, the selected frames keep their original
``times`` (``traj[10:].times`` starts at ``10``).

Checking the load
=================

``load`` never fails quietly on a suspicious input: it either raises, or warns.
The warnings fire on **every** load, regardless of ``verbose``:

* the **lattice round-off residual** - if the coordinates do not sit on the
  integer lattice at the given spacing (residual above 0.05, indicating a wrong
  or omitted ``LATTICE_TO_ANGSTROMS`` / ``spacing=``), a warning is raised, since
  the recovered lattice would be corrupted;
* the **bond check** - two bonded PIMMS beads are always on neighbouring lattice
  sites, so every bond of the first loaded frame must be a single lattice step
  (Chebyshev length 1 under the minimum-image convention). This catches what
  the residual cannot: a spacing that is a whole multiple of the assumed one -
  a ``LATTICE_TO_ANGSTROMS : 7.3`` run loaded without its keyfile, at the
  default 3.65 - puts every bead exactly on every second site, with no
  round-off at all, and silently doubles the positions and the box. The
  warning gives the spacing the shortest bond points to, to pass as
  ``spacing=``. It also fires for a PDB whose chains are not this trajectory's
  chains. A system made only of single-bead chains has no bonds to test, so
  nothing is checked and nothing is said there: pass the keyfile for those;
* the **incomplete last frame** - the trajectory of a running or killed run is
  loaded up to its last complete frame, and the warning says how many frames
  that is (see :ref:`lemonade-running-trajectory`);
* the **hardwall box check** - nothing crosses a hard wall, so a bead outside
  the box under ``hardwall`` was not put there by the simulation, and there are
  two ways it can have got there. For a *single-chain* system it is what PIMMS
  before 1.0.8 wrote under ``AUTOCENTER`` with ``HARDWALL``, when the centred
  chain could stick out through a wall. lemonade used to wrap those beads
  through the wall, which tore the chain in two; it now translates each such
  frame rigidly back into the box - the same shift the engine itself applies
  since 1.0.8 - which keeps the chain whole, and warns. Every frame of an
  ``AUTOCENTER`` run was re-centred when it was written, so positions relative
  to the walls are not meaningful in *any* frame of such a trajectory, whether
  or not it needed translating. With *several chains* (``AUTOCENTER`` never
  acted on those), or a frame too wide to fit in the box at all, the cause is
  a wrong box or boundary condition: the coordinates are wrapped, as they
  always were, and the warning says to check ``hardwall=`` and
  ``dimensions=``;
* the **box cross-check** - the given/keyfile ``DIMENSIONS`` are compared against
  the trajectory's own CRYST1/XTC box record, and a disagreement warns loudly.
  The *dimensionality* is compared first and separately, because that is the
  destructive case: 2D dimensions against a 3D trajectory throws away the ``z``
  coordinate of every bead, which merges clusters and shrinks every ``Rg``, and
  the warning says so explicitly. Both the keyfile and an explicit
  ``dimensions=`` argument are checked, since the resolved box is what the
  analysis uses;
* the **resized-equilibration adjustment** - loading ``eq_traj.xtc`` /
  ``eq_START.pdb`` reports that hard walls and the compact box have been applied
  in place of the keyfile's production values (see
  :ref:`lemonade-effective-keyfile`);
* the **keyfile/PDB chain mismatch** - a keyfile whose ``CHAIN``/``EXTRA_CHAIN``
  block describes different chain types from the PDB's (not merely the same types
  in a different order) cannot be from the same run, so its chain types are
  dropped and the PDB labels kept. This does not apply to a keyfile with a
  ``RESTART_FILE``, whose ``CHAIN`` block is never used in the first place;
  that keyfile is checked against its restart file instead, when the restart
  file can be found, and warns if the two describe different chains;
* the **62-identifier collision** - raised only when no keyfile resolved the
  types and the topology ends up with 62 or more of them. Nothing in the PDB
  says whether a merge actually happened, only that any type past the 62nd
  would have been folded into the last label, so lemonade says so and asks for
  the keyfile - or, when a keyfile was given and could not resolve them, for
  the run's own ``keyfile_used.kf``;
* the **quench temperature** - a ``QUENCH_RUN`` keyfile reports that
  ``QUENCH_END`` is being used in place of the (ignored) ``TEMPERATURE``; if it
  sets no ``QUENCH_END`` at all, it reports that no temperature could be
  recorded.

Errors, by contrast, are raised (all as ``ValueError``) for the inputs lemonade
cannot interpret at all: an ``xtc`` without a ``pdb``, neither file given, a
non-positive ``n_frames`` or ``spacing``, a box that is neither given, in the
keyfile, nor recorded in the trajectory, ``dimensions`` that are not 2 or 3
positive integers, a non-boolean ``hardwall``, a ``temperature`` that is not a
finite positive number, coordinates too large for int32, a PDB whose bead count
disagrees with the trajectory (mdtraj refuses that pair as it reads it), a
``start``/``stop``/``step`` selection that keeps no frames (a window past the
end of a short run), a box record smaller than one lattice site at the spacing
in use (the spacing is then certainly wrong), and a keyfile that sets
``RESTART_OVERRIDE_HARDWALL`` / ``RESTART_OVERRIDE_DIMENSIONS`` whose
``RESTART_FILE`` cannot be found or read. An XTC with no complete frame in it,
or one damaged before its end, fails with the XTC reader's own error.

A path with non-ASCII characters in it (an accented directory name) is fine: the
XTC reader itself only accepts ASCII paths, so lemonade opens such a file through
a temporary ASCII-named symbolic link. The link is made in a fresh directory
under the system temporary directory (``tempfile.gettempdir()``, so ``TMPDIR``
moves it), and both are removed again as soon as the file has been read. Where no
such link can be made, ``load`` raises a ``ValueError`` that says so; loading by a
relative path from inside the trajectory's own directory avoids the link
altogether.

Pass ``verbose=True`` additionally for a one-line summary. It gains a trailing
``(WARNING lattice round-off ...)`` when the residual is not float32 noise
(above ``1e-3``). The output shown here and below is from a periodic run of 250
eight-bead chains (``CHAIN : 250 SSSSSSSS``) in a 30 x 30 x 30 box at
``TEMPERATURE : 120``, written every 10 of 1000 steps:

.. code-block:: python

   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf",
                        verbose=True)
   # [lemonade] 101 frames, 250 chains, 2000 beads; box (30, 30, 30), spacing 3.65 A

The loaded object reports the basics directly:

.. code-block:: python

   traj.n_frames, traj.n_chains, traj.n_beads   # (101, 250, 2000)
   traj.dimensions          # (30, 30, 30)   - 2 entries for a 2D system
   traj.n_dim               # 3              - 2 for a 2D system
   traj.spacing             # 3.65
   traj.hardwall            # False
   traj.temperature         # 120.0  (None if unknown)
   traj.sequences           # ['SSSSSSSS', 'SSSSSSSS', ...] one per chain
   traj.chain_types         # (n_chains,) int32 type label per chain
   traj.times               # (n_frames,) float64 frame times from the XTC: 0.0, 1.0, ...
   traj.topology            # the chain/bead topology (see the hierarchy page)

``traj.n_atoms`` (and ``frame.n_atoms``) still work as deprecated aliases of
``n_beads`` and raise a ``DeprecationWarning``; PIMMS particles are beads.

**Which frames are in the file.** Frame 0 of ``traj.xtc`` is the starting
configuration (except in a ``RESIZED_EQUILIBRATION`` run - see below). With the default ``SAVE_EQ : True`` the equilibration frames are
included too, so a run with ``EQUILIBRATION : 1000`` and ``XTC_FREQ : 100`` has
ten equilibration frames before the first production one; set ``SAVE_EQ : False``
in the keyfile (or slice, ``traj[11:]``) to analyse production only. ``traj.times``
and ``frame.time`` are the *frame indices* ``0, 1, 2, ...`` as written by PIMMS,
not Monte Carlo steps: frame ``k`` is step ``k * XTC_FREQ`` (with ``SAVE_EQ :
False``, frame 0 is still step 0 and frame ``k >= 1`` is the ``k``-th production
frame).

A ``RESIZED_EQUILIBRATION`` run is different. Its equilibration frames go to
``eq_traj.xtc`` (only with ``SAVE_EQ : True``), and ``traj.xtc`` is opened at the
resize, so it holds **no** equilibration frames at all: frame 0 is the post-resize
configuration at step ``EQUILIBRATION``, and frame ``k >= 1`` is the ``k``-th
multiple of ``XTC_FREQ`` after it (step ``EQUILIBRATION + k * XTC_FREQ`` when
``EQUILIBRATION`` is a multiple of ``XTC_FREQ``). Every frame of it is production,
so the ``traj[11:]`` slice above would throw away the first ten production frames.
The most you would drop is frame 0, the compact configuration the resize started
from, which has not yet relaxed in the production box (``traj[1:]``). See
:doc:`/output_files`.

Coordinate, time and cached-analysis arrays are read-only. lemonade memoises
whole-chain coordinates and batched observables; immutability prevents an edit to
``traj.positions`` from leaving an already-cached Rg or COM describing a different
trajectory. Make an explicit ``.copy()`` when you need a mutable working array.
