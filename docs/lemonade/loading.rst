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
     - The full trajectory. The PDB provides the topology (it matches the XTC atom
       order exactly); the XTC provides the coordinates over time.
   * - ``xtc`` + ``pdb`` + ``keyfile``
     - As above, **plus** authoritative box dimensions, lattice spacing, hardwall
       flag, temperature and chain *types* taken from the keyfile.
   * - ``pdb`` only
     - A single frame (e.g. ``START.pdb``) - handy for inspecting a starting
       configuration.

A PDB is always required alongside an XTC (mdtraj needs a topology to read the
trajectory). Passing the keyfile is optional but recommended - without it, lemonade
*infers* what it can (see below).

Where the numbers come from
===========================

PIMMS writes coordinates in nanometres as ``lattice_index x spacing / 10``.
lemonade inverts that in a single vectorised step - ``round(nm / (spacing/10))`` -
recovering the exact integer lattice (the round-off is float32 noise). The
remaining metadata is resolved in this order:

* **spacing** - from the keyfile ``LATTICE_TO_ANGSTROMS``; otherwise PIMMS's default
  of ``3.65`` angstroms. Override with ``spacing=``.
* **dimensions** - from the keyfile ``DIMENSIONS``; otherwise inferred from the
  trajectory's box record. PIMMS writes a 2D system with a ``z`` period of exactly
  one lattice unit, and that is what marks a trajectory as 2D - a 3D
  configuration that happens to lie in the ``z = 0`` plane stays 3D. Override with
  ``dimensions=(x, y, z)``. A restart override or a resized-equilibration
  trajectory changes where this comes from - see
  :ref:`lemonade-effective-keyfile`.
* **hardwall** - from the keyfile ``HARDWALL``; otherwise ``False``. Override with
  ``hardwall=``. Again, a restart override or a resized-equilibration trajectory
  changes where this comes from - see :ref:`lemonade-effective-keyfile`.
* **temperature** - from the keyfile ``TEMPERATURE`` (needed only for surface
  tension). For a ``QUENCH_RUN`` keyfile PIMMS ignores ``TEMPERATURE`` and samples
  production at ``QUENCH_END``, so that is what is used (with a warning; frames
  written during the ramp were sampled at intermediate temperatures). Override with
  ``temperature=``.
* **topology** (chain lengths, sequences, bead types) - always from the PDB, since
  it is written in lockstep with the trajectory. Chains come from the PDB's
  ``TER`` blocks and beads from atom order; the atom serial and residue number
  columns are never read, so the duplicated serials of a 100,000+ bead system
  (PIMMS writes them modulo 100000) and the per-chain residue numbering make no
  difference. Without a keyfile, PIMMS PDB chain
  identifiers are preserved as chain-type labels (rather than guessing type from
  sequence); a PDB whose chain column is blank carries no type information at
  all, and there chains sharing a sequence are grouped into one type instead.
  PIMMS has 62 identifiers (``A-Z``, ``a-z``, ``0-9``); chain types past the
  62nd share the last one, so for such systems pass the keyfile. If a keyfile is
  given and its expanded ``CHAIN``/``EXTRA_CHAIN`` specification matches the PDB,
  those authoritative type labels are used; a mismatch emits a warning and
  retains the PDB labels. A *restart* keyfile is the exception - see below.

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

Selecting frames
================

You can subsample at load time (cheaper than loading everything and slicing after):

.. code-block:: python

   # every 5th frame between 100 and 400
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", start=100, stop=400, step=5)

   # thin an arbitrarily long run down to 200 evenly spaced frames
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", n_frames=200)

``start``/``stop``/``step`` are applied first and ``n_frames`` then thins what
survives, so the two compose. ``n_frames`` is a no-op if the trajectory is
already that short or shorter.

You can also slice *after* loading - ``traj[100:400:5]`` returns a new
``LatticeTrajectory`` over those frames that shares the parent's data, so it is
cheap and does not re-read anything. An integer index gives a ``Frame`` instead;
a slice or an index array gives a trajectory.

Checking the load
=================

``load`` never fails quietly on a suspicious input: it either raises, or warns.
The warnings fire on **every** load, regardless of ``verbose``:

* the **lattice round-off residual** - if the coordinates do not sit on the
  integer lattice at the given spacing (residual above 0.05, indicating a wrong
  or omitted ``LATTICE_TO_ANGSTROMS`` / ``spacing=``), a warning is raised, since
  the recovered lattice would be corrupted;
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
* the **keyfile/PDB chain mismatch** - a keyfile whose expanded
  ``CHAIN``/``EXTRA_CHAIN`` block does not reproduce the PDB's chains cannot be
  from the same run, so its chain types are dropped and the PDB labels kept.
  This does not apply to a keyfile with a ``RESTART_FILE``, whose ``CHAIN`` block
  is never used in the first place;
* the **62-identifier collision** - raised only when no keyfile resolved the
  types and the topology ends up with 62 or more of them. Nothing in the PDB
  says whether a merge actually happened, only that any type past the 62nd
  would have been folded into the last label, so lemonade says so and asks for
  the keyfile;
* the **quench temperature** - a ``QUENCH_RUN`` keyfile reports that
  ``QUENCH_END`` is being used in place of the (ignored) ``TEMPERATURE``; if it
  sets no ``QUENCH_END`` at all, it reports that no temperature could be
  recorded.

Errors, by contrast, are raised for the inputs lemonade cannot interpret at all:
an ``xtc`` without a ``pdb``, neither file given, a non-positive ``n_frames`` or
``spacing``, a box that is neither given, in the keyfile, nor recorded in the
trajectory, a non-boolean ``hardwall``, coordinates too large for int32, a
PDB whose bead count disagrees with the trajectory, a ``start``/``stop``/``step``
selection that keeps no frames (a window past the end of a short run), and a
keyfile that sets ``RESTART_OVERRIDE_HARDWALL`` / ``RESTART_OVERRIDE_DIMENSIONS``
whose ``RESTART_FILE`` cannot be found or read.

Pass ``verbose=True`` additionally for a one-line summary. It gains a trailing
``(WARNING lattice round-off ...)`` when the residual is not float32 noise
(above ``1e-3``):

.. code-block:: python

   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf",
                        verbose=True)
   # [lemonade] 101 frames, 250 chains, 2000 beads; box (30, 30, 30), spacing 3.65 A

The loaded object reports the basics directly:

.. code-block:: python

   traj.n_frames, traj.n_chains, traj.n_atoms
   traj.dimensions          # (30, 30, 30)   - 2 entries for a 2D system
   traj.n_dim               # 2 or 3
   traj.spacing             # 3.65
   traj.hardwall            # False
   traj.temperature         # 90.0  (or None if unknown)
   traj.sequences           # ['AAAA', 'AAAA', ...] one per chain
   traj.chain_types         # (n_chains,) int32 type label per chain
   traj.times               # frame times from the XTC, shape (n_frames,)

**Which frames are in the file.** Frame 0 of ``traj.xtc`` is always the starting
configuration. With the default ``SAVE_EQ : True`` the equilibration frames are
included too, so a run with ``EQUILIBRATION : 1000`` and ``XTC_FREQ : 100`` has
ten equilibration frames before the first production one; set ``SAVE_EQ : False``
in the keyfile (or slice, ``traj[11:]``) to analyse production only. ``traj.times``
and ``frame.time`` are the *frame indices* ``0, 1, 2, ...`` as written by PIMMS,
not Monte Carlo steps: frame ``k`` is step ``k * XTC_FREQ`` (with ``SAVE_EQ :
False``, frame 0 is still step 0 and frame ``k >= 1`` is the ``k``-th production
frame).

Coordinate, time and cached-analysis arrays are read-only. lemonade memoises
whole-chain coordinates and batched observables; immutability prevents an edit to
``traj.positions`` from leaving an already-cached Rg or COM describing a different
trajectory. Make an explicit ``.copy()`` when you need a mutable working array.
