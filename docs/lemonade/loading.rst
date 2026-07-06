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
  trajectory's box (and flattened to 2D if every ``z`` is zero). Override with
  ``dimensions=(x, y, z)``.
* **hardwall** - from the keyfile ``HARDWALL``; otherwise ``False``. Override with
  ``hardwall=``.
* **temperature** - from the keyfile ``TEMPERATURE`` (needed only for surface
  tension). Override with ``temperature=``.
* **topology** (chain lengths, sequences, bead types) - always from the PDB, since
  it is written in lockstep with the trajectory. If a keyfile is given, the per-chain
  *type* labels are taken from its ``CHAIN`` order.

Whatever the trajectory was written with - wrapped positions, or whole molecules via
``TRAJECTORY_PBC_UNWRAP`` - lemonade canonicalises positions back into the box on
load and re-derives whole chains itself, so results do not depend on how the run was
saved.

Selecting frames
================

You can subsample at load time (cheaper than loading everything and slicing after):

.. code-block:: python

   # every 5th frame between 100 and 400
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", start=100, stop=400, step=5)

   # thin an arbitrarily long run down to ~200 evenly spaced frames
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", n_frames=200)

You can also slice *after* loading - ``traj[100:400:5]`` returns a new
``LatticeTrajectory`` over those frames that shares the parent's data, so it is
cheap and does not re-read anything.

Checking the load
=================

Pass ``verbose=True`` for a one-line summary, including the lattice round-off
residual (a warning appears if it is not essentially zero, which would indicate a
spacing mismatch):

.. code-block:: python

   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf",
                        verbose=True)
   # [lemonade] 101 frames, 250 chains, 2000 beads; box (30, 30, 30), spacing 3.65 A

The loaded object reports the basics directly:

.. code-block:: python

   traj.n_frames, traj.n_chains, traj.n_atoms
   traj.dimensions          # (30, 30, 30)
   traj.n_dim               # 2 or 3
   traj.spacing             # 3.65
   traj.hardwall            # False
   traj.temperature         # 90.0  (or None if unknown)
   traj.sequences           # ['AAAA', 'AAAA', ...] one per chain
   traj.times               # mdtraj frame times, shape (n_frames,)
