.. _lemonade:

========================
Analysis with lemonade
========================

**lemonade** is the analysis backend that ships inside PIMMS
(``pimms.lemonade``). It loads a finished PIMMS trajectory - an ``XTC`` of
coordinates, a ``PDB`` of topology, and the run's *keyfile* - and turns it into a
navigable, object-oriented representation of the simulation that you can query for
conformational, cluster and phase-separation properties.

It is built for two things at once:

* **A clean object model.** A trajectory is a sequence of frames; a frame is a set
  of polymers (and clusters of polymers); a polymer is a chain of beads. You walk
  that hierarchy the way you think about it -
  ``trajectory[frame][chain].radius_of_gyration``.
* **Speed.** The whole trajectory is held in contiguous arrays, coordinates are
  converted back to the integer lattice in one vectorised step, periodic-boundary
  unwrapping runs in a compiled kernel, and the standard per-chain analyses are
  computed for every chain in every frame at once. A hundred-frame, few-hundred-chain
  trajectory loads and analyses in well under a second.

.. note::

   lemonade contains a compiled kernel, so the package must be **built** before use
   (:doc:`/installation` covers this; from the repo root it is ``./build.sh``). If
   ``import pimms.lemonade`` fails with a missing-extension error, rebuild.

Quickstart
==========

.. code-block:: python

   import pimms.lemonade as lemonade

   # load a finished run (keyfile is optional but recommended)
   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

   traj                                 # <LatticeTrajectory 101 frames, 250 chains, ...>

   # --- whole-trajectory analyses (vectorised: shape (n_frames, n_chains)) ---
   rg = traj.radius_of_gyration()       # radius of gyration of every chain, every frame
   com = traj.center_of_mass()          # (n_frames, n_chains, n_dim)

   # --- navigate to a single object ---
   frame   = traj[0]                    # a Frame
   polymer = frame[3]                   # chain 3 in frame 0 (a Polymer)
   polymer.radius_of_gyration           # a scalar
   polymer.sequence                     # e.g. 'AABBAABB'
   polymer.whole_positions              # the chain made contiguous across PBC

   # --- clusters and condensates ---
   for cluster in frame.clusters:       # connected-component clusters, most beads first
       cluster.n_beads, cluster.radius_of_gyration, cluster.volume

   # --- phase separation ---
   from pimms.lemonade import phase_separation as ps
   result = ps.analyze(traj)            # binodal, condensed fraction, droplet shape

The object model
================

Loading returns a :class:`~pimms.lemonade.LatticeTrajectory`. Everything else is a
lightweight *view* onto it:

.. code-block:: text

   LatticeTrajectory                 the whole run
    |
    |-- [i] ------------> Frame       one snapshot (time point)
    |                      |
    |                      |-- [c] --> Polymer     one chain in that frame
    |                      |
    |                      '-- clusters --> Cluster --> Polymer
    |                                        (a connected group of chains)
    '-- iterate ---------> Frame, Frame, ...

* A :class:`~pimms.lemonade.LatticeTrajectory` is indexable and iterable over
  **frames**; slicing it (``traj[::2]``) returns another trajectory over that range,
  sharing the underlying data.
* A :class:`~pimms.lemonade.Frame` is indexable and iterable over **polymers**, and
  exposes ``clusters`` and ``droplet`` (the largest cluster).
* A :class:`~pimms.lemonade.Polymer` is a chain within a frame - positions plus
  cached conformational properties (Rg, COM, asphericity, ...).
* A :class:`~pimms.lemonade.Cluster` is a connected group of polymers with condensate
  geometry (Rg, volume, density, radial profile, ...).

These views allocate nothing per bead: the positions live in one big array in the
trajectory, and a ``Frame`` / ``Polymer`` / ``Cluster`` is just a set of indices into
it. That is what keeps navigation cheap even for large trajectories.

Units and conventions
======================

* **Positions** are integer lattice coordinates (the box is ``traj.dimensions``).
  ``.positions`` are wrapped into the box; ``.whole_positions`` are unwrapped so each
  chain is spatially contiguous across periodic boundaries.
* **Lengths** (Rg, distances, radii, interface widths) are in **lattice units**. The
  physical spacing is ``traj.spacing`` (angstroms), so multiply by it for angstroms.
* **Densities** are occupied lattice-site fractions in ``[0, 1]``.
* **Temperature** is PIMMS's ``TEMPERATURE`` and equals :math:`k_B T` (PIMMS uses
  :math:`\exp(-\Delta E / T)` with :math:`k_B = 1`), so surface tension comes out in
  reduced units (interaction energy per lattice area).

Read on
=======

.. toctree::
   :maxdepth: 2

   loading
   hierarchy
   conformational
   phase_separation
   reference
