.. _lemonade-hierarchy:

======================
The object hierarchy
======================

lemonade mirrors the way you think about a simulation: a **trajectory** is a
sequence of **frames**, a frame contains **polymers** (and **clusters** of
polymers), and a polymer is a chain of beads. Each level is indexable, iterable and
has a small, predictable set of attributes.

.. code-block:: text

   LatticeTrajectory --[i]--> Frame --[c]--> Polymer
                                    '--clusters--> Cluster --> Polymer

LatticeTrajectory
=================

The top-level object. It is a sequence of frames.

.. code-block:: python

   len(traj)                 # number of frames
   traj[0]                   # first Frame
   traj[-1]                  # last Frame
   traj[10:50:2]             # a sub-trajectory (a new LatticeTrajectory)
   for frame in traj:        # iterate over frames
       ...

Besides the metadata from :doc:`loading`, it exposes the raw arrays and the
whole-trajectory analyses (see :doc:`conformational`):

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Member
     - Meaning
   * - ``positions``
     - ``(n_frames, n_atoms, 3)`` int array of wrapped lattice positions.
   * - ``whole_positions()``
     - the same, with every chain made contiguous across periodic boundaries.
   * - ``radius_of_gyration()``
     - ``(n_frames, n_chains)`` Rg of every chain in every frame.
   * - ``center_of_mass()``
     - ``(n_frames, n_chains, n_dim)``.
   * - ``asphericity()`` / ``end_to_end_distance()``
     - ``(n_frames, n_chains)`` each.
   * - ``topology``
     - the :class:`~pimms.lemonade.LatticeTrajectory` topology (offsets, sequences,
       types); rarely needed directly.

Frame
=====

One snapshot. It is a sequence of polymers, and also gives you the clustering.

.. code-block:: python

   frame = traj[0]
   len(frame)                # number of chains
   frame[3]                  # Polymer for chain 3
   for polymer in frame:     # iterate over chains
       ...
   frame.index               # 0
   frame.time                # frame time (from the XTC)
   frame.positions           # (n_atoms, 3) wrapped positions this frame

Clusters and the condensate:

.. code-block:: python

   frame.clusters            # list of Cluster, largest (most beads) first
   frame.droplet             # the largest cluster (or None if the frame is empty)
   frame.grid                # a dimensions-shaped int grid (site = chain index + 1)

``clusters`` and ``grid`` are computed lazily the first time you ask for them (and
only for that frame), so iterating over frames does not pay for clustering you do
not use.

Polymer
=======

A single chain within a single frame. Creating one is free (it just stores three
indices); its properties are computed on demand and cached.

.. code-block:: python

   p = traj[0][3]
   len(p)                     # number of beads
   p.sequence                 # 'AABBAABB'
   p.chain_type               # integer type label
   p.positions                # (L, 3) wrapped integer positions (a view)
   p.whole_positions          # (L, 3) made contiguous across PBC

   p.radius_of_gyration       # scalar
   p.center_of_mass           # (n_dim,)
   p.asphericity
   p.end_to_end_distance
   p.straddles_boundary       # True if the chain crosses a periodic face

   p.distance_map()           # (L, L) inter-bead distance matrix
   p.internal_scaling()       # (separations, mean_distance)

All the scalar conformational properties are read straight from the trajectory's
batched arrays, so ``traj[f][c].radius_of_gyration`` and
``traj.radius_of_gyration()[f, c]`` are the same number.

Cluster
=======

A connected group of polymers - the natural unit for condensate analysis. You get
clusters from a frame; they are sorted **largest first, by bead count**. (Bead count,
not chain count: in a multi-component system with chains of different lengths the
cluster with the most chains is not necessarily the one with the most material, and
it is the material that makes a condensate.)

.. code-block:: python

   cl = traj[-1].clusters[0]           # the biggest cluster
   cl.n_chains, cl.n_beads
   for polymer in cl:                  # iterate over member chains
       ...
   cl.chain_indices                    # the member chain indices in the frame

   cl.positions                        # raw positions of all beads (n_beads, 3)
   cl.single_image_positions()         # gathered into one periodic image

   cl.center_of_mass
   cl.radius_of_gyration
   cl.asphericity
   cl.sphericity                       # isoperimetric, ~1 for a sphere (3D)
   cl.volume, cl.surface_area, cl.density   # convex-hull based
   cl.radial_density_profile()

   cl.chain_type_composition           # {type: count}
   cl.bead_type_composition            # {'A': count, 'B': count}

.. note::

   ``single_image_positions`` and the convex-hull quantities assume a *compact*
   cluster. A cluster that **percolates** the box (e.g. a slab that spans the
   periodic plane) cannot be gathered into a single image; for those, work from the
   wrapped ``positions`` and use the slab tools in :doc:`phase_separation`. The
   convex-hull ``volume`` / ``surface_area`` return ``-1`` for degenerate (too
   small, coplanar) clusters, and ``sphericity`` is ``nan`` there.

A worked example
================

Mean radius of gyration over the second half of a run, and the size of the largest
cluster in the final frame:

.. code-block:: python

   import numpy as np
   import pimms.lemonade as lemonade

   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

   rg = traj.radius_of_gyration()             # (n_frames, n_chains)
   mean_rg = rg[traj.n_frames // 2:].mean()   # average over time and chains

   biggest = traj[-1].droplet
   print(f"<Rg> = {mean_rg:.2f} lattice units")
   print(f"largest cluster: {biggest.n_chains} chains, "
         f"{biggest.n_beads} beads, Rg = {biggest.radius_of_gyration:.2f}")
