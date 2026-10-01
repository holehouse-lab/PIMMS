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
     - read-only ``(n_frames, n_beads, 3)`` int32 array of wrapped lattice
       positions. The third column is always present and is zero in 2D.
   * - ``whole_positions()``
     - the same, with every chain made contiguous across periodic boundaries.
   * - ``radius_of_gyration()``
     - ``(n_frames, n_chains)`` float64 Rg of every chain in every frame.
   * - ``center_of_mass()``
     - ``(n_frames, n_chains, n_dim)`` float64 - note the last axis is ``n_dim``,
       not 3.
   * - ``asphericity()`` / ``end_to_end_distance()``
     - ``(n_frames, n_chains)`` float64 each.
   * - ``chain_types``
     - ``(n_chains,)`` int32 type label of each chain; ``sequences`` is the
       matching list of 1-letter bead sequences.
   * - ``topology``
     - the columnar chain/bead topology (CSR ``offsets``, ``lengths``,
       ``sequences``, ``chain_types``, ``alphabet``); rarely needed directly.
   * - ``store``
     - the backing ``TrajectoryStore``; every view above is an index into it.

Frame
=====

One snapshot. It is a sequence of polymers, and also gives you the clustering.

.. code-block:: python

   frame = traj[0]
   len(frame)                # number of chains  (== frame.n_chains)
   frame[3]                  # Polymer for chain 3
   frame.polymer(3)          # the same thing, spelled out
   frame.polymer(-1)         # negative indices count from the end
   frame.polymers            # list of every Polymer in the frame
   for polymer in frame:     # iterate over chains
       ...
   frame.index               # 0
   frame.time                # frame time (from the XTC)
   frame.n_chains, frame.n_beads
   frame.positions           # (n_beads, 3) wrapped positions this frame

Both ``frame[c]`` and ``frame.polymer(c)`` range-check the index and let negative
values count from the end; an out-of-range index raises ``IndexError``.
``frame.all_bead_positions`` is an alias of ``frame.positions``.

Clusters and the condensate:

.. code-block:: python

   frame.clusters            # list of Cluster, largest (most beads) first
   frame.droplet             # the largest cluster (or None if the frame is empty)
   frame.grid                # a dimensions-shaped int32 grid (site = chain index + 1)

``clusters`` is computed lazily the first time you ask for it and only for that
frame, so iterating over frames does not pay for clustering you do not use. The
expensive part - the connected-component search - is memoised on the
*trajectory*, not on the ``Frame``, because ``traj[f]`` mints a new ``Frame``
every time; a second analysis pass over the same frames therefore does not
repeat the search. The per-cluster geometry (single image, hull) stays on the
``Cluster`` objects, so it can still be garbage collected, and a fresh ``Frame``
recomputes it. ``grid`` is painted afresh on *every* access, so hold a
reference if you need it repeatedly.

Polymer
=======

A single chain within a single frame. Creating one is free (it just stores three
indices); its properties are computed on demand and cached.

.. code-block:: python

   p = traj[0][3]
   len(p)                     # number of beads
   p.chain_index, p.frame_index   # 3, 0
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

Note the shapes: position arrays always carry three columns (``z`` is zero in
2D), while ``center_of_mass`` is trimmed to ``n_dim``.

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

**Connectivity** is contact-based: two chains belong to the same cluster when any
bead of one is within Chebyshev distance 1 of any bead of the other (the full
26-site Moore shell in 3D, 8 sites in 2D), i.e. the same short-range shell the
Hamiltonian uses. Contacts are found across periodic boundaries, never through a
hardwall. Long-range (Chebyshev 2/3) pairs do **not** join clusters here; in the
``LR_*`` cluster files PIMMS writes during a run they do, when their ``LR`` /
``SLR`` interaction energy is nonzero.

**Centre-of-mass frames.** ``Cluster.center_of_mass`` is the mean of the cluster's
single-image (gathered) positions and can lie outside the box - it is congruent to
the in-box centre modulo the box length. ``Polymer.center_of_mass`` is likewise the
mean of the chain made whole, anchored at bead 0's wrapped position, so it too may
fall outside ``[0, L)``. Wrap with ``np.mod(com, dims)`` when an in-box point is
needed.

.. code-block:: python

   cl = traj[-1].clusters[0]           # the biggest cluster
   len(cl), cl.n_chains, cl.n_beads    # len() is the chain count
   for polymer in cl:                  # iterate over member chains
       ...
   cl.polymers                         # the same, as a list
   cl.chain_indices                    # the member chain indices in the frame

   cl.positions                        # raw positions of all beads (n_beads, 3) int32
   cl.single_image_positions()         # (n_beads, n_dim) float64, one periodic image
   cl.spanning_axes()                  # axes the cluster winds (or touches both walls of)

   cl.center_of_mass                   # (n_dim,)
   cl.radius_of_gyration
   cl.asphericity                      # unnormalised, as for a Polymer
   cl.sphericity                       # isoperimetric, ~1 for a sphere (3D)
   cl.volume, cl.surface_area, cl.density   # convex-hull based
   cl.radial_density_profile()         # fraction of each Chebyshev shell
                                       # holding this cluster's beads

   cl.chain_type_composition           # {type: count}
   cl.bead_type_composition            # {'A': count, 'B': count}

``cl.asphericity`` is the same unnormalised gyration-tensor quantity as a
polymer's (see :doc:`conformational`), so it is not the dimensionless number in
PIMMS's ``CLUSTER_ASPH.dat``; ``cl.radius_of_gyration`` is the ordinary Rg of
the gathered cluster.

``radial_density_profile()`` returns a plain list of shells **starting at shell
1**: entry ``k`` is the fraction of the lattice sites at Chebyshev distance ``k +
1`` from the cluster's centre of mass that are occupied by this cluster's beads.
The list always has ``min(dimensions) // 2 - 1`` entries (so every shell fits
inside the shortest box axis), zero-padded once every bead has been counted.
Beads of other clusters and dilute chains in the same shell are not counted, so
the profile falls to zero outside the cluster rather than to the dilute-phase
density; for the density of every bead about the condensate, use
:func:`~pimms.lemonade.phase_separation.radial_density_profile` (see
:doc:`phase_separation`). Shell 0 is the single
site the centre of mass rounds onto - one site out of one, which carries no
density information - so it is not returned; plot the profile against ``1, 2, 3,
...`` or the interface lands one lattice unit too close to the centre. This is
the same convention as the ``CLUSTER_RADIAL_DENSITY_PROFILE.dat`` file PIMMS
writes. Under ``HARDWALL`` a shell running through a wall is
normalised by the sites that actually lie inside the box, so a cluster wetting a
wall is not reported as artificially dilute. Pass
``minimum_cluster_size_in_beads=`` to skip clusters below a size, in which case
``None`` is returned instead of a profile.

.. note::

   ``single_image_positions`` and the convex-hull quantities assume a *compact*
   cluster. A cluster that **percolates** the box (e.g. a slab that spans the
   periodic plane) cannot be gathered into a single image; asking for one warns
   (``single-image gather: cluster percolates the periodic box on axis ...``) and
   the gathered coordinates that come back are search-order dependent.
   ``spanning_axes()`` lists the axes it winds, from the same test run once on
   the gathered image and cached with it; it is what the phase-separation and
   surface-tension routines use to leave such frames out. For those clusters,
   work from the wrapped ``positions`` and use the slab tools in
   :doc:`phase_separation`. The convex-hull ``volume`` / ``surface_area`` /
   ``density`` return ``-1`` for degenerate (too small, collinear, coplanar)
   clusters, and ``sphericity`` is ``nan`` there.

   Under ``HARDWALL`` the box is not periodic: clustering never connects chains
   through the walls, and ``single_image_positions()`` simply returns the raw
   positions unchanged (they already form a single Cartesian image);
   ``spanning_axes()`` then lists the axes on which the cluster touches both
   walls.

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
