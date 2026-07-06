.. _lemonade-conformational:

========================
Conformational analysis
========================

lemonade computes the standard single-chain observables - radius of gyration,
centre of mass, asphericity, end-to-end distance, distance maps and internal
scaling - either for the whole trajectory at once or for a single chain.

Whole trajectory, at once
=========================

The trajectory-level methods return an array covering every chain in every frame,
computed in a handful of vectorised operations (no Python loop over chains or
frames). This is the fast path and what you want for ensemble averages:

.. code-block:: python

   rg  = traj.radius_of_gyration()      # (n_frames, n_chains)
   com = traj.center_of_mass()          # (n_frames, n_chains, n_dim)
   asph = traj.asphericity()            # (n_frames, n_chains)
   ete  = traj.end_to_end_distance()    # (n_frames, n_chains)

From there, ordinary numpy gives you whatever average you need:

.. code-block:: python

   import numpy as np

   equil = traj.n_frames // 2                 # discard the first half
   rg_mean_per_chain = rg[equil:].mean(axis=0)     # (n_chains,)
   rg_mean = rg[equil:].mean()                     # scalar ensemble average
   rg_by_type = {t: rg[equil:][:, traj.chain_types == t].mean()
                 for t in np.unique(traj.chain_types)}

Single chain
============

Navigating to a :class:`~pimms.lemonade.Polymer` gives the same quantities as
scalars, plus per-chain matrices. The scalars are read from the batched arrays, so
they agree exactly with the trajectory-level results.

.. code-block:: python

   p = traj[0][3]
   p.radius_of_gyration          # == traj.radius_of_gyration()[0, 3]
   p.end_to_end_distance
   p.center_of_mass

   dm = p.distance_map()         # (L, L) inter-bead Euclidean distances
   sep, dist = p.internal_scaling()   # mean distance vs sequence separation |i - j|

``internal_scaling`` returns the separations ``1 .. L-1`` and the mean bead-bead
distance at each - the lattice analogue of the internal-scaling profile used to read
off polymer scaling exponents.

How the numbers are defined
===========================

All conformational quantities are computed on **whole** positions - each chain is
first made contiguous across periodic boundaries (a compiled kernel does this for
the whole trajectory at once), so a chain that straddles a box face is measured as
one connected object rather than being torn in two.

* **Radius of gyration** is :math:`R_g = \sqrt{\langle |{\bf r}_i - {\bf r}_{cm}|^2
  \rangle}` over the chain's beads (equivalently the square root of the trace of the
  gyration tensor), using the exact (floating-point) centre of mass.
* **Asphericity** comes from the gyration-tensor eigenvalues:
  :math:`\lambda_3 - \tfrac12(\lambda_1 + \lambda_2)` in 3D (zero for an isotropic
  coil), :math:`\lambda_2 - \lambda_1` in 2D.
* **Centre of mass** is the mean of the whole positions; ``distance_map`` and
  ``end_to_end_distance`` are ordinary Euclidean distances on those contiguous
  coordinates.

.. note::

   **Agreement with PIMMS's own Rg.** For chains that are small compared with the
   box, lemonade's Rg matches PIMMS's built-in ``get_polymeric_properties``. For a
   chain larger than roughly half the box the two *intentionally* differ: lemonade
   measures the whole (contiguous) chain, whereas PIMMS's on-the-fly analysis uses
   the minimum-image convention, which collapses such a chain (a finite-size
   artefact PIMMS warns about during the run). If you see a large discrepancy, your
   box is small relative to your chains.

Positions
=========

Two views of the coordinates are available at every level:

* ``positions`` - integer lattice coordinates **wrapped into the box** (``0 <= x <
  dimensions``). Use these for anything grid-based (occupancy, clustering, slab
  profiles).
* ``whole_positions`` - the chain **unwrapped** so consecutive beads never jump a
  boundary (coordinates may fall outside the box). Use these for shape and distance
  measurements. This is what all the conformational analyses use internally.

.. code-block:: python

   p = traj[0][3]
   p.straddles_boundary          # does this chain cross a periodic face?
   p.positions                   # wrapped
   p.whole_positions             # contiguous

At the trajectory level, ``traj.positions`` and ``traj.whole_positions()`` give the
same two views for the entire run as ``(n_frames, n_atoms, 3)`` arrays.
