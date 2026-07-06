.. _restart-files:

=============
Restart files
=============

A **restart file** is a snapshot of a simulation's configuration that can be used
to seed a new simulation. Restart files let you checkpoint long runs, continue a
simulation under different conditions, or build up a complex system in stages.

.. _restart-how:

How restart files work
======================

During a run PIMMS periodically writes its state to ``restart.pimms`` (a Python
pickle). How often is controlled by ``RESTART_FREQ``, which is either an integer
step frequency or the default sentinel ``"Every 10th-percentile"`` (write at 10%,
20%, ... 100% of ``N_STEPS``).

The file stores exactly what is needed to reconstruct the configuration:

.. code-block:: python

   {
     'DIMENSIONS': [x, y, z],          # box size (length 2 for 2D)
     'HARDWALL'  : True | False,       # boundary condition the snapshot was taken under
     'ENERGY'    : <float>,            # total energy at write time
     'CHAINS'    : {                   # one entry per chain
         chainID: [positions, sequence, chainType],
         ...
     },
   }

where, for each chain, ``positions`` is the list of bead coordinates (in N→C
order), ``sequence`` is the one-letter sequence string, and ``chainType`` is an
integer grouping identical chains.

What a restart file **does** capture is the full spatial configuration (every bead
of every chain), the box dimensions, and the boundary condition. What it **does
not** capture is the move statistics, accumulated analysis, temperature, or random
seed - a restarted run begins fresh in those respects.

.. _restart-using:

Using a restart file as a starting configuration
================================================

Point the ``RESTART_FILE`` keyword at a ``restart.pimms`` file:

.. code-block:: text

   RESTART_FILE : restart.pimms

When ``RESTART_FILE`` is set, the chains come from the restart file, so the
``CHAIN`` keyword is **not** required (any ``CHAIN`` lines are ignored). You still
provide the run-control keywords (``TEMPERATURE``, ``N_STEPS``, ``EQUILIBRATION``,
the ``MOVE_*`` set, analysis keywords, and a ``PARAMETER_FILE`` whose bead types
cover the restart's sequences).

This makes a simple continuation trivial: take the ``restart.pimms`` from one
simulation and start another from it, optionally at a different temperature or
with a different move mix.

.. _restart-complexities:

Complexities: dimensions, hardwall and when restarts are valid
==============================================================

Restarting into a *different* box or boundary condition is supported, but with
rules - PIMMS will refuse combinations that could place beads illegally. By default
the new run uses the ``DIMENSIONS`` and ``HARDWALL`` from your **keyfile** and
reconciles them against the snapshot. The two ``RESTART_OVERRIDE_*`` keywords do the
opposite: they tell PIMMS to **ignore those keyfile values and inherit them from the
restart file** instead.

**Dimensions.** By default the run uses the keyfile ``DIMENSIONS``:

* For a **hardwall** snapshot, the keyfile box may be **equal to or larger than**
  the snapshot's. If it is larger, PIMMS grows the box and re-centres the
  configuration inside it - so growing into a bigger box needs **no** override, just
  set ``DIMENSIONS`` to the larger box. The box can never be made *smaller* than the
  snapshot (that could force overlaps).
* For a **periodic (PBC)** snapshot, the keyfile ``DIMENSIONS`` must match the
  snapshot **exactly** (changing a periodic box would break the wrapping).
* The dimensionality must always match - you cannot turn a 2D restart into a 3D run.

Setting ``RESTART_OVERRIDE_DIMENSIONS : True`` **ignores the keyfile** ``DIMENSIONS``
and adopts the snapshot's box exactly as it was. Use it to continue in the original
box without having to repeat its size in the keyfile; it does *not* grow the box (and
it is incompatible with ``RESIZED_EQUILIBRATION``).

**Hardwall.** By default the run uses the keyfile ``HARDWALL``, and PIMMS checks the
transition is legal:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Snapshot taken under
     - New run requests
     - Allowed?
   * - ``HARDWALL : True``
     - ``HARDWALL : False`` (PBC) or ``True``
     - **Yes** - hardwall chains never cross a boundary, so they are valid either way.
   * - ``HARDWALL : False`` (PBC)
     - ``HARDWALL : False`` (PBC)
     - **Yes** - unchanged boundary.
   * - ``HARDWALL : False`` (PBC)
     - ``HARDWALL : True``
     - **No** - PBC chains may already wrap across a boundary, which a hard wall forbids.

Setting ``RESTART_OVERRIDE_HARDWALL : True`` **ignores the keyfile** ``HARDWALL`` and
adopts the snapshot's boundary condition - a convenience for continuing a run under
the same boundaries it was generated with.

**Box-size transitions at a glance:**

.. list-table::
   :header-rows: 1
   :widths: 45 20 35

   * - Transition
     - Allowed?
     - Notes
   * - 30³ → 50³ (grow)
     - **Yes**
     - Set ``DIMENSIONS : 50 50 50`` from a hardwall original run; **no** override
       needed (``RESTART_OVERRIDE_DIMENSIONS`` would instead force the box back to 30³).
   * - 30³ → 30³ (same)
     - **Yes**
     - The default; no override needed.
   * - 30³ → 20³ (shrink)
     - **No**
     - Smaller boxes are not supported.

Resizing on restart also interacts with ``RESIZED_EQUILIBRATION`` - a PBC restart
file is incompatible with ``RESIZED_EQUILIBRATION``/``RESTART_OVERRIDE_DIMENSIONS``
because those features assume a hardwall, growable box.

.. _restart-extra-chains:

Adding new chains: EXTRA_CHAIN
==============================

``EXTRA_CHAIN`` lets you add chains that were **not** present in the restart file -
for example, to titrate a second component into a pre-equilibrated condensate. It
uses the same syntax as ``CHAIN`` and may be repeated:

.. code-block:: text

   RESTART_FILE : restart.pimms
   EXTRA_CHAIN  : 50 EEEEEEEE      # add 50 copies of an 8-bead chain
   EXTRA_CHAIN  : 10 KKKK          # ...and 10 more of another type

The new chains are inserted at random positions that do not overlap the existing
configuration, on top of the restart chains. Because this can be repeated, you can
build a system up in stages - equilibrate component A, restart and add component
B, restart again and add component C, and so on. ``EXTRA_CHAIN`` requires a
``RESTART_FILE`` (there must be an existing configuration to add to).
