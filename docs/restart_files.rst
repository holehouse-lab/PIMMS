.. _restart-files:

=============
Restart files
=============

A **restart file** is a snapshot of a simulation's configuration that can be used
to seed a new simulation. Restart files let you checkpoint long runs, continue a
simulation under different conditions, build up a complex system in stages, or,
with ``RESTART_CONTINUE``, resume a stopped run *exactly*, as though it had never
stopped (:ref:`restart-continue`).

.. _restart-how:

How restart files work
======================

During a run PIMMS periodically writes its state to ``restart.pimms`` (a Python
pickle). How often is controlled by ``RESTART_FREQ``, which is either a positive
integer step frequency or the default sentinel ``"Every 10th-percentile"`` (write
every ``N_STEPS // 10`` steps, or every step for a run shorter than 10 steps,
**once production begins** - restart snapshots, like all
analysis, are suppressed during equilibration). A restart of the final state is
always written when the run completes, whatever the frequency. Each write goes to
a temporary file that is then renamed into place, so a crash part-way through a
checkpoint cannot destroy the previous good ``restart.pimms``.

.. warning::

   Restart files are Python **pickles**, and unpickling executes arbitrary code
   embedded in the file. Only load restart files you created yourself or
   obtained from a trusted source - never a file downloaded from an untrusted
   location. (The reader validates the *contents* - box bounds, overlaps,
   structure - but validation happens after unpickling, which is where the risk
   lives.)

The file stores what is needed to reconstruct the configuration, plus what is
needed to resume the run that wrote it:

.. code-block:: python

   {
     'DIMENSIONS'   : [x, y, z],          # box size (length 2 for 2D)
     'HARDWALL'     : True | False,       # boundary condition the snapshot was taken under
     'ENERGY'       : <float>,            # total energy at write time
     'CHAINS'       : {                   # one entry per chain
         chainID: [positions, sequence, chainType],
         ...
     },
     'STEP'         : <int>,              # the master step the snapshot was taken at
     'TEMPERATURE'  : <float>,            # the temperature in force at that step
     'RNG_PYTHON'   : <tuple>,            # random.getstate() at that moment
     'RNG_NUMPY'    : <tuple>,            # numpy.random.get_state() at that moment
     'PIMMS_VERSION': '<version>',
   }

where, for each chain, ``positions`` is the list of bead coordinates (in N→C
order), ``sequence`` is the one-letter sequence string, and ``chainType`` is an
integer grouping identical chains. ``ENERGY`` is recorded for reference only - the
new run recomputes its energy from scratch and never reads it.

The first four keys are present in every restart file, and the configuration
they hold (``ENERGY`` aside) is all that a plain restart (``RESTART_FILE``
without ``RESTART_CONTINUE``) uses: the new run begins fresh in every other
respect, with its own step count, seed, temperature and statistics. The last
five are the *continuation state*, written by every checkpoint since PIMMS 1.0.8
and acted on only by ``RESTART_CONTINUE``, which restores the step, the two
generators and the temperature so that the resumed run reproduces the
uninterrupted one (:ref:`restart-continue`); ``PIMMS_VERSION`` is kept for the
record, and the step is also quoted in ``keyfile_used.kf``. Restart files from
earlier versions lack them, still load as configurations, and are refused for
exact continuation with a message saying why. Move statistics and the
accumulated end-of-run analyses are not stored in either case: they are running
totals of the process that produced them.

The reader validates the whole file before it changes anything, so a bad restart
fails immediately with a message saying what is wrong (naming the offending
chain for most chain-level problems) rather than corrupting
a run: the top-level object must be a dictionary with the four configuration
keys (the continuation keys are optional, but if present must be well-formed:
a non-negative integer step, a positive finite temperature, both generator
states together and each a tuple), dimensions
must be 2 or 3 positive integers, ``HARDWALL`` a boolean and ``ENERGY`` finite,
every chainID a positive integer (0 is the solvent sentinel in the occupancy
grid), every chainType a non-negative integer, every sequence non-empty and the
same length as its position list, every coordinate an integer, every
bead inside the box, no two beads on the same site, and every pair of consecutive
beads a Moore (Chebyshev-1, so diagonals count) neighbour pair under the stored
boundary mode (so under periodic boundaries a bond may cross a face, and under a
hard wall it may not). Once the file has loaded, the keyfile checks also refuse
a snapshot in which two chains of the same chainType have different sequences,
and, as for any run, the parameter file must define every bead type the
snapshot uses. Chains are loaded in
ascending chainID order whatever order the pickle happens to hold them in, so the
per-chain analysis columns and the trajectory bead order always agree.

.. _restart-using:

Using a restart file as a starting configuration
================================================

Point the ``RESTART_FILE`` keyword at a ``restart.pimms`` file:

.. code-block:: text

   RESTART_FILE : restart.pimms

When ``RESTART_FILE`` is set, the chains come from the restart file, so the
``CHAIN`` keyword is **not** required (any ``CHAIN`` lines are ignored). It is the
only required keyword a restart file replaces: you still provide ``DIMENSIONS``,
``TEMPERATURE``, ``N_STEPS``, ``EQUILIBRATION``, the ``MOVE_*`` set, any analysis
keywords, and a ``PARAMETER_FILE`` whose bead types cover the restart's sequences.

Used this way a restart file is a *starting configuration*: the new run begins
from the snapshot's positions but is otherwise a new run. It numbers its steps
from 1, draws its own seed (or uses the ``SEED`` you give it), equilibrates for
its own ``EQUILIBRATION``, and in a quench starts the ramp from ``QUENCH_START``.
That is what you want for building a system in stages, changing the temperature
or the move mix, or fanning several independent replicas out from one
checkpoint, and it is the default. It is **not** the run that would have
happened had the original not stopped; for that, see
:ref:`Resuming a run exactly <restart-continue>` below.

.. warning::

   Run the new segment in a **fresh directory** (copying ``restart.pimms``
   and the input files across). Every run start wipes the previous run's
   outputs from its working directory, so restarting in place destroys the
   previous segment's ``.dat`` files and trajectory. Without
   ``RESTART_CONTINUE`` the new run's step numbers also restart from 1, so
   segments must be concatenated with that in mind.

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
  the snapshot's. If it is larger, PIMMS grows the box by placing the snapshot's
  box in the middle of the new one (every bead is shifted by half the difference
  in each axis, rounded down), so the configuration keeps its place within the
  old box rather than being centred itself - and growing into a bigger box needs
  **no** override, just set ``DIMENSIONS`` to the larger box. The box can never
  be made *smaller* than the snapshot (that could force overlaps).
* For a **periodic (PBC)** snapshot, the keyfile ``DIMENSIONS`` must match the
  snapshot **exactly** (changing a periodic box would break the wrapping).
* The dimensionality must always match - you cannot turn a 2D restart into a 3D run.

Setting ``RESTART_OVERRIDE_DIMENSIONS : True`` **ignores the keyfile** ``DIMENSIONS``
and adopts the snapshot's box exactly as it was. ``DIMENSIONS`` is a required keyword
and must still be present (any valid box will do - it is replaced wholesale,
number of dimensions included); the override guarantees the run
continues in the original box whatever the keyfile says. It does *not* grow the box
(and it is incompatible with ``RESIZED_EQUILIBRATION``). The box inherited from the
snapshot is re-checked against the usual rules - every axis must still be at least
7 lattice units - so an override cannot smuggle an unusable box past the keyfile
checks.

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
the same boundaries it was generated with. It is only a shortcut for setting
``HARDWALL`` yourself: it does not relax any of the other rules. A periodic
snapshot is still incompatible with ``RESIZED_EQUILIBRATION``, and the rule that
``MOVE_CLUSTER_ROTATE`` needs a cubic/square box under periodic boundaries (see
``DIMENSIONS``) is applied to whichever ``HARDWALL`` and ``DIMENSIONS`` the run
finally ends up with, not to the keyfile's values.

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
file is incompatible with ``RESIZED_EQUILIBRATION`` (that feature assumes a
hardwall, growable box), whether or not ``RESTART_OVERRIDE_HARDWALL`` is set.
``RESTART_OVERRIDE_DIMENSIONS`` is accepted for a PBC restart: it simply adopts the
restart file's box.

Combining a **hardwall** restart with ``RESIZED_EQUILIBRATION`` is supported: the
equilibration box, not ``DIMENSIONS``, is what the snapshot is reconciled against
for the equilibration phase, so ``RESIZED_EQUILIBRATION`` must be at least as large
as the snapshot's box in every axis (if it is larger, the snapshot's box is placed
in its middle in the same way), and ``DIMENSIONS`` is the box the run grows into
afterwards. So the ordering is snapshot box <= ``RESIZED_EQUILIBRATION`` <=
``DIMENSIONS``, and
``RESTART_OVERRIDE_DIMENSIONS`` must be False.

.. _restart-continue:

Resuming a run exactly: RESTART_CONTINUE
========================================

Everything above treats a restart file as a *configuration*: a place to start a
new run from. A run started that way numbers its steps from 1, draws a fresh seed,
and (in a quench) starts the ramp again from ``QUENCH_START``. That is the right
behaviour for building a system in stages or fanning out replicas from one
checkpoint, and it is the default. It is not a continuation: the segments of a
run that was stopped and restarted this way are not the run that would have
happened had it not been stopped.

For that there is ``RESTART_CONTINUE : True``. A restart file written by PIMMS
1.0.8 or later also records the step it was written at, the state of both global
random-number generators at that moment, and the temperature in force, and
``RESTART_CONTINUE`` restores all three, so the resumed segment makes exactly the
moves the uninterrupted run would have made. Concretely:

* the step counter continues from the checkpoint's step, so the rows the resumed
  segment writes carry the same step numbers the uninterrupted run would have
  given them and ``N_STEPS`` is the **total** length of the run, which must
  exceed the checkpoint's step;
* the generators are restored as the last thing before the master loop, so every
  random draw from the first resumed step on is the one the uninterrupted run
  would have made, including the per-megamove seeds handed to the compiled
  kernels (which are drawn from the Python generator);
* the temperature is restored to its value at the checkpoint and the TSMMC
  coordinator, if any, is rebuilt on it, so a quench ramp picks up where it was
  rather than starting again;
* frame 0 of the resumed segment's trajectory is the checkpoint configuration,
  and the frames after it fall on the same ``XTC_FREQ`` multiples as before;
* ``ACCEPTANCE.dat``, ``MOVE_FREQS.dat``, ``TOTAL_MOVES.dat`` and the accumulated
  end-of-run analyses (internal scaling, distance maps) count from the start of
  the segment, since they are running totals of *this* process.

A continuation has to be the same system, so the keyfile is checked: ``SEED``
must not be given (the generator state comes from the file, and a seed would
contradict it), ``RESIZED_EQUILIBRATION`` and ``EXTRA_CHAIN`` are refused, and
``DIMENSIONS`` and ``HARDWALL`` must equal the restart file's, which
``RESTART_OVERRIDE_DIMENSIONS : True`` and ``RESTART_OVERRIDE_HARDWALL : True``
guarantee without your having to copy the values across. ``TEMPERATURE`` must equal
the checkpoint's temperature; in a ``QUENCH_RUN``, where ``TEMPERATURE`` is
replaced by ``QUENCH_START`` anyway, the checkpoint's temperature must instead lie
between ``QUENCH_START`` and ``QUENCH_END`` (only that range is checked, so to
reproduce the uninterrupted run give exactly the quench settings of the run
being resumed). The reason is that the Hamiltonian's
``ANGLE_PENALTY_T_NORM`` scaling is
built from the keyfile, so a different value would give a run that is neither the
original nor the one asked for. A restart file with no
continuation state (one written by an earlier PIMMS, or by a
``RestartObject`` built outside a run) is refused with a message saying so; it
can still seed a new run in the ordinary way. Every problem found is listed in
the one refusal, each with what to change. (Some combinations never reach
these checks: a periodic checkpoint with ``HARDWALL : True`` or with
``RESIZED_EQUILIBRATION`` is refused first by the general restart rules above.)
Run each segment in its own directory, as always.

.. code-block:: text

   RESTART_FILE                 : restart.pimms
   RESTART_CONTINUE             : True
   RESTART_OVERRIDE_DIMENSIONS  : True
   RESTART_OVERRIDE_HARDWALL    : True
   N_STEPS                      : 200000      # the run's total length, not the segment's

The guarantee is tested by running a system to completion, resuming a second
copy from a mid-run checkpoint, and demanding identical ``ENERGY.dat``,
``RG.dat`` and ``QUENCH.dat`` rows and identical final positions, under periodic
and hardwall boundaries, through a quench ramp, with TSMMC excursions and with
``PARALLELIZE`` on.

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
configuration, on top of the restart chains, and take chainIDs that follow on from
the highest chainID in the restart file. An ``EXTRA_CHAIN`` whose sequence already
exists in the snapshot joins that existing chain type rather than creating a new
one (the sequence is upper-cased first unless ``CASE_INSENSITIVE_CHAINS : False``,
so ``b`` joins a type of ``B``); a new sequence gets the next free chain type. A
chain that joins an existing type keeps its new, higher chainID, so unless that
type already held the highest chainIDs its chains are no longer contiguous in
chainID order - and hence in the per-chain analysis columns,
``chain_to_chainid.txt`` and ``START.pdb`` - although they share the type's PDB
chain identifier and its ``CHAIN_<TYPE>_`` files. ``keyfile_used.kf`` lists the
merged composition (one ``CHAIN`` line per type) and, when the chains are no
longer grouped by type, warns that a re-run from it would number them
differently.
``EXTRA_CHAIN`` cannot be combined with ``RESTART_CONTINUE``, since adding
chains makes a different system. Because this can be repeated,
you can build a system up in stages - equilibrate component A, restart and add
component B, restart again and add component C, and so on. ``EXTRA_CHAIN`` requires
a ``RESTART_FILE`` (there must be an existing configuration to add to) and is an
error without one. Note also that the parameter file must define every bead type
used by the new chains, and that if the box is too crowded to place them the run
aborts at start-up with an "overcrowded lattice" error.
