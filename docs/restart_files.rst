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
a temporary file (``restart.pimms.tmp.<pid>``, named after the process so that two
runs can never share it) which is flushed and synced to disk and then renamed
into place, so neither a crash part-way through a checkpoint nor a power failure
just after one can destroy the previous good ``restart.pimms``. A temporary left
behind by a run that was killed inside the write is removed by a later run that
starts in the directory, once it is certain that no write is still in progress:
at the next start if its process is known to be gone from that machine,
otherwise when the file is more than an hour old (until then ``log.txt`` says it
was left in place).

The checkpoint is the **last** thing written for its step: the step's energy row,
its trajectory frame (flushed to disk first) and every analysis due on that step
are written before it. A checkpoint therefore never describes a step whose
output is missing, which is what lets a continuation start cleanly after it
(:ref:`restart-stopping`). The one exception is the trajectory under
``SAVE_AT_END : True``: that mode keeps every frame in memory until the run
ends, so there is nothing to flush, and a run killed outright leaves a
``traj.xtc`` holding frame 0 alone whatever step its checkpoint is from.

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
     'RNG_NUMPY'    : <tuple>,            # numpy.random.get_state() at that moment, with
                                          # its 624 key words as a list of Python ints
     'EQUILIBRIUM_TEMPERATURE' : <float>, # the temperature the Hamiltonian was built at
     'QUENCH'       : {                   # the run's quench settings
         'QUENCH_RUN': True | False,      # (the four below only when it is True)
         'QUENCH_START': <float>, 'QUENCH_END': <float>,
         'QUENCH_STEPSIZE': <float>, 'QUENCH_FREQ': <int>,
     },
     'PIMMS_VERSION': '<version>',
   }

where, for each chain, ``positions`` is the list of bead coordinates (in N→C
order), ``sequence`` is the one-letter sequence string, and ``chainType`` is an
integer grouping identical chains. A new run always recomputes its energy from
scratch; a plain restart never reads ``ENERGY``, and ``RESTART_CONTINUE`` reads it
only to check that the recomputed value agrees with it (:ref:`restart-continue`).
The file holds plain Python objects only (no numpy arrays or numpy scalars), so
it can be read under any numpy version, whichever one wrote it.

The first four keys are present in every restart file, and the configuration
they hold (``ENERGY`` aside) is all that a plain restart (``RESTART_FILE``
without ``RESTART_CONTINUE``) uses: the new run begins fresh in every other
respect, with its own step count, seed, temperature and statistics. The rest
are the *continuation state*, written by every checkpoint since PIMMS 1.0.8
and acted on only by ``RESTART_CONTINUE``, which restores the step, the two
generators and the temperature so that the resumed run reproduces the
uninterrupted one, and checks ``EQUILIBRIUM_TEMPERATURE`` and ``QUENCH`` against
the keyfile so that it is the same run (:ref:`restart-continue`);
``PIMMS_VERSION`` is kept for the record, and the step is also quoted in
``keyfile_used.kf``. Restart files from
earlier versions lack them, still load as configurations, and are refused for
exact continuation with a message saying why. Move statistics and the
accumulated end-of-run analyses are not stored in either case: they are running
totals of the process that produced them.

The reader validates the whole file before it changes anything, so a bad restart
fails immediately with a message saying what is wrong (naming the offending
chain for most chain-level problems) rather than corrupting
a run: the top-level object must be a dictionary with the four configuration
keys (the continuation keys are optional, but if present must be well-formed:
a non-negative integer step, a positive finite temperature, and both generator
states together, each of which is loaded into a throwaway generator so that a
state the generator would refuse is caught here and not at the start of the
master loop), dimensions
must be 2 or 3 positive integers, ``HARDWALL`` a boolean and ``ENERGY`` finite,
``CHAINS`` must hold at least one chain,
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
   segments must be concatenated with that in mind. With ``RESTART_CONTINUE``
   PIMMS refuses to start in a directory that still holds simulation output,
   because the rows and frames before the checkpoint would be deleted and never
   written again (:ref:`restart-continue`).

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
guarantee without your having to copy the values across.

The temperature schedule must be the one the stopped run had. The resumed run
takes its temperature from the checkpoint, but the Hamiltonian's
``ANGLE_PENALTY_T_NORM`` scaling is built from the keyfile (at ``QUENCH_END`` in
a ``QUENCH_RUN``, at ``TEMPERATURE`` otherwise) and so is the rest of a quench
ramp, so a keyfile that differs would give a run that is neither the original
nor the one asked for. The checkpoint records whether the run was a quench and,
if so, its ``QUENCH_START``, ``QUENCH_END``, ``QUENCH_STEPSIZE`` and
``QUENCH_FREQ``, and all of them must match; a run that was not a quench must be
continued without ``QUENCH_RUN`` and at the checkpoint's ``TEMPERATURE``. A
quench checkpoint cannot be continued as a plain run at the temperature it had
reached (the refusal lists the quench settings to use). A checkpoint written
before these settings were stored (by an earlier 1.0.8 development version) is
held to the older, weaker rule and the message says so: ``TEMPERATURE`` must
equal the checkpoint's temperature, or in a quench the checkpoint's temperature
must lie between ``QUENCH_START`` and ``QUENCH_END``.

Two further checks protect the run itself:

* **The energy function must be the same.** Once the lattice has been rebuilt,
  PIMMS computes the energy of the checkpoint's configuration and compares it
  with the ``ENERGY`` the checkpoint recorded. They are equal for a genuine
  continuation, in every case (a quench in progress, a run whose
  ``RESIZED_EQUILIBRATION`` phase is over, ``FREEZE_FILE``, hard walls,
  ``PARALLELIZE``). A difference means a different ``PARAMETER_FILE``, a
  different ``ANGLES_OFF`` or ``NON_INTERACTING``, or different quench settings,
  and the run is refused, naming those causes. Settings that leave the energy
  function alone are *not* checked: a different move set, different
  ``*_SUBSTEPS`` or a different ``PARALLELIZE`` setting gives a valid run that is
  simply not the one that was stopped.
* **The working directory must not hold simulation output** (any analysis
  ``.dat`` file, ``START.pdb`` or ``traj.xtc``). Every run start deletes the
  previous run's output, and a continuation never writes the steps before its
  checkpoint again, so resuming in the stopped segment's own directory used to
  delete that segment silently. The refusal lists what it found and leaves the
  directory exactly as it was. ``restart.pimms``, ``log.txt``,
  ``parameters_used.prm`` and ``keyfile_used.kf`` do not count.

A restart file with no continuation state (one written by an earlier PIMMS, or
by a ``RestartObject`` built outside a run) is refused with a message saying so,
as is one that has a step and generator states but no ``TEMPERATURE``; either
can still seed a new run in the ordinary way. Every problem the keyfile checks
find is listed in the one refusal, each with what to change. (Some combinations
never reach these checks: a periodic checkpoint with ``HARDWALL : True`` or with
``RESIZED_EQUILIBRATION`` is refused first by the general restart rules above.)

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

.. note::

   **Custom analysis and random numbers.** The generator states in a checkpoint
   are captured after every analysis due on that step has run. The built-in
   analyses draw no random numbers (this is tested), so for them the moment of
   capture makes no difference. A custom ``ANALYSIS_MODULE`` that draws from
   Python's ``random`` or from ``numpy.random`` *does* advance the stream the
   moves use: such a run is still resumed exactly, because its draws on the
   checkpoint step are made before the state is saved and are not repeated, but
   it is a different run from the same keyfile without the module. Give a
   custom analysis its own generator (``random.Random(seed)``,
   ``numpy.random.default_rng(seed)``) if it needs random numbers; note that the
   state of such a private generator is not stored in the restart file.

.. _restart-stopping:

Stopping a run, and what a stopped run leaves behind
====================================================

**A stop you asked for.** ``SIGTERM`` (what ``kill`` and a batch scheduler's
time limit send) and ``SIGINT`` (Ctrl-C) both end a run cleanly: the trajectory
is closed with every completed frame in it (under ``SAVE_AT_END`` the frames
buffered so far are written out), and one line is printed and written to
``log.txt`` saying after which step the run stopped and which step the
checkpoint on disk is from, for example::

   Run stopped by SIGTERM after step 41873 of 200000. The checkpoint on disk
   (restart.pimms) is from step 40000; RESTART_CONTINUE resumes the run from
   there, in a fresh directory.

"After step N" means that everything due on step N is on disk. If the signal
landed while that step's analyses were being written the line says "during the
output of step N" instead: the step's move, energy row and trajectory frame are
complete but its analysis rows may not be. A further ``SIGTERM`` or ``SIGINT``
sent while the run is cleaning up is ignored, so an impatient second ``kill`` or
Ctrl-C cannot cut the clean-up short.

No checkpoint is written at the moment of the stop: a signal can arrive in the
middle of a move, so the last periodic checkpoint is the state to resume from
(choose ``RESTART_FREQ`` with that in mind). The exit status of the ``PIMMS``
command is the conventional one, 143 after ``SIGTERM`` and 130 after ``SIGINT``
(a program that calls ``Simulation.run_simulation()`` itself sees ``SystemExit``
and ``KeyboardInterrupt`` respectively, after the clean-up). A signal that
arrives before the run has started, while the keyfile is being parsed or the
simulation built, still ends the process at once.
The ``SIGTERM`` handling is installed only while the run is going and only if
nothing else has installed a handler of its own, so a program that embeds
``Simulation`` keeps whatever behaviour it set up; off the main thread no
handler is installed at all.

**A stop you did not ask for** (``SIGKILL``, an out-of-memory kill, a power
failure) cannot be tidied up, but the files are still consistent with each
other in the following sense. ``restart.pimms`` is always a complete
checkpoint, and every row and frame up to and including its step is on disk -
except the frames of a ``SAVE_AT_END : True`` run, which are only in memory
until the run ends and are all lost with the process (``traj.xtc`` then holds
frame 0 alone; a ``SIGTERM`` or ``SIGINT`` stop does write them out).
The other files may run *past* it: the killed segment carried on for some steps
after its last checkpoint, so its ``.dat`` files and trajectory hold rows and
frames for steps the continuation, which starts after the checkpoint's step,
writes again. Concatenating the two segments as they are therefore duplicates
those steps. Before concatenating, cut the stopped segment back to the
checkpoint's step, which the continuation's ``log.txt`` and ``keyfile_used.kf``
both quote (and ``pickle.load(open('restart.pimms', 'rb'))['STEP']`` gives):

.. code-block:: python

   import numpy as np
   step = 40000                                    # the checkpoint's step
   rows = np.loadtxt('stopped/ENERGY.dat', ndmin=2)
   rows = rows[rows[:, 0] <= step]                 # every per-step file has the step first

Trajectory frames carry no step label, so count them: a segment holds frame 0
(its starting configuration) followed by one frame per ``XTC_FREQ`` multiple it
reached, so with the default ``SAVE_EQ : True`` a segment that started at step 0
should keep its first ``1 + step // XTC_FREQ`` frames (frame 0 and the multiples
at or before the checkpoint's step). Frame 0 of the continuation is the
checkpoint configuration itself, and is **always** dropped when joining: if
``step`` is a multiple of ``XTC_FREQ`` the stopped segment already has that
frame, and if it is not, the frame falls between two ``XTC_FREQ`` multiples and
the uninterrupted run would not have written it at all.

.. code-block:: python

   import mdtraj as md
   xtc_freq = 500                                  # the run's XTC_FREQ
   stopped = md.load('stopped/traj.xtc', top='stopped/START.pdb')
   resumed = md.load('resumed/traj.xtc', top='resumed/START.pdb')
   joined = stopped[:1 + step // xtc_freq].join(resumed[1:])

For example, with ``XTC_FREQ : 5``, a run of 400 steps stopped and resumed from
a checkpoint at step 91: the stopped segment keeps ``1 + 91 // 5 = 19`` frames
(steps 0, 5, ..., 90), the continuation holds 63 (the checkpoint, then steps 95,
100, ..., 400) of which 62 are kept, and 19 + 62 is the 81 frames of the
uninterrupted run.

**No checkpoint is written during equilibration.** The first one is written at
the first ``RESTART_FREQ`` multiple after ``EQUILIBRATION``, and until then any
``restart.pimms`` already in the directory - an earlier run's, or the file this
run started from - stays where it is (it is deliberately not removed at
start-up, so that a run which dies early does not also destroy the only
checkpoint there is). A run that stops before its first checkpoint therefore
leaves a ``restart.pimms`` that does **not** describe it, and
``RESTART_CONTINUE`` from that file resumes the earlier run. ``log.txt`` records
at start-up when such a file is present, and the line written when a run is
interrupted says whether the run had written a checkpoint of its own. Check the
``STEP`` in the file against the log before continuing from it.

**A start that fails.** ``log.txt`` is initialised when the keyfile is parsed
and ``parameters_used.prm`` when the Hamiltonian is built, so both can describe
a run that then failed to start, next to an earlier run's data. In that case
``log.txt`` ends with a line beginning ``START FAILED - NO SIMULATION WAS RUN``.
``keyfile_used.kf`` is written only once the run has really started (the
simulation was built and its trajectory opened), so a ``keyfile_used.kf`` in a
directory always belongs to the run whose trajectory is there, and a failed
start removes none of an earlier run's data.

**Two runs in one directory** overwrite each other's files. While a run is
going it keeps a line (its process number, host and start time) in a small file,
``pimms_running.pid``, in the working directory; it takes the line out when it
ends, and the file goes with its last line. A run that finds the line of another
run that may still be going prints and logs a warning and adds its own line
beside it, so a third run is warned about a first that is still going after a
second has come and gone. It does not refuse, because a process number can be
reused and a line left by a killed run on another machine cannot be told from a
live one. The line of a run that was killed on the same machine is dropped
without comment at the next start; one from another machine cannot be checked,
so it stays, with a warning at every start, until ``pimms_running.pid`` is
deleted by hand.

**Input files that are also output files.** PIMMS writes its output under fixed
names in the working directory, so an input file with one of those names (or a
link to one) would be overwritten or deleted by the run that reads it. Such a
run is refused at start-up, naming the keyword and the output file. Three
pairings are allowed because the file is read in full before it is rewritten
with equivalent content: running ``keyfile_used.kf`` as the keyfile, using
``parameters_used.prm`` as the ``PARAMETER_FILE`` (it gains another four header
lines each time), and a ``RESTART_FILE`` named ``restart.pimms`` (replaced at
the first checkpoint). The check is made while the keyfile is being parsed,
before ``log.txt`` is started, so the refusal leaves every file as it was.

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
