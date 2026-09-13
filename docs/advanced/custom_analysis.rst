.. _advanced-custom-analysis:

===============
Custom analysis
===============

PIMMS ships a set of built-in analyses (radius of gyration, internal scaling,
cluster properties, and so on - see :doc:`/output_files`), but you will often want
to measure something specific to your problem. The **custom analysis** hook lets you
supply your own Python code that PIMMS loads and runs *during* the simulation, with
direct access to the live lattice, without touching or rebuilding PIMMS itself.

How it works
============

You point PIMMS at a plain Python file with two keywords:

.. code-block:: text

   ANALYSIS_MODULE : my_analysis.py    # path to your Python file
   ANA_CUSTOM      : 500               # run it every 500 steps

* ``ANALYSIS_MODULE`` is the path to your file (relative paths are resolved against
  the working directory; ``~`` is expanded). Default ``False`` - no custom analysis.
  The file is loaded and validated while the keyfile is parsed, so any problem with
  it stops the run before any simulation work is done. A successful load prints two
  lines::

     [Module Analysis]: loading custom analysis from [/abs/path/my_analysis.py]
     [Module Analysis]: 'analysis_function' loaded and validated successfully

* ``ANA_CUSTOM`` is how often, in steps, your code is called: it runs on every step
  where ``step % ANA_CUSTOM == 0``. Like every analysis in PIMMS it only runs
  **after equilibration** - ``EQUILIBRATION`` is the last equilibration step, so the
  first possible call is the first multiple of ``ANA_CUSTOM`` strictly greater than
  ``EQUILIBRATION``. Setting ``ANA_CUSTOM`` without ``ANALYSIS_MODULE`` prints a
  warning and does nothing; setting ``ANALYSIS_MODULE`` with ``ANA_CUSTOM`` at
  0/unset is a **parse-time error**, because a module that loads and validates but
  never runs is almost certainly a mistake.

Your file must define a single top-level function called **exactly**
``analysis_function`` that takes two arguments:

.. code-block:: python

   def analysis_function(step, lattice):
       ...

PIMMS calls it as ``analysis_function(step, lattice)`` every ``ANA_CUSTOM`` steps.
The return value is ignored - a custom analysis works by *doing* something (writing
a file, updating an accumulator), not by returning a value.

The arguments
=============

``step`` : int
    The current simulation step, handy for labelling output.

``lattice`` : the live ``Lattice`` object
    The actual simulation state - **not** a copy - so you can read anything about
    the current configuration. The most useful attributes are:

    .. list-table::
       :header-rows: 1
       :widths: 32 68

       * - Attribute
         - What it is
       * - ``lattice.dimensions``
         - The box dimensions, a list of length 2 or 3.
       * - ``lattice.chains``
         - Dict mapping ``chainID`` (int) to the chain object. ChainIDs start
           at 1.
       * - ``lattice.get_number_of_chains()``
         - How many chains are in the system.
       * - ``lattice.hardwall``
         - ``True`` under hardwall boundaries, ``False`` under periodic ones.
       * - ``lattice.chainIDtoType``
         - Dict mapping ``chainID`` to its integer chain *type*.
       * - ``lattice.chainTypeList``
         - List of the distinct chain types present.
       * - ``lattice.grid``
         - The occupancy grid (a NumPy array; ``0`` = empty, otherwise the
           occupying chainID).
       * - ``lattice.type_grid``
         - Companion grid holding the bead *type* at each occupied site.

    Each chain object in ``lattice.chains`` exposes, among others:

    .. list-table::
       :header-rows: 1
       :widths: 40 60

       * - Chain attribute / method
         - What it is
       * - ``chain.get_analysis_positions()``
         - The positions **every intra-chain observable should be computed
           from**: the chain's beads in sequence order (N→C), bond-walked into a
           single periodic image so the chain is contiguous even when it crosses
           a box face. Coordinates may fall outside the box.
       * - ``chain.get_ordered_positions()``
         - The raw on-lattice positions in sequence order. These are wrapped
           back into the box, so a chain that straddles a face is **not** a
           contiguous object - use ``get_analysis_positions()`` for geometry.
       * - ``chain.sequence``
         - The chain's bead sequence (a string).
       * - ``chain.seq_len``
         - Number of beads in the chain.
       * - ``chain.chainID`` / ``chain.chainType``
         - The chain's integer ID and type.

.. important::

   ``lattice`` is the real, live object. **Read** from it freely, but do **not**
   mutate it (moving beads, editing the grid, adding/removing chains) - that would
   corrupt the simulation. If you need to transform positions, copy them first.

Reusing PIMMS' own analysis
===========================

Your module can import PIMMS and reuse the same routines the built-in analyses use.
The simplest and safest route is the chain object's own ``analysis_*`` methods,
because they already apply PIMMS' conventions - notably that intra-chain geometry is
measured on the chain **made whole** rather than on minimum-image distances:

* ``chain.analysis_get_radius_of_gyration()``
* ``chain.analysis_get_polymeric_properties()`` (``[rg, asphericity]``)
* ``chain.analysis_get_end_to_end_distance()``
* ``chain.analysis_get_instantaneous_internal_scaling()`` and
  ``chain.analysis_get_instantaneous_distance_map()``
* ``chain.analysis_get_residue_residue_distance(i, j)``

Underneath, :mod:`pimms.lattice_analysis_utils` provides the raw routines -
``get_polymeric_properties(positions, dimensions, pbc_correction=...)`` (radius of
gyration and asphericity), ``get_inter_position_distance(...)``,
``get_distance_matrix(...)``, ``get_cluster_distribution(...)`` and more. If you call
these directly, pass ``chain.get_analysis_positions()`` **and**
``pbc_correction=False``: the whole-chain positions are already contiguous, and
re-applying the minimum-image correction on top of them is what the ``analysis_*``
methods exist to avoid. You are free to use NumPy, SciPy, or anything else installed
in your environment.

A worked example
================

A minimal module that records the radius of gyration of the first chain each time it
is called:

.. code-block:: python

   # my_analysis.py

   def analysis_function(step, lattice):
       # grab the first chain (chainIDs start at 1)
       first_id = sorted(lattice.chains.keys())[0]
       chain = lattice.chains[first_id]

       # radius of gyration, computed on the chain made whole
       rg = chain.analysis_get_radius_of_gyration()

       # append it to our own output file (one row per call)
       with open("custom_rg.dat", "a") as fh:
           fh.write("%d\t%.4f\n" % (step, rg))

With ``ANALYSIS_MODULE : my_analysis.py`` and ``ANA_CUSTOM : 500`` in the keyfile,
PIMMS writes a ``custom_rg.dat`` row every 500 post-equilibration steps. Written out
the long way, without the convenience method, the same quantity is

.. code-block:: python

   from pimms import lattice_analysis_utils as lau

   rg = lau.get_polymeric_properties(chain.get_analysis_positions(),
                                     lattice.dimensions,
                                     pbc_correction=False)[0]

Practical notes
===============

* **Output files.** PIMMS does not manage your custom output - you open and write
  files yourself. Open in append mode (``"a"``) if you want one growing file across
  the run, and remember the working directory is wherever PIMMS was launched.
* **Helper modules.** Your module's own directory is placed first on the import path
  while the module loads and whenever ``analysis_function`` runs, so it may
  ``import`` helpers that sit alongside it (including imports made lazily inside the
  function). PIMMS restores the original import path afterwards, so loading custom
  analysis does not change module precedence for the rest of the simulation.
* **Keep it light.** Your function runs inside the simulation loop; expensive work
  every few steps will slow the run down. Prefer a modest ``ANA_CUSTOM`` frequency
  and cache anything you can.
* **State across calls.** Because the module is imported once and the same function
  object is reused, you can keep running state in module-level variables (e.g. an
  accumulator) between calls.

Validation and error handling
=============================

The custom-analysis hook is designed to fail early and clearly:

* **Load-time validation (at keyfile parse).** As soon as the keyfile is read, PIMMS
  loads your file and checks that it exists, imports without error, defines a
  ``analysis_function``, that it is callable, and that its signature can accept the
  ``(step, lattice)`` call. Any problem aborts immediately with a clear message that
  names the file and the issue - so a typo or a missing entry point is caught in
  seconds, before a long run starts, rather than part-way through. The same pass
  rejects a loaded module whose ``ANA_CUSTOM`` is 0/unset (the module would never
  run). For example, a file that defines the function under the wrong name fails
  with::

     The custom analysis module 'my_analysis.py' does not define an
     'analysis_function'. PIMMS calls 'analysis_function(step, lattice)', so the
     file must define a top-level function with exactly that name. Callables
     defined in the file: my_other_function.

* **Isolated import.** Your file is imported by path under a private module name, so
  it cannot collide with (or be shadowed by) a PIMMS or standard-library module, even
  if you call it something like ``random.py`` or ``energy.py``.

* **Runtime errors.** If your ``analysis_function`` raises once it is running (say it
  hits a configuration it did not expect), PIMMS stops the run and reports an
  ``AnalysisRoutineException`` that names the step and makes clear the fault is in
  your analysis code rather than in PIMMS - instead of surfacing as an opaque
  traceback deep inside the engine::

     The custom analysis function (from ANALYSIS_MODULE) raised ValueError at
     step 100: something unexpected. This is an error in your custom analysis
     code, not in PIMMS.
