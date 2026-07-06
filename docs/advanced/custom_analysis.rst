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
* ``ANA_CUSTOM`` is how often, in steps, your code is called. Like every analysis in
  PIMMS it only runs **after equilibration**, and a value below ``1`` disables it.
  If ``ANALYSIS_MODULE`` is not set, ``ANA_CUSTOM`` is ignored.

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
         - Dict mapping ``chainID`` (int) to the chain object.
       * - ``lattice.chainIDtoType``
         - Dict mapping ``chainID`` to its integer chain *type*.
       * - ``lattice.chainTypeList``
         - List of the distinct chain types present.
       * - ``lattice.grid``
         - The occupancy grid (a NumPy array; ``0`` = empty).
       * - ``lattice.type_grid``
         - Companion grid holding the bead *type* at each occupied site.

    Each chain object in ``lattice.chains`` exposes, among others:

    .. list-table::
       :header-rows: 1
       :widths: 40 60

       * - Chain attribute / method
         - What it is
       * - ``chain.get_ordered_positions()``
         - The chain's bead positions in sequence order (N→C), as a list of
           ``[x, y, z]`` (or ``[x, y]`` in 2D).
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
For example, :mod:`pimms.lattice_analysis_utils` provides
``get_polymeric_properties(positions, dimensions)`` (radius of gyration and
asphericity), ``get_inter_position_distance(...)``, cluster helpers, and more; and
each chain object also carries convenience methods such as
``analysis_get_polymeric_properties()`` and
``analysis_get_end_to_end_distance()``. You are free to use NumPy, SciPy, or anything
else installed in your environment.

A worked example
================

A minimal module that records the radius of gyration of the first chain each time it
is called:

.. code-block:: python

   # my_analysis.py
   from pimms import lattice_analysis_utils as lau

   def analysis_function(step, lattice):
       # grab the first chain
       first_id = sorted(lattice.chains.keys())[0]
       chain = lattice.chains[first_id]

       # radius of gyration from the chain's current positions
       rg = lau.get_polymeric_properties(chain.get_ordered_positions(),
                                         lattice.dimensions)[0]

       # append it to our own output file (one row per call)
       with open("custom_rg.dat", "a") as fh:
           fh.write("%d\t%.4f\n" % (step, rg))

With ``ANALYSIS_MODULE : my_analysis.py`` and ``ANA_CUSTOM : 500`` in the keyfile,
PIMMS writes a ``custom_rg.dat`` row every 500 post-equilibration steps.

Practical notes
===============

* **Output files.** PIMMS does not manage your custom output - you open and write
  files yourself. Open in append mode (``"a"``) if you want one growing file across
  the run, and remember the working directory is wherever PIMMS was launched.
* **Helper modules.** Your module's own directory is added to the import path, so it
  may ``import`` helper modules that sit alongside it.
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
  seconds, before a long run starts, rather than part-way through. For example, a
  file that defines the function under the wrong name fails with::

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
  traceback deep inside the engine.
