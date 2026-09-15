.. _development:

===========================
Development / API reference
===========================

This section documents PIMMS' internal Python API, generated automatically from
the NumPy-style docstrings in the source (Sphinx ``autodoc`` plus ``napoleon``,
configured in ``docs/conf.py``). It is aimed at developers extending PIMMS or
scripting against it; **users configuring simulations want the**
:doc:`keyword reference </keywords>` **instead**.

The package is organised around a few core objects: a :class:`~pimms.lattice.Lattice`
holds the grid and the :class:`~pimms.chain.Chain` objects; a
:class:`~pimms.energy.Hamiltonian` evaluates the energy; a
:class:`~pimms.moves.MoveObject` implements the Monte Carlo moves; an
:class:`~pimms.acceptance.AcceptanceCalculator` handles move selection and the
Metropolis criterion; and the :class:`~pimms.simulation.Simulation` ties them
together and drives the run, configured by the :class:`~pimms.keyfile_parser.KeyFileParser`.

Two things are documented elsewhere. The compiled Cython kernels (``mega_crank``,
``mega_crank_fast``, ``mega_crank_2D``, ``inner_loops``, ``inner_loops_hardwall``,
``hyperloop``, ``system_utils``, ``cluster_kernels``, and ``bookkeeping``, the
per-megamove copies between the Chain objects and the kernels' bead table) are
mocked during the docs build and so are not listed below. ``pimms/kernels.md`` in
the source tree describes what each of them does, how the Python objects are
turned into the arrays the kernels work on and back, and where the random numbers
come from; the ``.pyx`` sources carry the per-function detail. The analysis package
has its own page: :doc:`/lemonade/reference`.

Simulation engine
=================

.. automodule:: pimms.simulation
   :members:
   :show-inheritance:

.. automodule:: pimms.nonequilibrium_utils
   :members:
   :show-inheritance:

Monte Carlo moves
=================

.. automodule:: pimms.moves
   :members:
   :show-inheritance:

.. automodule:: pimms.acceptance
   :members:
   :show-inheritance:

.. automodule:: pimms.moveEvent
   :members:
   :show-inheritance:

.. automodule:: pimms.chainTSMMC
   :members:
   :show-inheritance:

Energy
======

.. automodule:: pimms.energy
   :members:
   :show-inheritance:

Lattice, chains and geometry
============================

.. automodule:: pimms.lattice
   :members:
   :show-inheritance:

.. automodule:: pimms.chain
   :members:
   :show-inheritance:

.. automodule:: pimms.lattice_utils
   :members:
   :show-inheritance:

.. automodule:: pimms.longrange_utils
   :members:
   :show-inheritance:

Input parsing & configuration
=============================

.. automodule:: pimms.keyfile_parser
   :members:
   :show-inheritance:

.. automodule:: pimms.parameterfile_parser
   :members:
   :show-inheritance:

.. automodule:: pimms.restart
   :members:
   :show-inheritance:

.. automodule:: pimms.data_structures
   :members:
   :show-inheritance:

Analysis & output
=================

.. automodule:: pimms.analysis_IO
   :members:
   :show-inheritance:

.. automodule:: pimms.analysis_general
   :members:
   :show-inheritance:

.. automodule:: pimms.analysis_structures
   :members:
   :show-inheritance:

.. automodule:: pimms.lattice_analysis_utils
   :members:
   :show-inheritance:

.. automodule:: pimms.cluster_utils
   :members:
   :show-inheritance:

Working on PIMMS
================

**Rebuilding after a kernel edit.** ``pip``/``uv`` compiles every extension on
install, so a normal install needs nothing extra. After editing a ``.pyx`` file
run ``./build.sh uv`` (or ``./build.sh pip``) from the repo root, which deletes the generated C and the
compiled extensions (in ``pimms/`` and in ``pimms/lemonade/kernels/``) along with
the ``build/`` cache before reinstalling - without that, Cython and ``build_ext``
can silently reuse stale artefacts. See :doc:`/installation` for the details.

**Tests.** ``pytest`` runs everything under ``pimms/tests`` and
``pimms/lemonade/tests``. Alongside the unit tests there are kernel-correctness
tests (serial vs optimised vs parallel kernels, 2D and 3D, hardwall and periodic),
detailed-balance tests for the megamove and collective kernels - slither, pull,
TSMMC, VMMC, jump-and-relax, the cluster moves and their parallel variants, all
marked ``slow`` - and an end-to-end regression suite in
``pimms/tests/simulation_tests/`` that runs ``scripts/PIMMS`` over 15 scenarios and
diffs the output against stored expected output. That harness forces the repo onto
the subprocess ``PYTHONPATH``, so it always tests the working tree rather than an
installed copy. ``pimms/tests/simulation_tests/readme.md`` explains how to add a
scenario and regenerate its expected output.

**Benchmarks and kernel validation.** ``pimms/fast_kernels/`` holds the harnesses
used to develop the optimised kernels, each run from the repo root:

.. code-block:: bash

   python pimms/fast_kernels/benchmark.py           # reference vs optimised crankshaft:
                                                    #   bit-for-bit equality + speedup
   python pimms/fast_kernels/benchmark_parallel.py  # parallel checkerboard kernel: energy
                                                    #   consistency, bead conservation, scaling
   python pimms/fast_kernels/benchmark_parallel_2d.py  # the 2D speed-up tables on the
                                                    #   parallelization page (wall vs kernel time)
   python pimms/fast_kernels/end_to_end.py          # full simulation, both kernels, compared

``python -m pimms.tests.benchmarks --scenarios all --repeats 3`` times full CLI runs
over the regression scenarios and writes a status log (see
``pimms/tests/benchmarks/README.md`` for the options).
