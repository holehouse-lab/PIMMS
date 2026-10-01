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
has its own page: :doc:`/lemonade/reference`. (The keyword tables in
``pimms/CONFIG.py`` - names, types, defaults and descriptions - are rendered as the
:doc:`keyword reference </keywords>` rather than here.)

The developer scripts that ship inside the package but are never imported by it
(``check_randomness``, ``cython_testing``, ``dev_megacrank_rng_check``,
``dev_randint_check``, ``dev_randneg_check`` and ``print_interaction_matrix``) are
not listed either; each does its work only under ``if __name__ == "__main__"``.

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

.. automodule:: pimms.crankshaft_list_functions
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

.. automodule:: pimms.numpy_utils
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

.. automodule:: pimms.file_utilities
   :members:
   :show-inheritance:

.. automodule:: pimms.CONFIG
   :members:
   :show-inheritance:

Analysis, output & logging
==========================

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

.. automodule:: pimms.pdb_utils
   :members:
   :show-inheritance:

.. automodule:: pimms.IO_utils
   :members:
   :show-inheritance:

.. automodule:: pimms.pimmslogger
   :members:
   :show-inheritance:

Exceptions
==========

.. automodule:: pimms.latticeExceptions
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
the detailed-balance tests in ``pimms/tests/test_detailed_balance.py`` (all marked
``slow``) - slither, pull, the parallel crankshaft, slither and pull kernels, the
single-chain translate, rotate, pivot and head-pivot moves, TSMMC, VMMC,
jump-and-relax and the
cluster moves - and an end-to-end regression suite in
``pimms/tests/simulation_tests/`` that runs ``scripts/PIMMS`` over 15 scenarios and
diffs the output against stored expected output. Each scenario's input files are
copied into a temporary directory made by pytest and the run happens there, so
the suite writes nothing into the source tree and two runs on one checkout do not
collide; a failure message gives the path of that scenario's log. From a checkout
the harness forces the repo onto the subprocess ``PYTHONPATH``, so it always tests
the working tree rather than an installed copy. The tests also ship in the wheel,
and run from an installed copy (``pytest --pyargs pimms``) the same harness runs
the installed ``PIMMS`` executable against the installed package instead - it
decides which case it is in from its own location only, never from a checkout
that happens to enclose the environment. A few tests need files that only a
checkout has: from an installed copy the packaging tests (which read
``MANIFEST.in`` and ``pyproject.toml``) are skipped, and the command-line tests
that look for ``scripts/PIMMS`` fail.
``pimms/tests/simulation_tests/readme.md`` explains how to add a
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
