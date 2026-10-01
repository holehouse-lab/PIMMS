.. _installation:

============
Installation
============

PIMMS is distributed on PyPI (as ``idptools-pimms``) and from GitHub, and is
installed with ``pip`` (or ``uv``). Because the performance-critical parts of PIMMS
are written in `Cython <https://cython.org/>`_, installing PIMMS **compiles native
C extensions** on your machine - this happens automatically, but it means you need
a working C compiler and the build dependencies described below.

Requirements
============

* **Python ≥ 3.10.** The full test suite - including the compiled Cython kernels - is
  verified on **3.10, 3.11, 3.12, 3.13 and 3.14**; the development environment is 3.12.
  There is no upper bound: PIMMS uses no version-gated syntax or standard library API,
  so a new Python release is expected to work as soon as its ``numpy``, ``scipy`` and
  ``mdtraj`` wheels are published - that, rather than PIMMS itself, is the practical
  gate on using a just-released interpreter.
* A **C compiler** (clang on macOS, gcc on Linux).
* The **runtime** dependencies, installed automatically with PIMMS: ``numpy``
  (≥ 1.21), ``scipy`` (≥ 1.9) and ``mdtraj`` (≥ 1.10; provides the XTC trajectory
  backend). These minimums are tested, not nominal - the suite is run against
  exactly these versions as well as against current releases.
  Note ``scipy`` ≥ 1.9 specifically: 1.7 does not expose ``scipy.spatial.QhullError``
  (which the analysis code imports, so no simulation can start), and the 1.8 macOS
  arm64 wheels crash inside their own LAPACK during the binodal ``curve_fit``.
* The **build-time** dependencies, which ``pip`` fetches into its isolated build
  environment unless you pass ``--no-build-isolation`` (see Step 1):
  ``setuptools`` (≥ 77), ``wheel``, ``cython``, ``numpy`` and ``versioningit`` (≥ 2;
  it derives the version from the git tags, which is why the package has no
  hardcoded version string).

We strongly recommend installing into a clean, dedicated environment.

.. code-block:: bash

   # with conda
   conda create -n pimms python=3.12 -y
   conda activate pimms

   # ...or with uv
   uv venv --python 3.12
   source .venv/bin/activate

Step 1 - install the dependencies
=================================

.. code-block:: bash

   pip install numpy scipy cython versioningit "setuptools>=77" wheel
   pip install mdtraj

(With ``uv``, use ``uv pip install ...`` instead of ``pip install ...``.)

Step 2 - install PIMMS
======================

Install from PyPI:

.. code-block:: bash

   pip install idptools-pimms


Install directly from GitHub:

.. code-block:: bash

   pip install --no-build-isolation git+https://github.com/holehouse-lab/PIMMS.git

.. note::

   The ``--no-build-isolation`` flag is optional. PIMMS' ``pyproject.toml``
   declares its build dependencies (setuptools, wheel, Cython, NumPy,
   versioningit), so pip's default isolated build already has what it needs.
   Passing ``--no-build-isolation`` simply tells pip to build against the packages
   you installed in Step 1 rather than fetching them again into a throwaway build
   environment - slightly faster, and the reason Step 1 installs the build tools up
   front.

Or clone and install from source (recommended if you intend to develop PIMMS):

.. code-block:: bash

   git clone https://github.com/holehouse-lab/PIMMS.git
   cd PIMMS
   pip install -e . --upgrade --force-reinstall          # editable install
   # ...or, with uv:
   uv pip install -e . --no-deps --reinstall

Do I need to run ``build.sh``?
==============================

**No - not for a normal install.** PIMMS' ``setup.py`` declares all of its Cython
modules as ``ext_modules`` via ``cythonize(...)``, so ``pip install`` (whether
from the GitHub URL or from source) compiles every extension automatically. There
is nothing extra to run.

``build.sh`` is a **developer convenience for rebuilding after you edit a**
``.pyx`` **file**. Cython skips regenerating a ``.c`` file that is newer than its
``.pyx`` and ``build_ext`` reuses cached object files, so a plain reinstall may
not pick up ``.pyx`` edits. ``build.sh`` forces a clean rebuild by deleting the
generated C (``pimms/*.c``), the compiled extensions (``pimms/*.so``) and the
``build/`` cache, then reinstalling:

.. code-block:: bash

   ./build.sh uv     # clean rebuild + editable reinstall with uv  (development only)
   ./build.sh pip    # the same, reinstalling with pip

It is a ``zsh`` script that takes one argument naming the installer: ``uv``
reinstalls with ``uv pip install -e . --no-deps --reinstall``, ``pip`` with
``python -m pip install -e . --force-reinstall --no-deps``. With no argument, or
an unknown one, it prints its usage and touches nothing, and it checks that the
chosen tool works before it deletes anything (a ``uv``-managed environment ships
without ``pip``). It sweeps the generated C and the
compiled extensions in ``pimms/`` and in ``pimms/lemonade/kernels/``, so every
kernel, including the ``lemonade`` kernel (``pimms/lemonade/kernels/_pbc.pyx``),
is recompiled from its ``.pyx``.

The compiled modules are the move kernels (``mega_crank_fast``, which holds every
production crankshaft, slither and pull kernel - serial and multi-threaded, 2D and
3D - and ``mega_crank`` / ``mega_crank_2D``, the original crankshaft kernels kept
as the bit-exactness reference for the fast ones), the energy inner loops
(``inner_loops``, ``inner_loops_hardwall``), the ``hyperloop``, ``system_utils``,
``cluster_kernels`` and ``bookkeeping`` utilities, and the ``lemonade``
periodic-boundary kernel (``pimms.lemonade.kernels._pbc``). What each kernel does
and how data flows between the Python objects and the compiled code is described
in ``pimms/kernels.md`` in the source tree. Only ``mega_crank_fast`` is built with
OpenMP: on Linux via the compiler's ``-fopenmp``, which is always passed there (so
a compiler without OpenMP support fails the build rather than falling back), and
on macOS via Homebrew ``libomp`` if present (looked for in
``/opt/homebrew/opt/libomp`` and ``/usr/local/opt/libomp``). A macOS build
without ``libomp`` still succeeds, and the parallel kernels then run their blocks
one after another - identical sampling, no speed-up. To check which you have:

.. code-block:: bash

   python -c "from pimms import mega_crank_fast; print(mega_crank_fast.openmp_info())"

which prints, for example, ``{'enabled': True, 'max_threads': 16}`` (``enabled``
is False for a build without OpenMP).

Verifying the installation
==========================

Open a **new terminal**, activate the environment, and check the CLI:

.. code-block:: bash

   PIMMS --help             # the four command-line flags
   PIMMS --version          # prints the installed version
   PIMMS --info             # lists every keyfile keyword, grouped by purpose
   PIMMS --info DIMENSIONS  # type, description and default for one keyword
   PIMMS --info ALL         # the full description of every keyword

That is the whole command-line interface: ``-k``/``-keyfile`` runs a simulation
from a keyfile (note the single dash - ``--keyfile`` is not accepted),
``-i``/``--info`` prints keyword documentation, ``-v``/``--version`` prints the
version, and ``-h``/``--help`` prints the usage. Keyword lookups are
case-insensitive (``--info dimensions`` works), and an unrecognised name prints a
short pointer back to ``--info``. Running ``PIMMS`` with no arguments prints a
reminder to try ``--help``. For scripting: a completed run, ``--help``,
``--version``, a successful ``--info`` lookup and a bare ``PIMMS`` exit with
status 0; a keyfile that cannot be opened (including ``-k ""``), a keyfile that
fails validation, any error during the run (printed with its traceback), or an
unrecognised ``--info`` keyword exits with status 1; and an unrecognised flag
exits with status 2 (argparse).

.. note::

   ``--help``, ``--version`` and ``--info`` deliberately avoid importing the
   simulation stack (SciPy, MDTraj and the compiled kernels), so they answer
   immediately and still work if an extension failed to build. Starting an actual
   simulation pays that import cost once, at startup.

You can also confirm the package imports from Python:

.. code-block:: bash

   python -c "import pimms; print(pimms.__version__)"
   python -c "import pimms.simulation"   # also loads the compiled kernels

``import pimms`` on its own is deliberately light and does not touch the compiled
extensions; importing ``pimms.simulation`` loads the simulation kernels and fails
with an ``ImportError`` if one of them did not build.

Finally, run one of the demos in the repository. These live under
``demo_keyfiles/`` in a git clone (they are not shipped inside the installed
package), and each directory contains a keyfile (usually ``KEYFILE.kf``; the
simulation configuration) and a parameter file (usually ``params.prm``; the force
field); most also carry a ``readme.md`` describing what the demo does:

.. code-block:: bash

   cd demo_keyfiles/single_chain_polymer
   PIMMS -k KEYFILE.kf

This writes ``ENERGY.dat``, a ``START.pdb``/``traj.xtc`` trajectory, and the
requested analysis files into the working directory (see :doc:`output_files`).

Running the tests
=================

PIMMS ships with an extensive test suite. It needs ``pytest``, which the ``test``
extra installs (``pip install -e ".[test]"``; there is a ``docs`` extra for the
Sphinx build in the same way). From a source checkout:

.. code-block:: bash

   pytest                                 # full suite (incl. slow detailed-balance tests)
   pytest -m "not slow"                   # fast suite
   pytest pimms/tests/test_moves.py       # one file
   pytest --cov=pimms                     # with coverage (needs pytest-cov)

A bare ``pytest`` picks up both test packages (``pimms/tests`` and the ``lemonade``
tests) via the ``testpaths`` setting in ``pyproject.toml``; pass a path explicitly to
run a subset. The ``slow`` marker (also defined in ``pyproject.toml``) tags the heavier
detailed-balance tests; deselect them with ``-m "not slow"`` for quick iteration.
The suite includes end-to-end regression runs that drive the ``PIMMS`` executable
over the scenarios in ``pimms/tests/simulation_tests/`` and compare the output files
against stored expected output, so a clean run of the full suite is the strongest
confirmation that the Cython extensions built correctly and PIMMS is behaving as
expected.
