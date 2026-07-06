.. _installation:

============
Installation
============

PIMMS is distributed from GitHub and installed with ``pip`` (or ``uv``). Because
the performance-critical parts of PIMMS are written in `Cython
<https://cython.org/>`_, installing PIMMS **compiles native C extensions** on your
machine - this happens automatically, but it means you need a working C compiler
and the build dependencies described below.

.. note::

   PIMMS was recently added to PyPI (1.0.x release); let us know if you have any issues.

Requirements
============

* **Python ≥ 3.8** (3.10+ recommended; the development/test environment is 3.12).
* A **C compiler** (clang on macOS, gcc on Linux).
* The build/runtime Python packages: ``numpy``, ``scipy``, ``cython``,
  ``versioningit`` and ``mdtraj`` (the last provides the XTC trajectory backend).

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

   pip install numpy scipy cython versioningit
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
   declares its build dependencies (Cython, NumPy, versioningit), so pip's default
   isolated build already has what it needs. Passing ``--no-build-isolation`` simply
   tells pip to build against the packages you installed in Step 1 rather than
   fetching them again into a throwaway build environment - slightly faster, and
   the reason Step 1 installs the build tools up front.

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

   ./build.sh        # clean rebuild + editable reinstall (development only)

The compiled modules include the serial and parallel move kernels
(``mega_crank``, ``mega_crank_2D`` and ``mega_crank_fast`` - the last holding the
multi-threaded crankshaft, slither and pull kernels), the energy inner loops
(``inner_loops``, ``inner_loops_hardwall``) and several utilities. On macOS the
multi-threaded kernels use OpenMP via Homebrew ``libomp`` if present, and degrade
gracefully (single-threaded) if not.

Verifying the installation
==========================

Open a **new terminal**, activate the environment, and check the CLI:

.. code-block:: bash

   PIMMS --version          # prints the installed version
   PIMMS --info             # lists every keyfile keyword, grouped
   PIMMS --info DIMENSIONS  # full details on one keyword

.. note::

   The *very first* ``PIMMS`` invocation may take 5-10 seconds while Python
   initialises; subsequent calls are fast.

You can also confirm the package imports from Python:

.. code-block:: bash

   python -c "import pimms; print(pimms.__version__)"

Finally, run a bundled demo. Each directory under ``demo_keyfiles/`` contains a
``KEYFILE.kf`` (simulation configuration) and a ``params.prm`` (force field):

.. code-block:: bash

   cd demo_keyfiles/single_chain_polymer
   PIMMS -k KEYFILE.kf

This writes ``ENERGY.dat``, a ``START.pdb``/``traj.xtc`` trajectory, and the
requested analysis files into the working directory (see :doc:`output_files`).

Running the tests
=================

PIMMS ships with an extensive test suite. From a source checkout:

.. code-block:: bash

   pytest -m "not slow" pimms/tests/      # fast suite
   pytest pimms/tests/                    # full suite (incl. slow detailed-balance tests)
   pytest --cov=pimms pimms/tests/        # with coverage

The ``slow`` marker (defined in ``pyproject.toml``) tags the heavier
detailed-balance tests; deselect them with ``-m "not slow"`` for quick iteration.
A clean run of the full suite is the strongest confirmation that the Cython
extensions built correctly and PIMMS is behaving as expected.
