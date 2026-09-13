# Compiling pimms's Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/). To compile the docs, first install the documentation requirements (Sphinx, the ReadTheDocs theme, plus `numpy` - `conf.py` imports `pimms.CONFIG` at build time to auto-generate the keyword reference - and `versioningit`, which resolves the version shown on the index page):

```bash
pip install -r requirements.txt
```

PIMMS itself does not need to be installed or compiled: `conf.py` mocks the compiled Cython extensions and the heavy runtime dependencies (`mdtraj`, `scipy`, `dateutil`) and puts the repository root on `sys.path`, so autodoc reads the docstrings straight from the source tree. This is the same setup Read the Docs uses (see `.readthedocs.yaml`).

Once installed, you can use the `Makefile` in this directory to compile static HTML pages by

```bash
make html
```

The compiled docs are written to `_build/html`, and can be viewed by opening `_build/html/index.html`. `make clean` removes the build directory, and `make help` lists the other builders.

Note that every build regenerates `keywords.rst` from `pimms/CONFIG.py` (the same source of truth that drives `PIMMS --info`), so the keyword reference cannot drift from the code. Do not edit `keywords.rst` by hand; edit the keyword descriptions in `pimms/CONFIG.py` instead.
