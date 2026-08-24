# Compiling pimms's Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/).
To compile the docs, first install the documentation requirements (Sphinx, the
ReadTheDocs theme, plus `numpy` - `conf.py` imports `pimms.CONFIG` at build time
to auto-generate the keyword reference - and `versioningit`, which resolves the
version shown on the index page):

```bash
pip install -r requirements.txt
```


Once installed, you can use the `Makefile` in this directory to compile static HTML pages by
```bash
make html
```

The compiled docs will be in the `_build` directory and can be viewed by opening `index.html` (which may itself 
be inside a directory called `html/` depending on what version of Sphinx is installed).