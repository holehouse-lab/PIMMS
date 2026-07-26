# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# Incase the project was not installed
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, '..')))   # repo root, so `import pimms` works
sys.path.insert(0, _HERE)                                        # so `import generate_keywords` works

# ---------------------------------------------------------------------------
# Mock compiled extensions and heavy runtime dependencies.
#
# PIMMS' hot loops are compiled Cython extensions, and a few modules import heavy
# runtime dependencies (mdtraj, scipy, dateutil). None of these are built or
# installed in the Read the Docs environment. Mocking them lets autodoc import the
# pure-Python modules for their docstrings without compiling anything, and (crucially)
# lets the keyword-reference generator below import ``pimms.CONFIG``, which imports the
# compiled ``get_randmax`` at module-load time. The compiled extensions must be mocked
# *now* (via ``sys.modules``) because that import happens while this file executes,
# before autodoc's own mocking (``autodoc_mock_imports``) is active.
from unittest.mock import MagicMock

_COMPILED_EXTENSIONS = [
    'pimms.get_randmax', 'pimms.hyperloop', 'pimms.inner_loops',
    'pimms.inner_loops_hardwall', 'pimms.lattice_tools', 'pimms.mega_crank',
    'pimms.mega_crank_fast', 'pimms.mega_crank_2D', 'pimms.random_number',
    'pimms.system_utils', 'pimms.cluster_kernels', 'pimms.lemonade.kernels._pbc',
]
for _name in _COMPILED_EXTENSIONS:
    sys.modules.setdefault(_name, MagicMock())

import pimms

# Regenerate the keyword-reference page from CONFIG (the single source of truth
# that also drives `PIMMS --info`) at the start of every build, so the docs never
# drift from the code.
import generate_keywords
generate_keywords.generate(os.path.join(_HERE, 'keywords.rst'))


# -- Project information -----------------------------------------------------

project = 'PIMMS'
copyright = "2016-2026, Alex Holehouse & Ryan Emenecker (www.holehouse.wustl.edu)"
author = 'Alex Holehouse'

# The version is read automatically so the docs always build with the current PIMMS
# version rather than a hardcoded string. We try, in order:
#   1. versioningit computed straight from the git tags. This works on Read the Docs
#      (versioningit is a docs dependency) even though RTD never installs the PIMMS
#      package, and it ignores any stale installed distribution.
#   2. the versioningit-written pimms/_version.py, for built/installed trees (e.g. an
#      sdist) that have no .git directory. NOTE this file is gitignored, so it is
#      absent on a fresh clone - hence versioningit is tried first.
#   3. the installed package metadata (the distribution is named "idptools-pimms";
#      the import package is "pimms").
#   4. "unknown".
import re

_REPO_ROOT = os.path.join(_HERE, "..")


def _get_pimms_version():
    # 1. versioningit from git
    try:
        import versioningit

        _v = versioningit.get_version(project_dir=_REPO_ROOT)
        # reject the pyproject default-version ("0+unknown") used when git/tags cannot
        # be resolved (e.g. a too-shallow clone with no reachable tag).
        if _v and "unknown" not in _v:
            return _v
    except Exception:
        pass

    # 2. versioningit-written _version.py
    try:
        with open(os.path.join(_REPO_ROOT, "pimms", "_version.py")) as _fh:
            _m = re.search(r"""__version__\s*=\s*['"]([^'"]+)['"]""", _fh.read())
            if _m:
                return _m.group(1)
    except OSError:
        pass

    # 3. installed package metadata
    _v = getattr(pimms, "__version__", "")
    if _v and "unknown" not in _v:
        return _v

    return "unknown"


def _get_release_date(rel):
    """Month/Year the given version was released, from changelog.md.

    The changelog headers carry the release month (e.g. ``## 1.0.0 (July 2026)``). We
    look up the entry matching ``rel`` and, failing that, fall back to the most recent
    versioned entry - so a development build (``1.0.2.post1+g47fe7be``) still reports the
    date of the release it is built on top of. Returns "" if nothing is found.
    """
    if rel == "unknown":
        return ""
    changelog = os.path.join(_REPO_ROOT, "changelog.md")
    try:
        with open(changelog) as _fh:
            _text = _fh.read()
    except OSError:
        return ""
    # exact match: "## <rel> (<Month Year>)"
    _m = re.search(
        r"^##\s+" + re.escape(rel) + r"\s+\(([^)]+)\)", _text, re.MULTILINE
    )
    if _m:
        return _m.group(1).strip()
    # fall back to the first "## X.Y.Z (<Month Year>)" header
    _m = re.search(r"^##\s+\d[\w.]*\s+\(([^)]+)\)", _text, re.MULTILINE)
    return _m.group(1).strip() if _m else ""


# The full version, including alpha/beta/rc tags (PEP 440 local segment, e.g.
# "+g47fe7be.d20260726", is dropped for display); the short X.Y version is derived
# from it.
release = _get_pimms_version().split("+")[0]
version = ".".join(release.split(".")[:2]) if release != "unknown" else release

# The release Month/Year (from changelog.md). Exposed to .rst as the |version_info|
# substitution below (version, optionally with the date).
release_date = _get_release_date(release)
if release_date:
    _version_info = f"{release} (released {release_date})"
else:
    _version_info = release
rst_prolog = f".. |version_info| replace:: {_version_info}\n"


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

# When autodoc imports each documented module it must not fail on the compiled Cython
# kernels (mocked above) or on heavy runtime dependencies that are not installed in the
# docs environment (mdtraj, scipy, dateutil). Mocking these keeps the docs build free of
# any compilation step; only the pure-Python modules' own docstrings are rendered.
autodoc_mock_imports = _COMPILED_EXTENSIONS + ['mdtraj', 'scipy', 'dateutil']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'default'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Project logo, shown at the top of the navigation sidebar (the contents list).
# Referenced from the repo's branding/ directory so there is a single source of truth.
html_logo = os.path.join('..', 'branding', 'logo.png')

html_theme_options = {
    'logo_only': False,        # show the project name under the logo as well
    'style_external_links': True,
}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Brand-colour overrides for the sphinx_rtd_theme accent (see _static/custom.css).
html_css_files = ['custom.css']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'pimmsdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'pimms.tex', 'pimms Documentation',
     'pimms', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'pimms', 'pimms Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'pimms', 'pimms Documentation',
     author, 'pimms', 'Lattice simulation package for biomolecule',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
