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
# runtime dependencies (mdtraj, scipy). None of these are built or
# installed in the Read the Docs environment. Mocking them lets autodoc import the
# pure-Python modules for their docstrings without compiling anything. The
# keyword-reference generator below imports ``pimms.CONFIG`` while this file
# executes, before autodoc's own mocking (``autodoc_mock_imports``) is active;
# CONFIG itself needs only numpy today, but the compiled extensions are mocked
# *now* (via ``sys.modules``) so that import stays safe whatever the package
# ``__init__`` or CONFIG pulls in later.
from unittest.mock import MagicMock

_COMPILED_EXTENSIONS = [
    'pimms.hyperloop', 'pimms.inner_loops', 'pimms.inner_loops_hardwall',
    'pimms.mega_crank', 'pimms.mega_crank_fast', 'pimms.mega_crank_2D',
    'pimms.system_utils', 'pimms.cluster_kernels', 'pimms.bookkeeping',
    'pimms.lemonade.kernels._pbc',
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
copyright = "2015-2026, Alex Holehouse & Ryan Emenecker (www.holehouse.wustl.edu)"
author = 'Alex Holehouse'

# The version is read automatically so the docs always build with the current PIMMS
# version rather than a hardcoded string. We try, in order:
#   1. versioningit computed straight from the git tags. This works on Read the Docs
#      (versioningit is a docs dependency) even though RTD never installs the PIMMS
#      package, and it ignores any stale installed distribution.
#   2. the installed package metadata (the distribution is named "idptools-pimms";
#      the import package is "pimms"), for a tree with no usable git history.
#   3. "unknown".
# NB: there is no pimms/_version.py to fall back on - pyproject.toml has no
# [tool.versioningit.write] table, so versioningit never writes one.
import re

_REPO_ROOT = os.path.join(_HERE, "..")


def _get_pimms_version():
    """Return the PIMMS version string the docs should display.

    Returns
    -------
    str
        The version computed by versioningit from the git tags, or failing that the
        installed package's ``pimms.__version__``, or ``"unknown"`` if neither gives
        a real version. The pyproject default version (``"0+unknown"``) counts as
        no version at all.
    """
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

    # 2. installed package metadata
    _v = getattr(pimms, "__version__", "")
    if _v and "unknown" not in _v:
        return _v

    return "unknown"


def _get_release_date(rel):
    """Return the Month/Year the given version was released, from changelog.md.

    The changelog headers carry the release month (e.g. ``## 1.0.0 (July 2026)``). We
    look up the header for ``rel`` itself first. Failing that, a post-release build
    (``1.0.7.post8``, ``1.0.7.post8+g1e71170``, or a ``.devN`` of one) is dated by the
    release it is built on top of, so we strip the ``.postN`` / ``.devN`` / ``+local``
    segments and look up that base release (``1.0.7``). We never fall back to the
    newest header: the top of the changelog is usually the next, not-yet-released
    version, and dating a build with it is exactly the wrong answer. A bare ``.devN``
    build (``1.0.8.dev3``) is a pre-release of a version that has not shipped, so it
    gets no date either.

    Parameters
    ----------
    rel : str
        The version string being documented (normally ``release`` below, which has
        already had any ``+local`` segment removed; one left on is ignored).

    Returns
    -------
    str
        The text inside the parentheses of the matching changelog header (e.g.
        ``"August 2026"``), or ``""`` if ``rel`` is ``"unknown"``, the changelog
        cannot be read, or neither ``rel`` nor its base release has a dated header.
    """
    if rel == "unknown":
        return ""
    changelog = os.path.join(_REPO_ROOT, "changelog.md")
    try:
        with open(changelog) as _fh:
            _text = _fh.read()
    except OSError:
        return ""

    def _header_date(ver):
        """Return the date of the ``## <ver> (<Month Year>)`` header, or "".

        Parameters
        ----------
        ver : str
            Version to look up. It is matched as a whole word, so ``1.0.7`` does not
            pick up a ``1.0.70`` header.

        Returns
        -------
        str
            The text inside the header's parentheses, or ``""`` if there is no such
            header.
        """
        _m = re.search(
            r"^##\s+" + re.escape(ver) + r"\s+\(([^)]+)\)", _text, re.MULTILINE
        )
        return _m.group(1).strip() if _m else ""

    _date = _header_date(rel)
    if _date:
        return _date

    # base release of a post-release build; a bare .devN (no .postN) is a pre-release
    # of an unreleased version and is deliberately not matched
    _m = re.fullmatch(r"(\d+(?:\.\d+)*)\.post\d+(?:\.dev\d+)?(?:\+.*)?", rel.strip())
    if _m:
        return _header_date(_m.group(1))
    _m = re.fullmatch(r"(\d+(?:\.\d+)*)\+.*", rel.strip())
    if _m:
        return _header_date(_m.group(1))
    return ""


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
# docs environment (mdtraj, scipy). Mocking these keeps the docs build free of
# any compilation step; only the pure-Python modules' own docstrings are rendered.
autodoc_mock_imports = _COMPILED_EXTENSIONS + ['mdtraj', 'scipy']

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
