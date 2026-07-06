"""
Generate the keyword-reference page (``keywords.rst``) from PIMMS' single source
of truth in ``pimms/CONFIG.py``.

The same ``CONFIG.KEYWORD_GROUPS`` / ``CONFIG.KEYWORDS_DESCRIPTION`` /
``CONFIG.DEFAULTS`` data structures drive the ``PIMMS --info`` command, so the
documentation can never drift from the CLI: regenerate this page whenever a
keyword changes. ``conf.py`` runs :func:`generate` automatically at the start of
every Sphinx build, but it can also be run by hand::

    python docs/generate_keywords.py

This file does NOT contain any keyword descriptions itself - edit the descriptions
in ``pimms/CONFIG.py`` (``KEYWORDS_DESCRIPTION``) and regenerate.
"""
import os

from pimms import CONFIG


_INTRO = """\
.. _keyword-reference:

=================
Keyword reference
=================

.. note::

   This page is auto-generated from ``pimms/CONFIG.py`` (the same source that
   powers ``PIMMS --info``). Do not edit it by hand - edit the descriptions in
   ``CONFIG.KEYWORDS_DESCRIPTION`` and rebuild the docs.

A PIMMS simulation is configured with a plain-text **keyfile**. Each line sets one
keyword::

    KEYWORD : value

Lines beginning with ``#`` are comments, and blank lines are ignored. Most
keywords may appear at most once; the ``CHAIN`` and ``EXTRA_CHAIN`` keywords are
the exception and may be repeated to define a multi-component system. The
required keywords are marked **required** below; everything else falls back to
the default shown.

You can query any keyword from the command line without opening this page::

    PIMMS --info                 # list every keyword, grouped
    PIMMS --info <KEYWORD>       # full details on one keyword
    PIMMS --info ALL             # print every keyword description

The keywords are grouped by purpose below.

"""


def _format_default(keyword):
    """
    Return a human-readable default string for ``keyword`` (or ``None``).

    Parameters
    ----------
    keyword : str
        The keyword to look up in ``CONFIG.DEFAULTS``.

    Returns
    -------
    str or None
        A display string for the default value, or ``None`` if the keyword has
        no entry in ``CONFIG.DEFAULTS``.
    """
    if keyword not in CONFIG.DEFAULTS:
        return None
    value = CONFIG.DEFAULTS[keyword]
    if isinstance(value, list) and len(value) == 0:
        return "none"
    if value == 'UNSET':
        return "unset"
    if value == 'N/A':
        return "none"
    return str(value)


def _term_line(keyword):
    """
    Build the reStructuredText definition-list term for a keyword.

    Parameters
    ----------
    keyword : str
        The keyword name.

    Returns
    -------
    str
        A term line of the form ````KEYWORD`` (*type*, default: ``value``)`` (or
        ``*required*`` when the description is flagged ``[REQUIRED]``).
    """
    type_string, description = CONFIG.KEYWORDS_DESCRIPTION[keyword]
    if "[REQUIRED]" in description:
        qualifier = "**required**"
    else:
        default = _format_default(keyword)
        qualifier = f"default: ``{default}``" if default is not None else "optional"
    return f"``{keyword}`` (*{type_string}*, {qualifier})"


def generate(output_path):
    """
    Write the keyword-reference reStructuredText page.

    Iterates ``CONFIG.KEYWORD_GROUPS`` (preserving order), emitting one section per
    group and a definition-list entry per keyword (term = name/type/default,
    definition = description). Any documented keyword not assigned to a group is
    collected under a trailing "Other" section so nothing is silently dropped.

    Parameters
    ----------
    output_path : str
        Filesystem path of the ``keywords.rst`` file to (over)write.

    Returns
    -------
    None
        The file is written to ``output_path``.
    """
    lines = [_INTRO]
    seen = set()

    groups = list(CONFIG.KEYWORD_GROUPS)
    leftover = [k for k in CONFIG.KEYWORDS_DESCRIPTION
                if k not in {kw for _, kws in groups for kw in kws}]
    if leftover:
        groups = groups + [("Other", leftover)]

    for heading, keywords in groups:
        present = [k for k in keywords if k in CONFIG.KEYWORDS_DESCRIPTION]
        if not present:
            continue
        lines.append(heading)
        lines.append("-" * len(heading))
        lines.append("")
        for keyword in present:
            # escape the reStructuredText substitution character so descriptions
            # containing a bare "|" (e.g. an absolute-value notation) do not trip
            # the parser.
            description = CONFIG.KEYWORDS_DESCRIPTION[keyword][1].replace("|", r"\|")
            lines.append(_term_line(keyword))
            # indent the description as the definition body (4 spaces)
            for paragraph in description.split("\n"):
                lines.append(f"    {paragraph.strip()}" if paragraph.strip() else "")
            lines.append("")
            seen.add(keyword)
        lines.append("")

    with open(output_path, "w") as fh:
        fh.write("\n".join(lines).rstrip() + "\n")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    target = os.path.join(here, "keywords.rst")
    generate(target)
    print(f"Wrote {target} ({len(CONFIG.KEYWORDS_DESCRIPTION)} keywords)")
