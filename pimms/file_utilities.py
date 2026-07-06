## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import os
import re
import sys
import inspect
import importlib
import importlib.util

from .latticeExceptions import KeyFileException


def is_comment_line(line):
    """
    Determine whether a line is a comment (or effectively blank) line.

    A comment line is any line where the first non-whitespace character
    is a ``'#'``. A fully blank / whitespace-only line also counts as a
    comment line, so that callers can skip it uniformly.

    Parameters
    ----------
    line : str
        The line of text to inspect.

    Returns
    -------
    bool
        True if the line is a comment or blank/whitespace-only line,
        False otherwise.

    """
    hash_index = line.find('#')

    # no comment character present: only blank/whitespace-only lines count.
    # (must special-case this because str.find returns -1 when not found, and
    # line[0:-1] would silently drop the final character.)
    if hash_index == -1:
        return len(line.strip()) == 0

    # otherwise it's a comment line iff everything before the '#' is whitespace
    return len(line[0:hash_index].strip()) == 0


def remove_comments(line):
    """
    Strip an inline comment from a line and return the remaining text.

    Removes everything from the first comment character (``'#'``)
    onwards and returns the whitespace/newline-stripped string to the
    left of it.

    Parameters
    ----------
    line : str
        The line of text to process.

    Returns
    -------
    str
        The portion of the line before the first ``'#'``, stripped of
        surrounding whitespace. Returns an empty string if the line is
        entirely a comment.

    """
    return line.split('#')[0].strip()



def custom_analysis_module_import(module_name):
    """
    Load, validate and return a user-supplied ``analysis_function``.

    Loads the Python file at ``module_name`` and returns its top-level
    ``analysis_function`` callable, which PIMMS invokes as
    ``analysis_function(step, lattice)`` during the run (see the ``ANALYSIS_MODULE``
    / ``ANA_CUSTOM`` keywords). The file is imported by path under a private,
    unique module name (rather than by bare basename) so that it cannot collide
    with, or be shadowed by, a PIMMS/standard-library module or a previously
    loaded custom module of the same name. The module's own directory is added to
    ``sys.path`` so it may import helper modules that live alongside it.

    Every failure mode is reported as a clear :class:`KeyFileException` rather than
    an opaque import error, so that problems in user code are caught up front (at
    keyfile-parse time) rather than part-way through a run:

    * the path is missing or is not a file,
    * the file raises while being imported (syntax error, bad import, etc.),
    * the file does not define ``analysis_function``,
    * ``analysis_function`` is not callable, or
    * ``analysis_function`` cannot accept the ``(step, lattice)`` call signature.

    Parameters
    ----------
    module_name : str
        Path to the Python file (relative paths are resolved against the current
        working directory; ``~`` is expanded) that defines an
        ``analysis_function(step, lattice)`` function.

    Returns
    -------
    callable
        The validated ``analysis_function`` callable defined in the module.

    Raises
    ------
    KeyFileException
        If the module cannot be found, imported, or does not expose a valid
        ``analysis_function(step, lattice)`` entry point.
    """

    if not isinstance(module_name, str) or len(module_name.strip()) == 0:
        raise KeyFileException(
            "ANALYSIS_MODULE must be a path to a Python file that defines an "
            "'analysis_function(step, lattice)' function."
        )

    module_path = os.path.abspath(os.path.expanduser(module_name.strip()))

    if not os.path.isfile(module_path):
        raise KeyFileException(
            f"Could not find the custom analysis module '{module_name}' "
            f"(resolved to '{module_path}'). Relative paths are resolved against "
            "the current working directory - check the ANALYSIS_MODULE path."
        )

    dirname = os.path.dirname(module_path)
    filename = os.path.basename(module_path)

    print("[Module Analysis]: loading custom analysis from [%s]" % module_path)

    # allow the user's module to import helper modules that sit alongside it
    if dirname not in sys.path:
        sys.path.append(dirname)

    # Load by explicit file path under a private, unique module name. This avoids
    # importing the user's bare basename into the global module namespace (where it
    # could shadow, or be shadowed by, a PIMMS/stdlib module or a previously loaded
    # custom module with the same filename).
    unique_name = "pimms_custom_analysis__" + re.sub(r"\W", "_", module_path)

    try:
        spec = importlib.util.spec_from_file_location(unique_name, module_path)
        if spec is None or spec.loader is None:
            raise KeyFileException(
                f"Could not create an import spec for the custom analysis module "
                f"'{module_path}'. Is it a valid .py file?"
            )
        module = importlib.util.module_from_spec(spec)
        # register before executing so the module can reference itself during
        # import (e.g. dataclasses, typing.get_type_hints, pickling helpers)
        sys.modules[unique_name] = module
        spec.loader.exec_module(module)
    except KeyFileException:
        raise
    except Exception as e:
        # syntax error / failed import / any error raised at import time
        sys.modules.pop(unique_name, None)
        raise KeyFileException(
            f"Failed to import the custom analysis module '{module_path}': it "
            f"raised {type(e).__name__} while being loaded ({e}). Fix the error in "
            "your analysis module and try again."
        ) from e

    # the module MUST expose a top-level 'analysis_function'
    if not hasattr(module, "analysis_function"):
        available = sorted(
            n for n in vars(module)
            if not n.startswith("_") and callable(getattr(module, n, None))
        )
        hint = (f" Callables defined in the file: {', '.join(available)}."
                if available else "")
        raise KeyFileException(
            f"The custom analysis module '{filename}' does not define an "
            "'analysis_function'. PIMMS calls 'analysis_function(step, lattice)', "
            "so the file must define a top-level function with exactly that name." + hint
        )

    analysis_function = module.analysis_function

    if not callable(analysis_function):
        raise KeyFileException(
            f"'analysis_function' in '{filename}' is not callable (it is a "
            f"{type(analysis_function).__name__}). It must be a function that "
            "accepts (step, lattice)."
        )

    # confirm the signature can accept the (step, lattice) call PIMMS makes
    try:
        signature = inspect.signature(analysis_function)
    except (ValueError, TypeError):
        signature = None  # some builtins/C callables expose no signature - skip check

    if signature is not None:
        try:
            signature.bind(0, None)
        except TypeError:
            raise KeyFileException(
                f"'analysis_function' in '{filename}' has signature "
                f"{analysis_function.__name__}{signature}, but PIMMS calls it as "
                "analysis_function(step, lattice). Define it to accept exactly those "
                "two positional arguments (any extra arguments must have defaults)."
            )

    print("[Module Analysis]: 'analysis_function' loaded and validated successfully")

    return analysis_function

