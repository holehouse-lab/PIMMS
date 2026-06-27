## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................


import os
import sys
import importlib


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
    Dynamically import a user module and return its ``analysis_function``.

    Splits the supplied path into a directory and filename, appends the
    directory to ``sys.path`` (defaulting to the current working
    directory if no directory component is present), imports the module
    by name, and returns its ``analysis_function`` attribute. This is
    used to load user-supplied custom analysis code at runtime.

    Parameters
    ----------
    module_name : str
        Path to the Python module file (with or without directory
        components and ``.py`` extension) that defines an
        ``analysis_function``.

    Returns
    -------
    callable
        The ``analysis_function`` callable defined within the imported
        module.

    """

    dirname, filename = os.path.split(module_name)
    
    print("[Module Analysis]: Reading [%s] from directory [%s]" %(filename, dirname))

    # if dirname empty assume relative path
    if len(dirname) == 0:
        dirname = os.getcwd()
    sys.path.append(dirname)
    mname = filename.split('.')[0]
    MOD = importlib.import_module(mname)
    return MOD.analysis_function
        
        
