## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................
# 

## Arbitrary non-specific I/O functions
##
##

import os

from .latticeExceptions import IOException
from .CONFIG import TERMINAL_WIDTH




# ............................................................
#
def wipe_file(filename):
    """
    Simple function to wipe the contents from a file.

    Parameters
    ----------

    filename : string
        filename as string (absolute or relative path)

    Returns
    -------
    None
        No return variable, but if the file exists it's contents
        will be erased. If the file does not exist an empty version
        of the file will be generated.

    """

    # wipe the file
    with open(filename, 'w') as fh:
        fh.write("")


# ............................................................
#
def remove_files(filenames):
    """
    Delete a set of files, ignoring any that are absent.

    This is what start-up uses to clear a previous run's output out of the
    working directory. Note that it deletes rather than truncates: analysis
    output files are created lazily (they exist only once something has been
    written to them), so leaving a truncated file behind would advertise an
    output the current run may never write.

    Failures to delete (a permissions problem, say) are ignored rather than
    raised - a file we cannot remove is not a reason to refuse to start a
    simulation.

    Parameters
    ----------
    filenames : iterable of str
        Filenames to delete (absolute or relative paths).

    Returns
    -------
    list
        The filenames that were actually deleted, in the order they were
        passed. Useful for reporting and for testing.
    """

    removed = []

    for filename in filenames:
        if os.path.exists(filename):
            try:
                os.remove(filename)
                removed.append(filename)
            except OSError:
                pass

    return removed


        
# ............................................................
#
def write_list_to_file(contents, filename, mode='w'):
    """
    Write the contents of a list out to a file.

    Each element of ``contents`` is written verbatim (no newlines are
    added), so callers are responsible for including any trailing
    newline characters they require.

    Parameters
    ----------
    contents : list of str
        The lines/strings to write to the file, in order.

    filename : str
        Name of the file to write to (absolute or relative path).

    mode : {'w', 'a'}, optional
        File open mode. ``'w'`` overwrites any existing content
        (default), ``'a'`` appends to it.

    Returns
    -------
    None
        No return value; the file is written to disk.

    Raises
    ------
    IOException
        If ``mode`` is not ``'w'`` or ``'a'``.

    """
    if mode not in ['w','a']:
        raise IOException("write_list_to_file requires either 'a' or 'w' to be massed as mode")
        

    with open(filename,mode) as fh:
        for line in contents:
            fh.write(line)

# ............................................................
#
def newline(nlines=1):
    """
    Print blank line(s) to stdout.

    Note that because ``print`` itself emits a trailing newline, calling
    this with the default ``nlines=1`` prints a single blank line.

    Parameters
    ----------
    nlines : int, optional
        Controls how many newline characters are emitted (``nlines - 1``
        explicit newlines plus the one ``print`` adds). Default is 1.

    Returns
    -------
    None
        No return value; output is written to stdout.

    """

    print('\n'*(nlines-1))


# ............................................................
#
def horizontal_line(hzlen=TERMINAL_WIDTH, linechar='-',leader=''):
    """
    Print a horizontal divider line to stdout.

    Parameters
    -------------

    hzlen : int, optional
        Line length, i.e. the number of times ``linechar`` is repeated.
        Defaults to ``TERMINAL_WIDTH`` (from CONFIG.py).

    linechar : str, optional
        The character used to draw the line. Default is ``'-'``.

    leader : str, optional
        String printed immediately before the line (e.g. for
        indentation). Default is an empty string.

    Returns
    -------
    None
        No return value; output is written to stdout.

    """

    print(leader+linechar*hzlen)


# ............................................................

# Module-level reduced-printing flag. Simulation.__init__ sets this from the
# REDUCED_PRINTING keyword so that deep call sites (the move functions, restart
# writer) can suppress their per-event chatter without threading the flag
# through every signature. status_message(..., allow_suppress=True) honours it.
REDUCED_PRINTING = False


def set_reduced_printing(value):
    """
    Record the run's REDUCED_PRINTING setting for per-event message gating.

    Sets the module-level ``REDUCED_PRINTING`` flag that
    ``status_message(..., allow_suppress=True)`` consults, so that deep
    call sites can drop their per-event chatter without the flag being
    threaded through their signatures.

    Parameters
    ----------
    value : bool
        The run's REDUCED_PRINTING setting. Cast to a bool before being
        stored, so any truthy/falsey value is accepted.

    Returns
    -------
    None
        No return value; the module-level flag is updated.

    """
    global REDUCED_PRINTING
    REDUCED_PRINTING = bool(value)


#
def status_message(msg, msg_type='info', allow_suppress=False):
    """
    Function that prints a status message to stdout with an
    associated header. Also ensures the width matches the
    TERMINAL_WIDTH which is hardcoded in CONFIG.py

    Parameters
    ------------

    msg : string
        message of interest

    msg_type : {'startup', 'info', 'warning', 'error', 'major', 'vanilla', 'update', 'null'}, optional
        Mode for message, which selects the header/prefix and formatting
        used when printing. Default is 'info'.

    allow_suppress : bool, optional
        Flag which, if set to True, marks this as a per-event message
        that is skipped entirely when the run was started with
        REDUCED_PRINTING (see set_reduced_printing()). Messages left at
        the default of False always print.

    Returns
    -------
    None
        No return value; the formatted message is printed to stdout.

    Raises
    ------
    IOException
        If ``msg_type`` is ``'major'`` but the message is longer than
        ``TERMINAL_WIDTH``, or if an unrecognised ``msg_type`` is passed.

    """

    # per-event messages opt in to suppression under REDUCED_PRINTING (the
    # per-move/per-restart chatter that previously bypassed the keyword)
    if allow_suppress and REDUCED_PRINTING:
        return
    leader='             '

    if msg_type == 'vanilla':
        stdout(msg)

    elif msg_type == 'null':
        leader='      '
        stdout(msg, multiline_leader=leader)
        
    elif msg_type == 'startup':

        s    = '  [STARTUP]: ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'info':
        s =    '  [INFO]:    ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'warning':
        s    = '  [WARNING]: ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'error':
        s =    '  [ERROR]:   ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'update':
        s =    '  [UPDATE]:  ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'major':
        
        # 'major' messages sandwhiches a message between two lines but MUST
        # be shorter than the TERMINAL_WIDTH.
        #
        # use for section headers

        sl=len(msg)
        if sl > TERMINAL_WIDTH:
            raise IOException("[THIS IS A BUG] 'major' messages must be shorter than %i"%(TERMINAL_WIDTH))

        line = '.'*TERMINAL_WIDTH

        print('')
        print(msg)
        print(line)
        print('')

    else:
        raise IOException("[THIS IS A BUG] Invalid msg_type passed to status_message")
        

# ............................................................
#
def stdout(string, maxlinelength=TERMINAL_WIDTH, multiline_leader='', print_to_stdout=True):
    """
    Function that prints a string to stdout, but ensures that
    the string is wrapped to a maximum line length.

    Parameters
    ------------
    string : str
        String to be printed

    maxlinelength : int, optional
        Maximum line length (default is TERMINAL_WIDTH, defined in CONFIG.py)

    multiline_leader : str, optional
        String to be printed at the start of each line for the 2nd line onwards,
        useful for indenting multiline strings. Default is an empty string.

    print_to_stdout : bool, optional
        If True (the default), the string is printed to stdout, if False the
        string is returned as a string.

    Returns
    ---------
    None or str
        If print_to_stdout is True, then the function returns None, otherwise
        it returns the string that would have been printed to stdout.


    """

    full_string=''
    newstring=''

    # count 
    c = -1


    # for each character in the string
    for s in string:

        # increment the counter
        c = c + 1

        # if a newline is encountered we add his to our
        # ever growing string and reset the newstring
        if repr(string[c:c+1]) == repr('\n'):
            full_string = full_string + newstring + '\n'
            newstring = multiline_leader
            continue
        
        # skip leading whitespace on a line
        if newstring == multiline_leader:
            if s == ' ':
                continue

        # add the current character to the newstring
        newstring = newstring + s

        # if newstring is as long as maxlinelength then
        if len(newstring) >= maxlinelength:

            try:
                
                # if next position is not whitespace we 
                # skip and iterate over each charcter so we
                # don't split words
                if string[c+1] != ' ':
                    continue

            # if we're at the end of the string we're done
            except IndexError:
                break

            # if we get here then the next character was
            # whitespace, so we can split on a new line
            
            full_string = full_string + newstring + '\n'

            # reset newstring
            newstring = multiline_leader

    # add the last newstring to the full_string
    full_string = full_string + newstring

    if print_to_stdout:
        print(full_string)
    else:
        return full_string
            
        
    

# ............................................................
