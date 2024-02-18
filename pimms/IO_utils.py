## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2024
## ...........................................................................
# 

## Arbitrary non-specific I/O functions
##
##

from .latticeExceptions import IOException
from .CONFIG import TERMINAL_WIDTH
from os import path


    
    
# ............................................................
#
def wipe_file(filename, header=None):
    """
    Simple function to wipe the contents from a file. 

    Parameters
    ----------

    filename : string
        filename as string (absolute or relative path)

    header : str
        String to write to the top of the file (optional).
        Default is None.

    Returns
    -------
    None
        No return variable, but if the file exists it's contents
        will be erased. If the file does not exist an empty version
        of the file will be generated.

    """


    """
    if os.path.exists("demofile.txt"):
        os.remove("demofile.txt")
    else:
        pass
    """

    # wipe the file
    with open(filename, 'w') as fh:
        if header is not None:
            fh.write(header)
        else:
            fh.write("")
            
        
            
        
    



        
# ............................................................
#
def write_list_to_file(contents, filename, mode='w'):
    """
    Function that writes the contents of a list out to a new file,
    overwriting any content that was in previously
    
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
    Function that writes a newline to output
    """

    print('\n'*(nlines-1))


# ............................................................
#
def horizontal_line(hzlen=TERMINAL_WIDTH, linechar='-',leader=''):
    """
    Function that prints a horizontal line to output

    Parameters
    -------------

    hzlen : int
        Line length
    
    """

    print(leader+linechar*hzlen)


# ............................................................
#
def status_message(msg, msg_type='info'):
    """
    Function that prints a status message to stdout with an
    associated header. Also ensures the width matches the
    TERMINAL_WIDTH which is hardcoded in CONFIG.py

    Parameters
    ------------

    msg : string
        message of interest

    msg_type : {'startup', 'info', 'warning', 'error','major', 'vanilla', 'update', 'null''}
        Mode for message

    """
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
            IOException("[THIS IS A BUG] 'major' messages must be shorter than %i"%(TERMINAL_WIDTH))

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

    maxlinelength : int
        Maximum line length (default is TERMINAL_WIDTH, defined in CONFIG.py)

    multiline_leader : str
        String to be printed at the start of each line for the 2nd line onwards, 
        useful for indenting multiline strings.

    print_to_stdout : bool
        If True, the string is printed to stdout, if False the string is returned
        as a string.

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
            
        
