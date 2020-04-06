## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2020
## ...........................................................................
# 

## Arbitrary non-specific I/O functions
##
##

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


    """
    if os.path.exists("demofile.txt"):
        os.remove("demofile.txt")
    else:
        pass
    """

    with open(filename, 'w') as fh:
        fh.write('')

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

    msg_type : {'startup', 'info', 'warning', 'error','major', 'vanilla', 'update''}
        Mode for message

    """

    if msg_type == 'vanilla':
        stdout(msg)
        

    elif msg_type == 'startup':
        leader='             '
        s    = '  [STARTUP]: ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'info':
        leader='             '
        s =    '  [INFO]:    ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'warning':
        leader='             '
        s    = '  [WARNING]: ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'error':
        leader='             '
        s =    '  [ERROR]:   ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'update':
        leader='             '
        s =    '  [UPDATE]:  ' + msg
        stdout(s, multiline_leader=leader)

    elif msg_type == 'major':
        
        # 'major' messages sandwitch a message between two lines but MUST
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
def stdout(string, maxlinelength=TERMINAL_WIDTH, multiline_leader=''):
    
    newstring=''
    for s in string:

        newstring = newstring + s
        if len(newstring) >= maxlinelength:
            print(newstring)
            newstring = multiline_leader

    print(newstring)
            
        
