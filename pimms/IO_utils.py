## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

## Arbitrary non-specific I/O functions
##
##
import importlib
from .latticeExceptions import IOException

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


def status_message(msg, msg_type='vanila'):

    # do no formatting
    if msg_type == 'vanila':
        print(msg)


def stdout(string, maxlinelength):
    
    newstring=''
    out=False
    for s in string:

        # ignore leading whitespace...
        if len(newstring) == 0:
            if s == ' ' :
                continue
            
        newstring =newstring + s
        if len(newstring) >= maxlinelength:
            out = True
            if s != ' ':
                continue
            else:
                print(newstring)
                newstring=''
                out=False
    print(newstring)
            
        
