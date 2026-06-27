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
    A comment line is any line where the first readable character is a '#'
    (a fully blank line also counts, so it is skipped by callers).

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
    Removes all the the content after a comment character
    and returns the whitespace/newline stripped string
    left of the comment charater ('#')

    """
    return line.split('#')[0].strip()



def custom_analysis_module_import(module_name):
    """
    Function that reads in and returns the 'analysis_function'
    function from the passed module

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
        
        
