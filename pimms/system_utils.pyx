## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2024
## 
## ...........................................................................

import numpy as np
from pimms import pimmslogger
from pimms import IO_utils
cimport numpy as cnp
cnp.import_array()

cimport cython 
import pimms.inner_loops as inner_loops

from pimms.cython_config cimport NUMPY_INT_TYPE
from pimms.CONFIG import NP_INT_TYPE as NUMPY_INT_TYPE_PYTHON


##
#################################################################################################
##
def check_dtype_consistency():
    """
    Function that ensures that the cython and python integers match one another.
    """

    try:
        __build_array()
    except:
        
        msg       = 'Error in PIMMS setup: Python and Cython intsize are not consistent. Please ensure that CONFIG.NP_INT_TYPE is set to the same value as cython_config.NUMPY_INT_TYPE. NOTE: You will need to recompile the cython code after making changes to cython_config.NUMPY_INT_TYPE'

        IO_utils.status_message(msg, msg_type="error")
        pimmslogger.log_error(IO_utils.stdout(msg, maxlinelength=115, multiline_leader="                                  ", print_to_stdout=False))
        exit(1)


def __build_array():
    """
    Function that builds a numpy array using the cython integer type.
    """
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] test = np.zeros(10, dtype=NUMPY_INT_TYPE_PYTHON)
    
        
def check_beads_to_grid_mapping(chains_list):
    """
    Function that checks that number of beads to be created can actually be accomodated by grids using
    the current integer type. 

    Parameters
    ----------
    chains_list : list
        List of chains to be created. Each element of the list is a tuple containing the number of beads
        in the chain and the chain sequence.
    """

    bits = NUMPY_INT_TYPE_PYTHON().itemsize*8

    # define largest possible possitive integer; note we could cram in more
    # if we REALLY wanted to by using unsigned integers, but let's worry abut
    # that if we need to....
    max_num_beads = 2**(bits-1) - 1
    
    bead_count = 0
    for x in chains_list:
        bead_count = bead_count + x[0]*len(x[1])

        
    if bead_count > max_num_beads:            
        msg = f"Error in PIMMS setup: The number of beads in the system ({bead_count} beads) exceeds the maximum number of beads that can be accomodated by the current integer type (signed int{bits}; max number of beads = {max_num_beads}). Please consider using a larger integer type."
        IO_utils.status_message(msg, msg_type="error")
        pimmslogger.log_error(IO_utils.stdout(msg, maxlinelength=115, multiline_leader="                                  ", print_to_stdout=False))
        exit(1)

    
    
