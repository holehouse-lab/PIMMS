## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2026
## 
## ...........................................................................

import numpy as np
from pimms import pimmslogger
from pimms import IO_utils
from pimms.latticeExceptions import SimulationException
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

    Builds a throwaway array using the Python-side integer dtype and asks Cython
    to bind it to a typed buffer of the compiled integer type. If the two do not
    agree the buffer acquisition fails, and we turn that into a readable setup
    error rather than an opaque ValueError deep inside a kernel.

    Returns
    -------
    None

    Raises
    ------
    SimulationException
        If CONFIG.NP_INT_TYPE and cython_config.NUMPY_INT_TYPE do not match.
    """

    try:
        __build_array()
    except Exception as exc:
        
        msg       = 'Error in PIMMS setup: Python and Cython intsize are not consistent. Please ensure that CONFIG.NP_INT_TYPE is set to the same value as cython_config.NUMPY_INT_TYPE. NOTE: You will need to recompile the cython code after making changes to cython_config.NUMPY_INT_TYPE'

        IO_utils.status_message(msg, msg_type="error")
        pimmslogger.log_error(IO_utils.stdout(msg, maxlinelength=115, multiline_leader="                                  ", print_to_stdout=False))
        raise SimulationException(msg) from exc


def __build_array():
    """
    Function that builds a numpy array using the cython integer type.

    The array itself is allocated with the Python-side dtype, so the typed
    buffer assignment only succeeds if the two integer types agree. Used solely
    as the probe behind check_dtype_consistency().

    Returns
    -------
    None

    Raises
    ------
    ValueError
        Raised by the buffer machinery if the Python dtype does not match the
        compiled Cython integer type.
    """
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] test = np.zeros(10, dtype=NUMPY_INT_TYPE_PYTHON)
    
        
def check_beads_to_grid_mapping(chains_list):
    """
    Check that every chain ID can be represented by the occupancy-grid dtype.

    The main grid stores ``chainID`` at each occupied site, not a globally
    unique bead ID.  The relevant upper bound is therefore the total number of
    chains, irrespective of how many beads each chain contains.

    Parameters
    ----------
    chains_list : list of [int, str]
        One entry per chain type, as parsed from the CHAIN (and, where relevant,
        EXTRA_CHAIN) keywords: element 0 is the number of copies of that type,
        element 1 is the sequence. Only the copy counts are read here.

    Returns
    -------
    None

    Raises
    ------
    SimulationException
        If the total number of chains exceeds the largest positive value the
        occupancy-grid integer type can hold.
    """

    bits = NUMPY_INT_TYPE_PYTHON().itemsize*8

    # define largest possible possitive integer; note we could cram in more
    # if we REALLY wanted to by using unsigned integers, but let's worry abut
    # that if we need to....
    max_chain_id = 2**(bits-1) - 1
    
    chain_count = 0
    for x in chains_list:
        chain_count = chain_count + x[0]

        
    if chain_count > max_chain_id:
        msg = f"Error in PIMMS setup: The number of chains in the system ({chain_count}) exceeds the maximum chainID that can be represented by the occupancy-grid integer type (signed int{bits}; max chainID = {max_chain_id}). Please consider using a larger integer type."
        IO_utils.status_message(msg, msg_type="error")
        pimmslogger.log_error(IO_utils.stdout(msg, maxlinelength=115, multiline_leader="                                  ", print_to_stdout=False))
        raise SimulationException(msg)

    
    
