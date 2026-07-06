## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

import datetime

def message_preprocess(msg):
    """
    Wrap an error message in a standardized PIMMS crash banner.

    Decorates the passed message with a timestamped header and footer
    (including bug-reporting instructions) so that exceptions raised
    throughout PIMMS present a consistent, easily-spotted crash report.

    Parameters
    ----------
    msg : str
        The raw error message to embed within the banner.

    Returns
    -------
    str
        The formatted multi-line message, ready to be passed to an
        exception constructor.

    """
    s = ''
    s = s + '\n\n#############################################\n%s\n'%((datetime.datetime.now()))
    s = s + 'OH NOOOOO! PIMMS HAS CRASHED!\n\n'
    s = s + "See above for the traceback for debugging information - please retain\nthis information if you think PIMMS crashed due to a bug.\n\n"    
    s = s + "Report bugs by contacting Alex directly (alex.holehouse@wustl.edu)\n\n"
    s = s + "Error message:\n"
    s = s + msg
    s = s + '\n############################################\n'

    return s
          


##
## SIMULATION EXCEPTIONS
##

class SimulationException(Exception):
    """Raised when a general, non-specific error occurs while running a simulation."""
    pass

class SimulationEnergyException(Exception):
    """Raised when an energy inconsistency is detected during a simulation (e.g. tracked vs. recomputed energy mismatch)."""
    pass

class AcceptanceException(Exception):
    """Raised when an error occurs in the Monte Carlo move acceptance/rejection logic."""
    pass

class RestartException(Exception):
    """Raised when a restart object/file is inconsistent with the requested simulation setup."""
    pass



##
## ENERGY EXCEPTIONS
##

class EnergyException(Exception):
    """Raised when an error occurs while constructing or evaluating the system Hamiltonian/energy."""
    pass

class TemperatureException(Exception):
    """Raised when an invalid or inconsistent temperature is encountered."""
    pass


##
## LATTICE EXCEPTIONS

class LatticeInitializationException(Exception):
    """Raised when a lattice cannot be constructed consistently (e.g. mismatched dimensions or grids)."""
    pass

class ChainInitializationException(Exception):
    """Raised when a chain cannot be successfully initialized onto the lattice."""
    pass

class ChainInsertionFailure(Exception):
    """Raised when a chain (or chain segment) cannot be inserted into the lattice, typically due to overlaps."""
    pass

class ChainDeletionFailure(Exception):
    """Raised when a chain (or chain segment) cannot be removed from the lattice as expected."""
    pass

class ChainAugmentFailure(Exception):
    """Raised when a chain cannot be grown/augmented (e.g. no valid free neighbouring sites)."""
    pass

class ChainConnectivityError(Exception):
    """Raised when a chain's beads are found to violate the expected nearest-neighbour connectivity."""
    pass

class TypeGridException(Exception):
    """Raised when an inconsistency is detected while updating the bead-type grid."""
    pass


class ClusterSizeThresholdException(Exception):
    """
    Exception raised by get_all_chains_in_connected_component iff a threshold
    cluster sized is reached.
    """
    pass



##
## FILE EXCEPTIONS
##

class KeyFileException(Exception):
    """Raised when an error occurs while parsing or validating a simulation keyfile."""
    pass

class ParameterFileException(Exception):
    """Raised when an error occurs while parsing or validating an energy parameter file."""
    pass

class PDBException(Exception):
    """Raised when an error occurs while constructing or writing a PDB file."""
    pass

class IOException(Exception):
    """Raised when a general input/output error occurs (e.g. invalid file write mode)."""
    pass


##
## MOVE EXCEPTIONS
## 

class MoveSetException(Exception):
    """Raised when an invalid or inconsistent set of Monte Carlo moves is configured."""
    pass

class MoveException(Exception):
    """Raised when an error occurs while constructing or executing a Monte Carlo move."""
    pass

class CustomInitializationException(Exception):
    """Raised when a user-supplied custom system initialization fails or is inconsistent."""
    pass

##
## ANALYSIS EXCEPTIONS
class AnalysisStructureException(Exception):
    """Raised when an error occurs in the structural analysis machinery."""
    pass

class AnalysisRoutineException(Exception):
    """Raised when an error occurs while executing an analysis routine."""
    pass


##
## LATTICE UTILS
class RotationException(Exception):
    """Raised when an invalid rotation operation is attempted on the lattice."""
    pass

class LatticeUtilsException(Exception):
    """Raised when a general error occurs within the low-level lattice utility functions."""
    pass


##
## Cython-code exceptions

class InnerLoopException(Exception):
    """Raised when an error is detected inside the performance-critical (Cython) inner loop."""
    pass

##
## MISC EXCEPTIONS
##

class ResidueAugmentException(Exception):
    """Raised when an individual residue cannot be added/augmented as expected."""
    pass

class UnfinishedCodeException(Exception):
    """Raised when execution reaches a code path that has not yet been implemented."""
    pass
