## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................

import datetime

def message_preprocess(msg):
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
    pass

class SimulationEnergyException(Exception):
    pass

class AcceptanceException(Exception):
    pass

class RestartException(Exception):
    pass



##
## ENERGY EXCEPTIONS
##

class EnergyException(Exception):
    pass

class TemperatureException(Exception):
    pass


##
## LATTICE EXCEPTIONS
##
class ParticleException(Exception):
    pass

class LatticeInitializationException(Exception):
    pass

class ChainInitializationException(Exception):
    pass

class ChainInsertionFailure(Exception):
    pass

class ChainDeletionFailure(Exception):
    pass

class ChainAugmentFailure(Exception):
    pass

class ChainConnectivityError(Exception):
    pass

class TypeGridException(Exception):
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
    pass

class ParameterFileException(Exception):
    pass

class PDBException(Exception):
    pass

class IOException(Exception):
    pass


##
## MOVE EXCEPTIONS
## 

class MoveSetException(Exception):
    pass

class MoveException(Exception):
    pass

class CustomInitializationException(Exception):
    pass

##
## ANALYSIS EXCEPTIONS
class AnalysisStructureException(Exception):
    pass

class AnalysisRoutineException(Exception):
    pass


##
## LATTICE UTILS
class RotationException(Exception):
    pass

class DebuggingException(Exception):
    pass

class LatticeUtilsException(Exception):
    pass


##
## Cython-code exceptions

class InnerLoopException(Exception):
    pass

##
## MISC EXCEPTIONS
##

class ResidueAugmentException(Exception):
    pass

class UnfinishedCodeException(Exception):
    pass
