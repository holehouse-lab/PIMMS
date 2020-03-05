## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................


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
