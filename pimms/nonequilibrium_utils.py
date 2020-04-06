## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2020
## ...........................................................................
# 



from .latticeExceptions import TemperatureException
from . import IO_utils


#-----------------------------------------------------------------
#    
def update_temperature_in_quench(QUENCH_STEPSIZE, QUENCH_START, QUENCH_END, temperature):
    """
    Quenching function which determines what the simulation's new temperature must be...
    


    """

    if temperature < 1:
        raise TemperatureException('Temperature is less than 1 [%i] - suggests something is very wrong...' % (temperature))
    
    # if a cooling run
    if QUENCH_END < QUENCH_START:


        # easy case - step gets us close to the target 
        if (temperature - QUENCH_STEPSIZE) >= QUENCH_END:                            
            IO_utils.status_message("QUENCH: Updating temperature from %0.3f to %0.3f" % (temperature, temperature - QUENCH_STEPSIZE))
            return (temperature - QUENCH_STEPSIZE)

        else:
            IO_utils.status_message("QUENCH: Trying to update the temperature from %i to %i, but this would skip the target temperature [%i]. Setting to target temperature now... " % (temperature, temperature-QUENCH_STEPSIZE, QUENCH_END))
            return QUENCH_END


    # if a heating run..
    else:
        # easy case - step gets us close to the target 
        if (temperature - QUENCH_STEPSIZE) <= QUENCH_END:                            
            IO_utils.status_message("QUENCH: Updating temperature from %0.3f to %0.3f" % (temperature, temperature - QUENCH_STEPSIZE))
            return (temperature - QUENCH_STEPSIZE)


        else:
            print("QUENCH: Trying to update the temperature from %i to %i, but this would skip the target temperature [%i]. Setting to target temperature now... " % (temperature, temperature-QUENCH_STEPSIZE, QUENCH_END))
            return QUENCH_END



