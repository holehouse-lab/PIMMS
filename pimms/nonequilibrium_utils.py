## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................
# 



from .latticeExceptions import TemperatureException
from . import IO_utils


#-----------------------------------------------------------------
#    
def update_temperature_in_quench(QUENCH_STEPSIZE, QUENCH_START, QUENCH_END, temperature, reduced_printing=False):
    """
    Determine the simulation's new temperature for the next step of a quench
    (a slow, monotonic temperature ramp).

    A single quench step of magnitude ``QUENCH_STEPSIZE`` is applied to the
    current temperature in the direction of ``QUENCH_END``. If the step would
    overshoot the target end temperature, the temperature is clamped exactly to
    ``QUENCH_END`` instead. A status message describing the change is printed
    unless ``reduced_printing`` is True.

    Note that whether the run is a cooling or heating run is inferred from the
    relationship between ``QUENCH_START`` and ``QUENCH_END``, and the update is
    always computed as ``temperature - QUENCH_STEPSIZE`` (so for a heating run
    ``QUENCH_STEPSIZE`` is expected to be negative).

    Parameters
    ----------
    QUENCH_STEPSIZE : float
        The magnitude (and sign, for heating runs) of the temperature change
        applied per quench step.

    QUENCH_START : float
        The temperature at which the quench began. Only used (together with
        ``QUENCH_END``) to decide whether this is a cooling or heating run.

    QUENCH_END : float
        The target temperature the quench is ramping towards. The returned
        temperature is clamped to this value if a step would overshoot it.

    temperature : float
        The simulation's current temperature.

    reduced_printing : bool, optional
        If True, suppress the informational status messages describing the
        temperature update. Default is False.

    Returns
    -------
    float
        The simulation's new temperature for the next step (either
        ``temperature - QUENCH_STEPSIZE`` or ``QUENCH_END`` if the step would
        overshoot the target).

    Raises
    ------
    TemperatureException
        If the current temperature is not positive, which indicates the
        simulation has reached a non-physical state.
    """

    # <= 0, not < 1: TEMPERATURE is a float and the acceptance calculator only
    # requires it to be positive, so a legitimate cooling quench can pass through
    # sub-unity temperatures (e.g. 1.5 -> 0.5 in 0.25 steps). Only a non-positive
    # temperature is genuinely unphysical.
    if temperature <= 0:
        raise TemperatureException('Temperature is <= 0 [%s] - suggests something is very wrong...' % (temperature))
    
    # A fractional step size accumulates floating-point error (1.0 -> 0.2 in
    # steps of 0.1 arrives at 0.20000000000000004), so a new temperature within
    # rounding of the target IS the target: snap to it so the "target reached"
    # bookkeeping (an exact equality test) fires on the right rung rather than
    # one rung late.
    proposed = temperature - QUENCH_STEPSIZE
    if abs(proposed - QUENCH_END) <= 1e-9 * max(1.0, abs(QUENCH_END)):
        if reduced_printing is False:
            IO_utils.status_message("QUENCH: Updating temperature from %0.3f to %0.3f (target reached)" % (temperature, QUENCH_END))
        return QUENCH_END

    # if a cooling run
    if QUENCH_END < QUENCH_START:


        # easy case - step gets us close to the target
        if (temperature - QUENCH_STEPSIZE) >= QUENCH_END:
            if reduced_printing is False:
                IO_utils.status_message("QUENCH: Updating temperature from %0.3f to %0.3f" % (temperature, temperature - QUENCH_STEPSIZE))
            return (temperature - QUENCH_STEPSIZE)

        else:
            if reduced_printing is False:
                IO_utils.status_message("QUENCH: Trying to update the temperature from %0.3f to %0.3f, but this would skip the target temperature [%0.3f]. Setting to target temperature now... " % (temperature, temperature-QUENCH_STEPSIZE, QUENCH_END))
            return QUENCH_END


    # if a heating run..
    else:
        # easy case - step gets us close to the target 
        if (temperature - QUENCH_STEPSIZE) <= QUENCH_END:

            if reduced_printing is False:
                IO_utils.status_message("QUENCH: Updating temperature from %0.3f to %0.3f" % (temperature, temperature - QUENCH_STEPSIZE))
            return (temperature - QUENCH_STEPSIZE)


        else:

            if reduced_printing is False:
                IO_utils.status_message("QUENCH: Trying to update the temperature from %0.3f to %0.3f, but this would skip the target temperature [%0.3f]. Setting to target temperature now... " % (temperature, temperature-QUENCH_STEPSIZE, QUENCH_END))
            return QUENCH_END



