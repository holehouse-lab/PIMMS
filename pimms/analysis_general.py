## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2023
## ...........................................................................

##
## analysis_general
##
## This file contains some general purpose analysis routines which are relevant for
## the entire simulation system rather than a specific chain. Kind of a catch-all 
## for analysis


from . import CONFIG
from . import analysis_IO
from datetime import datetime
import numpy as np
from . import pimmslogger

def evaluate_performance(step, start_time, total_steps, equilibration):
    """
    Function for analysing simulation performance and then writing to an
    appropriate file. This has been significantly revamped in 0.1.36 to 
    provide estimates of time remaining and elapsed time.

    Parameters
    ----------
    step : int
        Current step number

    start_time : datetime
        The time at which the simulation was started. This should 
        come from the Simulation.global_start_time variable.

    total_steps : int
        The total number of steps in the simulation

    equilibration : int
        The number of steps to be used for equilibration. This is used
        to determine whether the simulation is in the equilibration phase
        or the production phase.

    Returns
    -------
    None
        No return variable, but the performance information is written to
        the appropriate file.

    """

    # get current time

    now = datetime.now()
    time_elapsed = now - start_time

    passed_seconds = time_elapsed.seconds
    hours = passed_seconds // 3600
    minutes = (passed_seconds % 3600) // 60
    seconds = (passed_seconds % 3600) % 60

    time_elapsed_string = f'{hours:02d}:{minutes:02d}:{seconds:02d}'


    # get steps per second 
    seconds_elapsed = (now - start_time).total_seconds()    
    steps_per_second = step / seconds_elapsed

    
    # calculate anticipated time remaining assuming a constant step rate
    steps_remaining = total_steps - step       
    remaining_seconds = int(np.ceil(steps_remaining / steps_per_second))


    hours = remaining_seconds // 3600
    minutes = (remaining_seconds % 3600) // 60
    seconds = remaining_seconds % 60

    time_remaining_string = f'{hours:02d}:{minutes:02d}:{seconds:02d}'

    if step <= equilibration:
        eq_string = 'E'
    else:
        eq_string = 'P'
        
    analysis_IO.write_performance(step, eq_string, steps_per_second, time_elapsed_string, time_remaining_string)

    percentage_steps_left = (steps_remaining / total_steps) * 100

    # also log estimated time remaining to the logging system
    pimmslogger.log_status(f"Estimated remaining time: {time_remaining_string} (hh:mm:ss) | {percentage_steps_left:2.0f}% left")


