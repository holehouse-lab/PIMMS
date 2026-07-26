## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
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

def evaluate_performance(step, start_time, total_steps, equilibration, acceptanceObject):
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

    acceptanceObject : AcceptanceCalculator
        Object tracking move statistics. Used here via
        ``total_attempted_moves()`` to compute the overall MC move throughput
        (every accept/reject across all sub-loops, not just the outer master
        loop).

    Returns
    -------
    None
        No return variable, but the performance information is written to
        the appropriate file.

    """

    # get current time

    now = datetime.now()
    time_elapsed = now - start_time

    # NOTE: use total_seconds(), not timedelta.seconds, because the latter only
    # returns the sub-day component (0-86399) and would wrap after 24 hours.
    passed_seconds = int(time_elapsed.total_seconds())
    hours = passed_seconds // 3600
    minutes = (passed_seconds % 3600) // 60
    seconds = (passed_seconds % 3600) % 60

    time_elapsed_string = f'{hours:02d}:{minutes:02d}:{seconds:02d}'


    # get steps per second (guard against a zero elapsed time on the very first
    # call, which would otherwise raise ZeroDivisionError)
    seconds_elapsed = (now - start_time).total_seconds()
    if seconds_elapsed <= 0:
        steps_per_second = max(float(step), 1.0)
    else:
        steps_per_second = step / seconds_elapsed

    # OVERALL move rate: every individual accept/reject across ALL sub-loops
    # (megamove substeps + TSMMC excursion substeps), not just the outer master
    # loop. This is the true MC throughput - a single master-loop step can be a
    # megamove worth of thousands of sub-moves.
    total_attempted = acceptanceObject.total_attempted_moves()
    if seconds_elapsed <= 0:
        overall_moves_per_second = float(total_attempted)
    else:
        overall_moves_per_second = total_attempted / seconds_elapsed


    # calculate anticipated time remaining assuming a constant step rate. On step 0 no
    # step has been taken yet, so steps_per_second is exactly 0 and there is nothing to
    # extrapolate from - report 0 rather than dividing by zero. Any real rate is used as
    # is, so a genuinely slow run still gets an honest (long) estimate.
    steps_remaining = total_steps - step
    if steps_per_second > 0:
        remaining_seconds = int(np.ceil(steps_remaining / steps_per_second))
    else:
        remaining_seconds = 0


    hours = remaining_seconds // 3600
    minutes = (remaining_seconds % 3600) // 60
    seconds = remaining_seconds % 60

    time_remaining_string = f'{hours:02d}:{minutes:02d}:{seconds:02d}'

    if step <= equilibration:
        eq_string = 'E'
    else:
        eq_string = 'P'
        
    analysis_IO.write_performance(step, eq_string, steps_per_second, overall_moves_per_second, time_elapsed_string, time_remaining_string)

    percentage_steps_left = (steps_remaining / total_steps) * 100

    # also log estimated time remaining + both throughput rates to the logging system
    pimmslogger.log_status(f"Estimated remaining time: {time_remaining_string} (hh:mm:ss) | {percentage_steps_left:2.0f}% left "
                           f"| {steps_per_second:.1f} loop-steps/s | {overall_moves_per_second:,.0f} MC-moves/s")


