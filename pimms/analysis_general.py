## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

##
## analysis_general
##
## This file contains some general purpose analysis routines which are relevant for
## the entire simulation system rather than a specific chain. Kind of a catch-all 
## for analysis


from . import CONFIG
from . import analysis_IO

def evaluate_performance(step, step_interval, dt):
    """
    Function for analysing simulation performance and then writing to an
    appropriate file. Right now we're only evaluating the time per step.

    """

    analysis_IO.write_time_per_step(step, step_interval, 1000*(dt/step_interval))


