## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2021
## ...........................................................................


import random
import numpy as np

from . import CONFIG
from .latticeExceptions import MoveException

class SystemTSMMC:

    def __init__(self, target_temperature, jump_temp, interp_mode, step_multiplier, number_points, fixed_offset):
        """
        Class which allows the pre-computation of temperature-switch MMC moves, basically
        means we pre-calculated the switching schedule here once and then re-call it on
        each move rather than recalculating the same thing millions of times for no good
        reason.

        Some pre-conditions

        target_temperature : The main temperature the system is running at

        jump_temp          : Default temperature to jump to upon TSMMC moves

        interp_mode        : Interpolation mode through which temperature shift happens (right
                             now the only mode is 'linear'

        step_multiplier    : Multiplier for steps per temperature 

        number_of_points   : Number of points between the max and target temperature that
                             defines the temperature

        fixed_offset       : If a fixed offset is set then this offset is used to define the
                             jump temperature. This is specifically useful when a quench is 
                             being run, so that as the temperature changes the jump temperature
                             changes is a linearly proportional manner

        """

        self.mode = interp_mode
        self.steps_per_quench_multiplier = step_multiplier
        self.target_temperature = target_temperature
        self.inv_target_temperature = CONFIG.INVTEMP_FACTOR/(target_temperature)
        self.fixed_offset = fixed_offset
            
        
        if interp_mode == 'LINEAR':
            
            # get the schedule rounded to 2 decimalplaces going from the jump temperature to the target temperature
        
            # If we're using a fixed value to jump to
            if not fixed_offset:
                dT = jump_temp - target_temperature        
                step = dT/float(number_points)       

            # else if we're using a fixed offset from the target temperature
            else:
                dT = (target_temperature+fixed_offset) - target_temperature        
                step = dT/float(number_points)            
                jump_temp=target_temperature+fixed_offset
            
            self.true_temp_schedule = np.around(np.hstack((np.arange(target_temperature+step, jump_temp+(step), step), np.repeat(jump_temp, CONFIG.TOP_TEMP), np.flipud(np.arange(target_temperature+step, jump_temp+(step), step)))),5)

            self.inv_temperature_schedule = CONFIG.INVTEMP_FACTOR/self.true_temp_schedule
            

                

    def accept_TSMMC(self, new_energy, old_energy, inv_temp, prev_inv_temp):        

        A = (-inv_temp)*new_energy
        B = (-prev_inv_temp)*old_energy
        C = (-inv_temp)*old_energy
        D = (-prev_inv_temp)*new_energy
        
        to_exp = (A+B) - (C+D)

        #print "To EXP val: %10.10f " % to_exp

        if to_exp > 0.0:
            return True

        expterm = np.exp(to_exp)
        #print expterm

        if random.random() < expterm:
            #print "ACCEPTED"
            return True
             
        else:
            return False


    # this function might never be called...
    def start_system_TMMMC(self, backup_tuple, original_energy):
        
        self.system_move_count = 0
        self.system_move_temp_idx = 0
        self.original_info = (backup_tuple[0], backup_tuple[1], backup_tuple[2], original_energy)



    def check_in(self, ACC):
        
        # increment the general counters
        self.system_move_count = self.system_move_count+1
        
        # if the counter is mod-0 to the number of steps per temperature
        # then we update the temperature
        if self.system_move_count % self.steps_per_quench_multiplier == 0:

            ACC.update_temperature(self.true_temp_schedule[self.system_move_temp_idx])
            self.system_move_temp_idx = self.system_move_temp_idx+1

        return ACC
