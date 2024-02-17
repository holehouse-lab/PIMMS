## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2023
## ...........................................................................



##
## simulation
##
## This file represents the main() function for PIMMS. This is where the magic happens!
## 
##

import random
import sys
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta 
from copy import deepcopy

from .lattice import Lattice
from .chain import Chain
from .acceptance import AcceptanceCalculator
from .moves import MoveObject
from .latticeExceptions import SimulationEnergyException
from .latticeExceptions import SimulationException
from .moveEvent import MoveEvent
from .chainTSMMC import TSMMC
from . import pimmslogger
from . import data_structures       

from . import inner_loops

from . import energy
from . import acceptance
from . import initialized_systems
from . import analysis_IO
from . import analysis_general
from . import pdb_utils
from . import mega_crank # needed to set random seed...
from . import restart

# utility modules 
from . import lattice_utils
from . import lattice_analysis_utils
from . import longrange_utils
from . import cluster_utils
from . import IO_utils
from . import nonequilibrium_utils
from . import system_utils

from . import CONFIG

# if we want to check memory usage we can
# eidt this and then call hp.heap() as
# needed in the code
CHECK_MEMORY = False
if CHECK_MEMORY:
    from guppy import hpy
    hp = hpy()


class Simulation:
    """
    The Simulation object: this is the master object from which all simulation-related 
    events happen

    """

    #-----------------------------------------------------------------
    #
    def __init__(self, keyword_lookup):        
        """
        The simulation constructor is the the main object that consyructs the physical system. 

        
        keyword_lookup is a dictionary with a controlled vocabulary that will read in all of the
        information needed to run a simulation. This dictionary should be generated from the 
        keyfile parser, and expects to have the following key-value pairs
        
        CHAIN : list of lists, where each sublist is a tuple where element 0 is the number of chains
        and element 1 is the sequence of the chain.

        TEMPERATURE : int, temperature for simulation
        
        SEED : int, random seed for simulation
        
        PARAMETER_FILE : str, string that defines parameter file location
        
        NON_INTERACTING : bool, flag that sets if this is a non_interacting run or not
        
        ANGLES_OFF : bool, flag that sets if the angle energies should be considered or not
        
        ENERGY_CHECK : int, frequency with which global/local energy are compared
        
        PRINT_FREQ : int, frequency with which status is printed to STDOUT
        
        EN_FREQ : int, frequency with which energy is written to ENERGY.dat
        
        XTC_FREQ : int, frequency with which data are written to the output XTC file
        
        N_STEPS : int, total number of simulation steps
        
        EQUILIBRATION 
        
        ANALYSIS_FREQ
        
        CRANKSHAFT_SUBSTEPS
        
        CRANKSHAFT_MODE
        
        QUENCH_RUN
        
        QUENCH_START
        
        QUENCH_END
        
        QUENCH_FREQ
        
        QUENCH_STEPSIZE
        
        __TSMMC_USED
        
        TSMMC_INTERPOLATION_MODE
        
        TSMMC_JUMP_TEMP
        
        TSMMC_STEP_MULTIPLIER
        
        TSMMC_NUMBER_OF_POINTS
        
        TSMMC_FIXED_OFFSET
        
        HARDWALL
        
        Optional keywords are


        Parameters
        ----------------
        keyword_lookup : dict


        
         


        """

        
        ## SET UP THE LOGGER
        IO_utils.status_message('SETTING UP THE SIMULATION', 'major')

        ## CORE CONSISTENCY TESTS
        system_utils.check_dtype_consistency()

        system_utils.check_beads_to_grid_mapping(keyword_lookup['CHAIN'])
        

        ## SET LOCAL VARIABLES
        # set local variables for use in initialization
        chains                  = keyword_lookup['CHAIN'] 
        temperature             = keyword_lookup['TEMPERATURE']
        random_seed             = keyword_lookup['SEED']
        parameter_file          = keyword_lookup['PARAMETER_FILE']
        non_interacting         = keyword_lookup['NON_INTERACTING'] 
        angles_off              = keyword_lookup['ANGLES_OFF'] 
                
        # set simulation object variables to be used throughout the simulation    
        self.compare_energyfreq   = keyword_lookup['ENERGY_CHECK']
        self.printfreq            = keyword_lookup['PRINT_FREQ']
        self.reduced_printing     = keyword_lookup['REDUCED_PRINTING']
        self.enfreq               = keyword_lookup['EN_FREQ'] 
        self.xtcfreq              = keyword_lookup['XTC_FREQ']
        self.n_steps              = keyword_lookup['N_STEPS']
        self.equilibration        = keyword_lookup['EQUILIBRATION']
        self.anafreq              = keyword_lookup['ANALYSIS_FREQ']
        self.CS_substeps          = keyword_lookup['CRANKSHAFT_SUBSTEPS'] 
        self.CS_mode              = keyword_lookup['CRANKSHAFT_MODE'] 
        self.LATTICE_TO_ANGSTROMS = keyword_lookup['LATTICE_TO_ANGSTROMS']
        self.autocenter           = keyword_lookup['AUTOCENTER']

        # set quench keywords
        self.QUENCH_RUN         = keyword_lookup['QUENCH_RUN'] 
        self.QUENCH_START       = keyword_lookup['QUENCH_START'] 
        self.QUENCH_END         = keyword_lookup['QUENCH_END'] 
        self.QUENCH_FREQ        = keyword_lookup['QUENCH_FREQ'] 
        self.QUENCH_STEPSIZE    = keyword_lookup['QUENCH_STEPSIZE'] 
                    
        # set for updates to the TSMMC mode 
        self.TSMMC_USED                = keyword_lookup['__TSMMC_USED']        
        self.TSMMC_INTERPOLATION_MODE  = keyword_lookup['TSMMC_INTERPOLATION_MODE']
        self.TSMMC_JUMP_TEMP           = keyword_lookup['TSMMC_JUMP_TEMP']
        self.TSMMC_STEP_MULTIPLIER     = keyword_lookup['TSMMC_STEP_MULTIPLIER']
        self.TSMMC_NUMBER_OF_POINTS    = keyword_lookup['TSMMC_NUMBER_OF_POINTS']
        self.TSMMC_FIXED_OFFSET        = keyword_lookup['TSMMC_FIXED_OFFSET']
        self.production_hardwall       = keyword_lookup['HARDWALL']

        # set whether saving at end. 
        self.SAVE_AT_END       = keyword_lookup['SAVE_AT_END']
        # set whether saving equilibration steps
        self.SAVE_EQ           = keyword_lookup['SAVE_EQ']

        # set None as the mdtraj obj for now. This will be updated every time the coordinates of the system are saved
        # if we use set self.SAVE_AT_END=True. 
        self.master_traj_obj = None

        # analysis settings
        self.analysis_settings  = data_structures.AnalysisSettings(cluster_threshold=keyword_lookup['ANA_CLUSTER_THRESHOLD'])

        # set flags for auxillary chain MC moves (e.g. TSMMC). Set to False to start with
        self.auxillary_chain = False

        # set box size - this is a bit fiddly...
        if keyword_lookup['RESIZED_EQUILIBRATION']:
            
            dimensions = keyword_lookup['RESIZED_EQUILIBRATION']
            self.resize_eq = True
            self.current_xtc_filename = 'eq_traj.xtc'
            self.current_pdb_filename = 'eq_START.pdb'

            # regardless of what keyfile says, we must run initial compact sims with a hardwall
            # boundary to avoid the scenario in which we're re-sizing a system with chains crossing
            # a PBC
            self.hardwall = True 
                        
        else:
            dimensions = keyword_lookup['DIMENSIONS']
            self.resize_eq = False
            self.current_xtc_filename = 'traj.xtc'
            self.current_pdb_filename = 'START.pdb'
            self.hardwall = self.production_hardwall
            
        self.production_dims = keyword_lookup['DIMENSIONS']

        # set values for 10 and 5 percent of the simulation with over-ride values
        # in case we're running especially short simulations        
        if self.n_steps >= 10: 
            self.ten_percent = round(self.n_steps/10)
        else:
            self.ten_percent = 1

        if self.n_steps >= 20:            
            self.five_percent = round(self.n_steps/20)
        else:
            self.five_percent = 1
            
        self.global_start_time = None
        

        ## --------------------------------------------------------------------
        ## Part 1 - Randomization stuff
        ##        

        IO_utils.status_message("Using random seed   : %i" % (random_seed), 'startup')
        IO_utils.status_message("Using C random seed : %i" % (random_seed % CONFIG.C_RAND_MAX),'startup')
        IO_utils.status_message("System RAND_MAX     : %i" % (CONFIG.C_RAND_MAX),'startup')

        pimmslogger.log_status('Random Seed: %i'% (random_seed))
        pimmslogger.log_status('C random Seed: %i'% (random_seed%CONFIG.C_RAND_MAX))
        pimmslogger.log_status('C RAND_MAX (system): %i'% (CONFIG.C_RAND_MAX))

        random.seed(random_seed)
        np.random.seed(random_seed%CONFIG.C_RAND_MAX)
        mega_crank.seed_C_rand(random_seed%CONFIG.C_RAND_MAX)
            

        ## Part 2 - Build the Markov Chain Monte Carlo Metrpolis Acceptance
        #           object and set the various move probabilities therein
        self.ACC       = AcceptanceCalculator(temperature, keyword_lookup)


        ## Part 3 - Build the chain-mover object
        self.MOVER     = MoveObject()


        ## Part 4 - Build the system Hamiltonian based on the
        #           parameter file, or using an empty Hamiltonian
        #           for a non-interacting (Excluded volume) run. Note that non-interacting
        #           only gets used if set to True. Also note that we provide the equilibrium
        #           temperature which may be used to define the angle interaction energies if
        #           requested.
        self.Hamiltonian  = energy.Hamiltonian(parameter_file, len(dimensions), non_interacting, angles_off, hardwall = self.hardwall, temperature = keyword_lookup['EQUILIBRIUM_TEMPERATURE'], reduced_printing=self.reduced_printing)

        
        ## Part 5 - Build the actual simulation lattice!
        if keyword_lookup['RESTART_FILE']:


            # if  we passed a restart file then construct the lattice object using the restart file directly. Note             
            self.LATTICE   = Lattice(dimensions, chains, self.Hamiltonian, self.LATTICE_TO_ANGSTROMS, restart_object=keyword_lookup['RESTART_FILE'], hardwall=self.hardwall)
            
            # safety to ensure we don't break things when reading a restart file 
            if self.LATTICE.any_chains_straddle_boundary():
                self.hardwall = False
                self.Hamiltonian.set_hardwall(False)
                IO_utils.status_message("Restart-read file incompatible with hardwall simulation -> switching to PBC",'warning')
                pimmslogger.log_status("Restart-read file incompatible with hardwall simulation -> switching to PBC")
        else:
            self.LATTICE   = Lattice(dimensions, chains, self.Hamiltonian, self.LATTICE_TO_ANGSTROMS, hardwall = self.hardwall )


        ## Part 6 - Build the Chain Temperature Switch Metropolis Monte Carlo if 
        #           this is being used (if its not being used don't even try and create the 
        #           TSMMC_coordinator object - this is because if TSMMC moves are not being used we 
        #           don't want to force the user to have sane TSMMC parameters, which creating 
        #           a TSMMC_coordinator object would required
        #           
        #           
        if self.TSMMC_USED:
            self.TSMMC_coordinator = TSMMC(temperature, self.TSMMC_JUMP_TEMP, self.TSMMC_INTERPOLATION_MODE, self.TSMMC_STEP_MULTIPLIER, self.TSMMC_NUMBER_OF_POINTS, self.TSMMC_FIXED_OFFSET)
        else:
            self.TSMMC_coordinator = None

        ## Part 7 - Set all the custom analysis frequencies
        #
        #
        (self.non_default_freq_analysis, self.default_freq_analysis) = self.setup_analysis(keyword_lookup)
    
       
    #-----------------------------------------------------------------
    #       
    def run_simulation(self):
        """
        Run the simulation!

        Parameters:
        

        None

        Returns:

        None

        """


        # get the time everything kicks off...
        self.global_start_time = datetime.now()

        IO_utils.status_message("Simulation started at %s" % (str(self.global_start_time)),'startup')
        if CHECK_MEMORY:
            heap = hp.heap()
            print(heap)

        IO_utils.newline()

        # evaluate the initial energy of the system
        (old_energy, old_energy_local, old_energy_LR, old_energy_SLR, old_energy_angles) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)

        if self.QUENCH_RUN:
            with open('QUENCH.dat', 'w') as fh:
                fh.write('')

        # setup the initial trajectory and pdb files
        IO_utils.status_message("Building initial trajectory and pdb files...",'startup')
        lattice_utils.start_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename=self.current_pdb_filename, xtc_filename=self.current_xtc_filename)

        self.startup_analysis()

        IO_utils.status_message("Evaluating initial energy...",'startup')   
        IO_utils.newline()
        IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')        
        print("   ENERGY COMPARISON")   
        print("     STEP             : %i   " % 0)
        print("     GLOBAL           : %i" % old_energy)        
        print("     SHORT RANGE      : %i" % old_energy_local)
        print("     LONG RANGE       : %i" % old_energy_LR)
        print("     SUPER LONG RANGE : %i" % old_energy_SLR)
        print("     ANGLES           : %i" % old_energy_angles)

        IO_utils.newline()
        IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')
        print("   MEMORY USAGE")
        print(f"     GRID             : {sys.getsizeof(self.LATTICE.grid)/1048576:.1f} MB ")
        print(f"     TYPEGRID         : {sys.getsizeof(self.LATTICE.type_grid)/1048576:.1f} MB")
        
        
        IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')
        IO_utils.newline()

        IO_utils.status_message('STARTING SIMULATION','major')
        IO_utils.status_message('  Start time: %s'%(self.global_start_time), 'vanilla')
        

        # flush means we flush all the premable text to STDOUT - useful for running
        # jobs on clusters 
        sys.stdout.flush()

        ##==================================================================##
        ##                                                                  ##
        ##                   MASTER LOOP BEGINS HERE!                       ##
        ##                                                                  ##
        ##==================================================================##
        i = 0
        chain_selection_override=[]


        while i < self.n_steps:
            i = i + 1
            
            # if we're not using an auxillary chain (i.e. this is what happens
            # 99.9% of the time)
            if not self.auxillary_chain:

                # if we're doing a temperature quench..
                if self.QUENCH_RUN:
                    self.quench_update(i, old_energy)


                # if we're equilibrating in a different size box check what's goin' on there. If equilibration 
                # is done then update the energy as calculated via PBC 
                if self.resize_eq:

                    # note the chain_selection_override here will usually be an empty list UNLESS we get
                    # to the end of an equilibration period and there are chains that straddle the boundary,
                    # in which case we need to force those chains to move, so the chain_selection override
                    # ends up defining which chains are forced to move (i.e. all chains that don't straddle
                    # the boundary and frozen).
                    (chain_selection_override, old_energy) = self.update_dimensions(i, old_energy)
                            
                ## ***************************************************************
                ## Pre-move functionality
                ##
                ## Any analysis, IO or other things is done here...
                ##

                # run simulation I/O (write trajectory, energy, STDOUT, also uses reduced_printing)
                self.simulation_IO(i, old_energy)

                # run any/all analysis                
                self.run_all_analysis(i)

            # this is what happens if we're inside an auxillary chain
            else:

                # decrement the global counter, as auxillary chain moves don't count towards the 
                # global move count
                i=i-1                

                # this is where any/all updates happen to do with the TSMMC. The returned status tuple
                # tells us if the auxilary chain was complete and if the move was accepted or not
                tsmmc_move_status = self.auxillary_chain_update(old_energy)
                
                # IF the TSMMC move is finished!
                if tsmmc_move_status[0]:
                    i=i+1
                    
                    # if move was accepted 
                    if tsmmc_move_status[1]:
                        success=True
                        pass
                        
                    # if move was rejected
                    else:
                        # NOTE reverting back to the pre-move lattice is done in the auxillary_chain_update
                        # function, so we just have to revert the energy back to the pre-move value
                        old_energy = self.TSMMC_coordinator.system_move_original_energy
                        success=False
                        
                    # finally we reset the temperature to the system temperature and 
                    # zero out temporary information held during the TSMMMC move
                    self.ACC = self.TSMMC_coordinator.system_move_finalize(self.ACC)
                    self.auxillary_chain = False
                    self.ACC.auxillary_chain = False

                    # NOTE this has to come after we turn the ACC auxillary chain
                    # flag in in the AcceptanceCalculator object
                    self.ACC.update_move_logs(12, success)            
                    
                    ## Finally, now the move has happened do any analysis or IO necessary on this step
                    # run simulation I/O (write trajectory, energy, STDOUT)
                    self.simulation_IO(i, old_energy)

                    # run any/all analysis                    
                    self.run_all_analysis(i)
                    
                    # finally continue to the next real main-chain move
                    continue
                    
            #*************************************************************
            ## Move time! 
            ##
            ## First we select a random chain to 
            ##
            
            # select a random chain to perturb            
            chain_to_move   = self.LATTICE.get_random_chain(override=chain_selection_override)            
            chainID         = chain_to_move.chainID


            
            # get the currentposition of the chain we're going to move
            chain_length = len(chain_to_move.get_ordered_positions())

            ## MOVE SELECTION -------------------------------------------------------------------
            
            # reset the move accepted flag
            move_accepted = False

            # select a move to make (note the chain length matters - right now for single bead
            # chains we have a specific move set, though this should be changed in the future...)
            selection = self.ACC.move_selector(chain_length)
                        
            # if the chain is fixed skip that bad boy
            if chain_to_move.fixed:
                selection = 0
                success = False


            #
            #for chainID in self.LATTICE.chains:
            #    lattice_utils.check_chain_connectivity(chainID, self.LATTICE.chains[chainID].get_ordered_positions(), self.LATTICE.dimensions)

            # system shake
            if selection == 1:

                ## system_shake moves            
                (new_latticeObject, new_energy, total_proposed, total_accepted) = self.MOVER.system_shake(self.LATTICE, old_energy, self.ACC, self.Hamiltonian, self.CS_substeps, self.CS_mode, self.hardwall)

                ## Finally record moves for post-hoc analysis of movesets
                self.ACC.megastep_update_move_logs(1, total_accepted, total_proposed)

                # update energy
                old_energy = new_energy                                
                
                # skip everything else, all hail the megamove! NOTE that we have induvidual accept/rejects inside the system_shake() so this is still performing
                # Metropolis Monte Carlo ON THE SAME MARKOV CHAIN [important] - the place where the move is accepted/rejected has just moved, but we're evaluating
                # with the same Hamiltonian at the same temperature.
                continue
                                
            # translation
            elif selection == 2:                                
                (move_event, success) = self.MOVER.chain_translate(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)

            # rotation
            elif selection == 3:
                (move_event, success) = self.MOVER.chain_rotate(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                
            # chain pivot
            elif selection == 4:       
                (move_event, success) = self.MOVER.chain_pivot(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                
            # head pivoting
            elif selection == 5:                                
                (move_event, success) = self.MOVER.head_pivot(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                
            # chain slither
            elif selection == 6:
                (move_event, success) = self.MOVER.chain_slither(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                                
            # cluster translate
            elif selection == 7:                
                (move_event, success) = self.MOVER.cluster_translate(chain_to_move, 
                                                                     self.LATTICE, 
                                                                     cluster_move_threshold=None,
                                                                     cluster_size_threshold=self.LATTICE.get_number_of_chains()-1,
                                                                     hardwall=self.hardwall)
                
            # cluster rotation
            elif selection == 8:
                (move_event, success) = self.MOVER.cluster_rotate(chain_to_move, 
                                                                  self.LATTICE, 
                                                                  cluster_move_threshold=None,
                                                                  cluster_size_threshold=self.LATTICE.get_number_of_chains()-1,
                                                                  hardwall=self.hardwall)

            ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ##
            ## Note that in TSMMC and ratchet_pivot moves we create an alternative Markov chain and accept-reject on that chain
            ## before finally accepting/rejecting the final conformation *back* into the true system chain, where all this accept
            ## and rejection occurs inside the MOVER's move function (hence the 'continue' at the end of these moves).
            ##
            ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                
            # chain-based temperature sweep Metropolis Monte Carlo (TSMMC). 
            elif selection == 9:
                                
                
                (new_latticeObject, new_energy, total_moves, success) = self.MOVER.Chain_based_TSMMC(chainID, self.LATTICE, old_energy, self.Hamiltonian, self.TSMMC_coordinator, self.hardwall)

                # Update the lattice object!
                self.LATTICE = new_latticeObject

                old_energy = new_energy                                

                # Rinally record moves for post-hoc analysis of movesets
                self.ACC.update_move_logs(selection, success)

                # Finally update the alternative Markov chain move count - used for
                # for performance 
                self.ACC.alt_Markov_chain_update_move_logs(total_moves)

                continue 

            # multichain-based temperature sweep Metropolis Monte Carlo (TSMMC)
            elif selection == 10:
                (new_latticeObject, new_energy, total_moves, success) = self.MOVER.multichain_based_TSMMC(chainID, self.LATTICE, old_energy, self.Hamiltonian, self.TSMMC_coordinator, self.hardwall)
                

                self.LATTICE = new_latticeObject
                
                old_energy = new_energy                                

                # Record moves for post-hoc analysis of movesets
                self.ACC.update_move_logs(selection, success)

                # Finally update the alternative Markov chain move count - used for
                # for performance 
                self.ACC.alt_Markov_chain_update_move_logs(total_moves)

                continue 

            # ratchet pivot
            elif selection == 11:
                raise SimulationException('Ratchet pivot does not seem to maintain detailed balanace - do not use. This may be removed in later')
                (new_latticeObject, new_energy, total_moves, success) = self.MOVER.ratchet_pivot(chainID, self.LATTICE, old_energy, self.ACC, self.Hamiltonian, self.hardwall)
                
                old_energy = new_energy                                

                # Finally record moves for post-hoc analysis of movesets

                # Update the lattice object!
                self.LATTICE = new_latticeObject

                # Record moves for post-hoc analysis of movesets
                self.ACC.update_move_logs(selection, success)

                # Finally update the alternative Markov chain move count - used for
                # for performance 
                self.ACC.alt_Markov_chain_update_move_logs(total_moves)

                continue 

            # system-wide TSMMC
            elif selection == 12:

                if self.reduced_printing is False:
                    IO_utils.status_message("Performing System TSMMC...",'info')

                # create a backup and activate the auxillary chain flags
                self.TSMMC_coordinator.start_system_TSMMC(self.LATTICE.lattice_backupcopy(), old_energy, self.ACC)
                self.auxillary_chain = True
                self.ACC.auxillary_chain = True
            
                # from this moment we are *in* the system TSMMC
                continue


            elif selection == 13:
                
                # NEED to deepcopy here 
                old_chain_positions = deepcopy(chain_to_move.get_ordered_positions())
                
                # [STEP 1] Firstly we do an intial relaxation - this move updates everything implicitly
                (new_latticeObject, new_energy, total_proposed_part1, total_accepted) = self.MOVER.single_chain_shake(chainID, self.LATTICE, old_energy, self.ACC, self.Hamiltonian, self.CS_substeps, self.CS_mode, self.hardwall)
            
                # [STEP 2] next translate the chain to some random position (chain_to_move has updated positions AFTER the initial relaxation state 
                (move_event, success) = self.MOVER.chain_translate(chain_to_move, self.LATTICE.grid, self.hardwall)

                #success = False
                # if the chain movement worked...
                if success:

                    # calculate the new energy in the new position (after wiggling and then jumping)
                    local_dif = self.single_chain_move(move_event, chainID) 

                    # note that NOW the system (all lattices and the associated chains objects) are in an updated state with respect to
                    # the translation move that just happened
                    
                    # So, we now perform an additional relaxation. This move updates everything implicitly (note we use new_energy and the local_dif
                    # to define the starting energy)
                    
                    (new_latticeObject, new_energy, total_proposed_part2, total_accepted) = self.MOVER.single_chain_shake(chainID, self.LATTICE, new_energy+local_dif, self.ACC, self.Hamiltonian, self.CS_substeps, self.CS_mode, self.hardwall)
                    
                    # and get the chain's new positions :-)
                    new_chain_positions = chain_to_move.get_ordered_positions()

                    # update total proposed moves for book-keeping and delete the total_proposed_part1/2
                    # moves from the namespace (safety & sanity)
                    total_proposed = total_proposed_part1+total_proposed_part2
                    del total_proposed_part1 
                    del total_proposed_part2

                    # finally create the 'new' moveEvent which gets evaluated. Note that the energy change evaluated refects the full shake-jump-shake, and it is
                    # THIS that is accepted or rejected. Every move within this move respects microscopic reversibility, so the full move respects detailed balance                    
                    move_event = MoveEvent(original_positions        = old_chain_positions,
                                           moved_positions           = new_chain_positions,
                                           original_chain_positions  = old_chain_positions,
                                           moved_chain_positions     = new_chain_positions,
                                           moved_indices             = list(range(0,len(new_chain_positions))),
                                           move_type                 = 13)

                
                # if we get here the 'JUMP' part was rejected due to a hardsphere clash on the jump, so we just need to 
                # revert the initial single_chain_shake move
                else:

                    # get the new position of the chain (i.e. after the initial relaxation step)
                    new_chain_positions = chain_to_move.get_ordered_positions()

                    # construct a new MoveEvent
                    move_event = MoveEvent(original_positions        = old_chain_positions,
                                           moved_positions           = new_chain_positions,
                                           original_chain_positions  = old_chain_positions,
                                           moved_chain_positions     = new_chain_positions,
                                           moved_indices             = list(range(0,len(new_chain_positions))),
                                           move_type                 = 13)

                    # ..and reject in the usual manner
                    self.single_chain_revert(move_event, chainID)

                                    
            
            else:
                raise SimulationException(latticeException.message_preprocess('Invalid option passed... [%s]'%str(selection)))
                


            # If the hard-sphere energy allowed the move we then evaluate the change to the system energy
            if success:                

                # ..........................................................................................
                # SINGLE CHAIN MOVES! (1/2/3/4/5/6)
                if selection > 0 and selection < 7:
                                        
                    # determine the change in energy associated with this single chain move
                    local_dif             = self.single_chain_move(move_event, chainID) 

                    # Check if the move is accepted based on the Metropolis-Hasting's criterion
                    if self.ACC.boltzmann_acceptance(old_energy, old_energy + local_dif):

                        # accepted - update old_energy
                        old_energy = old_energy + local_dif
                        
                        # update the flag!
                        move_accepted = True
                    else:
                        # rejected - re-configure the system back to its former glory!
                        self.single_chain_revert(move_event, chainID)
                    
                # ..........................................................................................
                # CLUSTER rotation/translation (7/8)
                elif selection > 6 and selection < 9:
                    
                    local_dif  = self.rigid_cluster_move(move_event.moved_positions, move_event.original_positions)
                    
                    # Check if the move is accepted based on the Metropolis-Hasting's criterion
                    if self.ACC.boltzmann_acceptance(old_energy, old_energy + local_dif):
                        # accepted!
                        old_energy = old_energy + local_dif

                        # update the flag!
                        move_accepted = True

                    else:
                        # rejected!
                        self.rigid_cluster_revert(move_event.moved_positions, move_event.original_positions)


                # if we're doing a jump and relax move
                elif selection == 13:

                    # Check if the move is accepted based on the Metropolis-Hasting's criterion
                    if self.ACC.boltzmann_acceptance(old_energy, new_energy):
                        # accepted - update old_energy
                        old_energy = new_energy

                        # update the flag!
                        move_accepted = True

                        ## update the auxillary chain move information
                        self.ACC.alt_Markov_chain_update_move_logs(total_proposed)

                    else:                        
                        # rejected - re-configure the system back to its former glory!
                        self.single_chain_revert(move_event, chainID)

                    
                    

            # in the case of success being False the move caused a hard-sphere clash and is rejected out
            # of hand
            else:
                pass

            ## Finally record move for post-hoc analysis of movesets
            self.ACC.update_move_logs(selection, move_accepted)


        ###
        ### THE END IS NIGH!
        ### 

        # if we get here we have finished looping over the main simulation loop. Congrats?
            
        # save out the master traj if we are saving at end. Only do if True or we will overwrite the traj file. 
        if self.SAVE_AT_END == True:
            lattice_utils.save_out_sim(self.master_traj_obj, self.current_xtc_filename)

            
        global_end_time = datetime.now()
        IO_utils.newline()            
        IO_utils.status_message("Simulation complete", 'info')

        # extract time and build an easy to read string!
        diff = relativedelta(global_end_time, self.global_start_time)
        total_time_msg = "Simulation time:  %d hours, %d minutes, %d seconds" % (diff.hours, diff.minutes, diff.seconds)

        IO_utils.status_message("Simulation finished at %s" % (str(global_end_time)), 'info')
        IO_utils.status_message(total_time_msg, 'info')
        IO_utils.newline()
        IO_utils.status_message("Performing final analysis output...", 'info')

        self.end_of_simulation_analysis()
    
        ### Always (regardless of interval) save a restart file corresponding to the final state of the simulation. 
        self.ANAFUNCT_save_restart(i)
    
        IO_utils.status_message(".... done!\n\nWe hope the results are all you hoped for!", 'info')
        IO_utils.newline()



    #-----------------------------------------------------------------
    #               
    def auxillary_chain_update(self, old_energy):
        """
        Function that performs all the busywork associated with the system-wide
        TSMMC move. Whether or not this function should be here or a) inside the
        TSMMC object or inside the MOVER object I'm not sure. For now it can live
        here based on the logical that it's making global changes to the actual
        simulation system so should remain associated with the Simulation object
        but this might change in the future.

        Return:

        Two-place boolean tuple

        [Move complete, Move accepted]
        Move accepted is always false if the move has not complete (obviously)

        """

        #print "On move %i" %(self.TSMMC_coordinator.system_move_count)
        
        # if check to see if the temperature-sweep has finished
        if self.TSMMC_coordinator.system_move_complete():

            # check if move was accepted
            if self.TSMMC_coordinator.accept_system_TSMMC(old_energy):

                if self.reduced_printing is False:
                    IO_utils.status_message("System TSMMC: ACCEPTED [dE = %5.5f]" % (old_energy - self.TSMMC_coordinator.system_move_original_energy),'info')
                
                # if we get here the move was accepted!!
                # DO NOT RESET THE LATTICE!                
                success = True     

            else:                

                if self.reduced_printing is False:
                    IO_utils.status_message("System TSMMC: REJECTED [dE = %5.5f]" % (old_energy - self.TSMMC_coordinator.system_move_original_energy),'info')
                # RESET THE LATTICE
                self.LATTICE.lattice_restorefrombackup(self.TSMMC_coordinator.system_move_original_info[0], self.TSMMC_coordinator.system_move_original_info[1], self.TSMMC_coordinator.system_move_original_info[2])
                success = False

            # reset the ACC back to its pre TSMMC move status                                    
            self.auxillary_chain = False        
            return (True, success)

        else:
            # note this only changes the ACC object if the temperature
            # has changed, but updates various local chain parameters
            # for book-keeping
            self.ACC = self.TSMMC_coordinator.check_in_system_TSMMC(self.ACC)                    

            return (False, False)






    #-----------------------------------------------------------------
    #               
    def quench_update(self, i, old_energy):
        """
        Helper function to run quench update if a temperature quench
        run is being performed. Updates all relevant simulation variables
        appropriately. No return value.

        """

        # if the current step is requires a temperature update
        if i % self.QUENCH_FREQ == 0:
                    
            if self.ACC.temperature == self.QUENCH_END:


                IO_utils.status_message(f'Reached target temperature of [{self.ACC.temperature}] - no change', 'info')
                pimmslogger.log_status(f'Target temperature reached on step {i} (Target={self.ACC.temperature})')

                # turn off the quench run flag as we're no longer performing a quench run
                self.QUENCH_RUN = False
            else:
                        
                # update the temperature in an inteligent way
                self.ACC.update_temperature(nonequilibrium_utils.update_temperature_in_quench(self.QUENCH_STEPSIZE, self.QUENCH_START, self.QUENCH_END, self.ACC.temperature, self.reduced_printing))
                        
                # update the TSMMC_coordinator temperature if TSMMC is being used (specifically, the TSMMC_coordinator object needs to know the main Markov Chain temperature so it
                # knows what temperature to return to
                if self.TSMMC_USED:
                    self.TSMMC_coordinator = TSMMC(self.ACC.temperature, self.TSMMC_JUMP_TEMP, self.TSMMC_INTERPOLATION_MODE, self.TSMMC_STEP_MULTIPLIER, self.TSMMC_NUMBER_OF_POINTS, self.TSMMC_FIXED_OFFSET)
                else:
                    self.TSMMC_coordinator = None
                        
                # finally write out to the quench file reporting on the quench event
                analysis_IO.write_quench_file(i, self.ACC.temperature, old_energy)


    #-----------------------------------------------------------------
    #                               
    def simulation_IO(self, i, old_energy):
        """
        Helper function to run simulation IO (O) for various different things. Additional
        output should be added her as an if statement comparing against the appropriate
        keyword frequency. 

        NOTE: All analysis is dealt with seperatly and shouldn't be added here - this is
        for non-analysis IO (i.e. status IO).

        This includes

        1. Printing status of the simulation

        2. Writing out trajectory information

        3. Writing out energy information

        4. Performing global energy comparison

        Parameters
        -----------------
        i : int
            Current step that the simulation is on

        old_energy : int
            

        """

        ##
        # define a local function which we can then call in different
        # places. This just avoids us re-writing the same code in multiple
        # places
        def local_status():
            IO_utils.status_message("Step %i of %i [%2.3f %%] (Energy = %i)" %(i, self.n_steps, 100*(float(i)/float(self.n_steps)),old_energy), 'update')


        # flag that avoids this function re-printing the same information multiple times
        statusPrinted = False
        

        # first up we're going to do some performance analysis. This happens every 1/20th of the simulation AND 20
        # steps in so we get an initial estimate on how long this is gonna take quite quickly. 
        if i % self.five_percent == 0 or i == 20:
            analysis_general.evaluate_performance(i, self.global_start_time, self.n_steps, self.equilibration)
        
        # print status if we're at a printfreq interval of steps
        if i % self.printfreq == 0:
            
            if statusPrinted is False:
                local_status()            
                statusPrinted = True 

        # save coordinates
        if i % self.xtcfreq == 0:

            if statusPrinted is False:
                
                # if we're not doing reduced printing print!
                if self.reduced_printing is False:

                    # remark about saving coordinates only if we're saving coordinates
                    if self.SAVE_EQ==False:
                        if i > self.equilibration:
                            local_status()
                            statusPrinted = True
                            IO_utils.status_message("Saving coordinates...")
                    else:
                        local_status()
                        statusPrinted = True
                        IO_utils.status_message("Saving coordinates...")
                        

                # if we are doing reduced printing 
                else:

                    # is this 1/0th of the way through the simulation?
                    if i % self.ten_percent == 0:
                        local_status()
                        statusPrinted = True
                        # remark about saving coordinates only if we're saving coordinates
                        if self.SAVE_EQ==False:
                            if i > self.equilibration:
                                IO_utils.status_message("Saving coordinates [reduced printing mode]...")
                        else:
                            IO_utils.status_message("Saving coordinates [reduced printing mode]...")

            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=- # 
            # -=-=-=-=-=- SAVING THE traj.xtc FILE -=-=-=-=-=- #
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=- #
            # commenting out for now... This is the old saving version. Keeping in code for now in case we want to use it for debugging down the line. 
            #lattice_utils.append_to_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter)
            
            # this is the new version that doesn't write a frame.pdb file and then make the xtc from that. 
            # this is used by default in case we have memory issues with the approach of just updating the 
            # mdtraj Trajectory object.
            if self.SAVE_AT_END==False:
                # if we are saving eq, save regardless of eq step.  
                if self.SAVE_EQ==True:
                    # if we are not saving at the end, we need to append the new coordinates to the xtc file. 
                    lattice_utils.append_to_xtc_file_non_redundant(self.LATTICE, self.LATTICE.lattice_to_angstroms, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter) 
                else:
                    # check if we are passed the eq.
                    if i > self.equilibration:
                        lattice_utils.append_to_xtc_file_non_redundant(self.LATTICE, self.LATTICE.lattice_to_angstroms, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter) 
            else:
                # if we are saving the xtc file at the end, we need to update the master traj object. 
                # however, we don't want to do this if we aren't saving at the end because it will slow things
                # down and take up memory. 
                if self.SAVE_EQ==True:
                    self.master_traj_obj = lattice_utils.update_master_traj(self.LATTICE, 
                                            self.LATTICE.lattice_to_angstroms, self.master_traj_obj,
                                            self.current_pdb_filename, autocenter= self.autocenter)

                else:
                    if i > self.equilibration:
                        self.master_traj_obj = lattice_utils.update_master_traj(self.LATTICE, 
                                            self.LATTICE.lattice_to_angstroms, self.master_traj_obj,
                                            self.current_pdb_filename, autocenter= self.autocenter)



        # save energy
        if i % self.enfreq == 0:
            analysis_IO.write_energy(i, old_energy)
                
        # check global energy
        if i % self.compare_energyfreq == 0:

            # if we haven't printed the status message yet, print it
            if statusPrinted is False:
                local_status()
                statusPrinted = True

            IO_utils.newline()
            IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')

            # recalculate the energy using the full Hamiltonian from scratch
            (recalculated_energy, new_energy_local, new_energy_long_range, new_SLR_energy, new_energy_angles) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)

            # calculate the difference between our locally-tracked energy and the fully recalcalculated energy (these should be the same)
            current_diff = recalculated_energy - old_energy

            # print out the energy comparison and all current energy info
            print("   ENERGY COMPARISON")   
            print("     STEP             : %i   " % i)
            print("     GLOBAL           : %i" % recalculated_energy)
            print("     CURRENT          : %i" % old_energy)
            print("     DIFFERENCE       : %i" % current_diff)     
            print("     SHORT RANGE      : %i" % new_energy_local)
            print("     LONG RANGE       : %i" % new_energy_long_range)
            print("     SUPER LONG RANGE : %i" % new_SLR_energy)
            print("     ANGLES           : %i" % new_energy_angles)
            IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')

            # uncomment for memory info...
            if CHECK_MEMORY:
                heap = hp.heap()
                print(heap)

            IO_utils.newline()

            # if the energy comparison is off, raise an exception and write out the current configuration
            if not current_diff == 0:                    
                lattice_utils.start_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename='CONFIG_AT_ENERGY_FAIL.pdb', xtc_filename='CONFIG_AT_ENERGY_FAIL.xtc')
                print('Writing out abort trajectory to CONFIG_AT_ENERGY_FAIL.pdb/xtc') 
                raise SimulationEnergyException("ERROR: Something is wrong because energy comparisons were off...")

        # flush output
        sys.stdout.flush()
                

    #-----------------------------------------------------------------
    #               
    def single_chain_move(self, move_event, chainID):
        """
        Function which implements optimized energy calculations for moves which move a single chain

        Input Arguments:

        > move_event
        MoveEvent object containing all the move details necessary
    
        > chainID
        The ID of the single chain being removed (EACH CHAIN has a unique ID, starting at 1 and going up).
                
        """
        
        moved_positions       = move_event.moved_positions
        original_positions    = move_event.original_positions
        moved_chain_positions = move_event.moved_chain_positions
        moved_indices         = move_event.moved_indices        
        angle_indices         = move_event.get_angle_indice(self.LATTICE.chains[chainID].seq_len)
        dimensions            = self.LATTICE.dimensions
        num_moved             = len(moved_positions)
        binary_LR_array       = self.LATTICE.chains[chainID].get_LR_binary_array()[moved_indices]

            
        # We want to evaluate the energy with the chain in both positions - right now
        # 1) self.LATTICE.grid has the chain in it's new positoin
        # 2) The chain object in self.LATTICE.chains has it's position in the OLD position
        # 3) The self.LATTICE.type_grid  has the chain in its old position too
        
        # So we revert self.LATTICE.grid back to the original position to get the energy (note the
        # type grid was never changed so doesn't have to be 'reverted' back)
        lattice_utils.delete_chain_by_position(moved_positions, self.LATTICE.grid, chainID)
        lattice_utils.place_chain_by_position(original_positions, self.LATTICE.grid, chainID, safe=True)

        # extact out all the short-range and long-range inter-residue pairs
        (old_region_SR_pairs, old_region_LR_pairs, old_region_SLR_pairs) = lattice_utils.build_all_envelope_pairs(original_positions, binary_LR_array, self.LATTICE.type_grid, dimensions)

        ### get the energy of the area around the chain we're moving
        old_lattice_old_region     = self.Hamiltonian.evaluate_local_energy(self.LATTICE, old_region_SR_pairs)
        old_lattice_old_region_LR  = self.Hamiltonian.evaluate_local_energy_LR(self.LATTICE, old_region_LR_pairs)
        old_lattice_old_region_SLR = self.Hamiltonian.evaluate_local_energy_SLR(self.LATTICE, old_region_SLR_pairs)

        # old_restraint_energy = self.Hamiltonian.evaluate_restraints(self.LATTICE, chainID, moved_indices)
        
        ## evaluate the angle energy (NOTE that most of the moves only perturb a SMALL number of angles so the
        # number of iterations in the list comprehension is typically < 5 (i.e. super fast). This implementatoin
        # is ~20x faster than the old implementation, making the angle energy basically free :-)
        temporary_positions = self.LATTICE.chains[chainID].get_ordered_positions()
        intcode_seq         = self.LATTICE.chains[chainID].get_intcode_sequence()
        old_angle_energy = self.Hamiltonian.evaluate_angle_energy([temporary_positions[i] for i in angle_indices], [intcode_seq[i] for i in angle_indices], dimensions)


        #print "old_lattice_old_region   : %3.2F" % old_lattice_old_region
        #print ' ""        ""  LR        : %3.2F' % old_lattice_old_region_LR
        #print ' ""        ""  SLR       : %3.2F' % old_lattice_old_region_SLR
                
        ### delete the regions of the chains we're going to move from the grid                
        lattice_utils.delete_chain_by_position(original_positions, self.LATTICE.grid, chainID)                
        self.LATTICE.delete_chain_from_type_grid(chainID, original_positions, moved_indices, safe=True)

        # get the SR interactions of the new positions with the empty array (already have the SR interactions for the original position)
        # note we ensure that we get the SR interactions by defining the binary_LR array as all zero (np.zeroes(num_moved)), and we then
        # return the 0-th index to only return the SR interactions
        new_region_SR_pairs = lattice_utils.build_all_envelope_pairs(moved_positions, np.zeros(num_moved, dtype=int), self.LATTICE.type_grid, dimensions)[0] 

        # NOTE that *RIGHT NOW* we haven't deleted the chain from the self.LATTICE.chains list, however
        # the chain is overwritten when we insert a new chain (and the chains list is NOT used in the
        # energy calculations) so this is OK!

        # evaluate the energy of the old space after we've yanked the old chain out (we have to do this to capture the solvent
        # interaction changes at the two sites)
        empty_lattice_old_region      = self.Hamiltonian.evaluate_local_energy(self.LATTICE, old_region_SR_pairs) 
        empty_lattice_old_region_LR   = 0 # LR interactions must be zero
        empty_lattice_old_region_SLR  = 0 # LR interactions must be zero
               
        # evaluate the energy of the space the chain is going to fill
        empty_lattice_new_region      = self.Hamiltonian.evaluate_local_energy(self.LATTICE, new_region_SR_pairs)
        empty_lattice_new_region_LR   = 0 # LR interactions must be zero
        empty_lattice_new_region_SLR  = 0 # LR interactions must be zero

        #print "empty_lattice_old_region   : %3.2F" % empty_lattice_old_region
        #print ' ""           new region   : %3.2F' % empty_lattice_new_region


        ### insert chain into new position
        self.LATTICE.chains[chainID].set_ordered_positions(moved_chain_positions)
        lattice_utils.place_chain_by_position(moved_positions, self.LATTICE.grid, chainID, safe=True)                
        self.LATTICE.insert_chain_into_type_grid(chainID, moved_positions, moved_indices, safe=True)

        # get the LONG-RANGE interactions for the new position
        (new_region_LR_pairs, new_region_SLR_pairs) = longrange_utils.build_LR_envelope_pairs(moved_positions, binary_LR_array, self.LATTICE.type_grid, dimensions)
        
                                                            
        # get the energy of the local area around the chain we've just inserted
        new_lattice_new_region       = self.Hamiltonian.evaluate_local_energy(self.LATTICE, new_region_SR_pairs)
        new_lattice_new_region_LR    = self.Hamiltonian.evaluate_local_energy_LR(self.LATTICE, new_region_LR_pairs)
        new_lattice_new_region_SLR   = self.Hamiltonian.evaluate_local_energy_SLR(self.LATTICE, new_region_SLR_pairs)

        # new_restraint_energy = self.Hamiltonian.evaluate_restraints(self.LATTICE, chainID, moved_indices)
        
        #print "new_lattice_new_region   : %3.2F" % new_lattice_new_region
        #print ' ""        ""  LR        : %3.2F' % new_lattice_new_region_LR
        #print ' ""        "" SLR        : %3.2F' % new_lattice_new_region_SLR

        # and calculate angle changes for 
        temporary_positions = self.LATTICE.chains[chainID].get_ordered_positions()
        intcode_seq         = self.LATTICE.chains[chainID].get_intcode_sequence()
        new_angle_energy = self.Hamiltonian.evaluate_angle_energy([temporary_positions[i] for i in angle_indices], [intcode_seq[i] for i in angle_indices], dimensions)

        # Calculate the energy difference                    
        local_dif     =  ((new_lattice_new_region + new_lattice_new_region_LR + new_lattice_new_region_SLR) + (empty_lattice_old_region + empty_lattice_old_region_LR + empty_lattice_old_region_SLR)) - ((old_lattice_old_region + old_lattice_old_region_LR + old_lattice_old_region_SLR) + (empty_lattice_new_region + empty_lattice_new_region_LR + empty_lattice_new_region_SLR))
                                                                                                                                                                                                                          
        local_dif = local_dif + (new_angle_energy - old_angle_energy)

        
        #print "Short range : %3.2f" % ((new_lattice_new_region + empty_lattice_old_region) - (old_lattice_old_region + empty_lattice_new_region))
        #print "Long range  : %3.2f" % ((new_lattice_new_region_LR + empty_lattice_old_region_LR) - (old_lattice_old_region_LR + empty_lattice_new_region_LR))
        #print "Total range : %3.2f" % local_dif
        #print ""
        return local_dif


    #-----------------------------------------------------------------
    #       
    def single_chain_revert(self, move_event, chainID):
        """
        Function which reverts the system back to the original state after a singe chain
        move is rejected based on the energy difference
    
        """
        moved_positions            = move_event.moved_positions
        original_positions         = move_event.original_positions
        moved_chain_positions      = move_event.moved_chain_positions
        original_chain_positions   = move_event.original_chain_positions
        moved_indices              = move_event.moved_indices
                

        # revert the lattice to it's pre-move state 
        lattice_utils.delete_chain_by_position(moved_chain_positions, self.LATTICE.grid, chainID)
        lattice_utils.place_chain_by_position(original_chain_positions, self.LATTICE.grid, chainID, safe=True)
        
        self.LATTICE.chains[chainID].set_ordered_positions(original_chain_positions)
        
        # update the type_grid variable BACK 
        self.LATTICE.update_type_grid(chainID, moved_positions, original_positions, moved_indices, safe=True)


    #-----------------------------------------------------------------
    #       
    def rigid_cluster_move(self, new_chain_positions, old_chain_positions):
        """
        Function which implements optimized energy calculations for rigid body cluster moves. 

        ********************************************************************************************************

        NOTE that unlike the single chain moves we actually determine the set of moved pairs inside this function. There's a reason for
        this! So, when making a rigid cluster move we first determine the positions of all the chains in the cluster we're moving. This provides us with
        a useful set of prior information because we KNOW that in that clusters' original position the ONLY short range interactions we care about are between lattice
        sites occupied by the cluster and lattice sites occupied by the solvent. We know this because any sites between a cluster-component and a NON solvent 
        site would be an intra-cluster pair, and given rigid cluster movements cannot be changing intra-cluster interactions the change in energy associated with 
        intracluster sites must be zero. However, the key computational cost here is evaluating how the long range interactions contribute to the cost
        of moving the cluster.

        Because we know all the relative interactions WITHIN the cluster must be held fixed (both short and long-range) the only interactions we care about
        are between the cluster interface. If we have NO long range interactions we can actually just perform the move without worrying about the energy because 
        - by definition - we cannot be moving a cluster into direct contact with another solute molecule so the cluster-system interface is purely solute-solvent
        before and after - i.e. no change in energy. If we do have long range interactions (as would be usual) we have to compute their influence. 

        1) Determine all the long-range interactions goin' on
        2) Just TRANSLATE/ROTATE those positions to get the interfacial pairs in the clusters' new position

        However, to do this we need the cluster back in its original position - hence why we have to move the lattice BACK to its original position before
        we determine the interfacial residues
        
        ********************************************************************************************************

        Input Arguments:
        
        new_chain_positions [dictionary of chainIDs, where each new_chain_positions[chainID] represents a list of the new positions associated with that chain]
        Full set of new positions for the set of chains which make up the cluster

        old_chain_positions [dictionary of chainIDs, where each old_chain_positions[chainID] represents a list of the original positions associated with that chain]
        Full set of old positions for the set of chains which make up the cluster

        new_region_pairs
        The complete, redundant set of new interface residues associated with the cluster movement in the new position where the cluster is moving to
        
        old_region_pairs
        The complete, redundant set of old interface residues associated with the cluster's original position

        """
        
        dimensions = self.LATTICE.dimensions

        if not list(new_chain_positions.keys()) == list(old_chain_positions.keys()):
            raise Exception("I don't even care this should NOT HAPPEN")

        # ------------------------------------------------------------------------------------
        # shortcut incase we have no LR interactions then this move is automatically performed as it's
        # energy neutral, so no need to compute the <DELTA> energy as by definition it must be moving
        # a cluster from a fully solvated environment to a fully solvated environment
        if len(self.Hamiltonian.LR_residue_names) == 0:
                    
            # update the type grid (note we have to do this in two independent steps)
            for chainID in old_chain_positions:
                self.LATTICE.delete_chain_from_type_grid(chainID, old_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)

            # update the chain positions on the type grid and in the chains list
            for chainID in new_chain_positions:
                self.LATTICE.chains[chainID].set_ordered_positions(new_chain_positions[chainID])
                self.LATTICE.insert_chain_into_type_grid(chainID, new_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)

            # energy neutral move
            return 0.0
        # ------------------------------------------------------------------------------------

        ## If there are LR interactions...
        # We want to evaluate the energy with the chain in both positions - right now
        # 1) self.LATTICE.grid has the chain in it's new positoin
        # 2) The chain object in self.LATTICE.chains has it's position in the OLD position
        # 3) The self.LATTICE.type_grid  has the chain in its old position too
        
        # So we revert self.LATTICE.grid back to the original position to get the energy 
        # associated with the cluster in its original position
        # (note the type grid was never changed so doesn't have to be 'reverted' back)
        chainIDs = list(new_chain_positions.keys())

        old_region_LR_pairs = {}
        new_region_LR_pairs = {}

        ## xoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxo
        ##
        ## STAGE 1 - original configuration
        ##        
        # for each chain (note this *has* to be a two step process of removal followed by
        # insertion - otherwise you may get clashes if you use a single for-loop
               
        # delete and re-insert (need to do in two steps)
        for chainID in chainIDs:
            lattice_utils.delete_chain_by_position(new_chain_positions[chainID], self.LATTICE.grid, chainID)

        for chainID in chainIDs:
            lattice_utils.place_chain_by_position(old_chain_positions[chainID], self.LATTICE.grid, chainID, safe=True)
            
        # now **all** chains have been moved back to their original possitione we can calculate the reduced 
        # LR pairs (where both residues in the pair participate in LR interactions) - note many (most) of 
        # these interactions will be *within* the cluster, but for simplicity of code we just let this happen - 
        # the computational cost of finding which pairs have one member outside the cluster is greater than
        # just doing all the pairs..
        old_region_LR_pairs = {}
        non_redundant_LR_pairs_old_full = []

        old_region_SLR_pairs = {}
        non_redundant_SLR_pairs_old_full = []

        for chainID in chainIDs:
            
            # get the positions of LR interaction residues in the chain            
            #LR_original_positions_tmp   = longrange_utils.get_LR_positions(old_chain_positions[chainID], range(0,len(old_chain_positions[chainID])), self.LATTICE.chains[chainID].LR_IDX)
            
            # get all the LR
            #old_region_LR_pairs[chainID] = longrange_utils.build_LR_envelope_pairs(LR_original_positions_tmp, self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions)            
            (old_region_LR_pairs[chainID], old_region_SLR_pairs[chainID])  = longrange_utils.build_LR_envelope_pairs(old_chain_positions[chainID], self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions)

            # get all the pairs of LR interactions between for the chainID
            non_redundant_LR_pairs_old_full.extend(old_region_LR_pairs[chainID])
            non_redundant_SLR_pairs_old_full.extend(old_region_SLR_pairs[chainID])

        # perform energy evaluation
        ENERGY_old_lattice_old_region = self.Hamiltonian.evaluate_local_energy_LR(self.LATTICE, np.array(non_redundant_LR_pairs_old_full)) + self.Hamiltonian.evaluate_local_energy_SLR(self.LATTICE, np.array(non_redundant_SLR_pairs_old_full))
        
        ## xoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxo
        ##
        ## STAGE 2 - cluster chains removed (so a reduced lattice) - recall we do now have to delete from main 
        ## and type lattice

        # first delete all the original chains
        for chainID in chainIDs:
            lattice_utils.delete_chain_by_position(old_chain_positions[chainID], self.LATTICE.grid, chainID)                
            self.LATTICE.delete_chain_from_type_grid(chainID, old_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)

        # now we set the LR energy associated with the empty lattice to 0 - because it must be by definition. long-range interactoins ONLY occur when
        # there is a PAIR of beads which participate in long-range interactions. In the empty lattice, at the sites where the chain(s) of interest were 
        # and will be the site *must* be empty, so there cannot be a pair of LR interacting beads.
        ENERGY_empty_lattice_old_region = 0
        ENERGY_empty_lattice_new_region = 0
        
        ## xoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxo
        ##
        ## STAGE 3 - add in all the cluster chains back into their new positions

        # now insert all the chains into their new positions in the lattice
        for chainID in chainIDs:

            # update the new chain position
            self.LATTICE.chains[chainID].set_ordered_positions(new_chain_positions[chainID])

            # add chains into main and type grids
            lattice_utils.place_chain_by_position(new_chain_positions[chainID], self.LATTICE.grid, chainID, safe=True)                
            self.LATTICE.insert_chain_into_type_grid(chainID, new_chain_positions[chainID], list(range(0,len(new_chain_positions[chainID]))), safe=True)

        new_region_LR_pairs = {}
        non_redundant_LR_pairs_new_full = []

        new_region_SLR_pairs = {}
        non_redundant_SLR_pairs_new_full = []

        for chainID in chainIDs:
            
            # get the positions of LR interaction residues in the chain            
            #LR_original_positions_tmp   = longrange_utils.get_LR_positions(new_chain_positions[chainID], range(0,len(new_chain_positions[chainID])), self.LATTICE.chains[chainID].LR_IDX)

            # get the dumb full LR pairs 
            #new_region_LR_pairs[chainID] = longrange_utils.build_LR_envelope_pairs(LR_original_positions_tmp, self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions)
            (new_region_LR_pairs[chainID], new_region_SLR_pairs[chainID]) = longrange_utils.build_LR_envelope_pairs(new_chain_positions[chainID], self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions)

            # get all the pairs of LR interactions between for the chainID
            non_redundant_LR_pairs_new_full.extend(new_region_LR_pairs[chainID])
            non_redundant_SLR_pairs_new_full.extend(new_region_SLR_pairs[chainID])

        # finally  perform energy evaluation
        ENERGY_new_lattice_new_region = self.Hamiltonian.evaluate_local_energy_LR(self.LATTICE, np.array(non_redundant_LR_pairs_new_full)) + self.Hamiltonian.evaluate_local_energy_SLR(self.LATTICE, np.array(non_redundant_SLR_pairs_new_full))

        ## STAGE 4
        # Having now calculated all the relevant LR interactions we sum them up and use them to evaluate the energy of the move
        local_dif     = (ENERGY_new_lattice_new_region + ENERGY_empty_lattice_old_region) - (ENERGY_old_lattice_old_region + ENERGY_empty_lattice_new_region)    
                
        return local_dif


    #-----------------------------------------------------------------
    #       
    def rigid_cluster_revert(self, new_chain_positions, old_chain_positions):
        """
        Function which reverts the system back to the original state after a rigid cluster
        move is rejected
        
        """
        # revert the lattice to it's pre-move state  (delete everything)        
        for chainID in new_chain_positions:            
            lattice_utils.delete_chain_by_position(new_chain_positions[chainID], self.LATTICE.grid, chainID)
            self.LATTICE.delete_chain_from_type_grid(chainID, new_chain_positions[chainID], list(range(0,len(new_chain_positions[chainID]))), safe=True)

        # now re-insert everything
        for chainID in old_chain_positions:
            lattice_utils.place_chain_by_position(old_chain_positions[chainID], self.LATTICE.grid, chainID, safe=True)                    
            self.LATTICE.insert_chain_into_type_grid(chainID, old_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)
            self.LATTICE.chains[chainID].set_ordered_positions(old_chain_positions[chainID])

            
    #-----------------------------------------------------------------
    #           
    def update_dimensions(self, step, old_energy):
        """
        Function that updates the dimensions of the lattice.  This is done by


        """
        
        # if this step is the end of equilibration do all the fun jazz, else we simply return
        # an empty list and the old energy
        if step == self.equilibration:
                
            # for each chain is in a non-periodic configuration
            # if no - keep selecting a cluster and move move move until yess
            # also print warning - probably means equilibration is too short!
            # repeat


            # assess if each chain on the lattice straddles the boundary or not 
            offending_chains=[]
            for chainID in self.LATTICE.chains:
                if self.LATTICE.chains[chainID].does_chain_stradle_pbc_boundary():
                    offending_chains.append(chainID)
                    
            # if we found one or more offending chains, print a warning, increment the number of steps and  
            if len(offending_chains) > 0:

                IO_utils.status_message("%i chains are still crossing the periodic boundary despite this being a hardwall simulation...." % (len(offending_chains)), 'warning')
                IO_utils.status_message("Dynamically extending number of steps and equilibration",'info')
                pimmslogger.log_warning("%i chain(s) are still crossing the periodic boundary despite this being a hardwall simulation.\nOffending chains are:[%s]"% ( len(offending_chains), offending_chains))
                
                self.n_steps = self.n_steps+100
                self.equilibration = self.equilibration+100
                
                return (offending_chains, old_energy)

            # if we get here all the chains are valid, inasmuch as they are all within a non-periodic space, allowing us to change 
            # the lattice dimensions without fear of breaking everything! 
            pimmslogger.log_status("Resizing lattice dimensions from [%s] to [%s]" %(self.LATTICE.dimensions, self.production_dims))

            # create a new restart object, instantiate it with the current lattice, and then update the restart object's positions
            # to center 
            R = restart.RestartObject()
            R.build_from_lattice(self.LATTICE, self.production_hardwall)
            R.update_lattice_dimensions(self.production_dims)
                                             
            # use this restart object to construct a new lattice
            # the [] is the 'empty' chains list which would normally be passed from the keyfile, but we can disregard here,
            # but is a required parameter (ugly, but it's OK...)
            new_lattice = Lattice(self.production_dims, [], self.Hamiltonian, self.LATTICE_TO_ANGSTROMS, restart_object=R)
            # once that's done then 

            # finally assign this new lattice to the simulation object 
            self.LATTICE = new_lattice
            
            # turn off the resize flag and update the output file names
            # see if we need to save the output when 'save at end' is set to True. . 
            if self.SAVE_AT_END==True:
                # if saving EQ == True
                if self.SAVE_EQ == True:
                    # save the output
                    lattice_utils.save_out_sim(self.master_traj_obj, self.current_xtc_filename)
                    # reset master_traj_obj to None. 
                    self.master_traj_obj=None
            
            # set self.resize_eq to false, reset the namds of pdb and xtc files. 
            self.resize_eq = False
            self.current_pdb_filename='START.pdb'
            self.current_xtc_filename='traj.xtc'
            
            # initialize the xtc/pdb output files with these new names
            lattice_utils.start_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename=self.current_pdb_filename, xtc_filename=self.current_xtc_filename)

            # clean up if possible!
            import gc                
            gc.collect()

            # If we want to switch to PBC based on the keyfile HARDWALL variable, then do so.
            # Regardless, we now recalculate the new energy and return this
            if self.production_hardwall is False:
                self.Hamiltonian.set_hardwall(False)
                self.hardwall = self.production_hardwall
            
            (energy, _, _, _, _) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)

            # return an empty list which sets the chain_selection_override to empty. Note that the calling function is aware of success
            # because self.resize_eq has been switched from True to False
            return([], energy)
        else:
            return ([], old_energy)

        
    

    ######################################################################################
    ##                                                                                  ##
    ##                             ANALYSIS ROUTINES                                    ##
    ##                                                                                  ##
    ######################################################################################
    #
    # The functions below are general setup for running sytem-wide analysis. Note that
    # actual analysis logic should NOT be included here, and should be implemented in 
    # either the Chains class or in the analysis_general.py.
    #
    # ANAFUNCT functions are the functions called by run_all_analysis(), which calls
    # analysis functions at different frequencies depending on how often the analysis
    # is to be performed as defined by the keyfile.
    #
    # These functions must take a single argument (the step number) - they don't have
    # to use it but it will always be passed.
    #
    # Many of the functions write data to disk, while others just update internal running 
    # totals.
    #
    #

    #-----------------------------------------------------------------
    #   
    def startup_analysis(self):
        """
        Function for including all the analysis activity which should be run
        BEFORE the simulation starts. This includes wiping files we're going
        to progressively append to during the simulation, but one could imagine
        in the future additional code might be included here for functional
        purposes.

        """
        
        # wipe any existing files (or create them if they don't exist)
        IO_utils.wipe_file(CONFIG.OUTNAME_ENERGY)
        
        ## All writen by the cluster analysis routine
        # cluster analysis stuff
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTERS)
        IO_utils.wipe_file(CONFIG.OUTNAME_NUM_CLUSTERS)
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTER_RG)
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTER_ASPH)
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTER_VOL)
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTER_AREA)
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTER_DENSITY)
        IO_utils.wipe_file(CONFIG.OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE)

        # long-range clusters
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTERS)
        IO_utils.wipe_file(CONFIG.OUTNAME_NUM_LR_CLUSTERS)
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTER_RG)
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTER_ASPH)
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTER_VOL)
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTER_AREA)
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTER_DENSITY)
        IO_utils.wipe_file(CONFIG.OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE)

        # note both of these are written out by the 
        # polymeric properties analysis
        # 
        IO_utils.wipe_file(CONFIG.OUTNAME_RG)
        IO_utils.wipe_file(CONFIG.OUTNAME_ASPH)
        

        IO_utils.wipe_file(CONFIG.OUTNAME_ACCEPTANCE)
        IO_utils.wipe_file(CONFIG.OUTNAME_MOVES)
        IO_utils.wipe_file(CONFIG.OUTNAME_PERFORMANCE, header="Step\tE or P\tSteps-per-second\tElapsed time (hh:mm:ss)\tRemaining time (hh:mm:ss)\n")
        IO_utils.wipe_file(CONFIG.OUTNAME_TOTAL_MOVES)

        IO_utils.wipe_file(CONFIG.OUTNAME_E2E)
        IO_utils.wipe_file(CONFIG.OUTNAME_R2R)

        #IO_utils.wipe_file(CONFIG.OUTNAME_INTER_INTRA)            

        
        # if we have a multicomponent system initialize cluster heterogenity
        # output files
        if len(self.LATTICE.chainTypeList) > 1:
            for chainType in self.LATTICE.chainTypeList:
                IO_utils.wipe_file("CHAIN_%i_"%(chainType) + CONFIG.OUTNAME_CLUSTERS)
                IO_utils.wipe_file("CHAIN_%i_"%(chainType) + CONFIG.OUTNAME_LR_CLUSTERS)
                #IO_utils.wipe_file("CHAIN_%i_"%(chainType) + CONFIG.OUTNAME_INTER_INTRA)
                #IO_utils.wipe_file("CHAIN_%i_"%(chainType) + CONFIG.OUTNAME_MIXING)


                
    #-----------------------------------------------------------------
    #           
    def run_all_analysis(self, step):
        """
        Master analysis function - cycles over each type of analysis to
        assess if that analysis should be performed this step (or not), and
        launching the analysis if it should be done. 

        Updates various analysis state information (e.g. internal scaling
        distances associated with each Chain etc.) but does not change any
        lattice positions or anything like that

        """

        # do not perform analysis if we're still in equilibration
        if step < self.equilibration:
            return

        # for any analysis routines we've defined as occuring at a non-
        # default frequency (recall that the keys in self.[non_]default_freq_analysis
        # are actually the function signatures that are functions of the Simulation
        # class and expect a single parameter to be passed (step)
        for analysis_function in self.non_default_freq_analysis:                        
            if step % self.non_default_freq_analysis[analysis_function] == 0:
                analysis_function(step)

        # for all general analysis we haven't defined
        if step % self.anafreq == 0:
            for analysis_function in self.default_freq_analysis:
                analysis_function(step)



    #-----------------------------------------------------------------
    #           
    def setup_analysis(self, keyword_lookup):
        """
        This function constructs two dictionaries. Each key-value pair in the dictionary is a function-frequency 
        pair, where the function is an analysis routine of the format FXC(step) and the frequency is the frequency
        with which that analysis is performed.

        The two lists correspond to analysis which occurs with the same frequency as the general analysis and then
        the analysis which occurs at a frequency *different* to the general analysis.

        This is a bit of work at the start, but allows us to run bespoke, custom-frequency analysis in a very
        simple way during the simulation.

        keyword_lookup provides all the info needed.

        """
        

        non_default_freq_analysis = {}
        default_freq_analysis = {}

        # set the analysis names here
        all_ana_keywords = ['ANA_POL','ANA_INTSCAL', 'ANA_DISTMAP', 'ANA_ACCEPTANCE', 'ANA_CLUSTER', 'ANA_INTER_RESIDUE', 'ANA_END_TO_END', 'ANA_CUSTOM', 'RESTART_FREQ']

        # define the functions and initialze any closures needed
        analysis_keywords = {}
        analysis_keywords['ANA_POL']             = self.ANAFUNCT_polymeric_properties
        analysis_keywords['ANA_INTSCAL']         = self.ANAFUNCT_internal_scaling
        analysis_keywords['ANA_DISTMAP']         = self.ANAFUNCT_distance_map
        analysis_keywords['ANA_ACCEPTANCE']      = self.ANAFUNCT_acceptance
        analysis_keywords['ANA_CLUSTER']         = self.ANAFUNCT_cluster_analysis
        analysis_keywords['ANA_INTER_RESIDUE']   = self.build_R2R_distance_distribution_analysis(keyword_lookup['ANA_RESIDUE_PAIRS'])
        analysis_keywords['ANA_END_TO_END']      = self.ANAFUNCT_end_to_end
        analysis_keywords['RESTART_FREQ']        = self.ANAFUNCT_save_restart
        
    
        # if a side-loading module was provided
        if keyword_lookup['ANALYSIS_MODULE']:

            # define a closure that adds the LATTICE object to the
            # function call and then calls the custom analysis function
            # with the step and the self.LATTICE object passed as variables
            def fx(step):
                custom_analysis = keyword_lookup['ANALYSIS_MODULE']
                return custom_analysis(step, self.LATTICE)

            analysis_keywords['ANA_CUSTOM']          = fx

        else:
            analysis_keywords['ANA_CUSTOM']          = self.ANAFUNCT_custom_stubb
            
        
        # ------------------------------------------------------->>

        # quick check to ensure all our ducks are in a row...
        if not len(all_ana_keywords) == len(analysis_keywords):
            raise SimulationException('Bug in the the analysis setup routines. This was triggered by a failsafe check and indicates a software bug')

        for AKW in all_ana_keywords:
            if AKW not in analysis_keywords:
                raise SimulationException('Bug in the the analysis setup routines. This was triggered by a failsafe check and indicates a software bug')

        # ------------------------------------------------------->>

    
        # get the default analysis frequency
        anafreq = keyword_lookup['ANALYSIS_FREQ']

        # Having set up all that we now cycle through the analysis types as defined by the all_ana_keywords
        # list. This means that the default_freq_analysis and non_default_freq_analysis dictionaries have 
        # a key-value pairing where the _key_ is the actual function signature and the _value_ is the frequency
        # with which that analysis is done 

        for AKW in all_ana_keywords:
            if keyword_lookup[AKW] == anafreq:
                default_freq_analysis[analysis_keywords[AKW]] = anafreq
            else:
                non_default_freq_analysis[analysis_keywords[AKW]] = keyword_lookup[AKW]

        return (non_default_freq_analysis, default_freq_analysis)
                
                            

        
    #-----------------------------------------------------------------
    #       
    def end_of_simulation_analysis(self):
        """
        Final analysis routines run at the end of the simulation. For all analysis
        where a final average value makes sense this is going to be where the code
        to calculate and save that output is written. 

        """

        ### ()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
        ###
        ### SECTION: INTERNAL SCALING
        ###
        # right now we write out the chain average and chain STD on internal scaling for every 
        # chain probaly want to allow specific analysis groups but all in good time!
        # 
        
        # if the system contains only one type of chain...
        if len(self.LATTICE.chainTypeList) == 1:

            # get the internal scaling and distance map results for each individual chain 
            all_IS         = []
            all_IS_squared = []
            all_nu         = []
            all_R0         = []            
            all_dMap       = []

            for chain in self.LATTICE.chains:
                all_IS.append(self.LATTICE.chains[chain].analysis_get_cumulative_internal_scaling())
                all_IS_squared.append(self.LATTICE.chains[chain].analysis_get_internal_scaling_squared())
                                
                scaling_info = self.LATTICE.chains[chain].analysis_fit_scaling_exponent()                

                all_nu.append(scaling_info[0])
                all_R0.append(scaling_info[1])
                                
                all_dMap.append(self.LATTICE.chains[chain].analysis_get_cumulative_distance_map())

            # calculate the mean internal scaling and write to disk
            mean_IS         = np.array(all_IS).mean(0)        
            mean_IS_squared = np.array(all_IS_squared).mean(0)       
            
            analysis_IO.write_internal_scaling(mean_IS, mean_IS_squared)
            analysis_IO.write_scaling_information(all_nu, all_R0)

            # calculate the mean distance map and write to disk
            mean_dMap = np.asarray(all_dMap)
            analysis_IO.write_distance_map(mean_dMap.mean(axis=0))

        # if the system contains two or more different types of chains
        else:

            all_IS         = {}
            all_IS_squared = {}
            all_nu         = {}
            all_R0         = {}            
            all_dMap       = {}
            for chainTypeID in self.LATTICE.chainTypeList:
                all_IS[chainTypeID]         = []
                all_IS_squared[chainTypeID] = []
                all_nu[chainTypeID]         = []
                all_R0[chainTypeID]         = []
                all_dMap[chainTypeID]       = []
                

            # for each chain add the internal scaling for that chain to a list of IS for each specific chain
            # type
            for chain in self.LATTICE.chains:
                all_IS[self.LATTICE.chains[chain].chainType].append(self.LATTICE.chains[chain].analysis_get_cumulative_internal_scaling())
                all_IS_squared[self.LATTICE.chains[chain].chainType].append(self.LATTICE.chains[chain].analysis_get_internal_scaling_squared())

                # do scaling fitting..
                scaling_info = self.LATTICE.chains[chain].analysis_fit_scaling_exponent()                
                all_nu[self.LATTICE.chains[chain].chainType].append(scaling_info[0])
                all_R0[self.LATTICE.chains[chain].chainType].append(scaling_info[1])

                all_dMap[self.LATTICE.chains[chain].chainType].append(self.LATTICE.chains[chain].analysis_get_cumulative_distance_map())
                
                
            for chainTypeID in self.LATTICE.chainTypeList:
                mean_IS = np.array(all_IS[chainTypeID]).mean(0)
                mean_IS_squared = np.array(all_IS_squared[chainTypeID]).mean(0)        
                analysis_IO.write_internal_scaling(mean_IS, mean_IS_squared, prefix='CHAIN_%i_'%chainTypeID)

                analysis_IO.write_scaling_information(all_nu[chainTypeID], all_R0[chainTypeID], prefix='CHAIN_%i_'%chainTypeID)

                mean_dMap = np.asarray(all_dMap[chainTypeID])
                analysis_IO.write_distance_map(mean_dMap.mean(axis=0), prefix='CHAIN_%i_'%chainTypeID)
            
                
            
        ### END OF INTERNAL SCALING
        ### ()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()


    #-----------------------------------------------------------------
    #       
    def build_R2R_distance_distribution_analysis(self, R2R_info):
        """
        This function returns a function which is initialized by the variables pass
        in by the R2R_info - in essence generating a closure.

        Basically, if you're not familiar with functional programming, this creates
        a new function where the $R2R_info variable inside the ANAFUNCT_R2R_distance
        function is set by the build_R2R_distance_distribution_analysis function.

        This function (ANAFUNCTION_R2R_distance) is then returned, and next time
        its called the R2R_info variable IN THE FUNCTION BEING CALLED is already
        initialzed.

        """
        
        def ANAFUNCT_R2R_distance(step):
                    
            # just skip if no pairs defined...
            if len(R2R_info) == 0:
                pass
                                                                                                            
            all_data = []
            for pair in R2R_info:

                pair_data = []
                for chainID in sorted(self.LATTICE.chains.keys()):                
                    pair_data.append(self.LATTICE.chains[chainID].analysis_get_residue_residue_distance(pair[0], pair[1]))
                
                all_data.append(pair_data)
            
            # finally write the analysis to file
            analysis_IO.write_residue_residue_distance(step, R2R_info, all_data)


        # return the closure function for use
        return ANAFUNCT_R2R_distance
                
                                
    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_internal_scaling(self, step):
        """
        Run internal scaling analysis. Updates running counters associated
        with each chain
        """

        for chainID in self.LATTICE.chains:
            
            # note this updates both normal and squared internal scaling info
            self.LATTICE.chains[chainID].analysis_update_internal_scaling()
            


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_distance_map(self, step):
        """
        Run distance map analysis. Updates running counters associated 
        with each chain.
        """

        for chainID in self.LATTICE.chains:
            self.LATTICE.chains[chainID].analysis_update_distance_map()


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_cluster_analysis(self, step):
        """
        Run cluster analysis

        Writes data out on each call (I/O heavy)
        
        """

        # get clusters list - note this is really computationally expensive
        # so we try and only do this once and then perform any/all cluster analysis
        # subsequent to this!
        (clusters) = lattice_analysis_utils.get_cluster_distribution(self.LATTICE.grid, self.LATTICE.chains)        
        (LR_clusters) = lattice_analysis_utils.get_LR_cluster_distribution(self.LATTICE)        

        big_cluster_idx = []        
        for c_idx in range(0,len(clusters)):            
            if len(clusters[c_idx]) > self.analysis_settings.cluster_threshold:
                big_cluster_idx.append(c_idx)
            else:
                break
        

        big_LR_cluster_idx = []
        for c_idx in range(0,len(LR_clusters)):            
            if len(LR_clusters[c_idx]) > self.analysis_settings.cluster_threshold:
                big_LR_cluster_idx.append(c_idx)
            else:
                break

        # for each cluster extract the set of induvidual positions to get a list of positions
        cluster_positions    = lattice_analysis_utils.extract_positions_from_clusters(clusters, self.LATTICE.chains)
        LR_cluster_positions = lattice_analysis_utils.extract_positions_from_clusters(LR_clusters, self.LATTICE.chains)

        # for each cluster correct the cluster's positions such that each cluster lies in a single periodic image as best can be achieved. Note that when
        # we don't correct for this the cluster analysis ends up being confusing...
        corrected_cluster_positions    = lattice_analysis_utils.correct_cluster_positions_to_single_image(cluster_positions, self.LATTICE.dimensions)
        corrected_LR_cluster_positions = lattice_analysis_utils.correct_LR_cluster_positions_to_single_image(LR_cluster_positions, self.LATTICE.dimensions)

        ## subselect size-thresholded clusters for polymer/gross property/radial distribution analysis. The clusters are sorted by size, so we know that
        # once we find one cluster below the the threshold we've found all the big clusters, hence the 'break' statements
        big_clusters = [corrected_cluster_positions[i] for i in big_cluster_idx] 
        big_clusters_LR = [corrected_LR_cluster_positions[i] for i in big_LR_cluster_idx] 

        # for each set of positions get the polymeric properties associated with each whole cluster using the single image convention corrected values
        cluster_polymeric_properties_list      = lattice_analysis_utils.extract_cluster_polymeric_properties(big_clusters)
        LR_cluster_polymeric_properties_list   = lattice_analysis_utils.extract_cluster_polymeric_properties(big_clusters_LR)

        # for each cluster calculate the volume, surface area and density (requires corrected cluster positions)
        cluster_size_properties     = lattice_analysis_utils.compute_cluster_gross_properties(big_clusters)
        LR_cluster_size_properties  = lattice_analysis_utils.compute_cluster_gross_properties(big_clusters_LR)

        # for each cluster calculate the radial density profile IF the cluster contains more than 27 beads (3x3x3). We should probably make this number a keyfile
        # value
        
        cluster_radial_density     = lattice_analysis_utils.compute_cluster_radial_density_profile(big_clusters, self.LATTICE.dimensions, minimum_cluster_size_in_beads = CONFIG.RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD)
        LR_cluster_radial_density  = lattice_analysis_utils.compute_cluster_radial_density_profile(big_clusters_LR, self.LATTICE.dimensions, minimum_cluster_size_in_beads = CONFIG.RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD)

        # We'll leave the following in as a sanity check

        # remove soon - > for debugging
        """
        count=0
        for idx in range(0, len(corrected_cluster_positions)):
            pdb_utils.write_positions_to_file(corrected_cluster_positions[idx], 'clusters/%i_cluster_%i_CORR.pdb'%(step,count))
            pdb_utils.write_positions_to_file(cluster_positions[idx], 'clusters/%i_cluster_%i_UNCORR.pdb'%(step,count))
            count=count+1
        """            
        
        # write cluster list
        analysis_IO.write_clusters(step, clusters, self.LATTICE.chainIDtoType)
        analysis_IO.write_LR_clusters(step, LR_clusters, self.LATTICE.chainIDtoType)

        # write cluster size/shape analysis
        analysis_IO.write_cluster_properties(step, cluster_polymeric_properties_list, cluster_size_properties, cluster_radial_density)
        analysis_IO.write_LR_cluster_properties(step, cluster_polymeric_properties_list, LR_cluster_size_properties, LR_cluster_radial_density)


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_polymeric_properties(self, step):
        """
        Run radius of gyration analysis

        Writes data out on each call (I/O heavy)
        
        """

        RG_list = []
        asph_list = []
        
        for chainID in sorted(self.LATTICE.chains.keys()):            
            tmp  = self.LATTICE.chains[chainID].analysis_get_polymeric_properties()

            RG_list.append(tmp[0])
            asph_list.append(tmp[1])

        analysis_IO.write_radius_of_gyration(step, RG_list)
        analysis_IO.write_asphericity(step, asph_list)


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_end_to_end(self, step):
        """
        Run radius of gyration analysis

        Writes data out on each call (I/O heavy)
        
        """

        e2e_list = []
        
        for chainID in sorted(self.LATTICE.chains.keys()):
            e2e_list.append(self.LATTICE.chains[chainID].analysis_get_end_to_end_distance())

        analysis_IO.write_end_to_end(step, e2e_list)


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_acceptance(self, step):
        """
        Run acceptance criterion analysis
        
        Writes data out on each call (I/O heavy)

        """
        analysis_IO.write_acceptance_statistics(step, self.ACC)

        
        
    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_custom_stubb(self, step):
        """
        Stubb that is used if no custom analysis module is passed        

        """
        pass




    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_save_restart(self, step):
        """
        Function that outputs current system state to a file that can be used to initialze the system

        """
        IO_utils.status_message("Writing restart file on step %i..." %(step),'info')
        R = restart.RestartObject()

        # build using lattice, and also pass the hardwall status of the current simulation
        R.build_from_lattice(self.LATTICE, self.hardwall)

        # evaluate the total energy and provide this as well
        (energy, _, _, _, _) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)        
        R.set_energy(energy)

        # output restart file to disk
        R.write_to_file()

        
        


