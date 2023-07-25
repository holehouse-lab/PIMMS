## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2020
## 
## ...........................................................................


import numpy as np
cimport numpy as cnp
cnp.import_array()

cimport cython 
import random
from libc.math cimport exp


## Define high performance local min and max functions...
##

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

from libc.stdlib cimport rand, srand, RAND_MAX

ctypedef cnp.int64_t NUMPY_INT_TYPE


def seed_C_rand(int seedval):
    srand(seedval)


#-----------------------------------------------------------------
# 

@cython.wraparound(False)
@cython.boundscheck(False)
#def update_position(cnp.ndarray[NUMPY_INT_TYPE, ndim=1] old_position, cnp.ndarray[NUMPY_INT_TYPE, ndim=3] grid, NUMPY_INT_TYPE x_off, NUMPY_INT_TYPE y_off, NUMPY_INT_TYPE z_off, NUMPY_INT_TYPE XDIM, NUMPY_INT_TYPE YDIM, NUMPY_INT_TYPE ZDIM):

def update_position(NUMPY_INT_TYPE[:] old_position, NUMPY_INT_TYPE[:,:,:] grid, NUMPY_INT_TYPE x_off, NUMPY_INT_TYPE y_off, NUMPY_INT_TYPE z_off, NUMPY_INT_TYPE XDIM, NUMPY_INT_TYPE YDIM, NUMPY_INT_TYPE ZDIM):
    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position = np.zeros([3], dtype=int)

    cdef int local_x = pbc_correction(old_position[0] + x_off, XDIM)
    cdef int local_y = pbc_correction(old_position[1] + y_off, YDIM)
    cdef int local_z = pbc_correction(old_position[2] + z_off, ZDIM)

    #if grid[new_position[0], new_position[1], new_position[2]] > 0:
    if grid[local_x, local_y, local_z] > 0:
        # fail
        new_position[0] = -1

    # hard sphere clash
    else:
        # success        
        new_position[0] = local_x
        new_position[1] = local_y
        new_position[2] = local_z

    return (new_position)
    

#-----------------------------------------------------------------
# 
@cython.wraparound(False)
@cython.boundscheck(False)
def mega_crank(NUMPY_INT_TYPE[:,:,:] grid, 
               NUMPY_INT_TYPE[:,:,:] type_grid, 
               NUMPY_INT_TYPE[:,:] idx_to_bead,
               NUMPY_INT_TYPE[:,:] interaction_table, 
               NUMPY_INT_TYPE[:,:] LR_interaction_table, 
               NUMPY_INT_TYPE[:,:] SLR_interaction_table, 
               NUMPY_INT_TYPE[:,:,:,:,:,:,:] angle_lookup,
               long energy,
               float invtemp,g
               int nsteps,
               cnp.ndarray[NUMPY_INT_TYPE, ndim=1] bead_selector,
               int passed_seed,
               int hardwall):



#def mega_crank(cnp.ndarray[NUMPY_INT_TYPE, ndim=3] grid, 
#               cnp.ndarray[NUMPY_INT_TYPE, ndim=3] type_grid, 
#               cnp.ndarray[NUMPY_INT_TYPE, ndim=2] idx_to_bead,
#               cnp.ndarray[long, ndim=2] interaction_table, 
#               cnp.ndarray[long, ndim=2] LR_interaction_table, 
#               cnp.ndarray[long, ndim=2] SLR_interaction_table, 
#               cnp.ndarray[long, ndim=7] angle_lookup,
#               long energy,
#               float invtemp,
#               int nsteps,
#               cnp.ndarray[NUMPY_INT_TYPE, ndim=1] bead_selector,
#               int passed_seed,
#               int hardwall):


    """
    Some explanation is in order re: what the input variables here below.


    THE MOST IMPORT COMMENT: For some inexplicable reason, if we randomly select a bead and randomly select integers for
    x, y, and z pertubation, we get a correlation between the bead index and the value associated with the final dimension
    being perturbed. This pertains specifically to 

    grid - the main lattice grid

    type_grid - mirrors the lattice grid but the integers represent residue types not chain IDs

    idx_to_bead - this contains ALL the information we need to perturb the lattice. This a x by 6 matrix, where 
                        a is the TOTAL number of beads in the system. Each  index position reveals a vector of length 
                        5 that contains the following information:

                        0 - bead_flag (0 = single bead, 1 = N-terminal bead, 2 = fully central bead (OOXOO), 3 = C-terminal bead
                                      4 = central bead in a OXO configuration, 5 - N-terminal bead +1 from start)
                        1 - LR binary flag
                        2 - intcode value
                        3 - skip angles (1 = true, 0 = false)
                        4 - chainID
                        5 - X position
                        6 - Y position
                        7 - Z position (optional - depends on if we're in 3D or not)

    interaction_table - lookup table for the short range interactions of 2 beads
    
    LR_interaction_table - lookup table for the long range interactions of 2 beads

    SLR_interaction_table - lookup table for the super long range interactions of 2 beads

    energy - current energy

    invtemp - current inverse temperature

    nsteps - number of steps to perform

    bead_selection - pre-alloacted random selection of beads. Avoids a bug in PRNG

    passed_seeed - a random seed generated by the main program, ensures reproducibility

    """
    # set randomseed
    srand(passed_seed)

    cdef unsigned int i, bead_index;
    cdef int accepted_moves;
    cdef int XDIM, YDIM, ZDIM;
    cdef int num_beads
    cdef long delta_energy, delta_angle_energy

    accepted_moves = 0

    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=2] position_triptic = np.zeros([3, 3], dtype = int)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=2] three_position_holder = np.zeros([3, 3], dtype = int)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=2] two_position_holder = np.zeros([2, 3], dtype = int)    

    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] old_position = np.zeros([3], dtype = int)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] anchor_bead = np.zeros([3], dtype = int)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position;

    XDIM = grid.shape[0]
    YDIM = grid.shape[1]
    ZDIM = grid.shape[2]


    # get the number of beads
    num_beads = len(idx_to_bead)
    
    # angle short-circuit
    for i in range(nsteps):

        # select the bead from the bead_selector, a pre-allocated array of random numbers
        bead_index   = bead_selector[i]

        # get position
        old_position[0] = idx_to_bead[bead_index,5]
        old_position[1] = idx_to_bead[bead_index,6]
        old_position[2] = idx_to_bead[bead_index,7]


        # ------------------------------------------------------------
        # if single bead (beadflag == 0)
        if idx_to_bead[bead_index,0] == 0:
            new_position = single_bead_crank(old_position, grid, XDIM, YDIM, ZDIM)


        # ------------------------------------------------------------
        # if N-terminal bead (beadflag == 1)
        elif idx_to_bead[bead_index,0] == 1:
            anchor_bead[0] = idx_to_bead[bead_index+1,5]
            anchor_bead[1] = idx_to_bead[bead_index+1,6]
            anchor_bead[2] = idx_to_bead[bead_index+1,7]

            new_position = single_bead_crank(anchor_bead, grid, XDIM, YDIM, ZDIM)
            
            if hardwall == 1:
                two_position_holder[0,0] = new_position[0]
                two_position_holder[0,1] = new_position[1]
                two_position_holder[0,2] = new_position[2]

                two_position_holder[1,0] = anchor_bead[0]
                two_position_holder[1,1] = anchor_bead[1]
                two_position_holder[1,2] = anchor_bead[2]

                if do_positions_stradle_pbc_boundary(two_position_holder, 2) == 1:
                    continue

        # ------------------------------------------------------------
        # if C-terminal bead (beadflag == 3)
        elif idx_to_bead[bead_index,0] == 3:

            anchor_bead[0] = idx_to_bead[bead_index-1,5]
            anchor_bead[1] = idx_to_bead[bead_index-1,6]
            anchor_bead[2] = idx_to_bead[bead_index-1,7]

            new_position = single_bead_crank(anchor_bead, grid, XDIM, YDIM, ZDIM)

            if hardwall == 1:
                two_position_holder[0,0] = anchor_bead[0]
                two_position_holder[0,1] = anchor_bead[1]
                two_position_holder[0,2] = anchor_bead[2]

                two_position_holder[1,0] = new_position[0]
                two_position_holder[1,1] = new_position[1]
                two_position_holder[1,2] = new_position[2]

                if do_positions_stradle_pbc_boundary(two_position_holder,2) == 1:
                    continue

        # ------------------------------------------------------------
        # we're somewhere inside the chain, and beadflag is 2,4,5, or 6
        else:
            position_triptic[0,0] = idx_to_bead[bead_index-1,5]
            position_triptic[0,1] = idx_to_bead[bead_index-1,6]
            position_triptic[0,2] = idx_to_bead[bead_index-1,7]

            position_triptic[1,0] = idx_to_bead[bead_index,5]
            position_triptic[1,1] = idx_to_bead[bead_index,6]
            position_triptic[1,2] = idx_to_bead[bead_index,7]

            position_triptic[2,0] = idx_to_bead[bead_index+1,5]
            position_triptic[2,1] = idx_to_bead[bead_index+1,6]
            position_triptic[2,2] = idx_to_bead[bead_index+1,7]

            #position_triptic[1] = idx_to_bead[bead_index][5:8]
            #position_triptic[2] = idx_to_bead[bead_index+1][5:8]
        
            # normal crank move!
            new_position = crank_it(position_triptic, grid, XDIM, YDIM, ZDIM)

            if hardwall == 1:
                three_position_holder[0,0] = position_triptic[0,0]
                three_position_holder[0,1] = position_triptic[0,1]
                three_position_holder[0,2] = position_triptic[0,2]


                three_position_holder[1,0] = new_position[0]
                three_position_holder[1,1] = new_position[1]
                three_position_holder[1,2] = new_position[2]

                three_position_holder[2,0] = position_triptic[2,0]
                three_position_holder[2,1] = position_triptic[2,1]
                three_position_holder[2,2] = position_triptic[2,2]
                
                if do_positions_stradle_pbc_boundary(three_position_holder,3) == 1:
                    continue

            
        # if hardsphere success
        if not new_position[0] < 0:
            
            delta_energy = get_energy_change(grid, type_grid, old_position, new_position, idx_to_bead[bead_index,1],  interaction_table, LR_interaction_table, SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
            delta_angle_energy = get_angle_energy_change(bead_index, idx_to_bead, new_position, angle_lookup)
            
            if accept_or_reject(invtemp, energy, energy+delta_energy+delta_angle_energy) == 1:

                #with open('crankshaft_accepted.txt','a') as fh:
                #    fh.write('%i, %i, %i, %i, %i \n' % (old_position[0] - new_position[0], old_position[1] - new_position[1], old_position[2] - new_position[2], bead_index, chainID))

                # delete from main grid and insert into new position                                            
                grid[old_position[0], old_position[1], old_position[2]] = 0
                
                # inert the bead to its new position
                grid[new_position[0],new_position[1],new_position[2]] = idx_to_bead[bead_index,4]
                
                # set new type grid to old type grid value
                type_grid[new_position[0],new_position[1],new_position[2]] = type_grid[old_position[0], old_position[1], old_position[2]]

                # zero out old position on type grid                
                type_grid[old_position[0], old_position[1], old_position[2]] = 0

                # IMPORTANT!!! This change *MUST* be done last as the old_position is actually a reference
                # to chain_position[bead_index] so changing it BEFORE hand changes the old_position variable and
                # breaks everything
                idx_to_bead[bead_index,5] = new_position[0]
                idx_to_bead[bead_index,6] = new_position[1]
                idx_to_bead[bead_index,7] = new_position[2]

                energy = energy + delta_energy + delta_angle_energy
                accepted_moves = accepted_moves + 1

            else:
                #with open('crankshaft_rejected.txt','a') as fh:
                #    fh.write('%i, %i, %i, %i, %i \n' % (old_position[0] - new_position[0], old_position[1] - new_position[1], old_position[2] - new_position[2], bead_index, chainID))
                pass
        else:
            
            #with open('crankshaft_failed.txt','a') as fh:
            #        fh.write('%i, %i, %i, %i, %i \n' % (old_position[0] - new_position[0], old_position[1] - new_position[1], old_position[2] - new_position[2], bead_index, chainID))
            pass

                             

                
    return (energy, accepted_moves)


#-----------------------------------------------------------------
#                 
cdef int get_random_position(int chain_length):
    """
    Returns a random number between 0 and chain_length,
    where both 0 and chain_length are inclusive

    """
    
    cdef  int r;        
    r = randint(0, chain_length)

    return (r)

def get_random_position_python(int chain_length):
    return get_random_position(chain_length)

def randint_python(start,end):
    return randint(start, end)


#-----------------------------------------------------------------
# 
cdef int randint(int start, int end):
    """
    returns a random integer between start and end inclusive
    of ends.

    This mirrors random.randint's behaviour. NOTE 
    numpy.random.randint and random.randint HAVE DIFFERENT
    BEHAVIOUR. This is daft... We use random.randint here
    (where start and end are both incusive).
    
    e.g. start = 0 and end=5 gives one of 0,1,2,3,4,5

    """

    # this is inelegant but makes everything work out...
    if start == 0:
        end=end+1
    else:
        pass
         
    # ok so this is kind of inelegant too, but the rand()-1 is actually important
    # if we don'te have this we risk the situation where rand() == RAND_MAX
    # which would cause r = (start+end) which if start was 0 end has
    # become end+1 so we get a value of end+1 (i.e. outside the range)
    cdef  int r = start+int((float(rand()-1)/float(RAND_MAX))*(end))


    return r


#-----------------------------------------------------------------
# 
@cython.wraparound(False)
@cython.boundscheck(False)
cdef crank_it (NUMPY_INT_TYPE[:,:] position_triptic, NUMPY_INT_TYPE[:,:,:] grid, int XDIM, int YDIM, int ZDIM):
#cdef crank_it (cnp.ndarray[NUMPY_INT_TYPE, ndim=2] position_triptic, cnp.ndarray[NUMPY_INT_TYPE, ndim=3] grid, int XDIM, int YDIM, int ZDIM):
    """
    Perform crankshaft move!

    Returns a tuple of success and the PBC corrected posistion

    """    
   
    cdef int N_side_x, N_side_y, N_side_z, C_side_x, C_side_y, C_side_z;
    cdef int x_min, x_max, y_min, y_max, z_min, z_max
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position = np.zeros([3], dtype=int)

    
    # extract out x and y positions and perform a series of PBC corrections as needed
    N_side_x = position_triptic[0,0]
    N_side_y = position_triptic[0,1]
    N_side_z = position_triptic[0,2]

    C_side_x = position_triptic[2,0]
    C_side_y = position_triptic[2,1]
    C_side_z = position_triptic[2,2]


    #print "Before 'correction' "
    #print "Before positions (from var) [%i, %i, %i]]" % (N_side_x, N_side_y, N_side_z)
    #print "After positions  (from var) [%i, %i, %i]]" % (C_side_x, C_side_y, C_side_z)


    ## Periodicity fix 1
    # move into a non-periodic space just for now

    # if we're over 2 lattice sites away in any dimension
    # this must be straddling a PBC, so figure out which of 
    # the two dimensions is the biggest and move the smaller
    # into the 'fake' space so they're no longer in PBC
    if abs(N_side_x - C_side_x) > 2:
        if N_side_x > C_side_x:
            C_side_x = C_side_x + XDIM
        else:
            N_side_x = N_side_x + XDIM
        
    if abs(N_side_y - C_side_y) > 2:
        if N_side_y > C_side_y:
            C_side_y = C_side_y + YDIM
        else:
            N_side_y = N_side_y + YDIM

    if abs(N_side_z - C_side_z) > 2:
        if N_side_z > C_side_z:
            C_side_z = C_side_z + ZDIM
        else:
            N_side_z = N_side_z + ZDIM

    ## So we are now no longer in periodic space

    # get smallest/largest possible new x value
    x_min = int_max(N_side_x-1, C_side_x-1) 
    x_max = int_min(N_side_x+1, C_side_x+1) 

    # get smallest/largest possible new y value
    y_min = int_max(N_side_y-1, C_side_y-1) 
    y_max = int_min(N_side_y+1, C_side_y+1) 

    # get smallest/largest possible new z value
    z_min = int_max(N_side_z-1, C_side_z-1) 
    z_max = int_min(N_side_z+1, C_side_z+1) 
    

    # if a random bead is selected in the code, the ORDER in which these are defined introduce a bias
    # such that the system drives high number beads up in the Z dimension but low number beads dow
    
    cdef int local_x = pbc_correction((x_min + randint(1, (x_max - x_min + 1)) - 1 ) , XDIM)
    cdef int local_y = pbc_correction((y_min + randint(1, (y_max - y_min + 1)) - 1 ) , YDIM)
    cdef int local_z = pbc_correction((z_min + randint(1, (z_max - z_min + 1)) - 1 ) , ZDIM)
    

    # leave this in - can be used to check moves maintain detailed balance
    #with open('crankshaft_moves','a') as fh:
    #    fh.write('%i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i \n' % (x_min, x_max, y_min, y_max, z_min, z_max, local_x, local_y, local_z, position_triptic[1,0] - local_x, position_triptic[1,1] - local_y, position_triptic[1,2] - local_z))


    # hard sphere clash
    if grid[local_x, local_y, local_z] > 0:
        
        
        #with open('crankshaft_fail_move.txt','a') as fh:            
        #    fh.write('%i, %i, %i \n' % (position_triptic[1,0] - local_x, position_triptic[1,1] - local_y, position_triptic[1,2] - local_z))

        # fail
        new_position[0] = -1
        return (new_position)

        
    else:
        #with open('crankshaft_success_move.txt','a') as fh:            
        #    fh.write('%i, %i, %i \n' % (position_triptic[1,0] - local_x, position_triptic[1,1] - local_y, position_triptic[1,2] - local_z))

        # success        
        new_position[0] = local_x;
        new_position[1] = local_y;
        new_position[2] = local_z;
        return (new_position)



#-----------------------------------------------------------------
# 
#cdef crank_it_good (cnp.ndarray[NUMPY_INT_TYPE, ndim=2] position_triptic, cnp.ndarray[NUMPY_INT_TYPE, ndim=3] grid, int XDIM, int YDIM, int ZDIM):
cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] crank_it_good(NUMPY_INT_TYPE[:,:] position_triptic, NUMPY_INT_TYPE[:,:,:] grid, int XDIM, int YDIM, int ZDIM):
    """
    Perform crankshaft move!

    Returns a tuple of success and the PBC corrected posistion. This is just a less optimized and differently implemented version of crank_it. Was used in debugging, but is totally legit, so keeping around incase its needed again. 

    """    
   
    cdef int N_side_x, N_side_y, N_side_z, C_side_x, C_side_y, C_side_z;
    cdef int new_x, new_y, new_z
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position = np.zeros([3], dtype=int)

    # randomly perturb in 3D
    new_x = pbc_correction(position_triptic[1,0] + (randint(0,2)-1), XDIM)
    new_y = pbc_correction(position_triptic[1,1] + (randint(0,2)-1), YDIM)
    new_z = pbc_correction(position_triptic[1,2] + (randint(0,2)-1), ZDIM)
    
    # first check if the position is empty (note this rejects if we don't move, whic is good)
    if grid[new_x, new_y, new_z] > 0:
        new_position[0] = -1
        return new_position

    new_position[0] = -1

    #print "CRANK IT GOOOOOOOOOOOD"

    if grid[new_x, new_y, new_z ] > 0:        
        return new_position
        
        
    # next check we're within 1 in each direction including PBC considerations (x)
    if (abs(new_x - position_triptic[0,0]) < 2) or (abs(new_x - position_triptic[0,0]) == XDIM-1):
           if (abs(new_x - position_triptic[2,0]) < 2) or (abs(new_x - position_triptic[2,0]) == XDIM-1):

                  # in y
                  if (abs(new_y - position_triptic[0,1]) < 2) or (abs(new_y - position_triptic[0,1]) == YDIM-1):
                         if (abs(new_y - position_triptic[2,1]) < 2) or (abs(new_y - position_triptic[2,1]) == YDIM-1):
                                
                                # in z
                                if (abs(new_z - position_triptic[0,2]) < 2) or (abs(new_z - position_triptic[0,2]) == ZDIM-1):
                                       if (abs(new_z - position_triptic[2,2]) < 2) or (abs(new_z - position_triptic[2,2]) == ZDIM-1):
                                              
                                              new_position[0] = new_x
                                              new_position[1] = new_y
                                              new_position[2] = new_z
    return new_position
                                                      

#-----------------------------------------------------------------
#
#def single_bead_crank (cnp.ndarray[NUMPY_INT_TYPE, ndim=1] old_position, cnp.ndarray[NUMPY_INT_TYPE, ndim=3] grid, int XDIM, int YDIM, int ZDIM):
def single_bead_crank (NUMPY_INT_TYPE[:] old_position, NUMPY_INT_TYPE[:,:,:] grid, NUMPY_INT_TYPE XDIM, NUMPY_INT_TYPE YDIM, NUMPY_INT_TYPE ZDIM):
    """
    Perform crankshaft move!

    Single beads are easy as the new position is just a random pertubation 
    in all three coordinates. We then update the position, returning a new
    position if there were no hardsphere clashes or a [-1, 0, 0] position 
    if there were hardsphere clases
    
    """

    ### -----------------------------------------
    # delete all of this!
    #cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position = np.zeros([3], dtype=int)
    #new_position[0] = -1
    #return new_position
    ### -----------------------------------------
    
    cdef int x_off, y_off, z_off
    x_off = (randint(0,2)-1)
    y_off = (randint(0,2)-1)
    z_off = (randint(0,2)-1)

    RETURN = update_position(old_position, grid, x_off, y_off, z_off, XDIM, YDIM, ZDIM)
    
    
    #if RETURN[0] == -1:
    #    with open('single_bead_crank_fail.txt','a') as fh:        
    #        fh.write('%i, %i, %i \n' % (x_off, y_off, z_off))
    #else:
    #    with open('single_bead_crank_success.txt','a') as fh:        
    #        fh.write('%i, %i, %i \n' % (x_off, y_off, z_off))
    




    return RETURN


    
cdef int accept_or_reject(float invtemp, long old_energy, long new_energy):
    
    cdef float expterm
    cdef float randval

    # return
    if new_energy <= old_energy:
        return 1

    # note exp here is from libc and is imported at the start
    expterm = exp(-(new_energy-old_energy)*invtemp)
    randval = float(rand())/float(RAND_MAX)

    if randval < expterm:
        return 1
    else:
        return 0


#cdef long get_angle_energy_change(NUMPY_INT_TYPE bead_index,    
@cython.wraparound(False)
@cython.boundscheck(False)
cdef long get_angle_energy_change(int bead_index,
                                  NUMPY_INT_TYPE[:,:] idx_to_bead, 
                                  NUMPY_INT_TYPE[:] new_position, 
                                  NUMPY_INT_TYPE[:,:,:,:,:,:,:] angle_lookup):

    """
    cdef long get_angle_energy_change(NUMPY_INT_TYPE bead_index,
    cnp.ndarray[NUMPY_INT_TYPE, ndim=2] idx_to_bead, 
    cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position, 
    cnp.ndarray[long, ndim=7] angle_lookup):
    """
                        
    """
    Function that takes a pre-computed angle energy lookup table and returns the change
    in angle energy before and after the move, where the move must be moving a single bead. 
    Note this function is almost purely C (as measured by the output from cython -a)

    ## Parameters are as follows:
    
    bead_index
    Index position of the specific bead of interest

    idx_to_bead
    Full matrix containing the complete system information (see the description in
    the mega_crank function description for a full summary of this variable)

    new_position
    3D np.array ([x,y,z]) that defines the position that the moved bead 
    is moving to 

    angle_lookup
    6D lookup table that translates vector distances into an angle energy
    penalty based on precomputed values.

    """

    #print "Angle energy:"
    #print bead_index
    #print idx_to_bead
    #print idx_to_bead[bead_index]


    
    # if the skip angle is set to true
    if idx_to_bead[bead_index,3] == 1:
        return 0

    # initialize a bunch of values
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] a = np.zeros([3], dtype=int)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] b = np.zeros([3], dtype=int)    
    
    cdef long angle_penalty_new = 0;
    cdef long angle_penalty_old = 0;

    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=2] angle_positions  = np.zeros([5, 3], dtype=int)    
    cdef cnp.ndarray[NUMPY_INT_TYPE, ndim=1] intcode_lookup = np.zeros([5], dtype=int)    
    cdef int offset;
    cdef int offset_start, offset_end;
    cdef int i;
    cdef int angle_idx;

    cdef int local_move_idx = -1


    ## Firstly we assess what the type of bead is, as this determines over 
    #  what range we are assessing 
    
    # N-terminal bead
    if idx_to_bead[bead_index,0] == 1:
        offset_start = 0
        offset_end   = 3

    # central bead
    elif idx_to_bead[bead_index,0] == 2:
        offset_start = -2
        offset_end   = 3

    # central bead in an OXO configuration
    elif idx_to_bead[bead_index,0] == 4:
        offset_start = -1
        offset_end   = 2

    # C terminal bead
    elif idx_to_bead[bead_index,0] == 3:
        offset_start = -2
        offset_end   = 1

    # N terminal bead +1 from start
    elif idx_to_bead[bead_index,0] == 5:
        offset_start = -1
        offset_end   = 3

    # C terminal bead -1 from end
    elif idx_to_bead[bead_index,0] == 6:
        offset_start = -2
        offset_end   = 2

    # this is the index we use to keep track of the dynamic angle_positions array
    angle_idx = 0

    # next construct a dynamic list of positions over which the angles will be evaliated
    for offset in range(offset_start,offset_end):
        pos = bead_index + offset

        angle_positions[angle_idx, 0] = idx_to_bead[pos,5] # x pos
        angle_positions[angle_idx, 1] = idx_to_bead[pos,6] # y pos
        angle_positions[angle_idx, 2] = idx_to_bead[pos,7] # z pos
        intcode_lookup[angle_idx]     = idx_to_bead[pos,2] # intcode

        if pos == bead_index:
            local_move_idx = angle_idx

        angle_idx = angle_idx + 1

        
    for i in xrange(0, (angle_idx)-2):

        # compute the two vectors between the positions
        a[0] = fix_angle_pbc_issues(angle_positions[i+1, 0] - angle_positions[i, 0])
        a[1] = fix_angle_pbc_issues(angle_positions[i+1, 1] - angle_positions[i, 1])
        a[2] = fix_angle_pbc_issues(angle_positions[i+1, 2] - angle_positions[i, 2])

        #pritn "P2 to P1 [X] --> [%i] to [%i] = %i" % (angle_positions[i+1, 0], angle_positions[i, 0], a[0])
        #print "P2 to P1 [Y] --> [%i] to [%i] = %i" % (angle_positions[i+1, 1], angle_positions[i, 1], a[1])
        #print "P2 to P1 [Z] --> [%i] to [%i] = %i" % (angle_positions[i+1, 2], angle_positions[i, 2], a[2])

        b[0] = fix_angle_pbc_issues(angle_positions[i+1, 0] - angle_positions[i+2, 0])
        b[1] = fix_angle_pbc_issues(angle_positions[i+1, 1] - angle_positions[i+2, 1])
        b[2] = fix_angle_pbc_issues(angle_positions[i+1, 2] - angle_positions[i+2, 2])

        angle_penalty_old = angle_lookup[intcode_lookup[i+1], a[0]+1, a[1]+1, a[2]+1, b[0]+1, b[1]+1, b[2]+1] + angle_penalty_old

    angle_positions[local_move_idx,0] = new_position[0]
    angle_positions[local_move_idx,1] = new_position[1]
    angle_positions[local_move_idx,2] = new_position[2]

    for i in xrange(0, (angle_idx)-2):

        # compute the two vectors between the positions
        a[0] = fix_angle_pbc_issues(angle_positions[i+1, 0] - angle_positions[i, 0])
        a[1] = fix_angle_pbc_issues(angle_positions[i+1, 1] - angle_positions[i, 1])
        a[2] = fix_angle_pbc_issues(angle_positions[i+1, 2] - angle_positions[i, 2])

        b[0] = fix_angle_pbc_issues(angle_positions[i+1, 0] - angle_positions[i+2, 0])
        b[1] = fix_angle_pbc_issues(angle_positions[i+1, 1] - angle_positions[i+2, 1])
        b[2] = fix_angle_pbc_issues(angle_positions[i+1, 2] - angle_positions[i+2, 2])

        angle_penalty_new = angle_lookup[intcode_lookup[i+1], a[0]+1, a[1]+1, a[2]+1, b[0]+1, b[1]+1, b[2]+1] + angle_penalty_new
        
    return (angle_penalty_new - angle_penalty_old)    

    

@cython.cdivision(True)
cdef int fix_angle_pbc_issues(int distance):
    """
    Hack that takes advantage of the fact that the distances all must be -1,0, +1 and 
    we're always computing a X->O vector where O is the bend point

    """
    if distance < -1:
        return 1
    elif distance > 1:
        return -1
    else:
        return distance
        



# these directives made the function slower...
@cython.wraparound(False)
@cython.boundscheck(False)
cdef get_energy_change(NUMPY_INT_TYPE[:,:,:] grid, 
                       NUMPY_INT_TYPE[:,:,:] type_grid,
                       NUMPY_INT_TYPE[:] old_position,
                       NUMPY_INT_TYPE[:] new_position,
                       int LR_vs_SR,
                       NUMPY_INT_TYPE[:,:] interaction_table, 
                       NUMPY_INT_TYPE[:,:] LR_interaction_table,
                       NUMPY_INT_TYPE[:,:] SLR_interaction_table,
                       int XDIM, 
                       int YDIM, 
                       int ZDIM,
                       int hardwall):

    """
    cdef get_energy_change(cnp.ndarray[NUMPY_INT_TYPE, ndim=3] grid, 
    cnp.ndarray[NUMPY_INT_TYPE, ndim=3] type_grid,
    cnp.ndarray[NUMPY_INT_TYPE, ndim=1] old_position,
    cnp.ndarray[NUMPY_INT_TYPE, ndim=1] new_position,
    int LR_vs_SR,
    cnp.ndarray[long, ndim=2] interaction_table, 
    cnp.ndarray[long, ndim=2] LR_interaction_table,
    cnp.ndarray[long, ndim=2] SLR_interaction_table,
    int XDIM, 
    int YDIM, 
    int ZDIM,
    int hardwall):
    """
    
    # We define four energy variables which are going to be used to calcualte
    # the DELTA energy. 
    #
    # --> energy old is the interactions between the bead and its
    # partners in its current location
    #
    # --> energy_old_empy is the energy associated with the site after this
    # bead has been removed (i.e. solvent to local partners only)
    #
    # --> energy_new is the interaction between the bead and its partners in
    #     its NEW location!
    #
    # --> energy_new_empty is the interaction between the bead and its 
    
    cdef long energy_old       = 0
    cdef long energy_old_empty = 0
    cdef long energy_new       = 0
    cdef long energy_new_empty = 0 
    

    cdef int old_x, old_y, old_z;
    cdef int new_x, new_y, new_z;
    cdef int tmp_x, tmp_y, tmp_z;

    cdef int x, y, z

    cdef unsigned int site_bead_type, bead_type;

    # extract 
    new_x = new_position[0]
    new_y = new_position[1]
    new_z = new_position[2]

    old_x = old_position[0]
    old_y = old_position[1]
    old_z = old_position[2]

    # get the type of the current beads
    bead_type = type_grid[old_x, old_y, old_z]
    
    ## ************************************************************************************************
    ##  Long range energy
    ##
    # if we're looking at a bead which engages in longrange interactions 
    if LR_vs_SR == 1:

        ## First evaluate the old energy        
        for x in range(-3,4):
            for y in range(-3,4):
                for z in range(-3,4):

                    tmp_x = pbc_correction(old_x + x, XDIM)
                    tmp_y = pbc_correction(old_y + y, YDIM)
                    tmp_z = pbc_correction(old_z + z, ZDIM)
                                                
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]
                    
                    # short-range interactions
                    if abs(x) < 2 and abs(y) < 2 and abs(z) < 2:

                        # if harwall boundary is on and the site bead is jumping over the PBC then it is by
                        # definition solvent (type = 0)
                        if hardwall == 1:
                            if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                                site_bead_type=0
                        
                        # bead-site interaction                 
                        energy_old = energy_old + interaction_table[bead_type, site_bead_type]

                        # might be faster to just read this from memory every time - is if or a memory read and addition faster?
                        if site_bead_type > 0:
                            energy_old_empty = energy_old_empty + interaction_table[0, site_bead_type]

                    # long range interactions (maybe optimize this in future...) - note long range interactions cannot
                    # contribute to the energy_old_empty variable
                    elif abs(x) < 3 and abs(y) < 3 and abs(z) <3:

                        # if harwall boundary is on - don't feel LR or SLR interactions over boundary
                        if hardwall == 1:
                            if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                                continue

                        energy_old = energy_old + LR_interaction_table[bead_type, site_bead_type]

                    else:

                        # if harwall boundary is on - don't feel LR or SLR interactions over boundary
                        if hardwall == 1:
                            if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                                continue

                        energy_old = energy_old + SLR_interaction_table[bead_type, site_bead_type]



        # reset the typegrid to the new configuation
        type_grid[new_x,new_y,new_z] = type_grid[old_x,old_y,old_z]
        type_grid[old_x,old_y,old_z] = 0

        ## then evaluate the new energy        
        for x in range(-3,4):
            for y in range(-3,4):
                for z in range(-3,4):

                    # get PBC corrected positions of interest
                    tmp_x = pbc_correction(new_x + x, XDIM)
                    tmp_y = pbc_correction(new_y + y, YDIM)
                    tmp_z = pbc_correction(new_z + z, ZDIM)

                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]
                    
                    # short-range interactions
                    if abs(x) < 2 and abs(y) < 2 and abs(z) <2:
                        
                        # if harwall boundary is on... (see previous codeblock for an explanation)
                        if hardwall == 1:
                            if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                                site_bead_type=0
                        
                        # bead-site interaction                 
                        energy_new = energy_new + interaction_table[bead_type, site_bead_type]
                        
                        if site_bead_type > 0:
                            energy_new_empty = energy_new_empty + interaction_table[0, site_bead_type]

                    # long range interactions (maybe optimize this in future...) - note long range interactions cannot
                    # contribute to the energy_old_empty variable
                    elif abs(x) < 3 and abs(y) < 3 and abs(z) <3:
                        
                        # if harwall boundary is on - skip LR interactions
                        if hardwall == 1:
                            if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                                continue

                        energy_new = energy_new + LR_interaction_table[bead_type, site_bead_type]

                    else:

                        # if harwall boundary is on - skip LR interactions
                        if hardwall == 1:
                            if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                                continue

                        energy_new = energy_new + SLR_interaction_table[bead_type, site_bead_type]


    ## ************************************************************************************************
    ##  only short range brah!
    ##
    else:
        for x in range(-1,2):
            for y in range(-1,2):
                for z in range(-1,2):

                    # get PBC corrected positions of interest
                    tmp_x = pbc_correction(old_x + x, XDIM)
                    tmp_y = pbc_correction(old_y + y, YDIM)
                    tmp_z = pbc_correction(old_z + z, ZDIM)

             
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]
                    
                    # if harwall boundary is on... then interactions that would straddle the boundary
                    # must be solvent, so overwrite
                    if hardwall == 1:
                        if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                            site_bead_type = 0
                                                                         
                    # bead-site interaction                 
                    energy_old = energy_old + interaction_table[bead_type, site_bead_type]

                    # might be faster to just read this from memory every time - is if or a memory read and addition faster?
                    if site_bead_type > 0:
                        energy_old_empty = energy_old_empty + interaction_table[0, site_bead_type]
            
        # reset the typegrid to the new configuation
        type_grid[new_x,new_y,new_z] = type_grid[old_x,old_y,old_z]
        type_grid[old_x,old_y,old_z] = 0
                        
        for x in range(-1,2):
            for y in range(-1,2):
                for z in range(-1,2):
                    
                    # get PBC corrected positions of interest
                    tmp_x = pbc_correction(new_x + x, XDIM)
                    tmp_y = pbc_correction(new_y + y, YDIM)
                    tmp_z = pbc_correction(new_z + z, ZDIM)

                    #tmp_x = (new_x + x) % XDIM
                    #tmp_y = (new_y + y) % YDIM
                    #tmp_z = (new_z + z) % ZDIM

                    # if harwall boundary is on...
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]
                    if hardwall == 1:
                        if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                            site_bead_type = 0
                    
                    # bead-site interaction                 
                    energy_new = energy_new + interaction_table[bead_type, site_bead_type]

                    # might be faster to just read this from memory every time - is if or a memory read and addition faster?
                    if site_bead_type > 0:
                        energy_new_empty = energy_new_empty + interaction_table[0, site_bead_type]
                        
                        
    # finally correct for the self-self interactions we cover when x == y == z == 0 
    energy_old_empty = energy_old_empty - interaction_table[0, bead_type]
    energy_old       = energy_old - interaction_table[bead_type, bead_type]
    energy_new_empty = energy_new_empty - interaction_table[0, bead_type]
    energy_new       = energy_new - interaction_table[bead_type, bead_type]

    # reset the type grid again...
    type_grid[old_x,old_y,old_z] = type_grid[new_x,new_y,new_z] 
    type_grid[new_x,new_y,new_z] = 0


    return (energy_new + energy_old_empty) - (energy_old + energy_new_empty)

@cython.cdivision(True)
cdef int pbc_correction(int value, int DIM):    
    """
    Performs intelligent periodic boundary correction
    which FIRST checks to see if we have a negative
    value and IF NOT uses the % operator - this means
    we can use % without checking the sign giving
    a 35% speedup per call
    """
    
    if value < 0:
        return DIM+value
    else:
        return (value % DIM)



@cython.wraparound(False)
@cython.boundscheck(False)
cdef int do_positions_stradle_pbc_boundary(cnp.ndarray[NUMPY_INT_TYPE, ndim=2] chain_positions, int chain_length):
    """                                       
    For a set of positions returns true if the positions straddle a boundary
    else return false

    """

    cdef int pidx
    cdef int p1_x, p2_x, p1_y, p2_y, p1_z, p2_z
    
    

    for pidx in range(0, chain_length-1):
        p1_x = chain_positions[pidx,0]
        p2_x = chain_positions[pidx+1,0]

        p1_y = chain_positions[pidx,1]
        p2_y = chain_positions[pidx+1,1]

        p1_z = chain_positions[pidx,2]
        p2_z = chain_positions[pidx+1,2]
        
        if abs(p1_x - p2_x) > 1:
            return 1

        if abs(p1_y - p2_y) > 1:
            return 1

        if abs(p1_z - p2_z) > 1:
            return 1

    return 0



### 
###        EXTERNAL FUNCTIONS
###

def accept_or_reject_ext(float invtemp, long old_energy, long new_energy):
    return accept_or_reject(invtemp, old_energy, new_energy)


def get_random_position_ext(int chain_length):
    return get_random_position(chain_length)


def randint_ext(int start, int end):
    return randint(start,end)


##
## TESTER FUNCTIONS
##                        
            



def randint_tester(int start, int end, int randval):
    """
    returns a random integer between start and end-1 (exclusive
    of ends
    This mirrors random.randint's behaviour
    
    e.g. start = 0 and end=5 gives one of 0,1,2,3,4,5

    """
    if start == 0:
        end=end+1
    else:
        pass
         
    # ok so this is kind of inelegant, but the rand()-1 is actually important
    # if we don'te have this we risk the situation where 
    cdef  int r = start+int((float(randval-1)/float(RAND_MAX))*(end))

    return r

def crand_test():
    return rand()

def RAND_MAX_test():
    return RAND_MAX
