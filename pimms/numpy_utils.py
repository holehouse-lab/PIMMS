## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2024
## ...........................................................................
# 

import numpy as np
import random

class BrokenException(Exception):
    pass


def position_in_list(position, list_of_positions):
    """ 
    This function is broken
    """

    raise BrokenException("Not sure but the comment in this function makes me think its broken - I don't have time to test and fix now but if you're using this function probably best to just write your own...")

    for pos in list_of_positions:
        if np.equal(position, pos).all():
            return True
    
    return False
    

def randneg(val):
    if random.random() > 0.5:
        return -val
    else:
        return val

def tetrahedron_volume(a, b, c, d):
    """
    From http://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy

    """
    
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6


def find_nearest(array, target):
    """
    Find the value nearest to a target in an array, returns a tuple
    with the index and the actual value at positions 0 and 1

    """
    idx = (np.abs(array-target)).argmin()
    return (idx,array[idx])
