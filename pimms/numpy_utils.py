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
    """Return True if a position exists in a list of positions."""
    if len(list_of_positions) == 0:
        return False

    position_array = np.asarray(position)
    positions_array = np.asarray(list_of_positions)

    # Handle the case where a single position (not a list) is passed.
    if positions_array.ndim == position_array.ndim:
        return bool(np.array_equal(position_array, positions_array))

    return bool(np.any(np.all(positions_array == position_array, axis=1)))
    

def randneg(val):
    if random.random() > 0.5:
        return -val
    else:
        return val

def tetrahedron_volume(a, b, c, d):
    """
    From http://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy

    """
    
    a = np.asarray(a)
    b = np.asarray(b)
    c = np.asarray(c)
    d = np.asarray(d)

    if not (a.shape == b.shape == c.shape == d.shape):
        raise ValueError("All tetrahedron points must have the same shape")

    a = np.atleast_2d(a)
    b = np.atleast_2d(b)
    c = np.atleast_2d(c)
    d = np.atleast_2d(d)

    if a.shape[1] != 3:
        raise ValueError("Tetrahedron points must be 3D coordinates")

    return np.abs(np.einsum('ij,ij->i', a - d, np.cross(b - d, c - d))) / 6


def find_nearest(array, target):
    """
    Find the value nearest to a target in an array, returns a tuple
    with the index and the actual value at positions 0 and 1

    """
    array = np.asarray(array)
    if array.size == 0:
        raise ValueError("Cannot find nearest value in an empty array")

    flat = array.ravel()
    idx = int((np.abs(flat - target)).argmin())
    return (idx, flat[idx])
