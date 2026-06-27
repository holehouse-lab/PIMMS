## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................
# 

import numpy as np
import random


def position_in_list(position, list_of_positions):
    """
    Return True if a position exists in a list of positions.

    Performs an element-wise comparison of ``position`` against every
    entry in ``list_of_positions``. The single-position case (where
    ``list_of_positions`` is itself just one position rather than a
    list of positions) is handled explicitly.

    Parameters
    ----------
    position : array_like
        Coordinate to search for (e.g. a length-2 or length-3 sequence).

    list_of_positions : array_like
        Either a list/array of positions to search within, or a single
        position of the same dimensionality as ``position``.

    Returns
    -------
    bool
        True if ``position`` is found, False otherwise (including when
        ``list_of_positions`` is empty).

    """
    if len(list_of_positions) == 0:
        return False

    position_array = np.asarray(position)
    positions_array = np.asarray(list_of_positions)

    # Handle the case where a single position (not a list) is passed.
    if positions_array.ndim == position_array.ndim:
        return bool(np.array_equal(position_array, positions_array))

    return bool(np.any(np.all(positions_array == position_array, axis=1)))
    

def randneg(val):
    """
    Randomly return the value or its negation with equal probability.

    Parameters
    ----------
    val : numeric
        The value to (possibly) negate.

    Returns
    -------
    numeric
        ``-val`` with probability 0.5, otherwise ``val``.

    """
    if random.random() > 0.5:
        return -val
    else:
        return val

def tetrahedron_volume(a, b, c, d):
    """
    Compute the volume of a tetrahedron defined by four 3D vertices.

    Uses the scalar-triple-product formula. Inputs are broadcast to at
    least 2D so that stacks of tetrahedra (arrays of vertices) can be
    handled in a single vectorized call.

    From http://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy

    Parameters
    ----------
    a, b, c, d : array_like
        The four vertices of the tetrahedron. Each must be a 3D
        coordinate (or a stack of 3D coordinates) and all four must
        share the same shape.

    Returns
    -------
    numpy.ndarray
        The tetrahedron volume(s) as a 1D array (one entry per
        tetrahedron).

    Raises
    ------
    ValueError
        If the four points do not share the same shape, or if they are
        not 3D coordinates.

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
    Find the value nearest to a target in an array.

    The array is flattened before searching, so the returned index is
    into the raveled (1D) view of the input.

    Parameters
    ----------
    array : array_like
        The array to search. Must be non-empty.

    target : numeric
        The value to find the closest element to.

    Returns
    -------
    tuple
        A 2-tuple ``(idx, value)`` where ``idx`` is the integer index
        (into the flattened array) of the nearest element and ``value``
        is the element itself.

    Raises
    ------
    ValueError
        If ``array`` is empty.

    """
    array = np.asarray(array)
    if array.size == 0:
        raise ValueError("Cannot find nearest value in an empty array")

    flat = array.ravel()
    idx = int((np.abs(flat - target)).argmin())
    return (idx, flat[idx])
