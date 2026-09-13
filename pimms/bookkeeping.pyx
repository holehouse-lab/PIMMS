# cython: language_level=3
## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""Compiled bookkeeping between the Chain objects and the bead table.

Every megamove hands the compiled kernels one integer table with a row per bead
(the ``idx_to_bead`` matrix: static columns describing the bead, then its
coordinates), and reads the moved coordinates back out of it afterwards. The
Chain objects, which everything else in PIMMS reads, keep their positions as a
Python list of ``[x, y(, z)]`` lists. So each megamove has to copy every bead's
coordinates from lists into the table before the kernel runs and back into lists
after it returns.

Done in Python that copy is not a tax on the megamove, it *is* the megamove: at
12,000 beads it cost about 4 ms against 2 ms for the parallel crankshaft kernel
itself, because the natural way to write it, one ``np.array(chain.positions)``
and one ``table[a:b].tolist()`` per chain, creates a numpy array and a list
object for every one of the thousands of chains. That was what capped
``PARALLELIZE`` at three to five times when the kernels alone reached eight.

The two functions here do the same copies as tight C loops over the Python list
objects directly, with no temporaries. They change nothing about what is
copied, so a run's trajectory is byte-identical with or without this module;
:mod:`pimms.crankshaft_list_functions` falls back to the pure-Python loops if
the extension is not built.
"""

cimport cython
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_New, PyList_SET_ITEM
from cpython.long cimport PyLong_AsLong, PyLong_FromLong
from cpython.ref cimport Py_INCREF

import numpy as np
cimport numpy as cnp

cnp.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
def gather_positions(cnp.int64_t[:, ::1] table, list chains,
                     cnp.int64_t[::1] offsets, cnp.int64_t[::1] lengths, int n_dim):
    """Copy every chain's positions into the coordinate columns of the bead table.

    Parameters
    ----------
    table : numpy.ndarray
        The ``idx_to_bead`` matrix, ``int64`` and C-contiguous, with the
        coordinate columns starting at column 5. Modified in place; the static
        columns 0 to 4 are not touched.

    chains : list
        The Chain objects in the order their rows appear in ``table`` (ascending
        chainID). Each must expose ``positions`` as a list of lists of ints.

    offsets : numpy.ndarray
        ``int64`` row index in ``table`` of each chain's first bead, one entry per
        chain, in the same order as ``chains``.

    lengths : numpy.ndarray
        ``int64`` number of beads of each chain, in the same order.

    n_dim : int
        2 or 3, the number of coordinate columns to copy.

    Returns
    -------
    None
        The table is updated in place.

    Raises
    ------
    ValueError
        If a chain holds a different number of positions than ``lengths`` says,
        or a position has fewer than ``n_dim`` coordinates. Either means the
        table and the chains have gone out of step, which must never be papered
        over by a partial copy.
    """
    cdef Py_ssize_t n_chains = len(chains)
    cdef Py_ssize_t ci, i, d, off, L
    cdef object chain, positions, row

    if offsets.shape[0] != n_chains or lengths.shape[0] != n_chains:
        raise ValueError("gather_positions: offsets/lengths do not match the number of chains")
    if table.shape[1] < 5 + n_dim:
        raise ValueError("gather_positions: the table has too few columns for %d dimensions" % n_dim)

    for ci in range(n_chains):
        chain = chains[ci]
        positions = chain.positions
        off = offsets[ci]
        L = lengths[ci]
        if not isinstance(positions, list) or PyList_GET_SIZE(positions) != L:
            raise ValueError(
                "gather_positions: chain %r holds %d positions but the bead table has "
                "%d rows for it - the table and the chains have gone out of step"
                % (getattr(chain, "chainID", ci), len(positions), L))
        if off < 0 or off + L > table.shape[0]:
            raise ValueError("gather_positions: chain rows fall outside the table")
        for i in range(L):
            row = <object>PyList_GET_ITEM(positions, i)
            if not isinstance(row, list) or PyList_GET_SIZE(row) < n_dim:
                raise ValueError(
                    "gather_positions: bead %d of chain %r is not a %d-coordinate list"
                    % (i, getattr(chain, "chainID", ci), n_dim))
            for d in range(n_dim):
                table[off + i, 5 + d] = PyLong_AsLong(<object>PyList_GET_ITEM(row, d))


@cython.boundscheck(False)
@cython.wraparound(False)
def scatter_positions(cnp.int64_t[:, ::1] table, list chains,
                      cnp.int64_t[::1] offsets, cnp.int64_t[::1] lengths, int n_dim):
    """Copy the coordinate columns of the bead table back into every chain.

    Builds, for each chain, a fresh list of ``[x, y(, z)]`` lists from its rows
    of the table and hands it to ``chain.set_ordered_positions``, which is the
    same call the Python loop made and keeps that method's length check.

    Parameters
    ----------
    table : numpy.ndarray
        The ``idx_to_bead`` matrix as returned by a kernel, ``int64`` and
        C-contiguous, coordinates from column 5 on. Read only.

    chains : list
        The Chain objects in the order their rows appear in ``table``.

    offsets : numpy.ndarray
        ``int64`` row index of each chain's first bead, in that order.

    lengths : numpy.ndarray
        ``int64`` number of beads of each chain, in that order.

    n_dim : int
        2 or 3, the number of coordinate columns to copy.

    Returns
    -------
    None
        Every chain's positions are replaced.

    Raises
    ------
    ValueError
        If the offsets and lengths do not describe the table.
    """
    cdef Py_ssize_t n_chains = len(chains)
    cdef Py_ssize_t ci, i, d, off, L
    cdef list chain_positions, row
    cdef object value

    if offsets.shape[0] != n_chains or lengths.shape[0] != n_chains:
        raise ValueError("scatter_positions: offsets/lengths do not match the number of chains")
    if table.shape[1] < 5 + n_dim:
        raise ValueError("scatter_positions: the table has too few columns for %d dimensions" % n_dim)

    for ci in range(n_chains):
        off = offsets[ci]
        L = lengths[ci]
        if off < 0 or off + L > table.shape[0]:
            raise ValueError("scatter_positions: chain rows fall outside the table")
        chain_positions = PyList_New(L)
        for i in range(L):
            row = PyList_New(n_dim)
            for d in range(n_dim):
                value = PyLong_FromLong(table[off + i, 5 + d])
                Py_INCREF(value)
                PyList_SET_ITEM(row, d, value)
            Py_INCREF(row)
            PyList_SET_ITEM(chain_positions, i, row)
        chains[ci].set_ordered_positions(chain_positions)
