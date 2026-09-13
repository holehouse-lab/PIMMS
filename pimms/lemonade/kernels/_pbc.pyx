## ...........................................................................
##
## lemonade - analysis backend for PIMMS lattice trajectories
##
## Compiled kernels for the performance-critical trajectory operations: batched
## PBC unwrapping ("make whole") and per-frame grid painting. These replace the
## per-bead / per-frame pure-Python loops that dominated the original lemonade.
##
## Copyright 2015 - 2026
## ...........................................................................

import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def unwrap_chains(const cnp.int32_t[:, :, ::1] positions,
                  const cnp.int64_t[::1] offsets,
                  const cnp.int64_t[::1] dims,
                  int n_dim):
    """
    Make every chain whole across periodic boundaries, in every frame, at once.

    For each frame and each chain (the atoms ``offsets[c] .. offsets[c+1]``), the
    first bead is left in place and every subsequent bead is shifted by whole
    multiples of the box so that consecutive beads never jump across a boundary -
    i.e. the chain becomes spatially contiguous (coordinates may fall outside the
    box). This is the batched, typed-C equivalent of
    ``pimms.lattice_utils.make_chain_whole`` applied to the whole trajectory; it
    replaces the per-call ``copy.deepcopy`` + triple-nested Python loop that made
    single-image conversion the load-time bottleneck.

    Parameters
    ----------
    positions : (n_frames, n_atoms, 3) const int32 memoryview
        Wrapped, in-box integer lattice positions, C-contiguous (z is 0 for 2D
        systems).
    offsets : (n_chains + 1,) const int64 memoryview
        CSR-style atom index boundaries, one contiguous block of atoms per chain.
    dims : (3,) const int64 memoryview
        Box size per axis (dims[2] is ignored / may be 1 for 2D).
    n_dim : int (C int)
        2 or 3. Only the first n_dim columns are unwrapped.

    Returns
    -------
    (n_frames, n_atoms, 3) numpy.ndarray, int32
        Unwrapped ("whole") positions, first bead of each chain unchanged.

    Raises
    ------
    ValueError
        If a bond cannot be resolved to a unit step by shifting whole boxes,
        which means the input chain has a broken or non-unit bond.
    """
    cdef Py_ssize_t n_frames = positions.shape[0]
    cdef Py_ssize_t n_atoms = positions.shape[1]
    cdef Py_ssize_t n_chains = offsets.shape[0] - 1
    cdef Py_ssize_t f, c, a, a0, a1
    cdef int d
    cdef long dim, cur, v, guard

    out_np = np.asarray(positions).copy()
    cdef cnp.int32_t[:, :, ::1] out = out_np

    for f in range(n_frames):
        for c in range(n_chains):
            a0 = offsets[c]
            a1 = offsets[c + 1]
            if a0 == a1:
                # defensive: a zero-length chain has no anchor bead - without this
                # the unconditional anchor read below would index one row past the
                # buffer for an empty FINAL chain (boundscheck is off)
                continue
            for d in range(n_dim):
                dim = dims[d]
                cur = positions[f, a0, d]           # anchor: first bead unchanged
                for a in range(a0 + 1, a1):
                    v = positions[f, a, d]
                    if cur - v > 1:
                        # neighbour sits across the +boundary; walk it up
                        v += dim
                        guard = 0
                        while v - cur > 1 or cur - v > 1:
                            v += dim
                            guard += 1
                            if guard > 100000:
                                raise ValueError(
                                    "unwrap_chains: unresolvable bond (impossible/non-unit step in chain)")
                    elif cur - v < -1:
                        # neighbour sits across the -boundary; walk it down
                        v -= dim
                        guard = 0
                        while v - cur > 1 or cur - v > 1:
                            v -= dim
                            guard += 1
                            if guard > 100000:
                                raise ValueError(
                                    "unwrap_chains: unresolvable bond (impossible/non-unit step in chain)")
                    out[f, a, d] = v
                    cur = v
    return out_np


@cython.boundscheck(False)
@cython.wraparound(False)
def paint_frame_grid_3d(const cnp.int32_t[:, ::1] frame_positions,
                        const cnp.int32_t[::1] atom_chainid,
                        cnp.int32_t[:, :, ::1] grid):
    """
    Paint one frame's beads onto a 3D grid: ``grid[x, y, z] = chainID`` for every
    occupied site (0 = empty). ``grid`` must be pre-zeroed and sized to the box.
    Positions must already be wrapped into the box.

    Bounds checking is off, so an out-of-box position writes past the end of
    ``grid``; the caller is responsible for wrapping (``lemonade.load`` does this
    at load time).

    Parameters
    ----------
    frame_positions : (n_atoms, 3) const int32 memoryview
        In-box integer lattice positions for a single frame, C-contiguous.
    atom_chainid : (n_atoms,) const int32 memoryview
        Site value to write for each bead. Callers pass chain index + 1 so that
        0 is reserved for empty.
    grid : (XDIM, YDIM, ZDIM) int32 memoryview
        Pre-zeroed occupancy grid, C-contiguous, written in place.

    Returns
    -------
    None
        ``grid`` is modified in place.
    """
    cdef Py_ssize_t n = frame_positions.shape[0]
    cdef Py_ssize_t i
    for i in range(n):
        grid[frame_positions[i, 0], frame_positions[i, 1], frame_positions[i, 2]] = atom_chainid[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def paint_frame_grid_2d(const cnp.int32_t[:, ::1] frame_positions,
                        const cnp.int32_t[::1] atom_chainid,
                        cnp.int32_t[:, ::1] grid):
    """
    2D counterpart of :func:`paint_frame_grid_3d` (ignores the z column).

    Parameters
    ----------
    frame_positions : (n_atoms, 3) const int32 memoryview
        In-box integer lattice positions for a single frame, C-contiguous. Only
        columns 0 and 1 are read.
    atom_chainid : (n_atoms,) const int32 memoryview
        Site value to write for each bead. Callers pass chain index + 1 so that
        0 is reserved for empty.
    grid : (XDIM, YDIM) int32 memoryview
        Pre-zeroed occupancy grid, C-contiguous, written in place.

    Returns
    -------
    None
        ``grid`` is modified in place.
    """
    cdef Py_ssize_t n = frame_positions.shape[0]
    cdef Py_ssize_t i
    for i in range(n):
        grid[frame_positions[i, 0], frame_positions[i, 1]] = atom_chainid[i]
