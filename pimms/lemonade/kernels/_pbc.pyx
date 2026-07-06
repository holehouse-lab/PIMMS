## ...........................................................................
##
## lemonade - analysis backend for PIMMS lattice trajectories
##
## Compiled kernels for the performance-critical trajectory operations: batched
## PBC unwrapping ("make whole") and per-frame grid painting. These replace the
## per-bead / per-frame pure-Python loops that dominated the original lemonade.
##
## ...........................................................................

import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def unwrap_chains(cnp.int32_t[:, :, ::1] positions,
                  cnp.int64_t[::1] offsets,
                  cnp.int64_t[::1] dims,
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
    positions : (n_frames, n_atoms, 3) int32
        Wrapped, in-box integer lattice positions (z is 0 for 2D systems).
    offsets : (n_chains + 1,) int64
        CSR-style atom index boundaries, one contiguous block of atoms per chain.
    dims : (3,) int64
        Box size per axis (dims[2] is ignored / may be 1 for 2D).
    n_dim : int
        2 or 3.

    Returns
    -------
    (n_frames, n_atoms, 3) int32
        Unwrapped ("whole") positions, first bead of each chain unchanged.
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
def paint_frame_grid_3d(cnp.int32_t[:, ::1] frame_positions,
                        cnp.int32_t[::1] atom_chainid,
                        cnp.int32_t[:, :, ::1] grid):
    """
    Paint one frame's beads onto a 3D grid: ``grid[x, y, z] = chainID`` for every
    occupied site (0 = empty). ``grid`` must be pre-zeroed and sized to the box.
    Positions must already be wrapped into the box.
    """
    cdef Py_ssize_t n = frame_positions.shape[0]
    cdef Py_ssize_t i
    for i in range(n):
        grid[frame_positions[i, 0], frame_positions[i, 1], frame_positions[i, 2]] = atom_chainid[i]


@cython.boundscheck(False)
@cython.wraparound(False)
def paint_frame_grid_2d(cnp.int32_t[:, ::1] frame_positions,
                        cnp.int32_t[::1] atom_chainid,
                        cnp.int32_t[:, ::1] grid):
    """2D counterpart of :func:`paint_frame_grid_3d` (ignores the z column)."""
    cdef Py_ssize_t n = frame_positions.shape[0]
    cdef Py_ssize_t i
    for i in range(n):
        grid[frame_positions[i, 0], frame_positions[i, 1]] = atom_chainid[i]
