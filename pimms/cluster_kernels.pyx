## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2026
##
## ...........................................................................
##
## Compiled kernels for cluster-analysis hot paths.
##

import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def snakesearch_single_image(cnp.int64_t[:, ::1] positions,
                             cnp.int64_t[::1] dims,
                             Py_ssize_t seed_idx,
                             int space_threshold):
    """
    Breadth-first "snakesearch" reconstruction of a cluster into a single periodic
    image.

    This is the compiled equivalent of the pure-Python
    ``cluster_utils.convert_positions_to_single_image_snakesearch`` and produces
    byte-for-byte identical output for a given ``seed_idx``. Starting from the seed
    bead, each unvisited neighbour within ``space_threshold`` (per axis, PBC-aware)
    is placed into the same periodic image as the reference bead and enqueued; the
    result is finally shifted so every coordinate is >= 0.

    Neighbour discovery uses a flat occupancy-index grid of size ``prod(dims)``
    (position -> bead index, O(1) lookup) instead of a Python dict, and the whole
    walk runs in typed C with integer arithmetic - removing the millions of tuple
    constructions / dict lookups / per-dimension Python loops that dominated the
    original.

    Parameters
    ----------
    positions : (N, n_dim) int64 memoryview
        Bead positions in PBC space (each coordinate in ``[0, dims[d])``).
    dims : (n_dim,) int64 memoryview
        Box dimensions per axis.
    seed_idx : int
        Index of the bead to start the walk from (its image is preserved).
    space_threshold : int
        Maximum per-axis distance for two beads to count as neighbours.

    Returns
    -------
    numpy.ndarray
        ``(N, n_dim)`` int64 array of single-image positions, shifted so every
        coordinate is >= 0.

    Raises
    ------
    ValueError
        If the beads do not form a single connected cluster within
        ``space_threshold`` (i.e. not every bead was reached).
    """
    cdef Py_ssize_t N = positions.shape[0]
    cdef int n_dim = positions.shape[1]
    cdef Py_ssize_t i, j, ref, head, tail, found
    cdef int t = space_threshold
    cdef int d
    cdef long Dx = dims[0]
    cdef long Dy = dims[1]
    cdef long Dz = 1
    cdef long ox, oy, oz
    cdef long rx, ry, rz, rpx, rpy, rpz, nx, ny, nz, delta, mn

    si = np.empty((N, n_dim), dtype=np.int64)
    cdef cnp.int64_t[:, ::1] si_v = si

    # array-backed BFS queue + visited flags
    cdef cnp.int64_t[::1] q = np.empty(N, dtype=np.int64)
    cdef cnp.uint8_t[::1] vis = np.zeros(N, dtype=np.uint8)

    # flat occupancy-index grid: occ[encode(position)] = bead index (or -1)
    cdef cnp.int64_t[::1] occ

    if n_dim == 3:
        Dz = dims[2]
        occ = np.full(Dx * Dy * Dz, -1, dtype=np.int64)
        for i in range(N):
            occ[(positions[i, 0] * Dy + positions[i, 1]) * Dz + positions[i, 2]] = i

        si_v[seed_idx, 0] = positions[seed_idx, 0]
        si_v[seed_idx, 1] = positions[seed_idx, 1]
        si_v[seed_idx, 2] = positions[seed_idx, 2]
        vis[seed_idx] = 1
        head = 0
        tail = 0
        q[tail] = seed_idx
        tail += 1
        found = 1

        while head < tail:
            ref = q[head]
            head += 1
            rx = si_v[ref, 0]
            ry = si_v[ref, 1]
            rz = si_v[ref, 2]
            rpx = rx % Dx
            if rpx < 0:
                rpx += Dx
            rpy = ry % Dy
            if rpy < 0:
                rpy += Dy
            rpz = rz % Dz
            if rpz < 0:
                rpz += Dz

            for ox in range(-t, t + 1):
                nx = (rpx + ox) % Dx
                if nx < 0:
                    nx += Dx
                for oy in range(-t, t + 1):
                    ny = (rpy + oy) % Dy
                    if ny < 0:
                        ny += Dy
                    for oz in range(-t, t + 1):
                        nz = (rpz + oz) % Dz
                        if nz < 0:
                            nz += Dz
                        j = occ[(nx * Dy + ny) * Dz + nz]
                        if j < 0 or vis[j]:
                            continue
                        vis[j] = 1
                        found += 1

                        delta = positions[j, 0] - rpx
                        if 2 * delta > Dx:
                            delta -= Dx
                        elif 2 * delta < -Dx:
                            delta += Dx
                        si_v[j, 0] = rx + delta

                        delta = positions[j, 1] - rpy
                        if 2 * delta > Dy:
                            delta -= Dy
                        elif 2 * delta < -Dy:
                            delta += Dy
                        si_v[j, 1] = ry + delta

                        delta = positions[j, 2] - rpz
                        if 2 * delta > Dz:
                            delta -= Dz
                        elif 2 * delta < -Dz:
                            delta += Dz
                        si_v[j, 2] = rz + delta

                        q[tail] = j
                        tail += 1

    else:
        # 2D
        occ = np.full(Dx * Dy, -1, dtype=np.int64)
        for i in range(N):
            occ[positions[i, 0] * Dy + positions[i, 1]] = i

        si_v[seed_idx, 0] = positions[seed_idx, 0]
        si_v[seed_idx, 1] = positions[seed_idx, 1]
        vis[seed_idx] = 1
        head = 0
        tail = 0
        q[tail] = seed_idx
        tail += 1
        found = 1

        while head < tail:
            ref = q[head]
            head += 1
            rx = si_v[ref, 0]
            ry = si_v[ref, 1]
            rpx = rx % Dx
            if rpx < 0:
                rpx += Dx
            rpy = ry % Dy
            if rpy < 0:
                rpy += Dy

            for ox in range(-t, t + 1):
                nx = (rpx + ox) % Dx
                if nx < 0:
                    nx += Dx
                for oy in range(-t, t + 1):
                    ny = (rpy + oy) % Dy
                    if ny < 0:
                        ny += Dy
                    j = occ[nx * Dy + ny]
                    if j < 0 or vis[j]:
                        continue
                    vis[j] = 1
                    found += 1

                    delta = positions[j, 0] - rpx
                    if 2 * delta > Dx:
                        delta -= Dx
                    elif 2 * delta < -Dx:
                        delta += Dx
                    si_v[j, 0] = rx + delta

                    delta = positions[j, 1] - rpy
                    if 2 * delta > Dy:
                        delta -= Dy
                    elif 2 * delta < -Dy:
                        delta += Dy
                    si_v[j, 1] = ry + delta

                    q[tail] = j
                    tail += 1

    if found != N:
        raise ValueError(
            "Input positions must form a connected cluster within the provided space_threshold"
        )

    # shift so every coordinate is >= 0 (matches the Python reference)
    for d in range(n_dim):
        mn = si_v[0, d]
        for i in range(1, N):
            if si_v[i, d] < mn:
                mn = si_v[i, d]
        if mn < 0:
            for i in range(N):
                si_v[i, d] = si_v[i, d] - mn

    return si
