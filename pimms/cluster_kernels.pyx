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
from libc.math cimport pow

# Site-lookup selection. Filling a flat occupancy grid costs ~0.1 ns per site
# (a memset-like np.full plus N writes), but each of the walk's ~N * shell grid
# probes is a random read into a box-sized array (cache-missing in a big box)
# while a binary-search probe makes log2(N) comparisons in a small, cache-
# resident key array. Independent verifiers measured the crossover on an
# M-series laptop across 2D and 3D boxes of 12-300 sites per axis and
# thresholds 1-3: the search stays ahead until N reaches about
# 0.1 * volume^0.7 / shell^0.5 (e.g. 70 beads at t = 3 and ~300 at t = 1 in a
# 100^3 box, 850 at t = 3 and ~3000 at t = 1 in a 300^3 box; 20-100 in 2D),
# i.e. the grid wins once volume <= 26.5 * shell^0.72 * N^1.43. Earlier rules
# linear in N * shell were 2-14x slower than a plain grid on ordinary clusters
# or kept the search where the grid was 1.6x faster.
cdef double _GRID_FACTOR = 26.5
cdef double _GRID_SHELL_EXPONENT = 0.72
cdef double _GRID_EXPONENT = 1.43


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline Py_ssize_t _find_site(cnp.int64_t[::1] sorted_keys, cnp.int64_t[::1] sorted_idx,
                                  Py_ssize_t n, cnp.int64_t key) noexcept nogil:
    """Binary search ``key`` in ``sorted_keys``; return the bead index or -1.

    Parameters
    ----------
    sorted_keys : (N,) int64 memoryview
        Encoded site keys in ascending order, C-contiguous. A key packs a
        lattice position into a single integer as ``(x * Dy + y) * Dz + z`` in
        3D or ``x * Dy + y`` in 2D.
    sorted_idx : (N,) int64 memoryview
        Bead index that produced each entry of sorted_keys, i.e. the permutation
        that sorted the keys, C-contiguous.
    n : Py_ssize_t
        Number of entries to search over. This is the full bead count; sites are
        unique on an excluded-volume lattice, so there is one key per bead.
    key : int64 (cnp.int64_t)
        Encoded site key to look up.

    Returns
    -------
    Py_ssize_t
        Index of the bead sitting on that site, or -1 if the site is empty.
    """
    cdef Py_ssize_t lo = 0
    cdef Py_ssize_t hi = n - 1
    cdef Py_ssize_t mid
    while lo <= hi:
        mid = (lo + hi) >> 1
        if sorted_keys[mid] < key:
            lo = mid + 1
        elif sorted_keys[mid] > key:
            hi = mid - 1
        else:
            return sorted_idx[mid]
    return -1


cdef inline bint _use_grid(long volume, long n_beads, long shell) noexcept nogil:
    """Flat grid unless the box is large relative to the walk's work.

    Applies the measured crossover documented above: the grid wins while
    ``volume <= 26.5 * shell^0.72 * n_beads^1.43``.

    Parameters
    ----------
    volume : long
        Number of lattice sites in the box, i.e. the product of the box
        dimensions, and so the size of the grid that would have to be allocated
        and filled.
    n_beads : long
        Number of beads in the cluster being reconstructed.
    shell : long
        Number of sites probed around each bead, ``(2 * space_threshold + 1)``
        raised to the dimensionality.

    Returns
    -------
    bint
        True to use the flat occupancy grid, False to use the sorted key array
        with binary search.
    """
    return <double> volume <= (_GRID_FACTOR * pow(<double> shell, _GRID_SHELL_EXPONENT)
                               * pow(<double> n_beads, _GRID_EXPONENT))


def snakesearch_uses_grid(dims, n_beads, int space_threshold):
    """
    Report which site-lookup strategy :func:`snakesearch_single_image` takes.

    Parameters
    ----------
    dims : sequence of int
        Box dimensions (2 or 3 values). The length sets the dimensionality of
        the probe shell.
    n_beads : int
        Number of beads in the cluster that would be reconstructed.
    space_threshold : int (C int)
        Per-axis neighbour threshold of the gather, which sets how many sites
        are probed around each bead.

    Returns
    -------
    bool
        True if the flat box-sized occupancy grid is used, False if the sorted
        key array with binary search is used.
    """
    cdef long shell = 2 * space_threshold + 1
    cdef long volume = 1
    for d in dims:
        volume *= int(d)
    shell = shell * shell * (shell if len(dims) == 3 else 1)
    return bool(_use_grid(volume, int(n_beads), shell))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline long _mi_abs(long a, long D) noexcept nogil:
    """Absolute minimum-image separation from a raw coordinate difference.

    Uses exactly the wrap convention of :func:`_reach`, so the distance the
    connectivity predicate sees is the distance the placement uses.

    Parameters
    ----------
    a : long
        Raw per-axis difference between two in-box coordinates, so in (-D, D).
    D : long
        Box size on that axis.

    Returns
    -------
    long
        The absolute value of the minimum-image difference.
    """
    if 2 * a > D:
        a -= D
    elif 2 * a < -D:
        a += D
    return a if a >= 0 else -a


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline bint _linked(const cnp.int64_t* pos, cnp.int64_t[::1] typ,
                         cnp.int64_t[:, ::1] lr, cnp.int64_t[:, ::1] slr,
                         Py_ssize_t i, Py_ssize_t j, int n_dim,
                         long rpx, long rpy, long rpz,
                         long Dx, long Dy, long Dz) noexcept nogil:
    """Does the defining long-range cluster relation actually join beads i and j?

    A long-range cluster is the set of chains joined by a short-range contact
    (Chebyshev 1, any bead types), a Chebyshev-2 pair with a nonzero LR table
    entry, or a Chebyshev-3 pair with a nonzero SLR table entry. The gather that
    turns those beads into one periodic image has to walk that same relation:
    walking the looser "anything within Chebyshev 3" relation instead links the
    two extremes of a long cluster through the periodic face even when they do
    not interact, which throws part of the cluster a box-length away and tears
    it (Rg, hulls and radial profiles were then computed on the torn image).

    Parameters
    ----------
    pos : const int64 pointer
        First element of the (N, n_dim) in-box position array, row-major.
    typ : (N,) int64 memoryview
        Residue integer code of every bead, in the same order as ``pos``.
    lr : (n_res, n_res) int64 memoryview
        Long-range residue interaction table, indexed by those codes.
    slr : (n_res, n_res) int64 memoryview
        Super-long-range residue interaction table, indexed the same way.
    i : Py_ssize_t
        Index of the bead the walk is standing on.
    j : Py_ssize_t
        Index of the candidate neighbour.
    n_dim : int (C int)
        2 or 3; in 2D the z term is skipped.
    rpx : long
        In-box (wrapped) x of bead i.
    rpy : long
        In-box (wrapped) y of bead i.
    rpz : long
        In-box (wrapped) z of bead i (pass 0 in 2D).
    Dx : long
        Box size along x.
    Dy : long
        Box size along y.
    Dz : long
        Box size along z (pass 1 in 2D).

    Returns
    -------
    bint
        True if the pair carries one of the interactions that define long-range
        cluster membership, False otherwise.
    """
    cdef long d = _mi_abs(pos[j * n_dim] - rpx, Dx)
    cdef long e = _mi_abs(pos[j * n_dim + 1] - rpy, Dy)
    if e > d:
        d = e
    if n_dim == 3:
        e = _mi_abs(pos[j * n_dim + 2] - rpz, Dz)
        if e > d:
            d = e
    if d <= 1:
        return True
    if d == 2:
        return lr[typ[i], typ[j]] != 0
    if d == 3:
        return slr[typ[i], typ[j]] != 0
    return False


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _reach(const cnp.int64_t* pos, cnp.int64_t* si, cnp.uint8_t* vis,
                        cnp.int64_t* q, Py_ssize_t* tail, Py_ssize_t j, int n_dim,
                        long rx, long ry, long rz, long rpx, long rpy, long rpz,
                        long Dx, long Dy, long Dz) noexcept nogil:
    """Place newly reached bead ``j`` in the image nearest its parent and enqueue it
    (raw pointers into the C-contiguous arrays, so the call inlines).

    The parent's single-image coordinate is offset by the minimum-image
    separation between the two beads' in-box positions, which is what keeps the
    reconstructed cluster contiguous across the periodic boundary.

    Parameters
    ----------
    pos : const int64 pointer
        First element of the (N, n_dim) in-box position array, row-major, so
        bead j starts at ``pos[j * n_dim]``.
    si : int64 pointer
        First element of the (N, n_dim) single-image output array, row-major.
        The row for bead j is written here.
    vis : uint8 pointer
        First element of the (N,) visited-flag array. Bead j is marked visited.
    q : int64 pointer
        First element of the (N,) BFS queue. Bead j is appended.
    tail : Py_ssize_t pointer
        Current queue tail, read and incremented in place by the caller's
        walk loop.
    j : Py_ssize_t
        Index of the newly reached bead being placed and enqueued.
    n_dim : int (C int)
        2 or 3. In 2D the z branch is skipped and no third column is written.
    rx : long
        Single-image x of the parent bead the walk reached j from.
    ry : long
        Single-image y of the parent bead.
    rz : long
        Single-image z of the parent bead (pass 0 in 2D).
    rpx : long
        In-box (wrapped) x of the parent bead.
    rpy : long
        In-box (wrapped) y of the parent bead.
    rpz : long
        In-box (wrapped) z of the parent bead (pass 0 in 2D).
    Dx : long
        Box size along x, used for the minimum-image wrap.
    Dy : long
        Box size along y.
    Dz : long
        Box size along z (pass 1 in 2D).

    Returns
    -------
    None
        si, vis, q and tail are all modified through their pointers.
    """
    cdef long delta
    cdef Py_ssize_t base = j * n_dim
    vis[j] = 1
    delta = pos[base] - rpx
    if 2 * delta > Dx:
        delta -= Dx
    elif 2 * delta < -Dx:
        delta += Dx
    si[base] = rx + delta
    delta = pos[base + 1] - rpy
    if 2 * delta > Dy:
        delta -= Dy
    elif 2 * delta < -Dy:
        delta += Dy
    si[base + 1] = ry + delta
    if n_dim == 3:
        delta = pos[base + 2] - rpz
        if 2 * delta > Dz:
            delta -= Dz
        elif 2 * delta < -Dz:
            delta += Dz
        si[base + 2] = rz + delta
    q[tail[0]] = j
    tail[0] += 1


@cython.boundscheck(False)
@cython.wraparound(False)
def snakesearch_single_image(cnp.int64_t[:, ::1] positions,
                             cnp.int64_t[::1] dims,
                             Py_ssize_t seed_idx,
                             int space_threshold,
                             types=None,
                             LR_table=None,
                             SLR_table=None):
    """
    Breadth-first "snakesearch" reconstruction of a cluster into a single periodic
    image.

    This is the compiled equivalent of the pure-Python
    ``cluster_utils.convert_positions_to_single_image_snakesearch`` and produces
    byte-for-byte identical output for a given ``seed_idx`` on valid
    (singly-occupied) clusters - the only kind the excluded-volume lattice can
    produce. (With two or more beads stacked on ONE site the site lookup keeps
    only one bead per site - which one is unspecified and differs between the
    two lookup modes - so the kernel may report the cluster as disconnected
    where the dict-based fallback tolerates the duplicates.)
    Starting from the seed bead, each unvisited neighbour within
    ``space_threshold`` (per axis, PBC-aware) is placed into the same periodic
    image as the reference bead and enqueued; the result is finally shifted so
    every coordinate is >= 0.

    Neighbour discovery runs in typed C with integer arithmetic. Site lookups use
    a flat occupancy-index grid of size ``prod(dims)`` (O(1) per lookup) unless
    the box is large relative to the cluster (``volume > 26.5 * shell^0.72 *
    N^1.43``, the measured crossover), in which case a sorted array of encoded site keys
    with binary search is used instead, so a handful of beads in a very large
    box no longer pays to allocate and fill a box-sized grid on every call.
    :func:`snakesearch_uses_grid` reports the choice.

    Parameters
    ----------
    positions : (N, n_dim) int64 memoryview
        Bead positions in PBC space, C-contiguous, with each coordinate in
        ``[0, dims[d])``. n_dim must be 2 or 3.
    dims : (n_dim,) int64 memoryview
        Box dimensions per axis, C-contiguous.
    seed_idx : Py_ssize_t
        Index of the bead to start the walk from (its image is preserved).
    space_threshold : int (C int)
        Maximum per-axis distance for two beads to count as neighbours.
    types : (N,) int64 array or None, optional
        Residue integer codes for the beads, in the same order as
        ``positions``. Supplying these (together with both tables) switches the
        walk from the plain distance rule to the long-range cluster membership
        rule - Chebyshev 1, or Chebyshev 2 / 3 with a nonzero LR / SLR entry -
        so the gather walks the relation that defines the cluster. Default is
        None (plain distance rule).
    LR_table : (n_res, n_res) int64 array or None, optional
        Long-range residue interaction table, required when ``types`` is given.
        Default is None.
    SLR_table : (n_res, n_res) int64 array or None, optional
        Super-long-range residue interaction table, required when ``types`` is
        given. Default is None.

    Returns
    -------
    numpy.ndarray
        ``(N, n_dim)`` int64 array of single-image positions, shifted so every
        coordinate is >= 0.

    Raises
    ------
    ValueError
        If ``positions`` is not 2D/3D, if any coordinate lies outside the box, if
        ``types`` is supplied without both tables or with the wrong length, or if the
        beads do not form a single connected cluster within ``space_threshold`` (i.e. not
        every bead was reached).
    """
    cdef Py_ssize_t N = positions.shape[0]
    cdef int n_dim = positions.shape[1]
    cdef Py_ssize_t i, j, ref, head, tail, found
    cdef int t = space_threshold
    cdef int d
    cdef long Dx = 1
    cdef long Dy = 1
    cdef long Dz = 1
    cdef long ox, oy, oz
    cdef long rx, ry, rz, rpx, rpy, rpz, nx, ny, nz, mn
    cdef cnp.int64_t key
    cdef bint use_grid
    cdef long shell
    cdef bint use_pred
    cdef cnp.int64_t[::1] typ_v
    cdef cnp.int64_t[:, ::1] lr_v
    cdef cnp.int64_t[:, ::1] slr_v

    if n_dim != dims.shape[0] or n_dim not in (2, 3):
        raise ValueError(
            "snakesearch_single_image: positions must be (N, 2) or (N, 3) and match dims")
    # read the box only after the dimensionality check: with bounds checking off,
    # dims[1] on a length-1 buffer is an out-of-range read
    Dx = dims[0]
    Dy = dims[1]
    if n_dim == 3:
        Dz = dims[2]

    # seed_idx indexes si_v/vis/q with bounds checking off, so an out-of-range value
    # (or an empty position set, where the queue buffer has zero length) would write
    # past the end of the buffers rather than raise
    if N == 0 or seed_idx < 0 or seed_idx >= N:
        raise ValueError(
            "snakesearch_single_image: seed_idx %d out of range for %d positions" % (seed_idx, N))

    # The occupancy grid below is indexed by raw position with bounds checking off, so an
    # out-of-box coordinate would read/write past the end of the buffer rather than raise.
    # This O(N) guard is negligible next to the BFS and turns silent memory corruption into
    # a clear error. (The pure-Python fallback in cluster_utils uses a dict and so tolerates
    # arbitrary coordinates; callers must pass in-box positions to reach this kernel.)
    for i in range(N):
        for d in range(n_dim):
            if positions[i, d] < 0 or positions[i, d] >= dims[d]:
                raise ValueError(
                    "snakesearch_single_image: position %d is outside the box in axis %d "
                    "(positions must already be wrapped into [0, dims))" % (i, d))

    # optional connectivity predicate. The memoryviews are always bound (to a
    # dummy when unused) because an unbound memoryview would be an unchecked
    # read inside the nogil walk.
    if types is None:
        use_pred = 0
        typ_v = np.zeros(1, dtype=np.int64)
        lr_v = np.zeros((1, 1), dtype=np.int64)
        slr_v = np.zeros((1, 1), dtype=np.int64)
    else:
        if LR_table is None or SLR_table is None:
            raise ValueError(
                "snakesearch_single_image: types requires both LR_table and SLR_table")
        typ_v = np.ascontiguousarray(np.asarray(types, dtype=np.int64))
        lr_v = np.ascontiguousarray(np.asarray(LR_table, dtype=np.int64))
        slr_v = np.ascontiguousarray(np.asarray(SLR_table, dtype=np.int64))
        if typ_v.shape[0] != N:
            raise ValueError(
                "snakesearch_single_image: types must have one entry per position")
        use_pred = 1

    si = np.empty((N, n_dim), dtype=np.int64)
    cdef cnp.int64_t[:, ::1] si_v = si

    # array-backed BFS queue + visited flags
    cdef cnp.int64_t[::1] q = np.empty(N, dtype=np.int64)
    cdef cnp.uint8_t[::1] vis = np.zeros(N, dtype=np.uint8)

    shell = 2 * t + 1
    shell = shell * shell * (shell if n_dim == 3 else 1)
    use_grid = _use_grid(Dx * Dy * Dz, N, shell)

    # site lookup structures: exactly one of the two is built
    cdef cnp.int64_t[::1] occ
    cdef cnp.int64_t[::1] sorted_keys
    cdef cnp.int64_t[::1] sorted_idx
    cdef cnp.int64_t[::1] keys_v
    if use_grid:
        occ = np.full(Dx * Dy * Dz, -1, dtype=np.int64)
        if n_dim == 3:
            for i in range(N):
                occ[(positions[i, 0] * Dy + positions[i, 1]) * Dz + positions[i, 2]] = i
        else:
            for i in range(N):
                occ[positions[i, 0] * Dy + positions[i, 1]] = i
    else:
        keys = np.empty(N, dtype=np.int64)
        keys_v = keys
        if n_dim == 3:
            for i in range(N):
                keys_v[i] = (positions[i, 0] * Dy + positions[i, 1]) * Dz + positions[i, 2]
        else:
            for i in range(N):
                keys_v[i] = positions[i, 0] * Dy + positions[i, 1]
        order = np.argsort(keys, kind='stable')
        sorted_keys = keys[order]
        sorted_idx = order.astype(np.int64)

    for d in range(n_dim):
        si_v[seed_idx, d] = positions[seed_idx, d]
    vis[seed_idx] = 1
    head = 0
    tail = 0
    q[tail] = seed_idx
    tail += 1
    rz = 0
    rpz = 0

    # raw pointers for the placement helper (all arrays are C-contiguous)
    cdef const cnp.int64_t* pos_p = &positions[0, 0]
    cdef cnp.int64_t* si_p = &si_v[0, 0]
    cdef cnp.uint8_t* vis_p = &vis[0]
    cdef cnp.int64_t* q_p = &q[0]

    # One fully specialised walk per (dimensionality, lookup mode): the mode is
    # loop-invariant and a branch anywhere inside the walk measurably slows the
    # grid path (8-16 % on 2D clusters of a few hundred beads).
    if n_dim == 3 and use_grid:
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
                        # use_pred is loop-invariant, so this costs a perfectly
                        # predicted branch on the (default) distance-only path
                        if use_pred and not _linked(pos_p, typ_v, lr_v, slr_v, ref, j, 3,
                                                    rpx, rpy, rpz, Dx, Dy, Dz):
                            continue
                        _reach(pos_p, si_p, vis_p, q_p, &tail, j, 3,
                               rx, ry, rz, rpx, rpy, rpz, Dx, Dy, Dz)
    elif n_dim == 3:
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
                        j = _find_site(sorted_keys, sorted_idx, N, (nx * Dy + ny) * Dz + nz)
                        if j < 0 or vis[j]:
                            continue
                        if use_pred and not _linked(pos_p, typ_v, lr_v, slr_v, ref, j, 3,
                                                    rpx, rpy, rpz, Dx, Dy, Dz):
                            continue
                        _reach(pos_p, si_p, vis_p, q_p, &tail, j, 3,
                               rx, ry, rz, rpx, rpy, rpz, Dx, Dy, Dz)
    elif use_grid:
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
                    if use_pred and not _linked(pos_p, typ_v, lr_v, slr_v, ref, j, 2,
                                                rpx, rpy, 0, Dx, Dy, 1):
                        continue
                    _reach(pos_p, si_p, vis_p, q_p, &tail, j, 2,
                           rx, ry, 0, rpx, rpy, 0, Dx, Dy, 1)
    else:
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
                    j = _find_site(sorted_keys, sorted_idx, N, nx * Dy + ny)
                    if j < 0 or vis[j]:
                        continue
                    if use_pred and not _linked(pos_p, typ_v, lr_v, slr_v, ref, j, 2,
                                                rpx, rpy, 0, Dx, Dy, 1):
                        continue
                    _reach(pos_p, si_p, vis_p, q_p, &tail, j, 2,
                           rx, ry, 0, rpx, rpy, 0, Dx, Dy, 1)

    found = tail

    if found != N:
        raise ValueError(
            "Input positions must form a connected cluster within the provided space_threshold"
        )

    # shift so every coordinate is >= 0, by a WHOLE NUMBER OF BOX PERIODS
    # (matches the Python reference; congruence mod the box must be preserved -
    # see cluster_utils.convert_positions_to_single_image_snakesearch)
    cdef long shift
    for d in range(n_dim):
        mn = si_v[0, d]
        for i in range(1, N):
            if si_v[i, d] < mn:
                mn = si_v[i, d]
        if mn < 0:
            shift = dims[d] * ((-mn + dims[d] - 1) // dims[d])
            for i in range(N):
                si_v[i, d] = si_v[i, d] + shift

    return si
