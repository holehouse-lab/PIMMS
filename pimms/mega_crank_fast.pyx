## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
##
## pimms/mega_crank_fast.pyx
##
## Drop-in optimized replacement for pimms.mega_crank.mega_crank().
##
## GOAL
## ----
## The reference kernel (pimms/mega_crank.pyx) is correct but allocates several
## numpy arrays *inside* the per-substep hot loop:
##
##   * crank_it / update_position / single_bead_crank each `np.zeros([3])`
##   * get_angle_energy_change allocates 4 arrays (a, b, angle_positions,
##     intcode_lookup) on EVERY call
##
## With CRANKSHAFT_SUBSTEPS ~ 10^4 and thousands of megamoves this is hundreds
## of millions of heap allocations routed through the Python/numpy layer - which
## profiling showed to be ~85% of total runtime.
##
## This module performs the IDENTICAL algorithm with IDENTICAL random-number
## call ordering, but every transient buffer is a C stack array. Because the
## RNG sequence is preserved bit-for-bit, given the same (grid, type_grid,
## idx_to_bead, seed, ...) inputs this kernel returns the SAME (energy,
## accepted_moves) and leaves the SAME mutated grid/type_grid/idx_to_bead as the
## reference kernel. That makes it verifiable by exact comparison (see
## pimms/fast_kernels/benchmark.py).
##
## The public signature matches pimms.mega_crank.mega_crank exactly, so it can
## be swapped in at moves.py:184 with no other changes.
## ...........................................................................

import numpy as np
cimport numpy as cnp
cnp.import_array()

cimport cython
from libc.math cimport exp
from libc.stdlib cimport malloc, free
from cython.parallel cimport prange

from pimms.cython_config cimport NUMPY_INT_TYPE
from pimms.cython_config cimport NUMPY_INT_TYPE_long


# ---- Monte-Carlo PRNG (splitmix64) --------------------------------------
# Replaces libc rand()/srand(). On macOS rand() is the weak Park-Miller MINSTD
# LCG (x -> 16807*x mod 2^31-1: short period ~2.1e9 and pronounced lattice
# structure) and is platform-dependent (glibc differs, Windows RAND_MAX=32767).
# splitmix64 has period 2^64, passes BigCrush, and is bit-identical on every
# platform. The serial kernels are single-threaded, so a module-global state is
# fine; the parallel checkerboard kernel keeps its own per-block PRNG.
cdef unsigned long long _RNG_STATE[1]   # module-global serial PRNG state (zero-init)
cdef int PRNG_MAX = 2147483647          # 2^31 - 1 (fixed, platform-independent)

cdef inline void mc_seed(unsigned int seedval) noexcept nogil:
    _RNG_STATE[0] = <unsigned long long>seedval

cdef inline int mc_rand() noexcept nogil:
    # one splitmix64 step; return the top 31 bits -> [0, PRNG_MAX]
    _RNG_STATE[0] = _RNG_STATE[0] + <unsigned long long>0x9E3779B97F4A7C15
    cdef unsigned long long z = _RNG_STATE[0]
    z = (z ^ (z >> 30)) * <unsigned long long>0xBF58476D1CE4E5B9
    z = (z ^ (z >> 27)) * <unsigned long long>0x94D049BB133111EB
    z = z ^ (z >> 31)
    return <int>(z >> 33)


cdef inline int int_max(int a, int b) noexcept nogil: return a if a >= b else b
cdef inline int int_min(int a, int b) noexcept nogil: return a if a <= b else b


# -----------------------------------------------------------------
# Periodic boundary correction. Identical to the reference kernel:
# avoids the modulo on the common non-negative path.
@cython.cdivision(True)
cdef inline int pbc_correction(int value, int DIM) noexcept nogil:
    if value < 0:
        return DIM + value
    else:
        return value % DIM


# -----------------------------------------------------------------
# Integer RNG. Byte-for-byte identical to pimms.mega_crank.randint so the
# random stream (and therefore every accept/reject decision) matches exactly.
cdef inline int randint(int start, int end) noexcept nogil:
    if start == 0:
        end = end + 1
    # NOTE: must use double precision here. The reference kernel writes
    # float(rand()-1)/float(RAND_MAX), and Cython's float() is a C *double*
    # (Python float). Using a 32-bit C float (<float>) rounds differently at
    # truncation boundaries and produces different draws for the same RNG state.
    cdef int r = start + <int>((<double>(mc_rand() - 1) / <double>PRNG_MAX) * (end))
    return r


# -----------------------------------------------------------------
# Metropolis criterion. Identical to the reference: note rand() is consumed
# ONLY when new_energy > old_energy, which must be preserved to keep the
# stream aligned.
cdef inline int accept_or_reject(float invtemp, long old_energy, long new_energy) noexcept nogil:
    cdef float expterm
    cdef float randval

    if new_energy <= old_energy:
        return 1

    expterm = exp(-(new_energy - old_energy) * invtemp)
    # double-precision division to match the reference's float(rand())/float(RAND_MAX)
    randval = <double>mc_rand() / <double>PRNG_MAX

    if randval < expterm:
        return 1
    else:
        return 0


@cython.cdivision(True)
cdef inline int fix_angle_pbc_issues(int distance) noexcept nogil:
    if distance < -1:
        return 1
    elif distance > 1:
        return -1
    else:
        return distance


# -----------------------------------------------------------------
# Straddle check between two lattice points (used for hardwall rejection).
# Mirrors do_positions_stradle_pbc_boundary's per-consecutive-pair test.
cdef inline int straddle_pair(int ax, int ay, int az,
                              int bx, int by, int bz) noexcept nogil:
    if abs(ax - bx) > 1:
        return 1
    if abs(ay - by) > 1:
        return 1
    if abs(az - bz) > 1:
        return 1
    return 0


# -----------------------------------------------------------------
# Single-bead / terminal-bead proposal. Reproduces single_bead_crank ->
# update_position: three randint() draws (x, y, z order) then a hardsphere
# test. Writes the proposed position into out[3]; returns 1 on success, 0 if
# the target site is occupied (out[0] set to -1, matching the reference).
@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int single_bead_crank_c(int ox, int oy, int oz,
                                    NUMPY_INT_TYPE[:, :, :] grid,
                                    int XDIM, int YDIM, int ZDIM,
                                    int* out) noexcept nogil:
    cdef int x_off = randint(0, 2) - 1
    cdef int y_off = randint(0, 2) - 1
    cdef int z_off = randint(0, 2) - 1

    cdef int local_x = pbc_correction(ox + x_off, XDIM)
    cdef int local_y = pbc_correction(oy + y_off, YDIM)
    cdef int local_z = pbc_correction(oz + z_off, ZDIM)

    if grid[local_x, local_y, local_z] > 0:
        out[0] = -1
        out[1] = 0
        out[2] = 0
        return 0

    out[0] = local_x
    out[1] = local_y
    out[2] = local_z
    return 1


# -----------------------------------------------------------------
# Internal-bead crankshaft proposal. Reproduces crank_it exactly: the
# de-periodisation of the N/C anchors, the constrained box, then three
# randint() draws (x, y, z order). Only the N-side (index 0) and C-side
# (index 2) of the triptic are needed - the moved bead itself is unused.
@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int crank_it_c(int N_side_x, int N_side_y, int N_side_z,
                           int C_side_x, int C_side_y, int C_side_z,
                           NUMPY_INT_TYPE[:, :, :] grid,
                           int XDIM, int YDIM, int ZDIM,
                           int* out) noexcept nogil:
    cdef int x_min, x_max, y_min, y_max, z_min, z_max
    cdef int local_x, local_y, local_z

    # de-periodise: if anchors are >2 apart in a dim they straddle the PBC,
    # so lift the smaller into the 'fake' extended space.
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

    x_min = int_max(N_side_x - 1, C_side_x - 1)
    x_max = int_min(N_side_x + 1, C_side_x + 1)
    y_min = int_max(N_side_y - 1, C_side_y - 1)
    y_max = int_min(N_side_y + 1, C_side_y + 1)
    z_min = int_max(N_side_z - 1, C_side_z - 1)
    z_max = int_min(N_side_z + 1, C_side_z + 1)

    local_x = pbc_correction((x_min + randint(1, (x_max - x_min + 1)) - 1), XDIM)
    local_y = pbc_correction((y_min + randint(1, (y_max - y_min + 1)) - 1), YDIM)
    local_z = pbc_correction((z_min + randint(1, (z_max - z_min + 1)) - 1), ZDIM)

    if grid[local_x, local_y, local_z] > 0:
        out[0] = -1
        out[1] = 0
        out[2] = 0
        return 0

    out[0] = local_x
    out[1] = local_y
    out[2] = local_z
    return 1


# -----------------------------------------------------------------
# Angle-energy delta. Identical maths to get_angle_energy_change but with C
# stack arrays in place of the four per-call numpy allocations.
@cython.wraparound(False)
@cython.boundscheck(False)
cdef long get_angle_energy_change_c(int bead_index,
                                    NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                                    int* new_position,
                                    NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup) noexcept nogil:

    if idx_to_bead[bead_index, 3] == 1:
        return 0

    cdef int a0, a1, a2, b0, b1, b2
    cdef long angle_penalty_new = 0
    cdef long angle_penalty_old = 0

    cdef int angle_positions[5][3]
    cdef int intcode_lookup[5]

    cdef int offset, offset_start, offset_end, i, pos
    cdef int angle_idx
    cdef int local_move_idx = -1

    if idx_to_bead[bead_index, 0] == 1:
        offset_start = 0
        offset_end = 3
    elif idx_to_bead[bead_index, 0] == 2:
        offset_start = -2
        offset_end = 3
    elif idx_to_bead[bead_index, 0] == 4:
        offset_start = -1
        offset_end = 2
    elif idx_to_bead[bead_index, 0] == 3:
        offset_start = -2
        offset_end = 1
    elif idx_to_bead[bead_index, 0] == 5:
        offset_start = -1
        offset_end = 3
    elif idx_to_bead[bead_index, 0] == 6:
        offset_start = -2
        offset_end = 2
    else:
        # single bead (flag 0) never reaches here in the reference because the
        # skip-angle flag short-circuits; keep a safe default.
        offset_start = 0
        offset_end = 0

    angle_idx = 0
    for offset in range(offset_start, offset_end):
        pos = bead_index + offset
        angle_positions[angle_idx][0] = idx_to_bead[pos, 5]
        angle_positions[angle_idx][1] = idx_to_bead[pos, 6]
        angle_positions[angle_idx][2] = idx_to_bead[pos, 7]
        intcode_lookup[angle_idx] = idx_to_bead[pos, 2]
        if pos == bead_index:
            local_move_idx = angle_idx
        angle_idx = angle_idx + 1

    # old configuration
    for i in range(0, angle_idx - 2):
        a0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i][0])
        a1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i][1])
        a2 = fix_angle_pbc_issues(angle_positions[i + 1][2] - angle_positions[i][2])

        b0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i + 2][0])
        b1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i + 2][1])
        b2 = fix_angle_pbc_issues(angle_positions[i + 1][2] - angle_positions[i + 2][2])

        angle_penalty_old = angle_lookup[intcode_lookup[i + 1], a0 + 1, a1 + 1, a2 + 1, b0 + 1, b1 + 1, b2 + 1] + angle_penalty_old

    # move the central bead to the proposed position
    if local_move_idx >= 0:
        angle_positions[local_move_idx][0] = new_position[0]
        angle_positions[local_move_idx][1] = new_position[1]
        angle_positions[local_move_idx][2] = new_position[2]

    # new configuration
    for i in range(0, angle_idx - 2):
        a0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i][0])
        a1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i][1])
        a2 = fix_angle_pbc_issues(angle_positions[i + 1][2] - angle_positions[i][2])

        b0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i + 2][0])
        b1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i + 2][1])
        b2 = fix_angle_pbc_issues(angle_positions[i + 1][2] - angle_positions[i + 2][2])

        angle_penalty_new = angle_lookup[intcode_lookup[i + 1], a0 + 1, a1 + 1, a2 + 1, b0 + 1, b1 + 1, b2 + 1] + angle_penalty_new

    return angle_penalty_new - angle_penalty_old


# -----------------------------------------------------------------
# Short/long-range interaction-energy delta. Same logic as get_energy_change,
# operating on plain ints rather than numpy position arrays. type_grid is
# temporarily mutated and restored exactly as in the reference.
@cython.wraparound(False)
@cython.boundscheck(False)
cdef long get_energy_change_c(NUMPY_INT_TYPE[:, :, :] grid,
                              NUMPY_INT_TYPE[:, :, :] type_grid,
                              int old_x, int old_y, int old_z,
                              int new_x, int new_y, int new_z,
                              int LR_vs_SR,
                              NUMPY_INT_TYPE[:, :] interaction_table,
                              NUMPY_INT_TYPE[:, :] LR_interaction_table,
                              NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                              int XDIM, int YDIM, int ZDIM,
                              int hardwall) noexcept nogil:

    cdef long energy_old = 0
    cdef long energy_old_empty = 0
    cdef long energy_new = 0
    cdef long energy_new_empty = 0

    cdef int tmp_x, tmp_y, tmp_z
    cdef int x, y, z
    cdef unsigned int site_bead_type, bead_type

    bead_type = type_grid[old_x, old_y, old_z]

    if LR_vs_SR == 1:
        # ---- old energy over the 7x7x7 neighbourhood ----
        for x in range(-3, 4):
            for y in range(-3, 4):
                for z in range(-3, 4):
                    tmp_x = pbc_correction(old_x + x, XDIM)
                    tmp_y = pbc_correction(old_y + y, YDIM)
                    tmp_z = pbc_correction(old_z + z, ZDIM)
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]

                    if abs(x) < 2 and abs(y) < 2 and abs(z) < 2:
                        if hardwall == 1:
                            if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                                site_bead_type = 0
                        energy_old = energy_old + interaction_table[bead_type, site_bead_type]
                        if site_bead_type > 0:
                            energy_old_empty = energy_old_empty + interaction_table[0, site_bead_type]
                    elif abs(x) < 3 and abs(y) < 3 and abs(z) < 3:
                        if hardwall == 1:
                            if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                                continue
                        energy_old = energy_old + LR_interaction_table[bead_type, site_bead_type]
                    else:
                        if hardwall == 1:
                            if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                                continue
                        energy_old = energy_old + SLR_interaction_table[bead_type, site_bead_type]

        type_grid[new_x, new_y, new_z] = type_grid[old_x, old_y, old_z]
        type_grid[old_x, old_y, old_z] = 0

        # ---- new energy ----
        for x in range(-3, 4):
            for y in range(-3, 4):
                for z in range(-3, 4):
                    tmp_x = pbc_correction(new_x + x, XDIM)
                    tmp_y = pbc_correction(new_y + y, YDIM)
                    tmp_z = pbc_correction(new_z + z, ZDIM)
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]

                    if abs(x) < 2 and abs(y) < 2 and abs(z) < 2:
                        if hardwall == 1:
                            if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                                site_bead_type = 0
                        energy_new = energy_new + interaction_table[bead_type, site_bead_type]
                        if site_bead_type > 0:
                            energy_new_empty = energy_new_empty + interaction_table[0, site_bead_type]
                    elif abs(x) < 3 and abs(y) < 3 and abs(z) < 3:
                        if hardwall == 1:
                            if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                                continue
                        energy_new = energy_new + LR_interaction_table[bead_type, site_bead_type]
                    else:
                        if hardwall == 1:
                            if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                                continue
                        energy_new = energy_new + SLR_interaction_table[bead_type, site_bead_type]

    else:
        # ---- short-range only: 3x3x3 ----
        for x in range(-1, 2):
            for y in range(-1, 2):
                for z in range(-1, 2):
                    tmp_x = pbc_correction(old_x + x, XDIM)
                    tmp_y = pbc_correction(old_y + y, YDIM)
                    tmp_z = pbc_correction(old_z + z, ZDIM)
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]

                    if hardwall == 1:
                        if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3 or abs(tmp_z - old_z) > 3:
                            site_bead_type = 0

                    energy_old = energy_old + interaction_table[bead_type, site_bead_type]
                    if site_bead_type > 0:
                        energy_old_empty = energy_old_empty + interaction_table[0, site_bead_type]

        type_grid[new_x, new_y, new_z] = type_grid[old_x, old_y, old_z]
        type_grid[old_x, old_y, old_z] = 0

        for x in range(-1, 2):
            for y in range(-1, 2):
                for z in range(-1, 2):
                    tmp_x = pbc_correction(new_x + x, XDIM)
                    tmp_y = pbc_correction(new_y + y, YDIM)
                    tmp_z = pbc_correction(new_z + z, ZDIM)
                    site_bead_type = type_grid[tmp_x, tmp_y, tmp_z]
                    if hardwall == 1:
                        if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3 or abs(tmp_z - new_z) > 3:
                            site_bead_type = 0

                    energy_new = energy_new + interaction_table[bead_type, site_bead_type]
                    if site_bead_type > 0:
                        energy_new_empty = energy_new_empty + interaction_table[0, site_bead_type]

    # remove the self-interaction double counting at x==y==z==0
    energy_old_empty = energy_old_empty - interaction_table[0, bead_type]
    energy_old = energy_old - interaction_table[bead_type, bead_type]
    energy_new_empty = energy_new_empty - interaction_table[0, bead_type]
    energy_new = energy_new - interaction_table[bead_type, bead_type]

    # restore the type grid
    type_grid[old_x, old_y, old_z] = type_grid[new_x, new_y, new_z]
    type_grid[new_x, new_y, new_z] = 0

    return (energy_new + energy_old_empty) - (energy_old + energy_new_empty)


# -----------------------------------------------------------------
# Public entry point - signature identical to pimms.mega_crank.mega_crank.
@cython.wraparound(False)
@cython.boundscheck(False)
def mega_crank(NUMPY_INT_TYPE[:, :, :] grid,
               NUMPY_INT_TYPE[:, :, :] type_grid,
               NUMPY_INT_TYPE_long[:, :] idx_to_bead,
               NUMPY_INT_TYPE[:, :] interaction_table,
               NUMPY_INT_TYPE[:, :] LR_interaction_table,
               NUMPY_INT_TYPE[:, :] SLR_interaction_table,
               NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup,
               long energy,
               float invtemp,
               int nsteps,
               NUMPY_INT_TYPE_long[:] bead_selector,
               int passed_seed,
               int hardwall):
    """
    Optimized, allocation-free crankshaft Monte Carlo kernel.

    Behaviourally identical (same RNG stream) to pimms.mega_crank.mega_crank;
    see module docstring. Returns (energy, accepted_moves).
    """
    mc_seed(passed_seed)

    cdef unsigned int i
    cdef int bead_index
    cdef int accepted_moves = 0
    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int ZDIM = grid.shape[2]
    cdef long delta_energy, delta_angle_energy
    cdef int beadflag
    cdef int ok

    # reusable C buffers (no per-substep heap allocation)
    cdef int old_position[3]
    cdef int new_position[3]
    cdef int anchor_bead[3]
    cdef int Nx, Ny, Nz, Cx, Cy, Cz

    for i in range(nsteps):
        bead_index = bead_selector[i]

        old_position[0] = idx_to_bead[bead_index, 5]
        old_position[1] = idx_to_bead[bead_index, 6]
        old_position[2] = idx_to_bead[bead_index, 7]

        beadflag = idx_to_bead[bead_index, 0]

        # ---- single bead (flag 0) ----
        if beadflag == 0:
            single_bead_crank_c(old_position[0], old_position[1], old_position[2],
                                grid, XDIM, YDIM, ZDIM, new_position)

        # ---- N-terminal bead (flag 1) ----
        elif beadflag == 1:
            anchor_bead[0] = idx_to_bead[bead_index + 1, 5]
            anchor_bead[1] = idx_to_bead[bead_index + 1, 6]
            anchor_bead[2] = idx_to_bead[bead_index + 1, 7]

            single_bead_crank_c(anchor_bead[0], anchor_bead[1], anchor_bead[2],
                                grid, XDIM, YDIM, ZDIM, new_position)

            if hardwall == 1:
                if straddle_pair(new_position[0], new_position[1], new_position[2],
                                 anchor_bead[0], anchor_bead[1], anchor_bead[2]) == 1:
                    continue

        # ---- C-terminal bead (flag 3) ----
        elif beadflag == 3:
            anchor_bead[0] = idx_to_bead[bead_index - 1, 5]
            anchor_bead[1] = idx_to_bead[bead_index - 1, 6]
            anchor_bead[2] = idx_to_bead[bead_index - 1, 7]

            single_bead_crank_c(anchor_bead[0], anchor_bead[1], anchor_bead[2],
                                grid, XDIM, YDIM, ZDIM, new_position)

            if hardwall == 1:
                if straddle_pair(anchor_bead[0], anchor_bead[1], anchor_bead[2],
                                 new_position[0], new_position[1], new_position[2]) == 1:
                    continue

        # ---- internal bead (flags 2,4,5,6) ----
        else:
            Nx = idx_to_bead[bead_index - 1, 5]
            Ny = idx_to_bead[bead_index - 1, 6]
            Nz = idx_to_bead[bead_index - 1, 7]
            Cx = idx_to_bead[bead_index + 1, 5]
            Cy = idx_to_bead[bead_index + 1, 6]
            Cz = idx_to_bead[bead_index + 1, 7]

            crank_it_c(Nx, Ny, Nz, Cx, Cy, Cz, grid, XDIM, YDIM, ZDIM, new_position)

            if hardwall == 1:
                # three-position holder: [N, new, C]
                if (straddle_pair(Nx, Ny, Nz,
                                  new_position[0], new_position[1], new_position[2]) == 1
                    or straddle_pair(new_position[0], new_position[1], new_position[2],
                                     Cx, Cy, Cz) == 1):
                    continue

        # ---- accept/reject if the hardsphere proposal succeeded ----
        if not new_position[0] < 0:
            delta_energy = get_energy_change_c(grid, type_grid,
                                               old_position[0], old_position[1], old_position[2],
                                               new_position[0], new_position[1], new_position[2],
                                               idx_to_bead[bead_index, 1],
                                               interaction_table, LR_interaction_table,
                                               SLR_interaction_table,
                                               XDIM, YDIM, ZDIM, hardwall)
            delta_angle_energy = get_angle_energy_change_c(bead_index, idx_to_bead,
                                                           new_position, angle_lookup)

            if accept_or_reject(invtemp, energy, energy + delta_energy + delta_angle_energy) == 1:
                grid[old_position[0], old_position[1], old_position[2]] = 0
                grid[new_position[0], new_position[1], new_position[2]] = idx_to_bead[bead_index, 4]

                type_grid[new_position[0], new_position[1], new_position[2]] = type_grid[old_position[0], old_position[1], old_position[2]]
                type_grid[old_position[0], old_position[1], old_position[2]] = 0

                idx_to_bead[bead_index, 5] = new_position[0]
                idx_to_bead[bead_index, 6] = new_position[1]
                idx_to_bead[bead_index, 7] = new_position[2]

                energy = energy + delta_energy + delta_angle_energy
                accepted_moves = accepted_moves + 1

    return (energy, accepted_moves)


# =====================================================================
#
#   2D crankshaft kernel (allocation-free, bit-exact drop-in for
#   pimms.mega_crank_2D.mega_crank_2D)
#
#   The reference 2D kernel calls mega_crank.randint_ext /
#   accept_or_reject_ext and the Python-level (def) get_energy_change_2D on
#   every substep - each crossing the Python/C boundary - and allocates the
#   angle arrays per call. This version uses the same nogil C helpers as the 3D
#   fast kernel (randint / accept_or_reject - already double precision, matching
#   the reference) plus 2D C-array energy/angle/proposal functions, so it is bit
#   identical but avoids all that overhead.
# =====================================================================

# ---- 2D single-bead / terminal proposal ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int single_bead_crank_2D_c(int ox, int oy,
                                       NUMPY_INT_TYPE[:, :] grid,
                                       int XDIM, int YDIM,
                                       int* out) noexcept nogil:
    cdef int x_off = randint(0, 2) - 1
    cdef int y_off = randint(0, 2) - 1
    cdef int local_x = pbc_correction(ox + x_off, XDIM)
    cdef int local_y = pbc_correction(oy + y_off, YDIM)
    if grid[local_x, local_y] > 0:
        out[0] = -1
        out[1] = 0
        return 0
    out[0] = local_x
    out[1] = local_y
    return 1


# ---- 2D internal-bead crankshaft proposal ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int crank_it_2D_c(int N_side_x, int N_side_y,
                              int C_side_x, int C_side_y,
                              NUMPY_INT_TYPE[:, :] grid,
                              int XDIM, int YDIM,
                              int* out) noexcept nogil:
    cdef int x_min, x_max, y_min, y_max, local_x, local_y

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

    x_min = int_max(N_side_x - 1, C_side_x - 1)
    x_max = int_min(N_side_x + 1, C_side_x + 1)
    y_min = int_max(N_side_y - 1, C_side_y - 1)
    y_max = int_min(N_side_y + 1, C_side_y + 1)

    local_x = pbc_correction((x_min + randint(1, (x_max - x_min + 1)) - 1), XDIM)
    local_y = pbc_correction((y_min + randint(1, (y_max - y_min + 1)) - 1), YDIM)

    if grid[local_x, local_y] > 0:
        out[0] = -1
        out[1] = 0
        return 0
    out[0] = local_x
    out[1] = local_y
    return 1


# ---- 2D angle-energy delta (angle_lookup is 5D in 2D) ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef long get_angle_energy_change_2D_c(int bead_index,
                                       NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                                       int* new_position,
                                       NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup) noexcept nogil:

    if idx_to_bead[bead_index, 3] == 1:
        return 0

    cdef int a0, a1, b0, b1
    cdef long angle_penalty_new = 0
    cdef long angle_penalty_old = 0

    cdef int angle_positions[5][2]
    cdef int intcode_lookup[5]
    cdef int offset, offset_start, offset_end, i, pos
    cdef int angle_idx
    cdef int local_move_idx = -1
    cdef int bead_flag = idx_to_bead[bead_index, 0]

    offset_start = 0
    offset_end = 0
    if bead_flag == 1:
        offset_start = 0
        offset_end = 3
    elif bead_flag == 2:
        offset_start = -2
        offset_end = 3
    elif bead_flag == 4:
        offset_start = -1
        offset_end = 2
    elif bead_flag == 3:
        offset_start = -2
        offset_end = 1
    elif bead_flag == 5:
        offset_start = -1
        offset_end = 3
    elif bead_flag == 6:
        offset_start = -2
        offset_end = 2

    angle_idx = 0
    for offset in range(offset_start, offset_end):
        pos = bead_index + offset
        angle_positions[angle_idx][0] = idx_to_bead[pos, 5]
        angle_positions[angle_idx][1] = idx_to_bead[pos, 6]
        intcode_lookup[angle_idx] = idx_to_bead[pos, 2]
        if pos == bead_index:
            local_move_idx = angle_idx
        angle_idx = angle_idx + 1

    for i in range(0, angle_idx - 2):
        a0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i][0])
        a1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i][1])
        b0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i + 2][0])
        b1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i + 2][1])
        angle_penalty_old = angle_lookup[intcode_lookup[i + 1], a0 + 1, a1 + 1, b0 + 1, b1 + 1] + angle_penalty_old

    if local_move_idx >= 0:
        angle_positions[local_move_idx][0] = new_position[0]
        angle_positions[local_move_idx][1] = new_position[1]

    for i in range(0, angle_idx - 2):
        a0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i][0])
        a1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i][1])
        b0 = fix_angle_pbc_issues(angle_positions[i + 1][0] - angle_positions[i + 2][0])
        b1 = fix_angle_pbc_issues(angle_positions[i + 1][1] - angle_positions[i + 2][1])
        angle_penalty_new = angle_lookup[intcode_lookup[i + 1], a0 + 1, a1 + 1, b0 + 1, b1 + 1] + angle_penalty_new

    return angle_penalty_new - angle_penalty_old


# ---- 2D short/long-range interaction-energy delta ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef long get_energy_change_2D_c(NUMPY_INT_TYPE[:, :] grid,
                                 NUMPY_INT_TYPE[:, :] type_grid,
                                 int old_x, int old_y,
                                 int new_x, int new_y,
                                 int LR_vs_SR,
                                 NUMPY_INT_TYPE[:, :] interaction_table,
                                 NUMPY_INT_TYPE[:, :] LR_interaction_table,
                                 NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                                 int XDIM, int YDIM,
                                 int hardwall) noexcept nogil:

    cdef long energy_old = 0
    cdef long energy_old_empty = 0
    cdef long energy_new = 0
    cdef long energy_new_empty = 0

    cdef int tmp_x, tmp_y
    cdef int x, y
    cdef unsigned int site_bead_type, bead_type

    bead_type = type_grid[old_x, old_y]

    if LR_vs_SR == 1:
        for x in range(-3, 4):
            for y in range(-3, 4):
                tmp_x = pbc_correction(old_x + x, XDIM)
                tmp_y = pbc_correction(old_y + y, YDIM)
                site_bead_type = type_grid[tmp_x, tmp_y]

                if abs(x) < 2 and abs(y) < 2:
                    if hardwall == 1:
                        if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3:
                            site_bead_type = 0
                    energy_old = energy_old + interaction_table[bead_type, site_bead_type]
                    if site_bead_type > 0:
                        energy_old_empty = energy_old_empty + interaction_table[0, site_bead_type]
                elif abs(x) < 3 and abs(y) < 3:
                    if hardwall == 1:
                        if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3:
                            continue
                    energy_old = energy_old + LR_interaction_table[bead_type, site_bead_type]
                else:
                    if hardwall == 1:
                        if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3:
                            continue
                    energy_old = energy_old + SLR_interaction_table[bead_type, site_bead_type]

        type_grid[new_x, new_y] = type_grid[old_x, old_y]
        type_grid[old_x, old_y] = 0

        for x in range(-3, 4):
            for y in range(-3, 4):
                tmp_x = pbc_correction(new_x + x, XDIM)
                tmp_y = pbc_correction(new_y + y, YDIM)
                site_bead_type = type_grid[tmp_x, tmp_y]

                if abs(x) < 2 and abs(y) < 2:
                    if hardwall == 1:
                        if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3:
                            site_bead_type = 0
                    energy_new = energy_new + interaction_table[bead_type, site_bead_type]
                    if site_bead_type > 0:
                        energy_new_empty = energy_new_empty + interaction_table[0, site_bead_type]
                elif abs(x) < 3 and abs(y) < 3:
                    if hardwall == 1:
                        if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3:
                            continue
                    energy_new = energy_new + LR_interaction_table[bead_type, site_bead_type]
                else:
                    if hardwall == 1:
                        if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3:
                            continue
                    energy_new = energy_new + SLR_interaction_table[bead_type, site_bead_type]

    else:
        for x in range(-1, 2):
            for y in range(-1, 2):
                tmp_x = pbc_correction(old_x + x, XDIM)
                tmp_y = pbc_correction(old_y + y, YDIM)
                site_bead_type = type_grid[tmp_x, tmp_y]
                if hardwall == 1:
                    if abs(tmp_x - old_x) > 3 or abs(tmp_y - old_y) > 3:
                        site_bead_type = 0
                energy_old = energy_old + interaction_table[bead_type, site_bead_type]
                if site_bead_type > 0:
                    energy_old_empty = energy_old_empty + interaction_table[0, site_bead_type]

        type_grid[new_x, new_y] = type_grid[old_x, old_y]
        type_grid[old_x, old_y] = 0

        for x in range(-1, 2):
            for y in range(-1, 2):
                tmp_x = pbc_correction(new_x + x, XDIM)
                tmp_y = pbc_correction(new_y + y, YDIM)
                site_bead_type = type_grid[tmp_x, tmp_y]
                if hardwall == 1:
                    if abs(tmp_x - new_x) > 3 or abs(tmp_y - new_y) > 3:
                        site_bead_type = 0
                energy_new = energy_new + interaction_table[bead_type, site_bead_type]
                if site_bead_type > 0:
                    energy_new_empty = energy_new_empty + interaction_table[0, site_bead_type]

    energy_old_empty = energy_old_empty - interaction_table[0, bead_type]
    energy_old = energy_old - interaction_table[bead_type, bead_type]
    energy_new_empty = energy_new_empty - interaction_table[0, bead_type]
    energy_new = energy_new - interaction_table[bead_type, bead_type]

    type_grid[old_x, old_y] = type_grid[new_x, new_y]
    type_grid[new_x, new_y] = 0

    return (energy_new + energy_old_empty) - (energy_old + energy_new_empty)


# ---- public 2D entry point (signature matches mega_crank_2D.mega_crank_2D) ----
@cython.wraparound(False)
@cython.boundscheck(False)
def mega_crank_2D(NUMPY_INT_TYPE[:, :] grid,
                  NUMPY_INT_TYPE[:, :] type_grid,
                  NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                  NUMPY_INT_TYPE[:, :] interaction_table,
                  NUMPY_INT_TYPE[:, :] LR_interaction_table,
                  NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                  NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup,
                  long energy,
                  float invtemp,
                  int nsteps,
                  NUMPY_INT_TYPE_long[:] bead_selector,
                  int passed_seed,
                  int hardwall):
    """
    Optimized, allocation-free 2D crankshaft Monte Carlo kernel.

    Behaviourally identical (same RNG stream) to
    pimms.mega_crank_2D.mega_crank_2D. Returns (energy, accepted_moves).
    """
    mc_seed(passed_seed)

    cdef unsigned int i
    cdef int bead_index
    cdef int accepted_moves = 0
    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef long delta_energy, delta_angle_energy
    cdef int beadflag, ok

    cdef int old_x, old_y, lr_vs_sr, bead_id
    cdef int new_position[2]
    cdef int anchor0, anchor1
    cdef int Nx, Ny, Cx, Cy

    for i in range(nsteps):
        bead_index = bead_selector[i]

        beadflag = idx_to_bead[bead_index, 0]
        old_x = idx_to_bead[bead_index, 5]
        old_y = idx_to_bead[bead_index, 6]
        lr_vs_sr = idx_to_bead[bead_index, 1]
        bead_id = idx_to_bead[bead_index, 4]

        # ---- single bead (flag 0) ----
        if beadflag == 0:
            ok = single_bead_crank_2D_c(old_x, old_y, grid, XDIM, YDIM, new_position)

        # ---- N-terminal bead (flag 1) ----
        elif beadflag == 1:
            anchor0 = idx_to_bead[bead_index + 1, 5]
            anchor1 = idx_to_bead[bead_index + 1, 6]
            ok = single_bead_crank_2D_c(anchor0, anchor1, grid, XDIM, YDIM, new_position)
            if hardwall == 1:
                if straddle_pair(new_position[0], new_position[1], 0,
                                 anchor0, anchor1, 0) == 1:
                    continue

        # ---- C-terminal bead (flag 3) ----
        elif beadflag == 3:
            anchor0 = idx_to_bead[bead_index - 1, 5]
            anchor1 = idx_to_bead[bead_index - 1, 6]
            ok = single_bead_crank_2D_c(anchor0, anchor1, grid, XDIM, YDIM, new_position)
            if hardwall == 1:
                if straddle_pair(anchor0, anchor1, 0,
                                 new_position[0], new_position[1], 0) == 1:
                    continue

        # ---- internal bead (flags 2,4,5,6) ----
        else:
            Nx = idx_to_bead[bead_index - 1, 5]
            Ny = idx_to_bead[bead_index - 1, 6]
            Cx = idx_to_bead[bead_index + 1, 5]
            Cy = idx_to_bead[bead_index + 1, 6]
            ok = crank_it_2D_c(Nx, Ny, Cx, Cy, grid, XDIM, YDIM, new_position)
            if hardwall == 1:
                if (straddle_pair(Nx, Ny, 0,
                                  new_position[0], new_position[1], 0) == 1
                    or straddle_pair(new_position[0], new_position[1], 0,
                                     Cx, Cy, 0) == 1):
                    continue

        if ok == 1:
            delta_energy = get_energy_change_2D_c(grid, type_grid,
                                                  old_x, old_y,
                                                  new_position[0], new_position[1],
                                                  lr_vs_sr,
                                                  interaction_table, LR_interaction_table,
                                                  SLR_interaction_table,
                                                  XDIM, YDIM, hardwall)
            delta_angle_energy = get_angle_energy_change_2D_c(bead_index, idx_to_bead,
                                                              new_position, angle_lookup)

            if accept_or_reject(invtemp, energy, energy + delta_energy + delta_angle_energy) == 1:
                grid[old_x, old_y] = 0
                grid[new_position[0], new_position[1]] = bead_id
                type_grid[new_position[0], new_position[1]] = type_grid[old_x, old_y]
                type_grid[old_x, old_y] = 0
                idx_to_bead[bead_index, 5] = new_position[0]
                idx_to_bead[bead_index, 6] = new_position[1]
                energy = energy + delta_energy + delta_angle_energy
                accepted_moves = accepted_moves + 1

    return (energy, accepted_moves)


# =====================================================================
#
#   PARALLEL (checkerboard / domain-decomposition) crankshaft kernel
#
#   Correctness argument (see pimms/fast_kernels/CHECKERBOARD_DESIGN.md):
#
#     * A single crankshaft/terminal move displaces a bead by at most 2
#       lattice sites and the energy evaluation reads +/- R_int around the
#       old and new positions (R_int = 3 if any long-range beads exist, else
#       1). Chain-neighbour (anchor / angle) reads are <= 2 sites away.
#       So the entire read+write footprint of a move lies within
#       W = R_int + 2 sites of the bead's current position.
#
#     * The box is split into blocks; only beads at least W sites inside their
#       block (a frozen halo of width W on every blocked face) are moved.
#       Hence every block's footprint stays strictly inside that block, blocks
#       are pairwise disjoint, and they can be updated concurrently with NO
#       locks and NO colouring. (A dim that is not split, nb==1, has no halo -
#       PBC wrap is fine because that whole dimension belongs to one block.)
#
#     * Metropolis acceptance depends only on dE, so each block accumulates a
#       private energy delta; the global energy is base + sum(block deltas).
#       Per-block independent PRNG streams avoid the non-thread-safe libc rand.
#
#     * Beads in the frozen halo are skipped this sweep; a fresh random origin
#       shift each call restores ergodicity over many sweeps. The kernel is
#       therefore NOT bit-identical to the serial kernel - it is a different
#       (but equally valid) Markov chain targeting the same Boltzmann
#       distribution, and is validated by ensemble + energy-consistency checks
#       (pimms/fast_kernels/benchmark_parallel.py).
#
# =====================================================================

# ---- per-thread PRNG (splitmix64) -----------------------------------
cdef inline unsigned long long splitmix64(unsigned long long* state) noexcept nogil:
    state[0] = state[0] + <unsigned long long>0x9E3779B97F4A7C15
    cdef unsigned long long z = state[0]
    z = (z ^ (z >> 30)) * <unsigned long long>0xBF58476D1CE4E5B9
    z = (z ^ (z >> 27)) * <unsigned long long>0x94D049BB133111EB
    return z ^ (z >> 31)

cdef inline double rng_uniform(unsigned long long* state) noexcept nogil:
    # 53-bit uniform in [0, 1)
    return <double>(splitmix64(state) >> 11) * (1.0 / 9007199254740992.0)

cdef inline int rng_randint(unsigned long long* state, int start, int end) noexcept nogil:
    # mirrors the integer RANGE of the serial randint(start, end)
    if start == 0:
        end = end + 1
    return start + <int>(rng_uniform(state) * end)

cdef inline int accept_p(float invtemp, long delta, unsigned long long* rng) noexcept nogil:
    cdef double expterm
    if delta <= 0:
        return 1
    expterm = exp(-(<double>delta) * invtemp)
    if rng_uniform(rng) < expterm:
        return 1
    return 0


# ---- move proposals using a per-thread PRNG (mirror the _c versions) ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int single_bead_crank_cp(int ox, int oy, int oz,
                                     NUMPY_INT_TYPE[:, :, :] grid,
                                     int XDIM, int YDIM, int ZDIM,
                                     unsigned long long* rng, int* out) noexcept nogil:
    cdef int x_off = rng_randint(rng, 0, 2) - 1
    cdef int y_off = rng_randint(rng, 0, 2) - 1
    cdef int z_off = rng_randint(rng, 0, 2) - 1
    cdef int local_x = pbc_correction(ox + x_off, XDIM)
    cdef int local_y = pbc_correction(oy + y_off, YDIM)
    cdef int local_z = pbc_correction(oz + z_off, ZDIM)
    if grid[local_x, local_y, local_z] > 0:
        out[0] = -1
        out[1] = 0
        out[2] = 0
        return 0
    out[0] = local_x
    out[1] = local_y
    out[2] = local_z
    return 1


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int crank_it_cp(int N_side_x, int N_side_y, int N_side_z,
                            int C_side_x, int C_side_y, int C_side_z,
                            NUMPY_INT_TYPE[:, :, :] grid,
                            int XDIM, int YDIM, int ZDIM,
                            unsigned long long* rng, int* out) noexcept nogil:
    cdef int x_min, x_max, y_min, y_max, z_min, z_max
    cdef int local_x, local_y, local_z

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

    x_min = int_max(N_side_x - 1, C_side_x - 1)
    x_max = int_min(N_side_x + 1, C_side_x + 1)
    y_min = int_max(N_side_y - 1, C_side_y - 1)
    y_max = int_min(N_side_y + 1, C_side_y + 1)
    z_min = int_max(N_side_z - 1, C_side_z - 1)
    z_max = int_min(N_side_z + 1, C_side_z + 1)

    local_x = pbc_correction((x_min + rng_randint(rng, 1, (x_max - x_min + 1)) - 1), XDIM)
    local_y = pbc_correction((y_min + rng_randint(rng, 1, (y_max - y_min + 1)) - 1), YDIM)
    local_z = pbc_correction((z_min + rng_randint(rng, 1, (z_max - z_min + 1)) - 1), ZDIM)

    if grid[local_x, local_y, local_z] > 0:
        out[0] = -1
        out[1] = 0
        out[2] = 0
        return 0
    out[0] = local_x
    out[1] = local_y
    out[2] = local_z
    return 1


# ---- per-block worker (runs entirely without the GIL) ---------------
@cython.wraparound(False)
@cython.boundscheck(False)
cdef void run_block(int b,
                    int[::1] ids, int[::1] starts, int[::1] attempts,
                    unsigned long long[::1] seeds,
                    int[::1] blo_x, int[::1] blo_y, int[::1] blo_z,
                    int Lx, int Ly, int Lz,
                    int nbx, int nby, int nbz,
                    int shift_x, int shift_y, int shift_z, int W,
                    NUMPY_INT_TYPE[:, :, :] grid,
                    NUMPY_INT_TYPE[:, :, :] type_grid,
                    NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                    NUMPY_INT_TYPE[:, :] interaction_table,
                    NUMPY_INT_TYPE[:, :] LR_interaction_table,
                    NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                    NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup,
                    float invtemp, int XDIM, int YDIM, int ZDIM, int hardwall,
                    long[::1] out_delta, int[::1] out_accepted) noexcept nogil:

    cdef int lo = starts[b]
    cdef int hi = starts[b + 1]
    cdef int n_block = hi - lo
    out_delta[b] = 0
    out_accepted[b] = 0
    if n_block <= 0:
        return

    cdef unsigned long long rng = seeds[b]
    cdef long dsum = 0
    cdef int acc = 0
    cdef int k, pick, bead_index, beadflag, ok
    cdef int gx, gy, gz, sx, sy, sz, wx, wy, wz
    cdef int old0, old1, old2
    cdef int new_position[3]
    cdef int anchor0, anchor1, anchor2
    cdef int Nx, Ny, Nz, Cx, Cy, Cz
    cdef long delta_energy, delta_angle_energy

    for k in range(attempts[b]):
        # pick a random bead from this block's list
        pick = lo + <int>(rng_uniform(&rng) * n_block)
        if pick >= hi:
            pick = hi - 1
        bead_index = ids[pick]

        gx = idx_to_bead[bead_index, 5]
        gy = idx_to_bead[bead_index, 6]
        gz = idx_to_bead[bead_index, 7]

        # re-test interior membership against the CURRENT position (a bead may
        # have drifted to the frozen halo during this block's run). Skip if not
        # at least W inside every blocked dimension - guarantees the footprint
        # stays in-block.
        if nbx > 1:
            sx = (gx - shift_x + XDIM) % XDIM
            wx = sx - blo_x[b]
            if wx < W or wx >= Lx - W:
                continue
        if nby > 1:
            sy = (gy - shift_y + YDIM) % YDIM
            wy = sy - blo_y[b]
            if wy < W or wy >= Ly - W:
                continue
        if nbz > 1:
            sz = (gz - shift_z + ZDIM) % ZDIM
            wz = sz - blo_z[b]
            if wz < W or wz >= Lz - W:
                continue

        old0 = gx
        old1 = gy
        old2 = gz
        beadflag = idx_to_bead[bead_index, 0]

        if beadflag == 0:
            single_bead_crank_cp(old0, old1, old2, grid, XDIM, YDIM, ZDIM, &rng, new_position)

        elif beadflag == 1:
            anchor0 = idx_to_bead[bead_index + 1, 5]
            anchor1 = idx_to_bead[bead_index + 1, 6]
            anchor2 = idx_to_bead[bead_index + 1, 7]
            single_bead_crank_cp(anchor0, anchor1, anchor2, grid, XDIM, YDIM, ZDIM, &rng, new_position)
            if hardwall == 1:
                if straddle_pair(new_position[0], new_position[1], new_position[2],
                                 anchor0, anchor1, anchor2) == 1:
                    continue

        elif beadflag == 3:
            anchor0 = idx_to_bead[bead_index - 1, 5]
            anchor1 = idx_to_bead[bead_index - 1, 6]
            anchor2 = idx_to_bead[bead_index - 1, 7]
            single_bead_crank_cp(anchor0, anchor1, anchor2, grid, XDIM, YDIM, ZDIM, &rng, new_position)
            if hardwall == 1:
                if straddle_pair(anchor0, anchor1, anchor2,
                                 new_position[0], new_position[1], new_position[2]) == 1:
                    continue

        else:
            Nx = idx_to_bead[bead_index - 1, 5]
            Ny = idx_to_bead[bead_index - 1, 6]
            Nz = idx_to_bead[bead_index - 1, 7]
            Cx = idx_to_bead[bead_index + 1, 5]
            Cy = idx_to_bead[bead_index + 1, 6]
            Cz = idx_to_bead[bead_index + 1, 7]
            crank_it_cp(Nx, Ny, Nz, Cx, Cy, Cz, grid, XDIM, YDIM, ZDIM, &rng, new_position)
            if hardwall == 1:
                if (straddle_pair(Nx, Ny, Nz,
                                  new_position[0], new_position[1], new_position[2]) == 1
                    or straddle_pair(new_position[0], new_position[1], new_position[2],
                                     Cx, Cy, Cz) == 1):
                    continue

        if not new_position[0] < 0:
            # DETAILED BALANCE: reject any move whose NEW position would leave
            # the movable interior (i.e. land in the frozen halo). A halo bead is
            # never selected, so it could never make the reverse move; allowing
            # the forward move (interior -> halo) without the reverse breaks
            # detailed balance and biases the sampled ensemble. Rejecting these
            # out-of-interior proposals (the bead simply stays put) keeps the
            # movable region closed and the move set symmetric.
            if nbx > 1:
                sx = (new_position[0] - shift_x + XDIM) % XDIM
                wx = sx - blo_x[b]
                if wx < W or wx >= Lx - W:
                    continue
            if nby > 1:
                sy = (new_position[1] - shift_y + YDIM) % YDIM
                wy = sy - blo_y[b]
                if wy < W or wy >= Ly - W:
                    continue
            if nbz > 1:
                sz = (new_position[2] - shift_z + ZDIM) % ZDIM
                wz = sz - blo_z[b]
                if wz < W or wz >= Lz - W:
                    continue

            delta_energy = get_energy_change_c(grid, type_grid,
                                               old0, old1, old2,
                                               new_position[0], new_position[1], new_position[2],
                                               idx_to_bead[bead_index, 1],
                                               interaction_table, LR_interaction_table,
                                               SLR_interaction_table,
                                               XDIM, YDIM, ZDIM, hardwall)
            delta_angle_energy = get_angle_energy_change_c(bead_index, idx_to_bead,
                                                           new_position, angle_lookup)
            if accept_p(invtemp, delta_energy + delta_angle_energy, &rng) == 1:
                grid[old0, old1, old2] = 0
                grid[new_position[0], new_position[1], new_position[2]] = idx_to_bead[bead_index, 4]
                type_grid[new_position[0], new_position[1], new_position[2]] = type_grid[old0, old1, old2]
                type_grid[old0, old1, old2] = 0
                idx_to_bead[bead_index, 5] = new_position[0]
                idx_to_bead[bead_index, 6] = new_position[1]
                idx_to_bead[bead_index, 7] = new_position[2]
                dsum = dsum + delta_energy + delta_angle_energy
                acc = acc + 1

    out_delta[b] = dsum
    out_accepted[b] = acc


def _choose_nb(int DIM, int W, int cap):
    """Pick the number of blocks along one dimension: large enough blocks
    (>= 4W so the movable interior is a healthy fraction) and at most `cap`."""
    cdef int max_nb = DIM // (4 * W)
    if max_nb < 1:
        return 1
    if max_nb > cap:
        return cap
    return max_nb


@cython.wraparound(False)
@cython.boundscheck(False)
def mega_crank_parallel(NUMPY_INT_TYPE[:, :, :] grid,
                        NUMPY_INT_TYPE[:, :, :] type_grid,
                        NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                        NUMPY_INT_TYPE[:, :] interaction_table,
                        NUMPY_INT_TYPE[:, :] LR_interaction_table,
                        NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                        NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup,
                        long energy,
                        float invtemp,
                        int nsteps,
                        int passed_seed,
                        int hardwall,
                        int num_threads):
    """
    Parallel checkerboard crankshaft kernel.

    Same role as mega_crank() but distributes the substep work across
    `num_threads` OpenMP threads using a frozen-halo block decomposition.
    Does NOT take a pre-built bead_selector (it selects beads per block); does
    NOT preserve the serial RNG stream. Returns (energy, accepted_moves).

    Falls back to a single-block (effectively serial) sweep when the box is too
    small to decompose.
    """
    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int ZDIM = grid.shape[2]
    cdef int num_beads = idx_to_bead.shape[0]

    # interaction radius -> halo width W = R_int + 2 (2 covers the max bead
    # displacement of a terminal/internal move; R_int covers the energy read).
    idx_np = np.asarray(idx_to_bead)
    cdef int R_int = 3 if np.any(idx_np[:, 1] == 1) else 1
    cdef int W = R_int + 2

    # block counts per dim. IMPORTANT: this depends ONLY on box geometry (and W),
    # NOT on num_threads, so the decomposition - and therefore the result - is
    # identical for any thread count. Threads only change how fast the fixed set
    # of independent blocks is processed. cap=4 -> up to 4^3=64 blocks, which is
    # plenty for dynamic load-balancing across typical core counts.
    cdef int cap = 4
    cdef int nbx = _choose_nb(XDIM, W, cap)
    cdef int nby = _choose_nb(YDIM, W, cap)
    cdef int nbz = _choose_nb(ZDIM, W, cap)
    cdef int Lx = XDIM // nbx
    cdef int Ly = YDIM // nby
    cdef int Lz = ZDIM // nbz
    cdef int num_blocks = nbx * nby * nbz

    # reproducible shifts + per-block seeds derived from passed_seed
    rstate = np.random.RandomState(passed_seed & 0x7FFFFFFF)
    cdef int shift_x = int(rstate.randint(0, Lx)) if nbx > 1 else 0
    cdef int shift_y = int(rstate.randint(0, Ly)) if nby > 1 else 0
    cdef int shift_z = int(rstate.randint(0, Lz)) if nbz > 1 else 0

    # ---- bucket beads into blocks (interior beads only) -----------------
    gx = idx_np[:, 5].astype(np.int64)
    gy = idx_np[:, 6].astype(np.int64)
    gz = idx_np[:, 7].astype(np.int64)

    # per-dim block index (-1 => frozen: in halo or remainder or wrong region)
    def dim_block(g, nb, L, DIM, shift):
        if nb == 1:
            return np.zeros(num_beads, dtype=np.int64)        # whole dim = block 0
        s = (g - shift) % DIM
        bj = s // L
        within = s - bj * L
        # frozen if outside valid region (remainder) or inside the W halo
        bad = (s >= nb * L) | (within < W) | (within >= L - W)
        bj = bj.copy()
        bj[bad] = -1
        return bj

    bxj = dim_block(gx, nbx, Lx, XDIM, shift_x)
    byj = dim_block(gy, nby, Ly, YDIM, shift_y)
    bzj = dim_block(gz, nbz, Lz, ZDIM, shift_z)

    movable = (bxj >= 0) & (byj >= 0) & (bzj >= 0)
    block_of_bead = (bxj * nby * nbz + byj * nbz + bzj)
    block_of_bead[~movable] = -1

    sel = np.nonzero(movable)[0].astype(np.int32)
    if sel.shape[0] == 0:
        # nothing movable this sweep (tiny box / unlucky shift) - no-op
        return (energy, 0)

    sel_blocks = block_of_bead[sel]
    order = np.argsort(sel_blocks, kind="stable")
    ids = sel[order].astype(np.int32)
    sorted_blocks = sel_blocks[order]

    counts = np.bincount(sorted_blocks, minlength=num_blocks).astype(np.int32)
    starts = np.zeros(num_blocks + 1, dtype=np.int32)
    starts[1:] = np.cumsum(counts)

    # distribute nsteps proportional to each block's movable bead count
    total_interior = int(counts.sum())
    attempts = np.maximum(
        (nsteps * counts.astype(np.int64) // max(total_interior, 1)), 0).astype(np.int32)

    # shifted lower bound of each block per dim (for the in-worker interior test)
    bix = np.arange(num_blocks, dtype=np.int32) // (nby * nbz)
    biy = (np.arange(num_blocks, dtype=np.int32) // nbz) % nby
    biz = np.arange(num_blocks, dtype=np.int32) % nbz
    blo_x = (bix * Lx).astype(np.int32)
    blo_y = (biy * Ly).astype(np.int32)
    blo_z = (biz * Lz).astype(np.int32)

    # independent PRNG seed per block
    seeds = (np.arange(num_blocks, dtype=np.uint64) + np.uint64(passed_seed) + np.uint64(0x9E3779B9)) \
        * np.uint64(0x2545F4914F6CDD1D) + np.uint64(1)

    out_delta = np.zeros(num_blocks, dtype=np.int64)
    out_accepted = np.zeros(num_blocks, dtype=np.int32)

    # typed memoryviews for the nogil region
    cdef int[::1] ids_mv = ids
    cdef int[::1] starts_mv = starts
    cdef int[::1] attempts_mv = attempts
    cdef unsigned long long[::1] seeds_mv = seeds
    cdef int[::1] blo_x_mv = blo_x
    cdef int[::1] blo_y_mv = blo_y
    cdef int[::1] blo_z_mv = blo_z
    cdef long[::1] out_delta_mv = out_delta
    cdef int[::1] out_accepted_mv = out_accepted

    cdef int b
    cdef int nthreads = num_threads if num_threads > 0 else 1

    # ---- parallel region: each block is independent -----------------
    for b in prange(num_blocks, nogil=True, num_threads=nthreads, schedule='dynamic'):
        run_block(b, ids_mv, starts_mv, attempts_mv, seeds_mv,
                  blo_x_mv, blo_y_mv, blo_z_mv,
                  Lx, Ly, Lz, nbx, nby, nbz,
                  shift_x, shift_y, shift_z, W,
                  grid, type_grid, idx_to_bead,
                  interaction_table, LR_interaction_table, SLR_interaction_table,
                  angle_lookup, invtemp, XDIM, YDIM, ZDIM, hardwall,
                  out_delta_mv, out_accepted_mv)

    cdef long total_delta = 0
    cdef int total_accepted = 0
    for b in range(num_blocks):
        total_delta += out_delta_mv[b]
        total_accepted += out_accepted_mv[b]

    return (energy + total_delta, total_accepted)


# =====================================================================
#   PARALLEL checkerboard crankshaft kernel - 2D
#
#   Exact 2D analogue of the 3D kernel above: same frozen-halo block
#   decomposition (independent of num_threads), same per-block splitmix64
#   streams, same detailed-balance halo handling - just over a 2D grid using
#   the 2D move/energy helpers. Targets the same Boltzmann distribution as the
#   serial 2D kernel (mega_crank_2D); not bit-identical to it.
# =====================================================================

# ---- 2D move proposals using a per-thread PRNG (mirror the _c versions) ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int single_bead_crank_cp_2D(int ox, int oy,
                                        NUMPY_INT_TYPE[:, :] grid,
                                        int XDIM, int YDIM,
                                        unsigned long long* rng, int* out) noexcept nogil:
    cdef int x_off = rng_randint(rng, 0, 2) - 1
    cdef int y_off = rng_randint(rng, 0, 2) - 1
    cdef int local_x = pbc_correction(ox + x_off, XDIM)
    cdef int local_y = pbc_correction(oy + y_off, YDIM)
    if grid[local_x, local_y] > 0:
        out[0] = -1
        out[1] = 0
        return 0
    out[0] = local_x
    out[1] = local_y
    return 1


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int crank_it_cp_2D(int N_side_x, int N_side_y,
                               int C_side_x, int C_side_y,
                               NUMPY_INT_TYPE[:, :] grid,
                               int XDIM, int YDIM,
                               unsigned long long* rng, int* out) noexcept nogil:
    cdef int x_min, x_max, y_min, y_max, local_x, local_y

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

    x_min = int_max(N_side_x - 1, C_side_x - 1)
    x_max = int_min(N_side_x + 1, C_side_x + 1)
    y_min = int_max(N_side_y - 1, C_side_y - 1)
    y_max = int_min(N_side_y + 1, C_side_y + 1)

    local_x = pbc_correction((x_min + rng_randint(rng, 1, (x_max - x_min + 1)) - 1), XDIM)
    local_y = pbc_correction((y_min + rng_randint(rng, 1, (y_max - y_min + 1)) - 1), YDIM)

    if grid[local_x, local_y] > 0:
        out[0] = -1
        out[1] = 0
        return 0
    out[0] = local_x
    out[1] = local_y
    return 1


# ---- per-block worker (2D), runs entirely without the GIL ----
@cython.wraparound(False)
@cython.boundscheck(False)
cdef void run_block_2D(int b,
                       int[::1] ids, int[::1] starts, int[::1] attempts,
                       unsigned long long[::1] seeds,
                       int[::1] blo_x, int[::1] blo_y,
                       int Lx, int Ly,
                       int nbx, int nby,
                       int shift_x, int shift_y, int W,
                       NUMPY_INT_TYPE[:, :] grid,
                       NUMPY_INT_TYPE[:, :] type_grid,
                       NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                       NUMPY_INT_TYPE[:, :] interaction_table,
                       NUMPY_INT_TYPE[:, :] LR_interaction_table,
                       NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                       NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup,
                       float invtemp, int XDIM, int YDIM, int hardwall,
                       long[::1] out_delta, int[::1] out_accepted) noexcept nogil:

    cdef int lo = starts[b]
    cdef int hi = starts[b + 1]
    cdef int n_block = hi - lo
    out_delta[b] = 0
    out_accepted[b] = 0
    if n_block <= 0:
        return

    cdef unsigned long long rng = seeds[b]
    cdef long dsum = 0
    cdef int acc = 0
    cdef int k, pick, bead_index, beadflag
    cdef int gx, gy, sx, sy, wx, wy
    cdef int old0, old1
    cdef int new_position[2]
    cdef int anchor0, anchor1
    cdef int Nx, Ny, Cx, Cy
    cdef long delta_energy, delta_angle_energy

    for k in range(attempts[b]):
        # pick a random bead from this block's list
        pick = lo + <int>(rng_uniform(&rng) * n_block)
        if pick >= hi:
            pick = hi - 1
        bead_index = ids[pick]

        gx = idx_to_bead[bead_index, 5]
        gy = idx_to_bead[bead_index, 6]

        # re-test interior membership against the CURRENT position (a bead may
        # have drifted into the frozen halo during this block's run).
        if nbx > 1:
            sx = (gx - shift_x + XDIM) % XDIM
            wx = sx - blo_x[b]
            if wx < W or wx >= Lx - W:
                continue
        if nby > 1:
            sy = (gy - shift_y + YDIM) % YDIM
            wy = sy - blo_y[b]
            if wy < W or wy >= Ly - W:
                continue

        old0 = gx
        old1 = gy
        beadflag = idx_to_bead[bead_index, 0]

        if beadflag == 0:
            single_bead_crank_cp_2D(old0, old1, grid, XDIM, YDIM, &rng, new_position)

        elif beadflag == 1:
            anchor0 = idx_to_bead[bead_index + 1, 5]
            anchor1 = idx_to_bead[bead_index + 1, 6]
            single_bead_crank_cp_2D(anchor0, anchor1, grid, XDIM, YDIM, &rng, new_position)
            if hardwall == 1:
                if straddle_pair(new_position[0], new_position[1], 0,
                                 anchor0, anchor1, 0) == 1:
                    continue

        elif beadflag == 3:
            anchor0 = idx_to_bead[bead_index - 1, 5]
            anchor1 = idx_to_bead[bead_index - 1, 6]
            single_bead_crank_cp_2D(anchor0, anchor1, grid, XDIM, YDIM, &rng, new_position)
            if hardwall == 1:
                if straddle_pair(anchor0, anchor1, 0,
                                 new_position[0], new_position[1], 0) == 1:
                    continue

        else:
            Nx = idx_to_bead[bead_index - 1, 5]
            Ny = idx_to_bead[bead_index - 1, 6]
            Cx = idx_to_bead[bead_index + 1, 5]
            Cy = idx_to_bead[bead_index + 1, 6]
            crank_it_cp_2D(Nx, Ny, Cx, Cy, grid, XDIM, YDIM, &rng, new_position)
            if hardwall == 1:
                if (straddle_pair(Nx, Ny, 0,
                                  new_position[0], new_position[1], 0) == 1
                    or straddle_pair(new_position[0], new_position[1], 0,
                                     Cx, Cy, 0) == 1):
                    continue

        if not new_position[0] < 0:
            # DETAILED BALANCE: reject any move whose NEW position would leave the
            # movable interior (land in the frozen halo); see run_block for the
            # full argument.
            if nbx > 1:
                sx = (new_position[0] - shift_x + XDIM) % XDIM
                wx = sx - blo_x[b]
                if wx < W or wx >= Lx - W:
                    continue
            if nby > 1:
                sy = (new_position[1] - shift_y + YDIM) % YDIM
                wy = sy - blo_y[b]
                if wy < W or wy >= Ly - W:
                    continue

            delta_energy = get_energy_change_2D_c(grid, type_grid,
                                                  old0, old1,
                                                  new_position[0], new_position[1],
                                                  idx_to_bead[bead_index, 1],
                                                  interaction_table, LR_interaction_table,
                                                  SLR_interaction_table,
                                                  XDIM, YDIM, hardwall)
            delta_angle_energy = get_angle_energy_change_2D_c(bead_index, idx_to_bead,
                                                              new_position, angle_lookup)
            if accept_p(invtemp, delta_energy + delta_angle_energy, &rng) == 1:
                grid[old0, old1] = 0
                grid[new_position[0], new_position[1]] = idx_to_bead[bead_index, 4]
                type_grid[new_position[0], new_position[1]] = type_grid[old0, old1]
                type_grid[old0, old1] = 0
                idx_to_bead[bead_index, 5] = new_position[0]
                idx_to_bead[bead_index, 6] = new_position[1]
                dsum = dsum + delta_energy + delta_angle_energy
                acc = acc + 1

    out_delta[b] = dsum
    out_accepted[b] = acc


@cython.wraparound(False)
@cython.boundscheck(False)
def mega_crank_parallel_2D(NUMPY_INT_TYPE[:, :] grid,
                           NUMPY_INT_TYPE[:, :] type_grid,
                           NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                           NUMPY_INT_TYPE[:, :] interaction_table,
                           NUMPY_INT_TYPE[:, :] LR_interaction_table,
                           NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                           NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup,
                           long energy,
                           float invtemp,
                           int nsteps,
                           int passed_seed,
                           int hardwall,
                           int num_threads):
    """
    Parallel checkerboard crankshaft kernel (2D).

    The 2D analogue of mega_crank_parallel: distributes the substep work across
    `num_threads` OpenMP threads using the same frozen-halo block decomposition,
    over a 2D grid. The decomposition depends only on box geometry (and the halo
    width W), NOT on num_threads, so the result is identical for any thread count.
    Targets the same Boltzmann distribution as the serial 2D kernel
    (mega_crank_2D); it is not bit-identical to it. Returns (energy, accepted).
    """
    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int num_beads = idx_to_bead.shape[0]

    # interaction radius -> halo width W = R_int + 2 (2 covers the max bead
    # displacement of a terminal/internal move; R_int covers the energy read).
    idx_np = np.asarray(idx_to_bead)
    cdef int R_int = 3 if np.any(idx_np[:, 1] == 1) else 1
    cdef int W = R_int + 2

    cdef int cap = 4
    cdef int nbx = _choose_nb(XDIM, W, cap)
    cdef int nby = _choose_nb(YDIM, W, cap)
    cdef int Lx = XDIM // nbx
    cdef int Ly = YDIM // nby
    cdef int num_blocks = nbx * nby

    # reproducible shifts + per-block seeds derived from passed_seed
    rstate = np.random.RandomState(passed_seed & 0x7FFFFFFF)
    cdef int shift_x = int(rstate.randint(0, Lx)) if nbx > 1 else 0
    cdef int shift_y = int(rstate.randint(0, Ly)) if nby > 1 else 0

    # ---- bucket beads into blocks (interior beads only) -----------------
    gx = idx_np[:, 5].astype(np.int64)
    gy = idx_np[:, 6].astype(np.int64)

    # per-dim block index (-1 => frozen: in halo or remainder or wrong region)
    def dim_block(g, nb, L, DIM, shift):
        if nb == 1:
            return np.zeros(num_beads, dtype=np.int64)        # whole dim = block 0
        s = (g - shift) % DIM
        bj = s // L
        within = s - bj * L
        bad = (s >= nb * L) | (within < W) | (within >= L - W)
        bj = bj.copy()
        bj[bad] = -1
        return bj

    bxj = dim_block(gx, nbx, Lx, XDIM, shift_x)
    byj = dim_block(gy, nby, Ly, YDIM, shift_y)

    movable = (bxj >= 0) & (byj >= 0)
    block_of_bead = (bxj * nby + byj)
    block_of_bead[~movable] = -1

    sel = np.nonzero(movable)[0].astype(np.int32)
    if sel.shape[0] == 0:
        # nothing movable this sweep (tiny box / unlucky shift) - no-op
        return (energy, 0)

    sel_blocks = block_of_bead[sel]
    order = np.argsort(sel_blocks, kind="stable")
    ids = sel[order].astype(np.int32)
    sorted_blocks = sel_blocks[order]

    counts = np.bincount(sorted_blocks, minlength=num_blocks).astype(np.int32)
    starts = np.zeros(num_blocks + 1, dtype=np.int32)
    starts[1:] = np.cumsum(counts)

    # distribute nsteps proportional to each block's movable bead count
    total_interior = int(counts.sum())
    attempts = np.maximum(
        (nsteps * counts.astype(np.int64) // max(total_interior, 1)), 0).astype(np.int32)

    # shifted lower bound of each block per dim (for the in-worker interior test)
    bix = np.arange(num_blocks, dtype=np.int32) // nby
    biy = np.arange(num_blocks, dtype=np.int32) % nby
    blo_x = (bix * Lx).astype(np.int32)
    blo_y = (biy * Ly).astype(np.int32)

    # independent PRNG seed per block
    seeds = (np.arange(num_blocks, dtype=np.uint64) + np.uint64(passed_seed) + np.uint64(0x9E3779B9)) \
        * np.uint64(0x2545F4914F6CDD1D) + np.uint64(1)

    out_delta = np.zeros(num_blocks, dtype=np.int64)
    out_accepted = np.zeros(num_blocks, dtype=np.int32)

    # typed memoryviews for the nogil region
    cdef int[::1] ids_mv = ids
    cdef int[::1] starts_mv = starts
    cdef int[::1] attempts_mv = attempts
    cdef unsigned long long[::1] seeds_mv = seeds
    cdef int[::1] blo_x_mv = blo_x
    cdef int[::1] blo_y_mv = blo_y
    cdef long[::1] out_delta_mv = out_delta
    cdef int[::1] out_accepted_mv = out_accepted

    cdef int b
    cdef int nthreads = num_threads if num_threads > 0 else 1

    # ---- parallel region: each block is independent -----------------
    for b in prange(num_blocks, nogil=True, num_threads=nthreads, schedule='dynamic'):
        run_block_2D(b, ids_mv, starts_mv, attempts_mv, seeds_mv,
                     blo_x_mv, blo_y_mv,
                     Lx, Ly, nbx, nby,
                     shift_x, shift_y, W,
                     grid, type_grid, idx_to_bead,
                     interaction_table, LR_interaction_table, SLR_interaction_table,
                     angle_lookup, invtemp, XDIM, YDIM, hardwall,
                     out_delta_mv, out_accepted_mv)

    cdef long total_delta_2d = 0
    cdef int total_accepted_2d = 0
    for b in range(num_blocks):
        total_delta_2d += out_delta_mv[b]
        total_accepted_2d += out_accepted_mv[b]

    return (energy + total_delta_2d, total_accepted_2d)


def parallel_layout_info(int XDIM, int YDIM, int ZDIM, has_LR, int num_threads=1):
    """Introspection helper: returns the block decomposition the parallel kernel
    would use for a given box (decomposition is independent of num_threads)."""
    cdef int R_int = 3 if has_LR else 1
    cdef int W = R_int + 2
    cdef int cap = 4
    cdef int nbx = _choose_nb(XDIM, W, cap)
    cdef int nby = _choose_nb(YDIM, W, cap)
    cdef int nbz = _choose_nb(ZDIM, W, cap)
    return {
        "W": W, "blocks": (nbx, nby, nbz), "num_blocks": nbx * nby * nbz,
        "block_size": (XDIM // nbx, YDIM // nby, ZDIM // nbz),
    }



# =====================================================================
#
#   SLITHER (reptation) megamove kernel
#
#   A slither advances a chain forwards or backwards like a snake: one end
#   "grows" into a random adjacent empty site and the other end retracts,
#   preserving the chain's shape (path) shifted by one bead.
#
#   * Single-bead chains  -> a local translation (single_bead_crank).
#   * Homopolymers        -> O(1) interaction energy: only the vacated tail and
#                            the new head sites change type, so dE is exactly
#                            get_energy_change_c(vacated, new_end). Angle energy
#                            is the (cheap) before/after chain-angle difference.
#   * Heteropolymers      -> the reptation is decomposed into L single-bead moves
#                            (move each bead, from the leading end inwards, into
#                            the now-empty site ahead of it). dE is the sum of the
#                            already-validated get_energy_change_c +
#                            get_angle_energy_change_c over those L moves, which
#                            is exact and energy-consistent by construction.
#
#   Detailed balance: direction (head/tail grow) is 1/2 each and the new end is
#   proposed uniformly over the leading end's 27 neighbour offsets (the self
#   offset lands on the occupied end bead and is rejected), so the forward and
#   reverse proposals are symmetric and plain Metropolis on the total dE holds.
# =====================================================================

cdef inline long _triplet_angle_3D(int mx, int my, int mz,
                                   int lx, int ly, int lz,
                                   int rx, int ry, int rz,
                                   int ic,
                                   NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup) noexcept nogil:
    cdef int a0 = fix_angle_pbc_issues(mx - lx)
    cdef int a1 = fix_angle_pbc_issues(my - ly)
    cdef int a2 = fix_angle_pbc_issues(mz - lz)
    cdef int b0 = fix_angle_pbc_issues(mx - rx)
    cdef int b1 = fix_angle_pbc_issues(my - ry)
    cdef int b2 = fix_angle_pbc_issues(mz - rz)
    return angle_lookup[ic, a0 + 1, a1 + 1, a2 + 1, b0 + 1, b1 + 1, b2 + 1]


# position of chain-index k under a slither mode:
#   mode 0 = current (old) positions; 1 = head-grow proposal; 2 = tail-grow proposal
cdef inline void _smode_pos(NUMPY_INT_TYPE_long[:, :] idx, int off, int L, int mode,
                            int nex, int ney, int nez, int k, int* o) noexcept nogil:
    if mode == 0:
        o[0] = idx[off + k, 5]; o[1] = idx[off + k, 6]; o[2] = idx[off + k, 7]
    elif mode == 1:                 # head-grow: bead k takes bead (k+1)'s old site, last takes new end
        if k < L - 1:
            o[0] = idx[off + k + 1, 5]; o[1] = idx[off + k + 1, 6]; o[2] = idx[off + k + 1, 7]
        else:
            o[0] = nex; o[1] = ney; o[2] = nez
    else:                           # tail-grow: bead k takes bead (k-1)'s old site, first takes new end
        if k > 0:
            o[0] = idx[off + k - 1, 5]; o[1] = idx[off + k - 1, 6]; o[2] = idx[off + k - 1, 7]
        else:
            o[0] = nex; o[1] = ney; o[2] = nez


cdef long _chain_angle_mode(NUMPY_INT_TYPE_long[:, :] idx, int off, int L, int mode,
                            int nex, int ney, int nez,
                            NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup) noexcept nogil:
    # total angle energy of the chain under the given mode (reads idx; does not mutate)
    cdef long e = 0
    cdef int k
    cdef int pl[3]
    cdef int pm[3]
    cdef int pr[3]
    for k in range(1, L - 1):
        if idx[off + k, 3] == 1:    # skip-angle bead
            continue
        _smode_pos(idx, off, L, mode, nex, ney, nez, k - 1, pl)
        _smode_pos(idx, off, L, mode, nex, ney, nez, k, pm)
        _smode_pos(idx, off, L, mode, nex, ney, nez, k + 1, pr)
        e = e + _triplet_angle_3D(pm[0], pm[1], pm[2], pl[0], pl[1], pl[2],
                                  pr[0], pr[1], pr[2], idx[off + k, 2], angle_lookup)
    return e


cdef inline void _apply_bead_move(NUMPY_INT_TYPE[:, :, :] grid, NUMPY_INT_TYPE[:, :, :] type_grid,
                                  NUMPY_INT_TYPE_long[:, :] idx, int bead_index,
                                  int ox, int oy, int oz, int nx, int ny, int nz,
                                  int chainID) noexcept nogil:
    grid[ox, oy, oz] = 0
    grid[nx, ny, nz] = chainID
    type_grid[nx, ny, nz] = type_grid[ox, oy, oz]
    type_grid[ox, oy, oz] = 0
    idx[bead_index, 5] = nx
    idx[bead_index, 6] = ny
    idx[bead_index, 7] = nz


@cython.wraparound(False)
@cython.boundscheck(False)
def mega_slither(NUMPY_INT_TYPE[:, :, :] grid,
                 NUMPY_INT_TYPE[:, :, :] type_grid,
                 NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                 int[::1] chain_offset,
                 int[::1] chain_length,
                 int[::1] chain_homo,
                 int[::1] chain_selector,
                 NUMPY_INT_TYPE[:, :] interaction_table,
                 NUMPY_INT_TYPE[:, :] LR_interaction_table,
                 NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                 NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup,
                 long energy,
                 float invtemp,
                 int passed_seed,
                 int hardwall,
                 int max_chain_len):
    """
    3D slither (reptation) megamove. chain_selector lists chain indices to slither,
    one per entry (build it so every chain appears SLITHER_SUBSTEPS times, shuffled).
    Mutates grid / type_grid / idx_to_bead in place. Returns (energy, accepted).
    """
    mc_seed(passed_seed)

    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int ZDIM = grid.shape[2]
    cdef int n_steps = chain_selector.shape[0]
    cdef int accepted = 0

    cdef int i, c, off, L, homo, chainID, direction, k
    cdef int ex, ey, ez, nx, ny, nz, vx, vy, vz, ox, oy, oz, tx, ty, tz, cx, cy, cz
    cdef long de
    cdef int newp[3]

    # work buffers for the heteropolymer revert (old chain positions)
    cdef int* bx = <int*> malloc(max_chain_len * sizeof(int))
    cdef int* by = <int*> malloc(max_chain_len * sizeof(int))
    cdef int* bz = <int*> malloc(max_chain_len * sizeof(int))
    if bx == NULL or by == NULL or bz == NULL:
        if bx != NULL: free(bx)
        if by != NULL: free(by)
        if bz != NULL: free(bz)
        raise MemoryError()

    try:
        for i in range(n_steps):
            c = chain_selector[i]
            off = chain_offset[c]
            L = chain_length[c]
            homo = chain_homo[c]
            chainID = idx_to_bead[off, 4]

            # ----- single-bead chain: local translation -----
            if L == 1:
                ox = idx_to_bead[off, 5]; oy = idx_to_bead[off, 6]; oz = idx_to_bead[off, 7]
                if single_bead_crank_c(ox, oy, oz, grid, XDIM, YDIM, ZDIM, newp) == 1:
                    de = get_energy_change_c(grid, type_grid, ox, oy, oz,
                                             newp[0], newp[1], newp[2], idx_to_bead[off, 1],
                                             interaction_table, LR_interaction_table,
                                             SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
                    de = de + get_angle_energy_change_c(off, idx_to_bead, newp, angle_lookup)
                    if accept_or_reject(invtemp, energy, energy + de) == 1:
                        _apply_bead_move(grid, type_grid, idx_to_bead, off,
                                         ox, oy, oz, newp[0], newp[1], newp[2], chainID)
                        energy = energy + de
                        accepted = accepted + 1
                continue

            # ----- reptation: pick a growing end and propose a new end site -----
            direction = randint(0, 1)                  # 0 = head-grow (last bead), 1 = tail-grow (first bead)
            if direction == 0:
                ex = idx_to_bead[off + L - 1, 5]; ey = idx_to_bead[off + L - 1, 6]; ez = idx_to_bead[off + L - 1, 7]
            else:
                ex = idx_to_bead[off, 5]; ey = idx_to_bead[off, 6]; ez = idx_to_bead[off, 7]

            nx = pbc_correction(ex + randint(0, 2) - 1, XDIM)
            ny = pbc_correction(ey + randint(0, 2) - 1, YDIM)
            nz = pbc_correction(ez + randint(0, 2) - 1, ZDIM)

            if grid[nx, ny, nz] > 0:                   # occupied (incl. self when offset is 0) -> reject
                continue
            if hardwall == 1:
                if straddle_pair(ex, ey, ez, nx, ny, nz) == 1:
                    continue

            if homo == 1:
                # vacated end (opposite the growing end)
                if direction == 0:
                    vx = idx_to_bead[off, 5]; vy = idx_to_bead[off, 6]; vz = idx_to_bead[off, 7]
                else:
                    vx = idx_to_bead[off + L - 1, 5]; vy = idx_to_bead[off + L - 1, 6]; vz = idx_to_bead[off + L - 1, 7]

                de = get_energy_change_c(grid, type_grid, vx, vy, vz, nx, ny, nz,
                                         idx_to_bead[off, 1], interaction_table, LR_interaction_table,
                                         SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
                # angle change: new-config chain angle minus old-config chain angle
                de = de + (_chain_angle_mode(idx_to_bead, off, L, 1 if direction == 0 else 2, nx, ny, nz, angle_lookup)
                           - _chain_angle_mode(idx_to_bead, off, L, 0, 0, 0, 0, angle_lookup))

                if accept_or_reject(invtemp, energy, energy + de) == 1:
                    # apply: only the two end sites change in grid/type_grid; all positions shift
                    grid[vx, vy, vz] = 0; type_grid[vx, vy, vz] = 0
                    grid[nx, ny, nz] = chainID; type_grid[nx, ny, nz] = idx_to_bead[off, 2]
                    if direction == 0:                  # head-grow: bead k -> bead (k+1)'s old site
                        for k in range(0, L - 1):
                            idx_to_bead[off + k, 5] = idx_to_bead[off + k + 1, 5]
                            idx_to_bead[off + k, 6] = idx_to_bead[off + k + 1, 6]
                            idx_to_bead[off + k, 7] = idx_to_bead[off + k + 1, 7]
                        idx_to_bead[off + L - 1, 5] = nx; idx_to_bead[off + L - 1, 6] = ny; idx_to_bead[off + L - 1, 7] = nz
                    else:                               # tail-grow: bead k -> bead (k-1)'s old site
                        for k in range(L - 1, 0, -1):
                            idx_to_bead[off + k, 5] = idx_to_bead[off + k - 1, 5]
                            idx_to_bead[off + k, 6] = idx_to_bead[off + k - 1, 6]
                            idx_to_bead[off + k, 7] = idx_to_bead[off + k - 1, 7]
                        idx_to_bead[off, 5] = nx; idx_to_bead[off, 6] = ny; idx_to_bead[off, 7] = nz
                    energy = energy + de
                    accepted = accepted + 1
                # reject: get_energy_change_c self-reverted type_grid; idx untouched
                continue

            # ----- heteropolymer: decompose into L single-bead moves -----
            for k in range(L):
                bx[k] = idx_to_bead[off + k, 5]; by[k] = idx_to_bead[off + k, 6]; bz[k] = idx_to_bead[off + k, 7]

            de = 0
            if direction == 0:                          # head-grow: move bead L-1 -> new, then L-2 -> P[L-1], ...
                for k in range(L - 1, -1, -1):
                    ox = idx_to_bead[off + k, 5]; oy = idx_to_bead[off + k, 6]; oz = idx_to_bead[off + k, 7]
                    if k == L - 1:
                        tx = nx; ty = ny; tz = nz
                    else:
                        tx = bx[k + 1]; ty = by[k + 1]; tz = bz[k + 1]
                    de = de + get_energy_change_c(grid, type_grid, ox, oy, oz, tx, ty, tz,
                                                  idx_to_bead[off + k, 1], interaction_table,
                                                  LR_interaction_table, SLR_interaction_table,
                                                  XDIM, YDIM, ZDIM, hardwall)
                    newp[0] = tx; newp[1] = ty; newp[2] = tz
                    de = de + get_angle_energy_change_c(off + k, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move(grid, type_grid, idx_to_bead, off + k, ox, oy, oz, tx, ty, tz, chainID)
            else:                                       # tail-grow: move bead 0 -> new, then 1 -> P[0], ...
                for k in range(0, L):
                    ox = idx_to_bead[off + k, 5]; oy = idx_to_bead[off + k, 6]; oz = idx_to_bead[off + k, 7]
                    if k == 0:
                        tx = nx; ty = ny; tz = nz
                    else:
                        tx = bx[k - 1]; ty = by[k - 1]; tz = bz[k - 1]
                    de = de + get_energy_change_c(grid, type_grid, ox, oy, oz, tx, ty, tz,
                                                  idx_to_bead[off + k, 1], interaction_table,
                                                  LR_interaction_table, SLR_interaction_table,
                                                  XDIM, YDIM, ZDIM, hardwall)
                    newp[0] = tx; newp[1] = ty; newp[2] = tz
                    de = de + get_angle_energy_change_c(off + k, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move(grid, type_grid, idx_to_bead, off + k, ox, oy, oz, tx, ty, tz, chainID)

            if accept_or_reject(invtemp, energy, energy + de) == 1:
                energy = energy + de
                accepted = accepted + 1
            else:
                # revert: clear current (new) chain sites, then restore the buffered old chain
                for k in range(L):
                    cx = idx_to_bead[off + k, 5]; cy = idx_to_bead[off + k, 6]; cz = idx_to_bead[off + k, 7]
                    grid[cx, cy, cz] = 0; type_grid[cx, cy, cz] = 0
                for k in range(L):
                    grid[bx[k], by[k], bz[k]] = chainID
                    type_grid[bx[k], by[k], bz[k]] = idx_to_bead[off + k, 2]
                    idx_to_bead[off + k, 5] = bx[k]; idx_to_bead[off + k, 6] = by[k]; idx_to_bead[off + k, 7] = bz[k]
    finally:
        free(bx); free(by); free(bz)

    return (energy, accepted)


# =====================================================================
#   SLITHER (reptation) megamove kernel - 2D variant
#
#   Identical logic to mega_slither (3D) but on a 2D lattice, reusing the
#   bit-exact 2D primitives (get_energy_change_2D_c, get_angle_energy_change_2D_c,
#   single_bead_crank_2D_c). See mega_slither for the full design notes.
# =====================================================================

cdef inline long _triplet_angle_2D(int mx, int my,
                                   int lx, int ly,
                                   int rx, int ry,
                                   int ic,
                                   NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup) noexcept nogil:
    cdef int a0 = fix_angle_pbc_issues(mx - lx)
    cdef int a1 = fix_angle_pbc_issues(my - ly)
    cdef int b0 = fix_angle_pbc_issues(mx - rx)
    cdef int b1 = fix_angle_pbc_issues(my - ry)
    return angle_lookup[ic, a0 + 1, a1 + 1, b0 + 1, b1 + 1]


cdef inline void _smode_pos_2D(NUMPY_INT_TYPE_long[:, :] idx, int off, int L, int mode,
                               int nex, int ney, int k, int* o) noexcept nogil:
    if mode == 0:
        o[0] = idx[off + k, 5]; o[1] = idx[off + k, 6]
    elif mode == 1:                 # head-grow
        if k < L - 1:
            o[0] = idx[off + k + 1, 5]; o[1] = idx[off + k + 1, 6]
        else:
            o[0] = nex; o[1] = ney
    else:                           # tail-grow
        if k > 0:
            o[0] = idx[off + k - 1, 5]; o[1] = idx[off + k - 1, 6]
        else:
            o[0] = nex; o[1] = ney


cdef long _chain_angle_mode_2D(NUMPY_INT_TYPE_long[:, :] idx, int off, int L, int mode,
                               int nex, int ney,
                               NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup) noexcept nogil:
    cdef long e = 0
    cdef int k
    cdef int pl[2]
    cdef int pm[2]
    cdef int pr[2]
    for k in range(1, L - 1):
        if idx[off + k, 3] == 1:
            continue
        _smode_pos_2D(idx, off, L, mode, nex, ney, k - 1, pl)
        _smode_pos_2D(idx, off, L, mode, nex, ney, k, pm)
        _smode_pos_2D(idx, off, L, mode, nex, ney, k + 1, pr)
        e = e + _triplet_angle_2D(pm[0], pm[1], pl[0], pl[1], pr[0], pr[1], idx[off + k, 2], angle_lookup)
    return e


cdef inline void _apply_bead_move_2D(NUMPY_INT_TYPE[:, :] grid, NUMPY_INT_TYPE[:, :] type_grid,
                                     NUMPY_INT_TYPE_long[:, :] idx, int bead_index,
                                     int ox, int oy, int nx, int ny, int chainID) noexcept nogil:
    grid[ox, oy] = 0
    grid[nx, ny] = chainID
    type_grid[nx, ny] = type_grid[ox, oy]
    type_grid[ox, oy] = 0
    idx[bead_index, 5] = nx
    idx[bead_index, 6] = ny


@cython.wraparound(False)
@cython.boundscheck(False)
def mega_slither_2D(NUMPY_INT_TYPE[:, :] grid,
                    NUMPY_INT_TYPE[:, :] type_grid,
                    NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                    int[::1] chain_offset,
                    int[::1] chain_length,
                    int[::1] chain_homo,
                    int[::1] chain_selector,
                    NUMPY_INT_TYPE[:, :] interaction_table,
                    NUMPY_INT_TYPE[:, :] LR_interaction_table,
                    NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                    NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup,
                    long energy,
                    float invtemp,
                    int passed_seed,
                    int hardwall,
                    int max_chain_len):
    """2D slither (reptation) megamove. See mega_slither for semantics."""
    mc_seed(passed_seed)

    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int n_steps = chain_selector.shape[0]
    cdef int accepted = 0

    cdef int i, c, off, L, homo, chainID, direction, k
    cdef int ex, ey, nx, ny, vx, vy, ox, oy, tx, ty, cx, cy
    cdef long de
    cdef int newp[2]

    cdef int* bx = <int*> malloc(max_chain_len * sizeof(int))
    cdef int* by = <int*> malloc(max_chain_len * sizeof(int))
    if bx == NULL or by == NULL:
        if bx != NULL: free(bx)
        if by != NULL: free(by)
        raise MemoryError()

    try:
        for i in range(n_steps):
            c = chain_selector[i]
            off = chain_offset[c]
            L = chain_length[c]
            homo = chain_homo[c]
            chainID = idx_to_bead[off, 4]

            # ----- single-bead chain: local translation -----
            if L == 1:
                ox = idx_to_bead[off, 5]; oy = idx_to_bead[off, 6]
                if single_bead_crank_2D_c(ox, oy, grid, XDIM, YDIM, newp) == 1:
                    de = get_energy_change_2D_c(grid, type_grid, ox, oy, newp[0], newp[1],
                                                idx_to_bead[off, 1], interaction_table,
                                                LR_interaction_table, SLR_interaction_table,
                                                XDIM, YDIM, hardwall)
                    de = de + get_angle_energy_change_2D_c(off, idx_to_bead, newp, angle_lookup)
                    if accept_or_reject(invtemp, energy, energy + de) == 1:
                        _apply_bead_move_2D(grid, type_grid, idx_to_bead, off,
                                            ox, oy, newp[0], newp[1], chainID)
                        energy = energy + de
                        accepted = accepted + 1
                continue

            # ----- reptation -----
            direction = randint(0, 1)
            if direction == 0:
                ex = idx_to_bead[off + L - 1, 5]; ey = idx_to_bead[off + L - 1, 6]
            else:
                ex = idx_to_bead[off, 5]; ey = idx_to_bead[off, 6]

            nx = pbc_correction(ex + randint(0, 2) - 1, XDIM)
            ny = pbc_correction(ey + randint(0, 2) - 1, YDIM)

            if grid[nx, ny] > 0:
                continue
            if hardwall == 1:
                if straddle_pair(ex, ey, 0, nx, ny, 0) == 1:
                    continue

            if homo == 1:
                if direction == 0:
                    vx = idx_to_bead[off, 5]; vy = idx_to_bead[off, 6]
                else:
                    vx = idx_to_bead[off + L - 1, 5]; vy = idx_to_bead[off + L - 1, 6]

                de = get_energy_change_2D_c(grid, type_grid, vx, vy, nx, ny,
                                            idx_to_bead[off, 1], interaction_table,
                                            LR_interaction_table, SLR_interaction_table,
                                            XDIM, YDIM, hardwall)
                de = de + (_chain_angle_mode_2D(idx_to_bead, off, L, 1 if direction == 0 else 2, nx, ny, angle_lookup)
                           - _chain_angle_mode_2D(idx_to_bead, off, L, 0, 0, 0, angle_lookup))

                if accept_or_reject(invtemp, energy, energy + de) == 1:
                    grid[vx, vy] = 0; type_grid[vx, vy] = 0
                    grid[nx, ny] = chainID; type_grid[nx, ny] = idx_to_bead[off, 2]
                    if direction == 0:
                        for k in range(0, L - 1):
                            idx_to_bead[off + k, 5] = idx_to_bead[off + k + 1, 5]
                            idx_to_bead[off + k, 6] = idx_to_bead[off + k + 1, 6]
                        idx_to_bead[off + L - 1, 5] = nx; idx_to_bead[off + L - 1, 6] = ny
                    else:
                        for k in range(L - 1, 0, -1):
                            idx_to_bead[off + k, 5] = idx_to_bead[off + k - 1, 5]
                            idx_to_bead[off + k, 6] = idx_to_bead[off + k - 1, 6]
                        idx_to_bead[off, 5] = nx; idx_to_bead[off, 6] = ny
                    energy = energy + de
                    accepted = accepted + 1
                continue

            # ----- heteropolymer: decompose into L single-bead moves -----
            for k in range(L):
                bx[k] = idx_to_bead[off + k, 5]; by[k] = idx_to_bead[off + k, 6]

            de = 0
            if direction == 0:
                for k in range(L - 1, -1, -1):
                    ox = idx_to_bead[off + k, 5]; oy = idx_to_bead[off + k, 6]
                    if k == L - 1:
                        tx = nx; ty = ny
                    else:
                        tx = bx[k + 1]; ty = by[k + 1]
                    de = de + get_energy_change_2D_c(grid, type_grid, ox, oy, tx, ty,
                                                     idx_to_bead[off + k, 1], interaction_table,
                                                     LR_interaction_table, SLR_interaction_table,
                                                     XDIM, YDIM, hardwall)
                    newp[0] = tx; newp[1] = ty
                    de = de + get_angle_energy_change_2D_c(off + k, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move_2D(grid, type_grid, idx_to_bead, off + k, ox, oy, tx, ty, chainID)
            else:
                for k in range(0, L):
                    ox = idx_to_bead[off + k, 5]; oy = idx_to_bead[off + k, 6]
                    if k == 0:
                        tx = nx; ty = ny
                    else:
                        tx = bx[k - 1]; ty = by[k - 1]
                    de = de + get_energy_change_2D_c(grid, type_grid, ox, oy, tx, ty,
                                                     idx_to_bead[off + k, 1], interaction_table,
                                                     LR_interaction_table, SLR_interaction_table,
                                                     XDIM, YDIM, hardwall)
                    newp[0] = tx; newp[1] = ty
                    de = de + get_angle_energy_change_2D_c(off + k, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move_2D(grid, type_grid, idx_to_bead, off + k, ox, oy, tx, ty, chainID)

            if accept_or_reject(invtemp, energy, energy + de) == 1:
                energy = energy + de
                accepted = accepted + 1
            else:
                for k in range(L):
                    cx = idx_to_bead[off + k, 5]; cy = idx_to_bead[off + k, 6]
                    grid[cx, cy] = 0; type_grid[cx, cy] = 0
                for k in range(L):
                    grid[bx[k], by[k]] = chainID
                    type_grid[bx[k], by[k]] = idx_to_bead[off + k, 2]
                    idx_to_bead[off + k, 5] = bx[k]; idx_to_bead[off + k, 6] = by[k]
    finally:
        free(bx); free(by)

    return (energy, accepted)


# =====================================================================
#
#   PULL (cooperative reptation) megamove kernel
#
#   A pull move displaces an interior bead i to a neighbouring empty site and
#   then cooperatively "pulls" the following beads into the just-vacated sites
#   until chain connectivity is restored - local reptation of a sub-segment.
#   Because it threads through occupied space one bead at a time it rearranges
#   chains in DENSE systems where rigid translate/rotate/pivot would clash.
#
#   Detailed balance: the move is reversible (the reverse of a C-ward pull whose
#   cascade ends at bead k is an N-ward pull at bead k, and vice-versa - the
#   forward "did-not-stop" conditions force the reverse cascade to retrace the
#   same beads). Acceptance is Metropolis-Hastings with the proposal-multiplicity
#   ratio nF/nR, where nF/nR are the forward/reverse first-target counts computed
#   by the SAME predicate (pull_first_targets) so they can never diverge. Cascades
#   that reach a terminus are rejected (their reverse has no valid anchor).
#   Energy is the telescoping per-bead-move sum (the proven mega_slither approach).
# =====================================================================

cdef inline int accept_or_reject_ratio(float invtemp, long old_energy, long new_energy,
                                       int nF, int nR) noexcept nogil:
    # Metropolis-Hastings: accept iff rand < (nF/nR) * exp(-(dE)*invtemp)
    cdef double acc
    if nR <= 0:
        return 0
    acc = (<double>nF / <double>nR) * exp(-(<double>(new_energy - old_energy)) * invtemp)
    if acc >= 1.0:
        return 1
    if (<double>mc_rand() / <double>PRNG_MAX) < acc:
        return 1
    return 0


cdef inline int cheby_adjacent(int ax, int ay, int az, int bx, int by, int bz,
                               int XDIM, int YDIM, int ZDIM, int hardwall) noexcept nogil:
    # 1 iff the two sites are Chebyshev-1 neighbours (minimal-image under PBC;
    # literal coordinate distance under hardwall, so wrapped sites are NOT adjacent)
    cdef int dx = ax - bx
    cdef int dy = ay - by
    cdef int dz = az - bz
    if dx < 0: dx = -dx
    if dy < 0: dy = -dy
    if dz < 0: dz = -dz
    if hardwall == 0:
        if XDIM - dx < dx: dx = XDIM - dx
        if YDIM - dy < dy: dy = YDIM - dy
        if ZDIM - dz < dz: dz = ZDIM - dz
    if dx <= 1 and dy <= 1 and dz <= 1:
        return 1
    return 0


cdef int pull_first_targets(NUMPY_INT_TYPE[:, :, :] grid,
                            int anchor_x, int anchor_y, int anchor_z,
                            int mov_x, int mov_y, int mov_z,
                            int hardwall, int XDIM, int YDIM, int ZDIM,
                            int* out_buf) noexcept nogil:
    # empty sites that are Chebyshev-1 neighbours of BOTH the moving bead and the
    # anchor (so the moved bead stays bonded to the anchor and the next pulled
    # bead can occupy the moving bead's old site). Writes (x,y,z) triples into
    # out_buf, returns the count. The current moving-bead site is occupied so it
    # is never returned. (Assumes each box dim >= 3 so the 27 offsets are distinct.)
    cdef int count = 0
    cdef int dx, dy, dz, tx, ty, tz
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            for dz in range(-1, 2):
                tx = pbc_correction(mov_x + dx, XDIM)
                ty = pbc_correction(mov_y + dy, YDIM)
                tz = pbc_correction(mov_z + dz, ZDIM)
                if grid[tx, ty, tz] != 0:
                    continue
                if cheby_adjacent(tx, ty, tz, mov_x, mov_y, mov_z, XDIM, YDIM, ZDIM, hardwall) == 0:
                    continue
                if cheby_adjacent(tx, ty, tz, anchor_x, anchor_y, anchor_z, XDIM, YDIM, ZDIM, hardwall) == 0:
                    continue
                out_buf[3 * count] = tx
                out_buf[3 * count + 1] = ty
                out_buf[3 * count + 2] = tz
                count = count + 1
    return count


cdef inline void _revert_segment(NUMPY_INT_TYPE[:, :, :] grid, NUMPY_INT_TYPE[:, :, :] type_grid,
                                 NUMPY_INT_TYPE_long[:, :] idx, int off, int lo, int hi,
                                 int* bx, int* by, int* bz, int chainID) noexcept nogil:
    # restore beads [lo..hi] to their buffered old positions (clear current sites
    # first, then write old, to handle overlap) - mirrors the mega_slither revert
    cdef int j, ox, oy, oz
    for j in range(lo, hi + 1):
        ox = idx[off + j, 5]; oy = idx[off + j, 6]; oz = idx[off + j, 7]
        grid[ox, oy, oz] = 0; type_grid[ox, oy, oz] = 0
    for j in range(lo, hi + 1):
        grid[bx[j], by[j], bz[j]] = chainID
        type_grid[bx[j], by[j], bz[j]] = idx[off + j, 2]
        idx[off + j, 5] = bx[j]; idx[off + j, 6] = by[j]; idx[off + j, 7] = bz[j]


@cython.wraparound(False)
@cython.boundscheck(False)
def mega_pull(NUMPY_INT_TYPE[:, :, :] grid,
              NUMPY_INT_TYPE[:, :, :] type_grid,
              NUMPY_INT_TYPE_long[:, :] idx_to_bead,
              int[::1] chain_offset,
              int[::1] chain_length,
              int[::1] chain_homo,
              int[::1] chain_selector,
              NUMPY_INT_TYPE[:, :] interaction_table,
              NUMPY_INT_TYPE[:, :] LR_interaction_table,
              NUMPY_INT_TYPE[:, :] SLR_interaction_table,
              NUMPY_INT_TYPE[:, :, :, :, :, :, :] angle_lookup,
              long energy,
              float invtemp,
              int passed_seed,
              int hardwall,
              int max_chain_len):
    """
    3D pull (cooperative reptation) megamove. chain_selector lists chain indices
    to pull (each pull-eligible chain appears PULL_SUBSTEPS times, shuffled).
    Mutates grid / type_grid / idx_to_bead in place. Returns (energy, accepted).
    chain_homo is accepted for signature symmetry with mega_slither but unused.
    """
    mc_seed(passed_seed)

    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int ZDIM = grid.shape[2]
    cdef int n_steps = chain_selector.shape[0]
    cdef int accepted = 0

    cdef int step, c, off, L, chainID, i, j, k, direction, nF, nR, r, restored, lo, hi
    cdef long de
    cdef int sx, sy, sz, ox, oy, oz, tx, ty, tz, pnx, pny, pnz
    cdef int newp[3]
    cdef int tbuf[81]

    cdef int* bx = <int*> malloc(max_chain_len * sizeof(int))
    cdef int* by = <int*> malloc(max_chain_len * sizeof(int))
    cdef int* bz = <int*> malloc(max_chain_len * sizeof(int))
    if bx == NULL or by == NULL or bz == NULL:
        if bx != NULL: free(bx)
        if by != NULL: free(by)
        if bz != NULL: free(bz)
        raise MemoryError()

    try:
        for step in range(n_steps):
            c = chain_selector[step]
            off = chain_offset[c]
            L = chain_length[c]
            if L < 3:
                continue
            chainID = idx_to_bead[off, 4]

            for j in range(L):
                bx[j] = idx_to_bead[off + j, 5]
                by[j] = idx_to_bead[off + j, 6]
                bz[j] = idx_to_bead[off + j, 7]

            i = 1 + randint(0, L - 3)          # interior bead in [1, L-2]
            direction = randint(0, 1)          # 0 = C-ward (anchor i-1), 1 = N-ward (anchor i+1)
            restored = 0
            nR = 0

            if direction == 0:
                nF = pull_first_targets(grid, bx[i - 1], by[i - 1], bz[i - 1],
                                        bx[i], by[i], bz[i], hardwall, XDIM, YDIM, ZDIM, tbuf)
                if nF == 0:
                    continue
                r = randint(0, nF - 1)
                sx = tbuf[3 * r]; sy = tbuf[3 * r + 1]; sz = tbuf[3 * r + 2]
                de = get_energy_change_c(grid, type_grid, bx[i], by[i], bz[i], sx, sy, sz,
                                         idx_to_bead[off + i, 1], interaction_table, LR_interaction_table,
                                         SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
                newp[0] = sx; newp[1] = sy; newp[2] = sz
                de = de + get_angle_energy_change_c(off + i, idx_to_bead, newp, angle_lookup)
                _apply_bead_move(grid, type_grid, idx_to_bead, off + i, bx[i], by[i], bz[i], sx, sy, sz, chainID)
                pnx = sx; pny = sy; pnz = sz
                k = i
                for j in range(i + 1, L):
                    if cheby_adjacent(pnx, pny, pnz, bx[j], by[j], bz[j], XDIM, YDIM, ZDIM, hardwall) == 1:
                        k = j - 1; restored = 1; break
                    ox = bx[j]; oy = by[j]; oz = bz[j]
                    tx = bx[j - 1]; ty = by[j - 1]; tz = bz[j - 1]
                    de = de + get_energy_change_c(grid, type_grid, ox, oy, oz, tx, ty, tz,
                                                  idx_to_bead[off + j, 1], interaction_table, LR_interaction_table,
                                                  SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
                    newp[0] = tx; newp[1] = ty; newp[2] = tz
                    de = de + get_angle_energy_change_c(off + j, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move(grid, type_grid, idx_to_bead, off + j, ox, oy, oz, tx, ty, tz, chainID)
                    pnx = tx; pny = ty; pnz = tz
                    k = j
                lo = i; hi = k
                if restored == 1:
                    nR = pull_first_targets(grid, bx[k + 1], by[k + 1], bz[k + 1],
                                            idx_to_bead[off + k, 5], idx_to_bead[off + k, 6], idx_to_bead[off + k, 7],
                                            hardwall, XDIM, YDIM, ZDIM, tbuf)
            else:
                nF = pull_first_targets(grid, bx[i + 1], by[i + 1], bz[i + 1],
                                        bx[i], by[i], bz[i], hardwall, XDIM, YDIM, ZDIM, tbuf)
                if nF == 0:
                    continue
                r = randint(0, nF - 1)
                sx = tbuf[3 * r]; sy = tbuf[3 * r + 1]; sz = tbuf[3 * r + 2]
                de = get_energy_change_c(grid, type_grid, bx[i], by[i], bz[i], sx, sy, sz,
                                         idx_to_bead[off + i, 1], interaction_table, LR_interaction_table,
                                         SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
                newp[0] = sx; newp[1] = sy; newp[2] = sz
                de = de + get_angle_energy_change_c(off + i, idx_to_bead, newp, angle_lookup)
                _apply_bead_move(grid, type_grid, idx_to_bead, off + i, bx[i], by[i], bz[i], sx, sy, sz, chainID)
                pnx = sx; pny = sy; pnz = sz
                k = i
                for j in range(i - 1, -1, -1):
                    if cheby_adjacent(pnx, pny, pnz, bx[j], by[j], bz[j], XDIM, YDIM, ZDIM, hardwall) == 1:
                        k = j + 1; restored = 1; break
                    ox = bx[j]; oy = by[j]; oz = bz[j]
                    tx = bx[j + 1]; ty = by[j + 1]; tz = bz[j + 1]
                    de = de + get_energy_change_c(grid, type_grid, ox, oy, oz, tx, ty, tz,
                                                  idx_to_bead[off + j, 1], interaction_table, LR_interaction_table,
                                                  SLR_interaction_table, XDIM, YDIM, ZDIM, hardwall)
                    newp[0] = tx; newp[1] = ty; newp[2] = tz
                    de = de + get_angle_energy_change_c(off + j, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move(grid, type_grid, idx_to_bead, off + j, ox, oy, oz, tx, ty, tz, chainID)
                    pnx = tx; pny = ty; pnz = tz
                    k = j
                lo = k; hi = i
                if restored == 1:
                    nR = pull_first_targets(grid, bx[k - 1], by[k - 1], bz[k - 1],
                                            idx_to_bead[off + k, 5], idx_to_bead[off + k, 6], idx_to_bead[off + k, 7],
                                            hardwall, XDIM, YDIM, ZDIM, tbuf)

            if restored == 0:
                # cascade reached a terminus -> reject (no valid reverse anchor)
                _revert_segment(grid, type_grid, idx_to_bead, off, lo, hi, bx, by, bz, chainID)
                continue

            # Metropolis-Hastings with the nF/nR proposal-multiplicity correction
            if accept_or_reject_ratio(invtemp, energy, energy + de, nF, nR) == 1:
                energy = energy + de
                accepted = accepted + 1
            else:
                _revert_segment(grid, type_grid, idx_to_bead, off, lo, hi, bx, by, bz, chainID)
    finally:
        free(bx); free(by); free(bz)

    return (energy, accepted)


# =====================================================================
#   PULL megamove kernel - 2D variant (mirrors mega_pull; see its notes)
# =====================================================================

cdef inline int cheby_adjacent_2D(int ax, int ay, int bx, int by,
                                  int XDIM, int YDIM, int hardwall) noexcept nogil:
    cdef int dx = ax - bx
    cdef int dy = ay - by
    if dx < 0: dx = -dx
    if dy < 0: dy = -dy
    if hardwall == 0:
        if XDIM - dx < dx: dx = XDIM - dx
        if YDIM - dy < dy: dy = YDIM - dy
    if dx <= 1 and dy <= 1:
        return 1
    return 0


cdef int pull_first_targets_2D(NUMPY_INT_TYPE[:, :] grid,
                               int anchor_x, int anchor_y,
                               int mov_x, int mov_y,
                               int hardwall, int XDIM, int YDIM,
                               int* out_buf) noexcept nogil:
    cdef int count = 0
    cdef int dx, dy, tx, ty
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            tx = pbc_correction(mov_x + dx, XDIM)
            ty = pbc_correction(mov_y + dy, YDIM)
            if grid[tx, ty] != 0:
                continue
            if cheby_adjacent_2D(tx, ty, mov_x, mov_y, XDIM, YDIM, hardwall) == 0:
                continue
            if cheby_adjacent_2D(tx, ty, anchor_x, anchor_y, XDIM, YDIM, hardwall) == 0:
                continue
            out_buf[2 * count] = tx
            out_buf[2 * count + 1] = ty
            count = count + 1
    return count


cdef inline void _revert_segment_2D(NUMPY_INT_TYPE[:, :] grid, NUMPY_INT_TYPE[:, :] type_grid,
                                    NUMPY_INT_TYPE_long[:, :] idx, int off, int lo, int hi,
                                    int* bx, int* by, int chainID) noexcept nogil:
    cdef int j, ox, oy
    for j in range(lo, hi + 1):
        ox = idx[off + j, 5]; oy = idx[off + j, 6]
        grid[ox, oy] = 0; type_grid[ox, oy] = 0
    for j in range(lo, hi + 1):
        grid[bx[j], by[j]] = chainID
        type_grid[bx[j], by[j]] = idx[off + j, 2]
        idx[off + j, 5] = bx[j]; idx[off + j, 6] = by[j]


@cython.wraparound(False)
@cython.boundscheck(False)
def mega_pull_2D(NUMPY_INT_TYPE[:, :] grid,
                 NUMPY_INT_TYPE[:, :] type_grid,
                 NUMPY_INT_TYPE_long[:, :] idx_to_bead,
                 int[::1] chain_offset,
                 int[::1] chain_length,
                 int[::1] chain_homo,
                 int[::1] chain_selector,
                 NUMPY_INT_TYPE[:, :] interaction_table,
                 NUMPY_INT_TYPE[:, :] LR_interaction_table,
                 NUMPY_INT_TYPE[:, :] SLR_interaction_table,
                 NUMPY_INT_TYPE[:, :, :, :, :] angle_lookup,
                 long energy,
                 float invtemp,
                 int passed_seed,
                 int hardwall,
                 int max_chain_len):
    """2D pull megamove. See mega_pull for semantics."""
    mc_seed(passed_seed)

    cdef int XDIM = grid.shape[0]
    cdef int YDIM = grid.shape[1]
    cdef int n_steps = chain_selector.shape[0]
    cdef int accepted = 0

    cdef int step, c, off, L, chainID, i, j, k, direction, nF, nR, r, restored, lo, hi
    cdef long de
    cdef int sx, sy, ox, oy, tx, ty, pnx, pny
    cdef int newp[2]
    cdef int tbuf[18]

    cdef int* bx = <int*> malloc(max_chain_len * sizeof(int))
    cdef int* by = <int*> malloc(max_chain_len * sizeof(int))
    if bx == NULL or by == NULL:
        if bx != NULL: free(bx)
        if by != NULL: free(by)
        raise MemoryError()

    try:
        for step in range(n_steps):
            c = chain_selector[step]
            off = chain_offset[c]
            L = chain_length[c]
            if L < 3:
                continue
            chainID = idx_to_bead[off, 4]

            for j in range(L):
                bx[j] = idx_to_bead[off + j, 5]
                by[j] = idx_to_bead[off + j, 6]

            i = 1 + randint(0, L - 3)
            direction = randint(0, 1)
            restored = 0
            nR = 0

            if direction == 0:
                nF = pull_first_targets_2D(grid, bx[i - 1], by[i - 1], bx[i], by[i],
                                           hardwall, XDIM, YDIM, tbuf)
                if nF == 0:
                    continue
                r = randint(0, nF - 1)
                sx = tbuf[2 * r]; sy = tbuf[2 * r + 1]
                de = get_energy_change_2D_c(grid, type_grid, bx[i], by[i], sx, sy,
                                            idx_to_bead[off + i, 1], interaction_table, LR_interaction_table,
                                            SLR_interaction_table, XDIM, YDIM, hardwall)
                newp[0] = sx; newp[1] = sy
                de = de + get_angle_energy_change_2D_c(off + i, idx_to_bead, newp, angle_lookup)
                _apply_bead_move_2D(grid, type_grid, idx_to_bead, off + i, bx[i], by[i], sx, sy, chainID)
                pnx = sx; pny = sy
                k = i
                for j in range(i + 1, L):
                    if cheby_adjacent_2D(pnx, pny, bx[j], by[j], XDIM, YDIM, hardwall) == 1:
                        k = j - 1; restored = 1; break
                    ox = bx[j]; oy = by[j]
                    tx = bx[j - 1]; ty = by[j - 1]
                    de = de + get_energy_change_2D_c(grid, type_grid, ox, oy, tx, ty,
                                                     idx_to_bead[off + j, 1], interaction_table, LR_interaction_table,
                                                     SLR_interaction_table, XDIM, YDIM, hardwall)
                    newp[0] = tx; newp[1] = ty
                    de = de + get_angle_energy_change_2D_c(off + j, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move_2D(grid, type_grid, idx_to_bead, off + j, ox, oy, tx, ty, chainID)
                    pnx = tx; pny = ty
                    k = j
                lo = i; hi = k
                if restored == 1:
                    nR = pull_first_targets_2D(grid, bx[k + 1], by[k + 1],
                                               idx_to_bead[off + k, 5], idx_to_bead[off + k, 6],
                                               hardwall, XDIM, YDIM, tbuf)
            else:
                nF = pull_first_targets_2D(grid, bx[i + 1], by[i + 1], bx[i], by[i],
                                           hardwall, XDIM, YDIM, tbuf)
                if nF == 0:
                    continue
                r = randint(0, nF - 1)
                sx = tbuf[2 * r]; sy = tbuf[2 * r + 1]
                de = get_energy_change_2D_c(grid, type_grid, bx[i], by[i], sx, sy,
                                            idx_to_bead[off + i, 1], interaction_table, LR_interaction_table,
                                            SLR_interaction_table, XDIM, YDIM, hardwall)
                newp[0] = sx; newp[1] = sy
                de = de + get_angle_energy_change_2D_c(off + i, idx_to_bead, newp, angle_lookup)
                _apply_bead_move_2D(grid, type_grid, idx_to_bead, off + i, bx[i], by[i], sx, sy, chainID)
                pnx = sx; pny = sy
                k = i
                for j in range(i - 1, -1, -1):
                    if cheby_adjacent_2D(pnx, pny, bx[j], by[j], XDIM, YDIM, hardwall) == 1:
                        k = j + 1; restored = 1; break
                    ox = bx[j]; oy = by[j]
                    tx = bx[j + 1]; ty = by[j + 1]
                    de = de + get_energy_change_2D_c(grid, type_grid, ox, oy, tx, ty,
                                                     idx_to_bead[off + j, 1], interaction_table, LR_interaction_table,
                                                     SLR_interaction_table, XDIM, YDIM, hardwall)
                    newp[0] = tx; newp[1] = ty
                    de = de + get_angle_energy_change_2D_c(off + j, idx_to_bead, newp, angle_lookup)
                    _apply_bead_move_2D(grid, type_grid, idx_to_bead, off + j, ox, oy, tx, ty, chainID)
                    pnx = tx; pny = ty
                    k = j
                lo = k; hi = i
                if restored == 1:
                    nR = pull_first_targets_2D(grid, bx[k - 1], by[k - 1],
                                               idx_to_bead[off + k, 5], idx_to_bead[off + k, 6],
                                               hardwall, XDIM, YDIM, tbuf)

            if restored == 0:
                _revert_segment_2D(grid, type_grid, idx_to_bead, off, lo, hi, bx, by, chainID)
                continue

            # Metropolis-Hastings with the nF/nR proposal-multiplicity correction
            if accept_or_reject_ratio(invtemp, energy, energy + de, nF, nR) == 1:
                energy = energy + de
                accepted = accepted + 1
            else:
                _revert_segment_2D(grid, type_grid, idx_to_bead, off, lo, hi, bx, by, chainID)
    finally:
        free(bx); free(by)

    return (energy, accepted)
