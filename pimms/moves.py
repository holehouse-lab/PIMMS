## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import random
import math
import itertools
import bisect
import warnings
from functools import lru_cache
import numpy as np
import copy
import sys

from . import lattice_utils
from . import cluster_utils
from . import numpy_utils
from . import mega_crank_fast
from . import crankshaft_list_functions
from . import IO_utils

from .latticeExceptions import ClusterSizeThresholdException
from .moveEvent import MoveEvent


def _vmmc_link_probability(beta, delta_energy):
    """Return ``max(0, 1 - exp(-beta * delta_energy))`` stably.

    The direct expression overflows when ``delta_energy`` is sufficiently
    negative, even though the mathematically correct (clamped) answer is simply
    zero. ``expm1`` also retains precision for weak positive links where
    subtracting a value close to one would lose significant digits.

    Parameters
    ----------
    beta : float
        Inverse temperature (``acceptanceObject.invtemp``).

    delta_energy : float
        Energy change of the link under the virtual move (the energy that would
        be paid by moving one partner alone). Only positive values give a
        non-zero link probability.

    Returns
    -------
    float
        The VMMC link formation probability, in ``[0, 1)``; exactly ``0.0``
        whenever ``delta_energy`` is zero or negative.
    """
    if delta_energy <= 0.0:
        return 0.0
    return -math.expm1(-beta * delta_energy)


def _outside_box(raw_position, dimensions):
    """True if a raw (unwrapped) lattice coordinate lies outside ``[0, dim)`` on
    any axis - i.e. under HARDWALL the bead would pass through a wall.

    Parameters
    ----------
    raw_position : list of int
        An unwrapped lattice coordinate (one entry per dimension), i.e. a
        position BEFORE any periodic wrapping has been applied.

    dimensions : list of int
        The lattice box dimensions, one entry per dimension.

    Returns
    -------
    bool
        True if any component of ``raw_position`` lies outside ``[0, dim)``,
        False otherwise.
    """
    for dim in range(len(dimensions)):
        if raw_position[dim] < 0 or raw_position[dim] >= dimensions[dim]:
            return True
    return False


def _is_identity_proposal(new_positions, original_positions, num_dimensions):
    """True if a rigid move's proposal puts every bead back where it started.

    Some rigid-body draws map the body exactly onto itself - in 3D, the three
    rotations of an axis-aligned straight chain about its own axis, and every
    rotation of an isolated single-bead cluster. Those proposals are rejected as
    null moves rather than returned as successes; see the comment in
    ``chain_rotate`` for why that is ensemble-neutral and why it matters for
    ACCEPTANCE.dat.

    The comparison is over ORDERED positions on purpose. A 180 degree rotation of
    an axis-aligned rod about its centre bead reverses the bead order: the set of
    occupied sites is unchanged, but for a labelled (heteropolymer) chain the
    reversed chain is a different microstate and must stay a real move. Comparing
    occupied sites would delete that transition from the move set.

    Parameters
    ----------
    new_positions : list of list of int
        The proposed positions, in bead order.

    original_positions : list of list of int
        The positions the beads currently occupy, in the same bead order.

    num_dimensions : int
        Number of spatial dimensions (2 or 3).

    Returns
    -------
    bool
        True if the two position lists agree bead for bead and axis for axis.
    """
    for new, old in zip(new_positions, original_positions):
        for dim in range(num_dimensions):
            if int(new[dim]) != int(old[dim]):
                return False
    return True


def _normalized_frozen_chains(latticeObject, frozen_chains=()):
    """Return the valid, unique chain IDs requested by the freeze mechanism.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose ``chains`` dictionary defines which chainIDs actually
        exist; anything else in ``frozen_chains`` is discarded.

    frozen_chains : sequence of int, optional
        The chainIDs the caller wants frozen. Duplicates and IDs that are not on
        the lattice are dropped. Default is ``()``.

    Returns
    -------
    tuple of int
        The valid frozen chainIDs, deduplicated and sorted ascending, so the
        result is hashable and order-stable across calls.
    """
    return tuple(sorted({chainID for chainID in frozen_chains
                         if chainID in latticeObject.chains}))


@lru_cache(maxsize=8)
def _vmmc_offset_shell(num_dimensions, radius):
    """Cache VMMC neighbour offsets together with their Chebyshev radius.

    Parameters
    ----------
    num_dimensions : int
        Number of lattice dimensions (2 or 3).

    radius : int
        Chebyshev half-width of the shell to build (1 for the short-range shell,
        3 for the shell that also covers LR and SLR interactions).

    Returns
    -------
    tuple of (tuple of int, int)
        One ``(offset, chebyshev)`` pair per neighbour site within the shell,
        where ``offset`` is the per-dimension displacement and ``chebyshev`` its
        Chebyshev distance (1, 2 or 3). The zero offset is excluded.
    """
    shell = []
    for delta in itertools.product(range(-radius, radius + 1), repeat=num_dimensions):
        chebyshev = max(abs(value) for value in delta)
        if chebyshev:
            shell.append((delta, chebyshev))
    return tuple(shell)


def _vmmc_log_acceptance(beta, delta_energy, cluster, formed_links, failed_links):
    """
    Log of the Metropolis-Hastings ratio for one VMMC cluster realisation.

    The recruitment is the proposal, so the ratio multiplies the Boltzmann factor
    ``exp(-beta * dE)`` by the reverse/forward generation probability of the
    realisation: every tested link contributes exactly one factor, ``p_r / p_f``
    if it formed and ``(1 - p_r) / (1 - p_f)`` if it failed. Which reverse
    probability ``p_r`` applies to a failed link depends on where its partner
    ended up: a partner OUTSIDE the final cluster (a boundary link) does not move
    in the reverse move either, so the post-move link probability applies; a
    partner INSIDE the final cluster (a "frustrated" link - it failed here but was
    recruited through another chain) translates with the cluster in the reverse
    move, so the same relative-displacement probability as a formed link applies.
    Dropping the frustrated-link factor (the pre-1.0.8 "boundary-only" rule) left
    the generation probabilities of a realisation and of its mirror unbalanced;
    exact enumeration of three monomers gave forward/reverse transition ratios of
    up to 1.36 between equal-energy states.

    Parameters
    ----------
    beta : float
        Inverse temperature.

    delta_energy : float
        Exact total energy change of the rigid cluster translation.

    cluster : set of int
        chainIDs in the final cluster.

    formed_links : list of (float, float)
        ``(p_f, p_r)`` for every link that formed.

    failed_links : list of (int, float, float, float)
        ``(j, p_f, p_r_boundary, p_r_internal)`` for every tested link that did
        not form, ``j`` being the partner chainID.

    Returns
    -------
    float or None
        ``log`` of the acceptance ratio (before the ``min(1, .)``), or ``None``
        when the reverse realisation has zero probability (a formed link with
        ``p_r == 0`` or a failed link with ``p_r == 1``), i.e. the move must be
        rejected.
    """
    log_ratio = -beta * delta_energy
    for (p_f, p_r) in formed_links:
        if p_r <= 0.0:
            return None
        log_ratio += math.log(p_r) - math.log(p_f)
    for (j, p_f, p_r_boundary, p_r_internal) in failed_links:
        p_r = p_r_internal if j in cluster else p_r_boundary
        q_r = 1.0 - p_r
        if q_r <= 0.0:
            return None
        log_ratio += math.log(q_r) - math.log(1.0 - p_f)
    return log_ratio


@lru_cache(maxsize=128)
def _vmmc_cutoff_cdf(cap):
    """Cache the harmonic CDF used for VMMC cluster-size cutoffs.

    Parameters
    ----------
    cap : int
        Largest cluster-size cutoff that can be drawn (must be >= 1). The
        distribution ``Q(n_c)`` is proportional to ``1/n_c`` over ``[1, cap]``.

    Returns
    -------
    tuple of float
        The cumulative distribution of ``Q``, i.e. element ``k - 1`` is the
        probability that the drawn cutoff is <= ``k``. The final element is 1.0.
    """
    total = sum(1.0 / k for k in range(1, cap + 1))
    running = 0.0
    cdf = []
    for k in range(1, cap + 1):
        running += (1.0 / k) / total
        cdf.append(running)
    return tuple(cdf)


def parallel_chain_metadata(latticeObject):
    """Per-chain arrays in the sorted-chainID order used by the slither/pull kernels.

    A chain is *homo* (eligible for the kernels' O(1) energy path) only if EVERY bead
    shares one intcode AND one LR flag - the kernels read bead 0's values for the
    whole chain, so both must be uniform.

    Parameters
    ----------
    latticeObject : Lattice
        The lattice whose chains are being described. Its ``crankshaft_lists``
        matrix is refreshed from the current chain positions as a side effect.

    Returns
    -------
    tuple
        ``(idx_to_bead, sorted_chains, chain_offset, chain_length, chain_homo)``
        where ``idx_to_bead`` is the ``int64`` bead bookkeeping matrix of shape
        ``(num_beads, 7)`` in 2D or ``(num_beads, 8)`` in 3D, ``sorted_chains``
        the ascending list of chainIDs, and ``chain_offset``, ``chain_length``
        and ``chain_homo`` are ``int32`` arrays of length ``num_chains`` holding,
        per chain (in that same order), the row index of its first bead in
        ``idx_to_bead``, its number of beads, and 1 if it is homopolymeric
        (uniform intcode and LR flag) else 0.
    """
    idx_to_bead = crankshaft_list_functions.update_idx_to_bead(latticeObject)
    layout = crankshaft_list_functions.chain_layout(latticeObject)
    return (idx_to_bead, list(layout.sorted_ids), layout.offset32.copy(),
            layout.length32.copy(), layout.homo.copy())


def parallel_chain_partition(idx_to_bead, chain_offset, chain_length, dimensions, has_LR,
                             chain_homo=None, cap_mode='all', frozen_chains=()):
    """Decide, from chain LENGTHS alone, which chains the parallel kernel may move.

    The parallel whole-chain kernels only move a chain when every one of its beads
    lies inside a single block interior (the block minus a width-W frozen halo on
    each face), and they reject any sub-move whose new site leaves that interior.
    A chain that does not fit is silently frozen for the sweep. So the whole-chain
    megamoves have to decide which kernel moves which chain - and that decision
    MUST NOT look at the current configuration.

    That is what this function is for. A chain of ``n`` beads spans at most ``n``
    sites on any axis (consecutive beads differ by at most one per axis), so
    ``length <= interior`` is a sufficient condition for fitting that reads only
    the chain length. Everything else it consults - the 512-bead per-thread buffer
    cap, the frozen set, the box, ``has_LR`` - is likewise fixed for the run, so
    the returned partition is a run constant: it changes only when the box changes
    (a resized equilibration), which is a scheduled, state-independent event.

    Chains that come back True are handed to the parallel kernel; the rest are
    handed to the serial kernel in the same megamove, so nothing is ever skipped
    and nothing is ever gated on the current configuration.

    Parameters
    ----------
    idx_to_bead : numpy.ndarray
        The ``int64`` bead bookkeeping matrix, shape ``(num_beads, 7)`` in 2D or
        ``(num_beads, 8)`` in 3D. Only column 4 (the chainID of each chain's
        first bead) is read, to resolve ``frozen_chains``. The bead COORDINATES
        are deliberately not read - that is the whole point of this function.

    chain_offset : numpy.ndarray
        ``int32`` array of length ``num_chains`` holding the row of
        ``idx_to_bead`` at which each chain starts.

    chain_length : numpy.ndarray
        ``int32`` array of length ``num_chains`` holding the number of beads in
        each chain.

    dimensions : list of int
        The lattice box dimensions (2 or 3 entries), which set the block
        decomposition.

    has_LR : bool
        Whether any bead engages in long-range interactions. This widens the
        frozen halo W and so shrinks the interior.

    chain_homo : numpy.ndarray or None, optional
        ``int32`` array of length ``num_chains``, 1 where a chain is
        homopolymeric. Only consulted when ``cap_mode`` is ``'hetero'``. Default
        is None.

    cap_mode : str, optional
        Which chains the kernels' 512-bead per-thread buffer applies to:
        ``'all'`` (the pull kernels) or ``'hetero'`` (the slither kernels, whose
        homopolymer path needs no per-bead buffer). Default is ``'all'``.

    frozen_chains : sequence of int, optional
        chainIDs that are permanently frozen. These come back False: they are not
        meant to move at all, and reporting them on the serial side keeps them out
        of the parallel kernel's movable set through the same frozen mask that
        holds back the long chains. Default is ``()``.

    Returns
    -------
    numpy.ndarray
        Boolean array of length ``num_chains``, True where the chain may be moved
        by the parallel kernel. All False when the box does not split into more
        than one block (there is no parallel work to do), or when there are no
        chains.
    """
    n_chains = len(chain_offset)
    parallel_ok = np.zeros(n_chains, dtype=bool)
    if n_chains == 0:
        return parallel_ok

    info = mega_crank_fast.parallel_layout_info(
        dimensions[0], dimensions[1],
        dimensions[2] if len(dimensions) == 3 else 1, bool(has_LR))
    if info["num_blocks"] == 1:
        return parallel_ok

    n_dim = len(dimensions)
    W = info["W"]
    split_interiors = [info["block_size"][d] - 2 * W
                       for d in range(n_dim) if info["blocks"][d] > 1]
    if not split_interiors:
        return parallel_ok

    # the smallest interior over the SPLIT axes: an unsplit axis has no halo and
    # therefore no length limit, so it must not tighten the bound
    interior = min(split_interiors)
    lengths = np.asarray(chain_length, dtype=np.int64)
    parallel_ok = lengths <= interior

    over_cap = lengths > 512
    if cap_mode == 'hetero' and chain_homo is not None:
        # the slither kernels' homopolymer path is O(1) and needs no per-bead
        # stack buffer, so the cap only bites on heteropolymers there
        over_cap = over_cap & (np.asarray(chain_homo, dtype=np.int64) == 0)
    parallel_ok = parallel_ok & ~over_cap

    if frozen_chains:
        first_ids = np.asarray(idx_to_bead)[np.asarray(chain_offset, dtype=np.int64), 4]
        parallel_ok = parallel_ok & ~np.isin(first_ids, list(frozen_chains))

    return parallel_ok


def parallel_chain_fit_report(idx_to_bead, chain_offset, chain_length, dimensions, has_LR,
                              chain_homo=None, cap_mode='all', frozen_chains=(), min_length=1):
    """Explain whether the parallel whole-chain kernels could move every eligible chain.

    This is a DIAGNOSTIC, not the dispatch gate. It is what the startup summary
    reads to describe the system, and the ``"partition"`` / ``"n_parallel"`` /
    ``"n_serial"`` entries are the actual dispatch decision, taken by
    :func:`parallel_chain_partition` from chain LENGTHS alone. The ``"ok"``,
    ``"n_too_extended"`` and ``"max_extent*"`` entries describe the CURRENT
    configuration and must never be used to pick a kernel: doing so was the
    pre-1.0.8 behaviour and it broke stationarity (see the long comment in
    :meth:`MoveObject.system_slither`).

    ``n_over_cap`` counts chains longer than the kernels' 512-bead per-thread buffer
    (hetero chains only for slither, ``cap_mode='hetero'``; every chain for pull,
    ``cap_mode='all'``); ``n_too_extended`` counts chains whose periodic extent on a
    split axis exceeds a block interior (``block_size - 2W``) and so could never fit
    for any random block shift. Either count > 0 means ``ok`` is False.

    Parameters
    ----------
    idx_to_bead : numpy.ndarray
        The ``int64`` bead bookkeeping matrix, shape ``(num_beads, 7)`` in 2D or
        ``(num_beads, 8)`` in 3D. Column 4 holds the chainID and columns 5
        onwards the bead coordinates, which is what the extent check reads.

    chain_offset : numpy.ndarray
        ``int32`` array of length ``num_chains`` holding the row of
        ``idx_to_bead`` at which each chain starts.

    chain_length : numpy.ndarray
        ``int32`` array of length ``num_chains`` holding the number of beads in
        each chain.

    dimensions : list of int
        The lattice box dimensions (2 or 3 entries), which set the block
        decomposition and the periodic extent calculation.

    has_LR : bool
        Whether any bead engages in long-range interactions. This widens the
        frozen halo (W) used by the checkerboard decomposition.

    chain_homo : numpy.ndarray or None, optional
        ``int32`` array of length ``num_chains``, 1 where a chain is
        homopolymeric. Only consulted when ``cap_mode`` is ``'hetero'``, where
        the 512-bead cap applies to heteropolymers only. Default is None.

    cap_mode : str, optional
        Which chains the kernels' 512-bead per-thread buffer applies to:
        ``'all'`` (the pull kernels) or ``'hetero'`` (the slither kernels, whose
        homopolymer path needs no per-bead buffer). Default is ``'all'``.

    frozen_chains : sequence of int, optional
        chainIDs that are permanently frozen. These are counted but never
        constrain the decision, because the kernels are not meant to move them.
        Default is ``()``.

    min_length : int, optional
        The shortest chain the move can act on (3 for pull, whose kernels skip
        shorter chains; 1 for slither). Chains below it are excluded from
        ``n_parallel`` / ``n_serial`` and counted in ``n_too_short``. Default
        is 1.

    Returns
    -------
    dict
        ``{"ok": bool, "layout": <parallel_layout_info dict>, "n_chains": int,
        "n_frozen": int, "n_over_cap": int, "n_too_extended": int,
        "max_extent": int, "interior": int, "single_block": bool,
        "interiors": tuple of int (movable width per axis),
        "max_extent_by_axis": tuple of int (largest chain extent per axis),
        "partition": bool ndarray, "n_parallel": int, "n_serial": int,
        "n_too_short": int}``.
        ``ok`` is True only when the decomposition has more than one block and
        both ``n_over_cap`` and ``n_too_extended`` are zero. ``partition`` is
        :func:`parallel_chain_partition` evaluated on the same arguments (True
        where a chain goes to the parallel kernel), and ``n_parallel`` /
        ``n_serial`` count its True / False entries over the non-frozen chains
        long enough for the move (``min_length``) - these are the dispatch, ``ok``
        is only the description. ``n_too_short`` counts the non-frozen chains
        shorter than ``min_length``, which neither kernel ever moves.
    """
    info = mega_crank_fast.parallel_layout_info(
        dimensions[0], dimensions[1],
        dimensions[2] if len(dimensions) == 3 else 1, bool(has_LR))
    W = info["W"]
    block_size = info["block_size"]
    nblocks = info["blocks"]
    idx = np.asarray(idx_to_bead)
    n_dim = len(dimensions)
    frozen_set = set(frozen_chains)
    n_frozen = n_over_cap = n_too_extended = 0
    max_extent = 0
    interiors = tuple(block_size[d] - 2 * W if nblocks[d] > 1 else dimensions[d]
                      for d in range(n_dim))
    max_extent_by_axis = [0] * n_dim
    interior = min(block_size[d] - 2 * W for d in range(n_dim) if nblocks[d] > 1) \
        if any(nblocks[d] > 1 for d in range(n_dim)) else int(max(dimensions))
    n_chains = len(chain_offset)
    offsets = np.asarray(chain_offset, dtype=np.int64)
    lengths = np.asarray(chain_length, dtype=np.int64)
    if n_chains == 0:
        frozen = np.zeros(0, dtype=bool)
        eligible = frozen
        extents = np.zeros((0, n_dim), dtype=np.int64)
    else:
        # chain index of every bead (the bookkeeping matrix is CSR-ordered)
        bead_chain = np.repeat(np.arange(n_chains), lengths)
        first_ids = idx[offsets, 4]
        frozen = np.isin(first_ids, list(frozen_set)) if frozen_set else np.zeros(n_chains, dtype=bool)
        eligible = ~frozen

        # periodic (circular) extent of every chain on every axis, vectorised: sort
        # the coordinates within each chain, take the largest gap between
        # neighbours (including the wrap-around gap from the last back to the
        # first), and the chain occupies the box minus that gap. A PBC-straddling
        # chain is thereby measured by its true physical span, not by its raw
        # wrapped coordinates (the kernel's random per-sweep block shift re-wraps
        # coordinates, so a physically short straddling chain CAN fit a block
        # interior). Coincident coordinates give a gap of 0, so duplicates need no
        # special handling, and a chain sitting on one coordinate has extent 1.
        beads = idx[:, 5:5 + n_dim]
        last = offsets + lengths - 1
        extents = np.empty((n_chains, n_dim), dtype=np.int64)
        for d in range(n_dim):
            box = int(dimensions[d])
            coords = beads[:, d]
            order = np.lexsort((coords, bead_chain))
            cs = coords[order]
            gap = np.zeros(len(cs), dtype=np.int64)
            gap[:-1] = cs[1:] - cs[:-1]
            gap[last] = 0                       # never read across a chain boundary
            within = np.maximum.reduceat(gap, offsets)
            wrap = cs[offsets] + box - cs[last]
            extents[:, d] = box - np.maximum(within, wrap) + 1

    n_frozen = int(frozen.sum())
    over_cap = eligible & (lengths > 512)
    if cap_mode == 'hetero':
        # the parallel kernels also hard-skip chains longer than their fixed
        # per-thread stack buffers (512 beads): HETERO chains in slither
        # (cap_mode='hetero', using chain_homo), and ALL chains in pull
        # (cap_mode='all'). Such a chain would silently never slither/pull under
        # the parallel kernel, so the partition sends it to the serial side.
        if chain_homo is not None:
            over_cap &= (np.asarray(chain_homo, dtype=np.int64) == 0)
    n_over_cap = int(over_cap.sum())

    if n_chains and eligible.any():
        elig_ext = extents[eligible]
        max_extent = int(elig_ext.max())
        max_extent_by_axis = [int(v) for v in elig_ext.max(axis=0)]
        # an axis that is not split has no interior limit (the extent is still
        # measured so the report never shows an impossible 0 for that axis)
        split = np.array([nblocks[d] > 1 for d in range(n_dim)])
        axis_interior = np.array([block_size[d] - 2 * W for d in range(n_dim)])
        too_extended = ((elig_ext > axis_interior) & split).any(axis=1)
        n_too_extended = int(too_extended.sum())
    else:
        max_extent = 0
        max_extent_by_axis = [0] * n_dim
        n_too_extended = 0
    # A one-block checkerboard has no parallel work to distribute.  Selecting the
    # parallel wrapper there only adds setup/OpenMP overhead and contradicted the
    # startup report, which correctly promised the serial kernel for such boxes.
    partition = parallel_chain_partition(idx_to_bead, chain_offset, chain_length, dimensions,
                                         has_LR, chain_homo=chain_homo, cap_mode=cap_mode,
                                         frozen_chains=frozen_chains)
    # the partition itself does not know the move's minimum length (pull never
    # touches a chain of fewer than three beads), so the counts apply it here
    long_enough = lengths >= int(min_length)
    return {"ok": info["num_blocks"] > 1 and n_over_cap == 0 and n_too_extended == 0,
            "layout": info, "single_block": info["num_blocks"] == 1,
            "n_chains": len(chain_offset), "n_frozen": n_frozen, "n_over_cap": n_over_cap,
            "n_too_extended": n_too_extended, "max_extent": max_extent,
            "max_extent_by_axis": tuple(max_extent_by_axis), "interior": interior,
            "interiors": interiors, "partition": partition,
            "n_parallel": int((partition & long_enough).sum()),
            "n_serial": int((~partition & eligible & long_enough).sum()),
            "n_too_short": int((eligible & ~long_enough).sum())}


def _parallel_can_move_all_chains(idx_to_bead, chain_offset, chain_length, dimensions, has_LR,
                                  chain_homo=None, cap_mode='all', frozen_chains=()):
    """Would the parallel checkerboard kernel be able to move every eligible chain?

    DIAGNOSTIC ONLY. Before 1.0.8 this was the dispatch gate for ``system_slither``
    and ``system_pull``, and that was wrong: it reads the CURRENT bead coordinates,
    so it made the choice of kernel a function of the configuration, which destroys
    stationarity even though each kernel on its own is pi-invariant (see the long
    comment in :meth:`MoveObject.system_slither`). The dispatch is now
    :func:`parallel_chain_partition`, which reads chain lengths only. Keep this
    function for reporting and for tests; do not route a kernel choice through it.

    The parallel slither/pull kernels only move a chain when all of its beads fit
    inside one block's interior (the block minus a width-W frozen halo on each face).
    A chain whose spatial extent on any split axis exceeds that interior does not fit
    for the current configuration; a chain LONGER than the interior can never fit for
    ANY random block shift, and it is that stricter, state-independent statement that
    the partition acts on.

    Permanently frozen chains do not constrain this answer: the kernel is not
    supposed to move them, and they are independently excluded by its frozen-bead
    mask. See :func:`parallel_chain_fit_report` for the reasons behind a False.

    Parameters
    ----------
    idx_to_bead : numpy.ndarray
        The ``int64`` bead bookkeeping matrix, shape ``(num_beads, 7)`` in 2D or
        ``(num_beads, 8)`` in 3D.

    chain_offset : numpy.ndarray
        ``int32`` array of length ``num_chains`` holding the row of
        ``idx_to_bead`` at which each chain starts.

    chain_length : numpy.ndarray
        ``int32`` array of length ``num_chains`` holding the number of beads in
        each chain.

    dimensions : list of int
        The lattice box dimensions (2 or 3 entries).

    has_LR : bool
        Whether any bead engages in long-range interactions (this sets the width
        of the frozen halo between blocks).

    chain_homo : numpy.ndarray or None, optional
        ``int32`` array of length ``num_chains``, 1 where a chain is
        homopolymeric. Only used when ``cap_mode`` is ``'hetero'``. Default is
        None.

    cap_mode : str, optional
        ``'all'`` (pull kernels) or ``'hetero'`` (slither kernels) - which chains
        the kernels' 512-bead per-thread buffer applies to. Default is ``'all'``.

    frozen_chains : sequence of int, optional
        chainIDs that are permanently frozen and therefore excluded from the
        check. Default is ``()``.

    Returns
    -------
    bool
        True only if every non-frozen chain currently fits a block interior and
        is under the buffer cap; False otherwise.
    """
    return parallel_chain_fit_report(idx_to_bead, chain_offset, chain_length, dimensions,
                                     has_LR, chain_homo=chain_homo, cap_mode=cap_mode,
                                     frozen_chains=frozen_chains)["ok"]


def _frozen_bead_mask(idx_to_bead, frozen_chains):
    """
    Build the per-bead frozen mask passed to the parallel Cython kernels.

    The parallel checkerboard kernels select movable beads/chains purely by their
    spatial block, so they need an explicit list of which beads must never move.
    A frozen bead is excluded from the movable set but remains in the grid as a
    fixed, energy-contributing obstacle - reproducing the serial behaviour, where
    frozen chains are simply never offered to the bead/chain selector.

    Parameters
    ----------
    idx_to_bead : numpy.ndarray
        The ``int64`` bead bookkeeping matrix, one row per bead, shape
        ``(num_beads, 7)`` in 2D or ``(num_beads, 8)`` in 3D. Column 4 holds the
        chainID of each bead, which is what the mask is built from.

    frozen_chains : sequence of int
        The chainIDs that are frozen. An empty sequence gives an all-zero mask.

    Returns
    -------
    numpy.ndarray
        A C-contiguous ``int32`` array of length ``num_beads`` whose entries are
        1 where the bead belongs to a frozen chain and 0 otherwise (all zeros
        when ``frozen_chains`` is empty, which leaves the kernels' behaviour
        bit-identical to the no-freeze case).
    """
    num_beads = idx_to_bead.shape[0]
    if not frozen_chains:
        return np.zeros(num_beads, dtype=np.int32)
    mask = np.isin(np.asarray(idx_to_bead)[:, 4], list(frozen_chains))
    return np.ascontiguousarray(mask.astype(np.int32))


def _two_pass_whole_chain_megamove(parallel_chains, serial_chains, substeps,
                                   parallel_kernel, serial_kernel, head_args, table_args,
                                   current_energy, invtemp, hardwall, max_len, num_threads,
                                   idx_to_bead, sorted_chains):
    """Run one slither or pull megamove over a length-partitioned chain set.

    The chains ``parallel_chains`` are moved by the parallel kernel and the
    chains ``serial_chains`` by the serial kernel, each pass holding the other's
    chains fixed. Each pass is reversible on its own: the serial pass visits a
    uniformly shuffled multiset of chains, and within a parallel sweep the
    movable sets are invariant and every block picks its chains i.i.d. But the
    two passes do not commute - the chains interact - so running them in a
    fixed order gives a kernel that is pi-invariant, which is all the main loop
    needs, but not reversible, and a system-wide TSMMC excursion (tempered
    transitions, which runs the same kernel up and down its schedule) needs
    each rung's kernel to be reversible. When both passes have chains the
    order is therefore drawn with a fair coin, which makes the megamove the
    self-adjoint mixture ``(PS + SP) / 2``. With only one non-empty pass there
    is nothing to order and no coin is drawn, so those runs are unchanged.

    The attempts logged are the ones the kernels report making: the parallel
    kernel makes none on a sweep in which its random block shift leaves no
    chain inside a block interior.

    Parameters
    ----------
    parallel_chains : numpy.ndarray
        Chain indices (rows of the chain layout) handed to the parallel kernel.

    serial_chains : numpy.ndarray
        Chain indices handed to the serial kernel.

    substeps : int
        Sub-moves per selectable chain (``SLITHER_SUBSTEPS`` or
        ``PULL_SUBSTEPS``).

    parallel_kernel, serial_kernel : callable
        The compiled megamove kernels for this move and dimensionality.

    head_args : tuple
        ``(grid, type_grid, idx_to_bead, chain_offset, chain_length,
        chain_homo)``, shared by both kernels, which mutate them in place.

    table_args : tuple
        The short-range, long-range and super-long-range interaction tables and
        the angle lookup.

    current_energy : int
        Total energy entering the megamove.

    invtemp : float
        Inverse temperature for the Metropolis test.

    hardwall : bool
        Whether the box has hard walls.

    max_len : int
        Length of the longest chain (the serial kernels' buffer size).

    num_threads : int
        Threads for the parallel kernel.

    idx_to_bead : numpy.ndarray
        The bead table, used to build the parallel pass's frozen mask.

    sorted_chains : list of int
        chainIDs in chain-layout order.

    Returns
    -------
    tuple
        ``(energy, proposed, accepted)`` after both passes.
    """
    energy = current_energy
    total_proposed = 0
    total_accepted = 0

    def parallel_pass(energy):
        """Run the parallel kernel over ``parallel_chains``.

        Parameters
        ----------
        energy : int
            Total energy entering the pass.

        Returns
        -------
        tuple
            ``(energy, accepted, attempted)`` as the parallel kernel reports them.
        """
        selector = np.repeat(parallel_chains, substeps)
        np.random.shuffle(selector)
        local_seed = random.randint(1, sys.maxsize - 1)
        # the parallel kernel reads chain_selector only for its LENGTH (the
        # sub-move budget) and picks chains per block, so the chains held out of
        # this pass - the long ones, the frozen ones and, for pull, any chain too
        # short to pull - are excluded through the frozen mask. They stay in the
        # grid as fixed, energy-contributing obstacles.
        parallel_set = set(int(ci) for ci in parallel_chains)
        held_out = [chainID for ci, chainID in enumerate(sorted_chains)
                    if ci not in parallel_set]
        frozen_mask = _frozen_bead_mask(idx_to_bead, held_out)
        return parallel_kernel(
            *(head_args + (selector,) + table_args
              + (energy, invtemp, local_seed, 1 if hardwall else 0, max_len)),
            num_threads, frozen_mask)

    def serial_pass(energy):
        """Run the serial kernel over ``serial_chains``.

        Parameters
        ----------
        energy : int
            Total energy entering the pass.

        Returns
        -------
        tuple
            ``(energy, accepted, attempted)``; the serial kernel makes every
            attempt in its selector, so ``attempted`` is the selector length.
        """
        # each selectable chain appears `substeps` times, in random order
        selector = np.repeat(serial_chains, substeps)
        np.random.shuffle(selector)
        local_seed = random.randint(1, sys.maxsize - 1)
        (energy, accepted) = serial_kernel(
            *(head_args + (selector,) + table_args
              + (energy, invtemp, local_seed, 1 if hardwall else 0, max_len)))
        return (energy, accepted, len(selector))

    passes = []
    if len(parallel_chains) > 0:
        passes.append(parallel_pass)
    if len(serial_chains) > 0:
        passes.append(serial_pass)
    if len(passes) == 2 and random.random() < 0.5:
        passes.reverse()

    for run_pass in passes:
        (energy, accepted, attempted) = run_pass(energy)
        total_proposed += attempted
        total_accepted += accepted

    return (energy, total_proposed, total_accepted)

## A note on single chain MC moves (cluster moves are fundementally different...)
## MoveCodes 2 3 4 5 and 6 
##
## All single chain MC moves defined here must fulfill the following properties
##
## 1) Input should be the current 2D or 3D lattice grid and the Chain object which is to
##    be moved. Note that we're passing a Grid [e.g. an np.array]! NOT a Lattice object.
## 
## 2) The Chain object should be treated as READ-ONLY. It is not actually READ ONLY 
##    - it *can* be augmented but it should not be. This object is also stored in 
##    the calling Simulation object's Lattice object, where it is assumed to remain 
##    the same through a move.
## 
## 3) The lattice grid should be altered IF a move can be made, and should be 
##    returned as it came if the move cannot. 
##
## 4) Moves defined here DO NOT evaluate energy but DO evaluate hardsphere clashes
##    - i.e. a move should be rejected if after the move happening we find a site is 
##    occupied by two residues
##
## 5) Moves return two variables: 1) a MoveEvent object which holds all the information needed for 
##    further move acceptance/rejection in a well defined structures 2) True or False to
##    declare if the move should be evaluated or not
##
##
## There are three other types of moves which should be written up
##
## 1) Cluster moves, which are fundementally different but do need to be evaluated by
##    functions in Simulation.py (i.e. energy evaluation is done here)
##
## 2) Moves which take avantage of the mega_crank subchain moves - these moves keep 
##    track of their state and energy and so do not need to be subsequently re-
##    evaluated. These moves are also insanely fast - example being system_shake
##
## 3) TSMMC moves - there are three classes of TSMMC moves, one where a single
##    chain is perturbed, one where a random set of chains are perturbed, and 
##    one where the ENTIRE system is perturbed. For the single chain and set 
##    of chains all moves are crankshaft-like so energy evaluation is done in
##    hyperloop. For the system wide one it basically shifts the entire simulation
##    engine into an auxillary chain so subsequent moves are done as normal but 
##    but at a different temperature, and at the end the complete series of changes
##    are accepted or rejected. 

 
## General 'GOTCHAS'
## The compiled kernels reseed their own generator at the start of every call, so
## any Cython function that uses random numbers must be PASSED a seed, drawn here
## from Python's generator. Because that generator is seeded from SEED (and saved in
## restart files for RESTART_CONTINUE), the passed seeds make every run
## reproducible.
##
##
##
##
#


#-----------------------------------------------------------------
#    
class MoveObject:
    """The Monte Carlo move set: one method per move type, each proposing (and
    on acceptance applying) a change to the lattice.

    Moves either return a ``MoveEvent`` for the shared
    acceptance machinery (the local and rigid-body moves) or perform their own
    accept/reject internally on the same Markov chain (the megamoves: system
    crankshaft/slither/pull dispatched to the Cython kernels, VMMC, jump-and-
    relax and the TSMMC family). Every move is constructed to satisfy detailed
    balance; see the per-move documentation pages for the specific arguments.
    """

    def __init__(self):
        """
        MoveObject is a stateless class, who's objects implement chain movement
        functionality but do not have any state associated with themselves

        """
        pass



    #-----------------------------------------------------------------
    #    
    #
    def system_shake(self, latticeObject, current_energy, acceptanceObject, hamiltonianObject, number_of_steps, mode, hardwall=False, frozen_chains=(), parallelize=False, num_threads=1):
        """
        Perform a whole-system crankshaft megamove (MoveType code 1).

        The system_shake move performs a large number of very local single-bead
        perturbations. Each perturbation randomly selects any bead on the
        lattice, ensuring that complete detailed balance is maintained. The
        individual accept/reject decisions happen per-sub-move inside the
        optimized Cython kernel (the same Markov chain), so the move does not
        need to be re-evaluated afterwards. The chain positions are written back
        from the idx_to_bead matrix once the kernel returns.

        The appropriate kernel is selected automatically (see the in-body
        comments): a run with ``parallelize`` set uses the parallel checkerboard
        kernel (2D or 3D, with frozen chains honoured via a per-bead frozen mask);
        otherwise the serial fast kernel (2D or 3D) is used.

        Parameters
        ----------
        latticeObject : Lattice
            The full lattice object upon which the simulation is being
            performed. Its grids are mutated in place by the kernel.

        current_energy : int or float
            The current system energy value (before this megamove).

        acceptanceObject : AcceptanceCalculator
            Object containing the inverse temperature and all details needed to
            accept or reject a move.

        hamiltonianObject : Hamiltonian
            Self-contained object that allows the evaluation of energy functions
            and contains the interaction tables passed to the external (Cython)
            kernel for energy evaluation.

        number_of_steps : int
            Number of Monte Carlo sub-moves to perform across the system.

        mode : str
            Mode used for determining the final number of steps. Currently
            obsolete but kept in case bead selection is changed in the future.

        hardwall : bool, optional
            If True a hard-wall (impenetrable solvent) boundary is used;
            otherwise periodic boundary conditions are used. Default is False.

        frozen_chains : sequence of int, optional
            chainIDs that are frozen and cannot be moved. Default is ``()``.
            Frozen chains are honoured by both the serial and the parallel
            kernels (their beads are kept as fixed obstacles but never selected
            for a move).

        parallelize : bool, optional
            If True, use the multi-threaded checkerboard kernel (2D or 3D); frozen
            chains are fully supported. Default is False.

        num_threads : int, optional
            Number of threads to use when the parallel kernel is selected.
            Default is 1.

        Returns
        -------
        tuple
            ``(latticeObject, current_energy, total_proposed, total_accepted)``
            where ``latticeObject`` is the updated lattice, ``current_energy``
            the new system energy, ``total_proposed`` the number of sub-moves
            attempted, and ``total_accepted`` the number accepted.
        """
        
        # construct the idx_to_bead matrix. which gets passed into megacrank. This matrix contains position and identity information
        # for all beads on the lattice
        frozen_chains = _normalized_frozen_chains(latticeObject, frozen_chains)
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead(latticeObject)

        # get the number of beads, build a selection vector, and figure out the number of dimensions
        num_beads = len(idx_to_bead)        
        num_dims  = len(latticeObject.dimensions)
        
        total_accepted = 0
        total_proposed = 0 

        if len(frozen_chains) >= len(latticeObject.chains):
            return (latticeObject, current_energy, 0, 0)

        # set hardwall flag
        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0

        local_seed = random.randint(1, sys.maxsize - 1)


        #bead_selector = np.random.randint(0,num_beads,number_of_steps)
        bead_selector = crankshaft_list_functions.bead_selector_constructor(num_beads, number_of_steps, latticeObject, frozen_chains=frozen_chains, safecheck=True)

        # per-bead frozen mask for the parallel kernel (1 = bead belongs to a
        # frozen chain). idx_to_bead column 4 is the chainID. Frozen beads are
        # excluded from the movable set but stay as fixed energy-contributing
        # obstacles - so the parallel kernel honours freezing exactly like serial.
        frozen_mask = _frozen_bead_mask(idx_to_bead, frozen_chains)

        ##
        ## Both functions alter alter the grids on the back end and do not explicity
        ## reassign these as they're passed by reference as memoryviews (direct access to
        ## the memory)
        ## 

        
        # ------------------------------------------------------------------
        # Kernel selection. IMPORTANT SAFETY RULES (so enabling parallelize can
        # never silently produce wrong physics):
        #
        #   * The parallel checkerboard kernel exists for both 2D
        #     (mega_crank_parallel_2D) and 3D (mega_crank_parallel). Both use the
        #     same frozen-halo block decomposition, which is independent of the
        #     thread count, and target the same Boltzmann distribution as the
        #     serial kernels (they are NOT bit-identical to them).
        #   * The parallel kernel buckets ALL beads spatially; frozen chains are
        #     honoured by passing a per-bead frozen_mask (frozen beads are kept as
        #     fixed obstacles but never selected for a move), so it can be used
        #     even with frozen chains.
        #   * The serial fast kernels (mega_crank_fast.mega_crank /
        #     .mega_crank_2D) are bit-exact drop-ins for the reference kernels,
        #     so they are safe for every keyword combination the reference
        #     handled (NON_INTERACTING, ANGLES_OFF, HARDWALL, QUENCH, frozen
        #     chains, etc.).
        # ------------------------------------------------------------------

        # the serial kernels make exactly number_of_steps attempts; the parallel
        # ones report their own count, which is 0 on a sweep whose random block
        # shift leaves no bead inside a block interior
        attempted_moves = number_of_steps
        use_parallel = False
        if parallelize:
            crank_layout = mega_crank_fast.parallel_crank_layout_info(
                latticeObject.dimensions[0], latticeObject.dimensions[1],
                latticeObject.dimensions[2] if num_dims == 3 else 1,
                bool(np.any(np.asarray(idx_to_bead)[:, 1] == 1)))
            use_parallel = crank_layout["num_blocks"] > 1

        # 2D
        if num_dims == 2:

            # 2D + parallel requested -> 2D checkerboard kernel (frozen chains
            # honoured via frozen_mask)
            if use_parallel:
                (new_energy, accepted_moves, attempted_moves) = mega_crank_fast.mega_crank_parallel_2D(latticeObject.grid,
                                                                                      latticeObject.type_grid,
                                                                                      idx_to_bead,
                                                                                      hamiltonianObject.residue_interaction_table,
                                                                                      hamiltonianObject.LR_residue_interaction_table,
                                                                                      hamiltonianObject.SLR_residue_interaction_table,
                                                                                      hamiltonianObject.angle_lookup,
                                                                                      current_energy,
                                                                                      acceptanceObject.invtemp,
                                                                                      number_of_steps,
                                                                                      local_seed,
                                                                                      hardwall_int,
                                                                                      num_threads,
                                                                                      frozen_mask)

            # 2D serial (default)
            else:
                (new_energy, accepted_moves)= mega_crank_fast.mega_crank_2D(latticeObject.grid,
                                                                          latticeObject.type_grid,
                                                                          idx_to_bead,
                                                                          hamiltonianObject.residue_interaction_table,
                                                                          hamiltonianObject.LR_residue_interaction_table,
                                                                          hamiltonianObject.SLR_residue_interaction_table,
                                                                          hamiltonianObject.angle_lookup,
                                                                          current_energy,
                                                                          acceptanceObject.invtemp,
                                                                          number_of_steps,
                                                                          bead_selector,
                                                                          local_seed,
                                                                          hardwall_int)

        # 3D + parallel requested -> checkerboard kernel (frozen chains honoured
        # via frozen_mask)
        elif use_parallel:
            (new_energy, accepted_moves, attempted_moves) = mega_crank_fast.mega_crank_parallel(latticeObject.grid,
                                                                               latticeObject.type_grid,
                                                                               idx_to_bead,
                                                                               hamiltonianObject.residue_interaction_table,
                                                                               hamiltonianObject.LR_residue_interaction_table,
                                                                               hamiltonianObject.SLR_residue_interaction_table,
                                                                               hamiltonianObject.angle_lookup,
                                                                               current_energy,
                                                                               acceptanceObject.invtemp,
                                                                               number_of_steps,
                                                                               local_seed,
                                                                               hardwall_int,
                                                                               num_threads,
                                                                               frozen_mask)

        # 3D serial (default)
        else:
            (new_energy, accepted_moves) = mega_crank_fast.mega_crank(latticeObject.grid,
                                                                 latticeObject.type_grid,
                                                                 idx_to_bead,
                                                                 hamiltonianObject.residue_interaction_table,
                                                                 hamiltonianObject.LR_residue_interaction_table,
                                                                 hamiltonianObject.SLR_residue_interaction_table,
                                                                 hamiltonianObject.angle_lookup,
                                                                 current_energy,
                                                                 acceptanceObject.invtemp,
                                                                 number_of_steps,
                                                                 bead_selector,
                                                                 local_seed,
                                                                 hardwall_int)

            





        total_accepted = total_accepted + accepted_moves
        total_proposed = total_proposed + attempted_moves
        
        # finally push the positions the kernel left in idx_to_bead back into the
        # Chain objects (the grids were updated in the kernel itself). This is the
        # compiled write-back; done per chain in Python it was a large part of
        # the cost of a megamove at scale.
        crankshaft_list_functions.write_back_positions(latticeObject, idx_to_bead)

        current_energy = new_energy

        return (latticeObject, current_energy, total_proposed, total_accepted)


    #-----------------------------------------------------------------
    #
    def system_slither(self, latticeObject, current_energy, acceptanceObject, hamiltonianObject, slither_substeps, hardwall=False, frozen_chains=(), parallelize=False, num_threads=1):
        """
        Whole-system slither (reptation) megamove (2D and 3D). Every non-frozen
        chain is slithered ``slither_substeps`` times, in random order, by the
        optimized Cython kernel (mega_crank_fast.mega_slither /
        mega_crank_fast.mega_slither_2D). A slither advances a chain forwards or
        backwards like a snake:

          * homopolymers     -> O(1) interaction energy (only the moved end matters)
          * heteropolymers   -> every residue re-evaluated (decomposed into single
                                bead moves, reusing the validated energy primitives)
          * single-bead chains -> a local translation

        Like system_shake the accept/reject happens per-sub-move inside the kernel
        (same Markov chain), and the chain positions are written back from the
        idx_to_bead matrix afterwards.

        MoveType code: 6

        Parameters
        ----------
        latticeObject : Lattice
            The full lattice object being simulated; its grids and chain
            positions are mutated in place.

        current_energy : int or float
            The current system energy value (before this megamove).

        acceptanceObject : AcceptanceCalculator
            Object providing the inverse temperature used by the kernel.

        hamiltonianObject : Hamiltonian
            Object providing the interaction tables and angle lookup passed to
            the Cython kernel for energy evaluation.

        slither_substeps : int
            Number of times each non-frozen chain is slithered (every selectable
            chain appears this many times in the randomized selection order).

        hardwall : bool, optional
            If True a hard-wall boundary is used; otherwise periodic boundary
            conditions are used. Default is False.

        frozen_chains : sequence of int, optional
            chainIDs that are frozen and excluded from slithering (their beads
            remain as fixed obstacles). Default is ``()``.

        parallelize : bool, optional
            If True, the chains are split by LENGTH (see
            :func:`parallel_chain_partition`) and the megamove runs two passes:
            the multi-threaded checkerboard kernel over the chains short enough
            to be guaranteed to fit a block interior, and the serial kernel over
            the rest. The split depends only on chain lengths, the frozen set and
            the box, never on the current configuration, so both passes and their
            composition leave the equilibrium distribution alone. Default is
            False, which runs everything on the serial kernel.

        num_threads : int, optional
            Number of threads handed to the parallel pass. Ignored when no chain
            is short enough for it. Default is 1.

        Returns
        -------
        tuple
            ``(latticeObject, current_energy, total_proposed, total_accepted)``,
            with the two passes' proposal and acceptance counts summed. If every
            chain is frozen, ``(latticeObject, current_energy, 0, 0)`` is returned
            unchanged.
        """

        frozen_chains = _normalized_frozen_chains(latticeObject, frozen_chains)
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead(latticeObject)

        # per-chain metadata in the same (sorted-chainID) order as idx_to_bead:
        # row offsets, lengths and the homopolymer flags (uniform intcode AND
        # long-range flag, the precondition of the kernels' O(1) energy path).
        # None of it changes during a run, so it comes from the layout the
        # lattice built once rather than being rebuilt chain by chain here on
        # every megamove. The kernels get their own copies of the int32 arrays.
        layout = crankshaft_list_functions.chain_layout(latticeObject)
        sorted_chains = layout.sorted_ids
        chain_offset = layout.offset32.copy()
        chain_length = layout.length32.copy()
        chain_homo   = layout.homo.copy()

        frozen_set = set(frozen_chains)
        selectable = [ci for ci, chainID in enumerate(sorted_chains)
                      if chainID not in frozen_set]

        # nothing to do if every chain is frozen
        if len(selectable) == 0:
            return (latticeObject, current_energy, 0, 0)

        # ------------------------------------------------------------------
        # Kernel dispatch. The choice of kernel MUST NOT depend on the current
        # configuration. Each kernel on its own is pi-invariant, but a
        # state-dependent choice BETWEEN them is not, and here the asymmetry is
        # total: the parallel kernel's movable interior is closed (it rejects any
        # sub-move whose new site leaves the interior), so the transition rate out
        # of "every chain is compact enough to fit" is identically zero through it,
        # while the serial kernel crosses that boundary freely. Before 1.0.8
        # this dispatch asked whether every chain's CURRENT extent fits a
        # block interior and sent the whole megamove to one kernel or the other on
        # the answer; that pumped probability into compact conformations and biased
        # <Rg^2> low by a few percent, one-signed, with nothing to warn the user.
        # Do NOT reintroduce a per-configuration gate here.
        #
        # Instead the chains are partitioned ONCE, by LENGTH, which is a run
        # constant (a chain of n beads spans at most n sites, so length <= interior
        # is a state-independent sufficient condition for fitting). Short chains are
        # slithered by the parallel kernel and long chains by the serial kernel, in
        # the same megamove: both passes are reversible and the partition never
        # moves, so the composition is pi-invariant too, and running the two passes
        # in a random order (_two_pass_whole_chain_megamove) keeps it reversible,
        # which a system-wide TSMMC excursion needs. Long chains keep reptating
        # (no ergodicity hole) and short chains keep the parallel speed-up.
        # ------------------------------------------------------------------
        selectable = np.array(selectable, dtype=np.int32)
        if parallelize:
            parallel_ok = parallel_chain_partition(
                idx_to_bead, chain_offset, chain_length, latticeObject.dimensions,
                np.any(np.asarray(idx_to_bead)[:, 1] == 1),
                chain_homo=chain_homo, cap_mode='hetero',
                frozen_chains=frozen_chains)
            parallel_chains = selectable[parallel_ok[selectable]]
            serial_chains = selectable[~parallel_ok[selectable]]
        else:
            parallel_chains = selectable[:0]
            serial_chains = selectable

        if len(latticeObject.dimensions) == 2:
            parallel_kernel = mega_crank_fast.mega_slither_parallel_2D
            serial_kernel = mega_crank_fast.mega_slither_2D
        else:
            parallel_kernel = mega_crank_fast.mega_slither_parallel
            serial_kernel = mega_crank_fast.mega_slither

        # the kernels mutate grid / type_grid / idx_to_bead in place, so the two
        # passes chain through the same buffers and the running energy
        head_args = (latticeObject.grid,
                     latticeObject.type_grid,
                     idx_to_bead,
                     chain_offset,
                     chain_length,
                     chain_homo)
        table_args = (hamiltonianObject.residue_interaction_table,
                      hamiltonianObject.LR_residue_interaction_table,
                      hamiltonianObject.SLR_residue_interaction_table,
                      hamiltonianObject.angle_lookup)
        max_len = int(chain_length.max())

        # both passes, in a random order when both have chains (see the helper)
        (new_energy, total_proposed, total_accepted) = _two_pass_whole_chain_megamove(
            parallel_chains, serial_chains, slither_substeps, parallel_kernel, serial_kernel,
            head_args, table_args, current_energy, acceptanceObject.invtemp, hardwall,
            max_len, num_threads, idx_to_bead, sorted_chains)

        # write the moved chain positions back from idx_to_bead (compiled)
        crankshaft_list_functions.write_back_positions(latticeObject, idx_to_bead)

        return (latticeObject, new_energy, total_proposed, total_accepted)


    #-----------------------------------------------------------------
    #
    def system_pull(self, latticeObject, current_energy, acceptanceObject, hamiltonianObject, pull_substeps, hardwall=False, frozen_chains=(), parallelize=False, num_threads=1):
        """
        Whole-system pull (cooperative reptation) megamove (2D and 3D). Every
        non-frozen chain of length >= 3 is pulled ``pull_substeps`` times, in
        random order, by the optimized Cython kernel (mega_crank_fast.mega_pull /
        mega_crank_fast.mega_pull_2D).

        A pull move displaces an interior bead to a neighbouring empty site and
        cooperatively pulls the following beads into the vacated sites until chain
        connectivity is restored - local reptation of a sub-segment that lets
        chains rearrange in DENSE systems where rigid moves would clash. Detailed
        balance is maintained inside the kernel via a Metropolis-Hastings
        proposal-multiplicity correction.

        Like system_slither the accept/reject happens per-sub-move inside the
        kernel (same Markov chain) and chain positions are written back from the
        idx_to_bead matrix afterwards.

        MoveType code: 11

        Parameters
        ----------
        latticeObject : Lattice
            The full lattice object being simulated; its grids and chain
            positions are mutated in place.

        current_energy : int or float
            The current system energy value (before this megamove).

        acceptanceObject : AcceptanceCalculator
            Object providing the inverse temperature used by the kernel.

        hamiltonianObject : Hamiltonian
            Object providing the interaction tables and angle lookup passed to
            the Cython kernel for energy evaluation.

        pull_substeps : int
            Number of times each eligible chain (non-frozen, length >= 3) is
            pulled (every selectable chain appears this many times in the
            randomized selection order).

        hardwall : bool, optional
            If True a hard-wall boundary is used; otherwise periodic boundary
            conditions are used. Default is False.

        frozen_chains : sequence of int, optional
            chainIDs that are frozen and excluded from pulling (their beads
            remain as fixed obstacles). Default is ``()``.

        parallelize : bool, optional
            If True, the chains are split by LENGTH (see
            :func:`parallel_chain_partition`) and the megamove runs two passes:
            the multi-threaded checkerboard kernel over the chains short enough to
            be guaranteed to fit a block interior, and the serial kernel over the
            rest. The split depends only on chain lengths, the frozen set and the
            box, never on the current configuration, so both passes and their
            composition leave the equilibrium distribution alone. Default is
            False, which runs everything on the serial kernel.

        num_threads : int, optional
            Number of threads handed to the parallel pass. Ignored when no chain
            is short enough for it. Default is 1.

        Returns
        -------
        tuple
            ``(latticeObject, current_energy, total_proposed, total_accepted)``,
            with the two passes' proposal and acceptance counts summed. If no
            chain is long enough or all chains are frozen,
            ``(latticeObject, current_energy, 0, 0)`` is returned unchanged.
        """

        frozen_chains = _normalized_frozen_chains(latticeObject, frozen_chains)
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead(latticeObject)

        # per-chain metadata in the same (sorted-chainID) order as idx_to_bead:
        # row offsets, lengths and the homopolymer flags (uniform intcode AND
        # long-range flag, the precondition of the kernels' O(1) energy path).
        # None of it changes during a run, so it comes from the layout the
        # lattice built once rather than being rebuilt chain by chain here on
        # every megamove. The kernels get their own copies of the int32 arrays.
        layout = crankshaft_list_functions.chain_layout(latticeObject)
        sorted_chains = layout.sorted_ids
        chain_offset = layout.offset32.copy()
        chain_length = layout.length32.copy()
        chain_homo   = layout.homo.copy()

        frozen_set = set(frozen_chains)
        # a pull needs an interior bead with neighbours on both sides (L >= 3)
        selectable = [ci for ci, chainID in enumerate(sorted_chains)
                      if chainID not in frozen_set and chain_length[ci] >= 3]

        # nothing to do if no chain is long enough / all frozen
        if len(selectable) == 0:
            return (latticeObject, current_energy, 0, 0)

        # Kernel dispatch, exactly as in system_slither: the chains are split ONCE
        # by LENGTH into a parallel set and a serial set and BOTH kernels run, one
        # after the other, on their own chains. The partition depends only on chain
        # lengths, the frozen set, the box and has_LR, so it is a run constant and
        # the composition of the two pi-invariant passes (run in a random order, so
        # it is also reversible) is pi-invariant. Read the
        # long comment in system_slither before changing this: choosing one kernel
        # per megamove from the current configuration (the pre-1.0.8 behaviour) is not
        # pi-invariant and silently over-samples compact conformations.
        selectable = np.array(selectable, dtype=np.int32)
        if parallelize:
            parallel_ok = parallel_chain_partition(
                idx_to_bead, chain_offset, chain_length, latticeObject.dimensions,
                np.any(np.asarray(idx_to_bead)[:, 1] == 1),
                cap_mode='all', frozen_chains=frozen_chains)
            parallel_chains = selectable[parallel_ok[selectable]]
            serial_chains = selectable[~parallel_ok[selectable]]
        else:
            parallel_chains = selectable[:0]
            serial_chains = selectable

        if len(latticeObject.dimensions) == 2:
            parallel_kernel = mega_crank_fast.mega_pull_parallel_2D
            serial_kernel = mega_crank_fast.mega_pull_2D
        else:
            parallel_kernel = mega_crank_fast.mega_pull_parallel
            serial_kernel = mega_crank_fast.mega_pull

        head_args = (latticeObject.grid,
                     latticeObject.type_grid,
                     idx_to_bead,
                     chain_offset,
                     chain_length,
                     chain_homo)
        table_args = (hamiltonianObject.residue_interaction_table,
                      hamiltonianObject.LR_residue_interaction_table,
                      hamiltonianObject.SLR_residue_interaction_table,
                      hamiltonianObject.angle_lookup)
        max_len = int(chain_length.max())

        # both passes, in a random order when both have chains (see the helper)
        (new_energy, total_proposed, total_accepted) = _two_pass_whole_chain_megamove(
            parallel_chains, serial_chains, pull_substeps, parallel_kernel, serial_kernel,
            head_args, table_args, current_energy, acceptanceObject.invtemp, hardwall,
            max_len, num_threads, idx_to_bead, sorted_chains)

        # write the moved chain positions back from idx_to_bead (compiled)
        crankshaft_list_functions.write_back_positions(latticeObject, idx_to_bead)

        return (latticeObject, new_energy, total_proposed, total_accepted)


    #-----------------------------------------------------------------
    #     VIRTUAL-MOVE MONTE CARLO (collective move, code 14)
    #
    def _vmmc_draw_nc(self, n_chains, max_cluster):
        """
        Draw a VMMC cluster-size cutoff from ``Q(n_c)`` proportional to ``1/n_c``.

        The cutoff ``n_c`` is drawn over ``[1, cap]`` where
        ``cap = min(max_cluster, n_chains)``. This per-particle move-frequency
        correction is symmetric (independent of move direction) so it cancels
        between the forward and reverse VMMC proposals.

        Parameters
        ----------
        n_chains : int
            Number of chains currently in the system.

        max_cluster : int
            Configured upper bound on the cutoff (``VMMC_MAX_CLUSTER``); the cap is
            the smaller of this and ``n_chains``.

        Returns
        -------
        int
            The drawn cluster-size cutoff in ``[1, cap]`` (``1`` when ``cap <= 1``).
        """
        cap = int(min(max_cluster, n_chains))
        if cap <= 1:
            return 1

        index = bisect.bisect_left(_vmmc_cutoff_cdf(cap), random.random())
        return cap if index >= cap else index + 1


    def _vmmc_neighbour_energies(self, latticeObject, hamiltonianObject, hardwall, offsets, m_id, positions, intcodes, lr_flags, offset, dimensions):
        """
        Per-neighbour interaction energy of a (virtually shifted) chain.

        Computes the interaction energy between chain ``m_id`` - whose beads sit at
        ``positions`` shifted by ``offset`` - and every OTHER chain it touches,
        returned as a ``{chainID: energy}`` mapping.

        Neighbours are read from the UN-MUTATED grid/type_grid at their real
        positions (so this implements "move chain m alone"); m's own beads
        (``grid == m_id``) and solvent (0) are skipped. SR uses the Chebyshev-1
        shell, LR/SLR (only for LR beads) the Chebyshev-2/-3 shells - matching the
        energy model. Only the cross m-j interaction is needed and it merely shapes
        the recruitment proposal (detailed balance is enforced by the exact dE plus
        the consistently-computed forward/reverse proposal ratio), so it need not be
        bit-identical to ``evaluate_total_energy``.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice whose ``grid`` (chainIDs) and ``type_grid`` (intcodes) are
            scanned for neighbours.

        hamiltonianObject : Hamiltonian
            Provides the SR/LR/SLR residue interaction tables.

        hardwall : bool
            If True, neighbours across the box boundary do not interact; otherwise
            periodic boundary conditions are applied.

        offsets : dict of int to sequence of (tuple, int)
            Precomputed ``(neighbour offset, Chebyshev radius)`` pairs keyed by
            half-width (``1`` for the SR shell, ``3`` for LR/SLR shells).

        m_id : int
            chainID of the chain being virtually moved.

        positions : list of list of int
            Ordered bead positions of chain ``m_id``, unshifted (as returned by
            ``Chain.get_ordered_positions()``).

        intcodes : list of int
            Integer residue-type code for each bead of chain ``m_id``, indexing
            the interaction tables.

        lr_flags : numpy.ndarray
            Per-bead long-range flags for chain ``m_id`` (1 where the bead
            engages in LR/SLR interactions, 0 otherwise), as returned by
            ``Chain.get_LR_binary_array()``.

        offset : list of int
            Per-dimension translation applied to ``positions`` before scanning
            neighbours (``[0, 0, ...]`` gives the current configuration).

        dimensions : list of int
            Lattice box dimensions, used for boundary handling (hardwall vs PBC).

        Returns
        -------
        dict
            Mapping of neighbouring ``chainID`` to the summed cross interaction
            energy with chain ``m_id`` in the (shifted) configuration. Chains with
            zero net interaction are omitted.
        """
        grid = latticeObject.grid
        tg   = latticeObject.type_grid
        SRT  = hamiltonianObject.residue_interaction_table
        LRT  = hamiltonianObject.LR_residue_interaction_table
        SLRT = hamiltonianObject.SLR_residue_interaction_table
        nd   = len(dimensions)
        energies = {}

        for b in range(len(positions)):
            t_b   = int(intcodes[b])
            is_lr = bool(lr_flags[b])
            rng   = 3 if is_lr else 1
            base  = [positions[b][d] + offset[d] for d in range(nd)]

            for delta, cheb in offsets[rng]:
                npos = []
                straddle = False
                for d in range(nd):
                    coord = base[d] + delta[d]
                    if hardwall:
                        if coord < 0 or coord >= dimensions[d]:
                            straddle = True
                            break
                        npos.append(coord)
                    else:
                        npos.append(coord % dimensions[d])
                if straddle:
                    continue

                j = int(lattice_utils.get_gridvalue(npos, grid))
                if j == 0 or j == m_id:
                    continue
                t_n = int(lattice_utils.get_gridvalue(npos, tg))

                if cheb == 1:
                    e = SRT[t_b][t_n]
                elif cheb == 2:
                    e = LRT[t_b][t_n] if is_lr else 0.0
                else:
                    e = SLRT[t_b][t_n] if is_lr else 0.0

                if e != 0.0:
                    energies[j] = energies.get(j, 0.0) + e

        return energies


    def vmmc_move(self, seed_chain, latticeObject, current_energy, acceptanceObject, hamiltonianObject, max_displacement, max_cluster, hardwall=False, frozen_chains=()):
        """
        Virtual-Move Monte Carlo collective move (Whitelam & Geissler, J. Chem.
        Phys. 127, 154101, 2007). Translation-only. MoveType code: 14.

        A seed chain is given a trial rigid lattice translation; neighbouring chains
        are recruited into a moving cluster according to interaction-energy gradients
        (a neighbour is recruited when moving the seed alone would break their mutual
        attraction), and the whole cluster translates together. This lets correlated
        groups of chains move collectively, escaping the kinetic traps that single
        chain moves hit in strongly-attractive / condensed phases.

        Detailed balance is enforced as Metropolis-Hastings: the recruitment is the
        proposal, and acceptance multiplies the exact Boltzmann factor exp(-beta*dE)
        by the reverse/forward proposal ratio assembled from the link formation (p)
        and failure (q = 1 - p) probabilities of EVERY tested link - formed, boundary
        and frustrated (failed, partner recruited elsewhere) - see
        :func:`_vmmc_log_acceptance`. The move is self-contained - it
        applies, accepts or reverts the configuration in place and returns the
        resulting energy - mirroring the other whole-system moves (e.g.
        :meth:`system_pull`).

        Parameters
        ----------
        seed_chain : Chain
            The (uniformly selected) seed chain.

        latticeObject : Lattice
            The system lattice (mutated in place on acceptance).

        current_energy : float
            Current total system energy (maintained exactly by the master loop).

        acceptanceObject : AcceptanceCalculator
            Supplies the inverse temperature ``beta`` (``acceptanceObject.invtemp``).

        hamiltonianObject : Hamiltonian
            Used both for the cross-chain link energies and for the from-scratch
            total-energy recompute that gives the exact ``dE``.

        max_displacement : int
            Maximum magnitude (per dimension) of the trial translation
            (``VMMC_MAX_DISPLACEMENT``).

        max_cluster : int
            Upper bound on the cluster-size cutoff draw (``VMMC_MAX_CLUSTER``).

        hardwall : bool, optional
            Whether the lattice has hard-wall (non-periodic) boundaries. Default
            is False.

        frozen_chains : sequence of int, optional
            chainIDs that may not move; recruiting a frozen chain (or seeding on
            one) rejects the move. Default is ``()``.

        Returns
        -------
        (Lattice, float, bool, int)
            The (mutated) lattice, the new total energy (unchanged on rejection),
            whether the move was accepted, and the size of the recruited cluster.
        """
        dimensions = latticeObject.dimensions
        nd         = len(dimensions)
        beta       = acceptanceObject.invtemp
        seed_id    = int(seed_chain.chainID)
        frozen_set = set(_normalized_frozen_chains(latticeObject, frozen_chains))

        if seed_id in frozen_set:
            return (latticeObject, current_energy, False, 1)

        # --- trial translation dr (symmetric: dr and -dr are equiprobable) -------
        dr = []
        for d in range(nd):
            mag = random.randint(1, min(dimensions[d] - 1, max_displacement))
            dr.append(mag if random.random() < 0.5 else -mag)
        neg_dr = [-x for x in dr]

        # neighbour-offset tuples for the SR (Chebyshev-1) and LR/SLR (Chebyshev-3)
        # shells, computed once and reused by every link-energy scan.
        offsets = {1: _vmmc_offset_shell(nd, 1),
                   3: _vmmc_offset_shell(nd, 3)}

        # --- cluster-size cutoff, drawn BEFORE growth (1/n_c frequency factor) ----
        n_c = self._vmmc_draw_nc(latticeObject.get_number_of_chains(), max_cluster)

        meta = {}
        def get_meta(cid):
            """
            Return (and memoize) the (positions, intcodes, LR-flags) of a chain.

            Parameters
            ----------
            cid : int
                chainID whose cached metadata is requested.

            Returns
            -------
            tuple
                ``(ordered_positions, intcode_sequence, LR_binary_array)`` for the
                chain.
            """
            if cid not in meta:
                ch = latticeObject.chains[cid]
                meta[cid] = (ch.get_ordered_positions(),
                             ch.get_intcode_sequence(),
                             ch.get_LR_binary_array())
            return meta[cid]

        # --- recruit the cluster (BFS over chains) on the un-mutated lattice ------
        cluster      = {seed_id}
        queue        = [seed_id]
        formed_links = []   # (p_f, p_r) for links that formed (built the cluster)
        failed_links = []   # (j, p_f, p_r_boundary, p_r_internal) for tested links that did NOT form
        tested       = set()

        while queue:
            m = queue.pop()
            (P, ic, lr) = get_meta(m)
            E0 = self._vmmc_neighbour_energies(latticeObject, hamiltonianObject, hardwall, offsets, m, P, ic, lr, [0] * nd, dimensions)
            Ef = self._vmmc_neighbour_energies(latticeObject, hamiltonianObject, hardwall, offsets, m, P, ic, lr, dr,       dimensions)
            Er = self._vmmc_neighbour_energies(latticeObject, hamiltonianObject, hardwall, offsets, m, P, ic, lr, neg_dr,   dimensions)

            # sorted: the reverse move must test the same directed links in the same
            # order for its realisation to mirror this one.
            for j in sorted(set(E0) | set(Ef) | set(Er)):
                pair = frozenset((m, j))
                if pair in tested:
                    continue
                tested.add(pair)

                e0  = E0.get(j, 0.0)
                p_f = _vmmc_link_probability(beta, Ef.get(j, 0.0) - e0)

                # Two reverse link probabilities; which one applies is decided at
                # acceptance time by where the partner j ends up:
                #  * j INSIDE the final cluster - the link formed, or it failed but j
                #    was recruited through another chain (a "frustrated" link): in
                #    the reverse move m and j translate together, only their relative
                #    displacement matters, and the reverse probability comes from
                #    the energy with m shifted by -dr (Er);
                #  * j OUTSIDE the final cluster (a boundary link): j does not move
                #    in the reverse move either, so the reverse probability is the
                #    link probability of the post-move state, from the energy with
                #    m shifted by +dr (Ef): 1 - exp(-beta (E0 - Ef)).
                # Every tested link carries one factor in the acceptance ratio
                # (see _vmmc_log_acceptance); the pre-1.0.8 code dropped the
                # frustrated-link factor and violated detailed balance.
                p_r_internal = _vmmc_link_probability(beta, Er.get(j, 0.0) - e0)
                p_r_boundary = _vmmc_link_probability(beta, e0 - Ef.get(j, 0.0))

                if p_f > 0.0 and random.random() < p_f:
                    formed_links.append((p_f, p_r_internal))
                    if j not in cluster:
                        if j in frozen_set:
                            return (latticeObject, current_energy, False, len(cluster))   # cannot move a frozen chain
                        cluster.add(j)
                        if len(cluster) > n_c:
                            return (latticeObject, current_energy, False, len(cluster))   # exceeded the cutoff -> reject
                        queue.append(j)
                else:
                    failed_links.append((j, p_f, p_r_boundary, p_r_internal))

        # --- apply the rigid translation; reject on hard-core / hardwall clash ----
        old_positions = {}
        for c in cluster:
            old_positions[c] = latticeObject.chains[c].get_ordered_positions()
            lattice_utils.delete_chain_by_position(old_positions[c], latticeObject.grid, c)

        new_positions = {}
        placed = []
        for c in cluster:
            translated = []
            clash = False
            for pos in old_positions[c]:
                raw_tpos = [pos[d] + dr[d] for d in range(nd)]
                if hardwall and any(raw_tpos[d] < 0 or raw_tpos[d] >= dimensions[d]
                                    for d in range(nd)):
                    clash = True
                    break
                tpos = raw_tpos if hardwall else lattice_utils.pbc_convert(raw_tpos, dimensions)
                if lattice_utils.get_gridvalue(tpos, latticeObject.grid) != 0:
                    clash = True
                    break
                lattice_utils.set_gridvalue(tpos, c, latticeObject.grid)
                translated.append(tpos)

            if clash:
                lattice_utils.delete_chain_by_position(translated, latticeObject.grid, c)
                for cc in placed:
                    lattice_utils.delete_chain_by_position(new_positions[cc], latticeObject.grid, cc)
                for cc in cluster:
                    lattice_utils.place_chain_by_position(old_positions[cc], latticeObject.grid, cc, safe=True)
                return (latticeObject, current_energy, False, len(cluster))

            new_positions[c] = translated
            placed.append(c)

        # --- commit to the type_grid + chain objects, then get the exact dE -------
        for c in cluster:
            latticeObject.delete_chain_from_type_grid(c, old_positions[c], list(range(len(old_positions[c]))), safe=True)
        for c in cluster:
            latticeObject.chains[c].set_ordered_positions(new_positions[c])
            latticeObject.insert_chain_into_type_grid(c, new_positions[c], list(range(len(new_positions[c]))), safe=True)

        E_after = hamiltonianObject.evaluate_total_energy(latticeObject)[0]
        dE      = E_after - current_energy

        # --- VMMC (Metropolis-Hastings) acceptance --------------------------------
        #   acc = min(1, exp(-beta*dE) * PROD_formed (p_r/p_f)
        #                              * PROD_failed ((1-p_r)/(1-p_f)))
        # with p_r the boundary probability for partners left outside the cluster
        # and the internal one for partners inside it (frustrated links). A reverse
        # probability of 0 (formed) or 1 (failed) makes the reverse move impossible
        # -> reject.
        log_ratio = _vmmc_log_acceptance(beta, dE, cluster, formed_links, failed_links)

        accept = False
        if log_ratio is not None:
            if log_ratio >= 0.0 or random.random() < math.exp(log_ratio):
                accept = True

        if accept:
            return (latticeObject, E_after, True, len(cluster))

        # --- reject: revert to the original state ---------------------------------
        for c in cluster:
            lattice_utils.delete_chain_by_position(new_positions[c], latticeObject.grid, c)
            latticeObject.delete_chain_from_type_grid(c, new_positions[c], list(range(len(new_positions[c]))), safe=True)
        for c in cluster:
            latticeObject.chains[c].set_ordered_positions(old_positions[c])
            lattice_utils.place_chain_by_position(old_positions[c], latticeObject.grid, c, safe=True)
            latticeObject.insert_chain_into_type_grid(c, old_positions[c], list(range(len(old_positions[c]))), safe=True)
        return (latticeObject, current_energy, False, len(cluster))


    #-----------------------------------------------------------------
    #     JUMP-AND-RELAX (single-chain composite move, code 13)
    def jump_and_relax_move(self, chain_to_move, latticeObject, current_energy, acceptanceObject, hamiltonianObject, cs_substeps, cs_mode, hardwall=False):
        """
        Jump-and-relax single-chain move (move code 13). Self-contained.

        Relaxes a chain, attempts to relocate it, then relaxes it again. The move
        is built from three sub-steps that EACH individually preserve the Boltzmann
        distribution, so their composition does too (a sequence of pi-preserving
        Monte Carlo updates preserves pi):

        1. **relax** - a single-chain crankshaft shake (:meth:`single_chain_shake`;
           many local perturbations of this chain, each accepted/rejected by its own
           Metropolis criterion inside the kernel); always committed.
        2. **jump** - a rigid translation of the whole chain (:meth:`chain_translate`),
           accepted or rejected on its OWN Metropolis criterion (and reverted on a
           hard-sphere clash); a standard, detailed-balanced single-chain
           translation.
        3. **relax** - a second single-chain shake; always committed.

        The move concentrates sampling effort on relocating one chain and letting it
        settle into its (possibly new) environment.

        .. note::

           Earlier versions deferred a single accept/reject to the energy *after*
           both relaxations, which is an asymmetric proposal and broke detailed
           balance (it over-accepted downhill moves). Accepting/rejecting the jump
           on its own merit, between two pi-preserving relaxations, fixes this. For
           aggressive relocation through dense/condensed phases prefer
           :meth:`vmmc_move` or :meth:`system_pull`.

        Parameters
        ----------
        chain_to_move : Chain
            The (uniformly selected) chain object to move.

        latticeObject : Lattice
            The system lattice (mutated in place).

        current_energy : float
            The current total system energy.

        acceptanceObject : AcceptanceCalculator
            Provides the Metropolis acceptance test and the auxiliary-move logging.

        hamiltonianObject : Hamiltonian
            Used by the relaxations and for the from-scratch total-energy recompute
            that gives the exact jump dE.

        cs_substeps : int
            Number of crankshaft sub-moves per relaxation (``CRANKSHAFT_SUBSTEPS``).

        cs_mode : str
            Crankshaft mode passed through to :meth:`single_chain_shake`.

        hardwall : bool, optional
            Whether the lattice has hard-wall (non-periodic) boundaries. Default
            is False.

        Returns
        -------
        (Lattice, float, bool)
            The (mutated) lattice, the new total energy, and whether the jump
            (step 2) was accepted.
        """
        chainID = chain_to_move.chainID

        # [STEP 1] relax the chain in place (pi-preserving; committed unconditionally)
        (_, energy, proposed1, _) = self.single_chain_shake(chainID, latticeObject, current_energy,
                                                            acceptanceObject, hamiltonianObject,
                                                            cs_substeps, cs_mode, hardwall)

        # [STEP 2] propose a rigid jump and accept/reject it on its own Metropolis
        # criterion - a standard, detailed-balanced single-chain translation.
        # chain_translate leaves the chain in its new GRID position (or reverts on a
        # hard-sphere clash); we then sync the type_grid + chain object, evaluate the
        # exact dE from scratch, and either keep or fully revert it.
        jump_accepted = False
        (move_event, success) = self.chain_translate(chain_to_move, latticeObject.grid, hardwall)
        if success:
            old_positions = move_event.original_positions
            new_positions = move_event.moved_positions
            indices       = move_event.moved_indices

            latticeObject.update_type_grid(chainID, old_positions, new_positions, indices, safe=True)
            chain_to_move.set_ordered_positions(new_positions)

            E_after = hamiltonianObject.evaluate_total_energy(latticeObject)[0]
            if acceptanceObject.boltzmann_acceptance(energy, E_after):
                energy = E_after
                jump_accepted = True
            else:
                # revert the jump (grid, type_grid and chain object) back to pre-jump
                lattice_utils.delete_chain_by_position(new_positions, latticeObject.grid, chainID)
                lattice_utils.place_chain_by_position(old_positions, latticeObject.grid, chainID, safe=True)
                latticeObject.update_type_grid(chainID, new_positions, old_positions, indices, safe=True)
                chain_to_move.set_ordered_positions(old_positions)

        # [STEP 3] relax the chain again in its (possibly new) location (pi-preserving; committed)
        (_, energy, proposed2, _) = self.single_chain_shake(chainID, latticeObject, energy,
                                                            acceptanceObject, hamiltonianObject,
                                                            cs_substeps, cs_mode, hardwall)

        # the shake sub-moves are auxiliary-Markov-chain MC moves (throughput accounting)
        acceptanceObject.alt_Markov_chain_update_move_logs(proposed1 + proposed2)
        return (latticeObject, energy, jump_accepted)


    #-----------------------------------------------------------------
    #
    def chain_translate(self, ChainToMove, lattice, hardwall=False):
        """
        The chain_translate move allows the full chain to be translated in rigid body
        space around the lattice.

        The cost of this move will scale linearly with chain length (note cost comes 
        from the energy evaluation). However this will often be one of the most expensive
        moves which can be made as we have to evaluate both short range and (if relevant)
        long range interactions for EVERY residue in the chain when re-calculating the 
        new energy.
    
        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        MoveType code: 2

        Parameters
        ----------
        ChainToMove : Chain
            The chain object to be translated. Treated as read-only (its
            positions are read but not modified here).

        lattice : numpy.ndarray
            The occupancy grid itself (``latticeObject.grid``, not the Lattice
            object): an ``int32`` array of the box dimensions holding the chainID
            occupying each site, or 0 for solvent. Mutated in place to
            reflect the new positions if the move succeeds.

        hardwall : bool, optional
            If True the move is rejected when the translated chain straddles a
            periodic boundary (enforcing a hard wall). Default is False.

        Returns
        -------
        tuple
            ``(MoveEvent, True)`` if the move was made (the MoveEvent describes
            the change for downstream energy evaluation), or ``(False, False)``
            if the move was rejected (hard-sphere clash or hardwall violation),
            in which case the lattice is left unchanged.
        """

        chainID         = ChainToMove.chainID
        chain_positions = ChainToMove.get_ordered_positions()
        dimensions      = lattice_utils.get_dimensions(lattice)
        num_dims        = len(dimensions)

        # delete the chain from the lattice
        lattice_utils.delete_chain_by_position(chain_positions, lattice, chainID)

        # define translation operations through a translation vector which
        # we apply to each chain unit
        offset_vector = []
        for i in range(0, num_dims):
            offset_vector.append(random.randint(0, dimensions[i]-1))

        translated_positions = []

        # for the position of each residue
        for position in chain_positions:
            translated_pos = []

            # for each dimension incremement the position
            for dim in range(0, num_dims):
                translated_pos.append(position[dim] + offset_vector[dim] )

            # carry out periodic boundary correction on the new position
            translated_pos = lattice_utils.pbc_convert(translated_pos, dimensions)

            # check if that position is empty - as soon as we find a position
            # in the lattice which is not empty consider the move rejected.
            # This means we *only* carry out as many transation operations as
            # absolutely necessary
            if not lattice_utils.get_gridvalue(translated_pos, lattice) == 0:

                ## REVERT BACK !!
                
                # Delete the positions we insterted so far
                lattice_utils.delete_chain_by_position(translated_positions, lattice, chainID)

                # re-insert the chain back into its old positions
                lattice_utils.place_chain_by_position(chain_positions, lattice, chainID, safe=True)
                                
                return (False, False)
                #return (False, False, False, False, False)

            # if the position was free update the lattice copy object
            lattice_utils.set_gridvalue(translated_pos, chainID, lattice)
            
            # and add the position to the growing list
            translated_positions.append(translated_pos)

        # If we're outside the for-loop translation operation was a success!

        # Now check for hardwall rules
        if hardwall:
            if lattice_utils.do_positions_stradle_pbc_boundary(translated_positions):
                
                lattice_utils.delete_chain_by_position(translated_positions, lattice, chainID)
                lattice_utils.place_chain_by_position(chain_positions, lattice, chainID, safe=True)                                
                return (False, False)
                

        ## Create the MoveEvent object
        # We moved all the residues so moved_positions and full_moved_chain_positions
        # are the same.        
        ME = MoveEvent(original_positions        = chain_positions,
                       moved_positions           = translated_positions,
                       original_chain_positions  = chain_positions,
                       moved_chain_positions     = translated_positions,                
                       moved_indices             = list(range(0,len(chain_positions))),
                       move_type                 = 2)
                               
        return (ME, True)


    

    #-----------------------------------------------------------------
    #    
    def chain_rotate(self, ChainToMove, lattice, hardwall=False):
        """
        The chain_rotate move rotates the full chain in rigid-body space about
        one of its own beads - the bead nearest the chain's (single-image)
        centroid. The chain's displacement vectors relative to that pivot bead
        are rotated by a random cardinal rotation and re-anchored on the pivot
        bead's original in-box position. Anchoring on a bead (rather than the
        rounded centre of mass, the pre-1.0.8 behaviour) makes every rotation
        exactly invertible for any box shape, which detailed balance requires -
        see the in-body comment for the full rationale.
              
        The cost of this move will scale linearly with chain length (note cost comes 
        from the energy evaluation).

        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        A draw that maps the chain exactly onto itself - in 3D, any of the three
        rotations of an axis-aligned straight chain about its own axis - is
        rejected as a null move rather than returned as a successful proposal.
        The sampled ensemble is identical either way (an accepted identity
        proposal and a rejected one leave the same state behind), but rejecting
        it keeps ACCEPTANCE.dat a count of moves that actually changed the
        configuration, and skips an energy evaluation that cannot change
        anything.

        MoveType code: 3

        Parameters
        ----------
        ChainToMove : Chain
            The chain object to be rotated. Treated as read-only.

        lattice : numpy.ndarray
            The occupancy grid itself (``latticeObject.grid``, not the Lattice
            object): an ``int32`` array of the box dimensions holding the chainID
            occupying each site, or 0 for solvent. Mutated in place to
            reflect the rotated positions if the move
            succeeds.

        hardwall : bool, optional
            If True the move is rejected when the rotated chain straddles a
            periodic boundary. Default is False.

        Returns
        -------
        tuple
            ``(MoveEvent, True)`` if the rotation was made, or
            ``(False, False)`` if rejected (a singleton chain, a draw that maps
            the chain exactly onto itself, a hard-sphere clash or a hardwall
            violation), in which case the lattice is left unchanged.
        """
        ## A note on rotations and offset. The offset parameter is calculated here
        ## so the chain can be first converted into a single image and then rotated
        ## (rather than rotating something that exists in periodic image space). This
        ## is actually irrelevant if you're using a box or square and the simulation
        ## vessel, but when vertices are unequal in length PBC conditions break and
        ## so this ensures we can still rotate chains in rectangular boxes


        chainID                  = ChainToMove.chainID
        chain_positions_original = ChainToMove.get_ordered_positions()
        if len(chain_positions_original) < 2:
            return (False, False)
        chain_positions          = ChainToMove.get_single_image_positions()
        dimensions               = lattice_utils.get_dimensions(lattice)
        num_dims                 = len(dimensions)


        # delete the chain from the lattice
        lattice_utils.delete_chain_by_position(chain_positions_original, lattice, chainID)

        # Rotation centre and reversibility.
        #
        # We rotate the chain as a rigid body about one of its own BEADS - the bead
        # nearest the (single-image) centroid - rather than about the rounded
        # circular-mean COM. Concretely we rotate the DISPLACEMENT VECTORS of the
        # single-image chain relative to that pivot bead, then re-anchor them on the
        # pivot bead's ORIGINAL (in-box, on-lattice) position and re-wrap:
        #
        #     rotated[i] = pbc_convert( raw_pivot + R( si[i] - si[pivot] ) )
        #
        # This is exactly reversible for ANY box (cubic or not) and for
        # PBC-straddling chains: R is only ever applied to integer displacement
        # vectors, and the anchor raw_pivot is an in-box lattice point that maps to
        # itself, so the reverse move (which re-single-images, finds the same physical
        # pivot bead, and applies R^-1) returns the original configuration bit for
        # bit. The old code rotated about a rounded COM and added a bead-0 offset;
        # the rounded COM of the rotated chain is generally a different lattice point,
        # so the inverse rotation did not undo the move and detailed balance was
        # violated (verified by exact enumeration: an equilibrium bias in
        # translate+rotate sampling that vanishes once rotations are made reversible).
        # Pivot selection must be EXACT-INTEGER arithmetic: with a float centroid,
        # two beads exactly equidistant from the centre get their tie broken by
        # sub-ulp rounding noise that is NOT invariant under the rotation, so the
        # reverse move could pick the OTHER tied bead and fail to invert (a
        # complete detailed-balance violation on tied shapes). Comparing
        # ||n*p_i - sum(p)||^2 in integers is tie-stable: rotation preserves the
        # whole distance vector and the bead order, so first-index argmin picks
        # the same physical bead in both directions.
        _si = np.asarray(chain_positions, dtype=np.int64)
        _d2 = ((len(_si) * _si - _si.sum(axis=0)) ** 2).sum(axis=1)
        _pivot_idx = int(np.argmin(_d2))
        _si_pivot = _si[_pivot_idx]
        _raw_pivot = chain_positions_original[_pivot_idx]

        # displacement of every bead from the pivot, in the single (unwrapped) image
        OC_positions      = [[int(p[d] - _si_pivot[d]) for d in range(num_dims)]
                             for p in chain_positions]
        rotated_positions = []

        ## 2D rotation
        if num_dims == 2:
            OC_rotated_positions = lattice_utils.rotate_positions_2D(OC_positions, [90,180,270][random.randint(0,2)])
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert(
                    [position[0] + _raw_pivot[0], position[1] + _raw_pivot[1]], dimensions))

        ## 3D rotation
        if num_dims == 3:
            OC_rotated_positions = lattice_utils.rotate_positions_3D(OC_positions, ['x','y','z'][random.randint(0,2)], [90,180,270][random.randint(0,2)])
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert(
                    [position[0] + _raw_pivot[0], position[1] + _raw_pivot[1], position[2] + _raw_pivot[2]], dimensions))

        # Guaranteed-null draws are REJECTED, not accepted.
        #
        # Some draws map the chain exactly onto itself: in 3D, the three rotations
        # about the axis of an axis-aligned straight chain (3 of the 9 draws). These
        # used to be returned as successful proposals, sent through the full
        # energy evaluation, accepted with dE = 0 and counted in ACCEPTANCE.dat -
        # so that file reported moves that provably did nothing as successes and
        # the acceptance ratio a user computes from ACCEPTANCE.dat / MOVE_FREQS.dat
        # was not the fraction of steps that changed anything.
        #
        # Rejecting cannot change what is sampled. Detailed balance constrains only
        # transitions between DIFFERENT states (for y = x both sides of
        # pi(x)P(x->y) = pi(y)P(y->x) are identically pi(x)P(x->x)), and accepting
        # an identity proposal and rejecting it leave the system in the same state,
        # so the trajectory is bit-identical either way. Nor does it perturb the RNG
        # stream: the axis/angle draws have already been made above, and the
        # acceptance test that follows a null would take the dE <= 0 branch without
        # calling random().
        #
        # The comparison must be over ORDERED positions - see _is_identity_proposal.
        if _is_identity_proposal(rotated_positions, chain_positions_original, num_dims):
            # nothing has been inserted yet, so putting the chain back is the whole revert
            lattice_utils.place_chain_by_position(chain_positions_original, lattice, chainID, safe=True)
            return (False, False)

        # Now check for hardwall rules
        if hardwall:
            if lattice_utils.do_positions_stradle_pbc_boundary(rotated_positions):                

                # no need to delete anything because nothing inserted yet
                
                lattice_utils.place_chain_by_position(chain_positions_original, lattice, chainID, safe=True)                                
                return (False, False)
            else:
                pass


        # having built a new list of rotated positions let's see if any of them clash. Note the inserted_chain
        # keeps track of what's going on so if we find a clash we only have to cycle over a small number of filled
        # positions to delete the part of the chain we inserted
        inserted_chain = []
        for rotated_pos in rotated_positions:

            # if it turns out a position was already occupied
            if not lattice_utils.get_gridvalue(rotated_pos, lattice) == 0:

                # delete whatever part(s) of the chain we've already added
                lattice_utils.delete_chain_by_position(inserted_chain, lattice, chainID)

                # re-insert the chain back where it was
                lattice_utils.place_chain_by_position(chain_positions_original, lattice, chainID, safe=True)

                # return all the failure
                return (False, False)
            else:
                
                # if the position was free update the lattice copy object
                lattice_utils.set_gridvalue(rotated_pos, chainID, lattice)
                inserted_chain.append(rotated_pos)
            
        # if we get here we have succesfully added all the rotated positions to the lattice.
        # Assume all positions moved (maybe they didn't but determining this ends up being 
        # more computationally expensive
        ME = MoveEvent(original_positions        = chain_positions_original,
                       moved_positions           = rotated_positions,
                       original_chain_positions  = chain_positions_original,
                       moved_chain_positions     = rotated_positions,
                       moved_indices             = list(range(0,len(chain_positions))),
                       move_type                 = 3)

        return (ME, True)
    


    #-----------------------------------------------------------------
    #        
    def chain_pivot(self, ChainToMove, lattice, pivotPoint_range=None, hardwall=False):
        """
        The chain_pivot move allows part of the chain to perform a rigid
        pivot. Specifically, we randomly select (with uniform probability)
        some position through the chain, and then pivot the shorter half in
        a rigid-body manner using the selected position as an anchor point.

        NOTE: this is only applied to chains where are 3 residues or longer
        - the move is automatically rejected for shorter chains.
        
        The cost of this move will scale linearly with chain length (note cost 
        comes  from the energy evaluation).

        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        !!! WARNING !!!
        The pivotPoint_range allows for the user to define a range of positions
        which can be pivoted. THIS FUNCTIONALITY IS NOT GENERALIZABLE YET - it
        was implemented for a specific use case, and while generalizing it is
        not a major issue this has not yet been done - PLEASE DO NOT USE.

        MoveType code: 4

        Parameters
        ----------
        ChainToMove : Chain
            The chain object to be pivoted. Treated as read-only. Chains shorter
            than 3 residues are automatically rejected.

        lattice : numpy.ndarray
            The occupancy grid itself (``latticeObject.grid``, not the Lattice
            object): an ``int32`` array of the box dimensions holding the chainID
            occupying each site, or 0 for solvent. Mutated in place to
            reflect the pivoted positions if the move
            succeeds.

        pivotPoint_range : list or None, optional
            Optional list of candidate pivot indices to draw from (see the
            warning above - NOT generalizable, do not use). If None (default) a
            uniformly random interior position is selected and the shorter half
            of the chain is pivoted.

        hardwall : bool, optional
            If True the move is rejected when the pivoted segment straddles a
            periodic boundary. Default is False.

        Returns
        -------
        tuple
            ``(MoveEvent, True)`` if the pivot was made, or ``(False, False)``
            if rejected (chain too short, hard-sphere clash, or hardwall
            violation), in which case the lattice is left unchanged.

        Raises
        ------
        Exception
            If the reconstructed pivoted chain length does not match the
            original chain length (an internal consistency check).
        """
    
        chainID         = ChainToMove.chainID
        chain_positions = ChainToMove.get_ordered_positions()
        chain_length    = len(chain_positions)
        dimensions      = lattice_utils.get_dimensions(lattice)
        num_dims        = len(dimensions)

        # reject if we have a chain 2 or less in length 
        if chain_length < 3:
            return (False, False)
            
        # select a position along the chain to pivot
        pivot_point = random.randint(1,len(chain_positions)-2)
                
        # if no pivot point range was provided (as default)
        if pivotPoint_range is None:

            # Pivot the SHORTER arm about the bead at pivot_point (which stays
            # fixed). The left arm is beads [0, pivot_point-1] (pivot_point beads);
            # the right arm is beads [pivot_point+1, L-1] (L-1-pivot_point beads).
            # Anchoring on bead pivot_point in BOTH cases makes the two termini
            # symmetric and removes the pivot_point==1 null move: the old else branch
            # anchored on bead pivot_point-1, so for pivot_point==1 it "rotated" only
            # the anchor bead - a guaranteed no-op that was still fully energy
            # evaluated and logged as an accepted pivot. For L==3 that was the ONLY
            # possible pivot_point, so chain_pivot never actually moved a 3-mer, and
            # for even L the C-terminal bead could never move. Compare against
            # (L-1)/2 so the C-terminal bead is reachable for even L too.
            if pivot_point > (chain_length - 1) / 2:

                # set if we're going to to rotate positions and then
                # add them to the end of a fixed region
                add_to_end=True

                # positions which will be rotated (beads after the pivot bead)
                positions_to_rotate    = chain_positions[pivot_point:]

                # positions which will be held fixed (0 to pivot point)
                positions_held_fixed   = chain_positions[:pivot_point]

                # sequence indices of positions which will be rotated
                indices = list(range(pivot_point, len(chain_positions)))

            else:
                # variable roles given above; anchor is bead pivot_point (the LAST
                # element of positions_to_rotate), beads 0..pivot_point-1 swing.
                add_to_end=False
                positions_to_rotate  = chain_positions[:pivot_point + 1]
                positions_held_fixed = chain_positions[pivot_point + 1:]
                indices = list(range(0, pivot_point + 1))

        else:
            # randomly select a position from the pivot point range - right now we're always
            # going to treat the pivot point as defining a region C-terminal which is pivoted while
            # the N-terminal region remains fixed
            pivot_point = pivotPoint_range[random.randint(0,len(pivotPoint_range)-1)]
           
            # set if we're going to to rotate positions and then
            # add them to the end of a fixed region
            add_to_end=True

            # positions which will be rotated (pivot point to end)
            positions_to_rotate    = chain_positions[pivot_point:]

            # positions which will be held fixed (0 to pivot point)
            positions_held_fixed   = chain_positions[:pivot_point]

            # sequence indices of positions which will be rotated
            indices = list(range(pivot_point, len(chain_positions)))
                              
        
        # delete the positions we're going to rotate but keep the rest
        lattice_utils.delete_chain_by_position(positions_to_rotate, lattice, chainID)

        # get the head position as a separate copy - this is the pivot position which should
        # remain fixed as we perform lever arm rotation on the rest of the residues in the
        # positions_to_rotate - note that depending on which side we're rotating we set
        # first or last residue as head - e.g. see diagram below
        #
        # - : residue reminaing fixe
        # p : pivot point
        # x : residue to pivot
        #
        # if add_to_end
        #       -...---------PXXXXX...X 
        # else
        #   XXX...XXXXXXXP------...---



        if add_to_end:
            head_position = positions_to_rotate[0][:]
        else:
            head_position = positions_to_rotate[-1][:]

        # this offset is, again, so we can use PBC boxes with non-equal sides (see intro in the rotate
        # _chain move
        original_positions_to_rotate = positions_to_rotate
        positions_to_rotate = lattice_utils.convert_chain_to_single_image(positions_to_rotate, dimensions)

        # translate the shorter halve's positions to the origin
        if num_dims == 2:
            x_ref = head_position[0]
            y_ref = head_position[1]
            
            # build a list of origin centered positions
            OC_positions = []
            for position in positions_to_rotate:
                OC_positions.append([position[0] - x_ref, position[1] - y_ref])

            # carry out a random rotation in 2D
            OC_rotated_positions = lattice_utils.rotate_positions_2D(OC_positions, [90,180,270][random.randint(0,2)])

            #
            # determine rotation offset - basically when we do the rotationn we want the residue which connects BACK to the chain to be in exactly
            # the same position - i.e.
            #
            #  XXP
            #    XXXX
            #
            #    PXX
            #    XXXX
            #
            # We have to make sure P (the pivot point) remians in the same position - if we're pivoting the first half the pivot point is at the
            # end of the OC_rotated_positions, while if we're pivoting the second half the pivot point is at the front of the OC_rotated_positions        
            #

            if add_to_end:
                x_return = x_ref - OC_rotated_positions[0][0]
                y_return = y_ref - OC_rotated_positions[0][1]
            else:
                x_return = x_ref - OC_rotated_positions[-1][0]
                y_return = y_ref - OC_rotated_positions[-1][1]
            
            # now move all the positions back again from the origin such that they move back to link up with the chain
            rotated_positions = []
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert([position[0] + x_return, position[1] + y_return], dimensions))

        ## 3D rotation
        if num_dims == 3:
            x_ref = head_position[0]
            y_ref = head_position[1]
            z_ref = head_position[2]

            # center the pivot section at the origin
            OC_positions = []
            for position in positions_to_rotate:
                OC_positions.append([position[0] - x_ref, position[1] - y_ref, position[2] - z_ref])
                
            # carry out a random rotation in 3D
            OC_rotated_positions = lattice_utils.rotate_positions_3D(OC_positions, ['x','y','z'][random.randint(0,2)], [90,180, 270][random.randint(0,2)])
            
            # determine rotation offset (sometimes rotation around the 0 axis will still move the head position - see 2D description
            # for more details!
            if add_to_end:
                return_correction = [x_ref - OC_rotated_positions[0][0], y_ref - OC_rotated_positions[0][1], z_ref - OC_rotated_positions[0][2]]
            else:
                return_correction = [x_ref - OC_rotated_positions[-1][0], y_ref - OC_rotated_positions[-1][1], z_ref - OC_rotated_positions[-1][2]]

            # now move all the positions back again from the origin
            rotated_positions = []
            for position in OC_rotated_positions:
                rotated_positions.append(lattice_utils.pbc_convert([position[0] + return_correction[0], position[1] + return_correction[1], position[2] + return_correction[2]], dimensions))
                


        # Now check for hardwall rules
        if hardwall:
            if lattice_utils.do_positions_stradle_pbc_boundary(rotated_positions):                
                # no need to delete anything because nothing inserted yet
                lattice_utils.place_chain_by_position(original_positions_to_rotate, lattice, chainID, safe=True)                                
                return (False, False)


        # having built a new list of rotated positions let's see if any of them clash. Note the inserted_chain
        # keeps track of what's going on so if we find a clash we only have to cycle over a small number of filled
        # positions to delete the part of the chain we inserted
        inserted_chain = []
        for rotated_pos in rotated_positions:

            # if it turns out a position was already occupied
            if not lattice_utils.get_gridvalue(rotated_pos, lattice) == 0:

                # delete whatever part(s) of the chain we've already added
                lattice_utils.delete_chain_by_position(inserted_chain, lattice, chainID)

                # re-insert the section of chain we deleted previously
                lattice_utils.place_chain_by_position(original_positions_to_rotate, lattice, chainID, safe=True)

                # return all the failure
                return (False, False)
            else:                
                # if the position was free update the lattice copy object
                lattice_utils.set_gridvalue(rotated_pos, chainID, lattice)
                inserted_chain.append(rotated_pos)
        
        # if we get here we succesfully pivoted the chain!
        if add_to_end:
            # add the newly rotated positions to the end of the already
            # known positions
            fully_pivoted_chain = positions_held_fixed + inserted_chain
        else:
            fully_pivoted_chain = inserted_chain + positions_held_fixed

        if not len(fully_pivoted_chain) == len(chain_positions):
            raise Exception('Yeah stop right there...')
            
        ME = MoveEvent(original_positions        = original_positions_to_rotate,
                       moved_positions           = inserted_chain,
                       original_chain_positions  = chain_positions,
                       moved_chain_positions     = fully_pivoted_chain,
                       moved_indices             = indices,
                       move_type                 = 4,
                       pivot_point               = pivot_point)
                   

        return (ME, True)                

    #-----------------------------------------------------------------
    #    
    def head_pivot(self, ChainToMove, lattice, hardwall=False):
        """
        The head_pivot move allows the chain 'head' (either the first or last residue)
        to pivot in some direction.

        This is always going to be very cheap as it 'always' translates to a single
        position change irrespective of chain length. We arbitrarily pick the first 
        or last residue to pivot - i.e. chains have 2 heads. 

        This is - honestly - kind of a stupid move. It was one of the first moves I
        coded up as a very simple and easy to debug move, but is probably not going
        to add much. However, the chain_pivot move won't pivot the ends of chains
        if you have a very short chain so it does actually serve a relevant purpose!

        The move is rejected if there's a hard-sphere clash, else we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        MoveType code: 5

        Parameters
        ----------
        ChainToMove : Chain
            The chain object whose head (first or last residue, chosen at
            random) is to be pivoted. Treated as read-only.

        lattice : numpy.ndarray
            The occupancy grid itself (``latticeObject.grid``, not the Lattice
            object): an ``int32`` array of the box dimensions holding the chainID
            occupying each site, or 0 for solvent. Mutated in place to
            reflect the new head position if the move
            succeeds.

        hardwall : bool, optional
            If True the move is rejected when the moved head and its neighbour
            straddle a periodic boundary. Default is False.

        Returns
        -------
        tuple
            ``(MoveEvent, True)`` if the head pivot was made, or
            ``(False, False)`` if rejected (the chain has no pivotable head,
            the head landed on its original site, a hard-sphere clash, or a
            hardwall violation), in which case the lattice is left unchanged.
        """
        chainID           = ChainToMove.chainID
        chain_positions   = ChainToMove.get_ordered_positions()
        if len(chain_positions) < 2:
            return (False, False)
        dimensions        = lattice_utils.get_dimensions(lattice)
        num_dims          = len(dimensions)

        # copy because we want to create a new list of positions
        updated_positions = chain_positions[:] 

        # select one end of the chain as the head
        if random.random() > 0.5:
            
            ## **************************
            ## Working with first residue
        
            # Delete the first residue
            lattice_utils.delete_residue(chain_positions[0], lattice, chainID)

            # get possible new positions by building a list of the sites which are adajcent
            # to the second residue in the chain (having just deleted the first)
            if num_dims == 2:
                possible_positions = lattice_utils.get_adjacent_sites_2D(chain_positions[1][0], chain_positions[1][1], dimensions)
            else:
                possible_positions = lattice_utils.get_adjacent_sites_3D(chain_positions[1][0], chain_positions[1][1],chain_positions[1][2], dimensions)

            # randomly select one of the positions from this list
            possible_position = list(possible_positions[random.randint(0, len(possible_positions)-1)])
                

            # if hardwall boundary
            if hardwall:          

                # do the first and second positions now cross a PBC boundary?
                if lattice_utils.do_positions_stradle_pbc_boundary([possible_position, chain_positions[1]]):                

                    # revert back to original position
                    lattice_utils.insert_residue(chain_positions[0], lattice, chainID)
                    return (False, False)



            # if 'moved' the residue to the same position the original residue came from...
            # no move (same thing) - re insert and return false 
            # NOTE: Philiosphically should this be a rejection or not? I really don't know...
            if lattice_utils.same_sites(possible_position, chain_positions[0]):

                # revert to original position
                lattice_utils.insert_residue(possible_position, lattice, chainID)
                return (False, False)
                                
            # if moved into occupied site then reject
            elif not lattice_utils.get_gridvalue(possible_position, lattice) == 0.0:
            
                # revert to original position
                lattice_utils.insert_residue(chain_positions[0], lattice, chainID)
                return (False, False)

            # else move is A-OK!
            else:               
                
                # update the lattice 
                lattice_utils.insert_residue(possible_position, lattice, chainID)
                
                # set the first residue to the new position in the copy list
                updated_positions[0] = possible_position                

                # indicies represent the first residue
                ME = MoveEvent(original_positions        = [chain_positions[0]],
                               moved_positions           = [possible_position],
                               original_chain_positions  = chain_positions,
                               moved_chain_positions     = updated_positions,
                               moved_indices             = [0],
                               move_type                 = 5)                       
        
                return (ME, True)

        else:
            ## Working with last residue            

            # Delete last residue
            lattice_utils.delete_residue(chain_positions[-1], lattice, chainID)

            # get possible new positions 
            if num_dims == 2:
                possible_positions = lattice_utils.get_adjacent_sites_2D(chain_positions[-2][0], chain_positions[-2][1], dimensions)
            else:
                possible_positions = lattice_utils.get_adjacent_sites_3D(chain_positions[-2][0], chain_positions[-2][1], chain_positions[-2][2], dimensions)

            # randomly select one of the positions from this list
            possible_position = list(possible_positions[random.randint(0, len(possible_positions)-1)])


            # if hardwall boundary
            if hardwall:          

                # do the second to last and last positions now cross a PBC boundary?
                if lattice_utils.do_positions_stradle_pbc_boundary([chain_positions[-2],possible_position]):                

                    # revert back to original position
                    lattice_utils.insert_residue(chain_positions[-1], lattice, chainID)
                    return (False, False)
                
            # if 'moved' the same position head came from...
            # no move (same thing) - re insert and return false
            if lattice_utils.same_sites(possible_position, chain_positions[-1]):

                # revert to original position
                lattice_utils.insert_residue(possible_position, lattice, chainID)
                return (False, False)
                
            # if moved into occupied site then reject
            elif not lattice_utils.get_gridvalue(possible_position, lattice) == 0.0:
            
                # revert to original position
                lattice_utils.insert_residue(chain_positions[-1], lattice, chainID)
                return (False, False)

            else:

                lattice_utils.insert_residue(possible_position, lattice, chainID)

                # set the last residue to the new position in the copy list
                updated_positions[-1] = possible_position

                # indices represent the terminal residue
                ME = MoveEvent(original_positions        = [chain_positions[-1]],
                               moved_positions           = [possible_position],
                               original_chain_positions  = chain_positions,
                               moved_chain_positions     = updated_positions,                               
                               moved_indices             = [len(updated_positions)-1],
                               move_type                 = 5)
                       
        
                return (ME, True)





    #-----------------------------------------------------------------
    #    
    def cluster_translate(self, selected_chain, latticeObject, cluster_move_threshold=None, cluster_size_threshold=None, hardwall=False, frozen_chains=()):
        """
        The cluster_translate move allows a connected components (cluster) to be 
        translated in rigid body space around the lattice.

        The cost of this move grows linearly as clusters get big - the 
        cluster_threshold variable facilitates a soft threshold on the max cluster 
        size.

        To ensure detailed balance moves a cluster move MUST NOT lead to the 
        incorporation of new chains into the cluster being moved. Explicitly (thank you 
        Tyler!), if a cluster translation leads to two clusters merging, there is 
        no move in our which would allow that cluster to unmerge again in a single 
        move - i.e. a cluster-merging translation move is irreversible, hence breaking
        detailed balance.

        On the plus side, this also means that cluster moves which we can accept (i.e. 
        which does not lead to a cluster merger or clash) must be energy neutral, 
        so we don't have to run any short-range energy evaluations on it. Note that we
        must still run long-range energy calculations, meaning for charged systems 
        this becomes particularly expensive...

        The move is rejected if there's a hard-sphere clash, *or* we change the cluster
        size after the move (i.e. incorporate new residues in). If not rejected  we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        The cluster_size_threshold is default to None, but in the simulations.py file we
        set this to be such that in the case that a cluster contains ALL the chains it is not
        rotated, otherwise it can be.

        MoveType code: 7

        Parameters
        ----------
        selected_chain : Chain
            A chain belonging to the cluster to be moved; its connected component
            defines the cluster.

        latticeObject : Lattice
            The lattice object containing the chains. Its grid is mutated in
            place if the move succeeds.

        cluster_move_threshold : int or None, optional
            Maximum per-dimension step size of the translation. If None the
            offset in each dimension is drawn from the full box length. Default
            is None.

        cluster_size_threshold : int or None, optional
            Soft maximum cluster size (number of chains); if the connected
            component exceeds it the move is rejected. In simulation.py this is
            set so a cluster spanning all chains is not moved. Default is None.

        hardwall : bool, optional
            If True the move is rejected when a translated chain would cross the
            wall or straddle a periodic boundary; otherwise periodic boundary
            conditions are used. Default is False.

        frozen_chains : sequence of int, optional
            chainIDs which are frozen and cannot be moved. If a frozen chain ends
            up in the cluster the move is rejected. Default is ``()``.

        Returns
        -------
        tuple
            ``(MoveEvent, True)`` if the cluster was translated (the move is
            energy-neutral for short-range interactions by construction), or
            ``(False, False)`` if rejected (cluster exceeds the size threshold,
            includes a frozen chain, a hard-sphere clash, a hardwall violation,
            or the translation would merge/resize the cluster), in which case
            the lattice is left unchanged.
        """

        original_chainID  = selected_chain.chainID
        dimensions        = latticeObject.dimensions
        num_dims          = len(dimensions)
        lattice           = latticeObject.grid 
        
        # note that get_all_chains_in_connected_component returns all the chainIDs in the connected
        # component chainID is part of *including* chainID!                                
        try:
            list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                       lattice, 
                                                                                       latticeObject.chains, 
                                                                                       threshold=cluster_size_threshold,
                                                                                       useChains=True, 
                                                                                       hardwall=hardwall)

        # this 'exception' occurs if we're scanning a connected component and discover it's larger than the cluster_threshold.
        # It isn't really an exception - it implements the size threshold in an interrupt style (maximally efficient).
        # NOTE the semantics are strict: the component search checks the size after EVERY BFS wave, so any component
        # strictly larger than the threshold ALWAYS raises - a cluster larger than the threshold is never moved.
        # (Accepted moves preserve cluster size, so this constraint is symmetric and detailed-balance safe.)
        except ClusterSizeThresholdException:
            return (False, False)

        # exclude clusters where one of the chains is in the frozen list
        frozen_chain_set = set(_normalized_frozen_chains(latticeObject, frozen_chains))
        for chainID in list_of_chains_in_CC:
            if chainID in frozen_chain_set:
                return (False, False)

        # these dictionaries hold chainID indexed list of positions associated with a chain in their original
        # and new position
        old_chain_positions = {}
        new_chain_positions = {}

        # determine what the translation operation is gonna be
        offset_vector = []

        # cluster move threshold allows us to define translational movement as occuring in a maximum stepsize in each
        # dimension
        # if not included than we randomly move some distance
    
        if cluster_move_threshold is None:
            for i in range(0, num_dims):
                offset_vector.append(numpy_utils.randneg(random.randint(1, dimensions[i]-1)))
                      
        else:
            for i in range(0, num_dims):
                offset_vector.append(numpy_utils.randneg(random.randint(1, min(dimensions[i]-1, cluster_move_threshold))))
            
        # Delete EVERY cluster chain from the grid before placing a single
        # translated bead. Deleting lazily, chain by chain, made a translated bead
        # that landed on a cluster-mate's not-yet-vacated site a spurious clash;
        # because the chains are processed in a fixed order that rejection was
        # DIRECTIONAL (a shift by +v was refused while -v from the shifted state
        # was accepted: 0 vs 2151 successes in 200,000 trials between two states
        # of identical energy), breaking the proposal symmetry the move relies on.
        for chainID in list_of_chains_in_CC:
            old_chain_positions[chainID] = latticeObject.chains[chainID].get_ordered_positions()
            lattice_utils.delete_chain_by_position(old_chain_positions[chainID], lattice, chainID)

        def _revert(current_chain, current_positions):
            """Undo every placement made so far and restore the whole cluster.

            Parameters
            ----------
            current_chain : int
                chainID of the chain being placed when the rejection happened.

            current_positions : list of list of int
                The translated positions of ``current_chain`` that have already
                been written into the grid (a partial chain), which must be
                cleared before the original cluster is put back.

            Returns
            -------
            None
            """
            lattice_utils.delete_chain_by_position(current_positions, lattice, current_chain)
            for placed_id in new_chain_positions:
                lattice_utils.delete_chain_by_position(new_chain_positions[placed_id], lattice, placed_id)
            for cluster_id in list_of_chains_in_CC:
                lattice_utils.place_chain_by_position(old_chain_positions[cluster_id], lattice, cluster_id, safe=True)

        # now cycle through each chain in the connected commponent
        for chainID in list_of_chains_in_CC:

            # move to its new position
            translated_positions = []
            for position in old_chain_positions[chainID]:

                translated_pos = []
                
                # determine the translated position
                for dim in range(0, num_dims):
                    translated_pos.append(position[dim] + offset_vector[dim] )

                # Under HARDWALL a raw coordinate outside the box means the bead
                # would pass THROUGH the wall - the move must be rejected, never
                # periodically wrapped. (Wrapping was the old behaviour: the
                # per-chain bond-straddle check below cannot see a monomer or a
                # whole-chain wrap, so a box-spanning cluster could be "translated"
                # by permuting its chains through the wall - occupancy unchanged,
                # committed as energy-neutral, while the true hardwall SR/solvation
                # energy changed. Same convention as vmmc_move.)
                out_of_box = False
                if hardwall:
                    for dim in range(0, num_dims):
                        if translated_pos[dim] < 0 or translated_pos[dim] >= dimensions[dim]:
                            out_of_box = True
                            break
                translated_pos = lattice_utils.pbc_convert(translated_pos, dimensions)

                # if the proposed position is already occupied (or, under hardwall,
                # would have crossed the wall) back the f*ck up
                if out_of_box or not lattice_utils.get_gridvalue(translated_pos, lattice) == 0:
                    _revert(chainID, translated_positions)
                    return (False, False)

                # if the position was free update the lattice grid object and add the position to
                # the growing list of new positions for this chain
                lattice_utils.set_gridvalue(translated_pos, chainID, lattice)
                translated_positions.append(translated_pos)            

            # if we get here we succesfully inserted an entire chain into the grid so save 
            # the translated_positions as the chain's positions, BUT FIRST check for hardwall rules and reject 
            # the move IF we're applying a hardwall boundary and the chain breaks that hardwall
            if hardwall:

                if lattice_utils.do_positions_stradle_pbc_boundary(translated_positions):
                    _revert(chainID, translated_positions)
                    return (False, False)

            new_chain_positions[chainID] = translated_positions
            
        # >>>>
        # if we get here we moved the cluster! However we have to determine if the NEW cluster position is also 
        # the same size connected component - if it's larger this move would break detailed balance, 
        # if it's the same this move is fine but it is (by definition) energy neutral so no need to do 
        # short range energy calculations (need long-range ones though!)

        size_of_original_cluster = len(list_of_chains_in_CC)
        
        # we now build a new, bespoke positions dictionary which contains the new positions the moved chains
        # and the original positions of the chains which haven't moved (i.e. just a lattice up-to-date list
        # of chain positions
        chainPositionDict={}
        for chainID in latticeObject.chains:
            if chainID in new_chain_positions:
                chainPositionDict[chainID] = new_chain_positions[chainID]
            else:
                chainPositionDict[chainID] = latticeObject.chains[chainID].get_ordered_positions()
        
        try:
            new_list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                           lattice, 
                                                                                           chainPositionDict, 
                                                                                           threshold=size_of_original_cluster,
                                                                                           useChains=False,
                                                                                           hardwall=hardwall)


            # finally make sure that if we re-selected this chain and got a list of chains it's THE SAME list of chains!
            # This is actually important - same size is too lenient, as you could move across a pbc but keep number of chain
            # fixed - the implementation below is a necessary and sufficient check to ensure that:

            # a) All the new chain IDs were in the list of old IDs
            for new_id in new_list_of_chains_in_CC:
                if new_id not in list_of_chains_in_CC:
                    raise ClusterSizeThresholdException

            # b) All the old chain IDs are in the new of new IDs
            for old_id in list_of_chains_in_CC:
                if old_id not in new_list_of_chains_in_CC:
                    raise ClusterSizeThresholdException
               

        # if we find adding more chains than we had before in the CC then an exception is raised and we know our 
        # cluster move caused cluster merging        
        # ****************************************************************************************************
        except ClusterSizeThresholdException:

            # revert back by deleting the chains we insterted and then re-setting the old chain
            chains_reinserted = list(new_chain_positions.keys())
            for chainIDs_inserted in chains_reinserted:
                lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted)

            for chainIDs_inserted in chains_reinserted:
                lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted, safe=True)

            return (False, False)
        # ****************************************************************************************************

        # if we get here move is a go!
        ME = MoveEvent(original_positions        = old_chain_positions,
                       moved_positions           = new_chain_positions,
                       original_chain_positions  = old_chain_positions,
                       moved_chain_positions     = new_chain_positions,
                       moved_indices             = None,
                       move_type                 = 7)
        
        return (ME, True)                



    #-----------------------------------------------------------------
    #    
    def cluster_rotate(self, selected_chain, latticeObject, cluster_move_threshold=None, cluster_size_threshold=None, hardwall=False, frozen_chains=()):
        """
        The cluster_rotate move allows a connected components (cluster) to be 
        rotated in rigid body space around the lattice. Right now rotation occurs only
        over the cardinal directions (0/90/180/270) because anything other than this is 
        hard on a lattice...

        The cost of this move becomes massive as clusters get big - the cluster_threshold 
        variable facilitates a soft threshold on the max cluster size. 

        The cluster_size_threshold is default to None, but in the simulations.py file we
        set this to be such that in the case that a cluster contains ALL the chains it is not
        rotated, otherwise it can be.
        
        The move is rejected if there's a hard-sphere clash, *or* we change the cluster
        size after the move (i.e. incorporate new residues in). If not rejected  we pass
        back the relevant MoveEvent object. Note that like all move functions this
        updates the lattice to contain the chain in the new position.

        A draw that maps the cluster exactly onto itself is rejected as a null
        move rather than returned as a successful proposal. This matters most for
        an isolated single-bead cluster, which is invariant under every rotation
        (9 of 9 draws in 3D, 3 of 3 in 2D): in a box with free monomers those
        identity proposals used to be accepted with zero energy change and
        counted, roughly doubling the acceptance ratio reported for this move.
        The sampled ensemble is unaffected by the change.

        MoveType code: 8

        Parameters
        ----------
        selected_chain : Chain
            A chain belonging to the cluster to be rotated; its connected
            component defines the cluster.

        latticeObject : Lattice
            The lattice object containing the chains. Its grid is mutated in
            place if the move succeeds.

        cluster_move_threshold : int or None, optional
            Unused by the rotation itself but accepted for signature symmetry
            with cluster_translate. Default is None.

        cluster_size_threshold : int or None, optional
            Soft maximum cluster size (number of chains); if the connected
            component exceeds it the move is rejected. In simulation.py this is
            set so a cluster spanning all chains is not rotated. Default is None.

        hardwall : bool, optional
            If True the move is rejected when a rotated chain straddles a
            periodic boundary. Default is False.

        frozen_chains : sequence of int, optional
            chainIDs that are frozen; if any frozen chain is in the cluster the
            move is rejected. Default is ``()``.

        Returns
        -------
        tuple
            ``(MoveEvent, True)`` if the cluster was rotated (energy-neutral for
            short-range interactions by construction), or ``(False, False)`` if
            rejected (size threshold exceeded, frozen chain present, a draw that
            maps the cluster exactly onto itself, hard-sphere clash, hardwall
            violation, or the rotation would merge/resize the cluster), in which
            case the lattice is left unchanged.
        """

        original_chainID        = selected_chain.chainID
        dimensions              = latticeObject.dimensions
        num_dims                = len(dimensions)
        lattice                 = latticeObject.grid 

        old_chain_positions            = {}
        new_chain_positions_OC         = {}
        new_chain_positions            = {}
        
        # note that get_all_chains_in_connected_component returns all the chainIDs in the connected
        # component chainID is part of *including* chainID!                
        try:
            list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                       lattice, 
                                                                                       latticeObject.chains, 
                                                                                       threshold=cluster_size_threshold,
                                                                                       useChains=True,
                                                                                       hardwall=hardwall)

        # this 'exception' occurs if we're scanning a connected component and discover it's larger than the cluster_threshold.
        # It isn't really an exception - it implements the size threshold in an interrupt style (maximally efficient).
        # NOTE the semantics are strict: the component search checks the size after EVERY BFS wave, so any component
        # strictly larger than the threshold ALWAYS raises - a cluster larger than the threshold is never moved.
        # (Accepted moves preserve cluster size, so this constraint is symmetric and detailed-balance safe.)

            
        except ClusterSizeThresholdException:
            return (False, False)

        # exclude clusters where one of the chains is in the cluster list
        frozen_chain_set = set(_normalized_frozen_chains(latticeObject, frozen_chains))
        for chainID in list_of_chains_in_CC:
            if chainID in frozen_chain_set:
                return (False, False)
        
        # these dictionaries hold chainID indexed list of positions associated with a chain in their original
        # and new position - delete these chains from the lattice!
        # NOTE: iterate in sorted-chainID order so the concatenated bead order (and
        # therefore the first-index argmin pivot tie-break) is canonical - the
        # connected-component search returns a SET, whose iteration order is only
        # accidentally stable, and the forward and reverse moves must concatenate
        # identically for the pivot selection to invert.
        all_cluster_positions =[]
        _chain_lengths = []
        list_of_chains_in_CC = sorted(list_of_chains_in_CC)
        for chainID in list_of_chains_in_CC:

            _cp = latticeObject.chains[chainID].get_ordered_positions()
            all_cluster_positions.extend(_cp)
            old_chain_positions[chainID] = _cp
            _chain_lengths.append(len(_cp))

            lattice_utils.delete_chain_by_position(old_chain_positions[chainID], lattice, chainID)


        # so now ALL the chains in the cluster have been deleted from the lattice.
        #
        # We rotate in SINGLE-IMAGE space and about a fixed cluster BEAD, not the
        # rounded circular-mean COM of the raw (possibly PBC-straddling) positions.
        # Two things were wrong before:
        #   1. Rotating the RAW positions: if the cluster straddles a periodic
        #      boundary the raw coordinates are discontiguous, so the rotation is not
        #      a rigid-body move of the physical cluster at all.
        #   2. Rotating about a rounded COM: the rounded COM of the rotated cluster is
        #      generally a different lattice point, so the inverse rotation does not
        #      return the original state -> detailed balance is violated, and because
        #      energy-neutral cluster rotations are effectively always accepted nothing
        #      compensates the resulting bias.
        # Single-imaging first makes it a genuine rigid body, and rotating about the
        # bead nearest the (single-image) centroid - a rotation+translation-invariant
        # choice for a rigid body - makes the move exactly reversible: the reverse
        # move re-single-images, picks the SAME physical bead, and inverts the
        # rotation, with the periodic re-wrap cancelling exactly.
        if hardwall:
            # Hardwall coordinates are plain Cartesian: there is no periodic image
            # to reconstruct and no winding to guard against. Running the periodic
            # gather and the winding guard here refused the rotation of any cluster
            # spanning a box axis - a legal configuration - and, because the
            # extents permute under a cardinal rotation in a non-cubic box, refused
            # it in ONE direction only (a 7-mer along x in a 7x10 hardwall box could
            # rotate to y but never back), and emitted a spurious "percolates the
            # periodic box" warning in a run with no periodic box.
            si_all = np.asarray(all_cluster_positions, dtype=np.int64)
        else:
            # The gather warns when the cluster winds around the box, and that
            # warning is written for the ANALYSIS callers, where a shape computed
            # from a winding cluster is meaningless. Here nothing is computed from
            # it: a winding cluster is simply rejected by the guard just below, as
            # documented. Emitting the analysis warning from inside a Monte Carlo
            # move would be misleading, and in a condensed system every rotation
            # drawn on the percolating cluster would emit it, so it is silenced for
            # this call only.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", message="single-image gather: cluster percolates")
                si_all = np.asarray(
                    cluster_utils.convert_positions_to_single_image_snakesearch(
                        all_cluster_positions, dimensions), dtype=np.int64)

            # A cluster that WINDS around the box (is connected to its own periodic
            # image) has a single-image extent >= the box length on some axis. A
            # cardinal rotation of such a cluster is NOT a rigid motion of the
            # periodic system: the winding closure vector maps onto an axis with a
            # different period, so intra-cluster minimum-image LR/SLR (and in
            # principle SR) relations change while the move's dE assumes they are
            # invariant - silently corrupting the tracked energy. Reject outright
            # (rejection is symmetric: the winding property is preserved by the
            # move, so detailed balance is unaffected).
            for _d in range(num_dims):
                if int(si_all[:, _d].max() - si_all[:, _d].min()) + 1 >= dimensions[_d]:
                    for _cid in list_of_chains_in_CC:
                        lattice_utils.place_chain_by_position(old_chain_positions[_cid], lattice, _cid, safe=True)
                    return (False, False)

        # exact-integer pivot selection (see chain_rotate: a float centroid breaks
        # exact distance ties by rounding noise that is not rotation-invariant,
        # making tied configurations non-invertible - a detailed-balance violation)
        _d2 = ((len(si_all) * si_all - si_all.sum(axis=0)) ** 2).sum(axis=1)
        _cl_pivot = int(np.argmin(_d2))
        _si_pivot = si_all[_cl_pivot]
        # the pivot bead's ORIGINAL (in-box, on-lattice) position - the anchor the
        # rotated displacement vectors are re-hung on (see chain_rotate for why this
        # makes the move exactly reversible).
        _raw_pivot = list(all_cluster_positions[_cl_pivot])

        # per-chain SINGLE-IMAGE displacement vectors relative to the pivot bead
        si_chain_positions = {}
        _cursor = 0
        for chainID, _L in zip(list_of_chains_in_CC, _chain_lengths):
            si_chain_positions[chainID] = [
                [int(si_all[_cursor + k][d] - _si_pivot[d]) for d in range(num_dims)]
                for k in range(_L)]
            _cursor += _L

        # chains whose rotated RAW coordinates leave the box (hardwall only; see
        # the rejection in the insertion loop below)
        wall_crossing_chains = set()
        
        ## ----------------------------------------------------------------------------------------------------
        ## 2D CASE FIRST
        ##        
        if num_dims == 2:

            rotationFactor = [90,180,270][random.randint(0,2)]

            # now cycle through each chain in the connected commponent moving it such that it's centered on the
            # origin (OC = origin centered)
            for chainID in list_of_chains_in_CC:
                                             
                # si_chain_positions already holds displacement vectors from the pivot
                new_chain_positions_OC[chainID] = si_chain_positions[chainID]

                # rotate 2D displacement vectors by the rotation operation defined
                new_chain_positions_OC_rotated = lattice_utils.rotate_positions_2D(new_chain_positions_OC[chainID], rotationFactor)

                # re-anchor on the pivot bead's original in-box position, then re-wrap
                new_chain_positions[chainID] = []
                for position in new_chain_positions_OC_rotated:
                    raw = [position[0] + _raw_pivot[0], position[1] + _raw_pivot[1]]
                    if hardwall and _outside_box(raw, dimensions):
                        wall_crossing_chains.add(chainID)
                    new_chain_positions[chainID].append(lattice_utils.pbc_convert(raw, dimensions))



        ## ----------------------------------------------------------------------------------------------------
        ## 3D CASE SECOND
        ##
        else:
            rotationFactor = [90,180,270][random.randint(0,2)]
            rotationDim    = ['x','y','z'][random.randint(0,2)]

            # now cycle through each chain in the connected commponent moving it such that it's centered on the
            # origin (OC = origin centered)
            for chainID in list_of_chains_in_CC:
                                                
                # si_chain_positions already holds displacement vectors from the pivot
                new_chain_positions_OC[chainID] = si_chain_positions[chainID]

                # rotate 3D displacement vectors by the rotation operation defined
                new_chain_positions_OC_rotated = lattice_utils.rotate_positions_3D(new_chain_positions_OC[chainID], rotationDim, rotationFactor)

                # re-anchor on the pivot bead's original in-box position, then re-wrap
                new_chain_positions[chainID] = []
                for position in new_chain_positions_OC_rotated:
                    raw = [position[0] + _raw_pivot[0], position[1] + _raw_pivot[1], position[2] + _raw_pivot[2]]
                    if hardwall and _outside_box(raw, dimensions):
                        wall_crossing_chains.add(chainID)
                    new_chain_positions[chainID].append(lattice_utils.pbc_convert(raw, dimensions))


        # Guaranteed-null draws are REJECTED, not accepted (see chain_rotate for the
        # full argument). A cluster that is a symmetry axis of itself maps onto
        # itself: an ISOLATED single-bead cluster is invariant under every draw
        # (9 of 9 in 3D, 3 of 3 in 2D), and an axis-aligned rigid cluster under the
        # three rotations about that axis. Because a monomer-rich box makes cluster
        # rotation land on an isolated monomer very often, these used to dominate
        # the move's ACCEPTANCE.dat column - a measured 2.2x over-report of the
        # acceptance ratio in a box of free monomers, every one of those "accepted"
        # rotations having changed nothing.
        #
        # Rejecting is ensemble-neutral: detailed balance places no constraint on
        # P(x->x), and accepting or rejecting an identity proposal leaves the same
        # state behind, so the trajectory is bit-identical.
        if all(_is_identity_proposal(new_chain_positions[chainID],
                                     old_chain_positions[chainID], num_dims)
               for chainID in list_of_chains_in_CC):
            # nothing has been re-inserted yet, so putting the cluster back is the whole revert
            for _cid in list_of_chains_in_CC:
                lattice_utils.place_chain_by_position(old_chain_positions[_cid], lattice, _cid, safe=True)
            return (False, False)


        ## ----------------------------------------------------------------------------------------------------
        # having built a new list of rotated positions for each chain let's see if any of them clash. Note the inserted_chain
        # keeps track of what's going on so if we find a clash we only have to cycle over a small number of filled
        # positions to delete the part of the chain we inserted

        # for each position in each chain
        chains_reinserted = []
        for chainID in new_chain_positions:

            # Under HARDWALL a raw rotated coordinate outside the box means the
            # bead would pass THROUGH the wall: reject, never wrap. The per-chain
            # bond-straddle check further down cannot see this for a single-bead
            # chain, or for a chain that leaves the box entirely - both were
            # periodically wrapped to the opposite face, committed as SR-energy
            # neutral while the true energy changed, and the reverse rotation was
            # then rejected by the winding check (irreversible move). Same
            # convention as cluster_translate and vmmc_move.
            if chainID in wall_crossing_chains:
                for chainIDs_rotated in chains_reinserted:
                    lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_rotated], lattice, chainIDs_rotated)
                for chainIDs_org in old_chain_positions:
                    lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_org], lattice, chainIDs_org, safe=True)
                return (False, False)

            rotated_positions = []
            for position in new_chain_positions[chainID]:

                # if the position we're rotating into is CURRENTLY occupied
                if not lattice_utils.get_gridvalue(position, lattice) == 0:

                    IO_utils.status_message("Rejection because of clash", 'info', allow_suppress=True)
                    
                    # Delete the positions we insterted so far in the *current* chain and then
                    # delete all the other chains which were fully rotated
                    lattice_utils.delete_chain_by_position(rotated_positions, lattice, chainID)
                    for chainIDs_rotated in chains_reinserted:
                        lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_rotated], lattice, chainIDs_rotated)

                    # now reinsert ALL the chains back....
                    for chainIDs_org in old_chain_positions:
                         lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_org], lattice, chainIDs_org, safe=True)

                    # reject the move!
                    return (False, False)

                # else we're OK
                rotated_positions.append(position)
                lattice_utils.set_gridvalue(position, chainID, lattice)
                                    

            if hardwall:
                if lattice_utils.do_positions_stradle_pbc_boundary(rotated_positions):
                    
                    # this is exactly the same protocol as we use to reject the move in the case of the clash above, just not annotated
                    # as heavily...
                    lattice_utils.delete_chain_by_position(rotated_positions, lattice, chainID)
                    for chainIDs_rotated in chains_reinserted:
                        lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_rotated], lattice, chainIDs_rotated)

                    # now reinsert ALL the chains back....
                    for chainIDs_org in old_chain_positions:
                         lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_org], lattice, chainIDs_org, safe=True)

                    # reject the move!
                    return (False, False)

            # if we get here inserted the whole chain, so add it to the list of [succesfully] re-inserted chains
            chains_reinserted.append(chainID)

        # if we get here we moved the cluster! However we have to determine if the NEW cluster position is also 
        # the same size connected component - if it's larger this move would break detailed balance, 
        # if it's the same this move is fine but it is (by definition) energy neutral so no need to do 
        # energy calculations

        size_of_original_cluster = len(list_of_chains_in_CC)
        
        # we now build a new, bespoke positions dictionary which contains the new positions of the moved chains
        # and the original positions of the chains which haven't moved (i.e. just a lattice up-to-date list
        # of chain positions
        chainPositionDict={}
        for chainID in latticeObject.chains:
            if chainID in new_chain_positions:
                chainPositionDict[chainID] = new_chain_positions[chainID]
            else:
                chainPositionDict[chainID] = latticeObject.chains[chainID].get_ordered_positions()
        
        try:
            new_list_of_chains_in_CC = lattice_utils.get_all_chains_in_connected_component(original_chainID, 
                                                                                           lattice, 
                                                                                           chainPositionDict, 
                                                                                           threshold=size_of_original_cluster,
                                                                                           useChains=False,
                                                                                           hardwall=hardwall)

            # finally make sure the new cluster isn't SMALLER (could happen in hardwall mode) 
            for new_id in new_list_of_chains_in_CC:
                if new_id not in list_of_chains_in_CC:
                    raise ClusterSizeThresholdException

            for old_id in list_of_chains_in_CC:
                if old_id not in new_list_of_chains_in_CC:
                    raise ClusterSizeThresholdException


        # if we find adding more chains than we had before in the CC then an exception is raised and we know our cluster move caused cluster
        # merging
        # ****************************************************************************************************
        except ClusterSizeThresholdException:

            IO_utils.status_message("Cluster resize rejection", 'info', allow_suppress=True)

            # revert back by deleting the chains we insterted and then re-setting the old chain
            chains_reinserted = list(new_chain_positions.keys())
            for chainIDs_inserted in chains_reinserted:
                lattice_utils.delete_chain_by_position(new_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted)

            for chainIDs_inserted in chains_reinserted:
                lattice_utils.place_chain_by_position(old_chain_positions[chainIDs_inserted], lattice, chainIDs_inserted, safe=True)

            return (False, False)
        # ****************************************************************************************************

        # if we get here move is a go!
        ME = MoveEvent(original_positions        = old_chain_positions,
                       moved_positions           = new_chain_positions,
                       original_chain_positions  = old_chain_positions,
                       moved_chain_positions     = new_chain_positions,
                       moved_indices             = None,
                       move_type                 = 8)
        
        return (ME, True)                



    #-----------------------------------------------------------------
    #    
    def Chain_based_TSMMC(self, chainID, latticeObject, current_energy, hamiltonianObject, CTSMMC, hardwall=False):
        """
        The chain-based Temperature Sweep Metropolis Monte Carlo move heats a single chain along a
        temperature schedule and cools it back down again, relaxing it at every rung on the way.

        In terms of big picture - this move involves creating an alternative Monte Carlo chain. The
        schedule is palindromic: it ramps up from the simulation temperature to TSMMC_JUMP_TEMP, holds
        there, and ramps back down, with the same number of crankshaft sub-moves (each a Metropolis
        move at that rung's temperature) at every rung and none at the simulation temperature itself.

        The excursion as a whole is accepted or rejected with the tempered-transitions criterion
        (Neal 1996): the work accumulated as sum((beta_before - beta_after) * E) over every temperature
        change of the schedule, including the step off and back onto the simulation temperature, is
        passed to TSMMC.accept_tempered_transition. It is NOT a Metropolis test on the end-point
        energies - that would break detailed balance - see docs/moves/tsmmc.rst.

        HOWEVER, we evaluate this move here and then don't through the standard single chain energy evaluation
        because throughout the actual move we keep track of the system energy so don't have to re-evaluate 
        after. 

        [1] Mittal, A., Lyle, N., Harmon, T.S., and Pappu, R.V. (2014). Hamiltonian Switch Metropolis Monte Carlo 
            Simulations for Improved Conformational Sampling of Intrinsically Disordered Regions Tethered to Ordered 
            Domains of Proteins. J. Chem. Theory Comput. 10, 3550-3562.

        [2] Gelb, L.D. (2003). Monte Carlo simulations using sampling from an approximate potential. J. Chem. Phys. 118, 7747-7750.

        MoveType code: 9

        Parameters
        ----------
        chainID : int
            The ID of the single chain to be perturbed by the temperature
            excursion.

        latticeObject : Lattice
            The full lattice object being simulated; its grids and chain
            positions are mutated in place (and reverted if the excursion is
            rejected).

        current_energy : int or float
            The current system energy at the start of the excursion.

        hamiltonianObject : Hamiltonian
            Object providing the interaction tables and angle lookup passed to
            the Cython kernel for energy evaluation.

        CTSMMC : TSMMC
            The TSMMC coordinator providing the inverse-temperature schedule,
            steps-per-temperature multiplier, and the tempered-transitions
            acceptance test.

        hardwall : bool, optional
            If True a hard-wall boundary is used; otherwise periodic boundary
            conditions are used. Default is False.

        Returns
        -------
        tuple
            ``(latticeObject, current_energy, total_moves, accepted)`` where
            ``current_energy`` is the new energy (or the original energy if
            rejected), ``total_moves`` is the number of sub-moves proposed during
            the excursion, and ``accepted`` is True if the excursion was
            accepted.
        """

        idx_to_bead = crankshaft_list_functions.update_idx_to_bead_single_chain(latticeObject, chainID)
        chain_length = len(idx_to_bead)
        
        # this is a copy because its a list (if this was a numpy array would be by reference and would
        # not be a copy) 
        original_chain_positions = copy.deepcopy(latticeObject.chains[chainID].get_ordered_positions())

        # save old energy
        old_energy = current_energy
        num_dims = len(latticeObject.dimensions)
        num_temps = len(CTSMMC.inv_temperature_schedule)
        
        steps_per_temperature = chain_length * CTSMMC.steps_per_quench_multiplier

        total_moves = steps_per_temperature * num_temps
            
        # these are passed by reference, but we set to them variables so we can iteratively pass them
        # at different temperatures
        #tmp_grid            = latticeObject.grid
        #tmp_type_grid       = latticeObject.type_grid
        new_energy = current_energy

        # set new energy to current energy - this will be updated sequentially as we proceed


        # set hardwall flag
        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0

        # Tempered-transitions / NCMC bookkeeping. The excursion drives the
        # temperature off the target value, through the schedule, and back. To
        # preserve detailed balance we must accumulate the work
        # (beta_before - beta_after) * U(x) at EVERY temperature change, using
        # the energy at the instant of the change. We start at the target
        # temperature with energy `current_energy`.
        log_work = 0.0
        prev_inv = CTSMMC.inv_target_temperature

        for temp_idx in range(0, num_temps):

            # set previous and current inverse temperatures
            inv_temp = CTSMMC.inv_temperature_schedule[temp_idx]

            # work contribution of changing temperature prev_inv -> inv_temp,
            # evaluated at the current configuration energy (before propagating)
            log_work = log_work + (prev_inv - inv_temp) * new_energy

            local_seed = random.randint(1, sys.maxsize - 1)

            bead_selector = np.random.randint(0, chain_length, steps_per_temperature)


            ##
            ## Both functions alter alter the grids on the back end and do not explicity
            ## reassign these as they're passed by reference as memoryviews (direct access to
            ## the memory)
            ##

            if num_dims == 2:
                (new_energy, accepted_moves)= mega_crank_fast.mega_crank_2D(latticeObject.grid,
                                                                          latticeObject.type_grid,
                                                                          idx_to_bead,
                                                                          hamiltonianObject.residue_interaction_table,
                                                                          hamiltonianObject.LR_residue_interaction_table,
                                                                          hamiltonianObject.SLR_residue_interaction_table,
                                                                          hamiltonianObject.angle_lookup,
                                                                          new_energy,
                                                                          inv_temp,
                                                                          steps_per_temperature,
                                                                          bead_selector,
                                                                          local_seed,
                                                                          hardwall_int)

            else:


                (new_energy, accepted_moves) = mega_crank_fast.mega_crank(latticeObject.grid,
                                                                     latticeObject.type_grid,
                                                                     idx_to_bead,
                                                                     hamiltonianObject.residue_interaction_table,
                                                                     hamiltonianObject.LR_residue_interaction_table,
                                                                     hamiltonianObject.SLR_residue_interaction_table,
                                                                     hamiltonianObject.angle_lookup,
                                                                     new_energy,
                                                                     inv_temp,
                                                                     steps_per_temperature,
                                                                     bead_selector,
                                                                     local_seed,
                                                                     hardwall_int)

            prev_inv = inv_temp

        # final temperature change: schedule[-1] -> target temperature, at the
        # final configuration energy. This completes the path work sum.
        log_work = log_work + (prev_inv - CTSMMC.inv_target_temperature) * new_energy

        # if move is accepted update the grids, the energy, and the chain positions
        if CTSMMC.accept_tempered_transition(log_work):

            # udpate the chain positions
            current_energy = new_energy 

            # the beauty is this works for the 2D and 3D case
            latticeObject.chains[chainID].positions = idx_to_bead[:,5:].tolist()
            
            return (latticeObject, current_energy, total_moves, True)
            
        # reject the whole move
        else:
            
            # construct a new list of the chain's new positions based on the tmp_chain_positions matrix
            deletable_positions = idx_to_bead[:,5:].tolist()
            
            # revert the lattice to it's pre-move state 
            lattice_utils.delete_chain_by_position(deletable_positions, latticeObject.grid, chainID)                
            lattice_utils.place_chain_by_position(original_chain_positions, latticeObject.grid, chainID, safe=True)
            
            # set chain positions in the chain-list positions                        
            latticeObject.chains[chainID].set_ordered_positions(original_chain_positions)
                
            # update the type_grid variable BACK 
            latticeObject.update_type_grid(chainID, deletable_positions, original_chain_positions, list(range(0,len(original_chain_positions))), safe=True)
            
            # return everything!
            return (latticeObject, old_energy, total_moves, False)
                                    

                
    #-----------------------------------------------------------------
    #    
    def multichain_based_TSMMC(self, original_chainID, latticeObject, current_energy, hamiltonianObject, CTSMMC, hardwall=False, frozen_chains=()):
        """
        Same idea as Chain_based_TSMMC except here we randomly select some number of chains (currently this is 
        defined by the max_number_selectable function, which is set at 25% of the total number of chains on the
        lattice.

        Then, we sequentially raise and then lower the temperature, and at each different temperature cycle through
        the chains and update their positions. At the end the full move is accepted or rejected. See the 
        chain_based_TSMMC write up for more details on what's actually going on in terms of the TSMMC-ness.

        MoveType code: 10

        Parameters
        ----------
        original_chainID : int
            A chain ID associated with the move (kept for signature symmetry).
            The actual chains perturbed are chosen randomly from the
            non-frozen chains (between 1 and ~25% of them).

        latticeObject : Lattice
            The full lattice object being simulated; its grids and chain
            positions are mutated in place (and reverted if rejected).

        current_energy : int or float
            The current system energy at the start of the excursion.

        hamiltonianObject : Hamiltonian
            Object providing the interaction tables and angle lookup passed to
            the Cython kernel for energy evaluation.

        CTSMMC : TSMMC
            The TSMMC coordinator providing the inverse-temperature schedule,
            steps-per-temperature multiplier, and the tempered-transitions
            acceptance test.

        hardwall : bool, optional
            If True a hard-wall boundary is used; otherwise periodic boundary
            conditions are used. Default is False.

        frozen_chains : sequence of int, optional
            chainIDs excluded from selection. Default is ``()``.

        Returns
        -------
        tuple
            ``(latticeObject, current_energy, total_moves, accepted)`` where
            ``current_energy`` is the new energy (or the original energy if
            rejected), ``total_moves`` the number of sub-moves proposed, and
            ``accepted`` is True if the excursion was accepted. If all chains
            are frozen, ``(latticeObject, current_energy, 0, False)`` is
            returned unchanged.
        """
                    
        dimensions      = latticeObject.dimensions
        num_dims        = len(dimensions)
        num_temps       = len(CTSMMC.inv_temperature_schedule)
        old_energy      = current_energy

        
        ## in the current implementation we randomly select between 1 and 25% of the chains in the system
        # First figure out what 25% of the number of chains is
        tmp_all_chains      = list(latticeObject.chains.keys())

        # exclude frozen chains
        frozen_chain_set = set(_normalized_frozen_chains(latticeObject, frozen_chains))
        if len(frozen_chain_set) > 0:
            all_chains = []
            for c in tmp_all_chains:
                if c not in frozen_chain_set:
                    all_chains.append(c)
        else:
            all_chains = tmp_all_chains

        # get number of chains we might select from
        num_chains          = len(all_chains)

        # If all chains are frozen (or no chains exist), there is nothing to do.
        if num_chains == 0:
            return (latticeObject, current_energy, 0, False)

        # this works with 1 through n chains and give sensible values
        max_number_selectable = int(np.floor(0.25*num_chains) + 1)
        max_number_selectable = min(num_chains, max_number_selectable)
        number_selectable     = random.randint(1,max_number_selectable)
        list_of_chains = np.random.choice(all_chains, number_selectable,replace=False) # DO NOT REPLACE!

        # list_of_chains is now a list of chain IDs that we're going to perturb

        # this dictionary allows us to map the chainID to the old positions so we can rever if needed
        all_original_positions = {}        

        # for each chain, save the original positions via a deepcopy operation
        for chainID in list_of_chains:            
            # save the original positions in case we have to revert
            all_original_positions[chainID] = copy.deepcopy(latticeObject.chains[chainID].get_ordered_positions())

        # construct a specific idx_to_bead matrix that reflects the beads taken from this list of chains in the order
        # they appear in the list_of_chains
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead_multiple_chains(latticeObject, list_of_chains)

        # total number of beads
        num_beads = idx_to_bead.shape[0]

        # calculate the number of steps per temperature
        steps_per_temperature = len(idx_to_bead)*CTSMMC.steps_per_quench_multiplier
                
        # total proposed moves (keep track for peformance analysis post-factor)
        total_moves = steps_per_temperature * num_temps

        new_energy          = current_energy
        #tmp_grid            = latticeObject.grid
        #tmp_type_grid       = latticeObject.type_grid

        # set hardwall flag
        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0

        # Tempered-transitions / NCMC work accumulator (see Chain_based_TSMMC and
        # TSMMC.accept_tempered_transition). Detailed balance requires summing
        # (beta_before - beta_after) * U(x) over every temperature change.
        log_work = 0.0
        prev_inv = CTSMMC.inv_target_temperature

        for temp_idx in range(0, num_temps):

            # set previous and current inverse temperatures
            inv_temp = CTSMMC.inv_temperature_schedule[temp_idx]

            # work for the temperature change prev_inv -> inv_temp at the current energy
            log_work = log_work + (prev_inv - inv_temp) * new_energy

            local_seed = random.randint(1, sys.maxsize - 1)


            bead_selector = np.random.randint(0, num_beads, steps_per_temperature)

            ##
            ## Both functions alter alter the grids on the back end and do not explicity
            ## reassign these as they're passed by reference as memoryviews (direct access to
            ## the memory)
            ## 

            if num_dims == 2:
                
                (new_energy, accepted_moves)= mega_crank_fast.mega_crank_2D(latticeObject.grid,
                                                                          latticeObject.type_grid,
                                                                          idx_to_bead,
                                                                          hamiltonianObject.residue_interaction_table,
                                                                          hamiltonianObject.LR_residue_interaction_table,
                                                                          hamiltonianObject.SLR_residue_interaction_table,
                                                                          hamiltonianObject.angle_lookup,
                                                                          new_energy,
                                                                          inv_temp,
                                                                          steps_per_temperature,
                                                                          bead_selector,
                                                                          local_seed,
                                                                          hardwall_int)
                
            else:

                (new_energy, accepted_moves) = mega_crank_fast.mega_crank(latticeObject.grid,
                                                                     latticeObject.type_grid,
                                                                     idx_to_bead,
                                                                     hamiltonianObject.residue_interaction_table,
                                                                     hamiltonianObject.LR_residue_interaction_table,
                                                                     hamiltonianObject.SLR_residue_interaction_table,
                                                                     hamiltonianObject.angle_lookup,
                                                                     new_energy,
                                                                     inv_temp,
                                                                     steps_per_temperature,
                                                                     bead_selector,
                                                                     local_seed,
                                                                     hardwall_int)

            prev_inv = inv_temp

        # final temperature change back to the target temperature, completing the path work sum
        log_work = log_work + (prev_inv - CTSMMC.inv_target_temperature) * new_energy

        # if move was accepted
        if CTSMMC.accept_tempered_transition(log_work):
            current_energy = new_energy


            # cycle over each chain, and for each chain update the positions by exracting the updated
            # positions from the tmp_chain_positions matrix. The tmp_chain_poistions matrix contains ONLY
            # bead positions that were moved in a specific order that corresponds to the the order of beads
            # from the chains in the list_of_chains, so this ensures we update position correctly
            idx=0
            for chainID in list_of_chains:
                
                # chain length
                chain_len = len(latticeObject.chains[chainID].positions)

                latticeObject.chains[chainID].positions = idx_to_bead[idx:idx+chain_len,5:].tolist()
                idx = idx + chain_len

            IO_utils.status_message("Multichain re-arrangement accepted [dE = %i]  (number of chains: %i)" %(new_energy - old_energy, len(list_of_chains)), allow_suppress=True)
            return (latticeObject, current_energy, total_moves, True)
                
        else:
            all_new_positions={}

            # same logic as was used to update the chain positions (see text in the move-success branch
            # for an explanation)
            idx=0
            for chainID in list_of_chains:
                
                # chain length
                chain_len = len(latticeObject.chains[chainID].positions)

                # UP-2023
                all_new_positions[chainID] = idx_to_bead[idx:idx+chain_len,5:].tolist()
                idx = idx + chain_len
                
            ## Having obtained the positions of all the moved beads...

            # delete all the chains off the lattice
            for chainID in list_of_chains:
                
                # delete from main grind
                lattice_utils.delete_chain_by_position(all_new_positions[chainID], latticeObject.grid, chainID)

                # delete from type grid
                latticeObject.delete_chain_from_type_grid(chainID, all_new_positions[chainID], list(range(0,len(all_new_positions[chainID]))), safe=True)
                    
            # re-insert the chains back into their original position
            for chainID in list_of_chains:
                original_chain_positions = all_original_positions[chainID]

                # insert into main grid
                lattice_utils.place_chain_by_position(original_chain_positions, latticeObject.grid, chainID, safe=True)

                # insert into type grid
                latticeObject.insert_chain_into_type_grid(chainID, original_chain_positions, list(range(0,len(original_chain_positions))), safe=True)
                                    
                # set chain positions in the chain-list positions
                latticeObject.chains[chainID].set_ordered_positions(original_chain_positions)

                # reset the energy

            current_energy = old_energy

            return (latticeObject, current_energy, total_moves, False)
           


    # Code 11 is the pull megamove, dispatched in simulation.py via
    # MoveObject.system_pull (above) - no per-chain method here.

    # System_based_TSMMC
    # CODE 12
    #
    # Not actually a move implemented here, but implemented in the simulation object with functionality in the TSMMC object
    # too. This is here mainly to ensure that the next move added uses CODE 13 (OOO unlucky!!!!! Sucks to be you!)

    #-----------------------------------------------------------------
    #    
    def single_chain_shake(self, chainID, latticeObject, current_energy, acceptanceObject, hamiltonianObject, number_of_steps, mode, hardwall):
        """
        Perform a single-chain crankshaft shake (many local perturbations of one chain).

        Like system_shake, but restricted to a single chain: a large number of
        local single-bead perturbations are performed on the chain identified by
        ``chainID`` via the optimized Cython crankshaft kernel, with individual
        accept/reject decisions happening per-sub-move inside the kernel. The
        chain's positions are written back from the idx_to_bead matrix once the
        kernel returns.

        Parameters
        ----------
        chainID : int
            The ID of the chain to shake.

        latticeObject : Lattice
            The full lattice object upon which the simulation is being
            performed. Its grids and the chain's positions are mutated in place.

        current_energy : int or float
            The current system energy value (before this megamove).

        acceptanceObject : AcceptanceCalculator
            Object providing the inverse temperature used by the kernel.

        hamiltonianObject : Hamiltonian
            Self-contained object providing the interaction tables and angle
            lookup passed to the external (Cython) kernel for energy evaluation.

        number_of_steps : int
            Number of Monte Carlo sub-moves to perform on the chain.

        mode : str
            Mode used for determining the final number of steps. Currently
            obsolete but kept in case bead selection is changed in the future.

        hardwall : bool
            If True a hard-wall (impenetrable solvent) boundary is used;
            otherwise periodic boundary conditions are used.

        Returns
        -------
        tuple
            ``(latticeObject, current_energy, total_proposed, total_accepted)``
            where ``current_energy`` is the new system energy, ``total_proposed``
            the number of sub-moves attempted, and ``total_accepted`` the number
            accepted.
        """
        
        # get number of dimenisons and set various initial values
        num_dims = len(latticeObject.dimensions)
        idx_to_bead = crankshaft_list_functions.update_idx_to_bead_single_chain(latticeObject, chainID)
        chain_length = len(idx_to_bead)
           
        # set some initial values, the hardwall flag, and set the randoms seed 
        total_accepted = 0
        total_proposed = 0 

        if hardwall:
            hardwall_int = 1
        else:
            hardwall_int = 0
            
        local_seed = random.randint(1, sys.maxsize - 1)

        bead_selector = np.random.randint(0, chain_length, number_of_steps)

        ##
        ## Both functions alter alter the grids on the back end and do not explicity
        ## reassign these as they're passed by reference as memoryviews (direct access to
        ## the memory)
        ## 

        # 2D
        if num_dims == 2:
            (new_energy, accepted_moves) = mega_crank_fast.mega_crank_2D(latticeObject.grid, 
                                                                       latticeObject.type_grid, 
                                                                       idx_to_bead,
                                                                       hamiltonianObject.residue_interaction_table,
                                                                       hamiltonianObject.LR_residue_interaction_table,
                                                                       hamiltonianObject.SLR_residue_interaction_table,
                                                                       hamiltonianObject.angle_lookup,
                                                                       current_energy,
                                                                       acceptanceObject.invtemp,
                                                                       number_of_steps,
                                                                       bead_selector,
                                                                       local_seed,
                                                                       hardwall_int)
                
        else:
            (new_energy, accepted_moves) = mega_crank_fast.mega_crank(latticeObject.grid, 
                                                                 latticeObject.type_grid, 
                                                                 idx_to_bead,
                                                                 hamiltonianObject.residue_interaction_table,
                                                                 hamiltonianObject.LR_residue_interaction_table,
                                                                 hamiltonianObject.SLR_residue_interaction_table,
                                                                 hamiltonianObject.angle_lookup,
                                                                 current_energy,
                                                                 acceptanceObject.invtemp,
                                                                 number_of_steps,
                                                                 bead_selector,
                                                                 local_seed,
                                                                 hardwall_int)

        total_accepted = total_accepted + accepted_moves
        total_proposed = total_proposed + number_of_steps


        # set energy
        current_energy = new_energy 

        latticeObject.chains[chainID].positions = idx_to_bead[:,5:].tolist()

        return (latticeObject, current_energy, total_proposed, total_accepted)
