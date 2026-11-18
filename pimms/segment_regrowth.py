## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
The path tables, proposal and window list behind ``MOVE_SEGMENT_REGROWTH``.

``MOVE_SEGMENT_REGROWTH`` (move code 15) moves the free beads of partly frozen
chains (``FREEZE_FILE`` ``B`` lines) a short window at a time. For each partly
frozen chain it picks one window of consecutive free beads, takes those beads
off the lattice and puts them back on a path drawn uniformly from every lattice
path of that many beads that joins the window's fixed neighbours. The new
arrangement is then accepted or rejected on the usual Metropolis criterion.

Why this samples the Boltzmann distribution exactly: the windows a chain can
offer are fixed at start-up, and the set of paths a window can be given depends
only on how far apart its two fixed neighbours are - which the move never
changes - and not on where the window's own beads are. So the chance of
proposing arrangement B from arrangement A is the chance of proposing A from B,
and Metropolis acceptance gives detailed balance. A path that lands on an
occupied site or crosses a hard wall is a rejected proposal. (Dropping such
paths from the set before drawing would be exact too: the window is off the
lattice while it is regrown, so the set would depend only on the rest of the
system, which is the same before and after the move. But it means testing
every candidate against the grid, up to 49,288 of them for a four-bead window
in 3D, where drawing first tests one.) A draw of the path the window already
has is rejected too, as a null move. A one-bead window proposes the same new
sites as the crankshaft does for that bead (the crankshaft also draws, and
rejects, its neighbours' own sites).

Why it is needed: the crankshaft moves one bead at a time, and that cannot
reach every arrangement of two or more free beads held at both ends by frozen
beads. Two free beads closing a 2x2 square with their frozen neighbours, for
example, never move. A window as long as such a stretch can reach any
arrangement of it in one move.

One limit: in a periodic box no wider than 2 x (stretch + 1) sites, a stretch
of free beads can join the two frozen beads at its ends going either way round
the box. The offset between a window's fixed neighbours is read by walking the
window's bonds and the new path keeps it, so this move (like the crankshaft, and
any local move that keeps the frozen beads fixed) keeps a stretch going the way
it already goes. That is still exact sampling, but only of the arrangements that
go that way round.

This module holds the parts that do not touch the lattice: the path tables,
the drawing of a new path for one window, and the window list. The move itself
(placing the beads, the energy change and the accept/reject) is
``MoveObject.segment_regrowth`` in moves.py.
"""

from __future__ import annotations

import functools
import itertools
import random
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt

from pimms import CONFIG
from pimms.data_structures import contiguous_runs
from pimms.latticeExceptions import MoveException

# Two beads are bonded when their sites differ by at most one lattice unit along
# every axis (the Chebyshev-1, or Moore, neighbourhood: 8 neighbours in 2D and 26
# in 3D). Every PIMMS move keeps this bond rule; it is the same 3x3(x3)
# neighbourhood the crankshaft draws its proposals from.
BOND_COMPONENTS = (-1, 0, 1)

# Relative coordinates of a window path stay within MAX_SEGMENT_REGROWTH_BEADS + 1
# lattice units of the anchor, so they fit in an int8.
RELATIVE_COORDINATE_DTYPE = np.int8

# Indices into the open-walk table; the largest table (3D, four beads) holds
# 388,994 walks.
PATH_INDEX_DTYPE = np.int32

# Only 2D and 3D lattices exist.
SUPPORTED_DIMENSIONS = (2, 3)

IntArray = npt.NDArray[np.int_]
CoordinateArray = npt.NDArray[np.int8]
IndexArray = npt.NDArray[np.int32]


def _check_table_request(n_dim: int, n_beads: int) -> None:
    """
    Refuse a path-table request for an unsupported lattice or window size.

    Parameters
    ----------
    n_dim : int
        Number of lattice dimensions.
    n_beads : int
        Number of beads in the window.

    Returns
    -------
    None

    Raises
    ------
    MoveException
        If ``n_dim`` is not 2 or 3, or ``n_beads`` is not between 1 and
        ``CONFIG.MAX_SEGMENT_REGROWTH_BEADS``.
    """
    if n_dim not in SUPPORTED_DIMENSIONS:
        raise MoveException('segment regrowth: a lattice has 2 or 3 dimensions, not %r' % (n_dim,))
    if not 1 <= n_beads <= CONFIG.MAX_SEGMENT_REGROWTH_BEADS:
        raise MoveException('segment regrowth: a window holds 1 to %d beads, not %r'
                            % (CONFIG.MAX_SEGMENT_REGROWTH_BEADS, n_beads))


@functools.lru_cache(maxsize=None)
def bond_vectors(n_dim: int) -> CoordinateArray:
    """
    Every bond vector of the lattice.

    Parameters
    ----------
    n_dim : int
        Number of lattice dimensions (2 or 3).

    Returns
    -------
    numpy.ndarray
        Shape (8, 2) or (26, 3), int8, read-only: every vector whose components
        are -1, 0 or 1, except the zero vector.
    """
    _check_table_request(n_dim, 1)
    vectors = np.array([step for step in itertools.product(BOND_COMPONENTS, repeat=n_dim) if any(step)],
                       dtype=RELATIVE_COORDINATE_DTYPE)
    vectors.setflags(write=False)
    return vectors


@functools.lru_cache(maxsize=None)
def open_walks(n_dim: int, n_beads: int) -> CoordinateArray:
    """
    Every self-avoiding walk of ``n_beads`` sites from an anchor at the origin.

    Site 0 of a walk is bonded to the origin and each later site to the one
    before it. No site is the origin and no two sites are the same. These are
    the paths a window at a chain end can be given (its one fixed neighbour is
    the origin), and the closed paths of :func:`closed_path_indices` are a
    subset of them.

    Parameters
    ----------
    n_dim : int
        Number of lattice dimensions (2 or 3).
    n_beads : int
        Number of beads in the window, 1 to ``CONFIG.MAX_SEGMENT_REGROWTH_BEADS``.

    Returns
    -------
    numpy.ndarray
        Shape (n_walks, n_beads, n_dim), int8, read-only, coordinates relative
        to the origin. There are 8, 56, 368 and 2,336 walks of 1 to 4 beads in
        2D and 26, 650, 15,986 and 388,994 in 3D.
    """
    _check_table_request(n_dim, n_beads)
    steps = bond_vectors(n_dim)
    walks = steps[:, np.newaxis, :].copy()
    for _ in range(n_beads - 1):
        # extend every walk by every bond vector, then drop the extensions that
        # step back onto the origin or onto a site the walk already holds
        previous = np.repeat(walks, len(steps), axis=0)
        extended = (walks[:, -1, np.newaxis, :] + steps[np.newaxis, :, :]).reshape(-1, n_dim)
        on_origin = np.all(extended == 0, axis=1)
        on_walk = np.any(np.all(previous == extended[:, np.newaxis, :], axis=2), axis=1)
        keep = ~(on_origin | on_walk)
        walks = np.concatenate([previous[keep], extended[keep][:, np.newaxis, :]], axis=1)
    walks.setflags(write=False)
    return walks


@functools.lru_cache(maxsize=None)
def closed_path_indices(n_dim: int, n_beads: int, offset: Tuple[int, ...]) -> IndexArray:
    """
    The paths of ``n_beads`` sites joining an anchor at the origin to one at ``offset``.

    These are the paths an interior window can be given: the open walks whose
    last site is bonded to the second anchor and which do not pass through it.
    Cached per offset as the run meets them; in 3D a four-bead window has at
    most 49,288 paths for one offset and 9.4 million over all 1,330 offsets.

    Parameters
    ----------
    n_dim : int
        Number of lattice dimensions (2 or 3).
    n_beads : int
        Number of beads in the window.
    offset : tuple of int
        Position of the second anchor relative to the first, unwrapped (the
        sum of the bond vectors from the first anchor through the window to
        the second).

    Returns
    -------
    numpy.ndarray
        int32, read-only: indices into ``open_walks(n_dim, n_beads)``, in
        ascending order. Empty if no path of this length can join the anchors.

    Raises
    ------
    MoveException
        If ``offset`` has the wrong number of components or is the origin
        (two beads can never share a site).
    """
    if len(offset) != n_dim:
        raise MoveException('segment regrowth: anchor offset %r does not have %d components' % (offset, n_dim))
    if not any(offset):
        raise MoveException('segment regrowth: the two anchors of a window are on the same site')
    walks = open_walks(n_dim, n_beads)
    target = np.asarray(offset, dtype=np.int64)
    gap_to_target = np.abs(walks[:, -1, :].astype(np.int64) - target).max(axis=1)
    through_target = np.any(np.all(walks == target, axis=2), axis=1)
    indices = np.flatnonzero((gap_to_target == 1) & ~through_target).astype(PATH_INDEX_DTYPE)
    indices.setflags(write=False)
    return indices


def bond_vector(site_from: Sequence[int], site_to: Sequence[int], dimensions: Sequence[int]) -> Tuple[int, ...]:
    """
    The bond vector from one bonded site to the next, read through the periodic wrap.

    Parameters
    ----------
    site_from, site_to : sequence of int
        Two sites in the box, bonded to each other.
    dimensions : sequence of int
        The box. Every axis is at least 3 sites (PIMMS requires 7), so a bond
        vector is never ambiguous.

    Returns
    -------
    tuple of int
        The vector, every component -1, 0 or 1.

    Raises
    ------
    MoveException
        If the two sites are not bonded, which would mean the chain's stored
        positions are corrupt.
    """
    vector = []
    for start, end, length in zip(site_from, site_to, dimensions):
        component = int(end) - int(start)
        if component > 1:
            component -= int(length)
        elif component < -1:
            component += int(length)
        if component not in BOND_COMPONENTS:
            raise MoveException('segment regrowth: sites %s and %s are not bonded in a box of %s'
                                % (list(site_from), list(site_to), list(dimensions)))
        vector.append(component)
    if not any(vector):
        raise MoveException('segment regrowth: two consecutive beads share site %s' % (list(site_from),))
    return tuple(vector)


def window_anchor_offset(chain_positions: Sequence[Sequence[int]], first: int, last: int,
                         dimensions: Sequence[int]) -> Tuple[int, ...]:
    """
    How far the second fixed neighbour of an interior window is from the first.

    Read by walking the bonds from bead ``first - 1`` through the window to bead
    ``last + 1``, so a window that crosses a periodic boundary gets its true
    (unwrapped) offset. A new path drawn for the window keeps the offset, so the
    reverse move draws from the same set of paths.

    Parameters
    ----------
    chain_positions : sequence of sequence of int
        The chain's bead positions, in chain order.
    first, last : int
        The window, inclusive; ``first >= 1`` and ``last <= len - 2``.
    dimensions : sequence of int
        The box.

    Returns
    -------
    tuple of int
        Position of bead ``last + 1`` relative to bead ``first - 1``.
    """
    offset = [0] * len(dimensions)
    for bead in range(first - 1, last + 1):
        for axis, component in enumerate(bond_vector(chain_positions[bead], chain_positions[bead + 1], dimensions)):
            offset[axis] += component
    return tuple(offset)


def draw_window_path(chain_positions: Sequence[Sequence[int]], first: int, last: int,
                     dimensions: Sequence[int], hardwall: bool) -> Optional[List[List[int]]]:
    """
    Draw new positions for a window of beads, uniformly from every path it can take.

    An interior window (fixed neighbours on both sides) is given one of the
    paths joining its neighbours; a window at a chain end (one fixed neighbour)
    one of the walks leaving its neighbour. The path is drawn with one call to
    Python's ``random.randrange``, whatever the configuration. Whether its sites
    are free is not checked here - the caller does that against the grid.

    Parameters
    ----------
    chain_positions : sequence of sequence of int
        The chain's bead positions, in chain order, in the box.
    first, last : int
        The window, inclusive. It must leave at least one bead of the chain
        outside it.
    dimensions : sequence of int
        The box.
    hardwall : bool
        If True a path that leaves the box is refused (None is returned);
        otherwise its sites are wrapped back into the box.

    Returns
    -------
    list of list of int, or None
        The new positions of beads ``first`` to ``last``, in chain order and in
        the box, or None if the path drawn crosses a hard wall (a rejected
        proposal).

    Raises
    ------
    MoveException
        If the window is not inside the chain or covers all of it, or the
        current arrangement is not among the paths (which would mean the
        stored positions are corrupt).
    """
    chain_length = len(chain_positions)
    n_beads = last - first + 1
    n_dim = len(dimensions)
    if not (0 <= first <= last < chain_length) or n_beads == chain_length:
        raise MoveException('segment regrowth: window %d-%d does not fit a chain of %d beads with a bead to spare'
                            % (first, last, chain_length))

    walks = open_walks(n_dim, n_beads)
    if first > 0 and last < chain_length - 1:
        # interior window: paths from bead first-1 that end bonded to bead last+1
        anchor = chain_positions[first - 1]
        candidates = closed_path_indices(n_dim, n_beads, window_anchor_offset(chain_positions, first, last,
                                                                              dimensions))
        if len(candidates) == 0:
            raise MoveException('segment regrowth: no path of %d beads joins the neighbours of window %d-%d, '
                                'yet the chain holds one' % (n_beads, first, last))
        relative = walks[candidates[random.randrange(len(candidates))]]
    elif first > 0:
        # window at the C-terminal end: a walk away from bead first-1
        anchor = chain_positions[first - 1]
        relative = walks[random.randrange(len(walks))]
    else:
        # window at the N-terminal end: a walk away from bead last+1, whose
        # first site is bead last and whose last site is bead 0
        anchor = chain_positions[last + 1]
        relative = walks[random.randrange(len(walks))][::-1]

    unwrapped = np.asarray(anchor, dtype=np.int64) + relative.astype(np.int64)
    box = np.asarray(dimensions, dtype=np.int64)
    if hardwall:
        if np.any(unwrapped < 0) or np.any(unwrapped >= box):
            return None
        placed = unwrapped
    else:
        placed = np.mod(unwrapped, box)
    return [[int(component) for component in site] for site in placed]


def free_segment_windows(chain_length: int, frozen_bead_indices: Iterable[int],
                         max_beads: int) -> List[Tuple[int, int]]:
    """
    Every window of consecutive free beads that a partly frozen chain offers.

    A window is one to ``max_beads`` consecutive beads, none of them frozen.
    The list depends only on which beads are frozen, never on where any bead
    is, which is what keeps the window choice from biasing the sampling.

    Parameters
    ----------
    chain_length : int
        Number of beads in the chain.
    frozen_bead_indices : iterable of int
        The chain's frozen beads, counted from 0. At least one, and not all.
    max_beads : int
        The longest window, 1 to ``CONFIG.MAX_SEGMENT_REGROWTH_BEADS``
        (``SEGMENT_REGROWTH_MAX_BEADS``).

    Returns
    -------
    list of tuple of (int, int)
        The inclusive (first, last) windows, ordered by first bead and then
        by length. Every free bead is in at least one (its one-bead window).

    Raises
    ------
    MoveException
        If ``max_beads`` is out of range, a frozen bead is not in the chain, or
        the chain has no frozen bead or no free one.
    """
    if not 1 <= max_beads <= CONFIG.MAX_SEGMENT_REGROWTH_BEADS:
        raise MoveException('segment regrowth: windows hold 1 to %d beads, not %r'
                            % (CONFIG.MAX_SEGMENT_REGROWTH_BEADS, max_beads))
    frozen = set(int(index) for index in frozen_bead_indices)
    if not frozen or not frozen <= set(range(chain_length)) or len(frozen) == chain_length:
        raise MoveException('segment regrowth: a partly frozen chain of %d beads needs some, but not all, of its '
                            'beads frozen (got %s)' % (chain_length, sorted(frozen)))
    free = [bead for bead in range(chain_length) if bead not in frozen]
    windows = []
    for run_first, run_last in contiguous_runs(free):
        for first in range(run_first, run_last + 1):
            for last in range(first, min(first + max_beads, run_last + 1)):
                windows.append((first, last))
    return windows


def segment_windows(chain_lengths: Dict[int, int], frozen_beads: Dict[int, Sequence[int]],
                    max_beads: int) -> Dict[int, List[Tuple[int, int]]]:
    """
    The windows of every partly frozen chain, as ``MOVE_SEGMENT_REGROWTH`` uses them.

    Parameters
    ----------
    chain_lengths : dict of int to int
        chainID -> number of beads, for (at least) every partly frozen chain.
    frozen_beads : dict of int to sequence of int
        chainID -> frozen beads, for the partly frozen chains only (the
        ``frozen_beads`` of :meth:`FreezeFile.resolve`). Chains frozen whole are
        not in it and get no windows.
    max_beads : int
        The longest window (``SEGMENT_REGROWTH_MAX_BEADS``).

    Returns
    -------
    dict of int to list of tuple of (int, int)
        chainID -> its windows (see :func:`free_segment_windows`). Empty when
        no chain is partly frozen, and the move then does nothing.
    """
    return {chainID: free_segment_windows(chain_lengths[chainID], beads, max_beads)
            for chainID, beads in sorted(frozen_beads.items())}
