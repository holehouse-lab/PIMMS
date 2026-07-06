"""
Periodic-boundary edge-case tests.

These focus on the regime the box is *not* normally used in: chains large enough
that, when unwrapped, they span more than the box in one or several axes - i.e. a
chain that straddles the periodic boundary across multiple faces / wraps the box
more than once. The goal is to pin down where PIMMS is correct, where it degrades
gracefully (a warning / a bounded exception), and to guard against silent
corruption or hangs.

What these tests establish (verified empirically before they were written):

* **Single-image reconstruction is robust.** ``convert_chain_to_single_image`` and
  ``make_chain_whole`` rebuild a chain that spans an arbitrary number of periodic
  images - in one axis or several at once - recovering the true unwrapped shape
  exactly, and re-wrapping (mod box) returns the original positions bit-for-bit.

* **Rg/asphericity analysis degrades *gracefully, with a warning*.**
  ``Chain.analysis_get_polymeric_properties`` returns the minimum-image value,
  which COLLAPSES once a chain exceeds ~half the box (min-image folds far beads
  back). The single-image value stays exact, and the code prints a ``[WARNING]``
  about finite-size artefacts whenever the two disagree. So the wrong number is
  never returned silently.

* **Impossible / degenerate inputs fail loudly and quickly.** A non-wrappable bond
  raises ``LatticeUtilsException`` (bounded, no hang); a disconnected cluster
  raises a clear ``ValueError``; empty/single-bead inputs are handled.

* **Moves keep the energy self-consistent under heavy straddling.** A full
  simulation whose chains wrap a small box, run with ``ENERGY_CHECK`` every step
  across every move type, never desynchronises the tracked energy.
"""

import io
import time
import contextlib

import numpy as np
import pytest

from pimms import lattice_utils as lu
from pimms import lattice_analysis_utils as lau
from pimms import cluster_utils as cu
from pimms.latticeExceptions import LatticeUtilsException
from pimms.tests import kernel_test_utils as U


# ---------------------------------------------------------------------------
# constructors for chains that wrap the box a controlled number of times
# ---------------------------------------------------------------------------

def _rod_wrapping(box_len, n_periods, axis):
    """A straight rod of ``n_periods * box_len`` beads pointing along ``axis``.

    The rod is stepped by +1 on a perpendicular axis every ``box_len`` beads so the
    WRAPPED positions stay self-avoiding while the UNWRAPPED extent along ``axis``
    is ``n_periods * box_len`` - i.e. the rod crosses that boundary ``n_periods``
    times. Returns ``(dims, unwrapped, wrapped)``.
    """
    dims = [box_len, box_len, box_len]
    perp = (axis + 1) % 3
    assert n_periods <= dims[perp], "perp carry would self-collide"
    n_beads = n_periods * box_len
    unwrapped = []
    for i in range(n_beads):
        p = [0, 0, 0]
        p[axis] = i
        p[perp] = i // box_len
        unwrapped.append(p)
    wrapped = [[c % dims[d] for d, c in enumerate(p)] for p in unwrapped]
    return dims, unwrapped, wrapped


def _odometer_chain(base, n_beads):
    """A base-``base`` odometer walk: a self-avoiding chain (for ``n_beads <=
    base**3``) whose consecutive beads increment like an odometer, so it crosses
    the z-boundary every ``base`` beads and the x-boundary every ``base**2`` beads
    - i.e. it wraps the box in MORE THAN ONE axis. Returns ``(dims, unwrapped,
    wrapped)``.
    """
    dims = [base, base, base]
    unwrapped = [[i // base, i // (base * base), i] for i in range(n_beads)]
    wrapped = [[c % base for c in p] for p in unwrapped]
    return dims, unwrapped, wrapped


def _rel(points):
    a = np.asarray(points)
    return a - a[0]


def _axes_crossed(wrapped):
    crossed = set()
    for a, b in zip(wrapped, wrapped[1:]):
        for d in range(len(a)):
            if abs(a[d] - b[d]) > 1:
                crossed.add(d)
    return crossed


# ---------------------------------------------------------------------------
# 1. single-image reconstruction of multi-image chains is exact
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("n_periods", [1, 2, 3, 4])
def test_single_image_recovers_multiimage_chain(axis, n_periods):
    dims, unwrapped, wrapped = _rod_wrapping(6, n_periods, axis)

    # the constructed configuration is a valid (self-avoiding) lattice chain
    assert len(wrapped) == len(set(map(tuple, wrapped)))
    # it straddles the boundary exactly when it spans more than one period
    assert lu.do_positions_stradle_pbc_boundary(wrapped) is (n_periods > 1)

    # both reconstruction routines recover the true unwrapped shape (up to the
    # arbitrary translation each applies)
    for reconstruct in (lu.convert_chain_to_single_image, lu.make_chain_whole):
        rebuilt = reconstruct(wrapped, dims)
        assert np.array_equal(_rel(rebuilt), _rel(unwrapped))


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_unwrap_then_rewrap_is_identity(axis):
    # unwrapping only ever shifts a coordinate by whole multiples of the box, so
    # re-wrapping (mod box) must return the original positions exactly
    dims, _unwrapped, wrapped = _rod_wrapping(6, 4, axis)
    for reconstruct in (lu.convert_chain_to_single_image, lu.make_chain_whole):
        si = reconstruct(wrapped, dims)
        rewrapped = [[int(c) % dims[d] for d, c in enumerate(p)] for p in si]
        assert rewrapped == [list(p) for p in wrapped]


def test_multi_axis_wrap_is_reconstructed():
    # a chain that wraps the box in more than one axis at once
    dims, unwrapped, wrapped = _odometer_chain(4, 40)
    assert len(wrapped) == len(set(map(tuple, wrapped)))       # self-avoiding
    assert len(_axes_crossed(wrapped)) >= 2                    # genuinely multi-axis
    rebuilt = lu.convert_chain_to_single_image(wrapped, dims)
    assert np.array_equal(_rel(rebuilt), _rel(unwrapped))


# ---------------------------------------------------------------------------
# 2. Rg analysis: single-image stays exact, min-image collapses + warns
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_periods,diverges", [(1, False), (2, True), (3, True)])
def test_minimage_rg_collapses_while_single_image_stays_exact(n_periods, diverges):
    dims, unwrapped, wrapped = _rod_wrapping(6, n_periods, axis=2)

    rg_true = lau.get_polymeric_properties(unwrapped, dims)[0]
    rg_single = lau.get_polymeric_properties(lu.convert_chain_to_single_image(wrapped, dims), dims)[0]
    rg_minimage = lau.get_polymeric_properties(wrapped, dims)[0]

    # single-image Rg equals the free-chain Rg (shape is recovered exactly)
    assert rg_single == pytest.approx(rg_true, abs=1e-9)

    if diverges:
        # min-image folds the far beads back, underestimating Rg
        assert rg_minimage < rg_true - 0.01
    else:
        assert rg_minimage == pytest.approx(rg_true, abs=1e-9)


def test_finite_size_warning_is_emitted(tmp_path):
    # build a real Chain, then drive it into a >half-box configuration and confirm
    # the finite-size artefact warning fires (graceful degradation, not silence)
    state = U.build_state(
        tmp_path, 3, "SR", hardwall=False,
        moves={"MOVE_CRANKSHAFT": 1.0}, box=[8, 8, 8], chains=[(1, "A" * 24)],
    )
    chain = next(iter(state.lattice.chains.values()))

    _dims, _unw, spanning = _rod_wrapping(8, 3, axis=2)   # spans three box-periods
    chain.positions = spanning
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rg_minimage, _asph = chain.analysis_get_polymeric_properties()[:2]
    out = buf.getvalue()
    assert "WARNING" in out and "finite size" in out
    # the returned value is the (collapsed) min-image one - documented behaviour
    assert rg_minimage < lau.get_polymeric_properties(spanning, _dims)[0] + 1e-9

    # a compact, in-box configuration must NOT warn
    _dims2, _unw2, compact = _rod_wrapping(8, 1, axis=2)
    chain.positions = compact
    buf2 = io.StringIO()
    with contextlib.redirect_stdout(buf2):
        chain.analysis_get_polymeric_properties()
    assert "WARNING" not in buf2.getvalue()


# ---------------------------------------------------------------------------
# 3. impossible / degenerate inputs fail loudly and quickly
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("reconstruct", [lu.convert_chain_to_single_image, lu.make_chain_whole])
def test_impossible_bond_raises_and_does_not_hang(reconstruct):
    # a gap of 3 in a box of 10 is not resolvable by any whole-box shift
    impossible = [[0, 0, 0], [3, 0, 0]]
    start = time.time()
    with pytest.raises(LatticeUtilsException):
        reconstruct(impossible, [10, 10, 10])
    assert time.time() - start < 5.0            # bounded escape counter, no hang


def test_snakesearch_disconnected_cluster_raises_clearly():
    with pytest.raises(ValueError, match="connected"):
        cu.convert_positions_to_single_image_snakesearch(
            [[0, 0, 0], [5, 5, 5]], [10, 10, 10], space_threshold=1)


def test_degenerate_inputs_are_handled():
    assert cu.convert_positions_to_single_image_snakesearch([], [10, 10, 10]) == []
    assert lu.do_positions_stradle_pbc_boundary([[2, 2, 2]]) is False
    assert [list(p) for p in lu.make_chain_whole([[1, 1, 1]], [8, 8, 8])] == [[1, 1, 1]]


# ---------------------------------------------------------------------------
# 4. engine level: moves keep the energy consistent under heavy straddling
# ---------------------------------------------------------------------------

# every move that operates on whole chains / clusters, so the single-image
# reconstruction path is exercised inside the moves under boundary crossing
_ALL_MOVES = {
    "MOVE_CRANKSHAFT": 0.3,
    "MOVE_SLITHER": 0.2,
    "MOVE_PULL": 0.2,
    "MOVE_CHAIN_TRANSLATE": 0.1,
    "MOVE_CHAIN_ROTATE": 0.1,
    "MOVE_CLUSTER_TRANSLATE": 0.05,
    "MOVE_CLUSTER_ROTATE": 0.05,
}


@pytest.mark.parametrize("ff", ["SR", "LR", "SLR"])
def test_moves_preserve_energy_under_heavy_straddle(tmp_path, ff):
    # 16-mers packed into the smallest supported box (8, the parser minimum) cannot
    # avoid straddling; ENERGY_CHECK=1 recomputes the energy from scratch after every
    # step and raises on any desync, so a clean return proves every move handled the
    # boundary-crossing chains correctly.
    trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, ff, hardwall=False, moves=_ALL_MOVES,
        box=[8, 8, 8], chains=[(4, "AABB" * 4)],
        n_steps=200, energy_check=1, return_sim=True,
    )
    assert sim.LATTICE.any_chains_straddle_boundary(), (
        "test is vacuous - no chain crossed a boundary")


def test_hardwall_moves_consistent_and_no_straddle(tmp_path):
    # sanity mirror: with a hardwall there is no periodic image to cross into, so
    # no chain may straddle, and the energy must still stay consistent
    trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SLR", hardwall=True,
        moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_SLITHER": 0.4},
        box=[8, 8, 8], chains=[(2, "AABB" * 3)],
        n_steps=150, energy_check=1, return_sim=True,
    )
    assert not sim.LATTICE.any_chains_straddle_boundary()
