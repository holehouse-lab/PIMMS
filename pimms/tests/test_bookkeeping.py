"""The compiled bead-table bookkeeping must be bit-identical to the Python it replaced.

Every megamove copies each chain's positions from the Chain objects (a Python
list of ``[x, y(, z)]`` lists) into the compiled kernels' bead table, and copies
the moved positions back afterwards. Those two copies are now compiled
(``pimms.bookkeeping``), and the per-chain metadata the whole-chain megamoves
used to rebuild every call now comes from a layout the lattice builds once.

None of that may change a single number: the copies are pure data movement,
and the metadata is invariant during a run. So the pure-Python loops are kept
as the oracle here, and the decisive check at the end runs whole megamoves
both ways from the same seed and demands byte-identical grids, positions and
energies. The regression fixtures pin the same thing end to end.
"""

import contextlib
import copy
import os

import numpy as np
import pytest

from pimms import crankshaft_list_functions as clf
from pimms.tests import kernel_test_utils as U

_MIXED = [(6, "AABBA"), (5, "A"), (4, "BB"), (3, "ABABABAB")]


def _build(tmp_path, dim, hardwall, chains, seed=5):
    os.makedirs(tmp_path, exist_ok=True)      # callers pass fresh subdirectories
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        return U.build_state(tmp_path, dim, "SR", hardwall, {"MOVE_CRANKSHAFT": 1.0},
                             box=[14] * dim, chains=chains, seed=seed, temperature=60)


def test_the_compiled_bookkeeping_is_what_the_suite_exercises():
    """If the extension were missing every test below would compare Python with Python."""
    assert clf._HAVE_BOOKKEEPING, "pimms.bookkeeping is not built; rebuild the extensions"


@pytest.mark.parametrize("dim,hardwall", [(3, False), (2, False), (3, True), (2, True)])
def test_gather_refreshes_every_position_exactly_like_the_python_loop(tmp_path, dim, hardwall):
    state = _build(tmp_path, dim, hardwall, _MIXED)
    lattice = state.lattice

    # move a few chains by hand so the table is genuinely stale
    for cid in list(lattice.chains)[::2]:
        chain = lattice.chains[cid]
        chain.set_ordered_positions([[c + 1 for c in p] for p in chain.get_ordered_positions()])

    expected = clf._update_idx_to_bead_python(lattice)

    # scribble over the position columns so a gather that skipped any row shows
    lattice.crankshaft_lists[:, 5:] = -12345
    got = clf.update_idx_to_bead(lattice)

    assert got.dtype == np.int64 and got.flags["C_CONTIGUOUS"]
    np.testing.assert_array_equal(got, expected)
    # and it is a fresh copy, not the live table
    assert got is not lattice.crankshaft_lists
    got[0, 5] = 999
    assert lattice.crankshaft_lists[0, 5] != 999


@pytest.mark.parametrize("dim", [3, 2])
def test_scatter_writes_back_exactly_what_the_python_loop_writes(tmp_path, dim):
    state_a = _build(tmp_path / "a", dim, False, _MIXED)
    state_b = _build(tmp_path / "b", dim, False, _MIXED)
    table = clf.update_idx_to_bead(state_a.lattice)
    rng = np.random.default_rng(3)
    table[:, 5:] = rng.integers(0, 14, size=table[:, 5:].shape)

    clf._write_back_positions_python(state_a.lattice, table.copy())
    clf.write_back_positions(state_b.lattice, table.copy())

    for cid in state_a.lattice.chains:
        pa = state_a.lattice.chains[cid].get_ordered_positions()
        pb = state_b.lattice.chains[cid].get_ordered_positions()
        assert pa == pb
        # downstream code relies on plain Python list-of-lists of ints, not numpy
        assert type(pb) is list
        assert all(type(row) is list and all(type(c) is int for c in row) for row in pb)
        assert len(pb) == len(state_b.lattice.chains[cid])


def test_layout_is_built_once_and_matches_the_chains(tmp_path):
    state = _build(tmp_path, 3, False, _MIXED)
    lattice = state.lattice
    layout = clf.chain_layout(lattice)
    assert layout is clf.chain_layout(lattice)          # cached
    assert layout.sorted_ids == sorted(lattice.chains)
    lengths = [len(lattice.chains[c]) for c in layout.sorted_ids]
    assert list(layout.length) == lengths
    assert list(layout.offset) == list(np.cumsum([0] + lengths[:-1]))
    assert layout.n_beads == sum(lengths)
    # homopolymer flag: uniform intcode and LR flag over the chain
    table = lattice.crankshaft_lists
    for ci, (off, L) in enumerate(zip(layout.offset, layout.length)):
        seg = table[off:off + L]
        expect = int(np.all(seg[:, 2] == seg[0, 2]) and np.all(seg[:, 1] == seg[0, 1]))
        assert int(layout.homo[ci]) == expect
    assert layout.offset32.dtype == np.int32 and layout.length32.dtype == np.int32


def test_layout_survives_non_contiguous_chain_ids_and_a_rebuilt_table(tmp_path):
    """Deleting a chain and rebuilding the table must give a fresh, correct layout."""
    state = _build(tmp_path, 3, False, _MIXED)
    lattice = state.lattice
    old_layout = clf.chain_layout(lattice)
    victim = sorted(lattice.chains)[3]
    del lattice.chains[victim]
    lattice.crankshaft_lists = clf.initialize_idx_to_bead(lattice)
    lattice.chain_to_firstbead_lookup = clf.initialize_chain_to_firstbead_lookup(lattice)

    layout = clf.chain_layout(lattice)
    assert layout is not old_layout
    assert victim not in layout.sorted_ids
    assert layout.n_beads == len(lattice.crankshaft_lists)

    expected = clf._update_idx_to_bead_python(lattice)
    lattice.crankshaft_lists[:, 5:] = -1
    np.testing.assert_array_equal(clf.update_idx_to_bead(lattice), expected)


def test_bead_selector_random_stream_is_unchanged(tmp_path):
    """The selector must draw exactly what the old per-chain loop drew, frozen or not."""
    state = _build(tmp_path, 3, False, _MIXED)
    lattice = state.lattice
    n = len(lattice.crankshaft_lists)
    frozen = tuple(sorted(lattice.chains)[1::4])

    def old_selector(number_of_steps, frozen_chains):
        if len(frozen_chains) == 0:
            return np.random.randint(0, n, number_of_steps)
        c, sel = 0, []
        for cid in sorted(lattice.chains.keys()):
            if cid in frozen_chains:
                c += len(lattice.chains[cid])
            else:
                for _ in range(len(lattice.chains[cid])):
                    sel.append(c)
                    c += 1
        return np.random.choice(sel, number_of_steps, replace=True)

    for fz in ((), frozen):
        np.random.seed(77)
        new = clf.bead_selector_constructor(n, 5000, lattice, frozen_chains=fz, safecheck=True)
        np.random.seed(77)
        old = old_selector(5000, fz)
        np.testing.assert_array_equal(new, old)
        if fz:
            frozen_rows = {r for cid in fz for r in range(
                lattice.chain_to_firstbead_lookup[cid],
                lattice.chain_to_firstbead_lookup[cid] + len(lattice.chains[cid]))}
            assert not (set(new.tolist()) & frozen_rows)


def test_gather_refuses_a_chain_that_is_out_of_step_with_the_table(tmp_path):
    state = _build(tmp_path, 3, False, _MIXED)
    lattice = state.lattice
    cid = sorted(lattice.chains)[0]
    chain = lattice.chains[cid]
    chain.positions = chain.positions[:-1]          # bypass the length check on purpose
    with pytest.raises(ValueError, match="out of step"):
        clf.update_idx_to_bead(lattice)


@pytest.mark.parametrize("move", ["crank", "slither", "pull"])
@pytest.mark.parametrize("parallel", [False, True], ids=["serial", "parallel"])
def test_megamoves_are_byte_identical_with_and_without_compiled_bookkeeping(
        tmp_path, monkeypatch, move, parallel):
    """The decisive check: whole megamoves, both ways, same seed, identical state."""
    import random

    def run(root, compiled):
        monkeypatch.setattr(clf, "_HAVE_BOOKKEEPING", compiled)
        os.makedirs(root, exist_ok=True)
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            state = U.build_state(root, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                                  box=[24, 24, 24],
                                  chains=[(30, "AABBA"), (20, "A"), (10, "BBBBBBBB")],
                                  seed=9, temperature=60)
        lattice, ham, acc = state.lattice, state.ham, state.acc
        energy = int(ham.evaluate_total_energy(lattice)[0])
        mover = state.sim.MOVER
        before = {c: copy.deepcopy(lattice.chains[c].get_ordered_positions())
                  for c in lattice.chains}
        random.seed(11)
        np.random.seed(11)
        kw = dict(parallelize=parallel, num_threads=2)
        for _ in range(4):
            if move == "crank":
                lattice, energy, *_ = mover.system_shake(lattice, energy, acc, ham, 3000,
                                                         "UNIFORM", frozen_chains=(4, 9), **kw)
            elif move == "slither":
                lattice, energy, *_ = mover.system_slither(lattice, energy, acc, ham, 6,
                                                           frozen_chains=(4, 9), **kw)
            else:
                lattice, energy, *_ = mover.system_pull(lattice, energy, acc, ham, 6,
                                                        frozen_chains=(4, 9), **kw)
        positions = {c: copy.deepcopy(lattice.chains[c].get_ordered_positions())
                     for c in lattice.chains}
        moved = sum(positions[c] != before[c] for c in positions)
        frozen_moved = any(positions[c] != before[c] for c in (4, 9))
        return (energy, np.asarray(lattice.grid).copy(),
                np.asarray(lattice.type_grid).copy(), positions, moved, frozen_moved)

    e_py, g_py, t_py, p_py, moved_py, fz_py = run(tmp_path / "py", False)
    e_cy, g_cy, t_cy, p_cy, moved_cy, fz_cy = run(tmp_path / "cy", True)
    assert e_py == e_cy
    np.testing.assert_array_equal(g_py, g_cy)
    np.testing.assert_array_equal(t_py, t_cy)
    assert p_py == p_cy
    # the megamoves must have moved a fair number of chains, or the comparison
    # above would be between two untouched copies of the start configuration
    # (pull skips the 20 monomers and, at this temperature and substep count,
    # rearranges only a dozen or so of the remaining chains)
    assert moved_cy == moved_py and moved_cy >= 5
    # and the frozen mask, which the layout now feeds, must still hold
    assert not fz_py and not fz_cy
