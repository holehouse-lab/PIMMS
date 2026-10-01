"""
The PARALLELIZE startup report must describe the implementation that is actually
used for THIS box/system: thread budget, OpenMP availability, the crankshaft block
decomposition (or the single-block warning), and, per whole-chain move, whether it
runs on the parallel kernel or falls back to serial - and why.

``build_state`` silences the Simulation constructor's stdout, so the tests read the
report from its return value and from ``log.txt`` (every line is also logged).
"""

import contextlib
import io
import os

import pytest

import numpy as np

from pimms import mega_crank_fast, moves
from pimms.tests import kernel_test_utils as U


def _build(tmp_path, dim, ff, box, chains, moveset, extra, hardwall=False):
    base = {"PARALLELIZE": "True", "PARALLEL_THREADS": 4}
    base.update(extra)
    state = U.build_state(tmp_path, dim, ff, hardwall, moveset, box=box, chains=chains,
                          extra=base)
    log = (tmp_path / "log.txt").read_text()
    return state, log



@pytest.fixture(autouse=True)
def _restore_cwd():
    cwd = os.getcwd()
    yield
    os.chdir(cwd)

def test_openmp_info_shape():
    info = mega_crank_fast.openmp_info()
    assert set(info) == {"enabled", "max_threads"}
    assert isinstance(info["enabled"], bool) and info["max_threads"] >= 1


def test_report_absent_when_parallelize_off(tmp_path):
    U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                  extra={"PARALLELIZE": "False"})
    assert "PARALLELIZATION REPORT" not in (tmp_path / "log.txt").read_text()


def test_report_multiblock_crank_and_parallel_slither(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)     # the report is logged to log.txt in the working directory
    state, log = _build(tmp_path, 3, "LR", [40, 40, 40],
                        [(30, "AABB"), (30, "AAAA")],
                        {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5}, {})
    assert "PARALLELIZATION REPORT" in log          # emitted at construction
    # the method returns the lines it prints/logs; check the content there
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        lines = state.sim.report_parallelization(note="again")
    out = "\n".join(lines)
    assert lines[0] == "PARALLELIZATION REPORT - again"
    assert "PARALLELIZATION REPORT - again" in buf.getvalue()   # printed too
    assert "Threads: 4 OpenMP threads" in out
    assert "interaction radius 3 (long-range beads present)" in out
    # the crankshaft numbers must be the kernel's own layout
    lay = mega_crank_fast.parallel_crank_layout_info(40, 40, 40, True)
    assert "kernel mega_crank_parallel - per-bead" in out
    assert "halo W=%d; block grid %s = %d blocks" % (
        lay["W"], "x".join(map(str, lay["blocks"])), lay["num_blocks"]) in out
    assert "%.0f%% of the box movable per sweep" % (100 * lay["movable_fraction"]) in out
    # slither: the box splits for the chain-level halo and every one of the 60
    # chains is a 4-mer, so the whole system sits on the parallel side of the
    # partition. The report must describe a PARTITION - how many chains each
    # kernel takes and the length threshold - not a per-megamove choice of one
    # kernel for the move, which is what it used to claim ("-> PARALLEL") and
    # which was the bug: choosing a kernel from the current configuration broke
    # stationarity.
    chain_lay = mega_crank_fast.parallel_layout_info(40, 40, 40, True)
    interior = chain_lay["block_size"][0] - 2 * chain_lay["W"]
    assert "Slither (MOVE_SLITHER): kernel mega_slither_parallel" in out
    assert ("parallel kernel for 60 chain(s) of length <= %d, serial kernel for 0 "
            "longer chain(s) - both run every megamove" % interior) in out
    assert "the split is fixed by chain length for the whole run" in out
    assert "-> PARALLEL" not in out
    assert "falls back to the SERIAL kernel" not in out
    assert "Pull (MOVE_PULL): not in the move set" in out
    assert "relax more slowly per step" in out


def test_report_single_block_warning_for_small_box(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)     # the report is logged to log.txt in the working directory
    state, _ = _build(tmp_path, 3, "LR", [14, 14, 14], [(6, "AABB")],
                      {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5}, {})
    out = "\n".join(state.sim.report_parallelization())
    assert "ONE block -> runs single-threaded, equivalent to the serial kernel" in out
    assert "Slither (MOVE_SLITHER): box does not split for the chain-level halo" in out
    assert "runs on the SERIAL kernel" in out
    # Crank blocks are kept >= 8W long, so TWO blocks require at least 16W.
    assert "needs >= 32 sites in a dimension" in out


def test_single_block_gate_and_crank_dispatch_use_serial_kernel(tmp_path, monkeypatch):
    state = U.build_state(tmp_path, 3, "LR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[14, 14, 14], chains=[(2, "AABB")],
                          extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 4})
    idx, _chains, off, length, homo = moves.parallel_chain_metadata(state.lattice)
    report = moves.parallel_chain_fit_report(
        idx, off, length, state.lattice.dimensions, True,
        chain_homo=homo, cap_mode="hetero")
    assert report["single_block"] is True
    assert report["ok"] is False
    assert not moves._parallel_can_move_all_chains(
        idx, off, length, state.lattice.dimensions, True,
        chain_homo=homo, cap_mode="hetero")

    called = []

    def serial(*args):
        called.append("serial")
        return args[7], 0

    def parallel(*_args):
        raise AssertionError("single-block crank should not enter the parallel wrapper")

    monkeypatch.setattr(moves.mega_crank_fast, "mega_crank", serial)
    monkeypatch.setattr(moves.mega_crank_fast, "mega_crank_parallel", parallel)
    state.sim.MOVER.system_shake(
        state.lattice, state.energy, state.acc, state.ham, 3, "UNSET",
        parallelize=True, num_threads=4)
    assert called == ["serial"]


def test_frozen_chain_is_excluded_from_whole_chain_megamove(tmp_path, monkeypatch):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_SLITHER": 1.0},
                          box=[16, 16, 16], chains=[(2, "AABB")])
    original = [p[:] for p in state.lattice.chains[1].get_ordered_positions()]
    selectors = []

    def serial(*args):
        selectors.append(np.asarray(args[6]).copy())
        return args[11], 0

    monkeypatch.setattr(moves.mega_crank_fast, "mega_slither", serial)
    state.sim.MOVER.system_slither(
        state.lattice, state.energy, state.acc, state.ham, 2,
        frozen_chains=[1], parallelize=False)

    # sorted chain index 0 is frozen; only chain index 1 may be proposed.
    assert np.array_equal(selectors[0], np.array([1, 1], dtype=np.int32))
    assert state.lattice.chains[1].get_ordered_positions() == original


def test_report_splits_long_and_short_chains_between_the_two_kernels(tmp_path, monkeypatch):
    """A chain too long for a block interior goes to the serial side of the split.

    80-bead chains in a 32^3 SR box: the chain-level layout is 2 blocks of 16 with
    W=3, so the interior is 10 and an 80-mer can never fit one. It does NOT follow
    that the move runs serially - the report used to say the move "falls back to
    the SERIAL kernel", and the dispatch used to re-check the fit against the
    current configuration every megamove, which is what broke stationarity. The
    chains are now partitioned once, by length: the 80-mers are handed to the
    serial kernel and the 4-mers to the parallel one, and BOTH kernels run every
    megamove.
    """
    monkeypatch.chdir(tmp_path)     # the report is logged to log.txt in the working directory
    lay = mega_crank_fast.parallel_layout_info(32, 32, 32, False)
    assert lay["num_blocks"] > 1
    interior = lay["block_size"][0] - 2 * lay["W"]
    assert interior < 80

    state, _ = _build(tmp_path, 3, "SR", [32, 32, 32], [(3, "A" * 80), (5, "AABB")],
                      {"MOVE_CRANKSHAFT": 0.4, "MOVE_SLITHER": 0.3, "MOVE_PULL": 0.3}, {})
    out = "\n".join(state.sim.report_parallelization())

    # 5 four-mers fit an interior, the 3 eighty-mers cannot - and both moves say so
    split = ("parallel kernel for 5 chain(s) of length <= %d, serial kernel for 3 "
             "longer chain(s) - both run every megamove" % interior)
    assert out.count(split) == 2
    assert out.count("the split is fixed by chain length for the whole run") == 2
    # the old per-megamove-choice wording must not come back
    assert "falls back to the SERIAL kernel" not in out
    assert "span more than a block interior" not in out
    assert "the fit is re-checked every megamove" not in out


def test_report_lists_frozen_chains(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)     # the report is logged to log.txt in the working directory
    (tmp_path / "fz.in").write_text("C 1 2\n")
    state, _ = _build(tmp_path, 3, "SR", [32, 32, 32], [(10, "AABB")],
                      {"MOVE_CRANKSHAFT": 1.0}, {"FREEZE_FILE": "fz.in"})
    out = "\n".join(state.sim.report_parallelization())
    assert "Frozen chains: 2 (8 beads)" in out


def test_report_reissued_after_resized_equilibration(tmp_path):
    """The block decomposition depends on the box, so the report is re-issued when a
    resized equilibration swaps in the production box."""
    state = U.build_state(tmp_path, 3, "SR", True, {"MOVE_CRANKSHAFT": 1.0},
                          box=[36, 36, 36], chains=[(10, "AABB")], n_steps=12, equilibration=10,
                          extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 2,
                                 "RESIZED_EQUILIBRATION": "16 16 16",
                                 "ENERGY_CHECK": 0, "XTC_FREQ": 1000, "PRINT_FREQ": 1000})
    os.chdir(tmp_path)
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state.sim.run_simulation()
    log = (tmp_path / "log.txt").read_text()
    assert log.count("PARALLELIZATION REPORT") == 2
    assert "after resized equilibration (production box)" in log
    # both boxes split for the SR crankshaft halo (W=1: 16^3 -> 2 blocks/dim, 36^3 -> 4)
    assert log.count("block grid") == 2


def test_fit_report_reasons_match_gate(tmp_path):
    """parallel_chain_fit_report must agree with the boolean gate and explain it."""
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[32, 32, 32], chains=[(3, "A" * 80), (10, "AABB")])
    idx, chains, off, length, homo = moves.parallel_chain_metadata(state.lattice)
    rep = moves.parallel_chain_fit_report(idx, off, length, [32, 32, 32], False,
                                          chain_homo=homo, cap_mode="all")
    assert rep["ok"] is False and rep["n_too_extended"] >= 1 and rep["n_over_cap"] == 0
    assert rep["ok"] == moves._parallel_can_move_all_chains(idx, off, length, [32, 32, 32], False,
                                                            chain_homo=homo, cap_mode="all")
    # freezing the long chains removes the objection
    rep2 = moves.parallel_chain_fit_report(idx, off, length, [32, 32, 32], False,
                                           chain_homo=homo, cap_mode="all",
                                           frozen_chains=tuple(chains[:3]))
    assert rep2["ok"] is True and rep2["n_frozen"] == 3
    assert rep2["ok"] == moves._parallel_can_move_all_chains(idx, off, length, [32, 32, 32], False,
                                                             chain_homo=homo, cap_mode="all",
                                                             frozen_chains=tuple(chains[:3]))
