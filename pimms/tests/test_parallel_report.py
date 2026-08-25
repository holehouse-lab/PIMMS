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

from pimms import mega_crank_fast, moves
from pimms.tests import kernel_test_utils as U


def _build(tmp_path, dim, ff, box, chains, moveset, extra, hardwall=False):
    base = {"PARALLELIZE": "True", "PARALLEL_THREADS": 4}
    base.update(extra)
    state = U.build_state(tmp_path, dim, ff, hardwall, moveset, box=box, chains=chains,
                          extra=base)
    log = (tmp_path / "log.txt").read_text()
    return state, log


def test_openmp_info_shape():
    info = mega_crank_fast.openmp_info()
    assert set(info) == {"enabled", "max_threads"}
    assert isinstance(info["enabled"], bool) and info["max_threads"] >= 1


def test_report_absent_when_parallelize_off(tmp_path):
    U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                  extra={"PARALLELIZE": "False"})
    assert "PARALLELIZATION REPORT" not in (tmp_path / "log.txt").read_text()


def test_report_multiblock_crank_and_parallel_slither(tmp_path):
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
    # slither: box splits for the chain-level halo and short chains fit -> PARALLEL
    assert "Slither (MOVE_SLITHER): kernel mega_slither_parallel" in out
    assert "-> PARALLEL" in out
    assert "Pull (MOVE_PULL): not in the move set" in out
    assert "relax more slowly per step" in out


def test_report_single_block_warning_for_small_box(tmp_path):
    state, _ = _build(tmp_path, 3, "LR", [14, 14, 14], [(6, "AABB")],
                      {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5}, {})
    out = "\n".join(state.sim.report_parallelization())
    assert "ONE block -> runs single-threaded, equivalent to the serial kernel" in out
    assert "Slither (MOVE_SLITHER): box does not split for the chain-level halo" in out
    assert "runs on the SERIAL kernel" in out


def test_report_serial_fallback_for_over_extended_chain(tmp_path):
    # 80-bead chains in a 32^3 SR box: the chain-level layout is 2 blocks of 16 with
    # W=3 -> interior 10, which an 80-mer cannot fit on every axis -> slither and
    # pull both fall back to the serial kernel
    lay = mega_crank_fast.parallel_layout_info(32, 32, 32, False)
    assert lay["num_blocks"] > 1
    state, _ = _build(tmp_path, 3, "SR", [32, 32, 32], [(3, "A" * 80)],
                      {"MOVE_CRANKSHAFT": 0.4, "MOVE_SLITHER": 0.3, "MOVE_PULL": 0.3}, {})
    out = "\n".join(state.sim.report_parallelization())
    assert out.count("falls back to the SERIAL kernel") == 2
    assert "span more than a block interior" in out


def test_report_lists_frozen_chains(tmp_path):
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
