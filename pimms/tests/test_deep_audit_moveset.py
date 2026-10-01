"""
Deep-audit regression tests for the start-up checks and reports about the move set.

* M2-3: the null-move warning counts every chain a move cannot act on (not just
  monomers), and its percentages match what a real run measures.
* M2-4: mobile chains too short for every enabled move get a warning.
* M2-5 / M1-1: a rigid-only move set freezes a dimer's bond length, and says so.
* M2-6 / M1-2: a cluster-only move set freezes the contact graph (warning), and is
  refused when it is provable at start-up that no chain can ever move.
* PERF5: the parallelization report warns when a megamove is too small to pay for
  the parallel kernel's fixed cost.
* P2-5: a discarded Simulation is freed by reference counting alone.

Nothing here compares the code against itself. Null fractions are checked against
draws counted in a real run, immobile chains against positions before and after a
run, clusters against a flood fill written here, and the parallel-kernel rule
against the arithmetic stated in the report's docstring.
"""

from __future__ import annotations

import contextlib
import gc
import itertools
import math
import os
import pathlib
import re
import weakref
from typing import Any

import numpy as np
import pytest

from pimms import mega_crank_fast, simulation
from pimms.latticeExceptions import AnalysisRoutineException, SimulationException
from pimms.simulation import Simulation, conformation_freezing_warnings
from pimms.tests import kernel_test_utils as U

QUIET = {
    "PRINT_FREQ": 10**8,
    "XTC_FREQ": 10**8,
    "ANALYSIS_FREQ": 10**8,
    "RESTART_FREQ": 10**8,
    "EN_FREQ": 10**8,
    "ENERGY_CHECK": 0,
    "CRANKSHAFT_SUBSTEPS": 20,
}
TSMMC = {
    "TSMMC_JUMP_TEMP": 120,
    "TSMMC_STEP_MULTIPLIER": 2,
    "TSMMC_NUMBER_OF_POINTS": 2,
}

NULL_CASES = {
    # name: (chains, frozen, moves, expected {keyword: fraction of steps})
    "dimers_pivot": (
        [(6, "AB")],
        None,
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_CHAIN_PIVOT": 0.5},
        {"MOVE_CHAIN_PIVOT": 0.5},
    ),
    "monomers_and_dimers_pivot": (
        [(5, "A"), (5, "AB")],
        None,
        {"MOVE_CRANKSHAFT": 0.6, "MOVE_CHAIN_PIVOT": 0.4},
        {"MOVE_CHAIN_PIVOT": 0.4},
    ),
    # 4 monomers, 3 dimers, 3 pentamers: rotate and head pivot miss 4 of 10, pivot 7 of 10
    "three_lengths": (
        [(4, "A"), (3, "AB"), (3, "AABBA")],
        None,
        {
            "MOVE_CRANKSHAFT": 0.4,
            "MOVE_CHAIN_ROTATE": 0.2,
            "MOVE_CHAIN_PIVOT": 0.2,
            "MOVE_HEAD_PIVOT": 0.2,
        },
        {
            "MOVE_CHAIN_ROTATE": 0.2 * 4 / 10,
            "MOVE_CHAIN_PIVOT": 0.2 * 7 / 10,
            "MOVE_HEAD_PIVOT": 0.2 * 4 / 10,
        },
    ),
    # the same system with the pentamers (chains 8-10) frozen: 7 mobile chains
    "three_lengths_frozen": (
        [(4, "A"), (3, "AB"), (3, "AABBA")],
        [8, 9, 10],
        {
            "MOVE_CRANKSHAFT": 0.4,
            "MOVE_CHAIN_ROTATE": 0.2,
            "MOVE_CHAIN_PIVOT": 0.2,
            "MOVE_HEAD_PIVOT": 0.2,
        },
        {
            "MOVE_CHAIN_ROTATE": 0.2 * 4 / 7,
            "MOVE_CHAIN_PIVOT": 0.2,
            "MOVE_HEAD_PIVOT": 0.2 * 4 / 7,
        },
    ),
    "dimers_pull": (
        [(6, "AB")],
        None,
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_PULL": 0.5},
        {"MOVE_PULL": 0.5},
    ),
}

CLUSTER_ONLY = {"MOVE_CLUSTER_TRANSLATE": 0.5, "MOVE_CLUSTER_ROTATE": 0.5}


@pytest.fixture(autouse=True)
def _restore_cwd(tmp_path: pathlib.Path) -> Any:
    """Run each test in its own temporary directory and put the old one back.

    The start-up checks log their warnings to ``log.txt`` in the working
    directory, so a test left in the directory pytest was started from appends
    to (or creates) a ``log.txt`` there.

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytest's per-test temporary directory.

    Returns
    -------
    Any
        Generator fixture; nothing is yielded.
    """
    cwd = os.getcwd()
    os.chdir(tmp_path)
    yield
    os.chdir(cwd)


def _build(
    tmp_path: pathlib.Path,
    dim: int,
    moves: dict[str, float],
    chains: list[tuple[int, str]],
    box: list[int],
    frozen: list[int] | None = None,
    extra: dict[str, Any] | None = None,
    hardwall: bool = False,
    n_steps: int = 10,
    bypass: bool = False,
) -> U.State:
    """Build a Simulation in ``tmp_path``.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write the keyfile, parameter file and output into.
    dim : int
        2 or 3.
    moves : dict of str to float
        ``MOVE_*`` fractions.
    chains : list of (int, str)
        ``(count, sequence)`` pairs.
    box : list of int
        ``DIMENSIONS``.
    frozen : list of int, optional
        chainIDs to freeze through a ``FREEZE_FILE``.
    extra : dict, optional
        Further keywords.
    hardwall : bool, optional
        ``HARDWALL``.
    n_steps : int, optional
        ``N_STEPS``.
    bypass : bool, optional
        Skip ``check_moveset_applicability`` during construction, so that a
        move set the check would refuse can still be built and run.

    Returns
    -------
    kernel_test_utils.State
        The built state; ``state.sim`` is the Simulation.
    """
    keywords = dict(QUIET)
    keywords.update(extra or {})
    if frozen:
        (tmp_path / "freeze.in").write_text(
            "C " + " ".join(str(c) for c in frozen) + "\n"
        )
        keywords["FREEZE_FILE"] = "freeze.in"
    original = Simulation.check_moveset_applicability
    if bypass:
        Simulation.check_moveset_applicability = lambda self, keyword_lookup: None
    try:
        return U.build_state(
            tmp_path,
            dim,
            "SR",
            hardwall,
            moves,
            box=box,
            chains=chains,
            n_steps=n_steps,
            equilibration=1,
            extra=keywords,
        )
    finally:
        Simulation.check_moveset_applicability = original


def _check(sim: Simulation, moves: dict[str, float] | None = None) -> dict[str, Any]:
    """Run the applicability check and return its summary.

    Parameters
    ----------
    sim : Simulation
        A built simulation.
    moves : dict of str to float, optional
        ``MOVE_*`` fractions to check instead of the ones the simulation was
        built with (every other move is set to zero).

    Returns
    -------
    dict
        The summary ``check_moveset_applicability`` returns, plus ``"messages"``:
        the warnings it logged, captured here from the logger and not from the
        summary. A test asserts on the messages first, so that on code which
        returns no summary it is the message check that fails.
    """
    keyword_lookup = dict(sim.keyword_lookup)
    if moves is not None:
        for keyword in keyword_lookup:
            if keyword.startswith("MOVE_"):
                keyword_lookup[keyword] = float(moves.get(keyword, 0.0))
    messages: list[str] = []
    real = simulation.pimmslogger.log_warning
    simulation.pimmslogger.log_warning = lambda message: (
        messages.append(message),
        real(message),
    )[1]
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            summary = sim.check_moveset_applicability(keyword_lookup)
    finally:
        simulation.pimmslogger.log_warning = real
    checked = dict(summary or {})
    checked["messages"] = messages
    return checked


def _run(sim: Simulation, directory: pathlib.Path) -> None:
    """Run a built simulation with output and analysis switched off.

    Parameters
    ----------
    sim : Simulation
        The simulation to run.
    directory : pathlib.Path
        Its working directory.

    Returns
    -------
    None
    """
    sim.simulation_IO = lambda step, energy: None
    sim.run_all_analysis = lambda step: None
    os.chdir(directory)
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        sim.run_simulation()


def _positions(sim: Simulation) -> dict[int, list[tuple[int, ...]]]:
    """Snapshot every chain's ordered bead positions.

    Parameters
    ----------
    sim : Simulation
        The simulation to read.

    Returns
    -------
    dict of int to list of tuple
        chainID -> bead positions.
    """
    return {
        cid: [tuple(int(x) for x in p) for p in chain.get_ordered_positions()]
        for cid, chain in sim.LATTICE.chains.items()
    }


def _clusters(
    positions: dict[int, list[tuple[int, ...]]], dims: list[int]
) -> dict[int, frozenset]:
    """Clusters of chains by flood fill, written independently of PIMMS.

    Two chains touch when two of their beads are within Chebyshev distance 1
    under the minimum image (the boxes used here are periodic).

    Parameters
    ----------
    positions : dict of int to list of tuple
        chainID -> bead positions.
    dims : list of int
        Box dimensions.

    Returns
    -------
    dict of int to frozenset
        chainID -> the set of chainIDs in its cluster.
    """
    d = np.asarray(dims)
    touch = {c: set() for c in positions}
    for a, b in itertools.combinations(sorted(positions), 2):
        diff = np.abs(
            np.asarray(positions[a])[:, None, :] - np.asarray(positions[b])[None, :, :]
        )
        diff = np.minimum(diff, d - diff)
        if (diff.max(axis=2) <= 1).any():
            touch[a].add(b)
            touch[b].add(a)
    cluster: dict[int, frozenset] = {}
    for c in positions:
        if c in cluster:
            continue
        seen, stack = {c}, [c]
        while stack:
            for y in touch[stack.pop()]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        for x in seen:
            cluster[x] = frozenset(seen)
    return cluster


# ---------------------------------------------------------------------------
# M2-3: the null-move fraction is right, and checkable
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", sorted(NULL_CASES))
def test_null_move_fraction_matches_a_real_run(
    tmp_path: pathlib.Path, case: str
) -> None:
    """The fractions in the warning are the fractions a run actually wastes.

    Every per-chain move call is recorded with the length of the chain it was
    drawn for and whether the move was possible. A draw is counted as null when
    it landed on a chain length for which that move never once succeeded in the
    run; pull megamoves are counted as null when they propose nothing.
    """
    chains, frozen, moves, expected = NULL_CASES[case]
    n_steps = 2000
    state = _build(
        tmp_path, 3, moves, chains, [12, 12, 12], frozen=frozen, n_steps=n_steps
    )
    sim = state.sim
    # the message, as logged at construction. The prediction is from the
    # definition: MOVE_x times the fraction of mobile chains the move cannot act on
    log = (tmp_path / "log.txt").read_text()
    total = re.search(r"About (\d+)% of steps will be rejected null moves", log)
    assert total is not None, "no null-move warning was issued"
    assert int(total.group(1)) == round(100 * sum(expected.values()))
    for keyword, fraction in expected.items():
        clause = re.search(keyword + r"[^;]*?([\d.]+)% of steps", log)
        assert clause is not None, "%s is not named in the warning" % keyword
        assert float(clause.group(1)) == pytest.approx(100 * fraction, abs=0.05)

    # the returned summary carries the same numbers
    summary = _check(sim)
    assert set(summary["null_fractions"]) == set(expected)
    for keyword, fraction in expected.items():
        assert summary["null_fractions"][keyword] == pytest.approx(fraction, abs=1e-12)
    assert summary["null_fraction"] == pytest.approx(sum(expected.values()), abs=1e-12)

    # the measurement
    calls: dict[str, list[tuple[int, bool]]] = {
        k: [] for k in ("MOVE_CHAIN_ROTATE", "MOVE_CHAIN_PIVOT", "MOVE_HEAD_PIVOT")
    }
    for keyword, name in (
        ("MOVE_CHAIN_ROTATE", "chain_rotate"),
        ("MOVE_CHAIN_PIVOT", "chain_pivot"),
        ("MOVE_HEAD_PIVOT", "head_pivot"),
    ):
        real = getattr(sim.MOVER, name)

        def recorder(chain, *args, _real=real, _keyword=keyword, **kwargs):
            result = _real(chain, *args, **kwargs)
            calls[_keyword].append(
                (len(chain.get_ordered_positions()), bool(result[1]))
            )
            return result

        setattr(sim.MOVER, name, recorder)
    pulls: list[int] = []
    real_pull = sim.MOVER.system_pull

    def pull_recorder(*args, **kwargs):
        result = real_pull(*args, **kwargs)
        pulls.append(int(result[2]))
        return result

    sim.MOVER.system_pull = pull_recorder
    _run(sim, tmp_path)

    need = {"MOVE_CHAIN_ROTATE": 2, "MOVE_CHAIN_PIVOT": 3, "MOVE_HEAD_PIVOT": 2}
    for keyword, fraction in expected.items():
        if keyword == "MOVE_PULL":
            measured = sum(1 for proposed in pulls if proposed == 0) / n_steps
            assert pulls and all(proposed == 0 for proposed in pulls)
        else:
            lengths = {length for length, _ in calls[keyword]}
            never = {
                length
                for length in lengths
                if not any(ok for drawn, ok in calls[keyword] if drawn == length)
            }
            # the move is impossible for exactly the lengths below its minimum
            assert never == {length for length in lengths if length < need[keyword]}
            measured = (
                sum(1 for length, _ in calls[keyword] if length in never) / n_steps
            )
        sigma = math.sqrt(fraction * (1 - fraction) / n_steps)
        assert abs(measured - fraction) < 4 * sigma + 1e-9, (
            keyword,
            measured,
            fraction,
        )


def test_null_fraction_inside_system_tsmmc_excursions(tmp_path: pathlib.Path) -> None:
    """An excursion draws its sub-moves from the non-TSMMC moves, renormalised.

    Dimers with crankshaft 0.5, pivot 0.25 and system TSMMC 0.25: a quarter of
    the main steps are null pivots, and a third (0.25 / 0.75) of the sub-moves
    inside an excursion. The excursion itself is not null - its crankshaft
    sub-moves act on dimers - which is what taking the SHORTEST sub-move minimum
    length means.
    """
    moves = {
        "MOVE_CRANKSHAFT": 0.5,
        "MOVE_SYSTEM_TSMMC": 0.25,
        "MOVE_CHAIN_PIVOT": 0.25,
    }
    n_steps = 400
    state = _build(
        tmp_path, 3, moves, [(6, "AB")], [9, 9, 9], extra=TSMMC, n_steps=n_steps
    )
    sim = state.sim
    assert "33.3% of them are null too" in (tmp_path / "log.txt").read_text()
    summary = _check(sim)
    assert summary["null_fractions"] == {"MOVE_CHAIN_PIVOT": pytest.approx(0.25)}
    assert summary["excursion_null_fraction"] == pytest.approx(1.0 / 3.0)

    counts = {"main_pivot": 0, "sub_pivot": 0, "sub_crank": 0}
    real_pivot, real_shake = sim.MOVER.chain_pivot, sim.MOVER.system_shake

    def pivot(*args, **kwargs):
        counts["sub_pivot" if sim.auxillary_chain else "main_pivot"] += 1
        result = real_pivot(*args, **kwargs)
        assert result[1] is False
        return result

    def shake(*args, **kwargs):
        if sim.auxillary_chain:
            counts["sub_crank"] += 1
        return real_shake(*args, **kwargs)

    sim.MOVER.chain_pivot, sim.MOVER.system_shake = pivot, shake
    _run(sim, tmp_path)

    main = counts["main_pivot"] / n_steps
    assert abs(main - 0.25) < 4 * math.sqrt(0.25 * 0.75 / n_steps)
    n_sub = counts["sub_pivot"] + counts["sub_crank"]
    assert n_sub > 1000
    assert abs(counts["sub_pivot"] / n_sub - 1 / 3) < 4 * math.sqrt(
        (1 / 3) * (2 / 3) / n_sub
    )


def test_no_null_warning_when_every_move_fits_every_chain(
    tmp_path: pathlib.Path,
) -> None:
    state = _build(
        tmp_path,
        3,
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_CHAIN_PIVOT": 0.25, "MOVE_PULL": 0.25},
        [(4, "AAB"), (2, "AABBA")],
        [12, 12, 12],
    )
    assert "rejected null moves" not in (tmp_path / "log.txt").read_text()
    summary = _check(state.sim)
    assert summary["messages"] == []
    # a stay-quiet test: it holds before and after the fix, so the summary is
    # read leniently and no message check is hidden behind a missing key
    assert (
        summary.get("null_fraction", 0.0) == 0.0
        and summary.get("null_fractions", {}) == {}
    )


# ---------------------------------------------------------------------------
# M2-4: chains no enabled move can act on
# ---------------------------------------------------------------------------


def test_chains_too_short_for_every_move_are_reported_and_really_never_move(
    tmp_path: pathlib.Path,
) -> None:
    """Pull alone on monomers, dimers and pentamers: only the pentamers can move.

    The run is not refused (the longest mobile chain is what decides that), the
    warning names the five short chains, and a run confirms they stay put.
    """
    state = _build(
        tmp_path,
        3,
        {"MOVE_PULL": 1.0},
        [(3, "A"), (2, "AB"), (2, "AABBA")],
        [10, 10, 10],
        n_steps=200,
    )
    sim = state.sim
    log = (tmp_path / "log.txt").read_text()
    assert (
        "No enabled move can act on 5 of 7 mobile chains (3 of 1 bead, 2 of 2 beads)"
        in log
    )
    assert "MOVE_PULL needs 3 or more beads" in log
    assert _check(sim)["unmoved_lengths"] == [1, 2]

    before = _positions(sim)
    _run(sim, tmp_path)
    after = _positions(sim)
    moved = {len(before[c]) for c in before if before[c] != after[c]}
    assert moved == {5}


def test_no_unmoved_warning_when_a_move_reaches_the_short_chains(
    tmp_path: pathlib.Path,
) -> None:
    state = _build(
        tmp_path,
        3,
        {"MOVE_PULL": 0.5, "MOVE_CHAIN_TRANSLATE": 0.5},
        [(3, "A"), (2, "AB"), (2, "AABBA")],
        [10, 10, 10],
    )
    assert "No enabled move can act on" not in (tmp_path / "log.txt").read_text()
    # a stay-quiet test (true before and after the fix): read the summary leniently
    assert _check(state.sim).get("unmoved_lengths", []) == []


def test_excursion_sub_moves_count_for_the_short_chains(tmp_path: pathlib.Path) -> None:
    """A system excursion with no other move falls back to crankshaft sub-moves,
    which reach every chain - so pairing it with pull leaves nothing unmoved."""
    state = _build(
        tmp_path,
        3,
        {"MOVE_SYSTEM_TSMMC": 1.0},
        [(3, "A"), (2, "AABBA")],
        [10, 10, 10],
        extra=TSMMC,
    )
    alone = _check(state.sim)
    assert not any("No enabled move can act on" in m for m in alone["messages"])
    # with pull as the only sub-move the excursion can do no more than pull can
    summary = _check(state.sim, {"MOVE_SYSTEM_TSMMC": 0.5, "MOVE_PULL": 0.5})
    assert any(
        "No enabled move can act on 3 of 5 mobile chains (3 of 1 bead)" in m
        for m in summary["messages"]
    )
    assert alone["unmoved_lengths"] == []
    assert summary["unmoved_lengths"] == [1]


# ---------------------------------------------------------------------------
# M2-5 / M1-1: a dimer's bond is frozen by a rigid-only move set
# ---------------------------------------------------------------------------


def _bond_lengths_squared(sim: Simulation) -> dict[int, int]:
    """Squared minimum-image bond length of every two-bead chain.

    Parameters
    ----------
    sim : Simulation
        The simulation to read.

    Returns
    -------
    dict of int to int
        chainID -> squared bond length (1, 2 or 3 on a 3D lattice).
    """
    dims = np.asarray(sim.LATTICE.dimensions)
    out = {}
    for cid, beads in _positions(sim).items():
        bond = np.abs(np.asarray(beads[1]) - np.asarray(beads[0]))
        bond = np.minimum(bond, dims - bond)
        out[cid] = int((bond**2).sum())
    return out


def test_rigid_only_moves_on_dimers_warn_and_the_bonds_really_are_frozen(
    tmp_path: pathlib.Path,
) -> None:
    state = _build(
        tmp_path,
        3,
        {"MOVE_CHAIN_TRANSLATE": 0.5, "MOVE_CHAIN_ROTATE": 0.5},
        [(8, "AB")],
        [10, 10, 10],
        n_steps=400,
    )
    log = (tmp_path / "log.txt").read_text()
    assert "contains no move that can change a chain's shape" in log
    assert "bond length" in log
    before, bonds = _positions(state.sim), _bond_lengths_squared(state.sim)
    _run(state.sim, tmp_path)
    assert _positions(state.sim) != before  # the chains did move
    assert _bond_lengths_squared(state.sim) == bonds  # and no bond changed


def test_head_pivot_on_dimers_does_not_warn_and_the_bonds_change(
    tmp_path: pathlib.Path,
) -> None:
    state = _build(
        tmp_path, 3, {"MOVE_HEAD_PIVOT": 1.0}, [(8, "AB")], [10, 10, 10], n_steps=400
    )
    assert (
        "contains no move that can change a chain's shape"
        not in (tmp_path / "log.txt").read_text()
    )
    bonds = _bond_lengths_squared(state.sim)
    _run(state.sim, tmp_path)
    assert _bond_lengths_squared(state.sim) != bonds


def test_conformation_freezing_warnings_thresholds_and_wording() -> None:
    rigid = ["MOVE_CHAIN_TRANSLATE", "MOVE_VMMC"]
    assert conformation_freezing_warnings(rigid, 1) == []
    assert len(conformation_freezing_warnings(rigid, 2)) == 1
    assert len(conformation_freezing_warnings(rigid, 3)) == 1
    # a dimer has no interior or midpoint bead: the narrower warnings start at three
    assert conformation_freezing_warnings(["MOVE_HEAD_PIVOT"], 2) == []
    assert len(conformation_freezing_warnings(["MOVE_HEAD_PIVOT"], 3)) == 1
    assert conformation_freezing_warnings(["MOVE_CRANKSHAFT"] + rigid, 2) == []
    # the rigid moves named in the pull and pivot warnings are the non-conformational ones
    (pull,) = conformation_freezing_warnings(["MOVE_PULL"] + rigid, 5)
    assert "The rigid moves (MOVE_CHAIN_TRANSLATE, MOVE_VMMC) do not release it" in pull
    (pull_alone,) = conformation_freezing_warnings(["MOVE_PULL"], 5)
    assert "The rigid moves (none enabled)" in pull_alone
    (pivots,) = conformation_freezing_warnings(
        ["MOVE_CHAIN_PIVOT", "MOVE_CHAIN_ROTATE"], 5
    )
    assert "The rigid moves (MOVE_CHAIN_ROTATE) still carry chains" in pivots


def test_longest_mobile_chain_decides_the_refusal(tmp_path: pathlib.Path) -> None:
    """Monomers plus pentamers with only a pivot: the pentamers pivot, so no refusal."""
    state = _build(
        tmp_path, 3, {"MOVE_CHAIN_PIVOT": 1.0}, [(4, "A"), (2, "AABBA")], [10, 10, 10]
    )
    summary = _check(state.sim)
    assert any(
        "MOVE_CHAIN_PIVOT cannot move a chain of fewer than 3 beads: 4 of 6 mobile chains"
        in m
        for m in summary["messages"]
    )
    assert summary["null_fractions"] == {"MOVE_CHAIN_PIVOT": pytest.approx(4 / 6)}
    assert summary["unmoved_lengths"] == [1]
    # freeze the pentamers and the same move set has nothing left to act on
    refused = tmp_path / "refused"
    refused.mkdir()
    with pytest.raises(SimulationException, match="No enabled move can act"):
        _build(
            refused,
            3,
            {"MOVE_CHAIN_PIVOT": 1.0},
            [(4, "A"), (2, "AABBA")],
            [10, 10, 10],
            frozen=[5, 6],
        )


# ---------------------------------------------------------------------------
# M2-6 / M1-2: cluster-only move sets
# ---------------------------------------------------------------------------


def test_cluster_only_move_set_warns_that_the_contact_graph_is_frozen(
    tmp_path: pathlib.Path,
) -> None:
    """Cluster moves reject merges and nothing splits a cluster, so which chains
    touch which is a constant of the run. Checked on a run: the clusters found
    by an independent flood fill are the same before and after."""
    state = _build(
        tmp_path, 2, CLUSTER_ONLY, [(6, "A"), (3, "AB")], [14, 14], n_steps=600
    )
    sim = state.sim
    assert "contact graph is frozen" in (tmp_path / "log.txt").read_text()
    before = _positions(sim)
    clusters = _clusters(before, sim.LATTICE.dimensions)
    _run(sim, tmp_path)
    after = _positions(sim)
    assert after != before
    assert _clusters(after, sim.LATTICE.dimensions) == clusters


def test_no_contact_graph_warning_when_another_move_is_enabled(
    tmp_path: pathlib.Path,
) -> None:
    _build(
        tmp_path,
        2,
        {"MOVE_CLUSTER_TRANSLATE": 0.5, "MOVE_CRANKSHAFT": 0.5},
        [(6, "A")],
        [14, 14],
    )
    assert "contact graph is frozen" not in (tmp_path / "log.txt").read_text()


def test_isolated_monomers_with_only_cluster_rotation_are_refused(
    tmp_path: pathlib.Path,
) -> None:
    """Every rotation maps an isolated bead onto itself, and the beads can never
    meet, so nothing can move. Translation, or one touching pair, changes that."""
    state = _build(
        tmp_path, 3, {"MOVE_CRANKSHAFT": 1.0}, [(6, "A")], [15, 15, 15], n_steps=300
    )
    sim = state.sim
    clusters = _clusters(_positions(sim), sim.LATTICE.dimensions)
    assert all(len(members) == 1 for members in clusters.values()), (
        "fixture: monomers must start isolated"
    )

    with pytest.raises(
        SimulationException, match="No enabled move can ever change this configuration"
    ):
        _check(sim, {"MOVE_CLUSTER_ROTATE": 1.0})
    with pytest.raises(SimulationException, match="isolated single bead"):
        _check(sim, {"MOVE_CLUSTER_ROTATE": 0.5, "MOVE_SYSTEM_TSMMC": 0.5})
    # translation moves an isolated bead; so does any single-chain move
    assert (
        _check(sim, {"MOVE_CLUSTER_ROTATE": 0.5, "MOVE_CLUSTER_TRANSLATE": 0.5})[
            "immobile_chains"
        ]
        == []
    )
    assert (
        _check(sim, {"MOVE_CLUSTER_ROTATE": 0.5, "MOVE_CRANKSHAFT": 0.5})[
            "immobile_chains"
        ]
        == []
    )

    # the refusal is true: run the refused move set and nothing moves
    rotate_only = tmp_path / "rotate_only"
    rotate_only.mkdir()
    refused = _build(
        rotate_only,
        3,
        {"MOVE_CLUSTER_ROTATE": 1.0},
        [(6, "A")],
        [15, 15, 15],
        n_steps=300,
        bypass=True,
    ).sim
    before = _positions(refused)
    assert all(
        len(m) == 1 for m in _clusters(before, refused.LATTICE.dimensions).values()
    )
    _run(refused, rotate_only)
    assert _positions(refused) == before


def test_touching_monomers_with_only_cluster_rotation_are_warned_not_refused(
    tmp_path: pathlib.Path,
) -> None:
    """Six monomers in a 14 x 14 box, two of which start in contact: that pair
    rotates, so the run is legitimate; the four isolated beads are reported."""
    state = _build(tmp_path, 2, {"MOVE_CRANKSHAFT": 1.0}, [(6, "A")], [14, 14])
    sim = state.sim
    clusters = _clusters(_positions(sim), sim.LATTICE.dimensions)
    isolated = sorted(c for c, members in clusters.items() if len(members) == 1)
    assert 0 < len(isolated) < 6, "fixture: needs both isolated and touching monomers"
    summary = _check(sim, {"MOVE_CLUSTER_ROTATE": 1.0})
    assert any(
        "%d of 6 mobile chain(s) can never move" % len(isolated) in w
        for w in summary["messages"]
    )
    assert summary["immobile_chains"] == isolated


def test_cluster_only_with_frozen_neighbours(tmp_path: pathlib.Path) -> None:
    """A mobile chain in a cluster with a frozen chain is never moved by a cluster
    move. If that is every mobile chain the run is refused; a single-chain move
    lifts the refusal; and a resized equilibration turns it into a warning."""
    state = _build(
        tmp_path,
        3,
        {"MOVE_CRANKSHAFT": 1.0},
        [(4, "AABBA")],
        [9, 9, 9],
        frozen=[1, 2, 3],
    )
    sim = state.sim
    clusters = _clusters(_positions(sim), sim.LATTICE.dimensions)
    assert clusters[4] & {1, 2, 3}, (
        "fixture: the mobile chain must touch a frozen chain"
    )

    with pytest.raises(
        SimulationException, match=r"contains a frozen chain \(FREEZE_FILE\)"
    ):
        _check(sim, CLUSTER_ONLY)
    summary = _check(sim, {"MOVE_CLUSTER_TRANSLATE": 0.5, "MOVE_CRANKSHAFT": 0.5})
    assert summary["messages"] == [] and summary["immobile_chains"] == []

    sim.resize_eq = True
    try:
        summary = _check(sim, CLUSTER_ONLY)
    finally:
        sim.resize_eq = False
    assert any(
        "1 of 1 mobile chain(s) can never move" in m for m in summary["messages"]
    )
    assert summary["immobile_chains"] == [4]

    # the refusal is true: with the check bypassed the mobile chain never moves
    stuck = tmp_path / "stuck"
    stuck.mkdir()
    run = _build(
        stuck,
        3,
        CLUSTER_ONLY,
        [(4, "AABBA")],
        [9, 9, 9],
        frozen=[1, 2, 3],
        n_steps=300,
        bypass=True,
    ).sim
    before = _positions(run)
    assert _clusters(before, run.LATTICE.dimensions)[4] & {1, 2, 3}
    _run(run, stuck)
    assert _positions(run) == before


def test_cluster_only_is_not_refused_when_the_mobile_chain_is_free(
    tmp_path: pathlib.Path,
) -> None:
    """The same frozen system in a box large enough that the mobile chain starts
    alone: the cluster moves do move it, so there is a warning and no refusal."""
    state = _build(
        tmp_path,
        3,
        CLUSTER_ONLY,
        [(1, "AABBA"), (1, "AB")],
        [16, 16, 16],
        frozen=[1],
        n_steps=300,
    )
    sim = state.sim
    before = _positions(sim)
    assert _clusters(before, sim.LATTICE.dimensions)[2] == frozenset({2}), (
        "fixture: chain 2 must start alone"
    )
    assert "contact graph is frozen" in (tmp_path / "log.txt").read_text()
    _run(sim, tmp_path)
    after = _positions(sim)
    assert after[1] == before[1] and after[2] != before[2]


def test_one_cluster_holding_every_chain_is_refused_but_a_single_chain_is_not(
    tmp_path: pathlib.Path,
) -> None:
    """The cluster moves skip a cluster that holds every chain of a multi-chain
    system. A one-chain system is the exception - the size check never fires -
    and its chain does move."""
    state = _build(tmp_path, 2, {"MOVE_CRANKSHAFT": 1.0}, [(3, "AABBAABBAA")], [7, 7])
    sim = state.sim
    assert len(_clusters(_positions(sim), sim.LATTICE.dimensions)[1]) == 3, (
        "fixture: the chains must touch"
    )
    with pytest.raises(SimulationException, match="holds every chain in the system"):
        _check(sim, CLUSTER_ONLY)

    single = tmp_path / "single"
    single.mkdir()
    one = _build(
        single,
        2,
        {"MOVE_CLUSTER_TRANSLATE": 1.0},
        [(1, "AABBA")],
        [12, 12],
        n_steps=100,
    ).sim
    before = _positions(one)
    _run(one, single)
    assert _positions(one) != before


# ---------------------------------------------------------------------------
# PERF5: megamoves too small for the parallel kernel
# ---------------------------------------------------------------------------


def _parallel_report(
    tmp_path: pathlib.Path,
    name: str,
    moves: dict[str, float],
    chains: list[tuple[int, str]],
    box: list[int],
    extra: dict[str, Any],
) -> tuple[Simulation, str]:
    """Build a PARALLELIZE run and return it with its report as one string.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Parent directory; the run goes in ``tmp_path / name``.
    name : str
        Sub-directory name.
    moves : dict of str to float
        ``MOVE_*`` fractions.
    chains : list of (int, str)
        ``(count, sequence)`` pairs.
    box : list of int
        ``DIMENSIONS``.
    extra : dict
        Further keywords (``PARALLEL_THREADS`` and the substep counts).

    Returns
    -------
    tuple of (Simulation, str)
        The simulation and the report lines joined with newlines.
    """
    directory = tmp_path / name
    directory.mkdir()
    keywords = {"PARALLELIZE": "True", "PARALLEL_THREADS": 4}
    keywords.update(extra)
    sim = _build(directory, len(box), moves, chains, box, extra=keywords).sim
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        return sim, "\n".join(sim.report_parallelization())


def test_parallel_report_warns_when_the_megamove_cannot_pay_for_the_kernel(
    tmp_path: pathlib.Path,
) -> None:
    """The rule, from the docstring: fixed cost 0.1 ms + 0.1 us per bead; a
    crankshaft sub-move costs 0.1 us; with T threads the most a megamove can
    save is substeps x cost x (1 - 1/T). Below break-even the report warns."""
    if not mega_crank_fast.openmp_info()["enabled"]:
        pytest.skip("kernels built without OpenMP")
    chains, box = [(250, "AABB")], [32, 32, 32]
    n_beads = 1000
    fixed_ms = 0.1 + 1.0e-4 * n_beads

    def break_even(threads: int) -> float:
        return fixed_ms / (1.0e-4 * (1.0 - 1.0 / threads))

    assert (
        mega_crank_fast.parallel_crank_layout_info(32, 32, 32, False)["num_blocks"] >= 4
    )
    be4, be2 = break_even(4), break_even(2)
    assert be4 < be2

    # the keyword default (500) is far below break-even
    sim, default = _parallel_report(
        tmp_path,
        "default",
        {"MOVE_CRANKSHAFT": 1.0},
        chains,
        box,
        {"CRANKSHAFT_SUBSTEPS": 500},
    )
    assert (
        "megamove too small for the parallel kernel - at CRANKSHAFT_SUBSTEPS : 500, 4 threads "
        "can save at most about %.2g ms per megamove (75%% of %.2g ms of sampling)"
        % (500 * 1.0e-4 * 0.75, 500 * 1.0e-4)
        in default
    )
    assert "approximate" in default
    recommended = int(
        re.search(r"Raise CRANKSHAFT_SUBSTEPS to about (\d+) or more", default).group(1)
    )
    assert 10 * fixed_ms / 1.0e-4 <= recommended <= 2 * 10 * fixed_ms / 1.0e-4
    # it is a warning in the log, inside the report
    log = (tmp_path / "default" / "log.txt").read_text()
    warning = [line for line in log.splitlines() if "megamove too small" in line]
    assert len(warning) == 1 and "WARNING" in warning[0].upper()
    # run constants only: the same lines after the configuration has changed
    _run(sim, tmp_path / "default")
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        assert "\n".join(sim.report_parallelization()) == default

    # either side of break-even, and the dependence on the thread count
    _, below = _parallel_report(
        tmp_path,
        "below",
        {"MOVE_CRANKSHAFT": 1.0},
        chains,
        box,
        {"CRANKSHAFT_SUBSTEPS": int(be4 * 0.9)},
    )
    _, above = _parallel_report(
        tmp_path,
        "above",
        {"MOVE_CRANKSHAFT": 1.0},
        chains,
        box,
        {"CRANKSHAFT_SUBSTEPS": int(be4 * 1.1)},
    )
    _, two = _parallel_report(
        tmp_path,
        "two",
        {"MOVE_CRANKSHAFT": 1.0},
        chains,
        box,
        {"CRANKSHAFT_SUBSTEPS": int(be4 * 1.1), "PARALLEL_THREADS": 2},
    )
    assert "megamove too small" in below
    assert "megamove too small" not in above
    assert int(be4 * 1.1) < be2 and "megamove too small" in two


def test_parallel_report_megamove_rule_for_slither_counts_the_parallel_chains(
    tmp_path: pathlib.Path,
) -> None:
    """Slither's parallel budget is SLITHER_SUBSTEPS x (chains on the parallel
    side), at about 0.4 us a sub-move: 250 chains at the default 10 pay for the
    kernel, 6 chains do not."""
    if not mega_crank_fast.openmp_info()["enabled"]:
        pytest.skip("kernels built without OpenMP")
    moves = {"MOVE_SLITHER": 1.0}
    _, many = _parallel_report(
        tmp_path, "many", moves, [(250, "AABB")], [32, 32, 32], {}
    )
    _, few = _parallel_report(tmp_path, "few", moves, [(6, "AABB")], [32, 32, 32], {})
    assert (
        10 * 250 * 0.4e-3 * 0.75 > 0.1 + 1.0e-4 * 1000
        and "megamove too small" not in many
    )
    assert 10 * 6 * 0.4e-3 * 0.75 < 0.1 + 1.0e-4 * 24
    assert (
        "megamove too small for the parallel kernel - at SLITHER_SUBSTEPS : 10" in few
    )
    # a box that does not split runs the serial kernel: nothing to warn about
    _, small = _parallel_report(
        tmp_path,
        "small",
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5},
        [(6, "AABB")],
        [12, 12, 12],
        {},
    )
    assert "megamove too small" not in small


def test_parallel_report_warns_once_about_a_single_thread(
    tmp_path: pathlib.Path,
) -> None:
    """One thread pays the fixed cost for nothing: one warning for the whole
    report, however many parallel moves are enabled."""
    if not mega_crank_fast.openmp_info()["enabled"]:
        pytest.skip("kernels built without OpenMP")
    _, out = _parallel_report(
        tmp_path,
        "one",
        {"MOVE_CRANKSHAFT": 0.4, "MOVE_SLITHER": 0.3, "MOVE_PULL": 0.3},
        [(250, "AABB")],
        [32, 32, 32],
        {"PARALLEL_THREADS": 1, "CRANKSHAFT_SUBSTEPS": 10**6},
    )
    assert out.count("only one thread will run the parallel kernels") == 1
    assert "megamove too small" not in out


def test_parallel_report_says_where_the_thread_count_came_from(
    tmp_path: pathlib.Path,
) -> None:
    """``PARALLEL_THREADS : 0`` resolves to OMP_NUM_THREADS when set, otherwise to
    the CPUs available to the process; an explicit value is reported as such."""
    moves = {"MOVE_CRANKSHAFT": 1.0}
    _, explicit = _parallel_report(
        tmp_path,
        "explicit",
        moves,
        [(20, "AABB")],
        [32, 32, 32],
        {"PARALLEL_THREADS": 3},
    )
    assert (
        "Threads: 3 OpenMP threads per parallel megamove (from PARALLEL_THREADS; "
        in explicit
    )
    sim, auto = _parallel_report(
        tmp_path, "auto", moves, [(20, "AABB")], [32, 32, 32], {"PARALLEL_THREADS": 0}
    )
    assert (
        "Threads: %d OpenMP threads per parallel megamove (PARALLEL_THREADS : 0 -> OMP_NUM_THREADS if set, "
        "otherwise every available CPU; " % sim.parallel_threads
    ) in auto
    assert "all cores" not in auto and "CPUs available to this process" in auto


def test_parallel_report_on_an_anisotropic_box_and_chains_too_short_to_pull(
    tmp_path: pathlib.Path,
) -> None:
    """A 48 x 32 x 16 box splits differently along each axis; the report must
    print the kernel's own per-axis layout in x, y, z order. Monomers are below
    pull's three-bead minimum and are counted as never moved by it."""
    box = [48, 32, 16]
    _, out = _parallel_report(
        tmp_path,
        "aniso",
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_PULL": 0.5},
        [(20, "AABBA"), (7, "A")],
        box,
        {"CRANKSHAFT_SUBSTEPS": 50000},
    )
    crank = mega_crank_fast.parallel_crank_layout_info(48, 32, 16, False)
    assert len(set(crank["blocks"])) > 1, "fixture: the layout must differ between axes"
    assert (
        "block grid %s = %d blocks of %s sites"
        % (
            "x".join(map(str, crank["blocks"])),
            crank["num_blocks"],
            "x".join(map(str, crank["block_size"])),
        )
    ) in out
    assert "Box: 48x32x16 (3D, periodic); 27 chains, 107 beads" in out
    assert "(7 chain(s) shorter than 3 beads are never moved by this move)" in out
    _, no_short = _parallel_report(
        tmp_path,
        "noshort",
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_PULL": 0.5},
        [(20, "AABBA")],
        box,
        {"CRANKSHAFT_SUBSTEPS": 50000},
    )
    assert "are never moved by this move" not in no_short


# ---------------------------------------------------------------------------
# P2-5: no reference cycle through the analysis routines
# ---------------------------------------------------------------------------


def test_discarded_simulation_is_freed_without_the_cyclic_collector(
    tmp_path: pathlib.Path,
) -> None:
    """The analysis routines stored on the Simulation must not hold it strongly,
    or its lattice grids survive until the cyclic garbage collector runs."""
    gc.collect()
    gc.disable()
    try:
        state = _build(
            tmp_path,
            3,
            {"MOVE_CRANKSHAFT": 1.0},
            [(4, "AABB")],
            [10, 10, 10],
            extra={"ANA_RESIDUE_PAIRS": "0 3"},
        )
        sim_ref = weakref.ref(state.sim)
        lattice_ref = weakref.ref(state.sim.LATTICE)
        del state
        assert sim_ref() is None, (
            "the Simulation is still alive: a reference cycle holds it"
        )
        assert lattice_ref() is None
    finally:
        gc.enable()


def test_analysis_routines_run_in_the_same_order_at_the_same_frequencies(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Order and cadence are part of the output format: the order of the calls is
    the order rows are written in. The expected order is the keyword order in
    ``setup_analysis`` (ANA_POL, ANA_INTSCAL, ANA_DISTMAP, ANA_ACCEPTANCE,
    ANA_CLUSTER, ANA_INTER_RESIDUE, ANA_END_TO_END, ANA_CUSTOM, RESTART_FREQ),
    non-default frequencies first. ANA_CUSTOM is disabled without a module."""
    order = [
        "ANAFUNCT_polymeric_properties",
        "ANAFUNCT_internal_scaling",
        "ANAFUNCT_distance_map",
        "ANAFUNCT_acceptance",
        "ANAFUNCT_cluster_analysis",
        "ANAFUNCT_R2R_distance",
        "ANAFUNCT_end_to_end",
        "ANAFUNCT_save_restart",
    ]
    called: list[tuple[str, int]] = []

    def recorder(name: str) -> Any:
        def record(self: Simulation, step: int) -> None:
            called.append((name, step))

        record.__name__ = name
        return record

    for name in order:
        if name != "ANAFUNCT_R2R_distance":
            monkeypatch.setattr(Simulation, name, recorder(name))
    state = _build(
        tmp_path,
        3,
        {"MOVE_CRANKSHAFT": 1.0},
        [(4, "AABB")],
        [10, 10, 10],
        extra={
            "ANALYSIS_FREQ": 4,
            "ANA_CLUSTER": 6,
            "ANA_POL": 3,
            "RESTART_FREQ": 4,
            "ANA_INTSCAL": 4,
            "ANA_DISTMAP": 4,
            "ANA_ACCEPTANCE": 4,
            "ANA_INTER_RESIDUE": 4,
        },
    )
    sim = state.sim
    # the end-to-end analysis has no keyword of its own: it runs at the ANA_POL frequency
    off_default = (
        "ANAFUNCT_polymeric_properties",
        "ANAFUNCT_cluster_analysis",
        "ANAFUNCT_end_to_end",
    )

    def name(routine: Any) -> str:
        # a routine is stored as a method name or as a callable
        return routine if isinstance(routine, str) else routine.__name__

    assert [name(f) for f in sim.non_default_freq_analysis] == list(off_default)
    assert [name(f) for f in sim.default_freq_analysis] == [
        n for n in order if n not in off_default
    ]
    assert list(sim.non_default_freq_analysis.values()) == [3, 6, 3]

    sim.equilibration = 0
    for step in range(1, 13):
        sim.run_all_analysis(step)
    default = [n for n in order if n not in off_default + ("ANAFUNCT_R2R_distance",)]
    expected = []
    for step in range(1, 13):
        if step % 3 == 0:
            expected.append(("ANAFUNCT_polymeric_properties", step))
        if step % 6 == 0:
            expected.append(("ANAFUNCT_cluster_analysis", step))
        if step % 3 == 0:
            expected.append(("ANAFUNCT_end_to_end", step))
        if step % 4 == 0:
            expected.extend((n, step) for n in default)
    assert called == expected


# ---------------------------------------------------------------------------
# Review round: RF7-1 (cost), RF7-2 (dead moves), RF7-3 (deep copy), RF7-4..6
# ---------------------------------------------------------------------------


def test_cluster_only_check_searches_each_cluster_once(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """300 chains, nearly all in one cluster: the check must run one
    connected-component search per cluster, not one per chain (which is
    quadratic in a percolated box)."""
    state = _build(tmp_path, 3, {"MOVE_CRANKSHAFT": 1.0}, [(300, "AABB")], [14, 14, 14])
    sim = state.sim
    clusters = set(_clusters(_positions(sim), sim.LATTICE.dimensions).values())
    assert max(len(c) for c in clusters) > 250, (
        "fixture: one cluster must hold most chains"
    )

    seeds: list[int] = []
    real = simulation.lattice_utils.get_all_chains_in_connected_component

    def counting(chain_id: int, *args: Any, **kwargs: Any) -> Any:
        seeds.append(chain_id)
        return real(chain_id, *args, **kwargs)

    monkeypatch.setattr(
        simulation.lattice_utils, "get_all_chains_in_connected_component", counting
    )
    try:
        immobile = _check(sim, CLUSTER_ONLY)["immobile_chains"]
    except SimulationException:
        immobile = sorted(sim.LATTICE.chains)
    assert len(seeds) == len(clusters)
    # and the answer is the one the flood fill gives: only a cluster holding
    # every chain is stuck when nothing is frozen
    assert immobile == (sorted(sim.LATTICE.chains) if len(clusters) == 1 else [])


def test_a_move_that_cannot_act_does_not_hide_a_cluster_only_move_set(
    tmp_path: pathlib.Path,
) -> None:
    """Cluster-only is judged on the moves that can act. A chain rotation cannot
    act on monomers, a pull or a pivot cannot act on dimers: enabling one beside
    the cluster moves changes nothing about what can move."""
    # six isolated monomers, cluster rotation plus a rotation that cannot act
    monomers = _build(
        tmp_path, 3, {"MOVE_CRANKSHAFT": 1.0}, [(6, "A")], [15, 15, 15]
    ).sim
    assert all(
        len(m) == 1
        for m in _clusters(_positions(monomers), monomers.LATTICE.dimensions).values()
    )
    with pytest.raises(SimulationException) as refusal:
        _check(monomers, {"MOVE_CLUSTER_ROTATE": 0.5, "MOVE_CHAIN_ROTATE": 0.5})
    message = str(refusal.value)
    assert "MOVE_CHAIN_ROTATE cannot act on any mobile chain" in message
    # the refusal says it depends on the starting placement
    assert (
        "another SEED may start" in message
        and "contact graph is frozen either way" in message
    )

    # a mobile dimer touching a frozen chain: translate plus a pull that cannot act
    stuck_dir = tmp_path / "stuck"
    stuck_dir.mkdir()
    moves = {"MOVE_CLUSTER_TRANSLATE": 0.5, "MOVE_PULL": 0.5}
    stuck = _build(
        stuck_dir,
        3,
        moves,
        [(3, "AABBA"), (1, "AB")],
        [9, 9, 9],
        frozen=[1, 2, 3],
        n_steps=300,
        bypass=True,
    ).sim
    before = _positions(stuck)
    assert _clusters(before, stuck.LATTICE.dimensions)[4] & {1, 2, 3}, (
        "fixture: the dimer must touch a frozen chain"
    )
    with pytest.raises(
        SimulationException, match="No enabled move can ever change this configuration"
    ):
        _check(stuck)
    _run(stuck, stuck_dir)
    assert _positions(stuck) == before  # and it really is static

    # dimers with translate plus a pivot that cannot act: the contact graph is frozen
    dimers_dir = tmp_path / "dimers"
    dimers_dir.mkdir()
    dimers = _build(
        dimers_dir, 3, {"MOVE_CRANKSHAFT": 1.0}, [(12, "AB")], [12, 12, 12]
    ).sim
    checked = _check(dimers, {"MOVE_CLUSTER_TRANSLATE": 0.5, "MOVE_CHAIN_PIVOT": 0.5})
    assert any(
        "contact graph is frozen" in m
        and "MOVE_CHAIN_PIVOT cannot act on any mobile chain" in m
        for m in checked["messages"]
    )
    assert checked["null_fractions"] == {"MOVE_CHAIN_PIVOT": pytest.approx(0.5)}

    # but a pull that CAN act (there are pentamers) is not a dead move
    mixed_dir = tmp_path / "mixed"
    mixed_dir.mkdir()
    mixed = _build(
        mixed_dir, 3, {"MOVE_CRANKSHAFT": 1.0}, [(6, "AB"), (2, "AABBA")], [12, 12, 12]
    ).sim
    checked = _check(mixed, {"MOVE_CLUSTER_TRANSLATE": 0.5, "MOVE_PULL": 0.5})
    assert not any("contact graph is frozen" in m for m in checked["messages"])


def test_deep_copy_of_a_simulation_analyses_the_copy(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every analysis routine of a deep-copied Simulation must act on the copy:
    the built-in methods, the residue-pair routine and a custom module."""
    import copy

    names = [
        "ANAFUNCT_polymeric_properties",
        "ANAFUNCT_internal_scaling",
        "ANAFUNCT_distance_map",
        "ANAFUNCT_acceptance",
        "ANAFUNCT_cluster_analysis",
        "ANAFUNCT_end_to_end",
        "ANAFUNCT_save_restart",
    ]
    acted_on: list[tuple[str, Simulation]] = []

    def recorder(name: str) -> Any:
        def record(self: Simulation, step: int) -> None:
            acted_on.append((name, self))

        record.__name__ = name
        return record

    for name in names:
        monkeypatch.setattr(Simulation, name, recorder(name))
    written: list[Any] = []
    monkeypatch.setattr(
        simulation.analysis_IO,
        "write_residue_residue_distance",
        lambda step, pairs, data: written.append(data),
    )

    sim = _build(
        tmp_path,
        3,
        {"MOVE_CRANKSHAFT": 1.0},
        [(4, "AABB")],
        [10, 10, 10],
        extra={"ANA_RESIDUE_PAIRS": "0 3", "ANALYSIS_FREQ": 1, "RESTART_FREQ": 1},
    ).sim
    lattices: list[Any] = []
    keyword_lookup = dict(sim.keyword_lookup)
    keyword_lookup["ANALYSIS_MODULE"] = lambda step, lattice: lattices.append(lattice)
    keyword_lookup["ANA_CUSTOM"] = 1
    keyword_lookup["__DISABLED_FREQUENCIES"] = set(
        keyword_lookup.get("__DISABLED_FREQUENCIES", ())
    ) - {"ANA_CUSTOM"}
    sim.non_default_freq_analysis, sim.default_freq_analysis = sim.setup_analysis(
        keyword_lookup
    )

    twin = copy.deepcopy(sim)
    assert twin.LATTICE is not sim.LATTICE
    # mark the copy's chains, so a distance measured on the copy is recognisable
    for chain in twin.LATTICE.chains.values():
        chain.analysis_get_residue_residue_distance = lambda i, j, positions=None: (
            12345.0
        )
    twin.equilibration = 0
    os.chdir(tmp_path)
    twin.run_all_analysis(1)

    assert sorted(name for name, _ in acted_on) == sorted(names)
    assert all(target is twin for _, target in acted_on)
    assert acted_on[-1][0] == "ANAFUNCT_save_restart"  # the checkpoint still runs last
    assert lattices == [twin.LATTICE] and lattices[0] is twin.LATTICE
    assert written == [[[12345.0] * 4]]


def test_analysis_routine_called_after_its_simulation_is_gone_says_so(
    tmp_path: pathlib.Path,
) -> None:
    """A stored routine called with the step alone falls back on the Simulation
    that built it; once that is gone the error must say what happened."""
    state = _build(
        tmp_path,
        3,
        {"MOVE_CRANKSHAFT": 1.0},
        [(4, "AABB")],
        [10, 10, 10],
        extra={"ANA_RESIDUE_PAIRS": "0 3"},
    )
    routines = [
        r
        for r in list(state.sim.non_default_freq_analysis)
        + list(state.sim.default_freq_analysis)
        if callable(r)
    ]
    assert [r.__name__ for r in routines] == ["ANAFUNCT_R2R_distance"]
    gc.collect()
    gc.disable()
    try:
        del state
        with pytest.raises(AnalysisRoutineException, match="no longer exists"):
            routines[0](5)
    finally:
        gc.enable()


def test_parallel_report_without_openmp_does_not_suggest_more_threads(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """On a build without OpenMP no PARALLEL_THREADS value helps."""
    monkeypatch.setattr(
        simulation.mega_crank_fast,
        "openmp_info",
        lambda: {"enabled": False, "max_threads": 1},
    )
    _, out = _parallel_report(
        tmp_path,
        "noomp",
        {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5},
        [(250, "AABB")],
        [32, 32, 32],
        {"PARALLEL_THREADS": 4},
    )
    assert "set PARALLEL_THREADS above 1" not in out
    assert out.count("rebuild with OpenMP or drop PARALLELIZE") == 1
    assert "megamove too small" not in out


def test_megamove_warning_states_what_the_threads_can_save(
    tmp_path: pathlib.Path,
) -> None:
    """At two threads half the sampling time is the most that can be saved; the
    message must give that number, or 'X ms of sampling against a smaller fixed
    cost ... runs slower' reads as a contradiction."""
    if not mega_crank_fast.openmp_info()["enabled"]:
        pytest.skip("kernels built without OpenMP")
    n_beads = 1000
    fixed_ms = 0.1 + 1.0e-4 * n_beads
    substeps = 3000  # 0.3 ms of sampling: more than the 0.2 ms fixed cost, yet only 0.15 ms can be saved
    assert substeps * 1.0e-4 > fixed_ms > substeps * 1.0e-4 * 0.5
    _, out = _parallel_report(
        tmp_path,
        "two",
        {"MOVE_CRANKSHAFT": 1.0},
        [(250, "AABB")],
        [32, 32, 32],
        {"PARALLEL_THREADS": 2, "CRANKSHAFT_SUBSTEPS": substeps},
    )
    assert (
        "at CRANKSHAFT_SUBSTEPS : 3000, 2 threads can save at most about %.2g ms per megamove "
        "(50%% of %.2g ms of sampling), against a fixed cost of about %.2g ms"
        % (substeps * 1.0e-4 * 0.5, substeps * 1.0e-4, fixed_ms)
    ) in out
