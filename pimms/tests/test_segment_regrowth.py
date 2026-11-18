"""
Tests for MOVE_SEGMENT_REGROWTH (move code 15), the move that regrows short
windows of the free beads of partly frozen chains (FREEZE_FILE B lines).

What is checked:

* the path tables hold exactly the lattice paths a window can take, against an
  independent brute-force enumeration, and the window list holds exactly the
  windows of consecutive free beads;
* a window's new path is drawn uniformly, keeps the window bonded to its fixed
  neighbours through a periodic boundary, and is refused at a hard wall;
* the tracked energy stays the model energy (the oracle of test_energy_oracle)
  in 3D with SR, LR and SLR terms and angles, periodic and hardwall, and a
  rejected move leaves the grid, the type grid and the chains as they were;
* the free beads sample their exact Boltzmann distribution - including two free
  beads closing a 2x2 square with their frozen neighbours, which the crankshaft
  can never move - with a positive control that a 1.3x error in beta fails;
* the keywords are gated as experimental and range checked, the start-up checks
  and warnings account for the move, MOVE_FREQS.dat counts it in column 15, a
  RESTART_CONTINUE run with it reproduces the uninterrupted run, and with the
  move off the move draw is the one it was before the move existed.
"""

from __future__ import annotations

import contextlib
import io
import itertools
import os
import random
import types
from pathlib import Path
from typing import Any, Iterator, Sequence

import numpy as np
import pytest

from pimms import CONFIG, lattice_utils, segment_regrowth
from pimms.acceptance import AcceptanceCalculator
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import KeyFileException, MoveException, SimulationException
from pimms.moves import MoveObject
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U
from pimms.tests.test_energy_oracle import ChainSpec, build_hamiltonian, build_lattice, make_forcefield, oracle_energy
from pimms.tests.test_restart_continue import _assert_same_continuation, _stop_and_resume

# The fourteen move keywords that existed before MOVE_SEGMENT_REGROWTH, in code order.
OLDER_MOVE_KEYWORDS = ("MOVE_CRANKSHAFT", "MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE", "MOVE_CHAIN_PIVOT",
                       "MOVE_HEAD_PIVOT", "MOVE_SLITHER", "MOVE_CLUSTER_TRANSLATE", "MOVE_CLUSTER_ROTATE",
                       "MOVE_CTSMMC", "MOVE_MULTICHAIN_TSMMC", "MOVE_PULL", "MOVE_SYSTEM_TSMMC",
                       "MOVE_JUMP_AND_RELAX", "MOVE_VMMC")

# Path counts, measured by an independent enumeration when the move was written
# (and checked against brute force below where that is cheap): open walks of
# 1-4 beads, and closed paths summed over every anchor offset.
OPEN_WALK_COUNTS = {(2, 1): 8, (2, 2): 56, (2, 3): 368, (2, 4): 2336,
                    (3, 1): 26, (3, 2): 650, (3, 3): 15986, (3, 4): 388994}
CLOSED_PATH_TOTALS = {(2, 1): 56, (2, 2): 368, (2, 3): 2336, (2, 4): 14576,
                      (3, 1): 650, (3, 2): 15986, (3, 3): 388994}
# the largest number of closed paths for one anchor offset
CLOSED_PATH_MAXIMA = {(2, 1): 4, (2, 2): 12, (2, 3): 58, (2, 4): 254, (3, 1): 16, (3, 2): 192, (3, 3): 3008}

# Uniformity of the draw: draws per candidate path, and the chi-square bar
# (dof + K sigma, sigma = sqrt(2 dof)) for independent draws.
UNIFORM_DRAWS_PER_PATH = 400
UNIFORM_SIGMA = 5.0

# The exact-equilibrium fixtures (2D, 7x7). The temperature and force-field seed
# are those of the crankshaft's exact test in test_bead_freezing, where a 1.3x
# beta error is resolvable; the sample counts were set so the positive control
# clears its bar by 2.5x or more (measured: trap 3x at 20,000, tail 2.5x at 40,000).
EXACT_TEMPERATURE = 8.0
EXACT_FORCEFIELD_SEED = 20261109
EXACT_BOX = [7, 7]
EXACT_SWEEPS_PER_SAMPLE = 2
EXACT_SIGMA = 5.0
EXACT_CONTROL_BETA_SCALE = 1.3
EXACT_CONTROL_MARGIN = 2.0
# beads 0, 3, 4 and 5 frozen: beads 1 and 2 close a 2x2 square with beads 0 and 3
TRAP_CHAIN: ChainSpec = ("ABCABC", [[2, 2], [2, 3], [3, 3], [3, 2], [4, 2], [5, 2]])
TRAP_FROZEN = (0, 3, 4, 5)
TRAP_OBSTACLE: ChainSpec = ("CA", [[1, 4], [2, 5]])
TRAP_SAMPLES = 20000
# beads 0-3 frozen with bead 3 on the x = 6 face of a periodic box: the free tail
# (beads 4 and 5) starts across the boundary
TAIL_CHAIN: ChainSpec = ("ABCABC", [[3, 3], [4, 3], [5, 3], [6, 3], [0, 3], [0, 4]])
TAIL_FROZEN = (0, 1, 2, 3)
TAIL_OBSTACLE: ChainSpec = ("CA", [[1, 2], [1, 1]])
TAIL_SAMPLES = 40000

# The 3D bookkeeping fixture: a 9^3 box, a partly frozen chain with a four-bead
# free stretch held at both ends and a two-bead free tail (crossing the periodic
# boundary in the periodic layout), a chain with free tails at both ends, a
# frozen obstacle and a free dimer.
BOOKKEEPING_BOX = [9, 9, 9]
BOOKKEEPING_SWEEPS = 400
BOOKKEEPING_CHECK_EVERY = 40
BOOKKEEPING_FROZEN = {1: (0, 5), 2: (2,)}
BOOKKEEPING_TEMPERATURE = 30.0

# The chain order of a sweep: sweeps over three one-window chains, and the
# chi-square bar for the six orders being equally likely.
ORDER_SWEEPS = 6000
ORDER_SIGMA = 5.0

# Move selector check: random move sets, and draws per set.
SELECTOR_MOVE_SETS = 50
SELECTOR_DRAWS = 2000


def _bookkeeping_system(periodic: bool) -> list[ChainSpec]:
    """The 3D bookkeeping fixture.

    Parameters
    ----------
    periodic : bool
        If True, chain 1 runs across the x boundary; otherwise it starts on the
        x = 0 wall.

    Returns
    -------
    list of (str, list of list of int)
        Chains 1-4 in order.
    """
    xs = [6, 7, 8, 0, 1, 2, 3, 4] if periodic else list(range(8))
    return [("ABCABCAB", [[x, 4, 4] for x in xs]),
            ("BBAC", [[4, 6, z] for z in (4, 5, 6, 7)]),
            ("CA", [[2, 5, 4], [2, 5, 5]]),
            ("AB", [[6, 6, 6], [7, 6, 6]])]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _quiet() -> Iterator[None]:
    """Discard standard output.

    Yields
    ------
    None
    """
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def _acceptance(temperature: float, **moves: float) -> AcceptanceCalculator:
    """An AcceptanceCalculator for ``temperature`` with the given move fractions.

    Parameters
    ----------
    temperature : float
        Simulation temperature.
    **moves : float
        MOVE_* fractions; the others are 0. Defaults to the segment regrowth alone.

    Returns
    -------
    AcceptanceCalculator
        The calculator.
    """
    lookup = {keyword: 0.0 for keyword in OLDER_MOVE_KEYWORDS}
    lookup["MOVE_SEGMENT_REGROWTH"] = 0.0
    lookup.update(moves or {"MOVE_SEGMENT_REGROWTH": 1.0})
    return AcceptanceCalculator(temperature, lookup)


def _brute_force_paths(n_dim: int, n_beads: int) -> dict[tuple[int, ...], set[tuple[tuple[int, ...], ...]]]:
    """Every closed path, by anchor offset, enumerated without numpy.

    Parameters
    ----------
    n_dim : int
        Number of dimensions.
    n_beads : int
        Beads in the window.

    Returns
    -------
    dict of tuple to set of tuple
        Anchor offset -> the set of paths (each a tuple of sites) joining an
        anchor at the origin to one at that offset.
    """
    steps = [step for step in itertools.product((-1, 0, 1), repeat=n_dim) if any(step)]
    origin = (0,) * n_dim
    found: dict[tuple[int, ...], set[tuple[tuple[int, ...], ...]]] = {}
    for walk_steps in itertools.product(steps, repeat=n_beads):
        sites = list(itertools.accumulate(walk_steps, lambda a, b: tuple(x + y for x, y in zip(a, b)),
                                          initial=origin))[1:]
        if origin in sites or len(set(sites)) < n_beads:
            continue
        for last_step in steps:
            target = tuple(x + y for x, y in zip(sites[-1], last_step))
            if target != origin and target not in sites:
                found.setdefault(target, set()).add(tuple(sites))
    return found


def _bonded(a: Sequence[int], b: Sequence[int], dims: Sequence[int], periodic: bool) -> bool:
    """Whether two sites are bonded (Chebyshev distance 1).

    Parameters
    ----------
    a, b : sequence of int
        The sites.
    dims : sequence of int
        The box.
    periodic : bool
        Read distances through the periodic wrap.

    Returns
    -------
    bool
        True if bonded.
    """
    gaps = [abs(int(x) - int(y)) for x, y in zip(a, b)]
    if periodic:
        gaps = [min(gap, length - gap) for gap, length in zip(gaps, dims)]
    return max(gaps) == 1


def _lattice_positions(lattice: Any) -> dict[int, list[list[int]]]:
    """Every chain's positions, as plain lists.

    Parameters
    ----------
    lattice : pimms.lattice.Lattice
        The lattice.

    Returns
    -------
    dict of int to list of list of int
        chainID -> positions.
    """
    return {chain_id: [[int(x) for x in site] for site in lattice.chains[chain_id].get_ordered_positions()]
            for chain_id in lattice.chains}


def _assert_lattice_consistent(lattice: Any) -> None:
    """The grid and the type grid hold exactly what the chain objects say.

    Parameters
    ----------
    lattice : pimms.lattice.Lattice
        The lattice.

    Returns
    -------
    None
    """
    expected_grid = np.zeros_like(lattice.grid)
    expected_types = np.zeros_like(lattice.type_grid)
    for chain_id, chain in lattice.chains.items():
        codes = chain.get_intcode_sequence()
        for site, code in zip(chain.get_ordered_positions(), codes):
            assert expected_grid[tuple(site)] == 0, "two beads on site %s" % (list(site),)
            expected_grid[tuple(site)] = chain_id
            expected_types[tuple(site)] = code
    assert np.array_equal(lattice.grid, expected_grid)
    assert np.array_equal(lattice.type_grid, expected_types)


def _enumerate_states(chain: ChainSpec, frozen: tuple[int, ...], obstacle: ChainSpec, dims: list[int],
                      periodic: bool) -> list[tuple[tuple[int, ...], ...]]:
    """Every placement of the chain's free beads, by brute force over all sites.

    Parameters
    ----------
    chain : (str, list of list of int)
        The partly frozen chain.
    frozen : tuple of int
        Its frozen beads.
    obstacle : (str, list of list of int)
        A frozen chain.
    dims : list of int
        The box.
    periodic : bool
        Periodic boundaries.

    Returns
    -------
    list of tuple
        The sites of the free beads, in chain order, for every valid state.
    """
    free = [bead for bead in range(len(chain[1])) if bead not in frozen]
    occupied = {tuple(chain[1][bead]) for bead in frozen} | {tuple(site) for site in obstacle[1]}
    empty = [site for site in itertools.product(*(range(length) for length in dims)) if site not in occupied]
    states = []
    for combo in itertools.product(empty, repeat=len(free)):
        if len(set(combo)) < len(combo):
            continue
        positions = [tuple(site) for site in chain[1]]
        for bead, site in zip(free, combo):
            positions[bead] = site
        if all(_bonded(positions[i], positions[i + 1], dims, periodic) for i in range(len(positions) - 1)):
            states.append(combo)
    return states


def _boltzmann(energies: np.ndarray, beta: float) -> np.ndarray:
    """Normalised Boltzmann probabilities.

    Parameters
    ----------
    energies : numpy.ndarray
        State energies.
    beta : float
        Inverse temperature.

    Returns
    -------
    numpy.ndarray
        Probabilities.
    """
    weights = np.exp(-beta * (energies - energies.min()))
    return weights / weights.sum()


def _chi_square(counts: np.ndarray, probabilities: np.ndarray, effective_samples: float) -> float:
    """Pearson's chi-square of state counts against probabilities, for correlated samples.

    Parameters
    ----------
    counts : numpy.ndarray
        How often each state was sampled.
    probabilities : numpy.ndarray
        Expected probability of each state.
    effective_samples : float
        The number of independent samples the counts are worth.

    Returns
    -------
    float
        effective_samples x sum((observed - expected)^2 / expected).
    """
    observed = counts / counts.sum()
    return float(effective_samples * np.sum((observed - probabilities) ** 2 / probabilities))


def _run_dir(directory: Path, moves: dict[str, float], freeze_text: str | None, *,
             chains: Any = ((6, "AABBAABB"), (10, "A"), (4, "BBBB")), extra: dict[str, Any] | None = None,
             n_steps: int = 60, experimental: bool = True, box: Sequence[int] = (14, 14, 14)) -> None:
    """Write a 3D SR keyfile (and freeze file) in ``directory``.

    Parameters
    ----------
    directory : pathlib.Path
        The run directory (created).
    moves : dict of str to float
        MOVE_* fractions.
    freeze_text : str or None
        The freeze file, or None for no FREEZE_FILE.
    chains : sequence of (int, str), optional
        CHAIN lines.
    extra : dict, optional
        More keywords.
    n_steps : int, optional
        Steps to run.
    experimental : bool, optional
        Write EXPERIMENTAL_FEATURES : True (otherwise False).
    box : sequence of int, optional
        DIMENSIONS.

    Returns
    -------
    None
    """
    directory.mkdir(parents=True, exist_ok=True)
    U.write_param_file(str(directory / "params.prm"), "SR")
    keywords: dict[str, Any] = {"ENERGY_CHECK": 1, "XTC_FREQ": 1000, "EN_FREQ": 1, "ANALYSIS_FREQ": 1,
                                "RESTART_FREQ": 10**9, "CRANKSHAFT_SUBSTEPS": 200}
    if freeze_text is not None:
        keywords["FREEZE_FILE"] = "freeze.txt"
        (directory / "freeze.txt").write_text(freeze_text)
    keywords.update(extra or {})
    U.write_keyfile(str(directory / "KEYFILE.kf"), 3, False, moves, box=list(box), chains=list(chains),
                    n_steps=n_steps, equilibration=1, seed=5, temperature=40, extra=keywords)
    if not experimental:
        keyfile = directory / "KEYFILE.kf"
        keyfile.write_text(keyfile.read_text().replace("EXPERIMENTAL_FEATURES : True",
                                                       "EXPERIMENTAL_FEATURES : False"))


def _parse(directory: Path) -> dict[str, Any]:
    """Parse the keyfile in ``directory``.

    Parameters
    ----------
    directory : pathlib.Path
        The run directory.

    Returns
    -------
    dict
        The keyword lookup.
    """
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        with _quiet():
            return KeyFileParser("KEYFILE.kf").keyword_lookup
    finally:
        os.chdir(cwd)


def _build(directory: Path) -> tuple[Simulation, str]:
    """Build (but do not run) the simulation in ``directory``.

    Parameters
    ----------
    directory : pathlib.Path
        The run directory.

    Returns
    -------
    tuple of (Simulation, str)
        The simulation and what start-up printed.
    """
    printed = io.StringIO()
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        with contextlib.redirect_stdout(printed):
            simulation = Simulation(KeyFileParser("KEYFILE.kf").keyword_lookup)
    finally:
        os.chdir(cwd)
    return simulation, " ".join(printed.getvalue().split())


def _finish(simulation: Simulation, directory: Path) -> None:
    """Run a built simulation in ``directory``.

    Parameters
    ----------
    simulation : Simulation
        Built by :func:`_build`.
    directory : pathlib.Path
        Its run directory.

    Returns
    -------
    None
    """
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        with _quiet():
            simulation.run_simulation()
    finally:
        os.chdir(cwd)


# ---------------------------------------------------------------------------
# the path tables and the window list
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_dim, n_beads", sorted(OPEN_WALK_COUNTS))
def test_open_walk_tables_have_the_measured_sizes(n_dim: int, n_beads: int) -> None:
    walks = segment_regrowth.open_walks(n_dim, n_beads)
    assert walks.shape == (OPEN_WALK_COUNTS[(n_dim, n_beads)], n_beads, n_dim)
    assert not walks.flags.writeable


@pytest.mark.parametrize("n_dim, n_beads", sorted(CLOSED_PATH_TOTALS))
def test_closed_path_tables_have_the_measured_sizes(n_dim: int, n_beads: int) -> None:
    reach = n_beads + 1
    sizes = [len(segment_regrowth.closed_path_indices(n_dim, n_beads, offset))
             for offset in itertools.product(range(-reach, reach + 1), repeat=n_dim) if any(offset)]
    assert sum(sizes) == CLOSED_PATH_TOTALS[(n_dim, n_beads)]
    assert max(sizes) == CLOSED_PATH_MAXIMA[(n_dim, n_beads)]


@pytest.mark.parametrize("n_dim, n_beads", [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2)])
def test_closed_paths_are_exactly_the_brute_force_paths(n_dim: int, n_beads: int) -> None:
    """Every anchor offset gets exactly the paths an independent enumeration finds."""
    expected = _brute_force_paths(n_dim, n_beads)
    walks = segment_regrowth.open_walks(n_dim, n_beads)
    reach = n_beads + 1
    for offset in itertools.product(range(-reach, reach + 1), repeat=n_dim):
        if not any(offset):
            continue
        indices = segment_regrowth.closed_path_indices(n_dim, n_beads, offset)
        got = {tuple(tuple(int(x) for x in site) for site in walks[index]) for index in indices}
        assert len(got) == len(indices)
        assert got == expected.get(offset, set()), offset


def test_closed_paths_refuse_a_window_whose_anchors_coincide() -> None:
    with pytest.raises(MoveException, match="same site"):
        segment_regrowth.closed_path_indices(3, 2, (0, 0, 0))


@pytest.mark.parametrize("n_dim, n_beads", [(1, 2), (4, 2), (3, 0), (3, CONFIG.MAX_SEGMENT_REGROWTH_BEADS + 1)])
def test_path_tables_refuse_unsupported_requests(n_dim: int, n_beads: int) -> None:
    with pytest.raises(MoveException):
        segment_regrowth.open_walks(n_dim, n_beads)


def test_windows_are_every_run_of_up_to_max_free_beads() -> None:
    # beads 0 and 4 frozen in a chain of 8: free stretch 1-3 (enclosed) and tail 5-7
    windows = segment_regrowth.free_segment_windows(8, (4, 0), 2)
    assert windows == [(1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (5, 5), (5, 6), (6, 6), (6, 7), (7, 7)]
    assert segment_regrowth.free_segment_windows(8, (0, 4), 1) == [(b, b) for b in (1, 2, 3, 5, 6, 7)]
    # a window never covers a frozen bead, however long it may be
    assert segment_regrowth.free_segment_windows(5, (2,), 4) == [(0, 0), (0, 1), (1, 1), (3, 3), (3, 4), (4, 4)]


@pytest.mark.parametrize("length, frozen, max_beads", [(5, (), 2), (3, (0, 1, 2), 2), (4, (4,), 2),
                                                       (4, (0,), 0), (4, (0,), 5)])
def test_window_list_refuses_what_is_not_a_partly_frozen_chain(length: int, frozen: tuple[int, ...],
                                                               max_beads: int) -> None:
    with pytest.raises(MoveException):
        segment_regrowth.free_segment_windows(length, frozen, max_beads)


def test_only_partly_frozen_chains_get_windows() -> None:
    windows = segment_regrowth.segment_windows({1: 4, 2: 6, 3: 3}, {2: (0,), 3: (1,)}, 3)
    assert sorted(windows) == [2, 3]
    assert windows[3] == [(0, 0), (2, 2)]


# ---------------------------------------------------------------------------
# drawing a new path
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("first, last", [(2, 3), (4, 5), (0, 1)])
def test_a_window_path_is_drawn_uniformly_from_its_candidates(first: int, last: int) -> None:
    """Interior, C-terminal and N-terminal windows: every candidate equally often."""
    dims = [30, 30, 30]
    chain = [[10, 10, 10], [11, 10, 10], [12, 11, 10], [13, 11, 11], [14, 12, 11], [15, 12, 12]]
    n_beads = last - first + 1
    walks = segment_regrowth.open_walks(3, n_beads)
    if 0 < first and last < len(chain) - 1:
        offset = segment_regrowth.window_anchor_offset(chain, first, last, dims)
        candidates = [walks[index] for index in segment_regrowth.closed_path_indices(3, n_beads, offset)]
        anchor = chain[first - 1]
    elif first > 0:
        candidates, anchor = list(walks), chain[first - 1]
    else:
        candidates, anchor = [walk[::-1] for walk in walks], chain[last + 1]
    key = {tuple(tuple(int(a) + int(r) for a, r in zip(anchor, site)) for site in walk): i
           for i, walk in enumerate(candidates)}
    assert len(key) == len(candidates)

    random.seed(11)
    counts = np.zeros(len(candidates))
    for _ in range(UNIFORM_DRAWS_PER_PATH * len(candidates)):
        drawn = segment_regrowth.draw_window_path(chain, first, last, dims, hardwall=False)
        assert drawn is not None
        counts[key[tuple(tuple(site) for site in drawn)]] += 1
    expected = UNIFORM_DRAWS_PER_PATH
    dof = len(candidates) - 1
    assert float(np.sum((counts - expected) ** 2 / expected)) < dof + UNIFORM_SIGMA * np.sqrt(2.0 * dof)


def test_a_window_across_the_periodic_boundary_stays_bonded_to_its_neighbours() -> None:
    dims = [7, 7, 7]
    chain = [[4, 3, 3], [5, 3, 3], [6, 3, 3], [0, 3, 3], [1, 3, 3], [2, 3, 3]]
    assert segment_regrowth.window_anchor_offset(chain, 1, 4, dims) == (5, 0, 0)
    random.seed(3)
    for _ in range(500):
        drawn = segment_regrowth.draw_window_path(chain, 1, 4, dims, hardwall=False)
        new_chain = [chain[0]] + drawn + [chain[5]]
        assert all(0 <= x < 7 for site in drawn for x in site)
        assert all(_bonded(a, b, dims, True) for a, b in zip(new_chain, new_chain[1:]))
        # the window still reaches its far neighbour the long way round
        assert segment_regrowth.window_anchor_offset(new_chain, 1, 4, dims) == (5, 0, 0)


def test_a_path_through_a_hard_wall_is_refused_and_no_other_leaves_the_box() -> None:
    dims = [8, 8, 8]
    chain = [[0, 3, 3], [0, 4, 3], [1, 4, 3], [2, 4, 3]]
    random.seed(4)
    outcomes = [segment_regrowth.draw_window_path(chain, 0, 1, dims, hardwall=True) for _ in range(2000)]
    refused = sum(drawn is None for drawn in outcomes)
    assert 0 < refused < len(outcomes)
    for drawn in outcomes:
        if drawn is not None:
            assert all(0 <= x < 8 for site in drawn for x in site)
            assert _bonded(drawn[-1], chain[2], dims, False)


def test_a_corrupt_chain_is_refused() -> None:
    with pytest.raises(MoveException, match="not bonded"):
        segment_regrowth.draw_window_path([[0, 0, 0], [3, 0, 0], [4, 0, 0], [5, 0, 0]], 1, 2, [9, 9, 9], False)
    with pytest.raises(MoveException, match="with a bead to spare"):
        segment_regrowth.draw_window_path([[0, 0, 0], [1, 0, 0]], 0, 1, [9, 9, 9], False)


# ---------------------------------------------------------------------------
# the move: energy bookkeeping and undoing a rejection
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind, angle_mode", [("SR", "abs"), ("LR", "tnorm"), ("SLR", "abs")])
@pytest.mark.parametrize("hardwall", [False, True], ids=["periodic", "hardwall"])
def test_the_tracked_energy_stays_the_model_energy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str,
                                                   angle_mode: str, hardwall: bool) -> None:
    monkeypatch.chdir(tmp_path)
    system = _bookkeeping_system(periodic=not hardwall)
    forcefield = make_forcefield(np.random.default_rng(7), kind, angle_mode, BOOKKEEPING_TEMPERATURE)
    with _quiet():
        hamiltonian = build_hamiltonian(forcefield, 3, hardwall)
        lattice = build_lattice(hamiltonian, BOOKKEEPING_BOX, hardwall, system)
    windows = segment_regrowth.segment_windows({1: 8, 2: 4}, BOOKKEEPING_FROZEN, 4)
    acceptance = _acceptance(BOOKKEEPING_TEMPERATURE)
    mover = MoveObject()
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])
    assert energy == oracle_energy(system, BOOKKEEPING_BOX, hardwall, forcefield).terms()[0]
    start = _lattice_positions(lattice)

    random.seed(9)
    attempted = accepted = 0
    for sweep in range(1, BOOKKEEPING_SWEEPS + 1):
        lattice, energy, n_attempted, n_accepted = mover.segment_regrowth(lattice, energy, acceptance, hamiltonian,
                                                                          windows, hardwall)
        attempted += n_attempted
        accepted += n_accepted
        if sweep % BOOKKEEPING_CHECK_EVERY:
            continue
        now = _lattice_positions(lattice)
        chains = [(spec[0], now[chain_id]) for chain_id, spec in enumerate(system, start=1)]
        assert energy == oracle_energy(chains, BOOKKEEPING_BOX, hardwall, forcefield).terms()[0]
        assert energy == int(hamiltonian.evaluate_total_energy(lattice)[0])
        _assert_lattice_consistent(lattice)
        for chain_id, positions in now.items():
            assert all(_bonded(a, b, BOOKKEEPING_BOX, not hardwall) for a, b in zip(positions, positions[1:]))
            frozen = BOOKKEEPING_FROZEN.get(chain_id, range(len(positions)))
            assert [positions[bead] for bead in frozen] == [start[chain_id][bead] for bead in frozen]

    assert attempted == 2 * BOOKKEEPING_SWEEPS
    assert 0 < accepted < attempted
    moved = _lattice_positions(lattice)
    assert moved[1] != start[1] and moved[2] != start[2]


def test_a_rejected_window_leaves_everything_as_it_was(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    system = _bookkeeping_system(periodic=True)
    forcefield = make_forcefield(np.random.default_rng(7), "SLR", "abs", BOOKKEEPING_TEMPERATURE)
    with _quiet():
        hamiltonian = build_hamiltonian(forcefield, 3, False)
        lattice = build_lattice(hamiltonian, BOOKKEEPING_BOX, False, system)
    windows = segment_regrowth.segment_windows({1: 8, 2: 4}, BOOKKEEPING_FROZEN, 4)
    never = types.SimpleNamespace(boltzmann_acceptance=lambda old, new: False)
    start = _lattice_positions(lattice)
    grid, type_grid = lattice.grid.copy(), lattice.type_grid.copy()
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])

    random.seed(2)
    for _ in range(200):
        lattice, new_energy, _, n_accepted = MoveObject().segment_regrowth(lattice, energy, never, hamiltonian,
                                                                           windows, False)
        assert n_accepted == 0 and new_energy == energy
    assert _lattice_positions(lattice) == start
    assert np.array_equal(lattice.grid, grid) and np.array_equal(lattice.type_grid, type_grid)


def test_a_draw_of_the_current_path_is_a_rejected_null_move(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A four-bead window stretched along the body diagonal between its frozen ends
    (every bond +1 on every axis) has one path: its own."""
    monkeypatch.chdir(tmp_path)
    assert len(segment_regrowth.closed_path_indices(3, 4, (5, 5, 5))) == 1
    straight: ChainSpec = ("ABCABC", [[x, x, x] for x in range(1, 7)])
    forcefield = make_forcefield(np.random.default_rng(7), "SR", "abs", BOOKKEEPING_TEMPERATURE)
    with _quiet():
        hamiltonian = build_hamiltonian(forcefield, 3, False)
        lattice = build_lattice(hamiltonian, BOOKKEEPING_BOX, False, [straight])
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])
    random.seed(4)
    for _ in range(50):
        lattice, energy, n_attempted, n_accepted = MoveObject().segment_regrowth(
            lattice, energy, _acceptance(BOOKKEEPING_TEMPERATURE), hamiltonian, {1: [(1, 4)]}, False)
        assert (n_attempted, n_accepted) == (1, 0)
    assert _lattice_positions(lattice)[1] == straight[1]


def test_every_accepted_window_move_changes_the_configuration(tmp_path: Path,
                                                              monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    forcefield = make_forcefield(np.random.default_rng(EXACT_FORCEFIELD_SEED), "SR", "abs", EXACT_TEMPERATURE)
    with _quiet():
        hamiltonian = build_hamiltonian(forcefield, 2, True)
        lattice = build_lattice(hamiltonian, EXACT_BOX, True, [TRAP_CHAIN, TRAP_OBSTACLE])
    windows = segment_regrowth.segment_windows({1: 6}, {1: TRAP_FROZEN}, 2)
    acceptance = _acceptance(EXACT_TEMPERATURE)
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])
    random.seed(8)
    accepted = changed = 0
    for _ in range(2000):
        before = _lattice_positions(lattice)[1]
        lattice, energy, _, n_accepted = MoveObject().segment_regrowth(lattice, energy, acceptance, hamiltonian,
                                                                       windows, True)
        accepted += n_accepted
        changed += _lattice_positions(lattice)[1] != before
    assert accepted == changed > 0


def test_a_sweep_takes_the_chains_in_a_uniformly_random_order() -> None:
    """Every order equally likely, which is what makes the sweep reversible."""
    seen: list[int] = []
    mover = MoveObject()
    mover._place_regrown_window = lambda lattice, chain_id, first, last, hardwall: seen.append(chain_id)
    windows = {1: [(1, 1)], 2: [(1, 1)], 3: [(1, 1)]}
    random.seed(21)
    for _ in range(ORDER_SWEEPS):
        mover.segment_regrowth(None, 0, None, None, windows, False)
    orders = [tuple(seen[i:i + 3]) for i in range(0, len(seen), 3)]
    counts = np.array([orders.count(order) for order in itertools.permutations((1, 2, 3))])
    assert counts.sum() == ORDER_SWEEPS
    expected = ORDER_SWEEPS / 6.0
    dof = 5
    assert float(np.sum((counts - expected) ** 2 / expected)) < dof + ORDER_SIGMA * np.sqrt(2.0 * dof)


def test_a_sweep_with_no_windows_does_nothing_and_draws_nothing() -> None:
    random.seed(1)
    state = random.getstate()
    lattice = types.SimpleNamespace()
    result = MoveObject().segment_regrowth(lattice, 17, None, None, {}, False)
    assert result == (lattice, 17, 0, 0)
    assert random.getstate() == state


# ---------------------------------------------------------------------------
# the exact equilibrium of the free beads
# ---------------------------------------------------------------------------

def _sample(tmp_path: Path, chain: ChainSpec, frozen: tuple[int, ...], obstacle: ChainSpec, periodic: bool,
            max_beads: int, n_samples: int) -> tuple[list[tuple[tuple[int, ...], ...]], np.ndarray, np.ndarray,
                                                     np.ndarray]:
    """Run the move on a 2D fixture and count the states of its free beads.

    Parameters
    ----------
    tmp_path : pathlib.Path
        The test's directory; the caller has made it the working directory, where
        the parameter file is written.
    chain, obstacle : (str, list of list of int)
        The partly frozen chain and a frozen chain.
    frozen : tuple of int
        The chain's frozen beads.
    periodic : bool
        Periodic boundaries (otherwise hard walls).
    max_beads : int
        SEGMENT_REGROWTH_MAX_BEADS.
    n_samples : int
        Samples, EXACT_SWEEPS_PER_SAMPLE sweeps apart.

    Returns
    -------
    tuple
        The states, their model energies, the sample counts per state and the
        energy trace.
    """
    hardwall = not periodic
    states = _enumerate_states(chain, frozen, obstacle, EXACT_BOX, periodic)
    forcefield = make_forcefield(np.random.default_rng(EXACT_FORCEFIELD_SEED), "SR", "abs", EXACT_TEMPERATURE)
    free = [bead for bead in range(len(chain[1])) if bead not in frozen]
    energies = []
    for state in states:
        placed = [list(site) for site in chain[1]]
        for bead, site in zip(free, state):
            placed[bead] = list(site)
        energies.append(oracle_energy([(chain[0], placed), obstacle], EXACT_BOX, hardwall, forcefield).terms()[0])

    with _quiet():
        hamiltonian = build_hamiltonian(forcefield, 2, hardwall)
        lattice = build_lattice(hamiltonian, EXACT_BOX, hardwall, [chain, obstacle])
    windows = segment_regrowth.segment_windows({1: len(chain[1])}, {1: frozen}, max_beads)
    acceptance = _acceptance(EXACT_TEMPERATURE)
    mover = MoveObject()
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])
    random.seed(5)
    index = {state: i for i, state in enumerate(states)}
    counts = np.zeros(len(states))
    trace = np.zeros(n_samples)
    for sample in range(n_samples):
        for _ in range(EXACT_SWEEPS_PER_SAMPLE):
            lattice, energy, _, _ = mover.segment_regrowth(lattice, energy, acceptance, hamiltonian, windows,
                                                           hardwall)
        now = lattice.chains[1].get_ordered_positions()
        counts[index[tuple(tuple(int(x) for x in now[bead]) for bead in free)]] += 1
        trace[sample] = energy
    final = [(chain[0], [[int(x) for x in site] for site in lattice.chains[1].get_ordered_positions()]), obstacle]
    assert energy == oracle_energy(final, EXACT_BOX, hardwall, forcefield).terms()[0]
    return states, np.asarray(energies, dtype=float), counts, trace


@pytest.mark.parametrize("case", ["trap", "periodic-tail"])
def test_the_free_beads_sample_their_exact_distribution(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                                        case: str) -> None:
    """Sampled state frequencies match the enumerated Boltzmann weights.

    The trap is two free beads closing a 2x2 square with their frozen
    neighbours, which the crankshaft cannot move at all; the tail starts across
    a periodic boundary.
    """
    monkeypatch.chdir(tmp_path)
    if case == "trap":
        chain, frozen, obstacle, periodic, n_samples = TRAP_CHAIN, TRAP_FROZEN, TRAP_OBSTACLE, False, TRAP_SAMPLES
    else:
        chain, frozen, obstacle, periodic, n_samples = TAIL_CHAIN, TAIL_FROZEN, TAIL_OBSTACLE, True, TAIL_SAMPLES
    states, energies, counts, trace = _sample(tmp_path, chain, frozen, obstacle, periodic, 2, n_samples)
    beta = CONFIG.INVTEMP_FACTOR / EXACT_TEMPERATURE
    exact = _boltzmann(energies, beta)

    effective = n_samples / (2.0 * U.integrated_autocorr_time(trace))
    dof = len(states) - 1
    threshold = dof + EXACT_SIGMA * np.sqrt(2.0 * dof)
    assert _chi_square(counts, exact, effective) < threshold
    assert counts[exact > 1e-3].min() > 0, "a reachable state was never sampled"

    # positive control: the same samples against a beta 30% off must fail by a wide margin
    wrong = _boltzmann(energies, EXACT_CONTROL_BETA_SCALE * beta)
    assert _chi_square(counts, wrong, effective) > EXACT_CONTROL_MARGIN * threshold


def test_the_trap_holds_against_single_bead_moves(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Why the move exists: with one-bead windows, or the crankshaft, the square never changes."""
    monkeypatch.chdir(tmp_path)
    states, _, counts, _ = _sample(tmp_path, TRAP_CHAIN, TRAP_FROZEN, TRAP_OBSTACLE, False, 1, 500)
    assert len(states) > 1 and np.count_nonzero(counts) == 1

    forcefield = make_forcefield(np.random.default_rng(EXACT_FORCEFIELD_SEED), "SR", "abs", EXACT_TEMPERATURE)
    with _quiet():
        hamiltonian = build_hamiltonian(forcefield, 2, True)
        lattice = build_lattice(hamiltonian, EXACT_BOX, True, [TRAP_CHAIN, TRAP_OBSTACLE])
    acceptance = types.SimpleNamespace(invtemp=CONFIG.INVTEMP_FACTOR / EXACT_TEMPERATURE)
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])
    random.seed(5)
    np.random.seed(5)
    lattice, _, _, _ = MoveObject().system_shake(lattice, energy, acceptance, hamiltonian, 2000, "UNIFORM",
                                                 hardwall=True, frozen_chains=[1, 2],
                                                 frozen_beads={1: TRAP_FROZEN})
    assert _lattice_positions(lattice)[1] == TRAP_CHAIN[1]


# ---------------------------------------------------------------------------
# keywords and start-up
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("keyword, value", [("MOVE_SEGMENT_REGROWTH", None), ("SEGMENT_REGROWTH_MAX_BEADS", 2)])
def test_the_keywords_are_experimental(tmp_path: Path, keyword: str, value: Any) -> None:
    moves = {"MOVE_CRANKSHAFT": 0.5, "MOVE_SEGMENT_REGROWTH": 0.5} if value is None else {"MOVE_CRANKSHAFT": 1.0}
    extra = {} if value is None else {keyword: value}
    _run_dir(tmp_path / "run", moves, "B 1 0 4\n", extra=extra, experimental=False)
    with pytest.raises(KeyFileException, match="EXPERIMENTAL_FEATURES"):
        _parse(tmp_path / "run")


@pytest.mark.parametrize("value, match", [(0, "larger than 0"), (5, "may not exceed 4")])
def test_the_window_length_is_range_checked(tmp_path: Path, value: int, match: str) -> None:
    _run_dir(tmp_path / "run", {"MOVE_CRANKSHAFT": 0.5, "MOVE_SEGMENT_REGROWTH": 0.5}, "B 1 0 4\n",
             extra={"SEGMENT_REGROWTH_MAX_BEADS": value})
    with pytest.raises(KeyFileException, match=match):
        _parse(tmp_path / "run")


def test_the_window_length_without_the_move_is_warned_about(tmp_path: Path) -> None:
    _run_dir(tmp_path / "run", {"MOVE_CRANKSHAFT": 1.0}, None, extra={"SEGMENT_REGROWTH_MAX_BEADS": 2})
    printed = io.StringIO()
    cwd = os.getcwd()
    os.chdir(tmp_path / "run")
    try:
        with contextlib.redirect_stdout(printed):
            KeyFileParser("KEYFILE.kf")
    finally:
        os.chdir(cwd)
    assert "SEGMENT_REGROWTH_MAX_BEADS set for a move that is never drawn" in printed.getvalue()


def test_a_run_moves_the_free_beads_and_counts_the_attempts(tmp_path: Path) -> None:
    """No crankshaft: only the segment regrowth can move the free beads, and it does."""
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_SEGMENT_REGROWTH": 0.6, "MOVE_CHAIN_TRANSLATE": 0.4}, "B 1 0 4\nB 2 2 7\n")
    simulation, printed = _build(directory)
    assert sorted(simulation.segment_windows) == [1, 2]
    start = _lattice_positions(simulation.LATTICE)
    _finish(simulation, directory)
    end = _lattice_positions(simulation.LATTICE)
    for chain_id, frozen in ((1, (0, 4)), (2, (2, 7))):
        assert [end[chain_id][bead] for bead in frozen] == [start[chain_id][bead] for bead in frozen]
        assert end[chain_id] != start[chain_id]

    attempts = np.loadtxt(directory / "MOVE_FREQS.dat")[-1]
    accepted = np.loadtxt(directory / "ACCEPTANCE.dat")[-1]
    assert len(attempts) == 16
    assert attempts[15] > 0 and attempts[15] % 2 == 0, "one attempt per partly frozen chain per sweep"
    assert 0 < accepted[15] <= attempts[15]
    assert attempts[1] == 0
    assert "partly frozen (FREEZE_FILE B lines)" not in printed


def test_without_partly_frozen_chains_the_move_does_nothing_and_says_so(tmp_path: Path) -> None:
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_CRANKSHAFT": 0.75, "MOVE_SEGMENT_REGROWTH": 0.25}, "C 3\n")
    simulation, printed = _build(directory)
    assert simulation.segment_windows == {}
    assert "MOVE_SEGMENT_REGROWTH only regrows the free beads of partly frozen chains" in printed
    assert "the 25.0% of steps it takes change nothing" in printed
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        with _quiet():
            summary = simulation.check_moveset_applicability(simulation.keyword_lookup)
    finally:
        os.chdir(cwd)
    assert summary["null_fractions"] == {"MOVE_SEGMENT_REGROWTH": pytest.approx(0.25)}
    assert summary["null_fraction"] == pytest.approx(0.25)
    _finish(simulation, directory)
    assert np.loadtxt(directory / "MOVE_FREQS.dat")[-1][15] == 0


def test_the_move_alone_with_nothing_to_regrow_is_refused(tmp_path: Path) -> None:
    _run_dir(tmp_path / "run", {"MOVE_SEGMENT_REGROWTH": 1.0}, None)
    with pytest.raises(SimulationException, match="only regrows the free beads of partly frozen chains"):
        _build(tmp_path / "run")


def test_the_move_alone_keeps_a_run_whose_every_chain_is_partly_frozen(tmp_path: Path) -> None:
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_SEGMENT_REGROWTH": 0.5, "MOVE_CHAIN_PIVOT": 0.5}, "B 1 0\nB 2 5\nB 3 2-3\n",
             chains=[(3, "AABBAA")])
    simulation, printed = _build(directory)
    assert "only the crankshaft and MOVE_SEGMENT_REGROWTH move them" in printed
    assert "MOVE_CHAIN_PIVOT cannot act on a chain with a frozen bead" in printed
    start = _lattice_positions(simulation.LATTICE)
    _finish(simulation, directory)
    assert _lattice_positions(simulation.LATTICE) != start


def test_free_chains_nothing_can_move_are_reported_not_refused(tmp_path: Path) -> None:
    """The segment regrowth keeps the partly frozen chains moving; the free ones never move."""
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_SEGMENT_REGROWTH": 1.0}, "B 1 0 4\n")
    simulation, printed = _build(directory)
    assert ("No enabled move can act on 19 of 19 mobile chains (10 of 1 bead, 4 of 4 beads, 5 of 8 beads): "
            "MOVE_SEGMENT_REGROWTH only regrows the free beads of partly frozen chains") in printed
    assert "contains no move that can change a chain's shape" not in printed


def test_excursions_whose_only_sub_move_is_the_regrowth_are_accepted(tmp_path: Path,
                                                                     monkeypatch: pytest.MonkeyPatch) -> None:
    """The excursions really run the regrowth, never a crankshaft fallback."""
    inside_excursions = []
    crankshaft_calls = []
    real_sweep = MoveObject.segment_regrowth

    def recording_sweep(self: MoveObject, lattice: Any, energy: Any, acceptance: Any, *args: Any,
                        **kwargs: Any) -> Any:
        inside_excursions.append(bool(acceptance.auxillary_chain))
        return real_sweep(self, lattice, energy, acceptance, *args, **kwargs)

    monkeypatch.setattr(MoveObject, "segment_regrowth", recording_sweep)
    monkeypatch.setattr(MoveObject, "system_shake", lambda *a, **k: crankshaft_calls.append(1))
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_SYSTEM_TSMMC": 0.2, "MOVE_SEGMENT_REGROWTH": 0.8}, "B 1 0\nB 2 5\nB 3 2-3\n",
             chains=[(3, "AABBAA")], n_steps=40)
    simulation, printed = _build(directory)
    assert "cannot act" not in printed and "refusing" not in printed
    _finish(simulation, directory)
    attempts = np.loadtxt(directory / "MOVE_FREQS.dat")[-1]
    assert attempts[12] > 0 and attempts[15] > 0
    assert any(inside_excursions) and not all(inside_excursions)
    assert not crankshaft_calls


def test_cluster_moves_with_the_regrowth_are_not_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A free dimer touching a partly frozen chain can move once the regrowth pulls the chain away."""
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_SEGMENT_REGROWTH": 0.5, "MOVE_CLUSTER_TRANSLATE": 0.5}, "B 1 0 4\n",
             chains=[(1, "AABBAABB"), (1, "AB")], box=(7, 7, 7))
    # precondition: with this SEED the dimer starts touching the partly frozen chain
    with monkeypatch.context() as patched:
        patched.setattr(Simulation, "check_moveset_applicability", lambda self, keyword_lookup: {})
        unchecked, _ = _build(directory)
    assert 1 in lattice_utils.get_all_chains_in_connected_component(
        2, unchecked.LATTICE.grid, unchecked.LATTICE.chains, threshold=None, useChains=True, hardwall=False)

    simulation, printed = _build(directory)
    assert "only contacts with the partly frozen chains, which MOVE_SEGMENT_REGROWTH reshapes" in printed
    assert "contact graph is frozen" not in printed
    _finish(simulation, directory)


def test_excursions_running_the_regrowth_are_not_null_moves(tmp_path: Path) -> None:
    directory = tmp_path / "run"
    _run_dir(directory, {"MOVE_CHAIN_PIVOT": 0.4, "MOVE_SYSTEM_TSMMC": 0.2, "MOVE_SEGMENT_REGROWTH": 0.4},
             "B 1 0 4\n", chains=[(1, "AABBAABB"), (5, "A"), (5, "AB")])
    simulation, printed = _build(directory)
    assert "MOVE_SYSTEM_TSMMC excursions run only" not in printed
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        with _quiet():
            summary = simulation.check_moveset_applicability(simulation.keyword_lookup)
    finally:
        os.chdir(cwd)
    assert summary["null_fractions"] == {"MOVE_CHAIN_PIVOT": pytest.approx(0.4)}


def test_the_fixed_conformation_warning_says_it_is_about_the_free_chains(tmp_path: Path) -> None:
    _run_dir(tmp_path / "run", {"MOVE_VMMC": 0.5, "MOVE_SEGMENT_REGROWTH": 0.5}, "B 1 0 4\n")
    _, printed = _build(tmp_path / "run")
    assert "contains no move that can change a chain's shape" in printed
    assert "MOVE_SEGMENT_REGROWTH does reshape the free beads of the partly frozen chains" in printed


@pytest.mark.parametrize("freeze_line, max_beads, warned", [("B 1 0 4", 3, False), ("B 1 0 5", 3, True),
                                                            ("B 1 0 5", 4, False)])
def test_the_trap_warning_names_only_stretches_longer_than_a_window(tmp_path: Path, freeze_line: str,
                                                                    max_beads: int, warned: bool) -> None:
    _run_dir(tmp_path / "run", {"MOVE_CRANKSHAFT": 0.5, "MOVE_SEGMENT_REGROWTH": 0.5}, freeze_line + "\n",
             extra={"SEGMENT_REGROWTH_MAX_BEADS": max_beads})
    _, printed = _build(tmp_path / "run")
    assert ("held at both ends by frozen beads" in printed) is warned
    if warned:
        assert "have a stretch of more than %d free beads" % max_beads in printed
        assert "chain 1 beads 1-4" in printed


def test_a_continued_run_reproduces_the_uninterrupted_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    freeze = tmp_path / "freeze.txt"
    freeze.write_text("B 1 0 4\nB 12 1\n")
    extra = ("EXPERIMENTAL_FEATURES : True", "FREEZE_FILE : %s" % freeze)
    ref, cont, ref_sim, cont_sim = _stop_and_resume(str(tmp_path), 30, 12, monkeypatch,
                                                    moves={"MOVE_CRANKSHAFT": 0.4, "MOVE_SEGMENT_REGROWTH": 0.6},
                                                    extra=extra)
    _assert_same_continuation(ref, cont, ref_sim, cont_sim, 12)
    assert np.loadtxt(os.path.join(cont, "MOVE_FREQS.dat"))[-1][15] > 0


def test_with_the_move_off_every_move_is_drawn_as_before() -> None:
    """With the move at 0 (or absent), each draw is one uniform number mapped onto
    the fourteen older moves' intervals, laid out in code order as they were
    before code 15 existed - checked against that layout built independently here."""
    rng = np.random.default_rng(3)
    for _ in range(SELECTOR_MOVE_SETS):
        weights = rng.random(len(OLDER_MOVE_KEYWORDS)) * (rng.random(len(OLDER_MOVE_KEYWORDS)) < 0.6)
        weights[int(rng.integers(len(weights)))] += 0.1
        fractions = [float(w) for w in weights / weights.sum()]
        upper = list(itertools.accumulate(fractions))
        last_used = max(i for i, fraction in enumerate(fractions) if fraction > 0)
        upper[last_used] = 1.0
        for lookup in (dict(zip(OLDER_MOVE_KEYWORDS, fractions)),
                       {**dict(zip(OLDER_MOVE_KEYWORDS, fractions)), "MOVE_SEGMENT_REGROWTH": 0.0}):
            selector = AcceptanceCalculator(40.0, lookup)
            random.seed(13)
            drawn = [selector.move_selector() for _ in range(SELECTOR_DRAWS)]
            after = random.getstate()
            random.seed(13)
            uniforms = [random.random() for _ in range(SELECTOR_DRAWS)]
            assert random.getstate() == after, "the selector must draw exactly one number per move"
            expected = [next(code for code, top in enumerate(upper, start=1) if u < top) for u in uniforms]
            assert drawn == expected


def test_the_move_is_drawn_at_its_fraction() -> None:
    acceptance = _acceptance(40.0, MOVE_CRANKSHAFT=0.7, MOVE_SEGMENT_REGROWTH=0.3)
    random.seed(17)
    draws = [acceptance.move_selector() for _ in range(20000)]
    assert set(draws) == {1, 15}
    assert draws.count(15) / len(draws) == pytest.approx(0.3, abs=0.015)
