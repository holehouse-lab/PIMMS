"""
Freezing single beads of a chain (``FREEZE_FILE`` ``B`` lines).

A ``B`` line freezes some beads of one chain. The crankshaft, which moves one
bead per sub-move, is the only move that acts on such a chain: its free beads
are drawn for moves like any other bead and its frozen beads never are, while
every move that displaces a whole chain or a segment of one (translation,
rotation, the pivots, slither, pull, the cluster moves, VMMC, the TSMMC moves,
jump-and-relax) leaves the chain alone.

The tests check the reader, the per-bead selection the serial and the parallel
crankshaft kernels are handed, that frozen beads stay put and free beads move in
real runs that use every move, that a chain frozen bead by bead is the same run
as one on a ``C`` line, and what happens when every chain is partly frozen. The
decisive one checks that the crankshaft samples the exact Boltzmann distribution
of the free beads given the frozen ones, enumerated with the model energy of
``test_energy_oracle`` (which recomputes the energy without any PIMMS energy
routine), with a positive control showing the comparison can fail.
"""

from __future__ import annotations

import contextlib
import io
import itertools
import os
import random
import types
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import pytest

from pimms import CONFIG, crankshaft_list_functions, mega_crank_fast
from pimms.data_structures import FreezeFile, contiguous_runs, format_bead_runs
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import KeyFileException, MoveException, SimulationException
from pimms.lemonade import load
from pimms.moves import MoveObject
from pimms.simulation import CHAIN_TARGETED_MOVE_CODES, Simulation
from pimms.tests import kernel_test_utils as U
from pimms.tests.test_energy_oracle import (ChainSpec, build_hamiltonian, build_lattice, make_forcefield,
                                            oracle_energy)

# A 3D system for the end-to-end runs: chainIDs 1-6 are AABBAABB, 7-16 monomers,
# 17-20 BBBB. Chain 1 has its first four beads frozen, chain 2 two isolated
# interior beads, and chain 17 is frozen whole.
RUN_BOX = [24, 24, 24]
RUN_CHAINS = [(6, "AABBAABB"), (10, "A"), (4, "BBBB")]
RUN_FREEZE_FILE = "B 1 0-3\nB 2 2 5\nC 17\n"
RUN_FROZEN_BEADS = {1: (0, 1, 2, 3), 2: (2, 5), 17: (0, 1, 2, 3)}

# Every move type at once, so a move that ignored a partly frozen chain's frozen
# beads would show up in the trajectory
EVERY_MOVE = {"MOVE_CRANKSHAFT": 0.3, "MOVE_CHAIN_TRANSLATE": 0.1, "MOVE_CHAIN_ROTATE": 0.1,
              "MOVE_CHAIN_PIVOT": 0.1, "MOVE_HEAD_PIVOT": 0.05, "MOVE_SLITHER": 0.1, "MOVE_PULL": 0.05,
              "MOVE_CLUSTER_TRANSLATE": 0.05, "MOVE_CLUSTER_ROTATE": 0.05, "MOVE_CTSMMC": 0.03,
              "MOVE_MULTICHAIN_TSMMC": 0.02, "MOVE_JUMP_AND_RELAX": 0.03, "MOVE_VMMC": 0.02}
# their move codes (every code but 12, the system-wide TSMMC)
EVERY_MOVE_CODES = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14)
# long enough for every move in EVERY_MOVE to be drawn several times
EVERY_MOVE_STEPS = 400

# The exact-equilibrium fixture: a 2D hardwall box, a six-bead chain with some
# beads frozen and a frozen two-bead obstacle next to its tail. With beads 0, 1
# and 3 frozen, bead 2 is a free bead between two frozen ones and beads 4-5 a
# free tail (55 states); with beads 1 and 3 frozen, bead 0 is also a free chain
# end anchored on a frozen bead (314 states). At this temperature the Boltzmann
# weights are far enough from uniform for a 1.3x error in beta to be resolved
# with the samples below (measured: the control clears its bar 2.2-4.5x).
EXACT_TEMPERATURE = 8.0
EXACT_FORCEFIELD_SEED = 20261109
EXACT_SUBSTEPS = {False: 20, True: 2000}
EXACT_BOX = {False: [7, 7], True: [24, 24]}
EXACT_ORIGIN = {False: (0, 0), True: (10, 10)}
# (frozen beads, parallel kernel, samples); the parallel cases are slow
EXACT_CASES = [((0, 1, 3), False, 40000), ((1, 3), False, 40000),
               ((0, 1, 3), True, 40000), ((1, 3), True, 60000)]
# the beta error the positive control applies
EXACT_CONTROL_BETA_SCALE = 1.3
# the chi-square must stay below dof + K sigma, and the control must clear twice that
EXACT_SIGMA = 5.0
EXACT_CONTROL_MARGIN = 2.0


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _in(path: Path) -> Iterator[None]:
    """Run in ``path`` with standard output discarded.

    Parameters
    ----------
    path : pathlib.Path
        The directory.

    Yields
    ------
    None
    """
    cwd = os.getcwd()
    os.chdir(path)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            yield
    finally:
        os.chdir(cwd)


def _freeze_file(tmp_path: Path, text: str) -> FreezeFile:
    """Write a freeze file and read it.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Where to write it.
    text : str
        Its contents.

    Returns
    -------
    FreezeFile
        The file as PIMMS reads it.
    """
    path = tmp_path / "freeze.txt"
    path.write_text(text)
    return FreezeFile(str(path))


def _write_run(directory: Path, moves: dict[str, float], freeze_text: str, *, chains: Any = None,
               parallel: bool = False, n_steps: int = 60, seed: int = 5) -> None:
    """Write a 3D SR keyfile with a freeze file, saving a frame every step.

    Parameters
    ----------
    directory : pathlib.Path
        The run directory (created).
    moves : dict of str to float
        MOVE_* fractions.
    freeze_text : str
        The freeze file.
    chains : list of (int, str), optional
        CHAIN lines; RUN_CHAINS by default.
    parallel : bool, optional
        Run the megamoves on the parallel kernels.
    n_steps : int, optional
        Steps to run.
    seed : int, optional
        SEED.

    Returns
    -------
    None
    """
    directory.mkdir(parents=True, exist_ok=True)
    U.write_param_file(str(directory / "params.prm"), "SR")
    extra: dict[str, Any] = {"FREEZE_FILE": "freeze.txt", "ENERGY_CHECK": 5, "XTC_FREQ": 1, "EN_FREQ": 1,
                             "ANALYSIS_FREQ": 10**9, "RESTART_FREQ": 10**9, "CRANKSHAFT_SUBSTEPS": 400,
                             "SAVE_EQ": True}
    if parallel:
        extra.update({"PARALLELIZE": True, "PARALLEL_THREADS": 2})
    U.write_keyfile(str(directory / "KEYFILE.kf"), 3, False, moves, box=RUN_BOX,
                    chains=chains if chains is not None else RUN_CHAINS, n_steps=n_steps, equilibration=2,
                    seed=seed, temperature=40, extra=extra)
    (directory / "freeze.txt").write_text(freeze_text)


def _run(directory: Path) -> tuple[Simulation, str]:
    """Build and run the simulation in ``directory``.

    Parameters
    ----------
    directory : pathlib.Path
        The run directory.

    Returns
    -------
    tuple of (Simulation, str)
        The finished simulation and what it printed.
    """
    printed = io.StringIO()
    cwd = os.getcwd()
    os.chdir(directory)
    try:
        with contextlib.redirect_stdout(printed):
            simulation = Simulation(KeyFileParser("KEYFILE.kf").keyword_lookup)
            simulation.run_simulation()
    finally:
        os.chdir(cwd)
    return simulation, printed.getvalue()


def _frames(directory: Path) -> np.ndarray:
    """Every saved frame of a run, as lattice coordinates.

    Parameters
    ----------
    directory : pathlib.Path
        The run directory.

    Returns
    -------
    numpy.ndarray
        Shape (frames, beads, 3), beads in bead-table order (ascending chainID).
    """
    with _in(directory):
        store = load("traj.xtc", "START.pdb", "KEYFILE.kf")
    return np.asarray(store.positions)


def _rows_of(simulation: Simulation, frozen_beads: dict[int, tuple[int, ...]]) -> np.ndarray:
    """Mark the frozen beads of a run, in bead-table order.

    Parameters
    ----------
    simulation : Simulation
        The run.
    frozen_beads : dict of int to tuple of int
        chainID -> indices of its frozen beads.

    Returns
    -------
    numpy.ndarray
        Boolean, one entry per bead.
    """
    layout = crankshaft_list_functions.chain_layout(simulation.LATTICE)
    frozen = np.zeros(layout.n_beads, dtype=bool)
    for position, chain_id in enumerate(layout.sorted_ids):
        for bead in frozen_beads.get(chain_id, ()):
            frozen[int(layout.offset[position]) + bead] = True
    return frozen


# ---------------------------------------------------------------------------
# the reader
# ---------------------------------------------------------------------------

def test_bead_lines_read_single_beads_and_ranges(tmp_path: Path) -> None:
    """Ranges and single beads are merged per chain, across lines."""
    freeze = _freeze_file(tmp_path, "# tails\nB 3 0-4 7\nB 3 5   # one more\nB 1 2\nC 2\n")
    assert freeze.chains == [2]
    assert freeze.beads == [(1, 2, 2), (3, 0, 5), (3, 7, 7)]


def test_resolve_gives_the_frozen_chains_and_the_beads_of_partly_frozen_ones(tmp_path: Path) -> None:
    """B lines covering a whole chain make it a C chain; a B line on a C chain changes nothing."""
    freeze = _freeze_file(tmp_path, "C 2\nB 3 0-4 7\nB 1 0-2\nB 1 3-4\nB 2 1\n")
    frozen_chains, frozen_beads = freeze.resolve({1: 5, 2: 4, 3: 10})
    assert frozen_chains == [1, 2, 3]
    assert frozen_beads == {3: (0, 1, 2, 3, 4, 7)}


@pytest.mark.parametrize("line, message", [
    ("B 1\n", "names chain 1 but no beads"),
    ("B 1 4-2\n", "runs backwards"),
    ("B 1 -2\n", "neither a bead index"),
    ("B 1 +2\n", "neither a bead index"),
    ("B 1 2.0\n", "neither a bead index"),
    ("B 1 1_0\n", "neither a bead index"),
    ("B 1 0-\n", "neither a bead index"),
    ("B x 2\n", "starts with the ID"),
    ("B -1 2\n", "starts with the ID"),
])
def test_malformed_bead_lines_are_refused(tmp_path: Path, line: str, message: str) -> None:
    """Every malformed B line is refused with a message naming the problem."""
    with pytest.raises(KeyFileException, match=message):
        _freeze_file(tmp_path, line)


def test_a_bead_past_the_end_of_its_chain_or_an_absent_chain_is_refused(tmp_path: Path) -> None:
    """Indices count from 0, so bead N of an N-bead chain does not exist."""
    with pytest.raises(KeyFileException, match=r"freezes bead 5 of chain 1, which has 5 bead\(s\), numbered 0 to 4"):
        _freeze_file(tmp_path, "B 1 0-5\n").resolve({1: 5})
    with pytest.raises(KeyFileException, match="chain 9, which is NOT present"):
        _freeze_file(tmp_path, "B 9 0\n").resolve({1: 5})
    # a mistyped huge range is refused, not expanded
    with pytest.raises(KeyFileException, match="freezes bead 99999999999"):
        _freeze_file(tmp_path, "B 1 0-99999999999\n").resolve({1: 5})


def test_validate_checks_bead_lines_against_the_lattice(tmp_path: Path) -> None:
    """validate_freeze_file reads the chain lengths from the lattice."""
    lattice = types.SimpleNamespace(chains={1: types.SimpleNamespace(positions=[[0, 0, 0]] * 4),
                                            2: types.SimpleNamespace(positions=[[1, 0, 0]])})
    _freeze_file(tmp_path, "B 1 0-3\nC 2\n").validate_freeze_file(lattice)
    with pytest.raises(KeyFileException, match="freezes bead 4 of chain 1"):
        _freeze_file(tmp_path, "B 1 4\n").validate_freeze_file(lattice)
    with pytest.raises(KeyFileException, match="chain 3, which is NOT present"):
        _freeze_file(tmp_path, "B 3 0\n").validate_freeze_file(lattice)


def test_bead_runs_are_written_compactly() -> None:
    """The helpers used for logs, messages and pdb2pimms."""
    assert contiguous_runs([7, 0, 1, 2, 2, 5, 6]) == [(0, 2), (5, 7)]
    assert contiguous_runs([]) == []
    assert format_bead_runs([0, 1, 2, 9]) == "0-2, 9"
    assert format_bead_runs([0, 1, 2, 9], " ") == "0-2 9"


# ---------------------------------------------------------------------------
# which beads the crankshaft may draw
# ---------------------------------------------------------------------------

def _lattice(tmp_path: Path) -> Any:
    """A small 3D lattice built from a keyfile.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory.

    Returns
    -------
    Lattice
        Its chains are those of RUN_CHAINS.
    """
    (tmp_path / "state").mkdir()
    with contextlib.redirect_stdout(io.StringIO()):
        state = U.build_state(str(tmp_path / "state"), 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                              box=RUN_BOX, chains=RUN_CHAINS, seed=3)
    return state.lattice


def test_frozen_bead_rows_marks_exactly_the_frozen_beads(tmp_path: Path) -> None:
    """Rows follow ascending chainID, each chain from its first bead."""
    lattice = _lattice(tmp_path)
    rows = crankshaft_list_functions.frozen_bead_rows(lattice, [1, 2, 17], {1: (0, 1, 2, 3), 2: (2, 5)})
    expected = set()
    for chain_id, beads in RUN_FROZEN_BEADS.items():
        first = lattice.chain_to_firstbead_lookup[chain_id]
        expected.update(first + bead for bead in beads)
    assert set(np.flatnonzero(rows).tolist()) == expected


def test_the_selector_draws_every_free_bead_and_no_frozen_one(tmp_path: Path) -> None:
    """Partly frozen chains keep their free beads in the draw."""
    lattice = _lattice(tmp_path)
    n_beads = len(lattice.crankshaft_lists)
    rows = crankshaft_list_functions.frozen_bead_rows(lattice, [1, 2, 17], {1: (0, 1, 2, 3), 2: (2, 5)})
    np.random.seed(4)
    drawn = crankshaft_list_functions.bead_selector_constructor(
        n_beads, 50000, lattice, frozen_chains=[1, 2, 17], frozen_beads={1: (0, 1, 2, 3), 2: (2, 5)})
    assert set(drawn.tolist()) == set(np.flatnonzero(~rows).tolist())


def test_freezing_every_bead_of_a_chain_draws_like_freezing_the_chain(tmp_path: Path) -> None:
    """Same rows, same random stream: B lines covering a chain are a C line."""
    lattice = _lattice(tmp_path)
    n_beads = len(lattice.crankshaft_lists)
    np.random.seed(8)
    whole = crankshaft_list_functions.bead_selector_constructor(n_beads, 5000, lattice, frozen_chains=[17])
    np.random.seed(8)
    by_bead = crankshaft_list_functions.bead_selector_constructor(n_beads, 5000, lattice, frozen_chains=[17],
                                                                  frozen_beads={17: (0, 1, 2, 3)})
    np.testing.assert_array_equal(whole, by_bead)


def test_inconsistent_frozen_beads_are_a_bug_and_say_so(tmp_path: Path) -> None:
    """Frozen beads must belong to a frozen chain, inside it."""
    lattice = _lattice(tmp_path)
    with pytest.raises(MoveException, match="not among the frozen chains"):
        crankshaft_list_functions.frozen_bead_rows(lattice, [2], {1: (0,)})
    with pytest.raises(MoveException, match="without their chains"):
        crankshaft_list_functions.frozen_bead_rows(lattice, [], {1: (0,)})
    with pytest.raises(MoveException, match="outside the chain"):
        crankshaft_list_functions.frozen_bead_rows(lattice, [1], {1: (8,)})


# ---------------------------------------------------------------------------
# real runs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("parallel", [False, True], ids=["serial", "parallel"])
def test_frozen_beads_never_move_and_free_beads_do(tmp_path: Path, parallel: bool) -> None:
    """Every move type, a frame every step, and the energy checked every 5 steps."""
    if parallel:
        # the parallel crankshaft must really be dispatched, or this is the serial test again
        layout = mega_crank_fast.parallel_crank_layout_info(*RUN_BOX, False)
        assert layout["num_blocks"] > 1
    _write_run(tmp_path / "run", EVERY_MOVE, RUN_FREEZE_FILE, parallel=parallel, n_steps=EVERY_MOVE_STEPS)
    simulation, _ = _run(tmp_path / "run")
    assert simulation.frozen_chains == [1, 2, 17]
    assert simulation.frozen_beads == {1: (0, 1, 2, 3), 2: (2, 5)}
    for code in EVERY_MOVE_CODES:
        assert simulation.ACC.move_count[code] > 0, "move %d was never attempted" % code

    frames = _frames(tmp_path / "run")
    assert frames.shape[0] == EVERY_MOVE_STEPS + 1
    moved = np.any(frames != frames[0], axis=2).any(axis=0)
    frozen = _rows_of(simulation, RUN_FROZEN_BEADS)
    assert not moved[frozen].any(), "a frozen bead moved"

    # the free beads of the two partly frozen chains moved, and so did the free chains
    free_of_partly_frozen = _rows_of(simulation, {1: (4, 5, 6, 7), 2: (0, 1, 3, 4, 6, 7)})
    assert moved[free_of_partly_frozen].all()
    assert moved[~(frozen | free_of_partly_frozen)].mean() > 0.9


def test_a_chain_frozen_bead_by_bead_is_the_same_run_as_one_on_a_c_line(tmp_path: Path) -> None:
    """B lines that cover chain 17 give the C 17 run, step for step."""
    for name, text in (("by_chain", "B 1 0-3\nC 17\n"), ("by_bead", "B 1 0-3\nB 17 0-1\nB 17 2-3\n")):
        _write_run(tmp_path / name, EVERY_MOVE, text, n_steps=30)
        _run(tmp_path / name)
    assert (tmp_path / "by_chain" / "ENERGY.dat").read_text() == (tmp_path / "by_bead" / "ENERGY.dat").read_text()
    np.testing.assert_array_equal(_frames(tmp_path / "by_chain"), _frames(tmp_path / "by_bead"))


def test_a_system_whose_every_chain_is_partly_frozen_runs_on_the_crankshaft(tmp_path: Path) -> None:
    """Every move but the crankshaft is a null move, and the chain-picking ones are logged as rejected attempts."""
    # the fractions must add up to 1: the system-wide TSMMC takes 0.02 from the crankshaft
    moves = dict(EVERY_MOVE, MOVE_CRANKSHAFT=EVERY_MOVE["MOVE_CRANKSHAFT"] - 0.02, MOVE_SYSTEM_TSMMC=0.02)
    _write_run(tmp_path / "run", moves, "B 1 0\nB 2 5\nB 3 2-3\n", chains=[(3, "AABBAA")],
               n_steps=EVERY_MOVE_STEPS)
    simulation, printed = _run(tmp_path / "run")
    printed = " ".join(printed.split())       # the screen output is wrapped
    assert "Every chain in the system has frozen beads" in printed
    assert "cannot act on a chain with a frozen bead" in printed

    # the start-up summary (printed by the PIMMS script) lists the B lines
    summary = io.StringIO()
    with _in(tmp_path / "run"), contextlib.redirect_stdout(summary):
        KeyFileParser("KEYFILE.kf").print_summary()
    assert "Beads to freeze  : chain 1, bead 0\n" in summary.getvalue()
    assert "Beads to freeze  : chain 3, beads 2-3\n" in summary.getvalue()

    # every move that needs a chain was drawn and never accepted, and neither the
    # whole-chain megamoves nor the multichain TSMMC ever accepted anything
    for code in CHAIN_TARGETED_MOVE_CODES:
        assert simulation.ACC.move_count[code] > 0, "move %d was never drawn" % code
        assert simulation.ACC.accepted_count[code] == 0
    for code in (6, 10, 11):
        assert simulation.ACC.accepted_count[code] == 0

    frames = _frames(tmp_path / "run")
    moved = np.any(frames != frames[0], axis=2).any(axis=0)
    frozen = _rows_of(simulation, {1: (0,), 2: (5,), 3: (2, 3)})
    assert not moved[frozen].any()
    assert moved[~frozen].all()


def test_a_run_whose_every_chain_is_partly_frozen_needs_the_crankshaft(tmp_path: Path) -> None:
    """Without it nothing can move, so the run is refused."""
    _write_run(tmp_path / "run", {"MOVE_CHAIN_TRANSLATE": 1.0}, "B 1 0\nB 2 5\nB 3 2-3\n",
               chains=[(3, "AABBAA")])
    with pytest.raises(SimulationException, match="enables neither MOVE_CRANKSHAFT"):
        _run(tmp_path / "run")


def test_partly_frozen_chains_without_the_crankshaft_are_warned_about(tmp_path: Path) -> None:
    """Free chains keep the run going, but the free beads of chain 1 would never move."""
    _write_run(tmp_path / "run", {"MOVE_CHAIN_TRANSLATE": 1.0}, "B 1 0-3\n", n_steps=5)
    _, printed = _run(tmp_path / "run")
    assert "1 chain(s) are partly frozen (FREEZE_FILE B lines)" in " ".join(printed.split())


# ---------------------------------------------------------------------------
# the exact equilibrium of the free beads
# ---------------------------------------------------------------------------

def _neighbours(site: tuple[int, int], dims: list[int]) -> Iterator[tuple[int, int]]:
    """The Moore (Chebyshev-1) neighbours of a site inside a hardwall box.

    Parameters
    ----------
    site : tuple of int
        The site.
    dims : list of int
        The box.

    Yields
    ------
    tuple of int
        Each neighbour inside the box.
    """
    for step in itertools.product((-1, 0, 1), repeat=2):
        if step == (0, 0):
            continue
        neighbour = (site[0] + step[0], site[1] + step[1])
        if all(0 <= neighbour[axis] < dims[axis] for axis in range(2)):
            yield neighbour


def _exact_system(dims: list[int], origin: tuple[int, int]) -> list[ChainSpec]:
    """The fixture's chain and obstacle, placed at ``origin``.

    Parameters
    ----------
    dims : list of int
        The box (unused beyond documenting where the system sits).
    origin : tuple of int
        Offset of the system in the box.

    Returns
    -------
    list of (str, list of list of int)
        The partly frozen chain, then the frozen obstacle.
    """
    x, y = origin
    chain = ("ABCABC", [[x + 1, y + 1], [x + 2, y + 1], [x + 3, y + 1], [x + 3, y + 2], [x + 3, y + 3],
                        [x + 2, y + 3]])
    obstacle = ("CA", [[x + 4, y + 3], [x + 4, y + 4]])
    return [chain, obstacle]


def _free_states(system: list[ChainSpec], dims: list[int],
                 frozen: tuple[int, ...]) -> list[tuple[tuple[int, int], ...]]:
    """Every placement of the chain's free beads that its frozen beads allow.

    The beads are placed in chain order: each free bead on an empty site next to
    the bead before it (bead 0 next to bead 1, which must then be frozen), and
    next to the bead after it when that one is frozen.

    Parameters
    ----------
    system : list of (str, list of list of int)
        The fixture.
    dims : list of int
        The box.
    frozen : tuple of int
        The chain's frozen beads.

    Returns
    -------
    list of tuple
        The sites of the free beads, in chain order, for each state.
    """
    positions = [tuple(site) for site in system[0][1]]
    assert 0 in frozen or 1 in frozen, "bead 0 is placed next to bead 1, which must be frozen"
    occupied = {positions[bead] for bead in frozen} | {tuple(site) for site in system[1][1]}
    states: list[tuple[tuple[int, int], ...]] = []

    def place(bead: int, placed: list[tuple[int, int]]) -> None:
        if bead == len(positions):
            states.append(tuple(site for index, site in enumerate(placed) if index not in frozen))
            return
        if bead in frozen:
            place(bead + 1, placed + [positions[bead]])
            return
        anchor = placed[bead - 1] if bead > 0 else positions[1]
        following = positions[bead + 1] if bead + 1 in frozen else None
        for site in _neighbours(anchor, dims):
            if site in occupied or site in placed:
                continue
            if following is not None and max(abs(site[0] - following[0]), abs(site[1] - following[1])) != 1:
                continue
            place(bead + 1, placed + [site])

    place(0, [])
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
        The expected probability of each state.
    effective_samples : float
        The number of independent samples the counts are worth.

    Returns
    -------
    float
        effective_samples x sum((observed - expected)^2 / expected), over
        frequencies.
    """
    observed = counts / counts.sum()
    return float(effective_samples * np.sum((observed - probabilities) ** 2 / probabilities))


@pytest.mark.parametrize("frozen, parallel, n_samples",
                         [pytest.param(*case, marks=[pytest.mark.slow] if case[1] else [],
                                       id="%s-frozen-%s" % ("parallel" if case[1] else "serial",
                                                            "".join(str(bead) for bead in case[0])))
                          for case in EXACT_CASES])
def test_the_crankshaft_samples_the_exact_distribution_of_the_free_beads(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch, frozen: tuple[int, ...], parallel: bool,
        n_samples: int) -> None:
    """Sampled state frequencies of the free beads match the enumerated Boltzmann weights."""
    monkeypatch.chdir(tmp_path)
    dims = EXACT_BOX[parallel]
    if parallel:
        assert mega_crank_fast.parallel_crank_layout_info(dims[0], dims[1], 1, False)["num_blocks"] > 1
    forcefield = make_forcefield(np.random.default_rng(EXACT_FORCEFIELD_SEED), "SR", "abs", EXACT_TEMPERATURE)
    system = _exact_system(dims, EXACT_ORIGIN[parallel])
    states = _free_states(system, dims, frozen)
    sequence, positions = system[0]
    free = [bead for bead in range(len(positions)) if bead not in frozen]
    energies = []
    for state in states:
        placed = [list(site) for site in positions]
        for bead, site in zip(free, state):
            placed[bead] = list(site)
        energies.append(oracle_energy([(sequence, placed), system[1]], dims, True, forcefield).terms()[0])
    beta = CONFIG.INVTEMP_FACTOR / EXACT_TEMPERATURE
    exact = _boltzmann(np.asarray(energies, dtype=float), beta)

    with contextlib.redirect_stdout(io.StringIO()):
        hamiltonian = build_hamiltonian(forcefield, 2, True)
        lattice = build_lattice(hamiltonian, dims, True, system)
    acceptance = types.SimpleNamespace(invtemp=beta)
    mover = MoveObject()
    energy = int(hamiltonian.evaluate_total_energy(lattice)[0])
    random.seed(5)
    np.random.seed(5)

    index = {state: i for i, state in enumerate(states)}
    counts = np.zeros(len(states))
    trace = np.zeros(n_samples)
    held = [list(map(int, positions[bead])) for bead in frozen]
    for sample in range(n_samples):
        lattice, energy, _, _ = mover.system_shake(lattice, energy, acceptance, hamiltonian,
                                                   EXACT_SUBSTEPS[parallel], "UNIFORM", hardwall=True,
                                                   frozen_chains=[1, 2], frozen_beads={1: frozen},
                                                   parallelize=parallel, num_threads=2)
        now = lattice.chains[1].get_ordered_positions()
        assert [list(map(int, now[bead])) for bead in frozen] == held
        counts[index[tuple(tuple(map(int, now[bead])) for bead in free)]] += 1
        trace[sample] = energy

    # the tracked energy is the model energy of where the kernel left the beads
    final = [(sequence, [list(map(int, p)) for p in lattice.chains[1].get_ordered_positions()]), system[1]]
    assert energy == oracle_energy(final, dims, True, forcefield).terms()[0]

    effective = n_samples / (2.0 * U.integrated_autocorr_time(trace))
    dof = len(states) - 1
    threshold = dof + EXACT_SIGMA * np.sqrt(2.0 * dof)
    assert _chi_square(counts, exact, effective) < threshold
    assert counts[exact > 1e-3].min() > 0, "a reachable state was never sampled"
    assert len(states) == {(0, 1, 3): 55, (1, 3): 314}[frozen]

    # positive control: the same samples against the distribution at a beta that
    # is 30% off must fail by a wide margin, or the comparison could not fail
    wrong = _boltzmann(np.asarray(energies, dtype=float), EXACT_CONTROL_BETA_SCALE * beta)
    assert _chi_square(counts, wrong, effective) > EXACT_CONTROL_MARGIN * threshold
