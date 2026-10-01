"""
Regression tests for the deep-audit fixes to the Python-level move code, the
start-up chain placement and the whole-system energy evaluation.

Every fix here is a performance or memory fix that must not change a single
result, so each one is pinned twice:

* **the result is unchanged** - compared with a literal copy of the routine it
  replaces (kept in this file), or with a value worked out from the definition
  (the model energy from ``test_energy_oracle``'s plain-numpy oracle, the draws
  written out by hand);
* **the cost is gone** - by counting the work the old code did (calls to a
  per-site lookup, passes over a container, whole-grid comparisons, traced
  memory), never by timing, so the tests are deterministic.

The second kind of assertion is the one that fails on the code as it was.

Covered: the VMMC neighbour-energy scan (numpy instead of a Python triple loop,
including the HARDWALL bounds on the first and last layer of every axis), the
jump of the jump-and-relax move (local energy change instead of two
whole-system energies), ``build_all_envelope_pairs`` (per-chunk folding),
``get_empty_site`` (no whole-grid scan per chain), ``Lattice.get_random_chain``
(cached candidate list), the parallel pass's frozen mask, and the coordinates
the Python moves leave in ``Chain.positions`` (Python ints).
"""

from __future__ import annotations

import random
import tracemalloc
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

import numpy as np
import numpy.typing as npt
import pytest

from pimms import inner_loops, inner_loops_hardwall, lattice, lattice_utils, moves
from pimms.CONFIG import NP_INT_TYPE
from pimms.latticeExceptions import LatticeUtilsException
from pimms.tests import kernel_test_utils as U
from pimms.tests.test_energy_oracle import ForceField, lattice_chains, oracle_energy

IntArray = npt.NDArray[Any]


# ---------------------------------------------------------------------------
# VMMC neighbour energies
# ---------------------------------------------------------------------------

# box 7 with 8-mers / box 9 with 10-mers: a chain can lie against its own periodic
# image, which is where a local envelope and a whole-system sum could disagree
JUMP_CASES = [
    (3, "SLR", False, [7, 7, 7], [(4, "AABBAABB"), (6, "A"), (3, "BAB")]),
    (3, "SLR", True, [8, 7, 9], [(4, "AABBAABB"), (6, "A"), (3, "BAB")]),
    (3, "LR", False, [9, 8, 7], [(5, "ABBAB"), (5, "B"), (3, "AAAA")]),
    (3, "SR", True, [7, 7, 7], [(5, "ABBAB"), (5, "B"), (3, "AAAA")]),
    (2, "SLR", False, [9, 11], [(3, "AABBAABBAA"), (5, "A"), (3, "BAB")]),
    (2, "SLR", True, [10, 9], [(3, "AABBAAB"), (5, "A"), (3, "BAB")]),
    (2, "LR", False, [11, 8], [(4, "ABBAB"), (4, "B"), (2, "AAAA")]),
    (2, "SR", True, [9, 9], [(4, "ABBAB"), (4, "B"), (2, "AAAA")]),
]


def _reference_vmmc_neighbour_energies(
    grid: IntArray,
    tg: IntArray,
    SRT: IntArray,
    LRT: IntArray,
    SLRT: IntArray,
    hardwall: bool,
    offsets: dict,
    m_id: int,
    positions: list,
    intcodes: list,
    lr_flags: IntArray,
    offset: list,
    dimensions: list,
) -> dict:
    """The loop ``MoveObject._vmmc_neighbour_energies`` used before the numpy scan.

    A literal copy of the 1.0.8 body (commit 59163ce), with the grids and tables
    passed in rather than read from the lattice and Hamiltonian objects.

    Parameters
    ----------
    grid, tg : numpy.ndarray
        The chainID grid and the type grid.
    SRT, LRT, SLRT : numpy.ndarray
        The short-range, long-range and super-long-range interaction tables.
    hardwall : bool
        Whether a site beyond the box does not exist (True) or wraps (False).
    offsets : dict
        ``{1: shell, 3: shell}`` of ``(offset, chebyshev)`` pairs.
    m_id : int
        chainID of the chain being scanned.
    positions, intcodes, lr_flags : sequence
        The chain's bead positions, type codes and long-range flags.
    offset : list of int
        Translation applied to the positions before scanning.
    dimensions : list of int
        The box.

    Returns
    -------
    dict
        ``{chainID: energy}`` exactly as the loop built it.
    """
    nd = len(dimensions)
    energies: dict = {}

    for b in range(len(positions)):
        t_b = int(intcodes[b])
        is_lr = bool(lr_flags[b])
        rng = 3 if is_lr else 1
        base = [positions[b][d] + offset[d] for d in range(nd)]

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

            j = int(grid[tuple(npos)])
            if j == 0 or j == m_id:
                continue
            t_n = int(tg[tuple(npos)])

            if cheb == 1:
                e = SRT[t_b][t_n]
            elif cheb == 2:
                e = LRT[t_b][t_n] if is_lr else 0.0
            else:
                e = SLRT[t_b][t_n] if is_lr else 0.0

            if e != 0.0:
                energies[j] = energies.get(j, 0.0) + e

    return energies


def _random_tables(rng: np.random.Generator) -> tuple[IntArray, IntArray, IntArray]:
    """Random integer SR / LR / SLR tables over three bead types plus solvent.

    Types 1 and 2 are long-range, type 3 is not, so its LR and SLR rows and
    columns are zero - as the parameter-file parser guarantees.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of the entries.

    Returns
    -------
    tuple of numpy.ndarray
        Symmetric ``(4, 4)`` int32 tables ``(SR, LR, SLR)``; index 0 is solvent.
    """

    def symmetric(low: int, high: int) -> IntArray:
        half = rng.integers(low, high, size=(4, 4))
        return (np.triu(half) + np.triu(half, 1).T).astype(NP_INT_TYPE)

    SRT = symmetric(-9, 6)
    LRT = symmetric(-5, 4)
    SLRT = symmetric(-3, 3)
    SRT[1, 2] = SRT[2, 1] = 0  # one silent short-range pair
    for table in (LRT, SLRT):
        table[0, :] = table[:, 0] = 0
        table[3, :] = table[:, 3] = 0
    return (SRT, LRT, SLRT)


def _random_dense_system(
    rng: np.random.Generator, dims: list[int], n_chains: int, occupancy: float
) -> tuple[IntArray, IntArray]:
    """A random occupied lattice: a chainID grid and a matching type grid.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of the occupancy, the chainIDs and the types.
    dims : list of int
        The box.
    n_chains : int
        chainIDs are drawn from ``1..n_chains``.
    occupancy : float
        Fraction of sites that hold a bead.

    Returns
    -------
    tuple of numpy.ndarray
        ``(grid, type_grid)``, both int32; an empty site is 0 in both.
    """
    occupied = rng.random(dims) < occupancy
    grid = np.where(occupied, rng.integers(1, n_chains + 1, size=dims), 0).astype(
        NP_INT_TYPE
    )
    type_grid = np.where(occupied, rng.integers(1, 4, size=dims), 0).astype(NP_INT_TYPE)
    return (grid, type_grid)


def _vmmc_call(
    grid: IntArray,
    tg: IntArray,
    tables: tuple,
    hardwall: bool,
    m_id: int,
    positions: list,
    intcodes: list,
    lr_flags: IntArray,
    offset: list,
    dims: list[int],
) -> tuple[dict, dict]:
    """Run the routine under test and the reference loop on one input.

    Parameters
    ----------
    grid, tg : numpy.ndarray
        The chainID grid and the type grid.
    tables : tuple of numpy.ndarray
        ``(SR, LR, SLR)`` interaction tables.
    hardwall : bool
        Boundary convention.
    m_id : int
        chainID of the chain being scanned.
    positions, intcodes, lr_flags : sequence
        The chain's bead positions, type codes and long-range flags.
    offset : list of int
        Translation applied before scanning.
    dims : list of int
        The box.

    Returns
    -------
    tuple of dict
        ``(result of MoveObject._vmmc_neighbour_energies, result of the loop)``.
    """
    nd = len(dims)
    offsets = {1: moves._vmmc_offset_shell(nd, 1), 3: moves._vmmc_offset_shell(nd, 3)}
    lat = SimpleNamespace(grid=grid, type_grid=tg)
    ham = SimpleNamespace(
        residue_interaction_table=tables[0],
        LR_residue_interaction_table=tables[1],
        SLR_residue_interaction_table=tables[2],
    )
    arguments = (
        lat,
        ham,
        hardwall,
        offsets,
        m_id,
        positions,
        intcodes,
        lr_flags,
        offset,
        dims,
    )
    got = moves.MoveObject()._vmmc_neighbour_energies(*arguments)
    by_loop = moves.MoveObject()._vmmc_neighbour_energies_loop(*arguments)
    by_numpy = moves.MoveObject()._vmmc_neighbour_energies_numpy(*arguments)
    want = _reference_vmmc_neighbour_energies(
        grid,
        tg,
        tables[0],
        tables[1],
        tables[2],
        hardwall,
        offsets,
        m_id,
        positions,
        intcodes,
        lr_flags,
        offset,
        dims,
    )
    # whichever implementation the dispatcher picked, BOTH must give the
    # reference dictionary, so every comparison made with this helper pins the
    # loop and the numpy scan alike
    _assert_same_dict(by_loop, want, ("loop", m_id, offset))
    _assert_same_dict(by_numpy, want, ("numpy", m_id, offset))
    return (got, want)


def _assert_same_dict(got: dict, want: dict, label: object) -> None:
    """Two neighbour-energy dictionaries agree in keys, order, values and types.

    Parameters
    ----------
    got, want : dict
        The dictionary under test and the reference.
    label : object
        Shown when an assertion fails.

    Returns
    -------
    None
    """
    assert list(got) == list(want), label
    assert [type(k) for k in got] == [type(k) for k in want], label
    for key in want:
        assert got[key] == want[key], (label, key)
        assert type(got[key]) is type(want[key]), (
            label,
            key,
            type(got[key]),
            type(want[key]),
        )


@pytest.mark.parametrize("dims", [[9, 11], [8, 7, 9]], ids=["2D", "3D"])
@pytest.mark.parametrize("hardwall", [False, True], ids=["periodic", "hardwall"])
def test_vmmc_neighbour_energies_equal_the_loop_on_dense_systems(
    dims: list[int], hardwall: bool
) -> None:
    """The numpy scan returns what the loop returned, on random dense boxes.

    Half the sites are occupied, a third of the beads are not long-range, and
    the scanned chain is spread over the whole box, so the first and last layer
    of every axis hold both scanned beads and interacting neighbours.
    """
    rng = np.random.default_rng(1000 * len(dims) + int(hardwall))
    nd = len(dims)
    n_calls = 0
    n_nonempty = 0
    n_with_face_beads = 0
    for _ in range(6):
        (grid, tg) = _random_dense_system(rng, dims, n_chains=4, occupancy=0.5)
        tables = _random_tables(rng)
        for m_id in (1, 2, 4):
            positions = [[int(x) for x in p] for p in np.argwhere(grid == m_id)]
            order = rng.permutation(len(positions))
            positions = [positions[i] for i in order]
            intcodes = [int(tg[tuple(p)]) for p in positions]
            lr_flags = np.array([1 if code in (1, 2) else 0 for code in intcodes])
            n_with_face_beads += all(
                any(p[d] == face for p in positions)
                for d in range(nd)
                for face in (0, dims[d] - 1)
            )
            for offset in (
                [0] * nd,
                [1] * nd,
                [-2] + [3] * (nd - 1),
                [int(v) for v in rng.integers(-3, 4, nd)],
            ):
                (got, want) = _vmmc_call(
                    grid,
                    tg,
                    tables,
                    hardwall,
                    m_id,
                    positions,
                    intcodes,
                    lr_flags,
                    offset,
                    dims,
                )
                _assert_same_dict(got, want, (dims, hardwall, m_id, offset))
                n_calls += 1
                n_nonempty += bool(want)
    assert n_calls == 72 and n_nonempty > 60
    # at least one scanned chain had a bead on BOTH faces of EVERY axis (the
    # hand-built test below pins each face on its own)
    assert n_with_face_beads >= 1


@pytest.mark.parametrize("n_dim", [2, 3])
def test_vmmc_neighbour_energies_hardwall_bounds_on_every_face(n_dim: int) -> None:
    """A wall is closed at exactly coordinate 0 and L - 1, on every axis.

    One scanned bead and one partner bead, placed by hand; the expected energy
    is read straight from the table. A partner ON the first or last layer must
    be counted (the scan may not stop one layer early), and a partner that is
    only reachable through the wall must not be (the scan may not wrap), while
    the same geometry under periodic boundaries does interact.
    """
    dims = [8, 9, 10][:n_dim]
    rng = np.random.default_rng(7)
    (SRT, LRT, SLRT) = _random_tables(rng)
    SRT[1, 2] = SRT[2, 1] = -7
    LRT[1, 2] = LRT[2, 1] = 4
    SLRT[1, 2] = SLRT[2, 1] = -2
    middle = [d // 2 for d in dims]

    def energy_of(
        bead: int, partner: int, axis: int, hardwall: bool, shift: int = 0
    ) -> dict:
        """Scan one long-range bead of chain 1 against one bead of chain 2.

        Parameters
        ----------
        bead, partner : int
            Coordinates of the two beads along ``axis`` (box middle elsewhere).
        axis : int
            The axis under test.
        hardwall : bool
            Boundary convention.
        shift : int
            The bead is stored ``shift`` sites away and moved back by the
            ``offset`` argument, so the virtual translation is exercised too.

        Returns
        -------
        dict
            The routine's result.
        """
        grid = np.zeros(dims, dtype=NP_INT_TYPE)
        tg = np.zeros(dims, dtype=NP_INT_TYPE)
        site = list(middle)
        site[axis] = partner
        grid[tuple(site)] = 2
        tg[tuple(site)] = 2
        stored = list(middle)
        stored[axis] = bead - shift
        offset = [0] * n_dim
        offset[axis] = shift
        (got, want) = _vmmc_call(
            grid,
            tg,
            (SRT, LRT, SLRT),
            hardwall,
            1,
            [stored],
            [1],
            np.array([1]),
            offset,
            dims,
        )
        _assert_same_dict(got, want, (axis, bead, partner, hardwall, shift))
        return got

    for axis in range(n_dim):
        last = dims[axis] - 1
        for shift in (0, 1, -2):
            for hardwall in (True, False):
                # partner on the first / last layer, bead 1, 2 and 3 sites inside
                assert energy_of(1, 0, axis, hardwall, shift) == {2: np.float64(-7)}
                assert energy_of(2, 0, axis, hardwall, shift) == {2: np.float64(4)}
                assert energy_of(3, 0, axis, hardwall, shift) == {2: np.float64(-2)}
                assert energy_of(last - 1, last, axis, hardwall, shift) == {
                    2: np.float64(-7)
                }
                assert energy_of(last - 2, last, axis, hardwall, shift) == {
                    2: np.float64(4)
                }
                assert energy_of(last - 3, last, axis, hardwall, shift) == {
                    2: np.float64(-2)
                }
                # bead on the first / last layer, partner on the layer next to it
                assert energy_of(0, 1, axis, hardwall, shift) == {2: np.float64(-7)}
                assert energy_of(last, last - 1, axis, hardwall, shift) == {
                    2: np.float64(-7)
                }
            # only reachable through the wall: nothing under a hardwall ...
            assert energy_of(0, last, axis, True, shift) == {}
            assert energy_of(last, 0, axis, True, shift) == {}
            assert energy_of(1, last, axis, True, shift) == {}
            assert energy_of(last - 2, 0, axis, True, shift) == {}
            # ... and the wrapped contact under periodic boundaries
            assert energy_of(0, last, axis, False, shift) == {2: np.float64(-7)}
            assert energy_of(last, 0, axis, False, shift) == {2: np.float64(-7)}
            assert energy_of(1, last, axis, False, shift) == {2: np.float64(4)}
            assert energy_of(last - 2, 0, axis, False, shift) == {2: np.float64(-2)}


def test_vmmc_neighbour_energies_large_scans_make_no_per_site_python_lookups(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A large scan no longer calls ``get_gridvalue`` once or twice per shell site.

    That per-site call from a Python triple loop was the cost (2-6 ms per call,
    three calls per chain scanned); above the crossover the numpy scan is used
    and makes none.
    """
    rng = np.random.default_rng(3)
    dims = [8, 8, 8]
    (grid, tg) = _random_dense_system(rng, dims, n_chains=5, occupancy=0.5)
    tables = _random_tables(rng)
    positions = [[int(x) for x in p] for p in np.argwhere(grid == 2)]
    intcodes = [int(tg[tuple(p)]) for p in positions]
    lr_flags = np.array([1 if code in (1, 2) else 0 for code in intcodes])
    assert len(positions) > 20

    calls = []
    real = lattice_utils.get_gridvalue

    def counting(position: Any, lattice_grid: Any) -> Any:
        calls.append(1)
        return real(position, lattice_grid)

    monkeypatch.setattr(moves.lattice_utils, "get_gridvalue", counting)
    offsets = {1: moves._vmmc_offset_shell(3, 1), 3: moves._vmmc_offset_shell(3, 3)}
    lat = SimpleNamespace(grid=grid, type_grid=tg)
    ham = SimpleNamespace(
        residue_interaction_table=tables[0],
        LR_residue_interaction_table=tables[1],
        SLR_residue_interaction_table=tables[2],
    )
    got = moves.MoveObject()._vmmc_neighbour_energies(
        lat, ham, False, offsets, 2, positions, intcodes, lr_flags, [1, 0, -1], dims
    )
    assert len(calls) == 0
    want = _reference_vmmc_neighbour_energies(
        grid,
        tg,
        tables[0],
        tables[1],
        tables[2],
        False,
        offsets,
        2,
        positions,
        intcodes,
        lr_flags,
        [1, 0, -1],
        dims,
    )
    assert want
    _assert_same_dict(got, want, "large scan")


@pytest.mark.parametrize("n_dim", [2, 3])
@pytest.mark.parametrize("hardwall", [False, True], ids=["periodic", "hardwall"])
def test_vmmc_neighbour_energies_small_and_large_scans_agree(
    monkeypatch: pytest.MonkeyPatch, n_dim: int, hardwall: bool
) -> None:
    """Monomers and short chains stay on the loop; either way the answer is the same.

    A monomer with short-range beads visits 8 (2D) or 26 (3D) sites, where the
    fixed cost of the numpy scan is several times the whole loop, so scans of at
    most ``_VMMC_LOOP_MAX_SITES`` sites keep the loop and larger ones use numpy.
    The dispatcher, the loop and the numpy scan are each compared with the
    reference for chain lengths 1 to 20, short-range, long-range and mixed.
    """
    dims = [9, 10, 8][:n_dim]
    rng = np.random.default_rng(50 + n_dim + int(hardwall))
    (grid, tg) = _random_dense_system(rng, dims, n_chains=6, occupancy=0.45)
    tables = _random_tables(rng)
    shell = {1: 3**n_dim - 1, 3: 7**n_dim - 1}

    lookups = []
    real = lattice_utils.get_gridvalue

    def counting(position: Any, lattice_grid: Any) -> Any:
        lookups.append(1)
        return real(position, lattice_grid)

    used = {"loop": 0, "numpy": 0}
    offsets = {
        1: moves._vmmc_offset_shell(n_dim, 1),
        3: moves._vmmc_offset_shell(n_dim, 3),
    }
    lat = SimpleNamespace(grid=grid, type_grid=tg)
    ham = SimpleNamespace(
        residue_interaction_table=tables[0],
        LR_residue_interaction_table=tables[1],
        SLR_residue_interaction_table=tables[2],
    )
    for length in (1, 2, 3, 5, 10, 20):
        for kind in ("SR", "LR", "MIX"):
            start = [int(v) for v in rng.integers(0, 3, n_dim)]
            positions = [[(start[0] + i) % dims[0]] + start[1:] for i in range(length)]
            if hardwall:
                positions = [p for p in positions if p[0] >= start[0]]
            m_id = 3
            intcodes = [
                3 if kind == "SR" else (1 + i % 2 if kind == "LR" else 1 + i % 3)
                for i in range(len(positions))
            ]
            lr_flags = np.array([1 if code in (1, 2) else 0 for code in intcodes])
            for offset in ([0] * n_dim, [1] + [-1] * (n_dim - 1)):
                (got, want) = _vmmc_call(
                    grid,
                    tg,
                    tables,
                    hardwall,
                    m_id,
                    positions,
                    intcodes,
                    lr_flags,
                    offset,
                    dims,
                )
                _assert_same_dict(got, want, (length, kind, offset))

                # which implementation did the dispatcher use?
                sites = sum(shell[3 if flag else 1] for flag in lr_flags)
                monkeypatch.setattr(moves.lattice_utils, "get_gridvalue", counting)
                del lookups[:]
                moves.MoveObject()._vmmc_neighbour_energies(
                    lat,
                    ham,
                    hardwall,
                    offsets,
                    m_id,
                    positions,
                    intcodes,
                    lr_flags,
                    offset,
                    dims,
                )
                monkeypatch.undo()
                if sites <= moves._VMMC_LOOP_MAX_SITES:
                    used["loop"] += 1
                    # periodic: every visited site is looked up at least once
                    assert len(lookups) >= (0 if hardwall else sites), (
                        length,
                        kind,
                        sites,
                    )
                    assert hardwall or len(lookups) > 0
                else:
                    used["numpy"] += 1
                    assert len(lookups) == 0, (length, kind, sites)
    assert used["loop"] >= 4 and used["numpy"] >= 12, used
    # the cases the crossover exists for: a short-range monomer is a loop scan
    assert shell[1] <= moves._VMMC_LOOP_MAX_SITES < 7**3 - 1


def test_vmmc_shell_arrays_are_cached_and_faithful() -> None:
    """The numpy form of a shell is built once and matches the shell tuple."""
    for n_dim in (2, 3):
        for radius in (1, 3):
            shell = moves._vmmc_offset_shell(n_dim, radius)
            (deltas, chebyshev) = moves._vmmc_shell_arrays(shell, n_dim)
            assert deltas.shape == ((2 * radius + 1) ** n_dim - 1, n_dim)
            assert [tuple(row) for row in deltas.tolist()] == [
                delta for (delta, _) in shell
            ]
            assert chebyshev.tolist() == [
                max(abs(v) for v in delta) for (delta, _) in shell
            ]
            (again, _) = moves._vmmc_shell_arrays(shell, n_dim)
            assert again is deltas
    # a shell that could be edited in place is converted afresh, never cached
    editable = [((1, 0), 1), ((0, 2), 2)]
    (first, _) = moves._vmmc_shell_arrays(editable, 2)
    editable[0] = ((-1, 0), 1)
    (second, _) = moves._vmmc_shell_arrays(editable, 2)
    assert first.tolist() == [[1, 0], [0, 2]] and second.tolist() == [[-1, 0], [0, 2]]


# ---------------------------------------------------------------------------
# jump-and-relax: local energy change of the jump
# ---------------------------------------------------------------------------


def _oracle_forcefield(kind: str) -> ForceField:
    """The ``kernel_test_utils`` parameter set, in the form the model oracle reads.

    Parameters
    ----------
    kind : str
        ``"SR"``, ``"LR"`` or ``"SLR"`` - the columns ``write_param_file`` writes.

    Returns
    -------
    ForceField
        Pair, solvation and angle energies keyed by residue letter.
    """
    n_cols = {"SR": 1, "LR": 2, "SLR": 3}[kind]
    sr: dict = {}
    lr: dict = {}
    slr: dict = {}
    for (a, b), values in U._PAIR_VALUES.items():
        sr[(a, b)] = sr[(b, a)] = values[0]
        if n_cols >= 2:
            lr[(a, b)] = lr[(b, a)] = values[1]
            slr[(a, b)] = slr[(b, a)] = values[2] if n_cols == 3 else 0
    for residue, solvation in (("A", -2), ("B", -1)):
        sr[(residue, "0")] = sr[("0", residue)] = solvation
    return ForceField(
        text="",
        sr=sr,
        lr=lr,
        slr=slr,
        long_range=frozenset("AB") if n_cols >= 2 else frozenset(),
        angles={"A": (30, 10, 0), "B": (50, 20, 0)},
        temperature=0.0,
    )


def _model_energy(lat: lattice.Lattice, hardwall: bool, ff: ForceField) -> int:
    """Total energy of the configuration held by the chain objects, from the model.

    Parameters
    ----------
    lat : pimms.lattice.Lattice
        The lattice; only its chains' sequences and positions are read.
    hardwall : bool
        Boundary convention.
    ff : ForceField
        The parameter set.

    Returns
    -------
    int
        SR + LR + SLR + angle energy, computed by ``test_energy_oracle``'s
        plain-numpy oracle (no PIMMS energy routine, pair extractor or table).
    """
    terms = oracle_energy(
        lattice_chains(lat), [int(d) for d in lat.dimensions], hardwall, ff
    )
    return int(terms.sr + terms.lr + terms.slr + terms.angle)


@pytest.mark.parametrize(
    "dim,ff_kind,hardwall,box,chains",
    JUMP_CASES,
    ids=[
        "%dD-%s-%s" % (c[0], c[1], "hardwall" if c[2] else "periodic")
        for c in JUMP_CASES
    ],
)
def test_jump_and_relax_jump_energy_is_local_and_exact(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    dim: int,
    ff_kind: str,
    hardwall: bool,
    box: list[int],
    chains: list,
) -> None:
    """The jump's energy change equals the model's, without a whole-system sum.

    For every jump that survives the hard-sphere test - accepted or rejected -
    the energy change returned by ``commit_single_chain_move`` is compared with
    the difference of the model energy (angles included) before and after the
    translation, both computed by the independent oracle. The tracked energy is
    compared with the model after every move, so a rejected jump is also shown to
    be undone exactly. ``Hamiltonian.evaluate_total_energy`` must not be called:
    that whole-system evaluation per move is what the fix removes.
    """
    state = U.build_state(
        tmp_path,
        dim,
        ff_kind,
        hardwall,
        {"MOVE_JUMP_AND_RELAX": 0.5, "MOVE_CRANKSHAFT": 0.5},
        box=box,
        chains=chains,
        seed=31 + dim,
        temperature=14,
    )
    sim = state.sim
    lat = sim.LATTICE
    ham = sim.Hamiltonian
    ff = _oracle_forcefield(ff_kind)
    energy = state.energy
    assert type(energy) is int and energy == _model_energy(lat, hardwall, ff)

    def forbidden(*args: Any, **kwargs: Any) -> None:
        raise AssertionError(
            "jump_and_relax_move evaluated the energy of the whole system"
        )

    monkeypatch.setattr(ham, "evaluate_total_energy", forbidden)

    jumps: list[tuple[int, int]] = []
    real_commit = moves.commit_single_chain_move

    def checked_commit(
        latticeObject: Any,
        hamiltonianObject: Any,
        move_event: Any,
        chainID: int,
        hardwall_flag: bool,
    ) -> int:
        # on entry the chain objects still hold the pre-jump configuration
        before = _model_energy(latticeObject, hardwall, ff)
        local_dif = real_commit(
            latticeObject, hamiltonianObject, move_event, chainID, hardwall_flag
        )
        after = _model_energy(latticeObject, hardwall, ff)
        assert type(local_dif) is int
        assert move_event.move_type == 2 and hardwall_flag == hardwall
        jumps.append((local_dif, after - before))
        return local_dif

    monkeypatch.setattr(moves, "commit_single_chain_move", checked_commit)

    random.seed(5)
    np.random.seed(5)
    accepted = 0
    for _ in range(160):
        chain = lat.get_random_chain()
        (_, energy, jump_accepted) = sim.MOVER.jump_and_relax_move(
            chain, lat, energy, sim.ACC, ham, 20, sim.CS_mode, hardwall
        )
        accepted += bool(jump_accepted)
        assert type(energy) is int
        assert energy == _model_energy(lat, hardwall, ff)

    assert len(jumps) >= 40
    assert all(local == model for (local, model) in jumps), [
        j for j in jumps if j[0] != j[1]
    ][:5]
    assert any(local != 0 for (local, _) in jumps)
    # both branches were taken: some jumps kept, some undone
    assert 0 < accepted < len(jumps)


def test_single_chain_move_is_the_shared_helper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``Simulation.single_chain_move`` and the jump use one and the same routine."""
    state = U.build_state(tmp_path, 3, "SLR", False, {"MOVE_CHAIN_TRANSLATE": 1.0})
    sim = state.sim
    seen = []

    def fake(
        latticeObject: Any,
        hamiltonianObject: Any,
        move_event: Any,
        chainID: int,
        hardwall: bool,
    ) -> int:
        seen.append(
            (
                latticeObject is sim.LATTICE,
                hamiltonianObject is sim.Hamiltonian,
                move_event,
                chainID,
                hardwall,
            )
        )
        return 17

    monkeypatch.setattr(moves, "commit_single_chain_move", fake)
    assert sim.single_chain_move("event", 3) == 17
    assert seen == [(True, True, "event", 3, sim.hardwall)]


# ---------------------------------------------------------------------------
# build_all_envelope_pairs: per-chunk folding
# ---------------------------------------------------------------------------


def _reference_build_all_envelope_pairs(
    positions: list,
    LR_binary_array: IntArray,
    type_lattice: IntArray,
    dimensions: list[int],
    hardwall: bool = False,
    deduplicate: bool = True,
) -> tuple:
    """``lattice_utils.build_all_envelope_pairs`` as it was (commit 59163ce).

    Every bead's arrays are kept until one final concatenate per pair class.

    Parameters
    ----------
    positions : list
        Bead positions.
    LR_binary_array : numpy.ndarray
        1 where the bead is long-range.
    type_lattice : numpy.ndarray
        The type grid.
    dimensions : list of int
        The box.
    hardwall : bool
        Use the hardwall extractors.
    deduplicate : bool
        Remove duplicate pairs.

    Returns
    -------
    tuple of numpy.ndarray
        ``(SR_pairs, LR_pairs, SLR_pairs)``.
    """
    n_dim = len(dimensions)
    if n_dim == 2:
        extract = (
            inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_2D_hardwall
            if hardwall
            else inner_loops.extract_SR_and_LR_pairs_from_position_2D
        )
    else:
        extract = (
            inner_loops_hardwall.extract_SR_and_LR_pairs_from_position_3D_hardwall
            if hardwall
            else inner_loops.extract_SR_and_LR_pairs_from_position_3D
        )
    short_range_list = []
    long_range_list = []
    super_long_range_list = []
    for i in range(0, len(positions)):
        (SR_tmp, LR_tmp, SLR_tmp) = extract(
            np.array(positions[i], dtype=NP_INT_TYPE),
            LR_binary_array[i],
            type_lattice,
            *dimensions,
        )
        short_range_list.append(SR_tmp)
        if len(LR_tmp) > 0:
            long_range_list.append(LR_tmp)
        if len(SLR_tmp) > 0:
            super_long_range_list.append(SLR_tmp)

    short_range_pairs = np.concatenate(short_range_list)
    if len(long_range_list) > 0:
        long_range_pairs = np.concatenate(long_range_list)
    else:
        long_range_pairs = np.array([], dtype=NP_INT_TYPE)
    if len(super_long_range_list) > 0:
        super_long_range_pairs = np.concatenate(super_long_range_list)
    else:
        super_long_range_pairs = np.array([], dtype=NP_INT_TYPE)

    out = []
    for pairs in (short_range_pairs, long_range_pairs, super_long_range_pairs):
        reshaped = np.reshape(pairs, (len(pairs), 2 * n_dim))
        if deduplicate:
            reshaped = lattice_utils._unique_rows(reshaped)
        out.append(np.reshape(reshaped, (len(reshaped), 2, n_dim)))
    return tuple(out)


def _random_beads(
    rng: np.random.Generator, dims: list[int], occupancy: float, lr_mode: str
) -> tuple[list, IntArray, IntArray]:
    """Random beads on a lattice, as ``build_all_envelope_pairs`` takes them.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of the occupancy and the long-range flags.
    dims : list of int
        The box.
    occupancy : float
        Fraction of sites that hold a bead.
    lr_mode : str
        ``"SR"`` (no long-range bead), ``"LR"`` (every bead) or ``"MIX"``.

    Returns
    -------
    tuple
        ``(positions, LR_binary_array, type_grid)``; the positions are shuffled
        so neighbouring beads are not adjacent in the list.
    """
    type_grid = np.where(
        rng.random(dims) < occupancy, rng.integers(1, 4, size=dims), 0
    ).astype(NP_INT_TYPE)
    positions = [
        [int(x) for x in p] for p in rng.permutation(np.argwhere(type_grid > 0))
    ]
    if lr_mode == "SR":
        flags = np.zeros(len(positions), dtype=int)
    elif lr_mode == "LR":
        flags = np.ones(len(positions), dtype=int)
    else:
        flags = rng.integers(0, 2, size=len(positions))
    return (positions, flags, type_grid)


@pytest.mark.parametrize("dims", [[13, 9], [8, 7, 9]], ids=["2D", "3D"])
@pytest.mark.parametrize("lr_mode", ["SR", "LR", "MIX"])
@pytest.mark.parametrize("hardwall", [False, True], ids=["periodic", "hardwall"])
def test_envelope_pairs_equal_the_unchunked_routine(
    monkeypatch: pytest.MonkeyPatch, dims: list[int], lr_mode: str, hardwall: bool
) -> None:
    """Folding in chunks gives the same pair arrays, row for row.

    The energy is a sum over these rows, so they must come back equal AND in
    the same order, with and without de-duplication, whatever the chunk size.
    """
    rng = np.random.default_rng(len(dims) * 100 + int(hardwall) * 10 + len(lr_mode))
    (positions, flags, type_grid) = _random_beads(rng, dims, 0.35, lr_mode)
    assert len(positions) > 30
    for deduplicate in (True, False):
        want = _reference_build_all_envelope_pairs(
            positions,
            flags,
            type_grid,
            dims,
            hardwall=hardwall,
            deduplicate=deduplicate,
        )
        for chunk in (1, 7, 1024):
            monkeypatch.setattr(lattice_utils, "_ENVELOPE_PAIR_CHUNK", chunk)
            got = lattice_utils.build_all_envelope_pairs(
                positions,
                flags,
                type_grid,
                dims,
                hardwall=hardwall,
                deduplicate=deduplicate,
            )
            assert len(got) == 3
            for a, b in zip(got, want):
                assert a.dtype == b.dtype and a.shape == b.shape, (
                    chunk,
                    deduplicate,
                    a.shape,
                    b.shape,
                )
                assert np.array_equal(a, b), (chunk, deduplicate)
        if lr_mode != "SR":
            assert len(want[1]) > 0 and len(want[2]) > 0
        else:
            assert want[1].shape == (0, 2, len(dims)) and want[2].shape == (
                0,
                2,
                len(dims),
            )


def test_envelope_pairs_of_no_beads_are_three_empty_arrays() -> None:
    """The empty input keeps its ``(0, 2, n_dim)`` shape for all three classes."""
    for dims in ([9, 9], [8, 8, 8]):
        got = lattice_utils.build_all_envelope_pairs(
            [], np.array([], dtype=int), np.zeros(dims, dtype=NP_INT_TYPE), dims
        )
        assert [a.shape for a in got] == [(0, 2, len(dims))] * 3


def _traced_peak(function: Callable[[], Any]) -> tuple[Any, int]:
    """Run ``function`` and return its result and the peak traced memory it added.

    Parameters
    ----------
    function : callable
        Called with no arguments.

    Returns
    -------
    tuple
        ``(result, peak bytes above the level traced on entry)``.
    """
    tracemalloc.start()
    try:
        baseline = tracemalloc.get_traced_memory()[0]
        tracemalloc.reset_peak()
        result = function()
        peak = tracemalloc.get_traced_memory()[1] - baseline
    finally:
        tracemalloc.stop()
    return (result, peak)


def test_envelope_pairs_do_not_pin_a_padded_buffer_per_bead() -> None:
    """The transient is bounded per bead, and far below the unchunked routine's.

    3,000 long-range beads in a dilute 36^3 box: each bead has only a few
    long-range neighbours, but the extractor hands back views of padded
    (98, 2, 3) and (218, 2, 3) buffers, which the old routine kept alive for
    every bead at once (about 10 kB per bead).
    """
    rng = np.random.default_rng(11)
    dims = [36, 36, 36]
    (positions, flags, type_grid) = _random_beads(rng, dims, 3000 / 36**3, "LR")
    n_beads = len(positions)
    assert 2700 < n_beads < 3300

    (got, new_peak) = _traced_peak(
        lambda: lattice_utils.build_all_envelope_pairs(
            positions, flags, type_grid, dims
        )
    )
    (want, old_peak) = _traced_peak(
        lambda: _reference_build_all_envelope_pairs(positions, flags, type_grid, dims)
    )
    for a, b in zip(got, want):
        assert np.array_equal(a, b)

    assert old_peak / n_beads > 8000, (
        "the reference no longer shows the transient this test bounds"
    )
    assert new_peak / n_beads < 4500, new_peak / n_beads
    assert new_peak < 0.5 * old_peak


# ---------------------------------------------------------------------------
# get_empty_site: no whole-grid scan before the first draw
# ---------------------------------------------------------------------------


class _CountingGrid(np.ndarray):
    """An int32 grid that counts whole-array ``==`` comparisons made on it."""

    comparisons = 0

    def __eq__(self, other: Any) -> Any:  # type: ignore[override]
        type(self).comparisons += 1
        return np.asarray(self).__eq__(other)

    __hash__ = None  # type: ignore[assignment]


def _first_empty_site_by_hand(grid: IntArray, seed: int) -> tuple[list[int], tuple]:
    """The site rejection sampling must return, and the generator state it leaves.

    Parameters
    ----------
    grid : numpy.ndarray
        The occupancy grid.
    seed : int
        Seed for Python's generator.

    Returns
    -------
    tuple
        ``(site, random.getstate())`` after drawing one coordinate per axis,
        axis by axis, until the drawn site is empty.
    """
    random.seed(seed)
    while True:
        site = [random.randint(0, n - 1) for n in grid.shape]
        if grid[tuple(site)] == 0:
            return (site, random.getstate())


@pytest.mark.parametrize("dims", [[12, 9], [6, 7, 5]], ids=["2D", "3D"])
def test_get_empty_site_draws_are_unchanged_and_the_grid_is_not_scanned(
    dims: list[int],
) -> None:
    """Same site, same generator state, and no pass over the whole box.

    On a mostly empty lattice the site is found within a few draws, and the old
    up-front ``np.any(grid == 0)`` - one pass over the whole box per chain placed -
    is not made at all.
    """
    rng = np.random.default_rng(len(dims))
    plain = np.where(rng.random(dims) < 0.3, 1, 0).astype(NP_INT_TYPE)
    for seed in range(25):
        (want_site, want_state) = _first_empty_site_by_hand(plain, seed)
        grid = plain.copy().view(_CountingGrid)
        _CountingGrid.comparisons = 0
        random.seed(seed)
        site = lattice_utils.get_empty_site(grid)
        assert [int(x) for x in site] == want_site
        assert random.getstate() == want_state
        assert _CountingGrid.comparisons == 0


def test_get_empty_site_on_a_nearly_full_and_a_full_lattice() -> None:
    """The last empty site is still found by the same draws; a full box still raises."""
    dims = [6, 6, 6]
    plain = np.ones(dims, dtype=NP_INT_TYPE)
    plain[4, 1, 3] = 0
    found_after_100 = 0
    for seed in range(6):
        (want_site, want_state) = _first_empty_site_by_hand(plain, seed)
        random.seed(seed)
        draws = 0
        while True:
            draws += 1
            if [random.randint(0, n - 1) for n in dims] == [4, 1, 3]:
                break
        found_after_100 += draws > 100
        random.seed(seed)
        assert (
            [int(x) for x in lattice_utils.get_empty_site(plain)]
            == want_site
            == [4, 1, 3]
        )
        assert random.getstate() == want_state
    assert found_after_100 > 0, (
        "no seed exercised the every-100-draws full-lattice check"
    )

    with pytest.raises(LatticeUtilsException, match="fully occupied"):
        lattice_utils.get_empty_site(np.ones(dims, dtype=NP_INT_TYPE))


# ---------------------------------------------------------------------------
# Lattice.get_random_chain: cached candidate list
# ---------------------------------------------------------------------------


class _CountingDict(dict):
    """A dict that counts how many times it is iterated from the start."""

    passes = 0

    def __iter__(self) -> Any:
        type(self).passes += 1
        return super().__iter__()


def _bare_lattice(chain_ids: list[int]) -> lattice.Lattice:
    """A Lattice holding only a chains dictionary (no grids are needed here).

    Parameters
    ----------
    chain_ids : list of int
        The chainIDs, in insertion order.

    Returns
    -------
    pimms.lattice.Lattice
        An otherwise uninitialised lattice whose chains map each ID to a marker.
    """
    lat = lattice.Lattice.__new__(lattice.Lattice)
    lat.chains = _CountingDict((cid, "chain-%d" % cid) for cid in chain_ids)
    return lat


def test_get_random_chain_same_draws_without_a_pass_over_the_chains_per_call() -> None:
    """Nothing frozen: the draw is ``random.choice`` over the chainIDs, and the
    chains dictionary is walked once, not once per step."""
    ids = [3, 1, 7, 12, 5, 9]
    lat = _bare_lattice(ids)
    for frozen in (None, [], ()):
        random.seed(42)
        want = ["chain-%d" % random.choice(ids) for _ in range(200)]
        random.seed(42)
        _CountingDict.passes = 0
        got = [lat.get_random_chain(frozen_chains=frozen) for _ in range(200)]
        assert got == want
        assert _CountingDict.passes <= 1


def test_get_random_chain_follows_changes_to_the_chains() -> None:
    """Adding or removing a chain is seen by the next call; freezing still works."""
    lat = _bare_lattice([1, 2, 3])
    random.seed(1)
    assert {lat.get_random_chain() for _ in range(100)} == {
        "chain-1",
        "chain-2",
        "chain-3",
    }

    lat.chains[4] = "chain-4"
    assert {lat.get_random_chain() for _ in range(200)} == {
        "chain-1",
        "chain-2",
        "chain-3",
        "chain-4",
    }

    del lat.chains[2]
    assert {lat.get_random_chain() for _ in range(200)} == {
        "chain-1",
        "chain-3",
        "chain-4",
    }

    # same count, different membership: 4 leaves, 6 arrives
    del lat.chains[4]
    lat.chains[6] = "chain-6"
    assert {lat.get_random_chain() for _ in range(200)} == {
        "chain-1",
        "chain-3",
        "chain-6",
    }

    assert {lat.get_random_chain(frozen_chains=[3]) for _ in range(200)} == {
        "chain-1",
        "chain-6",
    }
    random.seed(9)
    want = ["chain-%d" % random.choice([1, 6]) for _ in range(50)]
    random.seed(9)
    assert [lat.get_random_chain(frozen_chains=[3]) for _ in range(50)] == want


# ---------------------------------------------------------------------------
# parallel slither / pull pass: frozen mask from the chain lengths
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("dim", [2, 3])
def test_parallel_pass_frozen_mask_marks_exactly_the_held_out_chains(
    tmp_path: Path, dim: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The mask handed to the parallel kernel is 1 on every bead of every chain
    that is not in the parallel set, and is built without the per-megamove
    ``np.isin`` over the whole bead table."""
    state = U.build_state(tmp_path, dim, "LR", False, {"MOVE_SLITHER": 1.0})
    lat = state.lattice
    (idx_to_bead, sorted_chains, chain_offset, chain_length, chain_homo) = (
        moves.parallel_chain_metadata(lat)
    )
    bead_chain = np.asarray(idx_to_bead)[:, 4].tolist()
    assert len(set(chain_length.tolist())) > 1, "chains of different lengths are needed"

    def no_isin(*args: Any, **kwargs: Any) -> None:
        raise AssertionError(
            "the parallel pass rebuilt its mask with _frozen_bead_mask"
        )

    monkeypatch.setattr(moves, "_frozen_bead_mask", no_isin)

    captured = []

    def fake_parallel_kernel(*args: Any) -> tuple[int, int, int]:
        captured.append(args[-1])
        return (0, 0, 0)

    rng = np.random.default_rng(dim)
    n_chains = len(sorted_chains)
    head_args = (
        lat.grid,
        lat.type_grid,
        idx_to_bead,
        chain_offset,
        chain_length,
        chain_homo,
    )
    for trial in range(12):
        size = int(rng.integers(1, n_chains + 1))
        parallel_chains = np.sort(rng.choice(n_chains, size=size, replace=False))
        moves._two_pass_whole_chain_megamove(
            parallel_chains,
            np.array([], dtype=parallel_chains.dtype),
            2,
            fake_parallel_kernel,
            None,
            head_args,
            (),
            0,
            1.0,
            False,
            int(chain_length.max()),
            1,
            idx_to_bead,
            sorted_chains,
        )
        mask = captured[-1]
        parallel_ids = {sorted_chains[ci] for ci in parallel_chains.tolist()}
        want = [0 if cid in parallel_ids else 1 for cid in bead_chain]
        assert mask.dtype == np.int32 and mask.flags["C_CONTIGUOUS"]
        assert mask.tolist() == want, trial
    assert len(captured) == 12


# ---------------------------------------------------------------------------
# Chain.positions holds Python ints after every Python-level move
# ---------------------------------------------------------------------------


def _coordinate_types(lat: lattice.Lattice) -> set:
    """The set of types found among every coordinate of every chain.

    Parameters
    ----------
    lat : pimms.lattice.Lattice
        The lattice.

    Returns
    -------
    set of type
        ``{int}`` when every coordinate is a Python int.
    """
    return {
        type(x)
        for chain in lat.chains.values()
        for p in chain.get_ordered_positions()
        for x in p
    }


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("hardwall", [False, True], ids=["periodic", "hardwall"])
def test_python_moves_leave_python_int_coordinates(
    tmp_path: Path, dim: int, hardwall: bool
) -> None:
    """Rotations and head pivots no longer leave numpy scalars in ``Chain.positions``.

    ``np.dot`` in ``run_rotation`` gave ``np.int64`` coordinates (chain rotate,
    chain pivot, cluster rotate) and the adjacent-site array gave ``np.int32``
    (head pivot). A numpy scalar wraps silently in later scalar arithmetic and
    is pickled into the restart file as a numpy object.
    """
    state = U.build_state(tmp_path, dim, "SLR", hardwall, {"MOVE_CHAIN_ROTATE": 1.0})
    sim = state.sim
    lat = sim.LATTICE
    mover = sim.MOVER
    assert _coordinate_types(lat) == {int}

    random.seed(8)
    np.random.seed(8)
    long_chains = [
        cid
        for cid, chain in lat.chains.items()
        if len(chain.get_ordered_positions()) >= 4
    ]
    single = {
        "chain_rotate": mover.chain_rotate,
        "chain_pivot": mover.chain_pivot,
        "head_pivot": mover.head_pivot,
    }
    done = {name: 0 for name in list(single) + ["cluster_rotate"]}
    for attempt in range(400):
        chain_id = long_chains[attempt % len(long_chains)]
        for name, move in single.items():
            (event, success) = move(lat.chains[chain_id], lat.grid, hardwall=hardwall)
            if success:
                sim.single_chain_move(event, chain_id)
                done[name] += 1
                assert _coordinate_types(lat) == {int}, name
        (event, success) = mover.cluster_rotate(
            lat.chains[chain_id],
            lat,
            cluster_move_threshold=None,
            cluster_size_threshold=len(lat.chains) - 1,
            hardwall=hardwall,
            frozen_chains=[],
        )
        if success:
            sim.rigid_cluster_move(event.moved_positions, event.original_positions)
            done["cluster_rotate"] += 1
            assert _coordinate_types(lat) == {int}, "cluster_rotate"
        if min(done.values()) >= 3:
            break
    assert min(done.values()) >= 3, done
    assert int(sim.Hamiltonian.evaluate_total_energy(lat)[0]) == _model_energy(
        lat, hardwall, _oracle_forcefield("SLR")
    )


def test_run_rotation_returns_python_ints() -> None:
    """The rotation helper hands back lists of Python ints with the rotated values."""
    quarter_turn = np.array([[0, -1], [1, 0]])
    rotated = lattice_utils.run_rotation(
        [[1, 0], [2, -3], np.array([0, 4])], quarter_turn
    )
    assert rotated == [[0, 1], [3, 2], [-4, 0]]
    assert {type(x) for p in rotated for x in p} == {int} and all(
        type(p) is list for p in rotated
    )
    for rotated_3d, want in (
        (lattice_utils.rotate_positions_3D([[1, 2, 3]], "z", 90), None),
        (lattice_utils.rotate_positions_2D([[1, 2]], 180), [[-1, -2]]),
    ):
        assert {type(x) for p in rotated_3d for x in p} == {int}
        if want is not None:
            assert rotated_3d == want
