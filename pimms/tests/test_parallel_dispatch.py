"""
The PARALLELIZE dispatch of the whole-chain megamoves must not depend on the state.

``system_slither`` and ``system_pull`` can each run on two kernels: the serial one,
which moves any chain anywhere, and the parallel checkerboard one, which only moves
a chain lying wholly inside one block interior and rejects any sub-move that would
leave that interior. Both kernels are individually pi-invariant - that is what the
kernel-level detailed-balance tests check - but a choice BETWEEN them that reads the
current configuration is not pi-invariant, and here the asymmetry is total. The
parallel kernel's interior is closed, so the rate out of "every chain is compact
enough to fit" is identically zero through it, while the serial kernel crosses that
boundary freely. Before 1.0.8 the dispatch made exactly that choice per megamove, and
the result was a silent one-signed pump into compact conformations.

The dispatch is now a partition of the chains by LENGTH, computed once from run
constants: short chains go to the parallel kernel, long chains to the serial kernel,
both in the same megamove. These tests pin the two things that has to mean.

  1. The partition is a function of chain lengths, the frozen set and the box only.
     Cheap, exact tests: the partition does not move over a run, and it is identical
     for a hand-built compact configuration and a hand-built extended one (where the
     old gate demonstrably differed).

  2. Going through ``MoveObject.system_slither`` / ``system_pull`` with
     ``parallelize=True`` samples the same equilibrium as going through them with
     ``parallelize=False``. The fixture is ATHERMAL - every pair, solvation and angle
     parameter is exactly zero, so the target measure is uniform on self-avoiding
     configurations and no relaxation-rate difference can masquerade as agreement or
     as a bias. Tolerances are built from the SERIAL reference arm alone, never from
     the arm under test.

The last test is a positive control: it re-runs the same comparison with the
partition replaced by the old state-dependent gate and asserts the comparison FAILS,
so a green suite means the comparison has power rather than being vacuous.

The fixture is 2D 24x24 short-range: the chain-level halo is W=3, the box splits into
2x2 blocks of 12 and the interior is 6. Three 8-mers are therefore always on the
serial side of the partition and their extent straddles the interior (which is what
used to flip the old gate); three 5-mers are always on the parallel side, so both
kernels really do run every megamove. Everything is seeded, so the numbers below are
reproducible rather than a fresh coin flip on every CI run.
"""

from __future__ import annotations

import contextlib
import os
from typing import Callable, Iterator

import numpy as np
import pytest

from pimms import mega_crank_fast, moves
from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U

# 2D 24x24 SR -> W=3, 2x2 blocks of 12, interior 6 (asserted in the fixture builder)
BOX: list[int] = [24, 24]
LONG_SEQUENCE: str = "AABBAABB"      # 8 beads: longer than the interior -> serial side
SHORT_SEQUENCE: str = "AABBA"        # 5 beads: fits the interior -> parallel side
N_PER_KIND: int = 3

# statistical arms
BURN_IN_CYCLES: int = 300
SAMPLE_CYCLES: int = 6000
CRANK_SUBSTEPS: int = 1000
WHOLE_CHAIN_SUBSTEPS: int = 10
K_SIGMA: float = 6.0


@pytest.fixture(autouse=True)
def _restore_cwd() -> Iterator[None]:
    """Put the working directory back after a test, since building a Simulation chdirs."""
    cwd = os.getcwd()
    yield
    os.chdir(cwd)


# ---------------------------------------------------------------------------
# athermal fixture
# ---------------------------------------------------------------------------

def _write_athermal_param_file(path: str) -> None:
    """Write a parameter file in which every energy term is exactly zero.

    An athermal system is the cleanest possible discriminator for this defect. The
    target measure becomes the uniform measure over self-avoiding configurations, so
    any difference between the two arms is a property of the move dispatch and cannot
    be blamed on the two arms relaxing an energy landscape at different rates. It also
    gives a free consistency check: the tracked energy must stay identically 0.

    Parameters
    ----------
    path : str
        Where to write the parameter file (the keyfile points at ``params.prm``).
    """
    lines = [
        "ANGLE_PENALTY\tA\t0\t0\t0",
        "ANGLE_PENALTY\tB\t0\t0\t0",
        "A\t0\t0",
        "B\t0\t0",
        "A  A\t0",
        "B  B\t0",
        "A  B\t0",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _build_athermal_state(tmpdir: str, moveset: dict[str, float], seed: int = 11) -> U.State:
    """Build the athermal 2D fixture: three 8-mers and three 5-mers in a 24x24 box.

    Parameters
    ----------
    tmpdir : str
        Directory to write ``KEYFILE.kf`` and ``params.prm`` into and to build in.

    moveset : dict
        ``{MOVE_KEYWORD: fraction}`` for the keyfile. The tests drive the megamoves
        directly, so this only has to be a legal move set.

    seed : int, optional
        The keyfile SEED. The Simulation constructor reseeds both ``random`` and
        ``numpy.random`` from it, so a freshly built state starts a reproducible
        stream. Default is 11.

    Returns
    -------
    pimms.tests.kernel_test_utils.State
        The built state, whose energy is asserted to be exactly zero.
    """
    os.makedirs(tmpdir, exist_ok=True)
    U.write_keyfile(os.path.join(tmpdir, "KEYFILE.kf"), 2, False, moveset,
                    box=BOX, chains=[(N_PER_KIND, LONG_SEQUENCE), (N_PER_KIND, SHORT_SEQUENCE)],
                    seed=seed, temperature=50,
                    extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 2})
    _write_athermal_param_file(os.path.join(tmpdir, "params.prm"))

    cwd = os.getcwd()
    os.chdir(tmpdir)
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(keyfile.keyword_lookup)
    finally:
        os.chdir(cwd)

    state = U.State(sim)
    assert state.energy_terms == (0, 0, 0, 0, 0), \
        "the fixture is meant to be athermal: %r" % (state.energy_terms,)
    return state


def _interior() -> int:
    """The movable block interior of the fixture box, from the kernel's own layout.

    Returns
    -------
    int
        The smallest per-axis interior (``block_size - 2W``) over the split axes.
    """
    info = mega_crank_fast.parallel_layout_info(BOX[0], BOX[1], 1, False)
    assert info["num_blocks"] > 1, "the fixture box must split, else there is no dispatch to test"
    return min(info["block_size"][d] - 2 * info["W"] for d in range(2) if info["blocks"][d] > 1)


# ---------------------------------------------------------------------------
# observables, computed here from scratch rather than borrowed from the code
# under test
# ---------------------------------------------------------------------------

def _periodic_extent(coords: np.ndarray, box: int) -> int:
    """Number of sites a set of periodic coordinates spans along one axis.

    The span is the box minus the largest circular gap between neighbouring
    coordinates (including the wrap-around gap), so a chain straddling a box face is
    measured by its true physical span rather than by its raw wrapped coordinates.

    Parameters
    ----------
    coords : numpy.ndarray
        The coordinates along one axis, in any order.

    box : int
        The box length along that axis.

    Returns
    -------
    int
        The occupied span in sites (1 for a single coordinate).
    """
    c = np.sort(np.asarray(coords, dtype=np.int64))
    gap = int(np.diff(c).max()) if len(c) > 1 else 0
    wrap = int(c[0] + box - c[-1])
    return int(box - max(gap, wrap) + 1)


def _max_chain_extent(lattice) -> int:
    """Largest periodic extent of any chain along any axis, in sites.

    Parameters
    ----------
    lattice : Lattice
        The lattice whose chains are measured.

    Returns
    -------
    int
        ``max`` over chains and axes of the per-axis periodic extent.
    """
    best = 0
    for chainID in sorted(lattice.chains.keys()):
        positions = np.asarray(lattice.chains[chainID].get_ordered_positions(), dtype=np.int64)
        for d in range(positions.shape[1]):
            best = max(best, _periodic_extent(positions[:, d], int(lattice.dimensions[d])))
    return best


def _mean_rg2(lattice) -> float:
    """Mean squared radius of gyration per chain, on the bond-walked chains.

    Parameters
    ----------
    lattice : Lattice
        The lattice whose chains are measured. ``Chain.get_analysis_positions()``
        returns each chain made whole in a single periodic image, which is what any
        intra-chain observable has to be computed on.

    Returns
    -------
    float
        The average over chains of ``<(r - r_com)^2>``.
    """
    values = []
    for chainID in sorted(lattice.chains.keys()):
        positions = np.asarray(lattice.chains[chainID].get_analysis_positions(), dtype=float)
        values.append(float(((positions - positions.mean(axis=0)) ** 2).sum(axis=1).mean()))
    return float(np.mean(values))


# ---------------------------------------------------------------------------
# the old (1.0.8) state-dependent gate, for the positive control
# ---------------------------------------------------------------------------

def _old_state_dependent_gate(idx_to_bead, chain_offset, chain_length, dimensions, has_LR,
                              chain_homo=None, cap_mode='all', frozen_chains=()):
    """Reproduce the 1.0.8 dispatch: all chains parallel, or all chains serial.

    The old rule measured every chain's CURRENT periodic extent and sent the whole
    megamove to the parallel kernel only if every one of them fitted a block interior.
    This is written out here from the current positions rather than called back into
    ``moves``, both to keep the control independent of the code under test and because
    ``parallel_chain_fit_report`` now reports the new partition (so calling it from a
    replacement for the partition would recurse).

    Parameters
    ----------
    idx_to_bead : numpy.ndarray
        The ``int64`` bead bookkeeping matrix; columns 5 onwards are the coordinates
        this rule reads, which is precisely what makes it state-dependent.

    chain_offset : numpy.ndarray
        Row of each chain's first bead in ``idx_to_bead``.

    chain_length : numpy.ndarray
        Number of beads in each chain.

    dimensions : list of int
        Lattice box dimensions.

    has_LR : bool
        Whether any bead is long-range. Accepted for signature compatibility with
        :func:`pimms.moves.parallel_chain_partition`; the fixture is short-range and
        the interior is taken from the module-level layout instead.

    chain_homo : numpy.ndarray or None, optional
        Accepted for signature compatibility; the fixture has no chain near the
        512-bead buffer cap, so it is not consulted. Default is None.

    cap_mode : str, optional
        Accepted for signature compatibility; see ``chain_homo``. Default is ``'all'``.

    frozen_chains : sequence of int, optional
        Accepted for signature compatibility; the fixture freezes nothing. Default
        is ``()``.

    Returns
    -------
    numpy.ndarray
        All-True when every chain currently fits, all-False otherwise.
    """
    idx = np.asarray(idx_to_bead)
    offsets = np.asarray(chain_offset, dtype=np.int64)
    lengths = np.asarray(chain_length, dtype=np.int64)
    n_dim = len(dimensions)
    interior = _interior()
    fits = True
    for ci in range(len(offsets)):
        beads = idx[offsets[ci]:offsets[ci] + lengths[ci], 5:5 + n_dim]
        for d in range(n_dim):
            if _periodic_extent(beads[:, d], int(dimensions[d])) > interior:
                fits = False
    return np.full(len(offsets), fits, dtype=bool)


# ---------------------------------------------------------------------------
# cheap, exact tests of the partition itself
# ---------------------------------------------------------------------------

def test_fixture_splits_the_chains_between_both_kernels(tmp_path) -> None:
    """The fixture must exercise BOTH kernels, else the statistical arms prove nothing."""
    state = _build_athermal_state(str(tmp_path), {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5})
    idx, _chains, offsets, lengths, homo = moves.parallel_chain_metadata(state.lattice)
    interior = _interior()
    assert interior == 6
    assert sorted(lengths.tolist()) == [5, 5, 5, 8, 8, 8]

    partition = moves.parallel_chain_partition(idx, offsets, lengths, BOX, False,
                                               chain_homo=homo, cap_mode="hetero")
    assert partition.tolist() == [length <= interior for length in lengths]
    assert int(partition.sum()) == N_PER_KIND          # the 5-mers -> parallel kernel
    assert int((~partition).sum()) == N_PER_KIND       # the 8-mers -> serial kernel


def test_partition_ignores_the_configuration_but_the_old_gate_did_not(tmp_path) -> None:
    """Same chains, compact vs extended: the partition is identical, the old gate is not.

    This is the invariant the fix establishes, asserted directly and exactly rather
    than statistically: the dispatch decision is a function of chain lengths and the
    box, so driving the system between the two extremes of the extent distribution
    cannot move it.
    """
    state = _build_athermal_state(str(tmp_path), {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5})
    idx, _chains, offsets, lengths, homo = moves.parallel_chain_metadata(state.lattice)
    interior = _interior()

    def _with_positions(placer: Callable[[int, int], list[int]]) -> np.ndarray:
        """Copy the bookkeeping matrix with chain ci's bead j moved to placer(ci, j)."""
        moved = np.asarray(idx).copy()
        for ci in range(len(offsets)):
            for j in range(int(lengths[ci])):
                moved[offsets[ci] + j, 5:7] = placer(ci, j)
        return moved

    # extended: every chain laid out as a straight rod, so an 8-mer spans 8 > 6
    extended = _with_positions(lambda ci, j: [2 * ci, j])
    # compact: every chain folded into a 3x3 patch, so nothing spans more than 3
    compact = _with_positions(lambda ci, j: [4 * ci + j // 3, j % 3])

    kwargs = dict(chain_homo=homo, cap_mode="hetero")
    part_extended = moves.parallel_chain_partition(extended, offsets, lengths, BOX, False, **kwargs)
    part_compact = moves.parallel_chain_partition(compact, offsets, lengths, BOX, False, **kwargs)
    assert part_extended.tolist() == part_compact.tolist()
    assert part_compact.tolist() == [length <= interior for length in lengths]

    # ... and the control really is a different rule: the old gate flips between the two
    assert not moves._parallel_can_move_all_chains(extended, offsets, lengths, BOX, False, **kwargs)
    assert moves._parallel_can_move_all_chains(compact, offsets, lengths, BOX, False, **kwargs)
    assert _old_state_dependent_gate(extended, offsets, lengths, BOX, False, **kwargs).tolist() \
        == [False] * len(lengths)
    assert _old_state_dependent_gate(compact, offsets, lengths, BOX, False, **kwargs).tolist() \
        == [True] * len(lengths)


def test_partition_does_not_move_over_a_run(tmp_path) -> None:
    """Recompute the partition after every megamove of a real run; it must never change."""
    state = _build_athermal_state(str(tmp_path),
                                  {"MOVE_CRANKSHAFT": 0.4, "MOVE_SLITHER": 0.3, "MOVE_PULL": 0.3})
    interior = _interior()

    def _partitions() -> tuple[list[bool], list[bool]]:
        idx, _chains, offsets, lengths, homo = moves.parallel_chain_metadata(state.lattice)
        slither = moves.parallel_chain_partition(idx, offsets, lengths, BOX, False,
                                                 chain_homo=homo, cap_mode="hetero")
        pull = moves.parallel_chain_partition(idx, offsets, lengths, BOX, False,
                                              cap_mode="all")
        return (slither.tolist(), pull.tolist())

    first = _partitions()
    energy = state.energy
    extents_seen = set()
    for _ in range(200):
        state.lattice, energy, _, _ = state.sim.MOVER.system_shake(
            state.lattice, energy, state.acc, state.ham, CRANK_SUBSTEPS, "UNSET", parallelize=False)
        state.lattice, energy, _, _ = state.sim.MOVER.system_slither(
            state.lattice, energy, state.acc, state.ham, WHOLE_CHAIN_SUBSTEPS,
            parallelize=True, num_threads=2)
        state.lattice, energy, _, _ = state.sim.MOVER.system_pull(
            state.lattice, energy, state.acc, state.ham, WHOLE_CHAIN_SUBSTEPS,
            parallelize=True, num_threads=2)
        extents_seen.add(_max_chain_extent(state.lattice))
        assert _partitions() == first

    assert energy == 0, "the athermal fixture must track exactly zero energy throughout"
    # the run has to actually visit both sides of the interior, otherwise "the partition
    # never changed" would be true for the old state-dependent gate as well
    assert min(extents_seen) <= interior < max(extents_seen), sorted(extents_seen)


@pytest.mark.parametrize("ff,box,hardwall", [("SR", [32, 32, 32], False),
                                             ("LR", [40, 40, 40], False),
                                             ("SLR", [48, 48, 48], True)])
def test_split_megamove_keeps_the_tracked_energy_exact_in_3d(tmp_path, ff, box, hardwall) -> None:
    """The two passes chain through the same buffers, so the energy must stay exact.

    This covers the 3D kernels and the hardwall / long-range branches, where the split
    dispatch has to hand the running energy from the parallel pass into the serial one
    and back. The system deliberately mixes 20-mers (always serial: longer than any of
    these interiors) with 4-mers and 8-mers (parallel), so both passes do real work.
    """
    state = U.build_state(str(tmp_path), 3, ff, hardwall,
                          {"MOVE_SLITHER": 0.5, "MOVE_PULL": 0.5}, box=box,
                          chains=[(4, "A" * 20), (6, "AABB"), (3, "AABBAABB")],
                          temperature=45)
    idx, _chains, offsets, lengths, homo = moves.parallel_chain_metadata(state.lattice)
    partition = moves.parallel_chain_partition(idx, offsets, lengths, box, state.has_LR(),
                                               chain_homo=homo, cap_mode="hetero")
    assert partition[:4].tolist() == [False] * 4          # the 20-mers
    assert partition[4:].all()                            # the 4-mers and 8-mers

    energy = state.energy
    for _ in range(15):
        state.lattice, energy, _, _ = state.sim.MOVER.system_slither(
            state.lattice, energy, state.acc, state.ham, 8, hardwall=hardwall,
            parallelize=True, num_threads=2)
        state.lattice, energy, _, _ = state.sim.MOVER.system_pull(
            state.lattice, energy, state.acc, state.ham, 8, hardwall=hardwall,
            parallelize=True, num_threads=2)
    assert energy == int(state.ham.evaluate_total_energy(state.lattice)[0])


def test_both_kernels_are_called_every_megamove_with_a_fixed_chain_split(tmp_path,
                                                                        monkeypatch) -> None:
    """Every parallelized whole-chain megamove runs a parallel pass AND a serial pass.

    The two passes run in a random order, so the megamove is reversible.

    The kernels are replaced by recorders so the selectors and the frozen mask handed
    to each pass can be read off exactly; the system is still driven between megamoves
    by a real crankshaft megamove, so the recorded split is being re-tested against
    genuinely different configurations.
    """
    state = _build_athermal_state(str(tmp_path), {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5})
    idx, _chains, offsets, lengths, homo = moves.parallel_chain_metadata(state.lattice)
    partition = moves.parallel_chain_partition(idx, offsets, lengths, BOX, False,
                                               chain_homo=homo, cap_mode="hetero")
    expected_parallel = sorted(np.nonzero(partition)[0].tolist())
    expected_serial = sorted(np.nonzero(~partition)[0].tolist())

    calls: list[tuple[str, list[int], int, int]] = []

    def _recorder(kind: str) -> Callable[..., tuple[int, ...]]:
        def kernel(*args):
            selector = np.asarray(args[6])
            n_masked = int(np.asarray(args[-1]).sum()) if kind == "parallel" else -1
            calls.append((kind, sorted(np.unique(selector).tolist()), len(selector), n_masked))
            # the parallel kernels also report the attempts they made
            return (args[11], 0, len(selector)) if kind == "parallel" else (args[11], 0)
        return kernel

    monkeypatch.setattr(moves.mega_crank_fast, "mega_slither_parallel_2D", _recorder("parallel"))
    monkeypatch.setattr(moves.mega_crank_fast, "mega_slither_2D", _recorder("serial"))

    energy = state.energy
    for _ in range(30):
        state.lattice, energy, _, _ = state.sim.MOVER.system_shake(
            state.lattice, energy, state.acc, state.ham, CRANK_SUBSTEPS, "UNSET", parallelize=False)
        _lat, _e, proposed, _accepted = state.sim.MOVER.system_slither(
            state.lattice, energy, state.acc, state.ham, WHOLE_CHAIN_SUBSTEPS,
            parallelize=True, num_threads=2)
        # the two passes' proposal counts must sum to what a single serial pass over
        # every selectable chain would have proposed
        assert proposed == 2 * N_PER_KIND * WHOLE_CHAIN_SUBSTEPS

    assert len(calls) == 60
    # one parallel pass and one serial pass every megamove, in an order drawn by a
    # fair coin (a fixed order is not reversible, which a system-wide TSMMC
    # excursion needs); over 30 megamoves both orders must turn up
    orders = [tuple(kind for kind, _sel, _n, _m in calls[k:k + 2]) for k in range(0, 60, 2)]
    assert set(orders) <= {("parallel", "serial"), ("serial", "parallel")}
    assert len(set(orders)) == 2, "the pass order never changed - the coin is not being drawn"
    for kind, selector, n_substeps, n_masked in calls:
        assert n_substeps == N_PER_KIND * WHOLE_CHAIN_SUBSTEPS
        if kind == "parallel":
            assert selector == expected_parallel
            # the long chains are held out of the parallel pass through the frozen mask,
            # not through the selector, because the parallel kernel reads the selector
            # only for its length and picks its chains per block
            assert n_masked == N_PER_KIND * len(LONG_SEQUENCE)
        else:
            assert selector == expected_serial


# ---------------------------------------------------------------------------
# equilibrium arms
# ---------------------------------------------------------------------------

def _sample_arm(state: U.State, move: str, parallelize: bool,
                cycles: int = SAMPLE_CYCLES) -> dict[str, np.ndarray]:
    """Run one arm and return its conformational observable traces.

    A cycle is a SERIAL crankshaft megamove plus one whole-chain megamove, so the two
    arms differ only in how the whole-chain megamove is dispatched.

    Parameters
    ----------
    state : pimms.tests.kernel_test_utils.State
        The athermal fixture to run, mutated in place.

    move : str
        ``'slither'`` or ``'pull'`` - which whole-chain megamove to drive.

    parallelize : bool
        Passed straight through to the megamove.

    cycles : int, optional
        Number of sampled cycles after the burn-in. Default is ``SAMPLE_CYCLES``.

    Returns
    -------
    dict
        ``{"max_extent": ndarray, "fraction_extended": ndarray, "rg2": ndarray}``,
        one entry per sampled cycle.
    """
    megamove = (state.sim.MOVER.system_slither if move == "slither"
                else state.sim.MOVER.system_pull)
    interior = _interior()
    energy = state.energy
    max_extent, fraction_extended, rg2 = [], [], []
    for cycle in range(BURN_IN_CYCLES + cycles):
        state.lattice, energy, _, _ = state.sim.MOVER.system_shake(
            state.lattice, energy, state.acc, state.ham, CRANK_SUBSTEPS, "UNSET",
            parallelize=False)
        state.lattice, energy, _, _ = megamove(
            state.lattice, energy, state.acc, state.ham, WHOLE_CHAIN_SUBSTEPS,
            parallelize=parallelize, num_threads=2)
        if cycle >= BURN_IN_CYCLES:
            extent = _max_chain_extent(state.lattice)
            max_extent.append(float(extent))
            fraction_extended.append(1.0 if extent > interior else 0.0)
            rg2.append(_mean_rg2(state.lattice))
    assert energy == 0, "athermal run drifted off zero energy: %r" % (energy,)
    return {"max_extent": np.array(max_extent),
            "fraction_extended": np.array(fraction_extended),
            "rg2": np.array(rg2)}


def _compare_arms(tmp_path, move: str, partition_override=None) -> None:
    """Assert a parallelized arm samples the serial arm's equilibrium.

    Both arms are built fresh (which reseeds ``random`` and ``numpy.random`` from the
    keyfile SEED), burned in under their own dispatch, and compared on three
    conformational observables with ``assert_same_equilibrium(..., ref_only=True)``,
    whose tolerance is built from the SERIAL reference trace alone. Letting the arm
    under test contribute to its own tolerance is exactly how a biased arm hides.

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytest temporary directory; each arm gets its own subdirectory.

    move : str
        ``'slither'`` or ``'pull'``.

    partition_override : callable or None, optional
        Replacement for ``moves.parallel_chain_partition`` in the test arm only, used
        by the positive control to reinstate the old state-dependent gate. Default is
        None (the production partition).
    """
    keyword = "MOVE_SLITHER" if move == "slither" else "MOVE_PULL"
    reference = _sample_arm(_build_athermal_state(str(tmp_path / "serial"),
                                                  {"MOVE_CRANKSHAFT": 0.5, keyword: 0.5}),
                            move, parallelize=False)

    test_state = _build_athermal_state(str(tmp_path / "parallel"),
                                       {"MOVE_CRANKSHAFT": 0.5, keyword: 0.5})
    original = moves.parallel_chain_partition
    if partition_override is not None:
        moves.parallel_chain_partition = partition_override
    try:
        test = _sample_arm(test_state, move, parallelize=True)
    finally:
        moves.parallel_chain_partition = original

    # a serial-vs-serial comparison would pass trivially; the reference arm must
    # genuinely sample both sides of the interior
    assert 0.05 < reference["fraction_extended"].mean() < 0.95

    for observable in ("fraction_extended", "max_extent", "rg2"):
        U.assert_same_equilibrium(reference[observable], test[observable],
                                  "%s parallelize=True, %s" % (move, observable),
                                  k_sigma=K_SIGMA, ref_only=True)


def test_slither_parallelize_samples_the_serial_equilibrium(tmp_path) -> None:
    """PARALLELIZE must not change the conformational equilibrium of MOVE_SLITHER."""
    _compare_arms(tmp_path, "slither")


def test_pull_parallelize_samples_the_serial_equilibrium(tmp_path) -> None:
    """PARALLELIZE must not change the conformational equilibrium of MOVE_PULL."""
    _compare_arms(tmp_path, "pull")


def test_the_old_state_dependent_gate_is_caught(tmp_path) -> None:
    """Positive control: the comparison above must FAIL for the 1.0.8 dispatch.

    Without this the two tests above could be green because the comparison has no
    power. With the old gate reinstated the parallelized arm over-samples compact
    conformations, and the fraction of samples with a chain extended beyond the block
    interior drops from about 0.28 to about 0.11 - a ~5x tolerance violation.
    """
    with pytest.raises(AssertionError, match="fraction_extended"):
        _compare_arms(tmp_path, "slither", partition_override=_old_state_dependent_gate)


# ---------------------------------------------------------------------------
# every start position is reachable when the box does not divide into blocks
# ---------------------------------------------------------------------------

def _straight_chain(box: tuple[int, ...], x0: int, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A straight short-range homopolymer along x from ``x0``, in an otherwise empty box.

    Parameters
    ----------
    box : tuple of int
        Box dimensions (2 or 3 entries).
    x0 : int
        x coordinate of the first bead (wrapped into the box).
    n : int
        Number of beads.

    Returns
    -------
    tuple of numpy.ndarray
        ``(grid, type_grid, idx_to_bead)`` for the one-chain system.
    """
    n_dim = len(box)
    grid = np.zeros(box, dtype=np.int32)
    type_grid = np.zeros(box, dtype=np.int32)
    idx = np.zeros((n, 5 + n_dim), dtype=np.int64)
    for k in range(n):
        flag = 1 if k == 0 else (3 if k == n - 1 else (5 if k == 1 else (6 if k == n - 2 else 2)))
        site = ((x0 + k) % box[0],) + (5,) * (n_dim - 1)
        idx[k, :5] = (flag, 0, 1, 0, 1)
        idx[k, 5:] = site
        grid[site] = 1
        type_grid[site] = 1
    return grid, type_grid, idx


@pytest.mark.parametrize("dim", (2, 3))
@pytest.mark.parametrize("move", ("slither", "pull"))
def test_an_interior_length_chain_is_movable_from_every_start_on_a_remainder_axis(dim: int, move: str) -> None:
    """A chain the length partition hands to the parallel kernel can move from anywhere.

    The partition sends a chain to the parallel kernel when its length is at most
    the block interior, which is only a sufficient condition for fitting if the
    random block shift can put an interior start at every coordinate. A 25-site axis
    splits into two 12-site blocks plus one remainder site; with the shift drawn from
    one block length (as it was) the interior starts covered only 24 of the 25
    coordinates, and a straight chain as long as the interior sitting across the
    remainder was never inside an interior on any sweep - neither kernel ever moved
    it. The kernel's attempt count says whether a sweep found the chain movable.
    """
    box = (25,) * dim
    info = mega_crank_fast.parallel_layout_info(box[0], box[1], box[2] if dim == 3 else 1, False)
    assert info["num_blocks"] > 1
    assert box[0] % info["blocks"][0] != 0, "the axis must leave a remainder, or this tests nothing"
    interior = info["block_size"][0] - 2 * info["W"]

    if dim == 3:
        kernel = mega_crank_fast.mega_slither_parallel if move == "slither" else mega_crank_fast.mega_pull_parallel
        angles = np.zeros((2, 3, 3, 3, 3, 3, 3), dtype=np.int32)
        n_sweeps = 1500
    else:
        kernel = mega_crank_fast.mega_slither_parallel_2D if move == "slither" else mega_crank_fast.mega_pull_parallel_2D
        angles = np.zeros((2, 3, 3, 3, 3), dtype=np.int32)
        n_sweeps = 600
    tables = (np.zeros((2, 2), dtype=np.int32),) * 3 + (angles,)
    offsets = np.array([0], dtype=np.int32)
    lengths = np.array([interior], dtype=np.int32)
    homo = np.array([1], dtype=np.int32)
    selector = np.zeros(1, dtype=np.int32)
    frozen = np.zeros(interior, dtype=np.int32)

    never_movable = []
    for x0 in range(box[0]):
        grid, type_grid, idx = _straight_chain(box, x0, interior)
        assert moves.parallel_chain_partition(idx, offsets, lengths, list(box), False)[0]
        for seed in range(1, n_sweeps + 1):
            _energy, _accepted, attempted = kernel(
                grid.copy(), type_grid.copy(), idx.copy(), offsets, lengths, homo, selector,
                *tables, 0, 0.0, seed, 0, interior, 1, frozen)
            if attempted > 0:
                break
        else:
            never_movable.append(x0)
    assert not never_movable, (
        f"a straight {interior}-mer starting at x = {never_movable} was never inside a block "
        f"interior in {n_sweeps} sweeps")


def test_logged_proposals_are_the_attempts_the_kernels_report(tmp_path, monkeypatch) -> None:
    """A parallel sweep that finds nothing movable is logged as no attempts.

    The parallel kernels return early, attempting nothing, when their random block
    shift leaves no bead (crankshaft) or chain (slither, pull) inside a block
    interior. The wrappers used to add the whole requested budget to the proposal
    count regardless, so MOVE_FREQS.dat counted attempts never made and
    ACCEPTANCE.dat understated the acceptance of every parallelized move. Here the
    parallel kernels are stood in for by ones that report such an empty sweep.
    """
    state = _build_athermal_state(str(tmp_path), {"MOVE_CRANKSHAFT": 0.5, "MOVE_SLITHER": 0.5})

    def empty_sweep(*args):
        # the entry energy is the eighth positional argument of the crankshaft
        # kernels and the twelfth of the whole-chain ones
        return (args[7] if len(args) == 14 else args[11], 0, 0)

    monkeypatch.setattr(moves.mega_crank_fast, "mega_crank_parallel_2D", empty_sweep)
    monkeypatch.setattr(moves.mega_crank_fast, "mega_slither_parallel_2D", empty_sweep)

    _lat, _e, proposed, accepted = state.sim.MOVER.system_shake(
        state.lattice, state.energy, state.acc, state.ham, CRANK_SUBSTEPS, "UNSET",
        parallelize=True, num_threads=2)
    assert (proposed, accepted) == (0, 0)

    _lat, _e, proposed, _accepted = state.sim.MOVER.system_slither(
        state.lattice, state.energy, state.acc, state.ham, WHOLE_CHAIN_SUBSTEPS,
        parallelize=True, num_threads=2)
    # only the serial pass (the long chains) made any attempts
    assert proposed == N_PER_KIND * WHOLE_CHAIN_SUBSTEPS
