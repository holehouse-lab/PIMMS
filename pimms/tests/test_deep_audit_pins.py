## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Pins for the places where the deep audit's mutation testing found that a wrong
program passed every test.

Each test here was written against a specific surviving mutant and was checked
to FAIL with that mutant applied and to pass on the real code. In order:

* the HARDWALL bound in the VMMC link-energy scan (a neighbour sitting on the
  coordinate-0 layer of an axis could be ignored without any test noticing),
  pinned twice - directly against a brute-force definition of the link energy,
  and through the equilibrium a hard-walled VMMC-only run samples;
* the absolute Boltzmann weight of the serial crankshaft kernels, against a
  fully enumerated ensemble (a common-mode beta x 1.15 in every Metropolis test
  used to be invisible to every crankshaft-specific fast test);
* the entry and exit work terms of the system-wide TSMMC excursion;
* input validation in the restart reader, the restart resize and the Chain
  constructor;
* the cadence of PERFORMANCE.dat and the per-axis arithmetic of
  ``parallel_chain_fit_report``.

The statistical tests compare against EXACT expectation values and judge the
difference in units of a standard error taken from independent blocks of the
run (never from the per-sample spread of a correlated trace), with a 5 sigma
threshold. Each one states the size of the effect it is there to resolve.
"""

import contextlib
import io
import itertools
import math
import os
import pickle
import random

import numpy as np
import pytest

from pimms import CONFIG, moves
from pimms.chain import Chain
from pimms.chainTSMMC import TSMMC
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import ChainInitializationException, RestartException
from pimms.restart import RestartObject
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U


# ---------------------------------------------------------------------------
# the model, written out from its definition
# ---------------------------------------------------------------------------

# bead-solvent energy per solvent site in a bead's Chebyshev-1 shell
_SOLVATION = {"A": -2, "B": -1}

# pair energies as (short range, long range, super-long range): the energy of a
# bead pair at Chebyshev distance 1, 2 and 3 respectively
_PAIRS = {("A", "A"): (-8, -4, 2), ("B", "B"): (-6, -3, 3), ("A", "B"): (-3, -2, 1)}

# how many of the three columns each pair gets in the parameter file. A bead
# type is long-range only if it appears on a line with a long-range column, so
# "MIXED" makes A long-range and leaves B short-range only.
_COLUMNS = {
    "SR": {("A", "A"): 1, ("B", "B"): 1, ("A", "B"): 1},
    "LR": {("A", "A"): 2, ("B", "B"): 2, ("A", "B"): 2},
    "SLR": {("A", "A"): 3, ("B", "B"): 3, ("A", "B"): 3},
    "MIXED": {("A", "A"): 3, ("B", "B"): 1, ("A", "B"): 1},
}

# every statistical assertion in this file uses the same threshold
_Z_MAX = 5.0

# Two A monomers in a 7x7 hard-walled box, bound strongly enough (the contact
# is worth 4 / 1.2 = 3.3 kT net of the two solvent contacts it removes) that a
# bound pair travels mostly as a two-chain VMMC cluster.
_VMMC_BOX = 7
_VMMC_TEMPERATURE = 1.2
_VMMC_BLOCKS = 50
_VMMC_PER_BLOCK = 4000
_VMMC_BURN_IN = 2000
_VMMC_EVENTS = (
    "P(contact)",
    "P(contact, a bead on a coordinate-0 layer)",
    "P(contact, a bead on a coordinate L-1 layer)",
    "P(a bead on a coordinate-0 layer)",
)

# An A monomer and an AB dimer in a periodic box of side 7, with all three
# interaction ranges switched on. Side 7 is the smallest box in which the
# Chebyshev-3 shell does not meet its own periodic image.
_CRANK_BOX = 7
_CRANK_SEQUENCE = "AAB"  # bead 0 is the monomer, beads 1-2 the dimer
_CRANK_TEMPERATURE = 3.0
_CRANK_BLOCKS = 50
# megamoves per block. The reference kernel costs several times more per call,
# so it gets a shorter run (and a correspondingly smaller, still ample, margin).
_CRANK_PER_BLOCK = {True: 1500, False: 400}
_CRANK_SUBSTEPS = 30  # bead moves per megamove (10 per bead)
_CRANK_CONTROL_SCALE = 1.15

_TSMMC_PROTOCOL = dict(target=50.0, jump=80.0, n_points=3, multiplier=4)


def _write_parameter_file(path, kind):
    """Write the parameter file for one of the force fields in ``_COLUMNS``.

    Parameters
    ----------
    path : str
        Where to write the file.

    kind : str
        A key of ``_COLUMNS``: which of the short, long and super-long range
        columns each pair is given.

    Returns
    -------
    None
    """
    lines = ["ANGLE_PENALTY\tA\t30\t10\t0", "ANGLE_PENALTY\tB\t50\t20\t0"]
    for bead, value in _SOLVATION.items():
        lines.append(f"{bead}\t0\t{value}")
    for pair, values in _PAIRS.items():
        columns = "\t".join(str(v) for v in values[: _COLUMNS[kind][pair]])
        lines.append(f"{pair[0]}  {pair[1]}\t{columns}")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


@contextlib.contextmanager
def _working_directory(path):
    """Run a block with ``path`` as the working directory, then change back.

    Parameters
    ----------
    path : str or pathlib.Path
        Directory to change into.

    Yields
    ------
    None
    """
    previous = os.getcwd()
    os.chdir(str(path))
    try:
        yield
    finally:
        os.chdir(previous)


def _build(
    directory,
    kind,
    hardwall,
    moveset,
    box,
    chains,
    temperature,
    *,
    seed=5,
    n_steps=10,
    equilibration=1,
    extra=None,
):
    """Build a Simulation in ``directory`` and return the kernel-test bundle for it.

    The same thing as ``kernel_test_utils.build_state``, except that the
    parameter file is the one written by this module, so the constants the
    exact enumerations below use are the ones the engine was given.

    Parameters
    ----------
    directory : str or pathlib.Path
        Directory the key file, the parameter file and every output go into.

    kind : str
        Force field, a key of ``_COLUMNS``.

    hardwall : bool
        Hard walls (True) or periodic boundaries (False).

    moveset : dict
        ``{MOVE_KEYWORD: fraction}``.

    box : list of int
        Box dimensions; their number sets the dimensionality.

    chains : list of (int, str)
        ``(count, sequence)`` for each chain type.

    temperature : float
        Simulation temperature.

    seed : int, optional
        The key file's SEED.

    n_steps, equilibration : int, optional
        Run length and equilibration length, for the tests that run the main loop.

    extra : dict, optional
        Extra key file keywords.

    Returns
    -------
    kernel_test_utils.State
        The built simulation and the objects the kernel drivers need from it.
    """
    directory = str(directory)
    _write_parameter_file(os.path.join(directory, "params.prm"), kind)
    U.write_keyfile(
        os.path.join(directory, "KEYFILE.kf"),
        len(box),
        hardwall,
        moveset,
        box=box,
        chains=chains,
        seed=seed,
        n_steps=n_steps,
        equilibration=equilibration,
        temperature=temperature,
        extra=extra,
    )
    with _working_directory(directory), contextlib.redirect_stdout(io.StringIO()):
        simulation = Simulation(KeyFileParser("KEYFILE.kf").keyword_lookup)
    return U.State(simulation)


def _pair_value(first, second, column, kind):
    """Energy of one pair of bead types at one range, zero if that range is off.

    Parameters
    ----------
    first, second : str
        The two bead types.

    column : int
        0 for the short range (Chebyshev distance 1), 1 for long range
        (distance 2), 2 for super-long range (distance 3).

    kind : str
        Force field, a key of ``_COLUMNS``.

    Returns
    -------
    int
        The pair energy, or 0 if the force field has no such column for the pair.
    """
    pair = (first, second) if (first, second) in _PAIRS else (second, first)
    return _PAIRS[pair][column] if column < _COLUMNS[kind][pair] else 0


def _chebyshev(first, second, box, hardwall):
    """Chebyshev distance between bead positions, under the boundary in force.

    Parameters
    ----------
    first, second : numpy.ndarray
        Integer positions, broadcastable against each other, last axis the
        spatial one.

    box : sequence of int
        Box dimensions.

    hardwall : bool
        With hard walls the raw separation is used; with periodic boundaries
        the minimum image.

    Returns
    -------
    numpy.ndarray
        The largest per-axis separation, with the spatial axis removed.
    """
    separation = np.abs(
        np.asarray(first, dtype=np.int64) - np.asarray(second, dtype=np.int64)
    )
    if not hardwall:
        separation = np.minimum(
            separation, np.asarray(box, dtype=np.int64) - separation
        )
    return separation.max(axis=-1)


def _oracle_energy(positions, sequence, box, hardwall, kind):
    """Total energy of bead configurations, from the definition of the model.

    Every pair of beads contributes its short, long or super-long range energy
    when it sits at Chebyshev distance 1, 2 or 3, and every bead contributes its
    solvation energy once for each site of its Chebyshev-1 shell that holds no
    bead (a site beyond a hard wall counts as solvent). Nothing here calls the
    PIMMS energy code; it is valid for chains of one or two beads (no angle
    term) and for uniform force fields (SR, LR, SLR).

    Parameters
    ----------
    positions : numpy.ndarray
        Integer array of shape (n_states, n_beads, n_dimensions).

    sequence : str
        Bead type of each bead, in the order of the second axis.

    box : sequence of int
        Box dimensions.

    hardwall : bool
        Boundary convention.

    kind : str
        Force field, one of "SR", "LR", "SLR".

    Returns
    -------
    numpy.ndarray
        Integer energy of each state.
    """
    positions = np.asarray(positions, dtype=np.int64)
    n_states, n_beads, n_dim = positions.shape
    energy = np.zeros(n_states, dtype=np.int64)
    contacts = np.zeros((n_states, n_beads), dtype=np.int64)
    for i, j in itertools.combinations(range(n_beads), 2):
        distance = _chebyshev(positions[:, i], positions[:, j], box, hardwall)
        for column in range(3):
            energy += (distance == column + 1) * _pair_value(
                sequence[i], sequence[j], column, kind
            )
        contacts[:, i] += distance == 1
        contacts[:, j] += distance == 1
    solvation = np.array([_SOLVATION[bead] for bead in sequence])
    energy += ((3**n_dim - 1 - contacts) * solvation).sum(axis=1)
    return energy


def _block_z(samples, exact, n_blocks):
    """z-score of a sample mean against an exact value, from independent blocks.

    The trace is cut into ``n_blocks`` consecutive blocks, each far longer than
    the correlation time of the run, and the standard error of the mean is the
    standard deviation of the block means over the square root of their number.
    The per-sample spread of the trace never enters, so a correlated trace
    cannot make the error look smaller than it is.

    Parameters
    ----------
    samples : numpy.ndarray
        The trace, of a length divisible by ``n_blocks``.

    exact : float
        The exact expectation value.

    n_blocks : int
        Number of blocks.

    Returns
    -------
    (float, float, float)
        The z-score, the sample mean and the block standard error.
    """
    block_means = np.asarray(samples, dtype=float).reshape(n_blocks, -1).mean(axis=1)
    mean = float(block_means.mean())
    sem = float(block_means.std(ddof=1) / math.sqrt(n_blocks))
    return (mean - exact) / sem, mean, sem


# ---------------------------------------------------------------------------
# VMMC link energies under HARDWALL: direct, against a brute-force definition
# ---------------------------------------------------------------------------


def _brute_force_link_energies(state, chain_id, offset, hardwall):
    """Link energies of one (virtually shifted) chain, by looping over all beads.

    The definition ``MoveObject._vmmc_neighbour_energies`` implements: every
    bead of the chain, moved by ``offset``, against every bead of every OTHER
    chain at its real position - the short-range table at Chebyshev distance 1
    and, if the moved bead is a long-range bead, the long and super-long range
    tables at distances 2 and 3. Separations are raw under HARDWALL and minimum
    image otherwise. Positions, bead types and long-range flags are read from
    the Chain objects, not from the lattice grids the routine under test scans.

    Parameters
    ----------
    state : kernel_test_utils.State
        The built system.

    chain_id : int
        The chain being virtually moved.

    offset : sequence of int
        The virtual translation.

    hardwall : bool
        Boundary convention.

    Returns
    -------
    (dict, list)
        ``{chainID: energy}`` holding the non-zero link energies, and a list of
        the positions of the partner beads that contributed a non-zero term.
    """
    lattice = state.lattice
    box = list(lattice.dimensions)
    tables = (
        state.ham.residue_interaction_table,
        state.ham.LR_residue_interaction_table,
        state.ham.SLR_residue_interaction_table,
    )
    moved = lattice.chains[chain_id]
    moved_positions = np.asarray(
        moved.get_ordered_positions(), dtype=np.int64
    ) + np.asarray(offset)
    if not hardwall:
        moved_positions = moved_positions % np.asarray(box)
    moved_types = moved.get_intcode_sequence()
    moved_long_range = moved.get_LR_binary_array()

    energies = {}
    contributing = []
    for other_id, other in lattice.chains.items():
        if other_id == chain_id:
            continue
        other_types = other.get_intcode_sequence()
        total = 0
        for partner, partner_type in zip(other.get_ordered_positions(), other_types):
            for bead in range(len(moved_positions)):
                distance = int(
                    _chebyshev(
                        moved_positions[bead], np.asarray(partner), box, hardwall
                    )
                )
                if distance == 1:
                    term = tables[0][moved_types[bead]][partner_type]
                elif distance in (2, 3) and moved_long_range[bead]:
                    term = tables[distance - 1][moved_types[bead]][partner_type]
                else:
                    term = 0
                if term != 0:
                    total += term
                    contributing.append(tuple(int(c) for c in partner))
        if total != 0:
            energies[int(other_id)] = float(total)
    return energies, contributing


@pytest.mark.parametrize("hardwall", [True, False], ids=["HW", "PBC"])
@pytest.mark.parametrize("kind", ["SR", "LR", "SLR", "MIXED"])
@pytest.mark.parametrize("box", [[9, 8], [7, 8, 7]], ids=["2D", "3D"])
def test_vmmc_link_energies_match_the_brute_force_definition(
    tmp_path, box, kind, hardwall
):
    """The VMMC link-energy scan must see every partner the definition sees.

    Mutation testing found that under HARDWALL the scan could skip every
    neighbour on the coordinate-0 layer of an axis (``coord < 0`` weakened to
    ``coord <= 0``) and no test failed. Here every chain of a dense box is
    virtually shifted by a set of translations and the routine's answer is
    compared with a loop over all other beads. The box is dense enough that
    interacting partners sit on the coordinate-0 and the coordinate L-1 layer of
    every axis, and that is asserted rather than assumed.
    """
    n_dim = len(box)
    if n_dim == 2:
        chains = [(4, "AB"), (6, "A"), (4, "B"), (2, "ABA")]
    else:
        chains = [(10, "AB"), (20, "A"), (20, "B"), (6, "ABA")]
    state = _build(
        tmp_path, kind, hardwall, {"MOVE_VMMC": 1.0}, box, chains, 10, seed=3
    )
    lattice, mover = state.lattice, state.sim.MOVER
    shells = {
        1: moves._vmmc_offset_shell(n_dim, 1),
        3: moves._vmmc_offset_shell(n_dim, 3),
    }

    rng = np.random.RandomState(20261001)
    on_layer = np.zeros(
        (n_dim, 2), dtype=int
    )  # partners seen on layer 0 / L-1 per axis
    compared = 0
    for chain_id in sorted(lattice.chains):
        chain = lattice.chains[chain_id]
        positions = chain.get_ordered_positions()
        offsets = [[0] * n_dim] + [
            list(map(int, rng.randint(-2, 3, n_dim))) for _ in range(8)
        ]
        for offset in offsets:
            shifted = np.asarray(positions) + np.asarray(offset)
            if hardwall and ((shifted < 0).any() or (shifted >= np.asarray(box)).any()):
                # a hard-walled proposal that leaves the box is rejected whatever
                # its link energies are, so they are not pinned
                continue
            expected, contributing = _brute_force_link_energies(
                state, chain_id, offset, hardwall
            )
            observed = mover._vmmc_neighbour_energies(
                lattice,
                state.ham,
                hardwall,
                shells,
                int(chain_id),
                positions,
                chain.get_intcode_sequence(),
                chain.get_LR_binary_array(),
                offset,
                list(lattice.dimensions),
            )
            observed = {int(k): float(v) for k, v in observed.items() if v != 0}
            assert observed == expected, (
                f"chain {chain_id} shifted by {offset}: link energies {observed} "
                f"but the definition gives {expected}"
            )
            compared += 1
            for partner in contributing:
                for axis in range(n_dim):
                    on_layer[axis, 0] += partner[axis] == 0
                    on_layer[axis, 1] += partner[axis] == box[axis] - 1

    assert compared > 50
    assert (on_layer > 0).all(), (
        f"the fixture never put an interacting partner on some boundary layer: {on_layer}"
    )


# ---------------------------------------------------------------------------
# VMMC under HARDWALL: the sampled equilibrium against exact enumeration
# ---------------------------------------------------------------------------


def _vmmc_hardwall_exact():
    """Exact expectation values for the two-monomer hard-walled VMMC system.

    Enumerates every ordered pair of distinct sites, weights it with the oracle
    energy, and returns the probabilities of the events the test samples.

    Returns
    -------
    numpy.ndarray
        The exact probability of each event of :func:`_vmmc_events`, in the
        order of ``_VMMC_EVENTS``.
    """
    sites = np.array(list(itertools.product(range(_VMMC_BOX), repeat=2)))
    first, second = np.meshgrid(
        np.arange(len(sites)), np.arange(len(sites)), indexing="ij"
    )
    keep = first != second
    states = np.stack([sites[first[keep]], sites[second[keep]]], axis=1)
    energy = _oracle_energy(states, "AA", [_VMMC_BOX] * 2, True, "SR")
    weight = np.exp(
        -(CONFIG.INVTEMP_FACTOR / _VMMC_TEMPERATURE) * (energy - energy.min())
    )
    weight /= weight.sum()
    events = _vmmc_events(states)
    return (weight[:, None] * events).sum(axis=0)


def _vmmc_events(states):
    """The indicator observables of the hard-walled VMMC test.

    Parameters
    ----------
    states : numpy.ndarray
        Integer array of shape (n_states, 2, 2): the two monomer positions.

    Returns
    -------
    numpy.ndarray
        Float array of shape (n_states, 4), in the order of ``_VMMC_EVENTS``:
        in contact; in contact with a bead on the coordinate-0 layer of either
        axis; in contact with a bead on the coordinate L-1 layer of either axis;
        a bead on the coordinate-0 layer of either axis, in contact or not.
    """
    contact = _chebyshev(states[:, 0], states[:, 1], [_VMMC_BOX] * 2, True) == 1
    low = (states == 0).any(axis=(1, 2))
    high = (states == _VMMC_BOX - 1).any(axis=(1, 2))
    return np.stack([contact, contact & low, contact & high, low], axis=1).astype(float)


def _vmmc_hardwall_z(directory, seed):
    """Sample the hard-walled two-monomer system with VMMC only; return z-scores.

    Parameters
    ----------
    directory : str or pathlib.Path
        Directory to build the system in.

    seed : int
        Seed of both random number generators the Python moves draw from.

    Returns
    -------
    (numpy.ndarray, numpy.ndarray, numpy.ndarray)
        z-scores, sample means and exact values of the events of
        :func:`_vmmc_events`.
    """
    L = _VMMC_BOX
    state = _build(
        directory,
        "SR",
        True,
        {"MOVE_VMMC": 1.0},
        [L, L],
        [(2, "A")],
        _VMMC_TEMPERATURE,
        seed=7,
        extra={"VMMC_MAX_DISPLACEMENT": 2, "VMMC_MAX_CLUSTER": 2},
    )
    lattice, ham, acc, mover = state.lattice, state.ham, state.acc, state.sim.MOVER
    assert acc.invtemp == pytest.approx(CONFIG.INVTEMP_FACTOR / _VMMC_TEMPERATURE)

    def positions():
        """Return the two monomer positions as a (1, 2, 2) integer array."""
        return np.array(
            [
                [
                    lattice.chains[1].get_ordered_positions()[0],
                    lattice.chains[2].get_ordered_positions()[0],
                ]
            ],
            dtype=np.int64,
        )

    random.seed(seed)
    np.random.seed(seed)
    energy = ham.evaluate_total_energy(lattice)[0]
    assert energy == _oracle_energy(positions(), "AA", [L, L], True, "SR")[0]

    n_moves = _VMMC_BLOCKS * _VMMC_PER_BLOCK
    trace = np.empty((n_moves, 2, 2), dtype=np.int64)
    for move in range(-_VMMC_BURN_IN, n_moves):
        seed_chain = lattice.chains[random.randint(1, 2)]
        _, energy, _accepted, _size = mover.vmmc_move(
            seed_chain, lattice, energy, acc, ham, 2, 2, hardwall=True, frozen_chains=[]
        )
        if move >= 0:
            trace[move] = positions()[0]

    # a hard-walled move never leaves the box, and the running energy is exact
    assert trace.min() >= 0 and trace.max() <= L - 1
    assert energy == _oracle_energy(positions(), "AA", [L, L], True, "SR")[0]
    assert ham.evaluate_total_energy(lattice)[0] == pytest.approx(energy)

    exact = _vmmc_hardwall_exact()
    events = _vmmc_events(trace)
    z = np.empty(len(_VMMC_EVENTS))
    mean = np.empty(len(_VMMC_EVENTS))
    for k in range(len(_VMMC_EVENTS)):
        z[k], mean[k], _sem = _block_z(events[:, k], exact[k], _VMMC_BLOCKS)
    return z, mean, exact


@pytest.mark.slow
def test_hardwall_vmmc_wall_layer_occupancy_matches_exact_boltzmann(tmp_path):
    """Hard-walled VMMC must populate the two wall layers of an axis equally.

    If the link-energy scan ignores a partner on the coordinate-0 layer, a
    chain next to such a partner never recruits it, while the same pair is
    recruited normally everywhere else. Cluster moves then carry bound pairs
    onto the coordinate-0 walls about twice as often as off them, and the
    equilibrium piles up there: measured with that bug, P(contact and a bead on
    a coordinate-0 layer) reads 0.33 against an exact 0.24 and the layer
    occupancy 0.43 against 0.33, while P(contact) alone moves by only 3 to 4
    standard errors - which is why the older P(contact)-only test at T = 4
    could not see it.

    The four probabilities are compared with their exactly enumerated values
    at 5 block standard errors (50 blocks of 4000 moves; the correlation time
    of the contact indicator is below 100 moves). Over 20 seeds the bug put the
    two coordinate-0 probabilities 10 to 17 and 15 to 22 standard errors out,
    and the correct scan stayed within 2.5 on all four.

    Marked slow: the direct test above is the fast pin, and a VMMC-only run
    long enough to resolve this from the equilibrium is 200000 Python moves.
    """
    z, mean, exact = _vmmc_hardwall_z(tmp_path, seed=13)
    for k, name in enumerate(_VMMC_EVENTS):
        assert abs(z[k]) < _Z_MAX, (
            f"hardwall VMMC {name}: sampled {mean[k]:.4f}, exact {exact[k]:.4f}, "
            f"z = {z[k]:+.1f}"
        )


# ---------------------------------------------------------------------------
# serial crankshaft kernels: absolute Boltzmann weight, by exact enumeration
# ---------------------------------------------------------------------------


def _crank_exact(n_dim, beta):
    """Exact <E> and P(contact) of the monomer + dimer system at one beta.

    The ensemble is enumerated modulo translation: the first dimer bead sits at
    the origin, the second on any of its Chebyshev-1 neighbours, and the monomer
    on any remaining site. Every such state stands for the same number of
    translated copies, so the weights need no multiplicity.

    Parameters
    ----------
    n_dim : int
        2 or 3.

    beta : float
        Inverse temperature.

    Returns
    -------
    (float, float, int)
        The mean energy, the probability that the monomer touches the dimer
        (Chebyshev distance 1 to either dimer bead), and the number of states
        enumerated.
    """
    L = _CRANK_BOX
    origin = (0,) * n_dim
    states = []
    for bond in itertools.product((-1, 0, 1), repeat=n_dim):
        if not any(bond):
            continue
        second = tuple(c % L for c in bond)
        for monomer in itertools.product(range(L), repeat=n_dim):
            if monomer != origin and monomer != second:
                states.append([monomer, origin, second])
    states = np.array(states, dtype=np.int64)
    energy = _oracle_energy(states, _CRANK_SEQUENCE, [L] * n_dim, False, "SLR")
    weight = np.exp(-beta * (energy - energy.min()))
    weight /= weight.sum()
    contact = _crank_contact(states)
    return float((weight * energy).sum()), float((weight * contact).sum()), len(states)


def _crank_contact(states):
    """Whether the monomer touches the dimer, for an array of states.

    Parameters
    ----------
    states : numpy.ndarray
        Integer array of shape (n_states, 3, n_dimensions): monomer, dimer
        bead, dimer bead.

    Returns
    -------
    numpy.ndarray
        Boolean array, True where the monomer is at Chebyshev distance 1 of
        either dimer bead.
    """
    box = [_CRANK_BOX] * states.shape[2]
    return _chebyshev(states[:, 0:1], states[:, 1:], box, False).min(axis=1) == 1


def _crank_trace(state, fast, seed, scale=1.0):
    """Energy and contact traces of a crankshaft-only run of the kernel itself.

    Parameters
    ----------
    state : kernel_test_utils.State
        The monomer + dimer system.

    fast : bool
        Drive the optimised serial kernel (True) or the reference kernel (False).

    seed : int
        Seeds the bead selection (numpy) and, offset per megamove, the kernel.

    scale : float, optional
        Multiplier on the inverse temperature the kernel is handed. 1.0 is the
        real run; anything else is the deliberate Metropolis error of the power
        control.

    Returns
    -------
    (numpy.ndarray, numpy.ndarray)
        Energy after each sampled megamove, and whether the monomer touched the
        dimer after it.
    """
    n_dim = state.dim
    np.random.seed(seed)
    grid, type_grid, idx = state.fresh()
    energy = state.energy
    n_samples = _CRANK_BLOCKS * _CRANK_PER_BLOCK[fast]
    energies = np.empty(n_samples)
    states = np.empty((n_samples, 3, n_dim), dtype=np.int64)
    kernel_seed = 7919 * (seed + 1)
    with U.scaled_invtemp(state, scale):
        for sample in range(-200, n_samples):  # 200 megamoves of burn-in
            kernel_seed += 1
            energy = U.crank_megastep(
                state,
                grid,
                type_grid,
                idx,
                energy,
                kernel_seed,
                substeps=_CRANK_SUBSTEPS,
                fast=fast,
            )
            if sample >= 0:
                energies[sample] = energy
                states[sample] = np.asarray(idx)[:, 5 : 5 + n_dim]
    # the energy the kernel tracked incrementally is the energy of where it ended up
    assert (
        energy
        == _oracle_energy(
            states[-1:], _CRANK_SEQUENCE, [_CRANK_BOX] * n_dim, False, "SLR"
        )[0]
    )
    return energies, _crank_contact(states).astype(float)


@pytest.mark.parametrize("fast", [True, False], ids=["fast", "reference"])
@pytest.mark.parametrize("n_dim", [2, 3], ids=["2D", "3D"])
def test_serial_crankshaft_kernels_sample_the_exact_boltzmann_distribution(
    tmp_path, n_dim, fast
):
    """Crankshaft megamoves alone must reproduce an exactly enumerated ensemble.

    This is the absolute pin on the Metropolis weight of the serial crankshaft
    kernels: no reference move, no second kernel - the 376 (2D) or 8866 (3D)
    states of a monomer and a dimer are enumerated, weighted with an energy
    written out from the model definition, and <E> and P(contact) of a
    crankshaft-only run are compared with the exact values.

    Power. With beta scaled by 1.15 inside the kernel the exact <E> moves from
    -44.337 to -44.586 in 2D and from -134.041 to -134.389 in 3D. Against the
    block standard errors of these runs (0.009 and 0.010 for the fast kernel's
    75000 megamoves, 0.017 and 0.020 for the reference kernel's 20000) that is
    28 and 34 standard errors for the fast kernel and 15 and 18 for the
    reference kernel. Measured over 20 seeds per case, the smallest |z| that
    error produced was 27 and 34 (fast), 13 and 16 (reference), and the largest
    |z| the correct kernels produced, on either observable, was 2.9. A kernel
    compiled with beta x 1.15 in its Metropolis tests fails all four cases, at
    z = -35, -17, -36 and -19.

    The same error is injected here as a positive control and must be thrown
    out at twice the threshold, so the fixture cannot silently lose its power.
    False failures: the block means give a t statistic with 49 degrees of
    freedom, for which |t| > 5 has probability 8e-6 per assertion; and the
    seeds are fixed, so a given build either passes or fails every time.
    """
    L = _CRANK_BOX
    temperature = _CRANK_TEMPERATURE
    beta = CONFIG.INVTEMP_FACTOR / temperature
    state = _build(
        tmp_path,
        "SLR",
        False,
        {"MOVE_CRANKSHAFT": 1.0},
        [L] * n_dim,
        [(1, "A"), (1, "AB")],
        temperature,
    )
    assert state.acc.invtemp == pytest.approx(beta, rel=1e-12)
    start = np.asarray(state.idx0)
    assert list(start[:, 4]) == [1, 2, 2]  # monomer first, then the dimer
    assert (
        state.energy
        == _oracle_energy(
            start[None, :, 5 : 5 + n_dim], _CRANK_SEQUENCE, [L] * n_dim, False, "SLR"
        )[0]
    )

    exact_energy, exact_contact, n_states = _crank_exact(n_dim, beta)
    assert n_states == (3**n_dim - 1) * (L**n_dim - 2)

    energies, contact = _crank_trace(state, fast, seed=1)
    z_energy, mean_energy, sem_energy = _block_z(energies, exact_energy, _CRANK_BLOCKS)
    z_contact, mean_contact, _ = _block_z(contact, exact_contact, _CRANK_BLOCKS)
    assert abs(z_energy) < _Z_MAX, (
        f"crankshaft <E> = {mean_energy:.4f} +/- {sem_energy:.4f}, exact {exact_energy:.4f}, "
        f"z = {z_energy:+.1f}"
    )
    assert abs(z_contact) < _Z_MAX, (
        f"crankshaft P(contact) = {mean_contact:.4f}, exact {exact_contact:.4f}, "
        f"z = {z_contact:+.1f}"
    )

    # positive control: the same run with the kernel handed beta x 1.15 must be
    # thrown out by the same criterion, and by a wide margin
    shifted_energy, _ = _crank_exact(n_dim, beta * _CRANK_CONTROL_SCALE)[:2]
    assert abs(shifted_energy - exact_energy) > 2 * _Z_MAX * sem_energy
    control, _ = _crank_trace(state, fast, seed=2, scale=_CRANK_CONTROL_SCALE)
    z_control = _block_z(control, exact_energy, _CRANK_BLOCKS)[0]
    assert abs(z_control) > 2 * _Z_MAX, (
        f"POSITIVE CONTROL FAILED: a kernel run at beta x {_CRANK_CONTROL_SCALE} gave "
        f"z = {z_control:+.1f} against the exact <E>; this fixture no longer resolves a "
        f"15 % Metropolis error"
    )


# ---------------------------------------------------------------------------
# system TSMMC: the accumulated work, against the tempered-transitions definition
# ---------------------------------------------------------------------------


class _RecordingAcceptance:
    """Stand-in for the AcceptanceCalculator a system TSMMC excursion drives.

    Records every temperature it is switched to, which is all the excursion
    bookkeeping asks of it.
    """

    def __init__(self, temperature):
        """
        Parameters
        ----------
        temperature : float
            The temperature the stand-in starts at (the target temperature).
        """
        self.temperature = float(temperature)
        self.history = []

    def update_temperature(self, temperature):
        """Record a temperature switch.

        Parameters
        ----------
        temperature : float
            The new temperature.

        Returns
        -------
        None
        """
        self.temperature = float(temperature)
        self.history.append(float(temperature))

    def get_total_aux_chain_moves(self):
        """Number of auxiliary-chain moves made so far (none in this script).

        Returns
        -------
        int
            Always 0.
        """
        return 0


def _scripted_system_excursion(energies, *, target, jump, n_points, multiplier):
    """Drive one system TSMMC excursion through a scripted energy sequence.

    The calls are the ones the main loop makes: ``start_system_TSMMC`` with the
    energy before the excursion, then ``check_in_system_TSMMC`` before every
    sub-move with the energy at that moment, then ``accept_system_TSMMC`` with
    the final energy.

    Parameters
    ----------
    energies : sequence of float
        ``energies[k]`` is the system energy after ``k`` sub-moves; its length
        must be one more than the number of sub-moves in the excursion.

    target, jump : float
        Target and jump temperatures.

    n_points : int
        Number of temperatures on the ramp (TSMMC_NUMBER_OF_POINTS).

    multiplier : int
        Sub-moves per temperature (TSMMC_STEP_MULTIPLIER).

    Returns
    -------
    (TSMMC, float, list, bool)
        The coordinator, the log-work it handed to the acceptance test, the
        temperature in force during each sub-move, and the accept decision.
    """
    coordinator = TSMMC(target, jump, "LINEAR", multiplier, n_points, False)
    acceptance = _RecordingAcceptance(target)
    handed_over = []
    decide = coordinator.accept_tempered_transition

    def spy(log_work):
        """Record the log-work handed to the acceptance test, then run it.

        Parameters
        ----------
        log_work : float
            The accumulated work of the excursion.

        Returns
        -------
        bool
            The real acceptance decision.
        """
        handed_over.append(float(log_work))
        return decide(log_work)

    coordinator.accept_tempered_transition = spy
    coordinator.start_system_TSMMC((None, None, None), energies[0], acceptance)
    during = []
    sub_move = 0
    while not coordinator.system_move_complete():
        acceptance = coordinator.check_in_system_TSMMC(acceptance, energies[sub_move])
        during.append(acceptance.temperature)
        sub_move += 1
    assert sub_move == len(energies) - 1
    accepted = coordinator.accept_system_TSMMC(energies[sub_move])
    assert len(handed_over) == 1
    return coordinator, handed_over[0], during, accepted


def _tempered_transition_work(energies, *, target, jump, n_points, multiplier):
    """Log-work of a system TSMMC excursion, from the tempered-transitions definition.

    The ladder is the target temperature, a linear ramp of ``n_points``
    temperatures ending on the jump temperature, a hold of ``CONFIG.TOP_TEMP``
    rungs at the jump temperature, the ramp mirrored, and the target again.
    ``multiplier`` sub-moves are made on every rung except the two target ends,
    and the log-work is the sum over every change of rung of
    ``(beta_before - beta_after) * E``, with E the energy at the instant of the
    change.

    Parameters
    ----------
    energies : sequence of float
        ``energies[k]`` is the energy after ``k`` sub-moves.

    target, jump : float
        Target and jump temperatures.

    n_points : int
        Number of temperatures on the ramp.

    multiplier : int
        Sub-moves per rung.

    Returns
    -------
    (float, list)
        The log-work and the ladder of temperatures between the two target ends.
    """
    ramp = [target + (jump - target) * k / n_points for k in range(1, n_points + 1)]
    rungs = ramp + [jump] * CONFIG.TOP_TEMP + ramp[::-1]
    betas = [CONFIG.INVTEMP_FACTOR / t for t in [target] + rungs + [target]]
    work = 0.0
    for change in range(len(betas) - 1):
        # change number `change` happens once `change * multiplier` sub-moves are done
        work += (betas[change] - betas[change + 1]) * energies[change * multiplier]
    return work, rungs


def _tsmmc_energies(seed, low, high):
    """A scripted energy for every point of the protocol in ``_TSMMC_PROTOCOL``.

    Parameters
    ----------
    seed : int
        Seed of the generator the energies are drawn from.

    low, high : float
        Range of the energies.

    Returns
    -------
    numpy.ndarray
        One energy per sub-move plus the starting energy.
    """
    n_rungs = 2 * _TSMMC_PROTOCOL["n_points"] + CONFIG.TOP_TEMP
    return np.random.RandomState(seed).uniform(
        low, high, n_rungs * _TSMMC_PROTOCOL["multiplier"] + 1
    )


def test_system_tsmmc_work_is_the_tempered_transitions_sum():
    """The log-work of a system TSMMC excursion must be sum (beta_k - beta_k+1) E_k.

    Mutation testing found that flipping the sign of the entry term (target ->
    first rung, at the starting energy) or dropping the exit term (last rung ->
    target, at the final energy) was caught only by a slow statistical test and
    by the regression fixtures. Here an energy sequence is scripted through the
    coordinator's own call sequence and the log-work it hands to the acceptance
    test is compared with the sum written out by hand.
    """
    energies = _tsmmc_energies(11, -900.0, -300.0)
    coordinator, log_work, during, _ = _scripted_system_excursion(
        energies, **_TSMMC_PROTOCOL
    )
    expected, rungs = _tempered_transition_work(energies, **_TSMMC_PROTOCOL)

    # the protocol is the one the definition assumes: every rung, in order, for
    # exactly `multiplier` sub-moves, and never the target temperature
    assert during == pytest.approx(
        [t for t in rungs for _ in range(_TSMMC_PROTOCOL["multiplier"])]
    )
    assert log_work == pytest.approx(expected, rel=1e-9, abs=1e-9)

    # each end term on its own is far larger than the tolerance above, so the
    # comparison does pin both of them
    beta_target = CONFIG.INVTEMP_FACTOR / _TSMMC_PROTOCOL["target"]
    beta_first = CONFIG.INVTEMP_FACTOR / rungs[0]
    assert abs((beta_target - beta_first) * energies[0]) > 0.1
    assert abs((beta_first - beta_target) * energies[-1]) > 0.1

    # a constant energy does no net work, whatever the ladder
    flat = np.full(len(energies), -512.0)
    assert _scripted_system_excursion(flat, **_TSMMC_PROTOCOL)[1] == pytest.approx(
        0.0, abs=1e-9
    )


@pytest.mark.parametrize(
    "start, end, accepted", [(-9000.0, -100.0, False), (-100.0, -9000.0, True)]
)
def test_system_tsmmc_accepts_downhill_and_rejects_uphill_excursions(
    start, end, accepted
):
    """The accept decision itself must follow the sign of the end terms.

    The energy is held at ``start`` until the last sub-move and then set to
    ``end``, so only the entry and the exit term of the work survive:
    ``(beta_target - beta_first) * (start - end)``. An excursion that ends
    8900 units below where it began has a log-work of +29.7 and is always
    accepted; the reverse has -29.7 and is never accepted (probability e^-29.7).
    A flipped entry term or a dropped exit term changes the sign or the size of
    that number and reverses one of the two decisions.
    """
    energies = np.full(len(_tsmmc_energies(0, 0.0, 1.0)), start)
    energies[-1] = end
    random.seed(3)
    _, log_work, _, decision = _scripted_system_excursion(energies, **_TSMMC_PROTOCOL)
    first_rung = (
        _TSMMC_PROTOCOL["target"]
        + (_TSMMC_PROTOCOL["jump"] - _TSMMC_PROTOCOL["target"]) / 3
    )
    by_hand = (
        CONFIG.INVTEMP_FACTOR / _TSMMC_PROTOCOL["target"]
        - CONFIG.INVTEMP_FACTOR / first_rung
    ) * (start - end)
    assert log_work == pytest.approx(by_hand, rel=1e-9)
    assert abs(by_hand) == pytest.approx(8900 * (1 / 50 - 1 / 60))
    assert bool(decision) is accepted


# ---------------------------------------------------------------------------
# input validation
# ---------------------------------------------------------------------------


def _restart_payload(box, position):
    """A minimal restart-file dictionary holding one single-bead chain.

    Parameters
    ----------
    box : list of int
        Box dimensions.

    position : list
        The bead's position.

    Returns
    -------
    dict
        The dictionary a restart file pickles.
    """
    return {
        "DIMENSIONS": list(box),
        "ENERGY": 0.0,
        "HARDWALL": False,
        "CHAINS": {1: [[list(position)], "A", 0]},
    }


@pytest.mark.parametrize(
    "box, position",
    [
        ([4, 5], [4, 1]),  # x == L
        ([4, 5], [1, 5]),  # y == L
        ([4, 5, 6], [1, 2, 6]),  # z == L
        ([4, 5], [5, 1]),  # beyond L
        ([4, 5], [-1, 1]),  # below 0
    ],
)
def test_restart_file_with_a_bead_outside_the_box_is_refused(tmp_path, box, position):
    """A bead at coordinate L (one past the last site) must be refused.

    Valid coordinates run from 0 to L-1. Mutation testing found the upper bound
    could be weakened from ``>= L`` to ``> L`` with no test failing, which
    would let a bead at exactly L through to die later with a raw IndexError
    (or, for a negative coordinate, to wrap silently onto a real site).
    """
    path = tmp_path / "restart.pimms"
    path.write_bytes(pickle.dumps(_restart_payload(box, position)))
    with pytest.raises(RestartException):
        RestartObject().build_from_file(str(path))


@pytest.mark.parametrize(
    "box, position", [([4, 5], [3, 4]), ([4, 5], [0, 0]), ([4, 5, 6], [3, 4, 5])]
)
def test_restart_file_with_a_bead_on_the_last_site_is_accepted(tmp_path, box, position):
    """The other side of the bound: coordinates 0 and L-1 are inside the box."""
    path = tmp_path / "restart.pimms"
    path.write_bytes(pickle.dumps(_restart_payload(box, position)))
    restart = RestartObject()
    restart.build_from_file(str(path))
    assert restart.chains[1][0] == [list(position)]


@pytest.mark.parametrize(
    "offset",
    [[1, 0.5], [0.5, 1], [1, True], [1, "1"], [1, None], [1.0, 1.0], [1], [1, 1, 1]],
)
def test_restart_resize_refuses_an_offset_that_is_not_all_integers(offset):
    """A manual resize offset must be one integer per dimension - every entry.

    Mutation testing found the per-entry check could be turned from "any entry
    is bad" into "all entries are bad" unnoticed, which accepts an offset with
    one good and one bad entry and truncates the bad one to an integer. A
    refused resize must also leave the restart object as it was.
    """
    restart = RestartObject()
    restart.dimensions = [3, 3]
    restart.chains = {1: [[[0, 0], [1, 0]], "AB", 0]}
    with pytest.raises(RestartException):
        restart.update_lattice_dimensions([6, 6], manual_offset=offset)
    assert restart.dimensions == [3, 3]
    assert restart.chains[1][0] == [[0, 0], [1, 0]]


@pytest.mark.parametrize(
    "specification", [[2, "AB", "CD"], [2, "AB", 5], [2], [], (1, "A", "A", "A")]
)
def test_extra_chains_specification_of_the_wrong_length_is_refused(specification):
    """EXTRA_CHAINS is a count and a sequence - exactly two entries.

    With the length check removed a three-entry specification was accepted and
    its trailing entry silently dropped; no test noticed.
    """
    restart = RestartObject()
    with pytest.raises(RestartException, match="EXTRA_CHAINS"):
        restart.add_extra_chains(specification)
    assert restart.extra_chains == {}


@pytest.mark.parametrize("argument", ["int_seq", "LR_int_seq"])
@pytest.mark.parametrize(
    "bad_code",
    [True, np.bool_(False), 1.5, 2.0, "1", None, 2**31, -(2**31) - 1],
    ids=[
        "bool",
        "numpy-bool",
        "float",
        "integral-float",
        "str",
        "None",
        "above-int32",
        "below-int32",
    ],
)
def test_chain_refuses_interaction_codes_that_are_not_int32_integers(
    argument, bad_code
):
    """Every interaction code of a chain must be an integer that fits in int32.

    The codes are written into the lattice type grid, so a boolean, a float or
    an integer outside the int32 range must be refused when the chain is built.
    Mutation testing found the whole check could be removed with no test
    failing.
    """
    arguments = dict(
        lattice_grid=np.zeros((8, 8), dtype=np.int32),
        dimensions=[8, 8],
        sequence="ABCD",
        int_seq=[1, 2, 3, 4],
        LR_int_seq=[1, 2, 3, 4],
        LR_IDX=[0, 2],
        chainID=3,
        chainType=1,
        chain_positions=[[0, 0], [0, 1], [0, 2], [0, 3]],
    )
    Chain(**arguments)  # the unmodified arguments are valid
    codes = list(arguments[argument])
    codes[2] = bad_code
    arguments[argument] = codes
    with pytest.raises(ChainInitializationException, match=argument):
        Chain(**arguments)


def test_chain_accepts_the_extreme_int32_interaction_codes():
    """The other side of the range check: the int32 limits themselves are valid."""
    limits = np.iinfo(np.int32)
    chain = Chain(
        lattice_grid=np.zeros((8, 8), dtype=np.int32),
        dimensions=[8, 8],
        sequence="AB",
        int_seq=[int(limits.max), np.int64(7)],
        LR_int_seq=[int(limits.min), 0],
        LR_IDX=[],
        chainID=1,
        chainType=0,
        chain_positions=[[0, 0], [0, 1]],
    )
    assert chain.get_intcode_sequence() == [int(limits.max), 7]


# ---------------------------------------------------------------------------
# reporting logic
# ---------------------------------------------------------------------------


def test_performance_file_is_written_every_twentieth_of_the_run_and_at_step_20(
    tmp_path,
):
    """PERFORMANCE.dat gets a row every 1/20th of the run, and one at step 20.

    With N_STEPS = 120 a twentieth is 6 steps, so the rows are at the multiples
    of 6 plus the early estimate at step 20 (which is not a multiple of 6, so
    both halves of the rule are visible). Mutation testing found the cadence
    could be shifted to one step after each multiple with no test failing.
    """
    n_steps = 120
    state = _build(
        tmp_path,
        "SR",
        False,
        {"MOVE_CRANKSHAFT": 1.0},
        [12, 12, 12],
        [(4, "AABBA")],
        50,
        n_steps=n_steps,
        equilibration=2,
        extra={
            "PRINT_FREQ": 100000,
            "XTC_FREQ": 100000,
            "ENERGY_CHECK": 0,
            "ANALYSIS_FREQ": 100000,
            "RESTART_FREQ": 100000,
            "EN_FREQ": 100000,
        },
    )
    with _working_directory(tmp_path), contextlib.redirect_stdout(io.StringIO()):
        state.sim.run_simulation()

    lines = (tmp_path / "PERFORMANCE.dat").read_text().splitlines()
    assert lines[0].startswith("Step\t")
    written = [int(line.split("\t")[0]) for line in lines[1:]]
    twentieth = n_steps // 20
    expected = sorted(set(range(twentieth, n_steps + 1, twentieth)) | {20})
    assert [step for step in written if step > 0] == expected


def _rod(chain_id, start, axis, length, box):
    """Bookkeeping rows of a straight chain laid along one axis.

    Parameters
    ----------
    chain_id : int
        chainID written into column 4.

    start : sequence of int
        Position of the first bead.

    axis : int
        Axis the rod runs along (wrapping through the periodic boundary).

    length : int
        Number of beads.

    box : sequence of int
        Box dimensions.

    Returns
    -------
    numpy.ndarray
        ``int64`` rows of the bead table: chainID in column 4, the position
        from column 5 on, zeros elsewhere.
    """
    rows = np.zeros((length, 5 + len(box)), dtype=np.int64)
    rows[:, 4] = chain_id
    rows[:, 5:] = start
    rows[:, 5 + axis] = (start[axis] + np.arange(length)) % box[axis]
    return rows


def test_parallel_fit_report_measures_extents_per_axis_and_counts_chains():
    """The fit report's per-axis extents and its too-extended count, by hand.

    A 60 x 60 x 20 short-range box splits into 4 x 4 x 1 blocks of 15 x 15 x 20
    with a halo of 3, so a chain fits a block interior if it spans at most
    15 - 2*3 = 9 sites on x and on y; z is not split and has no limit. Five
    straight rods are laid down:

    ======  ====  ======  ==================  ============
    chain   axis  length  extent (x, y, z)    too extended
    ======  ====  ======  ==================  ============
    1       x     12      (12, 1, 1)          yes
    2       x     11      (11, 1, 1)          yes
    3       x     10      (10, 1, 1)          yes
    4       z     15      (1, 1, 15)          no (z is not split)
    5       y     4       (1, 4, 1)           no (wraps through y = 0)
    ======  ====  ======  ==================  ============

    so the largest extent per axis is (12, 4, 15) and three chains are too
    extended. The numbers are deliberately different on every axis and the
    number of chains differs from the number of axes: mutation testing found
    that taking the maximum over the wrong axis of the (chain, axis) table, or
    counting axes instead of chains, passed every test.
    """
    box = [60, 60, 20]
    layout = moves.mega_crank_fast.parallel_layout_info(60, 60, 20, False)
    assert layout["blocks"] == (4, 4, 1) and layout["block_size"] == (15, 15, 20)
    assert layout["W"] == 3

    rods = [
        _rod(1, [5, 5, 3], 0, 12, box),
        _rod(2, [30, 9, 3], 0, 11, box),
        _rod(3, [5, 20, 7], 0, 10, box),
        _rod(4, [40, 40, 2], 2, 15, box),
        _rod(5, [50, 58, 10], 1, 4, box),
    ]
    idx = np.vstack(rods)
    lengths = np.array([len(rod) for rod in rods], dtype=np.int32)
    offsets = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int32)

    report = moves.parallel_chain_fit_report(idx, offsets, lengths, box, False)
    assert tuple(report["max_extent_by_axis"]) == (12, 4, 15)
    assert report["max_extent"] == 15
    assert report["n_too_extended"] == 3
    assert report["interiors"] == (9, 9, 20) and report["interior"] == 9
    assert (
        report["n_chains"] == 5
        and report["n_frozen"] == 0
        and report["n_over_cap"] == 0
    )
    assert report["ok"] is False

    # freezing the three long x rods removes them from the description
    frozen = moves.parallel_chain_fit_report(
        idx, offsets, lengths, box, False, frozen_chains=(1, 2, 3)
    )
    assert tuple(frozen["max_extent_by_axis"]) == (1, 4, 15)
    assert frozen["n_too_extended"] == 0 and frozen["n_frozen"] == 3
    assert frozen["ok"] is True
