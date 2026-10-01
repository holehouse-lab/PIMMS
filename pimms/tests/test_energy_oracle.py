"""
The Hamiltonian checked against the definition of the energy model.

Every other energy test in the suite compares one PIMMS energy with another -
usually a kernel's incremental energy with ``Hamiltonian.evaluate_total_energy``.
Those two share their conventions (walls as solvent, the Chebyshev shells, the
angle classification, the pair extractors), so a convention that was wrong in
both would pass every one of them. Here the energy is recomputed from the model
as the documentation states it (``docs/overview.rst``, ``docs/input_files.rst``)
in plain numpy, without calling any PIMMS energy routine, pair extractor or
lattice utility, and compared term by term with ``evaluate_total_energy``.

The model, as the oracle implements it:

* **Distances** are Chebyshev distances: minimum image under periodic boundaries,
  plain coordinate differences under a hardwall (nothing wraps).
* **Short range (SR).** Every pair of beads at distance 1 scores its SR energy
  once, bonded pairs included. Every one of a bead's ``3^d - 1`` neighbour sites
  that holds no bead scores that bead's solvation energy (``R 0``); under a
  hardwall this includes the neighbour sites beyond a wall, because a wall is
  solvent.
* **Long range (LR) and super-long range (SLR).** Every pair of long-range beads
  at distance exactly 2 (LR) or exactly 3 (SLR) scores the pair's LR / SLR
  energy. A residue is long-range if it appears on any 4- or 5-column line, a
  pair of long-range residues with no line of its own scores 0, a 4-column line
  has an SLR energy of 0, and there is no solvent term at either range.
* **Angles.** For every interior bead ``i`` the displacement ``d = p[i+1] -
  p[i-1]`` (minimum image under periodic boundaries) picks the class: A1 if
  every ``|d_k| <= 1``, A3 if every non-zero ``|d_k|`` is 2, A2 otherwise; the
  penalty is that class's value for bead ``i``'s residue. ``ANGLE_PENALTY_T_NORM``
  values are multiplied by the temperature and rounded to the nearest integer.

The sweep covers 2D and 3D, SR-only, LR and SLR parameter sets, periodic and
hardwall boundaries, the smallest legal box (7 per axis) and non-cubic boxes, a
chain longer than the box, short chains and monomers. Two deliberately perturbed
oracles (walls not counted as solvent; no minimum image) are shown to disagree,
so the comparison is able to fail.
"""

from __future__ import annotations

import itertools
import random
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

import numpy as np
import numpy.typing as npt
import pytest

from pimms.chain import Chain
from pimms.energy import Hamiltonian
from pimms.lattice import Lattice

IntArray = npt.NDArray[np.int64]
ChainSpec = tuple[str, list[list[int]]]
EnergyTerms = tuple[int, int, int, int, int]

RESIDUES: tuple[str, ...] = ("A", "B", "C")

# residues allowed to carry LR/SLR lines; "C" never does, so every LR sweep case
# mixes long-range and short-range-only beads
LONG_RANGE_CANDIDATES: frozenset[str] = frozenset({"A", "B"})

BOXES: dict[str, list[int]] = {
    "2D-7": [7, 7],
    "2D-7x11": [7, 11],
    "3D-7": [7, 7, 7],
    "3D-7x9x8": [7, 9, 8],
}


# ---------------------------------------------------------------------------
# parameter sets
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ForceField:
    """A parameter set, as parameter-file text and as the tables the oracle reads.

    Attributes
    ----------
    text : str
        The parameter file handed to :class:`pimms.energy.Hamiltonian`.
    sr : dict of (str, str) to int
        Symmetric SR energies keyed by residue letters, including the
        solvation energies under ``(residue, "0")`` and ``("0", residue)``.
    lr : dict of (str, str) to int
        Symmetric LR energies of the pairs that have their own 4- or 5-column
        line. Long-range pairs without a line are absent (they score 0).
    slr : dict of (str, str) to int
        Symmetric SLR energies, keyed like ``lr`` (0 for a 4-column line).
    long_range : frozenset of str
        The residues that appear on any 4- or 5-column line.
    angles : dict of str to tuple of int
        The applied ``(A1, A2, A3)`` penalties per residue, after temperature
        scaling and rounding for ``ANGLE_PENALTY_T_NORM`` lines.
    temperature : float
        The temperature the Hamiltonian is built at.
    """

    text: str
    sr: dict[tuple[str, str], int]
    lr: dict[tuple[str, str], int]
    slr: dict[tuple[str, str], int]
    long_range: frozenset[str]
    angles: dict[str, tuple[int, int, int]]
    temperature: float


def _nonzero(rng: np.random.Generator, low: int, high: int) -> int:
    """Draw an integer in [low, high] that is not zero.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    low, high : int
        Inclusive bounds, with low < 0 < high.

    Returns
    -------
    int
        A non-zero integer, so that a term the sweep relies on cannot vanish by
        drawing a zero energy.
    """
    while True:
        value = int(rng.integers(low, high + 1))
        if value != 0:
            return value


def _t_norm_penalty(rng: np.random.Generator, temperature: float) -> tuple[str, int]:
    """Draw a temperature-normalised angle penalty and its applied value.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    temperature : float
        The temperature the penalty is scaled by.

    Returns
    -------
    tuple of (str, int)
        The value as written in the parameter file and the integer penalty the
        model applies, ``round(value * temperature)``. Values whose scaled penalty
        sits on a half-integer are redrawn, so the rounding convention for ties
        never matters.
    """
    while True:
        token = "%.2f" % rng.uniform(-0.2, 0.5)
        scaled = float(token) * temperature
        if abs(abs(scaled - np.floor(scaled)) - 0.5) > 1e-6:
            return token, int(round(scaled))


def make_forcefield(rng: np.random.Generator, kind: str, angle_mode: str,
                    temperature: float) -> ForceField:
    """Draw a random parameter set over RESIDUES.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    kind : str
        ``"SR"`` (3-column lines only), ``"LR"`` (4-column lines for some
        long-range pairs) or ``"SLR"`` (5-column lines for them). For LR and SLR
        the A-A and B-B pairs always get a long-range line and A-B only
        sometimes, so a long-range pair without a line of its own is covered.
    angle_mode : str
        ``"abs"`` for ``ANGLE_PENALTY`` lines, ``"tnorm"`` for
        ``ANGLE_PENALTY_T_NORM`` lines.
    temperature : float
        Temperature used to scale T-normalised penalties.

    Returns
    -------
    ForceField
        The parameter-file text together with the oracle's tables.
    """
    lines: list[str] = []
    sr: dict[tuple[str, str], int] = {}
    lr: dict[tuple[str, str], int] = {}
    slr: dict[tuple[str, str], int] = {}
    angles: dict[str, tuple[int, int, int]] = {}

    for residue in RESIDUES:
        solvation = _nonzero(rng, -4, 4)
        sr[(residue, "0")] = sr[("0", residue)] = solvation
        lines.append("%s 0 %i" % (residue, solvation))

    for r1, r2 in itertools.combinations_with_replacement(RESIDUES, 2):
        e_sr = int(rng.integers(-6, 5))
        sr[(r1, r2)] = sr[(r2, r1)] = e_sr
        long_range_pair = (kind != "SR" and r1 in LONG_RANGE_CANDIDATES and
                           r2 in LONG_RANGE_CANDIDATES and
                           (r1 == r2 or rng.random() < 0.5))
        if not long_range_pair:
            lines.append("%s %s %i" % (r1, r2, e_sr))
            continue
        e_lr = _nonzero(rng, -5, 5)
        lr[(r1, r2)] = lr[(r2, r1)] = e_lr
        if kind == "SLR":
            e_slr = _nonzero(rng, -4, 4)
            lines.append("%s %s %i %i %i" % (r1, r2, e_sr, e_lr, e_slr))
        else:
            e_slr = 0
            lines.append("%s %s %i %i" % (r1, r2, e_sr, e_lr))
        slr[(r1, r2)] = slr[(r2, r1)] = e_slr

    for residue in RESIDUES:
        if angle_mode == "tnorm":
            drawn = [_t_norm_penalty(rng, temperature) for _ in range(3)]
            lines.append("ANGLE_PENALTY_T_NORM %s %s %s %s"
                         % (residue, drawn[0][0], drawn[1][0], drawn[2][0]))
            angles[residue] = (drawn[0][1], drawn[1][1], drawn[2][1])
        else:
            a1, a2, a3 = (int(rng.integers(-3, 11)) for _ in range(3))
            lines.append("ANGLE_PENALTY %s %i %i %i" % (residue, a1, a2, a3))
            angles[residue] = (a1, a2, a3)

    long_range = frozenset(residue for pair in lr for residue in pair)
    return ForceField("\n".join(lines) + "\n", sr, lr, slr, long_range, angles,
                      temperature)


# ---------------------------------------------------------------------------
# the oracle
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OracleEnergy:
    """The model energy of a configuration, with counts of what it scored.

    Attributes
    ----------
    sr, lr, slr, angle : int
        The four energy components, in the order evaluate_total_energy reports
        them.
    wall_contacts : int
        Number of (bead, neighbour site beyond a wall) contacts scored as solvent.
    lr_pairs, slr_pairs : int
        Number of long-range bead pairs at distance 2 and 3.
    angle_triplets : int
        Number of interior beads that carry an angle penalty.
    """

    sr: int
    lr: int
    slr: int
    angle: int
    wall_contacts: int
    lr_pairs: int
    slr_pairs: int
    angle_triplets: int

    def terms(self) -> EnergyTerms:
        """Return ``(total, SR, LR, SLR, angle)``, the evaluate_total_energy layout.

        Returns
        -------
        tuple of int
            The total followed by the four components.
        """
        return (self.sr + self.lr + self.slr + self.angle,
                self.sr, self.lr, self.slr, self.angle)


def _separation(delta: IntArray, dims: IntArray, periodic: bool) -> IntArray:
    """Map coordinate differences onto the separations the model uses.

    Parameters
    ----------
    delta : numpy.ndarray
        Integer coordinate differences, last axis over dimensions.
    dims : numpy.ndarray
        Box size per axis.
    periodic : bool
        If True, apply the minimum-image convention, giving each component in
        ``[-(L // 2), L - L // 2 - 1]``; otherwise return delta unchanged.

    Returns
    -------
    numpy.ndarray
        The separations, same shape as delta.
    """
    if not periodic:
        return delta
    return (delta + dims // 2) % dims - dims // 2


def _flatten(chains: list[ChainSpec], n_dim: int) -> tuple[list[str], IntArray]:
    """Concatenate every chain's residues and positions.

    Parameters
    ----------
    chains : list of (str, list of list of int)
        The configuration, one (sequence, positions) entry per chain.
    n_dim : int
        Number of dimensions (2 or 3).

    Returns
    -------
    tuple of (list of str, numpy.ndarray)
        The residue letter of every bead and an ``(n_beads, n_dim)`` array of
        their positions.
    """
    residues = [residue for sequence, _ in chains for residue in sequence]
    positions = np.array([p for _, chain_positions in chains for p in chain_positions],
                         dtype=np.int64).reshape(-1, n_dim)
    return residues, positions


def check_configuration(chains: list[ChainSpec], dims: list[int], hardwall: bool) -> None:
    """Assert the configuration is one the model is defined on.

    Parameters
    ----------
    chains : list of (str, list of list of int)
        The configuration, one (sequence, positions) entry per chain.
    dims : list of int
        Box size per axis.
    hardwall : bool
        Boundary convention; under a hardwall no bond may cross a face.

    Returns
    -------
    None
        Returns only if every bead is inside the box, no site holds two beads and
        every bond joins Chebyshev neighbours.
    """
    box = np.asarray(dims, dtype=np.int64)
    _, positions = _flatten(chains, len(dims))
    assert np.all(positions >= 0) and np.all(positions < box), "bead outside the box"
    assert len({tuple(p) for p in positions.tolist()}) == len(positions), "overlap"
    for sequence, chain_positions in chains:
        assert len(sequence) == len(chain_positions)
        if len(sequence) == 1:
            continue
        p = np.asarray(chain_positions, dtype=np.int64)
        bonds = _separation(p[1:] - p[:-1], box, not hardwall)
        assert np.all(np.abs(bonds).max(axis=1) == 1), "broken bond"


def oracle_energy(chains: list[ChainSpec], dims: list[int], hardwall: bool,
                  ff: ForceField, perturbation: Optional[str] = None) -> OracleEnergy:
    """Energy of a configuration straight from the model definition.

    Parameters
    ----------
    chains : list of (str, list of list of int)
        The configuration, one (sequence, positions) entry per chain.
    dims : list of int
        Box size per axis (each >= 7).
    hardwall : bool
        Boundary convention: False for periodic, True for hardwall.
    ff : ForceField
        The parameter set.
    perturbation : str, optional
        Deliberately wrong variants of the model, used as positive controls:
        ``"walls_not_solvent"`` scores only the in-box neighbour sites of a bead
        as solvent, ``"no_minimum_image"`` measures every separation (pairs and
        angles) without the periodic wrap. Default None (the model).

    Returns
    -------
    OracleEnergy
        The SR, LR, SLR and angle energies and the counts of what was scored.
    """
    n_dim = len(dims)
    box = np.asarray(dims, dtype=np.int64)
    periodic = not hardwall and perturbation != "no_minimum_image"
    residues, positions = _flatten(chains, n_dim)
    n_beads = len(residues)

    # Chebyshev distance between every pair of beads
    delta = positions[:, None, :] - positions[None, :, :]
    chebyshev = np.abs(_separation(delta, box, periodic)).max(axis=2)

    sr = lr = slr = 0
    lr_pairs = slr_pairs = 0
    first, second = np.nonzero(np.triu(chebyshev <= 3, k=1))
    for i, j in zip(first.tolist(), second.tolist()):
        pair = (residues[i], residues[j])
        distance = int(chebyshev[i, j])
        if distance == 1:
            sr += ff.sr[pair]
        elif pair[0] in ff.long_range and pair[1] in ff.long_range:
            if distance == 2:
                lr += ff.lr.get(pair, 0)
                lr_pairs += 1
            else:
                slr += ff.slr.get(pair, 0)
                slr_pairs += 1

    # Every neighbour site that holds no bead is solvent. Under a hardwall a site
    # beyond a wall never holds a bead, so it is counted here as solvent too; we
    # also count those sites separately, for the coverage checks and the
    # "walls_not_solvent" control.
    bead_neighbours = (chebyshev == 1).sum(axis=1)
    all_sites = 3 ** n_dim - 1
    wall_contacts = 0
    for i in range(n_beads):
        in_box = 1
        for k in range(n_dim):
            in_box *= 3 - int(positions[i, k] == 0) - int(positions[i, k] == dims[k] - 1)
        in_box -= 1
        beyond_wall = all_sites - in_box if hardwall else 0
        wall_contacts += beyond_wall
        solvent_sites = all_sites - int(bead_neighbours[i])
        if perturbation == "walls_not_solvent":
            solvent_sites -= beyond_wall
        sr += solvent_sites * ff.sr[(residues[i], "0")]

    angle = 0
    angle_triplets = 0
    for sequence, chain_positions in chains:
        p = np.asarray(chain_positions, dtype=np.int64).reshape(-1, n_dim)
        for i in range(1, len(sequence) - 1):
            d = np.abs(_separation(p[i + 1] - p[i - 1], box, periodic))
            if d.max() <= 1:
                penalty_class = 0
            elif np.all((d == 0) | (d == 2)):
                penalty_class = 2
            else:
                penalty_class = 1
            angle += ff.angles[sequence[i]][penalty_class]
            angle_triplets += 1

    return OracleEnergy(sr, lr, slr, angle, wall_contacts, lr_pairs, slr_pairs,
                        angle_triplets)


# ---------------------------------------------------------------------------
# building PIMMS objects for a configuration
# ---------------------------------------------------------------------------

@pytest.fixture
def workdir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Path]:
    """Run in a scratch directory with the global ``random`` state restored after.

    The parameter-file parser writes copies of the file into the working
    directory, and de-novo chain placement draws from Python's ``random``, which
    the tests seed.

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytest's per-test scratch directory.
    monkeypatch : pytest.MonkeyPatch
        Used to change into tmp_path for the duration of the test.

    Yields
    ------
    pathlib.Path
        The scratch directory, which is also the working directory.
    """
    state = random.getstate()
    monkeypatch.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        random.setstate(state)


def build_hamiltonian(ff: ForceField, n_dim: int, hardwall: bool) -> Hamiltonian:
    """Write the parameter file into the working directory and build a Hamiltonian.

    Parameters
    ----------
    ff : ForceField
        The parameter set.
    n_dim : int
        Number of dimensions (2 or 3).
    hardwall : bool
        Boundary convention.

    Returns
    -------
    pimms.energy.Hamiltonian
        The Hamiltonian under test, with angles on.
    """
    Path("params.prm").write_text(ff.text)
    return Hamiltonian("params.prm", n_dim, False, False, hardwall=hardwall,
                       temperature=ff.temperature, reduced_printing=True)


def build_lattice(ham: Hamiltonian, dims: list[int], hardwall: bool,
                  chains: list[ChainSpec]) -> Lattice:
    """Build a Lattice holding exactly the given chains.

    Parameters
    ----------
    ham : pimms.energy.Hamiltonian
        Provides the residue integer codes.
    dims : list of int
        Box size per axis.
    hardwall : bool
        Boundary convention.
    chains : list of (str, list of list of int)
        The configuration, one (sequence, positions) entry per chain.

    Returns
    -------
    pimms.lattice.Lattice
        A fully specified Lattice (chains, occupancy grid and type grid).
    """
    grid = np.zeros(dims, dtype=np.int32)
    type_grid = np.zeros(dims, dtype=np.int32)
    chain_objects: dict[int, Chain] = {}
    for chain_id, (sequence, positions) in enumerate(chains, start=1):
        codes = ham.convert_sequence_to_integer_sequence(sequence)
        chain_objects[chain_id] = Chain(
            grid, dims, sequence, codes,
            ham.convert_sequence_to_LR_integer_sequence(sequence),
            ham.get_indices_of_long_range_residues(sequence),
            chain_id, chain_id - 1, chain_positions=positions, hardwall=hardwall)
        for position, code in zip(positions, codes):
            grid[tuple(position)] = chain_id
            type_grid[tuple(position)] = code
    return Lattice(dims, [], ham, 3.6, chainsDict=chain_objects, lattice_grid=grid,
                   type_grid=type_grid, hardwall=hardwall)


def lattice_chains(lat: Lattice) -> list[ChainSpec]:
    """Read the configuration back out of a Lattice.

    Parameters
    ----------
    lat : pimms.lattice.Lattice
        The lattice.

    Returns
    -------
    list of (str, list of list of int)
        One (sequence, positions) entry per chain, in chainID order.
    """
    return [(lat.chains[c].sequence, [[int(x) for x in p] for p in lat.chains[c].positions])
            for c in sorted(lat.chains)]


def pimms_terms(ham: Hamiltonian, lat: Lattice) -> EnergyTerms:
    """Return evaluate_total_energy's ``(total, SR, LR, SLR, angle)`` as ints.

    Parameters
    ----------
    ham : pimms.energy.Hamiltonian
        The Hamiltonian under test.
    lat : pimms.lattice.Lattice
        The configuration.

    Returns
    -------
    tuple of int
        The five terms.
    """
    total, sr, lr, slr, angle = ham.evaluate_total_energy(lat)
    return (int(total), int(sr), int(lr), int(slr), int(angle))


def random_system(rng: np.random.Generator, dims: list[int]) -> list[tuple[int, str]]:
    """A keyfile-style chain list filling about a third of the box.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of randomness.
    dims : list of int
        Box size per axis.

    Returns
    -------
    list of (int, str)
        ``(count, sequence)`` entries: one chain twice as long as the longest
        axis (placed first), a few chains of 3-7 beads, and monomers of each
        residue type for the rest.
    """
    def sequence(length: int) -> str:
        """Draw a random sequence of the given length over RESIDUES.

        Parameters
        ----------
        length : int
            Number of residues.

        Returns
        -------
        str
            The sequence.
        """
        return "".join(rng.choice(list(RESIDUES), size=length).tolist())

    n_sites = int(np.prod(dims))
    system = [(1, sequence(2 * max(dims)))]
    for _ in range(max(2, n_sites // 60)):
        system.append((1, sequence(int(rng.integers(3, 8)))))
    n_monomers = max(3, int(0.35 * n_sites) - sum(len(s) for _, s in system))
    counts = rng.multinomial(n_monomers, [1.0 / len(RESIDUES)] * len(RESIDUES))
    system.extend((int(n), residue) for n, residue in zip(counts, RESIDUES) if n > 0)
    return system


# ---------------------------------------------------------------------------
# the sweep
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hardwall", [False, True], ids=["pbc", "hardwall"])
@pytest.mark.parametrize("kind", ["SR", "LR", "SLR"])
@pytest.mark.parametrize("box", sorted(BOXES))
def test_total_energy_matches_the_model_definition(workdir: Path, box: str, kind: str,
                                                   hardwall: bool) -> None:
    """evaluate_total_energy equals the model energy, term by term, on random systems.

    Three configurations per case are placed de novo by the Lattice (seeded), each
    with a chain twice as long as the longest axis, a few short chains and
    monomers, at about a third occupancy.

    Parameters
    ----------
    workdir : pathlib.Path
        Scratch working directory (fixture).
    box : str
        Key into BOXES.
    kind : str
        Parameter-set kind passed to make_forcefield ("SR", "LR" or "SLR").
    hardwall : bool
        Boundary convention.

    Returns
    -------
    None
    """
    dims = BOXES[box]
    seed = zlib.crc32(("%s-%s-%s" % (box, kind, hardwall)).encode())
    rng = np.random.default_rng(seed)
    angle_mode = "tnorm" if kind == "LR" else "abs"
    ff = make_forcefield(rng, kind, angle_mode, temperature=float(rng.choice([37.0, 55.0])))
    ham = build_hamiltonian(ff, len(dims), hardwall)

    scored = []
    for replicate in range(3):
        random.seed(seed + replicate)
        lat = Lattice(dims, random_system(rng, dims), ham, 3.6, hardwall=hardwall)
        chains = lattice_chains(lat)
        check_configuration(chains, dims, hardwall)
        expected = oracle_energy(chains, dims, hardwall, ff)
        assert pimms_terms(ham, lat) == expected.terms(), (box, kind, hardwall, replicate)
        scored.append(expected)

    # the sweep must actually exercise what it claims to
    assert all(e.angle_triplets > 0 for e in scored)
    if hardwall:
        assert all(e.wall_contacts > 0 for e in scored)
    if kind != "SR":
        assert all(e.lr_pairs > 0 for e in scored)
    if kind == "SLR":
        assert all(e.slr_pairs > 0 and e.slr != 0 for e in scored)


# ---------------------------------------------------------------------------
# hand-built edge cases
# ---------------------------------------------------------------------------

def _fixed_forcefield(kind: str) -> ForceField:
    """A deterministic SLR-capable parameter set for the hand-built cases.

    Parameters
    ----------
    kind : str
        Passed to make_forcefield.

    Returns
    -------
    ForceField
        Every solvation, LR and SLR energy is non-zero.
    """
    return make_forcefield(np.random.default_rng(20260926), kind, "abs", 40.0)


# In a box of 7, long-range monomers A/B placed so that most pairs are within
# range only through the periodic wrap: [0,0]-[4,0] and [4,0]-[0,5] are 3 apart
# only across the x face, [0,0]-[0,5] is 2 apart only across the y face, and
# [3,3] is 3 from each of the other three. The C at [6,6] (not long-range)
# touches [0,0] and [0,5] only across the corner. Under a hardwall only the three
# pairs with [3,3] remain in range.
_WRAP_SHELLS_2D: list[ChainSpec] = [("A", [[0, 0]]), ("B", [[4, 0]]), ("A", [[0, 5]]),
                                    ("B", [[3, 3]]), ("C", [[6, 6]])]
# a 3D chain with a bond that crosses all three faces at once, so two of its
# angle displacements are (2, 2, 2) only after the minimum image
_FACE_CROSSING_3D: list[ChainSpec] = [("ABCAB", [[5, 5, 5], [6, 6, 6], [0, 0, 0],
                                                 [1, 1, 1], [2, 1, 0]]),
                                      ("A", [[3, 5, 6]]), ("BA", [[6, 3, 3], [0, 3, 4]])]
# a hardwall chain lying against the z = 0 face of a non-cubic 3D box, and
# monomers against a face and in the far corner
_ALONG_WALLS_3D: list[ChainSpec] = [("ABABCA", [[0, 0, 0], [1, 0, 0], [2, 1, 0], [3, 2, 0],
                                                [3, 3, 1], [2, 4, 2]]),
                                    ("B", [[0, 2, 0]]), ("A", [[6, 8, 7]])]

_EDGE_CASES: dict[str, tuple[list[int], bool, list[ChainSpec]]] = {
    "wrap_shells_pbc": ([7, 7], False, _WRAP_SHELLS_2D),
    "wrap_shells_hardwall": ([7, 7], True, _WRAP_SHELLS_2D),
    "face_crossing_3d_pbc": ([7, 7, 7], False, _FACE_CROSSING_3D),
    "along_walls_3d_hardwall": ([7, 9, 8], True, _ALONG_WALLS_3D),
}


@pytest.mark.parametrize("case", sorted(_EDGE_CASES))
def test_hand_built_configurations_match_the_model_definition(workdir: Path,
                                                              case: str) -> None:
    """evaluate_total_energy equals the model energy on hand-placed edge cases.

    Parameters
    ----------
    workdir : pathlib.Path
        Scratch working directory (fixture).
    case : str
        Key into _EDGE_CASES.

    Returns
    -------
    None
    """
    dims, hardwall, chains = _EDGE_CASES[case]
    ff = _fixed_forcefield("SLR")
    ham = build_hamiltonian(ff, len(dims), hardwall)
    lat = build_lattice(ham, dims, hardwall, chains)
    check_configuration(lattice_chains(lat), dims, hardwall)

    expected = oracle_energy(chains, dims, hardwall, ff)
    assert pimms_terms(ham, lat) == expected.terms()

    # the oracle sees the pairs laid out in the comment above _WRAP_SHELLS_2D
    if case == "wrap_shells_pbc":
        assert (expected.lr_pairs, expected.slr_pairs) == (1, 5)
    if case == "wrap_shells_hardwall":
        assert (expected.lr_pairs, expected.slr_pairs) == (0, 3)


# ---------------------------------------------------------------------------
# positive controls: a wrong model must be caught
# ---------------------------------------------------------------------------

_CONTROLS: dict[str, tuple[str, list[int], bool, list[ChainSpec]]] = {
    # A monomer in a corner, one against a face and a chain along a wall. Every
    # bead is an A, so ignoring the walls shifts the SR energy by (number of
    # beyond-wall sites) x (A's solvation energy), which cannot vanish.
    "walls_not_solvent_2d": ("walls_not_solvent", [7, 7], True,
                             [("A", [[0, 0]]), ("A", [[3, 6]]), ("AAA", [[6, 1], [6, 2], [6, 3]])]),
    "walls_not_solvent_3d": ("walls_not_solvent", [7, 9, 8], True,
                             [("A", [[0, 0, 0]]), ("AA", [[3, 8, 4], [4, 8, 5]])]),
    "no_minimum_image_2d": ("no_minimum_image", [7, 7], False, _WRAP_SHELLS_2D),
    "no_minimum_image_3d": ("no_minimum_image", [7, 7, 7], False, _FACE_CROSSING_3D),
}


@pytest.mark.parametrize("control", sorted(_CONTROLS))
def test_perturbed_oracle_disagrees(workdir: Path, control: str) -> None:
    """A deliberately wrong model disagrees with PIMMS where the right one agrees.

    This is what shows the comparison can fail: ignoring the walls as solvent
    (hardwall cases) or dropping the minimum image (periodic cases) changes the
    oracle's answer on configurations where PIMMS and the model agree.

    Parameters
    ----------
    workdir : pathlib.Path
        Scratch working directory (fixture).
    control : str
        Key into _CONTROLS.

    Returns
    -------
    None
    """
    perturbation, dims, hardwall, chains = _CONTROLS[control]
    ff = _fixed_forcefield("SLR")
    ham = build_hamiltonian(ff, len(dims), hardwall)
    lat = build_lattice(ham, dims, hardwall, chains)

    observed = pimms_terms(ham, lat)
    assert observed == oracle_energy(chains, dims, hardwall, ff).terms()
    assert observed != oracle_energy(chains, dims, hardwall, ff,
                                     perturbation=perturbation).terms()
