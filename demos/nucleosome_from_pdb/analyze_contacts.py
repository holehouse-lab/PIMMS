"""
Where the disordered proteins touch the nucleosome, and each other.

Reads each run of the demo with lemonade and, in every production frame, finds
which nucleosome beads are touched by a bead of an added protein. Two beads
touch when they are Chebyshev neighbours (at most one lattice step apart along
every axis), the shell over which PIMMS's short-range interactions act. The
nucleosome is frozen except its histone tails, which the freeze file leaves
free, so the tails are found again in every frame. In a run with two added
proteins it also finds which copies touch a copy of the other protein, so that
each copy in each frame is free, on the nucleosome, on its partner, or on both
at once. And for each histone tail it finds how much of it lies on the DNA.

For each added protein of each run it writes, into the run's directory
(``<name>`` is ``H1_CTD`` or ``ProTa``):

* ``contacts_nucleosome_<name>.csv``: for every nucleosome residue, the fraction
  of frames in which a bead of the protein touches it, and the mean number of
  the protein's beads touching it;
* ``contacts_<name>.csv``: for every residue of the protein, the fraction of
  frames (over every copy) in which it touches the nucleosome;
* ``1kx5_contacts_<name>.pdb``: the heavy atoms of 1KX5 with each residue's
  contact fraction as its B-factor, for colouring in PyMOL or ChimeraX.

It also writes ``nucleosome_contacts.pdf`` here, comparing the runs, and prints
a summary. Run it from this directory once the simulations have finished::

    python analyze_contacts.py

The beads are mapped back to residues by converting the structure again in
memory, with pdb2pimms's default options; the result is checked against the
restart file, so a restart file made with other options is refused rather than
mislabelled. Which beads are frozen is read from the freeze file, and a run in
which a frozen bead moved is refused.
"""

import argparse
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import matplotlib
import numpy as np
import numpy.typing as npt

from pimms import lemonade
from pimms.data_structures import FreezeFile, contiguous_runs
from pimms.pdb2pimms import ConversionOptions, Residue, convert_structure

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402  (the backend has to be chosen before pyplot is imported)
from matplotlib.axes import Axes  # noqa: E402


IntArray = npt.NDArray[np.int64]
FloatArray = npt.NDArray[np.float64]
BoolArray = npt.NDArray[np.bool_]

# the demo's inputs, relative to this directory
STRUCTURE_FILE = Path('1kx5.cif.gz')
RESTART_FILE = Path('1kx5_restart.pimms')
FREEZE_FILE = Path('1kx5_freeze.txt')

# the files every run directory holds
KEYFILE = 'KEYFILE.kf'
TRAJECTORY = 'traj.xtc'
TOPOLOGY = 'START.pdb'

# The runs, and for each EXTRA_CHAIN line of a run's keyfile (in order): the
# protein's name in the figure, its name in file names, and the number of its
# first residue in its UniProt sequence. The H1.0 C-terminal domain starts
# after the H15 globular domain, residues 24-97 of P07305.
H1_CTD = ('H1.0 CTD', 'H1_CTD', 98)
PROTA = ('ProTα', 'ProTa', 1)
RUNS = {'H1_CTD': (H1_CTD,), 'ProTa': (PROTA,), 'H1_CTD_ProTa': (H1_CTD, PROTA)}

# which 1KX5 chain is which, from the file's _entity and _entity_poly records;
# in its DNA numbering nucleotide n of strand I pairs with nucleotide -n of J
DNA_CHAINS = ('I', 'J')
HISTONE_CHAINS = (('H3', ('A', 'E')), ('H4', ('B', 'F')), ('H2A', ('C', 'G')), ('H2B', ('D', 'H')))

# the 26 Chebyshev neighbours of a lattice site
NEIGHBOUR_OFFSETS = np.array([(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
                              if (dx, dy, dz) != (0, 0, 0)], dtype=np.int64)
EMPTY_SITE = -1

# outputs
NUCLEOSOME_CSV = 'contacts_nucleosome_%s.csv'
PROTEIN_CSV = 'contacts_%s.csv'
BFACTOR_PDB = '1kx5_contacts_%s.pdb'
FIGURE_FILE = Path('nucleosome_contacts.pdf')
TOP_RESIDUES = 10

# figure style: Arial throughout, text kept as text (TrueType) so it stays
# editable in Illustrator, and the xmgrace colour order
FONT = 'Arial'
FONT_SIZE = 8
EDITABLE_FONT_TYPE = 42
COLOURS = ('black', 'red', 'blue', 'green', 'purple')
FIGURE_SIZE_INCHES = (7.0, 6.5)
LINE_WIDTH = 0.8


@dataclass
class AddedProtein:
    """
    The copies of one protein a run added (one EXTRA_CHAIN line).

    Attributes
    ----------
    label : str
        The protein's name in the figure.
    file_name : str
        The protein's name in file names.
    first_residue : int
        The number of its first residue in its UniProt sequence.
    sequence : str
        Its sequence.
    copies : int
        How many copies were added.
    """
    label: str
    file_name: str
    first_residue: int
    sequence: str
    copies: int


@dataclass
class ProteinContacts:
    """
    What one added protein touched over one run.

    Attributes
    ----------
    run : str
        The run directory's name.
    protein : AddedProtein
        The protein.
    partners : list of str
        The other proteins added in the same run.
    n_frames : int
        The number of production frames analysed.
    nucleosome_touched : ndarray of float
        Per nucleosome bead: the fraction of frames in which a bead of the protein touches it.
    nucleosome_pairs : ndarray of float
        Per nucleosome bead: the mean number of the protein's beads touching it.
    residue_touching : ndarray of float
        Per residue of the protein: the fraction of frames, over every copy, in
        which it touches the nucleosome.
    bound : ndarray of float
        Per copy: the fraction of frames in which it touches the nucleosome at all.
    states : dict of str to float
        The fraction of copy-frames in which a copy touches neither the
        nucleosome nor a partner ('free'), only the nucleosome, only a partner,
        or both.
    """
    run: str
    protein: AddedProtein
    partners: list[str]
    n_frames: int
    nucleosome_touched: FloatArray
    nucleosome_pairs: FloatArray
    residue_touching: FloatArray
    bound: FloatArray
    states: dict[str, float]

    @property
    def label(self) -> str:
        """
        The protein's name, with its partners in the run if it had any.

        Returns
        -------
        str
            e.g. ``H1.0 CTD`` or ``H1.0 CTD (with ProTα)``.
        """
        if not self.partners:
            return self.protein.label
        return '%s (with %s)' % (self.protein.label, ', '.join(self.partners))


# .............................................................................
# Reading the inputs
# .............................................................................

def nucleosome_residues(structure: Path, restart: Path) -> list[Residue]:
    """
    Convert the structure again and check it gives the chains of the restart file.

    Parameters
    ----------
    structure : Path
        The structure the restart file was made from.
    restart : Path
        The restart file.

    Returns
    -------
    list of Residue
        One residue per nucleosome bead, in bead order (chainID 1 first).

    Raises
    ------
    SystemExit
        If the chains differ, which means the restart file was made with other options.
    """
    conversion = convert_structure(structure, ConversionOptions())
    with open(restart, 'rb') as handle:
        restart_chains = pickle.load(handle)['CHAINS']
    converted = [segment.sequence for segment in conversion.segments]
    stored = [restart_chains[chain_id][1] for chain_id in sorted(restart_chains)]
    if converted != stored:
        raise SystemExit('%s does not hold the chains a default pdb2pimms conversion of %s gives; convert it '
                         'again with the pdb2pimms command under "Running it" in README.md (which frees the '
                         'histone tails with --mobile), adding --overwrite' % (restart, structure))
    return [residue for segment in conversion.segments for residue in segment.residues]


def restart_sites(restart: Path) -> tuple[IntArray, IntArray]:
    """
    The sites of every bead in the restart file, and its box.

    Parameters
    ----------
    restart : Path
        The restart file.

    Returns
    -------
    tuple of (ndarray, ndarray)
        (n_beads, 3) sites in chainID order, and the box (3,), in lattice units.
    """
    with open(restart, 'rb') as handle:
        payload = pickle.load(handle)
    sites = np.concatenate([np.array(payload['CHAINS'][chain_id][0], dtype=np.int64)
                            for chain_id in sorted(payload['CHAINS'])])
    return sites, np.array(payload['DIMENSIONS'], dtype=np.int64)


def frozen_nucleosome_beads(freeze_file: Path, restart: Path) -> BoolArray:
    """
    Which beads of the restart file the freeze file freezes.

    PIMMS's own freeze-file reader resolves the file, so a bead past the end of
    its chain, or a B line for a chain the restart file does not hold, is an
    error here as it is in a run.

    Parameters
    ----------
    freeze_file : Path
        The freeze file the runs used.
    restart : Path
        The restart file they started from.

    Returns
    -------
    ndarray of bool
        One per restart-file bead, in chainID order: True where the bead is frozen.

    Raises
    ------
    SystemExit
        If the freeze file is missing.
    """
    if not freeze_file.is_file():
        raise SystemExit('%s is missing; convert the structure first (see README.md)' % freeze_file)
    with open(restart, 'rb') as handle:
        chains = pickle.load(handle)['CHAINS']
    lengths = {int(chain_id): len(chains[chain_id][0]) for chain_id in chains}
    frozen_chains, frozen_beads = FreezeFile(str(freeze_file)).resolve(lengths)
    masks = []
    for chain_id in sorted(lengths):
        mask = np.zeros(lengths[chain_id], dtype=bool)
        if chain_id in frozen_beads:
            mask[list(frozen_beads[chain_id])] = True
        elif chain_id in frozen_chains:
            mask[:] = True
        masks.append(mask)
    return np.concatenate(masks)


def added_proteins(run: Path) -> list[AddedProtein]:
    """
    Read the EXTRA_CHAIN lines of a run's keyfile, in order.

    PIMMS numbers the added chains after the restart file's, line by line, so
    the first line's copies come first in the trajectory.

    Parameters
    ----------
    run : Path
        The run directory.

    Returns
    -------
    list of AddedProtein
        One per EXTRA_CHAIN line, named from RUNS (or 'protein N' for a run it does not list).

    Raises
    ------
    SystemExit
        If the keyfile adds nothing, or RUNS names a different number of proteins.
    """
    entries = []
    for line in (run / KEYFILE).read_text().splitlines():
        keyword, _, value = line.split('#')[0].partition(':')
        if keyword.strip() == 'EXTRA_CHAIN':
            count, sequence = value.split()
            # PIMMS upper-cases EXTRA_CHAIN sequences unless CASE_INSENSITIVE_CHAINS is False
            entries.append((sequence.upper(), int(count)))
    if not entries:
        raise SystemExit('%s has no EXTRA_CHAIN line' % (run / KEYFILE))
    names = RUNS.get(run.name, tuple(('protein %d' % i, 'protein_%d' % i, 1) for i in range(1, len(entries) + 1)))
    if len(names) != len(entries):
        raise SystemExit('%s has %d EXTRA_CHAIN lines but %d proteins are named for it in RUNS'
                         % (run / KEYFILE, len(entries), len(names)))
    return [AddedProtein(label=label, file_name=file_name, first_residue=first, sequence=sequence, copies=copies)
            for (label, file_name, first), (sequence, copies) in zip(names, entries)]


def production_frames(run: Path, restart: Path, frozen: BoolArray,
                      n_added: int) -> tuple[IntArray, tuple[int, ...], bool]:
    """
    Load a run's production frames and check the frozen beads stayed where they were put.

    Frame 0 (the starting configuration) is left out; the keyfiles set SAVE_EQ
    to False, so every other frame is a production frame.

    Parameters
    ----------
    run : Path
        The run directory.
    restart : Path
        The restart file the run started from.
    frozen : ndarray of bool
        One per restart-file bead: True where the freeze file freezes it.
    n_added : int
        The number of beads the run added.

    Returns
    -------
    tuple of (ndarray, tuple of int, bool)
        (n_frames, n_beads, 3) lattice sites with the nucleosome first, the
        run's box in lattice units, and whether the run used hard walls.

    Raises
    ------
    SystemExit
        If there are no production frames, the bead count is wrong, a frozen
        nucleosome bead is not where the restart file put it in every frame
        (PIMMS centres the restart box in a larger run box), or none of the
        nucleosome beads the freeze file leaves free ever moves (a run made with
        a freeze file that also froze them).
    """
    trajectory = lemonade.load(xtc=str(run / TRAJECTORY), pdb=str(run / TOPOLOGY),  # type: ignore[no-untyped-call]
                               keyfile=str(run / KEYFILE))
    store = trajectory.store
    positions = np.asarray(store.positions, dtype=np.int64)[1:]
    if len(positions) == 0:
        raise SystemExit('%s holds no production frames' % (run / TRAJECTORY))

    sites, restart_box = restart_sites(restart)
    centring = (np.array(store.dimensions, dtype=np.int64) - restart_box) // 2
    n_nucleosome = len(sites)
    if positions.shape[1] != n_nucleosome + n_added:
        raise SystemExit('%s has %d beads, not %d (nucleosome) + %d (added)'
                         % (run / TRAJECTORY, positions.shape[1], n_nucleosome, n_added))
    if len(frozen) != n_nucleosome:
        raise SystemExit('the freeze file describes %d beads, but %s holds %d' % (len(frozen), restart, n_nucleosome))
    start = sites + centring
    if not np.all(positions[:, :n_nucleosome][:, frozen] == start[frozen]):
        raise SystemExit('the frozen beads of the nucleosome in %s are not where %s put them; was the run made '
                         'with this freeze file?' % (run / TRAJECTORY, restart))
    # the converse mismatch: a run made with the tails frozen too passes the test
    # above, and its tails would then be reported where the crystal has them
    free = ~frozen
    if np.any(free) and np.all(positions[:, :n_nucleosome][:, free] == start[free]):
        raise SystemExit('none of the %d nucleosome beads %s leaves free ever moves in %s; was the run made '
                         'with this freeze file?' % (int(np.sum(free)), FREEZE_FILE, run / TRAJECTORY))
    return positions, tuple(int(d) for d in store.dimensions), bool(store.hardwall)


# .............................................................................
# Counting contacts
# .............................................................................

def neighbour_values(grid: IntArray, sites: IntArray, hardwall: bool) -> IntArray:
    """
    The value the grid holds on each Chebyshev neighbour of each site.

    Parameters
    ----------
    grid : ndarray of int
        The box, holding a value (a bead index or a protein index) on each
        occupied site and EMPTY_SITE elsewhere.
    sites : ndarray of int
        (n, 3) lattice sites.
    hardwall : bool
        If True, neighbours outside the box are empty; otherwise the box is periodic.

    Returns
    -------
    ndarray of int
        (n, 26) values, EMPTY_SITE where a neighbour is empty.
    """
    dimensions = np.array(grid.shape, dtype=np.int64)
    neighbours = sites[:, None, :] + NEIGHBOUR_OFFSETS[None, :, :]
    if not hardwall:
        neighbours %= dimensions
        wrapped: IntArray = grid[neighbours[..., 0], neighbours[..., 1], neighbours[..., 2]]
        return wrapped
    inside = np.all((neighbours >= 0) & (neighbours < dimensions), axis=2)
    clipped = np.clip(neighbours, 0, dimensions - 1)
    values: IntArray = grid[clipped[..., 0], clipped[..., 1], clipped[..., 2]]
    return np.where(inside, values, EMPTY_SITE)


def touches_partner(sites: IntArray, protein_of_bead: IntArray, dimensions: tuple[int, ...],
                    hardwall: bool) -> BoolArray:
    """
    Which added beads touch a bead of a different added protein.

    Parameters
    ----------
    sites : ndarray of int
        (n, 3) sites of the added beads in one frame.
    protein_of_bead : ndarray of int
        (n,) which protein each bead belongs to.
    dimensions : tuple of int
        The box.
    hardwall : bool
        Whether the box has hard walls.

    Returns
    -------
    ndarray of bool
        (n,) True where the bead touches another protein's bead.
    """
    grid = np.full(dimensions, EMPTY_SITE, dtype=np.int64)
    grid[sites[:, 0], sites[:, 1], sites[:, 2]] = protein_of_bead
    neighbours = neighbour_values(grid, sites, hardwall)
    return np.any((neighbours != EMPTY_SITE) & (neighbours != protein_of_bead[:, None]), axis=1)


def count_contacts(run: Path, restart: Path, frozen: BoolArray,
                   dna: BoolArray) -> tuple[list[ProteinContacts], FloatArray]:
    """
    Count what each added protein touches over a run, and how often each nucleosome bead touches the DNA.

    Parameters
    ----------
    run : Path
        The run directory.
    restart : Path
        The restart file the run started from.
    frozen : ndarray of bool
        One per nucleosome bead: True where the freeze file freezes it.
    dna : ndarray of bool
        One per nucleosome bead: True for a DNA nucleotide.

    Returns
    -------
    tuple of (list of ProteinContacts, ndarray of float)
        One ProteinContacts per added protein, in keyfile order, and per
        nucleosome bead the fraction of frames in which it touches a DNA bead
        (other than itself). Only the free beads' fractions change from frame
        to frame.

    Raises
    ------
    SystemExit
        If the run's files are missing or do not match the restart file.
    """
    for name in (KEYFILE, TRAJECTORY, TOPOLOGY):
        if not (run / name).is_file():
            raise SystemExit('%s is missing; run PIMMS in %s first (see README.md)' % (run / name, run))
    proteins = added_proteins(run)
    lengths = [protein.copies * len(protein.sequence) for protein in proteins]
    positions, dimensions, hardwall = production_frames(run, restart, frozen, sum(lengths))
    n_nucleosome = positions.shape[1] - sum(lengths)
    protein_of_bead = np.repeat(np.arange(len(proteins)), lengths)
    starts = np.concatenate([[0], np.cumsum(lengths)])

    # The frozen beads sit on the same sites in every frame, so one grid of
    # their bead indices serves every frame. The free beads (the histone tails)
    # are written onto it frame by frame and wiped off again; they can never
    # share a site with a frozen bead, so wiping them leaves the frozen beads.
    grid = np.full(dimensions, EMPTY_SITE, dtype=np.int64)
    frozen_beads = np.flatnonzero(frozen)
    free_beads = np.flatnonzero(~frozen)
    fixed = positions[0, frozen_beads]
    grid[fixed[:, 0], fixed[:, 1], fixed[:, 2]] = frozen_beads
    # a lookup from a grid value to "is a DNA bead", with EMPTY_SITE (-1) landing on the trailing False
    is_dna = np.append(dna, False)

    touched = np.zeros((len(proteins), n_nucleosome), dtype=np.int64)
    pairs = np.zeros((len(proteins), n_nucleosome), dtype=np.int64)
    on_dna = np.zeros(n_nucleosome, dtype=np.int64)
    touching = [np.zeros((p.copies, len(p.sequence)), dtype=np.int64) for p in proteins]
    # per protein, per copy: frames on [neither, nucleosome only, partner only, both]
    state_frames = [np.zeros((p.copies, 4), dtype=np.int64) for p in proteins]
    for frame in positions:
        tails = frame[free_beads]
        grid[tails[:, 0], tails[:, 1], tails[:, 2]] = free_beads
        on_dna += np.any(is_dna[neighbour_values(grid, frame[:n_nucleosome], hardwall)], axis=1)
        added = frame[n_nucleosome:]
        beads = neighbour_values(grid, added, hardwall)
        grid[tails[:, 0], tails[:, 1], tails[:, 2]] = EMPTY_SITE
        on_nucleosome = np.any(beads != EMPTY_SITE, axis=1)
        on_partner = (touches_partner(added, protein_of_bead, dimensions, hardwall) if len(proteins) > 1
                      else np.zeros(len(added), dtype=bool))
        for index, protein in enumerate(proteins):
            block = slice(starts[index], starts[index + 1])
            hits = beads[block][beads[block] != EMPTY_SITE]
            touched[index, np.unique(hits)] += 1
            pairs[index] += np.bincount(hits, minlength=n_nucleosome)
            per_bead = on_nucleosome[block].reshape(protein.copies, -1)
            touching[index] += per_bead
            copy_on_nucleosome = np.any(per_bead, axis=1)
            copy_on_partner = np.any(on_partner[block].reshape(protein.copies, -1), axis=1)
            state = copy_on_nucleosome.astype(np.int64) + 2 * copy_on_partner.astype(np.int64)
            state_frames[index][np.arange(protein.copies), state] += 1

    n_frames = len(positions)
    results = []
    for index, protein in enumerate(proteins):
        fractions = state_frames[index].sum(axis=0) / (n_frames * protein.copies)
        results.append(ProteinContacts(
            run=run.name, protein=protein, partners=[p.label for p in proteins if p is not protein],
            n_frames=n_frames, nucleosome_touched=touched[index] / n_frames, nucleosome_pairs=pairs[index] / n_frames,
            residue_touching=touching[index].sum(axis=0) / (n_frames * protein.copies),
            bound=state_frames[index][:, [1, 3]].sum(axis=1) / n_frames,
            states=dict(zip(('free', 'nucleosome', 'partner', 'both'), (float(f) for f in fractions)))))
    return results, on_dna / n_frames


# .............................................................................
# Writing the results
# .............................................................................

def write_nucleosome_csv(path: Path, residues: list[Residue], contacts: ProteinContacts) -> None:
    """
    Write how often each nucleosome residue is touched.

    Parameters
    ----------
    path : Path
        The CSV file.
    residues : list of Residue
        The nucleosome residues in bead order.
    contacts : ProteinContacts
        The protein's contacts.

    Returns
    -------
    None
    """
    lines = ['chain,residue,resname,bead,contact_probability,mean_beads_touching']
    for residue, touched, pairs in zip(residues, contacts.nucleosome_touched, contacts.nucleosome_pairs):
        lines.append('%s,%s%s,%s,%s,%.4f,%.4f' % (residue.chain, residue.resseq, residue.icode, residue.resname,
                                                residue.code, touched, pairs))
    path.write_text('\n'.join(lines) + '\n')


def write_protein_csv(path: Path, contacts: ProteinContacts) -> None:
    """
    Write how often each residue of the added protein touches the nucleosome.

    Parameters
    ----------
    path : Path
        The CSV file.
    contacts : ProteinContacts
        The protein's contacts.

    Returns
    -------
    None
    """
    lines = ['residue,amino_acid,contact_probability']
    for index, (amino_acid, probability) in enumerate(zip(contacts.protein.sequence, contacts.residue_touching)):
        lines.append('%d,%s,%.4f' % (contacts.protein.first_residue + index, amino_acid, probability))
    path.write_text('\n'.join(lines) + '\n')


def _pdb_atom_name(name: str, element: str) -> str:
    """
    Place an atom name in the four columns the PDB format gives it.

    A name starts in column 14 unless it fills all four columns or its element
    has two letters (the PDB format's alignment rule).

    Parameters
    ----------
    name : str
        The atom name.
    element : str
        The element symbol.

    Returns
    -------
    str
        The four-character field.
    """
    if len(name) >= 4 or len(element) == 2:
        return '%-4s' % name[:4]
    return ' %-3s' % name


def write_bfactor_pdb(path: Path, residues: list[Residue], probabilities: FloatArray) -> None:
    """
    Write the structure's heavy atoms with each residue's contact probability as its B-factor.

    Parameters
    ----------
    path : Path
        The PDB file.
    residues : list of Residue
        The nucleosome residues in bead order.
    probabilities : ndarray of float
        One per residue.

    Returns
    -------
    None
    """
    lines = []
    serial = 0
    previous_chain: Optional[str] = None
    for residue, probability in zip(residues, probabilities):
        if previous_chain is not None and residue.chain != previous_chain:
            lines.append('TER')
        previous_chain = residue.chain
        for (name, xyz), element in zip(residue.atoms.items(), residue.heavy_elements):
            serial += 1
            lines.append('%-6s%5d %s %3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s'
                         % ('ATOM', serial, _pdb_atom_name(name, element), residue.resname, residue.chain,
                            residue.resseq, residue.icode or ' ', xyz[0], xyz[1], xyz[2], 1.0, probability,
                            element))
    lines += ['TER', 'END']
    path.write_text('\n'.join(lines) + '\n')


def summarise(contacts: ProteinContacts, residues: list[Residue]) -> None:
    """
    Print how often the protein binds, and where.

    Parameters
    ----------
    contacts : ProteinContacts
        The protein's contacts.
    residues : list of Residue
        The nucleosome residues in bead order.

    Returns
    -------
    None
    """
    protein = contacts.protein
    print('\n%s in %s: %d copies of %d residues, %d production frames'
          % (contacts.label, contacts.run, protein.copies, len(protein.sequence), contacts.n_frames))
    print('  fraction of frames each copy touches the nucleosome: %s'
          % ' '.join('%.2f' % value for value in contacts.bound))
    if contacts.partners:
        print('  copies free %.2f, on the nucleosome only %.2f, on %s only %.2f, on both %.2f'
              % (contacts.states['free'], contacts.states['nucleosome'], ' or '.join(contacts.partners),
                 contacts.states['partner'], contacts.states['both']))
    print('  residues of a copy touching the nucleosome, on average: %.1f' % contacts.residue_touching.sum())
    total_pairs = contacts.nucleosome_pairs.sum()
    if total_pairs == 0:
        return
    chains = np.array([residue.chain for residue in residues])
    groups = [('DNA', DNA_CHAINS)] + list(HISTONE_CHAINS)
    shares = ['%s %.0f%%' % (name, 100 * contacts.nucleosome_pairs[np.isin(chains, members)].sum() / total_pairs)
              for name, members in groups]
    print('  share of the contacts with the nucleosome: %s' % ', '.join(shares))
    top = np.argsort(contacts.nucleosome_touched)[::-1][:TOP_RESIDUES]
    print('  most often touched: %s'
          % ', '.join('%s:%s%s %.2f' % (residues[i].chain, residues[i].resname, residues[i].resseq,
                                        contacts.nucleosome_touched[i]) for i in top))


def free_tails(residues: list[Residue], frozen: BoolArray) -> dict[str, list[tuple[str, IntArray]]]:
    """
    The free stretches of the nucleosome's chains (its tails), grouped by histone and end.

    Parameters
    ----------
    residues : list of Residue
        The nucleosome residues in bead order.
    frozen : ndarray of bool
        One per nucleosome bead: True where the freeze file freezes it.

    Returns
    -------
    dict of str to list of (str, ndarray of int)
        e.g. ``'H3 N-terminal'`` -> one (``'A:1-26'``, bead indices) pair per
        copy, in bead order. A free stretch that reaches neither end of its
        chain is called internal.
    """
    histone_of_chain = {chain: name for name, members in HISTONE_CHAINS for chain in members}
    chains = np.array([residue.chain for residue in residues])
    tails: dict[str, list[tuple[str, IntArray]]] = {}
    # one chain at a time, so that the end of one chain's tail and the start of
    # the next chain's are never joined into one stretch
    for chain in dict.fromkeys(chains):
        chain_beads = np.flatnonzero(chains == chain)
        for first, last in contiguous_runs(np.flatnonzero(~frozen[chain_beads])):
            if first == 0:
                end = 'N-terminal'
            elif last == len(chain_beads) - 1:
                end = 'C-terminal'
            else:
                end = 'internal'
            beads = chain_beads[first:last + 1]
            label = '%s:%s-%s' % (chain, residues[beads[0]].resseq, residues[beads[-1]].resseq)
            tails.setdefault('%s %s' % (histone_of_chain.get(chain, 'chain ' + chain), end), []).append((label, beads))
    return tails


def summarise_tails(run: str, on_dna: FloatArray, residues: list[Residue], frozen: BoolArray) -> None:
    """
    Print how much of each free histone tail lies on the DNA.

    Parameters
    ----------
    run : str
        The run directory's name.
    on_dna : ndarray of float
        Per nucleosome bead: the fraction of frames in which it touches a DNA bead.
    residues : list of Residue
        The nucleosome residues in bead order.
    frozen : ndarray of bool
        One per nucleosome bead: True where the freeze file freezes it.

    Returns
    -------
    None
    """
    tails = free_tails(residues, frozen)
    if not tails:
        return
    print('\nHistone tails in %s: the fraction of their residues touching the DNA, per copy' % run)
    for name, copies in tails.items():
        print('  %-16s %s' % (name, ', '.join('%s %.2f' % (label, on_dna[beads].mean()) for label, beads in copies)))
    every_tail = np.concatenate([beads for copies in tails.values() for _, beads in copies])
    print('  %-16s %.2f' % ('every tail', on_dna[every_tail].mean()))


# .............................................................................
# The figure
# .............................................................................

def _per_chain(residues: list[Residue], values: FloatArray, members: Sequence[str]) -> list[dict[int, float]]:
    """
    Residue number -> value, for each of the named chains.

    Parameters
    ----------
    residues : list of Residue
        The nucleosome residues in bead order.
    values : ndarray of float
        One per residue.
    members : sequence of str
        The chains.

    Returns
    -------
    list of dict of int to float
        One mapping per chain.
    """
    return [{int(residue.resseq): float(value) for residue, value in zip(residues, values)
             if residue.chain == chain} for chain in members]


def _style_axes(axes: Axes, xlabel: str, ylabel: Optional[str]) -> None:
    """
    Label an axes and strip it of everything else.

    Parameters
    ----------
    axes : matplotlib Axes
        The axes.
    xlabel : str
        The x label.
    ylabel : str or None
        The y label, or None for none.

    Returns
    -------
    None
    """
    axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.set_ylim(bottom=0)


def plot(entries: list[ProteinContacts], residues: list[Residue], path: Path) -> None:
    """
    Plot where each added protein touches the nucleosome, and with which of its residues.

    The top row is the DNA (each base pair's two nucleotides averaged, against
    its position from the dyad), the middle row the histones (the two copies of
    each averaged, on one scale) and the bottom row each added protein.

    Parameters
    ----------
    entries : list of ProteinContacts
        Every added protein of every run, in plotting order.
    residues : list of Residue
        The nucleosome residues in bead order.
    path : Path
        The PDF to write.

    Returns
    -------
    None

    Raises
    ------
    SystemExit
        If there are more entries than colours.
    """
    if len(entries) > len(COLOURS):
        raise SystemExit('the figure has %d colours, for %d proteins' % (len(COLOURS), len(entries)))
    plt.rcParams.update({'font.family': FONT, 'font.size': FONT_SIZE, 'pdf.fonttype': EDITABLE_FONT_TYPE,
                         'ps.fonttype': EDITABLE_FONT_TYPE, 'axes.linewidth': LINE_WIDTH,
                         'lines.linewidth': LINE_WIDTH})
    proteins = list(dict.fromkeys(entry.protein.label for entry in entries))
    figure = plt.figure(figsize=FIGURE_SIZE_INCHES, constrained_layout=True)
    grid = figure.add_gridspec(3, max(len(HISTONE_CHAINS), len(proteins)))
    probability = 'Contact probability'

    dna_axes = figure.add_subplot(grid[0, :])
    for entry, colour in zip(entries, COLOURS):
        strand_i, strand_j = _per_chain(residues, entry.nucleosome_touched, DNA_CHAINS)
        positions = sorted(strand_i)
        dna_axes.plot(positions, [(strand_i[n] + strand_j[-n]) / 2 for n in positions], color=colour,
                      label=entry.label)
    _style_axes(dna_axes, 'Base pair (from the dyad)', probability)
    # above the panel, where it cannot hide a line
    dna_axes.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=len(entries))

    histone_axes: list[Axes] = []
    for column, (histone, members) in enumerate(HISTONE_CHAINS):
        axes = figure.add_subplot(grid[1, column], sharey=histone_axes[0] if histone_axes else None)
        for entry, colour in zip(entries, COLOURS):
            copies = _per_chain(residues, entry.nucleosome_touched, members)
            numbers = sorted(copies[0])
            axes.plot(numbers, [np.mean([copy[n] for copy in copies]) for n in numbers], color=colour)
        _style_axes(axes, '%s residue' % histone, probability if column == 0 else None)
        histone_axes.append(axes)

    for column, label in enumerate(proteins):
        axes = figure.add_subplot(grid[2, column])
        for entry, colour in zip(entries, COLOURS):
            if entry.protein.label == label:
                residue_numbers = np.arange(len(entry.protein.sequence)) + entry.protein.first_residue
                axes.plot(residue_numbers, entry.residue_touching, color=colour)
        _style_axes(axes, '%s residue' % label, probability if column == 0 else None)

    figure.savefig(path)
    plt.close(figure)


# .............................................................................
# Command line
# .............................................................................

def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Analyse the runs and write the results.

    Parameters
    ----------
    argv : sequence of str or None
        The arguments; None means sys.argv[1:].

    Returns
    -------
    int
        0 on success.
    """
    parser = argparse.ArgumentParser(description='Where the added disordered proteins touch the nucleosome.')
    parser.add_argument('runs', nargs='*', default=list(RUNS), metavar='RUN_DIRECTORY',
                        help='Run directories to analyse (default: %s)' % ' '.join(RUNS))
    args = parser.parse_args(argv)

    residues = nucleosome_residues(STRUCTURE_FILE, RESTART_FILE)
    frozen = frozen_nucleosome_beads(FREEZE_FILE, RESTART_FILE)
    dna = np.isin([residue.chain for residue in residues], DNA_CHAINS)
    entries = []
    for name in args.runs:
        run = Path(name)
        proteins, on_dna = count_contacts(run, RESTART_FILE, frozen, dna)
        for contacts in proteins:
            file_name = contacts.protein.file_name
            write_nucleosome_csv(run / (NUCLEOSOME_CSV % file_name), residues, contacts)
            write_protein_csv(run / (PROTEIN_CSV % file_name), contacts)
            write_bfactor_pdb(run / (BFACTOR_PDB % file_name), residues, contacts.nucleosome_touched)
            summarise(contacts, residues)
            entries.append(contacts)
        summarise_tails(run.name, on_dna, residues, frozen)
    plot(entries, residues, FIGURE_FILE)
    print('\nWrote %s, and the per-protein contact tables and B-factor structures in each run directory'
          % FIGURE_FILE)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
