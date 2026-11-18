## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

##
## pdb2pimms
##
## Turns a PDB or mmCIF structure into a PIMMS restart file, plus a Markdown
## report that says exactly what the conversion did. Every residue becomes one
## bead: amino acids keep their one-letter code, nucleotides get lower-case
## codes, and the report lists every bead type the parameter file of the run
## must define. The beads are then fitted onto the lattice. A polymer chain
## has to become a self-avoiding walk of Moore-neighbour (Chebyshev-1) steps,
## so each chain is placed by a dynamic programme that finds the walk closest
## to the structure (least squares over every bead), then repaired where the
## walk reuses a site, with a greedy depth-first search held in reserve.
##
## The restart file is a starting configuration (RESTART_FILE), written with
## HARDWALL True so that a run may use either boundary condition and any box
## at least as large as the one written. It is read back with PIMMS's own
## restart reader before it is put in place, so a file this tool writes is one
## PIMMS accepts. Unless told not to, it also writes a freeze file (FREEZE_FILE)
## listing every chain, so the structure can be held fixed while other chains
## move around it; that file is read back with PIMMS's own reader too. Residues
## named with --mobile (a disordered tail, a loop) are left out of the freeze
## file, which then freezes the rest of their chain bead by bead.
##
## This module is the whole of the tool; scripts/pdb2pimms only calls main().
##

import argparse
import datetime
import hashlib
import itertools
import math
import os
import pickle
import re
import shlex
import sys
import textwrap
import unicodedata
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence, cast

import numpy as np
import numpy.typing as npt

from . import CONFIG
from .data_structures import contiguous_runs, format_bead_runs
from .latticeExceptions import (KeyFileException, RestartException, StructureConversionException,
                                StructureFileException)
from .structure_reader import AtomRecord, StructureFile, read_structure


FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]
BoolArray = npt.NDArray[np.bool_]
Site = tuple[int, int, int]

# .............................................................................
# Defaults and limits
# .............................................................................

# The lattice spacing used when none is given. It is PIMMS's own default
# LATTICE_TO_ANGSTROMS, so the START.pdb and traj.xtc of a run started from the
# restart file come out on the same scale as the input structure.
# (CONFIG.DEFAULTS holds values of every type, hence the cast.)
DEFAULT_SPACING_ANGSTROMS = float(cast(float, CONFIG.DEFAULTS['LATTICE_TO_ANGSTROMS']))

# Empty sites left between the structure and each face of the box. Three is the
# reach of the longest-range interaction (super-long-range, Chebyshev distance
# 3), so even in a periodic box of exactly this size no bead can interact with
# a periodic image of the structure.
DEFAULT_PADDING_SITES = 3

# PIMMS refuses a box axis shorter than this (the super-long-range shell needs it)
MIN_BOX_SITES = 7

# START.pdb holds each coordinate in an 8-column field that ends at 9999.999 A,
# so PIMMS refuses a box whose longest axis times LATTICE_TO_ANGSTROMS reaches this
MAX_BOX_EDGE_ANGSTROMS = 10000.0

# The two ways a bead's position can be chosen (see representative_position)
BEAD_POSITION_BACKBONE = 'backbone'
BEAD_POSITION_CENTROID = 'centroid'

# A placeholder for the restart file's ENERGY. A plain restart never reads it
# (every run recomputes its energy from scratch), and RESTART_CONTINUE, the only
# reader of it, refuses this file anyway because it holds no continuation state.
PLACEHOLDER_ENERGY = 0.0

# How the report (and --chains) writes a blank chain identifier
BLANK_CHAIN_LABEL = '-'

# Bead characters PIMMS cannot use: 0 is the solvent type and # starts a
# comment in the parameter file. Whitespace and unprintable characters are
# refused separately.
RESERVED_BEAD_CHARACTERS = ('0', '#')

# .............................................................................
# Residue kinds and bead codes
# .............................................................................

KIND_AMINO_ACID = 'amino acid'
KIND_RNA = 'RNA nucleotide'
KIND_DNA = 'DNA nucleotide'
KIND_HETERO = 'hetero group'
KIND_WATER = 'water'
POLYMER_KINDS = (KIND_AMINO_ACID, KIND_RNA, KIND_DNA)

# The 20 standard amino acids, plus selenocysteine and pyrrolysine, by their
# one-letter codes (IUPAC-IUB Joint Commission on Biochemical Nomenclature)
AMINO_ACID_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
    'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O',
}

# Nucleotides get lower-case codes so that they never share a bead type with
# an amino acid (A, C, G, T and U are all amino-acid codes too). RNA and DNA
# share a, c and g; a structure holding both gets a warning saying so. The
# standard names are those of the PDB Chemical Component Dictionary (CCD): A,
# C, G and U for RNA, and DA, DC, DG and DT for DNA.
RNA_BASES = {'A': 'a', 'C': 'c', 'G': 'g', 'U': 'u'}
DNA_BASES = {'A': 'a', 'C': 'c', 'G': 'g', 'T': 't'}
DNA_PREFIX = 'D'

# The CCD names of an unknown amino acid, ribonucleotide and deoxyribonucleotide
UNKNOWN_AMINO_ACID_NAME = 'UNK'
UNKNOWN_RNA_NAME = 'N'
UNKNOWN_DNA_NAME = 'DN'

# The code of an unknown amino acid (UNK, or a residue of unknown name with an
# amino-acid backbone) and of an unknown nucleotide (IUPAC N = any nucleotide)
UNKNOWN_AMINO_ACID_CODE = 'X'
UNKNOWN_NUCLEOTIDE_CODE = 'n'

# Names that simulation packages (AMBER, CHARMM, GROMACS) give to protonation
# states of the standard amino acids. Several are also CCD codes for unrelated
# molecules (the CCD's HSD, CYX, ASH and HID are ligands), so these names, like
# every simulation-package name below, are trusted only for a residue that has
# the backbone of the residue they claim (see _has_amino_backbone).
AMINO_ACID_VARIANTS = {
    'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
    'HISD': 'HIS', 'HISE': 'HIS', 'HISH': 'HIS', 'HISA': 'HIS', 'HISB': 'HIS',
    'CYX': 'CYS', 'CYM': 'CYS', 'CYS2': 'CYS', 'CYSH': 'CYS',
    'ASH': 'ASP', 'ASPH': 'ASP', 'ASPP': 'ASP',
    'GLH': 'GLU', 'GLUH': 'GLU', 'GLUP': 'GLU',
    'LYN': 'LYS', 'LYSN': 'LYS', 'LYSH': 'LYS',
}

# AMBER names a terminal amino acid by prefixing N or C to its three-letter name
# (NALA, CLYS), and a terminal nucleotide by suffixing 5 or 3 (DA5, RA3, A5)
TERMINAL_PREFIXES = ('N', 'C')
TERMINAL_SUFFIXES = ('5', '3')

# The protonation states (AMINO_ACID_VARIANTS names) that AMBER also has N- and
# C-terminal residues of, which it names the same way: NHID, NHIE, NHIP, NCYX,
# CHID, CHIE, CHIP, CCYX. These are the terminal residue templates of ff99SB-ILDN,
# ff14SB and ff19SB (from AmberTools' leaprc.protein files, as converted in
# OpenMM 8.4's amber99sbildn.xml, amber14/protein.ff14SB.xml and
# amber19/protein.ff19SB.xml). None of them has a terminal form of ASH, GLH,
# LYN or CYM.
AMBER_TERMINAL_VARIANTS = ('HID', 'HIE', 'HIP', 'CYX')

# AMBER's older RNA names put an R before the base (RA, RU)
RNA_PREFIX = 'R'

# Version 2 of the PDB format named thymidine T
LEGACY_THYMIDINE = 'T'

# CHARMM names nucleotides by base alone, for DNA and RNA alike, so the sugar
# decides which (a 2'-hydroxyl means RNA)
BASE_ONLY_NUCLEOTIDES = {'ADE': 'A', 'CYT': 'C', 'GUA': 'G', 'THY': 'T', 'URA': 'U'}

# Modified residues mapped to the residue they derive from when the file itself
# does not say (PDB MODRES / mmCIF _pdbx_struct_mod_residue take precedence).
# These are CCD codes, and the parents are the CCD's.
FALLBACK_PARENTS = {
    # amino acids
    'MSE': 'MET', 'SEP': 'SER', 'TPO': 'THR', 'PTR': 'TYR', 'HYP': 'PRO', 'MLY': 'LYS',
    'M3L': 'LYS', 'ALY': 'LYS', 'KCX': 'LYS', 'CSO': 'CYS', 'CME': 'CYS', 'OCS': 'CYS',
    'CGU': 'GLU',
    # nucleotides
    'PSU': 'U', 'H2U': 'U', '5MU': 'U', '4SU': 'U', 'OMU': 'U', '5MC': 'C', 'OMC': 'C',
    'OMG': 'G', '1MG': 'G', '2MG': 'G', 'M2G': 'G', '7MG': 'G', 'YYG': 'G', '1MA': 'A',
}

# Codes handed out, in this order, to hetero groups kept with --hetero that have
# no --map entry: lower-case letters that no nucleotide code uses
HETERO_CODE_POOL = 'bdefhijklmopqrsvwxyz'

# Water as deposited structures name it (the CCD codes of water and heavy water)
WATER_NAMES = frozenset(('HOH', 'DOD'))

# Water as simulation packages name it. Some of these are also CCD codes of other
# molecules (T3P is thymidine 3'-phosphate, SOL is L-sorbose), so a residue with
# one of these names is water only if it is shaped like water (see _is_water)
SIMULATION_WATER_NAMES = frozenset(('WAT', 'H2O', 'D2O', 'TIP', 'TIP3', 'TIP4', 'TIP5', 'T3P', 'T4P',
                                    'T5P', 'SOL', 'SPC', 'SPCE'))

# Elements a water molecule never holds. Its heavy atoms are one oxygen, plus a
# massless virtual site in four- and five-site water models.
NON_WATER_ELEMENTS = ('C', 'N', 'P', 'S')

# Elements whose atoms are hydrogens (deuterium included)
HYDROGEN_ELEMENTS = ('H', 'D')

# The element a CA-only model's lone CA must have: carbon, or none given (a
# calcium ion is also a lone atom named CA, with element CA)
CARBON_ELEMENTS = ('C', '')

# The residue name of a calcium ion in the CCD and in AMBER's ion parameters
# (OpenMM 8.4's amber14/tip3p.xml), whose one atom is also named CA. It tells the
# ion from a CA-only residue when the file gives no element: the element is then
# guessed from the atom name's first letter, which makes calcium a carbon.
CALCIUM_ION_NAME = 'CA'

# Atom names (PDB version 3, then the asterisk spellings of older files)
CA_ATOM = 'CA'
C4_PRIME_ATOMS = ("C4'", 'C4*')
O3_PRIME_ATOMS = ("O3'", 'O3*')
O5_PRIME_ATOMS = ("O5'", 'O5*')
O2_PRIME_ATOMS = ("O2'", 'O2*')
C2_PRIME_ATOMS = ("C2'", 'C2*')
PHOSPHORUS_ATOM = 'P'
PEPTIDE_N_ATOM = 'N'
PEPTIDE_C_ATOM = 'C'

# The only atoms a residue of a coarse (trace) model has: CA for an amino acid,
# C4' and/or P for a nucleotide
TRACE_ATOMS = (CA_ATOM,) + C4_PRIME_ATOMS + (PHOSPHORUS_ATOM,)

# .............................................................................
# Chain connectivity
# .............................................................................

# Consecutive residues are bonded if the atoms of their linking bond are within
# this distance. The peptide C-N bond is 1.33 A (Engh & Huber, Acta Cryst. A47,
# 392, 1991) and the phosphodiester O3'-P bond 1.61 A (Parkinson et al., Acta
# Cryst. D52, 57, 1996); 2.0 A accepts every real bond and no missing residue.
MAX_LINK_BOND_ANGSTROMS = 2.0

# The test used when the linking atoms are missing (CA-only or C4'-only models).
# Consecutive CA atoms are 3.8 A apart (2.9 A in a cis peptide), and a single
# missing residue puts them further apart than this.
MAX_CA_CA_ANGSTROMS = 4.2

# Consecutive C4' atoms of a nucleic-acid backbone are about 5 to 7 A apart, and
# consecutive phosphorus atoms about 6 to 7 A (the P-only fallback)
MAX_C4_C4_ANGSTROMS = 7.5
MAX_P_P_ANGSTROMS = 7.5

# Human-readable bond lengths quoted in chain-break messages
PEPTIDE_BOND_ANGSTROMS = 1.33
PHOSPHODIESTER_BOND_ANGSTROMS = 1.61

# .............................................................................
# Lattice fitting
# .............................................................................

# The 26 Moore (Chebyshev-1) neighbour offsets, the only steps a PIMMS bond can take
MOORE_OFFSETS: IntArray = np.array(
    [offset for offset in itertools.product((-1, 0, 1), repeat=3) if offset != (0, 0, 0)],
    dtype=np.int64)

# The dynamic programme considers, for each bead, every site within this
# Chebyshev distance of the site nearest its target. Two is ample: a fit worse
# than a couple of sites is a sign the spacing is wrong. A wider radius is tried
# only when no walk exists inside the first.
DYNAMIC_PROGRAMME_RADII = (2, 3)

# How many rounds of repair a chain's walk may have for reusing its own sites
# (each round forbids every reused site to all but one of its beads). A chain
# that still clashes after this many is packed more densely than the lattice can
# hold near its targets, and the greedy search places it instead.
MAX_CLASH_REPAIR_ROUNDS = 40

# Sites are packed into one integer for fast membership tests: each coordinate
# is shifted by SITE_KEY_OFFSET and given SITE_KEY_BITS bits
SITE_KEY_BITS = 21
SITE_KEY_OFFSET = 1 << (SITE_KEY_BITS - 1)

# The greedy search places the first bead of a chain at the free site nearest
# its target, looking out to this Chebyshev distance
GREEDY_START_SEARCH_RADIUS = 30

# The greedy search gives up after this many placements per bead (plus a constant),
# which bounds its backtracking
GREEDY_PLACEMENTS_PER_BEAD = 50
GREEDY_PLACEMENTS_MINIMUM = 10000

# A fit whose RMSD exceeds this many lattice spacings gets a warning: 1.5 times
# the half spacing that rounding to a cubic lattice costs on its own
POOR_FIT_RMSD_SPACINGS = 0.75

# Bytes of memory PIMMS uses per lattice site for its two integer grids
LATTICE_BYTES_PER_SITE = 8

# The report names at most this many residues in any one list
REPORT_LIST_LIMIT = 30

# The default restart file adds RESTART_SUFFIX to the input file's stem, and the
# report and the freeze file replace the restart file's ending (either part of
# it) with REPORT_SUFFIX and FREEZE_SUFFIX
RESTART_EXTENSION = '.pimms'
RESTART_STEM_SUFFIX = '_restart'
RESTART_SUFFIX = RESTART_STEM_SUFFIX + RESTART_EXTENSION
REPORT_SUFFIX = '_report.md'
FREEZE_SUFFIX = '_freeze.txt'

# Input extensions removed to form that stem
STRUCTURE_EXTENSIONS = ('.gz', '.pdb', '.ent', '.cif', '.mmcif')

# The residue numbers of a --mobile selection: one author residue number or an
# inclusive range of them, either of which may be negative (DNA strands are often
# numbered through zero), e.g. 5, 1-37 or -73--60
MOBILE_RESIDUE_RANGE = re.compile(r'(-?[0-9]+)(?:-(-?[0-9]+))?', re.ASCII)
# an author residue number as a file writes it
RESIDUE_NUMBER = re.compile(r'-?[0-9]+', re.ASCII)

# A PDB-format residue number above 9999 in the hybrid-36 encoding (from cctbx,
# and written by programs built on it) that fits it into the four-column field:
# A000 to ZZZZ count on from 10000, then a000 to zzzz (Grosse-Kunstleve, Adams
# and Moriarty, https://cci.lbl.gov/hybrid_36/). The offsets are that scheme's
# for a four-character field.
HYBRID36_RESIDUE_NUMBER = re.compile(r'[A-Z][0-9A-Z]{3}|[a-z][0-9a-z]{3}', re.ASCII)
HYBRID36_BASE = 36
HYBRID36_FIRST_VALUE = 10 ** 4
HYBRID36_UPPER_OFFSET = -10 * HYBRID36_BASE ** 3 + HYBRID36_FIRST_VALUE
HYBRID36_LOWER_OFFSET = 16 * HYBRID36_BASE ** 3 + HYBRID36_FIRST_VALUE

# The width the --mobile selections in the header comment of a freeze file are
# wrapped to, so that a long list of them stays readable (the header's other
# lines are at most this long)
FREEZE_HEADER_WIDTH = 88


# .............................................................................
# Data structures
# .............................................................................

@dataclass
class Residue:
    """
    One residue of the converted model, and what the conversion decided about it.

    Attributes
    ----------
    chain, resseq, icode, resname : str
        The residue's identity as the file gives it.
    record : str
        'ATOM' or 'HETATM' (of its first atom).
    atoms : dict of str to ndarray
        Heavy (non-hydrogen) atom name -> coordinates (A), after alternate
        locations have been resolved.
    hydrogens : list of ndarray
        Coordinates of the hydrogen atoms, which only a residue with no heavy
        atom uses.
    heavy_elements : list of str
        The element of each heavy atom (``''`` where the file gives none).
    altloc : str
        The alternate location kept, ``''`` if the residue has none.
    altlocs_present : tuple of str
        Every alternate location the file gives for the residue.
    kind : str
        One of the KIND_* constants.
    parent : str
        The standard residue whose code was used, when that is not the residue itself.
    parent_source : str
        Where ``parent`` came from: 'file', 'built-in table', 'other name' (a
        simulation package's or an older PDB name for a standard residue) or
        'unknown'; ``''`` for a standard residue.
    code : str
        The bead character.
    position : ndarray or None
        The bead's position (A).
    position_source : str
        Which atom or average gave ``position``.
    split_by_repeated_name : bool
        True if the residue was started by an atom that repeats an atom name
        of the residue just before it, which has the same chain, number and
        insertion code (see build_residues).
    """
    chain: str
    resseq: str
    icode: str
    resname: str
    record: str
    atoms: dict[str, FloatArray]
    hydrogens: list[FloatArray]
    heavy_elements: list[str] = field(default_factory=list)
    altloc: str = ''
    altlocs_present: tuple[str, ...] = ()
    kind: str = ''
    parent: str = ''
    parent_source: str = ''
    code: str = ''
    position: Optional[FloatArray] = None
    position_source: str = ''
    split_by_repeated_name: bool = False

    @property
    def label(self) -> str:
        """
        Identify the residue the way the report names it, e.g. ``A:LYS48``.

        Returns
        -------
        str
            Chain, residue name, number and any insertion code.
        """
        return '%s:%s%s%s' % (self.chain or BLANK_CHAIN_LABEL, self.resname, self.resseq, self.icode)


@dataclass
class ChainBreak:
    """
    A place where consecutive residues of one chain were not bonded.

    Attributes
    ----------
    chain : str
        The chain identifier.
    before, after : str
        Labels of the residues either side of the break.
    reason : str
        Why they were not bonded.
    """
    chain: str
    before: str
    after: str
    reason: str


@dataclass
class Segment:
    """
    One PIMMS chain: a run of bonded residues (or a single hetero group).

    Attributes
    ----------
    residues : list of Residue
        The beads, N to C (5' to 3').
    chain : str
        The chain identifier of the source chain.
    """
    residues: list[Residue]
    chain: str

    @property
    def sequence(self) -> str:
        """
        The chain's PIMMS sequence, one bead character per residue.

        Returns
        -------
        str
            The sequence.
        """
        return ''.join(residue.code for residue in self.residues)

    @property
    def span(self) -> str:
        """
        The residues the chain covers, e.g. ``A:MET1-GLY76``.

        Returns
        -------
        str
            First and last residue labels.
        """
        first = self.residues[0]
        last = self.residues[-1]
        if first is last:
            return first.label
        return '%s-%s%s%s' % (first.label, last.resname, last.resseq, last.icode)


def _residue_number(residue: Residue) -> Optional[int]:
    """
    The author residue number of a residue, as an integer.

    A number above 9999 in the hybrid-36 encoding of PDB files (``A000`` is
    10000, ``A001`` 10001, ...) is decoded, so that it compares like the plain
    number an mmCIF file gives.

    Parameters
    ----------
    residue : Residue
        The residue.

    Returns
    -------
    int or None
        The number, or None if the file gives something that is neither an
        integer nor a hybrid-36 number.
    """
    text = residue.resseq.strip()
    if RESIDUE_NUMBER.fullmatch(text) is not None:
        return int(text)
    if HYBRID36_RESIDUE_NUMBER.fullmatch(text) is None:
        return None
    # int() reads base-36 digits in either case
    offset = HYBRID36_UPPER_OFFSET if text[0].isupper() else HYBRID36_LOWER_OFFSET
    return int(text, HYBRID36_BASE) + offset


@dataclass(frozen=True)
class MobileSelection:
    """
    Residues that --mobile leaves free to move.

    Attributes
    ----------
    chain : str
        The author chain identifier (``''`` for a blank one).
    first, last : int or None
        The inclusive range of author residue numbers, or None for the whole
        chain. Insertion codes are ignored, so 52A is inside 50-55.
    text : str
        The selection as given on the command line.
    """
    chain: str
    first: Optional[int]
    last: Optional[int]
    text: str

    def matches(self, residue: Residue) -> bool:
        """
        Whether a residue is in the selection.

        Parameters
        ----------
        residue : Residue
            The residue.

        Returns
        -------
        bool
            True if the residue is on the selected chain and, for a range, its
            number is inside it (a residue whose number is neither an integer
            nor a hybrid-36 number is inside no range).
        """
        if residue.chain != self.chain:
            return False
        if self.first is None or self.last is None:
            return True
        number = _residue_number(residue)
        return number is not None and self.first <= number <= self.last


@dataclass
class ConversionOptions:
    """
    Everything that controls a conversion.

    Attributes
    ----------
    spacing : float
        Lattice spacing in Angstroms per site.
    padding : int
        Empty sites between the structure and each box face (ignored when ``box`` is given).
    box : list of int or None
        An explicit box (three site counts) to centre the structure in.
    model : int or None
        Model number to convert; None means the first in the file.
    chains : list of str or None
        Author chain identifiers to convert; None means every chain.
    bead_position : str
        BEAD_POSITION_BACKBONE or BEAD_POSITION_CENTROID.
    keep_hetero : bool
        If True, every non-water hetero group becomes a one-bead chain.
    bead_map : dict of str to str
        Residue name -> bead character, overriding the built-in codes.
    mobile : list of MobileSelection
        Residues the freeze file leaves free to move (it does not change the
        restart file).
    """
    spacing: float = DEFAULT_SPACING_ANGSTROMS
    padding: int = DEFAULT_PADDING_SITES
    box: Optional[list[int]] = None
    model: Optional[int] = None
    chains: Optional[list[str]] = None
    bead_position: str = BEAD_POSITION_BACKBONE
    keep_hetero: bool = False
    bead_map: dict[str, str] = field(default_factory=dict)
    mobile: list[MobileSelection] = field(default_factory=list)


@dataclass
class ExclusionSummary:
    """
    What the conversion read but did not turn into beads.

    Attributes
    ----------
    waters : int
        Water residues dropped.
    hetero : Counter
        Hetero groups dropped, by residue name (empty with --hetero).
    hydrogens : int
        Hydrogen atoms, which never affect a bead's position (except in a
        residue with no heavy atom).
    altloc_residues : list of Residue
        Residues with alternate locations (one location kept for each).
    other_models : int
        Models in the file other than the one converted.
    unselected_chains : list of str
        Chains left out by --chains.
    """
    waters: int = 0
    hetero: Counter[str] = field(default_factory=Counter)
    hydrogens: int = 0
    altloc_residues: list[Residue] = field(default_factory=list)
    other_models: int = 0
    unselected_chains: list[str] = field(default_factory=list)


@dataclass
class FitStatistics:
    """
    How closely the lattice configuration reproduces the input structure.

    Deviations are measured after the optimal rigid-body superposition of the
    lattice beads onto the input bead positions (Kabsch, Acta Cryst. A32, 922, 1976).

    Attributes
    ----------
    rmsd : float
        Root-mean-square deviation over every bead (A).
    segment_rmsd : list of float
        The same, per PIMMS chain (A).
    max_deviation : float
        The largest single-bead deviation (A).
    max_deviation_label : str
        The residue it belongs to.
    bond_steps : Counter
        Lattice bonds by the number of axes they step along (1, 2 or 3).
    input_bond_lengths, lattice_bond_lengths : ndarray
        Lengths of every bond between bead positions in the input and on the lattice (A).
    input_rg, lattice_rg : list of float
        Radius of gyration of each PIMMS chain in the input and on the lattice (A).
    """
    rmsd: float
    segment_rmsd: list[float]
    max_deviation: float
    max_deviation_label: str
    bond_steps: Counter[int]
    input_bond_lengths: FloatArray
    lattice_bond_lengths: FloatArray
    input_rg: list[float]
    lattice_rg: list[float]


@dataclass
class Conversion:
    """
    The complete result of converting one structure, ready to be written out.

    Attributes
    ----------
    structure : StructureFile
        The file that was read.
    options : ConversionOptions
        The options used.
    model : int
        The model converted.
    residues : list of Residue
        Every residue of that model (after --chains), waters included.
    segments : list of Segment
        The PIMMS chains, in chainID order.
    breaks : list of ChainBreak
        Where source chains were split.
    chain_types : list of int
        The chainType of each segment (identical sequences share one).
    sites : list of list of tuple
        The lattice site of every bead of every segment.
    dimensions : list of int
        The box (DIMENSIONS).
    methods : list of str
        How each segment was placed.
    statistics : FitStatistics
        Fit quality.
    excluded : ExclusionSummary
        What was left out.
    position_fallbacks : list of str
        Labels of residues whose backbone atom was missing (placed at their centroid).
    isolated_hetatm : list of str
        Labels of HETATM residues with a polymer residue's name or backbone but
        bonded to neither neighbour, which were treated as hetero groups.
    warnings : list of str
        Everything the user should look at.
    input_sha256 : str
        SHA-256 digest of the input file.
    """
    structure: StructureFile
    options: ConversionOptions
    model: int
    residues: list[Residue]
    segments: list[Segment]
    breaks: list[ChainBreak]
    chain_types: list[int]
    sites: list[list[Site]]
    dimensions: list[int]
    methods: list[str]
    statistics: FitStatistics
    excluded: ExclusionSummary
    position_fallbacks: list[str]
    isolated_hetatm: list[str]
    warnings: list[str]
    input_sha256: str


# .............................................................................
# Building residues from atoms
# .............................................................................

def _is_hydrogen(atom: AtomRecord) -> bool:
    """
    Say whether an atom is a hydrogen (or deuterium).

    Parameters
    ----------
    atom : AtomRecord
        The atom.

    Returns
    -------
    bool
        True for a hydrogen.
    """
    return atom.element in HYDROGEN_ELEMENTS


def _chosen_altloc(atoms: list[AtomRecord]) -> tuple[str, tuple[str, ...]]:
    """
    Choose which alternate location of a residue to keep.

    The location whose atoms have the highest mean occupancy wins (the mean,
    not the sum, so that a location modelled for fewer atoms is not penalised),
    and the first one the file lists wins a tie (normally ``A``). Atoms with no
    alternate location are kept whichever is chosen.

    Parameters
    ----------
    atoms : list of AtomRecord
        The residue's atoms, in file order.

    Returns
    -------
    tuple of (str, tuple of str)
        The chosen location (``''`` if the residue has none) and every location present.
    """
    total: dict[str, float] = {}
    count: dict[str, int] = {}
    for atom in atoms:
        if atom.altloc:
            total[atom.altloc] = total.get(atom.altloc, 0.0) + atom.occupancy
            count[atom.altloc] = count.get(atom.altloc, 0) + 1
    if not total:
        return '', ()
    present = tuple(total)
    best = max(present, key=lambda location: (total[location] / count[location], -present.index(location)))
    return best, present


def _residue_from_atoms(atoms: list[AtomRecord]) -> Residue:
    """
    Build one residue from its atoms, keeping a single alternate location.

    Parameters
    ----------
    atoms : list of AtomRecord
        The residue's atoms, in file order (non-empty).

    Returns
    -------
    Residue
        The residue, with its kind, code and position still to be decided.
    """
    chosen, present = _chosen_altloc(atoms)
    kept = [atom for atom in atoms if atom.altloc in ('', chosen)]

    # A residue with microheterogeneity has a different name in each alternate
    # location, so the name is taken from the atoms that were kept
    resname = Counter(atom.resname for atom in kept).most_common(1)[0][0]

    heavy: dict[str, FloatArray] = {}
    elements: list[str] = []
    hydrogens: list[FloatArray] = []
    for atom in kept:
        xyz = np.array((atom.x, atom.y, atom.z), dtype=np.float64)
        if _is_hydrogen(atom):
            hydrogens.append(xyz)
        elif atom.name not in heavy:
            heavy[atom.name] = xyz
            elements.append(atom.element)

    first = kept[0]
    return Residue(chain=first.chain, resseq=first.resseq, icode=first.icode, resname=resname,
                   record=first.record, atoms=heavy, hydrogens=hydrogens, heavy_elements=elements,
                   altloc=chosen, altlocs_present=present)


def build_residues(atoms: list[AtomRecord]) -> list[Residue]:
    """
    Group the atoms of one model into residues, in the order the residues first appear.

    A hydrogen joins the latest residue with its chain, number and insertion
    code, whatever its name: hydrogens are only collected, never told apart by
    name, and some programs list every residue's hydrogens after the whole
    chain or name them all H (Open Babel does both). A hydrogen starts a residue
    only when no residue of its identity has been met.

    Any other atom joins the residue being read when it has that residue's
    chain, number and insertion code and does not repeat one of its atoms (the
    same name and alternate location). Otherwise it starts a new residue, with
    two exceptions for the same residue identity met again later:

    - an atom that carries an alternate location joins an earlier residue of
      the same identity that lacks it and has alternate locations itself,
      because some files write the alternate locations of several residues as
      separate blocks (4A, 5A, 4B, 5B);
    - any other atom starts a new residue, which keeps apart residues that only
      share a number (simulation files whose residue numbers wrap at 10000,
      waters of an mmCIF file that has no author numbering, a ligand given the
      number of a protein residue).

    A residue started by an atom that repeats a name of the residue just
    before it, of the same identity, is marked (``split_by_repeated_name``):
    that tells apart consecutive residues that share a number, but it also
    splits one molecule whose atoms are not uniquely named, so the conversion
    warns when such a residue becomes a bead.

    Parameters
    ----------
    atoms : list of AtomRecord
        The atoms of one model, in file order.

    Returns
    -------
    list of Residue
        One entry per residue.
    """
    groups: list[list[AtomRecord]] = []
    labels: list[set[tuple[str, str]]] = []
    latest_with_identity: dict[tuple[str, str, str], int] = {}
    split_by_repeated_name: set[int] = set()
    for atom in atoms:
        identity = (atom.chain, atom.resseq, atom.icode)
        label = (atom.name, atom.altloc)
        target: Optional[int] = None
        repeats_a_name = False
        if _is_hydrogen(atom) and identity in latest_with_identity:
            target = latest_with_identity[identity]
        if target is None and groups:
            last = groups[-1][0]
            if (last.chain, last.resseq, last.icode) == identity:
                repeats_a_name = label in labels[-1]
                if not repeats_a_name:
                    target = len(groups) - 1
        if target is None and atom.altloc and identity in latest_with_identity:
            earlier = latest_with_identity[identity]
            # a residue with no alternate location of its own is a different molecule
            # that shares the identity (a ligand numbered like a protein residue)
            earlier_has_altlocs = any(altloc for _, altloc in labels[earlier])
            if label not in labels[earlier] and earlier_has_altlocs:
                target = earlier
        if target is None:
            groups.append([])
            labels.append(set())
            target = len(groups) - 1
            latest_with_identity[identity] = target
            if repeats_a_name:
                split_by_repeated_name.add(target)
        groups[target].append(atom)
        labels[target].add(label)

    residues = [_residue_from_atoms(group) for group in groups]
    for index in split_by_repeated_name:
        residues[index].split_by_repeated_name = True
    return residues


# .............................................................................
# Classifying residues and choosing bead codes
# .............................................................................

def _has_atom(residue: Residue, names: Sequence[str]) -> bool:
    """
    Say whether a residue has any of several atoms.

    Parameters
    ----------
    residue : Residue
        The residue.
    names : sequence of str
        Atom names (usually one name in its modern and older spelling).

    Returns
    -------
    bool
        True if at least one is present.
    """
    return any(name in residue.atoms for name in names)


def _atom(residue: Residue, names: Sequence[str]) -> Optional[FloatArray]:
    """
    Return the coordinates of the first of several atoms a residue has.

    Parameters
    ----------
    residue : Residue
        The residue.
    names : sequence of str
        Atom names, most preferred first.

    Returns
    -------
    ndarray or None
        The coordinates (A), or None if the residue has none of the atoms.
    """
    for name in names:
        if name in residue.atoms:
            return residue.atoms[name]
    return None


def _has_amino_backbone(residue: Residue) -> bool:
    """
    Say whether a residue has the backbone of an amino acid.

    That is N, CA and C, or, in a CA-only model, a CA that is the residue's
    only heavy atom and is a carbon. The element matters there because a
    calcium ion is also a lone atom named CA, and it is a hetero group like
    any other ion. A file with no elements gives the ion away by its residue
    name, CA, which no amino acid has.

    Parameters
    ----------
    residue : Residue
        The residue.

    Returns
    -------
    bool
        True if it has an amino-acid backbone.
    """
    if not _has_atom(residue, (CA_ATOM,)):
        return False
    if len(residue.atoms) == 1:
        # heavy_elements lists the elements in the order of residue.atoms
        element = residue.heavy_elements[0] if residue.heavy_elements else ''
        return element in CARBON_ELEMENTS and residue.resname != CALCIUM_ION_NAME
    return _has_atom(residue, (PEPTIDE_N_ATOM,)) and _has_atom(residue, (PEPTIDE_C_ATOM,))


def _has_nucleotide_backbone(residue: Residue) -> bool:
    """
    Say whether a residue has the backbone of a nucleotide.

    That is a C4' with a P, O3' or O5', or, in a coarse model, a C4' or P that is
    the residue's only heavy atom.

    Parameters
    ----------
    residue : Residue
        The residue.

    Returns
    -------
    bool
        True if it has a nucleotide backbone.
    """
    if len(residue.atoms) == 1:
        return _has_atom(residue, C4_PRIME_ATOMS) or _has_atom(residue, (PHOSPHORUS_ATOM,))
    return _has_atom(residue, C4_PRIME_ATOMS) and (
        _has_atom(residue, (PHOSPHORUS_ATOM,)) or _has_atom(residue, O3_PRIME_ATOMS)
        or _has_atom(residue, O5_PRIME_ATOMS))


def _is_deoxyribose(residue: Residue) -> bool:
    """
    Say whether a nucleotide's sugar is clearly a deoxyribose (a C2' with no O2').

    Parameters
    ----------
    residue : Residue
        A nucleotide.

    Returns
    -------
    bool
        True only when the C2' atom is present and the O2' atom is not.
    """
    return _has_atom(residue, C2_PRIME_ATOMS) and not _has_atom(residue, O2_PRIME_ATOMS)


def _nucleotide_kind(residue: Residue) -> str:
    """
    Say whether a nucleotide is RNA or DNA from its sugar.

    Parameters
    ----------
    residue : Residue
        A nucleotide.

    Returns
    -------
    str
        KIND_RNA if it has a 2'-hydroxyl, otherwise KIND_DNA.
    """
    return KIND_RNA if _has_atom(residue, O2_PRIME_ATOMS) else KIND_DNA


def _is_water(residue: Residue) -> bool:
    """
    Say whether a residue is a water molecule.

    A CCD water name (HOH, DOD) is enough. A simulation package's water name
    also has to come with a water's shape: exactly one oxygen and no carbon,
    nitrogen, phosphorus or sulfur.

    Parameters
    ----------
    residue : Residue
        The residue.

    Returns
    -------
    bool
        True for a water molecule.
    """
    if residue.resname in WATER_NAMES:
        return True
    if residue.resname not in SIMULATION_WATER_NAMES:
        return False
    oxygens = sum(1 for element in residue.heavy_elements if element == 'O')
    others = sum(1 for element in residue.heavy_elements if element in NON_WATER_ELEMENTS)
    return oxygens == 1 and others == 0


def _ribonucleotide_or_deoxy(base: str, residue: Residue) -> tuple[str, str, str]:
    """
    Identify A, C or G (or the unknown nucleotide), which RNA and DNA may both be named.

    Parameters
    ----------
    base : str
        'A', 'C' or 'G', or '' for the unknown nucleotide.
    residue : Residue
        The residue, whose sugar decides.

    Returns
    -------
    tuple of (str, str, str)
        (kind, bead code, ``''``).
    """
    kind = KIND_DNA if _is_deoxyribose(residue) else KIND_RNA
    if not base:
        return kind, UNKNOWN_NUCLEOTIDE_CODE, ''
    table = DNA_BASES if kind == KIND_DNA else RNA_BASES
    return kind, table[base], ''


def _standard_identity(name: str, residue: Residue) -> Optional[tuple[str, str, str]]:
    """
    Recognise a standard residue name of the PDB Chemical Component Dictionary.

    These names mean the same in every file, so they are trusted whatever the
    residue's atoms. A, C and G are RNA in a current PDB file but DNA in files
    written before version 3 of the format, so their sugar decides: a
    deoxyribose makes them DNA.

    Parameters
    ----------
    name : str
        The residue name (the residue's own, or a parent's).
    residue : Residue
        The residue.

    Returns
    -------
    tuple of (str, str, str) or None
        (kind, bead code, ``''``), or None if the name is not a standard one.
    """
    if name in AMINO_ACID_CODES:
        return KIND_AMINO_ACID, AMINO_ACID_CODES[name], ''
    if name == UNKNOWN_AMINO_ACID_NAME:
        return KIND_AMINO_ACID, UNKNOWN_AMINO_ACID_CODE, ''
    if len(name) == 2 and name[0] == DNA_PREFIX and name[1] in DNA_BASES:
        return KIND_DNA, DNA_BASES[name[1]], ''
    if name == UNKNOWN_DNA_NAME:
        return KIND_DNA, UNKNOWN_NUCLEOTIDE_CODE, ''
    if name == 'U':
        return KIND_RNA, RNA_BASES['U'], ''
    if name in ('A', 'C', 'G'):
        return _ribonucleotide_or_deoxy(name, residue)
    if name == UNKNOWN_RNA_NAME:
        return _ribonucleotide_or_deoxy('', residue)
    return None


def _simulation_identity(name: str, residue: Residue) -> Optional[tuple[str, str, str]]:
    """
    Recognise a simulation package's name for a standard residue, if the residue has its backbone.

    The names are AMBER, CHARMM and GROMACS protonation states (HID, HSE,
    HISE, CYX, ...), AMBER terminal residues (NALA, NHIE, CCYX, DA5, RA3, A5), AMBER's
    R-prefixed RNA (RA), version-2 PDB thymidine (T) and CHARMM's base-only
    nucleotides (ADE, ...). Many collide with CCD codes of unrelated molecules
    (RU is a ruthenium ion, GUA glutaric acid, T3 a thyroid hormone), so a name
    is accepted only for a residue with the backbone of what it claims to be.

    Parameters
    ----------
    name : str
        The residue name.
    residue : Residue
        The residue.

    Returns
    -------
    tuple of (str, str, str) or None
        (kind, bead code, the standard name it stands for), or None.
    """
    if _has_amino_backbone(residue):
        if name in AMINO_ACID_VARIANTS:
            standard = AMINO_ACID_VARIANTS[name]
            return KIND_AMINO_ACID, AMINO_ACID_CODES[standard], standard
        if len(name) == 4 and name[0] in TERMINAL_PREFIXES and name[1:] in AMINO_ACID_CODES:
            return KIND_AMINO_ACID, AMINO_ACID_CODES[name[1:]], name[1:]
        if len(name) == 4 and name[0] in TERMINAL_PREFIXES and name[1:] in AMBER_TERMINAL_VARIANTS:
            standard = AMINO_ACID_VARIANTS[name[1:]]
            return KIND_AMINO_ACID, AMINO_ACID_CODES[standard], standard

    if not _has_nucleotide_backbone(residue):
        return None
    stem = name[:-1] if len(name) > 1 and name[-1] in TERMINAL_SUFFIXES else name
    if len(stem) == 2 and stem[0] == RNA_PREFIX and stem[1] in RNA_BASES:
        stem = stem[1:]
    if stem == LEGACY_THYMIDINE:
        return KIND_DNA, DNA_BASES['T'], stem
    if stem != name:
        identity = _standard_identity(stem, residue)
        if identity is not None and identity[0] != KIND_AMINO_ACID:
            return identity[0], identity[1], stem
    if name in BASE_ONLY_NUCLEOTIDES:
        base = BASE_ONLY_NUCLEOTIDES[name]
        kind = _nucleotide_kind(residue)
        table = RNA_BASES if kind == KIND_RNA else DNA_BASES
        if base in table:
            return kind, table[base], base
    return None


def classify_residue(residue: Residue, file_parents: dict[str, str]) -> None:
    """
    Decide a residue's kind, bead code and (for a modified residue) parent.

    In order: water; a standard CCD residue name; a modified residue whose
    parent the file gives; a simulation package's name for a standard residue
    (if the residue has that residue's backbone); a modified residue whose
    parent the built-in table gives; a residue of unknown name with an
    amino-acid backbone (code X) or a nucleotide backbone (code n); and
    anything else is a hetero group, whose code is chosen later if it is kept.

    Parameters
    ----------
    residue : Residue
        The residue; its kind, code, parent and parent_source are set in place.
    file_parents : dict of str to str
        Modified residue -> parent, from the structure file.

    Returns
    -------
    None
    """
    name = residue.resname
    if _is_water(residue):
        residue.kind = KIND_WATER
        return

    identity = _standard_identity(name, residue)
    if identity is not None:
        residue.kind, residue.code, _ = identity
        if residue.code in (UNKNOWN_AMINO_ACID_CODE, UNKNOWN_NUCLEOTIDE_CODE):
            residue.parent_source = 'unknown'
        return

    if name in file_parents:
        identity = _standard_identity(file_parents[name], residue)
        if identity is not None:
            residue.kind, residue.code, _ = identity
            residue.parent, residue.parent_source = file_parents[name], 'file'
            return

    identity = _simulation_identity(name, residue)
    if identity is not None:
        residue.kind, residue.code, residue.parent = identity
        residue.parent_source = 'other name'
        return

    if name in FALLBACK_PARENTS:
        identity = _standard_identity(FALLBACK_PARENTS[name], residue)
        if identity is not None:
            residue.kind, residue.code, _ = identity
            residue.parent, residue.parent_source = FALLBACK_PARENTS[name], 'built-in table'
            return

    if _has_amino_backbone(residue):
        residue.kind, residue.code, residue.parent_source = KIND_AMINO_ACID, UNKNOWN_AMINO_ACID_CODE, 'unknown'
        return
    if _has_nucleotide_backbone(residue):
        residue.kind, residue.code, residue.parent_source = (_nucleotide_kind(residue), UNKNOWN_NUCLEOTIDE_CODE,
                                                             'unknown')
        return
    residue.kind = KIND_HETERO


def check_bead_character(character: str, context: str) -> None:
    """
    Refuse a bead character that PIMMS cannot use.

    Parameters
    ----------
    character : str
        The proposed bead character.
    context : str
        Where it came from, for the message.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If it is not a single printable, non-space character other than 0 and #.
    """
    if len(character) != 1:
        raise ValueError("%s: a bead type must be exactly one character (PIMMS sequences have one "
                         "character per bead), got %r" % (context, character))
    if character.isspace() or not character.isprintable():
        raise ValueError("%s: %r cannot be a bead type (whitespace and unprintable characters "
                         "separate parameter-file columns)" % (context, character))
    if character in RESERVED_BEAD_CHARACTERS:
        raise ValueError("%s: %r cannot be a bead type (0 is the solvent and # starts a "
                         "parameter-file comment)" % (context, character))


def apply_bead_codes(residues: list[Residue], bead_map: dict[str, str], keep_hetero: bool) -> list[str]:
    """
    Apply --map overrides and give each kept hetero group a code.

    A hetero group with no --map entry gets the next unused letter of
    HETERO_CODE_POOL, one per residue name, in the order the names first appear.

    Parameters
    ----------
    residues : list of Residue
        Every residue, already classified; codes are set in place.
    bead_map : dict of str to str
        Residue name -> bead character, already checked with check_bead_character().
    keep_hetero : bool
        Whether hetero groups become beads.

    Returns
    -------
    list of str
        Warnings about --map entries that had no effect.

    Raises
    ------
    StructureConversionException
        If more distinct hetero groups are kept than the pool has letters.
    """
    warnings = []
    present = {residue.resname for residue in residues}
    for name in sorted(bead_map):
        if name not in present:
            warnings.append("--map %s=%s had no effect: the converted model has no residue named %s"
                            % (name, bead_map[name], name))
        elif name in WATER_NAMES:
            warnings.append("--map %s=%s had no effect: water is never converted" % (name, bead_map[name]))
        elif not keep_hetero and all(r.kind == KIND_HETERO for r in residues if r.resname == name):
            warnings.append("--map %s=%s had no effect: %s is a hetero group, and hetero groups are "
                            "only converted with --hetero" % (name, bead_map[name], name))

    for residue in residues:
        if residue.kind != KIND_WATER and residue.resname in bead_map:
            residue.code = bead_map[residue.resname]

    if not keep_hetero:
        return warnings

    used = {residue.code for residue in residues if residue.code}
    hetero_codes: dict[str, str] = {}
    pool = [letter for letter in HETERO_CODE_POOL if letter not in used]
    for residue in residues:
        if residue.kind != KIND_HETERO or residue.code:
            continue
        if residue.resname not in hetero_codes:
            if not pool:
                raise StructureConversionException(
                    "The model holds more distinct hetero groups than pdb2pimms has spare bead letters "
                    "for (%d). Give the rest a bead type each with --map RESNAME=CHARACTER."
                    % len(HETERO_CODE_POOL))
            hetero_codes[residue.resname] = pool.pop(0)
        residue.code = hetero_codes[residue.resname]
    return warnings


# .............................................................................
# Bead positions and chain connectivity
# .............................................................................

def _centroid(residue: Residue) -> FloatArray:
    """
    The mean position of a residue's heavy atoms (of all its atoms if it has no heavy atom).

    Parameters
    ----------
    residue : Residue
        The residue (it has at least one atom).

    Returns
    -------
    ndarray
        The centroid (A).
    """
    points = list(residue.atoms.values()) or residue.hydrogens
    return np.mean(np.array(points), axis=0)


def representative_position(residue: Residue, bead_position: str) -> None:
    """
    Choose the point a residue's bead stands for.

    With 'backbone' an amino acid is placed at its CA and a nucleotide at its
    C4' (the backbone atom every nucleotide has, including a 5'-terminal one
    with no phosphate); a residue missing that atom, and every hetero group, is
    placed at its heavy-atom centroid. With 'centroid' every residue is placed
    at its heavy-atom centroid.

    Parameters
    ----------
    residue : Residue
        A classified residue; position and position_source are set in place.
    bead_position : str
        BEAD_POSITION_BACKBONE or BEAD_POSITION_CENTROID.

    Returns
    -------
    None
    """
    if bead_position == BEAD_POSITION_BACKBONE and residue.kind in POLYMER_KINDS:
        names = (CA_ATOM,) if residue.kind == KIND_AMINO_ACID else C4_PRIME_ATOMS
        xyz = _atom(residue, names)
        if xyz is not None:
            residue.position = xyz
            residue.position_source = names[0]
            return
        residue.position = _centroid(residue)
        residue.position_source = 'heavy-atom centroid (no %s)' % names[0]
        return
    residue.position = _centroid(residue)
    residue.position_source = 'heavy-atom centroid'


def _bond_between(previous: Residue, current: Residue) -> tuple[bool, str]:
    """
    Decide whether two consecutive residues of a chain are bonded.

    The linking bond itself is tested when both of its atoms are present (C-N
    for a peptide, O3'-P for a phosphodiester); otherwise the distance between
    backbone atoms is (CA-CA, or C4'-C4' and then P-P for a nucleic acid).

    Parameters
    ----------
    previous, current : Residue
        Consecutive polymer residues of one chain.

    Returns
    -------
    tuple of (bool, str)
        Whether they are bonded, and (when not) why not.
    """
    previous_protein = previous.kind == KIND_AMINO_ACID
    if previous_protein != (current.kind == KIND_AMINO_ACID):
        return False, '%s followed by %s' % (previous.kind, current.kind)

    if previous_protein:
        link = (_atom(previous, (PEPTIDE_C_ATOM,)), _atom(current, (PEPTIDE_N_ATOM,)))
        link_name, link_max, link_length = 'C-N', MAX_LINK_BOND_ANGSTROMS, PEPTIDE_BOND_ANGSTROMS
        backbones: list[tuple[tuple[str, ...], str, float]] = [((CA_ATOM,), 'CA-CA', MAX_CA_CA_ANGSTROMS)]
    else:
        link = (_atom(previous, O3_PRIME_ATOMS), _atom(current, (PHOSPHORUS_ATOM,)))
        link_name, link_max, link_length = "O3'-P", MAX_LINK_BOND_ANGSTROMS, PHOSPHODIESTER_BOND_ANGSTROMS
        backbones = [(C4_PRIME_ATOMS, "C4'-C4'", MAX_C4_C4_ANGSTROMS),
                     ((PHOSPHORUS_ATOM,), 'P-P', MAX_P_P_ANGSTROMS)]

    first, second = link
    if first is not None and second is not None:
        distance = float(np.linalg.norm(second - first))
        if distance <= link_max:
            return True, ''
        return False, '%s distance %.2f A (a bonded %s is %.2f A)' % (link_name, distance, link_name, link_length)

    for names, backbone_name, backbone_max in backbones:
        first, second = _atom(previous, names), _atom(current, names)
        if first is None or second is None:
            continue
        distance = float(np.linalg.norm(second - first))
        if distance <= backbone_max:
            return True, ''
        return False, ('%s distance %.2f A (more than %.1f A; the %s bond itself could not be tested '
                       'because one of its atoms is missing)' % (backbone_name, distance, backbone_max, link_name))
    return False, 'neither the %s atoms nor any backbone atom pair is present to test the bond' % link_name


def _full_atom_residue_lacks_link_atom(previous: Residue, current: Residue) -> bool:
    """
    Say whether either residue of a pair lacks its atom of the linking bond although it is not a coarse-model residue.

    _bond_between() falls back to the backbone distance (CA-CA, C4'-C4', P-P)
    when a linking atom is missing, which is what a CA-only or C4'-only model
    needs. A residue that has more than those trace atoms but no linking atom
    cannot be linked on that side at all: a nucleoside, or a nucleoside
    triphosphate such as ATP, whose alpha phosphorus the CCD names PA rather
    than P. For such a pair the fallback would join a ligand to a chain end
    within 7.5 A.

    Parameters
    ----------
    previous, current : Residue
        Consecutive polymer residues of one chain, of the same kind.

    Returns
    -------
    bool
        True if a residue that is more than a trace lacks its linking atom
        (C or O3' for ``previous``, N or P for ``current``).
    """
    sides: tuple[tuple[Residue, tuple[str, ...]], ...]
    if previous.kind == KIND_AMINO_ACID:
        sides = ((previous, (PEPTIDE_C_ATOM,)), (current, (PEPTIDE_N_ATOM,)))
    else:
        sides = ((previous, O3_PRIME_ATOMS), (current, (PHOSPHORUS_ATOM,)))
    for residue, link_names in sides:
        is_trace = all(name in TRACE_ATOMS for name in residue.atoms)
        if not _has_atom(residue, link_names) and not is_trace:
            return True
    return False


def demote_isolated_hetatm(residues: list[Residue]) -> list[Residue]:
    """
    Treat a HETATM residue that is bonded to neither neighbour as a hetero group, whatever its name.

    A deposited structure writes a free amino acid or nucleotide bound as a
    ligand (an arginine in a binding pocket, an AMP, an ATP) as a HETATM
    residue with a standard name or a nucleotide's backbone, which would
    otherwise become a one-bead polymer chain or join the end of a chain. A
    modified residue inside a chain is also HETATM, but it is bonded to its
    neighbours and is kept. A HETATM residue recognised by its backbone alone
    (no name, MODRES record or table gives its parent: an ATP, a dNTP, a
    nucleoside) needs the linking atoms (C-N, O3'-P) of every residue of the
    pair that is more than a coarse-model trace, so such a ligand is never
    joined to a chain by the backbone distance alone. A named modified residue
    (MSE, PSU, 5MC) keeps the distance fallback, because at a chain end its
    linking atom, or its neighbour's, is often just not modelled.

    Parameters
    ----------
    residues : list of Residue
        Every residue, classified; demoted residues have their kind, code and
        parent reset in place.

    Returns
    -------
    list of Residue
        The residues that were demoted.
    """
    polymer = [residue for residue in residues if residue.kind in POLYMER_KINDS]
    demoted = []
    for index, residue in enumerate(polymer):
        if residue.record != 'HETATM':
            continue
        bonded = False
        for neighbour_index in (index - 1, index + 1):
            if 0 <= neighbour_index < len(polymer):
                neighbour = polymer[neighbour_index]
                if neighbour.chain == residue.chain:
                    pair = (neighbour, residue) if neighbour_index < index else (residue, neighbour)
                    # the linking-atom rule only for a residue known by its backbone
                    # alone: applied to every HETATM residue it also dropped named
                    # modified residues at chain ends whose link was not modelled
                    strict = (residue.parent_source == 'unknown'
                              and _full_atom_residue_lacks_link_atom(*pair))
                    linked = _bond_between(*pair)[0] and not strict
                    bonded = bonded or linked
        if not bonded:
            demoted.append(residue)
    for residue in demoted:
        residue.kind, residue.code, residue.parent, residue.parent_source = KIND_HETERO, '', '', ''
    return demoted


def split_into_segments(residues: list[Residue]) -> tuple[list[Segment], list[ChainBreak]]:
    """
    Turn the polymer residues into PIMMS chains.

    Consecutive polymer residues of one chain that are bonded go into the same
    PIMMS chain; a chain with a gap (missing residues) becomes one PIMMS chain
    per bonded run, and each split is recorded. Chains come in file order.

    Parameters
    ----------
    residues : list of Residue
        Every residue, classified (only the polymer residues are used).

    Returns
    -------
    tuple of (list of Segment, list of ChainBreak)
        The PIMMS chains and the breaks.
    """
    segments: list[Segment] = []
    breaks: list[ChainBreak] = []
    previous: Optional[Residue] = None
    for residue in residues:
        if residue.kind not in POLYMER_KINDS:
            continue
        if previous is not None and previous.chain == residue.chain:
            bonded, reason = _bond_between(previous, residue)
            if bonded:
                segments[-1].residues.append(residue)
                previous = residue
                continue
            breaks.append(ChainBreak(residue.chain, previous.label, residue.label, reason))
        segments.append(Segment([residue], residue.chain))
        previous = residue
    return segments, breaks


# .............................................................................
# Fitting chains onto the lattice
# .............................................................................

def _cube_offsets(radius: int) -> IntArray:
    """
    Every integer offset within a Chebyshev distance, x varying slowest.

    Parameters
    ----------
    radius : int
        The Chebyshev radius.

    Returns
    -------
    ndarray
        Shape ((2 radius + 1)^3, 3).
    """
    span = range(-radius, radius + 1)
    return np.array(list(itertools.product(span, span, span)), dtype=np.int64)


def _offset_index(offset: IntArray, radius: int) -> Optional[int]:
    """
    The row of _cube_offsets(radius) that holds an offset.

    Parameters
    ----------
    offset : ndarray
        Shape (3,).
    radius : int
        The Chebyshev radius.

    Returns
    -------
    int or None
        The row, or None if the offset lies outside the radius.
    """
    if np.any(np.abs(offset) > radius):
        return None
    width = 2 * radius + 1
    x, y, z = (int(value) + radius for value in offset)
    return (x * width + y) * width + z


def _site_keys(sites: IntArray) -> IntArray:
    """
    Pack lattice sites into single integers so they can be tested against a set with numpy.

    Parameters
    ----------
    sites : ndarray
        Integer sites, shape (..., 3), every coordinate within SITE_KEY_OFFSET of zero.

    Returns
    -------
    ndarray
        One int64 per site, shape (...).
    """
    shifted = sites + SITE_KEY_OFFSET
    return ((shifted[..., 0] << (2 * SITE_KEY_BITS)) | (shifted[..., 1] << SITE_KEY_BITS)
            | shifted[..., 2])


def _squared_distance(site: Site, target: FloatArray) -> float:
    """
    Squared distance between a site and a target position (lattice units).

    Parameters
    ----------
    site : tuple of int
        The site.
    target : ndarray
        The target, shape (3,).

    Returns
    -------
    float
        The squared distance.
    """
    return float(np.sum((np.array(site, dtype=np.float64) - target) ** 2))


def _optimal_walk(targets: FloatArray, occupied_keys: IntArray, banned: dict[int, set[Site]],
                  radius: int) -> Optional[list[Site]]:
    """
    Find the lattice walk closest to a chain's targets by dynamic programming.

    Each bead may sit at any site within ``radius`` of the site nearest its
    target, except sites already taken by earlier chains and any site banned
    for that bead. Consecutive beads must be Moore neighbours. Among all such
    walks the one with the smallest sum of squared bead-to-target distances is
    found exactly (a Viterbi recursion over beads). The recursion cannot see
    the walk's own history, so the walk may return to a site it already used;
    _self_avoiding_walk() repairs that.

    Parameters
    ----------
    targets : ndarray
        Target positions in lattice units, shape (n, 3).
    occupied_keys : ndarray
        _site_keys() of every site already occupied.
    banned : dict of int to set of tuple
        Bead index -> sites that bead may not use.
    radius : int
        The Chebyshev radius of each bead's candidate sites.

    Returns
    -------
    list of tuple or None
        The walk, or None if no walk exists within the radius.
    """
    offsets = _cube_offsets(radius)
    n_beads = len(targets)
    n_states = len(offsets)
    nearest = np.rint(targets).astype(np.int64)
    candidates = nearest[:, None, :] + offsets[None, :, :]
    cost = np.sum((candidates - targets[:, None, :]) ** 2, axis=2)
    if occupied_keys.size:
        cost[np.isin(_site_keys(candidates), occupied_keys)] = np.inf
    for bead, sites in banned.items():
        for site in sites:
            index = _offset_index(np.array(site, dtype=np.int64) - nearest[bead], radius)
            if index is not None:
                cost[bead, index] = np.inf

    # state a of bead j and state b of bead j-1 are a step of (nearest[j] -
    # nearest[j-1]) + pair_offsets[a, b] apart; the allowed-transition mask
    # depends only on that first term, which takes few distinct values
    pair_offsets = offsets[:, None, :] - offsets[None, :, :]
    masks: dict[Site, BoolArray] = {}
    accumulated = cost[0].copy()
    backpointers = np.zeros((max(n_beads - 1, 0), n_states), dtype=np.int32)
    every_state = np.arange(n_states)
    for bead in range(1, n_beads):
        step = nearest[bead] - nearest[bead - 1]
        key = (int(step[0]), int(step[1]), int(step[2]))
        if key not in masks:
            masks[key] = np.max(np.abs(pair_offsets + step), axis=2) == 1
        scores = np.where(masks[key], accumulated[None, :], np.inf)
        best_previous = np.argmin(scores, axis=1)
        accumulated = scores[every_state, best_previous] + cost[bead]
        backpointers[bead - 1] = best_previous
        if not np.isfinite(accumulated).any():
            return None

    state = int(np.argmin(accumulated))
    if not math.isfinite(float(accumulated[state])):
        return None
    states = [state]
    for bead in range(n_beads - 1, 0, -1):
        state = int(backpointers[bead - 1][state])
        states.append(state)
    states.reverse()
    walk = []
    for bead, state in enumerate(states):
        x, y, z = (int(value) for value in candidates[bead, state])
        walk.append((x, y, z))
    return walk


def _shared_sites(walk: list[Site]) -> dict[Site, list[int]]:
    """
    Find the sites a walk visits more than once.

    Parameters
    ----------
    walk : list of tuple
        The walk.

    Returns
    -------
    dict of tuple to list of int
        Each reused site -> the indices of the beads on it (empty if the walk is self-avoiding).
    """
    visits: dict[Site, list[int]] = {}
    for index, site in enumerate(walk):
        visits.setdefault(site, []).append(index)
    return {site: beads for site, beads in visits.items() if len(beads) > 1}


def _self_avoiding_walk(targets: FloatArray, occupied_keys: IntArray, radius: int) -> Optional[list[Site]]:
    """
    Find a close self-avoiding walk by repairing the optimal walk where it reuses sites.

    While the optimal walk puts several beads on one site, that site is kept for
    the bead whose target is nearest it (the earlier bead on a tie) and banned
    for the others, for every reused site at once, and the walk is found again.
    Each round bans at least one more (bead, site) pair, so the loop ends; it is
    capped at MAX_CLASH_REPAIR_ROUNDS rounds.

    Parameters
    ----------
    targets : ndarray
        Target positions in lattice units, shape (n, 3).
    occupied_keys : ndarray
        _site_keys() of every site already occupied.
    radius : int
        The Chebyshev radius of each bead's candidate sites.

    Returns
    -------
    list of tuple or None
        A self-avoiding walk, or None if none was found.
    """
    banned: dict[int, set[Site]] = {}
    for _ in range(MAX_CLASH_REPAIR_ROUNDS):
        walk = _optimal_walk(targets, occupied_keys, banned, radius)
        if walk is None:
            return None
        shared = _shared_sites(walk)
        if not shared:
            return walk
        for site, beads in shared.items():
            keeper = min(beads, key=lambda bead: (_squared_distance(site, targets[bead]), bead))
            for bead in beads:
                if bead != keeper:
                    banned.setdefault(bead, set()).add(site)
    return None


def _start_candidates(target: FloatArray, occupied: set[Site]) -> list[Site]:
    """
    Free sites for the first bead of a greedy walk, nearest its target first.

    Parameters
    ----------
    target : ndarray
        The bead's target (lattice units), shape (3,).
    occupied : set of tuple
        Sites already taken.

    Returns
    -------
    list of tuple
        Free sites within the nearest Chebyshev shell that has any.
    """
    nearest = np.rint(target).astype(np.int64)
    for radius in range(GREEDY_START_SEARCH_RADIUS + 1):
        sites = []
        for offset in _cube_offsets(radius):
            x, y, z = (int(value) for value in nearest + offset)
            if (x, y, z) not in occupied:
                sites.append((x, y, z))
        if sites:
            return sorted(sites, key=lambda site: (_squared_distance(site, target), site))
    return []


def _step_candidates(previous: Site, target: FloatArray, occupied: set[Site]) -> list[Site]:
    """
    Free Moore neighbours of a site, nearest a target first.

    Parameters
    ----------
    previous : tuple of int
        The site of the previous bead.
    target : ndarray
        The next bead's target (lattice units), shape (3,).
    occupied : set of tuple
        Sites already taken.

    Returns
    -------
    list of tuple
        The free neighbours, sorted.
    """
    sites = []
    for offset in MOORE_OFFSETS:
        site = (previous[0] + int(offset[0]), previous[1] + int(offset[1]), previous[2] + int(offset[2]))
        if site not in occupied:
            sites.append(site)
    return sorted(sites, key=lambda site: (_squared_distance(site, target), site))


def _greedy_walk(targets: FloatArray, occupied: set[Site]) -> Optional[list[Site]]:
    """
    Place a chain bead by bead at the free neighbour nearest each target, backtracking at dead ends.

    This is the fallback for a chain the dynamic programme cannot place. It
    always keeps the walk valid, but a choice is never revisited unless it leads
    to a dead end, so its fit is usually worse.

    Parameters
    ----------
    targets : ndarray
        Target positions in lattice units, shape (n, 3).
    occupied : set of tuple
        Sites already taken by earlier chains (not changed).

    Returns
    -------
    list of tuple or None
        The walk, or None if none was found within the placement budget.
    """
    n_beads = len(targets)
    budget = GREEDY_PLACEMENTS_PER_BEAD * n_beads + GREEDY_PLACEMENTS_MINIMUM
    taken = set(occupied)
    walk: list[Site] = []
    # pending[k] holds the untried sites for bead k
    pending = [_start_candidates(targets[0], taken)]
    placements = 0
    while len(walk) < n_beads:
        options = pending[-1]
        if options:
            site = options.pop(0)
            walk.append(site)
            taken.add(site)
            placements += 1
            if placements > budget:
                return None
            if len(walk) < n_beads:
                pending.append(_step_candidates(site, targets[len(walk)], taken))
            continue
        # no site left for this bead: take back the previous bead and try its next option
        pending.pop()
        if not walk:
            return None
        taken.discard(walk.pop())
    return walk


def fit_to_lattice(targets: list[FloatArray]) -> tuple[list[list[Site]], list[str]]:
    """
    Place every chain on the lattice, in order, without overlaps.

    Each chain is placed by the dynamic programme (with self-clash repair),
    widening its search radius once if needed, and by the greedy search if that
    fails. Earlier chains are fixed when later ones are placed.

    Parameters
    ----------
    targets : list of ndarray
        Per chain, target positions in lattice units, shape (n, 3).

    Returns
    -------
    tuple of (list of list of tuple, list of str)
        The sites of every chain (in an unshifted frame) and how each was placed.

    Raises
    ------
    StructureConversionException
        If a chain cannot be placed at all.
    """
    occupied: set[Site] = set()
    occupied_keys = np.zeros(0, dtype=np.int64)
    walks = []
    methods = []
    for index, chain_targets in enumerate(targets):
        walk: Optional[list[Site]] = None
        method = ''
        for radius in DYNAMIC_PROGRAMME_RADII:
            walk = _self_avoiding_walk(chain_targets, occupied_keys, radius)
            if walk is not None:
                method = 'dynamic programming (radius %d)' % radius
                break
        if walk is None:
            walk = _greedy_walk(chain_targets, occupied)
            method = 'greedy search'
        if walk is None:
            raise StructureConversionException(
                "PIMMS chain %d (%d beads) could not be placed on the lattice without overlapping "
                "itself or an earlier chain. The beads are packed more densely than the lattice can "
                "hold at this spacing: try a smaller --spacing." % (index + 1, len(chain_targets)))
        walks.append(walk)
        methods.append(method)
        occupied.update(walk)
        occupied_keys = np.concatenate([occupied_keys, _site_keys(np.array(walk, dtype=np.int64))])
    return walks, methods


def place_in_box(walks: list[list[Site]], padding: int,
                 box: Optional[list[int]]) -> tuple[list[list[Site]], list[int]]:
    """
    Choose the box and move the beads into it.

    Without an explicit box, each axis is the structure's extent plus
    ``padding`` empty sites on either side, and at least MIN_BOX_SITES. The
    structure is centred in the box either way.

    Parameters
    ----------
    walks : list of list of tuple
        The sites of every chain, in the fitting frame.
    padding : int
        Empty sites on each side (ignored when ``box`` is given).
    box : list of int or None
        An explicit box.

    Returns
    -------
    tuple of (list of list of tuple, list of int)
        The sites in the box frame (every coordinate in [0, box)) and the box.

    Raises
    ------
    StructureConversionException
        If an explicit box is too small for the structure or for PIMMS.
    """
    every_site = np.array([site for walk in walks for site in walk], dtype=np.int64)
    low = every_site.min(axis=0)
    extent = every_site.max(axis=0) - low + 1
    if box is None:
        dimensions = np.maximum(extent + 2 * padding, MIN_BOX_SITES)
    else:
        dimensions = np.array(box, dtype=np.int64)
        if np.any(dimensions < extent) or np.any(dimensions < MIN_BOX_SITES):
            raise StructureConversionException(
                "--box %s cannot be used: the structure spans %s lattice sites, and PIMMS needs every "
                "axis to be at least %d" % (' '.join(str(v) for v in box),
                                           ' x '.join(str(int(v)) for v in extent), MIN_BOX_SITES))

    shift = (dimensions - extent) // 2 - low
    placed = []
    for walk in walks:
        moved = []
        for site in walk:
            x, y, z = (int(value) for value in np.array(site, dtype=np.int64) + shift)
            moved.append((x, y, z))
        placed.append(moved)
    return placed, [int(value) for value in dimensions]


# .............................................................................
# Fit statistics
# .............................................................................

def _superposed_deviations(reference: FloatArray, mobile: FloatArray) -> FloatArray:
    """
    Per-point deviations after the optimal rigid superposition of one point set onto another.

    Kabsch (Acta Cryst. A32, 922, 1976): centre both sets, take the rotation from
    the singular value decomposition of their covariance, and correct it to a
    proper rotation (no reflection).

    Parameters
    ----------
    reference, mobile : ndarray
        Matching point sets, shape (n, 3).

    Returns
    -------
    ndarray
        The distance of each superposed mobile point from its reference point, shape (n,).
    """
    reference_centred = reference - reference.mean(axis=0)
    mobile_centred = mobile - mobile.mean(axis=0)
    if len(reference) < 2:
        return np.sqrt(np.sum((mobile_centred - reference_centred) ** 2, axis=1))
    u, _, vt = np.linalg.svd(mobile_centred.T @ reference_centred)
    handedness = float(np.sign(np.linalg.det(u @ vt))) or 1.0
    rotation = u @ np.diag([1.0, 1.0, handedness]) @ vt
    superposed = mobile_centred @ rotation
    return np.sqrt(np.sum((superposed - reference_centred) ** 2, axis=1))


def _radius_of_gyration(points: FloatArray) -> float:
    """
    Radius of gyration of a set of points (equal weights).

    Parameters
    ----------
    points : ndarray
        Shape (n, 3).

    Returns
    -------
    float
        The radius of gyration, in the units of ``points``.
    """
    centred = points - points.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centred ** 2, axis=1))))


def fit_statistics(segments: list[Segment], sites: list[list[Site]], spacing: float) -> FitStatistics:
    """
    Measure how closely the lattice configuration reproduces the input bead positions.

    Parameters
    ----------
    segments : list of Segment
        The PIMMS chains (with bead positions).
    sites : list of list of tuple
        The lattice site of every bead.
    spacing : float
        Lattice spacing (A per site).

    Returns
    -------
    FitStatistics
        The fit statistics.
    """
    input_points = []
    lattice_points = []
    labels = []
    for segment, walk in zip(segments, sites):
        for residue, site in zip(segment.residues, walk):
            assert residue.position is not None
            input_points.append(residue.position)
            lattice_points.append(np.array(site, dtype=np.float64) * spacing)
            labels.append(residue.label)
    reference = np.array(input_points)
    mobile = np.array(lattice_points)
    deviations = _superposed_deviations(reference, mobile)

    segment_rmsd = []
    input_rg = []
    lattice_rg = []
    bond_steps: Counter[int] = Counter()
    input_bonds = []
    lattice_bonds = []
    start = 0
    for segment, walk in zip(segments, sites):
        stop = start + len(walk)
        segment_rmsd.append(float(np.sqrt(np.mean(deviations[start:stop] ** 2))))
        input_rg.append(_radius_of_gyration(reference[start:stop]))
        lattice_rg.append(_radius_of_gyration(mobile[start:stop]))
        for bead in range(start, stop - 1):
            input_bonds.append(float(np.linalg.norm(reference[bead + 1] - reference[bead])))
            lattice_bonds.append(float(np.linalg.norm(mobile[bead + 1] - mobile[bead])))
            step = np.array(walk[bead - start + 1]) - np.array(walk[bead - start])
            bond_steps[int(np.count_nonzero(step))] += 1
        start = stop

    worst = int(np.argmax(deviations))
    return FitStatistics(
        rmsd=float(np.sqrt(np.mean(deviations ** 2))), segment_rmsd=segment_rmsd,
        max_deviation=float(deviations[worst]), max_deviation_label=labels[worst],
        bond_steps=bond_steps, input_bond_lengths=np.array(input_bonds, dtype=np.float64),
        lattice_bond_lengths=np.array(lattice_bonds, dtype=np.float64),
        input_rg=input_rg, lattice_rg=lattice_rg)


# .............................................................................
# The conversion
# .............................................................................

def assign_chain_types(segments: list[Segment]) -> list[int]:
    """
    Give chains with identical sequences the same chainType, numbered from 0 in order of appearance.

    PIMMS refuses a restart file in which one chainType holds two different
    sequences, and numbers the chain types of a keyfile from 0 the same way.

    Parameters
    ----------
    segments : list of Segment
        The PIMMS chains.

    Returns
    -------
    list of int
        The chainType of each chain.
    """
    type_of_sequence: dict[str, int] = {}
    chain_types = []
    for segment in segments:
        sequence = segment.sequence
        if sequence not in type_of_sequence:
            type_of_sequence[sequence] = len(type_of_sequence)
        chain_types.append(type_of_sequence[sequence])
    return chain_types


def _file_sha256(path: Path) -> str:
    """
    SHA-256 digest of a file, so the report pins exactly which input was converted.

    Parameters
    ----------
    path : Path
        The file.

    Returns
    -------
    str
        The hexadecimal digest.
    """
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _select_model(structure: StructureFile, model: Optional[int]) -> int:
    """
    Pick the model to convert.

    Parameters
    ----------
    structure : StructureFile
        The structure.
    model : int or None
        The requested model number, or None for the first in the file.

    Returns
    -------
    int
        The model number.

    Raises
    ------
    StructureConversionException
        If the requested model is not in the file.
    """
    if model is None:
        return structure.models[0]
    if model not in structure.models:
        shown = ', '.join(str(number) for number in structure.models[:20])
        more = ' ...' if len(structure.models) > 20 else ''
        raise StructureConversionException(
            "--model %d is not in %s, which holds model(s) %s%s"
            % (model, structure.path, shown, more))
    return model


def _select_chains(residues: list[Residue], chains: Optional[list[str]]) -> tuple[list[Residue], list[str]]:
    """
    Keep the residues of the requested chains.

    Parameters
    ----------
    residues : list of Residue
        Every residue of the model.
    chains : list of str or None
        The chain identifiers to keep, or None for all.

    Returns
    -------
    tuple of (list of Residue, list of str)
        The residues kept and the chains left out.

    Raises
    ------
    StructureConversionException
        If a requested chain is not in the model.
    """
    present = list(dict.fromkeys(residue.chain for residue in residues))
    if chains is None:
        return residues, []
    # a blank chain identifier is asked for as '-', the way the report shows it
    wanted = ['' if chain == BLANK_CHAIN_LABEL and chain not in present else chain for chain in chains]
    missing = [chain for chain, asked in zip(wanted, chains) if chain not in present]
    if missing:
        raise StructureConversionException(
            "--chains asks for chain(s) %s, which the model does not have (it has %s)"
            % (', '.join(chain or BLANK_CHAIN_LABEL for chain in missing),
               ', '.join(chain or '%s (blank)' % BLANK_CHAIN_LABEL for chain in present)))
    kept = [residue for residue in residues if residue.chain in wanted]
    return kept, [chain or BLANK_CHAIN_LABEL for chain in present if chain not in wanted]


def _summarise_exclusions(structure: StructureFile, model: int, residues: list[Residue],
                          unselected: list[str], keep_hetero: bool) -> ExclusionSummary:
    """
    Count what the conversion read but did not turn into beads.

    Parameters
    ----------
    structure : StructureFile
        The structure.
    model : int
        The model converted.
    residues : list of Residue
        The residues of the selected chains, classified.
    unselected : list of str
        Chains left out by --chains.
    keep_hetero : bool
        Whether hetero groups were converted.

    Returns
    -------
    ExclusionSummary
        The counts.
    """
    summary = ExclusionSummary()
    summary.waters = sum(1 for residue in residues if residue.kind == KIND_WATER)
    if not keep_hetero:
        summary.hetero = Counter(residue.resname for residue in residues if residue.kind == KIND_HETERO)
    summary.hydrogens = sum(len(residue.hydrogens) for residue in residues)
    converted_kinds = POLYMER_KINDS + ((KIND_HETERO,) if keep_hetero else ())
    summary.altloc_residues = [residue for residue in residues
                               if residue.altlocs_present and residue.kind in converted_kinds]
    summary.other_models = len(structure.models) - 1
    summary.unselected_chains = unselected
    return summary


def _conversion_warnings(segments: list[Segment], residues: list[Residue], statistics: FitStatistics,
                         dimensions: list[int], spacing: float) -> list[str]:
    """
    Collect what the user should check about a finished conversion.

    Parameters
    ----------
    segments : list of Segment
        The PIMMS chains.
    residues : list of Residue
        Every residue, classified.
    statistics : FitStatistics
        Fit quality.
    dimensions : list of int
        The box.
    spacing : float
        Lattice spacing (A per site).

    Returns
    -------
    list of str
        The warnings.
    """
    warnings = []
    beads = [residue for segment in segments for residue in segment.residues]

    kinds_of_code: dict[str, set[str]] = {}
    for residue in beads:
        kinds_of_code.setdefault(residue.code, set()).add(residue.kind)
    for code in sorted(kinds_of_code):
        kinds = kinds_of_code[code]
        if len(kinds) > 1:
            names = sorted({residue.resname for residue in beads if residue.code == code})
            warnings.append(
                "Bead type %s stands for %s (residues %s), so they share every interaction in the "
                "parameter file. Give them separate types with --map if they should differ."
                % (code, ' and '.join(sorted(kinds)), ', '.join(names)))

    unknown = sorted({residue.resname for residue in beads if residue.parent_source == 'unknown'})
    if unknown:
        warnings.append(
            "Residue name(s) %s are not standard and the file gives no parent for them; they have "
            "a polymer backbone, so they were converted as the unknown amino acid %s or the unknown "
            "nucleotide %s. Use --map to give them a bead type of their own."
            % (', '.join(unknown), UNKNOWN_AMINO_ACID_CODE, UNKNOWN_NUCLEOTIDE_CODE))

    # the file cannot say whether a repeated atom name starts a residue that shares
    # a number with the one before, or belongs to a molecule without unique names
    extra_beads = Counter(residue.label for residue in beads if residue.split_by_repeated_name)
    if extra_beads:
        listed = ['%s (%d extra bead(s))' % (label, count) for label, count in extra_beads.items()]
        warnings.append(
            "An atom name repeats within residue(s) %s, so each was read as several consecutive residues with "
            "the same chain, number and insertion code, one bead each. That is right for residues that share a "
            "number, but not for one molecule whose atoms are not uniquely named (Open Babel, for one, names a "
            "ligand's atoms by element alone): give such a molecule unique atom names." % _first_items(listed))

    if statistics.rmsd > POOR_FIT_RMSD_SPACINGS * spacing:
        warnings.append(
            "The lattice reproduces the structure less well than it could (RMSD %.2f A, against the "
            "%.2f A that rounding to a %.2f A lattice costs on its own). Bonds longer than a lattice step "
            "can span (see the bond lengths below) mean the spacing is too small for them; residues packed "
            "more densely than one per site mean it is too large." % (statistics.rmsd, 0.5 * spacing, spacing))

    longest_edge = max(dimensions) * spacing
    if longest_edge >= MAX_BOX_EDGE_ANGSTROMS:
        warnings.append(
            "The longest box axis is %.0f A at this spacing. PIMMS refuses a box whose longest axis "
            "times LATTICE_TO_ANGSTROMS reaches %.0f A (the PDB coordinate field), so the run will "
            "need a smaller LATTICE_TO_ANGSTROMS than the spacing used here."
            % (longest_edge, MAX_BOX_EDGE_ANGSTROMS))

    if not any(residue.kind in POLYMER_KINDS for residue in residues):
        warnings.append("The model has no amino acids or nucleotides; only hetero groups were converted.")
    return warnings


def convert_structure(path: str | Path, options: ConversionOptions) -> Conversion:
    """
    Convert a structure file into a PIMMS lattice configuration (without writing anything).

    Parameters
    ----------
    path : str or Path
        The PDB or mmCIF file.
    options : ConversionOptions
        What to convert and how.

    Returns
    -------
    Conversion
        Everything the restart file and the report need.

    Raises
    ------
    StructureFileException
        If the file cannot be read.
    StructureConversionException
        If the options cannot be applied or nothing can be converted.
    """
    structure = read_structure(path)
    model = _select_model(structure, options.model)
    residues = build_residues([atom for atom in structure.atoms if atom.model == model])
    residues, unselected = _select_chains(residues, options.chains)

    for residue in residues:
        classify_residue(residue, structure.parent_residues)
    demoted = demote_isolated_hetatm(residues)
    warnings = apply_bead_codes(residues, options.bead_map, options.keep_hetero)

    for residue in residues:
        if residue.kind in POLYMER_KINDS or (residue.kind == KIND_HETERO and options.keep_hetero):
            representative_position(residue, options.bead_position)

    # every polymer chain first, in file order, then one single-bead chain per kept hetero group
    segments, breaks = split_into_segments(residues)
    if options.keep_hetero:
        segments += [Segment([residue], residue.chain) for residue in residues if residue.kind == KIND_HETERO]
    if not segments:
        hint = '' if options.keep_hetero else ' (use --hetero to convert its hetero groups)'
        raise StructureConversionException(
            "Model %d of %s holds no amino acids or nucleotides%s, so there is nothing to convert"
            % (model, structure.path, hint))
    check_mobile_selections(segments, options.mobile)

    every_position = np.array([residue.position for segment in segments for residue in segment.residues])
    origin = every_position.min(axis=0)
    targets = []
    for segment in segments:
        positions = np.array([residue.position for residue in segment.residues], dtype=np.float64)
        targets.append((positions - origin) / options.spacing)
    if np.max(np.abs(np.concatenate(targets))) >= SITE_KEY_OFFSET // 2:
        raise StructureConversionException(
            "The structure spans more than %d lattice sites at a spacing of %g A; use a larger --spacing"
            % (SITE_KEY_OFFSET // 2, options.spacing))

    walks, methods = fit_to_lattice(targets)
    sites, dimensions = place_in_box(walks, options.padding, options.box)
    statistics = fit_statistics(segments, sites, options.spacing)
    warnings.extend(_conversion_warnings(segments, residues, statistics, dimensions, options.spacing))

    fallbacks = [residue.label for segment in segments for residue in segment.residues
                 if residue.position_source.startswith('heavy-atom centroid (no')]
    return Conversion(
        structure=structure, options=options, model=model, residues=residues, segments=segments,
        breaks=breaks, chain_types=assign_chain_types(segments), sites=sites, dimensions=dimensions,
        methods=methods, statistics=statistics,
        excluded=_summarise_exclusions(structure, model, residues, unselected, options.keep_hetero),
        position_fallbacks=fallbacks, isolated_hetatm=[residue.label for residue in demoted], warnings=warnings,
        input_sha256=_file_sha256(structure.path))


# .............................................................................
# Residues left free by --mobile
# .............................................................................

def check_mobile_selections(segments: list[Segment], selections: list[MobileSelection]) -> None:
    """
    Refuse --mobile selections that select nothing, or everything.

    A selection with a typo (a wrong chain, a range outside the chain) would
    otherwise leave the whole structure frozen without a word, and one that
    frees every residue would leave a freeze file that freezes nothing.

    Parameters
    ----------
    segments : list of Segment
        The PIMMS chains of the conversion.
    selections : list of MobileSelection
        The --mobile selections.

    Returns
    -------
    None

    Raises
    ------
    StructureConversionException
        If a selection matches no converted residue, or the selections together
        match every one.
    """
    if not selections:
        return
    residues = [residue for segment in segments for residue in segment.residues]
    for selection in selections:
        if not any(selection.matches(residue) for residue in residues):
            chains = sorted(set(residue.chain or BLANK_CHAIN_LABEL for residue in residues))
            raise StructureConversionException(
                "--mobile %s matches no converted residue (the converted chains are %s; residue numbers "
                "are the structure's own)" % (selection.text, ', '.join(chains)))
    if all(any(selection.matches(residue) for selection in selections) for residue in residues):
        raise StructureConversionException(
            "--mobile leaves every converted residue free, so the freeze file would freeze nothing; "
            "use --no-freeze instead")


def free_bead_flags(conversion: Conversion) -> list[list[bool]]:
    """
    Which beads of each chain the freeze file leaves free (--mobile).

    Parameters
    ----------
    conversion : Conversion
        The conversion; its options hold the selections.

    Returns
    -------
    list of list of bool
        One list per chain, in chainID order, with one entry per bead: True if
        a --mobile selection holds the bead's residue. All False without
        --mobile.
    """
    selections = conversion.options.mobile
    return [[any(selection.matches(residue) for selection in selections) for residue in segment.residues]
            for segment in conversion.segments]


def has_partly_frozen_chain(flags: list[list[bool]]) -> bool:
    """
    Whether --mobile leaves some chain partly frozen (some beads free, some not).

    Only such a chain gets a ``B`` line, and only its free beads are restricted to
    the crankshaft; a chain freed whole is an ordinary mobile chain.

    Parameters
    ----------
    flags : list of list of bool
        Whether each bead of each chain is free (see :func:`free_bead_flags`).

    Returns
    -------
    bool
        True if any chain has both free and frozen beads.
    """
    return any(any(free) and not all(free) for free in flags)


def _residue_run_label(residues: list[Residue], first: int, last: int) -> str:
    """
    Name a run of residues, e.g. ``A:SER1-LYS37``.

    Parameters
    ----------
    residues : list of Residue
        The residues of one chain.
    first, last : int
        Bead indices (inclusive) of the run.

    Returns
    -------
    str
        The first residue's label, and the last one's name and number if it is
        a different residue.
    """
    start = residues[first]
    if first == last:
        return start.label
    end = residues[last]
    return '%s-%s%s%s' % (start.label, end.resname, end.resseq, end.icode)


def free_residue_text(segment: Segment, free: list[bool]) -> str:
    """
    Describe the free beads of a chain, e.g. ``beads 0-36 (A:SER1-LYS37)``.

    Parameters
    ----------
    segment : Segment
        The chain.
    free : list of bool
        Whether each of its beads is free.

    Returns
    -------
    str
        Each run of free beads, by bead index (from 0) and by residue.
    """
    runs = contiguous_runs(index for index, is_free in enumerate(free) if is_free)
    return ', '.join('%s %s (%s)' % ('bead' if first == last else 'beads',
                                     format_bead_runs(range(first, last + 1)),
                                     _residue_run_label(segment.residues, first, last))
                     for first, last in runs)


# .............................................................................
# Writing the restart file
# .............................................................................

def restart_payload(conversion: Conversion) -> dict[str, object]:
    """
    Build the dictionary a PIMMS restart file holds.

    Only the four configuration keys are written. With no continuation state
    (step, generator states, temperature) the file is a starting
    configuration, which is what a structure is; RESTART_CONTINUE refuses it.

    Parameters
    ----------
    conversion : Conversion
        The conversion.

    Returns
    -------
    dict
        CHAINS, DIMENSIONS, ENERGY and HARDWALL, as plain Python objects.
    """
    chains: dict[int, list[object]] = {}
    for chain_id, (segment, walk, chain_type) in enumerate(
            zip(conversion.segments, conversion.sites, conversion.chain_types), start=1):
        chains[chain_id] = [[list(site) for site in walk], segment.sequence, chain_type]
    return {
        'CHAINS': chains,
        'DIMENSIONS': list(conversion.dimensions),
        'ENERGY': PLACEHOLDER_ENERGY,
        'HARDWALL': True,
    }


def _comment_text(text: str) -> str:
    """
    Make text from the input safe to put in a one-line comment.

    Residue and chain names come from the input file, and an mmCIF text field
    can hold a line break, which would end the comment and leave the rest of
    the name as a line PIMMS reads. Every run of whitespace (line breaks
    included) becomes one space.

    Parameters
    ----------
    text : str
        The text.

    Returns
    -------
    str
        The text on one line.
    """
    return ' '.join(text.split())


def freeze_file_text(conversion: Conversion, restart_path: Path) -> str:
    """
    Write a PIMMS freeze file that holds every chain of the restart file fixed.

    There is one ``C`` line per chain, with a comment naming the residues it
    came from, so that a chain can be released by deleting its line. The
    chainIDs are those of the restart file (1 to N, in the order of the
    report's chain table); chains added with EXTRA_CHAIN are numbered after
    them, so the file stays correct when mobile chains are added.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    restart_path : Path
        The restart file the chainIDs refer to.

    Returns
    -------
    str
        The freeze file.
    """
    n_chains = len(conversion.segments)
    flags = free_bead_flags(conversion)
    if any(any(free) for free in flags):
        return _partial_freeze_file_text(conversion, restart_path, flags)
    lines = ['# PIMMS freeze file written by pdb2pimms. It freezes every chain of %s,'
             % _comment_text(restart_path.name),
             '# converted from %s. PIMMS will not run with every chain frozen, so a run that'
             % _comment_text(conversion.structure.path.name),
             '# uses this file needs mobile chains as well: EXTRA_CHAIN adds them, numbered',
             '# from %d. Delete a line to let that chain move.' % (n_chains + 1)]
    width = len(str(n_chains))
    for chain_id, segment in enumerate(conversion.segments, start=1):
        lines.append('C %-*d  # %s, %d bead(s)'
                     % (width, chain_id, _comment_text(segment.span), len(segment.residues)))
    return '\n'.join(lines) + '\n'


def _partial_freeze_file_text(conversion: Conversion, restart_path: Path, flags: list[list[bool]]) -> str:
    """
    Write a freeze file that leaves the --mobile residues free.

    A chain with no free residue gets a ``C`` line, as in the default file. A
    chain with some gets a ``B`` line listing its frozen beads by index along
    the chain, counted from 0, and a chain whose residues are all free gets only
    a comment.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    restart_path : Path
        The restart file the chainIDs refer to.
    flags : list of list of bool
        Whether each bead of each chain is free (see :func:`free_bead_flags`).

    Returns
    -------
    str
        The freeze file.
    """
    n_chains = len(conversion.segments)
    selections = ' '.join(selection.text for selection in conversion.options.mobile)
    # the --mobile selections can be a long list, so their line is wrapped (never
    # splitting a selection such as C:118-128 at its hyphen); the rest of the
    # header keeps its layout
    selection_lines = textwrap.wrap('(%s), which are left free to move. A B line freezes the listed beads'
                                    % _comment_text(selections), width=FREEZE_HEADER_WIDTH, initial_indent='# ',
                                    subsequent_indent='# ', break_on_hyphens=False, break_long_words=False)
    lines = ['# PIMMS freeze file written by pdb2pimms. It freezes the chains of %s,'
             % _comment_text(restart_path.name),
             '# converted from %s, except the residues given with --mobile'
             % _comment_text(conversion.structure.path.name),
             *selection_lines,
             '# of one chain, by their index along it counted from 0. The free beads of',
             '# such a chain are moved by the crankshaft (MOVE_CRANKSHAFT) only. Chains',
             '# added with EXTRA_CHAIN are numbered from %d.' % (n_chains + 1)]
    width = len(str(n_chains))
    for chain_id, (segment, free) in enumerate(zip(conversion.segments, flags), start=1):
        n_beads = len(segment.residues)
        frozen = [index for index, is_free in enumerate(free) if not is_free]
        span = _comment_text(segment.span)
        if not frozen:
            lines.append('# chain %d (%s, %d bead(s)) is left free' % (chain_id, span, n_beads))
        elif len(frozen) == n_beads:
            lines.append('C %-*d  # %s, %d bead(s)' % (width, chain_id, span, n_beads))
        else:
            lines.append('B %-*d %s  # %s: %d of %d bead(s) frozen; free: %s'
                         % (width, chain_id, format_bead_runs(frozen, ' '), span, len(frozen), n_beads,
                            _comment_text(free_residue_text(segment, free))))
    return '\n'.join(lines) + '\n'


def _temporary_path(path: Path) -> Path:
    """
    A temporary name next to an output file, unique to this process.

    Parameters
    ----------
    path : Path
        The output file.

    Returns
    -------
    Path
        ``<name>.tmp.<pid>`` in the same directory.
    """
    return path.with_name('%s.tmp.%d' % (path.name, os.getpid()))


def _validate_with_pimms(path: Path) -> None:
    """
    Read a restart file back with PIMMS's own reader, which checks connectivity, overlaps and the box.

    Parameters
    ----------
    path : Path
        The restart file.

    Returns
    -------
    None

    Raises
    ------
    StructureConversionException
        If PIMMS refuses the file (which would be a bug in this tool).
    """
    from .restart import RestartObject

    try:
        RestartObject().build_from_file(str(path))  # type: ignore[no-untyped-call]
    except RestartException as error:
        raise StructureConversionException(
            "PIMMS refused the restart file pdb2pimms built, which is a bug in pdb2pimms; "
            "please report it with the input structure. PIMMS said: %s" % error)


def _validate_freeze_file(path: Path, conversion: Conversion) -> None:
    """
    Read a freeze file back with PIMMS's own reader and check it freezes exactly what it should.

    Without --mobile that is every chain, 1 to N, on ``C`` lines. With it, PIMMS
    must work out from the file the same frozen chains and frozen beads that
    :func:`free_bead_flags` gives.

    Parameters
    ----------
    path : Path
        The freeze file.
    conversion : Conversion
        The conversion it was written for.

    Returns
    -------
    None

    Raises
    ------
    StructureConversionException
        If PIMMS refuses the file or reads different frozen chains or beads
        from it (either would be a bug in this tool).
    """
    from .data_structures import FreezeFile

    n_chains = len(conversion.segments)
    lengths = {chain_id: len(segment.residues) for chain_id, segment in enumerate(conversion.segments, start=1)}
    try:
        freeze = FreezeFile(str(path))
        (frozen_chains, frozen_beads) = freeze.resolve(lengths)
    except KeyFileException as error:
        raise StructureConversionException(
            "PIMMS refused the freeze file pdb2pimms built, which is a bug in pdb2pimms; "
            "please report it with the input structure. PIMMS said: %s" % error)

    flags = free_bead_flags(conversion)
    if not any(any(free) for free in flags):
        frozen = list(freeze.chains)
        if frozen != list(range(1, n_chains + 1)) or freeze.beads:
            raise StructureConversionException(
                "PIMMS reads chainIDs %s from the freeze file pdb2pimms built, not 1 to %d, which is a bug in "
                "pdb2pimms; please report it with the input structure" % (frozen, n_chains))
        return

    expected_beads = {chain_id: tuple(index for index, is_free in enumerate(free) if not is_free)
                      for chain_id, free in enumerate(flags, start=1) if any(free) and not all(free)}
    expected_chains = [chain_id for chain_id, free in enumerate(flags, start=1) if not all(free)]
    if (frozen_chains, frozen_beads) != (expected_chains, expected_beads):
        raise StructureConversionException(
            "PIMMS reads frozen chains %s and frozen beads %s from the freeze file pdb2pimms built, not %s and "
            "%s, which is a bug in pdb2pimms; please report it with the input structure"
            % (frozen_chains, frozen_beads, expected_chains, expected_beads))


def write_outputs(conversion: Conversion, restart_path: Path, report_text: str, report_path: Path,
                  freeze_path: Optional[Path] = None) -> None:
    """
    Write the restart file, the report and the freeze file, each only once it is complete.

    All of them are written under temporary names. The restart file and the
    freeze file are read back with PIMMS's own readers, and only then is
    anything moved into place, so a failure while writing or checking them
    leaves no output behind (the final renames are separate operations, so a
    filesystem error between them would leave the restart file without the
    others).

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    restart_path : Path
        Where the restart file goes.
    report_text : str
        The report.
    report_path : Path
        Where the report goes.
    freeze_path : Path or None
        Where the freeze file goes; None writes no freeze file.

    Returns
    -------
    None

    Raises
    ------
    StructureConversionException
        If PIMMS refuses the restart file or the freeze file.
    OSError
        If a file cannot be written.
    """
    # (temporary file, final name) in the order they are moved into place
    renames = [(_temporary_path(restart_path), restart_path)]
    if freeze_path is not None:
        renames.append((_temporary_path(freeze_path), freeze_path))
    renames.append((_temporary_path(report_path), report_path))
    restart_temporary = renames[0][0]
    report_temporary = renames[-1][0]
    try:
        with open(restart_temporary, 'wb') as handle:
            pickle.dump(restart_payload(conversion), handle)
        _validate_with_pimms(restart_temporary)
        if freeze_path is not None:
            freeze_temporary = renames[1][0]
            with open(freeze_temporary, 'w', encoding='utf-8') as handle:
                handle.write(freeze_file_text(conversion, restart_path))
            _validate_freeze_file(freeze_temporary, conversion)
        with open(report_temporary, 'w', encoding='utf-8') as handle:
            handle.write(report_text)
        for temporary, final in renames:
            os.replace(temporary, final)
    finally:
        for temporary, _ in renames:
            if temporary.exists():
                temporary.unlink()


# .............................................................................
# The report
# .............................................................................

def _markdown_table(headers: list[str], rows: list[list[str]]) -> list[str]:
    """
    Lay out a Markdown table.

    Parameters
    ----------
    headers : list of str
        Column headings.
    rows : list of list of str
        The cells.

    Returns
    -------
    list of str
        The table's lines.
    """
    lines = ['| ' + ' | '.join(headers) + ' |', '|' + '|'.join('---' for _ in headers) + '|']
    for row in rows:
        lines.append('| ' + ' | '.join(cell.replace('|', '\\|') for cell in row) + ' |')
    return lines


def bead_type_table(conversion: Conversion) -> list[tuple[str, str, list[str], int]]:
    """
    Summarise the bead types the converted system uses.

    Parameters
    ----------
    conversion : Conversion
        The conversion.

    Returns
    -------
    list of tuple of (str, str, list of str, int)
        For each bead character, in the order used by the report (amino acids,
        nucleotides, then hetero groups, alphabetically within each): the
        character, the kind(s) it stands for, a description of each residue name
        it covers, and how many beads have it.
    """
    beads = [residue for segment in conversion.segments for residue in segment.residues]
    by_code: dict[str, list[Residue]] = {}
    for residue in beads:
        by_code.setdefault(residue.code, []).append(residue)

    kind_order = {KIND_AMINO_ACID: 0, KIND_RNA: 1, KIND_DNA: 1, KIND_HETERO: 2}
    table = []
    for code, members in by_code.items():
        kinds = sorted({member.kind for member in members}, key=lambda kind: kind_order.get(kind, 3))
        names: Counter[str] = Counter(member.resname for member in members)
        descriptions = []
        for name in sorted(names):
            example = next(member for member in members if member.resname == name)
            description = '%s (%d)' % (name, names[name])
            if example.parent_source == 'file':
                description += ', modified %s (parent from the file)' % example.parent
            elif example.parent_source == 'built-in table':
                description += ', modified %s (parent from the built-in table)' % example.parent
            elif example.parent_source == 'other name':
                description += ', another name for %s' % example.parent
            elif example.parent_source == 'unknown':
                description += ', not a known residue'
            if name in conversion.options.bead_map:
                description += ', set with --map'
            descriptions.append(description)
        table.append((code, ' and '.join(kinds), descriptions, len(members)))
    table.sort(key=lambda row: (min(kind_order.get(kind, 3) for kind in row[1].split(' and ')), row[0]))
    return table


def _length_summary(lengths: FloatArray) -> str:
    """
    Describe a set of bond lengths in one phrase.

    Parameters
    ----------
    lengths : ndarray
        Bond lengths (A).

    Returns
    -------
    str
        Mean, standard deviation and range, or 'none' for an empty set.
    """
    if lengths.size == 0:
        return 'none'
    return '%.2f +/- %.2f A (range %.2f to %.2f A)' % (
        float(lengths.mean()), float(lengths.std()), float(lengths.min()), float(lengths.max()))


def _box_text(dimensions: list[int], separator: str = ' ') -> str:
    """
    Write a box as its site counts.

    Parameters
    ----------
    dimensions : list of int
        The box.
    separator : str, optional
        What goes between the axes. Default is a space (the keyfile form).

    Returns
    -------
    str
        e.g. ``14 14 16`` or ``14 x 14 x 16``.
    """
    return separator.join(str(edge) for edge in dimensions)


def _first_items(items: list[str]) -> str:
    """
    List at most REPORT_LIST_LIMIT items, saying how many more there are.

    Parameters
    ----------
    items : list of str
        The items.

    Returns
    -------
    str
        The items, comma-separated, with ``and N more`` when some were left out.
    """
    shown = ', '.join(items[:REPORT_LIST_LIMIT])
    if len(items) > REPORT_LIST_LIMIT:
        shown += ' and %d more' % (len(items) - REPORT_LIST_LIMIT)
    return shown


def _summary_section(conversion: Conversion, codes: list[str]) -> list[str]:
    """
    The report's summary table.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    codes : list of str
        The bead types, in report order.

    Returns
    -------
    list of str
        The section's lines.
    """
    options = conversion.options
    statistics = conversion.statistics
    structure = conversion.structure
    n_beads = sum(len(segment.residues) for segment in conversion.segments)
    n_sites = int(np.prod(conversion.dimensions))
    n_chain_types = len(set(conversion.chain_types))
    entry = ', entry %s' % structure.entry_id if structure.entry_id else ''
    box_angstroms = ' x '.join('%.1f' % (edge * options.spacing) for edge in conversion.dimensions)
    if options.bead_position == BEAD_POSITION_BACKBONE:
        position = "CA (amino acids) and C4' (nucleotides)"
    else:
        position = 'heavy-atom centroid of every residue'
    rows = [
        ['Input file', '`%s` (%s%s)' % (structure.path.resolve(), structure.file_format, entry)],
        ['Input SHA-256', '`%s`' % conversion.input_sha256],
        ['Model converted', '%d (of %d in the file)' % (conversion.model, len(structure.models))],
        ['PIMMS chains', '%d (%d chain type%s)' % (len(conversion.segments), n_chain_types,
                                                  '' if n_chain_types == 1 else 's')],
        ['Beads', str(n_beads)],
        ['Bead types', '%d: `%s`' % (len(codes), ' '.join(codes))],
        ['Box (DIMENSIONS)', '%s lattice sites (%s A)' % (_box_text(conversion.dimensions), box_angstroms)],
        ['Lattice spacing', '%g A per site' % options.spacing],
        ['Bead position', position],
        ['Fit RMSD', '%.2f A (largest single-bead deviation %.2f A, at %s)'
         % (statistics.rmsd, statistics.max_deviation, statistics.max_deviation_label)],
        ['Occupied sites', '%d of %d (%.2f %%)' % (n_beads, n_sites, 100.0 * n_beads / n_sites)],
    ]
    return ['## Summary', ''] + _markdown_table(['', ''], rows) + ['']


def _bead_type_section(bead_types: list[tuple[str, str, list[str], int]]) -> list[str]:
    """
    The report's list of the bead types the parameter file must define.

    Parameters
    ----------
    bead_types : list of tuple
        The rows of bead_type_table().

    Returns
    -------
    list of str
        The section's lines.
    """
    codes = [row[0] for row in bead_types]
    n_types = len(codes)
    text = ('The restart file uses these %d bead types: `%s`. A PIMMS parameter file is case-sensitive, so '
            'each must be defined exactly as written here. Every bead type needs a solvation line (an '
            'interaction with the solvent type `0`), an interaction line with every bead type including '
            'itself, and an `ANGLE_PENALTY` or `ANGLE_PENALTY_T_NORM` line unless the keyfile sets '
            '`ANGLES_OFF : True`. For %d types that is %d solvation lines, %d pair lines and %d angle lines.'
            % (n_types, ' '.join(codes), n_types, n_types, n_types * (n_types + 1) // 2, n_types))
    rows = [['`%s`' % code, kinds, '; '.join(descriptions), str(count)]
            for code, kinds, descriptions, count in bead_types]
    return (['## Bead types the parameter file must define', '', text, '']
            + _markdown_table(['Bead', 'Kind', 'Residue names (count)', 'Beads'], rows) + [''])


def _freeze_paragraph(conversion: Conversion, freeze_path: Path) -> str:
    """
    The report's explanation of the freeze file.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    freeze_path : Path
        The freeze file.

    Returns
    -------
    str
        The paragraph.
    """
    n_chains = len(conversion.segments)
    if any(any(free) for free in free_bead_flags(conversion)):
        return _partial_freeze_paragraph(conversion, freeze_path)
    return ('The freeze file `%s` lists every chain of the restart file (chainID%s), one per line, so the '
            'converted structure stays exactly as it is while the rest of the system moves around it. Frozen '
            'chains still occupy their sites and still contribute to the energy, so the parameter file must '
            'still define their bead types. PIMMS will not run with every chain frozen, so a run that uses the '
            'freeze file needs mobile chains as well: `EXTRA_CHAIN` adds them, and numbers them from %d, so the '
            'freeze file stays correct. `EXTRA_CHAIN` sequences are upper-cased unless the keyfile sets '
            '`CASE_INSENSITIVE_CHAINS : False`, which an added chain with lower-case bead types (nucleotides) '
            'needs. Delete a chain\'s line from the freeze file to let it move, or leave out `FREEZE_FILE` to '
            'let everything move.'
            % (freeze_path.name, ' 1' if n_chains == 1 else 's 1 to %d' % n_chains, n_chains + 1))


def _partial_freeze_paragraph(conversion: Conversion, freeze_path: Path) -> str:
    """
    The report's explanation of a freeze file that leaves the --mobile residues free.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    freeze_path : Path
        The freeze file.

    Returns
    -------
    str
        The paragraph.
    """
    flags = free_bead_flags(conversion)
    n_free = sum(sum(free) for free in flags)
    n_beads = sum(len(free) for free in flags)
    selections = ', '.join('`%s`' % selection.text for selection in conversion.options.mobile)
    if has_partly_frozen_chain(flags):
        moves = ('Only the crankshaft move (`MOVE_CRANKSHAFT`) moves the free beads of a partly frozen chain - '
                 'every move that displaces a whole chain or a segment of one leaves such a chain alone - so a '
                 'run that uses this file needs `MOVE_CRANKSHAFT`. ')
    else:
        moves = ''
    if any(all(free) for free in flags):
        moves += 'A chain whose residues are all free is left out of the file and moves like any other chain. '
    return ('The freeze file `%s` freezes the converted structure except the residues given with `--mobile` '
            '(%s): %d of its %d beads are free, listed below. A chain with no free residue has a `C` line; a '
            'chain with some has a `B` line that lists its frozen beads by their index along the chain, '
            'counted from 0. %sFrozen beads still occupy their sites and still contribute to the energy, so '
            'the parameter file must still define their bead types. Mobile chains added with `EXTRA_CHAIN` are '
            'numbered from %d; they are optional here, because the free residues can move on their own. '
            '`EXTRA_CHAIN` sequences are upper-cased unless the keyfile sets `CASE_INSENSITIVE_CHAINS : False`, '
            'which an added chain with lower-case bead types (nucleotides) needs.'
            % (freeze_path.name, selections, n_free, n_beads, moves, len(conversion.segments) + 1))


def _free_residue_lines(conversion: Conversion) -> list[str]:
    """
    The report's list of the beads the freeze file leaves free, one line per chain that has any.

    Parameters
    ----------
    conversion : Conversion
        The conversion.

    Returns
    -------
    list of str
        Markdown list lines.
    """
    lines = []
    for chain_id, (segment, free) in enumerate(zip(conversion.segments, free_bead_flags(conversion)), start=1):
        if not any(free):
            continue
        status = 'every bead free' if all(free) else '%d of %d beads free' % (sum(free), len(free))
        lines.append('- chain %d (%s): %s - %s' % (chain_id, segment.span, status, free_residue_text(segment, free)))
    return lines


def _usage_section(conversion: Conversion, restart_path: Path, n_types: int,
                   freeze_path: Optional[Path] = None) -> list[str]:
    """
    The report's instructions for using the restart file (and the freeze file) in a run.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    restart_path : Path
        The restart file.
    n_types : int
        The number of bead types.
    freeze_path : Path or None
        The freeze file, or None if none was written.

    Returns
    -------
    list of str
        The section's lines.
    """
    n_sites = int(np.prod(conversion.dimensions))
    keyfile = ['```',
               'RESTART_FILE         : %s' % restart_path.name,
               'DIMENSIONS           : %s' % _box_text(conversion.dimensions),
               'LATTICE_TO_ANGSTROMS : %g' % conversion.options.spacing,
               'PARAMETER_FILE       : <a file defining the %d bead types above>' % n_types]
    flags = free_bead_flags(conversion)
    partial = any(any(free) for free in flags)
    if freeze_path is not None and partial:
        keyfile += ['FREEZE_FILE          : %s' % freeze_path.name]
        if has_partly_frozen_chain(flags):
            keyfile += ['MOVE_CRANKSHAFT      : <fraction>   # the only move for a partly frozen chain\'s free beads']
        keyfile += ['EXTRA_CHAIN          : <count> <sequence>   # mobile chains (optional here)']
    elif freeze_path is not None:
        keyfile += ['FREEZE_FILE          : %s' % freeze_path.name,
                    'EXTRA_CHAIN          : <count> <sequence>   # mobile chains, needed with FREEZE_FILE']
    keyfile += ['```']
    boundary = ('The file was written with `HARDWALL : True`, because no chain crosses a face of the box. '
                'That makes it valid for a run with either boundary condition, and lets `DIMENSIONS` be any box '
                'at least %s: PIMMS centres this box inside a larger one. `LATTICE_TO_ANGSTROMS` only scales '
                '`START.pdb` and `traj.xtc`; setting it to the spacing used here keeps them on the scale of the '
                'input structure.' % _box_text(conversion.dimensions, ' x '))
    contents = ('The restart file holds the configuration only. Its `ENERGY` is a placeholder (%g), which no '
                'run reads: PIMMS computes the energy from the parameter file at start-up. Because it holds no '
                'step, random-number state or temperature, it starts a new run and `RESTART_CONTINUE` refuses '
                'it. The run box needs %d bytes per site for the lattice, so this box needs %.1f MB before '
                'anything else.' % (PLACEHOLDER_ENERGY, LATTICE_BYTES_PER_SITE,
                                    n_sites * LATTICE_BYTES_PER_SITE / 1e6))
    lines = (['## Using the restart file', '',
              'Give the restart file to a run with `RESTART_FILE`. The chains come from it, so no `CHAIN` lines '
              'are needed (any are ignored). The lines that matter are:', '']
             + keyfile + ['', boundary, '', contents, ''])
    if freeze_path is not None:
        lines += [_freeze_paragraph(conversion, freeze_path), '']
        if partial:
            lines += _free_residue_lines(conversion) + ['']
    return lines


def _chain_section(conversion: Conversion) -> list[str]:
    """
    The report's table of PIMMS chains and their sequences.

    Parameters
    ----------
    conversion : Conversion
        The conversion.

    Returns
    -------
    list of str
        The section's lines.
    """
    statistics = conversion.statistics
    rows = []
    for chain_id, (segment, chain_type, rmsd, rg_in, rg_out, method) in enumerate(zip(
            conversion.segments, conversion.chain_types, statistics.segment_rmsd, statistics.input_rg,
            statistics.lattice_rg, conversion.methods), start=1):
        rows.append([str(chain_id), str(chain_type), segment.chain or '-', segment.span,
                     str(len(segment.residues)), '%.2f' % rmsd, '%.1f / %.1f' % (rg_in, rg_out), method])
    lines = ['## Chains', '']
    lines += _markdown_table(['chainID', 'chainType', 'Source chain', 'Residues', 'Beads', 'RMSD (A)',
                              'Rg input / lattice (A)', 'Placed by'], rows)
    lines += ['', 'Sequences, one character per bead:', '', '```']
    for chain_id, segment in enumerate(conversion.segments, start=1):
        lines += ['>%d %s' % (chain_id, segment.span), segment.sequence]
    return lines + ['```', '']


def _break_section(breaks: list[ChainBreak]) -> list[str]:
    """
    The report's list of chain breaks (nothing when there are none).

    Parameters
    ----------
    breaks : list of ChainBreak
        The breaks.

    Returns
    -------
    list of str
        The section's lines.
    """
    if not breaks:
        return []
    lines = ['## Chain breaks', '',
             'Each source chain was split wherever consecutive residues are not bonded, usually because '
             'residues are missing from the model. PIMMS chains must be connected, so each bonded run is its '
             'own chain.', '']
    lines += ['- Chain %s, between %s and %s: %s' % (item.chain or '-', item.before, item.after, item.reason)
              for item in breaks]
    return lines + ['']


def _exclusion_section(conversion: Conversion) -> list[str]:
    """
    The report's list of what was read but not converted.

    Parameters
    ----------
    conversion : Conversion
        The conversion.

    Returns
    -------
    list of str
        The section's lines.
    """
    excluded = conversion.excluded
    items = []
    if excluded.waters:
        items.append('%d water molecules (water is never converted).' % excluded.waters)
    if excluded.hetero:
        counts = ', '.join('%s x %d' % (name, count) for name, count in sorted(excluded.hetero.items()))
        items.append('Hetero groups (ligands, ions, cofactors): %s. Run with `--hetero` to convert each into a '
                     'one-bead chain.' % counts)
    if conversion.isolated_hetatm:
        fate = ('converted as one-bead hetero chains' if conversion.options.keep_hetero
                else 'left out with the other hetero groups')
        items.append('%d HETATM residue(s) with an amino-acid or nucleotide name or backbone are bonded to neither '
                     'neighbour, so they are ligands rather than part of a chain, and were %s: %s.'
                     % (len(conversion.isolated_hetatm), fate, _first_items(conversion.isolated_hetatm)))
    if excluded.hydrogens:
        items.append('%d hydrogen atoms, which do not affect any bead position (except in a residue with no '
                     'heavy atom).' % excluded.hydrogens)
    if excluded.altloc_residues:
        kept = ['%s (kept %s of %s)' % (residue.label, residue.altloc, '/'.join(residue.altlocs_present))
                for residue in excluded.altloc_residues]
        items.append('Alternate locations: for each of %d residues only the location with the highest '
                     'occupancy was kept: %s.' % (len(kept), _first_items(kept)))
    if excluded.other_models:
        items.append('%d other model(s) in the file; this conversion used model %d (choose another with '
                     '`--model`).' % (excluded.other_models, conversion.model))
    if excluded.unselected_chains:
        items.append('Chains not selected with `--chains`: %s.' % ', '.join(excluded.unselected_chains))
    if not items:
        items.append('Nothing: every atom of the model went into a bead.')
    return ['## What was not converted', ''] + ['- %s' % item for item in items] + ['']


def _method_section(conversion: Conversion) -> list[str]:
    """
    The report's account of how the conversion was done, with the bond statistics.

    Parameters
    ----------
    conversion : Conversion
        The conversion.

    Returns
    -------
    list of str
        The section's lines.
    """
    spacing = conversion.options.spacing
    if conversion.options.bead_position == BEAD_POSITION_BACKBONE:
        positions = ("Each residue became one bead, placed at its CA atom (amino acids) or C4' atom "
                     "(nucleotides); hetero groups were placed at their heavy-atom centroid.")
    else:
        positions = 'Each residue became one bead, placed at the centroid of its heavy atoms.'
    if conversion.position_fallbacks:
        positions += (' %d residue(s) lacked that atom and were placed at their heavy-atom centroid: %s.'
                      % (len(conversion.position_fallbacks), _first_items(conversion.position_fallbacks)))
    bonds = ("Consecutive residues of a chain were treated as bonded when the atoms of the linking bond (C-N "
             "for a peptide, O3'-P for a phosphodiester) are within %.1f A, or, where those atoms are missing, "
             "when the CA atoms are within %.1f A or the C4' atoms within %.1f A. Hetero groups are never "
             "bonded to anything." % (MAX_LINK_BOND_ANGSTROMS, MAX_CA_CA_ANGSTROMS, MAX_C4_C4_ANGSTROMS))
    fitting = ('A PIMMS chain is a self-avoiding walk in which consecutive beads are Moore neighbours (they '
               'differ by at most one site along every axis). Each chain was placed by dynamic programming, '
               'which finds the walk with the smallest total squared distance between each bead and its scaled '
               'input position, among walks that keep every bead within a couple of sites of its position and '
               'avoid the chains already placed. Where that walk visited a site twice, the site was forbidden '
               'to the bead further from it and the walk found again. A greedy search with backtracking is the '
               'fallback when no such walk exists. Chains were placed in the order of the table above.')
    # rounding each coordinate to the nearest site leaves an error uniform over one
    # spacing, of variance spacing^2 / 12 per axis, so an RMSD of spacing / 2 in 3D
    deviations = ('Deviations are measured after the best rigid-body superposition of the lattice onto the '
                  'input. Rounding each position to the nearest site of a cubic lattice gives an RMSD of about '
                  'half a spacing (%.2f A here) on its own, so a fit close to that is as good as the lattice '
                  'allows.' % (0.5 * spacing))
    statistics = conversion.statistics
    steps = statistics.bond_steps
    rows = [
        ['Bonds', str(sum(steps.values()))],
        ['Along one axis (1 spacing)', str(steps.get(1, 0))],
        ['Along a face diagonal (1.41 spacings)', str(steps.get(2, 0))],
        ['Along a body diagonal (1.73 spacings)', str(steps.get(3, 0))],
        ['Bond length in the input', _length_summary(statistics.input_bond_lengths)],
        ['Bond length on the lattice', _length_summary(statistics.lattice_bond_lengths)],
    ]
    return (['## How the conversion was done', '', positions, '', bonds, '', fitting, '', deviations, '']
            + _markdown_table(['Bonds between consecutive beads', ''], rows) + [''])


def render_report(conversion: Conversion, restart_path: Path, report_path: Path, command: str,
                  version: str, timestamp: str, freeze_path: Optional[Path] = None) -> str:
    """
    Write the Markdown report of a conversion.

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    restart_path : Path
        The restart file it describes.
    report_path : Path
        Where the report will be written (so the restart file can be named relative to it).
    command : str
        The command line that ran the conversion.
    version : str
        The PIMMS version.
    timestamp : str
        When the conversion ran.
    freeze_path : Path or None
        The freeze file written with the restart file, or None if none was.

    Returns
    -------
    str
        The report.
    """
    bead_types = bead_type_table(conversion)
    codes = [row[0] for row in bead_types]
    try:
        restart_name = os.path.relpath(restart_path, report_path.parent)
    except ValueError:
        restart_name = str(restart_path)

    lines = ['# pdb2pimms conversion report', '',
             '`%s` was converted into the PIMMS restart file `%s` on %s with PIMMS %s.'
             % (conversion.structure.path.name, restart_name, timestamp, version),
             '', '```', command, '```', '']
    lines += _summary_section(conversion, codes)
    if conversion.warnings:
        lines += ['## Warnings', ''] + ['- %s' % warning for warning in conversion.warnings] + ['']
    lines += _bead_type_section(bead_types)
    lines += _usage_section(conversion, restart_path, len(codes), freeze_path)
    lines += _chain_section(conversion)
    lines += _break_section(conversion.breaks)
    lines += _exclusion_section(conversion)
    lines += _method_section(conversion)
    return '\n'.join(lines)


# .............................................................................
# Command line
# .............................................................................

def pimms_version() -> str:
    """
    The installed PIMMS version.

    Returns
    -------
    str
        The version, or 'unknown' if it cannot be read.
    """
    try:
        from . import __version__
        return str(__version__)
    except ImportError:
        try:
            from importlib.metadata import version
            return version('idptools-pimms')
        except Exception:
            return 'unknown'


def parse_bead_map(entries: list[str]) -> dict[str, str]:
    """
    Read --map entries of the form RESNAME=CHARACTER.

    Parameters
    ----------
    entries : list of str
        The --map values.

    Returns
    -------
    dict of str to str
        Residue name -> bead character.

    Raises
    ------
    ValueError
        If an entry is malformed, repeats a residue name, or names a character PIMMS cannot use.
    """
    bead_map: dict[str, str] = {}
    for entry in entries:
        name, separator, character = entry.partition('=')
        name = name.strip()
        if not separator or not name:
            raise ValueError("--map %s: expected RESNAME=CHARACTER, e.g. --map MSE=M" % entry)
        check_bead_character(character, '--map %s' % entry)
        if name in bead_map and bead_map[name] != character:
            raise ValueError("--map gives %s two bead types (%s and %s)" % (name, bead_map[name], character))
        bead_map[name] = character
    return bead_map


def parse_mobile_selection(text: str) -> MobileSelection:
    """
    Read a --mobile selection: ``CHAIN``, ``CHAIN:NUMBER`` or ``CHAIN:FIRST-LAST``.

    Parameters
    ----------
    text : str
        The --mobile value. CHAIN is an author chain ID (``-`` for a blank one);
        the numbers are author residue numbers and may be negative, e.g.
        ``I:-73--60``.

    Returns
    -------
    MobileSelection
        The selection.

    Raises
    ------
    ValueError
        If the value is not of one of those forms, or a range runs backwards.
    """
    chain, colon, numbers = text.strip().partition(':')
    chain = chain.strip()
    if not chain:
        raise ValueError("--mobile %s: expected CHAIN, CHAIN:NUMBER or CHAIN:FIRST-LAST, e.g. --mobile A:1-37"
                         % text)
    if chain == BLANK_CHAIN_LABEL:
        chain = ''
    if not colon:
        return MobileSelection(chain=chain, first=None, last=None, text=text.strip())
    match = MOBILE_RESIDUE_RANGE.fullmatch(numbers.strip())
    if match is None:
        raise ValueError("--mobile %s: expected residue numbers after the colon, as NUMBER or FIRST-LAST, "
                         "e.g. --mobile A:1-37" % text)
    first = int(match.group(1))
    last = first if match.group(2) is None else int(match.group(2))
    if last < first:
        raise ValueError("--mobile %s: the range runs backwards" % text)
    return MobileSelection(chain=chain, first=first, last=last, text=text.strip())


def build_argument_parser() -> argparse.ArgumentParser:
    """
    The pdb2pimms command-line interface.

    Returns
    -------
    argparse.ArgumentParser
        The parser.
    """
    parser = argparse.ArgumentParser(
        prog='pdb2pimms',
        description=(
            "Convert a PDB or mmCIF structure (optionally gzipped) into a PIMMS restart file, one bead\n"
            "per residue, and write a Markdown report of the conversion that lists every bead type the\n"
            "parameter file must define. Amino acids keep their one-letter codes and nucleotides get\n"
            "lower-case ones (a, c, g, u, t). A freeze file listing every chain is written as well, so\n"
            "the structure can be held fixed while other chains move around it."),
        epilog=(
            "examples:\n"
            "  pdb2pimms 1ubq.cif                     writes 1ubq_restart.pimms, 1ubq_report.md and\n"
            "                                         1ubq_freeze.txt\n"
            "  pdb2pimms 4hhb.pdb --chains A,B        convert two chains only\n"
            "  pdb2pimms 1ehz.cif --spacing 5         a coarser lattice\n"
            "  pdb2pimms model.pdb --map MSE=Z        give selenomethionine its own bead type\n"
            "  pdb2pimms 4hhb.pdb --hetero            keep ligands and ions as one-bead chains\n"
            "  pdb2pimms 1kx5.cif --mobile A:1-26 --mobile E:1-26\n"
            "                                         freeze everything except the two H3 tails"),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('structure', metavar='STRUCTURE',
                        help='The PDB (.pdb, .ent) or mmCIF (.cif, .mmcif) file, optionally ending in .gz')
    parser.add_argument('-o', '--output', metavar='RESTART_FILE',
                        help='The restart file to write (default: <structure name>%s in the current '
                             'directory)' % RESTART_SUFFIX)
    parser.add_argument('--report', metavar='REPORT_FILE',
                        help='The report to write (default: next to the restart file, named after it with '
                             'its ending replaced by %s)' % REPORT_SUFFIX)
    freeze = parser.add_mutually_exclusive_group()
    freeze.add_argument('--freeze', metavar='FREEZE_FILE',
                        help='The freeze file to write, which freezes every chain (default: next to the restart '
                             'file, named after it with its ending replaced by %s)' % FREEZE_SUFFIX)
    freeze.add_argument('--no-freeze', action='store_true', help='Do not write a freeze file')
    parser.add_argument('--mobile', action='append', default=[], metavar='SELECTION',
                        help="Leave these residues free to move in the freeze file: CHAIN (a whole chain), "
                             "CHAIN:NUMBER or CHAIN:FIRST-LAST, with the structure's own chain IDs ('-' for a "
                             "blank one) and residue numbers; repeat for more. A chain ID selects every converted "
                             "residue with it, hetero groups kept with --hetero included. The rest of a chain is "
                             "frozen bead by bead, and its free beads move by MOVE_CRANKSHAFT only. Write a "
                             "selection that starts with '-' with an equals sign (--mobile=-:1-5), or it is taken "
                             "for an option")
    parser.add_argument('--spacing', type=float, default=DEFAULT_SPACING_ANGSTROMS, metavar='ANGSTROMS',
                        help='Lattice spacing in Angstroms per site (default %g, the PIMMS default '
                             'LATTICE_TO_ANGSTROMS)' % DEFAULT_SPACING_ANGSTROMS)
    parser.add_argument('--padding', type=int, default=DEFAULT_PADDING_SITES, metavar='SITES',
                        help='Empty lattice sites between the structure and each face of the box (default %d)'
                             % DEFAULT_PADDING_SITES)
    parser.add_argument('--box', type=int, nargs=3, metavar=('X', 'Y', 'Z'),
                        help='Use this box (in lattice sites) instead, with the structure centred in it')
    parser.add_argument('--model', type=int, metavar='N',
                        help='The model to convert (default: the first model in the file)')
    parser.add_argument('--chains', metavar='IDS',
                        help="Comma-separated author chain IDs to convert, with '-' for a blank chain ID "
                             "(default: every chain; a list that starts with '-' needs an equals sign, "
                             "--chains=-,A)")
    parser.add_argument('--bead-position', choices=(BEAD_POSITION_BACKBONE, BEAD_POSITION_CENTROID),
                        default=BEAD_POSITION_BACKBONE,
                        help="Where each bead goes: 'backbone' puts amino acids at CA and nucleotides at "
                             "C4' (default); 'centroid' puts every residue at its heavy-atom centroid")
    parser.add_argument('--hetero', action='store_true',
                        help='Also convert every non-water hetero group (ligand, ion) into a one-bead chain')
    parser.add_argument('--map', action='append', default=[], metavar='RESNAME=CHARACTER',
                        help='Use CHARACTER as the bead type of residue RESNAME (repeat for more)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Replace the restart file, report and freeze file if they already exist')
    parser.add_argument('--version', action='version', version='version %s' % pimms_version(),
                        help='Print the installed PIMMS version and exit')
    return parser


def structure_stem(path: Path) -> str:
    """
    The input file's name without its structure and compression extensions.

    Parameters
    ----------
    path : Path
        The structure file.

    Returns
    -------
    str
        e.g. ``1ubq`` for ``1ubq.cif.gz``.
    """
    name = path.name
    changed = True
    while changed:
        changed = False
        for extension in STRUCTURE_EXTENSIONS:
            if name.lower().endswith(extension) and len(name) > len(extension):
                name = name[:-len(extension)]
                changed = True
    return name


def output_paths(structure: Path, output: Optional[str], report: Optional[str],
                 freeze: Optional[str] = None) -> tuple[Path, Path, Path]:
    """
    Decide where the restart file, the report and the freeze file go.

    Parameters
    ----------
    structure : Path
        The input file.
    output : str or None
        The --output value.
    report : str or None
        The --report value.
    freeze : str or None
        The --freeze value.

    Returns
    -------
    tuple of (Path, Path, Path)
        The restart file, the report and the freeze file. By default the
        restart file is ``<structure name>_restart.pimms`` in the current
        directory, and the report and the freeze file sit next to it, named
        after it (``x_restart.pimms`` or ``x.pimms`` -> ``x_report.md`` and
        ``x_freeze.txt``).
    """
    restart_path = Path(output) if output else Path.cwd() / (structure_stem(structure) + RESTART_SUFFIX)

    # the companions are named after the restart file rather than the input, so
    # that converting one structure into two restart files (-o a.pimms, -o
    # b.pimms) gives two of each. Names that differ only by the _restart ending
    # (x.pimms, x_restart.pimms) share them; the second conversion is then
    # refused unless --overwrite is given.
    base = restart_path.name
    for suffix in (RESTART_EXTENSION, RESTART_STEM_SUFFIX):
        if base.endswith(suffix) and len(base) > len(suffix):
            base = base[:-len(suffix)]
    report_path = Path(report) if report else restart_path.parent / (base + REPORT_SUFFIX)
    freeze_path = Path(freeze) if freeze else restart_path.parent / (base + FREEZE_SUFFIX)
    return restart_path, report_path, freeze_path


def _same_file_key(path: Path) -> str:
    """
    A key that two names of one file share, whether or not the file exists yet.

    The path is resolved (so ``./x`` and symbolic links count as their target),
    put in Unicode normal form D and case-folded, because the default file
    systems of macOS and Windows ignore case (there ``notes.md`` and
    ``NOTES.md`` are one file) and macOS's also ignores the difference between
    composed and decomposed accented letters. On other file systems this
    refuses two names that differ only in those ways, which is the safe side
    to err on.

    Parameters
    ----------
    path : Path
        A file name.

    Returns
    -------
    str
        The key.
    """
    # NFD first: macOS's file system also treats the composed and decomposed
    # spellings of a name (café written with U+00E9 or with e + U+0301) as one file
    return unicodedata.normalize('NFD', str(path.resolve())).casefold()


def _output_problem(outputs: dict[str, Path], overwrite: bool, structure: Optional[Path] = None) -> Optional[str]:
    """
    Check the output paths before any work is done.

    Parameters
    ----------
    outputs : dict of str to Path
        Each output file, keyed by what it is ('restart file', 'report', ...).
    overwrite : bool
        Whether existing files may be replaced.
    structure : Path or None
        The input structure file, which no output may replace (even with
        --overwrite); None skips that check.

    Returns
    -------
    str or None
        What is wrong, or None.
    """
    seen: dict[str, str] = {}
    if structure is not None:
        seen[_same_file_key(structure)] = 'input structure'
    for role, path in outputs.items():
        earlier = seen.setdefault(_same_file_key(path), role)
        if earlier != role:
            return 'the %s and the %s cannot be the same file (%s)' % (earlier, role, path)
    for path in outputs.values():
        if not path.parent.is_dir():
            return 'the directory %s does not exist' % path.parent
        if path.is_dir():
            return '%s is a directory' % path
        if path.exists() and not overwrite:
            return '%s already exists; give --overwrite to replace it' % path
    return None


def _print_summary(conversion: Conversion, restart_path: Path, report_path: Path,
                   freeze_path: Optional[Path]) -> None:
    """
    Print what was written, and any warnings (to standard error).

    Parameters
    ----------
    conversion : Conversion
        The conversion.
    restart_path, report_path : Path
        The files written.
    freeze_path : Path or None
        The freeze file written, or None if none was.

    Returns
    -------
    None
    """
    n_beads = sum(len(segment.residues) for segment in conversion.segments)
    codes = [row[0] for row in bead_type_table(conversion)]
    print('Wrote %s: %d chain(s), %d beads, box %s, fit RMSD %.2f A'
          % (restart_path, len(conversion.segments), n_beads,
             ' '.join(str(v) for v in conversion.dimensions), conversion.statistics.rmsd))
    print('Wrote %s' % report_path)
    flags = free_bead_flags(conversion)
    if freeze_path is not None and any(any(free) for free in flags):
        n_beads = sum(len(free) for free in flags)
        n_free = sum(sum(free) for free in flags)
        n_whole = sum(1 for free in flags if not any(free))
        n_partly = sum(1 for free in flags if any(free) and not all(free))
        n_mobile = sum(1 for free in flags if all(free))
        crankshaft = ('; the free beads of a partly frozen chain move by MOVE_CRANKSHAFT only'
                      if n_partly else '')
        print('Wrote %s: freezes %d of %d beads (%d chain(s) entirely, %d partly, %d left free)%s'
              % (freeze_path, n_beads - n_free, n_beads, n_whole, n_partly, n_mobile, crankshaft))
    elif freeze_path is not None:
        print('Wrote %s: freezes all %d chain(s); a run using it needs mobile chains too (EXTRA_CHAIN)'
              % (freeze_path, len(conversion.segments)))
    print('Bead types the parameter file must define (%d): %s' % (len(codes), ' '.join(codes)))
    for warning in conversion.warnings:
        sys.stderr.write('WARNING: %s\n' % warning)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Run pdb2pimms.

    Parameters
    ----------
    argv : sequence of str or None
        The arguments (without the program name); None means sys.argv[1:].

    Returns
    -------
    int
        The exit status: 0 on success, 1 if the conversion failed. Usage
        errors exit with status 2 from inside argparse.
    """
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = build_argument_parser()
    args = parser.parse_args(arguments)

    if not (math.isfinite(args.spacing) and args.spacing >= CONFIG.MIN_LATTICE_TO_ANGSTROMS):
        parser.error('--spacing must be a number of at least %g A (the smallest LATTICE_TO_ANGSTROMS '
                     'PIMMS accepts)' % CONFIG.MIN_LATTICE_TO_ANGSTROMS)
    if args.padding < 0:
        parser.error('--padding cannot be negative')
    try:
        bead_map = parse_bead_map(args.map)
    except ValueError as error:
        parser.error(str(error))
    chains = None
    if args.chains is not None:
        chains = [chain.strip() for chain in args.chains.split(',') if chain.strip()]
        if not chains:
            parser.error('--chains needs at least one chain ID')
    if args.mobile and args.no_freeze:
        parser.error('--mobile chooses what the freeze file leaves free, so it cannot be used with --no-freeze')
    try:
        mobile = [parse_mobile_selection(selection) for selection in args.mobile]
    except ValueError as error:
        parser.error(str(error))

    structure = Path(args.structure)
    restart_path, report_path, default_freeze_path = output_paths(structure, args.output, args.report, args.freeze)
    outputs = {'restart file': restart_path, 'report': report_path}
    freeze_path: Optional[Path] = None
    if not args.no_freeze:
        freeze_path = default_freeze_path
        outputs['freeze file'] = freeze_path
    problem = _output_problem(outputs, args.overwrite, structure)
    if problem is not None:
        sys.stderr.write('ERROR: %s\n' % problem)
        return 1

    options = ConversionOptions(spacing=args.spacing, padding=args.padding, box=args.box, model=args.model,
                                chains=chains, bead_position=args.bead_position, keep_hetero=args.hetero,
                                bead_map=bead_map, mobile=mobile)
    try:
        conversion = convert_structure(structure, options)
        command = 'pdb2pimms ' + shlex.join(arguments)
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        report_text = render_report(conversion, restart_path, report_path, command, pimms_version(), timestamp,
                                    freeze_path)
        write_outputs(conversion, restart_path, report_text, report_path, freeze_path)
    except (StructureFileException, StructureConversionException, OSError) as error:
        sys.stderr.write('ERROR: %s\n' % error)
        return 1

    _print_summary(conversion, restart_path, report_path, freeze_path)
    return 0
