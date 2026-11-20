"""
Tests for pdb2pimms, which converts a PDB or mmCIF structure into a PIMMS restart
file and a Markdown report.

The structures are built here rather than shipped: an ideal alpha-helix (CA
radius 2.3 A, rise 1.5 A and 100 degrees per residue, so consecutive CA atoms are
3.83 A apart) with N and C atoms placed so that every peptide C-N distance is
1.3 A, and a nucleic-acid helix (C4' radius 9 A, rise 3.4 A, 36 degrees per
nucleotide) with O3' and P atoms placed so that every phosphodiester link is
0.3 A long. Neither is a real structure; both have exactly the connectivity
pdb2pimms tests for. The strongest tests run PIMMS itself from a converted file,
once with every chain free and once with the converted chains frozen.
"""

import gzip
import importlib.util
import itertools
import os
import pathlib
import pickle
import shutil
import subprocess
import sys

import numpy as np
import pytest

from pimms import pdb2pimms as p2p
from pimms.data_structures import FreezeFile
from pimms.latticeExceptions import StructureConversionException, StructureFileException
from pimms.restart import RestartObject
from pimms.structure_reader import read_structure


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
PDB2PIMMS_SCRIPT = REPO_ROOT / "scripts" / "pdb2pimms"
PIMMS_SCRIPT = REPO_ROOT / "scripts" / "PIMMS"

# geometry of the synthetic structures (see the module docstring)
HELIX_RADIUS = 2.3
HELIX_RISE = 1.5
HELIX_TWIST_DEGREES = 100.0
NUCLEIC_RADIUS = 9.0
NUCLEIC_RISE = 3.4
NUCLEIC_TWIST_DEGREES = 36.0
LINK_FRACTION = 0.33

PROTEIN_12 = ['MET', 'ALA', 'GLY', 'HIS', 'LYS', 'GLU', 'LEU', 'SER', 'TRP', 'ASP', 'VAL', 'ARG']

# the nucleosome demo: its inputs (not whatever an earlier run left there), its
# runs, and the settings that cut each run to a few steps that still write
# production frames
NUCLEOSOME_DEMO = REPO_ROOT / "demos" / "nucleosome_from_pdb"
NUCLEOSOME_DEMO_RUNS = ("H1_CTD", "ProTa", "H1_CTD_ProTa")
NUCLEOSOME_DEMO_INPUTS = ("1kx5.cif.gz", "params.prm", "analyze_contacts.py", "reset.sh")
# the histone tails the demo leaves free, as its README's conversion command names them
NUCLEOSOME_DEMO_TAILS = ("A:1-26", "E:1-26", "B:1-17", "F:1-17", "C:1-11", "C:118-128", "G:1-11", "G:118-128",
                         "D:1-20", "H:1-20")
DEMO_SMOKE_SETTINGS = {"N_STEPS": "40", "EQUILIBRATION": "10", "XTC_FREQ": "10", "EN_FREQ": "10",
                       "PRINT_FREQ": "40", "ENERGY_CHECK": "10", "ANALYSIS_FREQ": "10"}

# --protein-mode: human ubiquitin (UniProt P0CG48, residues 1-76), which metapredict
# predicts folded, between two Gly/Ser-rich tails it predicts disordered (metapredict
# 3.1 puts the folded domain at residues 32-100 of the 136)
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
DISORDERED_TAIL = "GSEGSSGSQSPGSDGSGSESGSPSGSGSQS"
THREE_LETTER = {code: name for name, code in p2p.AMINO_ACID_CODES.items()}


# ---------------------------------------------------------------------------
# building synthetic structures
# ---------------------------------------------------------------------------

def atom(name, resname, chain, resseq, xyz, *, record="ATOM", altloc="", icode="", occupancy=1.0,
         element=None, model=1):
    """One atom of a synthetic structure, as a dict."""
    return dict(record=record, name=name, resname=resname, chain=chain, resseq=resseq,
                xyz=np.asarray(xyz, dtype=float), altloc=altloc, icode=icode, occupancy=occupancy,
                element=element if element is not None else name.lstrip("0123456789")[0], model=model)


def helix_trace(n, origin=(0.0, 0.0, 0.0), radius=HELIX_RADIUS, rise=HELIX_RISE,
                twist=HELIX_TWIST_DEGREES):
    """Points on an ideal helix along z."""
    angles = np.radians(twist) * np.arange(n)
    trace = np.column_stack([radius * np.cos(angles), radius * np.sin(angles), rise * np.arange(n)])
    return trace + np.asarray(origin, dtype=float)


def _neighbour_directions(trace, i):
    """Directions from point i towards its previous and next points (extrapolated at the ends)."""
    n = len(trace)
    if n == 1:
        unit = np.array([1.0, 0.0, 0.0])
        return -unit, unit
    previous = trace[i - 1] - trace[i] if i > 0 else trace[i] - trace[i + 1]
    following = trace[i + 1] - trace[i] if i < n - 1 else trace[i] - trace[i - 1]
    return previous, following


def protein_atoms(chain, trace, resnames, *, start=1, record="ATOM", backbone=True, model=1):
    """N, CA and C atoms (CA only if not backbone) for a protein chain whose CA atoms follow trace."""
    atoms = []
    for i, (ca, resname) in enumerate(zip(trace, resnames)):
        previous, following = _neighbour_directions(trace, i)
        if backbone:
            atoms.append(atom("N", resname, chain, start + i, ca + LINK_FRACTION * previous,
                              record=record, model=model))
        atoms.append(atom("CA", resname, chain, start + i, ca, record=record, model=model))
        if backbone:
            atoms.append(atom("C", resname, chain, start + i, ca + LINK_FRACTION * following,
                              record=record, model=model))
    return atoms


def nucleic_atoms(chain, trace, resnames, *, rna=True, start=1, record="ATOM"):
    """P, O5', C4', C2', (O2'), O3' and N1 atoms for a nucleic-acid chain whose C4' atoms follow trace."""
    atoms = []
    for i, (c4, resname) in enumerate(zip(trace, resnames)):
        previous, following = _neighbour_directions(trace, i)
        resseq = start + i
        atoms.append(atom("P", resname, chain, resseq, c4 + 0.4 * previous, record=record))
        atoms.append(atom("O5'", resname, chain, resseq, c4 + 0.25 * previous, record=record))
        atoms.append(atom("C4'", resname, chain, resseq, c4, record=record))
        atoms.append(atom("C2'", resname, chain, resseq, c4 + np.array([0.0, 0.0, 1.5]), record=record))
        if rna:
            atoms.append(atom("O2'", resname, chain, resseq, c4 + np.array([0.0, 1.0, 1.8]), record=record))
        atoms.append(atom("O3'", resname, chain, resseq, c4 + 0.55 * following, record=record))
        atoms.append(atom("N1", resname, chain, resseq, c4 + np.array([1.5, 0.0, 0.0]), record=record))
    return atoms


def pdb_line(serial, a):
    """A PDB ATOM/HETATM line with the version 3 fixed columns."""
    name = a["name"]
    name_field = name if len(name) == 4 else " " + name.ljust(3)
    resname_field = a["resname"].rjust(3).ljust(4)
    x, y, z = a["xyz"]
    return ("%-6s%5d %-4s%1s%4s%1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
            % (a["record"], serial % 100000, name_field, a["altloc"], resname_field, a["chain"],
               a["resseq"], a["icode"], x, y, z, a["occupancy"], 0.0, a["element"]))


def write_pdb(path, atoms, *, header_lines=()):
    """Write atoms as a PDB file (MODEL records when more than one model is present)."""
    models = list(dict.fromkeys(a["model"] for a in atoms))
    lines = list(header_lines)
    serial = 1
    for model in models:
        if len(models) > 1:
            lines.append("MODEL     %4d" % model)
        for a in atoms:
            if a["model"] == model:
                lines.append(pdb_line(serial, a))
                serial += 1
        if len(models) > 1:
            lines.append("ENDMDL")
    lines.append("END")
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n")


def _cif_text(value):
    """Quote a value for mmCIF where needed."""
    if value == "":
        return "?"
    if "'" in value or " " in value:
        return '"%s"' % value
    return value


def write_cif(path, atoms, *, extra=""):
    """Write atoms as an mmCIF file with an _atom_site loop."""
    lines = ["data_TEST", "#", "_entry.id TEST", "#", "loop_"]
    for tag in ("group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
                "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
                "occupancy", "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num"):
        lines.append("_atom_site.%s" % tag)
    for serial, a in enumerate(atoms, start=1):
        x, y, z = a["xyz"]
        lines.append(" ".join([
            a["record"], str(serial), a["element"], _cif_text(a["name"]), a["altloc"] or ".",
            a["resname"], a["chain"], str(a["resseq"]), a["icode"] or "?", "%.3f" % x, "%.3f" % y,
            "%.3f" % z, "%.2f" % a["occupancy"], str(a["resseq"]), a["chain"], str(a["model"])]))
    lines.append("#")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n" + extra)


def convert(path, **options):
    """Run the conversion in memory with the given ConversionOptions fields."""
    return p2p.convert_structure(path, p2p.ConversionOptions(**options))


def assert_valid_pimms_configuration(conversion):
    """Every chain is a Moore walk, no two beads share a site, every bead is inside the box."""
    seen = set()
    for walk in conversion.sites:
        for site in walk:
            assert all(0 <= c < d for c, d in zip(site, conversion.dimensions))
            assert site not in seen
            seen.add(site)
        for a, b in zip(walk, walk[1:]):
            assert max(abs(i - j) for i, j in zip(a, b)) == 1


# ---------------------------------------------------------------------------
# reading PDB files
# ---------------------------------------------------------------------------

def test_pdb_reader_reads_every_field(tmp_path):
    lines = [
        "HEADER    %-40s%9s   %4s" % ("TEST PROTEIN", "01-JAN-26", "9XYZ"),
        "MODRES 9XYZ MSE A    2  MET  SELENOMETHIONINE",
        pdb_line(1, atom("CA", "ALA", "A", 1, (1.0, 2.0, 3.0))),
        pdb_line(2, atom("CA", "MSE", "A", 2, (4.0, 5.0, 6.0), record="HETATM", altloc="B",
                         occupancy=0.6, icode="X")),
    ]
    # a hydrogen with no element column, and an atom with a blank chain but a segment ID
    lines.append(pdb_line(3, atom("1HB", "ALA", "A", 1, (0.0, 0.0, 0.0)))[:76])
    segment_line = pdb_line(4, atom("CA", "GLY", " ", 3, (7.0, 8.0, 9.0)))
    lines.append(segment_line[:72] + "SEGA" + segment_line[76:])
    path = tmp_path / "x.pdb"
    path.write_text("\n".join(lines) + "\n")

    structure = read_structure(path)
    assert structure.file_format == "PDB"
    assert structure.entry_id == "9XYZ"
    assert structure.parent_residues == {"MSE": "MET"}
    first, second, hydrogen, segment = structure.atoms
    assert (first.record, first.name, first.resname, first.chain, first.resseq) == ("ATOM", "CA", "ALA", "A", "1")
    assert (first.x, first.y, first.z) == (1.0, 2.0, 3.0)
    assert (second.record, second.altloc, second.icode, second.occupancy) == ("HETATM", "B", "X", 0.6)
    assert hydrogen.element == "H"
    assert segment.chain == "SEGA"
    assert structure.chains_from_segment_ids
    assert structure.models == [1]


def test_pdb_reader_reads_models_in_order(tmp_path):
    trace = helix_trace(3)
    atoms = protein_atoms("A", trace, PROTEIN_12[:3], model=4) + protein_atoms("A", trace + 1.0, PROTEIN_12[:3], model=7)
    write_pdb(tmp_path / "m.pdb", atoms)
    structure = read_structure(tmp_path / "m.pdb")
    assert structure.models == [4, 7]
    assert sum(1 for a in structure.atoms if a.model == 7) == 9


def test_pdb_reader_refuses_a_truncated_atom_line_and_a_bad_coordinate(tmp_path):
    line = pdb_line(1, atom("CA", "ALA", "A", 1, (1.0, 2.0, 3.0)))
    (tmp_path / "short.pdb").write_text(line[:50] + "\n")
    with pytest.raises(StructureFileException, match="too short"):
        read_structure(tmp_path / "short.pdb")
    (tmp_path / "bad.pdb").write_text(line[:30] + "   abcde" + line[38:] + "\n")
    with pytest.raises(StructureFileException, match="x coordinate"):
        read_structure(tmp_path / "bad.pdb")


def test_files_without_atoms_or_with_an_unknown_format_are_refused(tmp_path):
    (tmp_path / "empty.pdb").write_text("HEADER    NOTHING\nEND\n")
    with pytest.raises(StructureFileException, match="no ATOM or HETATM"):
        read_structure(tmp_path / "empty.pdb")
    (tmp_path / "mystery.txt").write_text("hello\n")
    with pytest.raises(StructureFileException, match="Cannot tell"):
        read_structure(tmp_path / "mystery.txt")
    with pytest.raises(StructureFileException, match="no such file"):
        read_structure(tmp_path / "absent.pdb")


def test_format_is_sniffed_when_the_extension_says_nothing(tmp_path):
    atoms = protein_atoms("A", helix_trace(4), PROTEIN_12[:4])
    write_pdb(tmp_path / "a.pdb", atoms)
    write_cif(tmp_path / "a.cif", atoms)
    os.rename(tmp_path / "a.pdb", tmp_path / "a_pdb.txt")
    os.rename(tmp_path / "a.cif", tmp_path / "a_cif.txt")
    assert read_structure(tmp_path / "a_pdb.txt").file_format == "PDB"
    assert read_structure(tmp_path / "a_cif.txt").file_format == "mmCIF"


# ---------------------------------------------------------------------------
# reading mmCIF files
# ---------------------------------------------------------------------------

CIF_WITH_EVERYTHING = """\
# a comment before the block
data_FIRST
_entry.id   1ABC
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.parent_comp_id
_pdbx_struct_mod_residue.details
1 PSU U
;a text field
over two lines
;
2 5MC C 'quoted details'
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM C "C4'" . PSU X 1 ? 1.000 2.000 3.000 1.00 10 B 1
ATOM P P     A PSU X 1 A 4.000 5.000 6.000 0.50 10 B 1 # trailing comment
HETATM O O   . HOH Y . ? 7.000 8.000 9.000 ?    501 B 1
#
data_SECOND
loop_
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
1 1 1 CA ALA Z 1
"""


def test_cif_reader_handles_quotes_text_fields_missing_values_and_blocks(tmp_path):
    path = tmp_path / "x.cif"
    path.write_text(CIF_WITH_EVERYTHING)
    structure = read_structure(path)
    assert structure.file_format == "mmCIF"
    assert structure.entry_id == "1ABC"
    assert structure.parent_residues == {"PSU": "U", "5MC": "C"}
    # the second data block is not read
    assert len(structure.atoms) == 3
    first, second, water = structure.atoms
    assert (first.name, first.altloc, first.icode, first.chain, first.resseq) == ("C4'", "", "", "B", "10")
    assert (second.altloc, second.icode, second.occupancy) == ("A", "A", 0.5)
    # a missing occupancy is full occupancy; auth_seq_id wins over a missing label_seq_id
    assert (water.record, water.occupancy, water.resseq) == ("HETATM", 1.0, "501")


@pytest.mark.parametrize("text, message", [
    ("data_X\nloop_\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n1.0\n", "whole number of rows"),
    ("data_X\nloop_\n_atom_site.label_atom_id\nCA\n", "no _atom_site.cartn_x column"),
    ("data_X\n_entry.id 'unclosed\n", "not closed"),
    ("data_X\n_entry.id\n;never closed\n", "never closed"),
    ("data_X\n_entry.id ABC\n", "no _atom_site category"),
])
def test_cif_reader_refuses_malformed_files(tmp_path, text, message):
    path = tmp_path / "bad.cif"
    path.write_text(text)
    with pytest.raises(StructureFileException, match=message):
        read_structure(path)


def test_pdb_cif_and_gzip_give_identical_restart_files(tmp_path):
    atoms = (protein_atoms("A", helix_trace(15), (PROTEIN_12 + PROTEIN_12)[:15])
             + nucleic_atoms("B", helix_trace(8, origin=(30.0, 0.0, 0.0), radius=NUCLEIC_RADIUS,
                                              rise=NUCLEIC_RISE, twist=NUCLEIC_TWIST_DEGREES), ["G", "C", "A", "U"] * 2))
    write_pdb(tmp_path / "s.pdb", atoms)
    write_pdb(tmp_path / "s.pdb.gz", atoms)
    write_cif(tmp_path / "s.cif", atoms)
    payloads = [p2p.restart_payload(convert(tmp_path / name)) for name in ("s.pdb", "s.pdb.gz", "s.cif")]
    assert payloads[0] == payloads[1] == payloads[2]
    assert [chain[1] for chain in payloads[0]["CHAINS"].values()] == ["MAGHKELSWDVRMAG", "gcaugcau"]


# ---------------------------------------------------------------------------
# fitting onto the lattice
# ---------------------------------------------------------------------------

def test_a_trace_already_on_a_lattice_is_reproduced_exactly(tmp_path):
    # a self-avoiding walk of random Moore steps, scaled by a spacing small enough
    # (2 A, so a body diagonal is 3.46 A) that CA-only residues stay bonded
    rng = np.random.default_rng(3)
    spacing = 2.0
    walk = [np.zeros(3, dtype=int)]
    seen = {(0, 0, 0)}
    while len(walk) < 40:
        step = p2p.MOORE_OFFSETS[rng.integers(len(p2p.MOORE_OFFSETS))]
        site = walk[-1] + step
        if tuple(site) not in seen:
            seen.add(tuple(site))
            walk.append(site)
    walk = np.array(walk)
    atoms = protein_atoms("A", walk * spacing + 10.0, ["ALA"] * len(walk), backbone=False)
    write_pdb(tmp_path / "lattice.pdb", atoms)

    conversion = convert(tmp_path / "lattice.pdb", spacing=spacing)
    assert len(conversion.segments) == 1
    sites = np.array(conversion.sites[0])
    assert np.array_equal(sites - sites[0], walk - walk[0])
    assert conversion.statistics.rmsd < 1e-9


def test_optimal_walk_can_revisit_a_site_and_the_repair_removes_it():
    targets = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
    nothing = np.zeros(0, dtype=np.int64)
    unrepaired = p2p._optimal_walk(targets, nothing, {}, 2)
    assert unrepaired[0] == unrepaired[2]

    walk = p2p._self_avoiding_walk(targets, nothing, 2)
    assert len(set(walk)) == 3
    assert all(max(abs(i - j) for i, j in zip(a, b)) == 1 for a, b in zip(walk, walk[1:]))
    # the best self-avoiding walk moves one end bead by one site: total cost 1.01
    cost = sum(p2p._squared_distance(site, target) for site, target in zip(walk, targets))
    assert cost == pytest.approx(1.01)


def test_later_chains_avoid_earlier_ones():
    # two chains whose targets coincide: the second must find other sites
    targets = [np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])] * 2
    walks, methods = p2p.fit_to_lattice(targets)
    assert not set(walks[0]) & set(walks[1])
    assert all(method.startswith("dynamic programming") for method in methods)


def test_greedy_search_takes_over_when_the_dynamic_programme_fails(monkeypatch):
    monkeypatch.setattr(p2p, "_self_avoiding_walk", lambda *args: None)
    targets = [helix_trace(20) / 3.65, helix_trace(10, origin=(1.0, 0.0, 0.0)) / 3.65]
    walks, methods = p2p.fit_to_lattice(targets)
    assert methods == ["greedy search", "greedy search"]
    assert not set(walks[0]) & set(walks[1])
    for walk in walks:
        assert len(set(walk)) == len(walk)
        assert all(max(abs(i - j) for i, j in zip(a, b)) == 1 for a, b in zip(walk, walk[1:]))


def test_greedy_search_backtracks_and_gives_up_when_boxed_in():
    # every site around the origin is taken, so a two-bead chain starting there has nowhere to go
    occupied = {site for site in itertools.product(range(-2, 3), repeat=3) if site != (0, 0, 0)}
    assert p2p._greedy_walk(np.zeros((2, 3)), occupied) is None


def test_a_chain_that_cannot_be_placed_is_refused(monkeypatch):
    monkeypatch.setattr(p2p, "_self_avoiding_walk", lambda *args: None)
    monkeypatch.setattr(p2p, "_greedy_walk", lambda *args: None)
    with pytest.raises(StructureConversionException, match="could not be placed"):
        p2p.fit_to_lattice([np.zeros((3, 3))])


# ---------------------------------------------------------------------------
# chains, bonds and breaks
# ---------------------------------------------------------------------------

def test_a_gap_splits_a_chain_and_is_reported(tmp_path):
    trace = helix_trace(20)
    resnames = (PROTEIN_12 + PROTEIN_12)[:20]
    atoms = [a for a in protein_atoms("A", trace, resnames) if a["resseq"] not in (9, 10, 11)]
    write_pdb(tmp_path / "gap.pdb", atoms)
    conversion = convert(tmp_path / "gap.pdb")
    assert [len(segment.residues) for segment in conversion.segments] == [8, 9]
    (gap,) = conversion.breaks
    assert (gap.chain, gap.before, gap.after) == ("A", "A:SER8", "A:ARG12")
    assert "C-N distance" in gap.reason
    assert_valid_pimms_configuration(conversion)


def test_bonds_use_the_linking_atoms_then_the_backbone_distance():
    def residue(atoms, kind=p2p.KIND_AMINO_ACID):
        return p2p.Residue(chain="A", resseq="1", icode="", resname="ALA", record="ATOM",
                           atoms={name: np.array(xyz, dtype=float) for name, xyz in atoms.items()},
                           hydrogens=[], kind=kind)
    linked = p2p._bond_between(residue({"CA": (0, 0, 0), "C": (1, 0, 0)}),
                               residue({"N": (2.3, 0, 0), "CA": (3.8, 0, 0)}))
    assert linked == (True, "")
    # CA atoms close but the peptide link 2.5 A long: the link decides
    bonded, reason = p2p._bond_between(residue({"CA": (0, 0, 0), "C": (0.5, 0, 0)}),
                                       residue({"N": (3.0, 0, 0), "CA": (3.8, 0, 0)}))
    assert not bonded and "C-N distance 2.50" in reason
    # CA-only residues fall back to the CA-CA distance
    assert p2p._bond_between(residue({"CA": (0, 0, 0)}), residue({"CA": (3.8, 0, 0)}))[0]
    bonded, reason = p2p._bond_between(residue({"CA": (0, 0, 0)}), residue({"CA": (5.0, 0, 0)}))
    assert not bonded and "CA-CA distance 5.00" in reason
    # an amino acid is never bonded to a nucleotide
    bonded, reason = p2p._bond_between(residue({"CA": (0, 0, 0)}),
                                       residue({"C4'": (1, 0, 0)}, kind=p2p.KIND_RNA))
    assert not bonded and "followed by" in reason


def test_identical_sequences_share_a_chain_type(tmp_path):
    atoms = (protein_atoms("A", helix_trace(6), PROTEIN_12[:6])
             + protein_atoms("B", helix_trace(6, origin=(20, 0, 0)), PROTEIN_12[:6])
             + protein_atoms("C", helix_trace(5, origin=(0, 20, 0)), PROTEIN_12[:5]))
    write_pdb(tmp_path / "three.pdb", atoms)
    conversion = convert(tmp_path / "three.pdb")
    assert conversion.chain_types == [0, 0, 1]


# ---------------------------------------------------------------------------
# bead codes
# ---------------------------------------------------------------------------

def test_standard_variant_and_modified_amino_acids(tmp_path):
    resnames = ["NALA", "HID", "MSE", "XYZ", "SEP", "CYX"]
    atoms = protein_atoms("A", helix_trace(6), resnames)
    write_pdb(tmp_path / "mod.pdb", atoms,
              header_lines=["MODRES TEST XYZ A    4  TRP  SOMETHING MADE UP",
                            "MODRES TEST SEP A    5  THR  THE FILE OVERRIDES THE TABLE"])
    conversion = convert(tmp_path / "mod.pdb")
    residues = conversion.segments[0].residues
    assert conversion.segments[0].sequence == "AHMWTC"
    assert [r.parent_source for r in residues] == [
        "other name", "other name", "built-in table", "file", "file", "other name"]


def test_nucleotides_get_lower_case_codes_and_the_sugar_decides_rna_or_dna(tmp_path):
    trace = helix_trace(4, radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE, twist=NUCLEIC_TWIST_DEGREES)
    rna = nucleic_atoms("R", trace, ["A", "C", "G", "U"], rna=True)
    # old-style DNA (single-letter names, deoxyribose) and current DNA names
    old_dna = nucleic_atoms("D", trace + 40.0, ["A", "C", "G", "T"], rna=False)
    new_dna = nucleic_atoms("E", trace - 40.0, ["DA", "DC", "DG", "DT"], rna=False)
    write_pdb(tmp_path / "na.pdb", rna + old_dna + new_dna)
    conversion = convert(tmp_path / "na.pdb")
    assert [s.sequence for s in conversion.segments] == ["acgu", "acgt", "acgt"]
    kinds = [{r.kind for r in s.residues} for s in conversion.segments]
    assert kinds == [{p2p.KIND_RNA}, {p2p.KIND_DNA}, {p2p.KIND_DNA}]
    # RNA and DNA share a, c and g, which the user is told about
    shared = [w for w in conversion.warnings if "stands for DNA nucleotide and RNA nucleotide" in w]
    assert len(shared) == 3


def test_unknown_residues_with_a_backbone_become_x_or_n_and_others_are_hetero(tmp_path):
    protein = protein_atoms("A", helix_trace(3), ["ALA", "QQQ", "ALA"])
    trace = helix_trace(3, origin=(30, 0, 0), radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE,
                        twist=NUCLEIC_TWIST_DEGREES)
    nucleic = nucleic_atoms("B", trace, ["G", "ZZN", "G"])
    ligand = [atom("C1", "LIG", "C", 1, (0, 30, 0), record="HETATM"),
              atom("C2", "LIG", "C", 1, (1.5, 30, 0), record="HETATM")]
    write_pdb(tmp_path / "unknown.pdb", protein + nucleic + ligand)
    conversion = convert(tmp_path / "unknown.pdb")
    assert [s.sequence for s in conversion.segments] == ["AXA", "gng"]
    assert conversion.excluded.hetero == {"LIG": 1}
    assert any("QQQ, ZZN" in w for w in conversion.warnings)


def test_hetero_groups_become_one_bead_chains_with_their_own_codes(tmp_path):
    protein = protein_atoms("A", helix_trace(5), PROTEIN_12[:5])
    zinc = [atom("ZN", "ZN", "A", 101, (10, 0, 0), record="HETATM", element="ZN")]
    heme = [atom("FE", "HEM", "A", 102, (0, 10, 0), record="HETATM", element="FE"),
            atom("C1", "HEM", "A", 102, (2, 10, 0), record="HETATM")]
    water = [atom("O", "HOH", "A", 201, (0, 0, 10), record="HETATM")]
    write_pdb(tmp_path / "het.pdb", protein + zinc + heme + water)

    dropped = convert(tmp_path / "het.pdb")
    assert len(dropped.segments) == 1
    assert dropped.excluded.hetero == {"ZN": 1, "HEM": 1}
    assert dropped.excluded.waters == 1

    kept = convert(tmp_path / "het.pdb", keep_hetero=True, bead_map={"HEM": "h"})
    assert [s.sequence for s in kept.segments] == ["MAGHK", "b", "h"]
    heme_bead = kept.segments[2].residues[0]
    assert np.allclose(heme_bead.position, (1.0, 10.0, 0.0))
    assert heme_bead.position_source == "heavy-atom centroid"
    assert not kept.excluded.hetero


def test_bead_map_overrides_and_its_no_effect_warnings(tmp_path):
    atoms = protein_atoms("A", helix_trace(4), ["MET", "MSE", "ALA", "ALA"])
    atoms.append(atom("C1", "LIG", "A", 50, (0, 10, 0), record="HETATM"))
    atoms.append(atom("O", "HOH", "A", 60, (0, 0, 10), record="HETATM"))
    write_pdb(tmp_path / "map.pdb", atoms)
    conversion = convert(tmp_path / "map.pdb", bead_map={"MSE": "Z", "LIG": "l", "HOH": "w", "ABC": "q"})
    assert conversion.segments[0].sequence == "MZAA"
    warnings = " ".join(conversion.warnings)
    assert "--map LIG=l had no effect" in warnings and "only converted with --hetero" in warnings
    assert "--map HOH=w had no effect: water is never converted" in warnings
    assert "--map ABC=q had no effect" in warnings


@pytest.mark.parametrize("character, message", [
    ("0", "solvent"), ("#", "comment"), (" ", "whitespace"), ("\t", "whitespace"), ("AB", "exactly one"),
    ("", "exactly one"),
])
def test_unusable_bead_characters_are_refused(character, message):
    with pytest.raises(ValueError, match=message):
        p2p.check_bead_character(character, "test")


def test_parse_bead_map():
    assert p2p.parse_bead_map(["MSE=M", " SEP=s", "MSE=M"]) == {"MSE": "M", "SEP": "s"}
    for bad in (["MSE"], ["=M"], ["MSE=M", "MSE=Z"]):
        with pytest.raises(ValueError):
            p2p.parse_bead_map(bad)


# ---------------------------------------------------------------------------
# what is read and left out
# ---------------------------------------------------------------------------

def test_alternate_locations_keep_the_most_occupied_and_hydrogens_never_move_a_bead(tmp_path):
    atoms = protein_atoms("A", helix_trace(3), ["ALA", "SER", "ALA"])
    middle_ca = next(a for a in atoms if a["resseq"] == 2 and a["name"] == "CA")
    atoms.remove(middle_ca)
    position = list(atoms).index(next(a for a in atoms if a["resseq"] == 2 and a["name"] == "C"))
    atoms[position:position] = [
        atom("CA", "SER", "A", 2, middle_ca["xyz"] + 0.4, altloc="A", occupancy=0.4),
        atom("CA", "SER", "A", 2, middle_ca["xyz"], altloc="B", occupancy=0.6),
        atom("H", "SER", "A", 2, middle_ca["xyz"] + 5.0),
    ]
    write_pdb(tmp_path / "alt.pdb", atoms)
    conversion = convert(tmp_path / "alt.pdb", bead_position=p2p.BEAD_POSITION_BACKBONE)
    middle = conversion.segments[0].residues[1]
    assert (middle.altloc, middle.altlocs_present, middle.code) == ("B", ("A", "B"), "S")
    # the PDB file holds coordinates to 0.001 A
    assert np.allclose(middle.position, middle_ca["xyz"], atol=1e-3)
    assert conversion.excluded.hydrogens == 1
    assert [r.label for r in conversion.excluded.altloc_residues] == ["A:SER2"]


def test_centroid_positions_use_heavy_atoms(tmp_path):
    atoms = [atom("N", "GLY", "A", 1, (0, 0, 0)), atom("CA", "GLY", "A", 1, (1, 0, 0)),
             atom("C", "GLY", "A", 1, (2, 0, 0)), atom("H", "GLY", "A", 1, (50, 50, 50))]
    write_pdb(tmp_path / "c.pdb", atoms)
    conversion = convert(tmp_path / "c.pdb", bead_position=p2p.BEAD_POSITION_CENTROID)
    assert np.allclose(conversion.segments[0].residues[0].position, (1.0, 0.0, 0.0))


def test_a_residue_without_its_backbone_atom_is_placed_at_its_centroid(tmp_path):
    atoms = protein_atoms("A", helix_trace(3), PROTEIN_12[:3])
    atoms = [a for a in atoms if not (a["resseq"] == 3 and a["name"] == "CA")]
    write_pdb(tmp_path / "noca.pdb", atoms)
    conversion = convert(tmp_path / "noca.pdb")
    assert conversion.position_fallbacks == ["A:GLY3"]


def test_model_and_chain_selection(tmp_path):
    first = protein_atoms("A", helix_trace(4), PROTEIN_12[:4], model=1)
    first += protein_atoms("B", helix_trace(4, origin=(20, 0, 0)), PROTEIN_12[4:8], model=1)
    second = protein_atoms("A", helix_trace(4, origin=(1, 1, 1)), PROTEIN_12[:4], model=2)
    second += protein_atoms("B", helix_trace(4, origin=(21, 1, 1)), PROTEIN_12[4:8], model=2)
    write_pdb(tmp_path / "models.pdb", first + second)

    conversion = convert(tmp_path / "models.pdb", model=2, chains=["B"])
    assert conversion.model == 2
    assert [s.sequence for s in conversion.segments] == ["KELS"]
    assert conversion.excluded.unselected_chains == ["A"]
    assert conversion.excluded.other_models == 1
    with pytest.raises(StructureConversionException, match="--model 5 is not in"):
        convert(tmp_path / "models.pdb", model=5)
    with pytest.raises(StructureConversionException, match="chain\\(s\\) Q"):
        convert(tmp_path / "models.pdb", chains=["Q"])


def test_a_model_with_nothing_to_convert_is_refused(tmp_path):
    write_pdb(tmp_path / "water.pdb", [atom("O", "HOH", "A", 1, (0, 0, 0), record="HETATM")])
    with pytest.raises(StructureConversionException, match="nothing to convert"):
        convert(tmp_path / "water.pdb")


# ---------------------------------------------------------------------------
# the box and the restart file
# ---------------------------------------------------------------------------

def test_box_is_padded_or_given_and_the_structure_centred(tmp_path):
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(12), PROTEIN_12))
    padded = convert(tmp_path / "h.pdb", padding=4)
    every_site = np.array([site for walk in padded.sites for site in walk])
    assert np.array_equal(every_site.min(axis=0), [4, 4, 4])
    assert np.array_equal(np.array(padded.dimensions) - every_site.max(axis=0) - 1, [4, 4, 4])

    boxed = convert(tmp_path / "h.pdb", box=[30, 31, 32])
    assert boxed.dimensions == [30, 31, 32]
    sites = np.array([site for walk in boxed.sites for site in walk])
    low = sites.min(axis=0)
    high = np.array(boxed.dimensions) - sites.max(axis=0) - 1
    assert np.all(np.abs(low - high) <= 1)
    with pytest.raises(StructureConversionException, match="cannot be used"):
        convert(tmp_path / "h.pdb", box=[3, 30, 30])

    # a tiny structure still gets the 7-site minimum box PIMMS needs
    write_pdb(tmp_path / "one.pdb", protein_atoms("A", helix_trace(1), ["ALA"]))
    assert convert(tmp_path / "one.pdb", padding=0).dimensions == [7, 7, 7]


def test_written_restart_file_is_accepted_by_the_pimms_reader(tmp_path):
    atoms = (protein_atoms("A", helix_trace(30), (PROTEIN_12 * 3)[:30])
             + protein_atoms("B", helix_trace(30, origin=(6.0, 0.0, 0.0)), (PROTEIN_12 * 3)[:30]))
    write_pdb(tmp_path / "pair.pdb", atoms)
    conversion = convert(tmp_path / "pair.pdb")
    assert_valid_pimms_configuration(conversion)
    restart_path = tmp_path / "pair_restart.pimms"
    p2p.write_outputs(conversion, restart_path, "report", tmp_path / "pair_report.md")

    restart = RestartObject()
    restart.build_from_file(str(restart_path))
    assert restart.dimensions == conversion.dimensions
    assert restart.hardwall is True
    assert restart.step is None and restart.rng_python is None
    assert sorted(restart.chains) == [1, 2]
    assert restart.chains[1][1] == restart.chains[2][1] == conversion.segments[0].sequence
    assert restart.chains[1][2] == restart.chains[2][2] == 0
    assert (tmp_path / "pair_report.md").read_text() == "report"
    # nothing temporary is left behind
    assert sorted(p.name for p in tmp_path.iterdir()) == ["pair.pdb", "pair_report.md", "pair_restart.pimms"]


def test_a_restart_file_pimms_refuses_is_never_put_in_place(tmp_path):
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(5), PROTEIN_12[:5]))
    conversion = convert(tmp_path / "h.pdb")
    conversion.sites[0][1] = conversion.sites[0][0]   # two beads on one site
    with pytest.raises(StructureConversionException, match="bug in pdb2pimms"):
        p2p.write_outputs(conversion, tmp_path / "r.pimms", "report", tmp_path / "r_report.md")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb"]


# ---------------------------------------------------------------------------
# the freeze file
# ---------------------------------------------------------------------------

def _two_helices_and_an_ion(path):
    """Two identical helices and a zinc ion: three PIMMS chains with --hetero."""
    atoms = (protein_atoms("A", helix_trace(10), PROTEIN_12[:10])
             + protein_atoms("B", helix_trace(10, origin=(8.0, 0.0, 0.0)), PROTEIN_12[:10]))
    atoms.append(atom("ZN", "ZN", "A", 101, (0, 12, 0), record="HETATM", element="ZN"))
    write_pdb(path, atoms)


def test_the_freeze_file_freezes_every_chain_and_pimms_reads_it(tmp_path):
    _two_helices_and_an_ion(tmp_path / "pair.pdb")
    conversion = convert(tmp_path / "pair.pdb", keep_hetero=True)
    freeze_path = tmp_path / "pair_freeze.txt"
    p2p.write_outputs(conversion, tmp_path / "pair_restart.pimms", "report", tmp_path / "pair_report.md",
                      freeze_path)

    assert FreezeFile(str(freeze_path)).chains == [1, 2, 3]
    text = freeze_path.read_text()
    assert "every chain of pair_restart.pimms" in text and "numbered\n# from 4" in text
    assert "C 1  # A:MET1-ASP10, 10 bead(s)\n" in text
    assert "C 3  # A:ZN101, 1 bead(s)\n" in text
    assert sorted(p.name for p in tmp_path.iterdir()) == ["pair.pdb", "pair_freeze.txt", "pair_report.md",
                                                          "pair_restart.pimms"]


def test_the_report_explains_the_freeze_file_only_when_one_is_written(tmp_path):
    _two_helices_and_an_ion(tmp_path / "pair.pdb")
    conversion = convert(tmp_path / "pair.pdb", keep_hetero=True)
    arguments = (conversion, tmp_path / "pair_restart.pimms", tmp_path / "pair_report.md", "pdb2pimms pair.pdb",
                 "9.9.9", "2026-01-01 00:00")
    with_freeze = p2p.render_report(*arguments, tmp_path / "pair_freeze.txt")
    assert "FREEZE_FILE          : pair_freeze.txt" in with_freeze
    assert "EXTRA_CHAIN          : <count> <sequence>" in with_freeze
    assert "(chainIDs 1 to 3)" in with_freeze and "numbers them from 4" in with_freeze
    assert "FREEZE_FILE" not in p2p.render_report(*arguments)


def test_a_line_break_in_a_name_cannot_add_a_chain_to_the_freeze_file(tmp_path):
    # an mmCIF text field can put a line break in a chain name; unguarded, the
    # rest of the name would become a freeze line of its own
    _two_helices_and_an_ion(tmp_path / "pair.pdb")
    conversion = convert(tmp_path / "pair.pdb")
    conversion.segments[0].residues[0].chain = "A\nC 99"
    text = p2p.freeze_file_text(conversion, tmp_path / "pair_restart.pimms")
    assert "C 1  # A C 99:MET1-ASP10, 10 bead(s)\n" in text
    p2p.write_outputs(conversion, tmp_path / "pair_restart.pimms", "report", tmp_path / "pair_report.md",
                      tmp_path / "pair_freeze.txt")
    assert FreezeFile(str(tmp_path / "pair_freeze.txt")).chains == [1, 2]


def test_a_freeze_file_pimms_refuses_is_never_put_in_place(tmp_path, monkeypatch):
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(5), PROTEIN_12[:5]))
    conversion = convert(tmp_path / "h.pdb")
    monkeypatch.setattr(p2p, "freeze_file_text", lambda conversion, restart_path: "c 1\n")
    with pytest.raises(StructureConversionException, match="refused the freeze file"):
        p2p.write_outputs(conversion, tmp_path / "r.pimms", "report", tmp_path / "r_report.md",
                          tmp_path / "r_freeze.txt")
    monkeypatch.setattr(p2p, "freeze_file_text", lambda conversion, restart_path: "C 1 2\n")
    with pytest.raises(StructureConversionException, match="not 1 to 1"):
        p2p.write_outputs(conversion, tmp_path / "r.pimms", "report", tmp_path / "r_report.md",
                          tmp_path / "r_freeze.txt")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb"]


# ---------------------------------------------------------------------------
# the report
# ---------------------------------------------------------------------------

def test_report_lists_bead_types_keyfile_lines_and_what_was_left_out(tmp_path):
    atoms = protein_atoms("A", helix_trace(6), ["MET", "MSE", "ALA", "ALA", "GLY", "GLY"])
    atoms += nucleic_atoms("B", helix_trace(3, origin=(30, 0, 0), radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE,
                                            twist=NUCLEIC_TWIST_DEGREES), ["G", "PSU", "C"])
    atoms.append(atom("O", "HOH", "A", 300, (0, 0, 20), record="HETATM"))
    atoms.append(atom("MG", "MG", "A", 301, (0, 20, 0), record="HETATM", element="MG"))
    write_pdb(tmp_path / "rep.pdb", atoms)
    conversion = convert(tmp_path / "rep.pdb")
    text = p2p.render_report(conversion, tmp_path / "rep_restart.pimms", tmp_path / "rep_report.md",
                             "pdb2pimms rep.pdb", "9.9.9", "2026-01-01 00:00")
    assert "The restart file uses these 6 bead types: `A G M c g u`" in text
    assert "6 solvation lines, 21 pair lines and 6 angle lines" in text
    assert "MSE (1), modified MET (parent from the built-in table)" in text
    assert "PSU (1), modified U (parent from the built-in table)" in text
    assert "RESTART_FILE         : rep_restart.pimms" in text
    assert "DIMENSIONS           : %s" % " ".join(str(v) for v in conversion.dimensions) in text
    assert "LATTICE_TO_ANGSTROMS : 3.65" in text
    assert "1 water molecules" in text and "MG x 1" in text
    assert ">1 A:MET1-GLY6\nMMAAGG" in text and ">2 B:G1-C3\nguc" in text


# ---------------------------------------------------------------------------
# the command line
# ---------------------------------------------------------------------------

def test_output_names_follow_the_input_or_the_restart_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    restart, report, freeze = p2p.output_paths(pathlib.Path("dir/1ubq.cif.gz"), None, None)
    assert restart.resolve() == (tmp_path / "1ubq_restart.pimms").resolve()
    assert report.resolve() == (tmp_path / "1ubq_report.md").resolve()
    assert freeze.resolve() == (tmp_path / "1ubq_freeze.txt").resolve()
    assert p2p.output_paths(pathlib.Path("x.pdb"), "out/run.pimms", None)[1] == pathlib.Path("out/run_report.md")
    assert p2p.output_paths(pathlib.Path("x.pdb"), "a_restart.pimms", None)[1] == pathlib.Path("a_report.md")
    assert p2p.output_paths(pathlib.Path("x.pdb"), "a.bin", "notes.md")[1] == pathlib.Path("notes.md")
    assert p2p.output_paths(pathlib.Path("x.pdb"), "out/run.pimms", None)[2] == pathlib.Path("out/run_freeze.txt")
    assert p2p.output_paths(pathlib.Path("x.pdb"), None, None, "hold.txt")[2] == pathlib.Path("hold.txt")


def test_main_writes_every_file_and_refuses_to_overwrite_them(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    assert p2p.main(["h.pdb"]) == 0
    out = capsys.readouterr().out
    assert "Bead types the parameter file must define (8): A E G H K L M S" in out
    assert "h_freeze.txt: freezes all 1 chain(s)" in out
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb", "h_freeze.txt", "h_report.md",
                                                          "h_restart.pimms"]

    before = (tmp_path / "h_restart.pimms").read_bytes()
    assert p2p.main(["h.pdb", "--spacing", "2.5"]) == 1
    assert "already exists; give --overwrite" in capsys.readouterr().err
    assert (tmp_path / "h_restart.pimms").read_bytes() == before
    assert p2p.main(["h.pdb", "--spacing", "2.5", "--overwrite"]) == 0
    assert (tmp_path / "h_restart.pimms").read_bytes() != before


def test_main_freeze_options(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))

    # --no-freeze writes no freeze file, and the report does not mention one
    assert p2p.main(["h.pdb", "--no-freeze"]) == 0
    assert "freezes all" not in capsys.readouterr().out
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb", "h_report.md", "h_restart.pimms"]
    assert "FREEZE_FILE" not in (tmp_path / "h_report.md").read_text()

    # --freeze names it, and the report's keyfile lines follow
    assert p2p.main(["h.pdb", "-o", "b.pimms", "--freeze", "hold.txt"]) == 0
    assert FreezeFile("hold.txt").chains == [1]
    assert "FREEZE_FILE          : hold.txt" in (tmp_path / "b_report.md").read_text()

    # an existing freeze file is refused like the other outputs, before anything is written
    (tmp_path / "c_freeze.txt").write_text("keep me")
    assert p2p.main(["h.pdb", "-o", "c.pimms"]) == 1
    assert "c_freeze.txt already exists; give --overwrite" in capsys.readouterr().err
    assert not (tmp_path / "c.pimms").exists() and not (tmp_path / "c_report.md").exists()
    assert (tmp_path / "c_freeze.txt").read_text() == "keep me"

    # and it cannot be another output
    assert p2p.main(["h.pdb", "-o", "d.pimms", "--freeze", "d.pimms"]) == 1
    assert "the restart file and the freeze file cannot be the same file" in capsys.readouterr().err


def test_main_reports_a_failed_conversion_with_status_1(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "water.pdb", [atom("O", "HOH", "A", 1, (0, 0, 0), record="HETATM")])
    assert p2p.main(["water.pdb"]) == 1
    assert "nothing to convert" in capsys.readouterr().err
    assert sorted(p.name for p in tmp_path.iterdir()) == ["water.pdb"]


@pytest.mark.parametrize("arguments", [
    ["h.pdb", "--spacing", "0"], ["h.pdb", "--spacing", "nan"], ["h.pdb", "--padding", "-1"],
    ["h.pdb", "--map", "MSE=0"], ["h.pdb", "--chains", ","], [],
    ["h.pdb", "--freeze", "f.txt", "--no-freeze"],
])
def test_main_usage_errors_exit_with_status_2(tmp_path, monkeypatch, arguments):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as stop:
        p2p.main(arguments)
    assert stop.value.code == 2


def _clean_env():
    """An environment with no repo-specific import help."""
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONSAFEPATH", None)
    return env


def test_script_runs_from_a_foreign_directory(tmp_path):
    write_cif(tmp_path / "h.cif", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))

    def run(*arguments):
        return subprocess.run([sys.executable, str(PDB2PIMMS_SCRIPT), *arguments], cwd=str(tmp_path),
                              env=_clean_env(), capture_output=True, text=True, timeout=300, check=False)

    result = run("h.cif")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "h_restart.pimms").exists() and (tmp_path / "h_report.md").exists()
    assert (tmp_path / "h_freeze.txt").exists()
    missing = run("absent.pdb", "-o", "z.pimms")
    assert missing.returncode == 1 and "no such file" in missing.stderr
    usage = run()
    assert usage.returncode == 2 and usage.stdout == ""
    version = run("--version")
    assert version.returncode == 0 and version.stdout.startswith("version ")


def _write_parameter_file(path, bead_types):
    """A parameter file defining every bead type: an angle penalty, a solvation term and weak attractions."""
    lines = ["ANGLE_PENALTY\t%s\t10\t5\t0" % bead for bead in bead_types]
    lines += ["%s\t0\t-1" % bead for bead in bead_types]
    for index, (first, second) in enumerate(itertools.combinations_with_replacement(bead_types, 2)):
        lines.append("%s\t%s\t%d" % (first, second, -1 - index % 3))
    path.write_text("\n".join(lines) + "\n")


def test_pimms_runs_from_a_converted_structure(tmp_path):
    """The whole point: PIMMS accepts the restart file and simulates it, with every energy checked."""
    atoms = protein_atoms("A", helix_trace(14), (PROTEIN_12 + ["MSE", "HID"]))
    atoms += nucleic_atoms("B", helix_trace(6, origin=(20, 0, 0), radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE,
                                            twist=NUCLEIC_TWIST_DEGREES), ["G", "C", "PSU", "A", "U", "C"])
    atoms.append(atom("ZN", "ZN", "A", 101, (0, 12, 0), record="HETATM", element="ZN"))
    write_pdb(tmp_path / "complex.pdb", atoms)
    result = subprocess.run([sys.executable, str(PDB2PIMMS_SCRIPT), "complex.pdb", "--hetero", "-o", "start.pimms"],
                            cwd=str(tmp_path), env=_clean_env(), capture_output=True, text=True, timeout=300,
                            check=False)
    assert result.returncode == 0, result.stderr

    with open(tmp_path / "start.pimms", "rb") as handle:
        payload = pickle.load(handle)
    sequences = [chain[1] for chain in payload["CHAINS"].values()]
    assert sequences == ["MAGHKELSWDVRMH", "gcuauc", "b"]
    bead_types = sorted(set("".join(sequences)))
    _write_parameter_file(tmp_path / "params.prm", bead_types)

    # a larger periodic box than the restart file's: PIMMS centres the restart box in it
    box = [d + 6 for d in payload["DIMENSIONS"]]
    (tmp_path / "KEYFILE.kf").write_text("\n".join([
        "DIMENSIONS : %s" % " ".join(str(d) for d in box), "PARAMETER_FILE : params.prm",
        "RESTART_FILE : start.pimms", "SEED : 3", "TEMPERATURE : 50", "HARDWALL : False",
        "N_STEPS : 4", "EQUILIBRATION : 0", "MOVE_CRANKSHAFT : 0.7", "MOVE_CHAIN_TRANSLATE : 0.3",
        "CRANKSHAFT_SUBSTEPS : 500", "ENERGY_CHECK : 1", "EN_FREQ : 1"]) + "\n")
    run = subprocess.run([sys.executable, str(PIMMS_SCRIPT), "-k", "KEYFILE.kf"], cwd=str(tmp_path),
                         env=_clean_env(), capture_output=True, text=True, timeout=600, check=False)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    assert "Simulation complete" in run.stdout
    start_beads = [line for line in (tmp_path / "START.pdb").read_text().splitlines() if line.startswith("ATOM")]
    assert len(start_beads) == sum(len(s) for s in sequences)


def test_pimms_holds_the_frozen_structure_still_while_added_chains_move(tmp_path):
    """The freeze file pins every converted chain; EXTRA_CHAIN chains are numbered after them and move."""
    _two_helices_and_an_ion(tmp_path / "pair.pdb")
    result = subprocess.run([sys.executable, str(PDB2PIMMS_SCRIPT), "pair.pdb", "--hetero"], cwd=str(tmp_path),
                            env=_clean_env(), capture_output=True, text=True, timeout=300, check=False)
    assert result.returncode == 0, result.stderr
    with open(tmp_path / "pair_restart.pimms", "rb") as handle:
        start = pickle.load(handle)
    added_sequence = "KEKE"
    _write_parameter_file(tmp_path / "params.prm", sorted(set("".join(c[1] for c in start["CHAINS"].values()))
                                                          | set(added_sequence)))

    n_steps = 20
    box = [d + 10 for d in start["DIMENSIONS"]]
    (tmp_path / "KEYFILE.kf").write_text("\n".join([
        "DIMENSIONS : %s" % " ".join(str(d) for d in box), "PARAMETER_FILE : params.prm",
        "RESTART_FILE : pair_restart.pimms", "FREEZE_FILE : pair_freeze.txt",
        "EXTRA_CHAIN : 6 %s" % added_sequence, "SEED : 3", "TEMPERATURE : 50", "HARDWALL : False",
        "N_STEPS : %d" % n_steps, "EQUILIBRATION : 0", "MOVE_CRANKSHAFT : 0.6", "MOVE_CHAIN_TRANSLATE : 0.4",
        "CRANKSHAFT_SUBSTEPS : 200", "ENERGY_CHECK : 1", "EN_FREQ : 1", "RESTART_FREQ : %d" % n_steps]) + "\n")
    run = subprocess.run([sys.executable, str(PIMMS_SCRIPT), "-k", "KEYFILE.kf"], cwd=str(tmp_path),
                         env=_clean_env(), capture_output=True, text=True, timeout=600, check=False)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    assert "Chains to freeze : [1, 2, 3]" in run.stdout

    with open(tmp_path / "restart.pimms", "rb") as handle:
        end = pickle.load(handle)
    # the restart box is centred in the larger one, so every frozen bead moves by the same offset and no more
    offset = (np.array(box) - np.array(start["DIMENSIONS"])) // 2
    for chain_id, (positions, sequence, _) in start["CHAINS"].items():
        assert np.array_equal(np.array(end["CHAINS"][chain_id][0]), np.array(positions) + offset)
        assert end["CHAINS"][chain_id][1] == sequence
    assert sorted(end["CHAINS"]) == list(range(1, 10))
    assert all(end["CHAINS"][chain_id][1] == added_sequence for chain_id in range(4, 10))


def _shortened_keyfile(text):
    """A demo keyfile with DEMO_SMOKE_SETTINGS in place of its own values; every one must be replaced."""
    replaced = set()
    lines = []
    for line in text.splitlines():
        keyword = line.split(":")[0].strip()
        if keyword in DEMO_SMOKE_SETTINGS:
            line = "%s : %s" % (keyword, DEMO_SMOKE_SETTINGS[keyword])
            replaced.add(keyword)
        lines.append(line)
    assert replaced == set(DEMO_SMOKE_SETTINGS)
    return "\n".join(lines) + "\n"


def test_the_nucleosome_demo_runs_as_its_readme_says(tmp_path):
    """The demo's commands, each run cut to a few steps, get from the structure to the figure."""
    pytest.importorskip("matplotlib")
    for name in NUCLEOSOME_DEMO_INPUTS:
        shutil.copy(NUCLEOSOME_DEMO / name, tmp_path / name)

    def run(*command, cwd=tmp_path):
        result = subprocess.run(list(command), cwd=str(cwd), env=_clean_env(), capture_output=True, text=True,
                                timeout=600, check=False)
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
        return result

    # the conversion command is the README's, which spreads it over lines ending in a backslash
    readme = (NUCLEOSOME_DEMO / "README.md").read_text().replace("\\\n", " ")
    mobile = [argument for tail in NUCLEOSOME_DEMO_TAILS for argument in ("--mobile", tail)]
    assert " ".join(["pdb2pimms", "1kx5.cif.gz"] + mobile) in " ".join(readme.split())
    run(sys.executable, str(PDB2PIMMS_SCRIPT), "1kx5.cif.gz", *mobile)
    for name in NUCLEOSOME_DEMO_RUNS:
        (tmp_path / name).mkdir()
        (tmp_path / name / "KEYFILE.kf").write_text(_shortened_keyfile((NUCLEOSOME_DEMO / name / "KEYFILE.kf")
                                                                       .read_text()))
        run(sys.executable, str(PIMMS_SCRIPT), "-k", "KEYFILE.kf", cwd=tmp_path / name)
    analysis = run(sys.executable, "analyze_contacts.py")
    assert "ProT\u03b1 (with H1.0 CTD) in H1_CTD_ProTa: 4 copies of 111 residues, 3 production frames" in analysis.stdout
    assert "H2A C-terminal   C:118-128" in analysis.stdout
    assert (tmp_path / "nucleosome_contacts.pdf").stat().st_size > 0
    for name, protein in (("H1_CTD", "H1_CTD"), ("ProTa", "ProTa"), ("H1_CTD_ProTa", "H1_CTD"),
                          ("H1_CTD_ProTa", "ProTa")):
        table = (tmp_path / name / ("contacts_nucleosome_%s.csv" % protein)).read_text().splitlines()
        assert len(table) == 1 + 1268

    # reset.sh takes the directory back to its inputs, so the demo runs again from the start
    reset = run(str(tmp_path / "reset.sh"))
    assert "reset.sh: deleted" in reset.stdout and reset.stderr == ""
    left = sorted(str(path.relative_to(tmp_path)) for path in tmp_path.rglob("*") if path.is_file())
    assert left == sorted(list(NUCLEOSOME_DEMO_INPUTS) + ["%s/KEYFILE.kf" % name for name in NUCLEOSOME_DEMO_RUNS])
    run(sys.executable, str(PDB2PIMMS_SCRIPT), "1kx5.cif.gz", *mobile)


# ---------------------------------------------------------------------------
# regressions from the independent review
# ---------------------------------------------------------------------------

def test_simulation_names_that_are_other_molecules_codes_stay_hetero_groups(tmp_path):
    # RU (a ruthenium ion), GUA (glutaric acid), T3 (a thyroid hormone) and SOL
    # (L-sorbose) are PDB Chemical Component Dictionary codes; without a polymer
    # backbone they must not become nucleotides or water
    atoms = protein_atoms("A", helix_trace(5), PROTEIN_12[:5])
    atoms += [atom("RU", "RU", "A", 101, (10, 0, 0), record="HETATM", element="RU")]
    atoms += [atom("C1", "GUA", "A", 102, (0, 10, 0), record="HETATM"),
              atom("O1", "GUA", "A", 102, (1.2, 10, 0), record="HETATM")]
    atoms += [atom("C1", "T3", "A", 103, (0, 0, 15), record="HETATM"),
              atom("I1", "T3", "A", 103, (2, 0, 15), record="HETATM", element="I")]
    atoms += [atom("C1", "SOL", "A", 104, (10, 10, 0), record="HETATM"),
              atom("O1", "SOL", "A", 104, (11.4, 10, 0), record="HETATM")]
    # and a GROMACS water called SOL, which is water
    atoms += [atom("OW", "SOL", "W", 1, (20, 20, 20), element="O"),
              atom("HW1", "SOL", "W", 1, (20.9, 20, 20), element="H")]
    write_pdb(tmp_path / "ligands.pdb", atoms)
    conversion = convert(tmp_path / "ligands.pdb")
    assert [s.sequence for s in conversion.segments] == ["MAGHK"]
    assert conversion.excluded.hetero == {"RU": 1, "GUA": 1, "T3": 1, "SOL": 1}
    assert conversion.excluded.waters == 1


def test_simulation_names_with_their_backbone_are_standard_residues(tmp_path):
    trace = helix_trace(4, radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE, twist=NUCLEIC_TWIST_DEGREES)
    atoms = protein_atoms("A", helix_trace(4, origin=(0, 0, -40)), ["HSD", "HISE", "LYSN", "CYS2"])
    atoms += nucleic_atoms("B", trace, ["RA5", "ADE", "URA", "RU3"], rna=True)
    atoms += nucleic_atoms("C", trace + 40.0, ["DA5", "T", "THY", "DG3"], rna=False)
    write_pdb(tmp_path / "md.pdb", atoms)
    conversion = convert(tmp_path / "md.pdb")
    assert [s.sequence for s in conversion.segments] == ["HHKC", "aauu", "attg"]
    assert {r.parent_source for s in conversion.segments for r in s.residues} == {"other name"}


def test_the_file_parent_comes_before_a_simulation_name(tmp_path):
    atoms = protein_atoms("A", helix_trace(3), ["ALA", "HSE", "ALA"])
    write_pdb(tmp_path / "hse.pdb", atoms, header_lines=["MODRES TEST HSE A    2  SER  HOMOSERINE"])
    residue = convert(tmp_path / "hse.pdb").segments[0].residues[1]
    assert (residue.code, residue.parent, residue.parent_source) == ("S", "SER", "file")


def test_an_isolated_hetatm_amino_acid_is_a_ligand_but_a_bonded_one_is_kept(tmp_path):
    atoms = protein_atoms("A", helix_trace(5), ["MET", "MSE", "GLY", "HIS", "LYS"])
    for a in atoms:
        if a["resname"] == "MSE":
            a["record"] = "HETATM"
    atoms += protein_atoms("A", helix_trace(1, origin=(15, 0, 0)), ["ARG"], start=200, record="HETATM")
    write_pdb(tmp_path / "free.pdb", atoms)
    conversion = convert(tmp_path / "free.pdb")
    assert [s.sequence for s in conversion.segments] == ["MMGHK"]
    assert conversion.isolated_hetatm == ["A:ARG200"]
    assert conversion.excluded.hetero == {"ARG": 1}
    assert not conversion.breaks
    kept = convert(tmp_path / "free.pdb", keep_hetero=True)
    assert [s.sequence for s in kept.segments] == ["MMGHK", "b"]


def test_alternate_location_blocks_written_apart_still_make_one_residue_each(tmp_path):
    atoms = protein_atoms("A", helix_trace(5), PROTEIN_12[:5])
    block = {3: [], 4: []}
    for a in [a for a in atoms if a["resseq"] in (3, 4)]:
        block[a["resseq"]].append(a)
        atoms.remove(a)
    interleaved = []
    for altloc, shift in (("A", 0.0), ("B", 0.3)):
        for resseq in (3, 4):
            for a in block[resseq]:
                copy = dict(a, altloc=altloc, occupancy=0.5, xyz=a["xyz"] + shift)
                interleaved.append(copy)
    position = next(i for i, a in enumerate(atoms) if a["resseq"] == 5)
    atoms[position:position] = interleaved
    write_pdb(tmp_path / "blocks.pdb", atoms)
    conversion = convert(tmp_path / "blocks.pdb")
    assert [s.sequence for s in conversion.segments] == ["MAGHK"]
    assert [r.label for r in conversion.excluded.altloc_residues] == ["A:GLY3", "A:HIS4"]


def test_residues_that_only_share_a_number_stay_apart(tmp_path):
    # residue numbers that recur later in the file (wrapped numbering), with no alternate locations
    atoms = (protein_atoms("A", helix_trace(3), ["ALA", "GLY", "SER"])
             + protein_atoms("A", helix_trace(3, origin=(20, 0, 0)), ["LYS", "LYS", "LYS"]))
    write_pdb(tmp_path / "wrapped.pdb", atoms)
    conversion = convert(tmp_path / "wrapped.pdb")
    assert [s.sequence for s in conversion.segments] == ["AGS", "KKK"]


def test_frames_separated_by_end_records_are_separate_models(tmp_path):
    frame = [pdb_line(i + 1, a) for i, a in enumerate(protein_atoms("A", helix_trace(4), PROTEIN_12[:4]))]
    (tmp_path / "frames.pdb").write_text("\n".join(frame + ["END"] + frame + ["END"]) + "\n")
    assert read_structure(tmp_path / "frames.pdb").models == [1, 2]
    conversion = convert(tmp_path / "frames.pdb")
    assert [s.sequence for s in conversion.segments] == ["MAGH"]
    assert conversion.excluded.other_models == 1
    assert [s.sequence for s in convert(tmp_path / "frames.pdb", model=2).segments] == ["MAGH"]


def test_an_unknown_ca_only_model_is_converted(tmp_path):
    write_pdb(tmp_path / "unk.pdb", protein_atoms("A", helix_trace(6), ["UNK"] * 6, backbone=False))
    conversion = convert(tmp_path / "unk.pdb")
    assert [s.sequence for s in conversion.segments] == ["XXXXXX"]
    assert any("UNK" in warning for warning in conversion.warnings)


def test_a_phosphate_only_trace_stays_one_chain(tmp_path):
    trace = helix_trace(6, radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE, twist=NUCLEIC_TWIST_DEGREES)
    atoms = [atom("P", name, "A", i + 1, xyz) for i, (name, xyz) in enumerate(zip("GCAUGC", trace))]
    write_pdb(tmp_path / "ptrace.pdb", atoms)
    conversion = convert(tmp_path / "ptrace.pdb")
    assert [s.sequence for s in conversion.segments] == ["gcaugc"]


def test_mmcif_waters_without_author_numbering_are_counted_one_by_one(tmp_path):
    lines = ["data_X", "loop_"] + ["_atom_site.%s" % tag for tag in (
        "group_PDB", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id",
        "Cartn_x", "Cartn_y", "Cartn_z")]
    for i, a in enumerate(protein_atoms("A", helix_trace(3), PROTEIN_12[:3])):
        x, y, z = a["xyz"]
        lines.append("ATOM %s %s %s A %d %.3f %.3f %.3f" % (a["element"], a["name"], a["resname"], a["resseq"], x, y, z))
    for i in range(3):
        lines.append("HETATM O O HOH B . %.3f 20.000 20.000" % (3.0 * i))
        lines.append("HETATM H H1 HOH B . %.3f 20.900 20.000" % (3.0 * i))
    (tmp_path / "label.cif").write_text("\n".join(lines) + "\n")
    conversion = convert(tmp_path / "label.cif")
    assert conversion.excluded.waters == 3
    assert [s.sequence for s in conversion.segments] == ["MAG"]


def test_broken_gzip_and_non_finite_numbers_are_refused(tmp_path):
    write_pdb(tmp_path / "good.pdb.gz", protein_atoms("A", helix_trace(40), (PROTEIN_12 * 4)[:40]))
    data = (tmp_path / "good.pdb.gz").read_bytes()
    (tmp_path / "cut.pdb.gz").write_bytes(data[:len(data) // 2])
    with pytest.raises(StructureFileException, match="Cannot read structure file"):
        read_structure(tmp_path / "cut.pdb.gz")

    line = pdb_line(1, atom("CA", "ALA", "A", 1, (1.0, 2.0, 3.0)))
    (tmp_path / "nan.pdb").write_text(line[:30] + "     nan" + line[38:] + "\n")
    with pytest.raises(StructureFileException, match="not a finite number"):
        read_structure(tmp_path / "nan.pdb")
    (tmp_path / "inf.cif").write_text("data_X\nloop_\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n"
                                      "_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n"
                                      "_atom_site.Cartn_y\n_atom_site.Cartn_z\nCA ALA A 1 inf 0 0\n")
    with pytest.raises(StructureFileException, match="not a finite number"):
        read_structure(tmp_path / "inf.cif")


def test_a_byte_order_mark_does_not_hide_the_first_record(tmp_path):
    atoms = protein_atoms("A", helix_trace(3), PROTEIN_12[:3])
    lines = [pdb_line(i + 1, a) for i, a in enumerate(atoms)]
    (tmp_path / "bom.pdb").write_text("﻿" + "\n".join(lines) + "\n", encoding="utf-8")
    assert len(read_structure(tmp_path / "bom.pdb").atoms) == len(atoms)


def test_a_blank_chain_is_selected_with_a_dash(tmp_path):
    atoms = (protein_atoms("", helix_trace(3), PROTEIN_12[:3])
             + protein_atoms("A", helix_trace(3, origin=(20, 0, 0)), PROTEIN_12[3:6]))
    write_pdb(tmp_path / "blank.pdb", atoms)
    conversion = convert(tmp_path / "blank.pdb", chains=["-"])
    assert [s.sequence for s in conversion.segments] == ["MAG"]
    assert conversion.excluded.unselected_chains == ["A"]
    with pytest.raises(StructureConversionException, match=r"it has -  ?\(blank\), A"):
        convert(tmp_path / "blank.pdb", chains=["Q"])


def test_format_sniffing_looks_past_header_records(tmp_path):
    atoms = protein_atoms("A", helix_trace(3), PROTEIN_12[:3])
    header = ["HELIX    1   1 ALA A    1  GLY A    3  1", "SEQRES   1 A    3  MET ALA GLY"]
    (tmp_path / "plain").write_text("\n".join(header + [pdb_line(i + 1, a) for i, a in enumerate(atoms)]) + "\n")
    assert read_structure(tmp_path / "plain").file_format == "PDB"


def test_the_mean_occupancy_decides_the_alternate_location(tmp_path):
    atoms = protein_atoms("A", helix_trace(3), PROTEIN_12[:3])
    middle = [a for a in atoms if a["resseq"] == 2]
    for a in middle:
        a.update(altloc="A", occupancy=0.4)
    ca = next(a for a in middle if a["name"] == "CA")
    atoms.insert(atoms.index(ca) + 1, dict(ca, altloc="B", occupancy=0.6, xyz=ca["xyz"] + 0.5))
    write_pdb(tmp_path / "occ.pdb", atoms)
    residue = convert(tmp_path / "occ.pdb").segments[0].residues[1]
    # A covers three atoms at 0.4 (a larger sum), B one atom at 0.6 (the larger occupancy)
    assert residue.altloc == "B"


def test_poor_fit_warning_threshold_is_one_and_a_half_times_the_rounding_floor():
    def statistics(rmsd):
        return p2p.FitStatistics(rmsd=rmsd, segment_rmsd=[rmsd], max_deviation=rmsd, max_deviation_label="A:ALA1",
                                 bond_steps=p2p.Counter(), input_bond_lengths=np.zeros(0),
                                 lattice_bond_lengths=np.zeros(0), input_rg=[0.0], lattice_rg=[0.0])

    def warned(rmsd, spacing):
        warnings = p2p._conversion_warnings([], [p2p.Residue("A", "1", "", "ALA", "ATOM", {}, [],
                                                             kind=p2p.KIND_AMINO_ACID)],
                                            statistics(rmsd), [10, 10, 10], spacing)
        return any("less well than it could" in w for w in warnings)

    assert not warned(0.74 * 4.0, 4.0)
    assert warned(0.76 * 4.0, 4.0)


# ---------------------------------------------------------------------------
# --mobile: residues the freeze file leaves free
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text, chain, first, last", [
    ("A", "A", None, None), ("A:5", "A", 5, 5), ("A:1-37", "A", 1, 37), ("I:-73--60", "I", -73, -60),
    ("-:2-4", "", 2, 4), (" B:-5-3 ", "B", -5, 3)])
def test_mobile_selections_are_read(text, chain, first, last):
    selection = p2p.parse_mobile_selection(text)
    assert (selection.chain, selection.first, selection.last) == (chain, first, last)


@pytest.mark.parametrize("text", [":1-3", "A:", "A:x", "A:5-2", "A:1-2-3", "A:1.5", "A:+3"])
def test_malformed_mobile_selections_are_refused(text):
    with pytest.raises(ValueError, match="--mobile"):
        p2p.parse_mobile_selection(text)


def _mobile_conversion(tmp_path, *selections):
    """The two helices and the ion (three chains, with --hetero) with --mobile selections."""
    _two_helices_and_an_ion(tmp_path / "pair.pdb")
    return convert(tmp_path / "pair.pdb", keep_hetero=True,
                   mobile=[p2p.parse_mobile_selection(selection) for selection in selections])


def test_the_freeze_file_leaves_the_mobile_residues_free_and_pimms_reads_it(tmp_path):
    conversion = _mobile_conversion(tmp_path, "A:1-3", "B")
    freeze_path = tmp_path / "pair_freeze.txt"
    p2p.write_outputs(conversion, tmp_path / "pair_restart.pimms", "report", tmp_path / "pair_report.md",
                      freeze_path)

    text = freeze_path.read_text()
    assert "except the residues given with --mobile\n# (A:1-3 B)" in text
    assert "B 1 3-9  # A:MET1-ASP10: 7 of 10 bead(s) frozen; free: beads 0-2 (A:MET1-GLY3)\n" in text
    assert "# chain 2 (B:MET1-ASP10, 10 bead(s)) is left free\n" in text
    assert "C 3  # A:ZN101, 1 bead(s)\n" in text
    freeze = FreezeFile(str(freeze_path))
    assert freeze.resolve({1: 10, 2: 10, 3: 1}) == ([1, 3], {1: (3, 4, 5, 6, 7, 8, 9)})


def test_the_report_lists_the_free_residues_and_asks_for_the_crankshaft(tmp_path):
    conversion = _mobile_conversion(tmp_path, "A:1-3", "B")
    text = p2p.render_report(conversion, tmp_path / "pair_restart.pimms", tmp_path / "pair_report.md",
                             "pdb2pimms pair.pdb", "9.9.9", "2026-01-01 00:00", tmp_path / "pair_freeze.txt")
    assert "MOVE_CRANKSHAFT      : <fraction>   # the only move for a partly frozen chain's free beads" in text
    assert "A chain whose residues are all free is left out of the file and moves like any other chain." in text
    assert "EXTRA_CHAIN          : <count> <sequence>   # mobile chains (optional here)" in text
    assert "(`A:1-3`, `B`): 13 of its 21 beads are free" in text
    assert "- chain 1 (A:MET1-ASP10): 3 of 10 beads free - beads 0-2 (A:MET1-GLY3)" in text
    assert "- chain 2 (B:MET1-ASP10): every bead free - beads 0-9 (B:MET1-ASP10)" in text


def test_a_chain_freed_whole_is_an_ordinary_mobile_chain(tmp_path, capsys):
    """With no partly frozen chain, nothing asks for the crankshaft."""
    conversion = _mobile_conversion(tmp_path, "B")
    text = p2p.render_report(conversion, tmp_path / "pair_restart.pimms", tmp_path / "pair_report.md",
                             "pdb2pimms pair.pdb", "9.9.9", "2026-01-01 00:00", tmp_path / "pair_freeze.txt")
    assert "MOVE_CRANKSHAFT" not in text
    assert "A chain whose residues are all free is left out of the file and moves like any other chain." in text
    p2p._print_summary(conversion, tmp_path / "pair_restart.pimms", tmp_path / "pair_report.md",
                       tmp_path / "pair_freeze.txt")
    printed = capsys.readouterr().out
    assert "freezes 11 of 21 beads (2 chain(s) entirely, 0 partly, 1 left free)\n" in printed
    assert "MOVE_CRANKSHAFT" not in printed


def test_a_selection_on_a_blank_chain_is_given_with_an_equals_sign():
    """argparse takes a value that starts with '-' for an option, so the help says to write --mobile=-:1-5."""
    arguments = p2p.build_argument_parser().parse_args(["x.pdb", "--mobile=-:1-5", "--chains=-,A"])
    assert arguments.mobile == ["-:1-5"] and arguments.chains == "-,A"
    assert p2p.parse_mobile_selection(arguments.mobile[0]) == p2p.MobileSelection("", 1, 5, "-:1-5")


def test_mobile_selections_that_match_nothing_or_everything_are_refused(tmp_path):
    with pytest.raises(StructureConversionException, match="--mobile C:1-3 matches no converted residue"):
        _mobile_conversion(tmp_path, "C:1-3")
    with pytest.raises(StructureConversionException, match="--mobile A:50-60 matches no converted residue"):
        _mobile_conversion(tmp_path, "A:50-60")
    # A holds the zinc as well as the first helix, so A and B are everything
    with pytest.raises(StructureConversionException, match="leaves every converted residue free"):
        _mobile_conversion(tmp_path, "A", "B")


def test_a_wrong_freeze_file_with_mobile_residues_is_never_put_in_place(tmp_path, monkeypatch):
    conversion = _mobile_conversion(tmp_path, "A:1-3")
    # a writer that also freezes bead 2, which --mobile A:1-3 leaves free
    monkeypatch.setattr(p2p, "_partial_freeze_file_text", lambda conversion, restart_path, flags: "B 1 2-9\nC 2 3\n")
    with pytest.raises(StructureConversionException, match="reads frozen chains"):
        p2p.write_outputs(conversion, tmp_path / "r.pimms", "report", tmp_path / "r_report.md",
                          tmp_path / "r_freeze.txt")
    # and one PIMMS refuses
    monkeypatch.setattr(p2p, "_partial_freeze_file_text", lambda conversion, restart_path, flags: "B 1 9-3\n")
    with pytest.raises(StructureConversionException, match="refused the freeze file"):
        p2p.write_outputs(conversion, tmp_path / "r.pimms", "report", tmp_path / "r_report.md",
                          tmp_path / "r_freeze.txt")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["pair.pdb"]


def test_main_mobile_option(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    assert p2p.main(["h.pdb", "--mobile", "A:1-2"]) == 0
    assert ("freezes 6 of 8 beads (0 chain(s) entirely, 1 partly, 0 left free); the free beads of a partly "
            "frozen chain move by MOVE_CRANKSHAFT only") in capsys.readouterr().out
    assert FreezeFile("h_freeze.txt").beads == [(1, 2, 7)]

    # a selection that matches nothing fails the conversion and writes nothing
    assert p2p.main(["h.pdb", "-o", "z.pimms", "--mobile", "B:1-2"]) == 1
    assert "--mobile B:1-2 matches no converted residue" in capsys.readouterr().err
    assert not (tmp_path / "z.pimms").exists() and not (tmp_path / "z_freeze.txt").exists()


@pytest.mark.parametrize("arguments", [
    ["h.pdb", "--mobile", "A:3-1"], ["h.pdb", "--mobile", "A:"], ["h.pdb", "--mobile", "A:1", "--no-freeze"],
])
def test_main_mobile_usage_errors_exit_with_status_2(tmp_path, monkeypatch, arguments):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as stop:
        p2p.main(arguments)
    assert stop.value.code == 2


def test_pimms_moves_only_the_mobile_residues(tmp_path):
    """Every chain holds a frozen bead, so the crankshaft alone moves the three free residues."""
    _two_helices_and_an_ion(tmp_path / "pair.pdb")
    result = subprocess.run([sys.executable, str(PDB2PIMMS_SCRIPT), "pair.pdb", "--hetero", "--mobile", "A:1-3"],
                            cwd=str(tmp_path), env=_clean_env(), capture_output=True, text=True, timeout=300,
                            check=False)
    assert result.returncode == 0, result.stderr
    with open(tmp_path / "pair_restart.pimms", "rb") as handle:
        start = pickle.load(handle)
    _write_parameter_file(tmp_path / "params.prm", sorted(set("".join(c[1] for c in start["CHAINS"].values()))))

    n_steps = 40
    box = [d + 10 for d in start["DIMENSIONS"]]
    (tmp_path / "KEYFILE.kf").write_text("\n".join([
        "DIMENSIONS : %s" % " ".join(str(d) for d in box), "PARAMETER_FILE : params.prm",
        "RESTART_FILE : pair_restart.pimms", "FREEZE_FILE : pair_freeze.txt", "SEED : 3", "TEMPERATURE : 50",
        "HARDWALL : False", "N_STEPS : %d" % n_steps, "EQUILIBRATION : 0", "MOVE_CRANKSHAFT : 1.0",
        "CRANKSHAFT_SUBSTEPS : 200", "ENERGY_CHECK : 1", "EN_FREQ : 1", "RESTART_FREQ : %d" % n_steps]) + "\n")
    run = subprocess.run([sys.executable, str(PIMMS_SCRIPT), "-k", "KEYFILE.kf"], cwd=str(tmp_path),
                         env=_clean_env(), capture_output=True, text=True, timeout=600, check=False)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    assert "Beads to freeze  : chain 1, beads 3-9" in run.stdout

    with open(tmp_path / "restart.pimms", "rb") as handle:
        end = pickle.load(handle)
    offset = (np.array(box) - np.array(start["DIMENSIONS"])) // 2
    first_chain = np.array(start["CHAINS"][1][0]) + offset
    moved = np.any(np.array(end["CHAINS"][1][0]) != first_chain, axis=1)
    assert not moved[3:].any() and moved[:3].any()
    for chain_id in (2, 3):
        assert np.array_equal(np.array(end["CHAINS"][chain_id][0]), np.array(start["CHAINS"][chain_id][0]) + offset)
    assert end["FROZEN_BEADS"] == {1: [3, 4, 5, 6, 7, 8, 9]}


# ---------------------------------------------------------------------------
# regressions from the second audit
# ---------------------------------------------------------------------------

def test_calcium_ions_are_hetero_groups_not_ca_only_residues(tmp_path):
    """A calcium ion is a lone atom named CA (residue CA); it must not pass for a CA-only amino acid."""
    protein = protein_atoms("A", helix_trace(4), PROTEIN_12[:4])
    # a binuclear site (3.9 A apart, so the CA-CA test would bond them) and an ion written as an ATOM record
    ions = [atom("CA", "CA", "A", 201, (20.0, 0, 0), record="HETATM", element="CA"),
            atom("CA", "CA", "A", 202, (23.9, 0, 0), record="HETATM", element="CA"),
            atom("CA", "CA", "I", 1, (0, 25.0, 0), element="CA")]
    write_pdb(tmp_path / "ca.pdb", protein + ions)
    write_cif(tmp_path / "ca.cif", protein + ions)
    for name in ("ca.pdb", "ca.cif"):
        conversion = convert(tmp_path / name)
        assert [s.sequence for s in conversion.segments] == ["MAGH"]
        assert conversion.excluded.hetero == {"CA": 3}
        assert not any("not standard" in w for w in conversion.warnings)
    # with --hetero each is a one-bead hetero chain, as a zinc ion is
    assert [s.sequence for s in convert(tmp_path / "ca.pdb", keep_hetero=True).segments] == ["MAGH", "b", "b", "b"]

    # a file without element columns: the residue name CA gives the ion away
    lines = [pdb_line(i + 1, a)[:76] for i, a in enumerate(protein + ions)]
    (tmp_path / "noelement.pdb").write_text("\n".join(lines) + "\nEND\n")
    conversion = convert(tmp_path / "noelement.pdb")
    assert [s.sequence for s in conversion.segments] == ["MAGH"]
    assert conversion.excluded.hetero == {"CA": 3}
    # and a CA-only model without element columns is still converted
    lines = [pdb_line(i + 1, a)[:76] for i, a in enumerate(protein_atoms("A", helix_trace(4), ["UNK"] * 4,
                                                                         backbone=False))]
    (tmp_path / "unk.pdb").write_text("\n".join(lines) + "\nEND\n")
    assert [s.sequence for s in convert(tmp_path / "unk.pdb").segments] == ["XXXX"]


def test_hydrogens_listed_after_the_chain_or_all_named_h_never_make_beads(tmp_path):
    """Open Babel lists every residue's hydrogens after the whole chain and names them all H."""
    heavy = protein_atoms("A", helix_trace(4), PROTEIN_12[:4])
    hydrogens = [atom("H", resname, "A", i + 1, xyz, element="H")
                 for i, (resname, ca) in enumerate(zip(PROTEIN_12[:4], helix_trace(4)))
                 for xyz in (ca + (0.0, 0.0, 1.0), ca + (0.0, 1.0, 0.0), ca + (1.0, 0.0, 0.0))]
    write_pdb(tmp_path / "after.pdb", heavy + hydrogens)
    write_cif(tmp_path / "after.cif", heavy + hydrogens)
    # the same hydrogens interleaved with their residues, still all named H
    interleaved = []
    for i in range(4):
        interleaved += [a for a in heavy if a["resseq"] == i + 1] + hydrogens[3 * i:3 * i + 3]
    write_pdb(tmp_path / "interleaved.pdb", interleaved)

    for name in ("after.pdb", "after.cif", "interleaved.pdb"):
        conversion = convert(tmp_path / name)
        assert [s.sequence for s in conversion.segments] == ["MAGH"]
        assert len(conversion.residues) == 4
        assert conversion.excluded.hydrogens == 12
        assert not conversion.breaks and not conversion.warnings


def test_a_molecule_split_by_a_repeated_atom_name_is_warned_about(tmp_path):
    """A ligand whose atoms are named by element alone is split on each repeated name; the user is told."""
    protein = protein_atoms("A", helix_trace(4), PROTEIN_12[:4])
    ligand = [atom(name, "UNL", "A", 900, xyz, record="HETATM")
              for name, xyz in (("C", (0, 20, 0)), ("C", (1.5, 20, 0)), ("O", (2.2, 21, 0)), ("O", (2.2, 19, 0)))]
    write_pdb(tmp_path / "unl.pdb", protein + ligand)
    kept = convert(tmp_path / "unl.pdb", keep_hetero=True)
    assert any("An atom name repeats within residue(s) A:UNL900 (2 extra bead(s))" in w for w in kept.warnings)
    assert [s.sequence for s in kept.segments] == ["MAGH", "b", "b", "b"]
    # without --hetero no bead comes of it, so there is nothing to warn about
    assert not any("repeats" in w for w in convert(tmp_path / "unl.pdb").warnings)


def test_a_nucleotide_ligand_without_a_linking_phosphorus_is_not_joined_to_a_chain(tmp_path):
    """The CCD names ATP's alpha phosphorus PA, so only the generous C4'-C4' distance could bond it; it must not."""
    trace = helix_trace(5, radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE, twist=NUCLEIC_TWIST_DEGREES)
    rna = nucleic_atoms("B", trace[:4], ["G", "C", "A", "U"])

    def triphosphate(resname, chain, resseq, c4):
        names = [("PG", (-5.0, 0, 0)), ("PB", (-3.6, 0, 0)), ("PA", (-2.4, 0, 0)), ("O5'", (-1.4, 0, 0)),
                 ("C4'", (0, 0, 0)), ("C2'", (0, 1.5, 0)), ("O2'", (0, 2.5, 0.8)), ("O3'", (1.4, 0, 0)),
                 ("N9", (0.5, -1.5, 0))]
        return [atom(name, resname, chain, resseq, np.asarray(c4) + offset, record="HETATM") for name, offset in names]

    # an ATP whose C4' sits where a fifth nucleotide's would (6.5 A from the 3' end), and two ADPs 6 A apart
    atp = triphosphate("ATP", "B", 101, trace[4])
    protein = protein_atoms("A", helix_trace(4, origin=(0, 40, 0)), PROTEIN_12[:4])
    adps = triphosphate("ADP", "A", 401, (0, 60, 0)) + triphosphate("ADP", "A", 402, (6, 60, 0))
    write_pdb(tmp_path / "ntp.pdb", rna + atp + protein + adps)
    conversion = convert(tmp_path / "ntp.pdb")
    assert [s.sequence for s in conversion.segments] == ["gcau", "MAGH"]
    assert conversion.excluded.hetero == {"ATP": 1, "ADP": 2}
    assert conversion.isolated_hetatm == ["B:ATP101", "A:ADP401", "A:ADP402"]

    # HETATM residues of coarse models still bond by the backbone distance alone
    ca_trace = [atom("CA", name, "A", i + 1, xyz, record="HETATM" if name == "MSE" else "ATOM")
                for i, (name, xyz) in enumerate(zip(["ALA", "MSE", "GLY"], helix_trace(3)))]
    p_trace = [atom("P", name, "B", i + 1, xyz, record="HETATM" if name == "PSU" else "ATOM")
               for i, (name, xyz) in enumerate(zip(["G", "PSU", "C"], trace[:3] + 40.0))]
    write_pdb(tmp_path / "coarse.pdb", ca_trace + p_trace)
    coarse = convert(tmp_path / "coarse.pdb")
    assert [s.sequence for s in coarse.segments] == ["AMG", "guc"]
    assert not coarse.isolated_hetatm


def test_amber_terminal_protonation_states_are_recognised(tmp_path):
    """AMBER prefixes N or C to a terminal residue's own name, which for His and bridged Cys is a protonation state."""
    atoms = (protein_atoms("A", helix_trace(3), ["NHIE", "ALA", "CHID"])
             + protein_atoms("B", helix_trace(3, origin=(20, 0, 0)), ["NHIP", "GLY", "CHIE"])
             + protein_atoms("C", helix_trace(3, origin=(0, 20, 0)), ["NCYX", "ALA", "CCYX"])
             + protein_atoms("D", helix_trace(2, origin=(20, 20, 0)), ["NHID", "CHIP"]))
    write_pdb(tmp_path / "termini.pdb", atoms)
    conversion = convert(tmp_path / "termini.pdb")
    assert [s.sequence for s in conversion.segments] == ["HAH", "HGH", "CAC", "HH"]
    terminal = [r for s in conversion.segments for r in s.residues if len(r.resname) == 4]
    assert {(r.parent, r.parent_source) for r in terminal} == {("HIS", "other name"), ("CYS", "other name")}
    assert not any("not standard" in w for w in conversion.warnings)
    # like every simulation-package name, it counts only with the backbone it claims
    write_pdb(tmp_path / "ligand.pdb", protein_atoms("A", helix_trace(3), PROTEIN_12[:3])
              + [atom("C1", "NHIE", "A", 101, (0, 20, 0), record="HETATM")])
    assert convert(tmp_path / "ligand.pdb").excluded.hetero == {"NHIE": 1}


def test_a_ligand_with_alternate_locations_numbered_like_a_residue_stays_apart(tmp_path):
    """Alternate locations join an earlier residue of the same identity only if it has alternate locations too."""
    protein = protein_atoms("A", helix_trace(5), PROTEIN_12[:5])
    ligand = [atom(name, "LIG", "A", 3, np.asarray(xyz) + shift, record="HETATM", altloc=altloc, occupancy=0.5)
              for altloc, shift in (("A", 0.0), ("B", 0.4))
              for name, xyz in (("C1", (5.0, 15.0, 0.0)), ("C2", (6.5, 15.0, 0.0)), ("O1", (7.3, 16.0, 0.0)))]
    write_pdb(tmp_path / "lig.pdb", protein + ligand)
    conversion = convert(tmp_path / "lig.pdb", bead_position=p2p.BEAD_POSITION_CENTROID)
    assert conversion.excluded.hetero == {"LIG": 1}
    third = conversion.segments[0].residues[2]
    assert (third.resname, sorted(third.atoms), third.altlocs_present) == ("GLY", ["C", "CA", "N"], ())
    own = np.array([a["xyz"] for a in protein if a["resseq"] == 3])
    assert np.allclose(third.position, own.mean(axis=0), atol=1e-3)
    kept = convert(tmp_path / "lig.pdb", keep_hetero=True)
    assert [s.sequence for s in kept.segments] == ["MAGHK", "b"]
    assert kept.segments[1].residues[0].altlocs_present == ("A", "B")


def test_an_mmcif_file_without_type_symbol_finds_hydrogens_and_water_like_a_pdb_file(tmp_path):
    """Neither format gives elements here; both must guess them from the atom names the same way."""
    atoms = [("N", "GLY", 1, (0.0, 0.0, 0.0)), ("CA", "GLY", 1, (1.45, 0.0, 0.0)), ("C", "GLY", 1, (2.0, 1.4, 0.0)),
             ("O", "GLY", 1, (1.3, 2.4, 0.0)), ("H", "GLY", 1, (-0.5, -0.85, 0.0)), ("HA2", "GLY", 1, (1.8, -0.5, 0.9)),
             ("HA3", "GLY", 1, (1.8, -0.5, -0.9)), ("OH2", "TIP3", 2, (10.0, 10.0, 10.0)),
             ("H1", "TIP3", 2, (10.96, 10.0, 10.0)), ("H2", "TIP3", 2, (9.76, 10.93, 10.0))]
    lines = ["data_X", "loop_"] + ["_atom_site.%s" % tag for tag in (
        "group_PDB", "label_atom_id", "label_comp_id", "label_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y",
        "Cartn_z")]
    lines += ["ATOM %s %s A %d %.3f %.3f %.3f" % (name, resname, resseq, *xyz) for name, resname, resseq, xyz in atoms]
    (tmp_path / "g.cif").write_text("\n".join(lines) + "\n")
    pdb_lines = [pdb_line(i + 1, atom(name, resname, "A", resseq, xyz))[:76]
                 for i, (name, resname, resseq, xyz) in enumerate(atoms)]
    (tmp_path / "g.pdb").write_text("\n".join(pdb_lines) + "\nEND\n")

    heavy_centroid = np.mean([xyz for _, _, _, xyz in atoms[:4]], axis=0)
    for name in ("g.cif", "g.pdb"):
        conversion = convert(tmp_path / name, bead_position=p2p.BEAD_POSITION_CENTROID, keep_hetero=True)
        assert np.allclose(conversion.segments[0].residues[0].position, heavy_centroid, atol=1e-3)
        assert (conversion.excluded.hydrogens, conversion.excluded.waters) == (5, 1)
        assert [s.sequence for s in conversion.segments] == ["G"]
    assert [a.element for a in read_structure(tmp_path / "g.cif").atoms] == list("NCCOHHHOHH")


def test_a_model_number_used_twice_is_refused(tmp_path):
    """Two frames both labelled MODEL 1 (two single-model files joined) must not merge into one doubled model."""
    frame = [pdb_line(i + 1, a) for i, a in enumerate(protein_atoms("A", helix_trace(4), PROTEIN_12[:4]))]
    (tmp_path / "twice.pdb").write_text(
        "\n".join(["MODEL        1"] + frame + ["ENDMDL", "MODEL        1"] + frame + ["ENDMDL", "END"]) + "\n")
    with pytest.raises(StructureFileException, match="second block of atoms for model 1"):
        read_structure(tmp_path / "twice.pdb")
    # frames numbered apart, or separated by END with no MODEL records, are still read
    (tmp_path / "apart.pdb").write_text(
        "\n".join(["MODEL        1"] + frame + ["ENDMDL", "MODEL        2"] + frame + ["ENDMDL", "END"]) + "\n")
    assert read_structure(tmp_path / "apart.pdb").models == [1, 2]


def test_no_output_may_replace_the_input_or_differ_from_another_only_by_case(tmp_path, monkeypatch, capsys):
    """The input is protected even with --overwrite, and names that differ only by case are one file on macOS."""
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "s.pdb", protein_atoms("A", helix_trace(5), PROTEIN_12[:5]))
    before = (tmp_path / "s.pdb").read_bytes()
    for arguments, message in ((["s.pdb", "--report", "s.pdb", "--overwrite"],
                                "the input structure and the report cannot be the same file"),
                               (["s.pdb", "-o", "S.PDB", "--overwrite"],
                                "the input structure and the restart file cannot be the same file"),
                               (["s.pdb", "--report", "notes.md", "--freeze", "NOTES.md"],
                                "the report and the freeze file cannot be the same file")):
        assert p2p.main(arguments) == 1
        assert message in capsys.readouterr().err
        assert (tmp_path / "s.pdb").read_bytes() == before
        assert sorted(p.name for p in tmp_path.iterdir()) == ["s.pdb"]


def test_a_corrupt_gzip_file_named_only_gz_gives_an_error_not_a_traceback(tmp_path, monkeypatch, capsys):
    """A .gz name gives no format, so the file is sniffed first; a cut gzip stream must fail there cleanly too."""
    write_pdb(tmp_path / "good.pdb.gz", protein_atoms("A", helix_trace(5), PROTEIN_12[:5]))
    (tmp_path / "cut.gz").write_bytes((tmp_path / "good.pdb.gz").read_bytes()[:12])
    with pytest.raises(StructureFileException, match="Cannot read structure file"):
        read_structure(tmp_path / "cut.gz")
    monkeypatch.chdir(tmp_path)
    assert p2p.main(["cut.gz"]) == 1
    assert "ERROR: Cannot read structure file cut.gz" in capsys.readouterr().err


def test_hybrid_36_residue_numbers_are_numbers_for_mobile(tmp_path):
    """PDB files that number past 9999 in hybrid-36 (A000 is 10000) must not split a --mobile range there."""
    numbers = ["9996", "9997", "9998", "9999", "A000", "A001", "A002", "A003"]
    atoms = protein_atoms("A", helix_trace(8), PROTEIN_12[:8])
    for a in atoms:
        a["resseq"] = numbers[a["resseq"] - 1]
    write_pdb(tmp_path / "h36.pdb", atoms)
    conversion = convert(tmp_path / "h36.pdb", mobile=[p2p.parse_mobile_selection("A:9998-10001")])
    assert [index for index, free in enumerate(p2p.free_bead_flags(conversion)[0]) if free] == [2, 3, 4, 5]
    # the ends of the upper- and lower-case ranges (https://cci.lbl.gov/hybrid_36/)
    residue = conversion.segments[0].residues[0]
    for text, value in (("A000", 10000), ("ZZZZ", 1223055), ("a000", 1223056), ("zzzz", 2436111), ("Ab00", None)):
        residue.resseq = text
        assert p2p._residue_number(residue) == value


def test_a_calcium_ion_with_another_residue_name_is_known_by_its_element(tmp_path):
    """Only the element tells a calcium named CA in a residue that is not named CA from an alpha carbon.
    The ion is an ATOM record, as simulation packages write ions, so the rule that drops an unbonded
    HETATM residue cannot hide it."""
    protein = protein_atoms("A", helix_trace(4), PROTEIN_12[:4])
    ion = [atom("CA", "CAL", "I", 1, (20.0, 0, 0), element="CA")]
    write_pdb(tmp_path / "cal.pdb", protein + ion)
    conversion = convert(tmp_path / "cal.pdb")
    assert [s.sequence for s in conversion.segments] == ["MAGH"]
    assert conversion.excluded.hetero == {"CAL": 1}


def test_a_named_modified_residue_at_a_chain_end_keeps_the_distance_fallback(tmp_path):
    """A C-terminal HETATM MSE next to a residue whose C atom is not modelled is still part of the chain;
    the linking-atom rule is for residues known by their backbone alone (an ATP, a dNTP)."""
    trace = helix_trace(3)
    atoms = [a for a in protein_atoms("A", trace[:2], ["ARG", "GLU"])
             if not (a["resname"] == "GLU" and a["name"] == "C")]
    atoms += protein_atoms("A", trace[2:], ["MSE"], start=3, record="HETATM")
    write_pdb(tmp_path / "mse.pdb", atoms)
    assert [s.sequence for s in convert(tmp_path / "mse.pdb").segments] == ["REM"]


def test_output_names_that_differ_only_in_unicode_normal_form_are_one_file(tmp_path):
    """macOS treats café written with U+00E9 and with e + U+0301 as one file name."""
    composed, decomposed = tmp_path / "café.md", tmp_path / "café.md"
    problem = p2p._output_problem({"report": composed, "freeze file": decomposed}, overwrite=True)
    assert problem is not None and "cannot be the same file" in problem


# ---------------------------------------------------------------------------
# the keyfile, and --protein-mode
# ---------------------------------------------------------------------------

def _protein_structure(path, sequence, chain="A"):
    """A one-chain protein on an ideal helix, with the given one-letter sequence."""
    write_pdb(path, protein_atoms(chain, helix_trace(len(sequence)), [THREE_LETTER[c] for c in sequence]))


def _fake_metapredict(monkeypatch, disordered, folded):
    """Stand in for metapredict's child process: every sequence gets these boundaries, and each call is recorded."""
    calls = []

    def run_metapredict(sequences):
        calls.append(list(sequences))
        return "9.9-test", [p2p.MetapredictChain(folded=folded, disordered=disordered, warnings=[])
                            for _ in sequences]

    monkeypatch.setattr(p2p, "_run_metapredict", run_metapredict)
    return calls


def _skip_without_metapredict():
    """Skip unless metapredict is installed. find_spec does not import it, so PyTorch never loads here."""
    if importlib.util.find_spec("metapredict") is None:
        pytest.skip("metapredict is not installed")


def _run_pimms(directory, keyfile):
    """Run PIMMS on a keyfile in a directory, cut to a few steps with a fixed seed."""
    text = _shortened_keyfile((directory / keyfile).read_text()) + "SEED : 3\nRESTART_FREQ : 40\n"
    (directory / "short.kf").write_text(text)
    return subprocess.run([sys.executable, str(PIMMS_SCRIPT), "-k", "short.kf"], cwd=str(directory),
                          env=_clean_env(), capture_output=True, text=True, timeout=600, check=False)


def test_without_a_parameter_file_no_keyfile_is_written_and_the_report_says_why(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    assert p2p.main(["h.pdb"]) == 0
    assert "No keyfile written: give --parameter-file and --temperature, or --protein-mode" in capsys.readouterr().out
    assert not (tmp_path / "h_keyfile.kf").exists()
    assert "No keyfile was written, because no parameter file was given." in (tmp_path / "h_report.md").read_text()
    # --no-keyfile says nothing about it
    assert p2p.main(["h.pdb", "-o", "b.pimms", "--no-keyfile"]) == 0
    assert "keyfile" not in capsys.readouterr().out
    assert "No keyfile was written" not in (tmp_path / "b_report.md").read_text()


def test_a_keyfile_with_the_users_parameter_file_runs_in_pimms(tmp_path):
    """The default freeze file freezes everything, so the keyfile leaves it out and lets the whole chain move."""
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    _write_parameter_file(tmp_path / "params.prm", sorted(set("MAGHKELS")))
    result = subprocess.run([sys.executable, str(PDB2PIMMS_SCRIPT), "h.pdb", "--parameter-file", "params.prm",
                             "--temperature", "37.5"], cwd=str(tmp_path), env=_clean_env(), capture_output=True,
                            text=True, timeout=300, check=False)
    assert result.returncode == 0, result.stderr
    assert "h_keyfile.kf: run it from" in result.stdout and "with: PIMMS -k h_keyfile.kf\n" in result.stdout

    keyfile = (tmp_path / "h_keyfile.kf").read_text()
    with open(tmp_path / "h_restart.pimms", "rb") as handle:
        box = pickle.load(handle)["DIMENSIONS"]
    for line in ("RESTART_FILE         : h_restart.pimms", "PARAMETER_FILE       : params.prm",
                 "DIMENSIONS           : %s" % " ".join(str(d) for d in box), "HARDWALL             : True",
                 "TEMPERATURE          : 37.5", "MOVE_CRANKSHAFT      : 0.9", "MOVE_CHAIN_PIVOT     : 0.1",
                 "CRANKSHAFT_SUBSTEPS  : 500", "N_STEPS              : %d" % p2p.KEYFILE_N_STEPS):
        assert line in keyfile
    assert "FREEZE_FILE" not in keyfile
    assert not (tmp_path / p2p.BUNDLED_PARAMETER_FILE).exists()
    report = (tmp_path / "h_report.md").read_text()
    assert "## Running the keyfile" in report and "the whole structure moves" in report

    run = _run_pimms(tmp_path, "h_keyfile.kf")
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    assert "Simulation complete" in run.stdout


def test_a_parameter_file_that_misses_a_bead_type_writes_nothing(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    _write_parameter_file(tmp_path / "params.prm", sorted(set("MAGHKEL")))  # no S
    assert p2p.main(["h.pdb", "--parameter-file", "params.prm", "--temperature", "50"]) == 1
    error = capsys.readouterr().err
    assert "does not define every bead type the structure uses (A E G H K L M S)" in error
    assert "no interactions for S; no ANGLE_PENALTY line for S" in error
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb", "params.prm"]


def test_a_keyfile_pimms_refuses_is_never_put_in_place(tmp_path, monkeypatch, capsys):
    """The keyfile is checked as written: one with no TEMPERATURE fails PIMMS's parser, and nothing is kept."""
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    _write_parameter_file(tmp_path / "params.prm", sorted(set("MAGHKELS")))
    original = p2p.keyfile_text

    def without_temperature(*arguments, **settings):
        return "\n".join(line for line in original(*arguments, **settings).split("\n")
                         if not line.startswith("TEMPERATURE"))

    monkeypatch.setattr(p2p, "keyfile_text", without_temperature)
    assert p2p.main(["h.pdb", "--parameter-file", "params.prm", "--temperature", "50"]) == 1
    assert "PIMMS refused the keyfile pdb2pimms built" in capsys.readouterr().err
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb", "params.prm"]
    # the check ran in a temporary directory, and the working directory is back where it was
    assert pathlib.Path.cwd().resolve() == tmp_path.resolve()


@pytest.mark.parametrize("arguments", [
    ["h.pdb", "--parameter-file", "p.prm"],
    ["h.pdb", "--temperature", "50"],
    ["h.pdb", "--keyfile", "k.kf"],
    ["h.pdb", "--protein-mode", "--no-freeze"],
    ["h.pdb", "--protein-mode", "--no-keyfile", "--temperature", "50"],
    ["h.pdb", "--protein-mode", "--keyfile", "k.kf", "--no-keyfile"],
    ["h.pdb", "--parameter-file", "p.prm", "--temperature", "-1"],
    ["h.pdb", "--protein-mode", "--temperature", "nan"],
])
def test_main_keyfile_usage_errors_exit_with_status_2(tmp_path, monkeypatch, arguments):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as stop:
        p2p.main(arguments)
    assert stop.value.code == 2


def test_a_missing_parameter_file_or_one_named_like_an_output_is_refused(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    write_pdb(tmp_path / "h.pdb", protein_atoms("A", helix_trace(8), PROTEIN_12[:8]))
    assert p2p.main(["h.pdb", "--parameter-file", "absent.prm", "--temperature", "50"]) == 1
    assert "--parameter-file absent.prm: no such file" in capsys.readouterr().err
    _write_parameter_file(tmp_path / "params.prm", sorted(set("MAGHKELS")))
    assert p2p.main(["h.pdb", "--parameter-file", "params.prm", "--temperature", "50", "--keyfile",
                     "params.prm"]) == 1
    assert "the parameter file and the keyfile cannot be the same file" in capsys.readouterr().err
    assert sorted(p.name for p in tmp_path.iterdir()) == ["h.pdb", "params.prm"]


def test_protein_mode_maps_the_prediction_onto_every_chain(tmp_path, monkeypatch):
    """One sequence per source chain, joined across a gap; nucleotides and hetero groups stay frozen."""
    calls = _fake_metapredict(monkeypatch, disordered=[[0, 3], [9, 11]], folded=[[3, 9]])
    atoms = protein_atoms("A", helix_trace(5), PROTEIN_12[:5])
    atoms += protein_atoms("A", helix_trace(6, origin=(0, 30, 0)), PROTEIN_12[5:11], start=7)
    atoms += nucleic_atoms("B", helix_trace(4, origin=(30, 0, 0), radius=NUCLEIC_RADIUS, rise=NUCLEIC_RISE,
                                            twist=NUCLEIC_TWIST_DEGREES), ["G", "C", "A", "U"])
    atoms.append(atom("ZN", "ZN", "A", 101, (0, 15, 0), record="HETATM", element="ZN"))
    write_pdb(tmp_path / "mixed.pdb", atoms)
    conversion = convert(tmp_path / "mixed.pdb", keep_hetero=True, protein_mode=True)

    assert [s.sequence for s in conversion.segments] == ["MAGHK", "ELSWDV", "gcau", "b"]
    assert calls == [["MAGHKELSWDV"]]
    assert p2p.METAPREDICT_SETTINGS == {"disorder_threshold": 0.5, "minimum_IDR_size": 12,
                                        "minimum_folded_domain": 50, "gap_closure": 10, "version": "V3"}
    assert p2p.free_bead_flags(conversion) == [[True, True, True, False, False],
                                               [False, False, False, False, True, True],
                                               [False] * 4, [False]]
    chain = conversion.disorder.chains[0]
    assert (chain.chain, chain.folded, chain.disordered) == ("A", [(3, 8)], [(0, 2), (9, 10)])
    # --mobile adds to it
    mobile = convert(tmp_path / "mixed.pdb", keep_hetero=True, protein_mode=True,
                     mobile=[p2p.parse_mobile_selection("B:2")])
    assert p2p.free_bead_flags(mobile)[2] == [False, True, False, False]


def test_protein_mode_reads_the_amino_acid_not_the_bead_type(tmp_path, monkeypatch):
    calls = _fake_metapredict(monkeypatch, disordered=[], folded=[[0, 5]])
    _protein_structure(tmp_path / "p.pdb", "MAGHK")
    conversion = convert(tmp_path / "p.pdb", protein_mode=True, bead_map={"ALA": "Z"})
    assert conversion.segments[0].sequence == "MZGHK" and calls[-1] == ["MAGHK"]

    def residue(identity, code):
        return p2p.Residue(chain="A", resseq="7", icode="", resname="XYZ", record="ATOM", atoms={}, hydrogens=[],
                           kind=p2p.KIND_AMINO_ACID, code=code, identity_code=identity)

    assert p2p._metapredict_letter(residue("U", "U")) == ("C", True)
    assert p2p._metapredict_letter(residue("O", "O")) == ("K", True)
    assert p2p._metapredict_letter(residue("X", "A")) == ("A", True)
    with pytest.raises(StructureConversionException, match="metapredict reads only the 20 standard amino acids"):
        p2p._metapredict_letter(residue("X", "X"))


def test_protein_mode_without_metapredict_fails_cleanly(tmp_path, monkeypatch, capsys):
    """The child process reports the missing module, and the conversion stops before writing anything."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(p2p, "METAPREDICT_MODULE", "metapredict_is_not_installed_here")
    _protein_structure(tmp_path / "p.pdb", "MAGHKELSWD")
    assert p2p.main(["p.pdb", "--protein-mode"]) == 1
    assert "metapredict, which is not installed" in capsys.readouterr().err
    assert sorted(p.name for p in tmp_path.iterdir()) == ["p.pdb"]


def test_protein_mode_freezes_the_folded_domain_and_writes_a_keyfile(tmp_path, monkeypatch, capsys):
    _skip_without_metapredict()
    monkeypatch.chdir(tmp_path)
    sequence = DISORDERED_TAIL + UBIQUITIN + DISORDERED_TAIL
    _protein_structure(tmp_path / "ubq.pdb", sequence)
    assert p2p.main(["ubq.pdb", "--protein-mode"]) == 0
    printed = capsys.readouterr().out
    assert sorted(p.name for p in tmp_path.iterdir()) == [p2p.BUNDLED_PARAMETER_FILE, "ubq.pdb", "ubq_freeze.txt",
                                                          "ubq_keyfile.kf", "ubq_report.md", "ubq_restart.pimms"]
    assert "%s: a copy of the %s that ships with PIMMS" % (p2p.BUNDLED_PARAMETER_FILE, p2p.BUNDLED_PARAMETER_FILE) in printed

    # the freeze file freezes exactly metapredict's folded residues
    folded = p2p._run_metapredict([sequence])[1][0].folded
    expected = tuple(index for first, end in folded for index in range(first, end))
    assert 0 not in expected and 50 in expected and len(sequence) - 1 not in expected
    assert FreezeFile("ubq_freeze.txt").resolve({1: len(sequence)}) == ([1], {1: expected})

    keyfile = (tmp_path / "ubq_keyfile.kf").read_text()
    for line in ("PARAMETER_FILE       : %s" % p2p.BUNDLED_PARAMETER_FILE, "TEMPERATURE          : 50",
                 "FREEZE_FILE          : ubq_freeze.txt", "MOVE_CRANKSHAFT      : 1", "HARDWALL             : True"):
        assert line in keyfile
    assert "MOVE_CHAIN_PIVOT" not in keyfile
    assert (tmp_path / p2p.BUNDLED_PARAMETER_FILE).read_bytes() == p2p.bundled_parameter_file().read_bytes()
    report = (tmp_path / "ubq_report.md").read_text()
    assert "## Disorder prediction" in report and "| A | 136 | A:" in report

    # a second conversion next to an identical copy keeps it without --overwrite; a changed copy is refused
    assert p2p.main(["ubq.pdb", "--protein-mode", "-o", "second.pimms"]) == 0
    assert "Kept" in capsys.readouterr().out
    (tmp_path / p2p.BUNDLED_PARAMETER_FILE).write_text("edited\n")
    assert p2p.main(["ubq.pdb", "--protein-mode", "-o", "third.pimms"]) == 1
    assert "%s already exists; give --overwrite" % p2p.BUNDLED_PARAMETER_FILE in capsys.readouterr().err
    assert not (tmp_path / "third.pimms").exists()


def test_protein_mode_with_nothing_disordered_or_nothing_folded(tmp_path, monkeypatch, capsys):
    _skip_without_metapredict()
    monkeypatch.chdir(tmp_path)

    # ubiquitin alone is folded throughout: the freeze file freezes it all, so the keyfile lets it all move
    _protein_structure(tmp_path / "ubq.pdb", UBIQUITIN)
    assert p2p.main(["ubq.pdb", "--protein-mode"]) == 0
    captured = capsys.readouterr()
    assert "metapredict predicts no disordered residue" in captured.err
    assert FreezeFile("ubq_freeze.txt").chains == [1]
    keyfile = (tmp_path / "ubq_keyfile.kf").read_text()
    assert "FREEZE_FILE" not in keyfile and "MOVE_CHAIN_PIVOT     : 0.1" in keyfile

    # a short disordered peptide: nothing to freeze, so no freeze file, and metapredict's warning is kept
    _protein_structure(tmp_path / "gs.pdb", DISORDERED_TAIL[:10])
    assert p2p.main(["gs.pdb", "--protein-mode"]) == 0
    captured = capsys.readouterr()
    assert "No freeze file written: every residue is free" in captured.out
    assert "metapredict, on chain A (10 residues)" in captured.err
    assert not (tmp_path / "gs_freeze.txt").exists()
    assert "FREEZE_FILE" not in (tmp_path / "gs_keyfile.kf").read_text()


def test_pimms_runs_the_protein_mode_keyfile_and_moves_only_the_disordered_tails(tmp_path):
    """The whole point: the keyfile runs as written, the folded domain stays put and the tails move."""
    _skip_without_metapredict()
    sequence = DISORDERED_TAIL + UBIQUITIN + DISORDERED_TAIL
    _protein_structure(tmp_path / "ubq.pdb", sequence)
    result = subprocess.run([sys.executable, str(PDB2PIMMS_SCRIPT), "ubq.pdb", "--protein-mode"], cwd=str(tmp_path),
                            env=_clean_env(), capture_output=True, text=True, timeout=300, check=False)
    assert result.returncode == 0, result.stderr

    run = _run_pimms(tmp_path, "ubq_keyfile.kf")
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    assert "Simulation complete" in run.stdout
    with open(tmp_path / "ubq_restart.pimms", "rb") as handle:
        start = np.array(pickle.load(handle)["CHAINS"][1][0])
    with open(tmp_path / "restart.pimms", "rb") as handle:
        end = pickle.load(handle)
    frozen = np.array(end["FROZEN_BEADS"][1])
    free = np.setdiff1d(np.arange(len(sequence)), frozen)
    moved = np.any(np.array(end["CHAINS"][1][0]) != start, axis=1)
    assert not moved[frozen].any() and moved[free].any()


def test_protein_mode_runs_in_a_process_whose_openmp_runtime_has_started(tmp_path):
    """metapredict imports PyTorch, whose OpenMP runtime used to abort a process that had already started
    PIMMS's ("OMP: Error #15"), as a full test-suite run does. It runs in a child process now."""
    _skip_without_metapredict()
    from pimms import mega_crank_fast
    assert mega_crank_fast.openmp_info()["max_threads"] >= 1  # starts this process's OpenMP runtime
    _protein_structure(tmp_path / "ubq.pdb", DISORDERED_TAIL + UBIQUITIN + DISORDERED_TAIL)
    conversion = convert(tmp_path / "ubq.pdb", protein_mode=True)
    assert conversion.disorder.chains[0].folded and conversion.disorder.chains[0].disordered
    assert "torch" not in sys.modules
