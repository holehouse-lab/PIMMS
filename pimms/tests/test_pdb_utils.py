import numpy as np
import pytest

from pimms import pdb_utils
from pimms.latticeExceptions import PDBException


def test_one_to_three_known_and_fallback(monkeypatch):
    monkeypatch.setattr(pdb_utils.CONFIG, "ONE_TO_THREE", {"A": "ALA"})

    assert pdb_utils.one_to_three("A") == "ALA"
    assert pdb_utils.one_to_three("Q") == "XXQ"
    assert pdb_utils.one_to_three("AB") == "XAB"
    assert pdb_utils.one_to_three("LONG") == "LON"


def test_build_section_string_justifications_and_validation():
    assert pdb_utils.build_section_string("AB", 4, "L") == "AB  "
    assert pdb_utils.build_section_string("AB", 4, "R") == "  AB"
    assert pdb_utils.build_section_string("AB", 5, "C") == " AB  "

    with pytest.raises(PDBException, match="longer than the allowed column"):
        pdb_utils.build_section_string("TOOLONG", 3)

    with pytest.raises(PDBException, match="Invalid section justification"):
        pdb_utils.build_section_string("AB", 4, "X")


def test_build_line_has_fixed_length_and_overflow_guard():
    line = pdb_utils.build_line(["ABCDEF"], [[1, 6]])
    assert len(line) == 80
    assert line[:6] == "ABCDEF"

    with pytest.raises(PDBException, match="longer than 80"):
        pdb_utils.build_line(["A" * 81], [[1, 81]])

    with pytest.raises(PDBException, match="must have equal length"):
        pdb_utils.build_line(["A", "B"], [[1, 1]])


def test_build_atom_ter_model_conect_lines_shape():
    atom = pdb_utils.build_atom_line(1, "CA", "ALA", "A", 1, 1.0, 2.0, 3.0, 1)
    ter = pdb_utils.build_ter_line(2, "ALA", "A", 1)
    model = pdb_utils.build_model_line(1)
    conect = pdb_utils.build_conect_line(1, 2)

    assert atom.startswith("ATOM")
    assert ter.startswith("TER")
    assert model.startswith("MODEL")
    assert conect.startswith("CONECT")

    assert len(atom.rstrip("\n")) == 80
    assert len(ter.rstrip("\n")) == 80
    assert len(model) == 80
    assert len(conect.rstrip("\n")) == 80


def test_build_cryst_line_2d_3d_and_invalid_dimensions():
    c2 = pdb_utils.build_cryst_line([3, 4], spacing=4.0)
    c3 = pdb_utils.build_cryst_line([3, 4, 5], spacing=4.0)

    assert c2.startswith("CRYST1")
    assert c3.startswith("CRYST1")
    assert len(c2.rstrip("\n")) == 80
    assert len(c3.rstrip("\n")) == 80

    with pytest.raises(PDBException, match="only supports 2D/3D"):
        pdb_utils.build_cryst_line([3], spacing=4.0)


def test_initialize_build_finalize_pdb_file(tmp_path):
    fn = tmp_path / "x.pdb"

    pdb_utils.initialize_pdb_file([4, 4], spacing=4.0, filename=str(fn))
    pdb_utils.build_pdb_file(
        latticeObject=None,
        spacing=4.0,
        filename=str(fn),
        usePositionsOnly={
            "dimensions": [4, 4],
            "length": 2,
            "positions": [[0, 0], [1, 0]],
            "sequence": "AG",
        },
        write_connect=False,
    )
    pdb_utils.finalize_pdb_file(str(fn))

    text = fn.read_text()
    assert text.startswith("CRYST1")
    assert "MODEL" in text
    assert "ATOM" in text
    assert text.endswith("END\n")


def test_write_positions_to_file_valid_and_sequence_checks(tmp_path):
    fn = tmp_path / "positions.pdb"

    pdb_utils.write_positions_to_file([[0, 0], [1, 0]], str(fn), spacing=4.0, sequence="AG")
    text = fn.read_text()
    assert "CRYST1" in text
    assert "ATOM" in text

    with pytest.raises(PDBException, match="sequence which is not a string"):
        pdb_utils.write_positions_to_file([[0, 0]], str(tmp_path / "bad1.pdb"), spacing=4.0, sequence=["A"])

    with pytest.raises(PDBException, match="not the same length"):
        pdb_utils.write_positions_to_file([[0, 0], [1, 0]], str(tmp_path / "bad2.pdb"), spacing=4.0, sequence="A")


def test_write_positions_to_file_rejects_empty_and_dimensionality_mismatch(tmp_path):
    with pytest.raises(PDBException, match="positions is empty"):
        pdb_utils.write_positions_to_file([], str(tmp_path / "empty.pdb"), spacing=4.0)

    with pytest.raises(PDBException, match="does not match"):
        pdb_utils.write_positions_to_file([[0, 0]], str(tmp_path / "dim_mismatch.pdb"), spacing=4.0, dimensions=[5, 5, 5])


def test_write_positions_to_file_rejects_out_of_bounds_equal_dimension(tmp_path):
    # Coordinate equal to dimension length is out of bounds for 0-indexed lattice coordinates.
    with pytest.raises(PDBException, match="lies outside"):
        pdb_utils.write_positions_to_file([[2, 0]], str(tmp_path / "oob.pdb"), spacing=4.0, dimensions=[2, 2])


def test_write_positions_to_file_infers_nonzero_box_for_small_positions(tmp_path):
    fn = tmp_path / "small_box.pdb"
    pdb_utils.write_positions_to_file([[0, 0], [1, 0]], str(fn), spacing=4.0)

    first_line = fn.read_text().splitlines()[0]
    # Columns 7-15 encode box dimension a. The inferred lattice is 2 sites wide, and the
    # periodic cell of an L-site axis is L * spacing (see the CRYST1 tests below).
    a_val = float(first_line[6:15])
    assert a_val == pytest.approx(8.0)


# ---------------------------------------------------------------------------
# CRYST1 records the PERIODIC UNIT CELL
#
# Regression tests for a real bug. The CRYST1 line used to be written as
# (L - 1) * spacing, i.e. the extent spanned by the occupied sites rather than the
# period of the lattice. Sites L-1 and 0 are periodic neighbours one lattice unit
# apart, so the period is L * spacing. The consequences of getting this wrong were
# not cosmetic: the PDB disagreed with the box PIMMS writes into its own XTC, every
# PBC-aware mdtraj/VMD calculation on START.pdb was off by a lattice unit, and
# lemonade.load(pdb=...) (which has to infer the box from the file when no keyfile is
# given) inferred an L-1 box and then wrapped the coordinates into it.
# ---------------------------------------------------------------------------

def _cryst_abc(line):
    """Pull (a, b, c) out of the fixed CRYST1 columns."""
    return (float(line[6:15]), float(line[15:24]), float(line[24:33]))


def test_cryst_line_encodes_full_periodic_cell_3d():
    a, b, c = _cryst_abc(pdb_utils.build_cryst_line([10, 20, 30], spacing=3.65))

    assert a == pytest.approx(10 * 3.65)
    assert b == pytest.approx(20 * 3.65)
    assert c == pytest.approx(30 * 3.65)


def test_cryst_line_2d_uses_one_lattice_unit_in_z():
    a, b, c = _cryst_abc(pdb_utils.build_cryst_line([7, 9], spacing=4.0))

    assert a == pytest.approx(7 * 4.0)
    assert b == pytest.approx(9 * 4.0)
    # z is not periodic in a 2D system; one lattice unit, matching the XTC box
    assert c == pytest.approx(4.0)


def test_cryst_line_matches_the_xtc_box_convention():
    """The PDB topology and the XTC frames must describe the same box."""
    from pimms.lattice_utils import _lattice_frame_xyz_and_box

    class _StubChain:
        def get_output_positions(self, autocenter=False, unwrap=False):
            return [[0, 0, 0], [0, 0, 1]]

    class _StubLattice:
        dimensions = [10, 20, 30]
        chains = {1: _StubChain()}

    spacing = 3.65
    _xyz, box = _lattice_frame_xyz_and_box(_StubLattice(), spacing)
    a, b, c = _cryst_abc(pdb_utils.build_cryst_line([10, 20, 30], spacing))

    # the XTC box is in nm, the CRYST1 line in angstroms
    assert a == pytest.approx(box[0][0][0] * 10.0)
    assert b == pytest.approx(box[0][1][1] * 10.0)
    assert c == pytest.approx(box[0][2][2] * 10.0)


def test_build_pdb_file_rejects_bad_use_positions_only_dict(tmp_path):
    fn = tmp_path / "bad_usepos.pdb"
    pdb_utils.initialize_pdb_file([4, 4], spacing=4.0, filename=str(fn))

    with pytest.raises(PDBException, match="INVALID usePositionsOnly dictionary"):
        pdb_utils.build_pdb_file(None, spacing=4.0, filename=str(fn), usePositionsOnly={"bad": "dict"})

    with pytest.raises(PDBException, match="missing one of the keywords"):
        pdb_utils.build_pdb_file(
            None,
            spacing=4.0,
            filename=str(fn),
            usePositionsOnly={"dimensions": [4, 4], "length": 1, "positions": [[0, 0]], "oops": "x"},
        )


# ---------------------------------------------------------------------------
# CONECT records once the 5-column atom serials wrap (>= 100,000 serials)
#
# Regression test for a real bug. Serials run over every ATOM and TER record and
# are written modulo 100000, so in a file with N > 99999 serials the written serials
# 1 .. N - 100000 each belong to TWO records. The writer only dropped CONECT records
# whose unwrapped serial was >= 100000; a CONECT between two LOW serials was still
# written although each of them now named two atoms, and mdtraj resolves a
# duplicated serial to the last atom carrying it. A 1010 x 100-mer run lost ~2000
# bonds, gained 39 cross-chain ones, and find_molecules() stopped returning the
# chains. The fix writes a CONECT only when both of its serials are unique in the
# file, so a bond can go missing but can never join the wrong atoms.
# ---------------------------------------------------------------------------

def _snake_positions(n_beads: int, width: int, z0: int) -> list[list[int]]:
    """
    Lay ``n_beads`` out as a boustrophedon walk filling ``width x width`` layers.

    Consecutive beads are always lattice neighbours (the walk reverses direction
    at the end of every row and every layer), so every true bond has length 1 and
    any bond longer than sqrt(3) lattice units is unambiguously wrong.

    Parameters
    ----------
    n_beads : int
        Number of beads in the chain.

    width : int
        Edge length of each square layer, in lattice sites.

    z0 : int
        Layer index the walk starts in.

    Returns
    -------
    list of list of int
        ``n_beads`` ``[x, y, z]`` lattice positions in walk order.
    """
    idx = np.arange(n_beads)
    row = idx // width                     # global row index, continuous across layers
    col = idx % width
    layer = row // width
    row_in_layer = row % width
    x = np.where(row % 2 == 0, col, width - 1 - col)
    y = np.where(layer % 2 == 0, row_in_layer, width - 1 - row_in_layer)
    z = z0 + layer
    return np.stack([x, y, z], axis=1).tolist()


def test_conect_never_joins_wrong_atoms_after_serial_wrap(tmp_path):
    md = pytest.importorskip("mdtraj")

    width = 50
    spacing = 3.65
    # three long chains take the serials to just below the wrap, then short chains
    # straddle it: the wrapped serials 1..N-100000 are then carried both by beads in
    # the middle of chain 1 and by beads AND TER records of the short chains (a TER
    # twin is what used to produce the cross-chain bonds)
    lengths = [40000, 30000, 29900] + [37] * 20

    class _Chain:
        def __init__(self, n_beads, z0):
            self.chainType = 0
            self.sequence = "A" * n_beads
            self._positions = _snake_positions(n_beads, width, z0)

        def get_output_positions(self, autocenter=False, unwrap=False):
            return self._positions

    chains = {}
    z0 = 0
    for chain_id, n_beads in enumerate(lengths, start=1):
        chains[chain_id] = _Chain(n_beads, z0)
        z0 += -(-n_beads // (width * width)) + 1     # own block of layers, one-layer gap

    class _Lattice:
        dimensions = [width, width, z0]

    _Lattice.chains = chains

    fn = str(tmp_path / "big.pdb")
    pdb_utils.initialize_pdb_file(_Lattice.dimensions, spacing, fn)
    pdb_utils.build_pdb_file(_Lattice(), spacing, filename=fn, write_connect=True)
    pdb_utils.finalize_pdb_file(fn)

    # the serial layout from the PDB spec: one serial per ATOM, then one for the TER
    # closing each chain, counted from 1
    n_serials = sum(lengths) + len(lengths)
    assert n_serials > 100000

    traj = md.load(fn)
    top = traj.topology
    assert top.n_atoms == sum(lengths)
    assert top.n_chains == len(lengths)

    chain_of = np.array([atom.residue.chain.index for atom in top.atoms])
    lattice_xyz = traj.xyz[0] * 10.0 / spacing
    bonds = np.array(sorted((min(b.atom1.index, b.atom2.index), max(b.atom1.index, b.atom2.index))
                            for b in top.bonds))
    assert len(bonds) > 0
    a, b = bonds[:, 0], bonds[:, 1]

    # safety: every bond mdtraj built from the CONECT records joins two successive
    # beads of one chain, one lattice unit apart
    assert np.all(chain_of[a] == chain_of[b])
    assert np.all(b - a == 1)
    bond_length = np.linalg.norm(lattice_xyz[a] - lattice_xyz[b], axis=1)
    assert np.all(bond_length <= np.sqrt(3) + 1e-6)

    # completeness: every bond whose two serials are unique in the file (below the
    # wrap, and with no wrapped twin s + 100000 among the serials used) is kept, and
    # only those are dropped
    expected = set()
    atom_offset = 0
    for chain_idx, n_beads in enumerate(lengths):
        for k in range(n_beads - 1):
            s1 = atom_offset + k + 1 + chain_idx      # serial = atom number + TERs before it
            s2 = s1 + 1
            if all(s < 100000 and s + 100000 > n_serials for s in (s1, s2)):
                expected.add((atom_offset + k, atom_offset + k + 1))
        atom_offset += n_beads
    assert set(map(tuple, bonds.tolist())) == expected

    # and the wrap really did cost bonds at both ends, so the test is not vacuous
    all_bonds = sum(n - 1 for n in lengths)
    assert all_bonds - 2000 < len(expected) < all_bonds
