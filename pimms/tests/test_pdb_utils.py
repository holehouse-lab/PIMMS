from pathlib import Path

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
