import pytest

from pimms import parameterfile_parser
from pimms.latticeExceptions import ParameterFileException


def _write_file(path, text):
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_parse_energy_basic_short_and_long_range(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )

    paramfile = _write_file(
        tmp_path / "params.prm",
        """
A A 1
A B 2 3 4
B B 5 6
A 0 -1
B 0 -2
0 0 0
""".strip()
        + "\n",
    )

    energy, residues, lr_energy, lr_residues, slr_energy = parameterfile_parser.parse_energy(paramfile)

    assert energy["A"]["B"] == 2
    assert energy["B"]["A"] == 2
    assert energy["0"]["0"] == 0
    assert residues[0] == "0"

    assert lr_energy["A"]["B"] == 3
    assert lr_energy["B"]["A"] == 3
    assert slr_energy["A"]["B"] == 4
    assert slr_energy["B"]["A"] == 4
    assert "A" in lr_residues and "B" in lr_residues

    used_file = tmp_path / "parameters_used_test.prm"
    assert used_file.exists()
    saved_text = used_file.read_text(encoding="utf-8")
    assert "This is a copy of the parameter file used for the simulation" in saved_text


def test_parse_energy_skips_blank_comment_and_angle_lines(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )

    paramfile = _write_file(
        tmp_path / "params.prm",
        """
# pure comment

ANGLE_PENALTY A 1 2 3 # inline comment
ANGLE_PENALTY_T_NORM B 0.1 0.2 0.3
A A 1
A 0 0
0 0 0
""".strip()
        + "\n",
    )

    energy, residues, lr_energy, lr_residues, slr_energy = parameterfile_parser.parse_energy(paramfile)

    assert energy["A"]["A"] == 1
    assert "0" in residues
    assert lr_energy == {}
    assert lr_residues == []
    assert slr_energy == {}


def test_parse_energy_raises_on_float_short_range_value(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
A A 1.5
A 0 0
0 0 0
""".strip()
        + "\n",
    )

    with pytest.raises(ParameterFileException, match="Unable to use floats"):
        parameterfile_parser.parse_energy(paramfile)


def test_parse_energy_raises_on_non_numeric_long_range_value(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
A A 1
A B 2 X
A 0 0
B 0 0
0 0 0
""".strip()
        + "\n",
    )

    with pytest.raises(ParameterFileException, match="Unable to parse line"):
        parameterfile_parser.parse_energy(paramfile)


def test_parse_energy_raises_on_missing_solvent_interactions(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
A A 1
A B 2
B B 3
""".strip()
        + "\n",
    )

    with pytest.raises(ParameterFileException, match="solvation interactions"):
        parameterfile_parser.parse_energy(paramfile)


def test_parse_energy_raises_on_nonzero_solvent_solvent(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
A A 1
A 0 0
0 0 2
""".strip()
        + "\n",
    )

    with pytest.raises(ParameterFileException, match="SOLVENT-SOLVENT"):
        parameterfile_parser.parse_energy(paramfile)


def test_parse_energy_raises_on_long_range_solvent_entries(tmp_path, monkeypatch):
    monkeypatch.setattr(
        parameterfile_parser.CONFIG,
        "OUTPUT_USED_PARAMETER_FILE",
        str(tmp_path / "parameters_used_test.prm"),
    )
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
A A 1 2
A 0 0 1
0 0 0
""".strip()
        + "\n",
    )

    with pytest.raises(ParameterFileException, match="long range solvent-solute"):
        parameterfile_parser.parse_energy(paramfile)


def test_parse_angles_handles_comments_blank_lines_and_absolute_values(tmp_path):
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
# comment

ANGLE_PENALTY A 1 2 3  # inline comment
ANGLE_PENALTY B 0 0 0
""".strip()
        + "\n",
    )

    angle_dict = parameterfile_parser.parse_angles(paramfile)
    assert angle_dict["A"] == [1, 2, 3]
    assert angle_dict["B"] == [0, 0, 0]


def test_parse_angles_t_norm_requires_temperature(tmp_path):
    paramfile = _write_file(
        tmp_path / "params.prm",
        "ANGLE_PENALTY_T_NORM A 0.5 1.0 1.5\n",
    )

    with pytest.raises(ParameterFileException, match="no temperature provided"):
        parameterfile_parser.parse_angles(paramfile)


def test_parse_angles_t_norm_scales_by_temperature(tmp_path):
    paramfile = _write_file(
        tmp_path / "params.prm",
        "ANGLE_PENALTY_T_NORM A 0.5 1.0 1.5\n",
    )

    angle_dict = parameterfile_parser.parse_angles(paramfile, temperature=100)
    assert angle_dict["A"] == [50.0, 100.0, 150.0]


def test_parse_angles_rejects_duplicate_residue_definitions(tmp_path):
    paramfile = _write_file(
        tmp_path / "params.prm",
        """
ANGLE_PENALTY A 1 2 3
ANGLE_PENALTY A 1 2 3
""".strip()
        + "\n",
    )

    with pytest.raises(ParameterFileException, match="Multiple ANGLE_PENALTY definitions"):
        parameterfile_parser.parse_angles(paramfile)
