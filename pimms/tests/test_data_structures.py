import pytest

from pimms.data_structures import AnalysisSettings, FreezeFile
from pimms.latticeExceptions import KeyFileException, UnfinishedCodeException


def _write_freeze_file(tmp_path, text, name="freeze.txt"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


class _StubLattice:
    def __init__(self, chain_ids):
        self.chains = {c: None for c in chain_ids}


def test_analysis_settings_holds_cluster_threshold():
    assert AnalysisSettings(cluster_threshold=3).cluster_threshold == 3


def test_freeze_file_parses_chains_comments_and_duplicates(tmp_path):
    path = _write_freeze_file(tmp_path, """
# freeze two groups
C 0 1 2

C 10 11   # inline comment
C 1
""".lstrip())

    freeze = FreezeFile(path)

    assert sorted(freeze.chains) == [0, 1, 2, 10, 11]      # deduplicated
    assert freeze.beads == []
    assert freeze.filename == path


def test_freeze_file_missing_file_raises():
    with pytest.raises(KeyFileException, match="Unable to find FREEZE FILE"):
        FreezeFile("this-file-does-not-exist.txt")


def test_freeze_file_bead_directive_is_not_implemented(tmp_path):
    """A well-formed 'B' line must reach the not-implemented guard.

    The bead branch used to split on a single space rather than on whitespace, so
    the empty field left by the space after the 'B' made int() fail and the line was
    reported as a malformed parse error instead.
    """
    path = _write_freeze_file(tmp_path, "B 1 2\n")

    with pytest.raises(UnfinishedCodeException):
        FreezeFile(path)


def test_freeze_file_reports_genuinely_malformed_lines(tmp_path):
    with pytest.raises(ValueError, match="Error parsing chains"):
        FreezeFile(_write_freeze_file(tmp_path, "C 1 not_an_int\n"))


def test_validate_freeze_file_accepts_chains_that_exist(tmp_path):
    freeze = FreezeFile(_write_freeze_file(tmp_path, "C 1 2\n"))

    freeze.validate_freeze_file(_StubLattice([1, 2, 3]))     # must not raise


def test_validate_freeze_file_rejects_unknown_chain_and_says_so(tmp_path):
    """The error used to claim the chain WAS present - the opposite of the check."""
    freeze = FreezeFile(_write_freeze_file(tmp_path, "C 1 99\n"))

    with pytest.raises(KeyFileException) as excinfo:
        freeze.validate_freeze_file(_StubLattice([1, 2, 3]))

    message = str(excinfo.value)
    assert "99" in message
    assert "NOT present" in message
