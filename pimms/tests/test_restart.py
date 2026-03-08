import pickle

import pytest

from pimms.latticeExceptions import RestartException
from pimms.restart import RestartObject


class _DummyChain:
    def __init__(self, positions, sequence, chain_type):
        self.positions = positions
        self.sequence = sequence
        self.chainType = chain_type


class _DummyLattice:
    def __init__(self):
        self.dimensions = [6, 6]
        self.chains = {
            7: _DummyChain([[1, 1], [1, 2]], "AB", 3),
            11: _DummyChain([[3, 3]], "C", 5),
        }


def test_restart_object_initialization_defaults():
    r = RestartObject()

    assert r.energy == 0
    assert r.dimensions == []
    assert r.hardwall is False
    assert r.chains == {}
    assert r.seq2chainType == {}
    assert r.extra_chains == {}


def test_add_extra_chains_handles_empty_base_state():
    r = RestartObject()

    r.add_extra_chains([2, "AB"])

    assert sorted(r.extra_chains.keys()) == [1, 2]
    assert r.extra_chains[1][0] is None
    assert r.extra_chains[1][1] == "AB"
    assert r.extra_chains[1][2] == 0
    assert r.seq2chainType == {"AB": [0]}


def test_add_extra_chains_rejects_non_positive_count():
    r = RestartObject()

    with pytest.raises(RestartException, match="positive integer"):
        r.add_extra_chains([0, "AB"])


def test_add_extra_chains_rejects_malformed_payload():
    r = RestartObject()

    with pytest.raises(RestartException, match="EXTRA_CHAINS"):
        r.add_extra_chains(["not-an-int"])


def test_build_from_lattice_copies_positions_and_resets_extra_chains():
    r = RestartObject()
    r.extra_chains = {99: [None, "ZZ", 42]}
    lattice = _DummyLattice()

    r.build_from_lattice(lattice, hardwall=True)

    # Ensure stale extra chains are cleared on rebuild.
    assert r.extra_chains == {}

    # Ensure positions are copied, not aliased.
    lattice.chains[7].positions[0][0] = 99
    assert r.chains[7][0][0][0] == 1

    assert r.hardwall is True
    assert r.dimensions == [6, 6]


def test_update_lattice_dimensions_applies_center_offset_when_valid():
    r = RestartObject()
    r.dimensions = [3, 3]
    r.chains = {1: [[[0, 0], [1, 0]], "AB", 0]}

    r.update_lattice_dimensions([5, 5])

    assert r.dimensions == [5, 5]
    assert r.chains[1][0] == [[1, 1], [2, 1]]


def test_update_lattice_dimensions_manual_offset_failure_is_atomic():
    r = RestartObject()
    r.dimensions = [3, 3]
    r.chains = {1: [[[2, 1]], "A", 0]}

    with pytest.raises(RestartException):
        r.update_lattice_dimensions([3, 3], manual_offset=[1, 0])

    # Ensure both dimensions and positions are unchanged after failure.
    assert r.dimensions == [3, 3]
    assert r.chains[1][0] == [[2, 1]]


def test_build_from_file_invalid_pickle_raises_restart_exception(tmp_path):
    p = tmp_path / "bad_restart.pimms"
    p.write_bytes(b"not-a-pickle")

    r = RestartObject()
    with pytest.raises(RestartException, match="Error reading restart file"):
        r.build_from_file(str(p))


def test_build_from_file_rejects_non_dict_chains(tmp_path):
    p = tmp_path / "bad_schema.pimms"
    payload = {
        "DIMENSIONS": [4, 4],
        "ENERGY": 0.0,
        "HARDWALL": False,
        "CHAINS": [],
    }
    p.write_bytes(pickle.dumps(payload))

    r = RestartObject()
    with pytest.raises(RestartException, match="CHAINS entry must be a dictionary"):
        r.build_from_file(str(p))


def test_build_from_file_rejects_malformed_chain_entry(tmp_path):
    p = tmp_path / "bad_chain_entry.pimms"
    payload = {
        "DIMENSIONS": [4, 4],
        "ENERGY": 0.0,
        "HARDWALL": False,
        "CHAINS": {
            1: [[[1, 1]], "A"],
        },
    }
    p.write_bytes(pickle.dumps(payload))

    r = RestartObject()
    with pytest.raises(RestartException, match="malformed chain entry"):
        r.build_from_file(str(p))


def test_build_from_file_resets_extra_chains_and_updates_seq_map(tmp_path):
    p = tmp_path / "ok_restart.pimms"
    payload = {
        "DIMENSIONS": [5, 5],
        "ENERGY": -3.5,
        "HARDWALL": True,
        "CHAINS": {
            3: [[[1, 2], [1, 3]], "AB", 7],
        },
    }
    p.write_bytes(pickle.dumps(payload))

    r = RestartObject()
    r.extra_chains = {10: [None, "X", 0]}

    r.build_from_file(str(p))

    assert r.dimensions == [5, 5]
    assert r.energy == -3.5
    assert r.hardwall is True
    assert r.extra_chains == {}
    assert r.chains[3][1] == "AB"
    assert r.seq2chainType == {"AB": [7]}


def test_write_to_file_roundtrip(tmp_path, monkeypatch):
    out = tmp_path / "restart_roundtrip.pimms"
    monkeypatch.setattr("pimms.restart.CONFIG.RESTART_FILENAME", str(out))

    r = RestartObject()
    r.dimensions = [8, 8]
    r.energy = -12.0
    r.hardwall = True
    r.chains = {
        2: [[[1, 1], [1, 2]], "AB", 3],
        9: [[[4, 4]], "C", 5],
    }

    r.write_to_file()

    r2 = RestartObject()
    r2.build_from_file(str(out))

    assert r2.dimensions == [8, 8]
    assert r2.energy == -12.0
    assert r2.hardwall is True
    assert r2.chains == r.chains
