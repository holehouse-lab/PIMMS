"""
Regression tests for the output-writing findings of the deep audit (fixer F2).

Each test pins one finding and fails on the code as it stood at commit 59163ce:

* O1-1 / R1-1 - the streamed ``traj.xtc`` was never flushed, so completed frames
  sat in a 4 KiB buffer inside the process.
* O1-2 / A2-8 / K1-7 - a coordinate of 10000 A or more crashed the PDB writer
  with a message naming no keyword, and left a truncated ``START.pdb``.
* O1-8 / K1-3 / S1-3 - ``SAVE_AT_END`` took its unit cell from the PDB as mdtraj
  read it back (rounded, absent for dense sub-angstrom boxes, and with spurious
  off-diagonal terms in large boxes).
* O1-6 - ``SAVE_AT_END`` joined the whole buffer into a second copy to write it.
* O1-5 - ``remove_files`` swallowed deletion failures.
* O1-9 - elapsed time came from the wall clock.
* O1-4 - outputs opened with ``'w'`` wrote through a symbolic link.
* E1-3 - ``ENERGY.dat`` pushed the integer energy through a double.
* A2-6 - ``AUTOCENTER`` with ``HARDWALL`` wrote beads outside the box.
* K1-8 - ``AUTOCENTER`` with more than one chain was switched off silently.

Everything is written under ``tmp_path``.
"""

from __future__ import annotations

import os
import subprocess
import errno
import shutil
import sys
import tracemalloc
from datetime import datetime, timedelta
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

import pimms
from pimms import (
    CONFIG,
    IO_utils,
    analysis_general,
    analysis_IO,
    lattice_utils,
    pdb_utils,
)
from pimms.latticeExceptions import PDBException

SPACING = 3.65


@pytest.fixture(autouse=True)
def _work_in_tmp_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Run every test from inside its own ``tmp_path``.

    Several of the routines under test write to ``log.txt`` in the working
    directory when one exists there, so a test must never run in a directory it
    does not own.
    """
    monkeypatch.chdir(tmp_path)


# ---------------------------------------------------------------------------
# stand-ins: just enough of a Lattice / Chain for the writers
# ---------------------------------------------------------------------------
class _Chain:
    """A chain that reports fixed positions, centring them the way Chain does."""

    def __init__(
        self, positions: list[list[int]], dimensions: list[int], chain_type: int = 0
    ) -> None:
        self.positions = positions
        self.dimensions = dimensions
        self.chainType = chain_type
        self.sequence = "A" * len(positions)

    def get_output_positions(
        self, autocenter: bool = False, unwrap: bool = False
    ) -> list[list[int]]:
        if autocenter:
            return lattice_utils.center_positions(self.positions, self.dimensions)
        return self.positions


class _Lattice:
    """A lattice holding stand-in chains."""

    def __init__(
        self,
        dimensions: list[int],
        chains: list[list[list[int]]],
        hardwall: bool = False,
    ) -> None:
        self.dimensions = dimensions
        self.hardwall = hardwall
        self.chains = {
            i + 1: _Chain(p, dimensions, chain_type=i) for i, p in enumerate(chains)
        }


def _walk(n_beads: int, dimensions: list[int], seed: int) -> list[list[int]]:
    """n_beads distinct in-box sites (not a bonded walk - the writers do not care)."""
    rng = np.random.default_rng(seed)
    sites = rng.permutation(int(np.prod(dimensions)))[:n_beads]
    return [[int(v) for v in np.unravel_index(s, dimensions)] for s in sites]


def _read_xtc(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(xyz, time, step, box) of every complete frame in an XTC file."""
    with md.formats.XTCTrajectoryFile(str(path)) as fh:
        return fh.read()


def _expected_box(dimensions: list[int], spacing: float) -> np.ndarray:
    """The unit cell by definition: DIMENSIONS x LATTICE_TO_ANGSTROMS, in nm, orthorhombic.

    A 2D system has a single lattice unit along z.
    """
    edges = [d * spacing * 0.1 for d in dimensions] + (
        [spacing * 0.1] if len(dimensions) == 2 else []
    )
    return np.diag(np.array(edges, dtype=np.float32))


# ---------------------------------------------------------------------------
# O1-1 / R1-1: every streamed frame is on disk as soon as it has been written
# ---------------------------------------------------------------------------
def test_streamed_frames_are_readable_while_the_writer_is_still_open(
    tmp_path: Path,
) -> None:
    dims = [12, 12, 12]
    lattice = _Lattice(dims, [_walk(6, dims, 0)])
    pdb, xtc = tmp_path / "START.pdb", tmp_path / "traj.xtc"

    writer = lattice_utils.open_xtc_writer(
        lattice, SPACING, pdb_filename=str(pdb), xtc_filename=str(xtc)
    )
    written = [np.array(lattice.chains[1].positions)]
    try:
        for frame in range(1, 6):
            lattice.chains[1].positions = _walk(6, dims, frame)
            written.append(np.array(lattice.chains[1].positions))
            lattice_utils.write_xtc_frame(writer, lattice, SPACING)

            # read through a SECOND handle while the writer is still open: this is
            # what a monitoring script, or the disk after a SIGKILL, sees. Six beads
            # are ~128 bytes a frame, far below the 4096-byte stdio buffer the
            # frames used to sit in.
            xyz, time, step, _box = _read_xtc(xtc)
            assert len(xyz) == frame + 1
            assert list(step) == list(range(frame + 1))
            assert np.allclose(xyz[-1], written[-1] * SPACING * 0.1, atol=1e-3)
    finally:
        lattice_utils.close_xtc_writer(writer)

    xyz, _time, _step, _box = _read_xtc(xtc)
    assert np.allclose(xyz, np.array(written) * SPACING * 0.1, atol=1e-3)


def test_flush_is_public_idempotent_and_safe_on_anything_the_run_may_hold(
    tmp_path: Path,
) -> None:
    dims = [8, 8, 8]
    lattice = _Lattice(dims, [_walk(4, dims, 1)])
    xtc = tmp_path / "traj.xtc"
    writer = lattice_utils.open_xtc_writer(
        lattice,
        SPACING,
        pdb_filename=str(tmp_path / "START.pdb"),
        xtc_filename=str(xtc),
    )
    size = xtc.stat().st_size
    assert size > 0  # frame 0 is already on disk

    writer.flush()
    writer.flush()  # nothing new: nothing changes
    assert xtc.stat().st_size == size

    lattice_utils.close_xtc_writer(writer)
    writer.flush()  # closed writer: a no-op, not an error
    lattice_utils.flush_xtc_writer(writer)
    lattice_utils.flush_xtc_writer(None)  # nothing open at all

    # the SAVE_AT_END buffer accepts the same call and writes nothing
    buffer = lattice_utils.update_master_traj(
        lattice, SPACING, None, str(tmp_path / "START.pdb")
    )
    lattice_utils.flush_xtc_writer(buffer)
    assert xtc.stat().st_size == size
    assert len(_read_xtc(xtc)[0]) == 1


# ---------------------------------------------------------------------------
# O1-2 / A2-8 / K1-7: PDB coordinate overflow
# ---------------------------------------------------------------------------
def test_pdb_coordinate_limits_are_the_eight_column_field() -> None:
    assert pdb_utils.format_pdb_coordinate(9999.999) == "9999.999"
    assert pdb_utils.format_pdb_coordinate(-999.999) == "-999.999"
    assert pdb_utils.format_pdb_coordinate(3.65) == "   3.650"
    for too_wide in (10000.0, 9999.9996):
        with pytest.raises(PDBException) as err:
            pdb_utils.format_pdb_coordinate(too_wide)
        assert "DIMENSIONS" in str(err.value) and "LATTICE_TO_ANGSTROMS" in str(
            err.value
        )
    # a negative coordinate is never a site of the box, whatever the box
    with pytest.raises(PDBException) as err:
        pdb_utils.format_pdb_coordinate(-1000.0)
    assert "TRAJECTORY_PBC_UNWRAP" in str(err.value)


def test_overflowing_topology_leaves_no_truncated_pdb_and_names_the_keywords(
    tmp_path: Path,
) -> None:
    # site 2800 at 3.65 A is 10220 A: one column too wide for the PDB
    dims = [20, 20, 3000]
    lattice = _Lattice(dims, [[[0, 0, 0], [0, 0, 1]], [[5, 5, 2799], [5, 5, 2800]]])
    pdb, xtc = tmp_path / "START.pdb", tmp_path / "traj.xtc"

    for opener in (lattice_utils.open_xtc_writer, lattice_utils.start_xtc_file):
        with pytest.raises(PDBException) as err:
            opener(lattice, SPACING, pdb_filename=str(pdb), xtc_filename=str(xtc))
        assert "DIMENSIONS" in str(err.value) and "LATTICE_TO_ANGSTROMS" in str(
            err.value
        )
        assert (
            sorted(p.name for p in tmp_path.iterdir()) == []
        )  # nothing, not even a temp file

    # an earlier run's START.pdb is still whole after the failed start
    previous = "CRYST1 of the previous run\nEND\n"
    pdb.write_text(previous)
    with pytest.raises(PDBException):
        lattice_utils.open_xtc_writer(
            lattice, SPACING, pdb_filename=str(pdb), xtc_filename=str(xtc)
        )
    assert pdb.read_text() == previous
    assert sorted(p.name for p in tmp_path.iterdir()) == ["START.pdb"]


def test_write_positions_to_file_is_all_or_nothing(tmp_path: Path) -> None:
    out = tmp_path / "positions.pdb"
    with pytest.raises(PDBException):
        pdb_utils.write_positions_to_file([[0, 0, 0], [0, 0, 2800]], str(out), SPACING)
    assert list(tmp_path.iterdir()) == []

    pdb_utils.write_positions_to_file([[0, 0, 0], [0, 0, 1]], str(out), SPACING)
    lines = out.read_text().splitlines()
    assert lines[0].startswith("CRYST1") and lines[-1] == "END"
    assert sum(line.startswith("ATOM") for line in lines) == 2


# ---------------------------------------------------------------------------
# O1-8 / K1-3 / S1-3 / O1-6: the SAVE_AT_END trajectory
# ---------------------------------------------------------------------------
def _write_both_ways(
    tmp_path: Path, dims: list[int], spacing: float, n_beads: int, n_frames: int
) -> tuple[Path, Path]:
    """Write the same frames through the stream writer and through SAVE_AT_END."""
    frames = [_walk(n_beads, dims, seed) for seed in range(n_frames + 1)]

    lattice = _Lattice(dims, [frames[0]])
    stream = tmp_path / "stream.xtc"
    writer = lattice_utils.open_xtc_writer(
        lattice,
        spacing,
        pdb_filename=str(tmp_path / "stream.pdb"),
        xtc_filename=str(stream),
    )
    for positions in frames[1:]:
        lattice.chains[1].positions = positions
        lattice_utils.write_xtc_frame(writer, lattice, spacing)
    lattice_utils.close_xtc_writer(writer)

    lattice = _Lattice(dims, [frames[0]])
    buffered = tmp_path / "buffered.xtc"
    pdb = str(tmp_path / "buffered.pdb")
    lattice_utils.start_xtc_file(
        lattice, spacing, pdb_filename=pdb, xtc_filename=str(buffered)
    )
    master = None
    for positions in frames[1:]:
        lattice.chains[1].positions = positions
        master = lattice_utils.update_master_traj(lattice, spacing, master, pdb)
    if master is None:
        master = lattice_utils.start_master_traj(pdb)
    assert len(master) == n_frames + 1
    lattice_utils.save_out_sim(master, str(buffered))
    return stream, buffered


@pytest.mark.parametrize(
    "dims, spacing, n_beads, n_frames",
    [
        ([12, 12, 12], 3.65, 20, 4),  # the default spacing
        ([9, 9, 9], 3.6555, 20, 4),  # CRYST1 rounds 32.8995 A to 32.900 A
        (
            [8, 8, 8],
            0.5,
            65,
            4,
        ),  # > 1000 beads/nm^3: mdtraj drops the CRYST1 cell (K1-3)
        (
            [100, 100, 100],
            3.65,
            20,
            4,
        ),  # > 23 nm: mdtraj's cell round trip adds off-diagonals
        ([14, 14], 3.65, 20, 4),  # 2D
        ([12, 12, 12], 3.65, 20, 0),  # no frame ever buffered: frame 0 only
    ],
)
def test_save_at_end_writes_what_the_stream_writer_writes(
    tmp_path: Path, dims: list[int], spacing: float, n_beads: int, n_frames: int
) -> None:
    stream, buffered = _write_both_ways(tmp_path, dims, spacing, n_beads, n_frames)

    xyz, time, step, box = _read_xtc(buffered)
    assert len(xyz) == n_frames + 1
    assert list(step) == list(range(n_frames + 1))
    assert list(time) == [float(k) for k in range(n_frames + 1)]

    # the unit cell is DIMENSIONS x LATTICE_TO_ANGSTROMS, exactly, on every frame,
    # and exactly orthorhombic
    expected = _expected_box(dims, spacing)
    assert box.shape == (n_frames + 1, 3, 3)
    for frame_box in box:
        assert np.array_equal(frame_box, expected)

    # and the file is the one the stream writer produces, byte for byte
    assert buffered.read_bytes() == stream.read_bytes()
    assert not [p.name for p in tmp_path.iterdir() if ".tmp." in p.name]


def test_save_out_sim_writes_each_buffered_frame_exactly_once(tmp_path: Path) -> None:
    """The buffer can be written from an error path and again at the end."""
    dims, n_beads = [12, 12, 12], 10
    frames = [_walk(n_beads, dims, seed) for seed in range(7)]
    lattice = _Lattice(dims, [frames[0]])
    pdb, xtc = str(tmp_path / "START.pdb"), tmp_path / "traj.xtc"

    lattice_utils.start_xtc_file(
        lattice, SPACING, pdb_filename=pdb, xtc_filename=str(xtc)
    )
    master = None
    for positions in frames[1:4]:
        lattice.chains[1].positions = positions
        master = lattice_utils.update_master_traj(lattice, SPACING, master, pdb)

    lattice_utils.save_out_sim(master, str(xtc))
    assert len(_read_xtc(xtc)[0]) == 4
    after_first = xtc.read_bytes()

    lattice_utils.save_out_sim(master, str(xtc))  # nothing new: file untouched
    assert xtc.read_bytes() == after_first

    for positions in frames[4:]:
        lattice.chains[1].positions = positions
        master = lattice_utils.update_master_traj(lattice, SPACING, master, pdb)
    lattice_utils.save_out_sim(master, str(xtc))
    lattice_utils.save_out_sim(master, str(xtc))

    xyz, _time, step, _box = _read_xtc(xtc)
    assert list(step) == [0, 1, 2, 3, 4, 5, 6]
    assert np.allclose(xyz, np.array(frames) * SPACING * 0.1, atol=1e-3)
    assert xtc.read_bytes()[: len(after_first)] == after_first


def test_save_at_end_does_not_copy_the_buffer_to_write_it(tmp_path: Path) -> None:
    """Writing the buffer must not need a second copy of the trajectory in memory.

    The frames used to be joined into one new Trajectory (a full copy) and then
    converted again on save: 2.2-2.7 times the coordinate memory at the very end
    of the run.
    """
    dims, n_beads, n_frames = [40, 40, 40], 2000, 60
    lattice = _Lattice(dims, [_walk(n_beads, dims, 0)])
    pdb, xtc = str(tmp_path / "START.pdb"), str(tmp_path / "traj.xtc")
    lattice_utils.start_xtc_file(lattice, SPACING, pdb_filename=pdb, xtc_filename=xtc)

    master = None
    for _ in range(n_frames):
        master = lattice_utils.update_master_traj(lattice, SPACING, master, pdb)

    coordinate_bytes = n_frames * n_beads * 3 * 4  # float32, by definition
    tracemalloc.start()
    try:
        lattice_utils.save_out_sim(master, xtc)
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert len(_read_xtc(Path(xtc))[0]) == n_frames + 1
    assert peak < 0.5 * coordinate_bytes


def test_save_at_end_run_in_a_dense_sub_angstrom_box_completes(tmp_path: Path) -> None:
    """S1's minimal reproducer of K1-3: 65 beads in an 8^3 box at 0.5 A per site."""
    (tmp_path / "p.prm").write_text("ANGLE_PENALTY\tA\t0\t0\t0\nA\tA\t0\nA\t0\t0\n")
    (tmp_path / "KEYFILE.kf").write_text(
        "DIMENSIONS : 8 8 8\nPARAMETER_FILE : p.prm\nCHAIN : 65 A\nTEMPERATURE : 10\nN_STEPS : 4\n"
        "EQUILIBRATION : 0\nSEED : 1\nMOVE_CRANKSHAFT : 1.0\nXTC_FREQ : 1\nSAVE_AT_END : True\n"
        "LATTICE_TO_ANGSTROMS : 0.5\n"
    )
    package_root = Path(pimms.__file__).resolve().parents[1]
    script = package_root / "scripts" / "PIMMS"
    env = dict(os.environ, PYTHONPATH=str(package_root), OMP_NUM_THREADS="1")
    run = subprocess.run(
        [sys.executable, str(script), "-k", "KEYFILE.kf"],
        cwd=tmp_path,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
    )
    assert run.returncode == 0, run.stdout.decode()[-2000:]

    xyz, _time, step, box = _read_xtc(tmp_path / "traj.xtc")
    assert list(step) == [0, 1, 2, 3, 4]
    assert xyz.shape == (5, 65, 3)
    for frame_box in box:
        assert np.array_equal(frame_box, _expected_box([8, 8, 8], 0.5))


# ---------------------------------------------------------------------------
# O1-5: a stale output that cannot be removed is reported
# ---------------------------------------------------------------------------
def test_remove_files_reports_what_it_could_not_remove(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    # the warning also goes to the log.txt of the working directory when there
    # is one, so the working directory must be ours
    monkeypatch.chdir(tmp_path)
    log = tmp_path / "log.txt"
    log.write_text("PIMMS Simulation\n")

    removable = tmp_path / "ENERGY.dat"
    removable.write_text("0\t0\n")
    stuck = tmp_path / "CLUSTER_RG.dat"
    stuck.mkdir()  # a directory under an output's name
    absent = tmp_path / "RG.dat"

    removed = IO_utils.remove_files([str(removable), str(stuck), str(absent)])

    assert removed == [str(removable)]
    assert not removable.exists() and stuck.exists()
    out = " ".join(capsys.readouterr().out.split())
    assert "WARNING" in out
    assert "CLUSTER_RG.dat" in out
    assert "ENERGY.dat" not in out and "RG.dat]" not in out.replace(
        "CLUSTER_RG.dat]", ""
    )

    # ...and in the log, once
    logged = " ".join(log.read_text().split())
    assert logged.startswith("PIMMS Simulation")
    assert logged.count("Could not delete the existing output file") == 1
    assert "CLUSTER_RG.dat" in logged

    # with no log.txt, none is created as a side effect
    log.unlink()
    IO_utils.remove_files([str(stuck)])
    assert not log.exists()


def test_remove_files_clears_a_dangling_symlink(tmp_path: Path) -> None:
    link = tmp_path / "ENERGY.dat"
    link.symlink_to(tmp_path / "gone.dat")
    assert IO_utils.remove_files([str(link)]) == [str(link)]
    assert not os.path.lexists(link)


# ---------------------------------------------------------------------------
# O1-4: files created with 'w' replace a symbolic link, they do not write through it
# ---------------------------------------------------------------------------
def test_new_output_files_replace_a_symlink_instead_of_writing_through_it(
    tmp_path: Path,
) -> None:
    elsewhere = tmp_path / "previous_segment"
    elsewhere.mkdir()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    precious = "the previous segment's file\n"

    dims = [8, 8, 8]
    lattice = _Lattice(dims, [_walk(4, dims, 0)])

    def linked(name: str) -> tuple[Path, Path]:
        target = elsewhere / name
        target.write_text(precious)
        link = run_dir / name
        link.symlink_to(target)
        return link, target

    link, target = linked("START.pdb")
    writer = lattice_utils.open_xtc_writer(
        lattice, SPACING, pdb_filename=str(link), xtc_filename=str(run_dir / "traj.xtc")
    )
    lattice_utils.close_xtc_writer(writer)
    assert target.read_text() == precious
    assert not link.is_symlink() and link.read_text().startswith("CRYST1")

    link, target = linked("log.txt")
    IO_utils.wipe_file(str(link))
    assert target.read_text() == precious
    assert not link.is_symlink() and link.read_text() == ""

    link, target = linked("parameters_used.prm")
    IO_utils.write_list_to_file(["a\n", "b\n"], str(link))
    assert target.read_text() == precious
    assert not link.is_symlink() and link.read_text() == "a\nb\n"

    link, target = linked("lattice.pdb")
    pdb_utils.initialize_pdb_file(dims, SPACING, str(link))
    assert target.read_text() == precious
    assert not link.is_symlink() and link.read_text().startswith("CRYST1")


# ---------------------------------------------------------------------------
# O1-9: elapsed time is measured on the monotonic clock
# ---------------------------------------------------------------------------
class _SteppingClock:
    """A wall clock and a monotonic clock that we move by hand."""

    def __init__(self, wall: datetime, monotonic: float) -> None:
        self.wall = wall
        self.mono = monotonic

    def now(self) -> datetime:
        return self.wall

    def monotonic(self) -> float:
        return self.mono


class _Acceptance:
    def total_attempted_moves(self) -> int:
        return 8000


def test_performance_survives_the_wall_clock_going_backwards(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows: list[tuple] = []
    monkeypatch.setattr(
        analysis_general.analysis_IO,
        "write_performance",
        lambda *args: rows.append(args),
    )
    monkeypatch.setattr(
        analysis_general.pimmslogger, "log_status", lambda *a, **k: None
    )
    monkeypatch.setattr(analysis_general, "_MONOTONIC_ANCHOR", {}, raising=False)

    start = datetime(2026, 11, 1, 1, 30, 0)
    clock = _SteppingClock(wall=start + timedelta(seconds=10), monotonic=500.0)
    monkeypatch.setattr(analysis_general, "datetime", clock)
    monkeypatch.setattr(analysis_general, "time", clock, raising=False)

    # 10 s into the run
    analysis_general.evaluate_performance(20, start, 1000, 0, _Acceptance())
    # 30 s later the wall clock has been set back an hour (the end of daylight saving)
    clock.mono += 30.0
    clock.wall = clock.wall + timedelta(seconds=30) - timedelta(hours=1)
    analysis_general.evaluate_performance(80, start, 1000, 0, _Acceptance())

    assert rows[0][4] == "00:00:10"
    _step, _phase, steps_per_second, moves_per_second, elapsed, remaining = rows[1]
    assert elapsed == "00:00:40"  # 10 s + 30 s, whatever the wall clock says
    assert steps_per_second == pytest.approx(80 / 40.0)
    assert moves_per_second == pytest.approx(8000 / 40.0)
    assert remaining == "00:07:40"  # 920 steps at 2 steps/s


# ---------------------------------------------------------------------------
# E1-3: integer energies are written exactly, in the same format
# ---------------------------------------------------------------------------
def test_integer_energies_are_written_exactly_and_in_the_old_format(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    energy_file, quench_file = tmp_path / "ENERGY.dat", tmp_path / "QUENCH.dat"
    monkeypatch.setattr(CONFIG, "OUTNAME_ENERGY", str(energy_file))
    monkeypatch.setattr(CONFIG, "QUENCHFILE_NAME", str(quench_file))

    # beyond 2^53 a double cannot hold every integer; the file must
    big = [2**53 + 1, -(2**53 + 1), 2**63 - 1, 10**30 + 7]
    for step, energy in enumerate(big):
        analysis_IO.write_energy(step, energy)
        analysis_IO.write_quench_file(step, 50.0, energy)
    assert energy_file.read_text() == "".join(
        "%d\t%d.0000\n" % (s, e) for s, e in enumerate(big)
    )
    assert quench_file.read_text() == "".join(
        "%d\t50\t%d.0000\n" % (s, e) for s, e in enumerate(big)
    )

    # every integer that was exact before is written with the same characters as
    # before, whatever integer type it arrives as
    rng = np.random.default_rng(3)
    exact = [
        0,
        1,
        -1,
        9,
        -9,
        99999,
        -99999,
        123456789,
        -(2**31),
        2**31,
        2**53,
        -(2**53),
    ]
    exact += [int(v) for v in rng.integers(-(10**12), 10**12, size=200)]
    for value in exact:
        for typed in (value, np.int64(value)):
            assert analysis_IO.format_energy(typed) == "%10.4f" % value
    assert analysis_IO.format_energy(np.int32(-123)) == " -123.0000"

    # a non-integer energy is formatted as it always was
    assert analysis_IO.format_energy(-12.5) == "  -12.5000"
    assert analysis_IO.format_energy(np.float64(3.25)) == "    3.2500"


# ---------------------------------------------------------------------------
# A2-6: AUTOCENTER in a hardwall box keeps every bead inside the box
# ---------------------------------------------------------------------------
def _lopsided_chain() -> list[list[int]]:
    """Seven beads up the x = 0 wall, then a tail straight across to x = 9."""
    return [[0, y, 5] for y in range(7)] + [[x, 6, 5] for x in range(1, 10)]


def test_centring_shift_is_clamped_under_hardwall() -> None:
    dims = [10, 10, 10]
    chain = _lopsided_chain()

    # the centre of mass (x ~ 2.8) goes to the box centre, which carries the far
    # end of the tail two sites through the wall
    free = np.array(lattice_utils.center_positions(chain, dims))
    assert free[:, 0].max() == 11

    clamped = np.array(lattice_utils.center_positions(chain, dims, hardwall=True))
    assert clamped.min() >= 0 and (clamped < np.array(dims)).all()
    # a rigid whole-site translation: the conformation is untouched
    shift = clamped - np.array(chain)
    assert (shift == shift[0]).all()
    # y and z did fit, and are centred exactly as before
    assert np.array_equal(clamped[:, 1:], free[:, 1:])

    # a chain that fits when centred is not moved by the clamp at all
    compact = [[1, 1, 1], [2, 1, 1], [2, 2, 1], [3, 2, 1]]
    assert lattice_utils.center_positions(
        compact, dims, hardwall=True
    ) == lattice_utils.center_positions(compact, dims)

    # both walls, both directions
    assert lattice_utils.clamp_positions_to_box([[-2, 0], [-1, 0], [0, 0]], [5, 5]) == [
        [0, 0],
        [1, 0],
        [2, 0],
    ]
    assert lattice_utils.clamp_positions_to_box([[3, 4], [4, 4], [5, 4]], [5, 5]) == [
        [2, 4],
        [3, 4],
        [4, 4],
    ]


@pytest.mark.parametrize("save_at_end", [False, True])
def test_autocenter_with_hardwall_writes_no_bead_outside_the_box(
    tmp_path: Path, save_at_end: bool
) -> None:
    dims = [10, 10, 10]
    lattice = _Lattice(dims, [_lopsided_chain()], hardwall=True)
    pdb, xtc = str(tmp_path / "START.pdb"), str(tmp_path / "traj.xtc")

    if save_at_end:
        lattice_utils.start_xtc_file(
            lattice, SPACING, pdb_filename=pdb, xtc_filename=xtc, autocenter=True
        )
        master = lattice_utils.update_master_traj(
            lattice, SPACING, None, pdb, autocenter=True
        )
        lattice_utils.save_out_sim(master, xtc)
    else:
        writer = lattice_utils.open_xtc_writer(
            lattice, SPACING, pdb_filename=pdb, xtc_filename=xtc, autocenter=True
        )
        lattice_utils.write_xtc_frame(writer, lattice, SPACING, autocenter=True)
        lattice_utils.close_xtc_writer(writer)

    sites = np.rint(_read_xtc(Path(xtc))[0] * 10 / SPACING).astype(int)
    assert sites.shape == (2, 16, 3)
    assert sites.min() >= 0 and (sites < np.array(dims)).all()

    atom_lines = [
        line for line in Path(pdb).read_text().splitlines() if line.startswith("ATOM")
    ]
    pdb_sites = np.rint(
        np.array(
            [
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                for line in atom_lines
            ]
        )
        / SPACING
    ).astype(int)
    assert np.array_equal(pdb_sites, sites[0])  # START.pdb is frame 0

    # periodic boundaries: the same chain is centred exactly, tail outside the box
    periodic = _Lattice(dims, [_lopsided_chain()], hardwall=False)
    xyz, _box = lattice_utils._lattice_frame_xyz_and_box(
        periodic, SPACING, autocenter=True
    )
    assert np.rint(xyz * 10 / SPACING).astype(int)[0, :, 0].max() == 11


# ---------------------------------------------------------------------------
# K1-8: AUTOCENTER with more than one chain says that it is being ignored
# ---------------------------------------------------------------------------
def test_autocenter_with_several_chains_warns_once_and_names_the_keyword(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.chdir(tmp_path)  # no log.txt here, so nothing is logged
    monkeypatch.setattr(
        lattice_utils, "_AUTOCENTER_IGNORED_WARNED", False, raising=False
    )
    dims = [10, 10, 10]
    chain_a, chain_b = _lopsided_chain(), [[5, 0, 0], [5, 1, 0]]
    lattice = _Lattice(dims, [chain_a, chain_b])

    writer = lattice_utils.open_xtc_writer(
        lattice,
        SPACING,
        pdb_filename="START.pdb",
        xtc_filename="traj.xtc",
        autocenter=True,
    )
    for _ in range(3):
        lattice_utils.write_xtc_frame(writer, lattice, SPACING, autocenter=True)
    lattice_utils.close_xtc_writer(writer)

    out = " ".join(capsys.readouterr().out.split())
    assert out.count("AUTOCENTER : True is ignored") == 1
    assert "2 chains" in out

    # and the frames are the uncentred positions
    sites = np.rint(_read_xtc(tmp_path / "traj.xtc")[0] * 10 / SPACING).astype(int)
    assert np.array_equal(sites[0], np.array(chain_a + chain_b))

    # a single chain is centred, without a warning
    monkeypatch.setattr(
        lattice_utils, "_AUTOCENTER_IGNORED_WARNED", False, raising=False
    )
    single = _Lattice(dims, [[[1, 1, 1], [2, 1, 1]]])
    xyz, _box = lattice_utils._lattice_frame_xyz_and_box(
        single, SPACING, autocenter=True
    )
    assert "AUTOCENTER" not in capsys.readouterr().out
    assert np.rint(xyz * 10 / SPACING).astype(int)[0].tolist() == [[4, 5, 5], [5, 5, 5]]


# ---------------------------------------------------------------------------
# Review round (REV_F2): failure paths of the writers
# ---------------------------------------------------------------------------
def _buffered_run(
    tmp_path: Path, n_frames: int, n_beads: int = 30
) -> tuple[object, Path, Path, int]:
    """A SAVE_AT_END trajectory with n_frames buffered, and its streamed twin.

    Returns the buffer, the SAVE_AT_END file (frame 0 only so far), the file the
    stream writer produced from the same frames, and the bytes of one buffered
    frame.
    """
    dims = [12, 12, 12]
    frames = [_walk(n_beads, dims, seed) for seed in range(n_frames + 1)]

    lattice = _Lattice(dims, [frames[0]])
    stream = tmp_path / "stream.xtc"
    writer = lattice_utils.open_xtc_writer(
        lattice,
        SPACING,
        pdb_filename=str(tmp_path / "stream.pdb"),
        xtc_filename=str(stream),
    )
    for positions in frames[1:]:
        lattice.chains[1].positions = positions
        lattice_utils.write_xtc_frame(writer, lattice, SPACING)
    lattice_utils.close_xtc_writer(writer)

    lattice = _Lattice(dims, [frames[0]])
    pdb, xtc = str(tmp_path / "START.pdb"), tmp_path / "traj.xtc"
    lattice_utils.start_xtc_file(
        lattice, SPACING, pdb_filename=pdb, xtc_filename=str(xtc)
    )
    master = None
    for positions in frames[1:]:
        lattice.chains[1].positions = positions
        master = lattice_utils.update_master_traj(lattice, SPACING, master, pdb)
    return master, xtc, stream, n_beads * 3 * 4


def _scratch_files(directory: Path) -> list[str]:
    return sorted(p.name for p in directory.iterdir() if ".tmp." in p.name)


def test_failed_save_at_end_write_leaves_a_valid_file_and_a_clean_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RF2-1: the disk fills up half-way through adding the buffered frames."""
    master, xtc, stream, _frame_bytes = _buffered_run(tmp_path, n_frames=7)
    before = xtc.read_bytes()
    real_copy = shutil.copyfileobj

    def half_then_no_space(fin, fout, *args):  # noqa: ANN001, ANN002, ANN202
        data = fin.read()
        fout.write(data[: len(data) // 2])
        fout.flush()
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(lattice_utils.shutil, "copyfileobj", half_then_no_space)
    with pytest.raises(OSError):
        lattice_utils.save_out_sim(master, str(xtc))
    monkeypatch.setattr(lattice_utils.shutil, "copyfileobj", real_copy)

    # the half-written frames have been taken back out: the file is the valid
    # one-frame trajectory it was, and nothing is left lying around
    assert xtc.read_bytes() == before
    assert len(_read_xtc(xtc)[0]) == 1
    assert md.load(str(xtc), top=str(tmp_path / "START.pdb")).n_frames == 1
    assert _scratch_files(tmp_path) == []

    # nothing was lost from the buffer, so the retry (the run's own clean-up
    # makes one) gives exactly the file the stream writer wrote
    lattice_utils.save_out_sim(master, str(xtc))
    assert xtc.read_bytes() == stream.read_bytes()
    assert list(_read_xtc(xtc)[2]) == list(range(8))


def test_save_at_end_writes_in_batches_and_keeps_the_batches_that_succeeded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RF2-1: the scratch file holds one batch, and a failure costs one batch."""
    master, xtc, stream, frame_bytes = _buffered_run(tmp_path, n_frames=7)
    # two frames to a batch: 7 buffered frames -> batches of 2, 2, 2, 1
    monkeypatch.setattr(
        lattice_utils.TrajectoryAccumulator, "BATCH_BYTES", 2 * frame_bytes
    )
    real_copy = shutil.copyfileobj
    scratch_frames: list[int] = []

    def fail_on_second_batch(fin, fout, *args):  # noqa: ANN001, ANN002, ANN202
        with md.formats.XTCTrajectoryFile(fin.name) as fh:
            scratch_frames.append(len(fh))
        if len(scratch_frames) == 2:
            data = fin.read()
            fout.write(data[: len(data) // 2])
            fout.flush()
            raise OSError(errno.ENOSPC, "No space left on device")
        real_copy(fin, fout, *args)

    monkeypatch.setattr(lattice_utils.shutil, "copyfileobj", fail_on_second_batch)
    with pytest.raises(OSError):
        lattice_utils.save_out_sim(master, str(xtc))

    # frame 0 and the first batch are in the file, whole; the second is not
    assert scratch_frames == [2, 2]
    assert list(_read_xtc(xtc)[2]) == [0, 1, 2]
    assert xtc.read_bytes() == stream.read_bytes()[: xtc.stat().st_size]
    assert _scratch_files(tmp_path) == []

    # the retry writes the five frames still buffered, and only those
    lattice_utils.save_out_sim(master, str(xtc))
    monkeypatch.setattr(lattice_utils.shutil, "copyfileobj", real_copy)
    assert scratch_frames == [2, 2, 2, 2, 1]
    assert xtc.read_bytes() == stream.read_bytes()


@pytest.mark.parametrize("save_at_end", [False, True])
def test_dangling_trajectory_symlink_is_replaced_not_written_through(
    tmp_path: Path, save_at_end: bool
) -> None:
    """RF2-4: a link whose target does not exist is still a link."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    target = elsewhere / "gone.xtc"
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    link = run_dir / "traj.xtc"
    link.symlink_to(target)

    dims = [8, 8, 8]
    lattice = _Lattice(dims, [_walk(4, dims, 0)])
    pdb = str(run_dir / "START.pdb")
    if save_at_end:
        lattice_utils.start_xtc_file(
            lattice, SPACING, pdb_filename=pdb, xtc_filename=str(link)
        )
    else:
        writer = lattice_utils.open_xtc_writer(
            lattice, SPACING, pdb_filename=pdb, xtc_filename=str(link)
        )
        lattice_utils.close_xtc_writer(writer)

    assert not target.exists()
    assert list(elsewhere.iterdir()) == []
    assert not link.is_symlink() and len(_read_xtc(link)[0]) == 1


def test_a_pdb_that_cannot_be_moved_into_place_leaves_no_temporary_file(
    tmp_path: Path,
) -> None:
    """RF2-5: the rename is inside the clean-up too."""
    dims = [8, 8, 8]
    lattice = _Lattice(dims, [_walk(4, dims, 0)])
    blocked = tmp_path / "START.pdb"
    blocked.mkdir()  # a directory sits under the output's name

    with pytest.raises(OSError):
        lattice_utils.write_topology_pdb(lattice, SPACING, str(blocked))
    with pytest.raises(OSError):
        pdb_utils.write_positions_to_file([[0, 0, 0], [0, 0, 1]], str(blocked), SPACING)

    assert [p.name for p in tmp_path.iterdir()] == ["START.pdb"]
    assert blocked.is_dir() and list(blocked.iterdir()) == []


def _dead_pid() -> int:
    """The pid of a process that has just exited."""
    finished = subprocess.Popen([sys.executable, "-c", "pass"])
    finished.wait()
    return finished.pid


def test_scratch_files_of_killed_runs_are_listed_unless_their_process_is_alive(
    tmp_path: Path,
) -> None:
    """RF2-5: what a SIGKILL inside a topology or SAVE_AT_END write leaves behind."""
    dead, alive = _dead_pid(), os.getppid()
    stale = [
        "traj.xtc.tmp.%d" % dead,
        "START.pdb.tmp.%d" % dead,
        "eq_traj.xtc.tmp.%d" % dead,
        "eq_START.pdb.tmp.%d" % dead,
    ]
    kept = [
        "traj.xtc.tmp.%d" % alive,  # another run's write in progress
        "traj.xtc.tmp.backup",  # not ours: no pid
        "traj.xtc",
        "START.pdb",
        "notes.xtc.tmp.%d" % dead,  # not an output name
    ]
    for name in stale + kept:
        (tmp_path / name).write_text("x")

    assert sorted(lattice_utils.stale_writer_temporaries()) == sorted(stale)


def test_a_new_run_removes_the_scratch_files_of_killed_runs(tmp_path: Path) -> None:
    """RF2-5, end to end: start-up clears them, and only them."""
    dead = _dead_pid()
    stale = [
        "traj.xtc.tmp.%d" % dead,
        "START.pdb.tmp.%d" % dead,
        "eq_traj.xtc.tmp.%d" % dead,
    ]
    in_use = "traj.xtc.tmp.%d" % os.getpid()  # this test process is alive throughout
    for name in stale + [in_use]:
        (tmp_path / name).write_text("left by a killed run")

    (tmp_path / "p.prm").write_text("ANGLE_PENALTY\tA\t0\t0\t0\nA\tA\t0\nA\t0\t0\n")
    (tmp_path / "KEYFILE.kf").write_text(
        "DIMENSIONS : 10 10 10\nPARAMETER_FILE : p.prm\nCHAIN : 2 AAAA\nTEMPERATURE : 10\n"
        "N_STEPS : 4\nEQUILIBRATION : 0\nSEED : 1\nMOVE_CRANKSHAFT : 1.0\nXTC_FREQ : 1\n"
        "SAVE_AT_END : True\n"
    )
    package_root = Path(pimms.__file__).resolve().parents[1]
    env = dict(os.environ, PYTHONPATH=str(package_root), OMP_NUM_THREADS="1")
    run = subprocess.run(
        [sys.executable, str(package_root / "scripts" / "PIMMS"), "-k", "KEYFILE.kf"],
        cwd=tmp_path,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=300,
    )
    assert run.returncode == 0, run.stdout.decode()[-2000:]

    assert _scratch_files(tmp_path) == [in_use]
    assert len(_read_xtc(tmp_path / "traj.xtc")[0]) == 5


@pytest.mark.filterwarnings("ignore:Unlikely unit cell vectors")
def test_rebuilt_first_frame_of_a_dense_box_keeps_its_unit_cell(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RF2-6: traj.xtc vanished, nothing buffered, and mdtraj drops the CRYST1 cell."""
    # 65 beads in a (0.4 nm)^3 box: 1016 beads/nm^3, above the 1000 where mdtraj drops the cell
    dims, spacing = [8, 8, 8], 0.5
    lattice = _Lattice(dims, [_walk(65, dims, 0)])
    pdb, xtc = str(tmp_path / "START.pdb"), tmp_path / "traj.xtc"

    lattice_utils.start_xtc_file(
        lattice, spacing, pdb_filename=pdb, xtc_filename=str(xtc)
    )
    assert md.load(pdb).unitcell_vectors is None  # the condition this test is about
    expected = _expected_box(dims, spacing)

    # the cell the trajectory was started with is used, exactly
    xtc.unlink()
    lattice_utils.save_out_sim(lattice_utils.start_master_traj(pdb), str(xtc))
    xyz, _time, step, box = _read_xtc(xtc)
    assert list(step) == [0] and xyz.shape == (1, 65, 3)
    assert np.array_equal(box[0], expected)

    # and a buffer that knows nothing about the start still writes a cell: the
    # CRYST1 record itself, good to its three decimals in angstroms
    monkeypatch.setattr(lattice_utils, "_STARTED_TRAJECTORY_BOX", {}, raising=False)
    xtc.unlink()
    lattice_utils.save_out_sim(lattice_utils.start_master_traj(pdb), str(xtc))
    box = _read_xtc(xtc)[3]
    assert np.allclose(box[0], expected, rtol=0, atol=1e-4)
    assert np.array_equal(box[0] != 0, expected != 0)


def test_coordinate_outside_the_box_is_blamed_on_unwrap_not_on_the_box(
    tmp_path: Path,
) -> None:
    """RF2-7: a 9600 A box is fine; a chain written 1200 A outside it is not."""
    with pytest.raises(PDBException) as err:
        pdb_utils.format_pdb_coordinate(-1200.0, box_edge=9600.0)
    message = str(err.value)
    assert "TRAJECTORY_PBC_UNWRAP" in message and "AUTOCENTER" in message
    assert "Reduce DIMENSIONS" not in message

    with pytest.raises(PDBException) as err:
        pdb_utils.format_pdb_coordinate(10500.0, box_edge=9600.0)
    assert "TRAJECTORY_PBC_UNWRAP" in str(err.value)

    # inside the box, it is the box
    with pytest.raises(PDBException) as err:
        pdb_utils.format_pdb_coordinate(10500.0, box_edge=10950.0)
    message = str(err.value)
    assert "Reduce DIMENSIONS or LATTICE_TO_ANGSTROMS" in message
    assert "TRAJECTORY_PBC_UNWRAP" not in message

    # through the writer: a 7300 A box, with a whole chain reaching 330 sites
    # (1204.5 A) below the lower face
    dims = [2000, 10, 10]
    lattice = _Lattice(dims, [[[-330, 0, 0], [-329, 0, 0]], [[5, 5, 5], [6, 5, 5]]])
    with pytest.raises(PDBException) as err:
        lattice_utils.open_xtc_writer(
            lattice,
            SPACING,
            pdb_filename=str(tmp_path / "START.pdb"),
            xtc_filename=str(tmp_path / "traj.xtc"),
            unwrap=True,
        )
    message = str(err.value)
    assert "-1204.500" in message and "TRAJECTORY_PBC_UNWRAP" in message
    assert "Reduce DIMENSIONS" not in message
    assert list(tmp_path.iterdir()) == []


def test_replacing_an_output_ignores_the_old_files_permissions(tmp_path: Path) -> None:
    """Documented behaviour: a write-protected output of an earlier run is replaced."""
    dims = [8, 8, 8]
    lattice = _Lattice(dims, [_walk(4, dims, 0)])
    pdb, dat = tmp_path / "START.pdb", tmp_path / "ENERGY.dat"
    for old in (pdb, dat):
        old.write_text("an earlier run\n")
        old.chmod(0o444)

    lattice_utils.write_topology_pdb(lattice, SPACING, str(pdb))
    assert pdb.read_text().startswith("CRYST1")
    assert IO_utils.remove_files([str(dat)]) == [str(dat)]
