"""
Regression tests for the deep-audit findings on lemonade's loading, topology
and storage code (A2-5, A2-6 load side, A2-7, A2-11 tie, A2-12, A2-13, P2-4, D9).

The trajectories are written here from known integer lattice coordinates, so
every expectation is derived from those coordinates and never from the loader.
"""

import os
import re
import shutil
import tempfile
import tracemalloc
import types
import warnings

import mdtraj as md
import numpy as np
import pytest

from pimms import CONFIG
from pimms import lemonade
from pimms.lemonade import _load
from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology

_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _straight_chains(rng, n_frames, n_chains, length, dims, wrap=True):
    """Lattice coordinates of straight chains lying along x.

    Every chain keeps its own (y, z) row, so no two beads share a site, and is
    put at a fresh random x in every frame. Bonds are single steps along x.

    Parameters
    ----------
    rng : numpy.random.Generator
        Source of the per-frame x origins.
    n_frames, n_chains, length : int
        Size of the system.
    dims : tuple of int
        Box extent; two entries give a 2D system (z = 0).
    wrap : bool, optional
        Wrap x into the box (default), or leave chains whole across the
        boundary as ``TRAJECTORY_PBC_UNWRAP`` writes them.

    Returns
    -------
    numpy.ndarray
        ``(n_frames, n_chains * length, 3)`` int32 coordinates.
    """
    out = np.zeros((n_frames, n_chains * length, 3), dtype=np.int32)
    for f in range(n_frames):
        for c in range(n_chains):
            x = rng.integers(0, dims[0]) + np.arange(length)
            block = out[f, c * length : (c + 1) * length]
            block[:, 0] = x % dims[0] if wrap else x
            block[:, 1] = c % dims[1]
            block[:, 2] = (c // dims[1]) % dims[2] if len(dims) == 3 else 0
    return out


def _write(
    directory,
    positions,
    sequences,
    dims,
    spacing=3.65,
    labels=None,
    xtc=True,
    xtc_box=True,
):
    """Write a PIMMS-style ``START.pdb`` (and ``traj.xtc``) for known coordinates.

    Parameters
    ----------
    directory : path-like
        Where to write.
    positions : numpy.ndarray
        ``(n_frames, n_beads, 3)`` integer lattice coordinates.
    sequences : list of str
        One sequence per chain.
    dims : tuple of int
        Box extent in lattice units (2 or 3 entries).
    spacing : float, optional
        Lattice spacing in angstroms.
    labels : list of str, optional
        PDB chain identifier per chain (default: one per distinct sequence, in
        order of first appearance, as PIMMS assigns them).
    xtc : bool, optional
        Also write the trajectory (default ``True``).
    xtc_box : bool, optional
        Record the box in the XTC (default ``True``). ``False`` writes the
        trajectory without a unit cell, as PIMMS before 1.0.8 did under
        ``SAVE_AT_END`` at a sub-angstrom spacing.

    Returns
    -------
    tuple of str
        ``(xtc_path or None, pdb_path)``.
    """
    directory = str(directory)
    os.makedirs(directory, exist_ok=True)
    if labels is None:
        seen = {}
        labels = [_IDS[seen.setdefault(s, len(seen))] for s in sequences]
    box = [d * spacing for d in dims] + [spacing] * (3 - len(dims))
    lines = [
        "CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n" % tuple(box)
    ]
    bead = 0
    for seq, label in zip(sequences, labels):
        for r, ch in enumerate(seq):
            x, y, z = positions[0, bead] * spacing
            lines.append(
                "ATOM  %5d  CA  %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
                % (
                    (bead + 1) % 100000,
                    CONFIG.ONE_TO_THREE[ch],
                    label,
                    (r + 1) % 10000,
                    x,
                    y,
                    z,
                )
            )
            bead += 1
        lines.append("TER\n")
    lines.append("END\n")
    pdb_path = os.path.join(directory, "START.pdb")
    with open(pdb_path, "w") as fh:
        fh.writelines(lines)
    if not xtc:
        return None, pdb_path
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        top = md.load(pdb_path).topology
    n = positions.shape[0]
    traj = md.Trajectory(
        (positions * (spacing * 0.1)).astype(np.float32),
        top,
        time=np.arange(n, dtype=np.float32),
        unitcell_lengths=np.tile(np.array(box) * 0.1, (n, 1)) if xtc_box else None,
        unitcell_angles=np.full((n, 3), 90.0) if xtc_box else None,
    )
    xtc_path = os.path.join(directory, "traj.xtc")
    traj.save_xtc(xtc_path)
    return xtc_path, pdb_path


def _quiet_load(*args, **kwargs):
    """Call ``lemonade.load`` and fail on any lemonade warning.

    Parameters
    ----------
    *args, **kwargs
        Passed to ``lemonade.load``.

    Returns
    -------
    LatticeTrajectory
        The loaded trajectory.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(*args, **kwargs)
    ours = [str(w.message) for w in caught if "lemonade.load" in str(w.message)]
    assert ours == []
    return traj


# -- P2-4: only the selected frames are read and held -------------------------


@pytest.mark.parametrize("chunk_bytes", [12 * 30 * 3, 32 * 1024 * 1024])
def test_frame_selection_returns_the_frames_a_numpy_slice_selects(
    tmp_path, monkeypatch, chunk_bytes
):
    rng = np.random.default_rng(1)
    dims = (12, 6, 5)
    truth = _straight_chains(rng, 37, 6, 5, dims)
    xtc, pdb = _write(tmp_path, truth, ["AAAAA"] * 6, dims)
    # a tiny block size forces many separate reads; the default reads in one
    monkeypatch.setattr(_load, "_CHUNK_BYTES", chunk_bytes, raising=False)
    n = truth.shape[0]
    for sl in [
        (None, None, None),
        (-1, None, None),
        (-5, None, None),
        (None, -3, None),
        (None, None, 2),
        (None, None, -1),
        (None, None, n + 10),
        (5, 2, -1),
        (3, 30, 7),
        (-n - 50, None, 3),
        (None, n + 50, 4),
        (20, 4, -3),
    ]:
        for n_frames in (None, 1, 2, 5):
            expect = np.arange(n)[slice(*sl)]
            if n_frames is not None and n_frames < len(expect):
                expect = expect[
                    np.linspace(0, len(expect) - 1, n_frames).round().astype(int)
                ]
            traj = _quiet_load(
                xtc, pdb, start=sl[0], stop=sl[1], step=sl[2], n_frames=n_frames
            )
            assert traj.positions.dtype == np.int32
            assert np.array_equal(traj.positions, truth[expect]), (sl, n_frames)
            assert np.array_equal(traj.times, expect.astype(np.float64)), (sl, n_frames)
    with pytest.raises(ValueError, match="keeps no frames of a 37-frame"):
        lemonade.load(xtc, pdb, start=5, stop=2)


def test_load_memory_follows_the_frames_kept_not_the_file(tmp_path, monkeypatch):
    rng = np.random.default_rng(2)
    dims = (40, 10, 5)
    n_frames, n_beads = 2000, 500
    truth = _straight_chains(rng, n_frames, 50, 10, dims)
    xtc, pdb = _write(tmp_path, truth, ["AAAAAAAAAA"] * 50, dims)
    one_float32_copy = n_frames * n_beads * 12  # the file, decoded once
    # small working blocks, so the bound below is about the frames kept
    monkeypatch.setattr(_load, "_CHUNK_BYTES", 256 * 1024, raising=False)

    def peak(**kwargs):
        """Peak traced memory of one load, in bytes above the starting level."""
        tracemalloc.start()
        try:
            before = tracemalloc.get_traced_memory()[0]
            traj = _quiet_load(xtc, pdb, **kwargs)
            return tracemalloc.get_traced_memory()[1] - before, traj
        finally:
            tracemalloc.stop()

    # every tenth frame: the old loader decoded the whole file first, so its peak
    # could not be below one float32 copy of it
    used, traj = peak(step=10)
    assert np.array_equal(traj.positions, truth[::10])
    assert used < 0.5 * one_float32_copy
    # everything: the lattice that is kept (the same 12 bytes per bead per frame)
    # plus bounded working blocks, where the old loader held about five copies
    used, traj = peak()
    assert np.array_equal(traj.positions, truth)
    assert used < 1.5 * one_float32_copy


def test_store_adopts_the_loaders_array_only_when_asked():
    positions = np.zeros((2, 3, 3), dtype=np.int32)
    topology = Topology(["AAA"])
    copied = TrajectoryStore(positions, (4, 4, 4), 3.65, False, topology)
    assert not np.shares_memory(copied.positions, positions)
    assert positions.flags.writeable
    adopted = TrajectoryStore(positions, (4, 4, 4), 3.65, False, topology, copy=False)
    assert np.shares_memory(adopted.positions, positions)
    assert not adopted.positions.flags.writeable


# -- A2-7: the torn last frame of a running or killed run ---------------------


def test_torn_trailing_frame_is_dropped_with_a_warning(tmp_path):
    rng = np.random.default_rng(3)
    dims = (12, 6, 5)
    truth = _straight_chains(rng, 9, 6, 5, dims)
    xtc, pdb = _write(tmp_path, truth, ["AAAAA"] * 6, dims)
    with open(xtc, "rb") as fh:
        raw = fh.read()
    with md.formats.XTCTrajectoryFile(xtc) as fh:
        starts = [int(o) for o in fh.offsets] + [len(raw)]
    torn = os.path.join(str(tmp_path), "torn.xtc")

    # cut inside the body of a frame, a few bytes into its header, and one byte
    # short of the end: the complete frames are those that end at or before the cut
    for cut in (
        starts[9] - 1,
        starts[9] - 40,
        starts[8] + 6,
        starts[6] + 30,
        starts[1] + 5,
    ):
        complete = sum(1 for end in starts[1:] if end <= cut)
        assert 1 <= complete < 9
        with open(torn, "wb") as fh:
            fh.write(raw[:cut])
        with pytest.warns(
            UserWarning,
            match=rf"incomplete frame.*Recovered the {complete} "
            rf"complete frame",
        ):
            traj = lemonade.load(torn, pdb)
        assert np.array_equal(traj.positions, truth[:complete])
        # the selection counts from the last COMPLETE frame
        with pytest.warns(UserWarning, match="incomplete frame"):
            last = lemonade.load(torn, pdb, start=-1)
        assert np.array_equal(last.positions, truth[complete - 1 : complete])

    # a file cut exactly between two frames is a healthy, shorter file: no warning
    with open(torn, "wb") as fh:
        fh.write(raw[: starts[5]])
    assert np.array_equal(_quiet_load(torn, pdb).positions, truth[:5])
    # a stub shorter than the 4-byte frame marker reads as the end of the file in
    # the XTC library itself; nothing is lost, and the complete frames come back
    with open(torn, "wb") as fh:
        fh.write(raw[: starts[5] + 3])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert np.array_equal(lemonade.load(torn, pdb).positions, truth[:5])
    assert np.array_equal(_quiet_load(xtc, pdb).positions, truth)


@pytest.mark.parametrize("n_beads", [1, 5, 9, 10, 30])
def test_torn_file_is_recovered_for_both_xtc_frame_encodings(tmp_path, n_beads):
    # XTC stores nine beads or fewer uncompressed, in frames of a fixed size, and
    # mdtraj then takes the frame count from the file size - which fails outright
    # on a torn file. Ten beads or more are compressed and indexed by headers.
    rng = np.random.default_rng(n_beads)
    dims = (40, 6, 5)
    truth = _straight_chains(rng, 6, 1, n_beads, dims)
    xtc, pdb = _write(tmp_path, truth, ["A" * n_beads], dims)
    with open(xtc, "rb") as fh:
        raw = fh.read()
    with md.formats.XTCTrajectoryFile(xtc) as fh:
        starts = [int(o) for o in fh.offsets] + [len(raw)]
    torn = os.path.join(str(tmp_path), "torn.xtc")
    checked = 0
    for frame in (1, 3, 5):
        length = starts[frame + 1] - starts[frame]
        # the first bytes of the header, then every ninth byte of the frame
        for extra in sorted(
            set(range(4, 16)) | set(range(16, length, 9)) | {length - 1}
        ):
            with open(torn, "wb") as fh:
                fh.write(raw[: starts[frame] + extra])
            with pytest.warns(
                UserWarning,
                match=rf"incomplete frame.*Recovered the {frame} complete frame",
            ):
                traj = lemonade.load(torn, pdb)
            assert np.array_equal(traj.positions, truth[:frame])
            checked += 1
        # selections count within the complete frames, in any order
        with pytest.warns(UserWarning, match="incomplete frame"):
            back = lemonade.load(torn, pdb, start=-1, step=-2)
        assert np.array_equal(back.positions, truth[:frame][-1::-2])
    assert checked >= 40
    # the healthy file, and a healthy shorter one, take the direct path silently
    for a, b, c in [(None, None, None), (-2, None, None), (None, None, -1), (1, 5, 2)]:
        assert np.array_equal(
            _quiet_load(xtc, pdb, start=a, stop=b, step=c).positions,
            truth[slice(a, b, c)],
        )
    with open(torn, "wb") as fh:
        fh.write(raw[: starts[4]])
    assert np.array_equal(_quiet_load(torn, pdb).positions, truth[:4])


def test_read_plan_fetches_exactly_the_selection_in_bounded_reads():
    rng = np.random.default_rng(11)
    cases = [
        np.linspace(0, 4000, 2667).round().astype(int)
    ]  # dense, uneven: gaps 1 and 2
    for _ in range(60):
        n_total = int(rng.integers(1, 200))
        wanted = np.sort(
            rng.integers(0, n_total, int(rng.integers(1, 80)))
        )  # with repeats
        cases.append(wanted)
    cases.append(np.arange(3, 4000, 7))
    for wanted in cases:
        for chunk, max_gap, seekable in [
            (100, 2, True),
            (1, 1, True),
            (7, 3, True),
            (64, 1, False),
        ]:
            got = np.full(len(wanted), -1)
            file_position = 0
            n_reads = 0
            for first, n_read, stride, lo, hi in _load._read_plan(
                wanted, chunk, min(chunk, max_gap), seekable
            ):
                if not seekable:
                    assert first == file_position  # no seeking: strictly forward
                read = first + stride * np.arange(
                    n_read
                )  # the frames this read decodes
                assert n_read <= chunk and (n_read - 1) * stride + 1 <= max(
                    chunk, n_read * stride
                )
                assert stride == 1 or hi - lo == n_read
                got[lo:hi] = read[(wanted[lo:hi] - first) // stride]
                file_position = first + n_read
                n_reads += 1
            assert np.array_equal(got, wanted)
        # a dense uneven selection is read in spans, not a frame or two at a time
        if len(wanted) == 2667:
            reads = list(_load._read_plan(wanted, 100, 2, True))
            assert len(reads) <= 4001 // 100 + 2
    # an even stride is still one strided read per block: nothing in between is decoded
    reads = list(_load._read_plan(np.arange(3, 4000, 7), 100, 2, True))
    assert all(stride == 7 for _f, _n, stride, _lo, _hi in reads) and len(reads) == 6


def test_dense_uneven_selection_matches_the_known_frames(tmp_path, monkeypatch):
    rng = np.random.default_rng(12)
    dims = (12, 6, 5)
    truth = _straight_chains(rng, 61, 6, 5, dims)
    xtc, pdb = _write(tmp_path, truth, ["AAAAA"] * 6, dims)
    for per_read in (1, 5000, 10**9):  # never merge, default, always merge
        monkeypatch.setattr(_load, "_BEADS_PER_READ", per_read, raising=False)
        for n_frames in (2, 7, 40, 41, 55, 60):
            expect = np.linspace(0, 60, n_frames).round().astype(int)
            traj = _quiet_load(xtc, pdb, n_frames=n_frames)
            assert np.array_equal(traj.positions, truth[expect])
            assert np.array_equal(traj.times, expect.astype(np.float64))
            back = _quiet_load(xtc, pdb, step=-1, n_frames=n_frames)
            assert np.array_equal(back.positions, truth[::-1][expect])


# -- A2-5: a spacing that is a whole multiple of the assumed one --------------


@pytest.mark.parametrize("dims", [(12, 6, 5), (14, 9)])
def test_integer_multiple_spacing_is_caught_by_the_bond_check(tmp_path, dims):
    rng = np.random.default_rng(4)
    truth = _straight_chains(rng, 4, 6, 5, dims)
    xtc, pdb = _write(tmp_path, truth, ["AAAAA"] * 6, dims, spacing=7.3)
    # no keyfile: lemonade assumes 3.65 A, every coordinate lands exactly on an
    # even site, and the round-off residual is zero
    with pytest.warns(
        UserWarning,
        match=r"24 of 24 bonds.*lengths 2 to 2.*probably "
        r"7\.3 A.*spacing=7\.3",
    ):
        doubled = lemonade.load(xtc, pdb)
    assert np.array_equal(doubled.positions, 2 * truth)
    right = _quiet_load(xtc, pdb, spacing=7.3)
    assert right.dimensions == dims
    assert np.array_equal(right.positions, truth)


@pytest.mark.parametrize(
    "dims,spacing,wrap",
    [
        ((12, 6, 5), 3.65, True),
        ((12, 6, 5), 3.65, False),
        ((14, 9), 3.65, True),
        ((14, 9), 3.65, False),
        ((12, 6, 5), 0.05, True),
        ((12, 6, 5), 100.0, True),
        ((3, 6, 5), 3.65, True),
        ((2, 6, 5), 1.0, True),
    ],
)
def test_bond_check_is_silent_on_correct_loads(tmp_path, dims, spacing, wrap):
    rng = np.random.default_rng(5)
    length = min(5, dims[0])
    truth = _straight_chains(rng, 3, 6, length, dims, wrap=wrap)
    xtc, pdb = _write(tmp_path, truth, ["A" * length] * 6, dims, spacing=spacing)
    for hardwall in (False,) if not wrap or dims[0] < 5 else (False, True):
        frames = truth
        if hardwall:  # nothing may cross a wall: keep chains inside
            frames = truth.copy()
            frames[..., 0] = np.tile(np.arange(length), 6)
            xtc, pdb = _write(
                tmp_path / "hw", frames, ["A" * length] * 6, dims, spacing=spacing
            )
        traj = _quiet_load(xtc, pdb, spacing=spacing, hardwall=hardwall)
        assert np.array_equal(traj.positions[..., 0], frames[..., 0] % dims[0])


def test_bond_check_has_nothing_to_test_for_single_bead_chains(tmp_path):
    rng = np.random.default_rng(6)
    dims = (12, 6, 5)
    truth = _straight_chains(rng, 3, 20, 1, dims)
    xtc, pdb = _write(tmp_path, truth, ["A"] * 20, dims, spacing=7.3)
    # no bonds, so the doubled lattice cannot be seen: the load is silent (this
    # is the documented limit of the check, not a claim that the load is right)
    doubled = _quiet_load(xtc, pdb)
    assert np.array_equal(doubled.positions, 2 * truth)
    assert (
        _load._check_unit_bonds(truth[0], np.arange(21), dims, 3.65, None, 0.0) is False
    )


# -- A2-6: beads outside a hardwall box ----------------------------------------


def test_hardwall_frame_outside_the_box_is_translated_not_wrapped(tmp_path):
    dims = (10, 6, 5)
    chain = np.zeros((3, 6, 3), dtype=np.int32)
    chain[0, :, 0] = np.arange(6) + 2  # inside
    chain[1, :, 0] = np.arange(6) - 3  # sticks out through the low wall
    chain[2, :, 0] = np.arange(6) + 7  # sticks out through the high wall
    chain[..., 1] = 2
    chain[..., 2] = 1
    xtc, pdb = _write(tmp_path, chain, ["AAAAAA"], dims)
    with pytest.warns(
        UserWarning,
        match=r"2 of 3 frames have beads outside the hardwall box.*AUTOCENTER.*"
        r"translated rigidly.*not meaningful in any frame",
    ):
        traj = lemonade.load(xtc, pdb, hardwall=True)
    # the chain is whole in every frame: pushed back just far enough to fit
    expect = chain.copy()
    expect[1, :, 0] = np.arange(6)
    expect[2, :, 0] = np.arange(6) + 4
    assert np.array_equal(traj.positions, expect)
    assert np.allclose(traj.end_to_end_distance()[:, 0], 5.0)
    # the same file read as periodic is wrapped, as it always was, without a word
    periodic = _quiet_load(xtc, pdb, hardwall=False)
    assert np.array_equal(periodic.positions[..., 0], chain[..., 0] % 10)


def test_hardwall_frame_wider_than_the_box_is_wrapped_and_says_so(tmp_path):
    dims = (10, 6, 5)
    chain = np.zeros((1, 12, 3), dtype=np.int32)
    chain[0, :, 0] = np.arange(12) - 1
    xtc, pdb = _write(tmp_path, chain, ["A" * 12], dims)
    with pytest.warns(
        UserWarning,
        match=r"outside the hardwall box.*1 of them wider than the box.*"
        r"check hardwall= and dimensions=",
    ) as caught:
        traj = lemonade.load(xtc, pdb, hardwall=True)
    assert not any("AUTOCENTER" in str(w.message) for w in caught)
    assert np.array_equal(traj.positions[0, :, 0], (np.arange(12) - 1) % 10)


def test_hardwall_out_of_box_with_several_chains_is_a_wrong_box_not_autocenter(
    tmp_path,
):
    # AUTOCENTER only ever acted on a single chain. Three chains in a 12-wide box,
    # loaded with dimensions that are too small: the coordinates are wrapped, as
    # they always were, and the message points at the box, not at AUTOCENTER.
    true_dims = (12, 6, 5)
    truth = np.zeros((2, 12, 3), dtype=np.int32)
    for c in range(3):
        truth[:, 4 * c : 4 * c + 4, 0] = np.arange(4) + 3 * c + 2  # x up to 11
        truth[:, 4 * c : 4 * c + 4, 1] = c
    xtc, pdb = _write(tmp_path, truth, ["AAAA"] * 3, true_dims)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(xtc, pdb, dimensions=(8, 6, 5), hardwall=True)
    outside = [
        str(w.message) for w in caught if "outside the hardwall box" in str(w.message)
    ]
    assert len(outside) == 1
    assert "check hardwall= and dimensions=" in outside[0]
    assert "AUTOCENTER" not in outside[0] and "translated" not in outside[0]
    assert np.array_equal(traj.positions[..., 0], truth[..., 0] % 8)
    # the right box needs no warning of this kind
    assert np.array_equal(_quiet_load(xtc, pdb, hardwall=True).positions, truth)


# -- A2-12: 62 or more chain types ---------------------------------------------


def _sixty_two_labels():
    """A topology that uses all 62 PDB identifiers.

    Returns
    -------
    tuple
        ``(sequences, pdb_types)``: four chains of one sequence under the first
        label, 60 single chains, and three chains of three different sequences
        under the shared last label (type 61).
    """
    sequences = (
        ["AAAA"] * 4 + ["AG" + "A" * (k + 1) for k in range(60)] + ["GGA", "GGG", "GAG"]
    )
    pdb_types = [0] * 4 + list(range(1, 61)) + [61, 61, 61]
    return sequences, pdb_types


def test_only_the_shared_last_identifier_may_be_split_by_a_keyfile():
    sequences, pdb_types = _sixty_two_labels()
    topology = Topology(sequences, chain_types=pdb_types)
    singles = [[1, s] for s in sequences[4:]]
    # the run's own lines: one type for the first label, the last label split in three
    good = topology.with_keyfile_types([[4, "AAAA"]] + singles, labelled=True)
    assert good is not topology
    assert good.chain_types.tolist() == [0] * 4 + list(range(1, 64))
    # a keyfile from another run that splits the FIRST label 2 + 2 is not this run
    wrong = topology.with_keyfile_types(
        [[2, "AAAA"], [2, "AAAA"]] + singles, labelled=True
    )
    assert wrong is topology


def test_exhausted_identifiers_are_resolved_by_sequence_when_lines_are_out_of_order():
    sequences, pdb_types = _sixty_two_labels()
    # a restart run with an EXTRA_CHAIN: one more "AAAA" chain appended at the end,
    # which PIMMS files under the existing first type (label 0)
    sequences = sequences + ["AAAA"]
    pdb_types = pdb_types + [0]
    topology = Topology(sequences, chain_types=pdb_types)
    # keyfile_used.kf lists one line per type, so "5 AAAA" does not expand in order
    specs = [[5, "AAAA"]] + [[1, s] for s in sequences[4:67]]
    typed = topology.with_keyfile_types(specs, labelled=True)
    assert typed is not topology
    assert typed.chain_types.tolist() == [0] * 4 + list(range(1, 64)) + [0]
    assert topology.matches_keyfile_composition(specs)
    assert not topology.matches_keyfile_composition([[6, "AAAA"]] + specs[1:])


def test_restart_snapshot_types_split_the_shared_identifier(tmp_path, monkeypatch):
    sequences, pdb_types = _sixty_two_labels()
    # PIMMS's own chain types for the snapshot's chains: 63 of them, numbered from 5
    engine_types = [5] * 4 + list(range(6, 66)) + [66, 67, 68]
    snapshot = types.SimpleNamespace(
        chains={
            i + 1: [None, s, t] for i, (s, t) in enumerate(zip(sequences, engine_types))
        }
    )
    seqs, restart_types, n_snapshot = _load._restart_chain_types(
        snapshot, [[2, "AAAA"], [1, "GA"]]
    )
    assert n_snapshot == 67 and seqs == sequences + ["AAAA", "AAAA", "GA"]
    # an existing sequence joins its type; a new one takes the next unused type
    assert restart_types == engine_types + [5, 5, 69]

    dims = (30, 10, 10)
    positions = np.zeros((2, sum(len(s) for s in sequences), 3), dtype=np.int32)
    bead = 0
    for c, s in enumerate(sequences):
        positions[:, bead : bead + len(s), 0] = np.arange(len(s))
        positions[:, bead : bead + len(s), 1] = c % 10
        positions[:, bead : bead + len(s), 2] = c // 10
        bead += len(s)
    xtc, pdb = _write(
        tmp_path, positions, sequences, dims, labels=[_IDS[t] for t in pdb_types]
    )
    keyfile = os.path.join(str(tmp_path), "KEYFILE.kf")
    with open(keyfile, "w") as fh:
        fh.write(
            "DIMENSIONS : 30 10 10\nTEMPERATURE : 50\nRESTART_FILE : start.pimms\n"
        )
    monkeypatch.setattr(
        _load, "_restart_snapshot", lambda keydict, kf, keyword: snapshot
    )
    traj = _quiet_load(xtc, pdb, keyfile)
    assert traj.chain_types.tolist() == [0] * 4 + list(range(1, 64))

    # a restart keyfile from a different run used to load without a word
    other = types.SimpleNamespace(chains={1: [None, "AAAA", 0], 2: [None, "GG", 1]})
    monkeypatch.setattr(_load, "_restart_snapshot", lambda keydict, kf, keyword: other)
    # (two warnings, both wanted: the mismatch itself, and the consequence - the
    # merged chain types could not be recovered, with what to pass instead)
    with pytest.warns(UserWarning) as caught:
        traj = lemonade.load(xtc, pdb, keyfile)
    messages = [str(w.message) for w in caught]
    assert len(messages) == 2
    assert any(
        re.search(
            r"RESTART_FILE : start\.pimms\) holds 2 chains.*do not match the 67 chains",
            m,
            re.S,
        )
        for m in messages
    )
    assert any(
        "could not be read back" in m and "keyfile_used.kf" in m for m in messages
    )
    assert len(set(traj.chain_types.tolist())) == 62


# -- A2-13: small gaps in load() -----------------------------------------------


def test_non_ascii_xtc_path_loads(tmp_path, monkeypatch):
    # the temporary ASCII-named link is made in the system temporary directory:
    # point that inside tmp_path so the test writes nowhere else
    links = tmp_path / "links"
    links.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(links))
    rng = np.random.default_rng(7)
    dims = (12, 6, 5)
    truth = _straight_chains(rng, 4, 6, 5, dims)
    # mdtraj cannot write to such a path either, so write beside it and copy in
    plain_xtc, plain_pdb = _write(tmp_path / "plain", truth, ["AAAAA"] * 6, dims)
    directory = tmp_path / "dir_\u00e9\u4e2d"
    directory.mkdir()
    xtc = shutil.copy(plain_xtc, str(directory / "traj.xtc"))
    pdb = shutil.copy(plain_pdb, str(directory / "START.pdb"))
    assert np.array_equal(_quiet_load(xtc, pdb).positions, truth)
    assert np.array_equal(_quiet_load(xtc, pdb, start=1, step=2).positions, truth[1::2])
    # nothing is left behind, next to the trajectory or where the link was made
    assert sorted(os.listdir(str(directory))) == ["START.pdb", "traj.xtc"]
    assert os.listdir(str(links)) == []


def test_bool_spacing_is_refused(tmp_path):
    truth = _straight_chains(np.random.default_rng(8), 2, 2, 3, (6, 6, 6))
    xtc, pdb = _write(tmp_path, truth, ["AAA"] * 2, (6, 6, 6))
    for bad in (True, False, np.True_):
        with pytest.raises(
            ValueError, match="spacing must be a finite positive number"
        ):
            lemonade.load(xtc, pdb, spacing=bad)


@pytest.mark.parametrize("dims", [(8, 8, 8), (12, 12)])
def test_pdb_only_load_below_one_angstrom_reads_the_box_from_cryst1(tmp_path, dims):
    truth = _straight_chains(np.random.default_rng(9), 1, 4, 4, dims)
    _xtc, pdb = _write(tmp_path, truth, ["AAAA"] * 4, dims, spacing=0.1, xtc=False)
    traj = _quiet_load(pdb=pdb, spacing=0.1)
    assert traj.dimensions == dims
    assert np.array_equal(traj.positions, truth)
    # at the default spacing that box is smaller than one site: say so, and how to
    # fix it. The coordinates do not fit that lattice either, which is warned
    # about before the box is refused; nothing else may be warned (in particular
    # not mdtraj's notice that it discarded the box, since load() reads it itself).
    with pytest.warns(UserWarning) as caught:
        with pytest.raises(
            ValueError, match="smaller than one lattice site at spacing 3.65"
        ):
            lemonade.load(pdb=pdb)
    messages = [str(w.message) for w in caught]
    assert len(messages) == 1 and "round-off residual" in messages[0]


def test_xtc_without_a_box_takes_it_from_the_pdb_cryst1(tmp_path):
    dims = (8, 8, 8)
    truth = _straight_chains(np.random.default_rng(13), 3, 4, 4, dims)
    xtc, pdb = _write(tmp_path, truth, ["AAAA"] * 4, dims, spacing=0.5, xtc_box=False)
    with md.formats.XTCTrajectoryFile(xtc) as fh:
        assert not np.any(fh.read()[3])  # the file really records no box
    traj = _quiet_load(xtc, pdb, spacing=0.5)
    assert traj.dimensions == dims
    assert np.array_equal(traj.positions, truth)


def test_restart_types_are_recovered_from_a_real_restart_file(tmp_path, monkeypatch):
    from pimms.restart import RestartObject

    sequences, pdb_types = _sixty_two_labels()
    engine_types = [5] * 4 + list(range(6, 66)) + [66, 67, 68]
    dims = (
        70,
        10,
        10,
    )  # the longest chain is 62 beads; the restart reader checks the box
    positions = np.zeros((2, sum(len(s) for s in sequences), 3), dtype=np.int32)
    bead = 0
    chains = {}
    for c, (seq, chain_type) in enumerate(zip(sequences, engine_types)):
        block = positions[:, bead : bead + len(seq)]
        block[..., 0] = np.arange(len(seq))
        block[..., 1] = c % 10
        block[..., 2] = c // 10
        # the layout PIMMS writes: positions as lists of Python ints, sequence, chainType
        chains[c + 1] = [[[int(v) for v in p] for p in block[0]], seq, chain_type]
        bead += len(seq)
    # written by the tree's own restart code, into tmp_path (it writes to the cwd)
    monkeypatch.chdir(tmp_path)
    snapshot = RestartObject()
    snapshot.dimensions = list(dims)
    snapshot.hardwall = False
    snapshot.chains = chains
    snapshot.write_to_file()
    assert os.path.isfile(CONFIG.RESTART_FILENAME)

    # the run restarted from it with one more chain of an existing sequence
    # (it joins that sequence's type) and one of a new sequence (a new type)
    extra = np.zeros((2, 6, 3), dtype=np.int32)
    extra[:, :4, 0] = np.arange(4)
    extra[:, :4, 1] = 8
    extra[:, 4:, 0] = np.arange(2)
    extra[:, 4:, 1] = 9
    extra[..., 2] = 9
    all_positions = np.concatenate([positions, extra], axis=1)
    xtc, pdb = _write(
        tmp_path / "run",
        all_positions,
        sequences + ["AAAA", "GA"],
        dims,
        labels=[_IDS[t] for t in pdb_types] + [_IDS[0], _IDS[61]],
    )
    keyfile = os.path.join(str(tmp_path), "KEYFILE.kf")
    with open(keyfile, "w") as fh:
        fh.write(
            "DIMENSIONS : 70 10 10\nTEMPERATURE : 50\n"
            f"RESTART_FILE : {CONFIG.RESTART_FILENAME}\n"
            "EXTRA_CHAIN : 1 AAAA\nEXTRA_CHAIN : 1 GA\n"
        )
    traj = _quiet_load(xtc, pdb, keyfile)
    # 63 types in the snapshot, the extra AAAA chain in the first, GA a 64th
    assert traj.chain_types.tolist() == [0] * 4 + list(range(1, 64)) + [0, 64]


# -- A2-11: equal-size clusters ------------------------------------------------


@pytest.mark.parametrize("shift", [0, 3, 7])
def test_equal_size_clusters_are_ordered_by_lowest_chain_index(shift):
    dims = (12, 12, 12)
    # chains 0 and 3 each alone, chains 1 and 2 touching: sizes 3, 6, 3
    base = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [0, 5, 5],
            [1, 5, 5],
            [2, 5, 5],
            [0, 6, 5],
            [1, 6, 5],
            [2, 6, 5],
            [6, 9, 9],
            [7, 9, 9],
            [8, 9, 9],
        ],
        dtype=np.int32,
    )
    frame = ((base + shift) % 12)[None]
    store = TrajectoryStore(frame, dims, 3.65, False, Topology(["AAA"] * 4))
    assert [sorted(c) for c in store.cluster_membership(0)] == [[1, 2], [0], [3]]
    mirrored = ((11 - (base + shift)) % 12)[None]
    store = TrajectoryStore(mirrored, dims, 3.65, False, Topology(["AAA"] * 4))
    assert [sorted(c) for c in store.cluster_membership(0)] == [[1, 2], [0], [3]]


# -- D9: the topology and the store are on the reference page ------------------


def test_reference_page_documents_topology_and_store():
    here = os.path.dirname(os.path.abspath(__file__))
    page = os.path.join(here, "..", "..", "..", "docs", "lemonade", "reference.rst")
    if not os.path.isfile(page):
        pytest.skip("documentation sources are not installed alongside the package")
    with open(page) as fh:
        text = fh.read()
    for name in (
        "pimms.lemonade._topology.Topology",
        "pimms.lemonade._store.TrajectoryStore",
    ):
        assert f".. autoclass:: {name}" in text
