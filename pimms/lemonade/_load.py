## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Loading PIMMS trajectories into lemonade.

``load()`` accepts any of an XTC (coordinates), a PDB (topology, and/or a single
frame) and a PIMMS keyfile (box size, lattice spacing, hardwall, chain types), in
the combinations a user actually has to hand:

* ``xtc`` + ``pdb``            - the usual case (full trajectory).
* ``xtc`` + ``pdb`` + ``keyfile`` - adds authoritative spacing / dimensions /
  hardwall / chain types.
* ``pdb`` only                 - a single frame (e.g. START.pdb).

Coordinates are converted back to the integer lattice with vectorised arithmetic
(``round(nm / (spacing/10))``), a bounded block of frames at a time, and only for
the frames the ``start``/``stop``/``step``/``n_frames`` selection keeps; the
topology comes from the PDB (which matches the XTC bead order exactly) and is
refined with keyfile chain types when available.
"""

import contextlib
import math
import numbers
import os
import tempfile
import warnings
from typing import Any, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

from ._topology import Topology, pdb_chain_labels, _N_PDB_CHAIN_IDS
from ._store import TrajectoryStore
from .trajectory import LatticeTrajectory

DEFAULT_SPACING = 3.65   # PIMMS LATTICE_TO_ANGSTROMS default (v0.1.34+)

# Upper bound on the float32 coordinates decoded from the XTC in one read. The
# conversion to the integer lattice needs about three arrays of this size at
# once, plus the block being filed, so the transient memory of a load is about
# five times this (20 MB) however long the trajectory is - or about five frames,
# for a system so large (above about 350,000 beads) that one frame is bigger
# than the block. What is retained is the int32 lattice of the kept frames.
_CHUNK_BYTES = 4 * 1024 * 1024

# What one separate read call costs, expressed as the number of beads that take
# as long to decode (about 0.4 ms against about 80 ns per bead). Kept frames
# closer together than this are read as one span and picked out of it.
_BEADS_PER_READ = 5000

# PIMMS hardcodes these two names for the compact phase of a RESIZED_EQUILIBRATION
# run (simulation.py:216-217); nothing else in PIMMS writes them, so seeing one is
# sufficient evidence that the frames were produced in the resized box under the
# forced hardwall of that phase - including when no keyfile is passed at all.
_EQUILIBRATION_FILENAMES = ("eq_traj.xtc", "eq_START.pdb")

# The start of the warning mdtraj raises when it reads a PDB whose CRYST1 box it
# takes for a placeholder and discards (every PIMMS box at a sub-angstrom lattice
# spacing). load() reads that record itself in that case, so the warning would
# tell the user the box was lost when it was not.
_MDTRAJ_DROPPED_CELL_WARNING = "Unlikely unit cell vectors detected in PDB file"


def _is_equilibration_file(path):
    """Is ``path`` one of PIMMS's resized-equilibration filenames?

    Parameters
    ----------
    path : str or None
        A trajectory or topology path, or ``None``.

    Returns
    -------
    bool
        ``True`` if the basename is ``eq_traj.xtc`` or ``eq_START.pdb``, i.e. the
        file belongs to the compact equilibration phase of a
        ``RESIZED_EQUILIBRATION`` run.
    """
    if path is None:
        return False
    return os.path.basename(str(path)) in _EQUILIBRATION_FILENAMES


@contextlib.contextmanager
def _without_dropped_cell_warning() -> Iterator[None]:
    """Read a PDB through mdtraj without its "dummy CRYST1" warning.

    mdtraj discards the unit cell of a box it judges implausibly small for the
    number of atoms, with a warning, and returns no box. ``load()`` then takes
    the box from the PDB's ``CRYST1`` record directly (``_pdb_cryst1_nm``), so
    the warning describes a loss that does not happen. Only that one message is
    silenced; anything else mdtraj has to say still reaches the user.

    Returns
    -------
    Iterator[None]
        Context manager; nothing is yielded.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=_MDTRAJ_DROPPED_CELL_WARNING)
        yield


def _restart_snapshot(keydict, keyfile, keyword):
    """Read the restart file a keyfile points at and return the parsed snapshot.

    ``RESTART_OVERRIDE_HARDWALL`` and ``RESTART_OVERRIDE_DIMENSIONS`` tell PIMMS to
    discard the keyfile's ``HARDWALL`` / ``DIMENSIONS`` and adopt the snapshot's
    instead. That substitution happens inside restart processing, which only a full
    keyfile parse runs, so lemonade (which parses with ``parse_only=True``) has to
    do it itself by opening the restart pickle.

    Parameters
    ----------
    keydict : dict
        The ``parse_only=True`` keyword lookup, which must carry ``RESTART_FILE``.
    keyfile : str
        Path to the keyfile. A relative ``RESTART_FILE`` is resolved against the
        keyfile's own directory first, because that is where a run directory keeps
        it; the current working directory is only a fallback.
    keyword : str
        The override keyword that made this read necessary, quoted in any error.

    Returns
    -------
    pimms.restart.RestartObject
        The snapshot, from which ``.hardwall`` and ``.dimensions`` are read.

    Raises
    ------
    ValueError
        If the restart file cannot be found, or cannot be read as a restart file.
        Guessing is not an option here: the whole point of the override keyword is
        that the keyfile value is the wrong one.
    """
    raw = str(keydict["RESTART_FILE"])
    candidates = []
    if os.path.isabs(raw):
        candidates.append(raw)
    else:
        candidates.append(os.path.join(os.path.dirname(os.path.abspath(keyfile)), raw))
        candidates.append(os.path.abspath(raw))

    for path in candidates:
        if os.path.isfile(path):
            from pimms.restart import RestartObject, RestartException
            snapshot = RestartObject()
            try:
                snapshot.build_from_file(path)
            except RestartException as e:
                raise ValueError(
                    f"lemonade.load: the keyfile sets {keyword}, so the effective "
                    f"boundary condition/box must come from the restart file "
                    f"'{path}', but that file could not be read ({e}). Pass "
                    f"hardwall= / dimensions= explicitly.")
            return snapshot

    raise ValueError(
        f"lemonade.load: the keyfile sets {keyword} : True, which means PIMMS took "
        f"the run's boundary condition/box from the restart file rather than from "
        f"the keyfile - so the keyfile values are NOT the ones the run used. The "
        f"restart file 'RESTART_FILE : {raw}' could not be found (looked in "
        f"{', '.join(repr(c) for c in candidates)}). Put it beside the keyfile, or "
        f"pass hardwall= / dimensions= explicitly.")


def _select_frames(n_total: int, start: Optional[int], stop: Optional[int],
                   step: Optional[int], n_frames: Optional[int]) -> np.ndarray:
    """Indices of the frames a load keeps, in the order they are returned.

    This is the selection ``load()`` has always applied, written against the
    frame count alone so it can be made before any coordinate is read:
    ``start``/``stop``/``step`` slice the frames exactly as a numpy slice would
    (negative values count from the end, a negative ``step`` reverses), and
    ``n_frames`` then thins what survives to evenly spaced frames.

    Parameters
    ----------
    n_total : int
        Number of frames available in the file.
    start, stop, step : int or None
        The frame slice. Anything a numpy slice refuses (a float, a zero step) is
        refused here with the same exception.
    n_frames : int or None
        If given and smaller than the number of sliced frames, the evenly spaced
        subsample to keep.

    Returns
    -------
    numpy.ndarray
        1D integer array of file frame indices, in output order.

    Raises
    ------
    ValueError
        If the selection keeps no frames.
    """
    selected = np.arange(n_total)
    if start is not None or stop is not None or step is not None:
        selected = selected[slice(start, stop, step)]
    if n_frames is not None and n_frames < len(selected):
        idx = np.linspace(0, len(selected) - 1, int(n_frames)).round().astype(int)
        selected = selected[idx]
    if len(selected) == 0:
        raise ValueError(
            f"frame selection start={start}, stop={stop}, step={step} keeps no frames "
            f"of a {n_total}-frame trajectory")
    return selected


def _nm_to_lattice(xyz: np.ndarray, spacing: float) -> Tuple[np.ndarray, float]:
    """Convert coordinates in nanometres to integer lattice coordinates.

    Parameters
    ----------
    xyz : numpy.ndarray
        ``(n_frames, n_beads, 3)`` float32 coordinates in nanometres, as mdtraj
        returns them.
    spacing : float
        Lattice spacing in angstroms.

    Returns
    -------
    lattice : numpy.ndarray
        ``(n_frames, n_beads, 3)`` int32 lattice coordinates, not yet wrapped
        into the box.
    residual : float
        The largest distance of any coordinate from its lattice site, in lattice
        units (``0.0`` for an empty array). The subtraction is done in the
        coordinates' own float32: a float and its nearest integer differ by an
        exactly representable amount, so this is the same number the float64
        subtraction gives, without the float64 copy.

    Raises
    ------
    ValueError
        If a coordinate is not finite or does not fit in int32.
    """
    lattice_f = xyz / (spacing * 0.1)
    rounded = np.rint(lattice_f)
    if not rounded.size:
        return rounded.astype(np.int32), 0.0
    # a nan or an inf anywhere shows up in the extremes, so no mask is needed
    int32_info = np.iinfo(np.int32)
    lowest, highest = float(rounded.min()), float(rounded.max())
    if (not (math.isfinite(lowest) and math.isfinite(highest)) or
            lowest < int32_info.min or highest > int32_info.max):
        raise ValueError(
            "trajectory coordinates cannot be represented as int32 lattice coordinates")
    np.subtract(lattice_f, rounded, out=lattice_f)
    np.abs(lattice_f, out=lattice_f)
    residual = float(lattice_f.max())
    del lattice_f
    return rounded.astype(np.int32), residual


@contextlib.contextmanager
def _open_xtc(path: Union[str, os.PathLike]) -> Iterator[Any]:
    """Open an XTC for reading, whatever characters its path contains.

    mdtraj hands the path to the xdr library as ASCII bytes, so a path with any
    other character in it (an accented directory name, say) used to stop the load
    with a bare ``UnicodeEncodeError``. Such a file is opened through a
    temporary ASCII-named symbolic link instead, made in a fresh directory under
    the system temporary directory (``tempfile.gettempdir()``); the link and
    its directory are removed again when the file is closed.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the XTC file.

    Yields
    ------
    mdtraj.formats.XTCTrajectoryFile
        The open file, closed on exit.

    Raises
    ------
    ValueError
        If the path is not ASCII and no ASCII-named link to it can be made (no
        symbolic links on this platform, or a temporary directory that is not
        ASCII either). The message names the file and says what to do.
    """
    from mdtraj.formats import XTCTrajectoryFile

    path = os.fspath(path)
    link_dir = None
    link = None
    try:
        try:
            handle = XTCTrajectoryFile(path)
        except UnicodeEncodeError:
            if not os.path.isfile(path):
                raise OSError(f"The file '{path}' doesn't exist")
            try:
                link_dir = tempfile.mkdtemp(prefix="lemonade_")
                link = os.path.join(link_dir, "traj.xtc")
                os.symlink(os.path.abspath(path), link)
                handle = XTCTrajectoryFile(link)
            except (OSError, NotImplementedError, UnicodeEncodeError) as e:
                raise ValueError(
                    f"lemonade.load: the XTC path '{path}' contains non-ASCII "
                    f"characters, which the XTC reader cannot open, and a temporary "
                    f"ASCII-named link to it could not be made ({e}). Load it by a "
                    f"relative path from inside its own directory, or copy/link the "
                    f"file to a path made of ASCII characters only.")
        try:
            yield handle
        finally:
            handle.close()
    finally:
        if link is not None and os.path.islink(link):
            os.remove(link)
        if link_dir is not None and os.path.isdir(link_dir):
            os.rmdir(link_dir)


def _count_complete_frames(handle: Any) -> int:
    """Count the frames at the start of an XTC that can be read in full.

    The frames are read one at a time from wherever the file is positioned,
    which must be its start: the count never touches the reader's frame index
    or seeks, because for a torn file of uncompressed frames (nine beads or
    fewer) building that index is exactly what fails.

    Parameters
    ----------
    handle : mdtraj.formats.XTCTrajectoryFile
        A freshly opened XTC file. Its position is moved.

    Returns
    -------
    int
        The number of leading frames that decode without an error; the count
        stops at the end of the file or at the first frame the reader refuses.
    """
    n_complete = 0
    while True:
        try:
            if len(handle.read(n_frames=1)[0]) == 0:
                break
        except (RuntimeError, OSError, AssertionError):
            break
        n_complete += 1
    return n_complete


def _probe_xtc(xtc: Union[str, os.PathLike]) -> Tuple[int, bool]:
    """Find how many complete frames an XTC holds, and whether it can be seeked.

    A file that ends cleanly is answered from the reader's frame index, with one
    read at the tail to confirm it, and costs nothing more. A file that ends in
    an incomplete frame - the trajectory of a run that is still going, or of one
    that was killed - is counted frame by frame from the start, and a warning
    gives the number of complete frames found.

    Parameters
    ----------
    xtc : str or os.PathLike
        Path to the XTC file.

    Returns
    -------
    n_frames : int
        The number of complete frames at the start of the file.
    seekable : bool
        ``False`` if the reader cannot build its frame index for this file, so
        the frames have to be read in order from the first. That is the case for
        a torn file of uncompressed frames (systems of nine beads or fewer),
        where the reader works the frame count out from the file size and
        refuses a size that is not a whole number of frames.

    Raises
    ------
    RuntimeError or OSError
        The reader's own error if the file holds no complete frame at all.
    """
    n_listed = None
    seekable = True
    with _open_xtc(xtc) as handle:
        # The frame index mdtraj builds lists every frame whose header is on disk,
        # so a torn last frame is either listed (its body is cut short) or is a
        # stub of a header after the last listed one. Asking for two frames from
        # the last listed one finds both: a whole file gives back exactly one.
        try:
            n_listed = len(handle)
            handle.seek(max(n_listed - 1, 0))
            if len(handle.read(n_frames=2)[0]) == 1:
                return n_listed, True
            error = None
        except AssertionError as e:
            # uncompressed frames: the index is the file size over the frame size
            seekable = False
            error = e
        except (RuntimeError, OSError) as e:
            error = e

    with _open_xtc(xtc) as handle:
        n_complete = _count_complete_frames(handle)
    if n_complete == 0:
        if isinstance(error, (RuntimeError, OSError)):
            raise error
        raise RuntimeError(f"lemonade.load: no complete frame could be read from "
                           f"'{os.fspath(xtc)}'") from error
    warnings.warn(
        f"lemonade.load: '{os.fspath(xtc)}' ends in an incomplete frame (the "
        f"run is still writing it, or was killed part-way through a write). "
        f"Recovered the {n_complete} complete frame(s) before it"
        + (f"; the file's own index lists {n_listed}, so frames after the "
           f"first unreadable one have been left out too"
           if n_listed is not None and n_listed > n_complete + 1 else "")
        + ". Frame selection (start/stop/step/n_frames) applies to those "
          "complete frames.", stacklevel=4)
    return n_complete, seekable


def _read_plan(in_file_order: np.ndarray, chunk: int, max_gap: int, seekable: bool
               ) -> Iterator[Tuple[int, int, int, int, int]]:
    """Plan the reads that fetch a sorted set of frames from an XTC.

    Each read is ``n_read`` frames taken every ``stride`` frames from file frame
    ``first``; the frames wanted from it are those of
    ``in_file_order[lo:hi]``. Three kinds of read are planned:

    * a run of evenly spaced frames is one strided read, which skips the frames
      in between without decoding them;
    * frames that are close together but unevenly spaced (an ``n_frames``
      subsample that keeps most of the file) are read as one contiguous span
      and picked out of it, because a separate read for every frame or two
      costs more than decoding the few frames in between;
    * a file that cannot be seeked is read in contiguous spans from its start.

    No read spans more than ``chunk`` frames.

    Parameters
    ----------
    in_file_order : numpy.ndarray
        The file frame indices wanted, ascending (repeats allowed).
    chunk : int
        The largest number of frames one read may decode.
    max_gap : int
        The largest spacing between two wanted frames that is still read
        through rather than seeked over.
    seekable : bool
        Whether the reader can be positioned at an arbitrary frame.

    Yields
    ------
    tuple of int
        ``(first, n_read, stride, lo, hi)``. ``lo == hi`` marks a span of an
        unseekable file that holds no wanted frame and is read only to get past
        it.
    """
    n_kept = len(in_file_order)
    if not seekable:
        position = 0
        file_position = 0
        last = int(in_file_order[-1])
        while position < n_kept:
            n_read = min(chunk, last + 1 - file_position)
            end = int(np.searchsorted(in_file_order, file_position + n_read, side="left"))
            yield file_position, n_read, 1, position, end
            position = end
            file_position += n_read
        return

    gaps = np.diff(in_file_order)
    position = 0
    while position < n_kept:
        first = int(in_file_order[position])
        stride = int(gaps[position]) if position < n_kept - 1 else 1
        run_end = position + 1
        if stride >= 1:
            limit = min(n_kept, position + chunk)
            while run_end < limit and gaps[run_end - 1] == stride:
                run_end += 1
        if stride > max_gap or (stride > 1 and run_end - position >= 8):
            yield first, run_end - position, stride, position, run_end
            position = run_end
            continue
        end = position + 1
        while (end < n_kept and gaps[end - 1] <= max_gap
               and in_file_order[end] - first < chunk):
            end += 1
        yield first, int(in_file_order[end - 1]) - first + 1, 1, position, end
        position = end


def _read_xtc_lattice(xtc: Union[str, os.PathLike], md_topology: Any, spacing: float,
                      start: Optional[int], stop: Optional[int], step: Optional[int],
                      n_frames: Optional[int]
                      ) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], float, np.ndarray]:
    """Read the selected frames of an XTC straight onto the integer lattice.

    Only the frames the selection keeps are held: the file is read in blocks of
    at most ``_CHUNK_BYTES`` of float32 coordinates (see :func:`_read_plan`), and
    each block is converted to int32 before the next is read, so the memory in
    use never much exceeds the lattice that is returned.

    A file that ends in an incomplete frame - the trajectory of a run that is
    still going, or of one that was killed - is loaded up to its last complete
    frame, with a warning that says how many frames that is (see
    :func:`_probe_xtc`). The frame selection is then made against the complete
    frames only.

    Parameters
    ----------
    xtc : str or os.PathLike
        Path to the XTC file.
    md_topology : mdtraj.Topology
        Topology read from the PDB; its bead count must match the file's.
    spacing : float
        Lattice spacing in angstroms.
    start, stop, step : int or None
        Frame slice, as for :func:`load`.
    n_frames : int or None
        Even subsample, as for :func:`load`.

    Returns
    -------
    lattice : numpy.ndarray
        ``(n_kept, n_beads, 3)`` int32 lattice coordinates, not yet wrapped.
    times : numpy.ndarray
        ``(n_kept,)`` float64 frame times.
    box_nm : numpy.ndarray or None
        Box lengths in nanometres recorded with the first kept frame, or
        ``None`` if the file records no box.
    residual : float
        Largest lattice round-off over the kept frames.
    first_xyz : numpy.ndarray
        ``(n_beads, 3)`` float32 coordinates (nm) of the first kept frame.

    Raises
    ------
    ValueError
        If the selection keeps no frames, if the coordinates do not fit in int32
        lattice coordinates, or (from mdtraj) if the PDB and XTC bead counts
        differ.
    RuntimeError or OSError
        Whatever mdtraj raises for a file with no complete frame at all, or for
        a frame that is damaged before the end of the file.
    """
    n_beads = md_topology.n_atoms
    chunk = max(1, _CHUNK_BYTES // (12 * max(1, n_beads)))
    # A separate read costs about as much as decoding _BEADS_PER_READ beads, so
    # frames closer together than this are cheaper to read through than to seek to.
    max_gap = min(chunk, 1 + _BEADS_PER_READ // max(1, n_beads))

    n_total, seekable = _probe_xtc(xtc)
    selected = _select_frames(n_total, start, stop, step, n_frames)
    n_kept = len(selected)
    lattice = np.empty((n_kept, n_beads, 3), dtype=np.int32)
    times = np.empty(n_kept, dtype=np.float64)
    box_nm = None
    first_xyz = None
    residual = 0.0

    # Read in file order (one forward pass), whatever order the selection
    # returns the frames in, and put each block where the selection wants it.
    order = np.argsort(selected, kind="stable")
    in_file_order = selected[order]
    with _open_xtc(xtc) as handle:
        for first, n_read, stride, lo, hi in _read_plan(in_file_order, chunk, max_gap,
                                                        seekable):
            if seekable:
                handle.seek(first)
            block = handle.read_as_traj(md_topology, n_frames=n_read, stride=stride)
            if block.n_frames != n_read:
                raise RuntimeError(
                    f"lemonade.load: expected {n_read} frame(s) from frame "
                    f"{first} of '{os.fspath(xtc)}' but read {block.n_frames}")
            if lo == hi:
                continue
            # which frames of the block are wanted (all of them, for a strided read)
            pick = (in_file_order[lo:hi] - first) // stride
            whole = len(pick) == n_read and bool(np.all(pick == np.arange(n_read)))
            xyz = np.asarray(block.xyz)
            targets = order[lo:hi]
            block_lattice, block_residual = _nm_to_lattice(xyz if whole else xyz[pick],
                                                           spacing)
            lattice[targets] = block_lattice
            times[targets] = block.time if whole else block.time[pick]
            residual = max(residual, block_residual)
            is_first = np.flatnonzero(targets == 0)
            if is_first.size:
                in_block = int(pick[is_first[0]])
                first_xyz = np.array(xyz[in_block], dtype=np.float32)
                if block.unitcell_lengths is not None:
                    box_nm = np.array(block.unitcell_lengths[in_block])
    return lattice, times, box_nm, residual, first_xyz


def _pdb_cryst1_nm(pdb: Union[str, os.PathLike]) -> Optional[np.ndarray]:
    """Read the box lengths from a PDB's own ``CRYST1`` record.

    mdtraj discards a ``CRYST1`` record whose lengths look like a placeholder,
    which is every PIMMS box below about an angstrom of lattice spacing, so a
    PDB-only load of such a system had no box to infer the dimensions from.

    Parameters
    ----------
    pdb : str or os.PathLike
        Path to the PDB file.

    Returns
    -------
    numpy.ndarray or None
        The three box lengths in nanometres, or ``None`` if the file has no
        readable ``CRYST1`` record.
    """
    try:
        with open(pdb, "r") as fh:
            for line in fh:
                if line.startswith("CRYST1"):
                    lengths = [float(line[6:15]), float(line[15:24]), float(line[24:33])]
                    return np.array(lengths, dtype=np.float64) * 0.1
                if line.startswith(("ATOM", "HETATM")):
                    break
    except (OSError, ValueError, UnicodeDecodeError):
        return None
    return None


def _check_unit_bonds(frame_lattice: np.ndarray, offsets: np.ndarray,
                      dimensions: Tuple[int, ...], spacing: float,
                      frame_xyz: Optional[np.ndarray], residual: float) -> bool:
    """Warn if the bonds of a frame are not single lattice steps.

    Two bonded PIMMS beads always sit on neighbouring lattice sites (one step,
    diagonals included), so every bond of a correctly loaded frame has Chebyshev
    length exactly 1 under the minimum-image convention. A spacing that is a
    whole multiple of the one assumed leaves no round-off residual at all - the
    coordinates land exactly on every second (third, ...) site - and this is the
    check that sees it. The same test catches a PDB whose chains are not this
    trajectory's chains.

    A system with no bonds (every chain a single bead) has nothing to test, and
    is passed without a warning.

    Parameters
    ----------
    frame_lattice : numpy.ndarray
        ``(n_beads, 3)`` integer lattice coordinates of one frame (wrapped or
        not; the minimum image is taken here).
    offsets : numpy.ndarray
        ``(n_chains + 1,)`` bead-index boundaries of the chains.
    dimensions : tuple of int
        Box extent in lattice units (2 or 3 entries).
    spacing : float
        The lattice spacing, in angstroms, the frame was converted with.
    frame_xyz : numpy.ndarray or None
        ``(n_beads, 3)`` coordinates of the same frame in nanometres, used to
        estimate the true spacing from the shortest bond component; ``None``
        skips the estimate.
    residual : float
        The lattice round-off residual of the load. When it is small the
        estimate is snapped to a whole multiple of ``spacing``.

    Returns
    -------
    bool
        ``True`` if a warning was raised, ``False`` if every bond is a unit step
        or there are no bonds.
    """
    n_beads = frame_lattice.shape[0]
    if n_beads < 2:
        return False
    bonded = np.ones(n_beads - 1, dtype=bool)
    bonded[np.asarray(offsets[1:-1], dtype=np.int64) - 1] = False
    if not bonded.any():
        return False
    n_dim = len(dimensions)
    dims = np.array(dimensions, dtype=np.int64)
    steps = np.diff(frame_lattice[:, :n_dim].astype(np.int64), axis=0)[bonded]
    half = dims // 2
    steps = (steps + half) % dims - half
    lengths = np.abs(steps).max(axis=1)
    bad = lengths != 1
    if not bad.any():
        return False

    advice = ("the PDB topology or the box dimensions probably do not belong to "
              "this trajectory")
    if frame_xyz is not None:
        real = np.abs(np.diff(frame_xyz[:, :n_dim].astype(np.float64), axis=0)[bonded])
        real = real[real > 0]
        if real.size:
            estimate = 10.0 * float(real.min())
            ratio = estimate / spacing
            if residual <= 0.05 and round(ratio) >= 1 and abs(ratio - round(ratio)) < 0.05:
                estimate = spacing * round(ratio)
            if abs(estimate / spacing - 1.0) >= 0.05:
                advice = (f"from the shortest bond the spacing is probably "
                          f"{estimate:.4g} A, so pass spacing={estimate:.4g} or the "
                          f"keyfile (LATTICE_TO_ANGSTROMS)")
    warnings.warn(
        f"lemonade.load: {int(bad.sum())} of {int(bad.size)} bonds in the first loaded "
        f"frame are not single lattice steps at spacing {spacing} A (bond lengths "
        f"{int(lengths.min())} to {int(lengths.max())} lattice units, and PIMMS bonds "
        f"are always exactly 1). Positions, box and every contact, cluster and "
        f"density calculation are therefore on the wrong lattice: {advice}.",
        stacklevel=3)
    return True


def _restart_chain_types(snapshot: Any, extra_specs: Sequence[Sequence[Any]]
                         ) -> Tuple[List[str], List[int], int]:
    """Per-chain sequences and chain types of a run started from a restart file.

    PIMMS builds such a run from the snapshot's chains, in chain order, and then
    appends the ``EXTRA_CHAIN`` chains; an extra chain whose sequence the
    snapshot already holds joins that sequence's existing chain type, and a new
    sequence gets the next unused type.

    Parameters
    ----------
    snapshot : pimms.restart.RestartObject
        The parsed restart file. Its ``chains`` maps a chain ID onto
        ``[positions, sequence, chainType]``.
    extra_specs : list of [int, str]
        The keyfile ``EXTRA_CHAIN`` lines as ``[count, sequence]`` pairs (empty
        for none).

    Returns
    -------
    sequences : list of str
        One sequence per chain: the snapshot's chains, then the extra chains.
    chain_types : list of int
        PIMMS's chain type for each of those chains.
    n_snapshot : int
        How many of the chains came from the snapshot itself.
    """
    sequences = []
    chain_types = []
    first_type = {}
    for chain_id in sorted(snapshot.chains):
        sequence = str(snapshot.chains[chain_id][1])
        chain_type = int(snapshot.chains[chain_id][2])
        sequences.append(sequence)
        chain_types.append(chain_type)
        first_type.setdefault(sequence, chain_type)
    n_snapshot = len(sequences)
    for count, sequence in extra_specs:
        sequence = str(sequence)
        if sequence not in first_type:
            first_type[sequence] = max(chain_types) + 1 if chain_types else 0
        for _ in range(int(count)):
            sequences.append(sequence)
            chain_types.append(first_type[sequence])
    return sequences, chain_types, n_snapshot


def load(xtc=None, pdb=None, keyfile=None, *, spacing=None, dimensions=None,
         hardwall=None, temperature=None, start=None, stop=None, step=None,
         n_frames=None, verbose=False):
    """Load a PIMMS trajectory and return a :class:`LatticeTrajectory`.

    Parameters
    ----------
    xtc : str, optional
        Path to the XTC trajectory holding the coordinates. Requires ``pdb``,
        because mdtraj needs a topology to read it (default ``None``).
    pdb : str, optional
        Path to the PDB giving the topology, and on its own a single frame (e.g.
        ``START.pdb``). Its bead order matches the XTC exactly (default
        ``None``).
    keyfile : str, optional
        Path to the PIMMS keyfile. Optional, but authoritative for spacing,
        dimensions, hardwall and chain types (default ``None``).
    spacing : float, optional
        Lattice spacing in angstroms, overriding the keyfile
        ``LATTICE_TO_ANGSTROMS``. Default ``None``: taken from the keyfile, or
        ``DEFAULT_SPACING`` when there is none. A bool is refused.
    dimensions : sequence of int, optional
        Box extent in lattice units (2 or 3 positive integers), overriding
        everything below. Default ``None``: the resized-equilibration box for an
        ``eq_`` trajectory, else the restart file's box under
        ``RESTART_OVERRIDE_DIMENSIONS``, else the keyfile ``DIMENSIONS``, else
        inferred from the trajectory's own box record.
    hardwall : bool, optional
        Whether the run used hard walls, overriding everything below. Default
        ``None``: ``True`` for an ``eq_`` trajectory (that phase is always run
        under hard walls), else the restart file's flag under
        ``RESTART_OVERRIDE_HARDWALL``, else the keyfile ``HARDWALL``, else
        ``False``.
    temperature : float, optional
        Override (or supply, when no keyfile is given) the simulation
        temperature; must be a finite positive number. Only needed by the
        surface-tension estimators, which use it for :math:`k_BT`. Default
        ``None``: the keyfile ``TEMPERATURE``, or ``QUENCH_END`` for a
        ``QUENCH_RUN`` keyfile, or ``None`` (unknown) without a keyfile.
    start, stop, step : int, optional
        Frame slice applied at load time (default ``None``, i.e. every frame).
        Only the frames the selection keeps are read from the XTC and held in
        memory.
    n_frames : int, optional
        If given (and smaller), evenly subsample down to this many frames
        (default ``None``).
    verbose : bool, optional
        Print a one-line load summary, including the lattice round-off residual
        (default ``False``).

    Returns
    -------
    LatticeTrajectory
        The loaded trajectory, wrapping a
        :class:`~pimms.lemonade._store.TrajectoryStore`.

    Raises
    ------
    ValueError
        If neither ``xtc`` nor ``pdb`` is given, if ``xtc`` is given without
        ``pdb``, if ``n_frames`` is not a positive integer, if ``spacing`` is not
        finite and positive, if the box dimensions are neither given, in the
        keyfile, nor recorded in the trajectory (or are not 2 or 3 positive
        integers), if ``hardwall`` is not a bool, if ``temperature`` is not a
        finite positive number, if the coordinates do not fit in int32 lattice
        coordinates, if the topology's bead count disagrees with the trajectory
        (mdtraj itself refuses a PDB/XTC pair with different bead counts while
        reading them), if the ``start``/``stop``/``step`` selection
        keeps no frames, if the box inferred from the trajectory's own record is
        smaller than one lattice site at the spacing in use, if the XTC path is
        not ASCII and cannot be reached through a temporary link either, or if
        the keyfile sets ``RESTART_OVERRIDE_HARDWALL`` /
        ``RESTART_OVERRIDE_DIMENSIONS`` but its ``RESTART_FILE`` cannot be found
        or read (the keyfile values are then known to be the wrong ones, so
        guessing is not an option).

    Notes
    -----
    A large lattice round-off residual, bonds in the first loaded frame that
    are not single lattice steps (a spacing that is a whole multiple of the
    assumed one, which leaves no residual; the warning gives the spacing the
    bond lengths point to), a box or box *dimensionality* that disagrees with
    the trajectory's own record, a keyfile CHAIN block - or, for a restart
    keyfile, a restart file - that does not match the PDB, beads outside a
    hardwall box, an XTC that ends in an incomplete frame, the forced hardwall
    of a resized-equilibration trajectory, and a quench run's use of
    ``QUENCH_END`` as the temperature all raise a warning rather than an error.

    The bond check has nothing to test in a system made only of single-bead
    chains, and says nothing there.

    The XTC of a run that is still going, or that was killed, usually ends in a
    partly written frame. It is loaded up to its last complete frame, with a
    warning giving the number of frames recovered, and ``start`` / ``stop`` /
    ``step`` / ``n_frames`` then apply to those frames. A file with no torn
    frame is read directly and raises no such warning.

    Under hard walls a bead outside the box cannot have come from the
    simulation. For a single-chain system it is what ``AUTOCENTER`` wrote
    before PIMMS 1.0.8, and such a frame is translated rigidly back into the
    box (and warned about) rather than wrapped through the wall, which tore
    the chain. With several chains, or a frame wider than the box, the box or
    the boundary condition is wrong: the coordinates are wrapped as before, and
    the warning says to check ``hardwall=`` and ``dimensions=``.

    The keyfile is read with ``parse_only=True``, which does not run PIMMS's
    restart reconciliation, so the keywords that make the literal keyfile differ
    from the run PIMMS actually performed are resolved here instead: the
    ``RESTART_OVERRIDE_*`` keywords are honoured by reading the restart file, the
    ``eq_`` files of a ``RESIZED_EQUILIBRATION`` run are loaded in the compact
    box with hard walls, and under a ``RESTART_FILE`` the keyfile ``CHAIN`` lines
    (which PIMMS discards) are not applied.

    Keyfile chain types are only applied when the ``CHAIN`` lines expand onto
    the trajectory's chains in order and, if the PDB carries chain identifiers,
    reproduce the PDB's own partition of chains into types (PIMMS writes one
    identifier per chain type). A keyfile listing the same ``(count, sequence)``
    types in a different order - the ``keyfile_used.kf`` of a restart run with
    ``EXTRA_CHAIN`` chains, for example - keeps the PDB labels without a
    warning, because they are then the right ones.
    """
    import mdtraj as md

    if xtc is None and pdb is None:
        raise ValueError("load() needs at least an xtc (with a pdb topology) or a pdb")
    if n_frames is not None and (isinstance(n_frames, (bool, np.bool_)) or
                                 not isinstance(n_frames, (int, np.integer)) or
                                 n_frames < 1):
        raise ValueError("n_frames must be a positive integer")
    if xtc is not None:
        if pdb is None:
            raise ValueError("loading an xtc requires a pdb topology - pass pdb=...")
        # the topology only; the coordinates are read below, once the spacing and
        # the frame selection are known, so that only the kept frames are decoded
        traj = None
        with _without_dropped_cell_warning():
            md_topology = md.load_topology(pdb)
    else:
        with _without_dropped_cell_warning():
            traj = md.load(pdb)
        md_topology = traj.topology

    keydict = None
    if keyfile is not None:
        from pimms.keyfile_parser import KeyFileParser
        keydict = KeyFileParser(keyfile, parse_only=True).keyword_lookup

    # lattice spacing (angstroms). NB: KeyFileParser(parse_only=True) does not fill
    # defaults, so optional keys are read with .get(). A bool is refused: True used
    # to be read as a spacing of 1.0 A.
    if spacing is None:
        spacing = float(keydict.get("LATTICE_TO_ANGSTROMS", DEFAULT_SPACING)) if keydict else DEFAULT_SPACING
    if isinstance(spacing, (bool, np.bool_)):
        raise ValueError("spacing must be a finite positive number")
    try:
        spacing = float(spacing)
    except (TypeError, ValueError):
        raise ValueError("spacing must be a finite positive number")
    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be a finite positive number")
    int32_info = np.iinfo(np.int32)

    # frame selection, and coordinates (nm) -> integer lattice. The XTC is read a
    # block of kept frames at a time; a lone PDB is a frame or two and is
    # converted in one go.
    if traj is None:
        lattice, times, box_nm, residual, first_xyz = _read_xtc_lattice(
            xtc, md_topology, spacing, start, stop, step, n_frames)
    else:
        selected = _select_frames(traj.n_frames, start, stop, step, n_frames)
        lattice, residual = _nm_to_lattice(np.asarray(traj.xyz)[selected], spacing)
        times = np.asarray(traj.time, dtype=np.float64)[selected]
        first_xyz = np.asarray(traj.xyz)[selected[0]]
        box_nm = (np.array(traj.unitcell_lengths[selected[0]])
                  if traj.unitcell_lengths is not None else None)
    if box_nm is None:
        # mdtraj drops the CRYST1 record of a box it takes for a placeholder,
        # which is every PIMMS box at a sub-angstrom lattice spacing. Before 1.0.8
        # a SAVE_AT_END trajectory took its unit cell from that same mdtraj read,
        # so such an XTC carries no box either; the PDB's own record is the
        # run's box in both cases.
        box_nm = _pdb_cryst1_nm(pdb)

    # A large round-off residual means the coordinates do not sit on the integer
    # lattice at this spacing - almost always a wrong/omitted LATTICE_TO_ANGSTROMS
    # (or wrong spacing=). That silently corrupts the recovered lattice (non-unit
    # bonds, wrong inferred box), so warn ALWAYS, not only under verbose.
    if residual > 0.05:
        warnings.warn(
            f"lemonade.load: lattice round-off residual is {residual:.3g} (>0.05); the "
            f"coordinates do not fit the integer lattice at spacing {spacing} A. The "
            f"recovered lattice is probably corrupted - pass the right spacing=/keyfile "
            f"(LATTICE_TO_ANGSTROMS).", stacklevel=2)

    # Is this the compact equilibration phase of a RESIZED_EQUILIBRATION run? That
    # phase is ALWAYS hardwall and always in the resized box, whatever the keyfile
    # says (simulation.py:219-222), so its frames must not be analysed as periodic
    # in the production box. The filename is the discriminator rather than the box,
    # because RESIZED_EQUILIBRATION == DIMENSIONS is legal and a box-equality test
    # would then force hardwall onto a genuinely periodic production trajectory.
    is_eq_phase = _is_equilibration_file(xtc) or _is_equilibration_file(pdb)
    eq_notes = []

    # The restart snapshot is only opened if something actually needs it, and only
    # once; a missing restart file must not break a load that does not depend on it.
    _restart_cache = {}

    def _restart(keyword):
        """Memoised :func:`_restart_snapshot` for this load.

        Parameters
        ----------
        keyword : str
            The override keyword that needs the snapshot, quoted in any error.

        Returns
        -------
        pimms.restart.RestartObject
            The parsed restart file, read at most once per ``load()`` call.
        """
        if "snapshot" not in _restart_cache:
            _restart_cache["snapshot"] = _restart_snapshot(keydict, keyfile, keyword)
        return _restart_cache["snapshot"]

    _has_restart = bool(keydict is not None and keydict.get("RESTART_FILE"))

    # box dimensions
    if dimensions is None:
        if is_eq_phase and keydict is not None and keydict.get("RESIZED_EQUILIBRATION"):
            dimensions = tuple(int(d) for d in keydict["RESIZED_EQUILIBRATION"])
            eq_notes.append(f"box set to the RESIZED_EQUILIBRATION box {dimensions} "
                            f"rather than DIMENSIONS")
        elif _has_restart and keydict.get("RESTART_OVERRIDE_DIMENSIONS"):
            # RESTART_OVERRIDE_DIMENSIONS discards the keyfile DIMENSIONS; reading
            # them literally here put every minimum-image distance in the wrong box.
            dimensions = tuple(int(d) for d in _restart("RESTART_OVERRIDE_DIMENSIONS").dimensions)
        elif keydict is not None and keydict.get("DIMENSIONS"):
            dimensions = tuple(int(d) for d in keydict["DIMENSIONS"])
        else:
            if box_nm is None:
                raise ValueError("trajectory has no box and no keyfile/dimensions were given")
            dimensions = tuple(int(round(b / (spacing * 0.1))) for b in box_nm)
            # PIMMS writes a 2D system with a z period of exactly one lattice unit
            # (a real 3D box is never that thin), so the box record - not the
            # bead coordinates - decides the dimensionality: a 3D configuration
            # that happens to lie in the z = 0 plane must stay 3D.
            if len(dimensions) == 3 and dimensions[2] <= 1:
                dimensions = dimensions[:2]
            if any(d < 1 for d in dimensions):
                raise ValueError(
                    f"lemonade.load: the trajectory's own box record "
                    f"({', '.join(f'{10 * b:g}' for b in box_nm)} A) is smaller than one "
                    f"lattice site at spacing {spacing} A, so that spacing cannot be "
                    f"the one the run used - pass spacing= or the keyfile "
                    f"(LATTICE_TO_ANGSTROMS), or dimensions= explicitly.")
    try:
        dimensions = tuple(dimensions)
    except TypeError:
        raise ValueError("dimensions must be a 2D or 3D sequence of positive integers")
    if (len(dimensions) not in (2, 3) or
            any(isinstance(d, (bool, np.bool_)) or
                not isinstance(d, numbers.Integral) or d <= 0 or d > int32_info.max
                for d in dimensions)):
        raise ValueError("dimensions must be a 2D or 3D sequence of positive integers")
    dimensions = tuple(int(d) for d in dimensions)
    n_dim = len(dimensions)

    # cross-check against the trajectory's own box record where one exists: a
    # keyfile DIMENSIONS that disagrees with the XTC/CRYST1 box (wrong keyfile,
    # or the eq_ trajectory of a RESIZED_EQUILIBRATION run, written in the
    # smaller box) silently breaks every minimum-image / cluster / profile
    # calculation, so surface it loudly.
    if box_nm is not None:
        _record = tuple(int(round(b / (spacing * 0.1))) for b in box_nm)
        # Dimensionality first. The old check truncated the box record to the length
        # of the resolved `dimensions` before comparing, so 2D dimensions against a
        # 3D trajectory compared only x and y, matched, and stayed silent while the
        # whole z column was thrown away below. The z period of exactly one lattice
        # unit is the same 2D marker used when inferring dimensions above, and it is
        # in the one column the truncated comparison never read.
        _record_n_dim = 2 if (len(_record) == 3 and _record[2] <= 1) else len(_record)
        if _record_n_dim != n_dim:
            _discarded = ("the z coordinate of every bead is being DISCARDED"
                          if n_dim == 2 else
                          "a flat 2D trajectory is being analysed in a 3D box")
            warnings.warn(
                f"lemonade.load: {n_dim}D box dimensions {dimensions} "
                f"(keyfile/argument) disagree with the trajectory's own "
                f"{_record_n_dim}D box record {_record[:_record_n_dim]} - "
                f"{_discarded}, which changes every distance, cluster and profile. "
                f"Check you are loading the matching keyfile/trajectory pair.",
                stacklevel=2)
        else:
            _box_dims = _record[:n_dim]
            if _box_dims != tuple(dimensions[:len(_box_dims)]):
                warnings.warn(
                    f"lemonade.load: box dimensions {dimensions} (keyfile/argument) "
                    f"disagree with the trajectory's own box record {_box_dims}. All "
                    f"periodic-image and cluster calculations will use {dimensions} - "
                    f"check you are loading the matching keyfile/trajectory pair "
                    f"(eq_ trajectories from RESIZED_EQUILIBRATION runs use the "
                    f"smaller equilibration box).", stacklevel=2)

    # topology from the PDB (exact XTC bead order); keyfile refines chain types
    topology = Topology.from_mdtraj(md_topology)
    pdb_labelled = pdb_chain_labels(md_topology) is not None
    keyfile_types_applied = False
    # Under a RESTART_FILE the keyfile CHAIN lines are NOT what the run used: PIMMS
    # discards them and rebuilds the composition from the snapshot
    # (keyfile_parser.py:1950-2015), and an EXTRA_CHAIN whose sequence already
    # exists joins that existing chain type instead of defining a new one
    # (restart.py:328-348). Numbering each keyfile line as a fresh type therefore
    # split PIMMS's merged types apart (and, in the converse case, merged two real
    # types into one), silently mislabelling every per-type average. The PDB chain
    # identifiers are PIMMS's own chainType order, so they are the ones to keep.
    if keydict is not None and keydict.get("CHAIN") and not _has_restart:
        specs = list(keydict["CHAIN"]) + list(keydict.get("EXTRA_CHAIN") or [])
        # PIMMS upper-cases CHAIN sequences during sanitisation (the default
        # CASE_INSENSITIVE_CHAINS=True), but parse_only=True skips that step -
        # so a lower-case keyfile would silently fail to match the (upper-case)
        # PDB residue names and the keyfile types would be dropped.
        if keydict.get("CASE_INSENSITIVE_CHAINS", True):
            specs = [[n, str(seq).upper()] for (n, seq) in specs]
        # The PDB chain identifiers, when present, are PIMMS's own partition of
        # the chains into types, so the keyfile types are only taken if they
        # reproduce it. The keyfile_used.kf of a restart run lists one CHAIN line
        # per type while the trajectory has its EXTRA_CHAIN chains appended at
        # the end, so its lines do not expand onto the chains in order; when two
        # types share a sequence they expanded onto the wrong chains without a
        # word. A keyfile that holds the same (count, sequence) types in another
        # order describes this run, and the PDB labels are then already right,
        # so they are kept without a warning.
        typed_topology = topology.with_keyfile_types(specs, labelled=pdb_labelled)
        if typed_topology is topology:
            if not (pdb_labelled and topology.matches_keyfile_composition(specs)):
                warnings.warn(
                    "lemonade.load: the keyfile CHAIN/EXTRA_CHAIN specification does "
                    "not match the PDB topology, so keyfile chain types could not be "
                    "applied. The PDB chain identifiers will be used instead; check "
                    "that the keyfile and trajectory belong to the same run.",
                    stacklevel=2)
        else:
            keyfile_types_applied = True
        topology = typed_topology
    elif _has_restart:
        # The CHAIN lines say nothing here, but the restart file does: its chains,
        # followed by the EXTRA_CHAIN chains, ARE the run's composition. That is
        # the only thing a restart keyfile can be checked against (a keyfile from
        # another run used to load without a word), and the only place the types
        # merged under the 62nd PDB identifier can be recovered from. A restart
        # file that cannot be found is not an error for a load that does not
        # otherwise need it, so this is skipped silently then.
        try:
            _snapshot = _restart("RESTART_FILE")
        except ValueError:
            _snapshot = None
        if _snapshot is not None:
            _extra = list(keydict.get("EXTRA_CHAIN") or [])
            if keydict.get("CASE_INSENSITIVE_CHAINS", True):
                _extra = [[n, str(seq).upper()] for (n, seq) in _extra]
            _seqs, _types, _n_snapshot = _restart_chain_types(_snapshot, _extra)
            # the run's own restart.pimms, written over the file it started from,
            # already holds the extra chains: accept the snapshot on its own too
            if _seqs == topology.sequences:
                _restart_types = _types
            elif _seqs[:_n_snapshot] == topology.sequences:
                _restart_types = _types[:_n_snapshot]
            else:
                _restart_types = None
                warnings.warn(
                    f"lemonade.load: the keyfile's restart file (RESTART_FILE : "
                    f"{keydict['RESTART_FILE']}) holds {_n_snapshot} chains, "
                    f"{len(_seqs)} with the keyfile's EXTRA_CHAIN lines, and their "
                    f"sequences do not match the {topology.n_chains} chains of the PDB "
                    f"topology, so this keyfile does not describe this trajectory. The "
                    f"PDB chain identifiers will be used; check that the keyfile and "
                    f"trajectory belong to the same run.", stacklevel=2)
            if (_restart_types is not None and pdb_labelled and
                    len(set(topology.chain_types.tolist())) >= _N_PDB_CHAIN_IDS):
                typed_topology = topology.with_chain_types(_restart_types, labelled=True)
                if typed_topology is not topology:
                    keyfile_types_applied = True
                topology = typed_topology
    # PIMMS has 62 chain identifiers (A-Z, a-z, 0-9) and every chain type past
    # the 62nd shares the last one, so a PDB using all 62 MAY hide merged types;
    # only the keyfile CHAIN lines (or the restart file) can tell. (A PDB with a
    # blank chain column is typed by sequence instead, and cannot have run out of
    # identifiers.)
    if (not keyfile_types_applied and pdb_labelled
            and len(set(int(t) for t in topology.chain_types)) >= _N_PDB_CHAIN_IDS):
        if keydict is None:
            _how = "pass keyfile= to recover the true chain types."
        elif _has_restart:
            _how = ("the keyfile given starts from a RESTART_FILE, whose chain types "
                    "could not be read back (the restart file was not found beside the "
                    "keyfile, or does not match this trajectory); pass the run's own "
                    "keyfile_used.kf instead to recover the true chain types.")
        else:
            _how = ("the keyfile given could not be matched to the PDB chains, so it "
                    "did not resolve them; pass the run's own keyfile_used.kf to "
                    "recover the true chain types.")
        warnings.warn(
            "lemonade.load: the PDB uses all 62 PIMMS chain identifiers, so any chain "
            "type past the 62nd shares a label with another and would have been merged; "
            + _how, stacklevel=2)
    if topology.n_beads != lattice.shape[1]:
        raise ValueError(f"topology describes {topology.n_beads} beads but the "
                         f"trajectory has {lattice.shape[1]}")

    if hardwall is None:
        if is_eq_phase:
            hardwall = True
            eq_notes.append("hardwall set to True (the compact equilibration phase "
                            "is always run under hard walls)")
        elif _has_restart and keydict.get("RESTART_OVERRIDE_HARDWALL"):
            # RESTART_OVERRIDE_HARDWALL discards the keyfile HARDWALL. Taking the
            # keyfile value literally here was silent in both directions and made
            # the connected-component search join (or refuse to join) chains through
            # box faces the run never had.
            hardwall = bool(_restart("RESTART_OVERRIDE_HARDWALL").hardwall)
        else:
            hardwall = bool(keydict.get("HARDWALL", False)) if keydict else False
    elif not isinstance(hardwall, (bool, np.bool_)):
        raise ValueError("hardwall must be True or False")

    # canonicalise into the box (agnostic to whether the trajectory was written
    # wrapped or PBC-unwrapped); lemonade re-derives whole chains itself
    _dims = np.array(dimensions, dtype=np.int32)
    _in_box = lattice[..., :n_dim]
    if hardwall and lattice.size:
        # Nothing crosses a hard wall, so a bead outside the box was not put there
        # by the simulation.
        _lo = _in_box.min(axis=1).astype(np.int64)
        _hi = _in_box.max(axis=1).astype(np.int64)
        _outside = ((_lo < 0) | (_hi >= _dims)).any(axis=1)
        _wide = _outside & ~((_hi - _lo) < _dims).all(axis=1)
        if _outside.any() and topology.n_chains == 1 and not _wide.any():
            # One chain, and every such frame would fit: this is what AUTOCENTER
            # wrote before 1.0.8, when the shift that centres the chain was not
            # kept inside the walls. Wrapping such a bead "through" the wall tears
            # the chain in two; the frame is a rigid translation of a
            # configuration that did fit, so it is translated back until it fits
            # (the same shift the engine itself now applies when it writes).
            _shift = np.where(_lo < 0, -_lo, np.where(_hi >= _dims, _dims - 1 - _hi, 0))
            for _f in np.flatnonzero(_outside):
                _in_box[_f] += _shift[_f].astype(np.int32)
            warnings.warn(
                f"lemonade.load: {int(_outside.sum())} of {lattice.shape[0]} frames have "
                f"beads outside the hardwall box {dimensions}. For a single chain that "
                f"is what PIMMS before 1.0.8 wrote under AUTOCENTER with HARDWALL (the "
                f"centred chain could stick out through a wall), so those frames have "
                f"been translated rigidly back inside the box, which keeps the chain "
                f"whole. Every frame of an AUTOCENTER run was re-centred when it was "
                f"written, so positions relative to the walls (density profiles along "
                f"a wall normal, distances to a wall) are not meaningful in any frame "
                f"of this trajectory.", stacklevel=2)
        elif _outside.any():
            # Several chains (AUTOCENTER never acted on those), or a frame wider
            # than the box: the box or the boundary condition is not the run's.
            # The coordinates are wrapped, as they always were, but not silently.
            warnings.warn(
                f"lemonade.load: {int(_outside.sum())} of {lattice.shape[0]} frames have "
                f"beads outside the hardwall box {dimensions}"
                + (f" ({int(_wide.sum())} of them wider than the box)"
                   if _wide.any() else "")
                + ", which a hardwall run in that box cannot produce. The coordinates "
                  "have been wrapped into the box as if it were periodic, which tears "
                  "chains and clusters through the walls: check hardwall= and "
                  "dimensions= (or the keyfile HARDWALL / DIMENSIONS) against the run "
                  "that wrote this trajectory.", stacklevel=2)
    np.mod(_in_box, _dims, out=_in_box)
    if n_dim == 2:
        lattice[..., 2] = 0

    # A spacing that is a whole multiple of the right one leaves no round-off
    # residual, so the residual check above cannot see it; the bonds can.
    _check_unit_bonds(lattice[0], topology.offsets, dimensions, spacing, first_xyz,
                      residual)

    if eq_notes:
        warnings.warn(
            f"lemonade.load: this is the resized-equilibration trajectory of a "
            f"RESIZED_EQUILIBRATION run ({', '.join(_EQUILIBRATION_FILENAMES)}), "
            f"which PIMMS runs in the compact box under FORCED hard walls whatever "
            f"the keyfile HARDWALL says, so " + "; ".join(eq_notes) +
            ". Pass hardwall= / dimensions= explicitly to override this.",
            stacklevel=2)
    if temperature is None and keydict is not None:
        if keydict.get("QUENCH_RUN"):
            # a quench run ignores TEMPERATURE: production is sampled at
            # QUENCH_END (frames written during the ramp sit at intermediate
            # temperatures)
            temperature = keydict.get("QUENCH_END")
            if temperature is None:
                warnings.warn(
                    "lemonade.load: the keyfile describes a QUENCH_RUN but sets no "
                    "QUENCH_END, so no trajectory temperature was recorded (TEMPERATURE "
                    "is ignored by PIMMS in a quench). Pass temperature= explicitly for "
                    "temperature-dependent analyses (surface tension).", stacklevel=2)
            else:
                temperature = float(temperature)
                warnings.warn(
                    "lemonade.load: the keyfile describes a QUENCH_RUN, whose production "
                    f"phase runs at QUENCH_END = {temperature} (TEMPERATURE is ignored by "
                    "PIMMS in a quench). Using QUENCH_END as the trajectory temperature; "
                    "frames written during the ramp were sampled at intermediate "
                    "temperatures, so restrict temperature-dependent analyses (surface "
                    "tension) to post-ramp frames.", stacklevel=2)
        else:
            temperature = keydict.get("TEMPERATURE")

    store = TrajectoryStore(lattice, dimensions, spacing, bool(hardwall), topology,
                            times=times, temperature=temperature, copy=False)

    if verbose:
        print(f"[lemonade] {store.n_frames} frames, {store.n_chains} chains, "
              f"{store.n_beads} beads; box {dimensions}, spacing {spacing} A"
              f"{'' if residual < 1e-3 else f'  (WARNING lattice round-off {residual:.3g})'}")
    return LatticeTrajectory(store)
