"""
Regression tests for the lemonade defects found by the fifth 1.0.8 review.

The theme of F-I-1 to F-I-3 is that ``load()`` reads the keyfile with
``parse_only=True``, which does not run PIMMS's restart reconciliation, so the
*literal* keyfile can differ from the run PIMMS actually performed. Those three
need real two-stage PIMMS runs, and the clustering they produce is checked either
against PIMMS's own ``NUM_CLUSTERS.dat`` or against a union-find oracle written
here from the raw lattice positions - never against another lemonade call.

F-I-4 to F-I-6 are checked on hand-built systems with an analytic answer.
"""

import contextlib
import os
import shutil
import warnings

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import phase_separation as ps
from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology
from pimms.lemonade.trajectory import LatticeTrajectory


_PARAMS = ("ANGLE_PENALTY\tA\t0\t0\t0\n"
           "ANGLE_PENALTY\tB\t0\t0\t0\n"
           "A\tA\t-8\n"
           "B\tB\t-8\n"
           "A\tB\t2\n"
           "A\t0\t0\n"
           "B\t0\t0\n")

_BASE = ["PARAMETER_FILE : params.prm",
         "TEMPERATURE : 40",
         "PRINT_FREQ : 100000",
         "EN_FREQ : 100000",
         "ANALYSIS_FREQ : 100000",
         "MOVE_CRANKSHAFT : 0.6",
         "MOVE_SLITHER : 0.4",
         "TRAJECTORY_PBC_UNWRAP : False"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_RESTART_TYPE_CASES = [("merge", [0, 0, 0, 1, 1, 0, 0, 2]),
                       ("split", [0, 0, 0, 1, 1]),
                       ("extra_only", [0, 0, 0, 1, 1, 0, 0, 2]),
                       ("twin_extra", [0, 0, 0, 1, 1, 0, 0])]


def _run_pimms(dirpath, lines):
    """Run a short PIMMS simulation in ``dirpath`` and return it plus its keyfile.

    Parameters
    ----------
    dirpath : str or pathlib.Path
        Directory to run in. It is created if it does not exist, and the
        parameter file and keyfile are written into it.
    lines : list of str
        Keyfile lines to add to the shared ``_BASE`` block.

    Returns
    -------
    tuple
        ``(simulation, keyfile_path)``. The finished :class:`Simulation` carries
        the *effective* configuration (``sim.hardwall``,
        ``sim.LATTICE.dimensions``, the chains' ``chainType``), which is the
        ground truth these tests compare lemonade against.
    """
    dirpath = str(dirpath)
    os.makedirs(dirpath, exist_ok=True)
    with open(os.path.join(dirpath, "params.prm"), "w") as fh:
        fh.write(_PARAMS)
    keyfile = os.path.join(dirpath, "KEYFILE.kf")
    with open(keyfile, "w") as fh:
        fh.write("\n".join(_BASE + list(lines)) + "\n")

    cwd = os.getcwd()
    os.chdir(dirpath)
    try:
        from pimms.keyfile_parser import KeyFileParser
        from pimms.simulation import Simulation
        # these runs are dense enough for a cluster to wind the box, and the
        # engine's gather then raises its documented UserWarning (the run's own
        # log records it, which is the channel that matters). It is incidental
        # to every test built on this helper, so it is kept out of the output
        with contextlib.redirect_stdout(open(os.devnull, "w")), \
                warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="single-image gather: cluster percolates")
            sim = Simulation(KeyFileParser("KEYFILE.kf").keyword_lookup)
            sim.run_simulation()
    finally:
        os.chdir(cwd)
    return sim, keyfile


def _oracle_cluster_count(positions, chain_of_bead, dimensions, hardwall):
    """Count connected components of chains by union-find over Chebyshev-1 contacts.

    This is deliberately independent of every lemonade and PIMMS clustering
    routine: it works from the raw integer positions and a bead-to-chain map, and
    implements the boundary condition directly (minimum image under PBC, no
    contact through a wall under a hardwall).

    Parameters
    ----------
    positions : numpy.ndarray
        ``(n_beads, 3)`` integer lattice positions for one frame.
    chain_of_bead : numpy.ndarray
        ``(n_beads,)`` chain index of each bead.
    dimensions : tuple of int
        Box extent, 2 or 3 entries.
    hardwall : bool
        ``True`` for hard walls (no contact across a face), ``False`` for
        periodic boundaries.

    Returns
    -------
    int
        The number of connected components.
    """
    dims = np.asarray(dimensions, dtype=int)
    n_dim = len(dims)
    occupancy = {}
    for bead, pos in enumerate(positions):
        occupancy[tuple(int(v) for v in pos[:n_dim])] = int(chain_of_bead[bead])

    parent = {c: c for c in occupancy.values()}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    offsets = np.array(np.meshgrid(*[[-1, 0, 1]] * n_dim,
                                   indexing="ij")).reshape(n_dim, -1).T
    for site, chain in occupancy.items():
        for offset in offsets:
            if not offset.any():
                continue
            neighbour = np.asarray(site) + offset
            if hardwall:
                if np.any(neighbour < 0) or np.any(neighbour >= dims):
                    continue
            else:
                neighbour = np.mod(neighbour, dims)
            key = tuple(int(v) for v in neighbour)
            if key in occupancy:
                union(chain, occupancy[key])
    return len(set(find(c) for c in parent))


def _oracle_counts_for(traj, dimensions, hardwall):
    """Per-frame union-find cluster count of ``traj`` under a stated boundary."""
    chain_of_bead = np.repeat(np.arange(traj.n_chains),
                              [len(s) for s in traj.sequences])
    return [_oracle_cluster_count(traj.positions[f], chain_of_bead,
                                  dimensions, hardwall)
            for f in range(traj.n_frames)]


def _num_clusters_dat(dirpath, xtc_freq, n_frames):
    """Read PIMMS's own ``NUM_CLUSTERS.dat`` as ``{frame index: n_clusters}``."""
    out = {}
    with open(os.path.join(str(dirpath), "NUM_CLUSTERS.dat")) as fh:
        for line in fh:
            if not line.strip():
                continue
            step, n_clusters = line.split()
            frame = int(step) // xtc_freq
            if frame < n_frames:
                out[frame] = int(n_clusters)
    return out


def _pimms_chain_types(sim):
    """PIMMS's own per-chain ``chainType``, in the trajectory's chain order."""
    return [int(sim.LATTICE.chains[c].chainType) for c in sorted(sim.LATTICE.chains)]


def _make_traj(frames, sequences, dimensions, hardwall):
    """Wrap hand-built integer positions in a LatticeTrajectory."""
    store = TrajectoryStore(np.asarray(frames, dtype=np.int32), tuple(dimensions),
                            3.65, hardwall, Topology(list(sequences)))
    return LatticeTrajectory(store)


# ---------------------------------------------------------------------------
# F-I-1: RESTART_OVERRIDE_HARDWALL / _DIMENSIONS
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def restart_override_hardwall_run(tmp_path_factory):
    """A hardwall run, restarted with ``HARDWALL : False`` + the hardwall override.

    PIMMS therefore runs hard walls while the keyfile literally says periodic -
    the exact situation ``RESTART_OVERRIDE_HARDWALL`` exists to create.
    """
    root = tmp_path_factory.mktemp("r5_restart_hw")
    original = root / "original"
    _run_pimms(original, ["DIMENSIONS : 9 9 9", "HARDWALL : True", "SEED : 17",
                          "N_STEPS : 20", "EQUILIBRATION : 2",
                          "XTC_FREQ : 100000", "RESTART_FREQ : 10",
                          "CHAIN : 40 AA", "CHAIN : 4 AAA"])

    restarted = root / "restarted"
    os.makedirs(restarted, exist_ok=True)
    shutil.copy(original / "restart.pimms", restarted / "restart.pimms")
    sim, keyfile = _run_pimms(restarted, [
        "DIMENSIONS : 9 9 9", "HARDWALL : False",
        "RESTART_FILE : restart.pimms", "RESTART_OVERRIDE_HARDWALL : True",
        "SEED : 5", "N_STEPS : 24", "EQUILIBRATION : 2",
        "XTC_FREQ : 2", "ANA_CLUSTER : 2", "RESTART_FREQ : 100000"])
    return sim, str(restarted), keyfile


def test_restart_override_hardwall_is_taken_from_the_restart_file(
        restart_override_hardwall_run):
    """The boundary condition must be the one PIMMS ran, not the keyfile's.

    ``RESTART_OVERRIDE_HARDWALL`` is a no-op unless the keyfile and the snapshot
    disagree, so whenever it does anything the literal keyfile value is
    guaranteed wrong. Before the fix this was silent in both directions and
    lemonade joined chains through box faces the run never had.
    """
    sim, rundir, keyfile = restart_override_hardwall_run
    assert sim.hardwall is True                       # what PIMMS actually ran

    traj = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                         pdb=os.path.join(rundir, "START.pdb"), keyfile=keyfile)
    assert traj.hardwall is True

    # PIMMS's own cluster count, frame by frame - not another lemonade call
    expected = _num_clusters_dat(rundir, 2, traj.n_frames)
    assert len(expected) >= 8
    got = ps.number_of_clusters(traj, min_beads=1)
    assert {f: int(got[f]) for f in expected} == expected

    # and the old behaviour (the literal keyfile HARDWALL : False) does not agree
    literal = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                            pdb=os.path.join(rundir, "START.pdb"),
                            keyfile=keyfile, hardwall=False)
    was = ps.number_of_clusters(literal, min_beads=1)
    assert {f: int(was[f]) for f in expected} != expected


@pytest.fixture(scope="module")
def restart_override_dimensions_run(tmp_path_factory):
    """A 9^3 periodic run restarted with a 14^3 keyfile plus the box override."""
    root = tmp_path_factory.mktemp("r5_restart_dims")
    original = root / "original"
    _run_pimms(original, ["DIMENSIONS : 9 9 9", "HARDWALL : False", "SEED : 21",
                          "N_STEPS : 20", "EQUILIBRATION : 2",
                          "XTC_FREQ : 100000", "RESTART_FREQ : 10",
                          "CHAIN : 40 AA", "CHAIN : 4 AAA"])

    restarted = root / "restarted"
    os.makedirs(restarted, exist_ok=True)
    shutil.copy(original / "restart.pimms", restarted / "restart.pimms")
    sim, keyfile = _run_pimms(restarted, [
        "DIMENSIONS : 14 14 14", "HARDWALL : False",
        "RESTART_FILE : restart.pimms", "RESTART_OVERRIDE_DIMENSIONS : True",
        "SEED : 5", "N_STEPS : 24", "EQUILIBRATION : 2",
        "XTC_FREQ : 2", "ANA_CLUSTER : 2", "RESTART_FREQ : 100000"])
    return sim, str(restarted), keyfile


def test_restart_override_dimensions_is_taken_from_the_restart_file(
        restart_override_dimensions_run):
    """The box must be the restart file's, not the keyfile's 14^3."""
    sim, rundir, keyfile = restart_override_dimensions_run
    assert tuple(sim.LATTICE.dimensions) == (9, 9, 9)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                             pdb=os.path.join(rundir, "START.pdb"), keyfile=keyfile)
    assert traj.dimensions == (9, 9, 9)
    # the box now agrees with the trajectory's own record, so nothing is suspicious
    assert not [w for w in caught if "box dimensions" in str(w.message)]

    expected = _num_clusters_dat(rundir, 2, traj.n_frames)
    assert len(expected) >= 8
    got = ps.number_of_clusters(traj, min_beads=1)
    assert {f: int(got[f]) for f in expected} == expected

    # the literal keyfile box (14^3) invents periodic images and disagrees - and
    # the loader is expected to say so, which is asserted rather than tolerated
    with pytest.warns(UserWarning, match="box dimensions"):
        literal = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                                pdb=os.path.join(rundir, "START.pdb"),
                                keyfile=keyfile, dimensions=(14, 14, 14))
    was = ps.number_of_clusters(literal, min_beads=1)
    assert {f: int(was[f]) for f in expected} != expected


def test_restart_override_with_a_missing_restart_file_raises(
        restart_override_hardwall_run, tmp_path):
    """With the snapshot gone the keyfile value is known-wrong, so refuse to guess."""
    _, rundir, _ = restart_override_hardwall_run
    for name in ("traj.xtc", "START.pdb", "KEYFILE.kf"):
        shutil.copy(os.path.join(rundir, name), tmp_path / name)
    # deliberately do NOT copy restart.pimms
    with pytest.raises(ValueError, match="RESTART_OVERRIDE_HARDWALL"):
        lemonade.load(xtc=str(tmp_path / "traj.xtc"), pdb=str(tmp_path / "START.pdb"),
                      keyfile=str(tmp_path / "KEYFILE.kf"))
    # ... but an explicit hardwall= still loads, since nothing is being guessed
    traj = lemonade.load(xtc=str(tmp_path / "traj.xtc"), pdb=str(tmp_path / "START.pdb"),
                         keyfile=str(tmp_path / "KEYFILE.kf"), hardwall=True)
    assert traj.hardwall is True


def test_restart_file_is_resolved_against_the_keyfile_directory(
        restart_override_hardwall_run, tmp_path, monkeypatch):
    """A relative RESTART_FILE lives beside the keyfile, not in the analysis cwd."""
    _, rundir, keyfile = restart_override_hardwall_run
    monkeypatch.chdir(tmp_path)                    # analyse from somewhere else
    traj = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                         pdb=os.path.join(rundir, "START.pdb"), keyfile=keyfile)
    assert traj.hardwall is True


# ---------------------------------------------------------------------------
# F-I-2: the eq_ trajectory of a RESIZED_EQUILIBRATION run
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def resized_equilibration_run(tmp_path_factory):
    """A ``RESIZED_EQUILIBRATION`` run with ``HARDWALL : False`` and ``SAVE_EQ``."""
    d = tmp_path_factory.mktemp("r5_resized")
    sim, keyfile = _run_pimms(d, [
        "DIMENSIONS : 14 14 14", "RESIZED_EQUILIBRATION : 9 9 9",
        "HARDWALL : False", "SAVE_EQ : True", "SEED : 3",
        "N_STEPS : 30", "EQUILIBRATION : 16", "XTC_FREQ : 2",
        "RESTART_FREQ : 100000", "CHAIN : 40 AA", "CHAIN : 4 AAA"])
    return sim, str(d), keyfile


def test_resized_equilibration_trajectory_loads_hardwall_in_the_compact_box(
        resized_equilibration_run):
    """``eq_traj.xtc`` is hardwall in the compact box, with or without a keyfile.

    The compact phase is forced to hard walls by ``Simulation.__init__``
    regardless of the keyfile, so loading its frames as periodic joined chains
    through walls the run had. The detection keys on the FILENAME rather than on
    box equality, because ``RESIZED_EQUILIBRATION == DIMENSIONS`` is legal and a
    box-equality test would then force hard walls onto a genuinely periodic
    production trajectory (checked below).
    """
    sim, rundir, keyfile = resized_equilibration_run
    assert sim.hardwall is False and tuple(sim.LATTICE.dimensions) == (14, 14, 14)

    eq_xtc = os.path.join(rundir, "eq_traj.xtc")
    eq_pdb = os.path.join(rundir, "eq_START.pdb")

    for kwargs in ({"keyfile": keyfile}, {}):
        with pytest.warns(UserWarning, match="resized-equilibration"):
            traj = lemonade.load(xtc=eq_xtc, pdb=eq_pdb, **kwargs)
        assert traj.hardwall is True
        assert traj.dimensions == (9, 9, 9)
        got = [int(n) for n in ps.number_of_clusters(traj, min_beads=1)]
        assert got == _oracle_counts_for(traj, (9, 9, 9), hardwall=True)

    # the old behaviour - periodic in the compact box - is a different answer
    stale = lemonade.load(xtc=eq_xtc, pdb=eq_pdb, keyfile=keyfile,
                          hardwall=False, dimensions=(9, 9, 9))
    was = [int(n) for n in ps.number_of_clusters(stale, min_beads=1)]
    assert was != _oracle_counts_for(stale, (9, 9, 9), hardwall=True)

    # the PRODUCTION trajectory of the same run keeps its real (periodic) boundary
    production = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                               pdb=os.path.join(rundir, "START.pdb"), keyfile=keyfile)
    assert production.hardwall is False
    assert production.dimensions == (14, 14, 14)


def test_resized_equilibration_equal_to_dimensions_leaves_production_periodic(
        tmp_path):
    """The degenerate box-equality case a box-based detector would break."""
    sim, keyfile = _run_pimms(tmp_path, [
        "DIMENSIONS : 9 9 9", "RESIZED_EQUILIBRATION : 9 9 9", "HARDWALL : False",
        "SEED : 3", "N_STEPS : 20", "EQUILIBRATION : 8", "XTC_FREQ : 4",
        "RESTART_FREQ : 100000", "CHAIN : 20 AA"])
    assert sim.hardwall is False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(xtc=str(tmp_path / "traj.xtc"),
                             pdb=str(tmp_path / "START.pdb"), keyfile=keyfile)
    assert traj.hardwall is False
    assert not [w for w in caught if "resized-equilibration" in str(w.message)]


# ---------------------------------------------------------------------------
# F-I-3: keyfile CHAIN lines under a RESTART_FILE
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def restart_chain_type_runs(tmp_path_factory):
    """Four restarts whose keyfile CHAIN lines disagree with PIMMS's own types.

    Each run directory also holds the ``keyfile_used.kf`` PIMMS writes at
    start-up, which drops ``RESTART_FILE`` and lists the merged composition as
    one CHAIN line per chain type - so its lines are NOT in the trajectory's
    chain order whenever EXTRA_CHAIN chains joined an earlier type.
    """
    root = tmp_path_factory.mktemp("r5_types")
    common = ["DIMENSIONS : 9 9 9", "HARDWALL : False", "SEED : 4",
              "N_STEPS : 12", "EQUILIBRATION : 2", "XTC_FREQ : 6"]

    def restart_from(source, tag, lines):
        d = root / tag
        os.makedirs(d, exist_ok=True)
        shutil.copy(source / "restart.pimms", d / "restart.pimms")
        sim, keyfile = _run_pimms(d, common + ["RESTART_FREQ : 100000",
                                               "RESTART_FILE : restart.pimms"] + lines)
        return sim, str(d), keyfile

    mixed = root / "orig_mixed"
    _run_pimms(mixed, common + ["RESTART_FREQ : 6", "CHAIN : 3 AAAA", "CHAIN : 2 BBBB"])
    twinned = root / "orig_twinned"
    _run_pimms(twinned, common + ["RESTART_FREQ : 6", "CHAIN : 3 AAAA", "CHAIN : 2 AAAA"])

    return {
        # PIMMS merges the two extra AAAA into type 0 and calls AB type 2; the
        # literal keyfile numbering would call them types 2, 2 and 3
        "merge": restart_from(mixed, "merge",
                              ["CHAIN : 3 AAAA", "CHAIN : 2 BBBB",
                               "EXTRA_CHAIN : 2 AAAA", "EXTRA_CHAIN : 1 AB"]),
        # the converse: two real types share a sequence, and the literal keyfile
        # would collapse them into one
        "split": restart_from(twinned, "split", ["CHAIN : 5 AAAA"]),
        # the common idiom - EXTRA_CHAIN only - which must stay correct AND silent
        "extra_only": restart_from(mixed, "extra_only",
                                   ["EXTRA_CHAIN : 2 AAAA", "EXTRA_CHAIN : 1 AB"]),
        # two types share a sequence AND extra chains join the first: the
        # keyfile_used.kf lines "5 AAAA" / "2 AAAA" expand onto the chains in
        # order with every sequence matching, but give the wrong partition
        "twin_extra": restart_from(twinned, "twin_extra", ["EXTRA_CHAIN : 2 AAAA"]),
    }


@pytest.mark.parametrize("which_keyfile", ["KEYFILE.kf", "keyfile_used.kf"])
@pytest.mark.parametrize("case,expected", _RESTART_TYPE_CASES)
def test_restart_keyfile_chain_lines_are_not_applied(restart_chain_type_runs,
                                                     case, expected, which_keyfile):
    """Chain types must follow PIMMS's restart rules, not the literal CHAIN lines.

    PIMMS discards ``CHAIN`` under a ``RESTART_FILE`` and merges an
    ``EXTRA_CHAIN`` whose sequence already exists into the existing type. Taking
    the keyfile lines literally split PIMMS's type 0 into three ("merge") and
    collapsed two real types into one ("split"), which silently mislabels every
    per-type average and breaks the join against PIMMS's ``CHAIN_<type>_*``
    output files.

    The same runs are also loaded with their ``keyfile_used.kf``, which has no
    ``RESTART_FILE`` and one CHAIN line per merged type. Expanding those lines in
    order put the "twin_extra" run's two appended chains into type 0 and two of
    the snapshot's type-0 chains into type 1 with no warning, and gave the
    "merge" run a false "does not match" warning; the PDB chain identifiers are
    PIMMS's own partition and must win, silently, in both.
    """
    sim, rundir, keyfile = restart_chain_type_runs[case]
    assert _pimms_chain_types(sim) == expected
    if which_keyfile == "keyfile_used.kf":
        keyfile = os.path.join(rundir, "keyfile_used.kf")
        assert os.path.isfile(keyfile)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        traj = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                             pdb=os.path.join(rundir, "START.pdb"), keyfile=keyfile)
    assert [int(t) for t in traj.chain_types] == expected
    # the mismatch warning is what saves a user whose restart keyfile CHAIN block
    # has been reordered, so it must not have been blanket-suppressed; these
    # keyfiles are consistent with their runs, so nothing may fire
    assert [str(w.message) for w in caught] == []


def test_keyfile_used_twin_types_are_the_ones_in_the_keyfile(restart_chain_type_runs):
    """The twin_extra keyfile_used.kf really does list two types of one sequence.

    Guards the case above against going vacuous: if PIMMS ever wrote the merged
    composition differently, the in-order expansion would no longer be the trap
    it is meant to test.
    """
    from pimms.keyfile_parser import KeyFileParser
    _sim, rundir, _keyfile = restart_chain_type_runs["twin_extra"]
    keydict = KeyFileParser(os.path.join(rundir, "keyfile_used.kf"),
                            parse_only=True).keyword_lookup
    assert "RESTART_FILE" not in keydict or not keydict["RESTART_FILE"]
    specs = [[int(n), str(seq)] for n, seq in keydict["CHAIN"]]
    specs += [[int(n), str(seq)] for n, seq in (keydict.get("EXTRA_CHAIN") or [])]
    assert specs == [[5, "AAAA"], [2, "AAAA"]]


def test_keyfile_with_a_different_composition_still_warns(restart_chain_type_runs,
                                                           tmp_path):
    """A real composition mismatch must still warn and keep the PDB labels.

    The keyfile_used.kf of the twin_extra run is rewritten to split the seven
    AAAA chains 4 + 3 instead of PIMMS's 5 + 2. Every sequence still matches in
    order, so only the partition check can catch it, and because the (count,
    sequence) types differ it is not the harmless reordering case either.
    """
    _sim, rundir, _keyfile = restart_chain_type_runs["twin_extra"]
    with open(os.path.join(rundir, "keyfile_used.kf")) as fh:
        lines = fh.read().splitlines()
    rewritten = []
    for line in lines:
        if line.strip().startswith("CHAIN"):
            if "CHAIN : 4 AAAA" not in rewritten:
                rewritten.extend(["CHAIN : 4 AAAA", "CHAIN : 3 AAAA"])
            continue
        rewritten.append(line)
    keyfile = tmp_path / "wrong_split.kf"
    keyfile.write_text("\n".join(rewritten) + "\n")
    shutil.copy(os.path.join(rundir, "params.prm"), tmp_path / "params.prm")

    with pytest.warns(UserWarning, match="does not match the PDB topology"):
        traj = lemonade.load(xtc=os.path.join(rundir, "traj.xtc"),
                             pdb=os.path.join(rundir, "START.pdb"),
                             keyfile=str(keyfile))
    assert [int(t) for t in traj.chain_types] == [0, 0, 0, 1, 1, 0, 0]


# ---------------------------------------------------------------------------
# F-I-4: wetting-film thickness is measured from the wall
# ---------------------------------------------------------------------------

def _film(planes, dims, wall):
    """Hand-build a fully occupied film of ``planes`` planes along z.

    Parameters
    ----------
    planes : sequence of int
        The z planes to fill completely.
    dims : tuple of int
        Box extent.
    wall : str
        Unused label kept for readability at the call site.

    Returns
    -------
    LatticeTrajectory
        A one-frame hardwall trajectory, one chain per (x, y) column.
    """
    Lx, Ly, _ = dims
    positions, sequences = [], []
    for x in range(Lx):
        for y in range(Ly):
            column = [[x, y, z] for z in planes]
            positions.extend(column)
            sequences.append("A" * len(column))
    return _make_traj([positions], sequences, dims, hardwall=True)


@pytest.mark.parametrize("thickness", [3, 5, 10])
def test_wetting_film_thickness_is_measured_from_the_wall(thickness):
    """A film of t planes must report thickness t, whichever wall it wets.

    Lattice sites are unit cells centred on integer coordinates, so the wall is
    at ``-0.5`` and a film filling planes ``0 .. t-1`` fills a slab of thickness
    exactly ``t``. Measuring from the centre of plane 0 instead reported ``t -
    0.5`` - 5 % at ten planes, 17 % at three, systematic and one-signed, and
    inconsistent with the two-interface fit on identical material (asserted
    here).
    """
    dims = (12, 12, 40)
    tolerance = 0.05

    low = _film(range(thickness), dims, "low")
    z, rho = ps.slab_density_profile(low, axis=2)
    low_fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert low_fit.success is True
    assert 2.0 * low_fit.half_width == pytest.approx(thickness, abs=tolerance)

    high = _film(range(dims[2] - thickness, dims[2]), dims, "high")
    z, rho = ps.slab_density_profile(high, axis=2)
    high_fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert high_fit.success is True
    assert 2.0 * high_fit.half_width == pytest.approx(thickness, abs=tolerance)

    # the same material as a FREE slab, fit with the two-interface model, must
    # give the same number - the two conventions were 0.5 apart before the fix
    start = (dims[2] - thickness) // 2
    free = _film(range(start, start + thickness), dims, "free")
    z, rho = ps.slab_density_profile(free, axis=2)
    free_fit = ps.fit_slab_profile(z, rho, hardwall=True)
    assert free_fit.success is True
    assert 2.0 * free_fit.half_width == pytest.approx(thickness, abs=tolerance)
    assert (2.0 * low_fit.half_width ==
            pytest.approx(2.0 * free_fit.half_width, abs=tolerance))


def test_a_two_plane_wetting_film_is_accepted():
    """The intended side effect of measuring from the wall.

    ``fit_radial_profile`` rejects a fitted radius below ``_MIN_DROPLET_RADIUS =
    2.0``, a guard against calling a handful of sites at a cluster centre a dense
    phase. A two-plane film is a real film of thickness exactly 2.0, so it now
    sits on the floor rather than 0.5 below it and is accepted. A one-plane film
    is still rejected, which is the case the guard is for.
    """
    dims = (12, 12, 40)
    fit = ps.fit_slab_profile(*ps.slab_density_profile(_film(range(2), dims, "low"),
                                                       axis=2), hardwall=True)
    assert fit.success is True
    assert 2.0 * fit.half_width == pytest.approx(2.0, abs=0.05)

    thin = ps.fit_slab_profile(*ps.slab_density_profile(_film(range(1), dims, "low"),
                                                        axis=2), hardwall=True)
    assert thin.success is False


# ---------------------------------------------------------------------------
# F-I-5: the radial profile starts at shell 1
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("radius", [3, 4])
def test_radial_density_profile_starts_at_shell_one(radius):
    """Entry ``k`` of the profile is the shell at Chebyshev distance ``k + 1``.

    A solid Chebyshev cube of half-width R about the COM fills shells 1..R
    completely and leaves shell R+1 empty, so the profile must read 1.0 up to
    index R-1 and 0.0 at index R. Shell 0 - the single site the COM rounds onto -
    is never emitted; the lemonade docstring and hierarchy.rst used to claim
    entry k was shell k, which would put the fall-off at index R and the
    interface one lattice unit too close to the centre.
    """
    dims = (24, 24, 24)
    centre = np.array([12, 12, 12])
    grid = np.indices((2 * radius + 1,) * 3).reshape(3, -1).T - radius
    positions = (grid + centre).astype(np.int32)

    traj = _make_traj([positions], ["A"] * len(positions), dims, hardwall=False)
    cluster = traj[0].clusters[0]
    assert cluster.n_beads == len(positions)          # one connected cube

    profile = cluster.radial_density_profile()
    assert len(profile) == min(dims) // 2 - 1
    assert profile[:radius] == pytest.approx([1.0] * radius)
    assert profile[radius] == pytest.approx(0.0)

    # the fix for this finding was to the DOCUMENTATION (the shipped
    # CLUSTER_RADIAL_DENSITY_PROFILE.dat format is shell-1-first and must not
    # change), so pin the docstring too - it used to say entry k was shell k
    doc = lemonade.Cluster.radial_density_profile.__doc__
    assert "shell 1" in doc and "k + 1" in doc


# ---------------------------------------------------------------------------
# F-I-6: the box cross-check must compare the DIMENSIONALITY too
# ---------------------------------------------------------------------------

def _keyfile_with_2d_dimensions(keyfile, destination):
    """Copy a keyfile, rewriting its DIMENSIONS line to drop the z axis."""
    lines = []
    for line in open(keyfile):
        if line.split(":")[0].strip() == "DIMENSIONS":
            box = line.split(":", 1)[1].split()
            line = "DIMENSIONS : " + " ".join(box[:2]) + "\n"
        lines.append(line)
    with open(destination, "w") as fh:
        fh.writelines(lines)
    return str(destination)


def _assert_discarded_and_bond_warnings(caught):
    """A 2D box on a 3D trajectory must give exactly the two expected warnings.

    Parameters
    ----------
    caught : pytest.WarningsRecorder
        The warnings recorded around one ``lemonade.load`` call.

    Returns
    -------
    None
    """
    messages = [str(w.message) for w in caught]
    assert len(messages) == 2, messages
    assert any("DISCARDED" in m for m in messages)
    assert any("not single lattice steps" in m for m in messages)


def test_two_dimensional_dimensions_on_a_three_dimensional_trajectory_warns(
        traj3d_files, traj2d_files, tmp_path):
    """Both the keyfile path and the ``dimensions=`` path must warn.

    The old cross-check truncated the trajectory's box record to the length of
    the resolved dimensions before comparing, so it never read the one column
    that distinguishes 2D from 3D and stayed silent while the z coordinate of
    every bead was discarded (which merges clusters and shrinks every Rg). The
    check keys on the RESOLVED dimensions, so threading it through the keyfile
    alone would have missed the argument path.
    """
    xtc, pdb, keyfile = traj3d_files
    flat_keyfile = _keyfile_with_2d_dimensions(keyfile, tmp_path / "FLAT.kf")

    # Each load warns twice, and both are wanted: the box cross-check (the z
    # coordinate is DISCARDED), and the bond check, which sees the consequence
    # (a bond along z collapses to length 0 once z is thrown away).
    with pytest.warns(UserWarning) as caught:
        by_keyfile = lemonade.load(xtc=xtc, pdb=pdb, keyfile=flat_keyfile)
    _assert_discarded_and_bond_warnings(caught)
    with pytest.warns(UserWarning) as caught:
        by_argument = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile,
                                    dimensions=(8, 8))
    _assert_discarded_and_bond_warnings(caught)
    assert by_keyfile.n_dim == 2 and by_argument.n_dim == 2

    # negative control: a real 2D run with its own 2D keyfile warns about nothing
    xtc2, pdb2, keyfile2 = traj2d_files
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        flat = lemonade.load(xtc=xtc2, pdb=pdb2, keyfile=keyfile2)
    assert flat.n_dim == 2
    assert [str(w.message) for w in caught] == []

    # and a 3D trajectory with its own 3D keyfile is likewise clean
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert [str(w.message) for w in caught] == []
