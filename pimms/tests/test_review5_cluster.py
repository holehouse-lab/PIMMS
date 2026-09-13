"""
Regression tests for the fifth review's cluster-output findings.

Three defects, one test each:

* **F-H-1** - the long-range single-image gather walked a plain "anything within
  three sites" rule instead of the relation that DEFINES a long-range cluster, so
  an elongated cluster whose extent came within three sites of the box was linked
  to itself through the periodic face and torn apart. Every ``LR_CLUSTER_*``
  quantity then disagreed with the ``CLUSTER_*`` quantity for the very same set of
  chains.
* **F-H-2** - under periodic boundaries the radial density profile normalised by
  the periodic shell site count but binned beads by their RAW single-image
  distance, so a bead whose periodic image occupies a shell site was dropped.
* **F-H-3** - a cluster connected to its own periodic image has no shape, hull or
  radial profile, yet plausible numbers were written for it in every frame and the
  only diagnostic was a ``UserWarning`` the interpreter shows once per process.

Every expected value here comes from a plain-numpy oracle written from the
definitions (a BFS over the defining relation, an explicit enumeration of shell
sites reduced mod the box), never from calling the code under test a second way.
"""

import contextlib
import itertools
import math
import os
import random
import warnings
from collections import deque

import numpy as np
import pytest

from scipy.spatial import ConvexHull

from pimms import lattice_analysis_utils as lau
from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _restore_cwd_and_rng():
    """Tests here chdir into their tmp_path and the simulation reseeds the global
    RNG; put both back so later tests do not depend on this module having run."""
    cwd = os.getcwd()
    state = random.getstate()
    yield
    random.setstate(state)
    os.chdir(cwd)


def _inject(lattice, box, positions_by_chain):
    """Drop a hand-built configuration onto the lattice."""
    grid = np.zeros(tuple(box), dtype=lattice.grid.dtype)
    type_grid = np.zeros(tuple(box), dtype=lattice.type_grid.dtype)
    chain_dict = {}
    for cid, positions in positions_by_chain.items():
        chain = lattice.chains[cid]
        assert len(positions) == len(chain), (cid, len(positions), len(chain))
        wrapped = []
        for k, p in enumerate(positions):
            p = [int(v) % int(box[d]) for d, v in enumerate(p)]
            grid[tuple(p)] = cid
            type_grid[tuple(p)] = int(chain.int_sequence[k])
            wrapped.append(p)
        chain_dict[cid] = wrapped
    lattice.lattice_restorefrombackup(grid, type_grid, chain_dict)
    return chain_dict


def _rows(path):
    """Read a comma-separated PIMMS .dat file into a list of string fields."""
    with open(path) as fh:
        return [[cell.strip() for cell in line.strip().rstrip(',').split(',')]
                for line in fh if line.strip()]


# --- oracles ---------------------------------------------------------------

def _min_image(delta, box):
    d = np.asarray(delta, dtype=np.int64)
    L = np.asarray(box, dtype=np.int64)
    return (d + L // 2) % L - L // 2


def _gather_over_contacts(chains, box):
    """Gather a cluster into one image by walking bonds and Chebyshev-1 contacts.

    This is the oracle for F-H-1: it walks ONLY the relation that joins the
    chains (consecutive beads of a chain, and any two beads at minimum-image
    Chebyshev distance 1), placing each newly reached bead at its parent's
    position plus the minimum-image displacement of the edge that reached it.

    Parameters
    ----------
    chains : list of numpy.ndarray
        One ``(n_beads, n_dim)`` array of wrapped positions per chain, beads in
        bonded order.

    box : list of int
        The box dimensions.

    Returns
    -------
    numpy.ndarray
        ``(N, n_dim)`` array of gathered positions, in chain-concatenation order.

    """
    pos = np.concatenate([np.asarray(c, dtype=np.int64) for c in chains])
    n = len(pos)

    bonds = set()
    offset = 0
    for c in chains:
        for k in range(len(c) - 1):
            bonds.add((offset + k, offset + k + 1))
            bonds.add((offset + k + 1, offset + k))
        offset += len(c)

    adjacency = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = _min_image(pos[j] - pos[i], box)
            if (i, j) in bonds or int(np.abs(d).max()) == 1:
                adjacency[i].append((j, d))
                adjacency[j].append((i, -d))

    si = np.zeros_like(pos)
    seen = np.zeros(n, dtype=bool)
    seen[0] = True
    si[0] = pos[0]
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j, d in adjacency[i]:
            if not seen[j]:
                seen[j] = True
                si[j] = si[i] + d
                queue.append(j)
    assert seen.all(), "oracle: the chains are not one contact-connected cluster"
    return si


def _rg(points):
    p = np.asarray(points, dtype=float)
    return float(np.sqrt(((p - p.mean(axis=0)) ** 2).sum(axis=1).mean()))


def _asphericity_3d(points):
    p = np.asarray(points, dtype=float)
    c = p - p.mean(axis=0)
    eig = np.sort(np.linalg.eigvalsh(c.T @ c / len(c)))
    total = eig.sum()
    return float(1 - 3 * (eig[0] * eig[1] + eig[1] * eig[2] + eig[2] * eig[0]) / total ** 2)


def _radial_profile_by_site_enumeration(points, box, hardwall=False):
    """Radial density profile from an explicit enumeration of the shell sites.

    Shell ``k`` is the set of sites at Chebyshev distance ``k`` from the rounded
    centre of mass. Under periodic boundaries those sites are ``(COM + v) mod L``
    for every offset ``v`` with ``|v|_inf == k``, and a bead occupies the site
    ``p mod L``; under a hardwall no wrapping happens and sites outside the box do
    not exist. The occupancy is the fraction of the shell's sites that hold a bead.

    Parameters
    ----------
    points : numpy.ndarray
        ``(N, n_dim)`` array of the cluster's (single-image) positions.

    box : list of int
        The box dimensions.

    hardwall : bool, optional
        True for a hardwall box. Default is False.

    Returns
    -------
    list of float
        One occupancy per shell, out to ``min(box) // 2 - 1`` shells.

    """
    pts = np.asarray(points, dtype=np.int64)
    n_dim = pts.shape[1]
    com = np.rint(pts.mean(axis=0)).astype(int)
    k_max = int(min(box) / 2) - 1

    if hardwall:
        occupied_sites = {tuple(int(v) for v in p) for p in pts}
    else:
        occupied_sites = {tuple(int(p[d]) % int(box[d]) for d in range(n_dim)) for p in pts}

    profile = []
    for k in range(1, k_max + 1):
        shell = set()
        for offsets in itertools.product(range(-k, k + 1), repeat=n_dim):
            if max(abs(o) for o in offsets) != k:
                continue
            site = tuple(int(com[d]) + offsets[d] for d in range(n_dim))
            if hardwall:
                if any(site[d] < 0 or site[d] >= int(box[d]) for d in range(n_dim)):
                    continue
            else:
                site = tuple(site[d] % int(box[d]) for d in range(n_dim))
            shell.add(site)

        if not hardwall:
            # the shell cap is what licenses the periodic site count: shells must
            # not fold back onto themselves, or the denominator would be wrong
            assert len(shell) == (2 * k + 1) ** n_dim - (2 * k - 1) ** n_dim

        profile.append(len(shell & occupied_sites) / len(shell))

    return profile


# ---------------------------------------------------------------------------
# F-H-1: the long-range gather must walk the relation that defines the cluster
# ---------------------------------------------------------------------------

def _blob_and_tail(tail_length):
    """A 3x3x3 blob of nine 3-mers at x = 0..2 plus a tail chain along +x."""
    def build(chain_ids):
        by_chain = {}
        i = 0
        for y in (4, 5, 6):
            for z in (4, 5, 6):
                by_chain[chain_ids[i]] = [[x, y, z] for x in range(3)]
                i += 1
        by_chain[chain_ids[i]] = [[x, 5, 5] for x in range(3, 3 + tail_length)]
        return by_chain
    return build


def test_lr_cluster_shape_matches_sr_when_the_memberships_are_identical(tmp_path):
    """With a short-range-only parameter file the LR and SLR tables are all zero, so
    a long-range cluster IS the short-range cluster - the two file families describe
    the same beads in the same configuration and must agree row for row, and must
    agree with a gather over the bond+contact connectivity.

    The shape matters: a compact blob with a protruding tail tears where a straight
    rod may not. Before the fix the extent-11 blob+tail gave LR_CLUSTER_RG 2.3516
    against CLUSTER_RG 2.8420 (the truth), a fictitious solid first radial shell
    (1.0000 vs 0.6923) and a spurious percolation warning.
    """
    box = [12, 12, 12]
    tail = 8                                    # single-image extent 11 = L - 1
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(9, "AAA"), (1, "A" * tail)],
                          extra={"ANA_CLUSTER_THRESHOLD": 1})
    ham = state.ham

    # the identity this test rests on: no shell other than contact can connect
    assert not np.any(ham.LR_residue_interaction_table)
    assert not np.any(ham.SLR_residue_interaction_table)

    chain_ids = sorted(state.lattice.chains)
    injected = _inject(state.lattice, box, _blob_and_tail(tail)(chain_ids))

    # independent oracle: gather over bonds + Chebyshev-1 contacts only
    oracle_si = _gather_over_contacts([np.asarray(injected[c]) for c in chain_ids], box)
    oracle_hull = ConvexHull(oracle_si)
    oracle_radial = _radial_profile_by_site_enumeration(oracle_si, box)

    os.chdir(tmp_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            state.sim.ANAFUNCT_cluster_analysis(100)

    assert [w for w in caught if "percolat" in str(w.message)] == [], \
        "a cluster that does not touch its own image must not be called percolating"

    # (1) identical membership
    assert _rows("CLUSTERS.dat")[-1] == _rows("LR_CLUSTERS.dat")[-1] == ["100", "10"]

    # (2) identical shape descriptors, and both equal to the oracle
    expected = {
        "CLUSTER_RG.dat": _rg(oracle_si),
        "CLUSTER_ASPH.dat": _asphericity_3d(oracle_si),
        "CLUSTER_VOL.dat": oracle_hull.volume,
        "CLUSTER_AREA.dat": oracle_hull.area,
        "CLUSTER_DEN.dat": len(oracle_si) / oracle_hull.volume,
    }
    for name, truth in expected.items():
        sr = float(_rows(name)[-1][1])
        lr = float(_rows("LR_" + name)[-1][1])
        assert sr == pytest.approx(lr, abs=1e-12), name
        assert sr == pytest.approx(round(truth, 4), abs=5e-5), name

    # (3) identical radial profiles, and both equal to the site enumeration
    sr_profile = [float(v) for v in _rows("CLUSTER_RADIAL_DENSITY_PROFILE.dat")[-1][2:]]
    lr_profile = [float(v) for v in _rows("LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat")[-1][2:]]
    assert sr_profile == lr_profile
    assert sr_profile == [pytest.approx(round(v, 4), abs=5e-5) for v in oracle_radial]


def test_lr_cluster_winding_through_a_real_interaction_is_still_percolation(tmp_path):
    """The mirror of the test above, so the fix cannot be satisfied by simply
    lowering the long-range gather to contact distance: an LR-capable rod whose two
    ends sit at Chebyshev 2 through the face WITH a nonzero LR entry genuinely
    interacts with its own image, so the LR cluster has no single image and its
    shape quantities must be blanked - while the SR cluster (whose contacts do not
    wind) keeps a real number for the very same beads.
    """
    box = [12, 12, 12]
    state = U.build_state(tmp_path, 3, "LR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(1, "AAAAAA"), (1, "AAAAA")],
                          extra={"ANA_CLUSTER_THRESHOLD": 1})
    assert np.any(state.ham.LR_residue_interaction_table)

    chain_ids = sorted(state.lattice.chains)
    _inject(state.lattice, box, {chain_ids[0]: [[x, 5, 5] for x in range(6)],
                                 chain_ids[1]: [[x, 5, 5] for x in range(6, 11)]})

    os.chdir(tmp_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            state.sim.ANAFUNCT_cluster_analysis(100)

    assert [w for w in caught if "percolat" in str(w.message)], \
        "a cluster interacting with its own image through the face must be flagged"

    # the 11-bead rod: a straight line, so its hull is degenerate either way, but Rg
    # is well defined for the contact cluster and undefined for the LR cluster
    assert float(_rows("CLUSTER_RG.dat")[-1][1]) == pytest.approx(math.sqrt(10), abs=5e-5)
    assert math.isnan(float(_rows("LR_CLUSTER_RG.dat")[-1][1]))
    assert math.isnan(float(_rows("LR_CLUSTER_ASPH.dat")[-1][1]))


# ---------------------------------------------------------------------------
# F-H-2: the radial profile must bin beads with the metric it normalises with
# ---------------------------------------------------------------------------

def _blob_plus_tail_positions(tail_to_x):
    """A 3x3x3 blob at the origin plus a straight tail along +x out to tail_to_x."""
    pts = [[x, y, z] for x in range(3) for y in range(3) for z in range(3)]
    pts += [[x, 1, 1] for x in range(3, tail_to_x + 1)]
    return np.array(pts, dtype=np.int64)


@pytest.mark.parametrize("box,tail_to_x", [
    ([12, 12, 12], 10),      # shells 4 and 5 pick up the wrapped tail beads
    ([12, 12, 12], 9),
    ([11, 11, 11], 9),       # odd box
    ([12, 16, 16], 9),       # non-cubic, tail along the SHORT axis (sets the cap)
    ([12, 12, 12], 6),       # control: nothing wraps into a shell
    ([12, 16, 16], 13),      # control: tail along a LONG axis, images miss the shells
])
def test_periodic_radial_profile_counts_beads_whose_image_occupies_a_shell_site(box, tail_to_x):
    """Shell k is a set of SITES, and the code already normalises by the periodic
    site count, so the numerator has to be counted mod the box too. A bead sitting
    L-4 past the centre of mass in its single image occupies the shell-4 site at
    COM-4; binning it by its raw single-image distance dropped it from every shell
    and the outer shells read too dilute. For the 12-box blob+tail to x=10 the old
    code wrote shells 4 and 5 as 0.0026 and 0.0017 where the site enumeration gives
    0.0052 and 0.0033 (one bead counted instead of two, in each).
    """
    if tail_to_x == 13:
        pts = np.array([[x, y, z] for x in range(3) for y in range(3) for z in range(3)]
                       + [[1, y, 1] for y in range(3, 14)], dtype=np.int64)
    else:
        pts = _blob_plus_tail_positions(tail_to_x)

    profile = lau.compute_cluster_radial_density_profile([pts], box)[0]
    expected = _radial_profile_by_site_enumeration(pts, box)

    assert profile == [pytest.approx(v) for v in expected]


def test_hardwall_radial_profile_is_unchanged_by_the_periodic_fix():
    """A hardwall box does not wrap, so its shells are plain Cartesian cubes clipped
    at the walls and its beads must keep being binned by plain Cartesian distance.
    The minimum-image fix must not leak into this branch (applying the periodic
    metric to hardwall trajectories is a bug that has been fixed once already)."""
    box = [12, 12, 12]
    pts = _blob_plus_tail_positions(10)

    profile = lau.compute_cluster_radial_density_profile([pts], box, hardwall=True)[0]
    expected = _radial_profile_by_site_enumeration(pts, box, hardwall=True)

    assert profile == [pytest.approx(v) for v in expected]
    # and the hardwall answer is genuinely different from the periodic one, so this
    # is not passing by accident
    assert profile != lau.compute_cluster_radial_density_profile([pts], box)[0]


def test_symmetric_cluster_profile_is_untouched_by_the_minimum_image_fix():
    """The trigger is asymmetry about the centre of mass, not extent: a dumbbell of
    extent 11 in a 12-box has every bead within offset_max of its COM, so nothing
    wraps into a shell and the profile is the same before and after the fix."""
    box = [12, 12, 12]
    pts = np.array([[x, 6, 6] for x in range(0, 3)] + [[x, 6, 6] for x in range(8, 11)]
                   + [[x, 6, 6] for x in range(3, 8)], dtype=np.int64)

    profile = lau.compute_cluster_radial_density_profile([pts], box)[0]
    expected = _radial_profile_by_site_enumeration(pts, box)

    assert profile == [pytest.approx(v) for v in expected]


# ---------------------------------------------------------------------------
# F-H-3: a percolating cluster must be marked in the files and in the log
# ---------------------------------------------------------------------------

def test_percolating_cluster_is_sentinelled_in_every_row_and_logged_every_step(tmp_path):
    """Drive a real run whose largest cluster winds the box: a 3x3 tube of nine
    frozen decamers spanning x in a 10-box (90 beads, so it also clears the
    27-bead radial threshold), a compact non-percolating control cluster
    straddling a face, and one mobile monomer (a freeze file may not name every
    chain).

    Before the fix every row of CLUSTER_RG.dat carried a plausible number for the
    tube (3.0957, a property of the box rather than of the condensate, and
    BFS-order dependent), a radial profile row was written for it, and log.txt
    contained zero lines about it - the warning reached stderr twice for the whole
    run. Now the tube's columns carry nan / -1, its radial row is gone, and the log
    records every analysis step at which it happened.
    """
    box = [10, 10, 10]
    n_steps = 6
    equilibration = 1

    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    tube_chains = list(range(1, 10))
    with open(tmp_path / "freeze.txt", "w") as fh:
        fh.write("C " + " ".join(str(i) for i in tube_chains + [10, 11]) + "\n")

    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=box, chains=[(9, "A" * 10), (2, "AAAA"), (1, "A")],
                    seed=7, n_steps=n_steps, equilibration=equilibration, temperature=40,
                    extra={"ANA_CLUSTER": 1, "ANALYSIS_FREQ": 1, "ANA_CLUSTER_THRESHOLD": 1,
                           "PRINT_FREQ": 1000000, "XTC_FREQ": 1000000, "RESTART_FREQ": 1000000,
                           "EN_FREQ": 1000000, "ANA_POL": 1000000,
                           "FREEZE_FILE": "freeze.txt"})

    os.chdir(tmp_path)
    keyfile = KeyFileParser("KEYFILE.kf")
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        sim = Simulation(keyfile.keyword_lookup)

    chain_ids = sorted(sim.LATTICE.chains)
    layout = {}
    for cid, (y, z) in zip(chain_ids[:9], [(y, z) for y in (4, 5, 6) for z in (4, 5, 6)]):
        layout[cid] = [[x, y, z] for x in range(box[0])]           # winds x
    # control: a 2x2x2 cube straddling the x = 0 face, contact-connected, not winding
    layout[chain_ids[9]] = [[9, 0, 0], [9, 1, 0], [9, 1, 1], [9, 0, 1]]
    layout[chain_ids[10]] = [[0, 0, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]]
    layout[chain_ids[11]] = [[5, 0, 8]]                            # the mobile monomer
    _inject(sim.LATTICE, box, layout)

    # the tube winds the box on purpose, so the gather's documented UserWarning
    # is expected here and is asserted rather than left to leak into the test
    # output; the per-step LOG record is what the assertions below check
    with contextlib.redirect_stdout(open(os.devnull, "w")), \
            pytest.warns(UserWarning, match="percolates the periodic box"):
        sim.run_simulation()

    n_analysis = n_steps - equilibration

    # clusters are size sorted, so column 1 is the 9-chain tube and column 2 the
    # control (the mobile monomer may join the control, but can never outrank the tube)
    for name in ("CLUSTER_RG.dat", "CLUSTER_ASPH.dat",
                 "LR_CLUSTER_RG.dat", "LR_CLUSTER_ASPH.dat"):
        rows = _rows(name)
        assert len(rows) == n_analysis
        assert all(math.isnan(float(row[1])) for row in rows), name

    for name in ("CLUSTER_VOL.dat", "CLUSTER_AREA.dat", "CLUSTER_DEN.dat",
                 "LR_CLUSTER_VOL.dat", "LR_CLUSTER_AREA.dat", "LR_CLUSTER_DEN.dat"):
        rows = _rows(name)
        assert all(float(row[1]) == -1 for row in rows), name

    # the control cluster - which straddles a face but does not wind - keeps a real
    # measurement in exactly the same rows, so the sentinel is targeted rather than
    # a blanket blanking of the step (its exact value is pinned in the test below)
    for row in _rows("CLUSTER_RG.dat"):
        assert float(row[2]) > 0
    for row in _rows("CLUSTER_VOL.dat"):
        assert float(row[2]) > 0

    # the tube is the only cluster over the 27-bead radial threshold, and its row is
    # gone, so there is nothing to write and the file is never created at all
    # (output files come into existence on their first row)
    assert not os.path.exists("CLUSTER_RADIAL_DENSITY_PROFILE.dat")
    assert not os.path.exists("LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat")

    # the size distribution still reports the percolating component - that is how a
    # user diagnoses it in the first place
    # (9 chains, or 10 if the mobile monomer happens to have wandered onto the tube)
    assert all(int(row[1]) >= 9 for row in _rows("CLUSTERS.dat"))

    # EVERY occurrence reaches the run's own record. The count is the load-bearing
    # assertion: ">= 1" passed before the fix too, once the warning was re-armed.
    with open("log.txt") as fh:
        logged = [line for line in fh if "percolate" in line]
    assert len(logged) == n_analysis
    assert all("Rg/asphericity written as nan" in line for line in logged)


def test_only_the_percolating_cluster_is_blanked_in_a_frame(tmp_path):
    """The deterministic companion to the run above: one analysis frame holding a
    tube that winds the box AND a compact cluster that merely straddles a face.
    The winding cluster is blanked; the straddling one keeps values equal to a
    plain-numpy / ConvexHull oracle on a gather over its own contacts."""
    box = [10, 10, 10]
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(9, "A" * 10), (2, "AAAA")],
                          extra={"ANA_CLUSTER_THRESHOLD": 1})

    chain_ids = sorted(state.lattice.chains)
    layout = {}
    for cid, (y, z) in zip(chain_ids[:9], [(y, z) for y in (4, 5, 6) for z in (4, 5, 6)]):
        layout[cid] = [[x, y, z] for x in range(box[0])]
    layout[chain_ids[9]] = [[9, 0, 0], [9, 1, 0], [9, 1, 1], [9, 0, 1]]
    layout[chain_ids[10]] = [[0, 0, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]]
    _inject(state.lattice, box, layout)

    control_si = _gather_over_contacts([np.asarray(layout[chain_ids[9]]),
                                        np.asarray(layout[chain_ids[10]])], box)
    control_hull = ConvexHull(control_si)

    os.chdir(tmp_path)
    # one of the two clusters winds the box on purpose, so the gather's documented
    # UserWarning is expected and asserted here
    with contextlib.redirect_stdout(open(os.devnull, "w")), \
            pytest.warns(UserWarning, match="percolates the periodic box"):
        state.sim.ANAFUNCT_cluster_analysis(7)

    assert _rows("CLUSTERS.dat")[-1] == ["7", "9", "2"]

    assert math.isnan(float(_rows("CLUSTER_RG.dat")[-1][1]))
    assert float(_rows("CLUSTER_VOL.dat")[-1][1]) == -1

    assert float(_rows("CLUSTER_RG.dat")[-1][2]) == \
        pytest.approx(round(_rg(control_si), 4), abs=5e-5)
    assert float(_rows("CLUSTER_ASPH.dat")[-1][2]) == \
        pytest.approx(round(_asphericity_3d(control_si), 4), abs=5e-5)
    assert float(_rows("CLUSTER_VOL.dat")[-1][2]) == \
        pytest.approx(round(control_hull.volume, 4), abs=5e-5)
    assert float(_rows("CLUSTER_AREA.dat")[-1][2]) == \
        pytest.approx(round(control_hull.area, 4), abs=5e-5)
