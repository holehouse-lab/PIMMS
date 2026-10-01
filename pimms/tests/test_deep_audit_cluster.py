"""
Regression tests for the deep audit's in-simulation cluster-analysis findings.

* **A1-1** - the percolation test only recognised a closing pair of beads one
  box length apart in the single image, so a cluster whose every closing contact
  winds the box twice or more (a helix in a box with a short axis) was written
  out with ordinary-looking shape values.
* **A1-2 / P2-1** - the same test built an (N_high x N_low) separation matrix,
  quadratic in the number of beads for any cluster that spans the box, and the
  simulation ran it twice per cluster.
* **A1-3** - the radial-profile centre was rounded half-to-even, so translating
  the configuration by one site changed the profile.
* **A1-4 / PERF1** - the long-range cluster search rebuilt a box-sized LR-flag
  grid on every call although nothing reads it when both tables are passed.
* **P2-3 / S1-6** - every chain allocated a seqlen x seqlen distance-map
  accumulator at construction, and nothing warned about the O(seqlen^2) cost.

The percolation oracles are written from the definition and share no code with
PIMMS: one unfolds the cluster over periodic images by flood fill and asks
whether a site is reached at two different images; the other replicates the box
along one axis and counts the connected components of the copies.
"""

import contextlib
import itertools
import os
import random
import tracemalloc

import numpy as np
import pytest

from pimms import analysis_structures
from pimms import cluster_utils
from pimms import lattice_analysis_utils as lau
from pimms import lattice_utils
from pimms.tests import kernel_test_utils as U

_RANDOM_CASES = {
    # mode: (contact distance, boxes with occupancies near the percolation threshold)
    "SR": (1, [((5, 7), .42), ((7, 25), .42), ((12, 12), .40), ((3, 3, 3), .2), ((4, 5, 6), .15),
               ((7, 12, 12), .11), ((7, 16, 16), .11), ((9, 9, 9), .11)]),
    "untyped3": (3, [((7, 7), .08), ((8, 11), .07), ((7, 30), .07), ((7, 7, 7), .012), ((7, 8, 12), .012)]),
    "LR": (3, [((7, 7), .13), ((8, 11), .12), ((7, 30), .12), ((7, 7, 7), .020), ((7, 8, 12), .020)]),
    "SLR": (3, [((7, 7), .09), ((8, 11), .08), ((7, 30), .08), ((7, 7, 7), .014), ((7, 8, 12), .014)]),
}

# One of the clusters the audit found in a random-occupancy 7 x 12 x 12 box
# (71 beads): every contact that closes it joins beads two box lengths apart.
_MULTI_WINDING_CLUSTER = [
    [0, 0, 9], [6, 1, 10], [0, 1, 9], [1, 11, 8], [1, 0, 9], [1, 10, 7], [1, 11, 7], [2, 10, 9],
    [1, 9, 9], [0, 9, 6], [1, 9, 7], [0, 8, 7], [2, 8, 7], [2, 7, 6], [2, 6, 7], [3, 6, 5],
    [3, 5, 8], [3, 4, 9], [4, 5, 9], [4, 6, 9], [4, 7, 8], [5, 6, 8], [6, 6, 7], [6, 6, 9],
    [6, 5, 10], [0, 6, 8], [1, 7, 9], [5, 4, 11], [5, 4, 0], [6, 3, 1], [0, 4, 0], [0, 4, 2],
    [1, 5, 11], [2, 5, 11], [3, 4, 11], [2, 3, 0], [3, 2, 11], [3, 4, 1], [2, 1, 0], [4, 1, 11],
    [5, 0, 11], [2, 0, 11], [2, 0, 0], [1, 11, 1], [0, 10, 1], [6, 10, 1], [1, 9, 2], [0, 8, 1],
    [1, 8, 0], [2, 8, 0], [2, 9, 11], [3, 8, 1], [4, 7, 0], [4, 7, 11], [5, 9, 0], [4, 10, 11],
    [5, 9, 11], [6, 8, 5], [6, 2, 8], [0, 2, 8], [0, 2, 9], [6, 3, 9], [0, 3, 8], [1, 4, 7],
    [1, 3, 6], [2, 4, 6], [3, 3, 5], [2, 4, 4], [4, 3, 6], [1, 5, 4], [5, 2, 9]]


@pytest.fixture(autouse=True)
def _restore_cwd_and_rng():
    """Tests here chdir into their tmp_path and a Simulation reseeds the global
    RNG; put both back so later tests do not depend on this module having run."""
    cwd = os.getcwd()
    state = random.getstate()
    yield
    random.setstate(state)
    os.chdir(cwd)


# ---------------------------------------------------------------------------
# oracles (from the definition; no PIMMS code)
# ---------------------------------------------------------------------------

def _linked(offset, type_a, type_b, LR, SLR, t):
    """Are two beads separated by this lattice offset linked?

    Parameters
    ----------
    offset : tuple of int
        The lattice offset from one bead to the other.
    type_a, type_b : int
        Residue codes of the two beads (ignored when ``LR`` is None).
    LR, SLR : list of list of int or None
        The long-range / super-long-range tables, or None for the plain
        distance rule.
    t : int
        The contact distance.

    Returns
    -------
    bool
        True if the pair is linked: within Chebyshev ``t`` for the distance
        rule; a contact, a Chebyshev-2 pair with a nonzero LR entry or a
        Chebyshev-3 pair with a nonzero SLR entry for the typed rule.
    """
    cheb = max(abs(x) for x in offset)
    if cheb == 0 or cheb > t:
        return False
    if LR is None or cheb == 1:
        return True
    if cheb == 2:
        return LR[type_a][type_b] != 0
    if cheb == 3:
        return SLR[type_a][type_b] != 0
    return False


def _components_with_winding(site_types, dims, t=1, LR=None, SLR=None):
    """Connected components of a periodic configuration and the axes each winds.

    Flood fill over periodic IMAGES: every site is given the image (a vector of
    box counts) it is first reached in, and a link that reaches an already
    placed site in a different image shows the component is connected to its
    own image along every axis on which the two images differ.

    Parameters
    ----------
    site_types : dict
        Occupied site (tuple of int) -> residue code.
    dims : tuple of int
        The box.
    t : int, optional
        The contact distance.
    LR, SLR : list of list of int or None, optional
        The typed-rule tables, or None for the distance rule.

    Returns
    -------
    list of tuple
        ``(sites, winding_axes)`` per component: the sorted site list and the
        set of axes along which the component winds the box.
    """
    n_dim = len(dims)
    offsets = [o for o in itertools.product(range(-t, t + 1), repeat=n_dim) if any(o)]
    seen = set()
    out = []
    for start in sorted(site_types):
        if start in seen:
            continue
        image = {start: (0,) * n_dim}
        seen.add(start)
        stack = [start]
        wind = set()
        while stack:
            p = stack.pop()
            for o in offsets:
                raw = tuple(p[e] + o[e] for e in range(n_dim))
                q = tuple(raw[e] % dims[e] for e in range(n_dim))
                if q not in site_types or not _linked(o, site_types[p], site_types[q], LR, SLR, t):
                    continue
                target = tuple(image[p][e] + (raw[e] - q[e]) // dims[e] for e in range(n_dim))
                if q not in image:
                    image[q] = target
                    seen.add(q)
                    stack.append(q)
                else:
                    wind.update(e for e in range(n_dim) if image[q][e] != target[e])
        out.append((sorted(image), wind))
    return out


def _supercell_winds(sites, site_types, dims, axis, t=1, LR=None, SLR=None):
    """Second, independent oracle: does the cluster join its copies along ``axis``?

    The box is replicated ``m`` times along one axis and the components of the
    replicated cluster are counted. ``m`` separate copies means the cluster is
    not connected to its image along that axis. Links reach at most ``t``
    sites, so a cluster of N <= L_axis * (other sites) beads cannot need more
    than ``t * other`` box lengths to close; ``m = t * other + 1`` copies can
    therefore never be fooled by a winding number that is a multiple of ``m``.

    Parameters
    ----------
    sites : list of tuple
        The sites of one component.
    site_types : dict
        Occupied site -> residue code.
    dims : tuple of int
        The box.
    axis : int
        The axis to replicate along.
    t : int, optional
        The contact distance.
    LR, SLR : list of list of int or None, optional
        The typed-rule tables, or None for the distance rule.

    Returns
    -------
    bool
        True if the copies are joined, i.e. the cluster winds ``axis``.
    """
    n_dim = len(dims)
    other = 1
    for e in range(n_dim):
        if e != axis:
            other *= dims[e]
    m = t * other + 1
    big = list(dims)
    big[axis] = dims[axis] * m
    occupied = {}
    for s in sites:
        for k in range(m):
            q = list(s)
            q[axis] += k * dims[axis]
            occupied[tuple(q)] = site_types[s]
    offsets = [o for o in itertools.product(range(-t, t + 1), repeat=n_dim) if any(o)]
    seen = set()
    n_components = 0
    for start in occupied:
        if start in seen:
            continue
        n_components += 1
        seen.add(start)
        stack = [start]
        while stack:
            p = stack.pop()
            for o in offsets:
                q = tuple((p[e] + o[e]) % big[e] for e in range(n_dim))
                if q in occupied and q not in seen and _linked(o, occupied[p], occupied[q], LR, SLR, t):
                    seen.add(q)
                    stack.append(q)
    return n_components < m


def _pimms_axes(sites, site_types, dims, t=1, LR=None, SLR=None):
    """Gather one component with PIMMS and return the axes PIMMS says it winds."""
    pos = np.array(sites, dtype=np.int64)
    if LR is None:
        types = LR_arr = SLR_arr = None
    else:
        types = np.array([site_types[s] for s in sites], dtype=np.int64)
        LR_arr, SLR_arr = np.array(LR, dtype=np.int64), np.array(SLR, dtype=np.int64)
    si = cluster_utils.convert_positions_to_single_image_snakesearch(
        pos, list(dims), t, types=types, LR_table=LR_arr, SLR_table=SLR_arr,
        warn_if_percolating=False)
    axes = cluster_utils.percolating_axes(si, list(dims), t, types, LR_arr, SLR_arr)
    first = cluster_utils.cluster_percolates(si, list(dims), t, types, LR_arr, SLR_arr)
    assert first == (axes[0] if axes else None)
    return set(axes)


# ---------------------------------------------------------------------------
# A1-1: the exact criterion against the oracles
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dims, min_compact", [((3, 4), 500), ((2, 2, 3), 20)])
def test_percolation_matches_the_image_flood_fill_on_every_configuration_of_a_small_box(
        dims, min_compact):
    """Every one of the 4096 occupancy patterns of a 12-site box, every
    component: the axes PIMMS reports are the axes the flood fill over periodic
    images finds, and (on every 16th pattern) the axes the supercell count
    finds. The 2 x 2 x 3 box has axes on which a bead touches its neighbour
    through both faces, so nearly everything in it percolates; the 3 x 4 box
    has both answers in bulk."""
    all_sites = list(itertools.product(*[range(d) for d in dims]))
    n_percolating = n_clusters = 0
    for pattern in range(1, 2 ** len(all_sites)):
        site_types = {s: 0 for k, s in enumerate(all_sites) if (pattern >> k) & 1}
        for sites, wind in _components_with_winding(site_types, dims):
            n_clusters += 1
            n_percolating += bool(wind)
            assert _pimms_axes(sites, site_types, dims) == wind, (dims, sites)
            if pattern % 16 == 0:
                supercell = {e for e in range(len(dims))
                             if _supercell_winds(sites, site_types, dims, e)}
                assert supercell == wind, (dims, sites)
    # the enumeration exercises both answers heavily
    assert n_percolating > 500 and n_clusters - n_percolating >= min_compact


@pytest.mark.parametrize("mode", sorted(_RANDOM_CASES))
def test_percolation_matches_the_oracles_on_random_clusters(mode):
    """Random-occupancy boxes in 2D and 3D, cubic and with a short axis, under
    contact connectivity, the plain distance-3 rule, and the typed LR and
    LR+SLR rules (random symmetric tables with zeros, three residue types)."""
    t, boxes = _RANDOM_CASES[mode]
    rng = np.random.default_rng({"SR": 101, "untyped3": 102, "LR": 103, "SLR": 104}[mode])
    n_percolating = n_clusters = n_supercell = 0
    for it in range(260):
        dims, rho = boxes[it % len(boxes)]
        LR = SLR = None
        if mode in ("LR", "SLR"):
            a = rng.integers(0, 2, (3, 3)) * -3
            LR = (np.triu(a) + np.triu(a, 1).T).tolist()
            SLR = [[0] * 3 for _ in range(3)]
            if mode == "SLR":
                a = rng.integers(0, 2, (3, 3)) * 2
                SLR = (np.triu(a) + np.triu(a, 1).T).tolist()
        occ = rng.random(dims) < rho * rng.uniform(0.7, 1.3)
        site_types = {tuple(int(x) for x in s): int(rng.integers(3)) for s in np.argwhere(occ)}
        for sites, wind in _components_with_winding(site_types, dims, t, LR, SLR):
            if len(sites) < 2:
                continue
            n_clusters += 1
            n_percolating += bool(wind)
            assert _pimms_axes(sites, site_types, dims, t, LR, SLR) == wind, (mode, dims, sites)
            if int(np.prod(dims)) <= 130 and len(sites) <= 40:
                n_supercell += 1
                supercell = {e for e in range(len(dims))
                             if _supercell_winds(sites, site_types, dims, e, t, LR, SLR)}
                assert supercell == wind, (mode, dims, sites)
    assert n_percolating >= 25 and n_clusters - n_percolating >= 100, (n_percolating, n_clusters)
    assert n_supercell >= 20


def test_a_rod_one_bead_short_of_the_box_spans_it_without_percolating():
    """A straight rod of L-1 beads reaches from one face to within a site of the
    other and does not touch its image; the rod of L beads does."""
    for dims in ([9, 11], [8, 9, 10]):
        for axis, L in enumerate(dims):
            def rod(n):
                pts = np.full((n, len(dims)), 3, dtype=np.int64)
                pts[:, axis] = (np.arange(n) + 5) % L      # straddles the face
                return [tuple(int(v) for v in p) for p in pts]
            short, full = rod(L - 1), rod(L)
            assert _pimms_axes(short, dict.fromkeys(short, 0), dims) == set()
            assert _pimms_axes(full, dict.fromkeys(full, 0), dims) == {axis}


def test_diagonal_windings_list_every_axis_the_diagonal_advances_along():
    face = [(i, i) for i in range(8)]
    assert _pimms_axes(face, dict.fromkeys(face, 0), (8, 8)) == {0, 1}
    body = [(i, i, i) for i in range(7)]
    assert _pimms_axes(body, dict.fromkeys(body, 0), (7, 7, 7)) == {0, 1, 2}
    in_plane = [(i, i, 4) for i in range(9)]
    assert _pimms_axes(in_plane, dict.fromkeys(in_plane, 0), (9, 9, 9)) == {0, 1}


def test_clusters_that_only_close_after_winding_the_box_twice_are_flagged():
    """The two shapes the one-box-length pair search could not see."""
    # a staircase that advances two sites in x per site in y: it meets itself
    # through the corner of a 10-box after going round x twice and y once
    stairs = [(i % 10, i // 2, 0) for i in range(20)]
    types = dict.fromkeys(stairs, 0)
    (sites, wind), = _components_with_winding(types, (10, 10, 10))
    assert wind == {0, 1}
    assert {e for e in range(3) if _supercell_winds(sites, types, (10, 10, 10), e)} == {0, 1}
    assert _pimms_axes(sites, types, (10, 10, 10)) == {0, 1}

    cluster = [tuple(p) for p in _MULTI_WINDING_CLUSTER]
    types = dict.fromkeys(cluster, 0)
    (sites, wind), = _components_with_winding(types, (7, 12, 12))
    assert wind == {0, 1}
    assert _pimms_axes(sites, types, (7, 12, 12)) == wind


def test_the_typed_rule_decides_whether_a_gap_through_the_face_closes_the_loop():
    """A contact rod with a gap of two (three) sites through the face winds the
    box exactly when the LR (SLR) entry for its residue type is nonzero - the
    relation the long-range cluster search and gather walk."""
    dims = (12, 12, 12)
    zero = [[0, 0], [0, 0]]
    on = [[0, 0], [0, -4]]
    gap2 = [(x, 5, 5) for x in range(11)]       # x = 10 and x = 0 are 2 apart through the face
    gap3 = [(x, 5, 5) for x in range(10)]       # x = 9 and x = 0 are 3 apart
    # (with SLR on, the gap-2 rod also closes: x = 10 and x = 1 are 3 apart)
    for rod, LR, SLR, expected in ((gap2, on, zero, {0}), (gap2, zero, on, {0}),
                                   (gap2, zero, zero, set()), (gap3, on, on, {0}),
                                   (gap3, zero, on, {0}), (gap3, on, zero, set())):
        types = dict.fromkeys(rod, 1)
        (sites, wind), = _components_with_winding(types, dims, 3, LR, SLR)
        assert wind == expected
        assert _pimms_axes(sites, types, dims, 3, LR, SLR) == expected


# ---------------------------------------------------------------------------
# A1-2 / P2-1: memory, and one run of the test per cluster
# ---------------------------------------------------------------------------

def _slab_cluster(dims, z_lo, z_hi, rho, seed):
    """The largest contact cluster of a random-occupancy slab that fills the
    periodic xy plane of the box.

    Parameters
    ----------
    dims : tuple of int
        The box.
    z_lo, z_hi : int
        The slab occupies ``z_lo <= z < z_hi``.
    rho : float
        Site occupancy inside the slab.
    seed : int
        Seed for the occupancy.

    Returns
    -------
    tuple
        ``(sites, site_types)`` of the largest component.
    """
    rng = np.random.default_rng(seed)
    occ = np.zeros(dims, dtype=bool)
    occ[:, :, z_lo:z_hi] = rng.random((dims[0], dims[1], z_hi - z_lo)) < rho
    site_types = {tuple(int(x) for x in s): 0 for s in np.argwhere(occ)}
    sites = max((c for c, _ in _components_with_winding(site_types, dims)), key=len)
    return sites, site_types


def test_percolation_test_memory_is_linear_in_the_beads_of_a_box_spanning_slab():
    """A slab of ~4,000 beads that fills the periodic plane. The pairwise
    separation matrix this replaced needed tens of kilobytes per bead here (and
    grew with the square of the bead count); the exact test stays within a few
    hundred bytes per bead plus its fixed probe block."""
    dims = (24, 24, 40)
    sites, site_types = _slab_cluster(dims, 10, 30, 0.35, seed=5)
    assert len(sites) > 3500
    pos = np.array(sites, dtype=np.int64)
    si = cluster_utils.convert_positions_to_single_image_snakesearch(
        pos, list(dims), 1, warn_if_percolating=False)

    tracemalloc.start()
    try:
        axes = cluster_utils.percolating_axes(si, list(dims), 1)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()

    assert axes == [0, 1]
    assert peak < 1000 * len(sites) + 8 * 2 ** 20, (peak, len(sites))


def _inject(lattice, box, positions_by_chain):
    """Drop a hand-built configuration onto the lattice.

    Parameters
    ----------
    lattice : pimms.lattice.Lattice
        The lattice to overwrite.
    box : list of int
        The box dimensions.
    positions_by_chain : dict
        chainID -> list of positions, in bonded order.

    Returns
    -------
    dict
        chainID -> the wrapped positions that were written.
    """
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


def test_cluster_analysis_blanks_a_double_winding_cluster_and_tests_each_cluster_once(
        tmp_path, monkeypatch):
    """Through the real ``ANAFUNCT_cluster_analysis``: two decamers laid out as
    the double-winding staircase (plus a compact block of three dimers as a
    control) in a 10-box. The staircase is an unbounded object, so its Rg is nan, its hull is
    -1 and the step is logged; and the percolation test runs exactly once per
    short-range cluster and once per long-range cluster."""
    box = [10, 10, 10]
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(2, "A" * 10), (3, "AA")],
                          extra={"ANA_CLUSTER_THRESHOLD": 1})
    chain_ids = sorted(state.lattice.chains)
    stairs = [[i % 10, i // 2, 0] for i in range(20)]
    _inject(state.lattice, box, {chain_ids[0]: stairs[:10], chain_ids[1]: stairs[10:],
                                 chain_ids[2]: [[4, 4, 5], [5, 4, 5]],
                                 chain_ids[3]: [[4, 5, 5], [5, 5, 5]],
                                 chain_ids[4]: [[4, 6, 5], [5, 6, 5]]})

    tested = []
    real = cluster_utils.cluster_percolates

    def counting(positions, *args, **kwargs):
        tested.append(len(positions))
        return real(positions, *args, **kwargs)

    monkeypatch.setattr(cluster_utils, "cluster_percolates", counting)

    os.chdir(tmp_path)
    with contextlib.redirect_stdout(open(os.devnull, "w")), \
            pytest.warns(UserWarning, match="percolates the periodic box"):
        state.sim.ANAFUNCT_cluster_analysis(100)

    # clusters are sorted by chain count: the control (3 chains), then the staircase
    assert _rows("CLUSTERS.dat")[-1] == ["100", "3", "2"]
    # two short-range clusters and two long-range clusters, each tested once
    assert sorted(tested) == [6, 6, 20, 20]

    for name in ("CLUSTER_RG.dat", "CLUSTER_ASPH.dat", "LR_CLUSTER_RG.dat", "LR_CLUSTER_ASPH.dat"):
        assert np.isnan(float(_rows(name)[-1][2])), name
    for name in ("CLUSTER_VOL.dat", "LR_CLUSTER_VOL.dat"):
        assert float(_rows(name)[-1][2]) == -1, name
    # the 2 x 3 x 1 control block keeps its real radius of gyration:
    # Rg^2 = var(x) + var(y) = 1/4 + 2/3
    assert float(_rows("CLUSTER_RG.dat")[-1][1]) == pytest.approx(np.sqrt(0.25 + 2.0 / 3.0), abs=5e-5)

    with open("log.txt") as fh:
        logged = [line for line in fh if "percolate" in line]
    assert len(logged) == 1 and "Step 100" in logged[0]


# ---------------------------------------------------------------------------
# A1-3: the radial-profile centre
# ---------------------------------------------------------------------------

def _radial_profile_oracle(points, dims, centre):
    """Shell occupancies by explicit enumeration of the shell sites (periodic).

    Parameters
    ----------
    points : numpy.ndarray
        The cluster's single-image positions.
    dims : list of int
        The box.
    centre : tuple of int
        The profile centre.

    Returns
    -------
    list of float
        Fraction of the sites at Chebyshev distance k = 1 .. min(dims)//2 - 1
        from the centre that hold a bead of the cluster.
    """
    occupied = {tuple(int(v) % d for v, d in zip(p, dims)) for p in points}
    out = []
    for k in range(1, min(dims) // 2):
        shell = [o for o in itertools.product(range(-k, k + 1), repeat=len(dims))
                 if max(abs(x) for x in o) == k]
        hits = sum(tuple((c + x) % d for c, x, d in zip(centre, o, dims)) in occupied for o in shell)
        out.append(hits / len(shell))
    return out


def _nearest_site_half_up(points):
    """The profile centre from its definition, by search rather than by rounding.

    Parameters
    ----------
    points : numpy.ndarray
        The cluster's single-image positions.

    Returns
    -------
    tuple of int
        Per axis, the lattice coordinate closest to the arithmetic mean of the
        positions; when two coordinates are equally close (a half-integer
        mean) the larger one.
    """
    centre = []
    for axis in range(points.shape[1]):
        mean = sum(int(v) for v in points[:, axis]) / len(points)
        candidates = range(int(points[:, axis].min()), int(points[:, axis].max()) + 1)
        centre.append(max(candidates, key=lambda c: (-abs(c - mean), c)))
    return tuple(centre)


def test_radial_profile_is_unchanged_by_a_one_site_translation():
    """An asymmetric cluster (a 3 x 3 x 3 cube with a tail) whose mean x is a
    half-integer. Wherever it sits, its profile is the site-enumeration profile
    about the site nearest the mean with the half-way case rounded up, so a
    one-site translation along any axis cannot change it. Rounding half-to-even
    picked the other candidate site for every second position."""
    dims = [20, 20, 20]
    cube = [[x, y, z] for x in range(3) for y in range(3) for z in range(3)]
    tail = [[3, 1, 1], [4, 1, 1], [5, 1, 1], [6, 1, 1], [3, 0, 1]]
    base = np.array(cube + tail) + 6
    assert base.mean(axis=0)[0] % 1 == 0.5            # (27 + 21) / 32 = 1.5

    # the two candidate centres give different profiles, so the choice is visible:
    # about the cube's own centre the first shell is full, about the next site
    # along the tail it holds 19 of its 26 sites
    low = _radial_profile_oracle(base, dims, (7, 7, 7))
    high = _radial_profile_oracle(base, dims, (8, 7, 7))
    assert low[0] == 1.0 and high[0] == pytest.approx(19 / 26)
    assert _nearest_site_half_up(base) == (8, 7, 7)

    for shift in itertools.product(range(2), repeat=3):
        pts = base + np.array(shift)
        expected = _radial_profile_oracle(pts, dims, _nearest_site_half_up(pts))
        profile = lau.compute_cluster_radial_density_profile([pts], dims, 27)[0]
        assert profile == pytest.approx(expected), shift
        # ... and it is the same profile at every position
        assert profile == pytest.approx(high), shift


def test_radial_profile_of_a_mirrored_half_integer_cluster_uses_the_other_candidate_site():
    """Why mirror invariance cannot hold: the mirror image of the cluster above
    has its mean half-way between the mirror images of the same two sites, and
    rounding up now selects the image of the site that was NOT chosen before -
    so the mirrored cluster's profile is the other candidate's profile."""
    dims = [20, 20, 20]
    cube = [[x, y, z] for x in range(3) for y in range(3) for z in range(3)]
    tail = [[3, 1, 1], [4, 1, 1], [5, 1, 1], [6, 1, 1], [3, 0, 1]]
    base = np.array(cube + tail) + 6
    mirrored = base.copy()
    mirrored[:, 0] = 19 - mirrored[:, 0]

    # the cube's own centre (7, 7, 7) maps to (12, 7, 7), and that is the centre
    assert _nearest_site_half_up(mirrored) == (12, 7, 7)
    low = _radial_profile_oracle(base, dims, (7, 7, 7))
    high = _radial_profile_oracle(base, dims, (8, 7, 7))
    profile = lau.compute_cluster_radial_density_profile([mirrored], dims, 27)[0]
    assert profile == pytest.approx(_radial_profile_oracle(mirrored, dims, (12, 7, 7)))
    assert profile == pytest.approx(low) and profile != pytest.approx(high)


# ---------------------------------------------------------------------------
# A1-4 / PERF1: the LR-flag grid
# ---------------------------------------------------------------------------

def _chain_level_lr_clusters(lattice, LR, SLR):
    """Brute-force partition of the chains into long-range clusters.

    Parameters
    ----------
    lattice : pimms.lattice.Lattice
        The lattice; only chain positions and residue codes are read.
    LR, SLR : numpy.ndarray
        The interaction tables.

    Returns
    -------
    dict
        chainID -> frozenset of the chainIDs in its long-range cluster.
    """
    dims = tuple(int(d) for d in lattice.dimensions)
    site_types, site_chain = {}, {}
    for cid, chain in lattice.chains.items():
        for p, code in zip(chain.get_ordered_positions(), chain.int_sequence):
            site_types[tuple(int(v) for v in p)] = int(code)
            site_chain[tuple(int(v) for v in p)] = cid
    # bonded neighbours are in contact, so bead-level components are chain-level clusters
    out = {}
    for sites, _ in _components_with_winding(site_types, dims, 3, LR.tolist(), SLR.tolist()):
        members = frozenset(site_chain[s] for s in sites)
        for cid in members:
            out[cid] = out.get(cid, frozenset()) | members
    return out


@pytest.mark.parametrize("ff_kind", ["LR", "SLR"])
def test_long_range_cluster_search_reads_only_the_chains_it_reaches_when_tables_are_given(
        tmp_path, ff_kind):
    """With both tables supplied the search must give the brute-force partition
    and must not walk every chain of the system to build an LR-flag grid it
    never reads: only the chains of the cluster being built are asked for their
    LR flags."""
    box = [16, 16, 16]
    state = U.build_state(tmp_path, 3, ff_kind, False, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(6, "AABB"), (6, "AAAA"), (8, "B")], seed=3)
    lattice = state.lattice
    LR = np.asarray(state.ham.LR_residue_interaction_table)
    SLR = np.asarray(state.ham.SLR_residue_interaction_table)
    assert np.any(LR)
    expected = _chain_level_lr_clusters(lattice, LR, SLR)
    assert len(set(expected.values())) > 1          # several clusters, so the test can tell

    asked = []
    for cid, chain in lattice.chains.items():
        real = chain.get_LR_binary_array

        def recording(cid=cid, real=real):
            asked.append(cid)
            return real()

        chain.get_LR_binary_array = recording

    for seed in sorted(lattice.chains):
        asked.clear()
        members = lattice_utils.get_all_chains_in_long_range_cluster(
            seed, lattice, hardwall=False, LR_table=LR, SLR_table=SLR)
        assert frozenset(members) == expected[seed]
        assert set(asked) == set(expected[seed])

    # the structural rule (no tables) still needs the grid, and still works: in
    # this force field every LR-capable pair has a nonzero LR entry, so with an
    # all-zero SLR table swapped for "LR-capable" the answer can only grow
    asked.clear()
    members = lattice_utils.get_all_chains_in_long_range_cluster(sorted(lattice.chains)[0], lattice)
    assert set(asked) == set(lattice.chains)
    assert frozenset(members) >= expected[sorted(lattice.chains)[0]]


# ---------------------------------------------------------------------------
# P2-3 / S1-6: the distance-map accumulator
# ---------------------------------------------------------------------------

def test_distance_map_accumulator_is_not_allocated_until_it_is_first_updated():
    """Every chain owns a DistanceMap from construction; the seqlen x seqlen
    matrix (32 MB at 2000 beads) must only exist once the map is sampled."""
    tracemalloc.start()
    try:
        before = tracemalloc.get_traced_memory()[0]
        dmap = analysis_structures.DistanceMap(2000)
        held = tracemalloc.get_traced_memory()[0] - before
    finally:
        tracemalloc.stop()
    assert held < 100_000
    assert dmap.count == 0
    assert dmap.get_distance_map().shape == (2000, 2000)
    assert not dmap.get_distance_map().any()


def test_distance_map_running_mean_is_bit_identical_to_the_one_line_form():
    """The in-place update must reproduce mean += (new - mean) / (count + 1)
    exactly, whether or not it is allowed to use the caller's array as scratch,
    and must leave the caller's array alone unless told otherwise."""
    rng = np.random.default_rng(8)
    maps = [np.sqrt(rng.integers(0, 400, (17, 17)).astype(np.float64)) for _ in range(7)]

    expected = np.zeros((17, 17))
    for count, m in enumerate(maps):
        expected += (m - expected) / (count + 1)

    kept = analysis_structures.DistanceMap(17)
    consumed = analysis_structures.DistanceMap(17)
    for m in maps:
        original = m.copy()
        kept.update_distance_map(m)
        assert np.array_equal(m, original)
        consumed.update_distance_map(m.copy(), consume=True)

    assert kept.count == consumed.count == 7
    assert np.array_equal(kept.get_distance_map(), expected)
    assert np.array_equal(consumed.get_distance_map(), expected)


def test_distance_map_cost_estimate_matches_what_a_run_allocates_and_writes(tmp_path):
    """The numbers the ANA_DISTMAP warning quotes, against a real Simulation's
    own analysis routines under tracemalloc: two samples through
    ``ANAFUNCT_distance_map``, then ``end_of_simulation_analysis``, which
    averages over the chains and writes the file.

    One map is 8 x seqlen^2 bytes. The estimate is (chains + 2) maps: exactly
    what the end-of-run average holds (the accumulators, the running sum and
    the mean), while a sample holds the accumulators, one instantaneous map
    and the block temporaries of the distance calculation (at most ~54 MB)."""
    seqlen, n_chains = 400, 6
    one_map = 8 * seqlen * seqlen
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[40, 40, 40], chains=[(n_chains, "A" * seqlen)], seed=5)
    sim = state.sim
    os.chdir(tmp_path)

    tracemalloc.start()
    try:
        before = tracemalloc.get_traced_memory()[0]
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim.ANAFUNCT_distance_map(1)
            sim.ANAFUNCT_distance_map(2)
        sampling_peak = tracemalloc.get_traced_memory()[1] - before
        held = tracemalloc.get_traced_memory()[0] - before
        tracemalloc.reset_peak()
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim.end_of_simulation_analysis()
        end_peak = tracemalloc.get_traced_memory()[1] - before
    finally:
        tracemalloc.stop()

    est_memory, est_file = analysis_structures.distance_map_cost(seqlen, n_chains)

    # held for the whole run: one map per chain
    assert held == pytest.approx(n_chains * one_map, rel=0.02)
    # the end-of-run average: the estimate, to the percent
    assert end_peak == pytest.approx(est_memory, rel=0.02)
    assert end_peak == pytest.approx((n_chains + 2) * one_map, rel=0.02)
    # a sample: the accumulators, the instantaneous map, and block temporaries
    assert (n_chains + 1) * one_map <= sampling_peak <= est_memory + 56_000_000
    # and the file
    assert os.path.getsize("DISTANCE_MAP.dat") == pytest.approx(est_file, rel=0.05)


def test_distance_map_analysis_warns_once_about_a_chain_too_long_for_it(tmp_path, monkeypatch):
    """The first ANA_DISTMAP sample names the keyword, the chain length and the
    estimated memory and file size when they exceed the thresholds - once, in
    the log as well as on screen - and says nothing for an ordinary system."""
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[14, 14, 14], chains=[(3, "A" * 12), (2, "AB")])
    sim = state.sim
    os.chdir(tmp_path)

    def distmap_lines():
        if not os.path.exists("log.txt"):
            return []
        with open("log.txt") as fh:
            return [line for line in fh if "ANA_DISTMAP" in line]

    n_before = len(distmap_lines())
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        sim.ANAFUNCT_distance_map(1)
    assert len(distmap_lines()) == n_before          # 12-bead chains: nothing to say

    # the same system against a threshold its 12-mers exceed and its dimers do not:
    # (3 + 2) * 8 * 12**2 = 5760 bytes, against (2 + 2) * 8 * 2**2 = 128
    monkeypatch.setattr(analysis_structures, "DISTANCE_MAP_WARN_MEMORY_BYTES", 4000)
    sim._distance_map_cost_reported = False
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        sim.ANAFUNCT_distance_map(2)
        sim.ANAFUNCT_distance_map(3)
    new = distmap_lines()[n_before:]
    assert len(new) == 1
    assert "3 chain(s) of 12 beads" in new[0] and "ANA_DISTMAP : 0" in new[0]
    assert "memory" in new[0] and "DISTANCE_MAP.dat" in new[0]

    # and the samples were all folded in regardless
    for chain in sim.LATTICE.chains.values():
        assert chain.distance_map.count == 3
