"""
End-to-end tests for the lemonade analysis backend, run against real PIMMS output.

They establish that: loading recovers the exact integer lattice and topology; the
compiled PBC unwrap is bit-identical to PIMMS's reference; the batched analyses
and the per-polymer accessors agree with plain-numpy definitions and cross-validate
against PIMMS's own Rg; and the navigational hierarchy (trajectory -> frame ->
polymer / cluster), slicing and 2D all behave.
"""

import warnings

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms import lattice_utils as lu
from pimms import lattice_analysis_utils as lau
from pimms.lemonade._store import TrajectoryStore
from pimms.lemonade._topology import Topology


# ---------------------------------------------------------------------------
# loading / topology
# ---------------------------------------------------------------------------

def test_load_and_metadata(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)

    assert traj.dimensions == (8, 8, 8)
    assert traj.n_dim == 3
    assert traj.hardwall is False
    assert traj.spacing == pytest.approx(3.65)
    assert traj.n_chains == 4
    assert traj.n_beads == 4 * 12
    assert traj.n_frames == len(traj) >= 2
    assert traj.sequences == ["AABBAABBAABB"] * 4
    assert list(traj.chain_types) == [0, 0, 0, 0]           # keyfile: one CHAIN spec


def test_pdb_chain_ids_preserve_types_without_keyfile():
    import mdtraj as md
    from mdtraj.core.element import carbon

    top = md.Topology()
    specs = [("A", "ALA"), ("A", "GLY"), ("B", "ALA")]
    for chain_id, residue_name in specs:
        chain = top.add_chain(chain_id)
        residue = top.add_residue(residue_name, chain)
        top.add_atom("CA", carbon, residue)

    decoded = Topology.from_mdtraj(top)
    assert decoded.sequences == ["A", "G", "A"]
    # Different sequences with PDB ID A share a PIMMS type; the identical A
    # sequence under PDB ID B is a distinct type.
    assert decoded.chain_types.tolist() == [0, 0, 1]


def test_topology_rejects_misaligned_or_fractional_chain_types():
    with pytest.raises(ValueError, match="one entry per chain"):
        Topology(["AA", "BB"], chain_types=[0])
    with pytest.raises(ValueError, match="non-negative int32"):
        Topology(["AA"], chain_types=[1.5])


def test_trajectory_store_validates_shape_topology_and_bounds():
    topology = Topology(["AA"])
    good = np.array([[[0, 0, 0], [1, 0, 0]]], dtype=np.int32)
    TrajectoryStore(good, (4, 4), 3.65, False, topology)

    with pytest.raises(ValueError, match="shape"):
        TrajectoryStore(good[0], (4, 4), 3.65, False, topology)
    with pytest.raises(ValueError, match="outside dimension"):
        TrajectoryStore(np.array([[[0, 0, 0], [4, 0, 0]]]),
                        (4, 4), 3.65, False, topology)
    with pytest.raises(ValueError, match="topology describes"):
        TrajectoryStore(good[:, :1], (4, 4), 3.65, False, topology)
    with pytest.raises(ValueError, match="dimensions"):
        TrajectoryStore(good, (4.5, 4), 3.65, False, topology)
    with pytest.raises(ValueError, match="spacing"):
        TrajectoryStore(good, (4, 4), float("nan"), False, topology)
    with pytest.raises(ValueError, match="hardwall"):
        TrajectoryStore(good, (4, 4), 3.65, "False", topology)
    with pytest.raises(ValueError, match="temperature"):
        TrajectoryStore(good, (4, 4), 3.65, False, topology, temperature=0)
    with pytest.raises(ValueError, match="finite"):
        TrajectoryStore(good, (4, 4), 3.65, False, topology, times=[np.nan])


def test_trajectory_store_checks_bounds_before_int32_conversion():
    topology = Topology(["A"])
    wrapping_value = np.iinfo(np.int32).max + 2

    with pytest.raises(ValueError, match="outside dimension"):
        TrajectoryStore(
            np.array([[[wrapping_value, 0, 0]]], dtype=np.int64),
            (4, 4), 3.65, False, topology,
        )


def test_trajectory_arrays_are_immutable_so_cached_analyses_cannot_go_stale():
    topology = Topology(["AA"])
    store = TrajectoryStore(
        np.array([[[0, 0, 0], [1, 0, 0]]], dtype=np.int32),
        (4, 4), 3.65, False, topology, times=[0.0],
    )

    with pytest.raises(ValueError, match="read-only"):
        store.positions[0, 0, 0] = 2
    with pytest.raises(ValueError, match="read-only"):
        store.times[0] = 1.0
    with pytest.raises(ValueError, match="read-only"):
        store.whole_positions()[0, 0, 0] = 2
    with pytest.raises(ValueError, match="read-only"):
        store.radius_of_gyration()[0, 0] = 99.0

    # Read-only storage must not prevent the grid-painting kernel from running.
    assert store.frame_grid(0)[0, 0] == 1


@pytest.mark.parametrize("kwargs", [
    {"spacing": 0}, {"spacing": float("nan")}, {"dimensions": (8,)},
    {"dimensions": (8, -1)}, {"dimensions": (8.5, 8, 8)},
    {"n_frames": 0}, {"hardwall": "False"},
])
def test_load_rejects_invalid_geometry_options(traj3d_files, kwargs):
    xtc, pdb, _keyfile = traj3d_files
    with pytest.raises(ValueError):
        lemonade.load(xtc=xtc, pdb=pdb, **kwargs)


def test_lattice_roundtrip_is_exact_and_in_box(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    pos = traj.positions
    assert pos.dtype == np.int32
    # canonicalised into the box on every axis
    for axis, extent in enumerate(traj.dimensions):
        assert pos[..., axis].min() >= 0
        assert pos[..., axis].max() < extent


# ---------------------------------------------------------------------------
# PBC unwrap kernel
# ---------------------------------------------------------------------------

def test_unwrap_is_bit_identical_to_pimms_and_contiguous(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    dims = list(traj.dimensions)

    n_straddling = 0
    for f in range(traj.n_frames):
        for c in range(traj.n_chains):
            p = lemonade.Polymer(traj.store, f, c)
            raw = p.positions[:, :3].tolist()
            if lu.do_positions_stradle_pbc_boundary(raw):
                n_straddling += 1
            ref = np.asarray(lu.make_chain_whole(raw, dims))
            got = np.asarray(p.whole_positions)
            # same up to the (identical) anchor translation
            assert np.array_equal(got - got[0], ref - ref[0])
            # and genuinely whole: no bond jumps a boundary
            assert np.abs(np.diff(got, axis=0)).max() <= 1

    assert n_straddling > 0, "fixture is meant to have straddling chains"


# ---------------------------------------------------------------------------
# analyses
# ---------------------------------------------------------------------------

def test_rg_matches_direct_definition(traj3d_files):
    # Rg computed on the whole (contiguous) chain must equal its mathematical
    # definition sqrt(<|r - r_com|^2>). This validates the batched reduceat math
    # even for chains that heavily straddle the box.
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    whole = traj.whole_positions()[..., :traj.n_dim].astype(np.float64)
    off = traj.topology.offsets
    rg = traj.radius_of_gyration()
    assert rg.shape == (traj.n_frames, traj.n_chains)
    for f in range(traj.n_frames):
        for c in range(traj.n_chains):
            pts = whole[f, off[c]:off[c + 1]]
            ref = np.sqrt(np.mean(((pts - pts.mean(axis=0)) ** 2).sum(axis=1)))
            assert rg[f, c] == pytest.approx(ref)


def test_rg_agrees_with_pimms_in_dilute_regime(traj_dilute_files):
    # lemonade's Rg must match PIMMS's own get_polymeric_properties, which since
    # 1.0.8 is computed on the same whole (bond-walked) chain.
    xtc, pdb, keyfile = traj_dilute_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    dims = list(traj.dimensions)
    for f in range(traj.n_frames):
        for c in range(traj.n_chains):
            p = lemonade.Polymer(traj.store, f, c)
            pimms_rg = lau.get_polymeric_properties(p.positions[:, :3].tolist(), dims)[0]
            assert p.radius_of_gyration == pytest.approx(pimms_rg, abs=0.03)


def _bond_walk(raw: np.ndarray, dims: np.ndarray) -> np.ndarray:
    """Make one chain whole by walking its bonds under the minimum image.

    Written here from the definition, independently of lemonade's unwrap kernel:
    each bead is placed at its predecessor plus the minimum-image bond vector.

    Parameters
    ----------
    raw : numpy.ndarray
        ``(L, n_dim)`` wrapped lattice positions of one chain, in bonded order.
    dims : numpy.ndarray
        ``(n_dim,)`` box extent.

    Returns
    -------
    numpy.ndarray
        ``(L, n_dim)`` float64 whole positions, anchored on the first bead.
    """
    raw = np.asarray(raw, dtype=np.float64)
    whole = raw.copy()
    for i in range(1, len(raw)):
        bond = raw[i] - raw[i - 1]
        bond -= dims * np.round(bond / dims)
        whole[i] = whole[i - 1] + bond
    return whole


@pytest.mark.parametrize("fixture", ["traj3d_files", "traj2d_files"])
def test_polymer_properties_match_a_plain_numpy_oracle(fixture, request):
    """Rg, asphericity, end-to-end distance and COM against their definitions.

    Every chain of every frame is made whole by a bond walk written here and the
    four quantities are computed from it with plain numpy: Rg is the root mean
    square distance from the centre of mass; the asphericity is ``lam_3 -
    (lam_1 + lam_2) / 2`` of the gyration tensor in 3D and ``lam_2 - lam_1`` in 2D;
    the end-to-end distance is the distance between the first and last bead. Both
    the batched arrays and the per-polymer accessors must agree with it. (This
    replaces a test that compared the per-polymer accessors with the batched arrays
    they are read from, which could not fail.)
    """
    xtc, pdb, keyfile = request.getfixturevalue(fixture)
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    nd = traj.n_dim
    dims = np.asarray(traj.dimensions, dtype=np.float64)
    off = traj.topology.offsets
    batched = {"rg": traj.radius_of_gyration(), "asph": traj.asphericity(),
               "ete": traj.end_to_end_distance(), "com": traj.center_of_mass()}

    for f in range(traj.n_frames):
        for c in range(traj.n_chains):
            whole = _bond_walk(traj.positions[f, off[c]:off[c + 1], :nd], dims)
            com = whole.mean(axis=0)
            d = whole - com
            eig = np.linalg.eigvalsh(d.T @ d / len(d))
            expected = {
                "rg": float(np.sqrt((d * d).sum(axis=1).mean())),
                "asph": float(eig[2] - 0.5 * (eig[0] + eig[1]) if nd == 3 else eig[1] - eig[0]),
                "ete": float(np.linalg.norm(whole[-1] - whole[0])),
            }
            p = traj[f][c]
            per_polymer = {"rg": p.radius_of_gyration, "asph": p.asphericity,
                           "ete": p.end_to_end_distance}
            for key, value in expected.items():
                assert batched[key][f, c] == pytest.approx(value, abs=1e-9)
                assert per_polymer[key] == pytest.approx(value, abs=1e-9)
            # the COM is only defined up to the image the whole chain was built in
            for got in (batched["com"][f, c], p.center_of_mass):
                shift = np.asarray(got, dtype=np.float64) - com
                assert np.allclose(shift - dims * np.round(shift / dims), 0.0, atol=1e-9)


def test_distance_map_symmetric_and_zero_diagonal(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    dm = traj[0][0].distance_map()
    n = len(traj[0][0])
    assert dm.shape == (n, n)
    assert np.array_equal(dm, dm.T)
    assert np.allclose(np.diag(dm), 0.0)


# ---------------------------------------------------------------------------
# navigation
# ---------------------------------------------------------------------------

def test_navigation_and_indexing(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)

    frame = traj[0]
    assert isinstance(frame, lemonade.Frame)
    assert len(frame) == traj.n_chains
    assert isinstance(frame[0], lemonade.Polymer)
    assert frame[-1].chain_index == traj.n_chains - 1
    assert [p.chain_index for p in frame] == list(range(traj.n_chains))
    with pytest.raises(IndexError):
        frame[traj.n_chains]
    with pytest.raises(IndexError):
        traj[traj.n_frames]
    # frames iterate
    assert len(list(traj)) == traj.n_frames


def test_slicing_returns_subtrajectory(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    sub = traj[::2]
    assert isinstance(sub, lemonade.LatticeTrajectory)
    assert sub.n_frames == len(range(0, traj.n_frames, 2))
    assert sub.n_chains == traj.n_chains
    # sliced data matches the parent's strided frames
    assert np.array_equal(sub.positions[0], traj.positions[0])
    assert np.array_equal(sub.positions[1], traj.positions[2])


# ---------------------------------------------------------------------------
# clusters
# ---------------------------------------------------------------------------

def test_clusters_partition_all_chains(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    clusters = traj[-1].clusters
    # every chain lands in exactly one cluster
    members = [c for cl in clusters for c in cl.chain_indices]
    assert sorted(members) == list(range(traj.n_chains))
    # geometry is computable and sane for the largest cluster. In this fixture the
    # largest cluster can wind the box, in which case the gather correctly warns
    # that its shape is BFS-order dependent; this test only checks that the
    # geometry is computable and has the right shape, not its value, so the
    # warning is silenced here
    big = clusters[0]
    assert big.n_beads == sum(len(traj[-1][c]) for c in big.chain_indices)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="single-image gather: cluster percolates")
        assert big.radius_of_gyration > 0
        assert big.single_image_positions().shape == (big.n_beads, 3)
    assert big.bead_type_composition.get("A", 0) + big.bead_type_composition.get("B", 0) == big.n_beads


# ---------------------------------------------------------------------------
# 2D + flexible loading
# ---------------------------------------------------------------------------

def test_2d_trajectory(traj2d_files):
    xtc, pdb, keyfile = traj2d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert traj.dimensions == (10, 10)
    assert traj.n_dim == 2
    rg = traj.radius_of_gyration()
    assert rg.shape == (traj.n_frames, traj.n_chains)
    assert np.all(rg > 0)
    # z is flat
    assert int(traj.positions[..., 2].max()) == 0


def test_load_without_keyfile_infers_metadata(traj3d_files):
    xtc, pdb, _keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb)         # no keyfile
    assert traj.dimensions == (8, 8, 8)            # inferred from the box
    assert traj.spacing == pytest.approx(3.65)     # default
    assert traj.n_chains == 4


def test_load_pdb_only_single_frame(traj3d_files):
    _xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(pdb=pdb, keyfile=keyfile)
    assert traj.n_frames == 1
    assert traj.n_chains == 4
    assert traj[0][0].radius_of_gyration > 0


def test_load_requires_topology_for_xtc(traj3d_files):
    xtc, _pdb, _keyfile = traj3d_files
    with pytest.raises(ValueError, match="topology"):
        lemonade.load(xtc=xtc)


# ---------------------------------------------------------------------------
# batched-analysis internals: memory-bounded gyration tensor
# ---------------------------------------------------------------------------

def test_gyration_eigenvalues_match_the_full_outer_product_form():
    """Accumulating the tensor per component must be bit-identical.

    The previous form built an (n_frames, n_beads, k, k) intermediate before reducing -
    0.7 GB for a 1000-frame, 10k-bead trajectory - which put long trajectories out of
    reach. Reducing one component at a time walks each chain's beads in the same order,
    so the numbers are unchanged.
    """
    from pimms.lemonade import _analysis

    def reference(whole, offsets, lengths):
        d = _analysis._centered(whole, offsets, lengths, None)
        outer = d[:, :, :, np.newaxis] * d[:, :, np.newaxis, :]
        tensor = np.add.reduceat(outer, offsets[:-1], axis=1) / lengths[np.newaxis, :, np.newaxis, np.newaxis]
        return np.linalg.eigvalsh(tensor)

    rng = np.random.default_rng(5)
    for k in (2, 3):
        for lens in ([4, 4, 4], [1, 7, 3, 12], [2, 2]):
            lengths = np.array(lens, dtype=np.int64)
            offsets = np.zeros(len(lens) + 1, dtype=np.int64)
            np.cumsum(lengths, out=offsets[1:])
            whole = rng.integers(-20, 20, size=(6, int(offsets[-1]), k)).astype(np.float64)

            assert np.array_equal(_analysis.gyration_eigenvalues(whole, offsets, lengths),
                                  reference(whole, offsets, lengths))
