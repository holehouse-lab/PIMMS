"""
Regression tests for the third 1.0.8 audit pass, which hunted for bugs that
produce SILENTLY WRONG OUTPUT (wrong numbers in files, a move mix that
contradicts the keyfile, moves that corrupt state without tripping the energy
check) rather than crashes.

Every expected value here is hand-derived or computed with plain numpy from the
raw bead positions - never by calling the code under test a second way.
"""

import contextlib
import os
import random

import numpy as np
import pytest

import mdtraj as md

from pimms import chain as chain_module
from pimms import lattice_analysis_utils as lau
from pimms import lattice_utils
from pimms import nonequilibrium_utils
from pimms.keyfile_parser import KeyFileParser
from pimms.latticeExceptions import KeyFileException, SimulationException
from pimms.tests import kernel_test_utils as U


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# a 9-bead rod that crosses the x boundary of a 13-box: raw x = 9,10,11,12,0,1,2,3,4
_ROD = [[9, 0, 0], [10, 0, 0], [11, 0, 0], [12, 0, 0], [0, 0, 0],
        [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]]


def _chain(positions, dims, hardwall=False, chainID=1):
    n = len(positions)
    return chain_module.Chain(
        lattice_grid=np.zeros(dims, dtype=np.int32), dimensions=list(dims),
        sequence="A" * n, int_seq=[1] * n, LR_int_seq=[1] * n, LR_IDX=[],
        chainID=chainID, chainType=0, chain_positions=[list(p) for p in positions],
        hardwall=hardwall)


@pytest.fixture(autouse=True)
def _restore_cwd_and_rng():
    """Several tests chdir into their tmp_path and reseed the global RNG; put
    both back so later tests do not depend on this module having run."""
    cwd = os.getcwd()
    state = random.getstate()
    yield
    random.setstate(state)
    os.chdir(cwd)


def _run(state, tmp_path):
    os.chdir(tmp_path)
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state.sim.run_simulation()


# ---------------------------------------------------------------------------
# 1. intra-chain observables are computed on the chain made whole
# ---------------------------------------------------------------------------


def test_straddling_rod_end_to_end_and_internal_scaling_use_whole_chain_geometry():
    """The rod's true end-to-end distance is 8. Minimum image gave 13 - 8 = 5, and
    the internal-scaling profile at gaps 7 and 8 gave 6 and 5 instead of 7 and 8."""
    chain = _chain(_ROD, (13, 13, 13))
    assert chain.analysis_get_end_to_end_distance() == pytest.approx(8.0)
    assert chain.analysis_get_residue_residue_distance(0, 8) == pytest.approx(8.0)

    profile = chain.analysis_get_instantaneous_internal_scaling(mode="dict")
    for gap, value in profile.items():
        assert value == pytest.approx(float(gap)), gap

    dmap = chain.analysis_get_instantaneous_distance_map()
    idx = np.arange(9)
    assert np.allclose(dmap, np.abs(idx[:, None] - idx[None, :]))

    # cumulative accumulators see the same geometry
    chain.analysis_update_internal_scaling()
    # the cumulative accessor returns the mean distance per gap, ordered by gap 1..L-1
    values = np.asarray(chain.analysis_get_cumulative_internal_scaling(), dtype=float)
    assert np.allclose(values, np.arange(1, len(values) + 1))


def test_rg_and_asphericity_correct_for_chains_spanning_over_half_the_box():
    """An L-shaped 25-mer in a 13-box (13 beads along x from x=0, then 12 up y) never
    crosses a boundary, yet the old COM-relative image selection shifted the far
    beads by a box length and tore it. Expected values from plain numpy on the
    Cartesian coordinates."""
    positions = [[x, 0, 0] for x in range(13)] + [[12, y, 0] for y in range(1, 13)]
    chain = _chain(positions, (13, 13, 13))
    pos = np.asarray(positions, dtype=float)
    delta = pos - pos.mean(axis=0)
    tensor = delta.T @ delta / len(pos)
    ev = np.sort(np.linalg.eigvalsh(tensor))
    rg_expected = float(np.sqrt(ev.sum()))
    rg, asph = chain.analysis_get_polymeric_properties()
    assert rg == pytest.approx(rg_expected)
    # asphericity (3D convention): l3^2 - (l1^2 + l2^2)/2 ... normalised by Rg^4
    # - pin only that it is the value the Cartesian tensor gives, via the
    # public helper with the PBC correction disabled
    _rg_ref, asph_ref = lau.get_polymeric_properties(positions, [13, 13, 13], pbc_correction=False)
    assert asph == pytest.approx(asph_ref)
    assert rg == pytest.approx(_rg_ref)

    # the straddling rod: exact Rg of 9 collinear points, sqrt(60/9)
    rod = _chain(_ROD, (13, 13, 13))
    assert rod.analysis_get_radius_of_gyration() == pytest.approx(np.sqrt(60.0 / 9.0))


def test_random_walks_match_a_bond_walk_oracle_in_a_small_box():
    """Self-avoiding walks in a 10-box, many spanning more than half of it: Rg,
    end-to-end and every pair distance must equal the values computed on the
    bond-walked (unwrapped) coordinates."""
    rng = np.random.default_rng(5)
    L = 10
    steps = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
    checked = 0
    for _ in range(300):
        # grow an unwrapped SAW of 24 beads
        walk = [np.array([5, 5, 5])]
        occupied = {tuple(walk[0] % L)}
        ok = True
        for _b in range(23):
            for _try in range(20):
                nxt = walk[-1] + steps[rng.integers(6)]
                if tuple(nxt % L) not in occupied:
                    walk.append(nxt)
                    occupied.add(tuple(nxt % L))
                    break
            else:
                ok = False
                break
        if not ok:
            continue
        unwrapped = np.asarray(walk, dtype=float)
        wrapped = (np.asarray(walk) % L).tolist()
        chain = _chain(wrapped, (L, L, L))

        delta = unwrapped - unwrapped.mean(axis=0)
        rg_oracle = float(np.sqrt((delta ** 2).sum(axis=1).mean()))
        e2e_oracle = float(np.linalg.norm(unwrapped[-1] - unwrapped[0]))
        dmap_oracle = np.linalg.norm(unwrapped[:, None, :] - unwrapped[None, :, :], axis=2)

        assert chain.analysis_get_radius_of_gyration() == pytest.approx(rg_oracle)
        assert chain.analysis_get_end_to_end_distance() == pytest.approx(e2e_oracle)
        assert np.allclose(chain.analysis_get_instantaneous_distance_map(), dmap_oracle)
        checked += 1
    assert checked > 200


def test_hardwall_chain_observables_are_unchanged():
    positions = [[i, 4, 4] for i in range(9)]
    chain = _chain(positions, (12, 12, 12), hardwall=True)
    assert chain.get_analysis_positions() is chain.positions
    assert chain.analysis_get_end_to_end_distance() == pytest.approx(8.0)
    assert chain.analysis_get_radius_of_gyration() == pytest.approx(np.sqrt(60.0 / 9.0))


# ---------------------------------------------------------------------------
# 2. LR clusters connect only through pairs with nonzero interaction energy
# ---------------------------------------------------------------------------

def _place(lat, positions_by_chain):
    lat.grid[:] = 0
    lat.type_grid[:] = 0
    for cid, pos in positions_by_chain.items():
        lat.chains[cid].set_ordered_positions([list(p) for p in pos])
        lattice_utils.place_chain_by_position(pos, lat.grid, cid, safe=False)
    lat.initialize_type_grid()


def test_lr_clusters_ignore_zero_energy_shells(tmp_path):
    """With an LR-only parameter file the SLR table is all zeros, so two LR chains
    at Chebyshev distance 3 have ZERO interaction energy - they must not be one LR
    cluster. At Chebyshev distance 2 (nonzero LR entry) they must."""
    state = U.build_state(tmp_path, 3, "LR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[13, 13, 13], chains=[(2, "AAAA")])
    lat, ham = state.lattice, state.ham
    assert not np.any(ham.SLR_residue_interaction_table)
    tables = dict(LR_table=ham.LR_residue_interaction_table,
                  SLR_table=ham.SLR_residue_interaction_table)

    A = [[0, 0, 6], [1, 0, 6], [2, 0, 6], [3, 0, 6]]
    far = [[0, 6, 6], [1, 6, 6], [2, 6, 6], [3, 6, 6]]
    cheb3 = [[0, 3, 6], [1, 3, 6], [2, 3, 6], [3, 3, 6]]
    cheb2 = [[0, 2, 6], [1, 2, 6], [2, 2, 6], [3, 2, 6]]

    _place(lat, {1: A, 2: far})
    E_far = ham.evaluate_total_energy(lat)[0]
    _place(lat, {1: A, 2: cheb3})
    assert ham.evaluate_total_energy(lat)[0] == E_far          # zero interaction
    assert sorted(map(sorted, lau.get_LR_cluster_distribution(lat, hardwall=False, **tables))) == [[1], [2]]
    # without the tables the structural (both-LR-capable) rule is documented behaviour
    assert sorted(map(sorted, lau.get_LR_cluster_distribution(lat, hardwall=False))) == [[1, 2]]

    _place(lat, {1: A, 2: cheb2})
    assert ham.evaluate_total_energy(lat)[0] != E_far          # LR contact counts
    assert sorted(map(sorted, lau.get_LR_cluster_distribution(lat, hardwall=False, **tables))) == [[1, 2]]


def test_simulation_cluster_analysis_passes_the_interaction_tables(tmp_path, monkeypatch):
    state = U.build_state(tmp_path, 3, "LR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[13, 13, 13], chains=[(2, "AAAA")])
    seen = {}

    def recorder(latticeObject, hardwall=False, LR_table=None, SLR_table=None):
        seen["LR"] = LR_table
        seen["SLR"] = SLR_table
        return [[cid] for cid in latticeObject.chains]

    monkeypatch.setattr(lau, "get_LR_cluster_distribution", recorder)
    os.chdir(tmp_path)
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state.sim.ANAFUNCT_cluster_analysis(1)
    assert seen["LR"] is state.ham.LR_residue_interaction_table
    assert seen["SLR"] is state.ham.SLR_residue_interaction_table


# ---------------------------------------------------------------------------
# 3. never-sampled INTSCAL / DISTANCE_MAP are not written as zeros
# ---------------------------------------------------------------------------

def test_unsampled_internal_scaling_is_not_written_as_zeros(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[12, 12, 12], chains=[(4, "AABBA")], n_steps=30, equilibration=2,
                          extra={"ANA_INTSCAL": 1000, "ANA_DISTMAP": 1000, "ANA_POL": 10,
                                 "PRINT_FREQ": 1000, "XTC_FREQ": 1000, "ENERGY_CHECK": 0})
    _run(state, tmp_path)
    assert not (tmp_path / "INTSCAL.dat").exists()
    assert not (tmp_path / "INTSCAL_SQUARED.dat").exists()
    assert not (tmp_path / "DISTANCE_MAP.dat").exists()
    # SCALING_INFORMATION.dat goes too. It used to be written with its -1 rows on
    # the grounds that -1 is a sentinel rather than plausible data, but -1 is also
    # what a chain below the 26-bead fitting floor writes in a fully sampled run
    # (these 5-mers would write exactly the same file either way), so it never
    # distinguished the two cases. The warning does.
    assert not (tmp_path / "SCALING_INFORMATION.dat").exists()
    assert "never sampled" in (tmp_path / "log.txt").read_text()


def test_sampled_internal_scaling_is_still_written(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[12, 12, 12], chains=[(4, "AABBA")], n_steps=30, equilibration=2,
                          extra={"ANA_INTSCAL": 5, "ANA_DISTMAP": 5,
                                 "PRINT_FREQ": 1000, "XTC_FREQ": 1000, "ENERGY_CHECK": 0})
    _run(state, tmp_path)
    profile = np.loadtxt(tmp_path / "INTSCAL.dat")
    assert profile.shape[0] == 4 and np.all(profile[:, 1] > 0)
    assert (tmp_path / "DISTANCE_MAP.dat").exists()


# ---------------------------------------------------------------------------
# 4. move selection honours the keyfile in systems with monomers
# ---------------------------------------------------------------------------

def test_monomer_draws_never_become_crankshaft_megamoves(tmp_path):
    """MOVE_CRANKSHAFT : 0 with monomers and rotate/translate: the crankshaft column
    of MOVE_FREQS.dat must stay at zero (it used to fill with
    CRANKSHAFT_SUBSTEPS attempts per monomer rotate draw)."""
    state = U.build_state(tmp_path, 3, "SR", False,
                          {"MOVE_CRANKSHAFT": 0.0, "MOVE_CHAIN_ROTATE": 0.5, "MOVE_CHAIN_TRANSLATE": 0.5},
                          box=[12, 12, 12], chains=[(6, "A"), (2, "AAAAAA")], n_steps=30, equilibration=1,
                          extra={"ANA_ACCEPTANCE": 10, "PRINT_FREQ": 1000, "XTC_FREQ": 1000,
                                 "ENERGY_CHECK": 0, "CRANKSHAFT_SUBSTEPS": 500})
    log = (tmp_path / "log.txt").read_text()
    assert "cannot move a single bead" in log             # startup warning
    _run(state, tmp_path)
    freqs = np.loadtxt(tmp_path / "MOVE_FREQS.dat")
    last = freqs[-1] if freqs.ndim == 2 else freqs
    assert last[0] == 30
    assert last[1] == 0                                   # crankshaft never ran
    assert last[3] + last[2] == 30                        # every step was a rotate or a translate


def test_moveset_that_can_never_act_is_refused(tmp_path):
    with pytest.raises(SimulationException, match="No enabled move can act"):
        U.build_state(tmp_path, 3, "SR", False, {"MOVE_CHAIN_PIVOT": 1.0},
                      box=[12, 12, 12], chains=[(6, "AA")])


def test_all_chains_frozen_is_refused(tmp_path):
    (tmp_path / "fz.in").write_text("C 1 2 3\n")
    with pytest.raises(SimulationException, match="Every chain in the system is frozen"):
        U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                      box=[12, 12, 12], chains=[(3, "AABB")], extra={"FREEZE_FILE": "fz.in"})


# ---------------------------------------------------------------------------
# 5. input validation that used to pass silently
# ---------------------------------------------------------------------------

def test_solvent_symbol_in_a_chain_sequence_is_rejected(tmp_path):
    with pytest.raises(Exception, match="solvent symbol"):
        U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                      box=[12, 12, 12], chains=[(2, "A0A")])


def test_empty_path_keyword_is_rejected(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    extra={"FREEZE_FILE": ""})
    os.chdir(tmp_path)
    with pytest.raises(KeyFileException, match="expects a file path"):
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            KeyFileParser("KEYFILE.kf")


def test_tsmmc_jump_temperature_is_checked_against_the_hottest_quench_temperature(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False,
                    {"MOVE_CRANKSHAFT": 0.5, "MOVE_CTSMMC": 0.5}, n_steps=40, equilibration=4,
                    extra={"QUENCH_RUN": "True", "QUENCH_FREQ": 2, "QUENCH_STEPSIZE": 10,
                           "QUENCH_START": 50, "QUENCH_END": 100, "TSMMC_JUMP_TEMP": 70,
                           "QUENCH_AS_EQUILIBRATION": "False"})
    os.chdir(tmp_path)
    with pytest.raises(KeyFileException, match="highest simulation temperature"):
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            KeyFileParser("KEYFILE.kf")


def test_quench_snaps_to_the_target_within_rounding():
    T = 1.0
    visited = []
    for _ in range(8):
        T = nonequilibrium_utils.update_temperature_in_quench(0.1, 1.0, 0.2, T, True)
        visited.append(T)
    assert visited[-1] == 0.2                     # exact, not 0.20000000000000004
    assert all(t >= 0.2 for t in visited)
    # heating direction
    T = 0.2
    for _ in range(8):
        T = nonequilibrium_utils.update_temperature_in_quench(-0.1, 0.2, 1.0, T, True)
    assert T == 1.0


def test_angle_energy_of_short_chains_is_an_integer(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[12, 12, 12], chains=[(3, "AA"), (2, "AABB")])
    chain = next(c for c in state.lattice.chains.values() if len(c.positions) == 2)
    value = state.ham.evaluate_angle_energy(chain.get_ordered_positions(), chain.int_sequence,
                                            state.lattice.dimensions)
    assert isinstance(value, (int, np.integer)) and not isinstance(value, float)
    assert isinstance(state.ham.evaluate_total_energy(state.lattice)[0], (int, np.integer))


# ---------------------------------------------------------------------------
# 6. grid consistency check (the blind spot of ENERGY_CHECK)
# ---------------------------------------------------------------------------

def test_grid_consistency_check_sees_a_type_grid_corruption(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[12, 12, 12], chains=[(3, "AABB")])
    lat = state.lattice
    assert lat.check_grid_consistency() == []
    chain = lat.chains[1]
    site = tuple(chain.get_ordered_positions()[1])
    correct = int(lat.type_grid[site])
    lat.type_grid[site] = 0
    problems = lat.check_grid_consistency()
    assert any("type_grid holds 0" in p for p in problems)
    lat.type_grid[site] = correct
    lat.grid[site] = 99
    problems = lat.check_grid_consistency()
    assert any("grid holds chain 99" in p for p in problems)


def test_energy_check_raises_on_a_type_grid_corruption(tmp_path):
    from pimms.latticeExceptions import SimulationEnergyException
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CHAIN_TRANSLATE": 1.0},
                          box=[12, 12, 12], chains=[(3, "AABB")], n_steps=6, equilibration=1,
                          extra={"ENERGY_CHECK": 3, "PRINT_FREQ": 1000, "XTC_FREQ": 1000})
    chain = state.lattice.chains[1]
    site = tuple(chain.get_ordered_positions()[1])
    state.lattice.type_grid[site] = 0          # tracked and recomputed energies now agree - wrongly
    os.chdir(tmp_path)
    with pytest.raises(SimulationEnergyException, match="inconsistent"):
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            state.sim.run_simulation()


# ---------------------------------------------------------------------------
# 7. trajectory files
# ---------------------------------------------------------------------------

def test_save_at_end_with_no_qualifying_frame_writes_frame_zero_only(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[12, 12, 12], chains=[(3, "AABB")], n_steps=14, equilibration=2,
                          extra={"XTC_FREQ": 20, "SAVE_AT_END": "True", "PRINT_FREQ": 1000,
                                 "ENERGY_CHECK": 0})
    _run(state, tmp_path)
    traj = md.load(str(tmp_path / "traj.xtc"), top=str(tmp_path / "START.pdb"))
    assert traj.n_frames == 1
    assert traj.time.tolist() == [0.0]


def test_resized_equilibration_save_at_end_does_not_put_the_production_box_in_eq_traj(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", True, {"MOVE_CRANKSHAFT": 1.0},
                          box=[14, 14, 14], chains=[(3, "AABB")], n_steps=14, equilibration=6,
                          extra={"RESIZED_EQUILIBRATION": "8 8 8", "XTC_FREQ": 10,
                                 "SAVE_AT_END": "True", "SAVE_EQ": "True", "PRINT_FREQ": 1000,
                                 "ENERGY_CHECK": 0})
    _run(state, tmp_path)
    eq = md.load(str(tmp_path / "eq_traj.xtc"), top=str(tmp_path / "eq_START.pdb"))
    assert eq.n_frames == 1
    assert eq.unitcell_lengths[0][0] == pytest.approx(8 * 0.365, rel=1e-4)
    # every coordinate of the single frame lies inside the equilibration box
    assert float(eq.xyz.max()) < 8 * 0.365


def test_resized_equilibration_last_equilibration_frame_stays_in_the_small_box(tmp_path):
    """EQUILIBRATION is the last equilibration step: its frame belongs to eq_traj.xtc
    (small box), and traj.xtc holds the production frames only."""
    state = U.build_state(tmp_path, 3, "SR", True, {"MOVE_CRANKSHAFT": 1.0},
                          box=[24, 24, 24], chains=[(6, "AABB")], n_steps=30, equilibration=10,
                          extra={"RESIZED_EQUILIBRATION": "12 12 12", "XTC_FREQ": 2,
                                 "SAVE_EQ": "True", "PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run(state, tmp_path)
    eq = md.load(str(tmp_path / "eq_traj.xtc"), top=str(tmp_path / "eq_START.pdb"))
    prod = md.load(str(tmp_path / "traj.xtc"), top=str(tmp_path / "START.pdb"))
    assert eq.n_frames == 6                     # steps 0, 2, 4, 6, 8, 10
    assert prod.n_frames == 11                  # start + steps 12..30
    assert float(eq.xyz.max()) < 12 * 0.365
    assert eq.unitcell_lengths[0][0] == pytest.approx(12 * 0.365, rel=1e-4)
    assert prod.unitcell_lengths[0][0] == pytest.approx(24 * 0.365, rel=1e-4)


def test_stale_single_type_and_eq_outputs_are_removed_at_startup(tmp_path):
    first = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[12, 12, 12], chains=[(4, "AABBA")], n_steps=20, equilibration=2,
                          extra={"ANA_INTSCAL": 5, "ANA_DISTMAP": 5, "PRINT_FREQ": 1000,
                                 "XTC_FREQ": 1000, "ENERGY_CHECK": 0})
    _run(first, tmp_path)
    for name in ("INTSCAL.dat", "INTSCAL_SQUARED.dat", "DISTANCE_MAP.dat", "SCALING_INFORMATION.dat"):
        assert (tmp_path / name).exists()
    (tmp_path / "eq_traj.xtc").write_bytes(b"stale")
    (tmp_path / "CONFIG_AT_ENERGY_FAIL.pdb").write_text("stale")

    second = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                           box=[12, 12, 12], chains=[(3, "AABBA"), (3, "AAAA")], n_steps=20,
                           equilibration=2,
                           extra={"ANA_INTSCAL": 5, "ANA_DISTMAP": 5, "PRINT_FREQ": 1000,
                                  "XTC_FREQ": 1000, "ENERGY_CHECK": 0})
    _run(second, tmp_path)
    for name in ("INTSCAL.dat", "INTSCAL_SQUARED.dat", "DISTANCE_MAP.dat", "SCALING_INFORMATION.dat",
                 "eq_traj.xtc", "CONFIG_AT_ENERGY_FAIL.pdb"):
        assert not (tmp_path / name).exists(), name
    assert (tmp_path / "CHAIN_0_INTSCAL.dat").exists()
    assert (tmp_path / "CHAIN_1_INTSCAL.dat").exists()


# ---------------------------------------------------------------------------
# 8. cluster_rotate under HARDWALL never wraps a bead through the wall
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("placement,expect_accepted", [("corner", False), ("legal", True)])
def test_cluster_rotate_hardwall_rejects_wall_crossing_monomer(tmp_path, placement, expect_accepted):
    """7x9 hardwall box: a vertical 7-mer at x=3 with a monomer touching its top
    bead. 'corner': the monomer sits at (2,7), so every rotation about the middle
    bead carries a bead through a wall (a 90-degree turn maps the monomer to
    x = -1); it used to be periodically wrapped to x = 6 (a different contact
    count, committed as energy neutral, reverse move impossible), so nothing may
    be accepted. 'legal': the rod spans y = 1..7 with the monomer at (3,8), so
    the 180-degree turn stays in the box and must be accepted with an exact
    energy change - the positive control that separates "rejects illegal
    rotations" from "rejects everything"."""
    state = U.build_state(tmp_path, 2, "SR", True, {"MOVE_CLUSTER_ROTATE": 1.0}, box=[7, 9],
                          chains=[(1, "AAAAAAA"), (1, "B"), (1, "A")], seed=1)
    lat, ham, sim = state.lattice, state.ham, state.sim
    for cid in list(lat.chains):
        ch = lat.chains[cid]
        lattice_utils.delete_chain_by_position(ch.get_ordered_positions(), lat.grid, cid)
    if placement == "corner":
        _place(lat, {1: [[3, y] for y in range(7)], 2: [[2, 7]], 3: [[6, 8]]})
    else:
        _place(lat, {1: [[3, y] for y in range(1, 8)], 2: [[3, 8]], 3: [[6, 0]]})
    iA = lat.chains[1].int_sequence[0]
    iB = lat.chains[2].int_sequence[0]
    ham.residue_interaction_table[iA, iB] = -8
    ham.residue_interaction_table[iB, iA] = -8
    E0 = int(ham.evaluate_total_energy(lat)[0])
    dims = lat.dimensions

    accepted = 0
    for trial in range(120):
        random.seed(trial)
        me, ok = sim.MOVER.cluster_rotate(lat.chains[1], lat, cluster_move_threshold=None,
                                          cluster_size_threshold=2, hardwall=True, frozen_chains=[])
        if not ok:
            continue
        accepted += 1
        # every committed bead is inside the box and the claimed energy change is real
        for cid, positions in me.moved_positions.items():
            for p in positions:
                assert all(0 <= p[d] < dims[d] for d in range(2)), (cid, p)
        dE_claimed = sim.rigid_cluster_move(me.moved_positions, me.original_positions)
        E1 = int(ham.evaluate_total_energy(lat)[0])
        assert dE_claimed == E1 - E0
        sim.rigid_cluster_revert(me.moved_positions, me.original_positions)
        assert int(ham.evaluate_total_energy(lat)[0]) == E0
    if expect_accepted:
        assert accepted > 20, accepted          # the 180-degree turn (1 draw in 3)
    else:
        # with the monomer at the corner every rotation about (3,3) leaves the box
        assert accepted == 0
    assert lat.check_grid_consistency() == []
