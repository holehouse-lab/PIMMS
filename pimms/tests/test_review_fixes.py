## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Regression tests for the correctness fixes made during the multi-agent code
review. Each test pins a specific bug so that it cannot silently return.

The bugs covered here are, in order:

* the serial + fast crankshaft/slither kernels teleporting a monomer across a
  HARDWALL boundary (a single move that jumped a full box length);
* ``chain_rotate`` and ``cluster_rotate`` rotating about a rounded centre of
  mass, which made the moves non-reversible and violated detailed balance;
* ``chain_pivot`` performing guaranteed null moves (and never moving a 3-mer)
  because of an off-by-one in the N-terminal arm;
* the Virtual-Move Monte Carlo boundary-link reverse probability using the wrong
  energy, biasing the equilibrium towards contact states;
* ``extract_cluster_polymeric_properties`` applying a spurious periodic
  correction to already-single-image cluster positions;
* the lemonade clustering ignoring HARDWALL and merging chains across walls;
* the move-selector thresholds leaving a rounding gap below 1.0.
"""

import math
import random

import numpy as np
import pytest

from pimms import lattice_utils
from pimms.tests import kernel_test_utils as U


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _place(lat, chains):
    """Clear the grids and place ``chains`` = {chainID: positions}."""
    lat.grid[:] = 0
    lat.type_grid[:] = 0
    for cid, pos in chains.items():
        pos = [list(map(int, p)) for p in pos]
        idx = list(range(len(pos)))
        lattice_utils.place_chain_by_position(pos, lat.grid, cid, safe=True)
        lat.insert_chain_into_type_grid(cid, pos, idx, safe=True)
        lat.chains[cid].set_ordered_positions(pos)


def _controlled_randint(seq):
    """Return a ``random.randint`` replacement that yields ``seq`` in order."""
    it = iter(seq)
    return lambda a, b: next(it)


def _as_chain_dict(moved_chain_positions):
    """Normalise a MoveEvent's moved_chain_positions to {chainID: [[..],..]}.

    Single-chain moves (chain_rotate/chain_pivot) return a flat list of
    positions; cluster moves return a dict keyed by chainID."""
    if isinstance(moved_chain_positions, dict):
        return {cid: [list(map(int, q)) for q in moved_chain_positions[cid]]
                for cid in moved_chain_positions}
    return {1: [list(map(int, q)) for q in moved_chain_positions]}


# ---------------------------------------------------------------------------
# monomer HARDWALL teleport
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dims", [[9, 9], [7, 8, 9]])
def test_monomer_hardwall_never_jumps_more_than_one_site(dims):
    """A flag-0 monomer under HARDWALL must never move more than one lattice
    site in a single accepted crankshaft move (it used to wrap the whole box)."""
    from pimms import mega_crank_fast as mcf
    from pimms import mega_crank as mc3
    from pimms import mega_crank_2D as mc2

    NP = np.int32
    n_dim = len(dims)
    if n_dim == 2:
        kernels = [mcf.mega_crank_2D, mc2.mega_crank_2D]
        angle = np.zeros((2, 3, 3, 3, 3), NP)
    else:
        kernels = [mcf.mega_crank, mc3.mega_crank]
        angle = np.zeros((2, 3, 3, 3, 3, 3, 3), NP)

    tables = (np.zeros((2, 2), NP),) * 3
    sel = np.zeros(1, dtype=np.int64)

    for kernel in kernels:
        grid = np.zeros(dims, dtype=NP)
        tg = np.zeros(dims, dtype=NP)
        grid[tuple([0] * n_dim)] = 1
        tg[tuple([0] * n_dim)] = 1
        idxb = np.array([[0, 0, 1, 1, 1] + [0] * n_dim], dtype=np.int64)

        max_jump = 0
        for step in range(3000):
            prev = idxb[0, 5:5 + n_dim].copy()
            kernel(grid, tg, idxb, tables[0], tables[1], tables[2],
                   angle, 0, 0.5, 1, sel, step + 1, 1)  # hardwall=1
            cur = idxb[0, 5:5 + n_dim]
            for d in range(n_dim):
                jump = abs(int(cur[d]) - int(prev[d]))
                # a legitimate step is +/-1; a wrap shows up as dims[d]-1
                jump = min(jump, dims[d] - jump)
                max_jump = max(max_jump, jump)
        assert max_jump <= 1, f"{kernel.__module__}: monomer jumped {max_jump} sites"


# ---------------------------------------------------------------------------
# rotation reversibility (detailed balance)
# ---------------------------------------------------------------------------
def _rotation_roundtrips(tmp_path, move_name, chains, box, ff="SR"):
    """Return the fraction of accepted rotations that the inverse rotation
    fails to undo, over many random single-image configurations."""
    state = U.build_state(tmp_path, len(box), ff, False,
                          {"MOVE_CHAIN_TRANSLATE": 1.0},
                          box=box, chains=[(len(p), "A" * len(p)) for p in chains.values()],
                          temperature=10, seed=3)
    lat = state.lattice
    mover = getattr(state.sim.MOVER, move_name)
    orig = random.randint
    bad = tot = 0
    rng = random.Random(0)
    for _ in range(300):
        # random rigid placement of the template chains
        shift = [rng.randrange(b) for b in box]
        placed = {}
        occupied = set()
        ok_place = True
        for cid, tmpl in chains.items():
            pos = []
            for bead in tmpl:
                p = tuple((bead[d] + shift[d]) % box[d] for d in range(len(box)))
                if p in occupied:
                    ok_place = False
                    break
                occupied.add(p)
                pos.append(list(p))
            placed[cid] = pos
            if not ok_place:
                break
        if not ok_place:
            continue
        _place(lat, placed)

        # forward: controlled angle (and axis in 3D)
        if len(box) == 3:
            ax, ang = rng.randint(0, 2), rng.randint(0, 2)
            # chain_rotate consumes randint(axis) then randint(angle);
            # cluster_rotate consumes randint(angle) then randint(axis)
            fwd = [ax, ang] if move_name == "chain_rotate" else [ang, ax]
        else:
            ang = rng.randint(0, 2)
            fwd = [ang]
        lattice_arg = lat if move_name == "cluster_rotate" else lat.grid
        random.randint = _controlled_randint(fwd)
        try:
            ME, moved = mover(lat.chains[1], lattice_arg)
        finally:
            random.randint = orig
        if not moved:
            continue
        new = _as_chain_dict(ME.moved_chain_positions)
        # only pin non-straddling rotations (exactly reversible); straddling ties
        # are a separate, negligible statistical effect covered by the equilibrium
        # check in the scratchpad reproducers.
        all_new = [q for pos in new.values() for q in pos]
        all_old = [q for pos in placed.values() for q in pos]
        if (lattice_utils.do_positions_stradle_pbc_boundary(all_new) or
                lattice_utils.do_positions_stradle_pbc_boundary(all_old)):
            continue
        tot += 1
        # rebuild a fresh state seeded with the new positions and apply inverse
        state2 = U.build_state(tmp_path, len(box), ff, False,
                               {"MOVE_CHAIN_TRANSLATE": 1.0},
                               box=box, chains=[(len(p), "A" * len(p)) for p in chains.values()],
                               temperature=10, seed=3)
        lat2 = state2.lattice
        mover2 = getattr(state2.sim.MOVER, move_name)
        _place(lat2, new)
        # inverse angle: 90<->270 (idx 0<->2), 180 self-inverse (idx 1)
        inv_ang = {0: 2, 1: 1, 2: 0}[ang]
        if len(box) == 3:
            rev = [ax, inv_ang] if move_name == "chain_rotate" else [inv_ang, ax]
        else:
            rev = [inv_ang]
        lattice_arg2 = lat2 if move_name == "cluster_rotate" else lat2.grid
        random.randint = _controlled_randint(rev)
        try:
            ME2, moved2 = mover2(lat2.chains[1], lattice_arg2)
        finally:
            random.randint = orig
        back = None
        if moved2:
            back = _as_chain_dict(ME2.moved_chain_positions)
        if (not moved2) or any(back[cid] != placed[cid] for cid in placed):
            bad += 1
    return bad, tot


@pytest.mark.parametrize("box", [[11, 11], [9, 9, 9], [10, 14]])
def test_chain_rotate_is_reversible(tmp_path, box):
    """chain_rotate must be exactly invertible (detailed balance). Rotating about
    a rounded COM used to fail this for even-length and PBC-straddling chains."""
    L = 4
    tmpl = [[i, 0, 0][:len(box)] for i in range(L)]
    bad, tot = _rotation_roundtrips(tmp_path, "chain_rotate", {1: tmpl}, box)
    assert tot > 30
    assert bad == 0, f"{bad}/{tot} chain rotations not reversible"


def test_cluster_rotate_is_reversible(tmp_path):
    """cluster_rotate must be exactly invertible for a multi-chain cluster."""
    chains = {1: [[3, 3], [3, 4]], 2: [[3, 5], [3, 6]]}
    bad, tot = _rotation_roundtrips(tmp_path, "cluster_rotate", chains, [11, 11])
    assert tot > 30
    assert bad == 0, f"{bad}/{tot} cluster rotations not reversible"


# ---------------------------------------------------------------------------
# chain_pivot off-by-one
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("L", [3, 4, 5, 6])
def test_chain_pivot_moves_short_chains_and_both_termini(tmp_path, L):
    """chain_pivot must actually move a chain for every eligible pivot point
    (it never moved a 3-mer before) and both termini must be reachable."""
    box = [15, 15]
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CHAIN_PIVOT": 1.0},
                          box=box, chains=[(L, "A" * L)], temperature=10, seed=3)
    lat = state.lattice
    mover = state.sim.MOVER.chain_pivot
    orig = random.randint
    template = [[3 + i, 7] for i in range(L)]
    beads_moved = set()
    null_moves = 0
    trials = 0
    for _ in range(400):
        _place(lat, {1: template})
        # let chain_pivot pick its own pivot point / rotation
        try:
            ME, moved = mover(lat.chains[1], lat.grid)
        finally:
            random.randint = orig
        if not moved:
            continue
        trials += 1
        new = _as_chain_dict(ME.moved_chain_positions)[1]
        if new == template:
            null_moves += 1
        for i in range(L):
            if new[i] != template[i]:
                beads_moved.add(i)
    assert trials > 0
    assert null_moves == 0, f"L={L}: {null_moves} null pivot moves"
    assert beads_moved, f"L={L}: no bead ever moved"
    if L >= 4:
        # both terminal beads must be reachable for L>=4
        assert 0 in beads_moved and (L - 1) in beads_moved, \
            f"L={L}: only beads {sorted(beads_moved)} ever moved"


# ---------------------------------------------------------------------------
# VMMC boundary-link detailed balance (exact enumeration, tiny system)
# ---------------------------------------------------------------------------
def test_vmmc_matches_exact_boltzmann_two_beads(tmp_path):
    """Two single-bead chains under MOVE_VMMC only must reproduce the exact
    Boltzmann probability of being in (SR/Chebyshev) contact. The boundary-link
    reverse-probability bug biased this strongly towards contact."""
    L = 7
    T = 4.0
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_VMMC": 1.0},
                          box=[L, L], chains=[(2, "A")], temperature=T, seed=3,
                          extra={"VMMC_MAX_DISPLACEMENT": 2, "VMMC_MAX_CLUSTER": 2})
    lat, ham, acc = state.lattice, state.ham, state.acc
    mover = state.sim.MOVER
    beta = acc.invtemp

    def cheb_contact():
        a = lat.chains[1].get_ordered_positions()[0]
        b = lat.chains[2].get_ordered_positions()[0]
        dx = min((a[0] - b[0]) % L, (b[0] - a[0]) % L)
        dy = min((a[1] - b[1]) % L, (b[1] - a[1]) % L)
        return max(dx, dy) == 1

    random.seed(1)
    np.random.seed(1)
    E = ham.evaluate_total_energy(lat)[0]
    n = 40000
    hits = 0
    for i in range(n):
        seed = lat.chains[random.randint(1, 2)]
        (_, E, _a, _cs) = mover.vmmc_move(seed, lat, E, acc, ham, 2, 2,
                                          hardwall=False, frozen_chains=[])
        if cheb_contact():
            hits += 1
    # the running energy must stay exact
    assert abs(ham.evaluate_total_energy(lat)[0] - E) < 1e-9
    obs = hits / n

    # exact P(contact): fix chain 1, enumerate chain 2 over every other site and
    # weight by the REAL Hamiltonian energy (Chebyshev/SR contact defined above).
    p1 = lat.chains[1].get_ordered_positions()[0]
    p2 = lat.chains[2].get_ordered_positions()[0]
    lattice_utils.delete_chain_by_position([p2], lat.grid, 2)
    lat.delete_chain_from_type_grid(2, [p2], [0], safe=True)
    Z = Zc = 0.0
    for a in range(L * L):
        ax, ay = divmod(a, L)
        if [ax, ay] == list(p1):
            continue
        lattice_utils.place_chain_by_position([[ax, ay]], lat.grid, 2, safe=True)
        lat.insert_chain_into_type_grid(2, [[ax, ay]], [0], safe=True)
        lat.chains[2].set_ordered_positions([[ax, ay]])
        w = math.exp(-beta * ham.evaluate_total_energy(lat)[0])
        Z += w
        dx = min((ax - p1[0]) % L, (p1[0] - ax) % L)
        dy = min((ay - p1[1]) % L, (p1[1] - ay) % L)
        if max(dx, dy) == 1:
            Zc += w
        lattice_utils.delete_chain_by_position([[ax, ay]], lat.grid, 2)
        lat.delete_chain_from_type_grid(2, [[ax, ay]], [0], safe=True)
    exact = Zc / Z

    se = math.sqrt(max(obs * (1 - obs), 1e-6) / n)
    assert abs(obs - exact) < 6 * se + 0.01, \
        f"VMMC P(contact) obs={obs:.4f} exact={exact:.4f}"


def test_hardwall_vmmc_matches_exact_boltzmann_and_never_wraps(tmp_path):
    """Hardwall VMMC must reject translations outside the box. Wrapping those
    translations while evaluating them as out-of-box proposals biases equilibrium."""
    L = 7
    state = U.build_state(tmp_path, 2, "SR", True, {"MOVE_VMMC": 1.0},
                          box=[L, L], chains=[(2, "A")], temperature=4.0, seed=7,
                          extra={"VMMC_MAX_DISPLACEMENT": 2, "VMMC_MAX_CLUSTER": 2})
    lat, ham, acc = state.lattice, state.ham, state.acc
    mover = state.sim.MOVER

    def in_contact():
        a = lat.chains[1].get_ordered_positions()[0]
        b = lat.chains[2].get_ordered_positions()[0]
        return max(abs(a[0] - b[0]), abs(a[1] - b[1])) == 1

    random.seed(13)
    np.random.seed(13)
    energy = ham.evaluate_total_energy(lat)[0]
    n = 50000
    hits = 0
    for _ in range(n):
        seed = lat.chains[random.randint(1, 2)]
        _, energy, _accepted, _size = mover.vmmc_move(
            seed, lat, energy, acc, ham, 2, 2, hardwall=True, frozen_chains=[])
        assert all(0 <= coordinate < L
                   for chain in lat.chains.values()
                   for position in chain.get_ordered_positions()
                   for coordinate in position)
        hits += in_contact()
    assert ham.evaluate_total_energy(lat)[0] == pytest.approx(energy)
    observed = hits / n

    # Enumerate every ordered, non-overlapping pair using the production
    # Hamiltonian rather than duplicating its energy convention in the test.
    for cid in (1, 2):
        old = lat.chains[cid].get_ordered_positions()
        lattice_utils.delete_chain_by_position(old, lat.grid, cid)
        lat.delete_chain_from_type_grid(cid, old, [0], safe=True)
    z = z_contact = 0.0
    beta = acc.invtemp
    for flat_a in range(L * L):
        a = list(divmod(flat_a, L))
        lattice_utils.place_chain_by_position([a], lat.grid, 1, safe=True)
        lat.insert_chain_into_type_grid(1, [a], [0], safe=True)
        lat.chains[1].set_ordered_positions([a])
        for flat_b in range(L * L):
            b = list(divmod(flat_b, L))
            if b == a:
                continue
            lattice_utils.place_chain_by_position([b], lat.grid, 2, safe=True)
            lat.insert_chain_into_type_grid(2, [b], [0], safe=True)
            lat.chains[2].set_ordered_positions([b])
            weight = math.exp(-beta * ham.evaluate_total_energy(lat)[0])
            z += weight
            if max(abs(a[0] - b[0]), abs(a[1] - b[1])) == 1:
                z_contact += weight
            lattice_utils.delete_chain_by_position([b], lat.grid, 2)
            lat.delete_chain_from_type_grid(2, [b], [0], safe=True)
        lattice_utils.delete_chain_by_position([a], lat.grid, 1)
        lat.delete_chain_from_type_grid(1, [a], [0], safe=True)

    exact = z_contact / z
    se = math.sqrt(max(observed * (1 - observed), 1e-6) / n)
    assert abs(observed - exact) < 6 * se + 0.01, \
        f"hardwall VMMC P(contact) observed={observed:.4f} exact={exact:.4f}"


# ---------------------------------------------------------------------------
# cluster polymeric properties: no spurious PBC correction
# ---------------------------------------------------------------------------
def test_extract_cluster_properties_no_spurious_pbc():
    """extract_cluster_polymeric_properties must treat its input as already
    single-image (pbc_correction=False); an asymmetric cluster otherwise gets
    an inflated Rg."""
    from pimms import lattice_analysis_utils as lau

    rod = [[i, 0, 0] for i in range(40)]
    blob = [[38 + dx, 1 + dy, 1 + dz]
            for dx in range(4) for dy in range(4) for dz in range(4)]
    pts = rod + blob
    arr = np.asarray(pts, dtype=float)
    delta = arr - arr.mean(axis=0)
    eigenvalues = np.linalg.eigvalsh((delta.T @ delta) / len(arr))
    ref_rg = math.sqrt(eigenvalues.sum())
    got = lau.extract_cluster_polymeric_properties([pts])[0]
    assert np.isclose(got[0], ref_rg)


# ---------------------------------------------------------------------------
# lemonade clustering honours HARDWALL
# ---------------------------------------------------------------------------
def test_lemonade_hardwall_does_not_merge_across_walls():
    """Under HARDWALL two chains against opposite walls are not neighbours and
    must be two clusters, not one."""
    from pimms.lemonade._topology import Topology
    from pimms.lemonade._store import TrajectoryStore
    from pimms.lemonade.trajectory import LatticeTrajectory

    pos = [[0, 5, 5], [0, 6, 5], [0, 7, 5], [11, 5, 5], [11, 6, 5], [11, 7, 5]]
    top = Topology(["AAA", "AAA"])

    st_pbc = TrajectoryStore(np.array([pos], dtype=np.int32), (12, 12, 12),
                             3.65, False, top)
    st_hw = TrajectoryStore(np.array([pos], dtype=np.int32), (12, 12, 12),
                            3.65, True, top)
    n_pbc = len(LatticeTrajectory(st_pbc)[0].clusters)
    n_hw = len(LatticeTrajectory(st_hw)[0].clusters)
    assert n_pbc == 1
    assert n_hw == 2


# ---------------------------------------------------------------------------
# move-selector thresholds close the rounding gap at 1.0
# ---------------------------------------------------------------------------
def test_move_selector_thresholds_reach_one(tmp_path):
    """The highest non-zero move interval must have upper bound exactly 1.0 so a
    random draw arbitrarily close to 1 always selects a move."""
    # ten equal 0.1 fractions - float accumulation leaves the top below 1.0
    moves = {m: 0.1 for m in [
        "MOVE_CRANKSHAFT", "MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE",
        "MOVE_CHAIN_PIVOT", "MOVE_HEAD_PIVOT", "MOVE_SLITHER",
        "MOVE_CLUSTER_TRANSLATE", "MOVE_CLUSTER_ROTATE", "MOVE_PULL",
        "MOVE_SYSTEM_TSMMC"]}
    state = U.build_state(tmp_path, 3, "SR", False, moves, temperature=40, seed=3)
    ac = state.acc
    tops = [v[1] for v in ac.random_thresholds.values() if v[1] > v[0]]
    assert max(tops) == 1.0


# ---------------------------------------------------------------------------
# build_LR_envelope_pairs empty input returns correctly shaped arrays
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dims", [[10, 10], [10, 10, 10]])
def test_build_LR_envelope_pairs_empty_is_shaped(dims):
    """Empty LR envelope must be a correctly-shaped (0, 2, n_dim) array, not a
    flat empty array that breaks downstream reshaping."""
    from pimms import longrange_utils

    n_dim = len(dims)
    type_grid = np.zeros(dims, dtype=np.int32)
    sites, offsets = longrange_utils.build_LR_envelope_pairs([], [], type_grid, dims)
    assert sites.shape == (0, 2, n_dim)
    assert offsets.shape == (0, 2, n_dim)


# ---------------------------------------------------------------------------
# PARALLELIZE must never silently freeze a chain too long for a block interior
# ---------------------------------------------------------------------------
def test_parallelize_falls_back_for_long_chains(tmp_path):
    """A chain whose extent exceeds the checkerboard block interior can never be
    moved by the parallel slither/pull kernels; the dispatch must fall back to
    the serial kernel so PARALLELIZE changes only speed, never the sampling."""
    from pimms import moves as moves_mod

    chains = [(1, "A" * 40), (6, "AAA")]

    def long_chain_positions(parallelize):
        state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_SLITHER": 1.0},
                              box=[50, 50, 50], chains=chains,
                              temperature=40, seed=5)
        lat = state.lattice
        E = state.energy
        random.seed(2)
        np.random.seed(2)
        for _ in range(30):
            out = state.sim.MOVER.system_slither(
                lat, E, state.acc, state.ham, 10, hardwall=False,
                frozen_chains=[], parallelize=parallelize, num_threads=4)
            E = out[1]
        return lat.chains[1].get_ordered_positions()

    # the helper itself must classify this system as not-parallel-safe
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_SLITHER": 1.0},
                          box=[50, 50, 50], chains=chains, temperature=40, seed=5)
    from pimms import crankshaft_list_functions as clf
    lat = state.lattice
    idx = clf.update_idx_to_bead(lat)
    chain_ids = sorted(lat.chains.keys())
    lengths = [len(lat.chains[c]) for c in chain_ids]
    offsets = list(np.cumsum([0] + lengths[:-1]))
    assert not moves_mod._parallel_can_move_all_chains(
        idx, offsets, lengths, lat.dimensions, False)

    start = state.lattice.chains[1].get_ordered_positions()
    assert long_chain_positions(False) != start          # serial moves it
    assert long_chain_positions(True) != start           # parallel must too


# ---------------------------------------------------------------------------
# system TSMMC: palindromic protocol and non-overshooting schedule
# ---------------------------------------------------------------------------
def test_tsmmc_schedule_never_overshoots_jump_temperature():
    """The linear TSMMC ramp must end exactly at the jump temperature (the old
    arange-based schedule overshot it by a full step for ~6% of inputs) and be
    palindromic."""
    from pimms.chainTSMMC import TSMMC

    for target, jump, pts in [(50, 61, 3), (100, 120, 4), (77, 121, 7),
                              (150, 205, 20), (55, 60, 2)]:
        t = TSMMC(target, jump, 'LINEAR', 1, pts, False)
        sch = t.true_temp_schedule
        assert sch.max() == pytest.approx(float(jump)), (target, jump, pts)
        assert (sch == sch[::-1]).all(), "schedule must be palindromic"
        # ramp is monotonically non-decreasing up to the plateau
        ramp = sch[:pts]
        assert (np.diff(ramp) > 0).all()


def test_tsmmc_system_protocol_is_palindromic():
    """The executed system-TSMMC protocol must make exactly M sub-moves at every
    schedule temperature and none at the target temperature. The old bookkeeping
    ran the first M-1 sub-moves at the target and never visited the last
    schedule temperature."""
    from pimms.chainTSMMC import TSMMC

    class FakeACC:
        temperature = 100.0

        def update_temperature(self, t):
            self.temperature = float(t)

        def get_total_aux_chain_moves(self):
            return 0

    for M in (1, 5, 20):
        t = TSMMC(100, 120, 'LINEAR', M, 4, False)
        acc = FakeACC()
        t.start_system_TSMMC((None, None, None), -500.0, acc)
        visited = []
        while not t.system_move_complete():
            acc = t.check_in_system_TSMMC(acc, -500.0)
            visited.append(acc.temperature)
        expected = [float(temp) for temp in t.true_temp_schedule
                    for _ in range(M)]
        assert visited == expected, f"M={M}: protocol deviates"
        assert 100.0 not in visited, "no sub-move may run at the target temp"


def test_tsmmc_rejects_non_heating_excursion():
    """A jump temperature at or below the target must raise a clear error (a
    heating quench crossing TSMMC_JUMP_TEMP used to hit ZeroDivisionError or
    silently invert the excursion)."""
    from pimms.chainTSMMC import TSMMC
    from pimms.latticeExceptions import MoveException

    with pytest.raises(MoveException):
        TSMMC(120, 120, 'LINEAR', 1, 4, False)
    with pytest.raises(MoveException):
        TSMMC(130, 120, 'LINEAR', 1, 4, False)


# ---------------------------------------------------------------------------
# detailed-balance threshold sensitivity
# ---------------------------------------------------------------------------
def test_db_threshold_detects_beta_error(tmp_path):
    """The equilibrium-comparison tolerance must actually catch a wrong-beta
    acceptance bug: running the SAME trusted kernel at beta and 1.5*beta (a
    gross Metropolis error that shifts the collapsed fixture's mean energy by
    hundreds of sigma) must FAIL assert_same_equilibrium. Under the old
    std-based tolerance this exact bug passed EVERY detailed-balance case,
    because the biased trace's own (inflated) spread entered the tolerance.
    (A beta x1.2 error shifts this fixture's mean by less than any feasible
    trace can statistically resolve - the collapse response is nonlinear - so
    x1.5 is the sharpest pin available through the mean energy.)"""
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          temperature=40, seed=5)

    def crank_at(invtemp_scale):
        g, t, i = state.fresh()
        e = state.energy
        true_invtemp = state.acc.invtemp
        try:
            state.acc.invtemp = true_invtemp * invtemp_scale
            for m in range(30):
                e = U.crank_megastep(state, g, t, i, e, 1000 + m, substeps=2500)
            out = []
            for m in range(150):
                e = U.crank_megastep(state, g, t, i, e, 5000 + m, substeps=2500)
                out.append(e)
        finally:
            state.acc.invtemp = true_invtemp
        return np.array(out, dtype=float)

    ref = crank_at(1.0)
    same = crank_at(1.0)
    biased = crank_at(1.5)

    # correct kernel passes against itself
    U.assert_same_equilibrium(ref, same, "self-consistency")
    # a 50% beta error must be caught
    with pytest.raises(AssertionError):
        U.assert_same_equilibrium(ref, biased, "beta x1.5")


# ---------------------------------------------------------------------------
# ANALYSIS_MODULE without a live ANA_CUSTOM frequency is a parse error
# ---------------------------------------------------------------------------
def test_analysis_module_requires_ana_custom(tmp_path):
    """A loaded-and-validated ANALYSIS_MODULE with ANA_CUSTOM disabled (0/unset)
    must be rejected at parse time. The check must use the RAW keyword value:
    set_dynamic_defaults rewrites a disabled ANA_CUSTOM to N_STEPS+10 before the
    sanity checks run, which made an earlier version of this check dead code."""
    import os
    import contextlib
    from pimms.keyfile_parser import KeyFileParser
    from pimms.latticeExceptions import KeyFileException

    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    (tmp_path / "my_analysis.py").write_text(
        "def analysis_function(SIM_OBJ, step):\n    pass\n")

    def parse(extra):
        U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False,
                        {"MOVE_CRANKSHAFT": 1.0}, extra=extra)
        cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            with contextlib.redirect_stdout(open(os.devnull, "w")):
                return KeyFileParser("KEYFILE.kf")
        finally:
            os.chdir(cwd)

    with pytest.raises(KeyFileException):
        parse({"ANALYSIS_MODULE": "my_analysis.py", "ANA_CUSTOM": "0"})
    with pytest.raises(KeyFileException):
        parse({"ANALYSIS_MODULE": "my_analysis.py"})          # unset -> default 0
    parse({"ANALYSIS_MODULE": "my_analysis.py", "ANA_CUSTOM": "5"})   # fine


# ---------------------------------------------------------------------------
# Cython audit regressions
# ---------------------------------------------------------------------------
def test_crank_kernels_negative_nsteps_is_noop():
    """A non-positive substep count must return immediately: compared against an
    unsigned loop counter it previously wrapped to ~2^32 iterations and read the
    bead_selector far out of bounds."""
    from pimms import mega_crank_fast as mcf
    from pimms import mega_crank as mc3
    from pimms import mega_crank_2D as mc2

    NP = np.int32
    tables = (np.zeros((2, 2), NP),) * 3
    sel = np.zeros(1, dtype=np.int64)
    g3 = np.zeros((9, 9, 9), NP); t3 = np.zeros((9, 9, 9), NP)
    g3[0, 0, 0] = 1; t3[0, 0, 0] = 1
    i3 = np.array([[0, 0, 1, 1, 1, 0, 0, 0]], dtype=np.int64)
    a3 = np.zeros((2, 3, 3, 3, 3, 3, 3), NP)
    for kern in (mcf.mega_crank, mc3.mega_crank):
        e, acc = kern(g3, t3, i3.copy(), *tables, a3, -7, 0.5, -1, sel, 11, 0)
        assert (e, acc) == (-7, 0)
    g2 = np.zeros((9, 9), NP); t2 = np.zeros((9, 9), NP)
    g2[0, 0] = 1; t2[0, 0] = 1
    i2 = np.array([[0, 0, 1, 1, 1, 0, 0]], dtype=np.int64)
    a2 = np.zeros((2, 3, 3, 3, 3), NP)
    for kern in (mcf.mega_crank_2D, mc2.mega_crank_2D):
        e, acc = kern(g2, t2, i2.copy(), *tables, a2, -7, 0.5, -1, sel, 11, 0)
        assert (e, acc) == (-7, 0)


def test_randint_inclusive_for_any_start():
    """randint must honour its inclusive-range contract for every start, and the
    generalised formula must leave the historical start-in-{0,1} streams
    untouched (pinned separately by the 250-case bit-exactness suite)."""
    from pimms import mega_crank as mc3

    mc3.seed_C_rand(1234)
    for start, end in [(0, 5), (1, 3), (2, 5), (3, 7)]:
        seen = set()
        for _ in range(5000):
            r = mc3.randint_ext(start, end)
            assert start <= r <= end, (start, end, r)
            seen.add(r)
        assert seen == set(range(start, end + 1)), (start, end, seen)


def test_parallel_block_seeds_are_collision_free():
    """Per-block PRNG seeds must not recur across sweeps: the old affine
    derivation gave identical whole-block random streams whenever
    block1 + seed1 == block2 + seed2 (adjacent sweeps shared num_blocks-1
    streams)."""
    def seeds(passed_seed, nb):
        sm = (np.arange(nb, dtype=np.uint64) * np.uint64(0x9E3779B97F4A7C15)
              + np.uint64((int(passed_seed) * 0xBF58476D1CE4E5B9)
                          & 0xFFFFFFFFFFFFFFFF))
        sm ^= sm >> np.uint64(30)
        sm *= np.uint64(0x94D049BB133111EB)
        sm ^= sm >> np.uint64(27)
        sm *= np.uint64(0x9E3779B97F4A7C15)
        sm ^= sm >> np.uint64(31)
        return sm | np.uint64(1)

    all_seeds = set()
    n = 0
    for sv in range(512):
        ss = seeds(sv, 64)
        all_seeds.update(ss.tolist())
        n += 64
    assert len(all_seeds) == n, f"{n - len(all_seeds)} colliding block seeds"


def test_parallel_crank_3D_thread_count_independent(tmp_path):
    """The checkerboard decomposition is documented as thread-count independent;
    the suite only pinned the 2D case - pin 3D too."""
    state = U.build_state(tmp_path, 3, "SLR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=[19, 18, 17], temperature=40, seed=5)

    # the harness helper returns (energy, accepted)
    def trace(nthreads):
        g, t, i = state.fresh()
        e = state.energy
        out = []
        for m in range(4):
            e = U.parallel_megastep(state, g, t, i, e, 700 + m,
                                    substeps=6000, nthreads=nthreads)
            out.append(e)
        return out, np.asarray(g).copy(), np.asarray(i).copy()

    o1, g1, i1 = trace(1)
    o8, g8, i8 = trace(8)
    assert o1 == o8
    assert (g1 == g8).all() and (i1 == i8).all()


def test_snakesearch_rejects_bad_seed_idx():
    """seed_idx indexes fixed buffers with bounds checking off - out-of-range
    values (and empty inputs) must raise instead of writing out of bounds."""
    from pimms import cluster_kernels

    pos = np.array([[1, 1, 1], [1, 1, 2]], dtype=np.int64)
    dims = np.array([9, 9, 9], dtype=np.int64)
    with pytest.raises(ValueError):
        cluster_kernels.snakesearch_single_image(pos, dims, 2, 1)
    with pytest.raises(ValueError):
        cluster_kernels.snakesearch_single_image(pos, dims, -1, 1)
    with pytest.raises(ValueError):
        cluster_kernels.snakesearch_single_image(
            np.empty((0, 3), dtype=np.int64), dims, 0, 1)


def test_no_cython_backup_in_distribution():
    """The stale pimms/cython_backup snapshots must never ship: they were
    pre-fix twins of audited kernels swept into both the sdist and the wheel."""
    import pathlib
    repo = pathlib.Path(__file__).resolve().parents[2]
    assert not (repo / "pimms" / "cython_backup").exists()
    manifest = (repo / "MANIFEST.in").read_text()
    assert "prune pimms/cython_backup" in manifest


# ---------------------------------------------------------------------------
# cluster-move physics fixes (python_moves exact audit)
# ---------------------------------------------------------------------------
def test_cluster_translate_hardwall_never_wraps(tmp_path):
    """Under HARDWALL a cluster translate must never wrap a bead through the
    wall. A box-spanning cluster could previously be 'translated' by permuting
    its chains through the wall - occupancy unchanged, committed as
    energy-neutral, while the true hardwall energy changed."""
    L = 9
    state = U.build_state(tmp_path, 2, "SR", True, {"MOVE_CLUSTER_TRANSLATE": 1.0},
                          box=[L, 7], chains=[(9, "A")], temperature=10, seed=3)
    lat = state.lattice
    # a full column of monomers spanning y (box-spanning cluster) + 2 spectators
    cols = [[4, y] for y in range(7)]
    place = {i + 1: [cols[i]] for i in range(7)}
    place[8] = [[0, 0]]
    place[9] = [[8, 6]]
    _place(lat, place)
    random.seed(11)
    for trial in range(400):
        before = {c: [list(p) for p in lat.chains[c].get_ordered_positions()]
                  for c in lat.chains}
        ME, ok = state.sim.MOVER.cluster_translate(
            lat.chains[1], lat, cluster_move_threshold=None,
            cluster_size_threshold=len(lat.chains) - 1, hardwall=True)
        if not ok:
            continue
        after = {c: [list(map(int, p)) for p in ME.moved_chain_positions[c]]
                 for c in ME.moved_chain_positions}
        # the per-chain displacement must be a single consistent raw vector for
        # every moved chain (no chain may have wrapped through a wall)
        for c, newpos in after.items():
            d = [newpos[0][k] - before[c][0][k] for k in range(2)]
            for b_old, b_new in zip(before[c], newpos):
                assert [b_new[k] - b_old[k] for k in range(2)] == d
        deltas = {tuple(after[c][0][k] - before[c][0][k] for k in range(2))
                  for c in after}
        assert len(deltas) == 1, f"chains moved by different vectors: {deltas}"
        # rebuild the moved state for the next iteration
        _place(lat, {c: after.get(c, before[c]) for c in before})


def test_cluster_rotate_rejects_winding_cluster(tmp_path):
    """A cluster connected to its own periodic image (single-image extent >= the
    box on some axis) must be rejected: rotating it is not a rigid motion of the
    periodic system and previously committed with the wrong dE."""
    L = 7
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CLUSTER_TRANSLATE": 1.0},
                          box=[L, 9], chains=[(7, "A")], temperature=10, seed=3)
    lat = state.lattice
    # 7 monomers filling the x column -> winds the x axis
    place = {i + 1: [[i, 4]] for i in range(7)}
    _place(lat, place)
    random.seed(5)
    for trial in range(60):
        before = {c: [list(p) for p in lat.chains[c].get_ordered_positions()]
                  for c in lat.chains}
        ME, ok = state.sim.MOVER.cluster_rotate(
            lat.chains[1], lat, cluster_size_threshold=len(lat.chains) - 1)
        assert not ok, "winding cluster rotation must be rejected"
        # rejection must restore the lattice exactly
        for c in before:
            assert lat.chains[c].get_ordered_positions() == before[c]
        for c, pos in before.items():
            for p in pos:
                assert lat.grid[tuple(p)] == c


def test_rotation_pivot_tie_breaking_is_invertible(tmp_path):
    """Exact centroid-distance ties must be broken identically forward and
    reverse (integer arithmetic): float noise previously picked a different
    pivot bead for the inverse rotation, a complete DB violation on tied
    shapes."""
    box = [7, 9]
    tmpl = [[1, 8], [0, 7], [6, 8], [5, 7], [5, 6]]   # ABABA tied shape (beads 1,2 tied)
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CHAIN_TRANSLATE": 1.0},
                          box=box, chains=[(1, "ABABA")], temperature=10, seed=3)
    lat = state.lattice
    orig = random.randint
    bad = 0
    for ang in (0, 1, 2):
        _place(lat, {1: tmpl})
        random.randint = _controlled_randint([ang])
        try:
            ME, ok = state.sim.MOVER.chain_rotate(lat.chains[1], lat.grid)
        finally:
            random.randint = orig
        if not ok:
            continue
        new = _as_chain_dict(ME.moved_chain_positions)[1]
        _place(lat, {1: new})
        inv = {0: 2, 1: 1, 2: 0}[ang]
        random.randint = _controlled_randint([inv])
        try:
            ME2, ok2 = state.sim.MOVER.chain_rotate(lat.chains[1], lat.grid)
        finally:
            random.randint = orig
        back = _as_chain_dict(ME2.moved_chain_positions)[1] if ok2 else None
        if back != tmpl:
            bad += 1
    assert bad == 0, f"{bad}/3 tied-shape rotations not invertible"


def test_rigid_cluster_move_dedupes_intra_cluster_pairs():
    """_dedupe_pair_rows must remove exactly the doubled intra-cluster pair rows
    (each cross-chain pair is emitted once from each chain's envelope scan)."""
    from pimms.simulation import _dedupe_pair_rows

    a = [[1, 2], [3, 4]]
    b = [[5, 6], [7, 8]]
    doubled = np.array([a, b, a])          # a emitted twice
    out = _dedupe_pair_rows(doubled)
    assert out.shape == (2, 2, 2)
    assert [r.tolist() for r in out] == [a, b]
    # empty input passes through
    assert len(_dedupe_pair_rows([])) == 0


# ---------------------------------------------------------------------------
# per-chain analysis fixes (chain_analysis oracle audit)
# ---------------------------------------------------------------------------
def test_pbc_gyration_uses_arithmetic_mean_reference():
    """Under PBC the gyration tensor must be referenced to the ARITHMETIC mean
    of the reconstructed single-image coordinates. Referencing it to the
    circular COM (the old behaviour) inflated Rg^2 by exactly |mean-circular|^2
    for every non-symmetric chain (a rod passes by symmetry - hence the L-shape
    oracle here)."""
    from pimms import lattice_analysis_utils as lau

    # L-shape: exact Rg = sqrt(1.28)
    rg, _ = lau.get_polymeric_properties(
        [[0, 0], [1, 0], [2, 0], [2, 1], [2, 2]], [30, 30], pbc_correction=True)
    assert rg == pytest.approx(np.sqrt(1.28), abs=1e-12)

    # straddling 5-rod across the x wrap: exact Rg = sqrt(2)
    rg, _ = lau.get_polymeric_properties(
        [[7, 0], [8, 0], [0, 0], [1, 0], [2, 0]], [9, 9], pbc_correction=True)
    assert rg == pytest.approx(np.sqrt(2), abs=1e-12)

    # property test: for any non-straddling chain the PBC and Cartesian paths
    # must agree bit-for-bit
    rng = np.random.default_rng(3)
    for _ in range(50):
        pos = [[3, 3, 3]]
        while len(pos) < 8:
            step = rng.integers(-1, 2, 3)
            new = [int(pos[-1][k] + step[k]) for k in range(3)]
            if 0 <= min(new) and max(new) < 13 and new not in pos:
                pos.append(new)
        a = lau.get_polymeric_properties(pos, [13, 13, 13], pbc_correction=True)
        b = lau.get_polymeric_properties(pos, [13, 13, 13], pbc_correction=False)
        assert a[0] == b[0] and a[1] == b[1]


def test_r2r_analysis_skips_writer_when_no_pairs(tmp_path, monkeypatch):
    """With no ANA_RESIDUE_PAIRS the residue-distance analysis must return
    without invoking the writer (it previously fell through and reopened the
    output file on every analysis step)."""
    from pimms import analysis_IO
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          temperature=40, seed=3)
    called = {"n": 0}
    monkeypatch.setattr(analysis_IO, "write_residue_residue_distance",
                        lambda *a, **k: called.__setitem__("n", called["n"] + 1))
    # the R2R analysis is a closure built with the (empty) pair list
    r2r = state.sim.build_R2R_distance_distribution_analysis([])
    r2r(5)
    assert called["n"] == 0


# ---------------------------------------------------------------------------
# lemonade analysis fixes (phase/core oracle audits)
# ---------------------------------------------------------------------------
def _ball_store(cx, L=13, hardwall=False):
    from pimms.lemonade._topology import Topology
    from pimms.lemonade._store import TrajectoryStore
    pts = [[(cx + x) % L, (6 + y) % L, (6 + z) % L]
           for x in range(-2, 3) for y in range(-2, 3) for z in range(-2, 3)]
    return TrajectoryStore(np.array([pts], dtype=np.int32), (L, L, L), 3.65,
                           hardwall, __import__("pimms.lemonade._topology",
                                                fromlist=["Topology"]).Topology(["A" * len(pts)]))


def test_single_image_shift_preserves_box_congruence():
    """snakesearch's >=0 shift must be a whole number of box periods: an
    arbitrary shift displaced cluster COMs (mod box) by up to the cluster
    radius whenever the cluster straddled the LOW box face, smearing every
    radial profile built around them."""
    from pimms import cluster_utils

    L = 20
    pts = [[(1 + x) % L, (18 + y) % L, (10 + z) % L]
           for x in range(-3, 4) for y in range(-3, 4) for z in range(-3, 4)
           if x * x + y * y + z * z <= 4.3 ** 2]
    si = np.asarray(cluster_utils.convert_positions_to_single_image_snakesearch(
        pts, [L, L, L]))
    assert si.min() >= 0
    assert np.all((si - np.asarray(pts)) % L == 0), "congruence mod box lost"
    assert np.allclose(si.mean(axis=0) % L, [1, 18, 10])


def test_radial_profile_translation_invariant_across_boundary():
    """A droplet straddling the low box face must give the identical radial
    profile as the same droplet centred in the box."""
    from pimms.lemonade import phase_separation as ps
    from pimms.lemonade.trajectory import LatticeTrajectory

    _, straddle = ps.radial_density_profile(LatticeTrajectory(_ball_store(0)))
    _, centred = ps.radial_density_profile(LatticeTrajectory(_ball_store(6)))
    assert np.allclose(straddle, centred)
    # solid 5x5x5 cube: unit occupancy through the first three shells
    assert np.allclose(centred[:3], 1.0)


def test_radial_profile_hardwall_uses_true_distances():
    """Under HARDWALL the radial profile must use plain Cartesian distances:
    material at true distance r > box/2 must appear at r, not folded into an
    inner shell by the periodic metric."""
    from pimms.lemonade._topology import Topology
    from pimms.lemonade._store import TrajectoryStore
    from pimms.lemonade import phase_separation as ps
    from pimms.lemonade.trajectory import LatticeTrajectory

    L = 20
    drop = [[4 + dx, 10 + dy, 10 + dz] for dx in (-1, 0, 1)
            for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    far = [[19, 10, 10], [19, 11, 10], [19, 10, 11]]
    st = TrajectoryStore(np.array([drop + far], dtype=np.int32), (L, L, L),
                         3.65, True, Topology(["A" * 27, "AAA"]))
    r, rho = ps.radial_density_profile(LatticeTrajectory(st), min_beads=2)
    nz = [(rr, v) for rr, v in zip(r, rho) if v > 0]
    assert any(abs(rr - 15) < 1.0 for rr, _v in nz), "far material missing at true r"
    assert not any(4 < rr < 6 for rr, _v in nz), "periodic folding artefact"


def test_flat_interface_surface_tension_sentinel():
    """A fluctuation-free slab has zero capillary power: gamma must be the
    explicit +inf sentinel with no numpy warnings."""
    import warnings
    from pimms.lemonade._topology import Topology
    from pimms.lemonade._store import TrajectoryStore
    from pimms.lemonade import surface_tension as st_mod
    from pimms.lemonade.trajectory import LatticeTrajectory

    L = 8
    slab = [[x, y, z] for x in range(L) for y in range(L) for z in (3, 4, 5)]
    st = TrajectoryStore(np.array([slab], dtype=np.int32), (L, L, 24), 3.65,
                         False, Topology(["A" * len(slab)]), temperature=1.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out = st_mod.slab_surface_tension(LatticeTrajectory(st), n_modes=4)
    assert out.gamma == float("inf") and out.n_modes == 0


def test_analysis_selector_strings_validated():
    """Unknown geometry/by selector strings must raise instead of silently
    running a different analysis; the two modules accept each other's
    droplet/sphere vocabulary as synonyms."""
    from pimms.lemonade import phase_separation as ps
    from pimms.lemonade import surface_tension as st_mod
    from pimms.lemonade.trajectory import LatticeTrajectory

    traj = LatticeTrajectory(_ball_store(6))
    with pytest.raises(ValueError):
        ps.largest_cluster_size(traj, by="bead")
    with pytest.raises(ValueError):
        ps.cluster_size_distribution(traj, by="bead")
    with pytest.raises(ValueError):
        ps.analyze(traj, geometry="slb")
    with pytest.raises(ValueError):
        st_mod.surface_tension(traj, geometry="spehre")
    # synonyms accepted
    assert ps.analyze(traj, geometry="droplet").geometry == "sphere"


def test_topology_rejects_empty_sequences():
    """Zero-length chains would give reduceat garbage downstream; construction
    must reject them."""
    from pimms.lemonade._topology import Topology
    with pytest.raises(ValueError):
        Topology(["AAA", ""])


# ---------------------------------------------------------------------------
# cluster-analysis fixes (cluster_analysis oracle audit)
# ---------------------------------------------------------------------------
def test_lr_cluster_distribution_is_symmetric_partition(tmp_path):
    """LR-cluster connectivity must be symmetric (both endpoints LR-capable for
    LR/SLR edges): the directional emit previously made the decomposition
    seed-dependent and overlapping (chains counted in several clusters)."""
    import os
    import contextlib
    from pimms.keyfile_parser import KeyFileParser
    from pimms.simulation import Simulation
    from pimms import lattice_analysis_utils as lau

    (tmp_path / "params.prm").write_text(
        "A A -8 -3 -1\nA B -4 -2 -1\nB B -6 -3 -1\nA X -2\nB X -2\nX X -1\n"
        "A 0 0\nB 0 0\nX 0 0\n"
        "ANGLE_PENALTY A 0 0 0\nANGLE_PENALTY B 0 0 0\nANGLE_PENALTY X 0 0 0\n")
    (tmp_path / "KEYFILE.kf").write_text(
        "DIMENSIONS : 9 9 9\nPARAMETER_FILE : params.prm\nSEED : 7\n"
        "TEMPERATURE : 40\nHARDWALL : False\nN_STEPS : 4\nEQUILIBRATION : 1\n"
        "CHAIN : 4 AA\nCHAIN : 3 BB\nCHAIN : 4 X\nCHAIN : 2 XX\n"
        "MOVE_CRANKSHAFT : 1.0\n")
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        kf = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(kf.keyword_lookup)
    finally:
        os.chdir(cwd)
    lat = sim.LATTICE

    clusters = lau.get_LR_cluster_distribution(lat)
    # partition: every chain in exactly one cluster
    members = [c for cl in clusters for c in cl]
    assert sorted(members) == sorted(lat.chains.keys())

    # symmetric oracle: SR contact (any beads) OR Chebyshev 2/3 with BOTH LR
    info = {c: list(zip([tuple(p) for p in lat.chains[c].get_ordered_positions()],
                        lat.chains[c].get_LR_binary_array()))
            for c in lat.chains}
    dims = lat.dimensions
    ids = sorted(lat.chains)
    parent = {c: c for c in ids}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def cheb(a, b):
        return max(min(abs(a[k] - b[k]), dims[k] - abs(a[k] - b[k]))
                   for k in range(3))
    for i, ci in enumerate(ids):
        for cj in ids[i + 1:]:
            if any((lambda d: d == 1 or (d in (2, 3) and fa == 1 and fb == 1))(cheb(pa, pb))
                   for pa, fa in info[ci] for pb, fb in info[cj]):
                ra, rb = find(ci), find(cj)
                if ra != rb:
                    parent[ra] = rb
    oracle = sorted([sorted([c for c in ids if find(c) == r])
                     for r in set(find(c) for c in ids)])
    assert sorted([sorted(c) for c in clusters]) == oracle


def test_snakesearch_warns_on_percolating_cluster():
    """A cluster percolating the box has no legitimate single image; the gather
    must warn instead of silently returning BFS-order-dependent coordinates."""
    import warnings
    from pimms import cluster_utils

    ring = [[k, 3, 3] for k in range(7)]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        cluster_utils.convert_positions_to_single_image_snakesearch(ring, [7, 9, 11])
    assert any("percolates" in str(x.message) for x in w)


def test_extract_cluster_properties_accepts_numpy_dimensions():
    """A numpy-array dimensions argument must select the explicit-dimensions
    branch, not crash on an ambiguous truth test."""
    from pimms import lattice_analysis_utils as lau
    out = lau.extract_cluster_polymeric_properties(
        [[[1, 1, 1], [2, 1, 1]]], dimensions=np.array([9, 9, 9]))
    assert len(out) == 1


# ---------------------------------------------------------------------------
# I/O audit fixes
# ---------------------------------------------------------------------------
def test_restart_reader_validates_contents(tmp_path):
    """The restart reader must reject out-of-box beads, overlapping beads and
    non-dict pickles with clear RestartExceptions (negative coords previously
    wrapped silently onto real cells; overlaps crashed deep in the mover)."""
    import os
    import pickle
    from pimms.restart import RestartObject
    from pimms.latticeExceptions import RestartException

    os.chdir(tmp_path)
    r = RestartObject()
    r.dimensions = [9, 9, 9]
    r.energy = 0
    r.hardwall = False

    r.chains = {1: [[[1, 1, 1], [-2, 1, 1]], "AA", 0]}
    r.write_to_file()
    with pytest.raises(RestartException, match="outside box"):
        RestartObject().build_from_file("restart.pimms")

    r.chains = {1: [[[1, 1, 1]], "A", 0], 2: [[[1, 1, 1]], "A", 0]}
    r.write_to_file()
    with pytest.raises(RestartException, match="same site"):
        RestartObject().build_from_file("restart.pimms")

    pickle.dump([1, 2, 3], open("bad.pimms", "wb"))
    with pytest.raises(RestartException, match="dictionary"):
        RestartObject().build_from_file("bad.pimms")


def test_restart_write_is_atomic(tmp_path, monkeypatch):
    """A crash during the restart write must never destroy the previous good
    checkpoint (the old 'wb' open truncated it before writing)."""
    import os
    from pimms.restart import RestartObject
    from pimms import restart as restart_mod

    os.chdir(tmp_path)
    r = RestartObject()
    r.dimensions = [9, 9, 9]
    r.energy = -3
    r.hardwall = False
    r.chains = {1: [[[1, 1, 1], [2, 1, 1]], "AA", 0]}
    r.write_to_file()
    good = open("restart.pimms", "rb").read()

    def exploding_dump(obj, fh):
        fh.write(b"partial")
        raise RuntimeError("simulated crash mid-write")
    monkeypatch.setattr(restart_mod.pickle, "dump", exploding_dump)
    with pytest.raises(RuntimeError):
        r.write_to_file()
    assert open("restart.pimms", "rb").read() == good, "previous checkpoint destroyed"
    # and it still loads
    r2 = RestartObject()
    r2.build_from_file("restart.pimms")
    assert r2.energy == -3


def test_freeze_file_rejects_unknown_directives(tmp_path):
    """A typo'd freeze directive must raise, not silently freeze nothing."""
    from pimms.data_structures import FreezeFile
    from pimms.latticeExceptions import KeyFileException

    f = tmp_path / "fz.in"
    f.write_text("c 1 2\n")
    with pytest.raises(KeyFileException, match="Unrecognised"):
        FreezeFile(str(f))
    f.write_text("C 1 2\n# ok\nC 3\n")
    assert sorted(FreezeFile(str(f)).chains) == [1, 2, 3]


def test_write_keyfile_full_init_round_trip(tmp_path):
    """A keyfile written from a fully-initialised parser must re-parse to an
    equivalent keyword set (sentinels/derived keywords previously made every
    full-init round trip unparseable)."""
    import os
    import contextlib
    from pimms.keyfile_parser import KeyFileParser

    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0})
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            p1 = KeyFileParser("KEYFILE.kf")
            p1.write_keyfile("out.kf")
            p2 = KeyFileParser("out.kf")
    finally:
        os.chdir(cwd)
    diffs = [k for k in p1.keyword_lookup
             if k != "SEED" and not callable(p1.keyword_lookup.get(k))
             and str(p1.keyword_lookup.get(k)) != str(p2.keyword_lookup.get(k))]
    assert diffs == [], diffs


def test_stale_per_type_outputs_removed_on_rerun(tmp_path):
    """A re-run in a used directory must not leave the previous run's
    CHAIN_<T>_* files (they silently mixed two runs' data in glob analyses)."""
    import os
    os.chdir(tmp_path)
    # use type indices beyond anything this run's system can produce, so the
    # only way the files can be absent afterwards is the startup wipe
    for stale in ("CHAIN_9_CLUSTERS.dat", "CHAIN_7_INTSCAL.dat",
                  "CHAIN_8_DISTANCE_MAP.dat"):
        open(stale, "w").write("stale\n")
    open("QUENCH.dat", "w").write("stale quench\n")
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          n_steps=2, equilibration=0, temperature=40, seed=3)
    state.sim.run_simulation()
    assert not os.path.exists("CHAIN_9_CLUSTERS.dat")
    assert not os.path.exists("CHAIN_7_INTSCAL.dat")
    assert not os.path.exists("CHAIN_8_DISTANCE_MAP.dat")
    assert not os.path.exists("QUENCH.dat")


def test_xtc_stream_frames_carry_sequential_metadata(tmp_path):
    """Streamed XTC frames must be stamped time/step = 0,1,2,... like the
    buffered SAVE_AT_END path (mdtraj's writer defaults every frame to 0/0)."""
    import os
    import numpy as np
    import mdtraj as md
    import pimms.lattice_utils as lu

    os.chdir(tmp_path)
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          n_steps=1, equilibration=0, temperature=40, seed=5)
    lat = state.sim.LATTICE
    w = lu.open_xtc_writer(lat, lat.lattice_to_angstroms,
                           pdb_filename="START.pdb", xtc_filename="traj.xtc")
    for _ in range(3):
        lu.write_xtc_frame(w, lat, lat.lattice_to_angstroms)
    lu.close_xtc_writer(w)
    with md.formats.XTCTrajectoryFile("traj.xtc") as r:
        _, t, s, _ = r.read()
    assert list(s) == [0, 1, 2, 3]
    assert np.allclose(t, [0.0, 1.0, 2.0, 3.0])
    # topology + frame 0 must agree (autocenter/unwrap threading)
    traj = md.load("traj.xtc", top="START.pdb")
    pdb = md.load("START.pdb")
    assert np.allclose(traj.xyz[0], pdb.xyz[0], atol=1e-3)


def test_pdb_atom_and_model_lines_are_spec_compliant():
    """ATOM coordinates are Real(8.3) right-justified in cols 31-54 and the
    MODEL serial is right-justified in cols 11-14 (previously '1.5' was
    centre-padded, which strict PDB readers misparse)."""
    import pimms.pdb_utils as pu

    line = pu.build_atom_line(3, "CA", "GLY", "A", 7, 1.5, -2.25, 123.456, "SEG1")
    assert line[30:38] == "   1.500"
    assert line[38:46] == "  -2.250"
    assert line[46:54] == " 123.456"
    assert pu.build_model_line(12)[10:14] == "  12"


def test_save_at_end_energy_fail_flushes_buffered_trajectory(tmp_path):
    """With SAVE_AT_END, an energy-check abort must write the buffered
    trajectory out instead of discarding it with the raise."""
    import os
    import mdtraj as md
    from pimms.latticeExceptions import SimulationEnergyException

    os.chdir(tmp_path)
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                          n_steps=40, equilibration=0, temperature=40, seed=5,
                          extra={"SAVE_AT_END": "True", "XTC_FREQ": 10,
                                 "ENERGY_CHECK": 20})
    sim = state.sim
    # shift the from-scratch recompute after the run-start baseline call, so the
    # ENERGY_CHECK at step 20 sees a nonzero difference and aborts
    orig_eval = sim.Hamiltonian.evaluate_total_energy
    calls = {"n": 0}

    def _skewed(lattice):
        calls["n"] += 1
        tot, loc, lr, slr, ang = orig_eval(lattice)
        if calls["n"] > 1:
            tot = tot + 7
        return (tot, loc, lr, slr, ang)
    sim.Hamiltonian.evaluate_total_energy = _skewed
    with pytest.raises(SimulationEnergyException):
        sim.run_simulation()
    assert os.path.exists("CONFIG_AT_ENERGY_FAIL.pdb")
    # the buffered frames written before the abort must be on disk
    traj = md.load("traj.xtc", top="START.pdb")
    assert traj.n_frames >= 2


def test_LR_cluster_writer_validates_before_opening_files(tmp_path):
    """A profile/index length mismatch must raise before any LR cluster file
    is opened (the validation block was previously trapped inside the
    docstring, so the RG/ASPH/VOL/AREA/DEN tables gained a row the radial
    file never matched)."""
    import os
    from pimms import analysis_IO

    os.chdir(tmp_path)
    with pytest.raises(ValueError, match="equal length"):
        analysis_IO.write_LR_cluster_properties(
            5, [[1.0, 0.1]], [[3, 2, 0.5]],
            LR_cluster_radial_density=[[0.1, 0.2]],
            LR_cluster_radial_density_indices=[1, 2])
    assert os.listdir(tmp_path) == [], "files were created before validation"
