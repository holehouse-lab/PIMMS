"""
Regression tests for the fixes from the complete-review pass of 1.0.8.

Every test here reproduces a defect found by the review with an independent
oracle (an exact enumeration, a hand calculation or the documented behaviour)
and pins the corrected behaviour.
"""
from __future__ import annotations

import contextlib
import itertools
import math
import os
import random

import numpy as np
import pytest

from pimms import CONFIG
from pimms import lattice_analysis_utils as lau
from pimms import lattice_utils, moves
from pimms.keyfile_parser import KeyFileParser, _quench_rung_count
from pimms.latticeExceptions import KeyFileException
from pimms.simulation import Simulation
from pimms.tests import kernel_test_utils as U


@pytest.fixture(autouse=True)
def _restore_global_rng():
    """Some tests reseed the process-global RNG (the moves draw from it); put it
    back so later tests do not become order dependent."""
    state = random.getstate()
    yield
    random.setstate(state)


# ---------------------------------------------------------------------------
# VMMC: frustrated internal links carry a reverse-probability factor
# ---------------------------------------------------------------------------

def test_vmmc_log_acceptance_uses_internal_probability_for_frustrated_links():
    """Every tested link contributes one factor; a failed link whose partner ended
    up inside the cluster uses the internal (relative-displacement) reverse
    probability, a failed link whose partner stayed outside uses the boundary one."""
    beta, dE = 0.125, 4.0
    cluster = {1, 2}
    formed = [(0.5, 0.25)]
    failed = [(2, 0.3, 0.9, 0.4),     # partner 2 inside the cluster: frustrated -> p_r = 0.4
              (3, 0.3, 0.9, 0.4)]     # partner 3 outside: boundary -> p_r = 0.9
    expected = (-beta * dE + math.log(0.25 / 0.5)
                + math.log((1 - 0.4) / (1 - 0.3))
                + math.log((1 - 0.9) / (1 - 0.3)))
    got = moves._vmmc_log_acceptance(beta, dE, cluster, formed, failed)
    assert got == pytest.approx(expected)

    # a reverse probability of exactly one on a failed link (or zero on a formed
    # link) makes the reverse realisation impossible -> None (reject)
    assert moves._vmmc_log_acceptance(beta, 0.0, cluster, [], [(3, 0.3, 1.0, 0.4)]) is None
    assert moves._vmmc_log_acceptance(beta, 0.0, cluster, [], [(2, 0.3, 0.9, 1.0)]) is None
    assert moves._vmmc_log_acceptance(beta, 0.0, cluster, [(0.5, 0.0)], []) is None


class _VMMCModel:
    """Exact enumeration of the VMMC recruitment for SR monomers on a 2D periodic
    lattice, written independently of moves.py from the definition of the move:
    LIFO queue, each unordered pair tested once from the popped chain's side in
    ascending partner order, link probability max(0, 1 - exp(-beta dE_link))."""

    def __init__(self, L, eps, beta):
        self.L, self.eps, self.beta = L, eps, beta

    def _cheb(self, a, b):
        L = self.L
        return max(min((a[k] - b[k]) % L, (b[k] - a[k]) % L) for k in range(2))

    def _pair_e(self, x, y):
        return 0.0 if x == y else (self.eps if self._cheb(x, y) == 1 else 0.0)

    def _lp(self, d):
        return 0.0 if d <= 0 else 1.0 - math.exp(-self.beta * d)

    def shift(self, p, d):
        return ((p[0] + d[0]) % self.L, (p[1] + d[1]) % self.L)

    def _energies(self, pos, m, off):
        pm = self.shift(pos[m], off)
        return {j: e for j in pos if j != m for e in [self._pair_e(pm, pos[j])] if e != 0.0}

    def _link(self, pos, m, j, dr):
        neg = (-dr[0], -dr[1])
        e0 = self._energies(pos, m, (0, 0)).get(j, 0.0)
        ef = self._energies(pos, m, dr).get(j, 0.0)
        er = self._energies(pos, m, neg).get(j, 0.0)
        return self._lp(ef - e0), self._lp(er - e0), self._lp(e0 - ef)

    def realisations(self, pos, seed, dr):
        neg = (-dr[0], -dr[1])

        def rec(queue, cluster, tested, links, prob):
            if not queue:
                yield prob, frozenset(cluster), list(links)
                return
            queue = list(queue)
            m = queue.pop()
            cand = set(self._energies(pos, m, (0, 0))) | set(self._energies(pos, m, dr)) | set(self._energies(pos, m, neg))
            pending = [j for j in sorted(cand) if frozenset((m, j)) not in tested]

            def branch(k, queue, cluster, tested, links, prob):
                if k == len(pending):
                    yield from rec(queue, cluster, tested, links, prob)
                    return
                j = pending[k]
                p_f, p_ri, p_rb = self._link(pos, m, j, dr)
                t2 = tested | {frozenset((m, j))}
                if p_f > 0.0:
                    c2, q2 = set(cluster), list(queue)
                    if j not in c2:
                        c2.add(j)
                        q2.append(j)
                    yield from branch(k + 1, q2, c2, t2, links + [(j, True, p_f, p_ri, p_rb)], prob * p_f)
                if p_f < 1.0:
                    yield from branch(k + 1, queue, cluster, t2, links + [(j, False, p_f, p_ri, p_rb)], prob * (1.0 - p_f))

            yield from branch(0, queue, cluster, tested, links, prob)

        yield from rec([seed], {seed}, set(), [], 1.0)

    def rate(self, pos, target, dr):
        """T(mu -> nu) for translating exactly `target` by dr, averaged over the
        uniform seed choice; the dr and cluster-cutoff draws are direction
        independent and omitted. dE = 0 for a whole-system translation."""
        assert target == frozenset(pos)
        tot = 0.0
        for seed in target:
            for prob, cluster, links in self.realisations(pos, seed, dr):
                if cluster != target:
                    continue
                lr = 0.0
                ok = True
                for (j, formed, p_f, p_ri, p_rb) in links:
                    if formed:
                        if p_ri <= 0.0:
                            ok = False
                            break
                        lr += math.log(p_ri) - math.log(p_f)
                    else:
                        p_r = p_ri if j in cluster else p_rb
                        if p_r >= 1.0:
                            ok = False
                            break
                        lr += math.log(1.0 - p_r) - math.log(1.0 - p_f)
                if ok:
                    tot += prob * min(1.0, math.exp(lr))
        return tot / len(pos)


def test_vmmc_transition_flows_balance_between_equal_energy_states(tmp_path, monkeypatch):
    """Three single-bead chains in an L-shaped contact triangle and the same
    triangle rigidly translated have identical energy, so detailed balance demands
    T(mu -> nu) == T(nu -> mu). This geometry produces frustrated links (a tested
    link fails, the partner is recruited through the third bead); the pre-fix
    boundary-only rule gave a forward/reverse ratio of 1.36 here (exact
    enumeration, confirmed by direct measurement on the move). Both directions
    must now agree with each other and with the exact enumeration."""
    L, T, eps = 9, 8.0, -8.0
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_VMMC": 1.0},
                          box=[L, L], chains=[(3, "A")], temperature=T, seed=3,
                          extra={"VMMC_MAX_DISPLACEMENT": 2, "VMMC_MAX_CLUSTER": 3})
    lat, ham, acc, mover = state.lattice, state.ham, state.acc, state.sim.MOVER
    assert ham.residue_interaction_table[1][1] == eps

    mu = {1: (4, 4), 2: (5, 5), 3: (5, 4)}
    dr = (-2, -2)
    model = _VMMCModel(L, eps, acc.invtemp)
    nu = {c: model.shift(p, dr) for c, p in mu.items()}
    full = frozenset(mu)
    t_model = model.rate(mu, full, dr)
    assert model.rate(nu, full, (2, 2)) == pytest.approx(t_model)

    # force |dr| = 2 on both axes (4 equiprobable sign choices, one of which is
    # the wanted vector) and a cutoff that never rejects, so the measured rate is
    # t_model / 4 and the test has power at a few thousand trials
    monkeypatch.setattr(moves.random, "randint", lambda a, b: 2)
    monkeypatch.setattr(mover, "_vmmc_draw_nc", lambda n_chains, max_cluster: 3)

    def place(chains):
        lat.grid[:] = 0
        lat.type_grid[:] = 0
        for cid, p in chains.items():
            pos = [list(p)]
            lattice_utils.place_chain_by_position(pos, lat.grid, cid, safe=True)
            lat.insert_chain_into_type_grid(cid, pos, [0], safe=True)
            lat.chains[cid].set_ordered_positions(pos)

    def snapshot():
        return {c: tuple(int(v) for v in lat.chains[c].get_ordered_positions()[0]) for c in lat.chains}

    def measure(start, target, n, seed):
        rng = random.Random(seed)
        random.seed(seed + 1)
        place(start)
        E0 = ham.evaluate_total_energy(lat)[0]
        hits = 0
        for _ in range(n):
            place(start)
            chain = lat.chains[rng.randint(1, 3)]
            (_, E, accepted, _cs) = mover.vmmc_move(chain, lat, E0, acc, ham, 2, 3, hardwall=False, frozen_chains=[])
            if accepted and snapshot() == target:
                assert E == E0
                hits += 1
        return hits

    n = 6000
    place(mu)
    e_mu = ham.evaluate_total_energy(lat)[0]
    place(nu)
    assert ham.evaluate_total_energy(lat)[0] == e_mu

    h_f = measure(mu, nu, n, 101)
    h_r = measure(nu, mu, n, 202)
    expected = n * t_model / 4.0
    se = math.sqrt(expected * (1 - t_model / 4.0))
    assert expected > 150                       # the test must have power
    assert abs(h_f - expected) < 4.5 * se, (h_f, expected, se)
    assert abs(h_r - expected) < 4.5 * se, (h_r, expected, se)
    z = (h_f - h_r) / math.sqrt(h_f + h_r)
    assert abs(z) < 4.0, (h_f, h_r, z)


# ---------------------------------------------------------------------------
# helpers for simulation-level tests
# ---------------------------------------------------------------------------


def _run_dir(path, keyfile="KEYFILE.kf"):
    """Parse and run a PIMMS keyfile inside `path`, silently."""
    cwd = os.getcwd()
    os.chdir(str(path))
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            kf = KeyFileParser(keyfile)
            sim = Simulation(kf.keyword_lookup)
            sim.run_simulation()
        return kf, sim
    finally:
        os.chdir(cwd)


def _parse_dir(path, keyfile="KEYFILE.kf"):
    cwd = os.getcwd()
    os.chdir(str(path))
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            return KeyFileParser(keyfile)
    finally:
        os.chdir(cwd)


# ---------------------------------------------------------------------------
# quench rung count on fractional temperature ramps
# ---------------------------------------------------------------------------

def test_quench_rung_count_is_robust_to_binary_float_round_off():
    # (1.0 - 0.7) / 0.1 = 3.0000000000000004 -> three rungs, not four
    assert _quench_rung_count(1.0 - 0.7, 0.1) == 3
    # 0.3 - 0.2 = 0.09999999999999998 -> one rung, not an overshoot
    assert _quench_rung_count(0.3 - 0.2, 0.1) == 1
    # genuinely fractional ratios still round up
    assert _quench_rung_count(0.85, 0.2) == 5
    assert _quench_rung_count(0.5, 0.25) == 2


def test_fractional_quench_ramp_sets_equilibration_to_the_documented_window(tmp_path):
    """1.0 -> 0.7 by 0.1 makes three temperature changes; with QUENCH_FREQ 2 and
    QUENCH_AS_EQUILIBRATION the documented window is (1 + 3) * 2 = 8 steps (the
    ceil on binary floats used to give 10)."""
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[9, 9, 9], chains=[(2, "AAB")], n_steps=30, equilibration=1,
                    extra={"QUENCH_RUN": "True", "QUENCH_START": 1.0, "QUENCH_END": 0.7,
                           "QUENCH_STEPSIZE": 0.1, "QUENCH_FREQ": 2,
                           "QUENCH_AS_EQUILIBRATION": "True"})
    kf = _parse_dir(tmp_path)
    assert kf.keyword_lookup["EQUILIBRATION"] == 8
    assert kf.quench_rungs == 3

    # a one-rung ramp whose range is a hair under the step size is valid
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[9, 9, 9], chains=[(2, "AAB")], n_steps=30, equilibration=1,
                    extra={"QUENCH_RUN": "True", "QUENCH_START": 0.3, "QUENCH_END": 0.2,
                           "QUENCH_STEPSIZE": 0.1, "QUENCH_FREQ": 2,
                           "QUENCH_AS_EQUILIBRATION": "True"})
    kf = _parse_dir(tmp_path)
    assert kf.keyword_lookup["EQUILIBRATION"] == 4

    # a genuine overshoot is still refused
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[9, 9, 9], chains=[(2, "AAB")], n_steps=30, equilibration=1,
                    extra={"QUENCH_RUN": "True", "QUENCH_START": 0.3, "QUENCH_END": 0.25,
                           "QUENCH_STEPSIZE": 0.1, "QUENCH_FREQ": 2,
                           "QUENCH_AS_EQUILIBRATION": "True"})
    with pytest.raises(KeyFileException, match="overshoots"):
        _parse_dir(tmp_path)


# ---------------------------------------------------------------------------
# cluster-rotate box rule runs on the FINAL boundary mode / box of a restart run
# ---------------------------------------------------------------------------

def _make_restart(path, hardwall, box):
    U.write_param_file(str(path / "params.prm"), "SR")
    U.write_keyfile(str(path / "KEYFILE.kf"), 3, hardwall, {"MOVE_CRANKSHAFT": 1.0},
                    box=box, chains=[(3, "AAB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0, "XTC_FREQ": 2})
    _run_dir(path)
    assert (path / "restart.pimms").exists()


def test_restart_override_hardwall_cannot_bypass_the_periodic_cluster_rotate_rule(tmp_path):
    src = tmp_path / "pbc"
    src.mkdir()
    _make_restart(src, hardwall=False, box=[10, 10, 12])
    run = tmp_path / "run"
    run.mkdir()
    U.write_param_file(str(run / "params.prm"), "SR")
    # keyfile says HARDWALL, restart is periodic, override takes the restart's
    # (periodic) mode -> a periodic non-cubic run with cluster rotation
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, True,
                    {"MOVE_CRANKSHAFT": 0.5, "MOVE_CLUSTER_ROTATE": 0.5},
                    box=[10, 10, 12], chains=[(3, "AAB")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": str(src / "restart.pimms"),
                           "RESTART_OVERRIDE_HARDWALL": "True"})
    with pytest.raises(KeyFileException, match="MOVE_CLUSTER_ROTATE"):
        _parse_dir(run)


def test_periodic_keyfile_with_hardwall_restart_override_is_accepted(tmp_path):
    """The converse: the keyfile says periodic but the restart's hardwall mode is
    taken, so the production run IS hardwall and cluster rotation is fine."""
    src = tmp_path / "hw"
    src.mkdir()
    _make_restart(src, hardwall=True, box=[10, 10, 12])
    run = tmp_path / "run"
    run.mkdir()
    U.write_param_file(str(run / "params.prm"), "SR")
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, False,
                    {"MOVE_CRANKSHAFT": 0.5, "MOVE_CLUSTER_ROTATE": 0.5},
                    box=[10, 10, 12], chains=[(3, "AAB")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": str(src / "restart.pimms"),
                           "RESTART_OVERRIDE_HARDWALL": "True"})
    kf = _parse_dir(run)
    assert kf.keyword_lookup["HARDWALL"] is True


# ---------------------------------------------------------------------------
# never-sampled INTSCAL: every chain type gets its sentinel file
# ---------------------------------------------------------------------------

def test_never_sampled_internal_scaling_writes_nothing_and_warns_for_every_chain_type(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[12, 12, 12], chains=[(3, "AABB"), (3, "BBAA")], n_steps=30,
                    equilibration=5, extra={"ANA_INTSCAL": 1000, "ANA_DISTMAP": 1000,
                                            "PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(tmp_path)
    # An analysis that never sampled writes no file at all, for EVERY chain type -
    # the bug this pins is that the never-sampled branch used to clear a loop-wide
    # flag, so only the first chain type was handled. SCALING_INFORMATION.dat used
    # to be exempted and written with -1 rows; it no longer is, because an all--1
    # file is byte-identical to what a well-sampled run of sub-26-bead chains
    # writes, so it could never have told the two apart.
    for t in (0, 1):
        for name in (CONFIG.OUTNAME_SCALING_INFORMATION,
                     CONFIG.OUTNAME_INTERNAL_SCALING,
                     CONFIG.OUTNAME_DMAP):
            assert not (tmp_path / ("CHAIN_%i_%s" % (t, name))).exists()
    log = (tmp_path / "log.txt").read_text() if (tmp_path / "log.txt").exists() else ""
    assert log.count("never sampled") >= 4          # 2 types x (INTSCAL + DISTMAP)


# ---------------------------------------------------------------------------
# stale outputs: the angle summary and eq_* files of a previous run
# ---------------------------------------------------------------------------

def test_angles_off_rerun_removes_the_previous_angle_summary(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(2, "AAB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(tmp_path)
    assert (tmp_path / CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL).exists()
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(2, "AAB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0, "ANGLES_OFF": "True"})
    _run_dir(tmp_path)
    assert not (tmp_path / CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL).exists()


def test_resized_rerun_without_save_eq_removes_stale_eq_files(tmp_path):
    common = dict(box=[14, 14, 14], chains=[(3, "AABB")], n_steps=8, equilibration=4)
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, True, {"MOVE_CRANKSHAFT": 1.0}, **common,
                    extra={"RESIZED_EQUILIBRATION": "8 8 8", "SAVE_EQ": "True", "XTC_FREQ": 2,
                           "PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(tmp_path)
    assert (tmp_path / "eq_traj.xtc").exists() and (tmp_path / "eq_START.pdb").exists()
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, True, {"MOVE_CRANKSHAFT": 1.0}, **common,
                    extra={"RESIZED_EQUILIBRATION": "8 8 8", "SAVE_EQ": "False", "XTC_FREQ": 2,
                           "PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(tmp_path)
    assert not (tmp_path / "eq_traj.xtc").exists()
    assert not (tmp_path / "eq_START.pdb").exists()


# ---------------------------------------------------------------------------
# AUTOCENTER on an odd box axis stays on the lattice
# ---------------------------------------------------------------------------

def test_center_positions_stays_on_the_lattice_for_odd_boxes():
    out = lattice_utils.center_positions([[1, 0, 0], [2, 0, 0], [3, 0, 0]], [9, 9, 9])
    assert all(float(v).is_integer() for pos in out for v in pos)
    assert [list(map(int, p)) for p in out] == [[3, 4, 4], [4, 4, 4], [5, 4, 4]]


# ---------------------------------------------------------------------------
# hardwall radial density profiles normalise by the sites that exist
# ---------------------------------------------------------------------------

def _shell_oracle(pts, dims, hardwall):
    com = np.rint(pts.mean(axis=0)).astype(int)
    occ = {tuple(p) for p in pts.tolist()}
    out = []
    for k in range(1, min(dims) // 2):
        shell = [s for s in itertools.product(*[range(com[d] - k, com[d] + k + 1) for d in range(3)])
                 if max(abs(s[d] - com[d]) for d in range(3)) == k
                 and (not hardwall or all(0 <= s[d] < dims[d] for d in range(3)))]
        out.append(sum(1 for s in shell if s in occ) / len(shell))
    return out


def test_hardwall_radial_density_profile_matches_in_box_shell_enumeration():
    # a 2-thick 9x9 slab wetting the x = 0 wall of a 12^3 box
    pts = np.array([[x, y, z] for x in (0, 1) for y in range(1, 10) for z in range(1, 10)])
    dims = [12, 12, 12]
    hw = lau.compute_cluster_radial_density_profile([pts], dims, hardwall=True)[0]
    pbc = lau.compute_cluster_radial_density_profile([pts], dims)[0]
    assert np.allclose(hw, _shell_oracle(pts, dims, True))
    assert np.allclose(pbc, _shell_oracle(pts, dims, False))
    assert hw[0] == 1.0                      # the first shell is completely full
    assert pbc[0] < 0.7                      # the periodic normalisation reads it as ~65 %


# ---------------------------------------------------------------------------
# lemonade: slab routing under hardwall, quench temperature, 2D inference
# ---------------------------------------------------------------------------

def test_hardwall_slab_away_from_the_walls_keeps_the_two_interface_fit():
    from pimms.lemonade import phase_separation as ps
    z = np.arange(40, dtype=float)
    dens = np.full(40, 0.55)
    dens[14:26] = 1.0
    fit = ps.fit_slab_profile(z, dens, hardwall=True)
    assert fit.success, fit.reason
    assert fit.rho_dense == pytest.approx(1.0, abs=0.05)
    assert fit.rho_dilute == pytest.approx(0.55, abs=0.05)
    assert fit.half_width == pytest.approx(6.0, abs=0.5)
    # a genuine wetting film (dense against z = 0) still takes the single-interface fit
    wet = np.full(40, 0.05)
    wet[:12] = 1.0
    fit_w = ps.fit_slab_profile(z, wet, hardwall=True)
    assert fit_w.success, fit_w.reason
    assert fit_w.rho_dense == pytest.approx(1.0, abs=0.05)
    assert fit_w.rho_dilute == pytest.approx(0.05, abs=0.05)


def test_lemonade_load_uses_quench_end_as_the_temperature_of_a_quench_run(tmp_path):
    import warnings
    from pimms import lemonade
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(2, "AAB")], n_steps=12, equilibration=1,
                    temperature=90,
                    extra={"QUENCH_RUN": "True", "QUENCH_START": 200, "QUENCH_END": 40,
                           "QUENCH_STEPSIZE": 40, "QUENCH_FREQ": 1,
                           "QUENCH_AS_EQUILIBRATION": "True", "XTC_FREQ": 2,
                           "PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(tmp_path)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        traj = lemonade.load(str(tmp_path / "traj.xtc"), str(tmp_path / "START.pdb"),
                             keyfile=str(tmp_path / "KEYFILE.kf"))
    assert traj.temperature == pytest.approx(40.0)
    assert any("QUENCH_END" in str(w.message) for w in rec)


def test_lemonade_load_keeps_a_flat_3d_configuration_three_dimensional(tmp_path):
    import mdtraj as md
    from pimms import lemonade
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue("A", chain)
    for _ in range(3):
        top.add_atom("CA", md.element.carbon, res)
    spacing = 0.365
    xyz = np.array([[[1, 1, 0], [2, 1, 0], [3, 1, 0]]], dtype=float) * spacing
    t = md.Trajectory(xyz, top, unitcell_lengths=np.array([[9, 9, 9]]) * spacing,
                      unitcell_angles=np.array([[90.0, 90.0, 90.0]]))
    t.save_pdb(str(tmp_path / "flat.pdb"))
    traj = lemonade.load(pdb=str(tmp_path / "flat.pdb"))
    assert traj.dimensions == (9, 9, 9)
    # and a real 2D box (z period of one lattice unit) is still recognised as 2D
    t2 = md.Trajectory(xyz, top, unitcell_lengths=np.array([[9, 9, 1]]) * spacing,
                       unitcell_angles=np.array([[90.0, 90.0, 90.0]]))
    t2.save_pdb(str(tmp_path / "flat2d.pdb"))
    traj2 = lemonade.load(pdb=str(tmp_path / "flat2d.pdb"))
    assert traj2.dimensions == (9, 9)


# ---------------------------------------------------------------------------
# PDB chain identifiers: one per chain type, 62 available
# ---------------------------------------------------------------------------

def test_pdb_chain_identifiers_are_distinct_for_28_chain_types(tmp_path):
    import string
    from pimms import lemonade
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[14, 14, 14], chains=[(1, "AB")] * 28, n_steps=2, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0, "XTC_FREQ": 1})
    _run_dir(tmp_path)
    ids = []
    for line in (tmp_path / "START.pdb").read_text().splitlines():
        if line.startswith("ATOM"):
            ids.append(line[21])
    distinct = list(dict.fromkeys(ids))
    assert len(distinct) == 28
    assert distinct[:26] == list(string.ascii_uppercase)
    assert distinct[26:] == ["a", "b"]
    traj = lemonade.load(str(tmp_path / "traj.xtc"), str(tmp_path / "START.pdb"))
    assert len(set(int(t) for t in traj.topology.chain_types)) == 28


# ---------------------------------------------------------------------------
# keyfile parsing: '~' in paths, residue pairs under a restart, round trip
# ---------------------------------------------------------------------------

def test_tilde_is_expanded_for_every_path_keyword(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    # a restart file and a freeze file living under ~
    src = home / "src"
    src.mkdir()
    U.write_param_file(str(src / "params.prm"), "SR")
    U.write_keyfile(str(src / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(2, "AABB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(src)
    (home / "freeze.txt").write_text("C 1\n")
    U.write_param_file(str(home / "ff.prm"), "SR")
    run = tmp_path / "run"
    run.mkdir()
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(2, "AABB")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": "~/src/restart.pimms", "FREEZE_FILE": "~/freeze.txt"})
    text = (run / "KEYFILE.kf").read_text().replace("PARAMETER_FILE : params.prm",
                                                     "PARAMETER_FILE : ~/ff.prm")
    (run / "KEYFILE.kf").write_text(text)
    kf = _parse_dir(run)
    assert kf.keyword_lookup["PARAMETER_FILE"] == str(home / "ff.prm")
    assert kf.keyword_lookup["RESTART_FILE"].dimensions is not None      # the restart loaded
    assert sorted(kf.keyword_lookup["FREEZE_FILE"].chains) == [1]


def test_residue_pair_bounds_are_checked_against_the_restart_chains(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    U.write_param_file(str(src / "params.prm"), "SR")
    U.write_keyfile(str(src / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(2, "AABB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(src)
    run = tmp_path / "run"
    run.mkdir()
    U.write_param_file(str(run / "params.prm"), "SR")
    # a leftover, shorter CHAIN line that the restart run ignores
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(1, "AA")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": str(src / "restart.pimms"),
                           "ANA_RESIDUE_PAIRS": "0 3", "ANA_INTER_RESIDUE": 1})
    kf = _parse_dir(run)
    assert kf.keyword_lookup["ANA_RESIDUE_PAIRS"] == [[0, 3]]
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[10, 10, 10], chains=[(1, "AA")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": str(src / "restart.pimms"),
                           "ANA_RESIDUE_PAIRS": "0 5", "ANA_INTER_RESIDUE": 1})
    with pytest.raises(KeyFileException):
        _parse_dir(run)


def test_write_keyfile_round_trips_a_restart_run_with_extra_chains(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    U.write_param_file(str(src / "params.prm"), "SR")
    U.write_keyfile(str(src / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[12, 12, 12], chains=[(2, "AABB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(src)
    run = tmp_path / "run"
    run.mkdir()
    U.write_param_file(str(run / "params.prm"), "SR")
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[12, 12, 12], chains=[(2, "AABB")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": str(src / "restart.pimms"), "EXTRA_CHAIN": "2 BB"})
    kf = _parse_dir(run)
    composition = sorted(tuple(c) for c in kf.keyword_lookup["CHAIN"])
    cwd = os.getcwd()
    os.chdir(str(run))
    try:
        kf.write_keyfile("out.kf")
    finally:
        os.chdir(cwd)
    kf2 = _parse_dir(run, "out.kf")
    assert sorted(tuple(c) for c in kf2.keyword_lookup["CHAIN"]) == composition


def test_restart_chains_are_loaded_in_ascending_chainid_order(tmp_path):
    """The per-chain analysis columns are written in sorted chainID order while
    the trajectory follows dict order, so a restart pickle whose CHAINS dict is
    not ascending must still be loaded ascending or the two disagree."""
    import pickle
    src = tmp_path / "src"
    src.mkdir()
    U.write_param_file(str(src / "params.prm"), "SR")
    U.write_keyfile(str(src / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[12, 12, 12], chains=[(3, "AABB")], n_steps=4, equilibration=1,
                    extra={"PRINT_FREQ": 1000, "ENERGY_CHECK": 0})
    _run_dir(src)
    with open(src / "restart.pimms", "rb") as fh:
        payload = pickle.load(fh)
    assert list(payload["CHAINS"]) == [1, 2, 3]
    payload["CHAINS"] = dict(reversed(list(payload["CHAINS"].items())))
    with open(src / "restart.pimms", "wb") as fh:
        pickle.dump(payload, fh)

    run = tmp_path / "run"
    run.mkdir()
    U.write_param_file(str(run / "params.prm"), "SR")
    U.write_keyfile(str(run / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[12, 12, 12], chains=[(3, "AABB")], n_steps=4, equilibration=1,
                    extra={"RESTART_FILE": str(src / "restart.pimms"), "ENERGY_CHECK": 0})
    kf = _parse_dir(run)
    cwd = os.getcwd()
    os.chdir(str(run))
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(kf.keyword_lookup)
    finally:
        os.chdir(cwd)
    assert list(sim.LATTICE.chains) == [1, 2, 3]


def test_summary_reports_the_frame_split_of_a_resized_run(tmp_path, capsys):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, True, {"MOVE_CRANKSHAFT": 1.0},
                    box=[14, 14, 14], chains=[(3, "AABB")], n_steps=16, equilibration=4,
                    extra={"RESIZED_EQUILIBRATION": "8 8 8", "SAVE_EQ": "True", "XTC_FREQ": 3})
    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        kf = KeyFileParser("KEYFILE.kf")
        kf.print_summary()
    finally:
        os.chdir(cwd)
    out = capsys.readouterr().out
    assert "Expected number of frames : 2 (eq_traj.xtc) + 5 (traj.xtc)" in out


def test_summary_reports_the_quench_step_magnitude_for_a_heating_run(tmp_path, capsys):
    # the stored step is negated for a heating ramp; the summary used to echo
    # that internal sign rather than the magnitude the keyfile asked for
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[12, 12, 12], chains=[(2, "AABB")], n_steps=60, equilibration=4,
                    temperature=10,
                    extra={"QUENCH_RUN": "True", "QUENCH_START": 10, "QUENCH_END": 40,
                           "QUENCH_STEPSIZE": 10, "QUENCH_FREQ": 5,
                           "QUENCH_AS_EQUILIBRATION": "True"})
    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        kf = KeyFileParser("KEYFILE.kf")
        kf.print_summary()
    finally:
        os.chdir(cwd)
    out = capsys.readouterr().out
    assert kf.keyword_lookup["QUENCH_STEPSIZE"] == -10
    assert "QUENCH STEP  : 10.00" in out


def test_quench_file_records_fractional_temperatures_exactly(tmp_path, monkeypatch):
    from pimms import analysis_IO
    monkeypatch.chdir(tmp_path)
    analysis_IO.write_quench_file(5, 0.975, -12.0)
    analysis_IO.write_quench_file(10, 50.0, -13.0)
    lines = (tmp_path / CONFIG.QUENCHFILE_NAME).read_text().splitlines()
    assert lines[0].split("\t")[1] == "0.975"
    assert float(lines[1].split("\t")[1]) == 50.0


# ---------------------------------------------------------------------------
# system-wide TSMMC excursions draw their sub-moves from the non-TSMMC moves
# ---------------------------------------------------------------------------

def test_tsmmc_excursion_sub_moves_follow_the_keyfile_without_crankshaft(tmp_path):
    """MOVE_SLITHER 0.5 / MOVE_SYSTEM_TSMMC 0.5 / MOVE_CRANKSHAFT 0: inside the
    excursions every sub-move must be a slither; a nested TSMMC draw used to be
    turned into a crankshaft megamove the keyfile never enabled."""
    _trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SR", False, {"MOVE_SLITHER": 0.5, "MOVE_SYSTEM_TSMMC": 0.5},
        n_steps=40, energy_check=10, return_sim=True)
    counts = sim.ACC.aux_chain_move_count
    assert counts[6] > 0                         # slither sub-moves happened
    assert counts[1] == 0                        # no crankshaft was ever executed
    assert counts[9] == counts[10] == counts[12] == 0


def test_excursion_draw_renormalises_the_non_tsmmc_fractions():
    from pimms.acceptance import AcceptanceCalculator
    kw = {m: 0.0 for m in ['MOVE_CRANKSHAFT', 'MOVE_CHAIN_TRANSLATE', 'MOVE_CHAIN_ROTATE', 'MOVE_CHAIN_PIVOT',
                           'MOVE_HEAD_PIVOT', 'MOVE_SLITHER', 'MOVE_CLUSTER_TRANSLATE', 'MOVE_CLUSTER_ROTATE',
                           'MOVE_CTSMMC', 'MOVE_MULTICHAIN_TSMMC', 'MOVE_PULL', 'MOVE_SYSTEM_TSMMC',
                           'MOVE_JUMP_AND_RELAX', 'MOVE_VMMC']}
    kw.update({'MOVE_CHAIN_TRANSLATE': 0.3, 'MOVE_CHAIN_ROTATE': 0.1, 'MOVE_SYSTEM_TSMMC': 0.6})
    acc = AcceptanceCalculator(10.0, kw)
    acc.auxillary_chain = True
    random.seed(7)
    draws = [acc.move_selector() for _ in range(20000)]
    assert set(draws) == {2, 3}
    frac = draws.count(2) / len(draws)
    assert abs(frac - 0.75) < 0.02              # 0.3 / (0.3 + 0.1)


# ---------------------------------------------------------------------------
# console noise and warnings
# ---------------------------------------------------------------------------

def test_multichain_tsmmc_accept_message_is_silenced_by_reduced_printing(tmp_path, capsys):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_MULTICHAIN_TSMMC": 1.0},
                    box=[10, 10, 10], chains=[(6, "AAAAAA")], n_steps=30, equilibration=1,
                    extra={"REDUCED_PRINTING": "True", "PRINT_FREQ": 1000, "ENERGY_CHECK": 0,
                           "TSMMC_JUMP_TEMP": 120, "TSMMC_STEP_MULTIPLIER": 10,
                           "TSMMC_NUMBER_OF_POINTS": 6})
    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        kf = KeyFileParser("KEYFILE.kf")
        sim = Simulation(kf.keyword_lookup)
        capsys.readouterr()
        sim.run_simulation()
    finally:
        os.chdir(cwd)
    out = capsys.readouterr().out
    assert "Multichain re-arrangement accepted" not in out


def test_percolation_warning_honours_the_long_range_gather_threshold():
    """A ring of beads two sites apart in a 12-box winds through the LR (Chebyshev
    3) gather - its single image spans 11 sites, one short of the box - and must
    be flagged; the same gap is a legitimate non-winding contact cluster."""
    import warnings
    from pimms import cluster_utils
    ring = [[x, 5, 5] for x in range(0, 12, 2)]          # 0,2,4,6,8,10: 10 -> 0 is 2 apart
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        cluster_utils.convert_positions_to_single_image_snakesearch(ring, [12, 12, 12], space_threshold=3)
    assert any("percolates" in str(w.message) for w in rec)
    line = [[x, 5, 5] for x in range(0, 11)]             # 0..10 contact rod: extent 11 < 12
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        cluster_utils.convert_positions_to_single_image_snakesearch(line, [12, 12, 12], space_threshold=1)
    assert not any("percolates" in str(w.message) for w in rec)


# ---------------------------------------------------------------------------
# parallel crank: every requested substep is attempted
# ---------------------------------------------------------------------------

def test_parallel_crank_attempts_all_requested_substeps(tmp_path):
    """With fewer substeps than blocks the floored per-block split used to give
    every block zero attempts, so the sweep did nothing while system_shake logged
    the full number as proposed. At a very high temperature every non-clashing
    proposal is accepted, so a handful of substeps must change the configuration."""
    from pimms import mega_crank_fast as mcf
    box = [40, 40, 40]
    state = U.build_state(tmp_path, 3, "LR", False, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(40, "AABB"), (40, "AAAA"), (20, "A")],
                          temperature=1e6, seed=5)
    layout = mcf.parallel_crank_layout_info(*box, True)
    assert layout["num_blocks"] > 1
    n_changed = 0
    for seed in range(6):
        g, t, i = state.fresh()
        before = np.asarray(g).copy()
        U.parallel_megastep(state, g, t, i, state.energy, 900 + seed,
                            substeps=layout["num_blocks"] - 1, nthreads=2)
        n_changed += int((np.asarray(g) != before).any())
    assert n_changed >= 4, n_changed


# ---------------------------------------------------------------------------
# cluster moves: proposal symmetry
# ---------------------------------------------------------------------------

def _place_monomers(lat, chains):
    lat.grid[:] = 0
    lat.type_grid[:] = 0
    for cid, p in chains.items():
        pos = [list(p)]
        lattice_utils.place_chain_by_position(pos, lat.grid, cid, safe=True)
        lat.insert_chain_into_type_grid(cid, pos, [0], safe=True)
        lat.chains[cid].set_ordered_positions(pos)


def test_cluster_translate_does_not_clash_on_a_cluster_mates_old_site(tmp_path, monkeypatch):
    """Two touching monomers translated by (+1, +1): bead 1 lands on bead 2's old
    site. Deleting the cluster chain by chain made that a spurious clash, so +v was
    refused while -v from the shifted state was accepted (a directional proposal).
    Both directions must now succeed and be exact inverses."""
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CLUSTER_TRANSLATE": 1.0},
                          box=[9, 9], chains=[(3, "A")], temperature=10, seed=3)
    lat, mover = state.lattice, state.sim.MOVER
    start = {1: (2, 2), 2: (3, 3), 3: (7, 7)}
    monkeypatch.setattr(moves.random, "randint", lambda a, b: 1)

    _place_monomers(lat, start)
    monkeypatch.setattr(moves.numpy_utils, "randneg", lambda x: x)          # offset (+1, +1)
    me, ok = mover.cluster_translate(lat.chains[1], lat, cluster_move_threshold=1,
                                     cluster_size_threshold=3, hardwall=False, frozen_chains=[])
    assert ok, "translation onto a cluster-mate's vacated site was refused"
    moved = {cid: tuple(pos[0]) for cid, pos in me.moved_chain_positions.items()}
    assert moved == {1: (3, 3), 2: (4, 4)}

    # the reverse move from the shifted state
    _place_monomers(lat, {1: (3, 3), 2: (4, 4), 3: (7, 7)})
    monkeypatch.setattr(moves.numpy_utils, "randneg", lambda x: -x)         # offset (-1, -1)
    me2, ok2 = mover.cluster_translate(lat.chains[1], lat, cluster_move_threshold=1,
                                       cluster_size_threshold=3, hardwall=False, frozen_chains=[])
    assert ok2
    assert {cid: tuple(pos[0]) for cid, pos in me2.moved_chain_positions.items()} == {1: (2, 2), 2: (3, 3)}
    # the grid is consistent with the committed positions
    assert int(lat.grid[2, 2]) == 1 and int(lat.grid[3, 3]) == 2 and int(lat.grid[4, 4]) == 0


def test_cluster_translate_reverts_the_whole_cluster_on_a_genuine_clash(tmp_path, monkeypatch):
    state = U.build_state(tmp_path, 2, "SR", False, {"MOVE_CLUSTER_TRANSLATE": 1.0},
                          box=[9, 9], chains=[(3, "A")], temperature=10, seed=3)
    lat, mover = state.lattice, state.sim.MOVER
    start = {1: (2, 2), 2: (3, 3), 3: (4, 4)}            # 3 is NOT in the cluster (threshold 2)
    _place_monomers(lat, start)
    monkeypatch.setattr(moves.random, "randint", lambda a, b: 1)
    monkeypatch.setattr(moves.numpy_utils, "randneg", lambda x: x)
    before = np.asarray(lat.grid).copy()
    _me, ok = mover.cluster_translate(lat.chains[1], lat, cluster_move_threshold=1,
                                      cluster_size_threshold=2, hardwall=False, frozen_chains=[])
    assert not ok
    assert (np.asarray(lat.grid) == before).all()


def test_hardwall_cluster_rotation_is_invertible_in_a_non_cubic_box(tmp_path, monkeypatch):
    """A 7-mer lying along x spans a 7x10 hardwall box. The periodic winding guard
    used to refuse its rotation (extent 7 >= 7) while accepting the rotation back
    from the vertical state (extent 7 < 10): a one-way move. Under a hardwall the
    guard does not apply; both directions must be accepted and be exact inverses."""
    state = U.build_state(tmp_path, 2, "SR", True, {"MOVE_CLUSTER_ROTATE": 1.0},
                          box=[7, 10], chains=[(1, "AAAAAAA")], temperature=10, seed=3)
    lat, mover = state.lattice, state.sim.MOVER
    rod = [[x, 4] for x in range(7)]

    def place(pos):
        lat.grid[:] = 0
        lat.type_grid[:] = 0
        lattice_utils.place_chain_by_position(pos, lat.grid, 1, safe=True)
        lat.insert_chain_into_type_grid(1, pos, list(range(len(pos))), safe=True)
        lat.chains[1].set_ordered_positions(pos)

    place(rod)
    monkeypatch.setattr(moves.random, "randint", lambda a, b: 0)            # 90 degrees
    me, ok = mover.cluster_rotate(lat.chains[1], lat, cluster_move_threshold=None,
                                  cluster_size_threshold=2, hardwall=True, frozen_chains=[])
    assert ok, "rotation of a box-spanning hardwall rod was refused"
    vertical = [list(p) for p in me.moved_chain_positions[1]]
    xs = {p[0] for p in vertical}
    assert len(xs) == 1 and sorted(p[1] for p in vertical) == list(range(min(p[1] for p in vertical), min(p[1] for p in vertical) + 7))
    assert all(0 <= p[0] < 7 and 0 <= p[1] < 10 for p in vertical)

    place(vertical)
    monkeypatch.setattr(moves.random, "randint", lambda a, b: 2)            # 270 degrees
    me2, ok2 = mover.cluster_rotate(lat.chains[1], lat, cluster_move_threshold=None,
                                    cluster_size_threshold=2, hardwall=True, frozen_chains=[])
    assert ok2, "the inverse rotation was refused"
    assert [list(p) for p in me2.moved_chain_positions[1]] == rod


# ---------------------------------------------------------------------------
# write_keyfile round trips a quench keyfile with QUENCH_AS_EQUILIBRATION False
# ---------------------------------------------------------------------------

def test_write_keyfile_keeps_a_false_quench_as_equilibration(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[9, 9, 9], chains=[(2, "AAB")], n_steps=30, equilibration=2,
                    extra={"QUENCH_RUN": "True", "QUENCH_START": 1.0, "QUENCH_END": 0.5,
                           "QUENCH_STEPSIZE": 0.25, "QUENCH_FREQ": 2,
                           "QUENCH_AS_EQUILIBRATION": "False"})
    kf = _parse_dir(tmp_path)
    assert kf.keyword_lookup["QUENCH_AS_EQUILIBRATION"] is False
    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        kf.write_keyfile("out.kf")
    finally:
        os.chdir(cwd)
    kf2 = _parse_dir(tmp_path, "out.kf")           # used to raise: quench keyword missing
    assert kf2.keyword_lookup["QUENCH_AS_EQUILIBRATION"] is False
    assert kf2.keyword_lookup["EQUILIBRATION"] == kf.keyword_lookup["EQUILIBRATION"]


def test_write_keyfile_keeps_the_freeze_file(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    (tmp_path / "freeze.txt").write_text("C 1\n")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[9, 9, 9], chains=[(2, "AAB")], n_steps=4, equilibration=1,
                    extra={"FREEZE_FILE": "freeze.txt"})
    kf = _parse_dir(tmp_path)
    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        kf.write_keyfile("out.kf")
    finally:
        os.chdir(cwd)
    assert "FREEZE_FILE" in (tmp_path / "out.kf").read_text()
    kf2 = _parse_dir(tmp_path, "out.kf")
    assert sorted(kf2.keyword_lookup["FREEZE_FILE"].chains) == sorted(kf.keyword_lookup["FREEZE_FILE"].chains)


@pytest.mark.parametrize("move", ["slither", "pull"])
def test_parallel_chain_kernels_attempt_all_requested_submoves(tmp_path, move):
    """8 chains in a 40^3 SR box give 27 chain blocks: one sub-move per chain
    (8 < 27) floored to zero attempts per block, so the megamove did nothing. The
    kernels now report the attempts they make, so this is pinned exactly: a sweep
    either finds no chain inside a block interior (and attempts nothing) or makes
    every requested attempt. At a very high temperature almost every proposal is
    accepted, so the sweeps that attempt must also move things."""
    from pimms import mega_crank_fast as mcf
    box = [40, 40, 40]
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_" + move.upper(): 1.0},
                          box=box, chains=[(8, "AAAA")], temperature=1e6, seed=5)
    assert mcf.parallel_layout_info(*box, False, 4)["num_blocks"] > 8
    offsets, lengths, homo = U.chain_meta(state.idx0)
    selector = np.arange(len(offsets), dtype=np.int32)     # one sub-move per chain
    kernel = mcf.mega_slither_parallel if move == "slither" else mcf.mega_pull_parallel
    n_attempting = n_accepted = 0
    for seed in range(20):
        g, t, i = state.fresh()
        _e, accepted, attempted = kernel(g, t, i, offsets, lengths, homo, selector,
                                         *state.tables, state.energy, state.acc.invtemp,
                                         700 + seed, state.hardwall_int, int(lengths.max()),
                                         2, U.frozen_bead_mask(i))
        assert attempted in (0, len(selector)), attempted
        if attempted:
            n_attempting += 1
            n_accepted += accepted
    assert n_attempting >= 5, n_attempting
    assert n_accepted >= n_attempting, (n_accepted, n_attempting)


# ---------------------------------------------------------------------------
# compiled single-image gather: box-independent lookups, identical results
# ---------------------------------------------------------------------------

def _python_gather(positions, dims, seed_idx, threshold):
    """Independent dict-based BFS gather (the reference convention: each newly
    reached bead is placed in the image nearest its BFS parent, then the whole
    cluster is shifted by whole periods so every coordinate is >= 0)."""
    from collections import deque
    pos = [tuple(int(v) for v in p) for p in positions]
    n_dim = len(dims)
    site = {}
    for i, p in enumerate(pos):
        site[p] = i
    si = {seed_idx: pos[seed_idx]}
    queue = deque([seed_idx])
    while queue:
        ref = queue.popleft()
        r = si[ref]
        rp = tuple(r[d] % dims[d] for d in range(n_dim))
        for off in itertools.product(range(-threshold, threshold + 1), repeat=n_dim):
            nb = tuple((rp[d] + off[d]) % dims[d] for d in range(n_dim))
            j = site.get(nb)
            if j is None or j in si:
                continue
            new = []
            for d in range(n_dim):
                delta = pos[j][d] - rp[d]
                if 2 * delta > dims[d]:
                    delta -= dims[d]
                elif 2 * delta < -dims[d]:
                    delta += dims[d]
                new.append(r[d] + delta)
            si[j] = tuple(new)
            queue.append(j)
    assert len(si) == len(pos)
    out = np.array([si[i] for i in range(len(pos))], dtype=np.int64)
    for d in range(n_dim):
        mn = out[:, d].min()
        if mn < 0:
            out[:, d] += dims[d] * ((-mn + dims[d] - 1) // dims[d])
    return out


def _random_wrapped_cluster(rng, n, dims, threshold):
    """A connected cluster grown by random steps of at most `threshold`, wrapped."""
    n_dim = len(dims)
    pts = [tuple(int(rng.integers(0, dims[d])) for d in range(n_dim))]
    seen = set(pts)
    while len(pts) < n:
        base = pts[rng.integers(0, len(pts))]
        cand = tuple((base[d] + int(rng.integers(-threshold, threshold + 1))) % dims[d] for d in range(n_dim))
        if cand not in seen:
            seen.add(cand)
            pts.append(cand)
    return np.array(pts, dtype=np.int64)


def _random_snake(rng, n, dims, threshold):
    """A self-avoiding random walk of `n` beads with steps of at most `threshold`
    per axis, wrapped into the box: its single-image extent typically exceeds
    half the box, so a wrong bead-index mapping in the lookup (rank instead of
    index) changes the placed images and is caught."""
    n_dim = len(dims)
    pts = [tuple(int(rng.integers(0, dims[d])) for d in range(n_dim))]
    seen = set(pts)
    tries = 0
    while len(pts) < n and tries < 100000:
        tries += 1
        step = tuple(int(rng.integers(-threshold, threshold + 1)) for _ in range(n_dim))
        if all(v == 0 for v in step):
            continue
        cand = tuple((pts[-1][d] + step[d]) % dims[d] for d in range(n_dim))
        if cand not in seen:
            seen.add(cand)
            pts.append(cand)
    assert len(pts) == n
    return np.array(pts, dtype=np.int64)


def _drifting_snake(rng, n, dims, threshold):
    """A self-avoiding walk that mostly steps +threshold along x, so its
    single-image extent along x is well over half the box (and it usually
    crosses the x boundary): the case where a rank-for-index lookup bug changes
    the placed images."""
    n_dim = len(dims)
    pts = [tuple(int(rng.integers(0, dims[d])) for d in range(n_dim))]
    seen = set(pts)
    tries = 0
    while len(pts) < n and tries < 100000:
        tries += 1
        step = [int(rng.integers(-threshold, threshold + 1)) for _ in range(n_dim)]
        if rng.random() < 0.85:
            step[0] = threshold
        if all(v == 0 for v in step):
            continue
        cand = tuple((pts[-1][d] + step[d]) % dims[d] for d in range(n_dim))
        if cand not in seen:
            seen.add(cand)
            pts.append(cand)
    assert len(pts) == n
    return np.array(pts, dtype=np.int64)


@pytest.mark.parametrize("dims,n,threshold,grid,shape", [
    ([120, 130, 110], 12, 1, False, "blob"),   # small cluster, large box: binary search
    ([120, 130, 110], 12, 3, False, "blob"),
    ([80, 150, 150], 60, 1, False, "drift"),   # x-spanning snake on the search path
    ([120, 150, 150], 30, 3, False, "drift"),
    ([9, 9, 9], 300, 1, True, "blob"),         # large cluster, small box: flat grid
    ([9, 9, 9], 200, 3, True, "blob"),
    ([12, 12, 12], 60, 1, True, "snake"),
    ([140, 150], 10, 1, False, "blob"),        # 2D, both paths, both shapes
    ([100, 2000], 80, 1, False, "drift"),
    ([11, 13], 90, 1, True, "blob"),
    ([11, 13], 40, 3, True, "snake"),
])
def test_compiled_gather_matches_reference_on_both_lookup_paths(dims, n, threshold, grid, shape):
    from pimms import cluster_kernels
    # pin which path the case exercises, so a retuned selection rule cannot
    # silently drop coverage of either lookup
    assert cluster_kernels.snakesearch_uses_grid(dims, n, threshold) is grid
    rng = np.random.default_rng(7)
    makers = {"blob": _random_wrapped_cluster, "snake": _random_snake, "drift": _drifting_snake}
    for trial in range(6):
        pts = makers[shape](rng, n, dims, threshold)
        if shape == "drift":
            ref0 = _python_gather(pts, dims, 0, threshold)
            assert ref0[:, 0].max() - ref0[:, 0].min() + 1 > dims[0] // 2   # spans over half the box
        seed = int(rng.integers(0, n))
        got = cluster_kernels.snakesearch_single_image(pts, np.asarray(dims, dtype=np.int64), seed, threshold)
        ref = _python_gather(pts, dims, seed, threshold)
        assert np.array_equal(np.asarray(got), ref), (dims, n, threshold, shape, trial)


def test_compiled_gather_selection_rule_keeps_production_clusters_on_the_grid():
    """The flat grid is faster until the box is a few tens of times larger than
    the walk's work; the rule must not push ordinary clusters onto the search."""
    from pimms import cluster_kernels
    assert cluster_kernels.snakesearch_uses_grid([60, 60, 60], 150, 3)      # LR cluster, 60-box
    assert cluster_kernels.snakesearch_uses_grid([100, 100, 100], 700, 3)   # LR cluster, 100-box
    assert cluster_kernels.snakesearch_uses_grid([100, 100, 100], 3000, 1)  # big contact cluster
    assert cluster_kernels.snakesearch_uses_grid([200, 200, 200], 4500, 1)  # bigger box, bigger cluster
    assert cluster_kernels.snakesearch_uses_grid([12, 12, 12], 10, 1)       # tiny box
    assert not cluster_kernels.snakesearch_uses_grid([200, 200, 200], 10, 1)
    assert not cluster_kernels.snakesearch_uses_grid([100, 100, 100], 12, 1)
    assert not cluster_kernels.snakesearch_uses_grid([300, 300, 300], 40, 3)
    # measured crossovers (search still 3-4x faster below them at t = 3, the grid
    # 1.6x faster above them for compact contact clusters at t = 1)
    assert not cluster_kernels.snakesearch_uses_grid([100, 100, 100], 20, 3)
    assert not cluster_kernels.snakesearch_uses_grid([150, 150, 150], 57, 3)
    assert cluster_kernels.snakesearch_uses_grid([150, 150, 150], 240, 3)
    assert not cluster_kernels.snakesearch_uses_grid([300, 300, 300], 315, 3)
    assert cluster_kernels.snakesearch_uses_grid([300, 300, 300], 1000, 3)
    assert cluster_kernels.snakesearch_uses_grid([150, 150, 150], 940, 1)
    assert cluster_kernels.snakesearch_uses_grid([200, 200, 200], 1717, 1)
    # an absurd bead count must return, not hang (a 32-bit shift once did)
    assert cluster_kernels.snakesearch_uses_grid([10, 10, 10], 2 ** 31 + 1, 1)


def test_compiled_gather_cost_does_not_scale_with_the_box():
    """Ten beads in a 200^3 box must take roughly the same time as in a 12^3 box
    (previously a full 8,000,000-site grid was allocated and filled per call)."""
    import time
    from pimms import cluster_kernels
    rng = np.random.default_rng(3)
    small = _random_wrapped_cluster(rng, 10, [12, 12, 12], 1)
    big = _random_wrapped_cluster(rng, 10, [200, 200, 200], 1)
    def timed(pts, dims):
        d = np.asarray(dims, dtype=np.int64)
        cluster_kernels.snakesearch_single_image(pts, d, 0, 1)
        t0 = time.perf_counter()
        for _ in range(200):
            cluster_kernels.snakesearch_single_image(pts, d, 0, 1)
        return time.perf_counter() - t0
    t_small, t_big = timed(small, [12, 12, 12]), timed(big, [200, 200, 200])
    # generous: the old kernel was ~150x slower here; the selection rule itself is
    # pinned exactly by the test above
    assert t_big < 50 * t_small + 0.1, (t_small, t_big)


# ---------------------------------------------------------------------------
# hardwall long-range envelopes never contain a pair across a wall
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dim", [2, 3])
def test_hardwall_envelope_pairs_are_the_periodic_pairs_minus_wall_crossings(tmp_path, dim):
    from pimms import longrange_utils
    box = [9, 9, 9][:dim]
    state = U.build_state(tmp_path, dim, "SLR", True, {"MOVE_CRANKSHAFT": 1.0},
                          box=box, chains=[(4, "AABB"), (4, "AAAA"), (6, "A")], temperature=10, seed=2)
    lat = state.lattice
    n_pairs = 0
    n_dropped = 0
    for cid, chain in lat.chains.items():
        pos = chain.get_ordered_positions()
        lr = chain.get_LR_binary_array()
        hw_LR, hw_SLR = longrange_utils.build_LR_envelope_pairs(pos, lr, lat.type_grid, box, hardwall=True)
        pb_LR, pb_SLR = longrange_utils.build_LR_envelope_pairs(pos, lr, lat.type_grid, box, hardwall=False)
        for hw, pb in ((hw_LR, pb_LR), (hw_SLR, pb_SLR)):
            def key(pair):
                return tuple(sorted(tuple(int(v) for v in p) for p in pair))
            hw_set = {key(p) for p in np.asarray(hw).reshape(-1, 2, dim)}
            pb_set = {key(p) for p in np.asarray(pb).reshape(-1, 2, dim)}
            # a periodic pair that is not a raw (unwrapped) neighbour crosses a wall
            in_box = {k for k in pb_set if max(abs(k[0][d] - k[1][d]) for d in range(dim)) <= 3}
            assert hw_set == in_box, (cid, hw_set ^ in_box)
            n_pairs += len(hw_set)
            n_dropped += len(pb_set) - len(in_box)
    assert n_pairs > 0
    assert n_dropped > 0                       # the fixture actually has wall-crossing pairs


# ---------------------------------------------------------------------------
# final-review pins
# ---------------------------------------------------------------------------

def test_percolation_warning_at_contact_threshold_also_needs_a_touching_pair():
    """A contact staircase from (0,0) to (5,3) in a 6x6 box spans every x site
    but no pair meets through the x face (the ends differ by 3 in y), so its
    single image is unambiguous; the threshold-1 path used to flag it on extent
    alone. A straight rod spanning the box does meet itself through the face."""
    import warnings
    from pimms import cluster_utils
    stairs = [[0, 0], [1, 0], [1, 1], [2, 1], [2, 2], [3, 2], [3, 3], [4, 3], [5, 3]]
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        cluster_utils.convert_positions_to_single_image_snakesearch(stairs, [6, 6], space_threshold=1)
    assert not any("percolates" in str(w.message) for w in rec)
    rod = [[x, 2] for x in range(6)]
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        cluster_utils.convert_positions_to_single_image_snakesearch(rod, [6, 6], space_threshold=1)
    assert any("percolates" in str(w.message) for w in rec)


def test_percolating_axes_reports_every_confirmed_axis():
    """The per-axis form of the same test: the staircase wraps nothing, a rod
    wraps its own axis, a slab that fills the periodic plane wraps both of those
    axes and not the third; cluster_percolates is its first entry."""
    from pimms import cluster_utils
    stairs = [[0, 0], [1, 0], [1, 1], [2, 1], [2, 2], [3, 2], [3, 3], [4, 3], [5, 3]]
    assert cluster_utils.percolating_axes(stairs, [6, 6]) == []
    assert cluster_utils.cluster_percolates(stairs, [6, 6]) is None
    rod = [[x, 2] for x in range(6)]
    assert cluster_utils.percolating_axes(rod, [6, 6]) == [0]
    assert cluster_utils.cluster_percolates(rod, [6, 6]) == 0
    slab = [[x, y, z] for x in range(8) for y in range(8) for z in range(8, 16)]
    assert cluster_utils.percolating_axes(slab, [8, 8, 24]) == [0, 1]
    assert cluster_utils.cluster_percolates(slab, [8, 8, 24]) == 0
    assert cluster_utils.percolating_axes([], [6, 6]) == []
    # first_only stops at the first confirmed axis
    assert cluster_utils.percolating_axes(slab, [8, 8, 24], first_only=True) == [0]
    assert cluster_utils.percolating_axes(stairs, [6, 6], first_only=True) == []


def test_cluster_percolates_keeps_its_first_axis_early_exit(monkeypatch):
    """cluster_percolates only ever reports one axis, so it must not pay the
    per-axis pair search on the axes after the first confirmed one; the all-axes
    form visits every axis."""
    from pimms import cluster_utils
    calls = []
    real = cluster_utils._axis_percolates

    def counting(arr, dimensions, d, *args, **kwargs):
        calls.append(d)
        return real(arr, dimensions, d, *args, **kwargs)

    monkeypatch.setattr(cluster_utils, "_axis_percolates", counting)
    slab = [[x, y, z] for x in range(8) for y in range(8) for z in range(8, 16)]
    assert cluster_utils.cluster_percolates(slab, [8, 8, 24]) == 0
    assert calls == [0]
    calls.clear()
    assert cluster_utils.percolating_axes(slab, [8, 8, 24]) == [0, 1]
    assert calls == [0, 1, 2]


def test_percolation_warning_requires_a_pair_that_touches_through_the_face():
    """At threshold 3 the per-axis extent test is necessary but not sufficient:
    beads at x = 0, 3, 6, 9 with y = 0, 2, 4, 6 span 10 of a 12-box (>= 12 - 3 + 1)
    but no pair is within 3 on every axis through the boundary, so the cluster
    does not wind and must not be flagged."""
    import warnings
    from pimms import cluster_utils
    pts = [[0, 0, 5], [3, 2, 5], [6, 4, 5], [9, 6, 5]]
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        cluster_utils.convert_positions_to_single_image_snakesearch(pts, [12, 12, 12], space_threshold=3)
    assert not any("percolates" in str(w.message) for w in rec)


def test_fit_report_measures_the_extent_on_unsplit_axes(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_SLITHER": 1.0},
                          box=[60, 60, 20], chains=[(10, "AABBAABB")], temperature=10, seed=2)
    idx, chains, off, length, homo = moves.parallel_chain_metadata(state.lattice)
    rep = moves.parallel_chain_fit_report(idx, off, length, [60, 60, 20], False,
                                          chain_homo=homo, cap_mode="hetero")
    assert rep["layout"]["blocks"][2] == 1                  # z is not split
    assert all(e >= 1 for e in rep["max_extent_by_axis"])   # z extent is measured, not 0


def test_parallel_report_after_resize_names_the_production_boundary(tmp_path):
    U.write_param_file(str(tmp_path / "params.prm"), "SR")
    U.write_keyfile(str(tmp_path / "KEYFILE.kf"), 3, False, {"MOVE_CRANKSHAFT": 1.0},
                    box=[36, 36, 36], chains=[(10, "AABB")], n_steps=12, equilibration=10,
                    extra={"PARALLELIZE": "True", "PARALLEL_THREADS": 2,
                           "RESIZED_EQUILIBRATION": "16 16 16",
                           "ENERGY_CHECK": 0, "XTC_FREQ": 1000, "PRINT_FREQ": 1000})
    _run_dir(tmp_path)
    log = (tmp_path / "log.txt").read_text()
    second = log.split("after resized equilibration (production box)")[1]
    assert "36x36x36 (3D, periodic)" in second
    first = log.split("PARALLELIZATION REPORT")[1]
    assert "16x16x16 (3D, hardwall)" in first


def test_quench_overshoot_message_keeps_the_fractional_temperatures(capsys):
    from pimms import nonequilibrium_utils
    new = nonequilibrium_utils.update_temperature_in_quench(0.3, 0.4, 0.2, 0.4)
    assert new == 0.2
    out = capsys.readouterr().out
    assert "0.400" in out and "0.200" in out and "from 0 to 0" not in out


def test_energy_neutral_cluster_move_reports_an_integer_energy_change(tmp_path):
    state = U.build_state(tmp_path, 3, "SR", False, {"MOVE_CLUSTER_TRANSLATE": 1.0},
                          box=[12, 12, 12], chains=[(6, "AABB")], temperature=20, seed=3)
    sim, lat = state.sim, state.lattice
    random.seed(4)
    for _ in range(300):
        chain = lat.get_random_chain(frozen_chains=[])
        me, ok = sim.MOVER.cluster_translate(chain, lat, cluster_move_threshold=None,
                                             cluster_size_threshold=lat.get_number_of_chains() - 1,
                                             hardwall=False, frozen_chains=[])
        if ok:
            delta = sim.rigid_cluster_move(me.moved_positions, me.original_positions)
            assert isinstance(delta, (int, np.integer)) and not isinstance(delta, bool)
            break
    else:
        pytest.skip("no cluster translation was accepted in 300 attempts")
