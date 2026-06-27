"""Validation + speed harness for the PARALLEL checkerboard kernel.

The parallel kernel is NOT bit-identical to the serial kernel (different RNG
stream, frozen-halo decomposition), so it is validated differently:

  1. ENERGY CONSISTENCY (the key correctness gate)
     After a parallel megamove, the energy it reports incrementally
     (base + sum of per-block dE) must equal a FROM-SCRATCH recomputation of the
     total energy of the mutated configuration. Any cross-block race/overlap
     would corrupt the increment and break this equality. Checked across several
     thread counts.

  2. BEAD CONSERVATION - no bead created/destroyed.

  3. SPEED / SCALING - serial-fast vs parallel at 1/2/4/8 threads on a large,
     spatially DISPERSED system (domain decomposition needs beads spread across
     blocks; a collapsed single droplet concentrates them and will not scale).

Run from the repository root:

    python pimms/fast_kernels/benchmark_parallel.py
"""
import os
import sys
import copy
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))      # pimms/fast_kernels
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))     # repo root
sys.path.insert(0, REPO_ROOT)

from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation
from pimms import crankshaft_list_functions
import pimms.mega_crank as ref_kernel
import pimms.mega_crank_fast as fk


def build_state(demo_dir):
    cwd = os.getcwd()
    os.chdir(demo_dir)
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        sim = Simulation(keyfile.keyword_lookup)
    finally:
        os.chdir(cwd)
    lattice = sim.LATTICE
    ham = sim.Hamiltonian
    acc = sim.ACC
    energy = int(ham.evaluate_total_energy(lattice)[0])
    hardwall_int = 1 if sim.hardwall else 0
    return sim, lattice, ham, acc, energy, hardwall_int


def recompute_energy(lattice, ham, grid, type_grid, idx):
    """Total energy of the mutated state, computed from scratch."""
    lat = copy.deepcopy(lattice)
    lat.grid = grid
    lat.type_grid = type_grid
    local_idx = 0
    for chainID in sorted(lat.chains.keys()):
        n = len(lat.chains[chainID].get_ordered_positions())
        lat.chains[chainID].set_ordered_positions(idx[local_idx:local_idx + n, 5:].tolist())
        local_idx += n
    return int(ham.evaluate_total_energy(lat)[0])


def correctness(demo_dir, substeps=10000):
    print("=" * 72)
    print(f"CORRECTNESS (energy-consistency)  -  {os.path.basename(demo_dir)}")
    print("=" * 72)
    sim, lattice, ham, acc, energy, hardwall_int = build_state(demo_dir)
    has_LR = bool(np.any(np.asarray(crankshaft_list_functions.update_idx_to_bead(lattice))[:, 1] == 1))
    print(f"  E0={energy}  hardwall={hardwall_int}  dims={tuple(lattice.dimensions)}  "
          f"layout={fk.parallel_layout_info(*lattice.dimensions, has_LR, 8)}")

    idx0 = crankshaft_list_functions.update_idx_to_bead(lattice)
    all_ok = True
    for nthreads in (1, 2, 4, 8):
        grid = lattice.grid.copy()
        tg = lattice.type_grid.copy()
        ix = idx0.copy()
        beads_before = int(np.count_nonzero(grid))

        E_ret, accepted = fk.mega_crank_parallel(
            grid, tg, ix,
            ham.residue_interaction_table, ham.LR_residue_interaction_table,
            ham.SLR_residue_interaction_table, ham.angle_lookup,
            energy, acc.invtemp, substeps, 4242, hardwall_int, nthreads)

        beads_after = int(np.count_nonzero(grid))
        E_recompute = recompute_energy(lattice, ham, grid, tg, ix)

        e_ok = (E_ret == E_recompute)
        b_ok = (beads_before == beads_after)
        all_ok = all_ok and e_ok and b_ok
        print(f"  threads={nthreads}: accepted={accepted:6d}  "
              f"E_incremental={E_ret:8d}  E_recomputed={E_recompute:8d}  "
              f"[{'OK' if e_ok else 'ENERGY MISMATCH'}] [{'beads OK' if b_ok else 'BEAD LOSS'}]")
    print(f"\n  --> {'ALL CONSISTENT' if all_ok else 'FAILURE'}\n")
    return all_ok


def detailed_balance(demo_dir, equilibrate=300, compare=80, substeps=10000, nthreads=8):
    """Detailed-balance check: the PARALLEL kernel must sample the same
    equilibrium as the serial kernel.

    IMPORTANT: energy-consistency (correctness() above) and thread-independence
    do NOT detect a detailed-balance violation - they held while the parallel
    kernel was driving a condensing droplet apart. So we test the ensemble
    directly: equilibrate with the (correct) serial kernel, then from that SAME
    config run both kernels and confirm the parallel one HOLDS the equilibrium
    energy rather than drifting.
    """
    print("=" * 72)
    print(f"DETAILED BALANCE (ensemble agreement)  -  {os.path.basename(demo_dir)}")
    print("=" * 72)
    sim, lattice, ham, acc, energy, hw = build_state(demo_dir)
    tables = (ham.residue_interaction_table, ham.LR_residue_interaction_table,
              ham.SLR_residue_interaction_table, ham.angle_lookup)
    idx0 = crankshaft_list_functions.update_idx_to_bead(lattice)

    def serial_step(g, t, i, e, seed):
        bsel = crankshaft_list_functions.bead_selector_constructor(
            len(i), substeps, lattice, frozen_chains=[], safecheck=True)
        return fk.mega_crank(g, t, i, *tables, e, acc.invtemp, substeps, bsel, seed, hw)[0]

    # 1) equilibrate with the serial kernel
    g, t, i, e = lattice.grid.copy(), lattice.type_grid.copy(), idx0.copy(), energy
    for m in range(equilibrate):
        e = serial_step(g, t, i, e, 1000 + m)
    print(f"  serial-equilibrated E = {e}  (start was {energy})")

    # 2) from that config, run each kernel and collect the energy trace
    cfg = (g.copy(), t.copy(), i.copy())

    def trace(parallel):
        gg, tt, ii, ee, out = cfg[0].copy(), cfg[1].copy(), cfg[2].copy(), e, []
        for m in range(compare):
            if parallel:
                ee = fk.mega_crank_parallel(gg, tt, ii, *tables, ee, acc.invtemp,
                                            substeps, 5000 + m, hw, nthreads)[0]
            else:
                ee = serial_step(gg, tt, ii, ee, 5000 + m)
            out.append(ee)
        return np.array(out)

    es, ep = trace(False), trace(True)
    drift = ep.mean() - es.mean()
    tol = 3.0 * max(es.std(), ep.std()) + 0.02 * abs(es.mean())
    ok = abs(drift) < tol
    print(f"  serial   continues at E = {es.mean():9.1f} +/- {es.std():5.1f}")
    print(f"  parallel continues at E = {ep.mean():9.1f} +/- {ep.std():5.1f}")
    print(f"  parallel drift from serial: {drift:+.1f}  (tol ~{tol:.0f})")
    print(f"\n  --> {'DETAILED BALANCE OK' if ok else 'DETAILED BALANCE VIOLATED'}\n")
    return ok


def speed(demo_dir, substeps=20000, megamoves=15):
    print("=" * 72)
    print(f"SPEED / SCALING  -  {os.path.basename(demo_dir)}")
    print("=" * 72)
    sim, lattice, ham, acc, energy, hardwall_int = build_state(demo_dir)
    has_LR = bool(np.any(np.asarray(crankshaft_list_functions.update_idx_to_bead(lattice))[:, 1] == 1))
    info = fk.parallel_layout_info(*lattice.dimensions, has_LR, 8)
    print(f"  dims={tuple(lattice.dimensions)}  beads={lattice.get_number_of_chains()} chains  "
          f"layout={info}\n")
    idx0 = crankshaft_list_functions.update_idx_to_bead(lattice)

    def time_serial():
        total = 0.0
        for m in range(megamoves):
            grid = lattice.grid.copy(); tg = lattice.type_grid.copy(); ix = idx0.copy()
            bsel = crankshaft_list_functions.bead_selector_constructor(
                len(ix), substeps, lattice, frozen_chains=[], safecheck=True)
            t0 = time.perf_counter()
            ref_kernel.mega_crank(grid, tg, ix,
                ham.residue_interaction_table, ham.LR_residue_interaction_table,
                ham.SLR_residue_interaction_table, ham.angle_lookup,
                energy, acc.invtemp, substeps, bsel, 7 + m, hardwall_int)
            total += time.perf_counter() - t0
        return total

    def time_fast_serial():
        total = 0.0
        for m in range(megamoves):
            grid = lattice.grid.copy(); tg = lattice.type_grid.copy(); ix = idx0.copy()
            bsel = crankshaft_list_functions.bead_selector_constructor(
                len(ix), substeps, lattice, frozen_chains=[], safecheck=True)
            t0 = time.perf_counter()
            fk.mega_crank(grid, tg, ix,
                ham.residue_interaction_table, ham.LR_residue_interaction_table,
                ham.SLR_residue_interaction_table, ham.angle_lookup,
                energy, acc.invtemp, substeps, bsel, 7 + m, hardwall_int)
            total += time.perf_counter() - t0
        return total

    def time_parallel(nthreads):
        total = 0.0
        for m in range(megamoves):
            grid = lattice.grid.copy(); tg = lattice.type_grid.copy(); ix = idx0.copy()
            t0 = time.perf_counter()
            fk.mega_crank_parallel(grid, tg, ix,
                ham.residue_interaction_table, ham.LR_residue_interaction_table,
                ham.SLR_residue_interaction_table, ham.angle_lookup,
                energy, acc.invtemp, substeps, 7 + m, hardwall_int, nthreads)
            total += time.perf_counter() - t0
        return total

    time_fast_serial(); time_parallel(4)  # warm up

    t_ref = time_serial()
    t_fs = time_fast_serial()
    print(f"  reference serial kernel : {t_ref:7.3f} s")
    print(f"  fast serial kernel      : {t_fs:7.3f} s   ({t_ref / t_fs:.2f}x vs reference)")
    for nt in (1, 2, 4, 8):
        t = time_parallel(nt)
        print(f"  parallel  ({nt} thread{'s' if nt > 1 else ' '})     : {t:7.3f} s   "
              f"({t_fs / t:.2f}x vs fast-serial)")
    print()


def main():
    two_phase = os.path.join(REPO_ROOT, "demo_keyfiles", "two_phase_equilibrium_demo")
    ok1 = correctness(two_phase)
    ok2 = correctness(os.path.join(HERE, "validation_system"), substeps=4000)

    # detailed-balance: the two_phase droplet is the case that exposed the
    # frozen-halo DB bug, so it is the key regression check.
    ok3 = detailed_balance(two_phase)

    large = os.path.join(HERE, "validation_large")
    if os.path.isdir(large):
        correctness(large, substeps=20000)
        speed(large)
    else:
        print(f"(skip speed: {large} not found - create it for a dispersed large-box test)")

    return 0 if (ok1 and ok2 and ok3) else 1


if __name__ == "__main__":
    sys.exit(main())
