"""A/B comparison of the reference vs optimized crankshaft kernel.

Builds a REAL system state from the two_phase_equilibrium_demo keyfile (the same
inputs pimms.moves.system_shake feeds to the kernel), then:

  1. CORRECTNESS - runs both kernels on independent copies with the same seed
     and asserts the returned (energy, accepted_moves) AND the mutated
     grid / type_grid / idx_to_bead are bit-for-bit identical.

  2. SPEED - times each kernel over several large megamoves and reports speedup.

Run from the repository root:

    python pimms/fast_kernels/benchmark.py

Optionally pass a keyfile dir and substeps:

    python pimms/fast_kernels/benchmark.py <demo_dir> <substeps> <megamoves>
"""
import os
import sys
import time
import copy

import numpy as np

# --- make the demo's relative paths (PARAMETER_FILE) resolve -----------------
HERE = os.path.dirname(os.path.abspath(__file__))      # pimms/fast_kernels
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))     # repo root
sys.path.insert(0, REPO_ROOT)

DEMO_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    REPO_ROOT, "demo_keyfiles", "two_phase_equilibrium_demo")
SUBSTEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
MEGAMOVES = int(sys.argv[3]) if len(sys.argv) > 3 else 20

from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation
from pimms import crankshaft_list_functions
import pimms.mega_crank as ref_kernel
import pimms.mega_crank_fast as fast_kernel


def build_state():
    """Construct the lattice/Hamiltonian/acceptance objects from the demo."""
    cwd = os.getcwd()
    os.chdir(DEMO_DIR)
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        sim = Simulation(keyfile.keyword_lookup)
    finally:
        os.chdir(cwd)

    lattice = sim.LATTICE
    ham = sim.Hamiltonian
    acc = sim.ACC

    # evaluate_total_energy returns (total, local, LR, SLR, angle)
    energy = int(ham.evaluate_total_energy(lattice)[0])
    hardwall_int = 1 if sim.hardwall else 0
    return sim, lattice, ham, acc, energy, hardwall_int


def kernel_inputs(lattice, ham, substeps):
    """Reproduce exactly what moves.system_shake passes into the kernel."""
    idx_to_bead = crankshaft_list_functions.update_idx_to_bead(lattice)
    num_beads = len(idx_to_bead)
    bead_selector = crankshaft_list_functions.bead_selector_constructor(
        num_beads, substeps, lattice, frozen_chains=[], safecheck=True)
    return idx_to_bead, bead_selector


def run_kernel(kernel, lattice, ham, acc, idx_to_bead, bead_selector,
               energy, substeps, seed, hardwall_int):
    """Run one kernel on private copies of the mutable state; return results."""
    grid = lattice.grid.copy()
    type_grid = lattice.type_grid.copy()
    idx = idx_to_bead.copy()

    new_energy, accepted = kernel.mega_crank(
        grid,
        type_grid,
        idx,
        ham.residue_interaction_table,
        ham.LR_residue_interaction_table,
        ham.SLR_residue_interaction_table,
        ham.angle_lookup,
        energy,
        acc.invtemp,
        substeps,
        bead_selector,
        seed,
        hardwall_int,
    )
    return new_energy, accepted, grid, type_grid, idx


def main():
    print(f"Demo dir : {DEMO_DIR}")
    print(f"Substeps : {SUBSTEPS}   Megamoves (timing): {MEGAMOVES}\n")

    sim, lattice, ham, acc, energy, hardwall_int = build_state()
    print(f"chains={lattice.get_number_of_chains()}  "
          f"dims={tuple(lattice.dimensions)}  hardwall={hardwall_int}  "
          f"invtemp={acc.invtemp:.6f}  E0={energy}\n")

    # ---------------- CORRECTNESS (bit-exact) ----------------
    print("=" * 70)
    print("CORRECTNESS: reference vs fast on identical inputs + seed")
    print("=" * 70)
    seed = 12345
    idx_to_bead, bead_selector = kernel_inputs(lattice, ham, SUBSTEPS)

    r_E, r_acc, r_grid, r_tg, r_idx = run_kernel(
        ref_kernel, lattice, ham, acc, idx_to_bead, bead_selector,
        energy, SUBSTEPS, seed, hardwall_int)
    f_E, f_acc, f_grid, f_tg, f_idx = run_kernel(
        fast_kernel, lattice, ham, acc, idx_to_bead, bead_selector,
        energy, SUBSTEPS, seed, hardwall_int)

    checks = {
        "energy": r_E == f_E,
        "accepted_moves": r_acc == f_acc,
        "grid": np.array_equal(r_grid, f_grid),
        "type_grid": np.array_equal(r_tg, f_tg),
        "idx_to_bead": np.array_equal(r_idx, f_idx),
    }
    print(f"  reference: E={r_E}  accepted={r_acc}")
    print(f"  fast     : E={f_E}  accepted={f_acc}")
    for name, ok in checks.items():
        print(f"   [{'PASS' if ok else 'FAIL'}] {name}")
    all_ok = all(checks.values())
    print(f"\n  --> {'BIT-EXACT MATCH' if all_ok else 'MISMATCH!'}\n")

    # ---------------- SPEED ----------------
    print("=" * 70)
    print(f"SPEED: {MEGAMOVES} megamoves x {SUBSTEPS} substeps each")
    print("=" * 70)

    def timeit(kernel):
        # fresh inputs each timing run so both kernels see the same workload
        total = 0.0
        e = energy
        for m in range(MEGAMOVES):
            idx_tb, bsel = kernel_inputs(lattice, ham, SUBSTEPS)
            grid = lattice.grid.copy()
            tg = lattice.type_grid.copy()
            idx = idx_tb.copy()
            t0 = time.perf_counter()
            kernel.mega_crank(grid, tg, idx,
                              ham.residue_interaction_table,
                              ham.LR_residue_interaction_table,
                              ham.SLR_residue_interaction_table,
                              ham.angle_lookup, e, acc.invtemp, SUBSTEPS,
                              bsel, 7 + m, hardwall_int)
            total += time.perf_counter() - t0
        return total

    # warm up
    timeit(ref_kernel)
    timeit(fast_kernel)

    t_ref = timeit(ref_kernel)
    t_fast = timeit(fast_kernel)

    print(f"  reference kernel : {t_ref:8.3f} s  "
          f"({1e6 * t_ref / (MEGAMOVES * SUBSTEPS):.3f} us/substep)")
    print(f"  fast kernel      : {t_fast:8.3f} s  "
          f"({1e6 * t_fast / (MEGAMOVES * SUBSTEPS):.3f} us/substep)")
    if t_fast > 0:
        print(f"\n  --> SPEEDUP: {t_ref / t_fast:.2f}x")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
