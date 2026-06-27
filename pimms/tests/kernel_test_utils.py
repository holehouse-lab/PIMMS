"""
Shared helpers for the kernel correctness / detailed-balance test suite.

These tests exercise the optimized Cython kernels (fast serial crankshaft, the
parallel checkerboard kernel, the slither/reptation megamove and the TSMMC
moves) across the full matrix of:

  * dimensionality      : 2D and 3D
  * interaction range   : SR-only, SR+LR, SR+LR+SLR  (the three forcefield kinds)
  * boundary conditions : HARDWALL True and HARDWALL False (PBC)

Two correctness properties are checked:

  ENERGY CONSISTENCY  - after a megamove the incrementally-tracked energy must
                        equal a from-scratch recomputation of the mutated state.
                        This is the bit-exact correctness gate for every kernel.

  DETAILED BALANCE    - a move-under-test must sample the same Boltzmann
                        equilibrium as the trusted crankshaft reference (the fast
                        serial crankshaft kernel is bit-exact to the reference
                        mega_crank kernel, so it is the trusted yardstick).

The forcefield kinds map directly onto the parameter-file format:

    "A A -8"        -> short range only                       (len-3 line)
    "A A -8 -4"     -> short + long range  (SLR = 0)          (len-4 line)
    "A A -8 -4 2"   -> short + long + super-long range        (len-5 line)

A residue is treated as long-range iff it appears in a len>=4 line, and LR/SLR
interactions only act between two long-range residues.
"""
import os
import copy
import contextlib

import numpy as np

from pimms.keyfile_parser import KeyFileParser
from pimms.simulation import Simulation
from pimms import crankshaft_list_functions as clf
import pimms.mega_crank as ref_kernel_3D
import pimms.mega_crank_2D as ref_kernel_2D
import pimms.mega_crank_fast as fk


# ---------------------------------------------------------------------------
# forcefield + keyfile generation
# ---------------------------------------------------------------------------

# interaction values per pair: (SR, LR, SLR). All distinct and non-zero so that a
# kernel mishandling any single range is caught by the consistency check.
_PAIR_VALUES = {
    ("A", "A"): (-8, -4, 2),
    ("B", "B"): (-6, -3, 3),
    ("A", "B"): (-3, -2, 1),
}


def write_param_file(path, ff_kind):
    """Write a parameter file for ff_kind in {"SR", "LR", "SLR"}.

    SR  -> only the short-range column (no residue is long-range)
    LR  -> short + long range columns (residues become long-range)
    SLR -> short + long + super-long-range columns
    """
    if ff_kind not in ("SR", "LR", "SLR"):
        raise ValueError(ff_kind)
    n_cols = {"SR": 1, "LR": 2, "SLR": 3}[ff_kind]
    lines = [
        "ANGLE_PENALTY\tA\t30\t10\t0",
        "ANGLE_PENALTY\tB\t50\t20\t0",
        # solvation (bead-solvent) energies - required for every residue type
        "A\t0\t-2",
        "B\t0\t-1",
    ]
    for (r1, r2), vals in _PAIR_VALUES.items():
        cols = "\t".join(str(v) for v in vals[:n_cols])
        lines.append(f"{r1}  {r2}\t{cols}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# A standard mixed system: heteropolymers (exercise the slither O(N) path and
# every residue type), homopolymers (the slither O(1) path) and single beads
# (the slither -> translation path). Counts/box scale with dimensionality so the
# system is dense enough that SR, LR and SLR shells are all populated.
_SYSTEMS = {
    2: dict(box=[18, 18], chains=[(7, "AABB"), (7, "AAAA"), (8, "A"), (4, "AABBA")]),
    3: dict(box=[13, 13, 13], chains=[(8, "AABB"), (8, "AAAA"), (10, "A"), (5, "AABBA")]),
}


def write_keyfile(path, dim, hardwall, moves, *, box=None, chains=None, seed=11,
                  n_steps=10, equilibration=1, temperature=55, extra=None):
    """Write a KEYFILE.kf. `moves` is a dict of {MOVE_KEYWORD: fraction}."""
    spec = _SYSTEMS[dim]
    box = box if box is not None else spec["box"]
    chains = chains if chains is not None else spec["chains"]
    lines = [
        "DIMENSIONS : " + " ".join(str(b) for b in box),
        "PARAMETER_FILE : params.prm",
        f"SEED : {seed}",
        f"TEMPERATURE : {temperature}",
        f"HARDWALL : {'True' if hardwall else 'False'}",
        f"N_STEPS : {n_steps}",
        f"EQUILIBRATION : {equilibration}",
        "EXPERIMENTAL_FEATURES : True",
    ]
    for count, seq in chains:
        lines.append(f"CHAIN : {count} {seq}")
    for kw, frac in moves.items():
        lines.append(f"{kw} : {frac}")
    if extra:
        for k, v in extra.items():
            lines.append(f"{k} : {v}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# building / inspecting simulation state
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _chdir(path):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


class State:
    """Bundle of the objects a kernel test needs from a built Simulation."""

    def __init__(self, sim):
        self.sim = sim
        self.lattice = sim.LATTICE
        self.ham = sim.Hamiltonian
        self.acc = sim.ACC
        self.hardwall_int = 1 if sim.hardwall else 0
        decomposition = self.ham.evaluate_total_energy(self.lattice)
        self.energy = int(decomposition[0])
        # (total, local/SR, LR, SLR, angle)
        self.energy_terms = tuple(int(x) for x in decomposition)
        self.idx0 = clf.update_idx_to_bead(self.lattice)
        self.dim = len(self.lattice.dimensions)

    @property
    def tables(self):
        return (self.ham.residue_interaction_table,
                self.ham.LR_residue_interaction_table,
                self.ham.SLR_residue_interaction_table,
                self.ham.angle_lookup)

    def fresh(self):
        """A fresh (grid, type_grid, idx) copy of the initial state."""
        return (self.lattice.grid.copy(), self.lattice.type_grid.copy(), self.idx0.copy())

    def has_LR(self):
        return bool(np.any(np.asarray(self.idx0)[:, 1] == 1))


def build_state(tmpdir, dim, ff_kind, hardwall, moves, **kw):
    """Write a forcefield + keyfile into tmpdir and build the Simulation."""
    tmpdir = str(tmpdir)
    write_param_file(os.path.join(tmpdir, "params.prm"), ff_kind)
    write_keyfile(os.path.join(tmpdir, "KEYFILE.kf"), dim, hardwall, moves, **kw)
    with _chdir(tmpdir):
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(keyfile.keyword_lookup)
    return State(sim)


def run_sim_with_energy_check(tmpdir, dim, ff, hardwall, moves, *, n_steps=60,
                              energy_check=5, seed=11, extra=None, temperature=40,
                              box=None, chains=None, return_sim=False):
    """Run a short in-process simulation with ENERGY_CHECK enabled.

    ENERGY_CHECK recomputes the total energy from scratch every `energy_check`
    steps and raises SimulationEnergyException if it disagrees with the tracked
    energy, so a clean return proves the move(s) kept the energy consistent. IO
    frequencies are pushed high to keep the run fast. Returns the ENERGY.dat trace
    (or `(trace, sim)` when `return_sim` is set, e.g. to inspect move diagnostics).
    """
    base = {
        "ENERGY_CHECK": energy_check,
        "PRINT_FREQ": 1000000,
        "XTC_FREQ": 1000000,
        "ANALYSIS_FREQ": 1000000,
        "RESTART_FREQ": 1000000,
        "EN_FREQ": 10,
        "TSMMC_JUMP_TEMP": 120,
        "TSMMC_STEP_MULTIPLIER": 15,
        "TSMMC_NUMBER_OF_POINTS": 8,
    }
    if extra:
        base.update(extra)
    write_param_file(os.path.join(str(tmpdir), "params.prm"), ff)
    write_keyfile(os.path.join(str(tmpdir), "KEYFILE.kf"), dim, hardwall, moves,
                  temperature=temperature, n_steps=n_steps, equilibration=max(1, n_steps // 6),
                  seed=seed, extra=base, box=box, chains=chains)
    with _chdir(str(tmpdir)):
        keyfile = KeyFileParser("KEYFILE.kf")
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            sim = Simulation(keyfile.keyword_lookup)
            sim.run_simulation()
        trace = np.loadtxt("ENERGY.dat", delimiter="\t")
    if return_sim:
        return trace, sim
    return trace


def chain_meta(idx):
    """(offsets, lengths, homo-flags) for the contiguous chains in idx_to_bead."""
    cids = np.asarray(idx)[:, 4]
    offs, lens, homo = [], [], []
    i, n = 0, len(cids)
    while i < n:
        j = i
        while j < n and cids[j] == cids[i]:
            j += 1
        offs.append(i)
        lens.append(j - i)
        homo.append(1 if len(np.unique(np.asarray(idx)[i:j, 2])) == 1 else 0)
        i = j
    return (np.array(offs, np.int32), np.array(lens, np.int32), np.array(homo, np.int32))


def recompute_energy(state, grid, type_grid, idx):
    """Total energy of the mutated (grid, type_grid, idx) state, from scratch."""
    lat = copy.deepcopy(state.lattice)
    lat.grid = grid
    lat.type_grid = type_grid
    local = 0
    for chainID in sorted(lat.chains.keys()):
        n = len(lat.chains[chainID].get_ordered_positions())
        lat.chains[chainID].set_ordered_positions(idx[local:local + n, 5:].tolist())
        local += n
    return int(state.ham.evaluate_total_energy(lat)[0])


# ---------------------------------------------------------------------------
# kernel drivers - each runs ONE megamove on the supplied (grid, tg, idx),
# mutating them in place, and returns the new incremental energy.
# ---------------------------------------------------------------------------

def make_bead_selector(state, substeps, seed=None):
    """A bead selector for a crankshaft megastep (random unless reused)."""
    if seed is not None:
        np.random.seed(seed)
    return clf.bead_selector_constructor(len(state.idx0), substeps, state.lattice,
                                         frozen_chains=[], safecheck=True)


def crank_megastep(state, grid, tg, idx, energy, seed, *, substeps=4000, fast=True, bsel=None):
    if bsel is None:
        bsel = make_bead_selector(state, substeps)
    kernel = (fk.mega_crank if state.dim == 3 else fk.mega_crank_2D) if fast \
        else (ref_kernel_3D.mega_crank if state.dim == 3 else ref_kernel_2D.mega_crank_2D)
    out = kernel(grid, tg, idx, *state.tables, energy, state.acc.invtemp,
                 substeps, bsel, seed, state.hardwall_int)
    return out[0]


def slither_megastep(state, grid, tg, idx, energy, seed, *, substeps=8):
    offs, lens, homo = chain_meta(idx)
    sel = np.repeat(np.arange(len(offs), dtype=np.int32), substeps)
    np.random.RandomState(seed).shuffle(sel)
    kernel = fk.mega_slither if state.dim == 3 else fk.mega_slither_2D
    e, _ = kernel(grid, tg, idx, offs, lens, homo, sel, *state.tables,
                  energy, state.acc.invtemp, seed, state.hardwall_int, int(lens.max()))
    return e


def pull_megastep(state, grid, tg, idx, energy, seed, *, substeps=10):
    """One pull megamove (each L>=3 chain pulled `substeps` times). Returns (energy, accepted)."""
    offs, lens, homo = chain_meta(idx)
    sel = np.repeat(np.arange(len(offs), dtype=np.int32), substeps)
    np.random.RandomState(seed).shuffle(sel)
    kernel = fk.mega_pull if state.dim == 3 else fk.mega_pull_2D
    return kernel(grid, tg, idx, offs, lens, homo, sel, *state.tables,
                  energy, state.acc.invtemp, seed, state.hardwall_int, int(lens.max()))


def parallel_megastep(state, grid, tg, idx, energy, seed, *, substeps=8000, nthreads=4):
    e, _ = fk.mega_crank_parallel(grid, tg, idx, *state.tables, energy,
                                  state.acc.invtemp, substeps, seed,
                                  state.hardwall_int, nthreads)
    return e


def db_compare(state, test_step, *, equilibrate, sample, crank_substeps=2500,
               equil_seed=1000, sample_seed=5000):
    """Detailed-balance comparison driven directly by the kernels.

    The system is first equilibrated with the trusted crankshaft kernel, then -
    starting from that SAME configuration - both the crankshaft reference and the
    move-under-test are run for `sample` megamoves and their energy traces
    collected. A move that respects detailed balance holds the same equilibrium
    as crankshaft; a move that violates it drifts away.

    Everything is seeded, so the result is deterministic (no statistical flake):
    returns (reference_trace, test_trace) as numpy arrays.
    """
    g, t, i = state.fresh()
    e = state.energy
    for m in range(equilibrate):
        e = crank_megastep(state, g, t, i, e, equil_seed + m, substeps=crank_substeps)
    base = (np.asarray(g).copy(), np.asarray(t).copy(), np.asarray(i).copy(), e)

    def trace(step):
        gg, tt, ii, ee = base[0].copy(), base[1].copy(), base[2].copy(), base[3]
        out = []
        for m in range(sample):
            ee = step(state, gg, tt, ii, ee, sample_seed + m)
            out.append(ee)
        return np.array(out, dtype=float)

    ref = trace(lambda s, g, t, i, e, sd: crank_megastep(s, g, t, i, e, sd,
                                                         substeps=crank_substeps))
    test = trace(test_step)
    return ref, test


def assert_same_equilibrium(ref, test, label, rel_floor=0.03, k_sigma=2.5):
    """Assert two equilibrium energy traces agree within a statistical tolerance.

    Tolerance = k_sigma * max(per-sample std) + rel_floor * |mean|. The known
    detailed-balance bugs shifted the mean by many sigma / a large fraction of the
    energy, so they fail this comfortably while a correct move passes.
    """
    mr, mt = ref.mean(), test.mean()
    tol = k_sigma * max(ref.std(), test.std()) + rel_floor * abs(mr)
    assert abs(mt - mr) <= tol, (
        f"{label}: detailed balance violated - reference E={mr:.1f}+/-{ref.std():.1f}, "
        f"test E={mt:.1f}+/-{test.std():.1f}, |diff|={abs(mt - mr):.1f} > tol={tol:.1f}")
