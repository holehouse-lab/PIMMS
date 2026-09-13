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


def frozen_bead_mask(idx, frozen=None):
    """Per-bead int32 frozen mask (1 = bead's chainID is in `frozen`).

    `idx` column 4 is the chainID. With `frozen=None` (or empty) returns all
    zeros, which leaves the parallel kernels bit-identical to the no-freeze case.
    """
    nb = idx.shape[0]
    if not frozen:
        return np.zeros(nb, dtype=np.int32)
    return np.ascontiguousarray(np.isin(np.asarray(idx)[:, 4], list(frozen)).astype(np.int32))


def parallel_megastep(state, grid, tg, idx, energy, seed, *, substeps=8000, nthreads=4, frozen=None):
    e, _ = fk.mega_crank_parallel(grid, tg, idx, *state.tables, energy,
                                  state.acc.invtemp, substeps, seed,
                                  state.hardwall_int, nthreads, frozen_bead_mask(idx, frozen))
    return e


def parallel_megastep_2D(state, grid, tg, idx, energy, seed, *, substeps=8000, nthreads=4, frozen=None):
    """Drive one 2D parallel checkerboard megamove (mega_crank_parallel_2D)."""
    e, _ = fk.mega_crank_parallel_2D(grid, tg, idx, *state.tables, energy,
                                     state.acc.invtemp, substeps, seed,
                                     state.hardwall_int, nthreads, frozen_bead_mask(idx, frozen))
    return e


def slither_parallel_megastep(state, grid, tg, idx, energy, seed, *, substeps=8, nthreads=4, frozen=None):
    """Drive one parallel slither megamove (mega_slither_parallel / _2D)."""
    offs, lens, homo = chain_meta(idx)
    sel = np.repeat(np.arange(len(offs), dtype=np.int32), substeps)
    np.random.RandomState(seed).shuffle(sel)
    kernel = fk.mega_slither_parallel if state.dim == 3 else fk.mega_slither_parallel_2D
    e, _ = kernel(grid, tg, idx, offs, lens, homo, sel, *state.tables,
                  energy, state.acc.invtemp, seed, state.hardwall_int, int(lens.max()),
                  nthreads, frozen_bead_mask(idx, frozen))
    return e


def pull_parallel_megastep(state, grid, tg, idx, energy, seed, *, substeps=10, nthreads=4, frozen=None):
    """Drive one parallel pull megamove (mega_pull_parallel / _2D)."""
    offs, lens, homo = chain_meta(idx)
    sel = np.repeat(np.arange(len(offs), dtype=np.int32), substeps)
    np.random.RandomState(seed).shuffle(sel)
    kernel = fk.mega_pull_parallel if state.dim == 3 else fk.mega_pull_parallel_2D
    e, _ = kernel(grid, tg, idx, offs, lens, homo, sel, *state.tables,
                  energy, state.acc.invtemp, seed, state.hardwall_int, int(lens.max()),
                  nthreads, frozen_bead_mask(idx, frozen))
    return e


@contextlib.contextmanager
def scaled_invtemp(state, scale):
    """Temporarily scale the inverse temperature the kernels are handed.

    This is the fault injection used by the detailed-balance positive controls:
    running ONE move at ``invtemp * scale`` while the crankshaft reference stays
    at the true beta is exactly what a wrong Metropolis exponent, a wrong
    Hastings ratio or a units regression inside that kernel would look like. A
    fixture that cannot resolve such a scaling has no power against that class
    of bug, whatever its assertion says.

    Parameters
    ----------
    state : State
        The bundle returned by :func:`build_state`; ``state.acc.invtemp`` is the
        value every kernel driver in this module passes down to Cython.

    scale : float
        Multiplier applied to ``state.acc.invtemp`` for the duration of the
        block. ``1.0`` is a no-op.

    Yields
    ------
    State
        The same state object, with the scaled inverse temperature in place.
    """
    original = state.acc.invtemp
    state.acc.invtemp = original * scale
    try:
        yield state
    finally:
        state.acc.invtemp = original


def bond_walked_positions(idx, dims):
    """Whole-chain positions, unwrapped by walking the bonds.

    Wrapped lattice positions cannot be used for an intra-chain observable: a
    chain that straddles a periodic face is torn in two, and a minimum-image
    correction against the chain centre tears any chain longer than L/2 the same
    way. Because every bond in this model is a Chebyshev-1 step, walking the
    bonds is exact - the minimum image of each consecutive difference IS the
    bond vector - so this reconstructs each chain in a single image regardless of
    its extent.

    Parameters
    ----------
    idx : numpy.ndarray
        The idx_to_bead table; column 4 is the chainID and columns 5 onwards are
        the wrapped positions. Beads of a chain are contiguous.

    dims : sequence of int
        Box dimensions, one per axis.

    Returns
    -------
    numpy.ndarray
        Float array of shape (n_beads, dim) holding the unwrapped positions,
        each chain anchored on its own first bead's in-box position.
    """
    idx = np.asarray(idx)
    dim = len(dims)
    L = np.asarray(dims, dtype=np.int64)
    pos = idx[:, 5:5 + dim].astype(np.int64)
    cids = idx[:, 4]
    steps = np.zeros_like(pos)
    diff = pos[1:] - pos[:-1]
    steps[1:] = ((diff + L // 2) % L) - L // 2
    # a chain boundary is not a bond - restart the walk there
    starts = np.flatnonzero(np.r_[True, cids[1:] != cids[:-1]])
    steps[starts] = 0
    walked = np.cumsum(steps, axis=0)
    # re-anchor each chain on its own first bead
    anchor_row = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, len(pos)]))
    walked = walked - walked[starts][anchor_row] + pos[starts][anchor_row]
    return walked.astype(float)


def mean_chain_rg2(state, idx):
    """Mean squared radius of gyration over the chains of a state.

    A purely conformational observable: it is blind to the interaction tables
    and so it sees biases that leave the mean energy untouched (a proposal
    asymmetry in a rigid move, for instance). Single-bead chains contribute
    zero, which is their exact Rg^2, so they are included rather than dropped.

    Parameters
    ----------
    state : State
        The bundle returned by :func:`build_state` (used only for the box).

    idx : numpy.ndarray
        The idx_to_bead table of the configuration to measure.

    Returns
    -------
    float
        The mean over chains of Rg^2, in lattice units squared.
    """
    dims = list(state.lattice.dimensions)
    pos = bond_walked_positions(idx, dims)
    offs, lens, _ = chain_meta(idx)
    offs = offs.astype(np.int64)
    com = np.add.reduceat(pos, offs, axis=0) / lens[:, None]
    dev = pos - np.repeat(com, lens, axis=0)
    sq = np.add.reduceat((dev * dev).sum(axis=1), offs) / lens
    return float(sq.mean())


def _equilibrated_base(state, *, equilibrate, crank_substeps, equil_seed):
    """Crank-equilibrate a fresh state and return the (grid, tg, idx, energy) it reached."""
    g, t, i = state.fresh()
    e = state.energy
    for m in range(equilibrate):
        e = crank_megastep(state, g, t, i, e, equil_seed + m, substeps=crank_substeps)
    return (np.asarray(g).copy(), np.asarray(t).copy(), np.asarray(i).copy(), e)


def _run_trace(state, base, step, sample, sample_seed, *, with_rg=False):
    """Run `sample` megamoves of `step` from `base`, returning the energy (and Rg^2) trace."""
    gg, tt, ii, ee = base[0].copy(), base[1].copy(), base[2].copy(), base[3]
    energies = np.empty(sample, dtype=float)
    rg2 = np.empty(sample, dtype=float) if with_rg else None
    for m in range(sample):
        ee = step(state, gg, tt, ii, ee, sample_seed + m)
        energies[m] = ee
        if with_rg:
            rg2[m] = mean_chain_rg2(state, ii)
    return energies if not with_rg else (energies, rg2)


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
    base = _equilibrated_base(state, equilibrate=equilibrate,
                              crank_substeps=crank_substeps, equil_seed=equil_seed)
    ref = _run_trace(state, base,
                     lambda s, g, t, i, e, sd: crank_megastep(s, g, t, i, e, sd,
                                                              substeps=crank_substeps),
                     sample, sample_seed)
    test = _run_trace(state, base, test_step, sample, sample_seed)
    return ref, test


class DBResult:
    """The three traces of a detailed-balance comparison with a power control.

    Attributes
    ----------
    ref_energy, test_energy, control_energy : numpy.ndarray
        Energy traces of the crankshaft reference, the move under test, and the
        deliberately mis-tempered copy of the move under test.

    ref_rg2, test_rg2, control_rg2 : numpy.ndarray
        The matching mean-Rg^2 traces (a conformational observable that does not
        go through the energy at all).

    control_scale : float
        The inverse-temperature scaling that was injected into the control run.
    """

    __slots__ = ("ref_energy", "test_energy", "control_energy",
                 "ref_rg2", "test_rg2", "control_rg2", "control_scale")

    def __init__(self, ref, test, control, control_scale):
        self.ref_energy, self.ref_rg2 = ref
        self.test_energy, self.test_rg2 = test
        self.control_energy, self.control_rg2 = control
        self.control_scale = control_scale


def db_compare_with_control(state, test_step, control_step, *, equilibrate, sample,
                            control_scale, crank_substeps=2500, equil_seed=1000,
                            sample_seed=5000):
    """Detailed-balance comparison plus the positive control that measures its power.

    Three traces are grown from ONE crank-equilibrated configuration: the
    crankshaft reference, the move under test, and `control_step` - the same move
    with its inverse temperature scaled by `control_scale`, i.e. a deliberate
    Metropolis error of known size. The caller asserts that the first two agree
    AND that the third does not. Without that second half a fixture can decay
    into a vacuous green (the state this suite was in: no kernel-level fixture
    resolved a 20 % acceptance error in the move it was testing).

    Both an energy trace and a mean-Rg^2 trace are recorded, because the mean
    energy is blind by construction to any bias that leaves it unchanged - a
    rigid move with a one-way proposal, for instance.

    Parameters
    ----------
    state : State
        The system, as returned by :func:`build_state`.

    test_step : callable
        ``f(state, grid, type_grid, idx, energy, seed) -> energy``, one megamove
        of the move under test.

    control_step : callable
        The same signature, running the move under test inside
        ``scaled_invtemp(state, control_scale)``.

    equilibrate : int
        Crankshaft megamoves used to equilibrate before any trace is taken.

    sample : int
        Megamoves per trace.

    control_scale : float
        The inverse-temperature scaling injected into `control_step`, recorded on
        the result so the assertion message can name it.

    crank_substeps : int, optional
        Bead moves per crankshaft megamove, for the equilibration and reference.

    equil_seed, sample_seed : int, optional
        Base seeds; every megamove uses a distinct derived seed, so the whole
        comparison is deterministic.

    Returns
    -------
    DBResult
        The six traces and the control scale.
    """
    base = _equilibrated_base(state, equilibrate=equilibrate,
                              crank_substeps=crank_substeps, equil_seed=equil_seed)
    ref = _run_trace(state, base,
                     lambda s, g, t, i, e, sd: crank_megastep(s, g, t, i, e, sd,
                                                              substeps=crank_substeps),
                     sample, sample_seed, with_rg=True)
    test = _run_trace(state, base, test_step, sample, sample_seed, with_rg=True)
    control = _run_trace(state, base, control_step, sample, sample_seed, with_rg=True)
    return DBResult(ref, test, control, control_scale)


def integrated_autocorr_time(x):
    """Integrated autocorrelation time of a Monte Carlo trace, in samples.

    Sums the empirical autocorrelation function from lag 1 until it first drops
    below zero (an initial-positive-sequence style truncation) and adds the 1/2
    that the lag-0 term contributes, so an uncorrelated trace gives 0.5.

    Parameters
    ----------
    x : array_like
        The trace.

    Returns
    -------
    float
        tau, in units of samples. A trace whose tau approaches its own length
        carries almost no independent information, which is the leading
        indicator that a fixture has become underpowered.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 4:
        return 0.5
    d = x - x.mean()
    var = float(np.dot(d, d) / n)
    if var == 0.0:
        return 0.5
    tau = 0.5
    for lag in range(1, min(n // 4, 200)):
        rho = float(np.dot(d[:-lag], d[lag:]) / ((n - lag) * var))
        if rho <= 0.0:
            break
        tau += rho
    return tau


def _sem_with_autocorr(x):
    """Standard error of the mean of a (possibly correlated) MC energy trace.

    Estimates the integrated autocorrelation time tau with
    :func:`integrated_autocorr_time`, then inflates the naive SEM by
    sqrt(2 * tau). For the megamove traces used here tau is O(1), so this stays
    close to the naive SEM while remaining honest for slower-mixing moves.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 4:
        return float(x.std()) if n else 0.0
    var = float(np.dot(x - x.mean(), x - x.mean()) / n)
    if var == 0.0:
        return 0.0
    return float(np.sqrt(2.0 * integrated_autocorr_time(x) * var / n))


def assert_trace_is_well_sampled(trace, label, *, min_samples_per_tau=20.0):
    """Assert a trace holds enough independent information to test anything with.

    A trace whose integrated autocorrelation time is a sizeable fraction of its
    own length has only a handful of effective samples, and an equivalence test
    built on it certifies nothing however tight its arithmetic looks. Asserting
    tau explicitly stops a fixture from silently becoming underpowered when a
    box, a temperature or a substep count is changed later.

    Parameters
    ----------
    trace : array_like
        The trace to check (normally the trace of the move under test).

    label : str
        Fixture name for the failure message.

    min_samples_per_tau : float, optional
        Required ratio of trace length to tau. The default of 20 corresponds to
        at least ~10 statistically independent samples (tau counts the lag-0 1/2).

    Returns
    -------
    float
        The measured tau, so a caller can report it.
    """
    trace = np.asarray(trace, dtype=float)
    tau = integrated_autocorr_time(trace)
    assert tau * min_samples_per_tau <= len(trace), (
        f"{label}: trace is underpowered - integrated autocorrelation time tau={tau:.1f} "
        f"over only {len(trace)} samples (need tau <= {len(trace) / min_samples_per_tau:.1f}); "
        f"raise the sample count or the substeps per megamove")
    return tau


def assert_same_equilibrium(ref, test, label, rel_floor=0.005, k_sigma=4.0, ref_only=False):
    """Assert two equilibrium energy traces agree within a statistical tolerance.

    Tolerance = k_sigma * SEM_diff + rel_floor * |mean|, where SEM_diff combines
    the autocorrelation-aware standard errors OF THE MEANS of the two traces.

    With ``ref_only=True`` the tolerance is built from the REFERENCE trace's SEM
    alone and the test trace must additionally be stationary (its two halves
    must agree within the same tolerance). A broken move whose run never
    equilibrates has a drifting trace with a huge SEM; letting that SEM into the
    tolerance let a beta x1.5 jump-and-relax bug pass at the Simulation level
    even under the autocorrelation-aware criterion.

    This replaces the old tolerance of 2.5 * max(per-sample std) + 3% * |E|,
    which was ~30x the SEM on the actual fixtures: an injected acceptance bug
    with beta scaled by 1.5 (a gross Metropolis error) passed EVERY
    detailed-balance case under it. Two specific defects are addressed: the
    per-sample std is not the uncertainty of a mean (the traces have hundreds of
    samples), and using the TEST trace's own spread let a broken move inflate
    its own tolerance. Only the reference and test SEMs enter now, and the
    relative floor is 0.5%.

    ANY comparison of a move under test against a trusted reference must pass
    ``ref_only=True``. The default two-SEM form is the SEM of a difference of two
    means and is only meaningful when NEITHER trace is on trial - two relabelled
    copies of the same physical system, say. Used as a pass criterion for a move
    under test it converts absence of evidence into evidence of absence: a broken
    move drifts, its autocorrelation time balloons, its SEM balloons, and the
    tolerance balloons with it. Measured on the kernel fixtures, the two-SEM form
    missed a 20 % acceptance error in the parallel slither in 3 seeds of 3 even
    with seven times the samples, where the reference-only form caught it 3 of 3.
    """
    ref = np.asarray(ref, dtype=float)
    test = np.asarray(test, dtype=float)
    mr, mt = ref.mean(), test.mean()
    if ref_only:
        sem = float(_sem_with_autocorr(ref))
    else:
        sem = float(np.hypot(_sem_with_autocorr(ref), _sem_with_autocorr(test)))
    tol = k_sigma * sem + rel_floor * abs(mr)
    if ref_only:
        half = len(test) // 2
        drift = abs(test[:half].mean() - test[half:].mean())
        # The difference of two HALF-trace means has twice the standard error of
        # the full-trace mean (each half has half the samples, and the two errors
        # add in quadrature), so the drift must be judged against 2 * sem, not
        # sem. Judging it against the full-trace tolerance made the guard a ~2
        # sigma test that fired on correct-but-slowly-mixing moves at the few-%
        # level - "works on the seed we picked" rather than a real criterion.
        drift_tol = k_sigma * 2.0 * sem + rel_floor * abs(mr)
        assert drift <= drift_tol, (
            f"{label}: test trace is not stationary - first-half E={test[:half].mean():.1f}, "
            f"second-half E={test[half:].mean():.1f}, |drift|={drift:.1f} > tol={drift_tol:.1f}")
    assert abs(mt - mr) <= tol, (
        f"{label}: detailed balance violated - reference E={mr:.1f} (SEM {_sem_with_autocorr(ref):.2f}), "
        f"test E={mt:.1f} (SEM {_sem_with_autocorr(test):.2f}), |diff|={abs(mt - mr):.1f} > tol={tol:.1f}")


def equilibrium_margin(ref, test, rel_floor=0.005, k_sigma=4.0):
    """|mean difference| in units of the reference-only tolerance.

    The quantity :func:`assert_same_equilibrium` thresholds at 1.0 under
    ``ref_only=True``. Useful for reporting how much headroom (or how little) a
    fixture has, and for the positive controls, which need the margin of the
    deliberately biased run.

    Parameters
    ----------
    ref, test : array_like
        Reference and test traces of the same observable.

    rel_floor, k_sigma : float, optional
        The tolerance parameters, as in :func:`assert_same_equilibrium`.

    Returns
    -------
    float
        |mean(test) - mean(ref)| / tolerance.
    """
    ref = np.asarray(ref, dtype=float)
    test = np.asarray(test, dtype=float)
    tol = k_sigma * float(_sem_with_autocorr(ref)) + rel_floor * abs(ref.mean())
    return float(abs(test.mean() - ref.mean()) / tol) if tol > 0 else float("inf")


def assert_control_is_resolved(ref, control, label, scale, rel_floor=0.005, k_sigma=4.0):
    """Assert a fixture actually RESOLVES a known injected error (its positive control).

    The companion of :func:`assert_same_equilibrium`. The move under test is
    re-run with its inverse temperature scaled by `scale`, which is a Metropolis
    error of known size, and the comparison against the reference must FAIL. A
    detailed-balance fixture that passes its main assertion but cannot fail this
    one is not evidence of anything: it has no power, and every fixture in this
    file was in that state for a 20 % error before these controls were added.

    Parameters
    ----------
    ref : array_like
        The crankshaft reference trace.

    control : array_like
        The trace of the move under test run at the scaled inverse temperature.

    label : str
        Fixture name for the failure message.

    scale : float
        The injected inverse-temperature scaling, for the message.

    rel_floor, k_sigma : float, optional
        The same tolerance parameters used for the main assertion, so the control
        measures the power of the assertion that is actually made.
    """
    margin = equilibrium_margin(ref, control, rel_floor=rel_floor, k_sigma=k_sigma)
    assert margin > 1.0, (
        f"{label}: POSITIVE CONTROL FAILED - the move under test run at invtemp x {scale} "
        f"(a {abs(scale - 1) * 100:.0f} % Metropolis error) is still within tolerance "
        f"(margin {margin:.2f} <= 1). This fixture cannot resolve an error of that size, so "
        f"its detailed-balance assertion certifies nothing; lengthen the traces, raise the "
        f"weight of the move under test, or lower rel_floor until the control bites")
    return margin
