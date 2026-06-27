"""End-to-end full-simulation comparison: reference kernel vs fast kernel.

Runs the complete two_phase_equilibrium_demo simulation twice in isolated scratch
directories - once forcing the stock pimms.mega_crank kernel, once with the
production pimms.mega_crank_fast - by monkeypatching what pimms.moves.system_shake
dispatches to (moves.mega_crank_fast) - then compares the final ENERGY.dat and
reports wall-clock speedup.

Because the fast kernel preserves the exact RNG stream, the two runs should
produce identical trajectories/energies.

Run from the repository root:

    python pimms/fast_kernels/end_to_end.py
"""
import os
import sys
import time
import shutil
import tempfile
import contextlib
import io

HERE = os.path.dirname(os.path.abspath(__file__))      # pimms/fast_kernels
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))     # repo root
sys.path.insert(0, REPO_ROOT)

DEMO_DIR = os.path.join(REPO_ROOT, "demo_keyfiles", "two_phase_equilibrium_demo")
SCRATCH = os.environ.get("PIMMS_SCRATCH", os.path.join(tempfile.gettempdir(), "pimms_e2e"))


def run_once(tag, use_fast):
    """Run the full demo simulation once in an isolated scratch directory.

    Creates a clean ``e2e_<tag>`` directory under ``SCRATCH``, copies in the
    demo's parameter file and a copy of the keyfile with a fixed ``SEED``
    prepended (if the demo keyfile has none) so the run is deterministic, then
    monkeypatches ``pimms.moves.mega_crank_fast`` to the chosen kernel module
    before running the complete simulation with stdout suppressed. After the run
    it reads back the final line of ``ENERGY.dat``.

    Parameters
    ----------
    tag : str
        Short label used to name the run's scratch subdirectory (e.g.
        ``"ref"`` or ``"fast"``).
    use_fast : bool
        If ``True`` dispatch ``system_shake`` to the production fast kernel
        (``pimms.mega_crank_fast``); if ``False`` dispatch to the reference
        kernel (``pimms.mega_crank``).

    Returns
    -------
    elapsed : float
        Wall-clock seconds spent inside ``sim.run_simulation()``.
    final : str or None
        The last non-empty line of ``ENERGY.dat``, or ``None`` if the file is
        absent or empty.
    """
    rundir = os.path.join(SCRATCH, f"e2e_{tag}")
    if os.path.isdir(rundir):
        shutil.rmtree(rundir)
    os.makedirs(rundir)
    shutil.copy(os.path.join(DEMO_DIR, "params.prm"), rundir)
    # Copy the keyfile but pin a fixed SEED so both runs are deterministic
    # (the stock demo keyfile has none -> each run would diverge randomly).
    with open(os.path.join(DEMO_DIR, "KEYFILE.kf")) as fh:
        kf = fh.read()
    if "SEED" not in kf:
        kf = "SEED : 424242\n" + kf
    with open(os.path.join(rundir, "KEYFILE.kf"), "w") as fh:
        fh.write(kf)

    # system_shake dispatches to moves.mega_crank_fast.* , so swap THAT to force
    # either the production fast kernel or the original reference kernel. (Valid
    # for this 3D, non-parallel, crankshaft-only demo, where only .mega_crank is
    # reached; the reference module has no .mega_crank_2D/.mega_crank_parallel.)
    import pimms.moves as moves
    import pimms.mega_crank as ref_kernel
    import pimms.mega_crank_fast as fast_kernel
    from pimms.keyfile_parser import KeyFileParser
    from pimms.simulation import Simulation

    if use_fast:
        moves.mega_crank_fast = fast_kernel
    else:
        moves.mega_crank_fast = ref_kernel

    cwd = os.getcwd()
    os.chdir(rundir)
    devnull = io.StringIO()
    try:
        keyfile = KeyFileParser("KEYFILE.kf")
        sim = Simulation(keyfile.keyword_lookup)
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(devnull):
            sim.run_simulation()
        elapsed = time.perf_counter() - t0
    finally:
        os.chdir(cwd)

    energy_path = os.path.join(rundir, "ENERGY.dat")
    final = None
    if os.path.exists(energy_path):
        with open(energy_path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        final = lines[-1] if lines else None
    return elapsed, final


def main():
    """Run the demo end to end under both kernels and compare results.

    Executes the full ``two_phase_equilibrium_demo`` simulation twice (reference
    kernel then fast kernel) via :func:`run_once`, prints each run's wall-clock
    time and final ``ENERGY.dat`` line, and reports whether the final energies
    are identical (they should be, since the fast kernel preserves the RNG
    stream) along with the whole-simulation speedup.

    Returns
    -------
    int
        ``0`` if the two runs produced identical final energy lines, otherwise
        ``1`` (suitable as a process exit code).
    """
    print("Full-simulation end-to-end comparison (two_phase_equilibrium_demo)\n")

    t_ref, e_ref = run_once("ref", use_fast=False)
    print(f"  reference kernel : {t_ref:7.2f} s   final ENERGY.dat: {e_ref}")

    t_fast, e_fast = run_once("fast", use_fast=True)
    print(f"  fast kernel      : {t_fast:7.2f} s   final ENERGY.dat: {e_fast}")

    print()
    match = (e_ref == e_fast)
    print(f"  final energy identical : {'YES' if match else 'NO'}")
    if t_fast > 0:
        print(f"  whole-simulation speedup: {t_ref / t_fast:.2f}x")
    return 0 if match else 1


if __name__ == "__main__":
    sys.exit(main())
