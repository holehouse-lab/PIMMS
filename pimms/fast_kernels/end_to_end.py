## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""End-to-end full-simulation comparison: reference kernel vs fast kernel.

Runs the complete two_phase_equilibrium_demo simulation twice in isolated scratch
directories - once forcing the stock pimms.mega_crank crankshaft kernel, once
with the production pimms.mega_crank_fast - by swapping the ``mega_crank`` entry
point that pimms.moves.system_shake dispatches to - then compares the final
ENERGY.dat and reports wall-clock speedup. The demo keyfile is rewritten into a
serial, crankshaft-only run for the comparison (the reference module implements
nothing else).

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


def _pin_keywords(keyfile_text, pins):
    """
    Return ``keyfile_text`` with every keyword in ``pins`` set to the given value.

    A keyword already present is replaced on its own line (anchored at the line
    start, so a mention inside a comment cannot suppress the pin); a keyword
    that is absent is prepended.

    Parameters
    ----------
    keyfile_text : str
        Contents of a PIMMS keyfile.

    pins : dict of str to str
        ``keyword -> value`` pairs to enforce.

    Returns
    -------
    str
        The rewritten keyfile text.
    """
    import re
    out = keyfile_text
    for key, value in pins.items():
        pattern = re.compile(r'^\s*' + re.escape(key) + r'\s*:.*$', re.M)
        if pattern.search(out):
            out = pattern.sub(f"{key} : {value}", out, count=1)
        else:
            out = f"{key} : {value}\n" + out
    return out


def run_once(tag, use_fast):
    """Run the full demo simulation once in an isolated scratch directory.

    Creates a clean ``e2e_<tag>`` directory under ``SCRATCH``, copies in the
    demo's parameter file and a copy of the keyfile pinned to a serial,
    crankshaft-only run with a fixed ``SEED`` and ``EN_FREQ`` so the run is
    deterministic and writes energies, then swaps the crankshaft entry point of
    ``pimms.moves.mega_crank_fast`` to the chosen kernel before running the
    complete simulation with stdout suppressed. After the run it reads back the
    final line of ``ENERGY.dat``.

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
    # The reference module implements only the serial 3D crankshaft kernel, so
    # the demo keyfile (which ships as a PARALLELIZE + MOVE_SLITHER run) is
    # rewritten line by line into a serial crankshaft-only run with a fixed
    # SEED and an energy cadence that actually produces ENERGY.dat rows.
    kf = _pin_keywords(kf, {"SEED": "424242", "EN_FREQ": "10", "PARALLELIZE": "False",
                            "MOVE_CRANKSHAFT": "1.0", "MOVE_SLITHER": "0", "MOVE_PULL": "0",
                            "MOVE_VMMC": "0", "MOVE_MULTICHAIN_TSMMC": "0",
                            "MOVE_CTSMMC": "0", "MOVE_SYSTEM_TSMMC": "0"})
    with open(os.path.join(rundir, "KEYFILE.kf"), "w") as fh:
        fh.write(kf)

    # system_shake dispatches to moves.mega_crank_fast.mega_crank, but the same
    # module object also supplies the layout helpers that Simulation.__init__ and
    # system_shake call unconditionally (parallel_layout_info,
    # parallel_crank_layout_info) and the other megamoves. Replacing the whole
    # module with the reference kernel used to crash inside Simulation.__init__,
    # so only the crankshaft entry point is swapped.
    import types
    import pimms.moves as moves
    import pimms.mega_crank as ref_kernel
    import pimms.mega_crank_fast as fast_kernel
    from pimms.keyfile_parser import KeyFileParser
    from pimms.simulation import Simulation

    if use_fast:
        moves.mega_crank_fast = fast_kernel
    else:
        shim = types.SimpleNamespace(**{k: getattr(fast_kernel, k) for k in dir(fast_kernel)
                                        if not k.startswith('__')})
        shim.mega_crank = ref_kernel.mega_crank
        moves.mega_crank_fast = shim

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
    if e_ref is None or e_fast is None:
        # None == None must never count as agreement
        print("  ENERGY.dat was empty for at least one run - nothing was compared")
        return 1
    match = (e_ref == e_fast)
    print(f"  final energy identical : {'YES' if match else 'NO'}")
    if t_fast > 0:
        print(f"  whole-simulation speedup: {t_ref / t_fast:.2f}x")
    return 0 if match else 1


if __name__ == "__main__":
    sys.exit(main())
