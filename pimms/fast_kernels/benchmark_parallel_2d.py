## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""Speed harness for the 2D parallel megamove kernels.

This reproduces the tables on the parallelization page of the documentation:
2D, short-range interactions, square boxes uniformly filled to ~7.5% with short
chains (``AABB`` for the crankshaft and slither, ``AABBAB`` for the pull), timed
through the same ``MoveObject.system_*`` dispatch a simulation uses, so the
numbers are what ``PARALLELIZE`` actually buys per megamove. Each move is timed at
two megamove sizes (the defaults and a heavy megamove), and the time spent inside
the compiled kernel is captured separately from the Python bookkeeping around it.
After every configuration the incrementally tracked energy is checked against a
from-scratch recompute.

Run from the repository root::

    python pimms/fast_kernels/benchmark_parallel_2d.py
    python pimms/fast_kernels/benchmark_parallel_2d.py --boxes 64 160 --repeats 3
"""
from __future__ import annotations

import argparse
import contextlib
import os
import statistics
import sys
import tempfile
import time
from typing import Callable

HERE = os.path.dirname(os.path.abspath(__file__))      # pimms/fast_kernels
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))     # repo root
sys.path.insert(0, REPO_ROOT)

import pimms.mega_crank_fast as fk                      # noqa: E402
from pimms import moves                                 # noqa: E402
from pimms.tests import kernel_test_utils as U          # noqa: E402

BOXES: tuple[int, ...] = (64, 96, 160, 256, 400)
THREADS: tuple[int, ...] = (1, 2, 4, 8)
FILL: float = 0.075
TEMPERATURE: float = 500.0
SIZES: dict[str, tuple[int, int]] = {
    "crankshaft": (50_000, 500_000),   # attempts per megamove
    "slither": (10, 100),              # substeps per chain
    "pull": (10, 100),
}
SEQUENCES: dict[str, str] = {"crankshaft": "AABB", "slither": "AABB", "pull": "AABBAB"}
KERNELS: tuple[str, ...] = ("mega_crank_2D", "mega_crank_parallel_2D",
                            "mega_slither_2D", "mega_slither_parallel_2D",
                            "mega_pull_2D", "mega_pull_parallel_2D")

_kernel_seconds: dict[str, float] = {"t": 0.0}


def _wrap_kernels() -> None:
    """Wrap the 2D kernels so the time spent inside them is accumulated.

    The dispatch in :mod:`pimms.moves` looks the kernels up on the module at
    call time, so replacing the module attributes is enough to time them
    without touching the dispatch itself.

    Returns
    -------
    None
    """
    for name in KERNELS:
        real = getattr(fk, name)

        def timed(*args, _real=real, **kwargs):
            """Call the wrapped kernel and add its wall time to the running total.

            Parameters
            ----------
            args : tuple
                The kernel's positional arguments.

            _real : callable, optional
                The kernel being timed (bound at definition time).

            kwargs : dict
                The kernel's keyword arguments.

            Returns
            -------
            tuple
                Whatever the kernel returns.
            """
            t0 = time.perf_counter()
            result = _real(*args, **kwargs)
            _kernel_seconds["t"] += time.perf_counter() - t0
            return result

        setattr(moves.mega_crank_fast, name, timed)


def build_state(box: int, sequence: str) -> tuple[U.State, int]:
    """Build a dispersed 2D short-range system in a temporary directory.

    Parameters
    ----------
    box : int
        Side length of the square box in lattice sites.
    sequence : str
        The chain sequence; the chain count is chosen so the box is ~7.5% full.

    Returns
    -------
    tuple of (State, int)
        The built state (see :class:`pimms.tests.kernel_test_utils.State`) and
        the number of chains it holds.
    """
    n_chains = int(round(FILL * box * box / len(sequence)))
    workdir = tempfile.mkdtemp(prefix="pimms_bench2d_")
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state = U.build_state(workdir, 2, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                              box=[box, box], chains=[(n_chains, sequence)], seed=3,
                              temperature=TEMPERATURE)
    return state, n_chains


def timed_runs(fn: Callable[[], None], repeats: int) -> tuple[float, float]:
    """Time ``fn`` over ``repeats`` calls after one warm-up call.

    Parameters
    ----------
    fn : callable
        A zero-argument callable that performs one megamove.
    repeats : int
        Number of timed calls.

    Returns
    -------
    tuple of (float, float)
        Median wall-clock seconds per call, and median seconds spent inside the
        compiled kernel per call.
    """
    fn()
    walls: list[float] = []
    kernels: list[float] = []
    for _ in range(repeats):
        _kernel_seconds["t"] = 0.0
        t0 = time.perf_counter()
        fn()
        walls.append(time.perf_counter() - t0)
        kernels.append(_kernel_seconds["t"])
    return statistics.median(walls), statistics.median(kernels)


def bench_move(name: str, state: U.State, size: int, repeats: int) -> dict[str, float | bool]:
    """Time one move at one megamove size, serially and at every thread count.

    Parameters
    ----------
    name : str
        ``'crankshaft'``, ``'slither'`` or ``'pull'``.
    state : State
        The system to run the megamoves on; it is mutated in place.
    size : int
        Attempts per megamove (crankshaft) or substeps per chain (slither, pull).
    repeats : int
        Timed megamoves per configuration.

    Returns
    -------
    dict
        ``serial_wall``, ``serial_kernel`` and, for every thread count ``n`` in
        :data:`THREADS`, ``par<n>_wall`` and ``par<n>_kernel`` (seconds), plus
        ``energy_ok``: whether the tracked energy still matches a from-scratch
        recompute after every configuration ran.
    """
    mover = state.sim.MOVER
    current = {"lattice": state.lattice,
               "energy": int(state.ham.evaluate_total_energy(state.lattice)[0])}

    def run(parallel: bool, threads: int) -> None:
        """Run one megamove of this move on the shared state.

        Parameters
        ----------
        parallel : bool
            Whether to request the parallel kernels.

        threads : int
            OpenMP threads for the parallel kernels.

        Returns
        -------
        None
            The shared lattice and energy are updated in place.
        """
        if name == "crankshaft":
            result = mover.system_shake(current["lattice"], current["energy"], state.acc,
                                        state.ham, size, "UNIFORM", parallelize=parallel,
                                        num_threads=threads)
        elif name == "slither":
            result = mover.system_slither(current["lattice"], current["energy"], state.acc,
                                          state.ham, size, parallelize=parallel,
                                          num_threads=threads)
        else:
            result = mover.system_pull(current["lattice"], current["energy"], state.acc,
                                       state.ham, size, parallelize=parallel,
                                       num_threads=threads)
        current["lattice"], current["energy"] = result[0], result[1]

    out: dict[str, float | bool] = {}
    out["serial_wall"], out["serial_kernel"] = timed_runs(lambda: run(False, 1), repeats)
    for n_threads in THREADS:
        wall, kernel = timed_runs(lambda: run(True, n_threads), repeats)
        out[f"par{n_threads}_wall"] = wall
        out[f"par{n_threads}_kernel"] = kernel
    recomputed = int(state.ham.evaluate_total_energy(current["lattice"])[0])
    out["energy_ok"] = recomputed == int(current["energy"])
    return out


def main() -> int:
    """Run the harness and print one line per move, box and megamove size.

    Returns
    -------
    int
        0 if every energy check passed, 1 otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--boxes", type=int, nargs="+", default=list(BOXES),
                        help="square box side lengths to measure")
    parser.add_argument("--repeats", type=int, default=7,
                        help="timed megamoves per configuration (median is reported)")
    args = parser.parse_args()

    _wrap_kernels()
    print(f"OpenMP: {fk.openmp_info()}   threads timed: {THREADS}")
    all_ok = True
    for box in args.boxes:
        for name in ("crankshaft", "slither", "pull"):
            state, n_chains = build_state(box, SEQUENCES[name])
            for label, size in zip(("default", "heavy"), SIZES[name]):
                r = bench_move(name, state, size, args.repeats)
                all_ok &= bool(r["energy_ok"])
                serial_wall = float(r["serial_wall"])
                serial_kernel = float(r["serial_kernel"])
                cols = "  ".join(
                    f"{n}t {serial_wall / float(r[f'par{n}_wall']):.2f}x "
                    f"(kernel {serial_kernel / max(float(r[f'par{n}_kernel']), 1e-9):.2f}x)"
                    for n in THREADS)
                print(f"{name:10s} {box:4d}x{box:<4d} {n_chains:5d} chains  {label:7s} "
                      f"serial {serial_wall * 1e3:8.2f} ms (kernel {serial_kernel * 1e3:8.2f} ms)  "
                      f"{cols}  energy_ok={r['energy_ok']}", flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
