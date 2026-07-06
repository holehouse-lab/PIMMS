"""
Shared fixtures for the lemonade tests.

Each fixture runs a short PIMMS simulation into a temp directory (producing
``START.pdb`` + ``traj.xtc`` + the keyfile) so the tests exercise the real
end-to-end load path against genuine PIMMS output. Fixtures are session-scoped so
the simulations run once.
"""

import os
import contextlib

import pytest


_PARAMS = (
    "ANGLE_PENALTY\tA\t0\t0\t0\n"
    "ANGLE_PENALTY\tB\t0\t0\t0\n"
    "A\tA\t-8\n"
    "B\tB\t-8\n"
    "A\tB\t2\n"
    "A\t0\t0\n"
    "B\t0\t0\n"
)


_STICKY = (
    "ANGLE_PENALTY\tS\t0\t0\t0\n"
    "S\tS\t-16\n"
    "S\t0\t0\n"
)


def _make_trajectory(dirpath, box, chains, *, seed=7, n_steps=300, xtc_freq=15,
                     spacing=None, params=_PARAMS, temperature=55, equilibration=40,
                     resized_equilibration=None):
    """Run a small PIMMS sim in ``dirpath`` and return (xtc, pdb, keyfile) paths."""
    dirpath = str(dirpath)
    with open(os.path.join(dirpath, "params.prm"), "w") as fh:
        fh.write(params)
    lines = [
        "DIMENSIONS : " + " ".join(str(b) for b in box),
        "PARAMETER_FILE : params.prm",
        "HARDWALL : False",
        f"TEMPERATURE : {temperature}",
        f"N_STEPS : {n_steps}",
        f"EQUILIBRATION : {equilibration}",
        f"SEED : {seed}",
        f"XTC_FREQ : {xtc_freq}",
        "PRINT_FREQ : 100000",
        "EN_FREQ : 100000",
        "RESTART_FREQ : 100000",
        "ANALYSIS_FREQ : 100000",
        "MOVE_CRANKSHAFT : 0.6",
        "MOVE_SLITHER : 0.4",
        "TRAJECTORY_PBC_UNWRAP : False",
    ]
    if spacing is not None:
        lines.append(f"LATTICE_TO_ANGSTROMS : {spacing}")
    if resized_equilibration is not None:
        lines.append("RESIZED_EQUILIBRATION : " + " ".join(str(b) for b in resized_equilibration))
    for count, seq in chains:
        lines.append(f"CHAIN : {count} {seq}")
    keyfile = os.path.join(dirpath, "KEYFILE.kf")
    with open(keyfile, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    cwd = os.getcwd()
    os.chdir(dirpath)
    try:
        from pimms.keyfile_parser import KeyFileParser
        from pimms.simulation import Simulation
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            Simulation(KeyFileParser("KEYFILE.kf").keyword_lookup).run_simulation()
    finally:
        os.chdir(cwd)

    return (os.path.join(dirpath, "traj.xtc"),
            os.path.join(dirpath, "START.pdb"),
            keyfile)


@pytest.fixture(scope="session")
def traj3d_files(tmp_path_factory):
    """A 3D trajectory in a small box (8^3) - chains straddle the boundary heavily."""
    d = tmp_path_factory.mktemp("lemon3d")
    return _make_trajectory(d, [8, 8, 8], [(4, "AABBAABBAABB")], seed=7)


@pytest.fixture(scope="session")
def traj2d_files(tmp_path_factory):
    """A 2D trajectory."""
    d = tmp_path_factory.mktemp("lemon2d")
    return _make_trajectory(d, [10, 10], [(6, "AABBAA")], seed=11)


@pytest.fixture(scope="session")
def traj_dilute_files(tmp_path_factory):
    """Short chains in a large box: chains stay small vs the box, so there are no
    finite-size artefacts and lemonade must agree with PIMMS's own Rg."""
    d = tmp_path_factory.mktemp("lemondilute")
    return _make_trajectory(d, [16, 16, 16], [(5, "AABB")], seed=5)


@pytest.fixture(scope="session")
def traj_condensed_files(tmp_path_factory):
    """A strongly self-attracting homopolymer in a small cold box - it phase
    separates into a dominant condensate, for the phase-separation tests."""
    d = tmp_path_factory.mktemp("lemoncond")
    return _make_trajectory(d, [12, 12, 12], [(70, "SSSS")], seed=3,
                            params=_STICKY, temperature=22, n_steps=3000,
                            equilibration=1000, xtc_freq=100)


@pytest.fixture(scope="session")
def traj_slab_files(tmp_path_factory):
    """A slab condensate in an elongated box (equilibrated compact, then grown along
    z) - two flat interfaces for the capillary-wave surface-tension test."""
    d = tmp_path_factory.mktemp("lemonslab")
    return _make_trajectory(d, [12, 12, 34], [(150, "SSSS")], seed=2,
                            params=_STICKY, temperature=26, n_steps=5000,
                            equilibration=1800, xtc_freq=150,
                            resized_equilibration=[12, 12, 12])
