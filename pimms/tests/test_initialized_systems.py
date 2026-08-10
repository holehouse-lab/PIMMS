"""
Tests for the NeurofilamentDemo programmatic-construction example.

Regression tests for a real bug: the demo had rotted against the evolving Chain and
Lattice constructors (it predates the long-range / chain-type arguments and the
``lattice_to_angstroms`` parameter), so instantiating it raised ``TypeError`` before
building anything. It is documented as the worked example of programmatic system
construction, so it has to actually run; these tests exercise it at a small scale
(the class is now parametrised rather than hardcoded to a 500^3 box).
"""

import contextlib
import io
import os

import numpy as np
import pytest

from pimms import initialized_systems
from pimms.energy import Hamiltonian
from pimms.latticeExceptions import CustomInitializationException


@pytest.fixture
def hamiltonian(tmp_path, monkeypatch):
    """A minimal Hamiltonian defining the 'E' bead the demo builds with."""
    monkeypatch.chdir(tmp_path)          # parameters_used.prm etc. land in tmp
    prm = tmp_path / "params.prm"
    prm.write_text("ANGLE_PENALTY\tE\t0\t0\t0\nE\tE\t-1\nE\t0\t0\n")
    with contextlib.redirect_stdout(io.StringIO()):
        return Hamiltonian(str(prm), 3, non_interacting=False, angles_off=False,
                           temperature=50, reduced_printing=True)


def test_neurofilament_demo_builds_a_consistent_lattice(hamiltonian, tmp_path):
    with contextlib.redirect_stdout(io.StringIO()):
        demo = initialized_systems.NeurofilamentDemo(
            hamiltonian, dimensions=[48, 48, 48], sidearm_length=10,
            sidearm_z_spacing=8, write_pdb=True, verbose=False)

    lattice = demo.LATTICE

    # chain 1 is the filament tube; every 8th z-layer contributes one sidearm
    assert 1 in lattice.chains
    assert lattice.get_number_of_chains() == 1 + 48 // 8

    # the tube is a hollow 10x10 square running the whole z axis: 36 sites per layer
    filament = lattice.chains[1]
    assert len(filament) == 36 * 48
    assert filament.fixed is True

    # grid and chain bookkeeping must agree: every bead's site carries its chainID
    for chainID, chain in lattice.chains.items():
        for pos in chain.get_ordered_positions():
            assert lattice.grid[pos[0]][pos[1]][pos[2]] == chainID

    # total occupancy matches the sum of chain lengths (no overlaps, no strays)
    total_beads = sum(len(c) for c in lattice.chains.values())
    assert int(np.count_nonzero(lattice.grid)) == total_beads
    assert int(np.count_nonzero(lattice.type_grid)) == total_beads

    # the worked example writes its PDB
    assert os.path.exists(tmp_path / "NEUROFILAMENT.pdb")


def test_neurofilament_demo_rejects_sidearms_that_cannot_fit(hamiltonian):
    with pytest.raises(CustomInitializationException, match="do not fit"):
        initialized_systems.NeurofilamentDemo(
            hamiltonian, dimensions=[48, 48, 48], sidearm_length=100,
            write_pdb=False, verbose=False)
