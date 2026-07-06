"""
Tests for the TRAJECTORY_PBC_UNWRAP option.

When on, every chain is written to the trajectory as a single "whole" periodic
image (not torn across a box face); coordinates may fall outside the box. When off
(the default) the raw on-lattice positions are written, so a chain that crosses a
boundary is split.
"""

import os

import numpy as np
import pytest
import mdtraj as md

from pimms import CONFIG
from pimms import lattice_utils
from pimms.tests import kernel_test_utils as U
from pimms.keyfile_parser import KeyFileParser


def _max_within_chain_step(chain, unwrap):
    """Largest single-axis step between consecutive beads of a chain (1 if whole)."""
    p = np.array(chain.get_output_positions(unwrap=unwrap))
    return int(np.abs(np.diff(p, axis=0)).max()) if len(p) > 1 else 0


def _straddling_lattice(tmp_path):
    # box with a short axis (8) and an 8-bead chain under PBC forces many chains to
    # wrap across the short-axis boundary
    _trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SR", hardwall=False, moves={"MOVE_CRANKSHAFT": 0.6, "MOVE_SLITHER": 0.4},
        box=[8, 8, 30], chains=[(30, "AABBAABB")],
        n_steps=200, energy_check=50, seed=3, temperature=60, return_sim=True)
    return sim.LATTICE


def test_default_is_off():
    assert CONFIG.DEFAULTS["TRAJECTORY_PBC_UNWRAP"] is False


def test_keyword_parses_true(tmp_path):
    kf = tmp_path / "KEYFILE.kf"
    kf.write_text("N_STEPS : 10\nTRAJECTORY_PBC_UNWRAP : True\n")
    parser = KeyFileParser(str(kf), parse_only=True)
    assert parser.keyword_lookup["TRAJECTORY_PBC_UNWRAP"] is True


def test_unwrap_makes_all_chains_whole(tmp_path):
    lattice = _straddling_lattice(tmp_path)

    # the setup must actually exercise the boundary (otherwise the test is vacuous)
    assert lattice.any_chains_straddle_boundary()
    assert any(c.does_chain_stradle_pbc_boundary() for c in lattice.chains.values())

    # raw output leaves at least one chain split across the boundary...
    assert any(_max_within_chain_step(c, unwrap=False) > 1 for c in lattice.chains.values())

    # ...but with unwrap every chain is contiguous (consecutive beads differ by 1)
    assert all(_max_within_chain_step(c, unwrap=True) <= 1 for c in lattice.chains.values())


def test_unwrap_may_place_beads_outside_the_box(tmp_path):
    lattice = _straddling_lattice(tmp_path)
    dims = lattice.dimensions
    outside = False
    for c in lattice.chains.values():
        p = np.array(c.get_output_positions(unwrap=True))
        if (p < 0).any() or any((p[:, d] >= dims[d]).any() for d in range(len(dims))):
            outside = True
            break
    # a whole chain wrapped across a face necessarily extends past the box edge
    assert outside


def test_make_chain_whole_anchors_first_bead_both_directions():
    dims = [10, 10, 10]
    # chain crosses the +x face (9 -> 0): stays anchored at 9, spills ABOVE the box
    up = lattice_utils.make_chain_whole([[9, 5, 5], [0, 5, 5], [1, 5, 5]], dims)
    assert up[0] == [9, 5, 5]                      # first bead unchanged (in place)
    assert up == [[9, 5, 5], [10, 5, 5], [11, 5, 5]]
    assert int(np.abs(np.diff(np.array(up), axis=0)).max()) == 1   # contiguous

    # chain crosses the -x face (0 -> 9): stays anchored at 0, spills BELOW the box
    down = lattice_utils.make_chain_whole([[0, 5, 5], [9, 5, 5], [8, 5, 5]], dims)
    assert down[0] == [0, 5, 5]
    assert down == [[0, 5, 5], [-1, 5, 5], [-2, 5, 5]]
    # crucially: this is NOT shifted to be non-negative (that would move the chain to
    # one side and is exactly the old one-sided behaviour we fixed)
    assert min(p[0] for p in down) < 0


def test_persistent_xtc_writer_roundtrip(tmp_path):
    # a real lattice, then write several frames through the persistent writer and
    # confirm the trajectory loads back with the right number of frames/atoms
    state = U.build_state(tmp_path, 3, "SR", hardwall=False, moves={"MOVE_CRANKSHAFT": 1.0},
                          box=[8, 8, 20], chains=[(20, "AABB")])
    lattice = state.lattice
    n_atoms = sum(len(c.get_ordered_positions()) for c in lattice.chains.values())

    cwd = os.getcwd()
    os.chdir(str(tmp_path))
    try:
        writer = lattice_utils.open_xtc_writer(lattice, lattice.lattice_to_angstroms,
                                               pdb_filename="START.pdb", xtc_filename="traj.xtc",
                                               unwrap=True)
        for _ in range(9):
            lattice_utils.write_xtc_frame(writer, lattice, lattice.lattice_to_angstroms, unwrap=True)
        lattice_utils.close_xtc_writer(writer)
        traj = md.load("traj.xtc", top="START.pdb")
    finally:
        os.chdir(cwd)

    assert traj.n_frames == 10          # 1 initial (open) + 9 appended
    assert traj.n_atoms == n_atoms


def test_autocenter_takes_precedence_over_unwrap(tmp_path):
    # a single-chain lattice: autocenter should win (single-image AND centred)
    _trace, sim = U.run_sim_with_energy_check(
        tmp_path, 3, "SR", hardwall=False, moves={"MOVE_CRANKSHAFT": 1.0},
        box=[8, 8, 20], chains=[(1, "AABBAABB")],
        n_steps=100, energy_check=50, seed=1, temperature=60, return_sim=True)
    chain = list(sim.LATTICE.chains.values())[0]
    both = np.array(chain.get_output_positions(autocenter=True, unwrap=True))
    centred_only = np.array(chain.get_ordered_positions(center_positions=True))
    assert np.array_equal(both, centred_only)
