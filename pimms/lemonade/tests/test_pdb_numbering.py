"""lemonade must not depend on PDB atom serials or residue numbers.

PIMMS numbers residues from 1 within every chain and writes atom serials
modulo 100 000, so a large system's START.pdb carries duplicated serials (and
duplicated residue IDs in every file). lemonade reads chains from the PDB's
chain/TER structure and atom order, never from those columns, so a PDB whose
numbering restarts in every chain, or whose serials wrap mid-chain, must load
and analyse identically to a uniquely numbered one.
"""

import contextlib
import os
import re

import numpy as np
import pytest

from pimms import lemonade
from pimms.lemonade import phase_separation as ps
from pimms.tests import kernel_test_utils as U


def _renumber(pdb_text, mode):
    """Rewrite the serial and residue-number columns of a PIMMS PDB.

    ``per_chain``: serials and residue numbers restart at 1 in every chain (the
    user's scenario). ``wrapped``: serials follow PIMMS's own modulo rule from a
    start that crosses 99999 within the first chain, and residue numbers roll
    over past 9999 within the second chain.
    """
    out = []
    serial = 99995 if mode == "wrapped" else 0
    resid = 9995 if mode == "wrapped" else 0
    for line in pdb_text.splitlines(keepends=True):
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            serial += 1
            resid += 1
            if mode == "wrapped" and resid > 9999:
                resid = 1
            s = serial % 100000 if mode == "wrapped" else serial
            line = line[:6] + f"{s:5d}" + line[11:22] + f"{resid:4d}" + line[26:]
        elif rec == "TER   ":
            serial += 1
            s = serial % 100000 if mode == "wrapped" else serial
            line = line[:6] + f"{s:5d}" + line[11:22] + f"{resid:4d}" + line[26:]
            if mode == "per_chain":
                serial = 0
                resid = 0
        elif rec == "CONECT":
            # PIMMS skips CONECT records past the wrap; drop them all here so the
            # renumbered files stay self-consistent (lemonade never reads bonds)
            continue
        out.append(line)
    return "".join(out)


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    """A short PIMMS run: 5 chains of 10 beads, START.pdb + traj.xtc + keyfile."""
    d = tmp_path_factory.mktemp("pdbnum")
    st = U.build_state(d, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                       box=[16, 16, 16], chains=[(5, "AABBAABBAA")], n_steps=8,
                       equilibration=0, temperature=40,
                       extra={"XTC_FREQ": 2, "ANA_POL": 2, "PRINT_FREQ": 1000,
                              "ENERGY_CHECK": 0, "SAVE_EQ": "True"})
    cwd = os.getcwd()
    os.chdir(str(d))
    try:
        with contextlib.redirect_stdout(open(os.devnull, "w")):
            st.sim.run_simulation()
    finally:
        os.chdir(cwd)
    text = (d / "START.pdb").read_text()
    (d / "per_chain.pdb").write_text(_renumber(text, "per_chain"))
    (d / "wrapped.pdb").write_text(_renumber(text, "wrapped"))
    return d


def _summary(traj):
    return {
        "n_chains": traj.n_chains,
        "lengths": traj.topology.lengths.copy(),
        "sequences": list(traj.topology.sequences),
        "types": traj.chain_types.copy(),
        "positions": np.asarray(traj.positions).copy(),
        "rg": traj.radius_of_gyration().copy(),
        "e2e": traj.end_to_end_distance().copy(),
        "asph": traj.asphericity().copy(),
        "clusters": [sorted(sorted(int(i) for i in c.chain_indices) for c in frame.clusters)
                     for frame in traj],
    }


def _load(run_dir, pdb, keyfile):
    kw = {"keyfile": str(run_dir / "KEYFILE.kf")} if keyfile else {"spacing": 3.65}
    return lemonade.load(xtc=str(run_dir / "traj.xtc"), pdb=str(run_dir / pdb), **kw)


def test_renumbered_pdbs_differ_only_in_the_numbering_columns(run_dir):
    ref = (run_dir / "START.pdb").read_text()
    for name in ("per_chain.pdb", "wrapped.pdb"):
        alt = (run_dir / name).read_text()
        ref_atoms = [line for line in ref.splitlines() if line.startswith("ATOM")]
        alt_atoms = [line for line in alt.splitlines() if line.startswith("ATOM")]
        assert len(ref_atoms) == len(alt_atoms) == 50
        # coordinates, names, chain letters untouched
        assert [line[12:22] + line[26:] for line in ref_atoms] == [line[12:22] + line[26:] for line in alt_atoms]
        assert [line[6:11] for line in ref_atoms] != [line[6:11] for line in alt_atoms]
    per = [line for line in (run_dir / "per_chain.pdb").read_text().splitlines() if line.startswith("ATOM")]
    assert [int(line[6:11]) for line in per] == list(range(1, 11)) * 5
    assert [int(line[22:26]) for line in per] == list(range(1, 11)) * 5
    wrapped = [line for line in (run_dir / "wrapped.pdb").read_text().splitlines() if line.startswith("ATOM")]
    serials = [int(line[6:11]) for line in wrapped]
    assert serials[:6] == [99996, 99997, 99998, 99999, 0, 1]
    resids = [int(line[22:26]) for line in wrapped]
    assert 9999 in resids and resids[resids.index(9999) + 1] == 1


@pytest.mark.parametrize("keyfile", [True, False], ids=["with_keyfile", "no_keyfile"])
@pytest.mark.parametrize("pdb", ["per_chain.pdb", "wrapped.pdb"])
def test_lemonade_ignores_serials_and_residue_numbers(run_dir, pdb, keyfile):
    ref = _summary(_load(run_dir, "START.pdb", keyfile))
    alt = _summary(_load(run_dir, pdb, keyfile))
    assert alt["n_chains"] == ref["n_chains"] == 5
    np.testing.assert_array_equal(alt["lengths"], ref["lengths"])
    assert alt["sequences"] == ref["sequences"] == ["AABBAABBAA"] * 5
    np.testing.assert_array_equal(alt["types"], ref["types"])
    np.testing.assert_array_equal(alt["positions"], ref["positions"])
    for key in ("rg", "e2e", "asph"):
        np.testing.assert_array_equal(alt[key], ref[key])
    assert alt["clusters"] == ref["clusters"]


@pytest.mark.parametrize("pdb", ["per_chain.pdb", "wrapped.pdb"])
def test_phase_separation_analysis_is_unaffected(run_dir, pdb):
    ref = ps.analyze(_load(run_dir, "START.pdb", True))
    alt = ps.analyze(_load(run_dir, pdb, True))
    assert alt.is_phase_separated == ref.is_phase_separated
    assert alt.percolation_fraction == ref.percolation_fraction
    assert alt.n_clusters == ref.n_clusters
    np.testing.assert_array_equal(np.asarray(alt.condensed_fraction_series),
                                  np.asarray(ref.condensed_fraction_series))


def test_pdb_only_load_matches_the_trajectory_topology(run_dir):
    """A renumbered PDB on its own (no XTC) still gives one chain per TER block."""
    for pdb in ("per_chain.pdb", "wrapped.pdb"):
        single = lemonade.load(pdb=str(run_dir / pdb), keyfile=str(run_dir / "KEYFILE.kf"))
        assert single.n_chains == 5 and single.n_frames == 1
        assert list(single.topology.lengths) == [10] * 5
        full = _load(run_dir, "START.pdb", True)
        np.testing.assert_array_equal(np.asarray(single.positions)[0], np.asarray(full.positions)[0])


def test_keyfile_chain_block_is_matched_regardless_of_numbering(run_dir, recwarn):
    """The keyfile CHAIN block must still be recognised (no mismatch warning)."""
    for pdb in ("per_chain.pdb", "wrapped.pdb"):
        recwarn.clear()
        _load(run_dir, pdb, True)
        assert not [w for w in recwarn if re.search("does not reproduce|could not be applied", str(w.message))]
