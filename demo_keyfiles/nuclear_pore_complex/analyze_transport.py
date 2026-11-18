#!/usr/bin/env python
"""
Follow the transport receptors and the inert cargo through a nuclear pore run.

Reads ``START.pdb`` and ``traj.xtc`` (written by PIMMS) and ``npc.restart``
(written by ``build_npc.py``, to know which chains are which), and reports, frame
by frame, how many receptors and how many cargo chains are in the cytoplasm, in
the central channel and in the nucleus. A chain's place is set by its centre's
height: within the inner ring's height it's in the channel (the inner ring fills
the rest of the pore), above it in the cytoplasm, below it in the nucleus - so a
chain counts as through once it has crossed the inner ring.

Run from the directory holding the run's output:  ``python analyze_transport.py``
Writes the counts per frame to ``transport.csv`` as well.
"""

from __future__ import annotations

import csv
import pickle
from typing import Dict, List, Tuple

import mdtraj  # type: ignore[import-untyped]
import numpy as np

import build_npc

TRAJECTORY_FILENAME = "traj.xtc"
TOPOLOGY_FILENAME = "START.pdb"
CSV_FILENAME = "transport.csv"
NM_PER_ANGSTROM = 0.1
# the channel is the inner ring's height; a chain whose centre is below it
# has crossed the permeability barrier (in sites from the envelope mid-plane)
CHANNEL_HALF_HEIGHT = build_npc.IR_HALF_HEIGHT
REGIONS: Tuple[str, ...] = ("cytoplasm", "channel", "nucleus")
N_REPORT_ROWS = 11


def chain_atom_ranges(restart_filename: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Atom index ranges of the receptor and cargo chains in the PDB/XTC.

    PIMMS writes chains in ascending chainID order, so a chain's atoms start
    after those of every chain with a smaller ID.

    Parameters
    ----------
    restart_filename : str
        The restart file the run started from.

    Returns
    -------
    dict
        "receptor" and "cargo" -> list of (first atom, one past the last atom).
    """
    with open(restart_filename, "rb") as handle:
        chains = pickle.load(handle)["CHAINS"]
    ranges: Dict[str, List[Tuple[int, int]]] = {"receptor": [], "cargo": []}
    first = 0
    for chain_id in sorted(chains):
        sequence = chains[chain_id][1]
        if set(sequence) == {build_npc.RECEPTOR}:
            ranges["receptor"].append((first, first + len(sequence)))
        elif set(sequence) == {build_npc.CARGO}:
            ranges["cargo"].append((first, first + len(sequence)))
        first += len(sequence)
    return ranges


def check_ranges(topology: mdtraj.Topology, ranges: Dict[str, List[Tuple[int, int]]]) -> None:
    """
    Check the atom ranges hold the residues they should (LYS receptors, ASP cargo).

    Parameters
    ----------
    topology : mdtraj.Topology
        The run's topology.
    ranges : dict
        From :func:`chain_atom_ranges`.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If the topology doesn't match the restart file (a different build?).
    """
    expected = {"receptor": "LYS", "cargo": "ASP"}
    for kind, spans in ranges.items():
        for first, last in spans:
            names = {topology.atom(index).residue.name for index in range(first, last)}
            if names != {expected[kind]}:
                raise RuntimeError("atoms %d-%d should all be %s but are %s: does %s match this run?"
                                   % (first, last - 1, expected[kind], sorted(names), build_npc.RESTART_FILENAME))


def region_counts(z_sites: np.ndarray) -> np.ndarray:
    """
    Count chains in the cytoplasm, the channel and the nucleus.

    Parameters
    ----------
    z_sites : numpy.ndarray
        Shape (n_frames, n_chains): each chain's centre height, in sites from
        the envelope mid-plane.

    Returns
    -------
    numpy.ndarray
        Shape (n_frames, 3): counts in cytoplasm, channel, nucleus.
    """
    cytoplasm = (z_sites > CHANNEL_HALF_HEIGHT).sum(axis=1)
    nucleus = (z_sites < -CHANNEL_HALF_HEIGHT).sum(axis=1)
    channel = z_sites.shape[1] - cytoplasm - nucleus
    return np.stack([cytoplasm, channel, nucleus], axis=1)


def main() -> None:
    """
    Report where the receptors and the cargo are over the run.

    Returns
    -------
    None
    """
    trajectory = mdtraj.load(TRAJECTORY_FILENAME, top=TOPOLOGY_FILENAME)
    ranges = chain_atom_ranges(build_npc.RESTART_FILENAME)
    check_ranges(trajectory.topology, ranges)
    # PDB/XTC coordinates are lattice coordinates x LATTICE_TO_ANGSTROMS; mdtraj reads nm
    site_nm = build_npc.NM_PER_SITE
    centre_z: Dict[str, np.ndarray] = {}
    for kind, spans in ranges.items():
        heights = [trajectory.xyz[:, first:last, 2].mean(axis=1) / site_nm + build_npc.BOX_Z[0] for first, last in spans]
        centre_z[kind] = np.stack(heights, axis=1)
    counts = {kind: region_counts(z) for kind, z in centre_z.items()}

    n_frames = trajectory.n_frames
    print("%d frames; %d receptors and %d cargo chains" % (n_frames, len(ranges["receptor"]), len(ranges["cargo"])))
    print("chains in the cytoplasm / channel / nucleus:")
    print("  frame    receptors (cyt/chan/nuc)    cargo (cyt/chan/nuc)")
    for frame in np.unique(np.linspace(0, n_frames - 1, N_REPORT_ROWS).astype(int)):
        receptor = counts["receptor"][frame]
        cargo = counts["cargo"][frame]
        print("  %5d    %6d %6d %6d          %6d %6d %6d" % (frame, *receptor, *cargo))
    for kind in ("receptor", "cargo"):
        ever_in_nucleus = int((centre_z[kind] < -CHANNEL_HALF_HEIGHT).any(axis=0).sum())
        ever_in_channel = int((centre_z[kind] <= CHANNEL_HALF_HEIGHT).any(axis=0).sum())
        print("%-8s: %2d of %2d entered the channel, %2d reached the nucleus"
              % (kind, ever_in_channel, centre_z[kind].shape[1], ever_in_nucleus))

    with open(CSV_FILENAME, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frame"] + ["%s_%s" % (kind, region) for kind in ("receptor", "cargo") for region in REGIONS])
        for frame in range(n_frames):
            writer.writerow([frame] + list(counts["receptor"][frame]) + list(counts["cargo"][frame]))
    print("wrote %s" % CSV_FILENAME)


if __name__ == "__main__":
    main()
