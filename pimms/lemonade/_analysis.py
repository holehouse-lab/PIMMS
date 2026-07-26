"""
Vectorised, batched numeric core for lemonade.

Every function here operates on a whole trajectory at once - shape
``(n_frames, n_atoms, k)`` position arrays with a CSR-style ``offsets`` array
delimiting the atoms of each chain - and returns per-frame-per-chain results with
no Python-level loop over frames or chains. ``numpy.add.reduceat`` performs the
per-chain reductions in one C call.

These expect ``whole`` positions (each chain already made contiguous across
periodic boundaries), so intra-chain distances are ordinary Euclidean distances.
"""

import numpy as np


def centers_of_mass(whole, offsets, lengths):
    """Per-chain centre of mass. ``whole`` is ``(nf, na, k)`` -> ``(nf, nc, k)``."""
    sums = np.add.reduceat(whole, offsets[:-1], axis=1)
    return sums / lengths[np.newaxis, :, np.newaxis]


def _centered(whole, offsets, lengths, com):
    if com is None:
        com = centers_of_mass(whole, offsets, lengths)
    return whole - np.repeat(com, lengths, axis=1)


def radius_of_gyration(whole, offsets, lengths, com=None):
    """Per-chain radius of gyration, ``(nf, nc)``.

    ``Rg = sqrt( <|r_i - r_com|^2> )`` = ``sqrt(trace(gyration tensor))``.
    """
    d = _centered(whole, offsets, lengths, com)
    sq = np.einsum("fak,fak->fa", d, d)
    rg2 = np.add.reduceat(sq, offsets[:-1], axis=1) / lengths[np.newaxis, :]
    return np.sqrt(rg2)


def gyration_eigenvalues(whole, offsets, lengths, com=None):
    """Ascending eigenvalues of each chain's gyration tensor, ``(nf, nc, k)``.

    The tensor is accumulated one component at a time. Forming the full
    ``(nf, na, k, k)`` outer-product array first and reducing that costs ``k * k``
    times the size of the position array - 0.7 GB for a 1000-frame, 10k-bead
    trajectory and 5.8 GB at 2000 frames / 40k beads, which is enough to put loading a
    long trajectory out of reach on an ordinary machine. Reducing per component holds
    only one ``(nf, na)`` scratch array at a time, and symmetry means just
    ``k(k+1)/2`` of them need computing (6 rather than 9 in 3D).

    ``np.add.reduceat`` walks each chain's atoms in the same order either way, so the
    result is bit-identical to the full-array version.
    """
    d = _centered(whole, offsets, lengths, com)
    n_frames, _n_atoms, k = d.shape
    n_chains = len(lengths)

    tensor = np.empty((n_frames, n_chains, k, k), dtype=np.float64)
    norm = lengths[np.newaxis, :]

    for i in range(k):
        for j in range(i, k):
            component = np.add.reduceat(d[:, :, i] * d[:, :, j], offsets[:-1], axis=1) / norm
            tensor[:, :, i, j] = component
            if i != j:
                tensor[:, :, j, i] = component

    return np.linalg.eigvalsh(tensor)                                # ascending


def asphericity(eigenvalues):
    """Asphericity from ascending gyration eigenvalues.

    3D: ``lam_z - 0.5 (lam_x + lam_y)`` with ``lam_z`` the largest.
    2D: ``lam_max - lam_min``.
    Zero for a perfectly spherical/circular distribution.
    """
    if eigenvalues.shape[-1] == 3:
        return eigenvalues[..., 2] - 0.5 * (eigenvalues[..., 0] + eigenvalues[..., 1])
    return eigenvalues[..., -1] - eigenvalues[..., 0]


def end_to_end(whole, offsets, lengths):
    """Per-chain end-to-end distance, ``(nf, nc)``."""
    first = whole[:, offsets[:-1], :]
    last = whole[:, offsets[1:] - 1, :]
    d = last - first
    return np.sqrt(np.einsum("fck,fck->fc", d, d))


def distance_map(chain_positions):
    """Full inter-bead Euclidean distance matrix for one chain's whole positions.

    ``chain_positions`` is ``(L, k)`` -> ``(L, L)``.
    """
    diff = chain_positions[:, np.newaxis, :] - chain_positions[np.newaxis, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))


def internal_scaling(chain_positions):
    """Mean inter-bead distance as a function of sequence separation |i-j|.

    Returns ``(separations, mean_distance)`` where separation runs ``1 .. L-1``.
    """
    dmap = distance_map(chain_positions)
    n = dmap.shape[0]
    seps = np.arange(1, n)
    means = np.empty(n - 1, dtype=np.float64)
    for s in seps:
        means[s - 1] = np.mean(np.diagonal(dmap, offset=s))
    return seps, means
