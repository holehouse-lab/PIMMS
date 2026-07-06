"""
Tests for the phase-separation / droplet-physics analysis.

Machinery (shapes, bounds, fits) is checked on the general 3D fixture; the physics
(a real dense/dilute split, elevated condensed fraction) is checked on a strongly
self-attracting system that phase separates.
"""

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import phase_separation as ps


# ---------------------------------------------------------------------------
# order parameters
# ---------------------------------------------------------------------------

def test_order_parameter_shapes_and_bounds(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)

    cf = ps.condensed_fraction(traj)
    assert cf.shape == (traj.n_frames,)
    assert np.all((cf >= 0) & (cf <= 1))

    nc = ps.number_of_clusters(traj, min_beads=1)
    assert nc.shape == (traj.n_frames,)
    assert np.all(nc >= 1)                      # at least one cluster per frame

    lc = ps.largest_cluster_size(traj, by="beads")
    assert lc.shape == (traj.n_frames,)
    assert np.all(lc <= traj.n_atoms)

    sizes = ps.cluster_size_distribution(traj)
    assert sizes.ndim == 1 and sizes.sum() > 0


# ---------------------------------------------------------------------------
# density profiles are occupied fractions in [0, 1]
# ---------------------------------------------------------------------------

def test_radial_profile_is_occupied_fraction(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    r, rho = ps.radial_density_profile(traj)
    assert r.shape == rho.shape
    assert np.all(rho >= -1e-9) and np.all(rho <= 1.0 + 1e-9)
    # a condensate: dense near the COM, dilute far away
    assert rho[0] > rho[-1]


def test_slab_profile_bounds_and_length(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    z, rho = ps.slab_density_profile(traj, axis=2)
    assert len(z) == traj.dimensions[2]
    assert np.all(rho >= -1e-9) and np.all(rho <= 1.0 + 1e-9)


# ---------------------------------------------------------------------------
# tanh fits
# ---------------------------------------------------------------------------

def test_fits_return_ordered_binodal(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)

    r, rho = ps.radial_density_profile(traj)
    fit = ps.fit_radial_profile(r, rho)
    assert isinstance(fit, ps.BinodalFit)
    assert 0 <= fit.rho_dilute <= fit.rho_dense <= 1
    z, rhoz = ps.slab_density_profile(traj, axis=2)
    sfit = ps.fit_slab_profile(z, rhoz)
    assert 0 <= sfit.rho_dilute <= sfit.rho_dense <= 1


# ---------------------------------------------------------------------------
# full analysis on a phase-separated system
# ---------------------------------------------------------------------------

def test_analyze_detects_phase_separation(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    result = ps.analyze(traj)

    assert result.geometry == "sphere"                      # cubic box
    assert 0 <= result.rho_dilute <= result.rho_dense <= 1
    assert result.rho_dense > 2 * max(result.rho_dilute, 1e-6)  # a real density gap
    assert result.condensed_fraction > 0.3                  # most material condensed
    assert result.largest_cluster_beads > 0
    # droplet shape is computed
    assert np.isfinite(result.shape["radius_of_gyration"])
    assert result.is_phase_separated


def test_analyze_geometry_override(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    slab_result = ps.analyze(traj, geometry="slab")
    assert slab_result.geometry == "slab"
    assert 0 <= slab_result.rho_dilute <= slab_result.rho_dense <= 1


# ---------------------------------------------------------------------------
# droplet shape
# ---------------------------------------------------------------------------

def test_droplet_and_sphericity(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    droplet = traj[-1].droplet
    assert droplet is not None
    assert droplet.n_beads > 0
    sph = droplet.sphericity
    # sphericity is either a valid (0, 1.3] value or nan for a degenerate hull
    assert np.isnan(sph) or (0 < sph < 1.3)
