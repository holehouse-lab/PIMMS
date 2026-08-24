"""
Tests for capillary-wave surface-tension estimation.

Surface tension is inherently noisy for small lattice condensates, so these assert
the machinery, the units plumbing (temperature = kT), and a *positive, finite*
estimate on a genuinely phase-separated slab / droplet - not a precise value.
"""

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import surface_tension as st


def test_temperature_is_loaded(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    assert traj.temperature == pytest.approx(22.0)


def test_slab_surface_tension_positive(traj_slab_files):
    xtc, pdb, keyfile = traj_slab_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    result = st.slab_surface_tension(traj)
    assert isinstance(result, st.SurfaceTension)
    assert result.method == "slab"
    assert result.temperature == pytest.approx(26.0)
    assert np.isfinite(result.gamma) and result.gamma > 0
    assert result.n_modes > 0


def test_slab_counts_equal_power_wavevectors_as_distinct_modes(monkeypatch):
    class Cluster:
        n_beads = 100
        positions = np.zeros((100, 3), dtype=float)

    class Frame:
        clusters = [Cluster()]

    class Trajectory:
        dimensions = np.array([8, 8, 20])
        n_dim = 3
        n_frames = 1
        temperature = 2.0

        def __getitem__(self, index):
            return Frame()

    x = np.arange(8)[:, None]
    y = np.arange(8)[None, :]
    heights = np.cos(2 * np.pi * x / 8) + np.cos(2 * np.pi * y / 8)
    monkeypatch.setattr(st, "_interface_heights", lambda *args: (heights, None))

    result = st.slab_surface_tension(Trajectory(), n_modes=2)

    # The x and y modes happen to have identical power, but they are independent
    # wavevectors. Value-based np.unique used to collapse them to one mode.
    assert result.n_modes == 2
    assert len(result.spectrum[0]) == 2
    assert np.isfinite(result.gamma)


@pytest.mark.parametrize("n_modes", [0, -1, 1.5, True])
def test_slab_rejects_invalid_mode_count(traj_slab_files, n_modes):
    xtc, pdb, keyfile = traj_slab_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    with pytest.raises(ValueError, match="positive integer"):
        st.slab_surface_tension(traj, n_modes=n_modes)


def test_droplet_surface_tension_runs(traj_condensed_files):
    xtc, pdb, keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    result = st.droplet_surface_tension(traj, l_max=4, min_beads=20)
    assert result.method == "droplet"
    assert result.temperature == pytest.approx(22.0)
    assert np.isfinite(result.gamma)
    # spectrum is (l, <|u_l|^2>) for l = 2..l_max
    ls, u2 = result.spectrum
    assert list(ls) == [2, 3, 4]
    assert np.all(u2 > 0)


def test_dispatch_by_geometry(traj_condensed_files, traj_slab_files):
    cxtc, cpdb, ckey = traj_condensed_files
    sxtc, spdb, skey = traj_slab_files
    cubic = lemonade.load(xtc=cxtc, pdb=cpdb, keyfile=ckey)
    slab = lemonade.load(xtc=sxtc, pdb=spdb, keyfile=skey)
    assert st.surface_tension(cubic).method == "droplet"      # cubic box
    assert st.surface_tension(slab).method == "slab"          # elongated box


def test_requires_temperature(traj_condensed_files):
    xtc, pdb, _keyfile = traj_condensed_files
    traj = lemonade.load(xtc=xtc, pdb=pdb)                    # no keyfile -> no temperature
    assert traj.temperature is None
    with pytest.raises(ValueError, match="temperature"):
        st.droplet_surface_tension(traj)
    # explicit temperature override works
    result = st.droplet_surface_tension(traj, temperature=22.0, l_max=4, min_beads=20)
    assert result.temperature == pytest.approx(22.0)
