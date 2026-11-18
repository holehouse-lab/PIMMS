"""Fixes from the October 2026 survey of lemonade.

Each test here failed on the tree before its fix: the ``eq_START.pdb`` topology
forced hard walls onto the production trajectory it was paired with, a profile
with no qualifying cluster came back all-``nan`` without a word, a bool was a
valid frame or chain index, and the chain-extent sampler visited up to twice
the frames it said it would.
"""

import os
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

import pimms.lemonade as lemonade
from pimms.lemonade import phase_separation as ps
from pimms.lemonade.tests.test_review5_fixes import resized_equilibration_run  # noqa: F401 (fixture)


def test_eq_topology_with_the_production_trajectory_is_loaded_as_production(resized_equilibration_run):  # noqa: F811
    """``eq_START.pdb`` and ``START.pdb`` hold the same topology, so the former is a
    legitimate topology for ``traj.xtc``; the trajectory's name decides the phase."""
    sim, rundir, keyfile = resized_equilibration_run
    xtc = os.path.join(rundir, "traj.xtc")
    eq_pdb = os.path.join(rundir, "eq_START.pdb")
    reference = lemonade.load(xtc=xtc, pdb=os.path.join(rundir, "START.pdb"), keyfile=keyfile)
    assert reference.hardwall is False
    assert reference.dimensions == (14, 14, 14)

    for kwargs in ({"keyfile": keyfile}, {}):
        with pytest.warns(UserWarning, match="loaded as the PRODUCTION phase"):
            traj = lemonade.load(xtc=xtc, pdb=eq_pdb, **kwargs)
        assert traj.hardwall is False
        assert traj.dimensions == (14, 14, 14)
        assert np.array_equal(traj.positions, reference.positions)
        assert (list(ps.number_of_clusters(traj, min_beads=1))
                == list(ps.number_of_clusters(reference, min_beads=1)))

    # the eq_START.pdb on its own is still the compact, hardwall equilibration phase
    with pytest.warns(UserWarning, match="resized-equilibration"):
        alone = lemonade.load(pdb=eq_pdb, keyfile=keyfile)
    assert alone.hardwall is True
    assert alone.dimensions == (9, 9, 9)


def test_profiles_say_so_when_no_cluster_reaches_min_beads(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    with pytest.warns(UserWarning, match="no frame holds a cluster of at least 1000 beads"):
        _radii, density = ps.radial_density_profile(traj, min_beads=1000)
    assert np.all(np.isnan(density))
    with pytest.warns(UserWarning, match="droplet_shape: no frame holds a cluster"):
        shape = ps.droplet_shape(traj, min_beads=1000)
    assert all(np.isnan(value) for value in shape.values())

    # a threshold the clusters do reach raises no such warning
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ps.radial_density_profile(traj, min_beads=1)
        ps.droplet_shape(traj, min_beads=1)
    assert not [w for w in caught if "no frame holds" in str(w.message)]


def test_a_bool_is_not_a_frame_or_chain_index(traj3d_files):
    xtc, pdb, keyfile = traj3d_files
    traj = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    for key in (True, False, np.bool_(True)):
        with pytest.raises(TypeError, match="not a bool"):
            traj[key]
        with pytest.raises(TypeError, match="not a bool"):
            traj[0][key]
    # integers, numpy integers and negative indices are unchanged
    assert traj[np.int64(1)][0].radius_of_gyration == traj[1][0].radius_of_gyration
    assert traj[-1][np.int64(0)].radius_of_gyration == traj[traj.n_frames - 1][0].radius_of_gyration


def test_chain_extent_sampler_visits_at_most_max_frames():
    class CountingPositions:
        """A positions array that records which frames were read."""

        def __init__(self, array):
            self.array = array
            self.visited = []

        def __getitem__(self, index):
            self.visited.append(index)
            return self.array[index]

    for n_frames in (1, 9, 10, 11, 19, 20, 21, 39, 100):
        positions = CountingPositions(np.zeros((n_frames, 6, 3), dtype=np.int32))
        traj = SimpleNamespace(n_frames=n_frames, positions=positions, dimensions=(8, 8, 8),
                               hardwall=False, topology=SimpleNamespace(offsets=np.array([0, 3, 6])))
        ps._chain_extent_along(traj, 2, max_frames=10)
        assert 1 <= len(positions.visited) <= 10, (n_frames, positions.visited)
        if n_frames <= 10:
            assert len(positions.visited) == n_frames


def test_batched_rg_is_the_same_whichever_memoised_quantity_is_asked_for_first(traj3d_files):
    """Rg and the gyration eigenvalues hold ONE float64 copy of the whole trajectory.

    They need the centres of mass, and used to get them from centers_of_mass(),
    which made its own float64 copy of the whole trajectory while theirs was
    still alive: two full copies at once, for no difference in the result. The
    copies alive at once are counted; the results must not depend on the order
    in which the memoised quantities are first asked for.
    """
    import weakref

    xtc, pdb, keyfile = traj3d_files
    for quantity in ("radius_of_gyration", "gyration_eigenvalues"):
        store = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)._store
        make_copy = store._whole_k
        alive = [0]
        peak = [0]

        def counted_copy(make_copy=make_copy, alive=alive, peak=peak):
            """``_whole_k`` that counts how many of its copies are alive at once.

            Parameters
            ----------
            make_copy : callable
                The store's own ``_whole_k``.
            alive : list of int
                One-element counter of the copies alive now.
            peak : list of int
                One-element record of the most copies ever alive at once.

            Returns
            -------
            numpy.ndarray
                The float64 copy ``_whole_k`` made.
            """
            copy = make_copy()
            alive[0] += 1
            peak[0] = max(peak[0], alive[0])
            weakref.finalize(copy, lambda: alive.__setitem__(0, alive[0] - 1))
            return copy

        store._whole_k = counted_copy
        getattr(store, quantity)()
        assert peak[0] == 1, (quantity, peak[0])

    rg_first = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    rg = rg_first.radius_of_gyration()
    com_first = lemonade.load(xtc=xtc, pdb=pdb, keyfile=keyfile)
    com = com_first.center_of_mass()
    assert np.array_equal(com_first.radius_of_gyration(), rg)
    assert np.array_equal(rg_first.center_of_mass(), com)
    assert np.array_equal(rg_first._store.gyration_eigenvalues(), com_first._store.gyration_eigenvalues())
