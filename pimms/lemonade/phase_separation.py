## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Phase-separation and droplet-physics analysis for lemonade trajectories.

This module quantifies liquid-liquid phase separation of a PIMMS lattice system:
the coexistence (binodal) densities of the dense and dilute phases, the condensed
fraction and cluster-size order parameters, interfacial width, and droplet shape.

Two complementary geometries are supported:

* **Droplet (spherical)** - a radial density profile about the largest cluster's
  centre of mass, fit to ``rho(r) = 1/2(rho_d + rho_v) - 1/2(rho_d - rho_v)
  tanh((r - R)/w)`` to extract the dense density ``rho_d``, the dilute (vapour)
  density ``rho_v``, the droplet radius ``R`` and the interface width ``w``.
* **Slab** - a 1D density profile along the slab normal (the axis the condensate
  does not span, or the box's long axis when the clusters do not say; the
  geometry of the ``slab_phase_separation`` demo), slabs re-centred per frame and
  fit to a two-interface ``tanh`` to extract the same coexistence quantities.

Densities are volume fractions (occupied lattice sites per available lattice site),
so ``rho`` runs 0..1 and is directly comparable across box sizes.

Typical use::

    from pimms.lemonade import phase_separation as ps
    result = ps.analyze(traj)          # everything, auto-detecting the geometry
    r, rho, sites = ps.radial_density_profile_with_site_counts(traj)
    fit = ps.fit_radial_profile(r, rho, site_counts=sites)   # drops thin shells
"""

import warnings
from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# axis handling
# ---------------------------------------------------------------------------

#: A hardwall frame counts as holding a slab only when its dense block leaves at
#: least this many planes that are not dense on one side. A one-phase solution
#: between two walls is dense everywhere except in the depletion layer at each
#: wall, which is one to three planes of ramp, not a dilute phase.
_MIN_DILUTE_PLANES = 4

#: Fraction of its own asymptotic gap that a fitted profile must actually traverse
#: inside the data for the fit to count as a genuine two-phase description.
_MIN_REALISED_FRACTION = 0.5

#: Absolute floor on the coexistence gap, in occupied-site fraction.
_MIN_DENSITY_GAP = 1e-3

#: The gap must also stand clear of the scatter in the profile: a tanh will happily
#: fit a hump in the noise of a flat profile and report the hump as a coexistence gap.
_MIN_GAP_OVER_NOISE = 3.0

#: A radial shell with fewer lattice sites than this is too noisy to constrain a
#: fit: the ``[0, 1)`` shell is ONE site, always at the (rounded) centre of mass of
#: the largest cluster, so it reads ~0.6-0.9 even for a homogeneous solution. Fitted
#: with equal weight it turned that single site into a "droplet" of rho_dense = 1.
_MIN_SHELL_SITES = 20

#: The same floor for 2D profiles, whose shells hold ~2 pi r sites (see _min_shell_sites)
_MIN_SHELL_SITES_2D = 8

#: A fitted droplet radius below this (lattice units) is a handful of sites, not a
#: dense phase - the signature of the fit latching onto the central shells.
_MIN_DROPLET_RADIUS = 2.0

#: Fewest planes a slab may be thick, counted as the planes on which the fitted
#: profile is above the midpoint of the two densities. With one plane the dense
#: density and the slab thickness are a single number seen through two
#: parameters (a thinner, denser slab fits as well), so two planes are the least
#: that measures the density rather than extrapolating it. Counting planes
#: rather than comparing the fitted thickness with a length makes the answer the
#: same for a film and its mirror image: the fitted thickness of a sharp
#: two-plane film can land anywhere between 1.5 and 2.5.
_MIN_SLAB_PLANES = 2

#: The fitted slab profile must rise at least this fraction of the way from the
#: dilute density to the dense one somewhere in the data. For a free slab the
#: profile at its centre reaches tanh(half_width / interface_width) of the gap,
#: so 0.75 is the statement "the slab is thicker than its two interfaces"
#: (half_width > interface_width, tanh(1) = 0.76). Below that the slab has no
#: interior and its dense density is an extrapolation.
_MIN_DENSE_REACH = 0.75

#: A plane counts as part of the dilute plateau when the fitted slab profile is
#: within this fraction of the coexistence gap of the dilute density there.
_PLATEAU_TOLERANCE = 0.1

#: Fewest planes the dilute plateau may hold. Coexistence needs a dilute phase,
#: and a phase has an interior: a plane whose two neighbours are in the same
#: phase. The depletion layer of a one-phase solution at a hard wall is a ramp
#: of a few planes with no such plane, and a tanh fitted to it places its dilute
#: asymptote below every density actually observed.
_MIN_DILUTE_PLATEAU_PLANES = 3

#: The coexistence gap of a slab must be at least this many standard deviations
#: of the single-frame difference between one dense plane and one dilute plane:
#: a gap that is real is there in every single frame, while the hump that
#: re-centring builds out of the fluctuations of a one-phase solution is about
#: as tall as those fluctuations.
_MIN_GAP_OVER_FRAME_SCATTER = 3.0

#: The dense region of a slab must hold at least this many standard deviations
#: more chains than a uniform solution would put there (see
#: :func:`_chain_excess_sigma`).
_MIN_CHAIN_EXCESS_SIGMA = 3.0


def _normalise_axis(axis, n_dim, caller):
    """Turn a user-supplied axis into an index in ``0 .. n_dim - 1``.

    Negative values count from the last axis of the SYSTEM, as they do for a
    numpy array with ``n_dim`` axes. They used to be passed straight through:
    ``dims[-1]`` is the right box length, but the remaining axes were then
    selected with ``i != axis``, which is true for every ``i`` when ``axis`` is
    negative, so the cross-section came out as the whole box volume (densities
    too small by the box length), and in 2D ``positions[:, -1]`` is the padded
    zero z column rather than y.

    Parameters
    ----------
    axis : int
        The axis as the caller passed it.
    n_dim : int
        Dimensionality of the system, 2 or 3.
    caller : str
        Name of the public function, for the error message.

    Returns
    -------
    int
        The same axis as a non-negative index.

    Raises
    ------
    ValueError
        If ``axis`` is not an integer, or lies outside ``-n_dim .. n_dim - 1``.
    """
    if isinstance(axis, (bool, np.bool_)) or not isinstance(axis, (int, np.integer)):
        raise ValueError("%s: axis must be an integer in %d..%d, got %r"
                         % (caller, -n_dim, n_dim - 1, axis))
    if not -n_dim <= int(axis) < n_dim:
        raise ValueError("%s: axis %d is out of range for a %dD system (use %d..%d)"
                         % (caller, int(axis), n_dim, -n_dim, n_dim - 1))
    return int(axis) % n_dim


def _vote_slab_normal(spans, dims):
    """Pick the slab normal from the axes each frame's largest cluster spans.

    A slab spans every axis but one, and the one it does not span is its
    normal. Each frame whose largest cluster spans exactly ``n_dim - 1`` axes
    votes for the remaining one; an axis that collects the votes of more than
    half of the frames holding a cluster wins. Otherwise the clusters do not
    say (a droplet, a network, a one-phase solution) and the longest box axis
    is returned, the first of them on a tie, which is what was always used
    before.

    Parameters
    ----------
    spans : list of (tuple of int or None)
        Per frame, the axes the largest cluster spans, or ``None`` for a frame
        with no qualifying cluster.
    dims : sequence of int
        Box extent in lattice units, one entry per axis.

    Returns
    -------
    int
        Index of the slab normal.
    """
    n_dim = len(dims)
    votes = np.zeros(n_dim, dtype=np.int64)
    n_with_cluster = 0
    for axes in spans:
        if axes is None:
            continue
        n_with_cluster += 1
        if len(axes) == n_dim - 1:
            votes[[d for d in range(n_dim) if d not in axes][0]] += 1
    if n_with_cluster and 2 * int(votes.max()) > n_with_cluster:
        return int(np.argmax(votes))
    return int(np.argmax(dims))


# ---------------------------------------------------------------------------
# cluster / order-parameter helpers
# ---------------------------------------------------------------------------

def _largest_clusters(traj, min_beads=1):
    """Yield ``(frame_index, clusters)`` per frame, largest (most beads) first.

    :attr:`Frame.clusters <pimms.lemonade.Frame>` is ordered by bead count, so ``clusters[0]``
    here is the condensate. Frames with no qualifying cluster yield an empty list.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to walk, one frame at a time.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``1``, i.e. keep
        every cluster).

    Yields
    ------
    f : int
        The frame index.
    clusters : list of pimms.lemonade.Cluster
        That frame's qualifying clusters, largest (most beads) first.
    """
    for f in range(traj.n_frames):
        clusters = [c for c in traj[f].clusters if c.n_beads >= min_beads]
        yield f, clusters


def condensed_fraction(traj, min_beads=1):
    """Per-frame fraction of all beads that sit in the single largest cluster.

    This is the basic phase-separation order parameter: ~0 in a well-mixed
    phase, ->1 when most material is in one condensate.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``1``).

    Returns
    -------
    numpy.ndarray
        ``(n_frames,)`` float64 condensed fraction, ``0`` in frames with no
        qualifying cluster.
    """
    total = traj.n_beads
    out = np.zeros(traj.n_frames)
    for f, clusters in _largest_clusters(traj, min_beads):
        if clusters:
            out[f] = clusters[0].n_beads / total
    return out


def largest_cluster_size(traj, by="beads", min_beads=1):
    """Per-frame size of the largest cluster.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    by : str, optional
        ``'beads'`` (default) counts beads, ``'chains'`` counts chains.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``1``).

    Returns
    -------
    numpy.ndarray
        ``(n_frames,)`` int64 size of the largest cluster, ``0`` in frames with
        no qualifying cluster.

    Raises
    ------
    ValueError
        If ``by`` is neither ``'beads'`` nor ``'chains'``.
    """
    if by not in ("beads", "chains"):
        raise ValueError("largest_cluster_size: by must be 'beads' or 'chains', got %r" % (by,))
    out = np.zeros(traj.n_frames, dtype=np.int64)
    for f, clusters in _largest_clusters(traj, min_beads):
        if clusters:
            out[f] = clusters[0].n_beads if by == "beads" else clusters[0].n_chains
    return out


def number_of_clusters(traj, min_beads=2):
    """Per-frame count of clusters with at least ``min_beads`` beads.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    min_beads : int, optional
        Smallest cluster that counts, in beads (default ``2``, which drops
        single-bead clusters).

    Returns
    -------
    numpy.ndarray
        ``(n_frames,)`` int64 number of qualifying clusters per frame.
    """
    out = np.zeros(traj.n_frames, dtype=np.int64)
    for f, clusters in _largest_clusters(traj, min_beads):
        out[f] = len(clusters)
    return out


def cluster_size_distribution(traj, by="beads", min_beads=1):
    """All cluster sizes pooled across frames, as one flat array (for histograms).

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    by : str, optional
        ``'beads'`` (default) measures each cluster in beads, ``'chains'`` in
        chains.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``1``).

    Returns
    -------
    numpy.ndarray
        1D int64 array of every qualifying cluster's size, pooled over all
        frames (so it has no per-frame structure).

    Raises
    ------
    ValueError
        If ``by`` is neither ``'beads'`` nor ``'chains'``.
    """
    if by not in ("beads", "chains"):
        raise ValueError("cluster_size_distribution: by must be 'beads' or 'chains', got %r" % (by,))
    sizes = []
    for _f, clusters in _largest_clusters(traj, min_beads):
        sizes.extend(c.n_beads if by == "beads" else c.n_chains for c in clusters)
    return np.asarray(sizes, dtype=np.int64)


def _cluster_spans_box(cluster, all_axes):
    """Does this cluster span the box on any (or every) axis?

    The per-axis answer is :meth:`Cluster.spanning_axes
    <pimms.lemonade.Cluster.spanning_axes>`, computed once per cluster on its
    cached single image. Under periodic boundaries an axis counts when the
    cluster is connected to its own image through that face, which is PIMMS's
    own test, :func:`pimms.cluster_utils.percolating_axes`: a connected cluster
    that does not wind through the boundary always fits inside one image, so a
    per-axis extent reaching the box length is necessary for winding, but it is
    not sufficient - a contact staircase from one corner of the box to the other
    reaches the box length on an axis without any pair of beads meeting through
    that face, and has a perfectly good single image - so a candidate axis is
    confirmed by an actual pair of beads that touch through the face. This
    function used to apply the extent test alone, so it counted such clusters
    as spanning while PIMMS's gather (correctly) did not. Under HARDWALL
    nothing connects through a wall and an axis counts when the cluster
    touches both of its walls.

    Parameters
    ----------
    cluster : pimms.lemonade.Cluster
        The cluster to test.
    all_axes : bool
        ``True`` requires the cluster to span every axis; ``False`` accepts any
        one axis.

    Returns
    -------
    bool
        Whether the cluster spans the box under the chosen criterion.
    """
    # The first gather of a periodic cluster warns when it winds the box,
    # because a shape computed from such a gathering is meaningless. This
    # function IS the spanning detector: its callers (spanning_fraction,
    # analyze, the radial profile, droplet_shape, the surface-tension
    # estimators) use the answer to guard exactly those shape quantities, so
    # being warned that the thing it was asked to detect has been detected is
    # noise, and that one warning is silenced for this call only. The axes are
    # cached on the Cluster alongside the gathered image, so a later shape call
    # on the same object neither re-gathers nor re-tests.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="single-image gather: cluster percolates")
        axes = cluster.spanning_axes()
    spans = [d in axes for d in range(cluster._store.n_dim)]
    return all(spans) if all_axes else any(spans)


def _largest_cluster_spanning_axes(traj, min_beads=2):
    """Per frame, the axes the largest cluster spans.

    One gather per frame answers every spanning question :func:`analyze` asks
    (any axis, every axis, which axis is the slab normal), so it is done once
    here rather than once per question.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``2``).

    Returns
    -------
    list of (tuple of int or None)
        One entry per frame: the spanning axes of the largest cluster in axis
        order (empty for a compact cluster), or ``None`` when the frame holds
        no qualifying cluster.
    """
    spans = []
    for _f, clusters in _largest_clusters(traj, min_beads):
        if not clusters:
            spans.append(None)
            continue
        # same reasoning as in _cluster_spans_box: this IS the spanning
        # detector, so the gather's own percolation warning is noise here
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="single-image gather: cluster percolates")
            spans.append(tuple(clusters[0].spanning_axes()))
    return spans


def slab_normal(traj, min_beads=2):
    """The axis a slab analysis should profile along.

    A slab spans every box axis but one - it is connected to its own periodic
    image through those faces, or touches both of their walls under a hardwall
    - and the axis it does not span is its normal. When the largest cluster
    spans exactly ``n_dim - 1`` axes, leaving the same one free, in more than
    half of the frames that hold a cluster, that axis is returned. Otherwise
    the clusters do not single an axis out (a droplet, a network, a one-phase
    solution) and the longest box axis is returned, the first of them on a tie.

    The longest axis alone used to be used. A slab lying across a short axis,
    or in a box with two equal long axes, was then profiled along one of its
    own in-plane directions, where it is flat: a 10-plane slab normal to y in a
    ``(40, 40, 10)`` box came back as "not phase separated".

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    min_beads : int, optional
        Ignore clusters smaller than this many beads when picking the
        condensate (default ``2``).

    Returns
    -------
    int
        Index of the slab normal, ``0 .. n_dim - 1``.
    """
    return _vote_slab_normal(_largest_cluster_spanning_axes(traj, min_beads),
                             tuple(traj.dimensions)[:traj.n_dim])


def spanning_fraction(traj, min_beads=2, all_axes=False):
    """Fraction of frames in which the largest cluster spans the box.

    ``all_axes=False`` (default) counts a frame when the largest cluster spans
    **any** axis. Under periodic boundaries that means it is connected to its own
    periodic image through that face - a pair of its beads touch through it;
    merely reaching the box length is not enough - and so has no well-defined
    single image, and no droplet quantity (radial profile, hull, sphericity)
    computed from it means anything. Under a hardwall it means the cluster
    touches both walls of that axis: a wall-bounded film or network, not a
    droplet, although its image is exact. ``all_axes=True`` counts only
    frames where it spans **every** axis: a space-filling network, which is what the
    contact clustering of a *homogeneous* solution looks like at moderate volume
    fraction. A slab spans two axes but not the third, so it is caught by the first
    criterion and not the second.

    Frames with no qualifying cluster count as not spanning.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``2``).
    all_axes : bool, optional
        ``True`` counts only frames whose largest cluster spans every axis (a
        space-filling network); ``False`` (default) counts any axis.

    Returns
    -------
    float
        Fraction of frames that count as spanning, ``0.0`` for an empty
        trajectory.
    """
    n_span = 0
    for _f, clusters in _largest_clusters(traj, min_beads):
        if clusters and _cluster_spans_box(clusters[0], all_axes):
            n_span += 1
    return n_span / traj.n_frames if traj.n_frames else 0.0


# ---------------------------------------------------------------------------
# density profiles
# ---------------------------------------------------------------------------

def _min_image(delta, dims):
    """Wrap displacement(s) into ``[-L/2, L/2)`` per axis (broadcasts over dims).

    Parameters
    ----------
    delta : numpy.ndarray
        Float64 displacements, with the axis (dimension) index last.
    dims : numpy.ndarray
        Float64 box length per axis, broadcastable against ``delta``.

    Returns
    -------
    numpy.ndarray
        Float64 minimum-image displacements, same shape as ``delta``.
    """
    return delta - dims * np.round(delta / dims)


def _shell_site_counts(dimensions, edges):
    """Number of lattice sites whose PBC (min-image) distance from a point falls in
    each radial shell. Translationally invariant under PBC, so computed once.

    Built from separable per-axis offsets and accumulated a slab at a time. Listing
    every site explicitly (``np.indices(dims).reshape(nd, -1).T`` as float64) allocates
    ``n_dim`` float64 values per lattice site up front - about 200 MB for a 200^3 box,
    plus the distance array - which scales with box volume and has nothing to do with
    how many beads are actually being analysed. Histogram counts are additive, so
    accumulating over slabs gives exactly the same result with memory bounded by one
    slab.

    Parameters
    ----------
    dimensions : sequence of int
        Box extent in lattice units, 2 or 3 entries.
    edges : numpy.ndarray
        1D float64 radial bin edges, in lattice units and ascending.

    Returns
    -------
    numpy.ndarray
        ``(len(edges) - 1,)`` float64 number of lattice sites in each shell.
    """
    dims = np.asarray(dimensions, dtype=np.float64)

    # per-axis minimum-image offset of every coordinate from the origin
    offsets = [_min_image(np.arange(int(n), dtype=np.float64), d)
               for n, d in zip(dimensions, dims)]
    squared = [off ** 2 for off in offsets]

    counts = np.zeros(len(edges) - 1, dtype=np.int64)

    # sum the per-axis squared offsets by broadcasting, one x-slab at a time. The
    # addition order matches the axis-order reduction of the original.
    if len(dimensions) == 2:
        for x2 in squared[0]:
            r = np.sqrt(x2 + squared[1])
            counts += np.histogram(r, bins=edges)[0]
    else:
        plane = squared[1][:, np.newaxis] + squared[2][np.newaxis, :]
        for x2 in squared[0]:
            r = np.sqrt(x2 + plane)
            counts += np.histogram(r, bins=edges)[0]

    return counts.astype(np.float64)


def radial_density_profile(traj, bin_width=1.0, r_max=None, min_beads=2):
    """Spherically averaged density profile about the largest cluster's COM.

    For every frame the minimum-image distance of *every* bead from the condensate
    centre of mass is binned into radial shells and normalised by the number of
    lattice sites in each shell, giving a volume-fraction profile that falls from
    the dense core to the dilute background. Profiles are averaged over frames.

    A shell that contains **no lattice site** in any frame (under HARDWALL the
    outer shells are progressively cut off by the walls, and beyond the farthest
    in-box site there is nothing to be occupied) is returned as ``nan``, not
    ``0``: the two are different facts. Reporting such shells as zero density
    used to pad the tail of every hardwall profile with fake vacuum and pull the
    fitted dilute-phase density well below the true value.

    A frame whose largest cluster spans the box is left out: such a cluster is a
    network or a slab, not a droplet, so a profile about its centre is not a
    droplet profile. Under periodic boundaries it is connected to its own image
    and has no single image at all, so that centre is an arbitrary point of the
    window the gather hands back; under a hardwall it touches both walls, its
    centre of mass is exact and the profile about it is well defined, but it
    describes a wall-bounded condensate. The usual one is frame 0 of a run that
    saved its equilibration. A warning says how many frames were left out. If
    every frame spans there is no droplet to profile; the profile about that
    centre is then returned, with a different warning, because its percentiles
    still estimate the dense and dilute densities.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    bin_width : float, optional
        Radial shell width in lattice units (default ``1.0``).
    r_max : float, optional
        Outermost radius in lattice units. Default ``None``: half the shortest
        box axis under PBC (beyond that the minimum image folds back), or the
        full box diagonal under a hardwall.
    min_beads : int, optional
        Ignore clusters smaller than this many beads when picking the
        condensate (default ``2``).

    Returns
    -------
    radii : numpy.ndarray
        1D float64 shell-centre radii in lattice units.
    density : numpy.ndarray
        1D float64 frame-averaged occupied fraction ``rho(r)``, ``nan`` where
        the shell holds no lattice site.
    """
    centers, density, _sites, _scatter = _radial_profile_with_site_counts(
        traj, bin_width=bin_width, r_max=r_max, min_beads=min_beads)
    return centers, density


def radial_density_profile_with_site_counts(traj, bin_width=1.0, r_max=None, min_beads=2):
    """:func:`radial_density_profile` plus the mean number of lattice sites per shell.

    The site count is what a fit should weight by: the innermost shell ``[0, 1)``
    holds exactly one site, so its "density" is a single 0/1 draw per frame, while
    the 3D shell ``[10, 11)`` averages over about 1360 sites.
    :func:`fit_radial_profile` uses the counts to drop shells below a floor.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    bin_width : float, optional
        Radial shell width in lattice units (default ``1.0``).
    r_max : float, optional
        Outermost radius in lattice units (default ``None``, chosen from the box
        as in :func:`radial_density_profile`).
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``2``).

    Returns
    -------
    radii : numpy.ndarray
        1D float64 shell-centre radii in lattice units.
    density : numpy.ndarray
        1D float64 frame-averaged occupied fraction, ``nan`` where the shell
        holds no lattice site.
    site_counts : numpy.ndarray
        1D float64 mean number of lattice sites per shell over the frames that
        contributed, ``0`` for shells no frame could fill.
    """
    centers, density, site_counts, _scatter = _radial_profile_with_site_counts(
        traj, bin_width=bin_width, r_max=r_max, min_beads=min_beads)
    return centers, density, site_counts


def radial_density_profile_with_scatter(traj, bin_width=1.0, r_max=None, min_beads=2):
    """:func:`radial_density_profile_with_site_counts` plus the frame-to-frame
    scatter of each shell.

    The profile is centred, every frame, on the largest cluster. In a one-phase
    solution the largest cluster is a chain or a few touching chains, and a
    chain is denser than its surroundings, so the averaged profile of any
    solution shows a hump at the origin whether or not there is a droplet. What
    centring cannot do is make that hump steadier than the fluctuation it is: a
    droplet's core reads the same density in every frame, a coil's does not.
    :func:`fit_radial_profile` applies that test when it is handed the scatter,
    and :func:`analyze` always hands it over.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    bin_width : float, optional
        Radial shell width in lattice units (default ``1.0``).
    r_max : float, optional
        Outermost radius in lattice units (default ``None``, chosen from the box
        as in :func:`radial_density_profile`).
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``2``).

    Returns
    -------
    radii : numpy.ndarray
        1D float64 shell-centre radii in lattice units.
    density : numpy.ndarray
        1D float64 frame-averaged occupied fraction, ``nan`` where the shell
        holds no lattice site.
    site_counts : numpy.ndarray
        1D float64 mean number of lattice sites per shell over the frames that
        contributed, ``0`` for shells no frame could fill.
    scatter : numpy.ndarray
        1D float64 standard deviation over the contributing frames of each
        shell's occupied fraction; ``nan`` where fewer than two frames
        contributed (so all ``nan`` for a one-frame trajectory).
    """
    return _radial_profile_with_site_counts(
        traj, bin_width=bin_width, r_max=r_max, min_beads=min_beads)


def _radial_profile_with_site_counts(traj, bin_width=1.0, r_max=None, min_beads=2):
    """Bin every bead's distance from the condensate COM into radial shells.

    The single implementation behind :func:`radial_density_profile`,
    :func:`radial_density_profile_with_site_counts` and
    :func:`radial_density_profile_with_scatter`. Bead distances and shell
    site counts use the same metric (minimum image under PBC, plain Cartesian
    under a hardwall), so occupied never exceeds available and the density is a
    true occupied fraction in ``[0, 1]``.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    bin_width : float, optional
        Radial shell width in lattice units (default ``1.0``).
    r_max : float, optional
        Outermost radius in lattice units. Default ``None``: the full box
        diagonal under a hardwall, else half the shortest box axis.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``2``); frames
        with no qualifying cluster contribute nothing.

    Returns
    -------
    centers : numpy.ndarray
        1D float64 shell-centre radii in lattice units.
    density : numpy.ndarray
        1D float64 frame-averaged occupied fraction, ``nan`` in shells no frame
        could measure.
    site_counts : numpy.ndarray
        1D float64 mean number of lattice sites per shell over the contributing
        frames, ``0`` where there were none.
    scatter : numpy.ndarray
        1D float64 standard deviation over the contributing frames of each
        shell's occupied fraction, ``nan`` where fewer than two contributed.
    """
    dims = np.asarray(traj.dimensions, dtype=np.float64)
    nd = traj.n_dim
    hardwall = bool(getattr(traj, "hardwall", False))
    if r_max is None:
        if hardwall:
            # plain Cartesian distances can reach the full box diagonal
            r_max = float(np.sqrt((dims[:nd] ** 2).sum()))
        else:
            # beyond half the smallest box length the minimum image folds back
            r_max = float(min(traj.dimensions)) / 2.0
    edges = np.arange(0.0, r_max + bin_width, bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    if not hardwall:
        # PBC: shell counts are translation-invariant, compute them once
        site_counts = _shell_site_counts(traj.dimensions, edges)
        safe = site_counts.copy()
        safe[safe == 0] = np.nan
        _hw_safe_cache = None
    else:
        # HARDWALL: the box is not periodic, so distances are plain Cartesian
        # and the number of in-box sites per shell depends on where the COM
        # sits (shells get truncated by the walls). Compute per-COM shell
        # counts, cached by the rounded COM. Previously the periodic metric
        # was applied here too, silently folding beads more than half a box
        # from the condensate into inner shells.
        _hw_safe_cache = {}

    def _hardwall_safe(com):
        """Per-shell in-box site counts for a hardwall box, cached by rounded COM.

        Parameters
        ----------
        com : numpy.ndarray
            ``(n_dim,)`` float64 integer-valued centre of mass the shells are
            measured from.

        Returns
        -------
        numpy.ndarray
            ``(len(edges) - 1,)`` float64 site counts, with empty shells set to
            ``nan`` so dividing by them yields ``nan`` rather than a fake zero
            density.
        """
        key = tuple(int(v) for v in com)
        if key not in _hw_safe_cache:
            axes = [np.arange(int(n), dtype=np.float64) - c
                    for n, c in zip(traj.dimensions, com)]
            sq = [a ** 2 for a in axes]
            counts = np.zeros(len(edges) - 1, dtype=np.int64)
            if nd == 2:
                for x2 in sq[0]:
                    r = np.sqrt(x2 + sq[1])
                    counts += np.histogram(r, bins=edges)[0]
            else:
                plane = sq[1][:, np.newaxis] + sq[2][np.newaxis, :]
                for x2 in sq[0]:
                    r = np.sqrt(x2 + plane)
                    counts += np.histogram(r, bins=edges)[0]
            safe_local = counts.astype(np.float64)
            safe_local[safe_local == 0] = np.nan
            _hw_safe_cache[key] = safe_local
        return _hw_safe_cache[key]

    # Two sets of accumulators: frames whose largest cluster is a droplet, and
    # frames where it spans the box. A spanning cluster is a network or a slab,
    # not a droplet: under periodic boundaries it is connected to its own image
    # and the COM of the window the gather hands back is an arbitrary point;
    # under a hardwall it touches both walls and its COM is exact, but the
    # profile about it is that of a wall-bounded condensate. Either way it is
    # not a droplet profile. When the trajectory holds any droplet frame the
    # spanning frames are left out - frame 0 of a run that saved its
    # equilibration is the usual one, and it used to be averaged in. When EVERY
    # frame spans there is no droplet profile to give, and the profile about
    # that centre is returned instead with a warning: its percentiles still
    # estimate the two densities (which is how analyze() reads a network or a
    # slab it was asked to treat as a droplet), even though its shape describes
    # no droplet.
    acc = np.zeros(len(centers))
    n_valid = np.zeros(len(centers))
    sites_acc = np.zeros(len(centers))
    rows = []                              # per-frame profiles, for the scatter
    acc_span = np.zeros(len(centers))
    n_valid_span = np.zeros(len(centers))
    sites_acc_span = np.zeros(len(centers))
    rows_span = []
    n_droplet = 0
    n_spanning = 0
    positions = traj.positions
    for f, clusters in _largest_clusters(traj, min_beads):
        if not clusters:
            continue
        spanning = _cluster_spans_box(clusters[0], False)
        # bin bead distances from the INTEGER COM using the same metric as the
        # shell site counts, so occupied <= available in every shell and the
        # density is a true occupied fraction in [0, 1].
        # floor(x + 0.5), not np.round: numpy rounds halves to the nearest EVEN
        # integer, so a centre of mass at 11.5 went to 12 and its mirror image
        # at 12.5 also to 12, and the profile changed under a translation or a
        # reflection of the whole system. Halves now always go up. (Exact
        # mirror invariance is not attainable on a lattice: a centre halfway
        # between two sites has to be put on one of them.)
        com = np.mod(np.floor(np.asarray(clusters[0].center_of_mass, dtype=np.float64) + 0.5),
                     dims[:nd])
        delta = positions[f][:, :nd].astype(np.float64) - com
        if hardwall:
            d = delta                      # plain Cartesian
            safe_f = _hardwall_safe(com)
        else:
            d = _min_image(delta, dims[:nd])
            safe_f = safe
        r = np.sqrt((d * d).sum(axis=1))
        counts, _ = np.histogram(r, bins=edges)
        contribution = counts / safe_f        # nan where this frame's shell has no site
        finite = np.isfinite(contribution)
        if spanning:
            n_spanning += 1
            acc_span[finite] += contribution[finite]
            sites_acc_span[finite] += safe_f[finite]
            rows_span.append(contribution)
            n_valid_span[finite] += 1
        else:
            n_droplet += 1
            acc[finite] += contribution[finite]
            sites_acc[finite] += safe_f[finite]
            rows.append(contribution)
            n_valid[finite] += 1

    if n_spanning and n_droplet == 0:
        acc, sites_acc, n_valid = acc_span, sites_acc_span, n_valid_span
        rows = rows_span
        warnings.warn(
            "radial_density_profile: the largest cluster spans the box in every one of the "
            "%d frames with a cluster - a network or a slab, not a droplet - so the profile "
            "is centred on %s; its percentiles still estimate the two densities, its shape "
            "does not describe a droplet. See phase_separation.spanning_fraction."
            % (n_spanning, "the wall-bounded condensate" if hardwall
               else "an arbitrary point of its gathered image"), stacklevel=3)
    elif n_spanning:
        warnings.warn(
            "radial_density_profile: %d of %d frames with a cluster were left out because "
            "the largest cluster spans the box (%s) - a network or a slab, not a droplet, "
            "so a profile about its centre is not a droplet profile; see "
            "phase_separation.spanning_fraction."
            % (n_spanning, n_spanning + n_droplet, "touches both walls" if hardwall
               else "connected to its own periodic image"), stacklevel=3)

    with np.errstate(invalid="ignore", divide="ignore"):
        density = np.where(n_valid > 0, acc / np.maximum(n_valid, 1), np.nan)
        site_counts = np.where(n_valid > 0, sites_acc / np.maximum(n_valid, 1), 0.0)
        # the scatter is taken about the mean frame by frame (not as a difference
        # of two accumulated sums, which leaves 1e-9 where the answer is zero)
        per_frame = np.asarray(rows, dtype=np.float64).reshape(len(rows), len(centers))
        measured = np.isfinite(per_frame)
        mean = acc / np.maximum(n_valid, 1)
        squares = np.where(measured, (per_frame - mean) ** 2, 0.0).sum(axis=0)
        # sample standard deviation (n - 1): with a handful of frames the
        # population form reads low and made the single-frame test too easy
        scatter = np.where(n_valid > 1, np.sqrt(squares / np.maximum(n_valid - 1, 1)), np.nan)
    return centers, density, site_counts, scatter


def _hardwall_slab_frame_kind(counts):
    """Classify one frame's bead-count profile along the normal of a hardwall box.

    The dense bins are those holding at least half the frame's peak count, the
    same threshold :func:`slab_density_profile` aligns on.

    Parameters
    ----------
    counts : numpy.ndarray
        ``(L,)`` float64 number of beads in each plane along the slab normal.

    Returns
    -------
    kind : str
        ``'wall'`` if a dense bin is one of the two end planes (the condensate
        touches a wall, or there are no beads at all), ``'free'`` otherwise.
    clear : bool
        Whether the dense bins form one slab: a single block, allowing one-plane
        dips inside it, that does not reach from one wall to the other and
        leaves a dilute region at least ``_MIN_DILUTE_PLANES`` (4) planes wide
        on one side. A homogeneous frame does not pass: at low density its
        dense bins are scattered over the whole box, and at higher density they
        fill it apart from the depleted plane or two at each wall. The second
        case used to pass, because a block that misses only the two wall planes
        does not "reach from one wall to the other"; a one-phase solution was
        then treated as a slab trajectory. This is what tells a slab trajectory
        from a one-phase one.
    """
    length = counts.size
    peak = float(counts.max()) if length else 0.0
    if peak <= 0:
        return "wall", False
    idx = np.nonzero(counts >= 0.5 * peak)[0]
    kind = "wall" if (idx[0] == 0 or idx[-1] == length - 1) else "free"
    one_block = bool(np.all(np.diff(idx) <= 2)) and not (idx[0] == 0 and idx[-1] == length - 1)
    widest_margin = max(int(idx[0]), int(length - 1 - idx[-1]))
    clear = one_block and widest_margin >= _MIN_DILUTE_PLANES
    return kind, clear


def _slab_aligned_counts(traj, axis, stacklevel=3):
    """Per-frame bead counts along the slab normal, aligned on the dense slab.

    The single implementation behind :func:`slab_density_profile` and
    :func:`slab_density_profile_with_scatter`; see the first for what the
    alignment does and why.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    axis : int
        Index of the slab normal, already normalised to ``0 .. n_dim - 1``.
    stacklevel : int, optional
        ``stacklevel`` of the warning raised when hardwall frames are left out,
        so that it points at the caller of the public function (default ``3``).

    Returns
    -------
    aligned : numpy.ndarray
        ``(n_kept, L)`` float64 bead count per plane for every frame that is
        averaged, after alignment; ``nan`` in the planes a translated hardwall
        frame vacated.
    cross_section : int
        Number of lattice sites in one plane.
    """
    dims = traj.dimensions
    hardwall = bool(getattr(traj, "hardwall", False))
    length = dims[axis]
    cross_section = int(np.prod([d for i, d in enumerate(dims) if i != axis]))
    coord = np.arange(length)
    angle = 2.0 * np.pi * coord / length

    positions = traj.positions
    frame_counts = [np.bincount(positions[f][:, axis], minlength=length).astype(np.float64)
                    for f in range(traj.n_frames)]
    if not hardwall:
        aligned = np.zeros((traj.n_frames, length))
        for f, counts in enumerate(frame_counts):
            # circular centre of mass of the 1D density -> shift dense region to L/2
            cx = (counts * np.cos(angle)).sum()
            cy = (counts * np.sin(angle)).sum()
            com = (np.arctan2(-cy, -cx) + np.pi) / (2.0 * np.pi) * length
            shift = int(round(length / 2.0 - com))
            aligned[f] = np.roll(counts, shift)
        return aligned, cross_section

    # The box is not periodic: a slab cannot wrap, and rolling the profile would
    # move a condensate wetting a wall into the middle of the box, turning its
    # flat wall face into a fictitious second interface. A slab that does NOT
    # touch a wall still diffuses between the walls, so it is aligned by a plain
    # translation of its dense centroid to the window centre (frame-averaging the
    # raw counts smeared a wandering slab into a broad hump with no plateau).
    kinds, clear = [], []
    for counts in frame_counts:
        kind, is_slab = _hardwall_slab_frame_kind(counts)
        kinds.append(kind)
        clear.append(is_slab)
    n_wall = kinds.count("wall")
    n_free = kinds.count("free")
    # Only a trajectory that actually holds a slab (in at least half of its
    # frames) is aligned or split. In a one-phase system the dense bins are
    # noise: translating each frame by the centroid of its noise smears the
    # depletion layer at the walls into the bulk (the wall plane of an athermal
    # solution read 0.041 where the plain average is 0.010), and which of the
    # two ends the dense bins happen to reach is random, so splitting on it
    # throws frames away for nothing. Such a trajectory is averaged as it is.
    holds_slab = traj.n_frames > 0 and 2 * sum(clear) >= traj.n_frames
    keep_kind = None
    # Each frame was handled on its own, so a slab that wandered onto a wall in
    # one frame in six came back as a fake wetting film at the wall plus a slab
    # diluted by a sixth, and the fit still reported success.
    if holds_slab and n_wall and n_free:
        keep_kind = "free" if n_free >= n_wall else "wall"
        n_left = n_wall if keep_kind == "free" else n_free
        warnings.warn(
            "slab_density_profile: the condensate touches a wall of the hardwall box in "
            "%d frames and is free of both walls in %d; a wall-wetting film and a free "
            "slab are different profiles and cannot be averaged into one, so the %d "
            "%s frames were left out and the profile describes the %s."
            % (n_wall, n_free, n_left,
               "wall-touching" if keep_kind == "free" else "free-slab",
               "free slab" if keep_kind == "free" else "wall-wetting condensate"),
            stacklevel=stacklevel)

    # Vacated planes are not padded: the old padding with the frame's MEDIAN
    # dilute count sat well below the mean for Poisson-like dilute counts (0 for
    # half a bead per plane), and pulled the dilute density down by a quarter.
    # Each plane is instead averaged over the frames that actually cover it.
    rows = []
    for counts, kind in zip(frame_counts, kinds):
        if keep_kind is not None and kind != keep_kind:
            continue
        shift = 0
        if holds_slab and kind == "free":
            dense = counts >= 0.5 * counts.max()
            idx = np.nonzero(dense)[0]
            centroid = float((idx * counts[idx]).sum() / counts[idx].sum())
            shift = int(round(length / 2.0 - centroid))
        row = np.full(length, np.nan)
        if shift >= 0:
            row[shift:] = counts[:length - shift]
        else:
            row[:length + shift] = counts[-shift:]
        rows.append(row)
    aligned = np.asarray(rows, dtype=np.float64).reshape(len(rows), length)
    return aligned, cross_section


def _slab_profile_with_scatter(traj, axis, caller, stacklevel=4):
    """Frame-averaged slab profile and its frame-to-frame scatter.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    axis : int or None
        Index of the slab normal as the caller passed it; ``None`` picks it
        with :func:`slab_normal`.
    caller : str
        Name of the public function, for the axis error message.
    stacklevel : int, optional
        ``stacklevel`` handed to the hardwall "frames left out" warning
        (default ``4``: the caller of the public function).

    Returns
    -------
    coordinate : numpy.ndarray
        ``(L,)`` float64 lattice coordinate along the axis.
    density : numpy.ndarray
        ``(L,)`` float64 frame-averaged occupied fraction, ``nan`` where no
        kept frame covers the plane.
    scatter : numpy.ndarray
        ``(L,)`` float64 standard deviation over frames of that plane's
        occupied fraction, ``nan`` where fewer than two frames cover it.

    Raises
    ------
    ValueError
        If ``axis`` is not an integer in ``-n_dim .. n_dim - 1``.
    """
    if axis is None:
        axis = slab_normal(traj)
    else:
        axis = _normalise_axis(axis, traj.n_dim, caller)
    aligned, cross_section = _slab_aligned_counts(traj, axis, stacklevel=stacklevel)
    length = traj.dimensions[axis]
    covered = np.isfinite(aligned)
    coverage = covered.sum(axis=0).astype(np.float64)
    acc = np.where(covered, aligned, 0.0).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        density = np.where(coverage > 0, acc / (np.maximum(coverage, 1) * cross_section), np.nan)
        mean_counts = acc / np.maximum(coverage, 1)
        squares = np.where(covered, (aligned - mean_counts) ** 2, 0.0).sum(axis=0)
        scatter = np.where(coverage > 1,
                           np.sqrt(squares / np.maximum(coverage - 1, 1)) / cross_section, np.nan)
    return np.arange(length).astype(np.float64), density, scatter


def slab_density_profile(traj, axis=None):
    """1D density profile (volume fraction) along ``axis`` (default: the slab
    normal picked by :func:`slab_normal`), with the dense slab re-centred each
    frame so it does not smear out as the slab diffuses.

    Under periodic boundaries the profile is rolled so the circular centre of
    mass of the counts sits at the window centre. Under a hardwall a free slab is
    translated (no wrap) so that the centroid of its dense bins sits at the
    window centre, while a condensate touching a wall is left where it is (it is
    pinned, and moving it would turn its flat wall face into a fictitious second
    interface). A translated frame has no data in the planes it vacated, so it is
    simply not counted there: each plane is averaged over the frames that cover
    it, and a plane no frame covers is ``nan``.

    A wall-wetting film and a free slab are different profiles, and averaging
    one kind of frame with the other gives neither (a fake film at the wall plus
    a diluted slab). So when a hardwall trajectory holds a slab in most frames
    and the condensate touches a wall in some frames but not in others, only the
    majority kind is averaged (a tie goes to the free slab) and a warning says
    how many frames were left out.

    A one-phase hardwall trajectory is neither translated nor split: fewer than
    half of its frames hold a slab (one dense block with a dilute region at
    least four planes wide beside it), so every frame is averaged exactly where
    it is and the profile is the plain per-plane mean, depletion layers at the
    walls included. Under periodic boundaries there is no such test and every
    frame is re-centred, which lines the fluctuations of a one-phase solution up
    into a shallow central hump; :func:`slab_density_profile_with_scatter`
    returns the frame-to-frame scatter that :func:`fit_slab_profile` needs to
    tell that hump from a coexistence gap.

    Every bead in the box is binned, dense and dilute alike - that is what makes the
    result a *density* profile that a coexistence fit can be run against. (It therefore
    takes no ``min_beads``: unlike the radial profile, the profile itself does not
    need clusters. It previously accepted a ``min_beads`` argument that was never read.)

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    axis : int, optional
        Index of the axis to profile along; negative values count from the last
        axis. Default ``None``: the axis :func:`slab_normal` picks, i.e. the one
        the largest cluster does not span when it spans all the others in most
        frames, and the longest box axis otherwise.

    Returns
    -------
    coordinate : numpy.ndarray
        ``(L,)`` float64 lattice coordinate along the axis, ``0 .. L-1``.
    density : numpy.ndarray
        ``(L,)`` float64 frame-averaged occupied fraction at each coordinate;
        ``nan`` at a hardwall plane that no kept frame covers after alignment.

    Raises
    ------
    ValueError
        If ``axis`` is not an integer in ``-n_dim .. n_dim - 1``.
    """
    coord, density, _scatter = _slab_profile_with_scatter(traj, axis, "slab_density_profile")
    return coord, density


def slab_density_profile_with_scatter(traj, axis=None):
    """:func:`slab_density_profile` plus the frame-to-frame scatter of each plane.

    The scatter is what tells coexistence from a fluctuation. Re-centring every
    frame on its densest region lines the density fluctuations of a one-phase
    solution up into a hump, and the average over many frames is smooth, so a
    ``tanh`` fits it well and the gap stands clear of the residuals of the
    *averaged* profile. What re-centring cannot do is make the hump taller than
    the fluctuations it is built from: the fitted gap of a one-phase solution
    is of the order of the scatter of a single plane's density from frame to
    frame, while a coexistence gap is many times it. :func:`fit_slab_profile`
    applies that test when it is handed the scatter. The test compares single
    planes, whose scatter grows as the cross-section shrinks, so it also
    rejects real slabs in narrow boxes; :func:`analyze` does not use it and
    counts the chains in the slab instead (:func:`_chain_excess_sigma`).

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    axis : int, optional
        Index of the axis to profile along, as for :func:`slab_density_profile`
        (default ``None``, the axis :func:`slab_normal` picks).

    Returns
    -------
    coordinate : numpy.ndarray
        ``(L,)`` float64 lattice coordinate along the axis, ``0 .. L-1``.
    density : numpy.ndarray
        ``(L,)`` float64 frame-averaged occupied fraction at each coordinate;
        ``nan`` at a hardwall plane that no kept frame covers after alignment.
    scatter : numpy.ndarray
        ``(L,)`` float64 sample standard deviation (``n - 1``) over frames of
        the occupied fraction of each plane, after alignment; ``nan`` where
        fewer than two frames cover the plane (so all ``nan`` for a one-frame
        trajectory). Frames saved every step or two are strongly correlated and
        make it read low.

    Raises
    ------
    ValueError
        If ``axis`` is not an integer in ``-n_dim .. n_dim - 1``.
    """
    return _slab_profile_with_scatter(traj, axis, "slab_density_profile_with_scatter")


# ---------------------------------------------------------------------------
# tanh binodal fits
# ---------------------------------------------------------------------------

def _tanh_droplet(r, rho_d, rho_v, radius, width):
    """Single-interface ``tanh`` density profile of a spherical droplet.

    Parameters
    ----------
    r : numpy.ndarray
        Float64 radii from the droplet centre, in lattice units.
    rho_d : float
        Dense-phase (droplet) density, as an occupied-site fraction.
    rho_v : float
        Dilute (vapour) density, as an occupied-site fraction.
    radius : float
        Droplet radius: the radius at which the profile is halfway between the
        two phases.
    width : float
        Interface width in lattice units.

    Returns
    -------
    numpy.ndarray
        Float64 model density at each ``r``.
    """
    return 0.5 * (rho_d + rho_v) - 0.5 * (rho_d - rho_v) * np.tanh((r - radius) / width)


def _tanh_slab(z, rho_d, rho_v, half_width, width, center):
    """Two-interface ``tanh`` density profile of a slab.

    Parameters
    ----------
    z : numpy.ndarray
        Float64 coordinate along the slab normal, in lattice units.
    rho_d : float
        Dense-phase (slab) density, as an occupied-site fraction.
    rho_v : float
        Dilute density, as an occupied-site fraction.
    half_width : float
        Half the slab thickness, so the two interfaces sit at
        ``center -/+ half_width``.
    width : float
        Interface width in lattice units.
    center : float
        Position of the slab centre along the axis.

    Returns
    -------
    numpy.ndarray
        Float64 model density at each ``z``.
    """
    return rho_v + 0.5 * (rho_d - rho_v) * (
        np.tanh((z - (center - half_width)) / width) - np.tanh((z - (center + half_width)) / width))


@dataclass
class BinodalFit:
    """Result of a ``tanh`` fit to a density profile.

    ``success`` means **the fit is well posed** - the data actually constrains the
    parameters - not merely that the optimiser converged. A converged-but-degenerate
    fit (inverted, no density gap, asymptotes never reached inside the data, a gap
    smaller than three times the residual scatter, a slab that fills the box, a
    slab under two planes thick or thinner than its interfaces, a slab profile
    with no dilute plateau, a gap that is not resolved in single frames, or a
    droplet radius below two lattice units) is reported with ``success=False``,
    because trusting it silently is the more dangerous failure; so are a fit with
    fewer than four usable shells or planes and one whose optimiser failed. When
    ``success`` is ``False``, ``rho_dense`` and ``rho_dilute`` fall back to the 95th
    and 5th percentiles of the *observed* profile, so they stay bounded and physically meaningful (for a homogeneous system
    they simply coincide) instead of being unconstrained extrapolations.

    ``reason`` is empty on success and otherwise names the check that failed.

    .. important::

       ``success`` is a **numerical** guarantee, not a physical one. It says the fit is
       meaningful *as a fit*; it does **not** decide whether the system is phase
       separated. A well-posed fit can still describe a shallow density modulation in a
       one-phase system - and near the critical point the coexistence gap closes
       continuously, so no numerical check can draw that line for you.

       For the physical question, use :attr:`PhaseSeparationResult.is_phase_separated`,
       or apply your own density-contrast threshold to ``rho_dense`` / ``rho_dilute``.
       Be aware too that :func:`slab_density_profile` re-centres the slab every frame,
       which aligns the fluctuations of even a *homogeneous* system into a shallow
       central hump - so a small, well-fitted density gap is not by itself evidence of
       coexistence. The radial profile does the same thing by centring every frame
       on the largest cluster. A fit run on the averaged profile alone cannot tell
       that hump from a condensate; one that is also given the frame-to-frame
       scatter (the ``frame_scatter`` argument of :func:`fit_slab_profile` and
       :func:`fit_radial_profile`, which :func:`analyze` always passes) can, and
       rejects it.

    Attributes
    ----------
    rho_dense : float
        Dense-phase coexistence density, as an occupied-site fraction.
    rho_dilute : float
        Dilute-phase coexistence density, as an occupied-site fraction.
    interface_width : float
        Fitted interface width in lattice units (``nan`` if no fit was run or
        the optimiser failed).
    radius : float
        Droplet radius in lattice units; ``nan`` in slab geometry.
    half_width : float
        Slab half-thickness in lattice units; ``nan`` in droplet geometry. For
        a condensate wetting a hardwall it is half the film thickness measured
        from the wall face, so ``2 * half_width`` is comparable with a free
        slab's thickness.
    success : bool
        Whether the fit is well posed (see above), not merely converged.
    reason : str
        Empty on success, otherwise the check that failed.
    """
    rho_dense: float
    rho_dilute: float
    interface_width: float
    radius: float = float("nan")        # droplet radius (spherical geometry)
    half_width: float = float("nan")    # slab half-width (slab geometry)
    success: bool = True
    reason: str = ""


def _observed_binodal(density):
    """Robust dense/dilute estimates read straight off the profile.

    Used as the fallback when the tanh fit is unusable. For a homogeneous profile the
    two values coincide, which is the correct answer for a one-phase system.

    Parameters
    ----------
    density : numpy.ndarray
        1D float64 density profile, already stripped of any ``nan`` or otherwise
        unusable bins.

    Returns
    -------
    rho_dense : float
        The 95th percentile of the profile, clipped to ``[0, 1]``.
    rho_dilute : float
        The 5th percentile of the profile, clipped to ``[0, 1]``.
    """
    lo, hi = np.percentile(density, [5, 95])
    return float(np.clip(hi, 0.0, 1.0)), float(np.clip(lo, 0.0, 1.0))


def _fit_is_usable(density, model_values, rho_dense, rho_dilute):
    """Is this converged fit actually a two-phase description of the profile?

    Two independent ways for a ``tanh`` fit to converge on nonsense, and one check for
    each:

    **1. The asymptotes are never reached.** A ``tanh`` with a very wide interface is,
    over a finite window, almost a straight line. The optimiser can then place
    ``rho_dense`` and ``rho_dilute`` almost anywhere - typically pinned at the bounds,
    1 and 0 - because nothing in the data constrains them. The fit converges, the
    residual is small, and the reported coexistence densities are pure extrapolation,
    attained nowhere in the box. The signature is that the model, evaluated over the
    data, traverses only a small part of the gap between its own asymptotes. A genuine
    profile spans essentially all of it: it really does reach the dense plateau in the
    middle and the dilute plateau at the edges.

    **2. The gap is noise.** Even on a flat profile a ``tanh`` can find a small hump in
    the scatter and fit it as a very thin slab. The asymptotes *are* then realised, so
    check 1 passes - but the "coexistence gap" is the size of the noise. So the gap
    must also stand clear of the residual scatter.

    Parameters
    ----------
    density : numpy.ndarray
        1D float64 observed density profile that was fitted.
    model_values : numpy.ndarray
        1D float64 fitted model evaluated at the same coordinates as ``density``.
    rho_dense : float
        The fit's dense-phase asymptote.
    rho_dilute : float
        The fit's dilute-phase asymptote.

    Returns
    -------
    ok : bool
        Whether the fit is a genuine two-phase description of the profile.
    reason : str
        Empty when ``ok``, otherwise the check that failed and why.
    """
    if float(rho_dense) < float(rho_dilute):
        return False, ("inverted fit: the fitted dense density is below the dilute one - the "
                       "model has fitted a dilute slab in a dense background (under HARDWALL a "
                       "condensate wetting both walls can look like this)")
    gap = abs(float(rho_dense) - float(rho_dilute))
    if gap < _MIN_DENSITY_GAP:
        return False, "no density gap: the profile is homogeneous"

    span = float(np.max(model_values) - np.min(model_values))
    if span < _MIN_REALISED_FRACTION * gap:
        return False, (
            f"degenerate fit: the model spans only {span:.3g} of its own {gap:.3g} "
            f"asymptotic gap inside the box, so rho_dense/rho_dilute are "
            f"extrapolations that the data does not constrain"
        )

    noise = float(np.std(np.asarray(density, float) - np.asarray(model_values, float)))
    if gap < _MIN_GAP_OVER_NOISE * noise:
        return False, (
            f"density gap ({gap:.3g}) is not significant against the scatter in the "
            f"profile ({noise:.3g}): this is noise, not coexistence"
        )
    return True, ""


def _gap_resolved_in_single_frames(model_values, gap, frame_scatter, what):
    """Is the fitted density gap larger than the fluctuations of single frames?

    Both profiles are centred, frame by frame, on the densest part of the
    system: the slab profile on the centre of mass of the density, the radial
    profile on the largest cluster. In a one-phase solution that lines the
    density fluctuations up into a hump, smooth after averaging and well fitted
    by a ``tanh``. The hump is about as tall as the fluctuations it is built
    from, while a coexistence gap is there in every frame. With ``s_dense`` and
    ``s_dilute`` the root-mean-square frame-to-frame scatter of the bins at the
    top and at the bottom of the fitted profile (within 10% of the gap of its
    maximum and of its minimum), the difference between one dense bin and one
    dilute bin in one frame has standard deviation ``sqrt(s_dense^2 +
    s_dilute^2)``, and the gap must be at least three times that.

    Parameters
    ----------
    model_values : numpy.ndarray
        1D float64 fitted model at the fitted bins (planes or shells).
    gap : float
        The fitted coexistence gap, ``rho_dense - rho_dilute``.
    frame_scatter : numpy.ndarray or None
        1D float64 frame-to-frame standard deviation of the density in the same
        bins; ``nan`` entries are ignored. ``None`` skips the test, as does a
        scatter with no finite entry at either end of the profile (e.g. that of
        a one-frame trajectory).
    what : str
        ``'plane'`` or ``'shell'``, for the message.

    Returns
    -------
    ok : bool
        Whether the gap is resolved (or the test could not be made).
    reason : str
        Empty when ``ok``, otherwise why not.
    """
    if frame_scatter is None:
        return True, ""
    model_values = np.asarray(model_values, float)
    frame_scatter = np.asarray(frame_scatter, float)
    top, bottom = float(model_values.max()), float(model_values.min())
    variance = 0.0
    n_measured = 0
    for end in (model_values >= top - _PLATEAU_TOLERANCE * gap,
                model_values <= bottom + _PLATEAU_TOLERANCE * gap):
        values = frame_scatter[end]
        values = values[np.isfinite(values)]
        if values.size:
            variance += float(np.mean(values ** 2))
            n_measured += 1
    sigma = float(np.sqrt(variance))
    if n_measured and gap < _MIN_GAP_OVER_FRAME_SCATTER * sigma:
        return False, (
            f"density gap ({gap:.3g}) is not resolved in single frames: the difference "
            f"between a dense and a dilute {what} scatters by {sigma:.3g} from frame to "
            f"frame, so the gap cannot be told apart from a fluctuation lined up by centring "
            f"every frame on its densest region (more frames, leaving out equilibration frames "
            f"or a wider cross-section would help)"
        )
    return True, ""


def _slab_plateau_checks(model_values, rho_dense, rho_dilute, frame_scatter=None,
                         min_dilute_planes=None):
    """Does a converged slab fit describe two phases that are both really there?

    Run after :func:`_fit_is_usable`, on slab profiles only. Four ways a slab
    fit passes those generic checks and is still not coexistence:

    **1. A slab one plane thick.** Its density and its thickness cannot be told
    apart, and the fit returned a denser, thinner slab than the true one with
    ``success=True``. At least two planes must lie inside the slab (fitted
    profile above the midpoint of the two densities).

    **2. A slab thinner than its interfaces.** The profile never gets near the
    dense asymptote, so ``rho_dense`` is an extrapolation. The fitted profile
    must rise at least three quarters of the way from ``rho_dilute`` to
    ``rho_dense``; for a free slab that is ``half_width > interface_width``.

    **3. No dilute plateau.** A one-phase solution between two hard walls is
    depleted in the few planes next to each wall. The two-interface ``tanh``
    fits that ramp as the edge of a slab that almost fills the box, with a
    dilute density well below anything observed (0.03 against a bulk of 0.10).
    The model reaches the dense asymptote but comes near the dilute one on no
    plane, or on one or two. At least three planes must lie within 10% of the
    gap of ``rho_dilute``.

    **4. A gap that single frames do not resolve.** See
    :func:`_gap_resolved_in_single_frames`: re-centring lines fluctuations up
    into a hump that the averaged profile cannot be told from a thin slab by,
    and the frame-to-frame scatter is the yardstick.

    Parameters
    ----------
    model_values : numpy.ndarray
        1D float64 fitted model at the fitted planes.
    rho_dense : float
        The fit's dense-phase asymptote.
    rho_dilute : float
        The fit's dilute-phase asymptote.
    frame_scatter : numpy.ndarray, optional
        1D float64 frame-to-frame standard deviation of the density at the same
        planes; ``nan`` entries are ignored. Default ``None``, which skips the
        fourth check (as does a scatter with no finite entry, e.g. that of a
        one-frame trajectory).
    min_dilute_planes : int, optional
        Fewest planes the dilute plateau may hold (default ``None``, i.e. 3).
        :func:`analyze` passes the extent of a chain along the normal instead,
        ``max(2, ceil(2 Rg_z))``: a dilute phase must be wider than a chain.

    Returns
    -------
    ok : bool
        Whether the fit passes all four checks.
    reason : str
        Empty when ``ok``, otherwise the check that failed and why.
    """
    model_values = np.asarray(model_values, float)
    rho_dense, rho_dilute = float(rho_dense), float(rho_dilute)
    gap = rho_dense - rho_dilute
    top = float(model_values.max())

    n_inside =int((model_values > rho_dilute + 0.5 * gap).sum())
    if n_inside < _MIN_SLAB_PLANES:
        return False, (
            "the slab is too thin: the fitted profile is above the midpoint of the two "
            "densities on %d plane(s) and at least %d are needed - the density and the "
            "thickness of a one-plane film cannot be told apart"
            % (n_inside, _MIN_SLAB_PLANES))
    reach = (top - rho_dilute) / gap
    if reach < _MIN_DENSE_REACH:
        return False, (
            "the slab is thinner than its interfaces: the fitted profile rises only "
            "%.0f%% of the way from rho_dilute to rho_dense (at least %.0f%% is needed, "
            "i.e. half_width > interface_width for a free slab), so rho_dense is an "
            "extrapolation - a thicker slab (more material) is needed to measure it"
            % (100 * reach, 100 * _MIN_DENSE_REACH))
    n_dilute = int((model_values <= rho_dilute + _PLATEAU_TOLERANCE * gap).sum())
    need = _MIN_DILUTE_PLATEAU_PLANES if min_dilute_planes is None else int(min_dilute_planes)
    if n_dilute < need:
        return False, (
            "no dilute plateau: the fitted profile comes within %.0f%% of the gap of "
            "rho_dilute on %d plane(s) and at least %d are needed (a dilute phase must be "
            "wider than a chain) - this cannot be told apart from a depletion layer at a "
            "wall or a condensate that fills the box; a longer box would help"
            % (100 * _PLATEAU_TOLERANCE, n_dilute, need))
    return _gap_resolved_in_single_frames(model_values, gap, frame_scatter, "plane")


def _min_shell_sites(n_dim):
    """Smallest Chebyshev shell worth fitting: a shell holds ~2 pi r sites in 2D
    and ~24 r^2 in 3D, so the 3D floor of 20 would discard every shell inside r = 3
    in 2D and leave a small 2D box with fewer than four usable shells.

    Parameters
    ----------
    n_dim : int
        Dimensionality of the system, 2 or 3.

    Returns
    -------
    int
        Minimum number of lattice sites a shell must hold to be fitted.
    """
    return _MIN_SHELL_SITES if n_dim == 3 else _MIN_SHELL_SITES_2D


def fit_radial_profile(radii, density, site_counts=None, min_shell_sites=None,
                       frame_scatter=None):
    """Fit a spherical droplet profile; returns a :class:`BinodalFit`.

    A fit that converges but does not describe a real two-phase profile is returned
    with ``success=False`` - see :class:`BinodalFit`.

    Shells whose density is ``nan`` (no lattice site there - see
    :func:`radial_density_profile`) are ignored. When ``site_counts`` is given (from
    :func:`radial_density_profile_with_site_counts`) shells with fewer than 20
    sites (or ``min_shell_sites``) are ignored too: the innermost shell is a single
    site that sits at the largest cluster's own centre, and fitting it with the same
    weight as a 1000-site shell let a homogeneous solution pass as a droplet of
    ``rho_dense = 1`` and radius ``< 1``. A fitted radius below two lattice units
    is rejected for the same reason.

    **Single frames.** ``frame_scatter`` (from
    :func:`radial_density_profile_with_scatter`) adds the test the averaged
    profile cannot supply. The profile is centred on the largest cluster in
    every frame, and in a dilute one-phase solution that cluster is a chain or a
    few touching chains, denser than its surroundings: the averaged profile has
    a hump at the origin, a ``tanh`` fits it, and it used to be reported as a
    droplet of radius 3 and density 0.15 - 0.2. A coil's density fluctuates from
    frame to frame by about as much as the hump is tall, while a droplet's core
    reads the same in every frame. So the gap must be at least three standard
    deviations of the difference between one core shell and one outer shell in
    a single frame. Without ``frame_scatter`` that test is skipped.

    Parameters
    ----------
    radii : array_like
        1D shell-centre radii in lattice units, from
        :func:`radial_density_profile`.
    density : array_like
        1D density profile at those radii; ``nan`` marks a shell with no lattice
        site and is dropped.
    site_counts : array_like, optional
        1D number of lattice sites per shell, from
        :func:`radial_density_profile_with_site_counts`. Default ``None``, which
        fits every finite shell however few sites it holds.
    min_shell_sites : int, optional
        Override the site-count floor applied when ``site_counts`` is given
        (default ``None``, i.e. 20, the 3D floor). Pass 8 for a 2D profile,
        whose shells hold only about ``2 pi r`` sites; :func:`analyze` picks the
        floor from ``traj.n_dim`` itself.
    frame_scatter : array_like, optional
        1D frame-to-frame standard deviation of the density of each shell, from
        :func:`radial_density_profile_with_scatter`; ``nan`` entries are
        ignored. Default ``None``, which skips the single-frame test.

    Returns
    -------
    BinodalFit
        The fit, with ``success=False`` and a populated ``reason`` when fewer
        than four shells are usable, the optimiser fails, or the converged fit
        is degenerate (including a profile that ends before its dilute
        plateau, and a gap that single frames do not resolve). In those cases ``rho_dense`` / ``rho_dilute`` fall back to percentiles of
        the observed profile.

    Raises
    ------
    ValueError
        If ``frame_scatter`` is given and does not have one entry per shell.
    """
    from scipy.optimize import curve_fit
    radii = np.asarray(radii, float)
    density = np.asarray(density, float)
    scatter = None
    if frame_scatter is not None:
        scatter = np.asarray(frame_scatter, float)
        if scatter.shape != density.shape:
            raise ValueError(
                "fit_radial_profile: frame_scatter must have one entry per shell "
                "(%d), got shape %r" % (density.size, scatter.shape))
    keep = np.isfinite(density)
    if site_counts is not None:
        floor = _MIN_SHELL_SITES if min_shell_sites is None else int(min_shell_sites)
        keep &= np.asarray(site_counts, float) >= floor
    if keep.sum() < 4:
        # the fallback percentiles are read off the USABLE shells only: the one-site
        # shell at the cluster's own centre of mass reads ~0.7 even for a
        # homogeneous solution and would double the reported dense density
        usable = density[keep] if keep.any() else density[np.isfinite(density)]
        obs_d, obs_v = _observed_binodal(usable) if usable.size else (float("nan"), float("nan"))
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v, interface_width=float("nan"),
                          radius=float("nan"), success=False,
                          reason="fewer than four usable radial shells")
    radii = radii[keep]
    density = density[keep]
    if scatter is not None:
        scatter = scatter[keep]

    rho_d0 = float(density[:max(1, len(density) // 5)].mean())
    rho_v0 = float(density[-max(1, len(density) // 5):].mean())
    r0 = float(radii[np.argmin(np.abs(density - 0.5 * (rho_d0 + rho_v0)))])
    obs_d, obs_v = _observed_binodal(density)
    try:
        p, _ = curve_fit(_tanh_droplet, radii, density,
                         p0=[rho_d0, rho_v0, r0, 1.0],
                         bounds=([0, 0, 0, 0.1], [1, 1, radii.max(), radii.max()]),
                         maxfev=10000)
    except Exception:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float("nan"), radius=r0, success=False,
                          reason="curve_fit did not converge")

    ok, reason = _fit_is_usable(density, _tanh_droplet(radii, *p), p[0], p[1])
    if ok and float(p[2]) < _MIN_DROPLET_RADIUS:
        ok, reason = False, (
            f"degenerate fit: the fitted droplet radius ({float(p[2]):.2f} lattice units) "
            f"is a handful of sites at the cluster centre, not a dense phase")
    if ok:
        # the profile ends at half the shortest box axis; in a small box that is
        # before any dilute plateau, and rho_dilute is then read off the tail of
        # the tanh rather than off the data (it moved between 0.000 and 0.072
        # when single frames were left out of one 2D run)
        gap = float(p[0]) - float(p[1])
        tail = float(np.min(_tanh_droplet(radii, *p)))
        if tail > float(p[1]) + _PLATEAU_TOLERANCE * gap:
            ok, reason = False, (
                f"rho_dilute is an extrapolation: the fitted profile is still "
                f"{100 * (tail - float(p[1])) / gap:.0f}% of the gap above it at the outermost "
                f"usable shell (which reads {float(density[-1]):.3g}), so the profile ends "
                f"before any dilute plateau; a larger box would help")
    if ok:
        ok, reason = _gap_resolved_in_single_frames(
            _tanh_droplet(radii, *p), float(p[0]) - float(p[1]), scatter, "shell")
    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float(p[3]), radius=float(p[2]),
                          success=False, reason=reason)

    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      radius=float(p[2]), interface_width=float(p[3]))


def _fit_wetting_slab(coord, density, window=None, frame_scatter=None, min_dilute_planes=None):
    """Single-interface fit for a condensate that wets one wall of a HARDWALL box.

    The profile then has ONE interface (the wall face is not an interface at all),
    so it is fit with the droplet form in the distance from the wetted wall,
    ``rho(z') = 1/2(rho_d + rho_v) - 1/2(rho_d - rho_v) tanh((z' - t)/w)`` with
    ``t`` the slab thickness; ``half_width`` is reported as ``t/2`` so that
    ``2 * half_width`` is the thickness, as for the two-interface form.

    ``z'`` is measured from the **wall face**, not from the centre of the first
    site. Lattice sites are unit cells centred on integer coordinates, so the
    excluding face of the site at ``coord[0]`` sits at ``coord[0] - 0.5``, and a
    film filling ``t`` planes fills the volume from ``coord[0] - 0.5`` to
    ``coord[0] + t - 0.5``. That is the convention the two-interface fit already
    uses, so a wetting film and a free slab of the same material now report the
    same thickness.

    The film is checked as a slab, with :func:`_slab_plateau_checks`. It used to
    be handed to :func:`fit_radial_profile` and so inherited the droplet check
    (a fitted radius below two lattice units is a handful of sites): a film of
    one or two planes was rejected with a message about a droplet radius, and
    because the fitted thickness of a sharp two-plane film can land anywhere
    between 1.5 and 2.5, a two-plane film at the low wall and its mirror image
    at the high wall got opposite verdicts.

    Parameters
    ----------
    coord : numpy.ndarray
        1D float64 lattice coordinate along the axis.
    density : numpy.ndarray
        1D float64 density profile at those coordinates. Which wall is wetted is
        read off its two ends.
    window : tuple of float, optional
        ``(first, last)`` coordinate of the full window, whose outer faces are
        the walls (default ``None``, i.e. ``coord[0]`` and ``coord[-1]``). It
        differs from those only when planes with no data were dropped from the
        ends of the profile before fitting.
    frame_scatter : numpy.ndarray, optional
        1D float64 frame-to-frame scatter of the density at the same planes,
        for the single-frame check of :func:`_slab_plateau_checks` (default
        ``None``, which skips it).
    min_dilute_planes : int, optional
        Fewest planes the dilute plateau may hold (default ``None``, i.e. 3).

    Returns
    -------
    BinodalFit
        The fit, with ``half_width`` set to half the fitted film thickness and
        ``radius`` left at ``nan``. ``success`` is ``False`` with a populated
        ``reason`` when the optimiser fails, the fit is degenerate, or the film
        fails one of the slab checks.
    """
    from scipy.optimize import curve_fit
    # The +0.5 puts the origin ON THE WALL. Without it the distance was measured
    # from the CENTRE of the first site, half a lattice unit inside the wall, and
    # every wetting film came back exactly 0.5 units too thin - a 5% error at 10
    # planes but 17% at 3, systematic and one-signed, and inconsistent with the
    # two-interface fit, which reports 10.00 for the same 10-plane film.
    first, last = (float(coord[0]), float(coord[-1])) if window is None else window
    if density[0] >= density[-1]:
        z = coord - first + 0.5              # wets the low wall: profile falls with z
    else:
        z = last - coord + 0.5               # wets the high wall: mirror it
    order = np.argsort(z)
    z = z[order]
    rho = density[order]
    scatter = None if frame_scatter is None else np.asarray(frame_scatter, float)[order]

    # starting point and bounds are those of fit_radial_profile, which this
    # routine used to call, so a film that passed before is fitted identically
    rho_d0 = float(rho[:max(1, len(rho) // 5)].mean())
    rho_v0 = float(rho[-max(1, len(rho) // 5):].mean())
    t0 = float(z[np.argmin(np.abs(rho - 0.5 * (rho_d0 + rho_v0)))])
    obs_d, obs_v = _observed_binodal(rho)
    try:
        p, _ = curve_fit(_tanh_droplet, z, rho, p0=[rho_d0, rho_v0, t0, 1.0],
                         bounds=([0, 0, 0, 0.1], [1, 1, z.max(), z.max()]), maxfev=10000)
    except Exception:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v, interface_width=float("nan"),
                          half_width=0.5 * t0, success=False,
                          reason="curve_fit did not converge")

    model = _tanh_droplet(z, *p)
    ok, reason = _fit_is_usable(rho, model, p[0], p[1])
    if ok:
        ok, reason = _slab_plateau_checks(model, p[0], p[1], scatter, min_dilute_planes)
    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v, interface_width=float(p[3]),
                          half_width=0.5 * float(p[2]), success=False, reason=reason)
    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      interface_width=float(p[3]), half_width=0.5 * float(p[2]))


def _fit_both_walls(coord, density, window, frame_scatter=None, min_dilute_planes=None):
    """Fit a condensate that wets BOTH walls of a HARDWALL box along the normal.

    The profile is then dense at both ends and dilute in the middle: two films
    with a dilute region between them. That is a slab of the *dilute* phase in
    a dense background, so it is fit with the two-interface form with the roles
    of the two densities exchanged. It used to be sent to the single-film fit,
    which cannot describe a profile that rises again, and came back as "this is
    noise, not coexistence" for a perfectly clean pair of films.

    Parameters
    ----------
    coord : numpy.ndarray
        1D float64 lattice coordinate along the axis (planes with data only).
    density : numpy.ndarray
        1D float64 density profile at those coordinates.
    window : tuple of float
        ``(first, last)`` coordinate of the full window, whose outer faces are
        the two walls.
    frame_scatter : numpy.ndarray, optional
        1D float64 frame-to-frame scatter of the density at the same planes,
        for the single-frame check of :func:`_slab_plateau_checks` (default
        ``None``, which skips it).
    min_dilute_planes : int, optional
        Fewest planes the dilute plateau may hold (default ``None``, i.e. 3).

    Returns
    -------
    BinodalFit
        The fit. ``rho_dense`` is the density of the films and ``rho_dilute``
        that of the region between them. The two films need not be equally
        thick; ``half_width`` is half their MEAN thickness, each measured from
        its wall face, so that ``2 * half_width`` is comparable with the
        thickness of a single film or a free slab. ``success`` is ``False``
        with a populated ``reason`` when the optimiser fails, the fit is
        degenerate, a film has no thickness, or a plateau is missing.
    """
    from scipy.optimize import curve_fit
    first, last = window
    length = last - first + 1
    rho_d0 = float(max(density[0], density[-1]))
    rho_v0 = float(density.min())
    below = density < 0.5 * (rho_d0 + rho_v0)
    g0 = max(float(below.sum()) / 2.0, 1.0)
    weights = np.clip(rho_d0 - density, 0.0, None)
    c0 = float((coord * weights).sum() / weights.sum()) if weights.sum() > 0 else first + length / 2.0
    obs_d, obs_v = _observed_binodal(density)

    def _model(z, rd, rv, g, w, c):
        """Dilute region of half-width ``g`` centred at ``c`` between two films.

        Parameters
        ----------
        z : numpy.ndarray
            Float64 coordinate along the normal.
        rd : float
            Density of the two films.
        rv : float
            Density of the dilute region between them.
        g : float
            Half-width of the dilute region in lattice units.
        w : float
            Interface width in lattice units.
        c : float
            Position of the centre of the dilute region.

        Returns
        -------
        numpy.ndarray
            Float64 model density at each ``z``.
        """
        return _tanh_slab(z, rv, rd, g, w, c)

    try:
        p, _ = curve_fit(_model, coord, density, p0=[rho_d0, rho_v0, g0, 1.0, c0],
                         bounds=([0, 0, 0, 0.1, float(coord[0])],
                                 [1, 1, length, length, float(coord[-1])]), maxfev=10000)
    except Exception:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v, interface_width=float("nan"),
                          success=False, reason="curve_fit did not converge")

    model = _model(coord, *p)
    low_film = (float(p[4]) - float(p[2])) - (first - 0.5)
    high_film = (last + 0.5) - (float(p[4]) + float(p[2]))
    half_width = 0.25 * (low_film + high_film)
    ok, reason = _fit_is_usable(density, model, p[0], p[1])
    if ok and min(low_film, high_film) <= 0.0:
        ok, reason = False, ("degenerate fit: the dilute region between the two wall films "
                             "reaches a wall, so one of the films has no thickness")
    if ok:
        ok, reason = _slab_plateau_checks(model, p[0], p[1], frame_scatter, min_dilute_planes)
    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v, interface_width=float(p[3]),
                          half_width=half_width, success=False, reason=reason)
    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      interface_width=float(p[3]), half_width=half_width)


def fit_slab_profile(coord, density, hardwall=None, frame_scatter=None,
                     min_dilute_planes=None):
    """Fit a slab (two-interface) profile; returns a :class:`BinodalFit`.

    A fit that converges but does not describe a real two-phase profile is returned
    with ``success=False`` - see :class:`BinodalFit`. Above the critical temperature
    the profile is flat, and an unguarded ``tanh`` fit will happily report a large,
    entirely fictitious density gap for it.

    A profile whose dense phase sits against one **end** of the window can only
    come from a HARDWALL box (a periodic profile is re-centred by
    :func:`slab_density_profile`): the condensate wets a wall, its wall face is not
    an interface, and it is fit with a **single**-interface model - the
    two-interface form used to fit the wall step as a second, spuriously sharp
    interface. The wetting model is selected when an end of the window is dense
    (at least half the peak density) AND denser than the middle of the window;
    ``hardwall=False`` never uses it, ``hardwall=True`` and the default ``None``
    both decide from the profile shape (a hardwall slab that sits away from the
    walls, with a dilute phase at or above half the dense density, still has two
    interfaces and must not be routed to the single-interface fit). A profile
    dense at BOTH ends - a condensate wetting both walls, with the dilute phase
    between the two films - is fit as a dilute slab in a dense background, and
    ``half_width`` is then half the mean thickness of the two films.

    **Both phases must be there.** A slab fit also has to show the two phases
    it reports. The slab must be at least two planes thick (a one-plane film's
    density and thickness cannot be told apart) and thicker than its interfaces
    (the fitted profile rises at least three quarters of the way to
    ``rho_dense``, i.e. ``half_width > interface_width`` for a free slab;
    otherwise ``rho_dense`` is an extrapolation). And the fitted profile must
    come within 10% of the gap of ``rho_dilute`` on at least three planes: the
    depletion layer of a one-phase solution at a hard wall is a ramp of a few
    planes, which the ``tanh`` fits as the edge of a slab filling the box, with
    a "dilute density" below anything observed.

    **Single frames.** ``frame_scatter`` (from
    :func:`slab_density_profile_with_scatter`) adds the one test the averaged
    profile cannot supply. Re-centring every frame on its densest region lines
    the fluctuations of a one-phase solution up into a smooth hump, but cannot
    make it taller than those fluctuations; a coexistence gap is there in every
    frame. So the gap must be at least three standard deviations of the
    difference between one dense plane and one dilute plane in a single frame
    (the two planes' frame-to-frame scatters added in quadrature). Without
    ``frame_scatter`` that test is skipped and a well-fitted hump passes: pass
    it whenever the profile comes from a trajectory.

    Planes whose density is ``nan`` (a hardwall plane that no frame covered
    after :func:`slab_density_profile` aligned the slab) are left out of the fit.
    The window itself - its walls for the wetting model, its centre for a
    periodic profile - is still read from the full ``coord``.

    Parameters
    ----------
    coord : array_like
        1D lattice coordinate along the slab normal, from
        :func:`slab_density_profile`.
    density : array_like
        1D density profile at those coordinates; ``nan`` marks a plane with no
        data and is dropped.
    hardwall : bool, optional
        Whether the box has hard walls. ``False`` fixes the slab centre at the
        middle of the window (a periodic profile is re-centred per frame) and
        never uses the wetting model; ``True`` and the default ``None`` both fit
        the centre as a free parameter and pick the wetting model from the
        profile shape.
    frame_scatter : array_like, optional
        1D frame-to-frame standard deviation of the density at each coordinate,
        from :func:`slab_density_profile_with_scatter`; ``nan`` entries are
        ignored. Default ``None``, which skips the single-frame test. The test
        compares single planes, whose scatter grows as the cross-section
        shrinks: it rejects real slabs in narrow boxes (an 8-site-wide 2D
        stripe with a gap of 0.75), which is why :func:`analyze` does not use
        it and tests the number of chains in the slab instead.
    min_dilute_planes : int, optional
        Fewest planes the dilute plateau may hold (default ``None``, i.e. 3).
        :func:`analyze` passes ``max(2, ceil(2 Rg_z))`` with ``Rg_z`` the
        chains' radius of gyration along the normal: a dilute phase must be
        wider than a chain.

    Returns
    -------
    BinodalFit
        The fit, with ``success=False`` and a populated ``reason`` when fewer
        than four planes have data, the optimiser fails or the converged fit is
        degenerate (including a slab that fills the box, leaving no dilute
        phase, a slab under two planes thick or thinner than its interfaces, a
        profile with no dilute plateau, and a gap that single frames do not
        resolve).

    Raises
    ------
    ValueError
        If ``frame_scatter`` is given and does not have one entry per
        coordinate.
    """
    from scipy.optimize import curve_fit
    coord = np.asarray(coord, float)
    density = np.asarray(density, float)
    # the periodic slab centre is the middle of the WHOLE window, which is where
    # slab_density_profile put it, whatever planes are missing
    window_center = float(coord[0] + (coord[-1] - coord[0] + 1) / 2.0) if coord.size else 0.0
    window = (float(coord[0]), float(coord[-1])) if coord.size else (0.0, 0.0)
    keep = np.isfinite(density)
    scatter = None
    if frame_scatter is not None:
        scatter = np.asarray(frame_scatter, float)
        if scatter.shape != density.shape:
            raise ValueError(
                "fit_slab_profile: frame_scatter must have one entry per coordinate "
                "(%d), got shape %r" % (density.size, scatter.shape))
    if keep.sum() < 4:
        usable = density[keep]
        obs_d, obs_v = _observed_binodal(usable) if usable.size else (float("nan"), float("nan"))
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v, interface_width=float("nan"),
                          success=False, reason="fewer than four planes with data")
    coord = coord[keep]
    density = density[keep]
    if scatter is not None:
        scatter = scatter[keep]
    if hardwall is not False and density.size >= 4:
        peak = float(density.max())
        mid = float(density[density.size // 4: -max(1, density.size // 4)].mean())
        end_dense = peak > 0 and (density[0] >= 0.5 * peak or density[-1] >= 0.5 * peak)
        # dense at BOTH ends and dilute between: a condensate wetting both walls
        if (peak > 0 and min(density[0], density[-1]) >= 0.5 * peak
                and min(density[0], density[-1]) > mid):
            return _fit_both_walls(coord, density, window, frame_scatter=scatter,
                                   min_dilute_planes=min_dilute_planes)
        # the ends must also be denser than the middle: a re-centred periodic
        # slab, or a hardwall slab away from the walls whose dilute phase is at
        # least half the dense density, has two interfaces and stays on the
        # two-interface fit whatever the hardwall flag says
        if end_dense and max(density[0], density[-1]) > mid:
            return _fit_wetting_slab(coord, density, window=window, frame_scatter=scatter,
                                     min_dilute_planes=min_dilute_planes)
    length = coord[-1] - coord[0] + 1
    rho_d0 = float(density.max())
    rho_v0 = float(np.median(density[density < 0.5 * density.max()])) if np.any(density < 0.5 * density.max()) else 0.0
    hw0 = float(np.sum(density > 0.5 * (rho_d0 + rho_v0)) / 2.0)
    obs_d, obs_v = _observed_binodal(density)

    # The slab centre is FIXED at the middle of the window only for a periodic
    # profile, which slab_density_profile re-centres; a hardwall profile is not
    # re-centred (a wall-wetting film cannot be rolled, and a free slab is only
    # aligned when it clears the walls), so there the centre is a fitted
    # parameter, initialised from the density-weighted centroid. Anchor on
    # coord[0] rather than assuming the coordinate axis starts at zero (it does
    # for slab_density_profile, but not necessarily for a caller's own profile).
    if hardwall is False:
        center = window_center

        def _model(z, rd, rv, hw, w):
            """Two-interface slab with the centre pinned at the window centre.

            Parameters
            ----------
            z : numpy.ndarray
                Float64 coordinate along the slab normal.
            rd : float
                Dense-phase density.
            rv : float
                Dilute-phase density.
            hw : float
                Slab half-width in lattice units.
            w : float
                Interface width in lattice units.

            Returns
            -------
            numpy.ndarray
                Float64 model density at each ``z``.
            """
            return _tanh_slab(z, rd, rv, hw, w, center)

        p0 = [rho_d0, rho_v0, max(hw0, 1.0), 1.0]
        bounds = ([0, 0, 0, 0.1], [1, 1, length, length])
    else:
        weights = np.clip(density - rho_v0, 0.0, None)
        c0 = float((coord * weights).sum() / weights.sum()) if weights.sum() > 0 else float(coord[0] + length / 2.0)

        def _model(z, rd, rv, hw, w, c):
            """Two-interface slab with the centre fitted as a free parameter.

            Parameters
            ----------
            z : numpy.ndarray
                Float64 coordinate along the slab normal.
            rd : float
                Dense-phase density.
            rv : float
                Dilute-phase density.
            hw : float
                Slab half-width in lattice units.
            w : float
                Interface width in lattice units.
            c : float
                Position of the slab centre along the axis.

            Returns
            -------
            numpy.ndarray
                Float64 model density at each ``z``.
            """
            return _tanh_slab(z, rd, rv, hw, w, c)

        p0 = [rho_d0, rho_v0, max(hw0, 1.0), 1.0, c0]
        bounds = ([0, 0, 0, 0.1, float(coord[0])], [1, 1, length, length, float(coord[-1])])

    try:
        p, _ = curve_fit(_model, coord, density, p0=p0, bounds=bounds, maxfev=10000)
    except Exception:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float("nan"), half_width=hw0, success=False,
                          reason="curve_fit did not converge")

    ok, reason = _fit_is_usable(density, _model(coord, *p), p[0], p[1])

    # A "slab" that fills the box has no dilute phase to coexist with: the fitted
    # dilute density is then unconstrained, whatever the profile happens to look like.
    if ok and 2.0 * float(p[2]) >= length:
        ok, reason = False, "degenerate fit: the slab fills the box, leaving no dilute phase"

    # both phases must show a plateau, and the gap must be there in single frames
    if ok:
        ok, reason = _slab_plateau_checks(_model(coord, *p), p[0], p[1], scatter,
                                          min_dilute_planes)

    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float(p[3]), half_width=float(p[2]),
                          success=False, reason=reason)

    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      half_width=float(p[2]), interface_width=float(p[3]))


# ---------------------------------------------------------------------------
# droplet shape (frame-averaged over the largest cluster)
# ---------------------------------------------------------------------------

def droplet_shape(traj, min_beads=2):
    """Frame-averaged largest-cluster geometry.

    A frame whose largest cluster spans the box is left out of every average. Such
    a cluster is a network or a slab, not a droplet. Under periodic boundaries it
    is connected to its own periodic image: it has no single image, the gather
    hands back one search-order dependent window of it, and the radius of
    gyration, asphericity and hull of that window describe the search, not the
    cluster. Under a hardwall it touches both walls: its image is exact, but its
    shape is that of a wall-bounded condensate. The usual case is frame 0 of
    a run that saved its equilibration, where the random starting placement
    percolates as a contact network; averaged in with the droplet frames, an
    asphericity of 15 for that one frame moves a ten-frame mean from 0.0 to 1.5. A
    warning says how many frames were left out; :func:`spanning_fraction` reports
    the same count as a fraction.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    min_beads : int, optional
        Ignore clusters smaller than this many beads (default ``2``).

    Returns
    -------
    dict
        ``radius_of_gyration``, ``asphericity``, ``sphericity``, ``volume`` and
        ``density`` of the largest cluster, each a float averaged over the
        frames that contain one that does not span the box (``nan`` if none
        do). Degenerate convex-hull values (``-1``) are left out of the
        averages.
    """
    rg, asph, sph, vol, dens = [], [], [], [], []
    n_spanning = 0
    for _f, clusters in _largest_clusters(traj, min_beads):
        if not clusters:
            continue
        c = clusters[0]
        if _cluster_spans_box(c, False):
            n_spanning += 1
            continue
        rg.append(c.radius_of_gyration)
        asph.append(c.asphericity)
        sph.append(c.sphericity)
        vol.append(c.volume)
        dens.append(c.density)

    if n_spanning:
        warnings.warn(
            "droplet_shape: %d of %d frames with a cluster were left out of the averages "
            "because the largest cluster spans the box (%s) - a network or a slab, not a "
            "droplet, whose shape statistics do not describe a droplet; see "
            "phase_separation.spanning_fraction."
            % (n_spanning, n_spanning + len(rg), "touches both walls"
               if bool(getattr(traj, "hardwall", False)) else
               "connected to its own periodic image"), stacklevel=2)

    def _mean(x):
        """Mean of the usable entries of one per-frame series.

        Parameters
        ----------
        x : list of float
            One value per frame that held a cluster.

        Returns
        -------
        float
            Mean of the finite entries above ``-1`` (the degenerate convex-hull
            sentinel), or ``nan`` if none qualify.
        """
        x = np.asarray(x, float)
        x = x[np.isfinite(x) & (x > -1)]        # drop degenerate (-1) hull values
        return float(x.mean()) if x.size else float("nan")

    return {"radius_of_gyration": _mean(rg), "asphericity": _mean(asph),
            "sphericity": _mean(sph), "volume": _mean(vol), "density": _mean(dens)}


# ---------------------------------------------------------------------------
# chains as the unit of fluctuation
# ---------------------------------------------------------------------------


def _chain_extent_along(traj, axis, max_frames=10):
    """Mean radius of gyration of the chains along one axis, in lattice units.

    Each chain is made whole by walking its bonds (minimum image under periodic
    boundaries), and the standard deviation of its bead coordinates along
    ``axis`` is averaged over all chains and over up to ``max_frames`` evenly
    spaced frames.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    axis : int
        Index of the axis.
    max_frames : int, optional
        Number of frames sampled (default ``10``).

    Returns
    -------
    float
        The mean extent, ``0.0`` for a system of single beads.
    """
    offsets = np.asarray(traj.topology.offsets)
    length = traj.dimensions[axis]
    hardwall = bool(getattr(traj, "hardwall", False))
    positions = traj.positions
    step = max(1, traj.n_frames // max_frames)
    extents = []
    for f in range(0, traj.n_frames, step):
        z = positions[f][:, axis].astype(np.float64)
        for c in range(len(offsets) - 1):
            zz = z[offsets[c]:offsets[c + 1]]
            if zz.size < 2:
                extents.append(0.0)
                continue
            dz = np.diff(zz)
            if not hardwall:
                dz = dz - length * np.round(dz / length)
            extents.append(float(np.std(np.concatenate([[0.0], np.cumsum(dz)]))))
    return float(np.mean(extents)) if extents else 0.0


def _chain_excess_sigma(density, cross_section, rho_dense, rho_dilute, chain_lengths):
    """How many more chains the dense region holds than a uniform solution would.

    Coexistence is a statement about many chains. The independent unit of a
    polymer solution's density fluctuations is the chain, not the bead: beads
    arrive in the dense region a chain at a time. In a one-phase solution of
    ``N`` chains each chain is in the dense region (a fraction ``f`` of the
    planes) with probability ``f``, so the region holds ``N f`` chains give or
    take ``sqrt(N f (1 - f))``. Re-centring every frame on its densest region
    selects the upward fluctuation, which is why a one-phase profile shows a
    hump, but the hump is of that size: one to two standard deviations. A
    condensate holds most of the chains in a small part of the box, many
    standard deviations more.

    Counting chains rather than comparing plane densities makes the test
    independent of the cross-section in the right way: a narrow box has noisy
    planes but its slab still holds most of its chains, and a lump made of two
    or three long coils is not significant however smooth its averaged profile.

    Parameters
    ----------
    density : numpy.ndarray
        1D float64 frame-averaged slab profile (``nan`` planes are ignored).
    cross_section : int
        Lattice sites per plane.
    rho_dense : float
        Fitted dense-phase density.
    rho_dilute : float
        Fitted dilute-phase density.
    chain_lengths : numpy.ndarray
        1D beads per chain.

    Returns
    -------
    sigma : float
        ``(K - N f) / sqrt(N f (1 - f))``, with ``K`` the beads in the dense
        region (profile above the midpoint of the two densities) in units of
        the bead-weighted mean chain length; ``0.0`` if there is no dense
        region or it covers every plane.
    n_dense_chains : float
        ``K``.
    n_expected : float
        ``N f``, with ``N`` the beads of the system in the same units.
    """
    density = np.asarray(density, float)
    density = density[np.isfinite(density)]
    lengths = np.asarray(chain_lengths, float)
    mean_length = float((lengths ** 2).sum() / lengths.sum())
    n_chains = float(lengths.sum() / mean_length)
    inside = density > rho_dilute + 0.5 * (rho_dense - rho_dilute)
    fraction = float(inside.sum()) / density.size if density.size else 0.0
    if fraction <= 0.0 or fraction >= 1.0:
        return 0.0, 0.0, n_chains * fraction
    n_dense = float(density[inside].sum() * cross_section / mean_length)
    expected = n_chains * fraction
    return ((n_dense - expected) / np.sqrt(expected * (1.0 - fraction)), n_dense, expected)


# ---------------------------------------------------------------------------
# top-level summary
# ---------------------------------------------------------------------------

@dataclass
class PhaseSeparationResult:
    """Everything :func:`analyze` measured about one trajectory.

    Attributes
    ----------
    geometry : str
        The geometry actually analysed, ``'slab'`` or ``'sphere'`` (never
        ``'auto'``).
    condensed_fraction : float
        Trajectory-averaged fraction of beads in the largest cluster.
    condensed_fraction_series : numpy.ndarray
        ``(n_frames,)`` float64 per-frame condensed fraction behind that mean.
    n_clusters : float
        Trajectory-averaged number of clusters of at least ``min_beads`` beads.
    largest_cluster_beads : float
        Trajectory-averaged bead count of the largest cluster.
    binodal : BinodalFit
        The coexistence fit to the density profile.
    shape : dict or None
        :func:`droplet_shape` statistics over the frames whose largest cluster
        does not span the box, or ``None`` in slab geometry, where hull quantities
        of a box-spanning slab are meaningless.
    profile : tuple
        The fitted density profile as ``(coordinate, density)``: shell radii in
        droplet geometry, the axis coordinate in slab geometry.
    percolation_fraction : float
        Fraction of frames whose largest cluster spans the box on EVERY axis (a
        space-filling network: the contact clustering of a homogeneous
        solution).
    spanning_fraction : float
        Fraction of frames whose largest cluster spans the box on ANY axis:
        under periodic boundaries it then has no single periodic image, and
        under a hardwall it touches both walls of that axis.
    slab_axis : int or None
        Slab geometry only: index of the axis the profile was taken along (the
        slab normal). ``None`` in droplet geometry.
    chain_excess_sigma : float
        Slab geometry only: by how many standard deviations the number of
        chains in the dense region exceeds what a uniform solution would put
        there (see :func:`_chain_excess_sigma`); ``nan`` when the fit failed
        before that test, and in droplet geometry.
    """
    geometry: str
    condensed_fraction: float
    condensed_fraction_series: np.ndarray = field(repr=False)
    n_clusters: float
    largest_cluster_beads: float
    binodal: BinodalFit
    shape: dict
    profile: tuple = field(repr=False, default=None)
    #: fraction of frames in which the largest cluster spans the box on EVERY axis
    #: (a space-filling network: the contact clustering of a homogeneous solution)
    percolation_fraction: float = 0.0
    #: fraction of frames in which it spans the box on ANY axis
    spanning_fraction: float = 0.0
    #: slab geometry: the axis the profile was taken along (the slab normal)
    slab_axis: int = None
    #: slab geometry: standard deviations by which the chain count of the dense
    #: region exceeds that of a uniform solution (nan if no fit was tested)
    chain_excess_sigma: float = float("nan")

    @property
    def rho_dense(self):
        """float : Dense-phase coexistence density, from :attr:`binodal`."""
        return self.binodal.rho_dense

    @property
    def rho_dilute(self):
        """float : Dilute-phase coexistence density, from :attr:`binodal`."""
        return self.binodal.rho_dilute

    @property
    def is_phase_separated(self):
        """bool : Whether the trajectory looks phase separated.

        Heuristic: a *usable* binodal fit, a clear density gap, most material
        condensed, and the "condensate" not a box-filling network.

        In slab geometry "usable" includes the slab checks of
        :func:`fit_slab_profile` (a slab thicker than its interfaces, a dilute
        plateau wider than a chain) and the chain-number test of
        :func:`analyze`: the slab must hold at least three standard deviations
        more chains than a uniform solution would put in the same planes. They
        are what keeps a one-phase solution in an elongated box from passing: an
        athermal solution between hard walls has a depletion layer at each wall
        that fits as a "dilute phase", and under periodic boundaries re-centring
        builds a hump out of its fluctuations, and either used to be reported
        here as phase separated. Near a critical point the gap closes into the
        fluctuations, and such a system is reported as not phase separated.

        The ``binodal.success`` requirement is load-bearing. Without it, a homogeneous
        (supercritical) system passes: the unguarded tanh fit invents a large density
        gap, and ``condensed_fraction`` is no help either, because at moderate density
        the contact-based clustering *percolates* and reports nearly all the material
        as "condensed" even when there is no condensate at all. A spanning network is
        not a droplet - so the largest cluster must also not span the box on every
        axis in most frames (``percolation_fraction < 0.5``). A slab spans two axes
        but not the third and passes.
        """
        b = self.binodal
        return bool(b.success
                    and np.isfinite(b.rho_dense) and np.isfinite(b.rho_dilute)
                    and b.rho_dense > 2.0 * max(b.rho_dilute, 1e-6)
                    and self.condensed_fraction > 0.3
                    and self.percolation_fraction < 0.5)


def analyze(traj, geometry="auto", min_beads=2, axis=None):
    """Run the full phase-separation analysis and return a
    :class:`PhaseSeparationResult`.

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    geometry : str, optional
        ``'sphere'`` (or its synonym ``'droplet'``), ``'slab'``, or ``'auto'``
        (default), which picks slab when one box axis is at least 1.5 times the
        shortest and spherical otherwise.
    min_beads : int, optional
        Ignore clusters smaller than this many beads throughout (default ``2``).
    axis : int, optional
        Slab geometry only: index of the slab normal, the axis the density
        profile is taken along; negative values count from the last axis.
        Default ``None``: the axis the largest cluster does not span when it
        spans all the others in more than half of the frames (see
        :func:`slab_normal`), and the longest box axis otherwise. Ignored in
        droplet geometry.

    Returns
    -------
    PhaseSeparationResult
        The order parameters, the coexistence fit, the density profile and (in
        droplet geometry) the droplet shape statistics.

    Raises
    ------
    ValueError
        If ``geometry`` is not one of the accepted values, or if ``axis`` is
        not an integer in ``-n_dim .. n_dim - 1``.

    Notes
    -----
    In slab geometry the profile is :func:`slab_density_profile` along the slab
    normal (reported as ``slab_axis``), fit by :func:`fit_slab_profile` with
    ``hardwall=traj.hardwall`` and a dilute plateau at least as wide as a chain
    (``max(2, ceil(2 Rg_z))`` planes). A fit that succeeds must then pass the
    chain-number test: the dense region must hold at least three standard
    deviations more chains than a uniform solution would put there
    (``chain_excess_sigma``), or the fit is reported with ``success=False`` and
    a reason giving the numbers. ``shape`` is ``None``. In droplet geometry
    the profile is :func:`radial_density_profile_with_scatter`, fit by
    :func:`fit_radial_profile` with a site floor of 20 in 3D and 8 in 2D and the
    frame-to-frame scatter; if the
    largest cluster spans the box in more than half the frames the fit is
    reported with ``success=False`` (percentile densities, and a ``reason``
    saying so), and ``shape`` averages only the frames whose largest cluster
    does not span. The warnings the radial profile and :func:`droplet_shape`
    raise for spanning frames are silenced here, because the result reports the
    same thing as ``spanning_fraction``.
    """
    dims = traj.dimensions
    if geometry == "droplet":
        geometry = "sphere"   # synonym (surface_tension's vocabulary)
    if geometry == "auto":
        geometry = "slab" if max(dims) >= 1.5 * min(dims) else "sphere"
    if geometry not in ("slab", "sphere"):
        raise ValueError(
            "analyze: unknown geometry %r (use 'auto', 'slab', 'sphere' or the "
            "synonym 'droplet')" % (geometry,))

    if axis is not None:
        axis = _normalise_axis(axis, traj.n_dim, "analyze")

    cf = condensed_fraction(traj, min_beads=min_beads)
    n_clusters = number_of_clusters(traj, min_beads=min_beads)
    largest = largest_cluster_size(traj, by="beads", min_beads=min_beads)
    # one gather per frame answers all three spanning questions: any axis, every
    # axis (the two fractions below are what spanning_fraction returns), and
    # which axis a slab leaves free
    spans = _largest_cluster_spanning_axes(traj, min_beads=min_beads)
    n_frames = traj.n_frames
    span_any = (sum(1 for s in spans if s) / n_frames) if n_frames else 0.0
    span_all = (sum(1 for s in spans if s is not None and len(s) == traj.n_dim) / n_frames
                if n_frames else 0.0)
    hardwall = bool(getattr(traj, "hardwall", False))

    slab_axis = None
    chain_excess = float("nan")
    if geometry == "slab":
        # The longest box axis used to be taken as the slab normal, always. A
        # slab across a short axis, or in a box with two equal long axes, was
        # then profiled along a direction in which it is flat.
        slab_axis = axis if axis is not None else _vote_slab_normal(
            spans, tuple(traj.dimensions)[:traj.n_dim])
        coord, dens = slab_density_profile(traj, axis=slab_axis)
        # the chain is the unit: the dilute plateau must be wider than a chain,
        # and the slab must hold significantly more chains than a uniform
        # solution would put in the same planes
        chain_extent = _chain_extent_along(traj, slab_axis)
        fit = fit_slab_profile(coord, dens, hardwall=hardwall,
                               min_dilute_planes=max(2, int(np.ceil(2.0 * chain_extent))))
        profile = (coord, dens)
        if fit.success:
            cross_section = int(np.prod([d for i, d in enumerate(traj.dimensions)
                                         if i != slab_axis]))
            chain_excess, n_dense, n_expected = _chain_excess_sigma(
                dens, cross_section, fit.rho_dense, fit.rho_dilute, traj.topology.lengths)
            if chain_excess < _MIN_CHAIN_EXCESS_SIGMA:
                usable = dens[np.isfinite(dens)]
                obs_d, obs_v = _observed_binodal(usable)
                fit = BinodalFit(
                    rho_dense=obs_d, rho_dilute=obs_v, interface_width=fit.interface_width,
                    half_width=fit.half_width, success=False,
                    reason=("the dense region holds %.1f chains where a uniform solution "
                            "would put %.1f: %.1f standard deviations more, and at least "
                            "%.0f are needed - this cannot be told apart from the density "
                            "fluctuation of a one-phase solution (a lump of a few chains is "
                            "not a phase); more chains would help"
                            % (n_dense, n_expected, chain_excess, _MIN_CHAIN_EXCESS_SIGMA)))
    else:
        # spanning_fraction above has already measured how often the largest
        # cluster spans the box, and the block below acts on it (the fit is
        # marked unusable with a reason). The radial profile warns when it
        # leaves those frames out - and with SAVE_EQ on by default, frame 0 of
        # nearly every condensed run is a percolating random placement, so
        # analyze() would warn on essentially every trajectory about a case it
        # handles. That warning is silenced here; the standalone
        # radial_density_profile still raises it. The gather's own percolation
        # warning is silenced as well: the spanning test gathers each frame's
        # largest cluster first with it muted and caches the image, so it should
        # not surface here, but nothing downstream should depend on that.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="single-image gather: cluster percolates")
            warnings.filterwarnings("ignore", message="radial_density_profile:")
            coord, dens, sites, shell_scatter = radial_density_profile_with_scatter(
                traj, min_beads=min_beads)
        min_sites = _min_shell_sites(traj.n_dim)
        fit = fit_radial_profile(coord, dens, site_counts=sites, min_shell_sites=min_sites,
                                 frame_scatter=shell_scatter)
        profile = (coord, dens)
        if span_any > 0.5:
            # the largest cluster is connected to its own periodic image (or, under
            # hardwall, touches both walls) in most frames: it has no single image,
            # so the radial profile is centred on an arbitrary point of a spanning
            # network and cannot describe a droplet, however well a tanh fits it
            usable = np.isfinite(dens) & (np.asarray(sites, float) >= min_sites)
            finite = dens[usable] if usable.any() else dens[np.isfinite(dens)]
            obs_d, obs_v = _observed_binodal(finite) if finite.size else (float("nan"), float("nan"))
            fit = BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                             interface_width=fit.interface_width, radius=fit.radius,
                             success=False,
                             reason=("the largest cluster spans the box in %.0f%% of frames - a "
                                     "spanning network has no droplet profile to fit"
                                     % (100.0 * span_any)))

    # droplet_shape (snakesearch gather + convex hull + sphericity) is only
    # meaningful for a compact droplet: a box-spanning slab percolates through
    # the periodic boundaries and cannot be gathered into a single image (see
    # docs/lemonade/hierarchy.rst), so in slab geometry the shape statistics are
    # not computed rather than silently reporting hull quantities of an
    # artefactual gathering.
    if geometry == "slab":
        shape = None
    else:
        # same reasoning as for the radial profile above: the spanning fraction is
        # already measured and reported on the result, so neither the gather nor
        # droplet_shape need repeat it for the frames the shape statistics leave out
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="single-image gather: cluster percolates")
            warnings.filterwarnings("ignore", message="droplet_shape:")
            shape = droplet_shape(traj, min_beads=min_beads)

    return PhaseSeparationResult(
        geometry=geometry,
        condensed_fraction=float(cf.mean()),
        condensed_fraction_series=cf,
        n_clusters=float(n_clusters.mean()),
        largest_cluster_beads=float(largest.mean()),
        binodal=fit,
        shape=shape,
        profile=profile,
        percolation_fraction=float(span_all),
        spanning_fraction=float(span_any),
        slab_axis=slab_axis,
        chain_excess_sigma=float(chain_excess),
    )
