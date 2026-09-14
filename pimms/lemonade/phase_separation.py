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
* **Slab** - a 1D density profile along the box's long axis (the geometry of the
  ``slab_phase_separation`` demo), slabs re-centred per frame and fit to a
  two-interface ``tanh`` to extract the same coexistence quantities.

Densities are volume fractions (occupied lattice sites per available lattice site),
so ``rho`` runs 0..1 and is directly comparable across box sizes.

Typical use::

    from pimms.lemonade import phase_separation as ps
    result = ps.analyze(traj)          # everything, auto-detecting the geometry
    r, rho = ps.radial_density_profile(traj)
    fit = ps.fit_radial_profile(r, rho)
"""

import warnings
from dataclasses import dataclass, field

import numpy as np


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
    total = traj.n_atoms
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
    centers, density, _sites = _radial_profile_with_site_counts(
        traj, bin_width=bin_width, r_max=r_max, min_beads=min_beads)
    return centers, density


def radial_density_profile_with_site_counts(traj, bin_width=1.0, r_max=None, min_beads=2):
    """:func:`radial_density_profile` plus the mean number of lattice sites per shell.

    The site count is what a fit should weight by: the innermost shell ``[0, 1)``
    holds exactly one site, so its "density" is a single 0/1 draw per frame, while
    a shell at ``r = 10`` averages over ~1200 sites.

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
    return _radial_profile_with_site_counts(
        traj, bin_width=bin_width, r_max=r_max, min_beads=min_beads)


def _radial_profile_with_site_counts(traj, bin_width=1.0, r_max=None, min_beads=2):
    """Bin every bead's distance from the condensate COM into radial shells.

    The single implementation behind :func:`radial_density_profile` and
    :func:`radial_density_profile_with_site_counts`. Bead distances and shell
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
    acc_span = np.zeros(len(centers))
    n_valid_span = np.zeros(len(centers))
    sites_acc_span = np.zeros(len(centers))
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
        com = np.mod(np.round(np.asarray(clusters[0].center_of_mass, dtype=np.float64)), dims[:nd])
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
            n_valid_span[finite] += 1
        else:
            n_droplet += 1
            acc[finite] += contribution[finite]
            sites_acc[finite] += safe_f[finite]
            n_valid[finite] += 1

    if n_spanning and n_droplet == 0:
        acc, sites_acc, n_valid = acc_span, sites_acc_span, n_valid_span
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
    return centers, density, site_counts


def slab_density_profile(traj, axis=None):
    """1D density profile (volume fraction) along ``axis`` (default: the longest
    box axis), with the dense slab re-centred each frame so it does not smear out
    as the slab diffuses. Under periodic boundaries the profile is rolled;
    under a hardwall a free slab is translated (no wrap) so that the centroid of
    its dense bins sits at the window centre, with the vacated bins padded with
    that frame's dilute level, while a condensate touching a wall is left where
    it is (it is pinned, and rolling it would turn its flat wall face into a
    fictitious second interface).

    Every bead in the box is binned, dense and dilute alike - that is what makes the
    result a *density* profile that a coexistence fit can be run against. (It therefore
    takes no ``min_beads``: unlike the radial profile, this one does not need clusters
    at all. It previously accepted a ``min_beads`` argument that was never read.)

    Parameters
    ----------
    traj : pimms.lemonade.LatticeTrajectory
        The trajectory to analyse.
    axis : int, optional
        Index of the axis to profile along (default ``None``, i.e. the longest
        box axis, which is the slab normal in a slab geometry).

    Returns
    -------
    coordinate : numpy.ndarray
        ``(L,)`` float64 lattice coordinate along the axis, ``0 .. L-1``.
    density : numpy.ndarray
        ``(L,)`` float64 frame-averaged occupied fraction at each coordinate.
    """
    dims = traj.dimensions
    hardwall = bool(getattr(traj, "hardwall", False))
    if axis is None:
        axis = int(np.argmax(dims))
    length = dims[axis]
    cross_section = int(np.prod([d for i, d in enumerate(dims) if i != axis]))
    coord = np.arange(length)
    angle = 2.0 * np.pi * coord / length

    acc = np.zeros(length)
    positions = traj.positions
    for f in range(traj.n_frames):
        col = positions[f][:, axis]
        counts = np.bincount(col, minlength=length).astype(np.float64)
        if hardwall:
            # the box is not periodic: a slab cannot wrap, and rolling the profile
            # would move a condensate wetting a wall into the middle of the box,
            # turning its flat wall face into a fictitious second interface. A
            # slab that does NOT touch a wall still diffuses between the walls,
            # so it is aligned by a plain translation of its dense centroid to
            # the window centre (frame-averaging the raw counts smeared a
            # wandering slab into a broad hump with no plateau).
            peak = counts.max()
            dense = counts >= 0.5 * peak if peak > 0 else np.zeros(length, dtype=bool)
            if not dense.any() or dense[0] or dense[-1]:
                acc += counts
                continue
            idx = np.nonzero(dense)[0]
            centroid = float((idx * counts[idx]).sum() / counts[idx].sum())
            shift = int(round(length / 2.0 - centroid))
            if shift == 0:
                acc += counts
                continue
            dilute = float(np.median(counts[~dense])) if (~dense).any() else 0.0
            shifted = np.full(length, dilute)
            if shift > 0:
                shifted[shift:] = counts[:length - shift]
            else:
                shifted[:length + shift] = counts[-shift:]
            acc += shifted
            continue
        # circular centre of mass of the 1D density -> shift dense region to L/2
        cx = (counts * np.cos(angle)).sum()
        cy = (counts * np.sin(angle)).sum()
        com = (np.arctan2(-cy, -cx) + np.pi) / (2.0 * np.pi) * length
        shift = int(round(length / 2.0 - com))
        acc += np.roll(counts, shift)
    return coord.astype(np.float64), acc / (traj.n_frames * cross_section)


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
    fit (see :func:`_fit_is_usable`) is reported with ``success=False``, because
    trusting it silently is the more dangerous failure. When ``success`` is ``False``,
    ``rho_dense`` and ``rho_dilute`` fall back to robust percentiles of the *observed*
    profile, so they stay bounded and physically meaningful (for a homogeneous system
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
       coexistence.

    Attributes
    ----------
    rho_dense : float
        Dense-phase coexistence density, as an occupied-site fraction.
    rho_dilute : float
        Dilute-phase coexistence density, as an occupied-site fraction.
    interface_width : float
        Fitted interface width in lattice units (``nan`` if no fit converged).
    radius : float
        Droplet radius in lattice units; ``nan`` in slab geometry.
    half_width : float
        Slab half-thickness in lattice units; ``nan`` in droplet geometry.
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
                       "model has fitted a dilute slab in a dense background (a condensate "
                       "wetting both walls looks like this; fit it as a wetting film instead)")
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


def fit_radial_profile(radii, density, site_counts=None, min_shell_sites=None):
    """Fit a spherical droplet profile; returns a :class:`BinodalFit`.

    A fit that converges but does not describe a real two-phase profile is returned
    with ``success=False`` - see :class:`BinodalFit`.

    Shells whose density is ``nan`` (no lattice site there - see
    :func:`radial_density_profile`) are ignored. When ``site_counts`` is given (from
    :func:`radial_density_profile_with_site_counts`) shells with fewer than
    ``_MIN_SHELL_SITES`` sites are ignored too: the innermost shell is a single site
    that sits at the largest cluster's own centre, and fitting it with the same
    weight as a 1000-site shell let a homogeneous solution pass as a droplet of
    ``rho_dense = 1`` and radius ``< 1``. A fitted radius below
    ``_MIN_DROPLET_RADIUS`` is rejected for the same reason.

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
        (default ``None``, i.e. ``_MIN_SHELL_SITES``; use
        :func:`_min_shell_sites` for the dimension-aware value).

    Returns
    -------
    BinodalFit
        The fit, with ``success=False`` and a populated ``reason`` when fewer
        than four shells are usable, the optimiser fails, or the converged fit
        is degenerate. In those cases ``rho_dense`` / ``rho_dilute`` fall back to
        percentiles of the observed profile.
    """
    from scipy.optimize import curve_fit
    radii = np.asarray(radii, float)
    density = np.asarray(density, float)
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
    if not ok:
        return BinodalFit(rho_dense=obs_d, rho_dilute=obs_v,
                          interface_width=float(p[3]), radius=float(p[2]),
                          success=False, reason=reason)

    return BinodalFit(rho_dense=float(p[0]), rho_dilute=float(p[1]),
                      radius=float(p[2]), interface_width=float(p[3]))


def _fit_wetting_slab(coord, density):
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

    Parameters
    ----------
    coord : numpy.ndarray
        1D float64 lattice coordinate along the axis.
    density : numpy.ndarray
        1D float64 density profile at those coordinates. Which wall is wetted is
        read off its two ends.

    Returns
    -------
    BinodalFit
        The fit, with ``half_width`` set to half the fitted film thickness and
        ``radius`` left at ``nan``. ``success`` and ``reason`` are inherited from
        the underlying single-interface fit.
    """
    # The +0.5 puts the origin ON THE WALL. Without it the distance was measured
    # from the CENTRE of the first site, half a lattice unit inside the wall, and
    # every wetting film came back exactly 0.5 units too thin - a 5% error at 10
    # planes but 17% at 3, systematic and one-signed, and inconsistent with the
    # two-interface fit, which reports 10.00 for the same 10-plane film.
    if density[0] >= density[-1]:
        z = coord - coord[0] + 0.5           # wets the low wall: profile falls with z
    else:
        z = coord[-1] - coord + 0.5          # wets the high wall: mirror it
    order = np.argsort(z)
    fit = fit_radial_profile(z[order], density[order])
    return BinodalFit(rho_dense=fit.rho_dense, rho_dilute=fit.rho_dilute,
                      interface_width=fit.interface_width,
                      half_width=0.5 * fit.radius, success=fit.success,
                      reason=fit.reason)


def fit_slab_profile(coord, density, hardwall=None):
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
    interfaces and must not be routed to the single-interface fit).

    Parameters
    ----------
    coord : array_like
        1D lattice coordinate along the slab normal, from
        :func:`slab_density_profile`.
    density : array_like
        1D density profile at those coordinates.
    hardwall : bool, optional
        Whether the box has hard walls. ``False`` fixes the slab centre at the
        middle of the window (a periodic profile is re-centred per frame) and
        never uses the wetting model; ``True`` and the default ``None`` both fit
        the centre as a free parameter and pick the wetting model from the
        profile shape.

    Returns
    -------
    BinodalFit
        The fit, with ``success=False`` and a populated ``reason`` when the
        optimiser fails or the converged fit is degenerate (including a slab
        that fills the box, leaving no dilute phase).
    """
    from scipy.optimize import curve_fit
    coord = np.asarray(coord, float)
    density = np.asarray(density, float)
    if hardwall is not False and density.size >= 4:
        peak = float(density.max())
        mid = float(density[density.size // 4: -max(1, density.size // 4)].mean())
        end_dense = peak > 0 and (density[0] >= 0.5 * peak or density[-1] >= 0.5 * peak)
        # the ends must also be denser than the middle: a re-centred periodic
        # slab, or a hardwall slab away from the walls whose dilute phase is at
        # least half the dense density, has two interfaces and stays on the
        # two-interface fit whatever the hardwall flag says
        if end_dense and max(density[0], density[-1]) > mid:
            return _fit_wetting_slab(coord, density)
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
        center = float(coord[0] + length / 2.0)

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
        frames that contain one whose single image is unambiguous (``nan`` if
        none do). Degenerate convex-hull values (``-1``) are left out of the
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
        Trajectory-averaged number of clusters.
    largest_cluster_beads : float
        Trajectory-averaged bead count of the largest cluster.
    binodal : BinodalFit
        The coexistence fit to the density profile.
    shape : dict or None
        :func:`droplet_shape` statistics over the frames whose largest cluster
        has a single image, or ``None`` in slab geometry, where hull quantities
        of a box-spanning slab are meaningless.
    profile : tuple
        The fitted density profile as ``(coordinate, density)``: shell radii in
        droplet geometry, the axis coordinate in slab geometry.
    percolation_fraction : float
        Fraction of frames whose largest cluster spans the box on EVERY axis (a
        space-filling network: the contact clustering of a homogeneous
        solution).
    spanning_fraction : float
        Fraction of frames whose largest cluster spans the box on ANY axis, and
        so has no single periodic image.
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
    #: fraction of frames in which it spans the box on ANY axis (no single image)
    spanning_fraction: float = 0.0

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


def analyze(traj, geometry="auto", min_beads=2):
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

    Returns
    -------
    PhaseSeparationResult
        The order parameters, the coexistence fit, the density profile and (in
        droplet geometry) the droplet shape statistics.

    Raises
    ------
    ValueError
        If ``geometry`` is not one of the accepted values.
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

    cf = condensed_fraction(traj, min_beads=min_beads)
    n_clusters = number_of_clusters(traj, min_beads=min_beads)
    largest = largest_cluster_size(traj, by="beads", min_beads=min_beads)
    span_any = spanning_fraction(traj, min_beads=min_beads, all_axes=False)
    span_all = spanning_fraction(traj, min_beads=min_beads, all_axes=True)
    hardwall = bool(getattr(traj, "hardwall", False))

    if geometry == "slab":
        coord, dens = slab_density_profile(traj)
        fit = fit_slab_profile(coord, dens, hardwall=hardwall)
        profile = (coord, dens)
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
            coord, dens, sites = radial_density_profile_with_site_counts(traj, min_beads=min_beads)
        min_sites = _min_shell_sites(traj.n_dim)
        fit = fit_radial_profile(coord, dens, site_counts=sites, min_shell_sites=min_sites)
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
    )
