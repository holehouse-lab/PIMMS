.. _lemonade-phase-separation:

===================================
Phase separation & droplet physics
===================================

The ``pimms.lemonade.phase_separation`` and ``pimms.lemonade.surface_tension``
modules quantify liquid-liquid phase separation of a PIMMS system: the coexistence
(binodal) densities, the amount of condensed material, cluster-size order
parameters, droplet shape and the interfacial tension.

.. code-block:: python

   import pimms.lemonade as lemonade
   from pimms.lemonade import phase_separation as ps
   from pimms.lemonade import surface_tension as st

   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

Densities throughout are **occupied lattice-site fractions** in ``[0, 1]``, and
so are directly comparable across box sizes. ``shape['density']`` is the
exception: beads per convex-hull volume, which can exceed 1. In **slab geometry**
the ``shape`` statistics are not computed at all (``result.shape`` is ``None``),
because the convex hull of a box-spanning, percolating slab is meaningless.

Everything in ``phase_separation`` works in 2D and 3D. The two surface-tension
estimators are **3D only** and raise ``ValueError`` on a 2D trajectory.

Clusters and condensates
========================

Clusters are connected groups of chains, found per frame by
:attr:`Frame.clusters <pimms.lemonade.Frame.clusters>`, ordered largest first **by bead
count** (which is not the same as chain count once chains differ in length). The
largest cluster - the condensate - is ``frame.droplet``. Each cluster carries its
own geometry (:doc:`hierarchy`): ``radius_of_gyration``, ``asphericity``,
``sphericity``, ``volume``, ``surface_area``, ``density``,
``radial_density_profile()`` and its composition.

Order parameters
================

Simple, per-frame measures of *how* phase separated the system is:

.. code-block:: python

   ps.condensed_fraction(traj)          # (n_frames,) fraction of beads in the largest cluster
   ps.largest_cluster_size(traj)        # (n_frames,) beads in the largest cluster
   ps.largest_cluster_size(traj, by="chains")   # ... or chains
   ps.number_of_clusters(traj, min_beads=2)
   ps.cluster_size_distribution(traj)   # all cluster sizes pooled across frames
   ps.spanning_fraction(traj)           # fraction of frames whose largest cluster
                                        # spans the box on ANY axis (see below)
   ps.spanning_fraction(traj, all_axes=True)   # ... on EVERY axis
   ps.droplet_shape(traj)               # frame-averaged largest-cluster geometry

``condensed_fraction`` is the basic order parameter: near zero when the system is
well mixed, approaching one when most material collects into a single condensate.
Every one of these takes ``min_beads=`` to ignore clusters below a size (``1``
except for ``number_of_clusters``, ``spanning_fraction`` and ``droplet_shape``,
which default to ``2``), and ``by=`` must be ``'beads'`` or ``'chains'`` - any
other string raises rather than quietly running a different measurement.

A cluster *spans* an axis when, under periodic boundaries, it is connected to its
own periodic image through that face - a pair of its beads touch through it - and,
under a hardwall, when it touches both walls of that axis. Merely reaching the box
length is not enough: a contact staircase from one corner of the box to the other
does that on every axis without winding the box, and has a perfectly good single
image.
``droplet_shape`` returns the frame-averaged ``radius_of_gyration``,
``asphericity``, ``sphericity``, ``volume`` and ``density`` of the largest
cluster, leaving the degenerate ``-1`` convex-hull values out of the averages. A
frame whose largest cluster spans the box is left out too, with a warning saying
how many were: such a cluster is a network or a slab rather than a droplet. Under
periodic boundaries it has no single image, and the shape of the window the
gather hands back describes the search rather than the cluster; under a hardwall
it touches both walls, and its shape is that of a wall-bounded condensate (frame
0 of a run that saved its equilibration is the usual one, a random placement that
percolates as a contact network). ``droplet_shape`` assumes a compact droplet,
which is why ``analyze`` does not report it in slab geometry - call it directly
on a slab only if you know what you are asking for.

The binodal (coexistence densities)
===================================

The dense- and dilute-phase densities are read off a density profile fit to a
hyperbolic tangent. Two geometries are supported.

**Droplet (spherical).** A radial profile about the condensate centre of mass -
occupied fraction as a function of distance - falls from the dense core to the
dilute background:

.. code-block:: python

   r, rho, sites, scatter = ps.radial_density_profile_with_scatter(traj)
   fit = ps.fit_radial_profile(r, rho, site_counts=sites,    # sites: thin shells are ignored
                               frame_scatter=scatter)        # scatter: the single-frame test
   r, rho, sites = ps.radial_density_profile_with_site_counts(traj)   # without the scatter
   r, rho = ps.radial_density_profile(traj)                  # the profile alone, for plotting
   fit.rho_dense, fit.rho_dilute        # coexistence densities
   fit.radius, fit.interface_width      # droplet radius and interface width

The fit is
:math:`\rho(r) = \tfrac12(\rho_d + \rho_v) - \tfrac12(\rho_d - \rho_v)\tanh((r-R)/w)`.

Shells are one lattice unit wide by default (``bin_width=``). The profile runs
out to half the shortest box axis under periodic boundaries - beyond that the
minimum image folds back on itself - and to the full box diagonal under a
hardwall, where distances are plain Cartesian; pass ``r_max=`` to override.
Frames whose largest cluster spans the box are left out of the average, with a
warning saying how many: a spanning cluster is a network or a slab, so a profile
about its centre is not a droplet profile. Under periodic boundaries it has no
single image and that centre is an arbitrary point of the gathered window; under
a hardwall the centre of mass is exact and the profile is well defined, but it
describes a wall-bounded condensate. If *every* frame spans there is no droplet
to profile, and the profile about that centre is returned instead (with a
different warning) because its percentiles still estimate the two densities.

**Slab.** For a slab condensate that spans the periodic plane and is bounded along
one axis (the geometry of the ``slab_phase_separation`` demo), a 1D profile along
the slab normal - with the slab re-centred each frame so it does not smear - is fit to
a two-interface tanh:

.. code-block:: python

   z, rho, scatter = ps.slab_density_profile_with_scatter(traj)   # axis defaults to the slab normal
   fit = ps.fit_slab_profile(z, rho, hardwall=traj.hardwall, frame_scatter=scatter)
   fit.rho_dense, fit.rho_dilute, fit.interface_width
   fit.half_width                       # half the slab thickness

   z, rho = ps.slab_density_profile(traj)           # the same profile without the scatter
   ps.slab_normal(traj)                             # the axis both of them pick

``scatter`` is the sample standard deviation, over frames, of each plane's density
after alignment. ``frame_scatter=`` is an optional, strict test (see
:ref:`lemonade-binodal-success`): it compares single planes, whose scatter grows
as the cross-section shrinks, so it rejects real slabs in narrow boxes and
``analyze()`` does not use it in slab geometry. Frames saved every step or two
are strongly correlated and make the scatter read low.

``axis=`` picks the slab normal explicitly (``0``, ``1`` or, in 3D, ``2``;
negative values count from the last axis, and anything else raises). Left at
``None`` the normal is chosen by
:func:`~pimms.lemonade.phase_separation.slab_normal`: a slab spans every box axis
but one, so when the largest cluster spans all the axes except the same one in
more than half of the frames that hold a cluster, that one is the normal.
Otherwise the clusters do not single an axis out (a droplet, a network, a
one-phase solution) and the longest box axis is used, the first of them on a tie.
The longest axis alone used to be used, which profiled a slab lying across a short
axis, or in a box with two equal long axes, along one of its own in-plane
directions, where it is flat.

Passing ``hardwall=`` tells the fit which model to use: ``False`` pins the slab
centre at the middle of the window, which is right for a periodic profile because
:func:`~pimms.lemonade.phase_separation.slab_density_profile` re-centres it every
frame. ``True`` and the default ``None`` both fit the centre as a free fifth
parameter. ``analyze()`` passes ``traj.hardwall`` for you.

Under a hardwall the profile cannot be rolled. A slab that is free of both walls
is translated so its dense planes sit at the window centre; the planes a frame
vacates carry no data for that frame, so each plane is averaged over the frames
that actually cover it, and a plane that no frame covers is ``nan`` (the fits
skip it). A condensate that touches a wall is left where it is, because moving it
would turn its flat wall face into a second interface, and the fit then uses a
single-interface wetting model. A wetting film and a free slab are different
profiles, so when the condensate touches a wall in some frames and not in others
only the majority kind is averaged (a tie goes to the free slab) and a warning
says how many frames were left out.

A one-phase hardwall trajectory is neither translated nor split. A frame holds a
slab when its dense planes (at least half the frame's peak count) form one block
with a dilute region at least four planes wide on one side, and a trajectory is
treated as a slab trajectory only when at least half of its frames do. A
one-phase solution fails that test - at low density its dense planes are
scattered, and at higher density they fill the box apart from the depleted plane
or two next to each wall - so every frame is averaged exactly where it is and the
profile is the plain per-plane mean, depletion layers included. (The block test
alone used to pass such a box-filling frame, and each frame was then translated
by the centroid of its noise: the wall plane of an athermal solution read 0.041
where the plain mean is 0.010.)

Every bead in the box is binned, dense and dilute alike - that is what makes the
result a density profile a coexistence fit can be run against - so unlike the
radial profile the profile itself needs no clusters and takes no ``min_beads``
(only the default choice of ``axis`` looks at clusters).

.. _lemonade-binodal-success:

Always check ``fit.success``
============================

**A one-phase system does not make the fit fail - it makes it lie.** Above the
critical temperature the density profile is flat, and a ``tanh`` asked to find two
interfaces in a flat line does not raise: it converges to a very wide ``tanh``, which
over a finite box is almost a straight line, and then parks ``rho_dense`` and
``rho_dilute`` at whatever the data does not constrain (usually the bounds, 1 and 0).
The fit reports a large, entirely fictitious coexistence gap, with every appearance of
having succeeded.

Both fits therefore validate themselves, and set ``success = False`` when the fit is
**not well posed** - when the data does not constrain the parameters:

* the fit is **inverted** - the fitted dense density comes out *below* the dilute
  one, which means the model has fitted a dilute slab in a dense background. The
  densities are not silently swapped; or
* the density gap is **absent** - below ``1e-3`` in occupied fraction, i.e. the
  profile is homogeneous; or
* the fitted profile never actually **reaches its own asymptotes** inside the box, so
  the reported coexistence densities are extrapolation; or
* the density gap is **the size of the scatter** in the profile - noise, not signal; or
* the slab **fills the box**, leaving no dilute phase for the dense phase to coexist
  with; or
* (slab) the slab is **under two planes thick** - the fitted profile is above the
  midpoint of the two densities on a single plane, whose density and thickness
  cannot be told apart; or
* (slab) the slab is **thinner than its interfaces** - the fitted profile rises
  less than three quarters of the way from ``rho_dilute`` to ``rho_dense``
  (``half_width < interface_width`` for a free slab), so ``rho_dense`` is an
  extrapolation. A thicker slab (more material) is needed; or
* (slab) there is **no dilute plateau** - the fitted profile comes within 10% of
  the gap of ``rho_dilute`` on fewer than three planes (``min_dilute_planes=``;
  ``analyze()`` asks instead for a plateau wider than a chain,
  ``max(2, ceil(2 Rg_z))`` planes, with ``Rg_z`` the chains' radius of gyration
  along the normal). This is what the depletion layer of a one-phase solution next
  to a hard wall looks like: a ramp about as wide as a chain, which the ``tanh``
  fits as the edge of a slab filling the box, with a "dilute density" below
  anything observed; or
* (slab, in :func:`~pimms.lemonade.phase_separation.analyze`) the slab does not
  hold **significantly more chains than a uniform solution would** - see below; or
* (optional, when the fit is given ``frame_scatter``) the gap is **not resolved in
  single frames** - it is less than three standard deviations of the difference
  between one dense plane (or core shell) and one dilute plane (or outer shell) in
  a single frame. ``analyze()`` applies this in droplet geometry only; or
* (droplet) the fitted profile **ends before its dilute plateau** - at the
  outermost usable shell it is still more than 10% of the gap above
  ``rho_dilute``, which is then an extrapolation (the profile stops at half the
  shortest box axis; use a larger box); or
* (droplet) **fewer than four shells are usable** after the filters below, or
  (slab) fewer than four planes hold data; or
* (droplet) the fitted **radius is below two lattice units** - the fit has latched
  onto the handful of sites at the cluster centre rather than a dense phase; or
* (droplet, via :func:`~pimms.lemonade.phase_separation.analyze`) the largest cluster
  **spans the box** in more than half the frames - under periodic boundaries it is
  connected to its own periodic image, has no single image, and its "radial
  profile" is centred on an arbitrary point of a network; under a hardwall it
  touches both walls of an axis and is a wall-bounded film or network, not a
  droplet.

A fit whose optimiser did not converge at all is reported the same way, with
``reason = "curve_fit did not converge"``.

The radial fit also ignores shells that hold too few lattice sites to constrain
anything, and shells that hold none at all (returned as ``nan`` by the profile -
under ``HARDWALL`` the outer shells are cut off by the walls, and "no site" is not
"empty"). The innermost shell is a single site at the largest cluster's own centre
of mass and reads ~0.7 even for a homogeneous solution; given equal weight it used
to turn a one-phase system into a "droplet" of ``rho_dense = 1`` and radius ``< 1``,
with ``success = True``.

The site floor is dimension-aware: 20 sites in 3D, 8 in 2D, whose shells hold only
~2 pi r sites and would otherwise lose every shell inside ``r = 3``.
:func:`~pimms.lemonade.phase_separation.analyze` selects the right one from
``traj.n_dim``. Calling
:func:`~pimms.lemonade.phase_separation.fit_radial_profile` yourself applies the
3D floor of 20 whatever the dimensionality, so for a 2D profile pass
``min_shell_sites=8``. The counts themselves come from
:func:`~pimms.lemonade.phase_separation.radial_density_profile_with_site_counts`;
without a ``site_counts=`` argument no floor is applied at all.

When ``success`` is ``False``, ``reason`` names the check that failed, and
``rho_dense`` / ``rho_dilute`` fall back to the 95th and 5th percentiles of the
usable part of the observed profile - so they stay bounded and, for a homogeneous
system, simply coincide.

.. code-block:: python

   fit = ps.fit_slab_profile(z, rho)
   if not fit.success:
       print(f"fit is not well posed: {fit.reason}")

.. important::

   ``success`` is a **numerical** guarantee, not a physical one. It says the fit is
   meaningful *as a fit*. It does **not** tell you the system is phase separated.

   Two reasons. First,
   :func:`~pimms.lemonade.phase_separation.slab_density_profile` re-centres the
   slab every frame,
   which aligns the fluctuations of even a *homogeneous* system into a shallow central
   hump - and a ``tanh`` fits that hump perfectly well, giving a small but well-posed
   density gap. The radial profile does the same by centring every frame on the
   largest cluster, which in a dilute solution is a single coil. The averaged
   profile alone cannot tell that hump from a thin slab or a small droplet, and a
   fit called by hand has not been through the tests ``analyze()`` adds (the
   number of chains in the slab; the frame-to-frame scatter of a droplet's
   shells). Second, and
   more fundamentally, the coexistence gap **closes continuously** as the critical
   point is approached, so there is no numerical criterion that can draw the line
   for you: close to the critical point the gap sinks into the fluctuations and the
   single-frame test reports the system as not resolved, which is the cautious
   answer and not a measurement of the critical temperature.

   For the physical question use :attr:`~pimms.lemonade.phase_separation.PhaseSeparationResult.is_phase_separated`,
   or apply your own density-contrast threshold. Mapping a binodal across temperature
   needs both: ``fit.success`` to throw out the degenerate fits, and a contrast
   threshold to decide which of the survivors are really two-phase.

One call for everything
=======================

:func:`~pimms.lemonade.phase_separation.analyze` runs the whole pipeline and
auto-detects the geometry (slab when the longest box axis is at least 1.5 times
the shortest, else spherical):

.. code-block:: python

   result = ps.analyze(traj)                     # or geometry='slab' / 'sphere', min_beads=..., axis=...

   result.geometry                      # 'sphere' or 'slab' (never 'auto')
   result.rho_dense, result.rho_dilute  # binodal (shortcuts into result.binodal)
   result.binodal                       # the BinodalFit itself: .success, .reason, ...
   result.binodal.interface_width
   result.condensed_fraction            # time-averaged
   result.condensed_fraction_series     # (n_frames,) the per-frame values behind it
   result.n_clusters                    # time-averaged cluster count
   result.largest_cluster_beads         # time-averaged size of the largest cluster
   result.shape                         # {'radius_of_gyration', 'asphericity', 'sphericity',
                                        #  'volume', 'density'} - None in slab geometry
   result.is_phase_separated            # usable fit AND density gap AND most material condensed AND not a box-filling network
   result.percolation_fraction          # frames in which the largest cluster spans EVERY box axis
   result.spanning_fraction             # frames in which it spans ANY axis (not a droplet)
   result.profile                       # (coordinate, density) for plotting
   result.slab_axis                     # slab geometry: the axis profiled along; None otherwise
   result.chain_excess_sigma            # slab geometry: chains in the slab against a uniform solution

``geometry`` accepts ``'auto'`` (default), ``'slab'``, ``'sphere'`` and the
synonym ``'droplet'``; anything else raises. ``min_beads`` (default ``2``) is
applied to every cluster-based step. ``axis=`` sets the slab normal in slab
geometry (default: the axis :func:`~pimms.lemonade.phase_separation.slab_normal`
picks, reported back as ``slab_axis``) and is ignored in droplet geometry. In
slab geometry ``result.chain_excess_sigma`` reports the chain-number test described
below. The exact rule behind
``is_phase_separated`` is ``binodal.success`` **and** both densities finite
**and** ``rho_dense > 2 * max(rho_dilute, 1e-6)`` **and**
``condensed_fraction > 0.3`` **and** ``percolation_fraction < 0.5``.

``percolation_fraction`` is the guard against the classic false positive: at
moderate volume fraction the contact clustering of a perfectly homogeneous
solution *percolates*, so ``condensed_fraction`` is ~1 and the largest cluster
holds nearly every bead - yet there is no condensate. A cluster spanning every axis
of the box in most frames is a network, not a droplet, and
``is_phase_separated`` is ``False`` for it. A slab spans two axes and passes.

**One-phase solutions in an elongated box.** Coexistence is a statement about
many chains, and the independent unit of a polymer solution's density
fluctuations is the chain, not the bead. In a one-phase solution of ``N`` chains
each chain sits in the dense region (a fraction ``f`` of the planes) with
probability ``f``, so the region holds ``N f`` chains give or take
``sqrt(N f (1 - f))``. Re-centring every frame on its densest region selects the
upward fluctuation - that is why the averaged profile of a one-phase solution
shows a hump - but only by a standard deviation or two. ``analyze()`` therefore
asks, in slab geometry, that the dense region (profile above the midpoint of the
two fitted densities) hold at least **three standard deviations** more chains
than ``N f``, counting beads in units of the bead-weighted mean chain length, and
that the dilute plateau be **wider than a chain**. A fit that fails is reported
with ``success = False`` and a reason giving the numbers; ``chain_excess_sigma``
is on the result either way. A lump of two or three long coils is not a phase,
however clean its averaged profile, and a system of a handful of chains cannot
pass.

What this has been measured on: 352 real PIMMS runs in slab geometry, from three
independent sets (the auditor's, the fixer's and a reviewer's). The rules were
fixed by physical argument and the runs were split in two by a hash of their
names; the table gives both halves, as runs reported phase separated.

.. list-table::
   :header-rows: 1

   * - Class
     - n
     - before
     - now
   * - One-phase (athermal, or far above the critical temperature), first half
     - 104
     - 52
     - 0
   * - One-phase, second half
     - 98
     - 44
     - 0
   * - Two-phase slabs and stripes, first half
     - 41
     - 41
     - 41
   * - Two-phase slabs and stripes, second half
     - 34
     - 33
     - 33

Among fits that pass every other check the chain excess is at most 2.2 standard
deviations for the one-phase runs and at least 5.9 for those two-phase runs.
Three classes are **not** settled, and a verdict on them should not be trusted:

* **A narrow hardwall box close to the critical temperature** (8 x 8 x 40 at
  T = 130, whose periodic twin holds a slab of density 0.78): the condensate
  breaks into several lumps that wander between the walls, no frame shows a
  slab, and the slab profile has no plateau. 2 of 4 such runs are reported phase
  separated (with a "dense" density of 0.3 that is not a coexistence density)
  and 2 are not.
* **A hardwall condensate that fills the box** (12 x 12 x 24 with 450 to 600
  chains): its largest cluster touches every wall, so it is reported as a
  network, not phase separated (0 of 13, as before).
* **Near the critical point** (26 runs at the last temperatures at which a
  stripe or slab is still visible) the verdict varies from seed to seed: 13 of
  26, against 15 before.

Three further limits follow from the test itself, and an independent check on
runs outside the table above found each of them:

* **A system of only a few chains.** The chain-number test counts chains, so a
  fully condensed slab of ``N`` chains occupying a fraction ``f`` of the planes
  can score at most ``sqrt(N (1 - f) / f)`` standard deviations: it needs more
  than 9 chains when the slab fills half the box and more than 36 when it fills
  80% of it. A slab of twelve 100-bead chains (density gap 0.96) is rejected at
  2.9 standard deviations, and a single long chain collapsed among monomers is
  rejected too. Both are condensed, but neither can be told from a fluctuation
  by counting chains; use more chains or a longer box.
* **Short chains between hard walls.** For chains of four to six beads the
  required plateau is only two planes wide, and the wall depletion layer can
  then pass the slab fit: ``result.binodal.success`` is occasionally ``True``
  for an athermal solution, with a "dilute" density well below the bulk. In
  every such run seen ``is_phase_separated`` was still ``False`` (the cluster
  percolates), so read the verdict from ``is_phase_separated``, never from
  ``binodal.success`` alone.
* **A layer adsorbed on a wall.** A few dense planes against a hard wall over a
  uniform dilute bulk is reported as a wetting film. The density profile cannot
  distinguish adsorption from a thin wetting phase; that is a question about
  the thermodynamics, not about the profile.

The verdict needs frames that sample the equilibrium state. Leave equilibration
frames out. A handful of frames, or frames saved every step or two (which are
strongly correlated), give a profile that is one fluctuation rather than an
average, and a one-phase run can then pass; ``analyze()`` does not refuse a short
trajectory, so check ``traj.n_frames`` yourself.

In droplet geometry the test is the frame-to-frame scatter of the radial shells
(a coil's density fluctuates by about as much as its hump is tall; a droplet's
core reads the same in every frame). Profiled about their largest cluster, 5 of
54 athermal runs used to be reported as a "droplet" of radius 3 and density 0.15
to 0.2 - a single coil - and none is now; droplet geometry has been checked on far
fewer runs than slab geometry.

**Hardwall boxes.** Under ``HARDWALL`` a profile is never rolled (it cannot
wrap). A free slab that clears both walls is aligned frame by frame by a plain
translation of its dense centroid to the window centre (the planes it vacates get
no data from that frame, as described in the slab section above), so a slab
diffusing between the walls is not smeared into a broad hump; a condensate
touching a wall is left where it is and fit with a single-interface model, since
its wall face is not an interface -
``2 * half_width`` is then the slab thickness measured from the wall. The wall is
the *face* of the first lattice plane, half a lattice unit outside its centre, so
a film of ``t`` fully occupied planes reports a thickness of ``t`` - the same
number the two-interface fit gives for a free slab of ``t`` planes. The two
geometries are therefore directly comparable, which matters if you are plotting
film thickness against temperature or concentration.

The single-interface (wetting) model is chosen from the *shape* of the profile,
not from the flag: an end of the window must be at least half the peak density
**and** denser than the middle of the window. Both conditions matter. A hardwall
slab sitting away from the walls, whose dilute phase happens to be more than half
the dense density, still has two interfaces and stays on the two-interface fit.
``hardwall=False`` never uses the wetting model at all.

A condensate wetting **both** walls - the profile is dense at both ends, each at
least half the peak and denser than the middle - is two films with the dilute
phase between them. It is fit as a dilute slab in a dense background (the
two-interface form with the two densities exchanged), ``rho_dense`` is the
density of the films, and ``half_width`` is half the *mean* thickness of the two
films, each measured from its wall face.

Films are checked as slabs, by counting planes: a film of one plane is rejected
(its density and thickness cannot be told apart) with a reason that says so, and
a film of two planes is accepted at either wall. (The wetting fit used to borrow
the droplet check on the fitted radius, and since the fitted thickness of a sharp
two-plane film can land anywhere between 1.5 and 2.5, a film and its mirror image
could get opposite verdicts.)
Because only a periodic profile is guaranteed to be centred on the window, the
two-interface fit treats the slab centre as a fitted fifth parameter unless you
pass ``hardwall=False``.

The spanning guard is stated for periodic boxes, where it is PIMMS's own
percolation test on the gathered cluster: reaching the box length on an axis is
necessary for a cluster to be connected to its own image but not sufficient (a
contact staircase from corner to corner reaches the box length without any pair
of beads meeting through the face), so the axis is confirmed by a pair that
does. Under ``HARDWALL`` a cluster that touches both walls is still counted as
spanning: its single image is unambiguous, but it is a wall-bounded film or
network rather than a droplet, so it has no droplet geometry to fit, and the
radial profile and shape averages leave those frames out just as they leave out
periodic spanning frames. The per-axis answer is cached on the cluster as
``Cluster.spanning_axes()``, computed once per gather.

Surface tension from undulations
================================

Both surface-tension estimators use capillary-wave theory - the interface's
fluctuation spectrum. Because PIMMS uses :math:`\exp(-\Delta E/T)` with
:math:`k_B = 1`, the trajectory temperature *is* :math:`k_B T`, and the returned
:math:`\gamma` is in **reduced units** (interaction energy per lattice area).
Both are **3D only** and raise ``ValueError`` on a 2D trajectory.

.. code-block:: python

   st.surface_tension(traj)             # auto-dispatch by geometry
   st.slab_surface_tension(traj)        # planar capillary waves  (robust)
   st.droplet_surface_tension(traj)     # spherical-harmonic shape fluctuations

``surface_tension`` dispatches on ``geometry=``: ``'auto'`` (default) picks slab
when the longest box axis is at least 1.5 times the shortest and droplet
otherwise, and ``'slab'``, ``'droplet'`` and the synonym ``'sphere'`` force the
choice. Anything else raises. Remaining keyword arguments go straight through to
the chosen estimator.

* **Slab** - the two flat interfaces of a box-spanning condensate have a height
  field :math:`h(x,y)` obeying :math:`\langle|h(q)|^2\rangle = k_BT/(\gamma A q^2)`;
  averaging the low-:math:`q` spectrum gives :math:`\gamma`. On the lattice the
  estimator replaces the continuum :math:`q^2` with the exact lattice dispersion
  :math:`(2-2\cos q_x) + (2-2\cos q_y)` - identical in the continuum limit, but
  the continuum form under-estimates :math:`\gamma` by about 2% on a 20 x 20
  cross-section and 12% on an 8 x 8 one (with the default 8 modes).
  This is the reliable method. ``axis=`` sets the slab normal (negative values
  count from the last axis; default: the axis the largest cluster leaves
  unspanned in most frames, as for the density profile, and the longest box axis
  otherwise), ``n_modes=`` how many independent low-:math:`q` modes to
  average over (default 8) and ``min_beads=`` the smallest cluster that can be
  the condensate (default 2). Only frames whose largest cluster *is* a slab are
  used: it must span both in-plane axes and must not span the normal (see
  ``spanning_fraction`` above for what spanning means). A cluster that also spans
  the normal is a network - the contact clustering of a homogeneous solution at
  moderate volume fraction looks like this - and one that misses an in-plane axis
  is a droplet or a strip; both are skipped with a warning giving the count, and
  ``gamma`` is ``nan`` if no frame is left.
* **Droplet** - the radius :math:`R(\theta,\phi)` fluctuates in spherical-harmonic
  modes with :math:`\langle|u_{lm}|^2\rangle = k_BT/(\gamma R_0^2 (l-1)(l+2))` for
  :math:`l \ge 2` (up to ``l_max=``, default 5). Best-effort: it needs a single,
  compact, reasonably large droplet sampled over many frames, and it ignores
  clusters below ``min_beads=`` (default 30 here, not 2). Even then it reads
  **low** and depends on the angular grid (see below), so treat a droplet
  estimate as a rough number and prefer the slab estimator whenever the geometry
  allows it.

The droplet estimator takes the interface radius in each angular bin to be that
of the outermost bead in it. The grid is sized from the droplet itself, aiming
at about one surface bead per bin: ``n_polar = 2 R0`` rounded and clamped to
the range 8 to 64, and ``n_azim = 2 n_polar``, with :math:`R_0 = (3n/4\pi)^{1/3}`
and :math:`n` the **median** bead count of the largest cluster over every frame
that holds one of at least ``min_beads`` beads. Sizing it from the first frame
that held a cluster instead let a not-yet-condensed leading frame (frame 0 is
the start configuration, and ``SAVE_EQ`` keeps the equilibration frames) choose
a grid far too fine for the real droplets, which the per-frame coverage test
then rejected one by one, silently. Any frame whose droplet still fills fewer
than half the bins is skipped, and a warning names how many. Pass ``n_polar=``
/ ``n_azim=`` to size the grid yourself; either way the grid used is reported on
the result.

We measured the droplet estimator on lattice droplets filled from a known
capillary spectrum (modes :math:`l = 2` to :math:`1.5 R_0`, the droplet centre at
random positions relative to the lattice, 200 frames and 10 seeds per point,
15 at :math:`R_0 = 18`). The automatic grid read :math:`0.89\gamma`,
:math:`0.87\gamma` and :math:`0.83\gamma` for :math:`R_0 = 8`, 12 and 18 at :math:`\gamma = k_BT`, and
:math:`0.90\gamma` and :math:`0.84\gamma` at :math:`\gamma = 0.5` and
:math:`1.5\,k_BT` with :math:`R_0 = 12`, with a scatter of 2 to 3 % from one
seed to the next. The bias comes from the grid: near its poles the azimuthal
bins are narrower than a lattice site, so many of them hold no surface bead and
report an inner one. The radius there reads short in every frame, and because
that deficit is symmetric about the grid axis it lands in the even modes
(:math:`l = 2` and 4) as apparent fluctuation - which matters more the smaller
the true fluctuations are, so the estimate reads lower for a stiffer interface.
A droplet centred exactly on a lattice site, as in the oracle of the test
suite, is a special case that reads :math:`0.95` to :math:`0.98\gamma`. No other grid is reliably better: :math:`8 \times 16` read
:math:`0.97` to :math:`0.99\gamma` on these droplets but :math:`1.14` to
:math:`1.37\gamma` on smooth ones that carry only the fitted modes
:math:`l \le 5`, and grids twice as fine as the automatic one read anywhere from
:math:`0.5` to :math:`1.7\gamma`. Real PIMMS droplets are further off. The chains
of the ``slab_phase_separation`` demo at the same temperature, run as droplets
in cubic boxes, gave :math:`\gamma` = 36.5 and 37.4 (two seeds,
:math:`R_0 \approx 9`) and 37.9 (:math:`R_0 \approx 11`) on the automatic grid,
about three quarters of the 49.0 and 49.4 that the slab estimator gives for two
seeds of the demo itself (frames 20 onwards throughout), and the same droplet
frames gave 25 to 71 on other grids.

Under ``HARDWALL`` the slab estimator uses only faces that are not pressed against
a wall: a wall face is flat because the wall is, and counting its zero capillary
power halved the spectrum and doubled :math:`\gamma` for a wetting condensate.
Under ``HARDWALL`` the in-plane axes are walls too, so the height field is not
periodic in the plane and its capillary modes are not plane waves. They are the
modes of the lattice Laplacian with free ends,
:math:`\cos(\pi m (x + 1/2)/L_x)\,\cos(\pi n (y + 1/2)/L_y)`, with
:math:`q^2 = (2-2\cos(\pi m/L_x)) + (2-2\cos(\pi n/L_y))`, and the height field
is projected on those. (The periodic Fourier transform used to be applied here as
well, and read :math:`0.55\gamma` on synthetic walled slabs of known tension.)

.. warning::

   The hardwall slab estimator is validated **only on synthetic height fields**
   drawn from a known capillary spectrum, where it agrees with the periodic
   estimator. It is **not reliable on real hardwall slabs**: between walls the
   faces of the condensate are not flat on average but carry a static dome, which
   the projection reads as capillary power. One real run read
   :math:`2.7 \pm 56` where its periodic twin reads 41.9, and subtracting the
   time-averaged face does not repair it. ``slab_surface_tension`` raises a
   warning under ``HARDWALL``; measure the surface tension in a periodic box.

Each returns a :class:`~pimms.lemonade.surface_tension.SurfaceTension` with the
estimate ``gamma``, a per-mode spread ``gamma_std`` (an uncertainty proxy, not a
standard error), the number of modes used, the ``temperature`` used as
:math:`k_B T`, which ``method`` ran, and the raw ``spectrum`` for inspection.
For the slab the spectrum is ``(q, P(q))``, where ``q`` is the lattice
wavenumber :math:`\sqrt{(2-2\cos q_x) + (2-2\cos q_y)}` (``~|q|`` at low
:math:`q`) and ``P`` the frame- and face-averaged :math:`|\mathrm{FFT}(\delta
h)|^2` (under ``HARDWALL`` the cosine-mode power, scaled the same way), so that
:math:`\gamma = L_x L_y k_BT / \langle P q^2 \rangle`; for
the droplet it is ``(l, <|u_l|^2>)`` and ``n_polar`` / ``n_azim`` report the
angular grid. For the slab, ``n_modes`` counts *independent* Fourier
wavevectors: the conjugate :math:`+q` and :math:`-q` coefficients of a real
height field are identical and count once. For the droplet it is the number of
degrees :math:`l` that entered the fit:

.. code-block:: python

   result = st.slab_surface_tension(traj)
   result.gamma, result.gamma_std
   q, power = result.spectrum           # inspect the capillary spectrum

Two sentinels are worth knowing. ``gamma`` is ``nan`` (with ``n_modes = 0`` and
no ``spectrum``) when no frame yielded a usable interface, or, for the droplet,
when fewer than two modes survived. The slab estimator returns ``inf`` (again with
``n_modes = 0``) for a perfectly flat interface - zero capillary power is the
infinite-tension limit, returned explicitly rather than as a divide-by-zero.

.. warning::

   Surface tension is intrinsically noisy for small lattice condensates (few
   long-wavelength modes, rough interfaces). Always check ``gamma_std`` and plot the
   ``spectrum``; use a large interface (a big slab cross-section, or a single large
   droplet) and many frames for a precise value. If the temperature is unknown
   (loaded without a keyfile), pass ``temperature=`` explicitly - otherwise the
   estimators raise.

Worked example
==============

.. code-block:: python

   import pimms.lemonade as lemonade
   from pimms.lemonade import phase_separation as ps
   from pimms.lemonade import surface_tension as st

   traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

   result = ps.analyze(traj)
   if result.is_phase_separated:
       print(f"{result.geometry}: rho_dense = {result.rho_dense:.3f}, "
             f"rho_dilute = {result.rho_dilute:.3f}, "
             f"condensed fraction = {result.condensed_fraction:.2f}")
       if traj.n_dim == 3:
           gamma = st.surface_tension(traj)
           print(f"surface tension = {gamma.gamma:.2f} +/- {gamma.gamma_std:.2f} "
                 f"(reduced units, {gamma.method} estimator, {gamma.n_modes} modes)")
   else:
       print(f"not phase separated: binodal fit {result.binodal.reason or 'ok'}, "
             f"condensed fraction {result.condensed_fraction:.2f}, "
             f"percolating in {100 * result.percolation_fraction:.0f}% of frames")

Run on the ``traj.xtc`` of the ``slab_phase_separation`` demo
(``demo_keyfiles/slab_phase_separation``, 301 frames), one run printed:

.. code-block:: text

   slab: rho_dense = 0.924, rho_dilute = 0.001, condensed fraction = 0.99
   surface tension = 48.04 +/- 0.75 (reduced units, slab estimator, 8 modes)
