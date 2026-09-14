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
:attr:`Frame.clusters <pimms.lemonade.Frame>`, ordered largest first **by bead
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
                                        # reaches the box length on ANY axis
   ps.spanning_fraction(traj, all_axes=True)   # ... on EVERY axis
   ps.droplet_shape(traj)               # frame-averaged largest-cluster geometry

``condensed_fraction`` is the basic order parameter: near zero when the system is
well mixed, approaching one when most material collects into a single condensate.
Every one of these takes ``min_beads=`` to ignore clusters below a size (``1``
except for ``number_of_clusters``, ``spanning_fraction`` and ``droplet_shape``,
which default to ``2``), and ``by=`` must be ``'beads'`` or ``'chains'`` - any
other string raises rather than quietly running a different measurement.
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

   r, rho, sites = ps.radial_density_profile_with_site_counts(traj)
   fit = ps.fit_radial_profile(r, rho, site_counts=sites)   # sites: thin shells are ignored
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
the long axis - with the slab re-centred each frame so it does not smear - is fit to
a two-interface tanh:

.. code-block:: python

   z, rho = ps.slab_density_profile(traj)           # axis defaults to the longest
   fit = ps.fit_slab_profile(z, rho, hardwall=traj.hardwall)
   fit.rho_dense, fit.rho_dilute, fit.interface_width
   fit.half_width                       # half the slab thickness

``axis=`` picks the slab normal explicitly (``0``, ``1`` or, in 3D, ``2``).
Passing ``hardwall=`` tells the fit which model to use: ``False`` pins the slab
centre at the middle of the window, which is right for a periodic profile because
:func:`~pimms.lemonade.phase_separation.slab_density_profile` re-centres it every
frame. ``True`` and the default ``None`` both fit the centre as a free fifth
parameter. ``analyze()`` passes ``traj.hardwall`` for you.

Every bead in the box is binned, dense and dilute alike - that is what makes the
result a density profile a coexistence fit can be run against - so unlike the
radial profile this one needs no clusters and takes no ``min_beads``.

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
  one, which means the model has fitted a dilute slab in a dense background (a
  condensate wetting both walls looks like this). The densities are not silently
  swapped; or
* the density gap is **absent** - below ``1e-3`` in occupied fraction, i.e. the
  profile is homogeneous; or
* the fitted profile never actually **reaches its own asymptotes** inside the box, so
  the reported coexistence densities are extrapolation; or
* the density gap is **the size of the scatter** in the profile - noise, not signal; or
* the slab **fills the box**, leaving no dilute phase for the dense phase to coexist
  with; or
* (droplet) **fewer than four shells are usable** after the filters below; or
* (droplet) the fitted **radius is below two lattice units** - the fit has latched
  onto the handful of sites at the cluster centre rather than a dense phase; or
* (droplet, via :func:`~pimms.lemonade.phase_separation.analyze`) the largest cluster
  **spans the box** in most frames - it is connected to its own periodic image, has
  no single image, and its "radial profile" is centred on an arbitrary point of a
  network.

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
``rho_dense`` / ``rho_dilute`` fall back to robust percentiles of the observed profile
- so they stay bounded and, for a homogeneous system, simply coincide.

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
   density gap. Second, and more fundamentally, the coexistence gap **closes
   continuously** as the critical point is approached, so there is no numerical
   criterion that can draw the line for you.

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

   result = ps.analyze(traj)                     # or geometry='slab' / 'sphere', min_beads=...

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
   result.spanning_fraction             # frames in which it spans ANY axis (no single image)
   result.profile                       # (coordinate, density) for plotting

``geometry`` accepts ``'auto'`` (default), ``'slab'``, ``'sphere'`` and the
synonym ``'droplet'``; anything else raises. ``min_beads`` (default ``2``) is
applied to every cluster-based step. The exact rule behind
``is_phase_separated`` is ``binodal.success`` **and** both densities finite
**and** ``rho_dense > 2 * max(rho_dilute, 1e-6)`` **and**
``condensed_fraction > 0.3`` **and** ``percolation_fraction < 0.5``.

``percolation_fraction`` is the guard against the classic false positive: at
moderate volume fraction the contact clustering of a perfectly homogeneous
solution *percolates*, so ``condensed_fraction`` is ~1 and the largest cluster
holds nearly every bead - yet there is no condensate. A cluster spanning every axis
of the box in most frames is a network, not a droplet, and
``is_phase_separated`` is ``False`` for it. A slab spans two axes and passes.

**Hardwall boxes.** Under ``HARDWALL`` a profile is never rolled (it cannot
wrap). A free slab that clears both walls is aligned frame by frame by a plain
translation of its dense centroid to the window centre (the vacated bins take
that frame's dilute level), so a slab diffusing between the walls is not
smeared into a broad hump; a condensate touching a wall is left where it is and
fit with a single-interface model, since its wall face is not an interface -
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
Because a hardwall profile is not re-centred, the two-interface fit treats the
slab centre as a fitted fifth parameter there; it is fixed at the window centre
only when you pass ``hardwall=False``.

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
  avoiding a 2-10% under-estimate of :math:`\gamma` at typical PIMMS box sizes.
  This is the reliable method. ``axis=`` sets the slab normal (default: the
  longest box axis) and ``n_modes=`` how many independent low-:math:`q` modes to
  average over (default 8).
* **Droplet** - the radius :math:`R(\theta,\phi)` fluctuates in spherical-harmonic
  modes with :math:`\langle|u_{lm}|^2\rangle = k_BT/(\gamma R_0^2 (l-1)(l+2))` for
  :math:`l \ge 2` (up to ``l_max=``, default 5). Best-effort: it needs a single,
  compact, reasonably large droplet sampled over many frames, and it ignores
  clusters below ``min_beads=`` (default 30 here, not 2). The interface radius in
  each angular bin is that of the outermost bead, so a grid much coarser than one
  bead per bin *over*-estimates :math:`\gamma` (the maximum over a wide bin sits
  above the surface, and the bin-to-bin scatter of those maxima is counted as
  fluctuation): on deformed spheres of known :math:`\gamma` with :math:`R_0 = 12`
  a fixed :math:`8 \times 16` grid gave :math:`1.27\gamma`, :math:`16 \times 32`
  gave :math:`1.04\gamma` and :math:`24 \times 48` gave :math:`1.005\gamma`.
  Those spheres had continuous radii. On a droplet filled on the lattice the
  outermost radius in each bin is a whole-site quantity, and that rounding is
  white noise across the bins: it projects onto every mode, is read as extra
  fluctuation, and pulls the estimate *low*. Lattice droplets of known
  :math:`\gamma` sampled from the same spectrum gave :math:`0.92\gamma` at
  :math:`\gamma = 0.5\,k_BT` and :math:`0.85\gamma` at :math:`\gamma = 1.5\,k_BT`
  for :math:`R_0 = 12`, and :math:`0.95`, :math:`0.89` and :math:`0.84\gamma` for
  :math:`R_0 = 8, 12, 18` at :math:`\gamma = k_BT`; a finer grid makes this
  worse. Treat a droplet estimate as good to 10 to 20 % and prefer the slab
  estimator whenever the geometry allows it.

The droplet grid is therefore sized from the droplet itself:
``n_polar = 2 R0`` clamped to the range 8 to 64, and ``n_azim = 2 n_polar``,
with :math:`R_0` estimated from the **median** largest-cluster bead count over the frames
analysed. Sizing it from the first frame that held a cluster instead let a
not-yet-condensed leading frame (frame 0 is the start configuration, and
``SAVE_EQ`` keeps the equilibration frames) choose a grid far too fine for the
real droplets, which the per-frame coverage test then rejected one by one,
silently. Any frame whose droplet still fills fewer than half the bins is
skipped, and a warning names how many. Pass ``n_polar=`` / ``n_azim=`` to size
the grid yourself; either way the grid used is reported on the result.

Under ``HARDWALL`` the slab estimator uses only faces that are not pressed against
a wall: a wall face is flat because the wall is, and counting its zero capillary
power halved the spectrum and doubled :math:`\gamma` for a wetting condensate.
The capillary spectrum is a Fourier transform over the in-plane axes, which
assumes the height field is periodic in the plane; under ``HARDWALL`` the
in-plane axes are walls too, so treat the hardwall slab estimate as
approximate.

Each returns a :class:`~pimms.lemonade.surface_tension.SurfaceTension` with the
estimate ``gamma``, a per-mode spread ``gamma_std`` (an uncertainty proxy), the
number of modes used, the ``temperature`` used as :math:`k_B T`, which ``method``
ran, and the raw ``spectrum`` for inspection (``(q, P(q))`` for the slab,
``(l, <|u_l|^2>)`` for the droplet). For the droplet, ``n_polar`` and
``n_azim`` report the angular grid. ``n_modes`` counts *independent* Fourier
wavevectors: the conjugate :math:`+q` and :math:`-q` coefficients of a real
height field are identical and count once:

.. code-block:: python

   result = st.slab_surface_tension(traj)
   result.gamma, result.gamma_std
   q, power = result.spectrum           # inspect the capillary spectrum

Two sentinels are worth knowing. ``gamma`` is ``nan`` (with ``n_modes = 0``) when
no frame yielded a usable interface, and ``inf`` for a perfectly flat one - zero
capillary power is the infinite-tension limit, returned explicitly rather than as
a divide-by-zero.

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
