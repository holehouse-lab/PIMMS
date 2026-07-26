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

Densities throughout are **occupied lattice-site fractions** in ``[0, 1]``, so they
are directly comparable across box sizes.

Clusters and condensates
========================

Clusters are connected groups of chains, found per frame by
:attr:`Frame.clusters <pimms.lemonade.Frame>`, ordered largest first **by bead
count** (which is not the same as chain count once chains differ in length). The
largest cluster - the condensate - is ``frame.droplet``. Each cluster carries its own geometry
(:doc:`hierarchy`): ``radius_of_gyration``, ``asphericity``, ``sphericity``,
``volume``, ``surface_area``, ``density``, ``radial_density_profile()`` and its
composition.

Order parameters
================

Simple, per-frame measures of *how* phase separated the system is:

.. code-block:: python

   ps.condensed_fraction(traj)          # (n_frames,) fraction of beads in the largest cluster
   ps.largest_cluster_size(traj)        # (n_frames,) beads in the largest cluster
   ps.number_of_clusters(traj, min_beads=2)
   ps.cluster_size_distribution(traj)   # all cluster sizes pooled across frames

``condensed_fraction`` is the basic order parameter: near zero when the system is
well mixed, approaching one when most material collects into a single condensate.

The binodal (coexistence densities)
===================================

The dense- and dilute-phase densities are read off a density profile fit to a
hyperbolic tangent. Two geometries are supported.

**Droplet (spherical).** A radial profile about the condensate centre of mass -
occupied fraction as a function of distance - falls from the dense core to the
dilute background:

.. code-block:: python

   r, rho = ps.radial_density_profile(traj)
   fit = ps.fit_radial_profile(r, rho)
   fit.rho_dense, fit.rho_dilute        # coexistence densities
   fit.radius, fit.interface_width      # droplet radius and interface width

The fit is
:math:`\rho(r) = \tfrac12(\rho_d + \rho_v) - \tfrac12(\rho_d - \rho_v)\tanh((r-R)/w)`.

**Slab.** For a slab condensate that spans the periodic plane and is bounded along
one axis (the geometry of the ``slab_phase_separation`` demo), a 1D profile along
the long axis - with the slab re-centred each frame so it does not smear - is fit to
a two-interface tanh:

.. code-block:: python

   z, rho = ps.slab_density_profile(traj, axis=2)   # axis defaults to the longest
   fit = ps.fit_slab_profile(z, rho)
   fit.rho_dense, fit.rho_dilute, fit.interface_width

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

* the fitted profile never actually **reaches its own asymptotes** inside the box, so
  the reported coexistence densities are extrapolation; or
* the density gap is **the size of the scatter** in the profile - noise, not signal; or
* the slab **fills the box**, leaving no dilute phase for the dense phase to coexist
  with.

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

   Two reasons. First, :func:`slab_density_profile` re-centres the slab every frame,
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
auto-detects the geometry (slab if one box axis is much longer than the others,
else spherical):

.. code-block:: python

   result = ps.analyze(traj)

   result.geometry                      # 'sphere' or 'slab'
   result.rho_dense, result.rho_dilute  # binodal
   result.condensed_fraction            # time-averaged
   result.binodal.interface_width
   result.shape                         # {'radius_of_gyration', 'sphericity', ...}
   result.is_phase_separated            # usable fit AND density gap AND most material condensed
   result.profile                       # (coordinate, density) for plotting

Surface tension from undulations
================================

Both surface-tension estimators use capillary-wave theory - the interface's
fluctuation spectrum. Because PIMMS uses :math:`\exp(-\Delta E/T)` with
:math:`k_B = 1`, the trajectory temperature *is* :math:`k_B T`, and the returned
:math:`\gamma` is in **reduced units** (interaction energy per lattice area).

.. code-block:: python

   st.surface_tension(traj)             # auto-dispatch by geometry
   st.slab_surface_tension(traj)        # planar capillary waves  (robust)
   st.droplet_surface_tension(traj)     # spherical-harmonic shape fluctuations

* **Slab** - the two flat interfaces of a box-spanning condensate have a height
  field :math:`h(x,y)` obeying :math:`\langle|h(q)|^2\rangle = k_BT/(\gamma A q^2)`;
  fitting the low-:math:`q` spectrum gives :math:`\gamma`. This is the reliable
  method.
* **Droplet** - the radius :math:`R(\theta,\phi)` fluctuates in spherical-harmonic
  modes with :math:`\langle|u_{lm}|^2\rangle = k_BT/(\gamma R_0^2 (l-1)(l+2))` for
  :math:`l \ge 2`. Best-effort: it needs a single, compact, reasonably large droplet
  sampled over many frames.

Each returns a :class:`~pimms.lemonade.surface_tension.SurfaceTension` with the
estimate ``gamma``, a per-mode spread ``gamma_std`` (an uncertainty proxy), the
number of modes used, and the raw ``spectrum`` for inspection:

.. code-block:: python

   result = st.slab_surface_tension(traj)
   result.gamma, result.gamma_std
   q, power = result.spectrum           # inspect the capillary spectrum

.. warning::

   Surface tension is intrinsically noisy for small lattice condensates (few
   long-wavelength modes, rough interfaces). Always check ``gamma_std`` and plot the
   ``spectrum``; use a large interface (a big slab cross-section, or a single large
   droplet) and many frames for a precise value. If the temperature is unknown
   (loaded without a keyfile), pass ``temperature=`` explicitly.

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
       gamma = st.surface_tension(traj)
       print(f"surface tension = {gamma.gamma:.2f} +/- {gamma.gamma_std:.2f} (kT units)")
