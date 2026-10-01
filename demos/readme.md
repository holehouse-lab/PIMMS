# Demos

This directory contains worked examples of simulations with a fixed surface in PIMMS. Each subdirectory holds a `KEYFILE.kf` and is run from inside that subdirectory with

	PIMMS -k KEYFILE.kf

The executable is `PIMMS` in capitals; on a case-sensitive file system (e.g. Linux) `pimms` will not be found.

### Content

* `surface/` - a surface made from one long frozen chain of `A` beads that tiles the bottom face (z = 0) of a 60 x 60 x 60 hardwall box, with 2000 free `B` beads and 60 24-bead `C` chains simulated above it. **Ready to run**: the restart file the keyfile reads (`new_restart.pimms`) is included, and `build_surface.ipynb` shows how it was built.
* `chains-on-a-surface_sims/` - 50 copies of a 30-residue chain, each tethered by its first bead to a frozen surface on the bottom face of a 50 x 50 x 50 hardwall box. **Needs a preparation step**: run the notebook `build_surface_attach.ipynb` first. It writes `surface_restart_attached.pimms`, which the keyfile reads and which is not included, so `PIMMS -k KEYFILE.kf` fails until the notebook has been run. See the README in that directory for the full walk-through.

More example keyfiles (single chains, two-phase and multiphase systems, slabs, bilayers and others) live in `demo_keyfiles/` at the top of the repository.
