.. PIMMS documentation master file.

=======================================================
PIMMS: Polymer Interactions in Multi-component MixtureS
=======================================================

**PIMMS** is a lattice-based, coarse-grained Monte Carlo simulation engine for
exploring the phase behaviour and conformational properties of polymer systems -
single homo- or hetero-polymers, many-chain mixtures, and biomolecular
condensates - in both 2D and 3D.

Highlights
==========

* **Easy to use.** Installation provides a command-line executable (``PIMMS``)
  that runs a simulation directly from a plain-text *keyfile*.
* **Fast.** The hot loops are written in optimised Cython that compiles to native
  C, with an optional multi-threaded (OpenMP) kernel.
* **Flexible interactions.** Drive interactions over three distinct length scales
  (short / long / super-long range) plus solvation and backbone-angle terms,
  all defined in a simple parameter file.
* **Many components.** Simulate arbitrary mixtures of distinct polymer species to
  study co-assembly and phase separation.
* **Rich move set.** Local and collective Monte Carlo moves - crankshaft,
  reptation (slither), cooperative pull, rigid-body cluster moves, virtual-move
  Monte Carlo (VMMC) and temperature-switch (TSMMC) moves - to sample efficiently
  and escape kinetic traps.

If you are new to PIMMS, start with :doc:`installation` and then work through the
:doc:`overview`.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   overview
   input_files
   restart_files
   output_files
   keywords

.. toctree::
   :maxdepth: 2
   :caption: Moves

   moves/index

.. toctree::
   :maxdepth: 2
   :caption: Advanced Features

   advanced/index
   advanced/quench
   advanced/freeze
   advanced/parallelization
   advanced/tsmmc
   advanced/custom_analysis
   advanced/reference_controls

.. toctree::
   :maxdepth: 2
   :caption: Analysis (lemonade)

   lemonade/index
   lemonade/loading
   lemonade/hierarchy
   lemonade/conformational
   lemonade/phase_separation
   lemonade/reference

.. toctree::
   :maxdepth: 1
   :caption: Development

   development/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
