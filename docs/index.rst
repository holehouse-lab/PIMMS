.. pimms documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PIMMS!
=========================================================
PIMMS is a lattice-based coarse-grained Monte Carlo simulation engine for exploring complex polymer systems.

What is PIMMS?
-----------------


PIMMS is a lattice-based simulation engine that allows both 2D and 3D simulations to be performed. Useful features include:

Easy to use! Upon installation a command-line executable (PIMMS) is available, should be in you $PATH variable, and can be used to run simulations. No messing around, it (should) just work!
Easy to define interaction parameters through a simple parameter file (example included in /demo_keyfiles/demo_1/params.prm)
Easily run fast 2D or 3D lattice based simulations
Run simulations with many distinct components
Run simulations of a single homo or heteropolymer
Run simulations of many copies of polymers to explore phase behaviour
Drive interactions over three distinct length-scales
Various other things



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   output_files
   code/simulation   
   code/lattice
   code/chain
   code/energy



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
