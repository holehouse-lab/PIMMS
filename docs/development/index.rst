.. _development:

===========================
Development / API reference
===========================

This section documents PIMMS' internal Python API, generated automatically from
the NumPy-style docstrings in the source. It is aimed at developers extending
PIMMS or scripting against it; **users configuring simulations want the**
:doc:`keyword reference </keywords>` **instead**.

The package is organised around a few core objects: a :class:`~pimms.lattice.Lattice`
holds the grid and the :class:`~pimms.chain.Chain` objects; a
:class:`~pimms.energy.Hamiltonian` evaluates the energy; a
:class:`~pimms.moves.MoveObject` implements the Monte Carlo moves; an
:class:`~pimms.acceptance.AcceptanceCalculator` handles move selection and the
Metropolis criterion; and the :class:`~pimms.simulation.Simulation` ties them
together and drives the run, configured by the :class:`~pimms.keyfile_parser.KeyFileParser`.

Simulation engine
=================

.. automodule:: pimms.simulation
   :members:
   :show-inheritance:

Monte Carlo moves
=================

.. automodule:: pimms.moves
   :members:
   :show-inheritance:

.. automodule:: pimms.acceptance
   :members:
   :show-inheritance:

.. automodule:: pimms.chainTSMMC
   :members:
   :show-inheritance:

Energy
======

.. automodule:: pimms.energy
   :members:
   :show-inheritance:

Lattice, chains and geometry
============================

.. automodule:: pimms.lattice
   :members:
   :show-inheritance:

.. automodule:: pimms.chain
   :members:
   :show-inheritance:

.. automodule:: pimms.lattice_utils
   :members:
   :show-inheritance:

Input parsing & configuration
=============================

.. automodule:: pimms.keyfile_parser
   :members:
   :show-inheritance:

.. automodule:: pimms.parameterfile_parser
   :members:
   :show-inheritance:

.. automodule:: pimms.restart
   :members:
   :show-inheritance:

.. automodule:: pimms.data_structures
   :members:
   :show-inheritance:

Analysis & output
=================

.. automodule:: pimms.analysis_IO
   :members:
   :show-inheritance:

.. automodule:: pimms.analysis_structures
   :members:
   :show-inheritance:

.. automodule:: pimms.lattice_analysis_utils
   :members:
   :show-inheritance:
