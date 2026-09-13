.. _lemonade-reference:

=================
API reference
=================

Auto-generated from the ``pimms.lemonade`` docstrings. The public surface of the
package is exactly ``load``, ``LatticeTrajectory``, ``Frame``, ``Polymer``,
``Cluster`` and the ``phase_separation`` and ``surface_tension`` modules;
everything else (the topology, the backing store, the batched numeric core and
the compiled PBC kernel) is an implementation detail.

Loading
=======

.. autofunction:: pimms.lemonade.load

The object hierarchy
====================

.. autoclass:: pimms.lemonade.LatticeTrajectory
   :members:
   :special-members: __getitem__, __len__, __iter__

.. autoclass:: pimms.lemonade.Frame
   :members:
   :special-members: __getitem__, __len__, __iter__

.. autoclass:: pimms.lemonade.Polymer
   :members:
   :special-members: __len__

.. autoclass:: pimms.lemonade.Cluster
   :members:
   :special-members: __len__, __iter__

Phase separation
================

.. automodule:: pimms.lemonade.phase_separation
   :members:

Surface tension
===============

.. automodule:: pimms.lemonade.surface_tension
   :members:
