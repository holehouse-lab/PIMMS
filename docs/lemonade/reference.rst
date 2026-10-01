.. _lemonade-reference:

=================
API reference
=================

Auto-generated from the ``pimms.lemonade`` docstrings. The names exported by the
package are exactly ``load``, ``LatticeTrajectory``, ``Frame``, ``Polymer``,
``Cluster`` and the ``phase_separation`` and ``surface_tension`` modules. Two
further classes are reachable from a loaded trajectory - ``traj.topology`` is a
``Topology`` and ``traj.store`` a ``TrajectoryStore`` - and are documented at the
end of this page for that reason; they live in private modules and are not
constructed by hand in normal use. Everything else (the batched numeric core and
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

Topology and backing store
==========================

What ``traj.topology`` and ``traj.store`` return. ``load`` builds both; the
methods here are what the hierarchy above is implemented with.

.. autoclass:: pimms.lemonade._topology.Topology
   :members:

.. autoclass:: pimms.lemonade._store.TrajectoryStore
   :members:
