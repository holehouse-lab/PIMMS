"""
lemonade - a fast, hierarchical analysis backend for PIMMS lattice trajectories.

Load a simulation from any combination of an XTC, a PDB and a PIMMS keyfile, then
navigate a lazy object hierarchy::

    import pimms.lemonade as lemonade

    traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

    traj.radius_of_gyration()          # (n_frames, n_chains) - whole trajectory, vectorised
    frame = traj[0]                    # a Frame view
    polymer = frame[3]                 # a Polymer view (chain 3 in frame 0)
    polymer.radius_of_gyration         # scalar
    polymer.whole_positions            # chain made whole across PBC
    for cluster in frame.clusters:     # connected-component clusters
        cluster.radius_of_gyration, cluster.volume

The whole trajectory is stored columnar (one ``(n_frames, n_atoms, 3)`` int array),
positions are converted back to the integer lattice in one vectorised step, and PBC
unwrapping runs in a compiled kernel - so loading and analysis stay fast even for
large systems.
"""

from ._load import load
from .trajectory import LatticeTrajectory
from .frame import Frame
from .polymer import Polymer
from .cluster import Cluster
from . import phase_separation
from . import surface_tension

__all__ = ["load", "LatticeTrajectory", "Frame", "Polymer", "Cluster",
           "phase_separation", "surface_tension"]
