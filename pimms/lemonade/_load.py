"""
Loading PIMMS trajectories into lemonade.

``load()`` accepts any of an XTC (coordinates), a PDB (topology, and/or a single
frame) and a PIMMS keyfile (box size, lattice spacing, hardwall, chain types), in
the combinations a user actually has to hand:

* ``xtc`` + ``pdb``            - the usual case (full trajectory).
* ``xtc`` + ``pdb`` + ``keyfile`` - adds authoritative spacing / dimensions /
  hardwall / chain types.
* ``pdb`` only                 - a single frame (e.g. START.pdb).

Coordinates are converted back to the integer lattice in one vectorised step
(``round(nm / (spacing/10))``); the topology comes from the PDB (which matches the
XTC atom order exactly) and is refined with keyfile chain types when available.
"""

import numpy as np

from ._topology import Topology
from ._store import TrajectoryStore
from .trajectory import LatticeTrajectory

DEFAULT_SPACING = 3.65   # PIMMS LATTICE_TO_ANGSTROMS default (v0.1.34+)


def load(xtc=None, pdb=None, keyfile=None, *, spacing=None, dimensions=None,
         hardwall=None, temperature=None, start=None, stop=None, step=None,
         n_frames=None, verbose=False):
    """Load a PIMMS trajectory and return a :class:`LatticeTrajectory`.

    Parameters
    ----------
    xtc, pdb, keyfile : str or None
        Paths. ``xtc`` requires ``pdb`` (mdtraj needs a topology). ``pdb`` alone
        loads a single frame. ``keyfile`` is optional but authoritative for
        spacing / dimensions / hardwall / chain types.
    spacing, dimensions, hardwall : optional overrides
        Bypass the keyfile / inference for the lattice spacing (angstroms), the box
        dimensions (2- or 3-tuple), and the hardwall flag.
    start, stop, step : int, optional
        Frame slice applied at load time.
    n_frames : int, optional
        If given (and smaller), evenly subsample down to this many frames.
    verbose : bool
        Print a one-line load summary (including the lattice round-off residual).
    """
    import mdtraj as md

    if xtc is None and pdb is None:
        raise ValueError("load() needs at least an xtc (with a pdb topology) or a pdb")
    if xtc is not None:
        if pdb is None:
            raise ValueError("loading an xtc requires a pdb topology - pass pdb=...")
        traj = md.load(xtc, top=pdb)
    else:
        traj = md.load(pdb)

    keydict = None
    if keyfile is not None:
        from pimms.keyfile_parser import KeyFileParser
        keydict = KeyFileParser(keyfile, parse_only=True).keyword_lookup

    # frame selection
    if start is not None or stop is not None or step is not None:
        traj = traj[slice(start, stop, step)]
    if n_frames is not None and n_frames < traj.n_frames:
        idx = np.linspace(0, traj.n_frames - 1, int(n_frames)).round().astype(int)
        traj = traj[idx]

    # lattice spacing (angstroms). NB: KeyFileParser(parse_only=True) does not fill
    # defaults, so optional keys are read with .get().
    if spacing is None:
        spacing = float(keydict.get("LATTICE_TO_ANGSTROMS", DEFAULT_SPACING)) if keydict else DEFAULT_SPACING
    spacing = float(spacing)

    # coordinates (nm) -> integer lattice, vectorised
    xyz = np.asarray(traj.xyz)                       # (nf, na, 3), nm
    lattice_f = xyz / (spacing * 0.1)
    lattice = np.rint(lattice_f).astype(np.int32)
    residual = float(np.abs(lattice_f - lattice).max()) if lattice_f.size else 0.0

    # box dimensions
    if dimensions is None:
        if keydict is not None and keydict.get("DIMENSIONS"):
            dimensions = tuple(int(d) for d in keydict["DIMENSIONS"])
        else:
            box = traj.unitcell_lengths
            if box is None:
                raise ValueError("trajectory has no box and no keyfile/dimensions were given")
            dimensions = tuple(int(round(b / (spacing * 0.1))) for b in box[0])
            if lattice.shape[1] and int(lattice[..., 2].max()) == 0 and int(lattice[..., 2].min()) == 0:
                dimensions = dimensions[:2]           # flat in z -> 2D
    dimensions = tuple(int(d) for d in dimensions)
    n_dim = len(dimensions)

    # canonicalise into the box (agnostic to whether the trajectory was written
    # wrapped or PBC-unwrapped); lemonade re-derives whole chains itself
    lattice[..., :n_dim] = np.mod(lattice[..., :n_dim], np.array(dimensions, dtype=np.int32))
    if n_dim == 2:
        lattice[..., 2] = 0

    # topology from the PDB (exact XTC atom order); keyfile refines chain types
    topology = Topology.from_mdtraj(traj.topology)
    if keydict is not None and keydict.get("CHAIN"):
        specs = list(keydict["CHAIN"]) + list(keydict.get("EXTRA_CHAIN") or [])
        topology = topology.with_keyfile_types(specs)
    if topology.n_atoms != lattice.shape[1]:
        raise ValueError(f"topology describes {topology.n_atoms} beads but the "
                         f"trajectory has {lattice.shape[1]}")

    if hardwall is None:
        hardwall = bool(keydict.get("HARDWALL", False)) if keydict else False
    if temperature is None and keydict is not None:
        temperature = keydict.get("TEMPERATURE")

    store = TrajectoryStore(lattice, dimensions, spacing, bool(hardwall), topology,
                            times=np.asarray(traj.time, dtype=np.float64),
                            temperature=temperature)

    if verbose:
        print(f"[lemonade] {store.n_frames} frames, {store.n_chains} chains, "
              f"{store.n_atoms} beads; box {dimensions}, spacing {spacing} A"
              f"{'' if residual < 1e-3 else f'  (WARNING lattice round-off {residual:.3g})'}")
    return LatticeTrajectory(store)
