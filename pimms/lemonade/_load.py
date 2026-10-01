## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

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
XTC bead order exactly) and is refined with keyfile chain types when available.
"""

import math
import numbers
import os
import warnings

import numpy as np

from ._topology import Topology, pdb_chain_labels
from ._store import TrajectoryStore
from .trajectory import LatticeTrajectory

DEFAULT_SPACING = 3.65   # PIMMS LATTICE_TO_ANGSTROMS default (v0.1.34+)

# PIMMS hardcodes these two names for the compact phase of a RESIZED_EQUILIBRATION
# run (simulation.py:216-217); nothing else in PIMMS writes them, so seeing one is
# sufficient evidence that the frames were produced in the resized box under the
# forced hardwall of that phase - including when no keyfile is passed at all.
_EQUILIBRATION_FILENAMES = ("eq_traj.xtc", "eq_START.pdb")


def _is_equilibration_file(path):
    """Is ``path`` one of PIMMS's resized-equilibration filenames?

    Parameters
    ----------
    path : str or None
        A trajectory or topology path, or ``None``.

    Returns
    -------
    bool
        ``True`` if the basename is ``eq_traj.xtc`` or ``eq_START.pdb``, i.e. the
        file belongs to the compact equilibration phase of a
        ``RESIZED_EQUILIBRATION`` run.
    """
    if path is None:
        return False
    return os.path.basename(str(path)) in _EQUILIBRATION_FILENAMES


def _restart_snapshot(keydict, keyfile, keyword):
    """Read the restart file a keyfile points at and return the parsed snapshot.

    ``RESTART_OVERRIDE_HARDWALL`` and ``RESTART_OVERRIDE_DIMENSIONS`` tell PIMMS to
    discard the keyfile's ``HARDWALL`` / ``DIMENSIONS`` and adopt the snapshot's
    instead. That substitution happens inside restart processing, which only a full
    keyfile parse runs, so lemonade (which parses with ``parse_only=True``) has to
    do it itself by opening the restart pickle.

    Parameters
    ----------
    keydict : dict
        The ``parse_only=True`` keyword lookup, which must carry ``RESTART_FILE``.
    keyfile : str
        Path to the keyfile. A relative ``RESTART_FILE`` is resolved against the
        keyfile's own directory first, because that is where a run directory keeps
        it; the current working directory is only a fallback.
    keyword : str
        The override keyword that made this read necessary, quoted in any error.

    Returns
    -------
    pimms.restart.RestartObject
        The snapshot, from which ``.hardwall`` and ``.dimensions`` are read.

    Raises
    ------
    ValueError
        If the restart file cannot be found, or cannot be read as a restart file.
        Guessing is not an option here: the whole point of the override keyword is
        that the keyfile value is the wrong one.
    """
    raw = str(keydict["RESTART_FILE"])
    candidates = []
    if os.path.isabs(raw):
        candidates.append(raw)
    else:
        candidates.append(os.path.join(os.path.dirname(os.path.abspath(keyfile)), raw))
        candidates.append(os.path.abspath(raw))

    for path in candidates:
        if os.path.isfile(path):
            from pimms.restart import RestartObject, RestartException
            snapshot = RestartObject()
            try:
                snapshot.build_from_file(path)
            except RestartException as e:
                raise ValueError(
                    f"lemonade.load: the keyfile sets {keyword}, so the effective "
                    f"boundary condition/box must come from the restart file "
                    f"'{path}', but that file could not be read ({e}). Pass "
                    f"hardwall= / dimensions= explicitly.")
            return snapshot

    raise ValueError(
        f"lemonade.load: the keyfile sets {keyword} : True, which means PIMMS took "
        f"the run's boundary condition/box from the restart file rather than from "
        f"the keyfile - so the keyfile values are NOT the ones the run used. The "
        f"restart file 'RESTART_FILE : {raw}' could not be found (looked in "
        f"{', '.join(repr(c) for c in candidates)}). Put it beside the keyfile, or "
        f"pass hardwall= / dimensions= explicitly.")


def load(xtc=None, pdb=None, keyfile=None, *, spacing=None, dimensions=None,
         hardwall=None, temperature=None, start=None, stop=None, step=None,
         n_frames=None, verbose=False):
    """Load a PIMMS trajectory and return a :class:`LatticeTrajectory`.

    Parameters
    ----------
    xtc : str, optional
        Path to the XTC trajectory holding the coordinates. Requires ``pdb``,
        because mdtraj needs a topology to read it (default ``None``).
    pdb : str, optional
        Path to the PDB giving the topology, and on its own a single frame (e.g.
        ``START.pdb``). Its bead order matches the XTC exactly (default
        ``None``).
    keyfile : str, optional
        Path to the PIMMS keyfile. Optional, but authoritative for spacing,
        dimensions, hardwall and chain types (default ``None``).
    spacing : float, optional
        Lattice spacing in angstroms, overriding the keyfile
        ``LATTICE_TO_ANGSTROMS``. Default ``None``: taken from the keyfile, or
        ``DEFAULT_SPACING`` when there is none.
    dimensions : sequence of int, optional
        Box extent in lattice units (2 or 3 positive integers), overriding
        everything below. Default ``None``: the resized-equilibration box for an
        ``eq_`` trajectory, else the restart file's box under
        ``RESTART_OVERRIDE_DIMENSIONS``, else the keyfile ``DIMENSIONS``, else
        inferred from the trajectory's own box record.
    hardwall : bool, optional
        Whether the run used hard walls, overriding everything below. Default
        ``None``: ``True`` for an ``eq_`` trajectory (that phase is always run
        under hard walls), else the restart file's flag under
        ``RESTART_OVERRIDE_HARDWALL``, else the keyfile ``HARDWALL``, else
        ``False``.
    temperature : float, optional
        Override (or supply, when no keyfile is given) the simulation
        temperature; must be a finite positive number. Only needed by the
        surface-tension estimators, which use it for :math:`k_BT`. Default
        ``None``: the keyfile ``TEMPERATURE``, or ``QUENCH_END`` for a
        ``QUENCH_RUN`` keyfile, or ``None`` (unknown) without a keyfile.
    start, stop, step : int, optional
        Frame slice applied at load time (default ``None``, i.e. every frame).
    n_frames : int, optional
        If given (and smaller), evenly subsample down to this many frames
        (default ``None``).
    verbose : bool, optional
        Print a one-line load summary, including the lattice round-off residual
        (default ``False``).

    Returns
    -------
    LatticeTrajectory
        The loaded trajectory, wrapping a
        :class:`~pimms.lemonade._store.TrajectoryStore`.

    Raises
    ------
    ValueError
        If neither ``xtc`` nor ``pdb`` is given, if ``xtc`` is given without
        ``pdb``, if ``n_frames`` is not a positive integer, if ``spacing`` is not
        finite and positive, if the box dimensions are neither given, in the
        keyfile, nor recorded in the trajectory (or are not 2 or 3 positive
        integers), if ``hardwall`` is not a bool, if ``temperature`` is not a
        finite positive number, if the coordinates do not fit in int32 lattice
        coordinates, if the topology's bead count disagrees with the trajectory
        (mdtraj itself refuses a PDB/XTC pair with different bead counts while
        reading them), if the ``start``/``stop``/``step`` selection
        keeps no frames, or if the keyfile sets ``RESTART_OVERRIDE_HARDWALL`` /
        ``RESTART_OVERRIDE_DIMENSIONS`` but its ``RESTART_FILE`` cannot be found
        or read (the keyfile values are then known to be the wrong ones, so
        guessing is not an option).

    Notes
    -----
    A large lattice round-off residual, a box or box *dimensionality* that
    disagrees with the trajectory's own record, a keyfile CHAIN block that does
    not match the PDB, the forced hardwall of a resized-equilibration
    trajectory, and a quench run's use of ``QUENCH_END`` as the temperature all
    raise a warning rather than an error.

    The keyfile is read with ``parse_only=True``, which does not run PIMMS's
    restart reconciliation, so the keywords that make the literal keyfile differ
    from the run PIMMS actually performed are resolved here instead: the
    ``RESTART_OVERRIDE_*`` keywords are honoured by reading the restart file, the
    ``eq_`` files of a ``RESIZED_EQUILIBRATION`` run are loaded in the compact
    box with hard walls, and under a ``RESTART_FILE`` the keyfile ``CHAIN`` lines
    (which PIMMS discards) are not applied.

    Keyfile chain types are only applied when the ``CHAIN`` lines expand onto
    the trajectory's chains in order and, if the PDB carries chain identifiers,
    reproduce the PDB's own partition of chains into types (PIMMS writes one
    identifier per chain type). A keyfile listing the same ``(count, sequence)``
    types in a different order - the ``keyfile_used.kf`` of a restart run with
    ``EXTRA_CHAIN`` chains, for example - keeps the PDB labels without a
    warning, because they are then the right ones.
    """
    import mdtraj as md

    if xtc is None and pdb is None:
        raise ValueError("load() needs at least an xtc (with a pdb topology) or a pdb")
    if n_frames is not None and (isinstance(n_frames, (bool, np.bool_)) or
                                 not isinstance(n_frames, (int, np.integer)) or
                                 n_frames < 1):
        raise ValueError("n_frames must be a positive integer")
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
    n_loaded = traj.n_frames
    if start is not None or stop is not None or step is not None:
        traj = traj[slice(start, stop, step)]
    if n_frames is not None and n_frames < traj.n_frames:
        idx = np.linspace(0, traj.n_frames - 1, int(n_frames)).round().astype(int)
        traj = traj[idx]
    if traj.n_frames == 0:
        raise ValueError(
            f"frame selection start={start}, stop={stop}, step={step} keeps no frames "
            f"of a {n_loaded}-frame trajectory")

    # lattice spacing (angstroms). NB: KeyFileParser(parse_only=True) does not fill
    # defaults, so optional keys are read with .get().
    if spacing is None:
        spacing = float(keydict.get("LATTICE_TO_ANGSTROMS", DEFAULT_SPACING)) if keydict else DEFAULT_SPACING
    try:
        spacing = float(spacing)
    except (TypeError, ValueError):
        raise ValueError("spacing must be a finite positive number")
    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be a finite positive number")

    # coordinates (nm) -> integer lattice, vectorised
    xyz = np.asarray(traj.xyz)                       # (nf, na, 3), nm
    lattice_f = xyz / (spacing * 0.1)
    rounded_lattice = np.rint(lattice_f)
    int32_info = np.iinfo(np.int32)
    if (not np.all(np.isfinite(rounded_lattice)) or
            (rounded_lattice.size and
             (rounded_lattice.min() < int32_info.min or
              rounded_lattice.max() > int32_info.max))):
        raise ValueError(
            "trajectory coordinates cannot be represented as int32 lattice coordinates")
    lattice = rounded_lattice.astype(np.int32)
    residual = float(np.abs(lattice_f - lattice).max()) if lattice_f.size else 0.0

    # A large round-off residual means the coordinates do not sit on the integer
    # lattice at this spacing - almost always a wrong/omitted LATTICE_TO_ANGSTROMS
    # (or wrong spacing=). That silently corrupts the recovered lattice (non-unit
    # bonds, wrong inferred box), so warn ALWAYS, not only under verbose.
    if residual > 0.05:
        warnings.warn(
            f"lemonade.load: lattice round-off residual is {residual:.3g} (>0.05); the "
            f"coordinates do not fit the integer lattice at spacing {spacing} A. The "
            f"recovered lattice is probably corrupted - pass the right spacing=/keyfile "
            f"(LATTICE_TO_ANGSTROMS).", stacklevel=2)

    # Is this the compact equilibration phase of a RESIZED_EQUILIBRATION run? That
    # phase is ALWAYS hardwall and always in the resized box, whatever the keyfile
    # says (simulation.py:219-222), so its frames must not be analysed as periodic
    # in the production box. The filename is the discriminator rather than the box,
    # because RESIZED_EQUILIBRATION == DIMENSIONS is legal and a box-equality test
    # would then force hardwall onto a genuinely periodic production trajectory.
    is_eq_phase = _is_equilibration_file(xtc) or _is_equilibration_file(pdb)
    eq_notes = []

    # The restart snapshot is only opened if something actually needs it, and only
    # once; a missing restart file must not break a load that does not depend on it.
    _restart_cache = {}

    def _restart(keyword):
        """Memoised :func:`_restart_snapshot` for this load.

        Parameters
        ----------
        keyword : str
            The override keyword that needs the snapshot, quoted in any error.

        Returns
        -------
        pimms.restart.RestartObject
            The parsed restart file, read at most once per ``load()`` call.
        """
        if "snapshot" not in _restart_cache:
            _restart_cache["snapshot"] = _restart_snapshot(keydict, keyfile, keyword)
        return _restart_cache["snapshot"]

    _has_restart = bool(keydict is not None and keydict.get("RESTART_FILE"))

    # box dimensions
    if dimensions is None:
        if is_eq_phase and keydict is not None and keydict.get("RESIZED_EQUILIBRATION"):
            dimensions = tuple(int(d) for d in keydict["RESIZED_EQUILIBRATION"])
            eq_notes.append(f"box set to the RESIZED_EQUILIBRATION box {dimensions} "
                            f"rather than DIMENSIONS")
        elif _has_restart and keydict.get("RESTART_OVERRIDE_DIMENSIONS"):
            # RESTART_OVERRIDE_DIMENSIONS discards the keyfile DIMENSIONS; reading
            # them literally here put every minimum-image distance in the wrong box.
            dimensions = tuple(int(d) for d in _restart("RESTART_OVERRIDE_DIMENSIONS").dimensions)
        elif keydict is not None and keydict.get("DIMENSIONS"):
            dimensions = tuple(int(d) for d in keydict["DIMENSIONS"])
        else:
            box = traj.unitcell_lengths
            if box is None:
                raise ValueError("trajectory has no box and no keyfile/dimensions were given")
            dimensions = tuple(int(round(b / (spacing * 0.1))) for b in box[0])
            # PIMMS writes a 2D system with a z period of exactly one lattice unit
            # (a real 3D box is never that thin), so the box record - not the
            # bead coordinates - decides the dimensionality: a 3D configuration
            # that happens to lie in the z = 0 plane must stay 3D.
            if len(dimensions) == 3 and dimensions[2] <= 1:
                dimensions = dimensions[:2]
    try:
        dimensions = tuple(dimensions)
    except TypeError:
        raise ValueError("dimensions must be a 2D or 3D sequence of positive integers")
    if (len(dimensions) not in (2, 3) or
            any(isinstance(d, (bool, np.bool_)) or
                not isinstance(d, numbers.Integral) or d <= 0 or d > int32_info.max
                for d in dimensions)):
        raise ValueError("dimensions must be a 2D or 3D sequence of positive integers")
    dimensions = tuple(int(d) for d in dimensions)
    n_dim = len(dimensions)

    # cross-check against the trajectory's own box record where one exists: a
    # keyfile DIMENSIONS that disagrees with the XTC/CRYST1 box (wrong keyfile,
    # or the eq_ trajectory of a RESIZED_EQUILIBRATION run, written in the
    # smaller box) silently breaks every minimum-image / cluster / profile
    # calculation, so surface it loudly.
    if traj.unitcell_lengths is not None:
        _record = tuple(int(round(b / (spacing * 0.1)))
                        for b in traj.unitcell_lengths[0])
        # Dimensionality first. The old check truncated the box record to the length
        # of the resolved `dimensions` before comparing, so 2D dimensions against a
        # 3D trajectory compared only x and y, matched, and stayed silent while the
        # whole z column was thrown away below. The z period of exactly one lattice
        # unit is the same 2D marker used when inferring dimensions above, and it is
        # in the one column the truncated comparison never read.
        _record_n_dim = 2 if (len(_record) == 3 and _record[2] <= 1) else len(_record)
        if _record_n_dim != n_dim:
            _discarded = ("the z coordinate of every bead is being DISCARDED"
                          if n_dim == 2 else
                          "a flat 2D trajectory is being analysed in a 3D box")
            warnings.warn(
                f"lemonade.load: {n_dim}D box dimensions {dimensions} "
                f"(keyfile/argument) disagree with the trajectory's own "
                f"{_record_n_dim}D box record {_record[:_record_n_dim]} - "
                f"{_discarded}, which changes every distance, cluster and profile. "
                f"Check you are loading the matching keyfile/trajectory pair.",
                stacklevel=2)
        else:
            _box_dims = _record[:n_dim]
            if _box_dims != tuple(dimensions[:len(_box_dims)]):
                warnings.warn(
                    f"lemonade.load: box dimensions {dimensions} (keyfile/argument) "
                    f"disagree with the trajectory's own box record {_box_dims}. All "
                    f"periodic-image and cluster calculations will use {dimensions} - "
                    f"check you are loading the matching keyfile/trajectory pair "
                    f"(eq_ trajectories from RESIZED_EQUILIBRATION runs use the "
                    f"smaller equilibration box).", stacklevel=2)

    # canonicalise into the box (agnostic to whether the trajectory was written
    # wrapped or PBC-unwrapped); lemonade re-derives whole chains itself
    lattice[..., :n_dim] = np.mod(lattice[..., :n_dim], np.array(dimensions, dtype=np.int32))
    if n_dim == 2:
        lattice[..., 2] = 0

    # topology from the PDB (exact XTC bead order); keyfile refines chain types
    topology = Topology.from_mdtraj(traj.topology)
    pdb_labelled = pdb_chain_labels(traj.topology) is not None
    keyfile_types_applied = False
    # Under a RESTART_FILE the keyfile CHAIN lines are NOT what the run used: PIMMS
    # discards them and rebuilds the composition from the snapshot
    # (keyfile_parser.py:1950-2015), and an EXTRA_CHAIN whose sequence already
    # exists joins that existing chain type instead of defining a new one
    # (restart.py:328-348). Numbering each keyfile line as a fresh type therefore
    # split PIMMS's merged types apart (and, in the converse case, merged two real
    # types into one), silently mislabelling every per-type average. The PDB chain
    # identifiers are PIMMS's own chainType order, so they are the ones to keep.
    if keydict is not None and keydict.get("CHAIN") and not _has_restart:
        specs = list(keydict["CHAIN"]) + list(keydict.get("EXTRA_CHAIN") or [])
        # PIMMS upper-cases CHAIN sequences during sanitisation (the default
        # CASE_INSENSITIVE_CHAINS=True), but parse_only=True skips that step -
        # so a lower-case keyfile would silently fail to match the (upper-case)
        # PDB residue names and the keyfile types would be dropped.
        if keydict.get("CASE_INSENSITIVE_CHAINS", True):
            specs = [[n, str(seq).upper()] for (n, seq) in specs]
        # The PDB chain identifiers, when present, are PIMMS's own partition of
        # the chains into types, so the keyfile types are only taken if they
        # reproduce it. The keyfile_used.kf of a restart run lists one CHAIN line
        # per type while the trajectory has its EXTRA_CHAIN chains appended at
        # the end, so its lines do not expand onto the chains in order; when two
        # types share a sequence they expanded onto the wrong chains without a
        # word. A keyfile that holds the same (count, sequence) types in another
        # order describes this run, and the PDB labels are then already right,
        # so they are kept without a warning.
        typed_topology = topology.with_keyfile_types(specs, labelled=pdb_labelled)
        if typed_topology is topology:
            if not (pdb_labelled and topology.matches_keyfile_composition(specs)):
                warnings.warn(
                    "lemonade.load: the keyfile CHAIN/EXTRA_CHAIN specification does "
                    "not match the PDB topology, so keyfile chain types could not be "
                    "applied. The PDB chain identifiers will be used instead; check "
                    "that the keyfile and trajectory belong to the same run.",
                    stacklevel=2)
        else:
            keyfile_types_applied = True
        topology = typed_topology
    # PIMMS has 62 chain identifiers (A-Z, a-z, 0-9) and every chain type past
    # the 62nd shares the last one, so a PDB using all 62 MAY hide merged types;
    # only the keyfile CHAIN lines can tell. (A PDB with a blank chain column is
    # typed by sequence instead, and cannot have run out of identifiers.)
    if (not keyfile_types_applied and pdb_labelled
            and len(set(int(t) for t in topology.chain_types)) >= 62):
        warnings.warn(
            "lemonade.load: the PDB uses all 62 PIMMS chain identifiers, so any chain "
            "type past the 62nd shares a label with another and would have been merged; "
            "pass keyfile= to recover the true chain types.", stacklevel=2)
    if topology.n_beads != lattice.shape[1]:
        raise ValueError(f"topology describes {topology.n_beads} beads but the "
                         f"trajectory has {lattice.shape[1]}")

    if hardwall is None:
        if is_eq_phase:
            hardwall = True
            eq_notes.append("hardwall set to True (the compact equilibration phase "
                            "is always run under hard walls)")
        elif _has_restart and keydict.get("RESTART_OVERRIDE_HARDWALL"):
            # RESTART_OVERRIDE_HARDWALL discards the keyfile HARDWALL. Taking the
            # keyfile value literally here was silent in both directions and made
            # the connected-component search join (or refuse to join) chains through
            # box faces the run never had.
            hardwall = bool(_restart("RESTART_OVERRIDE_HARDWALL").hardwall)
        else:
            hardwall = bool(keydict.get("HARDWALL", False)) if keydict else False
    elif not isinstance(hardwall, (bool, np.bool_)):
        raise ValueError("hardwall must be True or False")

    if eq_notes:
        warnings.warn(
            f"lemonade.load: this is the resized-equilibration trajectory of a "
            f"RESIZED_EQUILIBRATION run ({', '.join(_EQUILIBRATION_FILENAMES)}), "
            f"which PIMMS runs in the compact box under FORCED hard walls whatever "
            f"the keyfile HARDWALL says, so " + "; ".join(eq_notes) +
            ". Pass hardwall= / dimensions= explicitly to override this.",
            stacklevel=2)
    if temperature is None and keydict is not None:
        if keydict.get("QUENCH_RUN"):
            # a quench run ignores TEMPERATURE: production is sampled at
            # QUENCH_END (frames written during the ramp sit at intermediate
            # temperatures)
            temperature = keydict.get("QUENCH_END")
            if temperature is None:
                warnings.warn(
                    "lemonade.load: the keyfile describes a QUENCH_RUN but sets no "
                    "QUENCH_END, so no trajectory temperature was recorded (TEMPERATURE "
                    "is ignored by PIMMS in a quench). Pass temperature= explicitly for "
                    "temperature-dependent analyses (surface tension).", stacklevel=2)
            else:
                temperature = float(temperature)
                warnings.warn(
                    "lemonade.load: the keyfile describes a QUENCH_RUN, whose production "
                    f"phase runs at QUENCH_END = {temperature} (TEMPERATURE is ignored by "
                    "PIMMS in a quench). Using QUENCH_END as the trajectory temperature; "
                    "frames written during the ramp were sampled at intermediate "
                    "temperatures, so restrict temperature-dependent analyses (surface "
                    "tension) to post-ramp frames.", stacklevel=2)
        else:
            temperature = keydict.get("TEMPERATURE")

    store = TrajectoryStore(lattice, dimensions, spacing, bool(hardwall), topology,
                            times=np.asarray(traj.time, dtype=np.float64),
                            temperature=temperature)

    if verbose:
        print(f"[lemonade] {store.n_frames} frames, {store.n_chains} chains, "
              f"{store.n_beads} beads; box {dimensions}, spacing {spacing} A"
              f"{'' if residual < 1e-3 else f'  (WARNING lattice round-off {residual:.3g})'}")
    return LatticeTrajectory(store)
