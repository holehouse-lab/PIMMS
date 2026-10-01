## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

##
## simulation
##
## This file represents the main() function for PIMMS. This is where the magic happens!
## 
##

import glob
import os
import random
import sys
import time
import numpy as np
from datetime import datetime

from .lattice import Lattice
from .acceptance import AcceptanceCalculator
from .moves import MoveObject
from . import moves
from . import mega_crank_fast
from .latticeExceptions import SimulationEnergyException
from .latticeExceptions import SimulationException
from .latticeExceptions import AnalysisRoutineException
from .chainTSMMC import TSMMC
from . import pimmslogger
from . import data_structures

from . import energy
from . import analysis_IO
from . import analysis_general
from . import mega_crank # needed to set random seed...
from . import restart

# utility modules 
from . import lattice_utils
from . import lattice_analysis_utils
from . import longrange_utils
from . import IO_utils
from . import nonequilibrium_utils
from . import system_utils
from . import numpy_utils

from . import CONFIG

# if we want to check memory usage we can
# edit this and then call hp.heap() as
# needed in the code
CHECK_MEMORY = False
if CHECK_MEMORY:
    from guppy import hpy
    hp = hpy()

# Moves that can change a chain's INTERNAL conformation. Everything not listed
# here (chain translate, chain rotate, cluster translate, cluster rotate, VMMC) is
# a rigid-body motion of one or more whole chains: it moves a chain around the box
# but leaves its shape, bond for bond, exactly as it was. The jump-and-relax, chain
# TSMMC and multichain TSMMC moves are in the list because their relaxations are
# crankshaft megamoves. A system-wide TSMMC excursion runs the other enabled moves
# as its sub-moves, so Simulation replaces it by those before using this list.
CONFORMATIONAL_MOVES = frozenset(('MOVE_CRANKSHAFT', 'MOVE_SLITHER', 'MOVE_PULL',
                                  'MOVE_CHAIN_PIVOT', 'MOVE_HEAD_PIVOT',
                                  'MOVE_JUMP_AND_RELAX', 'MOVE_CTSMMC',
                                  'MOVE_MULTICHAIN_TSMMC', 'MOVE_SYSTEM_TSMMC'))

# The conformational observables that become run constants when the move set
# cannot reshape a chain. Named in the warnings so the user knows exactly which
# numbers not to trust.
_CONFORMATIONAL_OUTPUT_FILES = ('RG.dat, ASPH.dat, END_TO_END_DIST.dat, '
                                'CHAIN_*_INTSCAL*.dat, CHAIN_*_DISTANCE_MAP.dat and '
                                'CHAIN_*_SCALING_INFORMATION.dat')


def _pimms_version():
    """The installed PIMMS version, or 'unknown' outside a built package."""
    try:
        from . import __version__ as _v
        return str(_v)
    except Exception:
        try:
            from importlib.metadata import version as _metadata_version
            return _metadata_version('idptools-pimms')
        except Exception:
            return 'unknown'


def _mean_over_chains(maps):
    """
    Average a sequence of equally shaped per-chain matrices without stacking them.

    The end-of-run distance map of a chain type is the plain mean, over the
    chains of that type, of each chain's own running-mean map. It used to be
    taken as ``np.asarray(list_of_maps).mean(axis=0)``, which builds an
    ``(n_chains, L, L)`` copy of every map next to the accumulators that
    already hold them, so the last thing a run did was double its distance-map
    memory. Here the maps are added one at a time into a single ``(L, L)``
    array instead.

    The result is bit-identical to the stacked mean. numpy reduces a stack
    along its first axis by adding the slices in order, which is exactly what
    this loop does; the one exception is a ``1 x 1`` map (a single-bead chain),
    where the stack is contiguous along the reduced axis and numpy switches to
    pairwise summation, so those are still stacked - they cost nothing.

    Parameters
    ----------
    maps : iterable of numpy.ndarray
        The per-chain matrices, all of one shape. Consumed once, in order.

    Returns
    -------
    numpy.ndarray
        The element-wise mean, as a new ``float64`` array of the common shape.
        An empty iterable gives what the stacked mean gave (``nan``, with
        numpy's own warning), so the caller's behaviour is unchanged.
    """
    summed = None
    count = 0
    small = []
    for chain_map in maps:
        chain_map = np.asarray(chain_map)
        if chain_map.size <= 1:
            small.append(chain_map)
            continue
        if summed is None:
            summed = np.array(chain_map, dtype=np.float64)
        else:
            summed += chain_map
        count += 1
    if summed is None:
        return np.asarray(small).mean(axis=0)
    return summed / count


def _blank_percolating_clusters(percolating, polymeric_properties, size_properties,
                                radial_density, radial_density_indices):
    """
    Replace the undefined quantities of percolating clusters with sentinels.

    A cluster connected to its own periodic image has no radius of gyration,
    asphericity, convex hull or radial density profile: those are properties of
    a finite object, and in the infinite periodic system the simulation
    represents this cluster is unbounded. The gather still hands back
    coordinates - an arbitrary finite window whose shape depends on which bead
    the walk started from - so every one of these numbers has to be marked
    rather than written as if it were a measurement.

    Parameters
    ----------
    percolating : list of bool
        One flag per size-thresholded cluster, True where the cluster
        percolates the box.

    polymeric_properties : list of list
        One ``[Rg, asphericity]`` entry per cluster; percolating entries are
        replaced by ``[nan, nan]``.

    size_properties : list of list
        One ``[volume, area, density]`` entry per cluster; percolating entries
        are replaced by ``[-1, -1, -1]``, the existing degenerate-hull
        convention for "undefined".

    radial_density : list of list of float
        The emitted radial profiles, one per cluster above the bead threshold.

    radial_density_indices : list of int
        The one-based cluster numbers those profiles belong to.

    Returns
    -------
    tuple
        ``(polymeric_properties, size_properties, radial_density,
        radial_density_indices)`` with the percolating clusters blanked and
        their radial rows removed.

    """
    if not any(percolating):
        return (polymeric_properties, size_properties, radial_density, radial_density_indices)

    polymeric_properties = [
        [float('nan'), float('nan')] if percolating[i] else values
        for i, values in enumerate(polymeric_properties)]

    size_properties = [
        [-1, -1, -1] if percolating[i] else values
        for i, values in enumerate(size_properties)]

    keep = [k for k, idx in enumerate(radial_density_indices) if not percolating[idx - 1]]
    radial_density = [radial_density[k] for k in keep]
    radial_density_indices = [radial_density_indices[k] for k in keep]

    return (polymeric_properties, size_properties, radial_density, radial_density_indices)


def conformation_freezing_warnings(usable_moves, longest_mobile):
    """
    Identify move sets under which a chain's conformation is a conserved quantity.

    Detailed balance is only half of what Boltzmann sampling needs: the move set
    must also be irreducible. PIMMS happily accepts move sets that are not. A set
    built only from the rigid moves (chain translate, chain rotate, cluster
    translate, cluster rotate, VMMC) moves chains around the box while leaving
    every chain's shape exactly as ``insert_chain()`` drew it, so the run samples
    the Boltzmann distribution restricted to the conformational class the seed
    happened to produce. The restriction is invisible from ENERGY.dat or the
    trajectory - the chains really are moving - but every conformational
    observable is a constant, and the bias in the energy has neither a bound nor
    a fixed sign, since it depends entirely on which shapes were drawn.

    Three weaker versions of the same problem are also caught here: pull as the
    only conformational move (pull never displaces a terminus, so the end-to-end
    DISTANCE, which is what END_TO_END_DIST.dat records, is frozen - and rotating
    or translating the chain does not change a distance, so mixing pull with the
    rigid moves does not help); head pivot alone (only the two terminal beads ever
    move); and pivots alone (chain pivot always rotates the SHORTER arm, the
    N-terminal one on a tie, so the two beads at the chain's midpoint are never in
    the rotated arm for any legal pivot point and can never be displaced - which
    also means such a chain can never diffuse).

    This is a warning and not a refusal on purpose. With ``FREEZE_FILE`` pinning
    position as well, a rigid-only move set is the only way PIMMS can express
    rigid-body assembly of fixed-conformation objects, which is legitimate science
    and something the shipped demos already do.

    Parameters
    ----------
    usable_moves : iterable of str
        The ``MOVE_*`` keywords that are enabled AND long enough to act on the
        system (i.e. whose minimum chain length is met by some mobile chain).

    longest_mobile : int
        Number of beads in the longest unfrozen chain. A single bead has no
        shape, so nothing is returned below two. A two-bead chain has exactly
        one internal degree of freedom, its bond length (1, sqrt(2) or, in 3D,
        sqrt(3) lattice units), and a rigid-only move set freezes it just as it
        freezes a longer chain's conformation - so the rigid-only warning is
        issued from two beads up. The three narrower warnings (pull only, head
        pivot only, pivots only) describe interior or midpoint beads a dimer
        does not have, and need three beads or more.

    Returns
    -------
    list of str
        Warning messages, empty if the move set can reshape the chains. At most
        one message is returned; the list form is so a caller can simply extend
        its warning list.
    """
    if longest_mobile < 2:
        return []

    usable = set(usable_moves)
    shape_movers = usable & CONFORMATIONAL_MOVES
    rigid = sorted(usable - CONFORMATIONAL_MOVES)

    if not shape_movers and longest_mobile == 2:
        return [("The enabled move set (%s) contains no move that can change a chain's shape - "
                 "these are all rigid-body moves of whole chains. The longest mobile chain has "
                 "two beads, so the shape at stake is the bond: every dimer keeps the bond length "
                 "(1, sqrt(2) or, in 3D, sqrt(3) lattice units) it was built with for the whole "
                 "run. END_TO_END_DIST.dat and RG.dat will be identical rows, and the energy is "
                 "that of the bond lengths the SEED happened to draw. Add crankshaft, slither, "
                 "head pivot, jump-and-relax or a chain / multi-chain TSMMC move unless fixed "
                 "bonds are what you intend." % ', '.join(sorted(usable)))]

    if longest_mobile < 3:
        return []

    if not shape_movers:
        return [("The enabled move set (%s) contains no move that can change a chain's shape - "
                 "these are all rigid-body moves of whole chains. Every chain will keep the "
                 "conformation it was built with at startup for the whole run. This is still a "
                 "valid Monte Carlo run, but it samples the Boltzmann distribution RESTRICTED to "
                 "the conformations the SEED happened to draw, so the energy and every "
                 "conformational average are conditional on that arbitrary starting set - the "
                 "difference from the true ensemble has no bound and no fixed sign. %s will be "
                 "identical rows describing the initialisation rather than the ensemble (a Flory "
                 "fit on constant data still succeeds and still looks plausible). Add a "
                 "conformational move - crankshaft, slither, pull, chain pivot, head pivot, "
                 "jump-and-relax or a TSMMC move - unless rigid-body assembly of fixed "
                 "conformations is exactly what you intend."
                 % (', '.join(sorted(usable)), _CONFORMATIONAL_OUTPUT_FILES))]

    if shape_movers == {'MOVE_PULL'}:
        return [("MOVE_PULL is the only enabled move that can change a chain's shape, and pull "
                 "never displaces a terminus, so every chain's end-to-end DISTANCE is a constant "
                 "of this run - END_TO_END_DIST.dat will be identical rows. The rigid moves "
                 "(%s) do not release it either: translating or rotating a chain does not change "
                 "a distance. Pair MOVE_PULL with crankshaft or slither if the end-to-end "
                 "distribution is something you intend to measure."
                 % (', '.join(rigid) if rigid else 'none enabled'))]

    if shape_movers == {'MOVE_HEAD_PIVOT'}:
        return [("MOVE_HEAD_PIVOT is the only enabled move that can change a chain's shape, and "
                 "it only ever moves the two terminal beads. Beads 1 to L-2 of every chain will "
                 "sit on their startup sites for the whole run, so %s describe a nearly frozen "
                 "conformation. Add crankshaft, slither, pull or chain pivot to sample the "
                 "interior." % _CONFORMATIONAL_OUTPUT_FILES)]

    if shape_movers <= {'MOVE_CHAIN_PIVOT', 'MOVE_HEAD_PIVOT'}:
        return [("The only enabled shape-changing move%s %s. Chain pivot always rotates the "
                 "SHORTER arm about the drawn pivot bead (the N-terminal arm on a tie), so the two "
                 "beads at each chain's midpoint are never in the rotated arm for any legal pivot "
                 "point and never move relative to the rest of the chain%s. %s Add crankshaft, "
                 "slither or pull."
                 % ('s are pivots' if len(shape_movers) > 1 else ' is a pivot',
                    '(%s)' % ', '.join(sorted(shape_movers)),
                    ', and head pivot only moves the termini'
                    if 'MOVE_HEAD_PIVOT' in shape_movers else '',
                    ("The rigid moves (%s) still carry chains around the box, but the shapes "
                     "reachable are only those pivots about fixed midpoint beads can make."
                     % ', '.join(rigid)) if rigid else
                    ("With no rigid-body move enabled those beads are nailed to their startup "
                     "sites, which also means a chain cannot diffuse and so cannot sample its "
                     "position relative to the other chains.")))]

    return []


class Simulation:
    """
    The master object coordinating a PIMMS simulation.

    A ``Simulation`` is constructed from the completed ``keyword_lookup``
    dictionary created by :class:`pimms.keyfile_parser.KeyFileParser`. The
    parser initializes the dictionary while reading a keyfile, then populates
    it with parsed values, defaults, dynamic defaults, and derived keywords.
    ``Simulation`` uses this dictionary as the configuration snapshot for
    constructing the lattice, chains, interaction model, move machinery, and
    analysis settings. It does not parse the keyfile or apply missing-value
    defaults itself.

    """

    #-----------------------------------------------------------------
    #
    def __init__(self, keyword_lookup):        
        """
        Construct a simulation from a parsed keyfile configuration.

        ``keyword_lookup`` is not built in this constructor. It is the
        ``keyword_lookup`` attribute of a :class:`~pimms.keyfile_parser.KeyFileParser`
        instance after that parser has read the keyfile and applied its
        standard, dynamic, and derived defaults. The dictionary is passed in
        as a completed configuration and is used throughout initialization to
        construct and configure the simulation. Its controlled vocabulary and
        validation are defined by ``KeyFileParser``.

        Parameters
        ----------------
        keyword_lookup : dict
            Completed configuration dictionary from
            :class:`~pimms.keyfile_parser.KeyFileParser.keyword_lookup`.
            Required and optional entries are defined and validated by the
            parser; this constructor consumes the resulting values without
            parsing the keyfile or filling in defaults.

        Raises
        ----------------
        SimulationException
            If an input file is also one of the files the run writes (see
            ``_input_output_collisions``), or if ``RESTART_CONTINUE`` is set and
            the restart file's configuration does not have, under this run's
            Hamiltonian, the energy the checkpoint recorded. Any exception
            raised while the simulation is being set up is first recorded in
            ``log.txt`` as a failed start and then re-raised unchanged.

        """
        # log.txt was (re)initialised by the keyfile parser, so by the time we
        # get here it describes THIS start. If the start fails, the log has to
        # say so in a line that cannot be missed: the directory may well hold a
        # previous run's data, which a failed start does not remove.
        self._last_completed_step = None
        self._last_analysed_step = None
        self._last_checkpoint_step = None
        self._stopped_by_signal = None
        self._stop_in_progress = False
        self._effective_keyfile_written = False
        try:
            self._construct(keyword_lookup)
        except BaseException as error:
            self._log_failed_start(error, 'while the simulation was being set up')
            raise


    #-----------------------------------------------------------------
    #
    def _log_failed_start(self, error, stage):
        """
        Record in ``log.txt`` that this start failed and no simulation was run.

        The parser wipes ``log.txt`` and the Hamiltonian rewrites
        ``parameters_used.prm`` (and, when angle penalties are on,
        ``absolute_energies_of_angles.txt``) before the simulation is known to
        be startable, so after a failed start those files describe a run that
        never happened, next to whatever data an earlier run left in the
        directory.
        This writes the one line that says so. It never raises: a log that
        cannot be written must not hide the error being reported.

        Parameters
        ----------
        error : BaseException
            The exception that stopped the start.

        stage : str
            Where the start had got to, completing the sentence "PIMMS stopped
            ...".

        Returns
        -------
        None
            One ERROR line is appended to ``log.txt``.
        """
        try:
            detail = ' '.join(str(error).split())
            pimmslogger.log_error(
                "START FAILED - NO SIMULATION WAS RUN. PIMMS stopped %s (%s%s). This log, %s and "
                "(if angle penalties are on) %s describe a start that did not happen; any .dat or "
                "restart files in this directory, and any %s, were left by an earlier run and "
                "have not been touched."
                % (stage, type(error).__name__, (': ' + detail[:400]) if detail else '',
                   CONFIG.OUTPUT_USED_PARAMETER_FILE, CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL,
                   CONFIG.EFFECTIVE_KEYFILE_NAME))
        except Exception:
            pass


    #-----------------------------------------------------------------
    #
    @staticmethod
    def _input_output_collisions(keyword_lookup):
        """
        Find input files that are also files this run writes or deletes.

        Every run writes its outputs under fixed names in the working directory
        and deletes a previous run's copies at start-up. An input file that IS
        one of those files (the same file on disk, so a symbolic link or a name
        differing only in case on a case-insensitive filesystem counts) used to
        be overwritten or deleted without a word.

        Three pairings are left alone because the file is read in full before it
        is rewritten with equivalent content: the keyfile and
        ``keyfile_used.kf`` (re-running a ``keyfile_used.kf`` in place is a
        supported way to repeat a run), ``PARAMETER_FILE`` and
        ``parameters_used.prm``, and ``RESTART_FILE`` and ``restart.pimms`` (a
        restart in place, replaced only at the first checkpoint).

        Parameters
        ----------
        keyword_lookup : dict
            The resolved keyword dictionary. ``__KEYFILE``, ``PARAMETER_FILE``,
            ``FREEZE_FILE`` and ``RESTART_FILE`` are the inputs examined.

        Returns
        -------
        list of tuple
            ``(keyword, input path, output name)`` for every collision, where
            ``keyword`` is the keyfile keyword naming the input (``'the
            keyfile'`` for the keyfile itself).
        """
        inputs = [('the keyfile', keyword_lookup.get('__KEYFILE'), CONFIG.EFFECTIVE_KEYFILE_NAME),
                  ('PARAMETER_FILE', keyword_lookup.get('PARAMETER_FILE'), CONFIG.OUTPUT_USED_PARAMETER_FILE),
                  ('FREEZE_FILE', getattr(keyword_lookup.get('FREEZE_FILE'), 'filename', None), None),
                  ('RESTART_FILE', getattr(keyword_lookup.get('RESTART_FILE'), 'filename', None), CONFIG.RESTART_FILENAME)]
        outputs = list(CONFIG.analysis_output_files()) + [
            CONFIG.OUTNAME_LOGFILE, CONFIG.OUTPUT_USED_PARAMETER_FILE, CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL,
            CONFIG.OUTPUT_CHAIN_TO_CHAINID, CONFIG.EFFECTIVE_KEYFILE_NAME, CONFIG.RESTART_FILENAME,
            'START.pdb', 'traj.xtc', 'eq_START.pdb', 'eq_traj.xtc',
            'CONFIG_AT_ENERGY_FAIL.pdb', 'CONFIG_AT_ENERGY_FAIL.xtc']
        collisions = []
        for keyword, path, allowed in inputs:
            if not isinstance(path, str):
                continue
            path = os.path.expanduser(path)
            if not os.path.exists(path):
                continue
            for output in outputs:
                if output == allowed or not os.path.exists(output):
                    continue
                try:
                    same = os.path.samefile(path, output)
                except OSError:
                    same = False
                if same:
                    collisions.append((keyword, path, output))
        return collisions


    #-----------------------------------------------------------------
    #
    def _construct(self, keyword_lookup):
        """
        Build every part of the simulation from the resolved keywords.

        This is the body of the constructor (see ``__init__`` for what
        ``keyword_lookup`` is); it is a separate method only so that
        ``__init__`` can record a failed start in ``log.txt``.

        Parameters
        ----------------
        keyword_lookup : dict
            Completed configuration dictionary from
            :class:`~pimms.keyfile_parser.KeyFileParser.keyword_lookup`.

        Returns
        ----------------
        None
            The simulation's attributes are set in place.

        Raises
        ----------------
        SimulationException
            See ``__init__``.

        """

        ## SET UP THE LOGGER
        IO_utils.status_message('SETTING UP THE SIMULATION', 'major')

        # An input file that is also an output of this run would be overwritten
        # or deleted by it. Refuse before anything else is written (log.txt is
        # the one exception: the parser has already initialised it).
        _collisions = self._input_output_collisions(keyword_lookup)
        if _collisions:
            msg = ("An input file is also a file this run writes, so the run would destroy it: "
                   + "; ".join("%s (%s) is the output file %s" % c for c in _collisions)
                   + ". PIMMS writes its output under fixed names in the working directory; rename "
                     "the input file, or keep it in another directory, and start again."
                   + (" (%s is initialised when the keyfile is parsed, so an input with that name "
                      "has already been overwritten.)" % CONFIG.OUTNAME_LOGFILE
                      if any(c[2] == CONFIG.OUTNAME_LOGFILE for c in _collisions) else ""))
            raise SimulationException(msg)

        ## CORE CONSISTENCY TESTS

        # check int-types in cython/python are consistent
        system_utils.check_dtype_consistency()

        # check we can accomodate 
        # under a restart the parser has already folded EXTRA_CHAIN into CHAIN, so
        # adding the extras again would count them twice
        _capacity_chains = list(keyword_lookup['CHAIN'])
        if not keyword_lookup.get('RESTART_FILE'):
            _capacity_chains += list(keyword_lookup.get('EXTRA_CHAIN') or [])
        system_utils.check_beads_to_grid_mapping(_capacity_chains)
        

        ## SET LOCAL VARIABLES
        # set local variables for use in initialization
        chains                  = keyword_lookup['CHAIN'] 
        temperature             = keyword_lookup['TEMPERATURE']
        random_seed             = keyword_lookup['SEED']
        parameter_file          = keyword_lookup['PARAMETER_FILE']
        non_interacting         = keyword_lookup['NON_INTERACTING'] 
        angles_off              = keyword_lookup['ANGLES_OFF'] 
                
        # set simulation object variables to be used throughout the simulation    
        self.compare_energyfreq   = keyword_lookup['ENERGY_CHECK']
        self.printfreq            = keyword_lookup['PRINT_FREQ']
        self.reduced_printing     = keyword_lookup['REDUCED_PRINTING']
        IO_utils.set_reduced_printing(self.reduced_printing)
        self.enfreq               = keyword_lookup['EN_FREQ'] 
        self.xtcfreq              = keyword_lookup['XTC_FREQ']
        self.n_steps              = keyword_lookup['N_STEPS']
        self.equilibration        = keyword_lookup['EQUILIBRATION']
        self.anafreq              = keyword_lookup['ANALYSIS_FREQ']
        # Frequencies below one mean "disabled".  Keep that state separate
        # from their legacy N_STEPS+10 display sentinel, so nothing that reads
        # a frequency can mistake the sentinel for a real cadence.
        self.disabled_frequencies = frozenset(
            keyword_lookup.get('__DISABLED_FREQUENCIES', ()))
        self.CS_substeps          = keyword_lookup['CRANKSHAFT_SUBSTEPS']
        self.CS_mode              = keyword_lookup['CRANKSHAFT_MODE']
        self.slither_substeps     = keyword_lookup['SLITHER_SUBSTEPS']   # number of slithers applied to each chain per slither megamove
        self.pull_substeps        = keyword_lookup['PULL_SUBSTEPS']      # number of pull moves applied to each chain per pull megamove
        self.vmmc_max_displacement = keyword_lookup['VMMC_MAX_DISPLACEMENT']  # max |translation| per dimension for a VMMC collective move
        self.vmmc_max_cluster      = keyword_lookup['VMMC_MAX_CLUSTER']       # cap on the VMMC cluster-size cutoff draw (clamped to n_chains at runtime)
        self.vmmc_accepted_multichain = 0    # diagnostics: accepted VMMC moves whose cluster had >1 chain
        self.vmmc_max_accepted_cluster = 0   # diagnostics: largest accepted VMMC cluster
        self.LATTICE_TO_ANGSTROMS = keyword_lookup['LATTICE_TO_ANGSTROMS']
        self.autocenter           = keyword_lookup['AUTOCENTER']
        self.trajectory_pbc_unwrap = keyword_lookup['TRAJECTORY_PBC_UNWRAP']

        # set quench keywords
        self.QUENCH_RUN         = keyword_lookup['QUENCH_RUN'] 
        self.QUENCH_START       = keyword_lookup['QUENCH_START'] 
        self.QUENCH_END         = keyword_lookup['QUENCH_END'] 
        self.QUENCH_FREQ        = keyword_lookup['QUENCH_FREQ'] 
        self.QUENCH_STEPSIZE    = keyword_lookup['QUENCH_STEPSIZE'] 
                    
        # set for updates to the TSMMC mode 
        self.TSMMC_USED                = keyword_lookup['__TSMMC_USED']        
        self.TSMMC_INTERPOLATION_MODE  = keyword_lookup['TSMMC_INTERPOLATION_MODE']
        self.TSMMC_JUMP_TEMP           = keyword_lookup['TSMMC_JUMP_TEMP']
        self.TSMMC_STEP_MULTIPLIER     = keyword_lookup['TSMMC_STEP_MULTIPLIER']
        self.TSMMC_NUMBER_OF_POINTS    = keyword_lookup['TSMMC_NUMBER_OF_POINTS']
        self.TSMMC_FIXED_OFFSET        = keyword_lookup['TSMMC_FIXED_OFFSET']
        self.production_hardwall       = keyword_lookup['HARDWALL']

        # initialize freezefile stuff
        self.frozen_chains = []

        # set whether saving at end. 
        self.keyword_lookup    = keyword_lookup
        self.SAVE_AT_END       = keyword_lookup['SAVE_AT_END']

        # set whether saving equilibration steps
        self.SAVE_EQ           = keyword_lookup['SAVE_EQ']

        # parallelization of the crankshaft, slither and pull megamoves (also when
        # they run as system-TSMMC sub-moves). PARALLEL_THREADS of 0 means "use all
        # available cores". The parallel checkerboard kernels are used in both 2D
        # and 3D and honour frozen chains via a per-bead frozen mask; they target
        # the same Boltzmann distribution as the serial kernels, though they follow
        # a different, per-step slower-relaxing Markov chain.
        self.parallelize       = keyword_lookup['PARALLELIZE']
        # (0 resolves to OMP_NUM_THREADS if that is set, otherwise to the CPUs this
        # process may actually run on - the affinity mask, not the machine's core
        # count - so a 4-core batch allocation on a 128-core node starts 4 threads)
        from . import keyfile_parser as _keyfile_parser
        self.parallel_threads  = _keyfile_parser.resolve_parallel_threads(keyword_lookup['PARALLEL_THREADS'])

        # set equilibration offset
        self.EQ_OFFSET = keyword_lookup['EQUILIBRATION_OFFSET']

        # set None as the mdtraj obj for now. This will be updated every time the coordinates of the system are saved
        # if we use set self.SAVE_AT_END=True.
        self.master_traj_obj = None

        # persistent XTC write handle used for the (default) incremental-save path.
        # Kept open for the whole run so each frame is an O(1) append rather than a
        # full reload+resave of the growing trajectory.
        self.xtc_writer = None

        # analysis settings
        self.analysis_settings  = data_structures.AnalysisSettings(cluster_threshold=keyword_lookup['ANA_CLUSTER_THRESHOLD'])

        # set flags for auxiliary chain MC moves (e.g. TSMMC). Set to False to start with
        self.auxillary_chain = False

        # set box size - this is a bit fiddly...
        if keyword_lookup['RESIZED_EQUILIBRATION']:
            
            dimensions = keyword_lookup['RESIZED_EQUILIBRATION']
            self.resize_eq = True
            self.current_xtc_filename = 'eq_traj.xtc'
            self.current_pdb_filename = 'eq_START.pdb'

            # regardless of what keyfile says, we must run initial compact sims with a hardwall
            # boundary to avoid the scenario in which we're re-sizing a system with chains crossing
            # a PBC
            self.hardwall = True 
                        
        else:
            dimensions = keyword_lookup['DIMENSIONS']
            self.resize_eq = False
            self.current_xtc_filename = 'traj.xtc'
            self.current_pdb_filename = 'START.pdb'
            self.hardwall = self.production_hardwall
            
        self.production_dims = keyword_lookup['DIMENSIONS']

        # set values for 10 and 5 percent of the simulation with over-ride values
        # in case we're running especially short simulations        
        if self.n_steps >= 10: 
            self.ten_percent = round(self.n_steps/10)
        else:
            self.ten_percent = 1

        if self.n_steps >= 20:            
            self.five_percent = round(self.n_steps/20)
        else:
            self.five_percent = 1

        self.global_start_time = None
        # the same instant on the monotonic clock: every elapsed time is measured
        # from this, so a wall-clock change during the run cannot corrupt one
        self.global_start_monotonic = None


        ## --------------------------------------------------------------------
        ## Part 1 - Randomization stuff
        ##        

        IO_utils.status_message("Using random seed   : %i" % (random_seed), 'startup')
        IO_utils.status_message("Reference-kernel seed: %i" % (random_seed % CONFIG.C_RAND_MAX),'startup')
        IO_utils.status_message("Kernel PRNG range   : %i" % (CONFIG.C_RAND_MAX),'startup')

        random.seed(random_seed)
        # numpy takes any seed below 2**32; the C_RAND_MAX reduction belongs to the
        # reference kernel's C int seed alone. Reducing numpy's seed by it too gave
        # seeds s and s + 2**31 - 1 identical bead-selector and shuffle streams (only
        # their Python and kernel streams differed). The production kernels are
        # seeded per megamove from Python's generator, not from this reduced value.
        np.random.seed(random_seed % 2**32)
        mega_crank.seed_C_rand(random_seed%CONFIG.C_RAND_MAX)
            

        ## Part 2 - Build the Markov Chain Monte Carlo Metrpolis Acceptance
        #           object and set the various move probabilities therein
        self.ACC       = AcceptanceCalculator(temperature, keyword_lookup)


        ## Part 3 - Build the chain-mover object
        self.MOVER     = MoveObject()


        ## Part 4 - Build the system Hamiltonian based on the
        #           parameter file, or using an empty Hamiltonian
        #           for a non-interacting (Excluded volume) run. Note that non-interacting
        #           only gets used if set to True. Also note that we provide the equilibrium
        #           temperature which may be used to define the angle interaction energies if
        #           requested.
        self.Hamiltonian  = energy.Hamiltonian(parameter_file, 
                                               len(dimensions), 
                                               non_interacting, 
                                               angles_off, 
                                               hardwall = self.hardwall, 
                                               temperature = keyword_lookup['EQUILIBRIUM_TEMPERATURE'], 
                                               reduced_printing=self.reduced_printing)

        
        ## Part 5 - Build the actual simulation lattice!
        if keyword_lookup['RESTART_FILE']:

            # if  we passed a restart file then construct the lattice object using the restart file directly. Note             
            self.LATTICE   = Lattice(dimensions, chains, self.Hamiltonian, self.LATTICE_TO_ANGSTROMS, restart_object=keyword_lookup['RESTART_FILE'], hardwall=self.hardwall)
            
            # a chain crossing a periodic face cannot exist under hard walls. The
            # parser refuses a periodic restart file for a hardwall run and a
            # hardwall file is validated bond by bond, so this should never fire;
            # it used to switch the run to periodic boundaries silently, which
            # contradicted both the keyfile and keyfile_used.kf
            if self.hardwall and self.LATTICE.any_chains_straddle_boundary():
                msg = ("The restart file has a chain crossing a periodic face, which cannot exist in "
                       "a hardwall box. Run it with HARDWALL : False (or RESTART_OVERRIDE_HARDWALL : "
                       "True if the file itself is periodic).")
                pimmslogger.log_error(msg)
                raise SimulationException(msg)
        else:
            self.LATTICE   = Lattice(dimensions, chains, self.Hamiltonian, self.LATTICE_TO_ANGSTROMS, hardwall = self.hardwall )

        # Continuation (RESTART_CONTINUE): the master step counter starts at the
        # checkpoint's step, and the generator states and temperature recorded
        # in the restart file are restored just before the master loop - not
        # here, because nothing between here and the loop may draw a random
        # number after the restore. The parser has already checked the file
        # carries the state and that this run is the same system.
        self.continue_from_step = 0
        self._continue_rng = None
        self._continue_temperature = None
        if keyword_lookup.get('RESTART_CONTINUE') and keyword_lookup['RESTART_FILE']:
            _restart = keyword_lookup['RESTART_FILE']
            self.continue_from_step = int(_restart.step)
            self._continue_rng = (_restart.rng_python, _restart.rng_numpy)
            self._continue_temperature = float(_restart.temperature) if _restart.temperature is not None else None

            # A continuation has to run under the Hamiltonian of the run it
            # resumes. The checkpoint records the total energy of its own
            # configuration, so evaluating that configuration here must give
            # the same number; if it does not, the energy function differs and
            # the resumed run would be neither the original run nor a new one.
            # (Evaluating the energy draws no random numbers.)
            _start_energy = self.Hamiltonian.evaluate_total_energy(self.LATTICE)[0]
            if _start_energy != _restart.energy:
                msg = ("RESTART_CONTINUE : True cannot resume this run: under this keyfile's energy "
                       "function the restart file's configuration has energy %s, but the run that "
                       "wrote the file (%s) recorded ENERGY %s for it at step %d. The energy function "
                       "is therefore not the one that run used. The likely causes are a different "
                       "PARAMETER_FILE, a different ANGLES_OFF or NON_INTERACTING, or quench settings "
                       "that differ (QUENCH_END - or TEMPERATURE when there is no quench - is the "
                       "temperature ANGLE_PENALTY_T_NORM penalties are scaled by). Use the inputs of "
                       "the run being resumed. To start a NEW run from this configuration under a "
                       "different energy function, remove RESTART_CONTINUE."
                       % (_start_energy, getattr(_restart, 'filename', 'restart file'),
                          _restart.energy, self.continue_from_step))
                pimmslogger.log_error(msg)
                raise SimulationException(msg)

            IO_utils.status_message(
                "Continuing the run that wrote the restart file: resuming at step %d of %d "
                "with its random-number generators restored (the seed announced above is not used)"
                % (self.continue_from_step, self.n_steps), 'startup')
            pimmslogger.log_status("RESTART_CONTINUE: resuming at step %d of %d from %s"
                                   % (self.continue_from_step, self.n_steps,
                                      getattr(_restart, 'filename', 'restart file')))


        ## Part 6 - Build the Chain Temperature Switch Metropolis Monte Carlo if 
        #           this is being used (if its not being used don't even try and create the 
        #           TSMMC_coordinator object - this is because if TSMMC moves are not being used we 
        #           don't want to force the user to have sane TSMMC parameters, which creating 
        #           a TSMMC_coordinator object would required
        #           
        #           
        if self.TSMMC_USED:
            self.TSMMC_coordinator = TSMMC(temperature, 
                                           self.TSMMC_JUMP_TEMP,
                                           self.TSMMC_INTERPOLATION_MODE,
                                           self.TSMMC_STEP_MULTIPLIER,
                                           self.TSMMC_NUMBER_OF_POINTS,
                                           self.TSMMC_FIXED_OFFSET)
        else:
            self.TSMMC_coordinator = None

        ## Part 7 - Set all the custom analysis frequencies
        #
        (self.non_default_freq_analysis, self.default_freq_analysis) = self.setup_analysis(keyword_lookup)


        ## Part 8 - Finalize any special output files we want to write once all initialization has been complete
        #
        self.write_chain_to_chainid = bool(keyword_lookup['WRITE_CHAIN_TO_CHAINID'])
        if self.write_chain_to_chainid:
            # written with a plain 'w' open, which would follow a link left in the
            # directory and overwrite another run's file: replace the link instead
            IO_utils.replace_symlink(CONFIG.OUTPUT_CHAIN_TO_CHAINID)
            self.LATTICE.write_chain_to_chainid_file()

        # check freeze file,  log status, and assign frozen chains
        if keyword_lookup['FREEZE_FILE']:
            keyword_lookup['FREEZE_FILE'].validate_freeze_file(self.LATTICE)
            keyword_lookup['FREEZE_FILE'].log_freeze_file()
            self.frozen_chains = keyword_lookup['FREEZE_FILE'].chains

        # can the configured move set act on this system at all?
        self.check_moveset_applicability(keyword_lookup)

        ## Part 9 - Final logging
        pimmslogger.log_status(f'Random Seed: {random_seed}')
        pimmslogger.log_status(f'Reference-kernel seed: {random_seed%CONFIG.C_RAND_MAX}')
        pimmslogger.log_status(f'Kernel PRNG range: {CONFIG.C_RAND_MAX}')

        # describe exactly which parallel implementation this box/system gets
        if self.parallelize:
            self.report_parallelization()

        # NB: keyfile_used.kf (the configuration this run actually uses) is NOT
        # written here. It is written by startup_analysis, once the run has
        # opened its trajectory: a start that fails before that point must not
        # leave "the one to trust" describing a run that never happened beside
        # an earlier run's data.


            
            
    
       
    #-----------------------------------------------------------------
    #
    def write_effective_keyfile(self, keyword_lookup, filename=None):
        """
        Write ``keyfile_used.kf``: the configuration this run actually uses.

        The keyfile a user writes is not always what runs. ``RESTART_OVERRIDE_*``
        replace ``HARDWALL`` or ``DIMENSIONS`` with the restart file's values, a
        ``RESIZED_EQUILIBRATION`` forces hard walls for its first phase whatever
        ``HARDWALL`` says, ``EXTRA_CHAIN`` lines merge into existing chain types,
        and a missing ``SEED`` is generated at start-up. This writes the resolved
        keyword set, with a comment header recording where it came from, so a
        post-processing tool can read the truth rather than re-derive it.

        The body is written by the same writer as ``KeyFileParser.write_keyfile``,
        so it re-parses. For a restarted run it describes the equivalent system
        started afresh (the chains as a plain ``CHAIN`` list, since a restart
        object has no keyfile form); the restart file it actually started from
        is named in the header.

        Parameters
        ----------
        keyword_lookup : dict
            The resolved keyword dictionary the Simulation was built from.

        filename : str, optional
            Where to write. Default ``CONFIG.EFFECTIVE_KEYFILE_NAME``.

        Returns
        -------
        None
            The file is written. If ``filename`` is a symbolic link the link is
            replaced by a regular file, so the file it pointed at (another
            run's ``keyfile_used.kf``, say) is left as it was.
        """
        from . import keyfile_parser as _kfp
        filename = filename or CONFIG.EFFECTIVE_KEYFILE_NAME
        restart_object = keyword_lookup.get('RESTART_FILE')
        if keyword_lookup.get('__SEED_GIVEN'):
            seed_note = 'from the keyfile'
        elif restart_object and keyword_lookup.get('RESTART_CONTINUE'):
            seed_note = ('generated at start-up but unused: the generators were restored from '
                         'the restart file')
        elif restart_object:
            seed_note = ('generated at start-up; this run began from the restart file, so the '
                         'value below does not reproduce it')
        else:
            seed_note = 'generated at start-up; the value below reproduces this run'
        header = [
            "Effective configuration of this PIMMS run, written at start-up (%s, PIMMS %s)."
            % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), _pimms_version()),
            "This is what the run actually used, after every start-up resolution, and it",
            "re-parses as a keyfile. Where it differs from the keyfile you wrote, this is",
            "the one to trust.",
            "",
            "source keyfile     : %s" % keyword_lookup.get('__KEYFILE', 'unknown'),
            "SEED               : %d (%s)" % (keyword_lookup['SEED'], seed_note),
            "DIMENSIONS         : %s" % ' '.join(str(d) for d in keyword_lookup['DIMENSIONS']),
            "HARDWALL           : %s" % keyword_lookup['HARDWALL'],
        ]
        if restart_object:
            header += [
                "started from       : restart file %s%s"
                % (getattr(restart_object, 'filename', 'restart.pimms'),
                   (', written at step %d' % restart_object.step) if getattr(restart_object, 'step', None) is not None else ''),
                "                     the CHAIN lines below give that file's chain composition (EXTRA_CHAIN",
                "                     already merged, by sequence, into the existing chain types); re-running",
                "                     this keyfile starts the same system afresh, not from the snapshot",
            ]
            # the CHAIN lines are grouped by chain type; when the snapshot's chain IDs
            # are not (an EXTRA_CHAIN that joined an existing type, or a later restart
            # from such a run) a re-run numbers the chains differently
            snapshot = {}
            snapshot.update(getattr(restart_object, 'chains', None) or {})
            snapshot.update(getattr(restart_object, 'extra_chains', None) or {})
            types_in_id_order = [snapshot[cid][2] for cid in sorted(snapshot)]
            if types_in_id_order != sorted(types_in_id_order):
                header += [
                    "                     NB: this run's chain IDs are not grouped by chain type, so a re-run",
                    "                     numbers the chains differently (chain IDs, and a FREEZE_FILE keyed",
                    "                     on them, refer to different chains; chain_to_chainid.txt has the map)",
                ]
            if keyword_lookup.get('RESTART_OVERRIDE_DIMENSIONS') or keyword_lookup.get('RESTART_OVERRIDE_HARDWALL'):
                header.append("                     DIMENSIONS / HARDWALL above are the restart file's, as the overrides asked")
            if keyword_lookup.get('RESTART_CONTINUE'):
                header.append("continuation       : RESTART_CONTINUE - resumed at step %d with that file's generator state"
                              % restart_object.step)
        if keyword_lookup.get('RESIZED_EQUILIBRATION'):
            header.append("resized equilibrn. : steps 1 to %d ran in box %s under HARDWALL : True (always, whatever"
                          % (keyword_lookup['EQUILIBRATION'],
                             ' '.join(str(d) for d in keyword_lookup['RESIZED_EQUILIBRATION'])))
            header.append("                     HARDWALL says); the box above and HARDWALL above apply from step %d"
                          % (keyword_lookup['EQUILIBRATION'] + 1))
        if keyword_lookup.get('FREEZE_FILE'):
            header.append("FREEZE_FILE        : %s" % getattr(keyword_lookup['FREEZE_FILE'], 'filename', keyword_lookup['FREEZE_FILE']))
        header.append("")
        # the body must be literal: the overrides have been applied to the values
        # above, so writing them as True would send a re-parse back to a file
        # that may no longer exist
        resolved = dict(keyword_lookup)
        resolved['RESTART_OVERRIDE_DIMENSIONS'] = False
        resolved['RESTART_OVERRIDE_HARDWALL'] = False
        resolved['RESTART_CONTINUE'] = False
        # a keyfile_used.kf linked in from another run's directory must not be
        # written through: replace the link, not the file it points at
        IO_utils.replace_symlink(filename)
        _kfp.write_keyword_lookup(resolved, filename, header_lines=header)

    #-----------------------------------------------------------------
    #
    def report_parallelization(self, note=''):
        """
        Print (and log) a detailed description of the parallel implementation this
        run will actually use.

        Emitted at startup when ``PARALLELIZE`` is set, and again after a resized
        equilibration swaps in the production box (the block decomposition depends
        on the box). Covers the thread budget and whether the compiled kernels have
        OpenMP at all; the per-bead crankshaft decomposition (halo width, block
        grid, block size, fraction of the box movable per sweep, single-block
        warning); and, for each enabled whole-chain move (slither, pull), the
        chain-level decomposition and how the chains are split between the two
        kernels. That split is a run constant taken from chain LENGTHS alone -
        both kernels run on every megamove, the parallel one over the chains short
        enough to fit a block interior and the serial one over the rest - so the
        report describes a partition, not a per-move choice.

        The report also warns, per move, when the megamove is too small to pay for
        the parallel kernel. Entering a parallel kernel has a fixed cost that the
        serial kernel does not have (about 0.1 ms plus 0.1 microseconds per bead;
        0.6 ms at 5,000 beads), so at the default ``CRANKSHAFT_SUBSTEPS : 500``
        (about 0.05 ms of sampling) ``PARALLELIZE`` makes the crankshaft several
        times slower, not faster. The warning is issued when the most the
        threads could save, ``sub-moves x cost per sub-move x (1 - 1/threads)``,
        is below that fixed cost; it names the ``*_SUBSTEPS`` keyword and the
        value to raise it to (ten times the fixed cost in sampling work). The
        rule uses only run constants - substeps, bead count, block layout and
        thread count - and the cost figures are approximate. It is advice only:
        the kernel that runs is never changed by it.

        Parameters
        ----------
        note : str, optional
            Context appended to the report header (e.g. why the report is being
            re-issued). Default is an empty string, which prints the bare header.

        Returns
        -------
        list of str
            The report lines (also printed via ``IO_utils.status_message`` and
            written to the log). A line that starts with ``WARNING:`` after its
            indent is issued as a warning rather than as a status line.
        """
        dims = list(self.LATTICE.dimensions)
        n_dim = len(dims)
        omp = mega_crank_fast.openmp_info()
        idx_to_bead, sorted_chains, chain_offset, chain_length, chain_homo = \
            moves.parallel_chain_metadata(self.LATTICE)
        has_LR = bool(np.any(np.asarray(idx_to_bead)[:, 1] == 1))
        frozen = list(self.frozen_chains) if self.frozen_chains else []
        n_frozen_beads = int(np.isin(np.asarray(idx_to_bead)[:, 4], frozen).sum()) if frozen else 0
        n_beads = int(len(idx_to_bead))
        # the CPUs this process may run on (affinity-aware), which is what
        # PARALLEL_THREADS : 0 resolves against - not the machine's core count
        from . import keyfile_parser
        cores = keyfile_parser.available_cpu_count()

        L = []
        head = 'PARALLELIZATION REPORT' + (' - ' + note if note else '')
        L.append(head)
        L.append('-' * len(head))
        L.append('Box: %s (%dD, %s); %d chains, %d beads; interaction radius %d (%s)'
                 % ('x'.join(str(d) for d in dims), n_dim,
                    'hardwall' if self.hardwall else 'periodic',
                    len(sorted_chains), n_beads, 3 if has_LR else 1,
                    'long-range beads present' if has_LR else 'short-range only'))
        if frozen:
            L.append('Frozen chains: %d (%d beads) - excluded from every parallel move, kept as fixed obstacles'
                     % (len(frozen), n_frozen_beads))
        req = self.keyword_lookup.get('PARALLEL_THREADS', 0) if hasattr(self, 'keyword_lookup') else None
        # PARALLEL_THREADS : 0 means OMP_NUM_THREADS when that is set, otherwise
        # every CPU available to the process (keyfile_parser.resolve_parallel_threads)
        L.append('Threads: %d OpenMP threads per parallel megamove (%s; %d CPUs available to this process)'
                 % (self.parallel_threads,
                    'PARALLEL_THREADS : 0 -> OMP_NUM_THREADS if set, otherwise every available CPU'
                    if (req is None or int(req) <= 0) else 'from PARALLEL_THREADS', cores))
        if omp['enabled']:
            L.append('OpenMP: compiled in (runtime default thread budget %d)' % omp['max_threads'])
        else:
            L.append('OpenMP: NOT compiled into pimms.mega_crank_fast - the "parallel" kernels will '
                     'execute their blocks serially (no speed-up; sampling unaffected). Rebuild with '
                     'OpenMP (macOS: brew install libomp) to use the threads.')

        # --- crankshaft: per-bead frozen-halo decomposition ---
        # the parallel crankshaft runs for MOVE_CRANKSHAFT, and for the sub-moves of
        # a system-wide TSMMC excursion when no non-TSMMC move is enabled (the
        # excursion then falls back to crankshaft megamoves)
        kl = self.keyword_lookup if hasattr(self, 'keyword_lookup') else {}

        def _enabled(keyword):
            """Whether a move keyword has a positive frequency in this run.

            Parameters
            ----------
            keyword : str
                A ``MOVE_*`` keyword.

            Returns
            -------
            bool
                True when the keyword is present and above zero.
            """
            return float(kl.get(keyword, 0) or 0) > 0

        # --- is a megamove big enough to pay for the parallel kernel? ---
        # Entering a parallel kernel costs a fixed amount before the first sub-move
        # (bucketing every bead into its block, drawing the block shift, starting
        # the threads): about 0.1 ms plus 0.1 us per bead, i.e. 0.6 ms at 5,000
        # beads, where the serial kernel costs next to nothing. A sub-move costs
        # about 0.1 us (crankshaft) or 0.4 us (slither or pull of a short chain) in
        # either kernel. These are measured, approximate and machine dependent.
        # With T threads the most a megamove of A sub-moves can save is
        # A x cost x (1 - 1/T), so below A = fixed / (cost x (1 - 1/T)) the
        # parallel kernel is slower than the serial one however the beads are
        # arranged. Everything in the rule is a run constant (substeps, bead count,
        # block layout, threads) - it only ever produces a line in this report and
        # never changes which kernel runs.
        fixed_ms = 0.1 + 1.0e-4 * n_beads
        one_thread_warned = []

        def _megamove_advice(keyword, substeps, attempts_per_substep, attempt_us, num_blocks):
            """Warn when a megamove is too small to pay for the parallel kernel.

            Parameters
            ----------
            keyword : str
                The substep keyword to name (``CRANKSHAFT_SUBSTEPS``,
                ``SLITHER_SUBSTEPS`` or ``PULL_SUBSTEPS``).

            substeps : int
                Its value in this run.

            attempts_per_substep : int
                Sub-moves the parallel kernel is given per unit of ``substeps``:
                1 for the crankshaft, the number of parallel-side chains for
                slither and pull.

            attempt_us : float
                Approximate cost of one sub-move in microseconds.

            num_blocks : int
                Blocks in the decomposition this kernel uses; no more threads
                than blocks can work at once.

            Returns
            -------
            list of str
                No lines when the megamove is large enough, otherwise one
                indented ``WARNING`` line naming the keyword and a value to
                raise it to.
            """
            import math
            threads = min(int(self.parallel_threads), int(num_blocks)) if omp['enabled'] else 1
            if attempts_per_substep <= 0 or num_blocks <= 1:
                return []
            if threads <= 1:
                # one condition, one warning: said under the first move it applies to
                if one_thread_warned:
                    return []
                one_thread_warned.append(keyword)
                # without OpenMP no PARALLEL_THREADS value helps: the blocks run
                # one after another whatever is asked for
                return ['    WARNING: %s, so every parallel megamove pays the parallel kernels\' '
                        'fixed cost (about %.2g ms here; approximate) for no gain - %s or drop '
                        'PARALLELIZE.'
                        % (('only one thread will run the parallel kernels', fixed_ms,
                            'set PARALLEL_THREADS above 1') if omp['enabled'] else
                           ('the kernels were built without OpenMP and run their blocks one after '
                            'another', fixed_ms, 'rebuild with OpenMP'))]
            per_substep_ms = 1.0e-3 * attempt_us * attempts_per_substep
            # the most the threads can take off a megamove's sampling time
            saving_fraction = 1.0 - 1.0 / threads
            break_even = fixed_ms / (per_substep_ms * saving_fraction)
            if substeps >= break_even:
                return []

            def _round_up(value):
                """Round a positive number up to one significant figure.

                Parameters
                ----------
                value : float
                    The number to round.

                Returns
                -------
                int
                    The rounded value, at least 1.
                """
                if value <= 1:
                    return 1
                scale = 10 ** int(math.floor(math.log10(value)))
                return int(math.ceil(value / scale - 1e-9) * scale)

            return ['    WARNING: megamove too small for the parallel kernel - at %s : %d, %d threads '
                    'can save at most about %.2g ms per megamove (%d%% of %.2g ms of sampling), '
                    'against a fixed cost of about %.2g ms for entering the parallel kernel (0.1 ms '
                    '+ 0.1 us per bead, %d beads), so this move runs slower than it would without '
                    'PARALLELIZE. Raise %s to about %d or more (break-even is near %d), or drop '
                    'PARALLELIZE. These costs are approximate and machine dependent.'
                    % (keyword, substeps, threads, per_substep_ms * substeps * saving_fraction,
                       int(round(100.0 * saving_fraction)), per_substep_ms * substeps, fixed_ms,
                       n_beads, keyword, _round_up(10.0 * fixed_ms / per_substep_ms),
                       _round_up(break_even))]

        non_tsmmc = [kw for kw in ('MOVE_CRANKSHAFT', 'MOVE_CHAIN_TRANSLATE', 'MOVE_CHAIN_ROTATE',
                                   'MOVE_CHAIN_PIVOT', 'MOVE_HEAD_PIVOT', 'MOVE_SLITHER',
                                   'MOVE_CLUSTER_TRANSLATE', 'MOVE_CLUSTER_ROTATE', 'MOVE_PULL',
                                   'MOVE_JUMP_AND_RELAX', 'MOVE_VMMC') if _enabled(kw)]
        crank_used = _enabled('MOVE_CRANKSHAFT') or (_enabled('MOVE_SYSTEM_TSMMC') and not non_tsmmc)
        ci = mega_crank_fast.parallel_crank_layout_info(dims[0], dims[1], dims[2] if n_dim == 3 else 1, has_LR)
        kname = 'mega_crank_parallel' + ('' if n_dim == 3 else '_2D')
        if not crank_used:
            L.append('Crankshaft (MOVE_CRANKSHAFT): not in the move set')
        else:
            L.append('Crankshaft (MOVE_CRANKSHAFT): kernel %s - per-bead frozen-halo checkerboard' % kname)
            if ci['num_blocks'] == 1:
                L.append('    halo W=%d; box does not split (needs >= %d sites in a dimension): ONE block -> '
                         'runs single-threaded, equivalent to the serial kernel (no speed-up)'
                         % (ci['W'], 16 * ci['W']))
            else:
                L.append('    halo W=%d; block grid %s = %d blocks of %s sites; %.0f%% of the box movable per sweep '
                         '(random block shift every sweep)'
                         % (ci['W'], 'x'.join(str(b) for b in ci['blocks'][:n_dim]), ci['num_blocks'],
                            'x'.join(str(b) for b in ci['block_size'][:n_dim]), 100.0 * ci['movable_fraction']))
                L.extend(_megamove_advice('CRANKSHAFT_SUBSTEPS', int(self.CS_substeps), 1, 0.1,
                                          ci['num_blocks']))

        # --- whole-chain moves: chain-level decomposition + fit gate ---
        for keyword, label, cap_mode, kernel, min_length in (
                ('MOVE_SLITHER', 'Slither (MOVE_SLITHER)', 'hetero', 'mega_slither_parallel', 1),
                ('MOVE_PULL', 'Pull (MOVE_PULL)', 'all', 'mega_pull_parallel', 3)):
            freq = self.keyword_lookup.get(keyword, 0) if hasattr(self, 'keyword_lookup') else 0
            if not freq or float(freq) <= 0:
                L.append('%s: not in the move set' % label)
                continue
            rep = moves.parallel_chain_fit_report(idx_to_bead, chain_offset, chain_length, dims, has_LR,
                                                  chain_homo=chain_homo, cap_mode=cap_mode,
                                                  frozen_chains=frozen, min_length=min_length)
            lay = rep['layout']
            kname = kernel + ('' if n_dim == 3 else '_2D')
            if lay['num_blocks'] == 1:
                # one block means there is no parallel work to hand out, so the
                # partition is empty and every chain runs on the serial kernel
                L.append('%s: box does not split for the chain-level halo (W=%d, needs >= %d sites): '
                         'every chain runs on the SERIAL kernel' % (label, lay['W'], 8 * lay['W']))
                continue
            # Both kernels run on every megamove: the chains are partitioned ONCE,
            # by length, and the parallel kernel takes the short ones while the
            # serial kernel takes the rest. Nothing is chosen per megamove from the
            # current configuration - doing that broke stationarity, so this report
            # must not describe a per-move choice either.
            L.append('%s: kernel %s - chain-level halo W=%d, block grid %s (%s sites, interiors %s); '
                     'parallel kernel for %d chain(s) of length <= %d, serial kernel for %d longer '
                     'chain(s) - both run every megamove'
                     % (label, kname, lay['W'], 'x'.join(str(b) for b in lay['blocks'][:n_dim]),
                        'x'.join(str(b) for b in lay['block_size'][:n_dim]),
                        'x'.join(str(v) for v in rep['interiors']),
                        rep['n_parallel'], rep['interior'], rep['n_serial']))
            if rep['n_over_cap']:
                L.append('    (%d of the serial-side chain(s) exceed the 512-bead kernel buffer)'
                         % rep['n_over_cap'])
            if rep['n_too_short']:
                L.append('    (%d chain(s) shorter than %d beads are never moved by this move)'
                         % (rep['n_too_short'], min_length))
            L.append('    (the split is fixed by chain length for the whole run - it is never '
                     're-decided from the configuration)')
            # the parallel kernel's budget is substeps x (chains on the parallel side)
            L.extend(_megamove_advice(
                keyword.replace('MOVE_', '') + '_SUBSTEPS',
                int(self.slither_substeps if keyword == 'MOVE_SLITHER' else self.pull_substeps),
                rep['n_parallel'], 0.4, lay['num_blocks']))

        L.append('Note: the parallel kernels sample the SAME equilibrium as the serial ones but relax more '
                 'slowly per step (only block interiors move each sweep) - judge equilibration by the '
                 'observable plateau, not by step count.')

        IO_utils.newline()
        for line in L:
            # a WARNING line goes out as a warning, in place, so it carries the
            # usual prefix on screen and in log.txt and still reads as part of
            # the report
            if line.lstrip().startswith('WARNING: '):
                text = 'PARALLELIZE: ' + line.strip()[len('WARNING: '):]
                IO_utils.status_message(text, 'warning')
                pimmslogger.log_warning(text)
            else:
                IO_utils.status_message(line, 'startup')
                pimmslogger.log_status(line)
        IO_utils.newline()
        return L


    #-----------------------------------------------------------------
    #
    def check_moveset_applicability(self, keyword_lookup):
        """
        Refuse to start a run whose move set can never move anything, and warn
        about draws that are guaranteed null moves.

        A run in which every chain is frozen, or in which every enabled move
        needs longer chains than the system contains (e.g. dimers with only
        ``MOVE_CHAIN_PIVOT``, or a pull-only move set with no chain of three or
        more beads), used to complete normally: ENERGY.dat, RG.dat and traj.xtc
        were all written from a configuration that never changed, with no
        warning anywhere. Both cases now raise at construction.

        A third refusal covers move sets in which only cluster moves can act -
        judged on the moves long enough for some mobile chain, so an enabled
        pivot or pull that no chain is long enough for does not count. A
        cluster move is rejected when the cluster holds a frozen chain, when it
        holds every chain of a multi-chain system, and (rotation) when the
        cluster is one isolated bead; and because cluster moves reject merges
        and nothing splits a cluster, which chains touch which never changes. So
        if every mobile chain starts in such a cluster, nothing can ever move and
        the run is refused. The clusters are labelled once, one search per
        cluster. The verdict depends on the starting placement, and the message
        says so. Under ``RESIZED_EQUILIBRATION`` the box and its boundary change
        mid-run, so this is only a warning there.

        Everything else is a warning, one per condition:

        * the move draw is independent of the chain draw, so a rotation or head
          pivot (two beads or more) or a pivot (three or more) drawn for a
          shorter chain is a rejected null move, and a pull megamove proposes
          nothing when no mobile chain has three beads. The warning names each
          such move, the chain lengths it cannot act on and the fraction of
          steps wasted, which is ``MOVE_x`` times the fraction of mobile chains
          that are too short (the whole of ``MOVE_PULL``). Since 1.0.8 such
          draws are no longer remapped to a crankshaft megamove - that made the
          executed move mix composition dependent;
        * mobile chains too short for EVERY enabled move never move at all, even
          though longer chains keep the run alive;
        * a cluster-only move set freezes the contact graph, and may leave some
          chains unable to move (see above);
        * ``MOVE_CLUSTER_ROTATE`` gets its own sentence: it CAN move a monomer,
          it just cannot move an isolated one;
        * the move set leaves a chain's conformation invariant - see
          :func:`conformation_freezing_warnings`. That is a failure of
          irreducibility rather than of detailed balance, and it is silent in
          ENERGY.dat and in the trajectory, so it needs saying out loud.

        Nothing here touches the lattice or draws a random number.

        Parameters
        ----------
        keyword_lookup : dict
            The completed keyfile configuration dictionary. Only the ``MOVE_*``
            frequencies are read from it here.

        Returns
        -------
        dict
            What was found, so that the numbers in the warnings can be checked:
            ``null_fractions`` (``MOVE_*`` keyword -> fraction of steps that are
            guaranteed null moves), ``null_fraction`` (their sum),
            ``excursion_null_fraction`` (the same for the sub-moves of a
            ``MOVE_SYSTEM_TSMMC`` excursion, 0.0 without one),
            ``unmoved_lengths`` (mobile chain lengths no enabled move can act
            on), ``immobile_chains`` (chainIDs a cluster-only move set can never
            move) and ``warnings`` (the messages issued, in order).

        Raises
        ------
        SimulationException
            If every chain is frozen, if no unfrozen chain is long enough for
            any enabled move, or if only cluster moves are enabled and none of
            them can ever move any mobile chain.
        """
        # minimum chain length each move needs to do anything (see the guards in
        # moves.py: chain_rotate < 2, chain_pivot < 3, head_pivot < 2 reject;
        # pull needs an interior bead). Every other move can act on a monomer.
        min_length = {'MOVE_CRANKSHAFT': 1, 'MOVE_CHAIN_TRANSLATE': 1, 'MOVE_CHAIN_ROTATE': 2,
                      'MOVE_CHAIN_PIVOT': 3, 'MOVE_HEAD_PIVOT': 2, 'MOVE_SLITHER': 1,
                      'MOVE_CLUSTER_TRANSLATE': 1, 'MOVE_CLUSTER_ROTATE': 1, 'MOVE_CTSMMC': 1,
                      'MOVE_MULTICHAIN_TSMMC': 1, 'MOVE_PULL': 3, 'MOVE_SYSTEM_TSMMC': 1,
                      'MOVE_JUMP_AND_RELAX': 1, 'MOVE_VMMC': 1}
        enabled = [kw for kw in min_length if keyword_lookup.get(kw, 0) > 0]

        # a system-wide TSMMC excursion runs the OTHER enabled non-TSMMC moves as
        # its sub-moves (crankshaft only when there are none), so it can do exactly
        # what they can and no more: it needs their shortest minimum length, and it
        # adds no conformational freedom of its own
        tsmmc_moves = ('MOVE_CTSMMC', 'MOVE_MULTICHAIN_TSMMC', 'MOVE_SYSTEM_TSMMC')
        excursion_moves = [kw for kw in enabled if kw not in tsmmc_moves] or ['MOVE_CRANKSHAFT']
        min_length['MOVE_SYSTEM_TSMMC'] = min(min_length[kw] for kw in excursion_moves)

        frozen = set(self.frozen_chains) if self.frozen_chains else set()
        mobile_lengths = [len(chain.get_ordered_positions())
                          for chainID, chain in self.LATTICE.chains.items()
                          if chainID not in frozen]

        if len(mobile_lengths) == 0:
            msg = ("Every chain in the system is frozen (FREEZE_FILE), so no move can change "
                   "the configuration - refusing to run a simulation whose output would be a "
                   "constant configuration repeated %i times" % self.n_steps)
            pimmslogger.log_error(msg)
            raise SimulationException(msg)

        longest = max(mobile_lengths)
        usable = [kw for kw in enabled if min_length[kw] <= longest]
        if enabled and not usable:
            msg = ("No enabled move can act on this system: the longest unfrozen chain has %i "
                   "bead(s) but %s need(s) at least %s. Every step would be a rejected null "
                   "move and the output would describe an unchanging configuration."
                   % (longest, ', '.join(enabled),
                      ', '.join('%i (%s)' % (min_length[kw], kw) for kw in enabled)))
            pimmslogger.log_error(msg)
            raise SimulationException(msg)

        summary = {'null_fractions': {}, 'null_fraction': 0.0, 'excursion_null_fraction': 0.0,
                   'unmoved_lengths': [], 'immobile_chains': [], 'warnings': []}

        def _warn(msg):
            """Issue one start-up warning on screen and in the log, and record it.

            Parameters
            ----------
            msg : str
                The warning text.

            Returns
            -------
            None
            """
            IO_utils.status_message(msg, 'warning')
            pimmslogger.log_warning(msg)
            summary['warnings'].append(msg)

        n_mobile = len(mobile_lengths)
        length_counts = {}
        for L in mobile_lengths:
            length_counts[L] = length_counts.get(L, 0) + 1

        def _shorter_than(need):
            """Describe the mobile chains with fewer than ``need`` beads.

            Parameters
            ----------
            need : int
                The minimum chain length a move needs.

            Returns
            -------
            tuple of (int, str)
                How many mobile chains are shorter than ``need``, and their
                lengths spelled out, e.g. ``'3 of 1 bead, 2 of 2 beads'``.
            """
            short = [L for L in sorted(length_counts) if L < need]
            return (sum(length_counts[L] for L in short),
                    ', '.join('%i of %i bead%s' % (length_counts[L], L, '' if L == 1 else 's')
                              for L in short))

        # every move that can change the configuration. A system excursion only
        # ever runs its sub-moves, so it stands for them and nothing more
        movers = set(kw for kw in enabled if kw != 'MOVE_SYSTEM_TSMMC')
        if 'MOVE_SYSTEM_TSMMC' in enabled:
            movers.update(excursion_moves)

        # --- cluster-only move sets ---------------------------------------
        # cluster_translate / cluster_rotate (moves.py) reject a cluster that holds
        # a frozen chain, a cluster the component search refuses at the size
        # threshold the main loop passes (number of chains - 1: every chain of a
        # multi-chain system), and - rotation only - a cluster that is one isolated
        # bead, which every rotation maps onto itself. An accepted cluster move
        # must leave the cluster's membership unchanged and nothing else can act to
        # split one, so the clusters found here are the clusters of the whole run.
        #
        # "Cluster-only" is judged on the moves that can ACT, not on the moves that
        # are enabled: a move too short-handed for every mobile chain (a pivot or a
        # pull on dimers, a rotation on monomers) changes nothing, so it must not
        # hide a frozen contact graph or a system that can never move.
        cluster_moves = ('MOVE_CLUSTER_TRANSLATE', 'MOVE_CLUSTER_ROTATE')
        acting = set(kw for kw in movers if min_length[kw] <= longest)
        if acting and acting <= set(cluster_moves):
            n_chains = self.LATTICE.get_number_of_chains()
            stuck = {'frozen': [], 'whole': [], 'isolated': []}
            # label each cluster ONCE: a search is started only from a chain no
            # earlier search reached, so the cost is one pass over the system
            # rather than one search per chain (quadratic in a percolated box).
            # No size threshold is passed; the main loop's threshold (number of
            # chains - 1) refuses exactly the cluster that holds every chain of a
            # system of two or more, which is tested for directly below.
            component_of = {}
            for chainID, chain in self.LATTICE.chains.items():
                if chainID in frozen:
                    continue
                if chainID not in component_of:
                    component = frozenset(lattice_utils.get_all_chains_in_connected_component(
                        chainID, self.LATTICE.grid, self.LATTICE.chains, threshold=None,
                        useChains=True, hardwall=self.hardwall))
                    for member in component:
                        component_of[member] = component
                component = component_of[chainID]
                if n_chains > 1 and len(component) == n_chains:
                    stuck['whole'].append(chainID)
                elif frozen.intersection(component):
                    stuck['frozen'].append(chainID)
                elif ('MOVE_CLUSTER_TRANSLATE' not in acting and len(component) == 1
                      and len(chain.get_ordered_positions()) == 1):
                    stuck['isolated'].append(chainID)

            reasons = []
            if stuck['frozen']:
                reasons.append('%i in a cluster that contains a frozen chain (FREEZE_FILE)'
                               % len(stuck['frozen']))
            if stuck['whole']:
                reasons.append('%i in a cluster that holds every chain in the system'
                               % len(stuck['whole']))
            if stuck['isolated']:
                reasons.append('%i isolated single bead(s), which every rotation maps onto itself'
                               % len(stuck['isolated']))
            immobile = sorted(stuck['frozen'] + stuck['whole'] + stuck['isolated'])
            summary['immobile_chains'] = immobile
            only = ', '.join(kw for kw in cluster_moves if kw in acting)
            # enabled moves that are too short-handed for every mobile chain
            dead = [kw for kw in min_length if kw in movers and kw not in acting]
            if dead:
                only += '; %s cannot act on any mobile chain' % ', '.join(dead)
            remedy = ("Add a move that acts on single chains (MOVE_CRANKSHAFT, MOVE_SLITHER or "
                      "MOVE_CHAIN_TRANSLATE).")

            if len(immobile) == n_mobile and not getattr(self, 'resize_eq', False):
                msg = ("No enabled move can ever change this configuration. Only cluster moves can "
                       "act on this system (%s) and they cannot move any of the %i mobile "
                       "chain(s): %s. Cluster moves never merge or split clusters, so this holds "
                       "for the whole run and the output would describe an unchanging "
                       "configuration. This is decided from the starting configuration, so "
                       "another SEED may start - but the contact graph is frozen either way. %s"
                       % (only, n_mobile, '; '.join(reasons), remedy))
                pimmslogger.log_error(msg)
                raise SimulationException(msg)

            msg = ("Only cluster moves can act on this system (%s). A cluster move is rejected if "
                   "it would merge two clusters and nothing here splits one, so the contact graph "
                   "is frozen: which chains touch which at start-up is fixed for the whole run, "
                   "and no inter-chain contact will ever form or break." % only)
            if immobile:
                msg += (" In addition %i of %i mobile chain(s) can never move: %s."
                        % (len(immobile), n_mobile, '; '.join(reasons)))
            _warn(msg + ' ' + remedy)

        # --- guaranteed null moves ----------------------------------------
        # the chain and the move are drawn independently (uniform over the mobile
        # chains; MOVE_* fractions), so a per-chain move lands on a chain too short
        # for it with probability (fraction of mobile chains that are too short).
        # A pull megamove sweeps every mobile chain of three or more beads and
        # proposes nothing when there is none; an excursion whose sub-moves cannot
        # act on any mobile chain is a whole excursion of null sub-moves.
        total_freq = float(sum(keyword_lookup[kw] for kw in enabled))
        non_tsmmc_freq = float(sum(keyword_lookup[kw] for kw in enabled if kw not in tsmmc_moves))
        clauses = []
        excursion_null = 0.0
        for kw in enabled:
            share = float(keyword_lookup[kw]) / total_freq
            if kw in ('MOVE_CHAIN_ROTATE', 'MOVE_CHAIN_PIVOT', 'MOVE_HEAD_PIVOT'):
                (n_short, described) = _shorter_than(min_length[kw])
                if n_short == 0:
                    continue
                null_share = float(n_short) / n_mobile
                clauses.append('%s cannot move %s: %i of %i mobile chains (%s), %.1f%% of steps'
                               % (kw, 'a single bead' if min_length[kw] == 2
                                  else 'a chain of fewer than %i beads' % min_length[kw],
                                  n_short, n_mobile, described, 100.0 * share * null_share))
            elif kw == 'MOVE_PULL' and longest < min_length[kw]:
                null_share = 1.0
                clauses.append('MOVE_PULL needs a chain of %i or more beads and no mobile chain has '
                               'one, so every pull megamove proposes nothing: %.1f%% of steps'
                               % (min_length[kw], 100.0 * share))
            elif kw == 'MOVE_SYSTEM_TSMMC' and longest < min_length[kw]:
                summary['null_fractions'][kw] = share
                clauses.append('MOVE_SYSTEM_TSMMC excursions run only %s, which cannot act on any '
                               'mobile chain: %.1f%% of steps'
                               % ('/'.join(excursion_moves), 100.0 * share))
                continue
            else:
                continue
            summary['null_fractions'][kw] = share * null_share
            excursion_null += null_share * float(keyword_lookup[kw]) / non_tsmmc_freq

        if clauses:
            summary['null_fraction'] = sum(summary['null_fractions'].values())
            msg = ("About %.0f%% of steps will be rejected null moves, because a move is drawn "
                   "independently of the chain it lands on. %s. Sampling stays correct; lower "
                   "those MOVE_* fractions to stop wasting the steps."
                   % (100.0 * summary['null_fraction'], '; '.join(clauses)))
            if 'MOVE_SYSTEM_TSMMC' in usable and excursion_null > 0:
                summary['excursion_null_fraction'] = excursion_null
                msg += (" The sub-moves of a MOVE_SYSTEM_TSMMC excursion are drawn the same way: "
                        "%.1f%% of them are null too." % (100.0 * excursion_null))
            _warn(msg)

        # --- mobile chains nothing can act on ------------------------------
        # the refusal above looks at the longest mobile chain only; shorter ones
        # can still be out of reach of every enabled move
        unmoved = [L for L in sorted(length_counts)
                   if not any(min_length[kw] <= L for kw in movers)]
        if unmoved:
            summary['unmoved_lengths'] = unmoved
            (n_short, described) = _shorter_than(max(unmoved) + 1)
            _warn("No enabled move can act on %i of %i mobile chains (%s): %s. Those chains will "
                  "sit where they were placed at start-up for the whole run, as fixed obstacles. "
                  "MOVE_CRANKSHAFT, MOVE_SLITHER and MOVE_CHAIN_TRANSLATE act on a chain of any "
                  "length."
                  % (n_short, n_mobile, described,
                     ', '.join('%s needs %i or more beads' % (kw, min_length[kw])
                               for kw in min_length if kw in movers)))

        n_monomers = length_counts.get(1, 0)

        # MOVE_CLUSTER_ROTATE needs its own sentence rather than a place in the
        # null-move warning above: the message there says a move "cannot move a
        # single bead", which is not what happens here. A cluster rotation of an
        # ISOLATED monomer is a perfectly legal rotation that maps the bead onto
        # itself, and only draws landing on an isolated monomer are null - a
        # monomer with any neighbour belongs to a bigger cluster that rotates
        # normally, so the null fraction is not the monomer fraction.
        if n_monomers and 'MOVE_CLUSTER_ROTATE' in enabled:
            msg = ("%i of %i mobile chains are single beads and MOVE_CLUSTER_ROTATE is enabled; "
                   "rotating a cluster that is a single ISOLATED bead maps it exactly onto itself, "
                   "so those draws are rejected as null moves (they used to be accepted with zero "
                   "energy change, which roughly doubled this move's column in ACCEPTANCE.dat in a "
                   "monomer-rich box). A monomer that touches another chain is part of a "
                   "larger cluster and rotates normally, so this is not every monomer draw, and "
                   "MOVE_FREQS.dat is unaffected - a null draw is still a genuine attempt."
                   % (n_monomers, len(mobile_lengths)))
            _warn(msg)

        # move sets under which a chain's conformation is a conserved quantity;
        # a system excursion reshapes chains only through its sub-moves
        shape_check = [kw for kw in usable if kw != 'MOVE_SYSTEM_TSMMC']
        if 'MOVE_SYSTEM_TSMMC' in usable:
            shape_check += [kw for kw in excursion_moves
                            if min_length[kw] <= longest and kw not in shape_check]
        for msg in conformation_freezing_warnings(shape_check, longest):
            _warn(msg)

        return summary


    #-----------------------------------------------------------------
    #
    def run_simulation(self):
        """Run the simulation, and leave the directory in a known state however it ends.

        The implementation lives in :meth:`_run_simulation`; this guard is the
        run's lifecycle. It

        * marks the working directory as in use for the duration of the run and
          warns if another run already seems to be using it
          (``restart.claim_run_directory``);
        * turns ``SIGTERM`` (what a scheduler or ``kill`` sends) into a clean
          stop, so that the clean-up below runs rather than the process dying
          with trajectory frames still buffered. The handler is installed on
          entry and the previous one restored on exit, only on the main thread
          (Python allows signal handlers nowhere else) and only if ``SIGTERM``
          still had its default action - a handler an embedding program
          installed is left alone. ``SIGINT`` goes through the same handler
          (again only where it still had Python's standard one), which raises
          ``KeyboardInterrupt`` as before; sharing the handler is what lets the
          first of two signals that arrive together be the only one acted on;
        * always releases the incremental XTC writer, so exceptions from a move,
          analysis routine or output path do not leave the persistent trajectory
          handle open, and under ``SAVE_AT_END`` writes the frames buffered so
          far if the run did not finish;
        * records in ``log.txt`` (and, for an interrupt, on screen) how the run
          ended when it did not finish: at which step it stopped and which
          checkpoint is on disk;
        * ignores ``SIGINT`` and ``SIGTERM`` while that clean-up is going on
          (main thread only), so a second signal cannot cut it short and lose
          the stop line or the buffered trajectory, and puts both dispositions
          back as they were before returning, on every path.

        No checkpoint is written on the way out. A signal can arrive in the
        middle of a move, when the lattice is not a state of the Markov chain;
        the last periodic checkpoint is the one to resume from.

        Returns
        -------
        None

        Raises
        ------
        KeyboardInterrupt
            On ``SIGINT`` (Ctrl-C), re-raised unchanged after the clean-up.

        SystemExit
            On ``SIGTERM``, with the conventional status 143 (128 + 15), after
            the clean-up.
        """
        import signal
        import threading

        active_exception = False
        buffer_already_written = False
        self._last_completed_step = None
        self._last_checkpoint_step = None
        self._stopped_by_signal = None

        self._last_analysed_step = None
        self._stop_in_progress = False

        # --- SIGTERM -> clean stop, for the duration of the run only; and SIGINT
        # through a handler of ours that does what Python's own does (raise
        # KeyboardInterrupt) but takes part in the "one stop only" rule below
        installed_handlers = {}
        on_main_thread = threading.current_thread() is threading.main_thread()
        if on_main_thread:
            def _stop_on_signal(signum, frame):
                """
                Turn SIGTERM or SIGINT into an exception so the run unwinds and
                cleans up - once.

                The first signal records that a stop is under way and raises.
                Any later call returns without doing anything: two signals that
                arrive together (a ``kill`` and a Ctrl-C, a scheduler that sends
                both) are handled one after the other, and the second used to
                raise a fresh exception in the first lines of the clean-up,
                before the stop had been reported.

                Parameters
                ----------
                signum : int
                    The signal number delivered (``signal.SIGTERM`` or
                    ``signal.SIGINT``).

                frame : frame or None
                    The interrupted stack frame (unused).

                Raises
                ------
                SystemExit
                    For the first signal, if it is ``SIGTERM``, with status
                    ``128 + signum``.

                KeyboardInterrupt
                    For the first signal, if it is ``SIGINT`` (exactly what
                    Python's default handler does).
                """
                if self._stop_in_progress:
                    return
                self._stop_in_progress = True
                if signum == signal.SIGINT:
                    raise KeyboardInterrupt
                self._stopped_by_signal = 'SIGTERM'
                raise SystemExit(128 + signum)

            # only where the signal still has its standard behaviour: a handler an
            # embedding program installed (or a signal it ignores) is left alone
            for _signum, _standard in ((signal.SIGTERM, signal.SIG_DFL),
                                       (signal.SIGINT, signal.default_int_handler)):
                try:
                    _previous = signal.getsignal(_signum)
                    if _previous is _standard:
                        signal.signal(_signum, _stop_on_signal)
                        installed_handlers[_signum] = _previous
                except (ValueError, OSError):
                    pass

        # --- no second stop during the clean-up
        held_signals = {}

        def _hold_stop_signals():
            """
            Ignore SIGINT and SIGTERM until the clean-up is over.

            The clean-up closes the trajectory, writes the SAVE_AT_END buffer
            and records where the run stopped. A second signal arriving in the
            middle of it (two ``kill`` commands in a row, a scheduler that
            repeats its SIGTERM, an impatient second Ctrl-C) used to raise a
            fresh exception there, and the stop line - or the buffered
            trajectory - was lost. From the first call until
            ``run_simulation`` returns both signals are ignored; the
            dispositions found here are put back at the very end. Does nothing
            off the main thread, where signal dispositions cannot be changed
            (and where no signal is ever delivered).

            Returns
            -------
            None
                ``held_signals`` records, for each signal now ignored, the
                disposition to restore.
            """
            self._stop_in_progress = True
            if not on_main_thread:
                return
            for signum in (signal.SIGINT, signal.SIGTERM):
                if signum in held_signals:
                    continue
                try:
                    current = signal.getsignal(signum)
                    if current is None:
                        # installed from outside Python: it could not be put back
                        continue
                    signal.signal(signum, signal.SIG_IGN)
                    held_signals[signum] = current
                except (ValueError, OSError):
                    pass

        marker_warning = restart.claim_run_directory()
        if marker_warning:
            IO_utils.status_message(marker_warning, 'warning')
            pimmslogger.log_warning(marker_warning)

        try:
            return self._run_simulation()
        except BaseException as error:
            active_exception = True
            _hold_stop_signals()
            # the ENERGY_CHECK abort writes the SAVE_AT_END buffer itself
            buffer_already_written = isinstance(error, SimulationEnergyException)
            self._report_unfinished_run(error)
            raise
        finally:
            # (again, for the paths that did not come through the except clause;
            # it does nothing the second time)
            _hold_stop_signals()
            try:
                writer = getattr(self, 'xtc_writer', None)
                if writer is not None:
                    # Clear first so an error raised by close can never lead to a
                    # second close attempt against the same native handle.
                    self.xtc_writer = None
                    try:
                        lattice_utils.close_xtc_writer(writer)
                    except Exception as cleanup_error:
                        if not active_exception:
                            raise
                        # Preserve the original simulation failure; cleanup errors
                        # are secondary but still useful in the log.
                        pimmslogger.log_warning(
                            'Unable to close XTC writer while handling another error: %s'
                            % cleanup_error)

                # SAVE_AT_END keeps every frame in memory until the end. A run
                # that does not reach the end would lose them all, so write what
                # there is (the ENERGY_CHECK abort has always done this itself).
                if (active_exception and not buffer_already_written
                        and getattr(self, 'SAVE_AT_END', False)
                        and getattr(self, 'master_traj_obj', None) is not None):
                    try:
                        lattice_utils.save_out_sim(self.master_traj_obj, self.current_xtc_filename)
                    except Exception as cleanup_error:
                        pimmslogger.log_warning(
                            'Unable to write the SAVE_AT_END trajectory buffer while handling '
                            'another error: %s' % cleanup_error)
            finally:
                try:
                    restart.release_run_directory()
                finally:
                    # put every signal disposition back as it was before the run:
                    # as the clean-up found it, or, where our own handler was in
                    # place (which is what the clean-up then found), as it was
                    # before that handler went in
                    restore = dict(held_signals)
                    restore.update(installed_handlers)
                    for signum, disposition in restore.items():
                        try:
                            signal.signal(signum, disposition)
                        except (ValueError, OSError, TypeError):
                            pass


    #-----------------------------------------------------------------
    #
    def _report_unfinished_run(self, error):
        """
        Say, in ``log.txt`` and on screen, how and where an unfinished run stopped.

        A run stopped by Ctrl-C used to end in a bare ``KeyboardInterrupt``
        traceback and one stopped by ``SIGTERM`` in nothing at all, and in both
        cases ``log.txt`` simply ended, indistinguishable from a run still
        going. This writes one line giving the last step that completed and the
        step of the checkpoint on disk, which is what is needed to decide how to
        resume. A run that fails before its first step is recorded as a failed
        start instead. It never raises: reporting must not replace the error
        being reported.

        Parameters
        ----------
        error : BaseException
            What ended the run: ``KeyboardInterrupt`` (``SIGINT``), the
            ``SystemExit`` raised by the ``SIGTERM`` handler, or any other
            exception.

        Returns
        -------
        None
            One line is appended to ``log.txt``; for an interrupt it is also
            printed.
        """
        try:
            done = getattr(self, '_last_completed_step', None)
            saved = getattr(self, '_last_checkpoint_step', None)
            if isinstance(error, KeyboardInterrupt):
                cause = 'interrupted by SIGINT (Ctrl-C)'
            elif getattr(self, '_stopped_by_signal', None) and isinstance(error, SystemExit):
                cause = 'stopped by %s' % self._stopped_by_signal
            else:
                cause = None

            # keyfile_used.kf is written by startup_analysis just before the
            # previous run's output is cleared; until then nothing of an earlier
            # run has been removed and the failure is a failed start
            if cause is None and done is None and not getattr(self, '_effective_keyfile_written', False):
                self._log_failed_start(error, 'before the first step, while its output files were being opened')
                return

            if saved is not None:
                checkpoint = ("The checkpoint on disk (%s) is from step %d; RESTART_CONTINUE resumes "
                              "the run from there, in a fresh directory."
                              % (CONFIG.RESTART_FILENAME, saved))
            elif os.path.exists(CONFIG.RESTART_FILENAME):
                checkpoint = ("This run had not yet written a checkpoint (none is written during "
                              "equilibration); the %s in this directory was left by an earlier run "
                              "or is the file this run started from." % CONFIG.RESTART_FILENAME)
            else:
                checkpoint = ("This run had not yet written a checkpoint (none is written during "
                              "equilibration), so there is nothing to resume from.")
            # "after step N" is only true once the analyses due on step N have
            # run too (the checkpoint among them, last); a stop that lands inside
            # them is reported as such
            if done is None:
                where = 'before the first step was completed'
            elif getattr(self, '_last_analysed_step', None) == done:
                where = 'after step %d of %d' % (done, self.n_steps)
            else:
                where = ('during the output of step %d of %d (its move, energy row and trajectory '
                         'frame were complete; the analyses due on that step may not be)'
                         % (done, self.n_steps))

            if cause is not None:
                msg = "Run %s %s. %s" % (cause, where, checkpoint)
                IO_utils.newline()
                IO_utils.status_message(msg, 'warning')
                sys.stdout.flush()
                pimmslogger.log_warning(msg)
            else:
                detail = ' '.join(str(error).split())
                pimmslogger.log_error("Run stopped by an error %s (%s%s). %s"
                                      % (where, type(error).__name__,
                                         (': ' + detail[:400]) if detail else '', checkpoint))
        except BaseException:
            # not even an interrupt may replace the error being reported: the
            # caller re-raises that error as soon as we return
            pass


    #-----------------------------------------------------------------
    #
    def _run_simulation(self):
        """
        Execute the full Monte Carlo simulation workflow.

        This method drives a complete production run from the current
        :class:`Simulation` state. It performs one-time startup tasks,
        iterates over ``self.n_steps`` Monte Carlo steps, dispatches move
        proposals, applies Metropolis acceptance/rejection logic, runs
        analysis and trajectory/energy I/O at configured frequencies, and
        performs final output and cleanup.

        High-level flow
        ---------------
        1. Record the global start time and print startup status.
        2. Compute initial system energy with the full Hamiltonian.
        3. Initialize quench and trajectory output files when enabled.
        4. Run startup analysis/file initialization routines.
        5. Enter the main simulation loop and repeat until ``n_steps`` is
             reached:

             - Skip move proposals if all chains are frozen, while still
               running scheduled quench, resize, output and analysis work.
             - If not in an auxiliary TSMMC chain:
                 - apply quench updates (if ``QUENCH_RUN``),
                 - handle box-resize equilibration logic (if ``resize_eq``),
             - If inside an auxiliary TSMMC chain, update/complete that chain
                 and only advance the global step counter once the TSMMC cycle is
                 complete.
             - Choose a random movable chain and sample a move type via
                 ``AcceptanceCalculator.move_selector``.
             - Execute the selected move implementation (single-chain,
                 cluster, TSMMC, system-shake, etc.).
             - For standard move families, compute energy deltas and apply
                 Boltzmann acceptance; on rejection, revert lattice state.
             - Update move statistics used for post-hoc diagnostics.
             - Write scheduled trajectory/energy output and run analysis
               against the resulting post-step state.

        6. After the loop, optionally flush an in-memory trajectory (when
             ``SAVE_AT_END`` is active), run final analysis, and always write a
             final restart snapshot.

        Side effects
        ------------
        - Mutates simulation state in place, including lattice coordinates,
            chain positions, energies, counters, and acceptance statistics.
        - Writes multiple output artifacts (trajectory, energy, quench,
            analysis files, restart file), depending on runtime options.
        - Emits progress and diagnostic messages to stdout/loggers.

        Notes
        -----
        - Invalid move-selection codes raise :class:`SimulationException`.
        - Energy consistency checks may raise
            :class:`SimulationEnergyException` via ``simulation_IO``.
        - The method is intentionally monolithic because step ordering is
            coupled to detailed-balance constraints and output semantics.

        Returns
        -------
        None

        """
        # get the time everything kicks off... The datetime is what we print; the
        # monotonic reading taken beside it is what elapsed times are measured
        # from, since the wall clock can be set (daylight saving, an NTP step)
        # while the run is going
        self.global_start_time = datetime.now()
        self.global_start_monotonic = time.monotonic()

        IO_utils.status_message("Simulation started at %s" % (str(self.global_start_time)),'startup')
        if CHECK_MEMORY:
            heap = hp.heap()
            print(heap)

        IO_utils.newline()

        # evaluate the initial energy of the system
        (old_energy, old_energy_local, old_energy_LR, old_energy_SLR, old_energy_angles) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)

        # NB: QUENCH.dat needs nothing here. It is in the output manifest, so
        # startup_analysis (below) deletes any stale copy, and write_quench_file
        # creates it if and when the run actually quenches.

        # setup the initial trajectory and pdb files
        if self.resize_eq is True and self.SAVE_EQ is False:
            # NB if we wanna (1) resize and (2) NOT save the equilibration stage then
            # we do NOT at this stage initialize an eq_start and eq_traj file because
            # it won't be written to
            pass
        else:
            # otherwise we initilize the eq_start and eq_traj files if self.resize_eq is True and SAVE_EQ is True, otherwise
            # initialize the START.pdb and traj.xtc files if self.resize_eq is False
            IO_utils.status_message("Building initial trajectory and pdb files...",'startup')
            if self.SAVE_AT_END:
                # SAVE_AT_END buffers the whole trajectory in memory (O(N)); just
                # write the topology PDB + an initial xtc (overwritten at the end).
                lattice_utils.start_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename=self.current_pdb_filename, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter, unwrap=self.trajectory_pbc_unwrap)
            else:
                # default incremental path: open a persistent XTC writer (O(1) per frame)
                self.xtc_writer = lattice_utils.open_xtc_writer(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename=self.current_pdb_filename, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter, unwrap=self.trajectory_pbc_unwrap)

        self.startup_analysis()

        IO_utils.status_message("Evaluating initial energy...",'startup')   
        IO_utils.newline()
        IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')        
        print("   ENERGY COMPARISON")   
        print("     STEP             : %i   " % 0)
        print("     GLOBAL           : %i" % old_energy)        
        print("     SHORT RANGE      : %i" % old_energy_local)
        print("     LONG RANGE       : %i" % old_energy_LR)
        print("     SUPER LONG RANGE : %i" % old_energy_SLR)
        print("     ANGLES           : %i" % old_energy_angles)

        IO_utils.newline()
        IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')
        print("   MEMORY USAGE")

        # Use .nbytes (the true size of the underlying data buffer) rather than
        # sys.getsizeof (which returns the Python wrapper size and is misleading
        # for numpy arrays - e.g. ~0 for views). The two lattice grids dominate
        # and scale as XDIM*YDIM*ZDIM * itemsize.
        _MB = 1048576.0
        # getattr(..., 'nbytes', 0) so a non-numpy grid (e.g. a test stub) never
        # crashes this cosmetic startup print; real grids are always numpy arrays.
        grids_mb = (getattr(self.LATTICE.grid, 'nbytes', 0) + getattr(self.LATTICE.type_grid, 'nbytes', 0)) / _MB

        # energy lookup tables: the SR/LR/SLR interaction matrices and the angle
        # lookup (the latter can be sizeable: n_residues x 3^6 in 3D)
        _tables = [getattr(self.Hamiltonian, _n, None) for _n in
                   ('residue_interaction_table', 'LR_residue_interaction_table',
                    'SLR_residue_interaction_table', 'angle_lookup')]
        tables_mb = sum(_t.nbytes for _t in _tables if hasattr(_t, 'nbytes')) / _MB

        print(f"     LATTICE GRIDS    : {grids_mb:8.1f} MB   (grid + type_grid)")
        print(f"     ENERGY TABLES    : {tables_mb:8.1f} MB   (interaction + angle lookup)")
        print(f"     DATA SUBTOTAL    : {grids_mb + tables_mb:8.1f} MB")

        # actual resident memory of the whole process (data + Python + numpy +
        # code) - the true footprint. ru_maxrss is bytes on macOS, KiB on Linux.
        try:
            import resource
            _rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            rss_mb = (_rss / _MB) if sys.platform == 'darwin' else (_rss / 1024.0)
            print(f"     PROCESS RESIDENT : {rss_mb:8.1f} MB   (actual; incl. Python + numpy)")
        except Exception:
            pass


        IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')
        IO_utils.newline()

        IO_utils.status_message('STARTING SIMULATION','major')
        IO_utils.status_message('  Start time: %s'%(self.global_start_time), 'vanilla')
        
        # flush means we flush all the premable text to STDOUT - useful for running
        # jobs on clusters 
        sys.stdout.flush()

        ##==================================================================##
        ##                                                                  ##
        ##                   MASTER LOOP BEGINS HERE!                       ##
        ##                                                                  ##
        ##==================================================================##
        i = getattr(self, 'continue_from_step', 0)
        pending_quench_output = {}
        if getattr(self, '_continue_rng', None) is not None:
            # Restore the generators the checkpoint recorded, as the last thing
            # before the loop: from here on every draw is the one the
            # uninterrupted run would have made. The temperature is restored
            # too (a quench ramp's value at the checkpoint), and the TSMMC
            # coordinator is rebuilt on it exactly as a quench rung does.
            random.setstate(self._continue_rng[0])
            np.random.set_state(self._continue_rng[1])
            if self._continue_temperature is not None and self._continue_temperature != self.ACC.temperature:
                self.ACC.update_temperature(self._continue_temperature)
                if self.TSMMC_USED:
                    self.TSMMC_coordinator = TSMMC(self.ACC.temperature, self.TSMMC_JUMP_TEMP,
                                                   self.TSMMC_INTERPOLATION_MODE, self.TSMMC_STEP_MULTIPLIER,
                                                   self.TSMMC_NUMBER_OF_POINTS, self.TSMMC_FIXED_OFFSET)
            IO_utils.status_message("Generator state and temperature (%.6g) restored from the restart file; "
                                    "continuing from step %d" % (self.ACC.temperature, i), 'info')
            pimmslogger.log_status("RESTART_CONTINUE: generator state restored, temperature %.6g, step %d"
                                   % (self.ACC.temperature, i))

        def run_post_step_output(step, energy_value):
            """
            Emit scheduled output only for completed master-chain steps.

            Writes any quench line that was deferred for this step, then runs the
            standard per-step IO and the scheduled analysis routines. Steps taken
            inside an auxiliary (TSMMC) chain are not part of the master Markov
            chain, so nothing is written for them - the output for the master step
            that owns the excursion is emitted once the excursion completes.

            Parameters
            ----------
            step : int
                The master-chain step number the output is being written for.

            energy_value : int or float
                The total system energy at the end of that step, written to the
                energy/quench files and passed to the analysis routines.

            Returns
            -------
            None
            """
            if self.auxillary_chain:
                return
            quench_temperature = pending_quench_output.pop(step, None)
            if quench_temperature is not None:
                analysis_IO.write_quench_file(step, quench_temperature, energy_value)
            self.simulation_IO(step, energy_value)
            self.run_all_analysis(step)

        # A system-wide TSMMC proposal selected on the final master step must be
        # allowed to finish its auxiliary sweep. Stopping solely because i reached
        # n_steps would serialize a half-completed Markov move and omit the final
        # scheduled state.
        while i < self.n_steps or self.auxillary_chain:
            i = i + 1
            
            # if we're not using an auxiliary chain (i.e. this is what happens
            # 99.9% of the time)
            if not self.auxillary_chain:

                # if we're doing a temperature quench..
                if self.QUENCH_RUN:
                    if self.quench_update(i, old_energy, write_output=False):
                        pending_quench_output[i] = self.ACC.temperature


                # if we're equilibrating in a different size box check what's goin' on there. If equilibration 
                # is done then update the energy as calculated via PBC 
                if self.resize_eq:

                    # the returned chain_selection_override is always empty: no
                    # forced-move mechanism exists, and a chain straddling a face at
                    # the resize step raises inside update_dimensions
                    (chain_selection_override, old_energy) = self.update_dimensions(i, old_energy)
                            
                # A fully frozen system still advances logical time. Quenches,
                # resizing and scheduled reporting above/below must continue even
                # though no proposal can be made. (With nothing frozen the count
                # is zero without building a set on every step.)
                if (len(set(self.frozen_chains)) if len(self.frozen_chains) > 0 else 0) >= self.LATTICE.get_number_of_chains():
                    run_post_step_output(i, old_energy)
                    continue

            # this is what happens if we're inside an auxiliary chain
            else:

                # decrement the global counter, as auxiliary chain moves don't count towards the 
                # global move count
                i=i-1                

                # this is where any/all updates happen to do with the TSMMC. The returned status tuple
                # tells us if the auxilary chain was complete and if the move was accepted or not
                tsmmc_move_status = self.auxillary_chain_update(old_energy)
                
                # IF the TSMMC move is finished!
                if tsmmc_move_status[0]:
                    # if move was accepted 
                    if tsmmc_move_status[1]:
                        success=True
                        pass
                        
                    # if move was rejected
                    else:
                        # NOTE reverting back to the pre-move lattice is done in the auxillary_chain_update
                        # function, so we just have to revert the energy back to the pre-move value
                        old_energy = self.TSMMC_coordinator.system_move_original_energy
                        success=False
                        
                    # finally we reset the temperature to the system temperature and 
                    # zero out temporary information held during the TSMMMC move
                    self.ACC = self.TSMMC_coordinator.system_move_finalize(self.ACC)
                    self.auxillary_chain = False
                    self.ACC.auxillary_chain = False

                    # NOTE this has to come after we turn the ACC auxiliary chain
                    # flag in in the AcceptanceCalculator object
                    self.ACC.update_move_logs(12, success)            
                    
                    # The TSMMC proposal was selected at master step i. Auxiliary
                    # updates hold i fixed; completion reports the resulting state
                    # at that same step. Quench/resize updates already ran before
                    # the proposal began and must not be applied twice.
                    run_post_step_output(i, old_energy)
                    
                    # finally continue to the next real main-chain move
                    continue
                    
            #*************************************************************
            ## Move time! 
            ##
            ## First we select a random chain to 
            ##
            
            # select a random chain to perturb            
            chain_to_move   = self.LATTICE.get_random_chain(frozen_chains=self.frozen_chains)
            chainID         = chain_to_move.chainID


            
            # get the currentposition of the chain we're going to move

            ## MOVE SELECTION -------------------------------------------------------------------
            
            # reset the move accepted flag
            move_accepted = False

            # select a move to make. The selection depends only on the MOVE_*
            # frequencies; a per-chain move that is undefined for the selected
            # chain (rotate/pivot of a single bead) is rejected as a null move
            # rather than being remapped to a different move type.
            selection = self.ACC.move_selector()
                        
            #
            #for chainID in self.LATTICE.chains:
            #    lattice_utils.check_chain_connectivity(chainID, self.LATTICE.chains[chainID].get_ordered_positions(), self.LATTICE.dimensions)

            # system shake
            if selection == 1:

                ## system_shake moves            
                (new_latticeObject, new_energy, total_proposed, total_accepted) = self.MOVER.system_shake(self.LATTICE,
                                                                                                          old_energy,
                                                                                                          self.ACC,
                                                                                                          self.Hamiltonian,
                                                                                                          self.CS_substeps,
                                                                                                          self.CS_mode,
                                                                                                          self.hardwall,
                                                                                                          self.frozen_chains,
                                                                                                          parallelize=self.parallelize,
                                                                                                          num_threads=self.parallel_threads)

                ## Finally record moves for post-hoc analysis of movesets
                self.ACC.megastep_update_move_logs(1, total_accepted, total_proposed)

                # update energy
                old_energy = new_energy                                
                
                # skip everything else, all hail the megamove! NOTE that we have induvidual accept/rejects inside the system_shake() so this is still performing
                # Metropolis Monte Carlo ON THE SAME MARKOV CHAIN [important] - the place where the move is accepted/rejected has just moved, but we're evaluating
                # with the same Hamiltonian at the same temperature.
                run_post_step_output(i, old_energy)
                continue
                                
            # translation
            elif selection == 2:                                
                (move_event, success) = self.MOVER.chain_translate(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)

            # rotation
            elif selection == 3:
                (move_event, success) = self.MOVER.chain_rotate(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                
            # chain pivot
            elif selection == 4:       
                (move_event, success) = self.MOVER.chain_pivot(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                
            # head pivoting
            elif selection == 5:                                
                (move_event, success) = self.MOVER.head_pivot(chain_to_move, self.LATTICE.grid, hardwall=self.hardwall)
                
            # chain slither
            elif selection == 6:

                # optimized whole-system slither megamove (every chain slithers
                # slither_substeps times, in random order). The 2D and 3D fast
                # kernels are selected inside system_slither.
                (new_latticeObject, new_energy, total_proposed, total_accepted) = self.MOVER.system_slither(self.LATTICE,
                                                                                                           old_energy,
                                                                                                           self.ACC,
                                                                                                           self.Hamiltonian,
                                                                                                           self.slither_substeps,
                                                                                                           self.hardwall,
                                                                                                           self.frozen_chains,
                                                                                                           parallelize=self.parallelize,
                                                                                                           num_threads=self.parallel_threads)

                self.ACC.megastep_update_move_logs(6, total_accepted, total_proposed)
                old_energy = new_energy

                # megamove: individual accept/rejects happen inside system_slither
                # on the SAME Markov chain, so skip the rest of the loop body.
                run_post_step_output(i, old_energy)
                continue


            # cluster translate
            elif selection == 7:
                (move_event, success) = self.MOVER.cluster_translate(chain_to_move, 
                                                                     self.LATTICE, 
                                                                     cluster_move_threshold=None,
                                                                     cluster_size_threshold=self.LATTICE.get_number_of_chains()-1,
                                                                     hardwall=self.hardwall,
                                                                     frozen_chains=self.frozen_chains)
                
            # cluster rotation
            elif selection == 8:
                (move_event, success) = self.MOVER.cluster_rotate(chain_to_move, 
                                                                  self.LATTICE, 
                                                                  cluster_move_threshold=None,
                                                                  cluster_size_threshold=self.LATTICE.get_number_of_chains()-1,
                                                                  hardwall=self.hardwall,
                                                                  frozen_chains=self.frozen_chains)

            ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ##
            ## Note that in TSMMC and ratchet_pivot moves we create an alternative Markov chain and accept-reject on that chain
            ## before finally accepting/rejecting the final conformation *back* into the true system chain, where all this accept
            ## and rejection occurs inside the MOVER's move function (hence the 'continue' at the end of these moves).
            ##
            ## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                
            # chain-based temperature sweep Metropolis Monte Carlo (TSMMC). 
            elif selection == 9:
                                                
                (new_latticeObject, new_energy, total_moves, success) = self.MOVER.Chain_based_TSMMC(chainID, self.LATTICE, old_energy, self.Hamiltonian, self.TSMMC_coordinator, self.hardwall)

                # Update the lattice object!
                self.LATTICE = new_latticeObject

                old_energy = new_energy                                

                # Rinally record moves for post-hoc analysis of movesets
                self.ACC.update_move_logs(selection, success)

                # Finally update the alternative Markov chain move count - used for
                # for performance 
                self.ACC.alt_Markov_chain_update_move_logs(total_moves)

                run_post_step_output(i, old_energy)
                continue 

            # multichain-based temperature sweep Metropolis Monte Carlo (TSMMC)
            elif selection == 10:
                (new_latticeObject, new_energy, total_moves, success) = self.MOVER.multichain_based_TSMMC(chainID, self.LATTICE, old_energy, self.Hamiltonian, self.TSMMC_coordinator, self.hardwall, self.frozen_chains)
                
                self.LATTICE = new_latticeObject
                
                old_energy = new_energy                                

                # Record moves for post-hoc analysis of movesets
                self.ACC.update_move_logs(selection, success)

                # Finally update the alternative Markov chain move count - used for
                # for performance 
                self.ACC.alt_Markov_chain_update_move_logs(total_moves)

                run_post_step_output(i, old_energy)
                continue 

            # pull
            elif selection == 11:

                # optimized whole-system pull (cooperative reptation) megamove -
                # every chain (length >= 3) is pulled pull_substeps times, in random
                # order. The 2D and 3D fast kernels are selected inside system_pull,
                # which maintains detailed balance internally.
                (new_latticeObject, new_energy, total_proposed, total_accepted) = self.MOVER.system_pull(self.LATTICE,
                                                                                                        old_energy,
                                                                                                        self.ACC,
                                                                                                        self.Hamiltonian,
                                                                                                        self.pull_substeps,
                                                                                                        self.hardwall,
                                                                                                        self.frozen_chains,
                                                                                                        parallelize=self.parallelize,
                                                                                                        num_threads=self.parallel_threads)

                self.ACC.megastep_update_move_logs(11, total_accepted, total_proposed)
                old_energy = new_energy

                # megamove: individual accept/rejects happen inside system_pull on
                # the SAME Markov chain, so skip the rest of the loop body.
                run_post_step_output(i, old_energy)
                continue

            # system-wide TSMMC
            elif selection == 12:

                if self.reduced_printing is False:
                    IO_utils.status_message("Performing System TSMMC...",'info')

                # create a backup and activate the auxiliary chain flags
                self.TSMMC_coordinator.start_system_TSMMC(self.LATTICE.lattice_backupcopy(), old_energy, self.ACC)
                self.auxillary_chain = True
                self.ACC.auxillary_chain = True
            
                # from this moment we are *in* the system TSMMC
                continue

            # jump and relax move (relax -> translate -> relax). Self-contained in
            # MoveObject: each of the three sub-steps preserves the Boltzmann
            # distribution, so the composite does too. Like the megamoves it does
            # its own accept/reject on the SAME Markov chain and skips the shared
            # acceptance block below.
            elif selection == 13:
                (new_latticeObject, new_energy, accepted) = self.MOVER.jump_and_relax_move(chain_to_move,
                                                                                           self.LATTICE,
                                                                                           old_energy,
                                                                                           self.ACC,
                                                                                           self.Hamiltonian,
                                                                                           self.CS_substeps,
                                                                                           self.CS_mode,
                                                                                           self.hardwall)
                old_energy = new_energy
                self.ACC.update_move_logs(13, accepted)
                run_post_step_output(i, old_energy)
                continue


            # VMMC (virtual-move Monte Carlo collective cluster move)
            elif selection == 14:

                # self-contained collective move: recruits a cluster of chains by
                # interaction-energy gradients (Whitelam & Geissler 2007) and
                # translates it rigidly, doing its own accept/reject on the SAME
                # Markov chain - so we skip the rest of the loop body, exactly like
                # the megamoves above.
                (new_latticeObject, new_energy, accepted, cluster_size) = self.MOVER.vmmc_move(chain_to_move,
                                                                                               self.LATTICE,
                                                                                               old_energy,
                                                                                               self.ACC,
                                                                                               self.Hamiltonian,
                                                                                               self.vmmc_max_displacement,
                                                                                               self.vmmc_max_cluster,
                                                                                               self.hardwall,
                                                                                               self.frozen_chains)
                old_energy = new_energy
                self.ACC.update_move_logs(14, accepted)

                # diagnostics: track accepted collective (multi-chain) moves
                if accepted and cluster_size > 1:
                    self.vmmc_accepted_multichain += 1
                    if cluster_size > self.vmmc_max_accepted_cluster:
                        self.vmmc_max_accepted_cluster = cluster_size
                run_post_step_output(i, old_energy)
                continue


            else:
                raise SimulationException('Invalid option passed... [%s]' % str(selection))
                

            ##
            ## MOVE SELECTION OVER, NOW WE DECIDED WHAT TO DO NEXT
            ##

            # If the hard-sphere energy allowed the move we then evaluate the change to the system energy
            if success:                

                # ..........................................................................................
                # SINGLE CHAIN MOVES! (1/2/3/4/5/6)
                if selection > 0 and selection < 7:
                                        
                    # determine the change in energy associated with this single chain move
                    local_dif = self.single_chain_move(move_event, chainID) 

                    # Check if the move is accepted based on the Metropolis-Hasting's criterion
                    if self.ACC.boltzmann_acceptance(old_energy, old_energy + local_dif):

                        # accepted - update old_energy
                        old_energy = old_energy + local_dif
                        
                        # update the flag!
                        move_accepted = True
                    else:
                        # rejected - re-configure the system back to its former glory!
                        self.single_chain_revert(move_event, chainID)
                    
                # ..........................................................................................
                # CLUSTER rotation/translation (7/8)
                elif selection > 6 and selection < 9:
                    
                    local_dif  = self.rigid_cluster_move(move_event.moved_positions, move_event.original_positions)
                    
                    # Check if the move is accepted based on the Metropolis-Hasting's criterion
                    if self.ACC.boltzmann_acceptance(old_energy, old_energy + local_dif):
                        # accepted!
                        old_energy = old_energy + local_dif

                        # update the flag!
                        move_accepted = True

                    else:
                        # rejected!
                        self.rigid_cluster_revert(move_event.moved_positions, move_event.original_positions)


            # in the case of success being False the move caused a hard-sphere clash and is rejected out
            # of hand
            else:
                pass

            ## Finally record move for post-hoc analysis of movesets
            self.ACC.update_move_logs(selection, move_accepted)

            run_post_step_output(i, old_energy)


        ###
        ### THE END IS NIGH!
        ### 

        # if we get here we have finished looping over the main simulation loop. Congrats?
            
        # save out the master traj if we are saving at end. Only do if True or we will overwrite the traj file.
        if self.SAVE_AT_END == True:
            if self.master_traj_obj is None:
                # no step qualified for a frame: the trajectory is frame 0 only,
                # exactly as the incremental writer would have left it. Do NOT
                # append the final state - it would be labelled as frame 1, i.e.
                # a step that never happened.
                self.master_traj_obj = lattice_utils.start_master_traj(self.current_pdb_filename)
            lattice_utils.save_out_sim(self.master_traj_obj, self.current_xtc_filename)
        else:
            # incremental path: flush and close the persistent XTC writer
            lattice_utils.close_xtc_writer(self.xtc_writer)
            self.xtc_writer = None

            
        global_end_time = datetime.now()
        IO_utils.newline()            
        IO_utils.status_message("Simulation complete", 'info')
        # record clean completion in log.txt too (the docs point users at the log
        # to see whether a run finished; stdout-only messages never got there)
        pimmslogger.log_status("Simulation complete (all %i steps finished)" % self.n_steps)

        # extract time and build an easy to read string! Use total_seconds rather
        # than relativedelta fields: relativedelta normalises >24 h into .days, so
        # printing only hours/minutes/seconds reported a 50-hour run as "2 hours".
        # Measured on the monotonic clock: the difference of two wall-clock
        # readings is wrong (even negative) if the clock was set during the run.
        total_secs = int(max(time.monotonic() - self.global_start_monotonic, 0.0))
        total_time_msg = "Simulation time:  %d hours, %d minutes, %d seconds" % (
            total_secs // 3600, (total_secs % 3600) // 60, total_secs % 60)

        IO_utils.status_message("Simulation finished at %s" % (str(global_end_time)), 'info')
        IO_utils.status_message(total_time_msg, 'info')
        IO_utils.newline()
        IO_utils.status_message("Performing final analysis output...", 'info')

        self.end_of_simulation_analysis()
    
        ### Always (regardless of interval) save a restart file corresponding to the final state of the simulation. 
        self.ANAFUNCT_save_restart(i)
    
        IO_utils.status_message(".... done!\n\nWe hope the results are all you hoped for!", 'info')
        IO_utils.newline()


    #-----------------------------------------------------------------
    #               
    def auxillary_chain_update(self, old_energy):
        """
        Advance (and possibly finalize) an in-progress system-wide TSMMC move.

        Performs all the busywork associated with the system-wide temperature
        switch Metropolis Monte Carlo (TSMMC) move. When the temperature sweep
        held by ``self.TSMMC_coordinator`` is not yet complete this simply checks
        the auxiliary chain in (updating the acceptance object only if the
        temperature changed). When the sweep is complete it applies the
        accept/reject decision: on rejection the lattice is restored from the
        coordinator's backup, on acceptance the (already-updated) lattice is kept.

        Whether this logic should live here, inside the TSMMC object, or inside the
        MOVER object is undecided; for now it lives here on the logic that it makes
        global changes to the simulation system and so belongs with the
        :class:`Simulation` object.

        Parameters
        ----------
        old_energy : int or float
            Current total system energy, used both for the acceptance test and for
            reporting the energy change of a completed move.

        Returns
        -------
        tuple of (bool, bool)
            A two-place tuple ``(move_complete, move_accepted)``. ``move_accepted``
            is always ``False`` while the move has not yet completed.
        """

        #print "On move %i" %(self.TSMMC_coordinator.system_move_count)
        
        # if check to see if the temperature-sweep has finished
        if self.TSMMC_coordinator.system_move_complete():

            # check if move was accepted
            if self.TSMMC_coordinator.accept_system_TSMMC(old_energy):

                if self.reduced_printing is False:
                    IO_utils.status_message("System TSMMC: ACCEPTED [dE = %5.5f]" % (old_energy - self.TSMMC_coordinator.system_move_original_energy),'info')
                
                # if we get here the move was accepted!!
                # DO NOT RESET THE LATTICE!                
                success = True     

            else:                

                if self.reduced_printing is False:
                    IO_utils.status_message("System TSMMC: REJECTED [dE = %5.5f]" % (old_energy - self.TSMMC_coordinator.system_move_original_energy),'info')
                # RESET THE LATTICE
                self.LATTICE.lattice_restorefrombackup(self.TSMMC_coordinator.system_move_original_info[0], self.TSMMC_coordinator.system_move_original_info[1], self.TSMMC_coordinator.system_move_original_info[2])
                success = False

            # reset the ACC back to its pre TSMMC move status                                    
            self.auxillary_chain = False        
            return (True, success)

        else:
            # note this only changes the ACC object if the temperature
            # has changed, but updates various local chain parameters
            # for book-keeping
            self.ACC = self.TSMMC_coordinator.check_in_system_TSMMC(self.ACC, old_energy)

            return (False, False)


    #-----------------------------------------------------------------
    #               
    def quench_update(self, i, old_energy, write_output=True):
        """
        Apply a temperature-quench update on quench steps.

        Helper that runs a quench update if a temperature quench run is being
        performed. On steps that are multiples of ``self.QUENCH_FREQ`` it either
        reports (and disables the quench) when the target temperature has been
        reached, or advances the temperature one ``QUENCH_STEPSIZE`` toward the
        target (negating the step for heating quenches). When TSMMC is in use the
        ``TSMMC_coordinator`` is rebuilt at the new temperature, and the quench
        event is optionally written to ``QUENCH.dat``. Updates all relevant
        simulation variables in place.

        Parameters
        ----------
        i : int
            Current simulation step number.

        old_energy : int or float
            Current total system energy, written to the quench output file when a
            temperature change occurs.

        write_output : bool, optional
            Write ``QUENCH.dat`` immediately. Default True. The master loop passes
            False so it can record the energy after the step's move; direct callers
            retain the historical immediate-write behaviour.

        Returns
        -------
        bool
            True if the temperature changed on this call, otherwise False.

        Raises
        ------
        SimulationException
            If ``self.QUENCH_FREQ`` is not a positive integer.
        """

        # if the current step is requires a temperature update
        if self.QUENCH_FREQ <= 0:
            raise SimulationException('QUENCH_FREQ must be a positive integer')

        if i % self.QUENCH_FREQ == 0:
                    
            if self.ACC.temperature == self.QUENCH_END:


                IO_utils.status_message(f'Reached target temperature of [{self.ACC.temperature}] - no change', 'info')
                pimmslogger.log_status(f'Target temperature reached on step {i} (Target={self.ACC.temperature})')

                # turn off the quench run flag as we're no longer performing a quench run
                self.QUENCH_RUN = False
            else:
                        
                # QUENCH_STEPSIZE already carries the correct sign from the keyfile
                # parser (negative for a heating run, positive for cooling), and the
                # update is always computed as `temperature - QUENCH_STEPSIZE`. Do NOT
                # re-negate here: doing so cancels the parser's sign flip and turns a
                # heating quench back into a cooling one.
                quench_step = self.QUENCH_STEPSIZE

                # update the temperature in an inteligent way
                self.ACC.update_temperature(nonequilibrium_utils.update_temperature_in_quench(quench_step, self.QUENCH_START, self.QUENCH_END, self.ACC.temperature, self.reduced_printing))
                        
                # update the TSMMC_coordinator temperature if TSMMC is being used (specifically, the TSMMC_coordinator object needs to know the main Markov Chain temperature so it
                # knows what temperature to return to
                if self.TSMMC_USED:
                    self.TSMMC_coordinator = TSMMC(self.ACC.temperature, self.TSMMC_JUMP_TEMP, self.TSMMC_INTERPOLATION_MODE, self.TSMMC_STEP_MULTIPLIER, self.TSMMC_NUMBER_OF_POINTS, self.TSMMC_FIXED_OFFSET)
                else:
                    self.TSMMC_coordinator = None
                        
                # finally write out to the quench file reporting on the quench event
                if write_output:
                    analysis_IO.write_quench_file(i, self.ACC.temperature, old_energy)
                return True

        return False


    #-----------------------------------------------------------------
    #                               
    def simulation_IO(self, i, old_energy):
        """
        Helper function to run simulation IO (O) for various different things. Additional
        output should be added her as an if statement comparing against the appropriate
        keyword frequency. 

        NOTE: All analysis is dealt with separately and shouldn't be added here - this is
        for non-analysis IO (i.e. status IO).

        This includes

        1. Printing status of the simulation

        2. Writing out trajectory information

        3. Writing out energy information

        4. Performing global energy comparison

        Parameters
        -----------------
        i : int
            Current step that the simulation is on

        old_energy : int or float
            Current (locally-tracked) total system energy. Used for status
            printing, written to the energy file, and compared against the
            from-scratch Hamiltonian recalculation during the periodic energy
            consistency check.

        Returns
        -------
        None

        Raises
        ------
        SimulationEnergyException
            If, on an energy-comparison step, the locally-tracked energy differs
            from the fully recalculated energy (a configuration snapshot is written
            to ``CONFIG_AT_ENERGY_FAIL.pdb``/``.xtc`` before raising).
        """

        ##
        # define a local function which we can then call in different
        # places. This just avoids us re-writing the same code in multiple
        # places
        def local_status():
            """
            Print a one-line step/progress/energy status message.

            Returns
            -------
            None
            """
            IO_utils.status_message("Step %i of %i [%2.3f %%] (Energy = %i)" %(i, self.n_steps, 100*(float(i)/float(self.n_steps)),old_energy), 'update')


        # flag that avoids this function re-printing the same information multiple times
        statusPrinted = False
        

        # first up we're going to do some performance analysis. This happens every 1/20th of the simulation AND 20
        # steps in so we get an initial estimate on how long this is gonna take quite quickly. 
        if i % self.five_percent == 0 or i == 20:
            # the monotonic start, so no elapsed time or rate ever touches the wall
            # clock (the datetime is the fallback for a caller that set only that)
            _start = getattr(self, 'global_start_monotonic', None)
            if _start is None:
                _start = self.global_start_time
            analysis_general.evaluate_performance(i, _start, self.n_steps, self.equilibration, self.ACC, start_step=getattr(self, 'continue_from_step', 0))
        
        # print status if we're at a printfreq interval of steps
        if i % self.printfreq == 0:
            
            if statusPrinted is False:
                local_status()            
                statusPrinted = True 

        # save coordinates
        if i % self.xtcfreq == 0:

            if statusPrinted is False:
                
                # if we're not doing reduced printing print!
                if self.reduced_printing is False:

                    # remark about saving coordinates only if we're saving coordinates
                    if self.SAVE_EQ == False:
                        if i > self.equilibration:
                            local_status()
                            statusPrinted = True
                            IO_utils.status_message("Saving coordinates...")
                    else:
                        local_status()
                        statusPrinted = True
                        IO_utils.status_message("Saving coordinates...")
                        

                # if we are doing reduced printing 
                else:

                    # is this 1/0th of the way through the simulation?
                    if i % self.ten_percent == 0:
                        local_status()
                        statusPrinted = True
                        # remark about saving coordinates only if we're saving coordinates
                        if self.SAVE_EQ==False:
                            if i > self.equilibration:
                                IO_utils.status_message("Saving coordinates [reduced printing mode]...")
                        else:
                            IO_utils.status_message("Saving coordinates [reduced printing mode]...")

            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=- # 
            # -=-=-=-=-=- SAVING THE traj.xtc FILE -=-=-=-=-=- #
            # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=- #
            
            # frames are buffered in the master trajectory object and written out once at
            # the end (no per-frame scratch PDB or XTC rewrite)
            if self.SAVE_AT_END==False:
                # default incremental path: O(1) append to the open XTC writer.
                # if we are saving eq, save regardless of eq step.
                if self.SAVE_EQ==True:
                    lattice_utils.write_xtc_frame(self.xtc_writer, self.LATTICE,
                                                  self.LATTICE.lattice_to_angstroms,
                                                  autocenter = self.autocenter, unwrap = self.trajectory_pbc_unwrap)
                else:
                    # check if we are passed the eq.
                    if i > self.equilibration:
                        lattice_utils.write_xtc_frame(self.xtc_writer, self.LATTICE,
                                                      self.LATTICE.lattice_to_angstroms,
                                                      autocenter = self.autocenter, unwrap = self.trajectory_pbc_unwrap)
            else:
                # if we are saving the xtc file at the end, we need to update the master traj object. 
                # however, we don't want to do this if we aren't saving at the end because it will slow things
                # down and take up memory. 
                if self.SAVE_EQ==True:
                    self.master_traj_obj = lattice_utils.update_master_traj(self.LATTICE, 
                                                                            self.LATTICE.lattice_to_angstroms,
                                                                            self.master_traj_obj,
                                                                            self.current_pdb_filename,
                                                                            autocenter = self.autocenter, unwrap = self.trajectory_pbc_unwrap)

                else:
                    if i > self.equilibration:
                        self.master_traj_obj = lattice_utils.update_master_traj(self.LATTICE, 
                                                                                self.LATTICE.lattice_to_angstroms,
                                                                                self.master_traj_obj,
                                                                                self.current_pdb_filename,
                                                                                autocenter = self.autocenter, unwrap = self.trajectory_pbc_unwrap)


        # save energy
        if i % self.enfreq == 0:
            analysis_IO.write_energy(i, old_energy)
                
        # check global energy
        if ('ENERGY_CHECK' not in getattr(self, 'disabled_frequencies', ()) and
                i % self.compare_energyfreq == 0):

            # if we haven't printed the status message yet, print it
            if statusPrinted is False:
                local_status()
                statusPrinted = True

            IO_utils.newline()
            IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')

            # recalculate the energy using the full Hamiltonian from scratch
            (recalculated_energy, new_energy_local, new_energy_long_range, new_SLR_energy, new_energy_angles) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)

            # calculate the difference between our locally-tracked energy and the fully recalcalculated energy (these should be the same)
            current_diff = recalculated_energy - old_energy

            # The recompute takes bead TYPES from type_grid - the same array the
            # kernels read - so a corrupted type grid would make tracked and
            # recomputed energies wrong identically and the comparison above could
            # not see it. Cross-check the grids against the chains' own sequences.
            grid_problems = self.LATTICE.check_grid_consistency()
            if grid_problems:
                for problem in grid_problems[:10]:
                    IO_utils.status_message(problem, 'error')
                    pimmslogger.log_error(problem)
                self._abort_energy_check(
                    "ERROR: occupancy/type grid is inconsistent with the chain objects "
                    "(%i problem(s) - see log)" % len(grid_problems))

            # print out the energy comparison and all current energy info
            print("   ENERGY COMPARISON")   
            print("     STEP             : %i   " % i)
            print("     GLOBAL           : %i" % recalculated_energy)
            print("     CURRENT          : %i" % old_energy)
            print("     DIFFERENCE       : %i" % current_diff)     
            print("     SHORT RANGE      : %i" % new_energy_local)
            print("     LONG RANGE       : %i" % new_energy_long_range)
            print("     SUPER LONG RANGE : %i" % new_SLR_energy)
            print("     ANGLES           : %i" % new_energy_angles)
            IO_utils.horizontal_line(hzlen=40, linechar='*', leader='  ')

            # uncomment for memory info...
            if CHECK_MEMORY:
                heap = hp.heap()
                print(heap)

            IO_utils.newline()

            # if the energy comparison is off, raise an exception and write out the current configuration
            if not current_diff == 0:
                self._abort_energy_check("ERROR: Something is wrong because energy comparisons were off...")

        # flush output
        sys.stdout.flush()

    #-----------------------------------------------------------------
    #
    def _abort_energy_check(self, message):
        """
        Save what can be saved and abort a run whose ``ENERGY_CHECK`` failed.

        Both failure modes (a tracked energy that disagrees with a from-scratch
        recompute, and occupancy/type grids that disagree with the chain objects)
        come through here, so both keep the trajectory up to the failure and
        leave a snapshot to inspect: the main XTC writer is closed, a
        ``SAVE_AT_END`` buffer is flushed to disk rather than discarded, and the
        configuration held by the chain objects is written to
        ``CONFIG_AT_ENERGY_FAIL.pdb`` / ``.xtc``.

        Parameters
        ----------
        message : str
            The message carried by the exception.

        Raises
        ------
        SimulationEnergyException
            Always, after the files are written.
        """
        # flush/close the main trajectory writer so traj.xtc is valid up to the
        # last frame before we abort
        lattice_utils.close_xtc_writer(self.xtc_writer)
        self.xtc_writer = None

        # under SAVE_AT_END the whole trajectory so far is buffered in memory;
        # write it out rather than discarding it with the abort
        if self.SAVE_AT_END and self.master_traj_obj is not None:
            lattice_utils.save_out_sim(self.master_traj_obj, self.current_xtc_filename)
            print('Writing out buffered trajectory to %s' % self.current_xtc_filename)

        lattice_utils.start_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename='CONFIG_AT_ENERGY_FAIL.pdb', xtc_filename='CONFIG_AT_ENERGY_FAIL.xtc')
        print('Writing out abort trajectory to CONFIG_AT_ENERGY_FAIL.pdb/xtc')
        # leave a record in log.txt whichever check failed
        pimmslogger.log_error(message + ' (configuration written to CONFIG_AT_ENERGY_FAIL.pdb/xtc)')
        raise SimulationEnergyException(message)
                

    #-----------------------------------------------------------------
    #               
    def single_chain_move(self, move_event, chainID):
        """
        Compute the energy change of a single-chain move.

        Implements the optimized local energy calculation for moves that perturb a
        single chain. The method evaluates the short-range, long-range, super
        long-range and angle contributions for the chain in both its original and
        moved positions (operating only on the local interaction envelope rather
        than the whole lattice) and commits the chain to its new position on the
        grid, type_grid and chain object. The returned value is the resulting
        change in total system energy, which the caller passes to the Metropolis
        acceptance test.

        Parameters
        ----------
        move_event : MoveEvent
            Object containing all the move details (original/moved positions, moved
            indices, move type, etc.).

        chainID : int
            ID of the single chain being moved. Each chain has a unique ID starting
            at 1 and increasing.

        Returns
        -------
        int
            The change in total system energy (``local_dif``) produced by the move,
            including short-range, long-range, super long-range and angle terms.
            Lattice energies are integer valued.
        """
        
        # the evaluation itself lives in moves.commit_single_chain_move, so that
        # the jump of the jump-and-relax move (a chain translation made inside
        # MoveObject) can use exactly the same local energy change
        return moves.commit_single_chain_move(self.LATTICE, self.Hamiltonian, move_event, chainID, self.hardwall)


    #-----------------------------------------------------------------
    #       
    def single_chain_revert(self, move_event, chainID):
        """
        Revert a rejected single-chain move.

        Restores the system back to its pre-move state after a single-chain move is
        rejected by the Metropolis criterion. The chain is removed from its moved
        positions and placed back at its original positions on the grid, the chain
        object's ordered positions are reset, and the type_grid is updated back.

        Parameters
        ----------
        move_event : MoveEvent
            Object containing the move details (moved/original positions and chain
            positions, and moved indices) used to undo the move.

        chainID : int
            ID of the single chain whose move is being reverted.

        Returns
        -------
        None
        """
        moved_positions            = move_event.moved_positions
        original_positions         = move_event.original_positions
        moved_chain_positions      = move_event.moved_chain_positions
        original_chain_positions   = move_event.original_chain_positions
        moved_indices              = move_event.moved_indices
                

        # revert the lattice to it's pre-move state 
        lattice_utils.delete_chain_by_position(moved_chain_positions, self.LATTICE.grid, chainID)
        lattice_utils.place_chain_by_position(original_chain_positions, self.LATTICE.grid, chainID, safe=True)
        
        self.LATTICE.chains[chainID].set_ordered_positions(original_chain_positions)
        
        # update the type_grid variable BACK
        self.LATTICE.update_type_grid(chainID, moved_positions, original_positions, moved_indices, safe=True)


    #-----------------------------------------------------------------
    #
    def rigid_cluster_move(self, new_chain_positions, old_chain_positions):
        """
        Function which implements optimized energy calculations for rigid body cluster moves.

        **Why the moved pairs are determined here**

        NOTE that unlike the single chain moves we actually determine the set of moved pairs inside this function. There's a reason for
        this! So, when making a rigid cluster move we first determine the positions of all the chains in the cluster we're moving. This provides us with
        a useful set of prior information because we KNOW that in that clusters' original position the ONLY short range interactions we care about are between lattice
        sites occupied by the cluster and lattice sites occupied by the solvent. We know this because any sites between a cluster-component and a NON solvent 
        site would be an intra-cluster pair, and given rigid cluster movements cannot be changing intra-cluster interactions the change in energy associated with 
        intracluster sites must be zero. However, the key computational cost here is evaluating how the long range interactions contribute to the cost
        of moving the cluster.

        Because we know all the relative interactions WITHIN the cluster must be held fixed (both short and long-range) the only interactions we care about
        are between the cluster interface. If we have NO long range interactions we can actually just perform the move without worrying about the energy because 
        - by definition - we cannot be moving a cluster into direct contact with another solute molecule so the cluster-system interface is purely solute-solvent
        before and after - i.e. no change in energy. If we do have long range interactions (as would be usual) we have to compute their influence. 

        1) Determine all the long-range interactions goin' on
        2) Just TRANSLATE/ROTATE those positions to get the interfacial pairs in the clusters' new position

        However, to do this we need the cluster back in its original position - hence why we have to move the lattice BACK to its original position before
        we determine the interfacial residues

        Parameters
        ----------
        new_chain_positions : dict
            Mapping of chainID to the list of new positions for that chain. The
            full set of new positions for the chains making up the cluster.

        old_chain_positions : dict
            Mapping of chainID to the list of original positions for that chain.
            The full set of old positions for the chains making up the cluster.

        Returns
        -------
        int
            The change in total system energy produced by the rigid cluster move.
            Returns ``0`` immediately for the energy-neutral case where the
            Hamiltonian has no long-range interactions (the move is still committed
            to the grids in that case).

        Raises
        ------
        Exception
            If ``new_chain_positions`` and ``old_chain_positions`` do not describe
            the same set of chainIDs.
        """
        
        dimensions = self.LATTICE.dimensions

        if not list(new_chain_positions.keys()) == list(old_chain_positions.keys()):
            raise Exception("I don't even care this should NOT HAPPEN")

        # ------------------------------------------------------------------------------------
        # shortcut incase we have no LR interactions then this move is automatically performed as it's
        # energy neutral, so no need to compute the <DELTA> energy as by definition it must be moving
        # a cluster from a fully solvated environment to a fully solvated environment
        if len(self.Hamiltonian.LR_residue_names) == 0:
                    
            # update the type grid (note we have to do this in two independent steps)
            for chainID in old_chain_positions:
                self.LATTICE.delete_chain_from_type_grid(chainID, old_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)

            # update the chain positions on the type grid and in the chains list
            for chainID in new_chain_positions:
                self.LATTICE.chains[chainID].set_ordered_positions(new_chain_positions[chainID])
                self.LATTICE.insert_chain_into_type_grid(chainID, new_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)

            # energy neutral move (an int: the tracked energy is an integer and a
            # float 0.0 here turned it into a float from the first accepted
            # cluster move on)
            return 0
        # ------------------------------------------------------------------------------------

        ## If there are LR interactions...
        # We want to evaluate the energy with the chain in both positions - right now
        # 1) self.LATTICE.grid has the chain in it's new positoin
        # 2) The chain object in self.LATTICE.chains has it's position in the OLD position
        # 3) The self.LATTICE.type_grid  has the chain in its old position too
        
        # So we revert self.LATTICE.grid back to the original position to get the energy 
        # associated with the cluster in its original position
        # (note the type grid was never changed so doesn't have to be 'reverted' back)
        chainIDs = list(new_chain_positions.keys())

        old_region_LR_pairs = {}
        new_region_LR_pairs = {}

        ## xoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxo
        ##
        ## STAGE 1 - original configuration
        ##        
        # for each chain (note this *has* to be a two step process of removal followed by
        # insertion - otherwise you may get clashes if you use a single for-loop
               
        # delete and re-insert (need to do in two steps)
        for chainID in chainIDs:
            lattice_utils.delete_chain_by_position(new_chain_positions[chainID], self.LATTICE.grid, chainID)

        for chainID in chainIDs:
            lattice_utils.place_chain_by_position(old_chain_positions[chainID], self.LATTICE.grid, chainID, safe=True)
            
        # now **all** chains have been moved back to their original possitione we can calculate the reduced 
        # LR pairs (where both residues in the pair participate in LR interactions) - note many (most) of 
        # these interactions will be *within* the cluster, but for simplicity of code we just let this happen - 
        # the computational cost of finding which pairs have one member outside the cluster is greater than
        # just doing all the pairs..
        old_region_LR_pairs = {}
        non_redundant_LR_pairs_old_full = []

        old_region_SLR_pairs = {}
        non_redundant_SLR_pairs_old_full = []

        for chainID in chainIDs:
            
            # get the positions of LR interaction residues in the chain            
            #LR_original_positions_tmp   = longrange_utils.get_LR_positions(old_chain_positions[chainID], range(0,len(old_chain_positions[chainID])), self.LATTICE.chains[chainID].LR_IDX)
            
            # get all the LR
            #old_region_LR_pairs[chainID] = longrange_utils.build_LR_envelope_pairs(LR_original_positions_tmp, self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions)            
            (old_region_LR_pairs[chainID], old_region_SLR_pairs[chainID])  = longrange_utils.build_LR_envelope_pairs(old_chain_positions[chainID], self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions, hardwall=self.hardwall)

            # get all the pairs of LR interactions between for the chainID
            non_redundant_LR_pairs_old_full.extend(old_region_LR_pairs[chainID])
            non_redundant_SLR_pairs_old_full.extend(old_region_SLR_pairs[chainID])

        # perform energy evaluation (pair lists deduped)
        ENERGY_old_lattice_old_region = self.Hamiltonian.evaluate_local_energy_LR(self.LATTICE, numpy_utils._dedupe_pair_rows(np.array(non_redundant_LR_pairs_old_full))) + self.Hamiltonian.evaluate_local_energy_SLR(self.LATTICE, numpy_utils._dedupe_pair_rows(np.array(non_redundant_SLR_pairs_old_full)))
        
        ## xoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxo
        ##
        ## STAGE 2 - cluster chains removed (so a reduced lattice) - recall we do now have to delete from main 
        ## and type lattice

        # first delete all the original chains
        for chainID in chainIDs:
            lattice_utils.delete_chain_by_position(old_chain_positions[chainID], self.LATTICE.grid, chainID)                
            self.LATTICE.delete_chain_from_type_grid(chainID, old_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)

        # now we set the LR energy associated with the empty lattice to 0 - because it must be by definition. long-range interactoins ONLY occur when
        # there is a PAIR of beads which participate in long-range interactions. In the empty lattice, at the sites where the chain(s) of interest were 
        # and will be the site *must* be empty, so there cannot be a pair of LR interacting beads.
        ENERGY_empty_lattice_old_region = 0
        ENERGY_empty_lattice_new_region = 0
        
        ## xoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxoxo
        ##
        ## STAGE 3 - add in all the cluster chains back into their new positions

        # now insert all the chains into their new positions in the lattice
        for chainID in chainIDs:

            # update the new chain position
            self.LATTICE.chains[chainID].set_ordered_positions(new_chain_positions[chainID])

            # add chains into main and type grids
            lattice_utils.place_chain_by_position(new_chain_positions[chainID], self.LATTICE.grid, chainID, safe=True)                
            self.LATTICE.insert_chain_into_type_grid(chainID, new_chain_positions[chainID], list(range(0,len(new_chain_positions[chainID]))), safe=True)

        new_region_LR_pairs = {}
        non_redundant_LR_pairs_new_full = []

        new_region_SLR_pairs = {}
        non_redundant_SLR_pairs_new_full = []

        for chainID in chainIDs:
            
            # get the positions of LR interaction residues in the chain            
            #LR_original_positions_tmp   = longrange_utils.get_LR_positions(new_chain_positions[chainID], range(0,len(new_chain_positions[chainID])), self.LATTICE.chains[chainID].LR_IDX)

            # get the dumb full LR pairs 
            #new_region_LR_pairs[chainID] = longrange_utils.build_LR_envelope_pairs(LR_original_positions_tmp, self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions)
            (new_region_LR_pairs[chainID], new_region_SLR_pairs[chainID]) = longrange_utils.build_LR_envelope_pairs(new_chain_positions[chainID], self.LATTICE.chains[chainID].get_LR_binary_array(), self.LATTICE.type_grid, dimensions, hardwall=self.hardwall)

            # get all the pairs of LR interactions between for the chainID
            non_redundant_LR_pairs_new_full.extend(new_region_LR_pairs[chainID])
            non_redundant_SLR_pairs_new_full.extend(new_region_SLR_pairs[chainID])

        # finally  perform energy evaluation (pair lists deduped)
        ENERGY_new_lattice_new_region = self.Hamiltonian.evaluate_local_energy_LR(self.LATTICE, numpy_utils._dedupe_pair_rows(np.array(non_redundant_LR_pairs_new_full))) + self.Hamiltonian.evaluate_local_energy_SLR(self.LATTICE, numpy_utils._dedupe_pair_rows(np.array(non_redundant_SLR_pairs_new_full)))

        ## STAGE 4
        # Having now calculated all the relevant LR interactions we sum them up and use them to evaluate the energy of the move
        local_dif     = (ENERGY_new_lattice_new_region + ENERGY_empty_lattice_old_region) - (ENERGY_old_lattice_old_region + ENERGY_empty_lattice_new_region)    
                
        return local_dif


    #-----------------------------------------------------------------
    #       
    def rigid_cluster_revert(self, new_chain_positions, old_chain_positions):
        """
        Revert a rejected rigid cluster move.

        Restores the system back to its pre-move state after a rigid cluster move
        is rejected. Every chain in the cluster is deleted from its new positions
        (on both the grid and the type_grid) and then re-inserted at its original
        positions, and each chain object's ordered positions are reset.

        Parameters
        ----------
        new_chain_positions : dict
            Mapping of chainID to the list of new (rejected) positions for that
            chain.

        old_chain_positions : dict
            Mapping of chainID to the list of original positions to restore for
            that chain.

        Returns
        -------
        None
        """
        # revert the lattice to it's pre-move state  (delete everything)        
        for chainID in new_chain_positions:            
            lattice_utils.delete_chain_by_position(new_chain_positions[chainID], self.LATTICE.grid, chainID)
            self.LATTICE.delete_chain_from_type_grid(chainID, new_chain_positions[chainID], list(range(0,len(new_chain_positions[chainID]))), safe=True)

        # now re-insert everything
        for chainID in old_chain_positions:
            lattice_utils.place_chain_by_position(old_chain_positions[chainID], self.LATTICE.grid, chainID, safe=True)
            self.LATTICE.insert_chain_into_type_grid(chainID, old_chain_positions[chainID], list(range(0,len(old_chain_positions[chainID]))), safe=True)
            self.LATTICE.chains[chainID].set_ordered_positions(old_chain_positions[chainID])


    #-----------------------------------------------------------------
    #          CHANGE ME
    def update_dimensions(self, step, old_energy):
        """
        Resize the lattice box at the end of resize-equilibration.

        Handles the box-resize-equilibration logic. On any step other than the end
        of equilibration this is a no-op that returns the energy unchanged. On the
        final equilibration step it checks whether any chain still straddles the
        periodic boundary:

        - A chain straddling a face is impossible there (the equilibration box is
          always hardwall and a periodic restart file is refused with
          ``RESIZED_EQUILIBRATION``), so finding one means the lattice state is
          corrupted and a :class:`latticeExceptions.SimulationException` is
          raised. Earlier versions extended the run by 100 steps and retried, a
          branch that could never fire.
        - Otherwise the lattice is rebuilt at the production dimensions via a
          :class:`restart.RestartObject` (optionally applying ``EQ_OFFSET``),
          output trajectory/PDB files are (re)initialized, the resize flag is
          cleared, hardwall is switched off if the production run uses PBC, and the
          energy is recomputed from scratch.

        Parameters
        ----------
        step : int
            Current simulation step number, compared against
            ``self.equilibration``.

        old_energy : int or float
            Current total system energy, returned unchanged except when the box is
            actually resized (in which case the energy is recomputed).

        Returns
        -------
        tuple
            ``(chain_selection_override, energy)``. ``chain_selection_override`` is
            always an empty list (kept so the caller's unpacking is unchanged; no
            forced-move mechanism exists). ``energy`` is the (possibly recomputed)
            total system energy.

        Raises
        ------
        latticeExceptions.SimulationException
            If any chain straddles a periodic face at the resize step.
        """
        
        # if this step is the end of equilibration do all the fun jazz, else we simply return
        # an empty list and the old energy
        # EQUILIBRATION is the LAST equilibration step everywhere else in PIMMS
        # (production analysis starts at EQUILIBRATION + 1), so the swap happens
        # before the first production move, i.e. once step > equilibration - the
        # step-EQUILIBRATION move and its trajectory frame stay in the equilibration
        # box. > (not ==) so that if the exact boundary step is consumed by something
        # that bypasses the pre-move block (historically the system-TSMMC completion
        # step) the resize still fires at the next opportunity; self.resize_eq is
        # cleared once the resize completes, so this cannot fire twice.
        if step > self.equilibration:
                
            # for each chain is in a non-periodic configuration
            # if no - keep selecting a cluster and move move move until yess
            # also print warning - probably means equilibration is too short!
            # repeat

            # assess if each chain on the lattice straddles the boundary or not 
            offending_chains=[]
            for chainID in self.LATTICE.chains:
                if self.LATTICE.chains[chainID].does_chain_stradle_pbc_boundary():
                    offending_chains.append(chainID)
                    
            # The equilibration box is hardwall whatever the keyfile says (set at
            # start-up) and a periodic restart file is refused with
            # RESIZED_EQUILIBRATION, so no move can leave a chain across a face
            # here. This used to extend N_STEPS and EQUILIBRATION by 100 and retry:
            # a branch that could never fire and, had it fired, would have changed
            # the run length silently. A chain across a face at this point is a
            # corrupted state, so stop rather than resize on top of it.
            if len(offending_chains) > 0:
                pimmslogger.log_error("%i chain(s) straddle a periodic face at the resize step of a hardwall equilibration. Offending chains: [%s]" % (len(offending_chains), offending_chains))
                raise SimulationException("RESIZED_EQUILIBRATION: %i chain(s) straddle a periodic face at the resize step, which cannot happen in a hardwall equilibration box - the lattice state is corrupted (chains %s)" % (len(offending_chains), offending_chains))

            # if we get here all the chains are valid, inasmuch as they are all within a non-periodic space, allowing us to change 
            # the lattice dimensions without fear of breaking everything! 
            pimmslogger.log_status("Resizing lattice dimensions from [%s] to [%s]" %(self.LATTICE.dimensions, self.production_dims))

            # create a new restart object, instantiate it with the current lattice, and then update the restart object's positions
            # to center 
            R = restart.RestartObject()
            R.build_from_lattice(self.LATTICE, self.production_hardwall)

            if self.EQ_OFFSET:
                R.update_lattice_dimensions(self.production_dims, manual_offset=self.EQ_OFFSET)
            else:
                R.update_lattice_dimensions(self.production_dims)
                                             
            # use this restart object to construct a new lattice
            # the [] is the 'empty' chains list which would normally be passed from the keyfile, but we can disregard here,
            # but is a required parameter (ugly, but it's OK...)
            new_lattice = Lattice(self.production_dims, [], self.Hamiltonian,
                                  self.LATTICE_TO_ANGSTROMS, restart_object=R,
                                  hardwall=self.production_hardwall)
            # once that's done then 

            # finally assign this new lattice to the simulation object 
            self.LATTICE = new_lattice

            # turn off the resize flag and update the output file names
            # see if we need to save the output when 'save at end' is set to True. . 
            if self.SAVE_AT_END == True:

                # if saving EQ == True
                if self.SAVE_EQ == True:

                    if self.master_traj_obj is None:
                        # no equilibration step qualified for a frame: save frame 0
                        # only. self.LATTICE is ALREADY the production lattice here,
                        # so appending it would write production coordinates into
                        # eq_traj.xtc under the equilibration unit cell.
                        self.master_traj_obj = lattice_utils.start_master_traj(self.current_pdb_filename)

                    # save the output

                    lattice_utils.save_out_sim(self.master_traj_obj, self.current_xtc_filename)

                    # reset master_traj_obj to None. 
                    self.master_traj_obj = None
            
            # set self.resize_eq to false, reset the namds of pdb and xtc files. 
            self.resize_eq = False
            self.current_pdb_filename = 'START.pdb'
            self.current_xtc_filename = 'traj.xtc'
            
            # initialize the xtc/pdb output files with these new names
            if self.SAVE_AT_END:
                lattice_utils.start_xtc_file(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename=self.current_pdb_filename, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter, unwrap=self.trajectory_pbc_unwrap)
            else:
                # close the equilibration writer (if any) and open a fresh persistent
                # writer for the production trajectory
                lattice_utils.close_xtc_writer(self.xtc_writer)
                self.xtc_writer = lattice_utils.open_xtc_writer(self.LATTICE, self.LATTICE.lattice_to_angstroms, pdb_filename=self.current_pdb_filename, xtc_filename=self.current_xtc_filename, autocenter=self.autocenter, unwrap=self.trajectory_pbc_unwrap)

            # clean up if possible!
            import gc                
            gc.collect()

            # If we want to switch to PBC based on the keyfile HARDWALL variable, then do so.
            # Regardless, we now recalculate the new energy and return this
            if self.production_hardwall is False:
                self.Hamiltonian.set_hardwall(False)
                self.hardwall = self.production_hardwall

            # the block decomposition depends on the box, so re-describe it - after
            # the boundary mode has switched to the production box's, so the report
            # does not label a periodic production box as hardwall
            if self.parallelize:
                self.report_parallelization(note='after resized equilibration (production box)')
            
            (energy, _, _, _, _) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)

            # return an empty list which sets the chain_selection_override to empty. Note that the calling function is aware of success
            # because self.resize_eq has been switched from True to False
            return([], energy)
        else:
            return ([], old_energy)

        
    

    ######################################################################################
    ##                                                                                  ##
    ##                             ANALYSIS ROUTINES                                    ##
    ##                                                                                  ##
    ######################################################################################
    #
    # The functions below are general setup for running sytem-wide analysis. Note that
    # actual analysis logic should NOT be included here, and should be implemented in 
    # either the Chains class or in the analysis_general.py.
    #
    # ANAFUNCT functions are the functions called by run_all_analysis(), which calls
    # analysis functions at different frequencies depending on how often the analysis
    # is to be performed as defined by the keyfile.
    #
    # These functions must take a single argument (the step number) - they don't have
    # to use it but it will always be passed.
    #
    # Many of the functions write data to disk, while others just update internal running 
    # totals.
    #
    #

    #-----------------------------------------------------------------
    #   
    def stale_output_files(self):
        """
        Build the list of output files a previous run may have left behind.

        This is the single canonical answer to "what could a PIMMS run write
        into this directory", and it is what start-up deletes. It is built from
        the manifest in CONFIG (``CONFIG.analysis_output_files``), every
        ``CHAIN_<type>_`` prefixed file matching that manifest, and the handful
        of files whose staleness depends on what this particular run will do.

        Two things are deliberately NOT in the returned list. The first is
        anything the current run has already written or opened by the time
        start-up analysis runs - the trajectory pair it has opened (for a
        resized-equilibration run that is the ``eq_`` pair, or nothing, so the
        production pair is stale and IS listed), ``parameters_used.prm``,
        ``log.txt``, and the angle/chain-to-chainID summaries when this run does
        write them. The second is ``restart.pimms``, which is overwritten at the
        first checkpoint rather than at start-up (so a run that dies early
        leaves the previous restart file usable; ``startup_analysis`` says so in
        ``log.txt``). The temporary files of interrupted checkpoint writes
        (``restart.pimms.tmp.<pid>``) ARE listed when they cannot belong to a
        write still in progress (the process is known to be gone from this
        machine, or the file is over an hour old: see
        ``restart.checkpoint_temporaries``), and so, when the process that owns
        one is no longer running, are the scratch files
        of the topology and ``SAVE_AT_END`` trajectory writers
        (``START.pdb.tmp.<pid>``, ``traj.xtc.tmp.<pid>`` and their ``eq_`` and
        ``CONFIG_AT_ENERGY_FAIL`` counterparts).

        Returns
        -------
        list
            Filenames to delete. Files that do not exist are included; the
            caller is expected to ignore absent files.
        """

        # every analysis output PIMMS knows how to write
        stale = CONFIG.analysis_output_files()

        # the per-chain-type variants are matched by glob rather than built from
        # this run's chain types, so that EVERY stale CHAIN_<T>_* file goes: a
        # re-run with fewer (or different) chain types otherwise left the old
        # run's files sitting beside the new outputs, silently mixing two runs'
        # data in any glob-based analysis
        for name in CONFIG.PER_CHAIN_TYPE_OUTPUT_NAMES:
            base = getattr(CONFIG, name)
            pattern = os.path.join(os.path.dirname(base),
                                   "CHAIN_*_" + os.path.basename(base))
            stale.extend(glob.glob(pattern))

        # the abort dump is written only by a run whose ENERGY_CHECK failed, so
        # a clean re-run must not leave the previous failure's snapshot behind
        stale.extend(['CONFIG_AT_ENERGY_FAIL.pdb', 'CONFIG_AT_ENERGY_FAIL.xtc'])

        # absolute_energies_of_angles.txt is (re)written by the Hamiltonian only
        # when angle penalties are on; an ANGLES_OFF re-run must not inherit the
        # previous run's summary of penalties it does not apply.
        if getattr(self, 'keyword_lookup', {}).get('ANGLES_OFF', False):
            stale.append(CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL)

        # likewise chain_to_chainid.txt, which is written at construction time
        # when WRITE_CHAIN_TO_CHAINID is on
        if not self.write_chain_to_chainid:
            stale.append(CONFIG.OUTPUT_CHAIN_TO_CHAINID)

        # a resized-equilibration run that saves its equilibration has ALREADY
        # opened its own eq_* files by the time this runs (the trajectory is
        # started before startup_analysis); with SAVE_EQ off it never opens them,
        # so a previous run's eq_* files are stale and must go
        if not (getattr(self, 'resize_eq', False) and getattr(self, 'SAVE_EQ', True)):
            stale.extend(['eq_traj.xtc', 'eq_START.pdb'])

        # conversely a resized-equilibration run opens the production pair only at
        # the resize, so at start-up any traj.xtc / START.pdb is a previous run's;
        # left in place, a run that died during equilibration sat beside another
        # run's production trajectory
        if getattr(self, 'resize_eq', False):
            stale.extend(['traj.xtc', 'START.pdb'])

        # temporary checkpoint files (restart.pimms.tmp.<pid>) left by runs that
        # were killed inside a checkpoint write; restart.pimms itself stays
        stale.extend(restart.stale_checkpoint_temporaries())

        # likewise the scratch files of the topology and SAVE_AT_END trajectory
        # writers (START.pdb.tmp.<pid>, traj.xtc.tmp.<pid>, eq_traj.xtc.tmp.<pid>,
        # ...), which a run killed inside one of those writes cannot remove
        stale.extend(lattice_utils.stale_writer_temporaries())

        return stale

    #-----------------------------------------------------------------
    #
    def startup_analysis(self):
        """
        Function for including all the analysis activity which should be run
        BEFORE the simulation starts.

        Output files are created lazily - every writer opens its file in append
        mode (or, for the end-of-run files, write mode) at the moment it has a
        row to put in it - so a file exists if and only if the run wrote to it.
        Nothing is created here. That is deliberate: which files a run produces
        cannot honestly be predicted at start-up (whether
        CLUSTER_RADIAL_DENSITY_PROFILE.dat gets a row depends on whether any
        cluster ever reaches the bead threshold), and before 1.0.8 start-up
        created about 25 files whether or not anything would ever be written to
        them.

        What DOES happen here is the reverse. Every output this run could write
        is deleted if a previous run in the same directory left a copy, since
        otherwise a re-run that does not write a given file silently inherits
        the previous run's data for it. See ``stale_output_files`` for the list.

        Two more things happen here because this is the first point at which
        the run is known to have started (the simulation was built and its
        trajectory opened). ``keyfile_used.kf`` is written, before anything is
        deleted, so that a start which fails earlier never leaves the effective
        keyfile of a run that did not happen beside an earlier run's data. And
        if the directory holds a ``restart.pimms``, ``log.txt`` is told that it
        stays there until this run writes its first checkpoint: no checkpoint is
        written during equilibration, so a run that stops before then leaves a
        restart file that does not describe it.

        Returns
        -------
        None
        """

        # Record the configuration this run actually uses. The keyfile the user
        # wrote is not always what ran: a restart override replaces HARDWALL or
        # DIMENSIONS, a resized equilibration forces hard walls for its first
        # phase, EXTRA_CHAIN merges into existing chain types, and a missing SEED
        # is generated at start-up. Post-processing tools (lemonade included)
        # used to have to re-derive all of that; now they can read it.
        keyword_lookup = getattr(self, 'keyword_lookup', None)
        if keyword_lookup is not None and not getattr(self, '_effective_keyfile_written', False):
            self.write_effective_keyfile(keyword_lookup)
            self._effective_keyfile_written = True

        if os.path.exists(CONFIG.RESTART_FILENAME):
            pimmslogger.log_status(
                "%s is present in this directory from before this run started (an earlier run's "
                "checkpoint, or the file this run started from). It stays until this run writes "
                "its first checkpoint, at the first RESTART_FREQ multiple after equilibration "
                "ends (step %s); if this run stops before then, that file does not describe "
                "this run." % (CONFIG.RESTART_FILENAME, getattr(self, 'equilibration', '?')))

        # a temporary checkpoint file that a running process may still own is
        # not removed (stale_output_files lists only the ones that are safe to)
        _kept = restart.checkpoint_temporaries()[1]
        if _kept:
            pimmslogger.log_status(
                "Left %s in place: a temporary checkpoint file less than an hour old whose "
                "process may still be running (it exists, cannot be checked from here, or may be "
                "on another machine sharing this directory). If no other run is using this "
                "directory it is the remnant of a killed run; a later start removes it once it "
                "is more than an hour old, or delete it by hand." % ", ".join(_kept))

        IO_utils.remove_files(self.stale_output_files())

                
    #-----------------------------------------------------------------
    #           
    def run_all_analysis(self, step):
        """
        Master analysis function - cycles over each type of analysis to
        assess if that analysis should be performed this step (or not), and
        launching the analysis if it should be done. 

        Updates various analysis state information (e.g. internal scaling
        distances associated with each Chain etc.) but does not change any
        lattice positions or anything like that

        Analysis is skipped entirely while the simulation is still in
        equilibration (``step <= self.equilibration``; ``EQUILIBRATION`` is the last equilibration step). Non-default-frequency
        analysis routines run on their own per-routine frequencies, while
        default-frequency routines run together every ``self.anafreq`` steps.

        The restart checkpoint (``ANAFUNCT_save_restart``) is always the LAST
        thing written for a step, whichever group it is in. A continuation
        starts after the checkpoint's step and never writes that step's rows, so
        everything the step is due to write has to be on disk before the
        checkpoint that says the step is done: written first, a run killed
        between the checkpoint and the remaining analyses left a step whose rows
        were in neither segment. This also means the generator states stored in
        the checkpoint are captured after the step's analyses rather than
        between them. The built-in analyses draw no random numbers, so the
        states are the same either way; a custom ``ANALYSIS_MODULE`` that draws
        from Python's or numpy's global generator is now resumed exactly when
        it runs on a checkpoint step (see docs/restart_files.rst).

        Parameters
        ----------
        step : int
            Current simulation step number, used to decide which analysis routines
            (if any) fire this step.

        Returns
        -------
        None
        """

        # the last master step whose move, energy row and trajectory frame are
        # complete (its analyses run below, and _last_analysed_step is set once
        # they have); reported if the run is interrupted (see
        # _report_unfinished_run)
        self._last_completed_step = step

        # do not perform analysis if we're still in equilibration
        # one boundary convention everywhere: step == EQUILIBRATION is the LAST
        # equilibration step - no analysis/restart snapshot (this gate), no
        # trajectory frame when SAVE_EQ is off (i > equilibration), and an 'E'
        # label in PERFORMANCE.dat (step <= equilibration). Analysis starts at
        # the first production step, equilibration + 1.
        if step <= self.equilibration:
            self._last_analysed_step = step
            return

        # for any analysis routines we've defined as occuring at a non-
        # default frequency (recall that the keys in self.[non_]default_freq_analysis
        # are the routines themselves: the NAME of a method of the Simulation class
        # that expects a single parameter (step), or a callable - see setup_analysis)
        # the checkpoint is held back until every other routine due this step has run
        checkpoint_functions = []

        def _routine_name(routine):
            """Name of an analysis routine, whichever way it is stored.

            Parameters
            ----------
            routine : str or callable
                A key of one of the two analysis dictionaries: the name of a
                method of this Simulation, or a callable.

            Returns
            -------
            str
                The method name, or the callable's ``__name__`` ('' if none).
            """
            if isinstance(routine, str):
                return routine
            return getattr(routine, '__name__', '')

        def _run_routine(routine):
            """Run one analysis routine for this step, on THIS Simulation.

            A routine stored by name is looked up on ``self`` here, at the call,
            so it is always this object's method - also for a deep copy of a
            Simulation, whose dictionaries hold the same names. A callable marked
            ``takes_simulation`` is handed ``self`` for the same reason; any
            other callable is called with the step alone.

            Parameters
            ----------
            routine : str or callable
                A key of one of the two analysis dictionaries.

            Returns
            -------
            None
            """
            if isinstance(routine, str):
                getattr(self, routine)(step)
            elif getattr(routine, 'takes_simulation', False):
                routine(step, self)
            else:
                routine(step)

        for analysis_function in self.non_default_freq_analysis:
            if step % self.non_default_freq_analysis[analysis_function] == 0:
                if _routine_name(analysis_function) == 'ANAFUNCT_save_restart':
                    checkpoint_functions.append(analysis_function)
                    continue
                _run_routine(analysis_function)

        # for all general analysis we haven't defined
        if self.default_freq_analysis and step % self.anafreq == 0:
            for analysis_function in self.default_freq_analysis:
                if _routine_name(analysis_function) == 'ANAFUNCT_save_restart':
                    checkpoint_functions.append(analysis_function)
                    continue
                _run_routine(analysis_function)

        for analysis_function in checkpoint_functions:
            _run_routine(analysis_function)

        # everything due on this step, checkpoint included, is now on disk
        self._last_analysed_step = step


    #-----------------------------------------------------------------
    #           
    def setup_analysis(self, keyword_lookup):
        """
        This function constructs two dictionaries. Each key-value pair in the dictionary is a function-frequency 
        pair, where the function is an analysis routine of the format FXC(step) and the frequency is the frequency
        with which that analysis is performed.

        The two lists correspond to analysis which occurs with the same frequency as the general analysis and then
        the analysis which occurs at a frequency *different* to the general analysis.

        This is a bit of work at the start, but allows us to run bespoke, custom-frequency analysis in a very
        simple way during the simulation.

        keyword_lookup provides all the info needed.

        Parameters
        ----------
        keyword_lookup : dict
            The controlled-vocabulary keyword dictionary. Provides each analysis
            keyword's frequency, the residue pairs for R2R analysis
            (``ANA_RESIDUE_PAIRS``), the default frequency (``ANALYSIS_FREQ``), and
            an optional side-loaded custom ``ANALYSIS_MODULE``.

        Returns
        -------
        tuple of (dict, dict)
            ``(non_default_freq_analysis, default_freq_analysis)``. Each maps an
            analysis routine to the integer frequency at which it should run.
            The first holds routines whose frequency differs from the default
            analysis frequency; the second holds those that match it. A routine
            is stored as the NAME of a ``Simulation`` method (the eight
            ``ANAFUNCT_*`` methods), or as a callable marked ``takes_simulation``
            that ``run_all_analysis`` calls as ``routine(step, simulation)`` (the
            residue-pair and custom-analysis routines). Neither keeps the
            Simulation alive - see the note in the body. The order is unchanged,
            and it is the order the routines run in.

        Raises
        ------
        SimulationException
            If the internal failsafe consistency checks on the analysis keyword
            tables fail (indicates a software bug).
        """
        import weakref

        non_default_freq_analysis = {}
        default_freq_analysis = {}

        # The two dictionaries built here are stored on the Simulation itself. A
        # bound method (or a closure over ``self``) used as a key therefore made a
        # reference cycle, Simulation -> dict -> method -> Simulation, and a
        # discarded Simulation - with its lattice grids - stayed alive until the
        # cyclic garbage collector happened to run: 51 dead lattices in a loop that
        # built 300 simulations. So nothing stored here keeps the Simulation alive.
        # The built-in routines are stored by method name and looked up on the
        # Simulation that is running them (``run_all_analysis``), and the two
        # routines that carry data of their own are plain functions that are
        # handed that Simulation as an argument. Looking the routine up at the call
        # is also what makes a ``copy.deepcopy`` of a Simulation analyse the copy
        # and not the original.

        # set the analysis names here
        all_ana_keywords = ['ANA_POL','ANA_INTSCAL', 'ANA_DISTMAP', 'ANA_ACCEPTANCE', 'ANA_CLUSTER', 'ANA_INTER_RESIDUE', 'ANA_END_TO_END', 'ANA_CUSTOM', 'RESTART_FREQ']

        # define the functions and initialze any closures needed
        analysis_keywords = {}
        analysis_keywords['ANA_POL']             = 'ANAFUNCT_polymeric_properties'
        analysis_keywords['ANA_INTSCAL']         = 'ANAFUNCT_internal_scaling'
        analysis_keywords['ANA_DISTMAP']         = 'ANAFUNCT_distance_map'
        analysis_keywords['ANA_ACCEPTANCE']      = 'ANAFUNCT_acceptance'
        analysis_keywords['ANA_CLUSTER']         = 'ANAFUNCT_cluster_analysis'
        analysis_keywords['ANA_INTER_RESIDUE']   = self.build_R2R_distance_distribution_analysis(keyword_lookup['ANA_RESIDUE_PAIRS'])
        analysis_keywords['ANA_END_TO_END']      = 'ANAFUNCT_end_to_end'
        analysis_keywords['RESTART_FREQ']        = 'ANAFUNCT_save_restart'

    
        # if a side-loading module was provided
        if keyword_lookup['ANALYSIS_MODULE']:

            # define a function that adds the LATTICE object to the call and then
            # calls the custom analysis function with the step and that lattice.
            # It is handed the Simulation by run_all_analysis (takes_simulation)
            # rather than closing over self. The weak reference is only the
            # fallback for a caller that invokes the routine directly with the
            # step alone, as could always be done: it names the Simulation that
            # built the routine, and it is never used by the run loop.
            built_by = weakref.ref(self)

            def fx(step, simulation=None):
                """
                Call the side-loaded custom analysis module with the live lattice.

                The user's ``analysis_function`` is validated at load time, but a
                runtime error can still occur once it sees real data. Any such
                exception is wrapped in an :class:`AnalysisRoutineException` that
                names the offending step and makes clear the fault is in the
                user-supplied analysis code, not in PIMMS itself, rather than
                surfacing as an opaque traceback deep inside the run loop.

                Parameters
                ----------
                step : int
                    Current simulation step number.

                simulation : Simulation, optional
                    The Simulation whose lattice is analysed. ``run_all_analysis``
                    always passes the Simulation it is running on. Left out, the
                    Simulation that built this routine is used.

                Returns
                -------
                object
                    Whatever the custom analysis module returns.

                Raises
                ------
                AnalysisRoutineException
                    If the custom ``analysis_function`` raises at runtime, or if
                    the routine is called without a Simulation after the one that
                    built it has been discarded.
                """
                if simulation is None:
                    simulation = built_by()
                if simulation is None:
                    raise AnalysisRoutineException(
                        f"The custom analysis routine was called at step {step} with no "
                        "Simulation, and the Simulation that built it no longer exists. Call "
                        "it as routine(step, simulation), or through "
                        "Simulation.run_all_analysis().")
                custom_analysis = keyword_lookup['ANALYSIS_MODULE']
                try:
                    return custom_analysis(step, simulation.LATTICE)
                except Exception as e:
                    raise AnalysisRoutineException(
                        f"The custom analysis function (from ANALYSIS_MODULE) raised "
                        f"{type(e).__name__} at step {step}: {e}. This is an error in "
                        "your custom analysis code, not in PIMMS."
                    ) from e

            fx.takes_simulation = True
            analysis_keywords['ANA_CUSTOM']          = fx

        else:
            analysis_keywords['ANA_CUSTOM']          = 'ANAFUNCT_custom_stubb'
            
        
        # ------------------------------------------------------->>

        # quick check to ensure all our ducks are in a row...
        if not len(all_ana_keywords) == len(analysis_keywords):
            raise SimulationException('Bug in the the analysis setup routines. This was triggered by a failsafe check and indicates a software bug')

        for AKW in all_ana_keywords:
            if AKW not in analysis_keywords:
                raise SimulationException('Bug in the the analysis setup routines. This was triggered by a failsafe check and indicates a software bug')

        # ------------------------------------------------------->>

    
        # get the default analysis frequency
        anafreq = keyword_lookup['ANALYSIS_FREQ']

        # Having set up all that we now cycle through the analysis types as defined by the all_ana_keywords
        # list. This means that the default_freq_analysis and non_default_freq_analysis dictionaries have 
        # a key-value pairing where the _key_ is the actual function signature and the _value_ is the frequency
        # with which that analysis is done 

        disabled_frequencies = set(keyword_lookup.get(
            '__DISABLED_FREQUENCIES', getattr(self, 'disabled_frequencies', ())))

        for AKW in all_ana_keywords:
            if AKW in disabled_frequencies:
                continue
            if keyword_lookup[AKW] == anafreq:
                default_freq_analysis[analysis_keywords[AKW]] = anafreq
            else:
                non_default_freq_analysis[analysis_keywords[AKW]] = keyword_lookup[AKW]

        return (non_default_freq_analysis, default_freq_analysis)
                
                            

        
    #-----------------------------------------------------------------
    #       
    def end_of_simulation_analysis(self):
        """
        Final analysis routines run at the end of the simulation. For all analysis
        where a final average value makes sense this is going to be where the code
        to calculate and save that output is written.

        Currently computes and writes the chain-averaged internal scaling, scaling
        exponents (nu, R0) and distance maps, handling both single-chain-type and
        multicomponent (per-chain-type) systems.

        Returns
        -------
        None
        """

        disabled = set(getattr(self, 'disabled_frequencies', ()))
        do_internal_scaling = 'ANA_INTSCAL' not in disabled
        do_distance_map = 'ANA_DISTMAP' not in disabled

        # A disabled analysis must neither execute nor produce a plausible-
        # looking zero/(-1) final file.  Previously these finalizers ran
        # unconditionally even though their per-step accumulators were disabled.
        if not do_internal_scaling and not do_distance_map:
            return

        # Group once, then calculate only the requested final products.  This
        # avoids distance-map allocation/fitting work for disabled analyses and
        # removes duplicate single-/multi-component implementations.
        chains_by_type = {chain_type: [] for chain_type in self.LATTICE.chainTypeList}
        for chain_object in self.LATTICE.chains.values():
            chains_by_type[chain_object.chainType].append(chain_object)

        single_type = len(self.LATTICE.chainTypeList) == 1
        for chain_type in self.LATTICE.chainTypeList:
            prefix = False if single_type else 'CHAIN_%i_' % chain_type
            chain_objects = chains_by_type[chain_type]

            # An analysis none of whose multiples fell in the production window was
            # never sampled. Its accumulators are all zero, and writing them out produced
            # a plausible-looking profile of exact zeros - so write nothing at all.
            #
            # SCALING_INFORMATION.dat used to be carved out of this and written with
            # its -1 -1 rows, on the grounds that -1 is a sentinel rather than
            # plausible data. That carve-out did not do what it claimed: -1 is also
            # what a chain shorter than the 26-bead fitting floor writes in a
            # perfectly well sampled run, so an all--1 file does not distinguish
            # "never sampled" from "sampled, nothing long enough to fit" - the two
            # cases are byte-identical for any system whose chains are short. The
            # never-sampled case is reported by the warning below, in the run's log,
            # and by the absence of INTSCAL.dat beside it.
            # Per chain type: the flags must not be cleared for the loop as a whole,
            # or every type after the first silently loses its warning (each type is
            # checked and warned on its own).
            write_internal_scaling = do_internal_scaling
            write_distance_map = do_distance_map
            if chain_objects and chain_objects[0].internal_scaling.count == 0:
                if do_internal_scaling:
                    msg = ("Internal scaling was never sampled (no production step is a "
                           "multiple of ANA_INTSCAL) - not writing %sINTSCAL.dat or "
                           "%sSCALING_INFORMATION.dat"
                           % (prefix or '', prefix or ''))
                    IO_utils.status_message(msg, 'warning')
                    pimmslogger.log_warning(msg)
                write_internal_scaling = False
            if chain_objects and chain_objects[0].distance_map.count == 0:
                if do_distance_map:
                    msg = ("Distance maps were never sampled (no production step is a "
                           "multiple of ANA_DISTMAP) - not writing %sDISTANCE_MAP.dat" % (prefix or ''))
                    IO_utils.status_message(msg, 'warning')
                    pimmslogger.log_warning(msg)
                write_distance_map = False

            if write_internal_scaling:
                all_is = [
                    chain.analysis_get_cumulative_internal_scaling()
                    for chain in chain_objects
                ]
                all_is_squared = [
                    chain.analysis_get_internal_scaling_squared()
                    for chain in chain_objects
                ]
                scaling_info = [
                    chain.analysis_fit_scaling_exponent()
                    for chain in chain_objects
                ]

                analysis_IO.write_internal_scaling(
                    np.asarray(all_is).mean(axis=0),
                    np.asarray(all_is_squared).mean(axis=0),
                    prefix=prefix)
                analysis_IO.write_scaling_information(
                    [values[0] for values in scaling_info],
                    [values[1] for values in scaling_info],
                    prefix=prefix)

            if write_distance_map:
                analysis_IO.write_distance_map(
                    _mean_over_chains(chain.analysis_get_cumulative_distance_map()
                                      for chain in chain_objects),
                    prefix=prefix)


    #-----------------------------------------------------------------
    #       
    def build_R2R_distance_distribution_analysis(self, R2R_info):
        """
        This function returns a function which is initialized by the variables pass
        in by the R2R_info - in essence generating a closure.

        Basically, if you're not familiar with functional programming, this creates
        a new function where the $R2R_info variable inside the ANAFUNCT_R2R_distance
        function is set by the build_R2R_distance_distribution_analysis function.

        This function (ANAFUNCTION_R2R_distance) is then returned, and next time
        its called the R2R_info variable IN THE FUNCTION BEING CALLED is already
        initialized.

        Parameters
        ----------
        R2R_info : list of tuple of int
            Residue-index pairs ``(i, j)``, taken from the ANA_RESIDUE_PAIRS
            keyword, for which the residue-residue distance should be computed and
            written on every analysis step. May be empty, in which case the
            returned closure does nothing.

        Returns
        -------
        callable
            A function ``ANAFUNCT_R2R_distance(step, simulation=None)`` that
            computes the requested residue-residue distances across all chains of
            ``simulation`` for the given step and writes them to disk. It carries
            the pairs and nothing else: the function is stored on the Simulation,
            so closing over ``self`` would form a reference cycle that keeps a
            discarded Simulation (and its lattice) alive until the cyclic garbage
            collector runs. It is marked ``takes_simulation``, and
            ``run_all_analysis`` passes the Simulation it is running on - which is
            what makes a deep copy of a Simulation measure the copy. Called with
            the step alone, it falls back on the Simulation that built it, through
            a weak reference the run loop never uses.
        """
        import weakref
        built_by = weakref.ref(self)

        def ANAFUNCT_R2R_distance(step, simulation=None):
            """
            Compute and write residue-residue distances for the captured pairs.

            Parameters
            ----------
            step : int
                Current simulation step number, written alongside the data.

            simulation : Simulation, optional
                The Simulation whose chains are measured. ``run_all_analysis``
                always passes the Simulation it is running on. Left out, the
                Simulation that built this routine is used.

            Returns
            -------
            None

            Raises
            ------
            AnalysisRoutineException
                If the routine is called without a Simulation after the one that
                built it has been discarded.
            """

            # just skip if no pairs defined...
            if len(R2R_info) == 0:
                # nothing to measure - return instead of falling through to the writer,
                # which would pointlessly reopen RES_TO_RES_DIST.dat every analysis step
                return

            if simulation is None:
                simulation = built_by()
            if simulation is None:
                raise AnalysisRoutineException(
                    "The residue-pair analysis was called at step %s with no Simulation, and the "
                    "Simulation that built it no longer exists. Call it as "
                    "routine(step, simulation), or through Simulation.run_all_analysis()." % step)

            all_data = []
            # whole-chain positions once per chain, not once per (chain, pair)
            chains = simulation.LATTICE.chains
            chain_ids = sorted(chains.keys())
            whole = {chainID: chains[chainID].get_analysis_positions() for chainID in chain_ids}
            for pair in R2R_info:

                pair_data = []
                for chainID in chain_ids:
                    pair_data.append(chains[chainID].analysis_get_residue_residue_distance(pair[0], pair[1], positions=whole[chainID]))
                
                all_data.append(pair_data)
            
            # finally write the analysis to file
            analysis_IO.write_residue_residue_distance(step, R2R_info, all_data)


        # return the function for use; run_all_analysis hands it the Simulation
        ANAFUNCT_R2R_distance.takes_simulation = True
        return ANAFUNCT_R2R_distance
                
                                
    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_internal_scaling(self, step):
        """
        Run internal scaling analysis.

        Updates the running internal-scaling counters (both normal and squared)
        associated with each chain on the lattice. No data is written to disk on
        each call; the accumulated averages are written at the end of the
        simulation.

        Parameters
        ----------
        step : int
            Current simulation step number (accepted for interface uniformity; not
            used directly).

        Returns
        -------
        None
        """

        for chainID in self.LATTICE.chains:
            
            # note this updates both normal and squared internal scaling info
            self.LATTICE.chains[chainID].analysis_update_internal_scaling()
            


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_distance_map(self, step):
        """
        Run distance map analysis.

        Updates the running distance-map counters associated with each chain on
        the lattice. No data is written to disk on each call; the accumulated map
        is written at the end of the simulation.

        The map is a ``seqlen x seqlen`` float64 matrix PER CHAIN, and the
        output file has ``seqlen ** 2`` entries, so the cost grows with the
        square of the chain length. On the first call, before any accumulator
        is allocated, the cost is estimated for every chain type
        (``analysis_structures.distance_map_cost``) and a warning naming
        ``ANA_DISTMAP`` is printed and logged for each type whose estimated
        peak memory exceeds 1 GiB or whose output file exceeds 256 MB.

        Parameters
        ----------
        step : int
            Current simulation step number (accepted for interface uniformity; not
            used directly).

        Returns
        -------
        None
        """

        if not getattr(self, '_distance_map_cost_reported', False):
            self._distance_map_cost_reported = True

            # imported here because this is the only place the simulation needs it
            from . import analysis_structures

            chains_by_type = {}
            for chain_object in self.LATTICE.chains.values():
                chains_by_type.setdefault(chain_object.chainType, []).append(chain_object)

            for chain_type in sorted(chains_by_type):
                n_chains = len(chains_by_type[chain_type])
                seqlen = max(len(chain_object) for chain_object in chains_by_type[chain_type])
                (peak_memory, file_bytes) = analysis_structures.distance_map_cost(seqlen, n_chains)
                if (peak_memory > analysis_structures.DISTANCE_MAP_WARN_MEMORY_BYTES or
                        file_bytes > analysis_structures.DISTANCE_MAP_WARN_FILE_BYTES):
                    msg = ("ANA_DISTMAP: the distance map of chain type %i (%i chain(s) of %i beads) "
                           "keeps one %i x %i matrix per chain for the whole run. Estimated peak "
                           "memory %.1f GB (those matrices plus about two more of working memory) "
                           "and an output file "
                           "(DISTANCE_MAP.dat, or CHAIN_%i_DISTANCE_MAP.dat in a multi-component "
                           "system) of about %.0f MB. Both grow with the square of the chain length; "
                           "set ANA_DISTMAP : 0 to switch the distance map off if you do not need it"
                           % (chain_type, n_chains, seqlen, seqlen, seqlen, peak_memory / 1e9,
                              chain_type, file_bytes / 1e6))
                    IO_utils.status_message(msg, 'warning')
                    pimmslogger.log_warning(msg)

        for chainID in self.LATTICE.chains:
            self.LATTICE.chains[chainID].analysis_update_distance_map()


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_cluster_analysis(self, step):
        """
        Run cluster analysis.

        Computes both the contact (short-range) and long-range cluster
        distributions for the current configuration, corrects cluster positions
        into a single periodic image, and derives polymeric properties, gross
        size/shape properties (volume, surface area, density) and radial density
        profiles for the size-thresholded clusters. All results are written to
        disk on each call, making this routine I/O heavy.

        Parameters
        ----------
        step : int
            Current simulation step number, written alongside the cluster data.

        Returns
        -------
        None
        """

        # get clusters list - note this is really computationally expensive
        # so we try and only do this once and then perform any/all cluster analysis
        # subsequent to this!
        (clusters) = lattice_analysis_utils.get_cluster_distribution(
            self.LATTICE.grid, self.LATTICE.chains, hardwall=self.hardwall)
        # pass the interaction tables so LR clusters are joined only through
        # pairs with nonzero LR/SLR energy (the documented definition)
        (LR_clusters) = lattice_analysis_utils.get_LR_cluster_distribution(
            self.LATTICE, hardwall=self.hardwall,
            LR_table=self.Hamiltonian.LR_residue_interaction_table,
            SLR_table=self.Hamiltonian.SLR_residue_interaction_table)

        big_cluster_idx = []        
        for c_idx in range(0,len(clusters)):            
            if len(clusters[c_idx]) > self.analysis_settings.cluster_threshold:
                big_cluster_idx.append(c_idx)
            else:
                break
        

        big_LR_cluster_idx = []
        for c_idx in range(0,len(LR_clusters)):            
            if len(LR_clusters[c_idx]) > self.analysis_settings.cluster_threshold:
                big_LR_cluster_idx.append(c_idx)
            else:
                break

        # for each cluster extract the set of induvidual positions to get a list of positions
        cluster_positions    = lattice_analysis_utils.extract_positions_from_clusters(clusters, self.LATTICE.chains)
        LR_cluster_positions = lattice_analysis_utils.extract_positions_from_clusters(LR_clusters, self.LATTICE.chains)

        # for each cluster correct the cluster's positions such that each cluster lies in a single periodic image as best can be achieved. Note that when
        # we don't correct for this the cluster analysis ends up being confusing...
        if self.hardwall:
            # Hardwall coordinates already occupy one Cartesian image; applying a
            # periodic snakesearch can move beads by a full box and corrupt Rg,
            # hull and radial-density measurements.
            corrected_cluster_positions = cluster_positions
            corrected_LR_cluster_positions = LR_cluster_positions
        else:
            # the gather's own percolation test is switched off: every cluster is
            # tested exactly once, below, and the same documented warning is
            # raised from there. Leaving it on ran the test twice per cluster.
            corrected_cluster_positions = lattice_analysis_utils.correct_cluster_positions_to_single_image(
                cluster_positions, self.LATTICE.dimensions, warn_if_percolating=False)
            # the LR gather walks the relation that DEFINES an LR cluster (contact,
            # or Chebyshev 2/3 with a nonzero table entry), so it needs the types
            # and the tables; walking a plain distance-3 rule instead tore
            # elongated clusters against the periodic face
            corrected_LR_cluster_positions = lattice_analysis_utils.correct_LR_cluster_positions_to_single_image(
                LR_cluster_positions, self.LATTICE.dimensions,
                type_grid=self.LATTICE.type_grid,
                LR_table=self.Hamiltonian.LR_residue_interaction_table,
                SLR_table=self.Hamiltonian.SLR_residue_interaction_table,
                warn_if_percolating=False)

        ## subselect size-thresholded clusters for polymer/gross property/radial distribution analysis. The clusters are sorted by size, so we know that
        # once we find one cluster below the the threshold we've found all the big clusters, hence the 'break' statements
        big_clusters = [corrected_cluster_positions[i] for i in big_cluster_idx]
        big_clusters_LR = [corrected_LR_cluster_positions[i] for i in big_LR_cluster_idx]

        # A cluster connected to its own periodic image is an unbounded object, so
        # its shape, hull and radial quantities are simply undefined - the gather
        # still returns coordinates, but they are an arbitrary window cut out of an
        # infinite object whose Rg is set by the box, not the condensate (rooting
        # the walk at different beads of one fixed configuration moves Rg by ~30%).
        # These used to be written as plausible numbers with no marker at all, the
        # only signal being a UserWarning the interpreter shows once per process.
        # Hardwall boxes cannot percolate and take no gather at all.
        #
        # The test is run once per cluster, on every cluster (the gather used to
        # test - and warn about - every cluster, and then the size-thresholded
        # ones were tested a second time here), with the connectivity that
        # cluster was built and gathered with: contact for the short-range
        # clusters, the typed LR/SLR relation for the long-range ones.
        if self.hardwall:
            cluster_percolating = [False] * len(big_clusters)
            LR_cluster_percolating = [False] * len(big_clusters_LR)
        else:
            all_percolating = lattice_analysis_utils.flag_percolating_clusters(
                corrected_cluster_positions, self.LATTICE.dimensions, space_threshold=1,
                warn=True)
            all_LR_percolating = lattice_analysis_utils.flag_percolating_clusters(
                corrected_LR_cluster_positions, self.LATTICE.dimensions, space_threshold=3,
                type_grid=self.LATTICE.type_grid,
                LR_table=self.Hamiltonian.LR_residue_interaction_table,
                SLR_table=self.Hamiltonian.SLR_residue_interaction_table,
                warn=True)
            cluster_percolating = [all_percolating[i] for i in big_cluster_idx]
            LR_cluster_percolating = [all_LR_percolating[i] for i in big_LR_cluster_idx]

        # for each set of positions get the polymeric properties associated with each whole cluster using the single image convention corrected values
        cluster_polymeric_properties_list      = lattice_analysis_utils.extract_cluster_polymeric_properties(big_clusters)
        LR_cluster_polymeric_properties_list   = lattice_analysis_utils.extract_cluster_polymeric_properties(big_clusters_LR)

        # for each cluster calculate the volume, surface area and density (requires corrected cluster positions)
        cluster_size_properties     = lattice_analysis_utils.compute_cluster_gross_properties(big_clusters)
        LR_cluster_size_properties  = lattice_analysis_utils.compute_cluster_gross_properties(big_clusters_LR)

        # for each cluster calculate the radial density profile IF the cluster contains more than 27 beads (3x3x3). We should probably make this number a keyfile
        # value
        
        cluster_radial_density     = lattice_analysis_utils.compute_cluster_radial_density_profile(big_clusters, self.LATTICE.dimensions, minimum_cluster_size_in_beads = CONFIG.RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD, hardwall=self.hardwall)
        LR_cluster_radial_density  = lattice_analysis_utils.compute_cluster_radial_density_profile(big_clusters_LR, self.LATTICE.dimensions, minimum_cluster_size_in_beads = CONFIG.RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD, hardwall=self.hardwall)
        cluster_radial_density_indices = [
            idx + 1 for idx, cluster in enumerate(big_clusters)
            if len(cluster) >= CONFIG.RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD]
        LR_cluster_radial_density_indices = [
            idx + 1 for idx, cluster in enumerate(big_clusters_LR)
            if len(cluster) >= CONFIG.RADIAL_DENSITY_PROFILE_BEAD_THRESHOLD]

        # blank out the undefined quantities of any percolating cluster. Rg and
        # asphericity become nan (self-describing, and it poisons a naive mean()
        # loudly rather than dragging it down quietly); the hull trio becomes -1,
        # which is already the documented "undefined" convention for a degenerate
        # hull; the radial profile row is dropped entirely, which the C<n> labels
        # already tolerate because they index the CLUSTERS.dat column explicitly.
        (cluster_polymeric_properties_list, cluster_size_properties,
         cluster_radial_density, cluster_radial_density_indices) = _blank_percolating_clusters(
            cluster_percolating, cluster_polymeric_properties_list, cluster_size_properties,
            cluster_radial_density, cluster_radial_density_indices)

        (LR_cluster_polymeric_properties_list, LR_cluster_size_properties,
         LR_cluster_radial_density, LR_cluster_radial_density_indices) = _blank_percolating_clusters(
            LR_cluster_percolating, LR_cluster_polymeric_properties_list, LR_cluster_size_properties,
            LR_cluster_radial_density, LR_cluster_radial_density_indices)

        # ... and put a line in the run's own record for EVERY analysis step at
        # which it happened. The UserWarning raised down in cluster_utils is
        # documented lemonade API and stays, but Python's default filter shows it
        # once per (message, location) per process, so a long run that percolates
        # in every frame produced at most a couple of stderr lines and nothing at
        # all in log.txt. The logging lives here rather than in cluster_utils
        # because that module is imported by lemonade, where writing a log.txt
        # into the working directory would be a surprise.
        percolating_report = []
        if any(cluster_percolating):
            percolating_report.append("short-range cluster(s) %s" % (
                ', '.join('C%i' % (i + 1) for i, f in enumerate(cluster_percolating) if f)))
        if any(LR_cluster_percolating):
            percolating_report.append("long-range cluster(s) %s" % (
                ', '.join('C%i' % (i + 1) for i, f in enumerate(LR_cluster_percolating) if f)))

        if percolating_report:
            msg = ("Step %i: %s percolate the periodic box - Rg/asphericity written as nan, "
                   "hull volume/area/density as -1, radial profile omitted" % (
                       step, ' and '.join(percolating_report)))
            IO_utils.status_message(msg, 'warning')
            pimmslogger.log_warning(msg)

        # We'll leave the following in as a sanity check

        # remove soon - > for debugging
        """
        count=0
        for idx in range(0, len(corrected_cluster_positions)):
            pdb_utils.write_positions_to_file(corrected_cluster_positions[idx], 'clusters/%i_cluster_%i_CORR.pdb'%(step,count))
            pdb_utils.write_positions_to_file(cluster_positions[idx], 'clusters/%i_cluster_%i_UNCORR.pdb'%(step,count))
            count=count+1
        """            
        
        # write cluster list
        analysis_IO.write_clusters(step, clusters, self.LATTICE.chainIDtoType)
        analysis_IO.write_LR_clusters(step, LR_clusters, self.LATTICE.chainIDtoType)

        # write cluster size/shape analysis
        analysis_IO.write_cluster_properties(
            step, cluster_polymeric_properties_list, cluster_size_properties,
            cluster_radial_density, cluster_radial_density_indices)
        analysis_IO.write_LR_cluster_properties(
            step, LR_cluster_polymeric_properties_list, LR_cluster_size_properties,
            LR_cluster_radial_density, LR_cluster_radial_density_indices)


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_polymeric_properties(self, step):
        """
        Run polymeric-properties (radius of gyration and asphericity) analysis.

        Computes the radius of gyration and asphericity for every chain and writes
        both lists to disk on each call, making this routine I/O heavy.

        Parameters
        ----------
        step : int
            Current simulation step number, written alongside the data.

        Returns
        -------
        None
        """

        RG_list = []
        asph_list = []
        
        for chainID in sorted(self.LATTICE.chains.keys()):            
            tmp  = self.LATTICE.chains[chainID].analysis_get_polymeric_properties()

            RG_list.append(tmp[0])
            asph_list.append(tmp[1])

        analysis_IO.write_radius_of_gyration(step, RG_list)
        analysis_IO.write_asphericity(step, asph_list)


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_end_to_end(self, step):
        """
        Run end-to-end distance analysis.

        Computes the end-to-end distance for every chain and writes the list to
        disk on each call, making this routine I/O heavy.

        Parameters
        ----------
        step : int
            Current simulation step number, written alongside the data.

        Returns
        -------
        None
        """

        e2e_list = []
        
        for chainID in sorted(self.LATTICE.chains.keys()):
            e2e_list.append(self.LATTICE.chains[chainID].analysis_get_end_to_end_distance())

        analysis_IO.write_end_to_end(step, e2e_list)


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_acceptance(self, step):
        """
        Run acceptance-criterion analysis.

        Writes the current move acceptance statistics (held by the
        :class:`AcceptanceCalculator`, ``self.ACC``) to disk on each call.

        Parameters
        ----------
        step : int
            Current simulation step number, written alongside the statistics.

        Returns
        -------
        None
        """
        analysis_IO.write_acceptance_statistics(step, self.ACC)

        
        
    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_custom_stubb(self, step):
        """
        No-op custom-analysis stub.

        Used as the ``ANA_CUSTOM`` analysis routine when no custom analysis module
        is side-loaded via the keyfile. Does nothing.

        Parameters
        ----------
        step : int
            Current simulation step number (ignored).

        Returns
        -------
        None
        """
        pass


    #-----------------------------------------------------------------
    #       
    def ANAFUNCT_save_restart(self, step):
        """
        Write a restart file capturing the current system state.

        Builds a :class:`restart.RestartObject` from the current lattice (recording
        the hardwall status), stores the freshly evaluated total energy in it, and
        writes the restart file to disk. The resulting file can be used to
        re-initialize the system in a later simulation, and holds everything
        ``RESTART_CONTINUE`` needs to resume this run exactly: the step, both
        global generator states, the temperature in force, and the settings a
        continuation must not change (the temperature the Hamiltonian was built
        at and the quench settings). ``RESTART_CONTINUE`` also checks the stored
        energy against the energy it computes for the configuration.

        The trajectory is flushed first, so the checkpoint never records a step
        whose frame is not yet on disk. ``run_all_analysis`` calls this after
        every other analysis due on the step, for the same reason.

        Parameters
        ----------
        step : int
            The master step the checkpoint is taken at (after that step's move
            and its scheduled output). Stored in the file.

        Returns
        -------
        None
        """
        IO_utils.status_message("Writing restart file on step %i..." %(step),'info', allow_suppress=True)

        # everything the checkpoint says has happened must be on disk before it.
        # lattice_utils.flush_xtc_writer accepts whatever the run holds (None, a
        # closed writer, the SAVE_AT_END buffer) and is idempotent; fall back to
        # the writer's own flush() where that helper is not available.
        _writer = getattr(self, 'xtc_writer', None)
        _flush_writer = getattr(lattice_utils, 'flush_xtc_writer', None)
        if callable(_flush_writer):
            _flush_writer(_writer)
        else:
            _flush = getattr(_writer, 'flush', None)
            if callable(_flush):
                _flush()

        R = restart.RestartObject()

        # build using lattice, and also pass the hardwall status of the current simulation
        R.build_from_lattice(self.LATTICE, self.hardwall)

        # evaluate the total energy and provide this as well (as a plain Python
        # number: a numpy scalar would tie the pickle to the numpy that wrote it)
        (energy, _, _, _, _) = self.Hamiltonian.evaluate_total_energy(self.LATTICE)
        if isinstance(energy, np.generic):
            energy = energy.item()
        R.set_energy(energy)

        # what a continuation must not change: the temperature the Hamiltonian's
        # ANGLE_PENALTY_T_NORM scaling was built at, and the quench settings
        # (read from the resolved keywords, not from self.QUENCH_RUN: that flag
        # is switched off once the ramp has reached QUENCH_END, and a checkpoint
        # written after that is still a checkpoint of a quench run)
        _keywords = getattr(self, 'keyword_lookup', None) or {}
        _quench = {'QUENCH_RUN': bool(_keywords.get('QUENCH_RUN', getattr(self, 'QUENCH_RUN', False)))}
        if _quench['QUENCH_RUN']:
            for _key in ('QUENCH_START', 'QUENCH_END', 'QUENCH_STEPSIZE', 'QUENCH_FREQ'):
                _value = _keywords[_key] if _key in _keywords else getattr(self, _key)
                _quench[_key] = _value.item() if isinstance(_value, np.generic) else _value

        # and everything a later run needs to resume from exactly here: the step,
        # both global generators (the compiled kernels are reseeded from the
        # Python generator every megamove, so its state is all that matters) and
        # the temperature in force
        R.set_continuation_state(step, random.getstate(), np.random.get_state(),
                                 self.ACC.temperature, pimms_version=_pimms_version(),
                                 equilibrium_temperature=_keywords.get('EQUILIBRIUM_TEMPERATURE'),
                                 quench=_quench)

        # output restart file to disk
        R.write_to_file()

        # the step of the checkpoint now on disk; reported if the run is interrupted
        self._last_checkpoint_step = step

        
        
