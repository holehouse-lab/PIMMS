## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................

## analysis_structures
##
## Objects defined in this analysis_structures.py file are analytical
## objects where a system-averaged value is useful at the end of the
## simulation. For such analyses there is typically a subset of logic
## involved in dealing with the associated underlying data structure.
## This file allows that data-structure and the associated code to be
## defined independently of anything else.
##

import numpy as np
from . import numpy_utils
from .latticeExceptions import AnalysisStructureException

# Thresholds above which the distance-map analysis warns before it allocates
# anything (see distance_map_cost and Simulation.ANAFUNCT_distance_map): an
# estimated peak memory of 1 GiB, or an output file of 256 MB.
DISTANCE_MAP_WARN_MEMORY_BYTES = 1 << 30
DISTANCE_MAP_WARN_FILE_BYTES = 256 * 1000 * 1000


class InternalScaling:
    """
    InternalScaling analysis is an analysis with provides insight into 
    the degree of expansion of the chain, but only makes sense in the
    the context of a full simulation (i.e. an instantaneous Internal
    Scaling value is not particularly useful).

    """

    def __init__(self, seqlen):
        """
        Initialize the running internal scaling accumulator.

        Parameters
        ----------
        seqlen : int
            Length (number of residues) of the chain being analysed. One
            internal-scaling bin is created for each sequence-separation gap
            from 1 to ``seqlen - 1`` inclusive.

        Returns
        -------
        None
        """

        self.internal_scaling = {}
        self.initialized = False
        self.count = 0

        for i in range(1, seqlen):
            self.internal_scaling[i] = 0


    def update_internal_scaling(self, IS):
        """
        Fold an instantaneous internal scaling profile into the running mean.

        Each per-gap value is updated as a running average so that
        ``self.internal_scaling`` always holds the mean over all snapshots
        seen so far.

        Parameters
        ----------
        IS : dict
            Instantaneous internal scaling values keyed by sequence-separation
            gap. Must have the same number of entries as the accumulator.

        Returns
        -------
        None

        Raises
        ------
        AnalysisStructureException
            If ``IS`` does not have the same length as the internal accumulator.
        """

        if not len(IS) == len(self.internal_scaling):
            raise AnalysisStructureException('ERROR: INTERNAL SCALING UPDATE')

        # update mean internal scaling to include current values (note if count = 0 this just
        # initializes the self.internal_scaling to the passed data
        for i in self.internal_scaling:
            self.internal_scaling[i] += (
                IS[i] - self.internal_scaling[i]) / (self.count + 1)

        # increment the count
        self.count = self.count+1


    def print_status(self):
        """
        Print the current mean internal scaling profile to stdout.

        Returns
        -------
        None
        """
        for i in self.internal_scaling:
            print('%i\t%4.4f' %(i, self.internal_scaling[i]))

    def get_internal_scaling_array(self):
        """
        Return the mean internal scaling values ordered by sequence separation.

        Returns
        -------
        list of float
            Mean internal scaling values ordered by increasing gap. The
            dictionary keys are sorted explicitly because dictionary iteration
            order is not guaranteed to be numerical.
        """
        ISGaps = list(self.internal_scaling.keys())
        
        # cannot assume the dictionary will return in numerical order
        # - it almost certainly will be that's not a fair assumption and
        # is not specified in the language        
        ISGaps.sort()

        ISArray = []
        for i in ISGaps:
            ISArray.append(self.internal_scaling[i])

        return ISArray


class InternalScalingSquared:
    """
    InternalScalingSquared analysis is an analysis with provides insight into 
    the degree of expansion of the chain, but only makes sense in the
    the context of a full simulation (i.e. an instantaneous Internal
    Scaling value is not particularly useful).

    """

    def __init__(self, seqlen):
        """
        Initialize the running internal scaling squared accumulator.

        Parameters
        ----------
        seqlen : int
            Length (number of residues) of the chain being analysed. One bin is
            created for each sequence-separation gap from 1 to ``seqlen - 1``
            inclusive.

        Returns
        -------
        None
        """

        self.internal_scaling_squared = {}
        self.initialized = False
        self.count = 0

        for i in range(1, seqlen):
            self.internal_scaling_squared[i] = 0


    def update_internal_scaling(self, IS):
        """
        Fold an instantaneous internal scaling profile into the running mean of
        the squared distances.

        Note that ``IS`` is expected to contain the instantaneous internal
        scaling distances; the value accumulated into the running mean is the
        square of each distance (``IS[i] * IS[i]``).

        Parameters
        ----------
        IS : dict
            Instantaneous internal scaling distances keyed by sequence-separation
            gap. Must have the same number of entries as the accumulator.

        Returns
        -------
        None

        Raises
        ------
        AnalysisStructureException
            If ``IS`` does not have the same length as the internal accumulator.
        """

        if not len(IS) == len(self.internal_scaling_squared):
            raise AnalysisStructureException('ERROR: INTERNAL SCALING UPDATE')

        # update mean internal scaling to include current values (note if count = 0 this just
        # initializes the self.internal_scaling to the passed data
        for i in self.internal_scaling_squared:

            # NOTE that the value we're adding is IS[i]*IS[i] - i.e. internal scaling squared
            squared = IS[i] * IS[i]
            self.internal_scaling_squared[i] += (
                squared - self.internal_scaling_squared[i]) / (self.count + 1)

        # increment the count
        self.count = self.count+1


    def update_internal_scaling_squared(self, IS_squared):
        """Fold already-averaged squared pair distances into the running mean.

        This is the scientifically correct update for an RMS internal-scaling
        profile. If the pair distances at one sequence gap are ``r_1 .. r_n``, the
        required instantaneous value is ``mean(r_i**2)``. Passing ``mean(r_i)`` to
        :meth:`update_internal_scaling` instead computes ``mean(r_i)**2``, which loses
        the within-snapshot variance. The legacy method remains available for callers
        that genuinely have one distance per gap; PIMMS's chain analysis uses this
        method with the pairwise second moment computed by
        :func:`pimms.lattice_analysis_utils.get_internal_scaling_profile`.

        Parameters
        ----------
        IS_squared : dict
            Instantaneous mean squared distances keyed by sequence-separation gap.
            Must have the same number of entries as the accumulator.

        Returns
        -------
        None
            No return value; the running mean and sample count are updated in place.

        Raises
        ------
        AnalysisStructureException
            If the profile length does not match the accumulator.
        """
        if len(IS_squared) != len(self.internal_scaling_squared):
            raise AnalysisStructureException('ERROR: INTERNAL SCALING UPDATE')

        for i in self.internal_scaling_squared:
            self.internal_scaling_squared[i] += (
                IS_squared[i] - self.internal_scaling_squared[i]) / (self.count + 1)

        self.count += 1


    def print_status(self):
        """
        Print the current mean internal scaling squared profile to stdout.

        Returns
        -------
        None
        """
        for i in self.internal_scaling_squared:
            print('%i\t%4.4f' %(i, self.internal_scaling_squared[i]))

    def get_internal_scaling_array(self):
        """
        Return the mean internal scaling squared values ordered by separation.

        Returns
        -------
        list of float
            Mean internal scaling squared values ordered by increasing gap. The
            dictionary keys are sorted explicitly because dictionary iteration
            order is not guaranteed to be numerical.
        """
        ISGaps = list(self.internal_scaling_squared.keys())
        
        # cannot assume the dictionary will return in numerical order
        # - it almost certainly will be that's not a fair assumption and
        # is not specified in the language        
        ISGaps.sort()

        ISArray = []
        for i in ISGaps:
            ISArray.append(self.internal_scaling_squared[i])

        return ISArray


    def fit_scaling_exponent(self):
        """
        Fit the polymer scaling exponent and prefactor from the mean profile.

        This method for extracting scaling relationships was developed to
        avoid the bias introduced by the fact that on a log scale, most
        inter-residue distances occupy the top-right part of the fitting
        regime; the idea is to shift to approximately evenly spaced points in
        log space for the linear fit. The first 15 sequence-separation gaps are
        always discarded, and the fit uses at most 41 log-spaced points (41
        rather than 40 because the largest sequence separation, the most
        informative point for the fit, is always included).

        Returns
        -------
        tuple of float
            ``(nu, R0)`` where ``nu`` is the fitted scaling exponent and ``R0``
            the prefactor. Returns ``(-1, -1)`` if the chain is too short
            (fewer than 25 internal-scaling gaps) to fit meaningfully, or if no
            internal-scaling sample has been accumulated yet.
        """

        # if the chain is shorter than 25 residues then don't bother doing
        # any kind of scaling analysis - too finite
        if len(self.get_internal_scaling_array()) < 25:
            return (-1,-1)

        # if no internal-scaling sample was ever folded in (e.g. ANA_INTSCAL exceeded
        # the production length) the profile is all zeros; log(0) -> -inf would make
        # polyfit return (nan, nan) and write NaNs to SCALING_INFORMATION.dat.
        if self.count == 0:
            return (-1, -1)
    
        # always discard for 15 residues!
        scaling_array = self.get_internal_scaling_array()[15:]
        seq_sep_vals   = np.arange(1,len(self.get_internal_scaling_array())+1)[15:]
    
        # next find indices for evenly spaced points in logspace
        if len(seq_sep_vals) > 40:
            num_fitting_points = 40
        else:
            num_fitting_points = len(seq_sep_vals)

        # this section basically identifies the indices that provide
        # a linearly spaces dataset in logspace
        y_data = np.log(seq_sep_vals)
        y_data_offset = y_data - y_data[0]
        interval = y_data_offset[-1]/num_fitting_points
        integer_vals = y_data_offset/interval

        # finally, identfy the indices that are used for fitting. Note the range is
        # inclusive of num_fitting_points so the LARGEST sequence separation - the
        # most informative point for a scaling fit - is always part of the fit set
        # (the exclusive range dropped the end of the fitting window).
        logspaced_idx = []
        for i in range(0, num_fitting_points + 1):
            [local_ix,_] = numpy_utils.find_nearest(integer_vals, i) 

            # if we already found this point then skip...
            if local_ix in logspaced_idx:
                continue
            else:
                logspaced_idx.append(local_ix)

        # defines the x and y values used for log linear fitting
        fitting_separation = [seq_sep_vals[i] for i in logspaced_idx]
        fitting_distances  = [np.sqrt(scaling_array[i]) for i in logspaced_idx]

        # do fitting and extract value
        out = np.polyfit(np.log(fitting_separation), np.log(fitting_distances), 1)
        nu = out[0]
        R0 = np.exp(out[1])
    

        return (nu, R0)

                              
def distance_map_cost(seqlen: int, n_chains: int) -> tuple[int, int]:
    """Estimate what the distance-map analysis costs for one chain type.

    Every chain keeps a ``seqlen x seqlen`` float64 running mean, so one "map"
    is ``8 * seqlen ** 2`` bytes and ``n_chains`` of them are held from the
    first sample to the end of the run. On top of those the analysis needs
    about two more maps of working memory, at two moments:

    * while a chain is sampled, its instantaneous map plus the temporaries of
      the distance calculation, which works in blocks of at most about 54 MB;
    * when the run ends and the per-chain means of the type are averaged, the
      running sum and the averaged map (the maps are added one at a time and
      are not stacked).

    The estimate returned is therefore ``(n_chains + 2)`` maps. Measured
    through the real analysis routines it is exact for the end-of-run average,
    and the sampling peak is ``(n_chains + 1)`` maps plus the block temporaries,
    so for chains shorter than about 2,600 beads (one map below 54 MB) the true
    peak is reached while sampling and exceeds the estimate by at most that
    block, which is small against the 1 GiB the warning starts at. The output
    file holds ``seqlen ** 2`` entries written as ``%4.4f`` plus a tab, about 8
    bytes each for a long chain.

    Parameters
    ----------
    seqlen : int
        Number of beads in a chain of this type.

    n_chains : int
        Number of chains of this type.

    Returns
    -------
    tuple of int
        ``(peak_memory_bytes, file_bytes)``, both estimates.

    """
    one_map = 8 * int(seqlen) * int(seqlen)
    peak_memory = (int(n_chains) + 2) * one_map
    file_bytes = one_map + int(seqlen)
    return (peak_memory, file_bytes)


class DistanceMap:
    """
    Distance map analysis is an analysis with provides insight into 
    the long-range interaction on a residue-by-residue level - essentially
    can be considered a contactmap which lacks cutoffs and instead computes the
    average distance between two residues

    """

    def __init__(self, seqlen):
        """
        Initialize a running square inter-residue distance map.

        Nothing of size ``seqlen x seqlen`` is allocated here. Every chain owns
        one of these objects from the moment it is built, whether or not the
        distance map is ever sampled, and the matrix is 8 bytes x seqlen^2: the
        accumulator is created by the first :meth:`update_distance_map`, so a
        run with ``ANA_DISTMAP`` switched off (or one that never reaches a
        sampling step) does not pay for it.

        Parameters
        ----------
        seqlen : int
            Length (number of residues) of the chain being analysed. The
            running mean is a full ``(seqlen, seqlen)`` symmetric matrix.

        Returns
        -------
        None
        """

        # allocated on the first update (see above)
        self.distance_map = None
        self.initialized = False
        self.seqlen = seqlen
        self.count = 0
        

    def update_distance_map(self, dMap, consume=False):
        """
        Fold an instantaneous distance map into the running mean distance map.

        Each matrix element is updated as a running average so that the stored
        matrix always holds the full symmetric mean over all snapshots seen so far.

        Parameters
        ----------
        dMap : numpy.ndarray
            Square instantaneous distance map with the same shape as the stored
            matrix.

        consume : bool, optional
            If True, ``dMap`` is used as scratch space and holds garbage on
            return. A caller that built the instantaneous map only to pass it
            here (as :meth:`pimms.chain.Chain.analysis_update_distance_map`
            does) saves one ``seqlen x seqlen`` temporary that way. It needs a
            float64 map; anything else is left untouched and a temporary is
            used. Default is False (``dMap`` is not modified).

        Returns
        -------
        None

        Raises
        ------
        AnalysisStructureException
            If ``dMap`` is not a numpy array, or its shape does not match the
            stored distance map.
        """

        if type(dMap) is not np.ndarray:
            raise AnalysisStructureException('ERROR: Passed the update distance map function a matrix but was not a numpy array')
            
        if not dMap.shape == (self.seqlen, self.seqlen):
            raise AnalysisStructureException('ERROR: Distance map to update and newly generated distance maps do not match in size')

        if self.distance_map is None:
            self.distance_map = np.zeros((self.seqlen, self.seqlen), dtype=float)
            
        # Numerically stable in-place running mean,
        #     mean += (new - mean) / (count + 1),
        # evaluated one operation at a time into a single work array. This
        # avoids both the old O(seqlen^2) Python loop and the two full temporary
        # matrices of the one-line form, and gives bit-for-bit the same mean.
        if consume and dMap.dtype == self.distance_map.dtype:
            work = dMap
            np.subtract(dMap, self.distance_map, out=work)
        else:
            work = dMap - self.distance_map
        work /= (self.count + 1)
        self.distance_map += work

        # increment the count
        self.count = self.count+1


    def get_distance_map(self):
        """
        Return the current system-average distance map.

        Returns
        -------
        numpy.ndarray
            The full symmetric ``(seqlen, seqlen)`` running-mean distance map.
            Before the first update this is a matrix of zeros (built on
            request, not stored).
        """
        if self.distance_map is None:
            return np.zeros((self.seqlen, self.seqlen), dtype=float)
        return self.distance_map
