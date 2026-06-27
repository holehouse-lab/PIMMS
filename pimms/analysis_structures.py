## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2024
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
            from 1 to ``seqlen - 2`` inclusive.

        Returns
        -------
        None
        """

        self.internal_scaling = {}
        self.initialized = False
        self.count = 0

        for i in range(1,seqlen-1):
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
            self.internal_scaling[i] = (self.internal_scaling[i]*self.count + IS[i])/(self.count+1)

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

    def write_status(self, filename='INTSCAL.dat'):
        """
        Write the current mean internal scaling profile to file.

        Parameters
        ----------
        filename : str, optional
            Output filename (default ``'INTSCAL.dat'``). Overwritten if it
            already exists.

        Returns
        -------
        None
        """
        with open(filename, 'w') as fh:
            for i in self.internal_scaling:
                fh.write('%i\t%4.4f \n' %(i, self.internal_scaling[i]))

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
            created for each sequence-separation gap from 1 to ``seqlen - 2``
            inclusive.

        Returns
        -------
        None
        """

        self.internal_scaling_squared = {}
        self.initialized = False
        self.count = 0

        for i in range(1,seqlen-1):
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
            self.internal_scaling_squared[i] = (self.internal_scaling_squared[i]*self.count + (IS[i]*IS[i]))/(self.count+1)

        # increment the count
        self.count = self.count+1


    def print_status(self):
        """
        Print the current mean internal scaling squared profile to stdout.

        Returns
        -------
        None
        """
        for i in self.internal_scaling_squared:
            print('%i\t%4.4f' %(i, self.internal_scaling_squared[i]))

    def write_status(self, filename='INTSCAL.dat'):
        """
        Write the current mean internal scaling squared profile to file.

        Parameters
        ----------
        filename : str, optional
            Output filename (default ``'INTSCAL.dat'``). Overwritten if it
            already exists.

        Returns
        -------
        None
        """
        with open(filename, 'w') as fh:
            for i in self.internal_scaling_squared:
                fh.write('%i\t%4.4f \n' %(i, self.internal_scaling_squared[i]))

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
        always discarded, and at most 40 log-spaced points are used for the fit.

        Returns
        -------
        tuple of float
            ``(nu, R0)`` where ``nu`` is the fitted scaling exponent and ``R0``
            the prefactor. Returns ``(-1, -1)`` if the chain is too short
            (fewer than 25 internal-scaling gaps) to fit meaningfully.
        """

        # if the chain is shorter than 25 residues then don't bother doing
        # any kind of scaling analysis - too finite
        if len(self.get_internal_scaling_array()) < 25:
            return (-1,-1)
    
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

        # finally, identfy the indices that are used for fitting
        logspaced_idx = []
        for i in range(0,num_fitting_points):
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

        Parameters
        ----------
        seqlen : int
            Length (number of residues) of the chain being analysed. A
            ``(seqlen, seqlen)`` matrix is allocated; only the upper-right
            triangle is ever populated.

        Returns
        -------
        None
        """

        # create a square matrix for the distance map. Note we'll only
        # populate the upper right triangle
        #
        # O O O O O O O O O O
        # * O O O O O O O O O
        # * * O O O O O O O O
        # * * * O O O O O O O
        # * * * * O O O O O O
        # * * * * * O O O O O
        # * * * * * * O O O O
        # * * * * * * * O O O
        # * * * * * * * * O O
        # * * * * * * * * * O
        # * * * * * * * * * *
        #
        # O filled 
        # * empty
        #

        self.distance_map = np.zeros((seqlen,seqlen),dtype=float)
        self.initialized = False
        self.seqlen = seqlen
        self.count = 0
        

    def update_distance_map(self, dMap):
        """
        Fold an instantaneous distance map into the running mean distance map.

        Each matrix element is updated as a running average so that the stored
        matrix always holds the mean over all snapshots seen so far. Only the
        upper-right triangle is ever populated (the lower triangle stays zero).

        Parameters
        ----------
        dMap : numpy.ndarray
            Square instantaneous distance map with the same shape as the stored
            matrix.

        Returns
        -------
        None

        Raises
        ------
        AnalysisStructureException
            If ``dMap`` is not a numpy array, or its shape does not match the
            stored distance map.
        """

        if not type(dMap) == np.ndarray:
            raise AnalysisStructureException('ERROR: Passed the update distance map function a matrix but was not a numpy array')
            
        if not dMap.shape == self.distance_map.shape:
            raise AnalysisStructureException('ERROR: Distance map to update and newly generated distance maps do not match in size')
            
        # update over the full square, but updating 0 with 0 is still zero so only the upper
        # right triangle will ever get filled
        for i in range(0, self.seqlen):
            for j in range(0, self.seqlen):
                self.distance_map[i][j] =(self.distance_map[i][j]*self.count + dMap[i][j])/(self.count+1)
                
        # increment the count
        self.count = self.count+1


    def write_status(self, filename='DISTANCE_MAP.dat'):
        """
        Write the current mean distance map out to file.

        Each row is written as comma-separated values, one row per source
        residue.

        Parameters
        ----------
        filename : str, optional
            Output filename (default ``'DISTANCE_MAP.dat'``). Overwritten if it
            already exists.

        Returns
        -------
        None
        """
        with open(filename, 'w') as fh:
            for i in range(0, self.seqlen):

                for j in range(0, self.seqlen-1):
                    fh.write('%4.4f, ' % self.distance_map[i][j])
                fh.write('%4.4f\n' % self.distance_map[i][self.seqlen-1])


    def get_distance_map(self):
        """
        Return the current system-average distance map.

        Returns
        -------
        numpy.ndarray
            The ``(seqlen, seqlen)`` running-mean distance map array (only the
            upper-right triangle is populated).
        """
        return self.distance_map
                            

                    
                


            
            
