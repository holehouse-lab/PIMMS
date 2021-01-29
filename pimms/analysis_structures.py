## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2021
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
        
        self.internal_scaling = {}
        self.initialized = False
        self.count = 0

        for i in range(1,seqlen-1):
            self.internal_scaling[i] = 0


    def update_internal_scaling(self, IS):

        if not len(IS) == len(self.internal_scaling):
            raise AnalysisStructureException('ERROR: INTERNAL SCALING UPDATE')

        # update mean internal scaling to include current values (note if count = 0 this just 
        # initializes the self.internal_scaling to the passed data
        for i in self.internal_scaling:
            self.internal_scaling[i] = (self.internal_scaling[i]*self.count + IS[i])/(self.count+1)

        # increment the count
        self.count = self.count+1


    def print_status(self):
        for i in self.internal_scaling:
            print('%i\t%4.4f' %(i, self.internal_scaling[i]))

    def write_status(self, filename='INTSCAL.dat'):
        with open(filename, 'w') as fh:
            for i in self.internal_scaling:
                fh.write('%i\t%4.4f \n' %(i, self.internal_scaling[i]))

    def get_internal_scaling_array(self):
        """

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
        
        self.internal_scaling_squared = {}
        self.initialized = False
        self.count = 0

        for i in range(1,seqlen-1):
            self.internal_scaling_squared[i] = 0


    def update_internal_scaling(self, IS):
        """
        NOTE that this expects IS to the the instantaneous internal scaling distance 
        but 

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
        for i in self.internal_scaling_squared:
            print('%i\t%4.4f' %(i, self.internal_scaling_squared[i]))

    def write_status(self, filename='INTSCAL.dat'):
        with open(filename, 'w') as fh:
            for i in self.internal_scaling_squared:
                fh.write('%i\t%4.4f \n' %(i, self.internal_scaling_squared[i]))

    def get_internal_scaling_array(self):
        """

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
        This methdod for extracting scaling relationships was developed to 
        avoid the bias introduced by the fact that on a log scale, most inter-residue
        distances occupy the top-right part of the fitting regime, so the
        idea is the shift to approximately evenly spaced points in log space for
        the linear fitting.
        
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
        Accepts a square distance map matrix which is then used to update the existing matrix

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
        Write this distance map out to file

        """
        with open(filename, 'w') as fh:
            for i in range(0, self.seqlen):

                for j in range(0, self.seqlen-1):
                    fh.write('%4.4f, ' % self.distance_map[i][j])
                fh.write('%4.4f\n' % self.distance_map[i][self.seqlen-1])


    def get_distance_map(self):
        """
        Returns the numpy array containing the current system average
        distance map information

        """
        return self.distance_map
                            

                    
                


            
            
