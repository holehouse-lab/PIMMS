## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Author: Alex Holehouse
## Developed by the Holehouse and Pappu labs
## Copyright 2015 - 2024
## 
## ...........................................................................

# Config file to define the types used in the cython code. This allows us to
# toggle the memory usage of the code by changing the types used. Note that
# an equivalent change must be made in the CONFIG file as well.
#
#


cimport numpy as cnp
import numpy as np

#from numpy cimport int16_t as NUMPY_INT16_TYPE
#ctypedef NUMPY_INT16_TYPE  NUMPY_INT_TYPE

## TOGGLE THIS TO MATCH THE TYPE IN THE CONFIG FILE
#ctypedef cnp.int16_t NUMPY_INT_TYPE
ctypedef cnp.int64_t NUMPY_INT_TYPE


# this should not be changed
ctypedef cnp.int64_t NUMPY_INT_TYPE_long
