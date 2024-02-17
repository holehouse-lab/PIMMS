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

# this should not be changed
ctypedef cnp.int64_t NUMPY_INT_TYPE_long

## TOGGLE THIS TO MATCH THE TYPE IN THE CONFIG FILE
## NOTE that if you change this value you must change
## the corresponding PYTHON config in CONFIG.py

# valid values are cnp.int8_t cnp.int16_t cnp.int32_t cnp.int64_t
ctypedef cnp.int32_t NUMPY_INT_TYPE


