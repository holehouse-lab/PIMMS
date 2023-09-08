## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2023
## ...........................................................................

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy



setup(
    ext_modules = cythonize("*.pyx"), 
    include_dirs=[numpy.get_include()],
    extra_compile_args=["-O3"],
)
