"""
PIMMS (Polymer Interactions in Multicomponent Mixtures)

Build script for the compiled Cython extensions. All package metadata (name,
version, dependencies, classifiers, ...) lives in ``pyproject.toml``; this file
exists only to declare the native extension modules - which cannot yet be expressed
in ``pyproject.toml`` - and the ``PIMMS`` command-line script. The version is
supplied automatically at build time by the versioningit plugin.

Author: Alex Holehouse
Developed by the Holehouse and Pappu labs
Copyright 2015 - 2026
"""

import os
import sys
from setuptools import setup, find_packages

# ................................
# added for cython construction (Nov 2018)
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
# ................................


def openmp_flags():
    """Return (compile_args, link_args, include_dirs, library_dirs) for OpenMP.

    The optimized parallel kernel (pimms.mega_crank_fast) uses
    cython.parallel.prange (OpenMP). Apple clang needs Homebrew libomp; Linux
    gcc/clang use -fopenmp. If no OpenMP is found the extension still builds and
    runs - prange simply executes serially - so this degrades gracefully and
    never breaks the install.
    """
    if sys.platform == "darwin":
        for prefix in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"):
            if os.path.isdir(prefix):
                return (["-Xpreprocessor", "-fopenmp", f"-I{prefix}/include"],
                        [f"-L{prefix}/lib", "-lomp"],
                        [f"{prefix}/include"], [f"{prefix}/lib"])
        return ([], [], [], [])
    return (["-fopenmp"], ["-fopenmp"], [], [])


_omp_c, _omp_l, _omp_inc, _omp_libdirs = openmp_flags()


# nb -O2 and -O3 do nothing or make performance worse
# on ARM CPU at least...
#
extensions = [
    Extension(
        "pimms.get_randmax",
        ["pimms/get_randmax.pyx"],
        include_dirs=[numpy.get_include()],
        
    ),

    Extension(
        "pimms.hyperloop",
        ["pimms/hyperloop.pyx"],
        include_dirs=[numpy.get_include()],
        
    ),

    Extension(
        "pimms.inner_loops",
        ["pimms/inner_loops.pyx"],
        include_dirs=[numpy.get_include()], 
        
    ),

    Extension(
        "pimms.inner_loops_hardwall",
        ["pimms/inner_loops_hardwall.pyx"],
        include_dirs=[numpy.get_include()], 
        
    ),

    Extension(
        "pimms.lattice_tools",
        ["pimms/lattice_tools.pyx"],
        include_dirs=[numpy.get_include()], 
        
    ),

    Extension(
        "pimms.mega_crank",
        ["pimms/mega_crank.pyx"],
        include_dirs=[numpy.get_include()],

    ),

    # optimized drop-in: allocation-free serial kernel + OpenMP parallel kernel
    Extension(
        "pimms.mega_crank_fast",
        ["pimms/mega_crank_fast.pyx"],
        include_dirs=[numpy.get_include()] + _omp_inc,
        library_dirs=_omp_libdirs,
        extra_compile_args=_omp_c,
        extra_link_args=_omp_l,
    ),

    Extension(
        "pimms.mega_crank_2D",
        ["pimms/mega_crank_2D.pyx"],
        include_dirs=[numpy.get_include()], 
        
    ),

    Extension(
        "pimms.random_number",
        ["pimms/random_number.pyx"],
        include_dirs=[numpy.get_include()],
        
    ),
    Extension(
        "pimms.system_utils",
        ["pimms/system_utils.pyx"],
        include_dirs=[numpy.get_include()],

    ),

    Extension(
        "pimms.cluster_kernels",
        ["pimms/cluster_kernels.pyx"],
        include_dirs=[numpy.get_include()],

    ),

    # lemonade analysis backend
    Extension(
        "pimms.lemonade.kernels._pbc",
        ["pimms/lemonade/kernels/_pbc.pyx"],
        include_dirs=[numpy.get_include()],

    )

]

setup(
    # Python importable packages (pimms, pimms.lemonade, ...) discovered automatically.
    packages=find_packages(),

    # Compiled Cython extension modules (the one thing that must live in setup.py).
    ext_modules=cythonize(extensions, language_level="3"),

    # Ship the non-Python package data listed in MANIFEST.in.
    include_package_data=True,

    # The PIMMS command-line executable (a script file, not a console entry point).
    scripts=['scripts/PIMMS'],

    # Compiled extensions cannot run from a zipped egg.
    zip_safe=False,
)
