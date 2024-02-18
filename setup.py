"""
PIMMS (Polymer Interactions in Multicomponent Mixtures)
Lattice simulation package for biomolecule
Author: Alex Holehouse
Developed by the Holehouse and Pappu labs
Copyright 2015 - 2024
"""

import sys
from setuptools import setup, find_packages
import versioneer

# ................................
# added for cython construction (Nov 2018)
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
# ................................


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
        
    )
    
]

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='pimms',
    author='Alex Holehouse',
    author_email='alex.holehouse@wustl.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='LGPLv3',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # external modules
    ext_modules = cythonize(extensions, language_level="3"),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,
    scripts=['scripts/PIMMS'], 

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.8",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    zip_safe=False,

)
