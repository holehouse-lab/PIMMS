#!/usr/bin/env python
## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2024
## ...........................................................................


"""
pimms
Lattice simulation package for biomolecules

"""

from importlib.metadata import PackageNotFoundError, version

try:
	__version__ = version("pimms")
except PackageNotFoundError:
	# Local source-tree import before install.
	__version__ = "0+unknown"

# Git revision is no longer injected by versioneer.
__git_revision__ = "unknown"
