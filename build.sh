#!/bin/zsh
set +e

# Simplified build script for PIMMS
#


# delete any existing build files 
for i in pimms/*so; do rm $i; done
for i in pimms/*c; do rm $i; done

# this command will install locally 
pip install -e . --upgrade --force-reinstall
