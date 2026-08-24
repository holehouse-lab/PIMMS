#!/bin/sh
# Remove generated simulation outputs from every test_* fixture directory.
# (conftest.py does this automatically before each run; this script is for
# manual cleanup.)
for d in test_*/; do
    (
        cd "$d" || exit 1
        rm -f -- *.dat *.pdb traj.xtc eq_traj.xtc restart.pimms log.txt \
            pytest_*_log.txt run_log.txt
    )
done
