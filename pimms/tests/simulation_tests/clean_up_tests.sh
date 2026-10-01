#!/bin/sh
# Remove generated simulation outputs from every test_* fixture directory.
# Neither the pytest harness (conftest.py) nor run_up_test_sims.sh runs inside
# these directories any more, so this is only needed after running a fixture
# by hand, or to clear what an older version of the harness left behind.
for d in test_*/; do
    (
        cd "$d" || exit 1
        rm -f -- *.dat *.pdb traj.xtc eq_traj.xtc restart.pimms log.txt \
            pytest_*_log.txt run_log.txt keyfile_used.kf
    )
done
