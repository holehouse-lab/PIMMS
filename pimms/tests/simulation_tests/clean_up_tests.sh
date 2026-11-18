#!/bin/sh
# Remove generated simulation outputs from every test_* fixture directory.
# Neither the pytest harness (conftest.py) nor run_up_test_sims.sh runs inside
# these directories any more, so this is only needed after running a fixture
# by hand, or to clear what an older version of the harness left behind.
#
# The names are what a PIMMS run writes next to its keyfile: the .dat tables,
# structures, trajectories, the restart snapshot, the status log, the echoed
# keyfile and parameter file (keyfile_used.kf, parameters_used.prm), the
# angle-energy table and the per-test logs the old harness kept. A file git
# tracks is never removed, whatever its name: test_11 carries a tracked copy of
# absolute_energies_of_angles.txt, and deleting it would only leave a deletion
# in the working tree for the next commit to pick up.
for d in test_*/; do
    (
        cd "$d" || exit 1
        for f in *.dat *.pdb traj.xtc eq_traj.xtc restart.pimms log.txt \
                 pytest_*_log.txt run_log.txt keyfile_used.kf \
                 parameters_used.prm absolute_energies_of_angles.txt; do
            [ -f "$f" ] || continue
            if git ls-files --error-unmatch -- "$f" >/dev/null 2>&1; then
                continue
            fi
            rm -f -- "$f"
        done
    )
done
