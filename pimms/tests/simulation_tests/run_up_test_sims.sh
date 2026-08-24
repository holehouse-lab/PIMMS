#!/bin/sh
# Clean and re-run every test_* fixture simulation with the PIMMS from THIS
# working tree, then regenerate the expected-output baselines. Run from this
# directory. (The old version of this script was an accidental copy of
# clean_up_tests.sh - it never ran anything - and both stopped at test_13.)
set -e
REPO=$(cd ../../.. && pwd)
export PYTHONPATH="$REPO"
sh clean_up_tests.sh
for d in test_*/; do
    d=${d%/}
    [ -f "$d/KEYFILE.kf" ] || continue
    echo "--- running $d ---"
    (cd "$d" && python "$REPO/scripts/PIMMS" -k KEYFILE.kf > run_log.txt 2>&1)
done
python generate_expected_output.py --allow-missing
