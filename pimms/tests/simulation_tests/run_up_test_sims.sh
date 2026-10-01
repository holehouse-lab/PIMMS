#!/bin/sh
# Re-run every test_* fixture simulation with the PIMMS from THIS working tree,
# then regenerate the expected-output baselines from those runs. Run from this
# directory:
#
#     sh run_up_test_sims.sh [RUN_ROOT [EXPECTED_ROOT]]
#
# Like the pytest harness (conftest.py), this does not run inside the fixture
# directories: the input files of each test_N are copied to RUN_ROOT/test_N and
# PIMMS runs there, so the tree is left as it was. RUN_ROOT defaults to a new
# temporary directory (its path is printed) and must not already hold test_N
# directories. EXPECTED_ROOT is where the baselines are written; it defaults to
# expected_output/ here, which is what the tests read. Give another directory to
# look at the regenerated baselines before replacing the tracked ones.
set -e
HERE=$(pwd)
REPO=$(cd ../../.. && pwd)
export PYTHONPATH="$REPO"
RUN_ROOT=${1:-$(mktemp -d)}
EXPECTED_ROOT=${2:-$HERE/expected_output}
mkdir -p "$RUN_ROOT"
RUN_ROOT=$(cd "$RUN_ROOT" && pwd)
echo "running the fixtures in $RUN_ROOT"
for d in test_*/; do
    d=${d%/}
    [ -f "$d/KEYFILE.kf" ] || continue
    if [ -e "$RUN_ROOT/$d" ]; then
        echo "$RUN_ROOT/$d already exists - give an empty RUN_ROOT" >&2
        exit 1
    fi
    mkdir "$RUN_ROOT/$d"
    # Copy the inputs and leave behind any output of an earlier run made by hand
    # in the fixture directory: the rule of _is_generated_output in conftest.py.
    for f in "$d"/*; do
        [ -f "$f" ] || continue
        case "${f##*/}" in
            KEYFILE*|*.kf|*.prm) cp "$f" "$RUN_ROOT/$d/" ;;
            restart.pimms|*.dat|*.xtc|*.pdb|*.txt|*.lat) ;;
            *) cp "$f" "$RUN_ROOT/$d/" ;;
        esac
    done
    echo "--- running $d ---"
    (cd "$RUN_ROOT/$d" && python "$REPO/scripts/PIMMS" -k KEYFILE.kf > run_log.txt 2>&1)
done
python generate_expected_output.py --allow-missing --tests-root "$RUN_ROOT" --expected-root "$EXPECTED_ROOT"
