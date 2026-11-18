#!/bin/sh
# Delete everything the demo's commands write, so that it can be run again
# from the start:
#
#   * what pdb2pimms writes here (the restart file, the report and the freeze
#     file) and the figure analyze_contacts.py writes here;
#   * what PIMMS writes in each run directory;
#   * the contact tables and B-factor structures analyze_contacts.py writes in
#     each run directory.
#
# The inputs (1kx5.cif.gz, params.prm, the keyfiles, the scripts, README.md and
# assets/) are kept. Every file is deleted by its exact name, never by a
# pattern, so a file of your own in a run directory is left alone: anything
# other than KEYFILE.kf still in a run directory at the end is listed, and the
# script then exits with status 1.
#
# Usage, from anywhere:
#
#     ./reset.sh        delete the outputs
#     ./reset.sh -n     list what would be deleted, and delete nothing

set -u

# the run directories, and the added proteins as analyze_contacts.py names
# them in file names (a run directory holds files for the ones it adds)
RUNS='H1_CTD ProTa H1_CTD_ProTa'
PROTEINS='H1_CTD ProTa'

# written in this directory by pdb2pimms and analyze_contacts.py
TOP_LEVEL_OUTPUTS='1kx5_restart.pimms 1kx5_report.md 1kx5_freeze.txt nucleosome_contacts.pdf'

# written in a run directory by PIMMS with the demo keyfiles' settings
RUN_OUTPUTS='START.pdb traj.xtc restart.pimms log.txt keyfile_used.kf parameters_used.prm
absolute_energies_of_angles.txt ACCEPTANCE.dat MOVE_FREQS.dat TOTAL_MOVES.dat PERFORMANCE.dat
ENERGY.dat RG.dat ASPH.dat END_TO_END_DIST.dat
CLUSTERS.dat NUM_CLUSTERS.dat CLUSTER_AREA.dat CLUSTER_ASPH.dat CLUSTER_DEN.dat
CLUSTER_RADIAL_DENSITY_PROFILE.dat CLUSTER_RG.dat CLUSTER_VOL.dat
LR_CLUSTERS.dat NUM_LR_CLUSTERS.dat LR_CLUSTER_AREA.dat LR_CLUSTER_ASPH.dat LR_CLUSTER_DEN.dat
LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat LR_CLUSTER_RG.dat LR_CLUSTER_VOL.dat'

# written in a run directory by PIMMS once per chain type, as CHAIN_0_<name>,
# CHAIN_1_<name>, ... (seven chain types with one added protein, eight with two)
PER_CHAIN_TYPE_OUTPUTS='CLUSTERS.dat LR_CLUSTERS.dat DISTANCE_MAP.dat INTSCAL.dat INTSCAL_SQUARED.dat
SCALING_INFORMATION.dat'

# files in a run directory that are not outputs and are never reported
RUN_INPUTS='KEYFILE.kf .DS_Store'

USAGE='usage: reset.sh [-n]'


# The outputs of one run directory that exist, one path per line. The
# per-chain-type files are counted up from chain type 0 until a chain type has
# none, so a run with more added proteins is covered too.
run_outputs() {
    run=$1
    for name in $RUN_OUTPUTS; do
        [ -f "$run/$name" ] && echo "$run/$name"
    done
    for protein in $PROTEINS; do
        for name in "contacts_nucleosome_$protein.csv" "contacts_$protein.csv" "1kx5_contacts_$protein.pdb"; do
            [ -f "$run/$name" ] && echo "$run/$name"
        done
    done
    chain_type=0
    while :; do
        found=0
        for name in $PER_CHAIN_TYPE_OUTPUTS; do
            if [ -f "$run/CHAIN_${chain_type}_$name" ]; then
                echo "$run/CHAIN_${chain_type}_$name"
                found=1
            fi
        done
        [ "$found" -eq 1 ] || break
        chain_type=$((chain_type + 1))
    done
}


# Every output that exists, one path per line.
all_outputs() {
    for name in $TOP_LEVEL_OUTPUTS; do
        [ -f "$name" ] && echo "$name"
    done
    for run in $RUNS; do
        [ -d "$run" ] && run_outputs "$run"
    done
    return 0
}


# What a run directory holds that is neither an input nor in the list of
# outputs given on standard input, one path per line.
unexpected_files() {
    outputs=$(cat)
    for run in $RUNS; do
        [ -d "$run" ] || continue
        ls -A "$run" | while IFS= read -r entry; do
            case " $RUN_INPUTS " in
                *" $entry "*) continue ;;
            esac
            if ! printf '%s\n' "$outputs" | grep -Fxq -- "$run/$entry"; then
                echo "$run/$entry"
            fi
        done
    done
}


dry_run=0
case $# in
    0) ;;
    1) if [ "$1" = '-n' ]; then dry_run=1; else echo "$USAGE" >&2; exit 2; fi ;;
    *) echo "$USAGE" >&2; exit 2 ;;
esac

# work in the demo directory wherever the script is run from, so that the
# relative names above can never point anywhere else. CDPATH is cleared for
# this cd: with an exported CDPATH, a relative directory name is looked up
# there first, and the script would then delete another directory's outputs
CDPATH= cd -- "$(dirname -- "$0")" > /dev/null || exit 1

outputs=$(all_outputs)
if [ -z "$outputs" ]; then
    count=0
else
    count=$(printf '%s\n' "$outputs" | wc -l | tr -d ' ')
fi

if [ "$dry_run" -eq 1 ]; then
    [ "$count" -gt 0 ] && printf '%s\n' "$outputs"
    echo "reset.sh -n: $count file(s) would be deleted"
else
    if [ "$count" -gt 0 ]; then
        printf '%s\n' "$outputs" | while IFS= read -r path; do
            rm -- "$path" || exit 1
        done || exit 1
    fi
    echo "reset.sh: deleted $count file(s)"
fi

# with -n the outputs are still there, so they are left out of this by name
leftovers=$(printf '%s\n' "$outputs" | unexpected_files)
if [ -n "$leftovers" ]; then
    if [ "$dry_run" -eq 1 ]; then
        echo "reset.sh -n: these files are not outputs this script knows, so they would be left in place:" >&2
    else
        echo "reset.sh: these files are not outputs this script knows, so they were left in place:" >&2
    fi
    printf '%s\n' "$leftovers" | sed 's/^/  /' >&2
    exit 1
fi
exit 0
