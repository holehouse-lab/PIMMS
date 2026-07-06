#!/bin/zsh
#
# Remove PIMMS run-output files (trajectories, energies, analysis .dat files,
# restart snapshots, status logs, and the echoed parameter / angle files) from
# the repo root and every demo / validation directory. Input keyfiles
# (KEYFILE.kf), parameter files (params.prm) and any *new_restart.pimms inputs
# are left in place. Safe to run repeatedly.

setopt NULL_GLOB        # demo_keyfiles/*/ etc. expand to nothing if absent (no error)

# every directory that may hold a PIMMS run
run_dirs=( . demo_keyfiles/*/ demos/*/ pimms/fast_kernels/validation_*/ )

for d in $run_dirs; do
    [[ -d "$d" ]] || continue
    # -maxdepth 1 keeps this to files directly in the run dir; an exact
    # 'restart.pimms' match leaves input files like new_restart.pimms untouched.
    find "$d" -maxdepth 1 -type f \( \
            -name '*.dat' \
        -o  -name '*.xtc' \
        -o  -name '*.pdb' \
        -o  -name 'restart.pimms' \
        -o  -name 'log.txt' \
        -o  -name 'parameters_used.prm' \
        -o  -name 'absolute_energies_of_angles.txt' \
    \) -delete
done

echo "Cleaned PIMMS output files from the repo root and all demo / validation directories."
