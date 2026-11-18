#!/bin/zsh
#
# Remove PIMMS run-output files (trajectories, energies, analysis .dat files,
# restart snapshots, status logs, the PDB files a run writes, and the echoed
# parameter / angle files) from the repo root and every demo / validation
# directory. Input keyfiles (KEYFILE.kf), parameter files (params.prm), any
# *new_restart.pimms inputs and any other PDB file (a structure for pdb2pimms,
# say) are left in place. Safe to run repeatedly, from any directory.
#
# demos/nucleosome_from_pdb runs in its own subdirectories and is not swept here;
# it has its own reset.sh.

setopt NULL_GLOB        # demo_keyfiles/*/ etc. expand to nothing if absent (no error)

# the directories below are relative to the repo root, whatever directory the
# script was launched from; without this, '.' was wherever it was called from
cd "${0:A:h}" || exit 1

# every directory that may hold a PIMMS run
run_dirs=( . demo_keyfiles/*/ demos/*/ pimms/fast_kernels/validation_*/ )

for d in $run_dirs; do
    [[ -d "$d" ]] || continue
    # -maxdepth 1 keeps this to files directly in the run dir; an exact
    # 'restart.pimms' match leaves input files like new_restart.pimms untouched,
    # and only the PDB names PIMMS writes are matched, never '*.pdb'.
    find "$d" -maxdepth 1 -type f \( \
            -name '*.dat' \
        -o  -name '*.xtc' \
        -o  -name 'START.pdb' \
        -o  -name 'eq_START.pdb' \
        -o  -name 'CONFIG_AT_ENERGY_FAIL.pdb' \
        -o  -name 'restart.pimms' \
        -o  -name 'log.txt' \
        -o  -name 'parameters_used.prm' \
        -o  -name 'keyfile_used.kf' \
        -o  -name 'absolute_energies_of_angles.txt' \
    \) -delete
done

echo "Cleaned PIMMS output files from the repo root and all demo / validation directories."
