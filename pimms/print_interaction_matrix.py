## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

import sys

from . import parameterfile_parser


# Developer script: dump a parsed interaction table as CSV. Run by hand with
# `python -m pimms.print_interaction_matrix` from a directory containing the
# hardcoded parameter file below. Guarded so importing the module (it ships in
# the package) does not try to read a file that only exists on one machine.
if __name__ == "__main__":
    FILENAME='GCF_parameters/CULSAC_GCF_5pC.prm'
    EF = parameterfile_parser.parse_energy(FILENAME)

    TABLE=EF[0]

    AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','0']

    for aa1 in AA:
        for aa2 in AA:
            sys.stdout.write('%s, ' %str(TABLE[aa1][aa2]))

        print("")
