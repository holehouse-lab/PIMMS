## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2021
## ...........................................................................

import sys

from . import parameterfile_parser

FILENAME='GCF_parameters/CULSAC_GCF_5pC.prm'
EF = parameterfile_parser.parse_energy(FILENAME)

TABLE=EF[0]

AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','0']

for aa1 in AA:
    for aa2 in AA:
        sys.stdout.write('%s, ' %str(TABLE[aa1][aa2]))

    print("")
