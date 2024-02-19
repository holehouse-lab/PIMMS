## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2024
## ...........................................................................



##
## Not used yet, but might be a good model for organization going forward
##

import os
from pimms.latticeExceptions import UnfinishedCodeException, KeyFileException
from pimms import pimmslogger

class AnalysisSettings:
    
    def __init__(self, cluster_threshold):
        self.cluster_threshold = cluster_threshold


class FreezeFile:


    # ...........................................................................
    #
    def __init__(self, filename):
        """
        Class for reading and storing information from a freeze file. This is a file
        which can be used to specify which chains or beads are to be frozen in the
        simulation. 

        Parameters
        ----------
        filename : str
            The name of the file to be read 


        """

        # check that the file exists
        if not os.path.isfile(filename):
            raise KeyFileException(f'Unable to find FREEZE FILE. Passed filename is: {filename}. Please verify the file actually exists')

                    
        # initialize the chains and beads lists
        chains = []
        beads = []

        # open the file and read the contents
        with open(filename, 'r') as f:
            content = f.readlines()


        # cycle through each line
        for idx, line in enumerate(content):


            # strip the line of whitespace
            sline = line.strip()

            # skip empty lines
            if len(sline) == 0:
                continue

            # skip comment lines
            if sline[0] == '#':
                continue

            # discard comments at the end of the line as well, if they exist
            sline = sline.split('#')[0].strip()


            # if this line is reporting on chains (C)
            if sline[0] == 'C':

                try:
                    local_chains = [int(i) for i in sline[1:].split()]
                except:
                    raise ValueError(f'Error parsing chains in freeze file on line {idx}: {line}')

                chains.extend(local_chains)
            
            # if this line is reporting on beads [NOT YET IMPLEMENTED]
            if sline[0] == 'B':

                try:
                    local_beads = [int(i) for i in sline[1:].split(' ')]
                except:
                    raise ValueError(f'Error parsing chains in freeze file on line {idx}: {line}')

                beads.extend(local_chains)
                raise UnfinishedCodeException('Beads not yet implemented for freezeing')

        # remove duplicates
        self._chains = list(set(chains))
        self._beads = list(set(beads))
        self._filename = filename

    # ...........................................................................
    #
    @property
    def chains(self):
        return self._chains

    # ...........................................................................
    #
    @property
    def beads(self):
        return self._beads

    # ...........................................................................
    #    
    @property
    def filename(self):
        return self._filename


    # ...........................................................................
    #
    def validate_freeze_file(self, latticeObject):
        """
        Function to validate that the chains and beads specified in the freeze file
        are actually present in the lattice object.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object to be validated against

        Returns
        -------
        None
            No return variable, but an exception is raised if the freeze file is not valid

        """

        # for each chain specified in the freeze file
        for chainID in self.chains:

            # if the chain is not present in the lattice object, raise an exception
            if chainID not in latticeObject.chains:
                
                raise KeyFileException(f"\n\nFreeze file {self.filename} specifies chain {chainID} which is present in the lattice object. Lattice object chains are {list(latticeObject.chains.keys())} while freeze file chains are {self.chains}.")
                

    # ...........................................................................
    #
    def log_freeze_file(self):
        """
        Function to log the chains and beads specified in the freeze file

        """

        pimmslogger.log_status(f"Freeze file             : {self.filename}", timestamp=False)
        pimmslogger.log_status(f"Number of frozen chains : {len(self.chains)}", timestamp=False)
        pimmslogger.log_status(f"Frozen chainIDs         : {str(self.chains)}", timestamp=False)

        
    # ...........................................................................
    #            
    def __str__(self):
        s = 'Freeze file:\n'
        s = s + 'Chains: %s\n'%self.chains
        s = s + 'Beads: %s\n'%self.beads
        return s

    # ...........................................................................
    #
    def __repr__(self):
        return self.__str__()
    
