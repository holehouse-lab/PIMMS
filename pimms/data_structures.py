## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................



##
## Not used yet, but might be a good model for organization going forward
##

import os
from pimms.latticeExceptions import UnfinishedCodeException, KeyFileException
from pimms import pimmslogger

class AnalysisSettings:
    """
    Lightweight container for on-the-fly analysis configuration.

    Holds settings used by the analysis machinery during a simulation,
    currently just the cluster-size threshold.
    """

    def __init__(self, cluster_threshold):
        """
        Initialize the analysis settings container.

        Parameters
        ----------
        cluster_threshold : int
            Threshold cluster size used by clustering analysis routines.

        """
        self.cluster_threshold = cluster_threshold


class FreezeFile:
    """
    Reader/container for a simulation freeze file.

    Parses a freeze file (see :meth:`__init__` for the file format) and
    exposes the set of frozen chain IDs (and, in future, bead IDs) so
    that the simulation can hold those chains fixed during sampling.
    """


    # ...........................................................................
    #
    def __init__(self, filename):
        """
        Class for reading and storing information from a freeze file. This is a file
        which can be used to specify which chains or beads are to be frozen in the
        simulation. 

                Expected freeze file structure
                ------------------------------
                The freeze file is plain text with one directive per line.

                Supported directives:
                - `C <chain_id> <chain_id> ...`
                    Freeze one or more chain IDs. Multiple `C` lines are allowed.

                Notes:
                - Empty lines are ignored.
                - Lines beginning with `#` are ignored.
                - Inline comments are allowed after `#`.
                - Duplicate chain IDs are allowed in the file and are deduplicated internally.

                Current limitation:
                - `B ...` bead-level directives are recognized but not implemented and will
                    raise `UnfinishedCodeException`.

                Example:
                # Freeze two groups of chains
                C 0 1 2
                C 10 11

                # Inline comments are also supported
                C 42  # freeze chain 42

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
                except ValueError:
                    raise ValueError(f'Error parsing chains in freeze file on line {idx}: {line}')

                chains.extend(local_chains)
            
            # if this line is reporting on beads [NOT YET IMPLEMENTED]
            if sline[0] == 'B':

                try:
                    # split() with no argument (as the 'C' branch above does), NOT
                    # split(' '): the latter keeps the empty field produced by the
                    # space right after the 'B', so int('') blew up on well-formed
                    # input and reported it as a parse error
                    local_beads = [int(i) for i in sline[1:].split()]
                except ValueError:
                    raise ValueError(f'Error parsing beads in freeze file on line {idx}: {line}')

                beads.extend(local_beads)
                raise UnfinishedCodeException('Beads not yet implemented for freezeing')

        # remove duplicates
        self._chains = list(set(chains))
        self._beads = list(set(beads))
        self._filename = filename

    # ...........................................................................
    #
    @property
    def chains(self):
        """
        list of int : The deduplicated chain IDs to be frozen.
        """
        return self._chains

    # ...........................................................................
    #
    @property
    def beads(self):
        """
        list of int : The deduplicated bead IDs to be frozen (bead-level
        freezing is not yet implemented, so this is typically empty).
        """
        return self._beads

    # ...........................................................................
    #    
    @property
    def filename(self):
        """
        str : The path of the freeze file that was read.
        """
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
                
                raise KeyFileException(f"\n\nFreeze file {self.filename} specifies chain {chainID}, which is NOT present in the lattice object. Lattice object chains are {list(latticeObject.chains.keys())} while freeze file chains are {self.chains}.")
                

    # ...........................................................................
    #
    def log_freeze_file(self):
        """
        Write a summary of the freeze file to the simulation log.

        Logs the freeze file path, the number of frozen chains, and the
        list of frozen chain IDs as STATUS entries.

        Returns
        -------
        None
            No return value; entries are appended to the simulation log.

        """

        pimmslogger.log_status(f"Freeze file             : {self.filename}", timestamp=False)
        pimmslogger.log_status(f"Number of frozen chains : {len(self.chains)}", timestamp=False)
        pimmslogger.log_status(f"Frozen chainIDs         : {str(self.chains)}", timestamp=False)

        
    # ...........................................................................
    #            
    def __str__(self):
        """
        Return a human-readable summary of the frozen chains and beads.

        Returns
        -------
        str
            Multi-line string listing the frozen chain and bead IDs.

        """
        s = 'Freeze file:\n'
        s = s + 'Chains: %s\n'%self.chains
        s = s + 'Beads: %s\n'%self.beads
        return s

    # ...........................................................................
    #
    def __repr__(self):
        """
        Return the same summary string as :meth:`__str__`.

        Returns
        -------
        str
            Multi-line string listing the frozen chain and bead IDs.

        """
        return self.__str__()
    
