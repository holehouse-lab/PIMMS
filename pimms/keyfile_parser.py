## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2022
## ...........................................................................

##
## keyfile_parser
##
## The KeyFileParser object includes all the functionality necessary for parsing keyfiles
## for running PIMMS simulation
##
##
### Derived keywords
##
## EQUILIBRIUM_TEMPERATURE : Is set to the temperature that the simulation treats as the final equilibrium temperature.


import random
import sys
import os.path

from . import IO_utils
from . import latticeExceptions
from .latticeExceptions import KeyFileException, RestartException
from . import file_utilities
from . import restart
from . import pimmslogger
from . import CONFIG


def print_keyword_info():
    maxlen = 25

    for d in CONFIG.KEYWORD_DESCRIPTION:

        spacer = " "*(maxlen - len(d))

        print("%s%s | %s" % (d, spacer, CONFIG.KEYWORD_DESCRIPTION[d]))


# ===================================================================================================
#
#
class KeyFileParser:
    """
    KeyFileParser is essentially where all the logic that deals with input information is defined.

    Specifically, a KeyFileParser object can read a keyfile and extract any/all pertinant information
    there in. It then sets any default values which can be estimated but weren't defined explicitly.
    Finally, it sanity checks all the input to ensure we're doing something sensible.

    """



    #-----------------------------------------------------------------
    #    
    def __init__(self, filename, parse_only=False):
        """
        Function which initializes the keyfile parser object by defining the expected keywords and required keywords.

        These words serve different purposes:

        EXPECTED_KEYWORDS are words which we define default values for, and are expected to be included in the keyfile.
        
        REQUIRED_KEYWORDS are words which MUST be defined in the keyfile - i.e. we cannot make a reasonable guess without
                          essentially making a design decision about the simulation being run

        The REQUIRED_KEYWORDS are a subset of the EXPECTED_KEYWORDS. ONLY EXPECTED_KEYWORDS are correctly parsed - i.e. if
        an unexpected keyword is identified an exception is raised. This is a slightly over the top behaviour, but is 
        deliberate to avoid the scenario where you mis-type a keyword, don't realize, and the system over-writes with
        a default without you knowing.

        Parameters
        -----------------

        filename : str
            The only argument required for this object is the location of a keyfile which is to be parsed by the KeyFileParser
            object. 

        parse_only : bool 
            Optional keyword which - if set to true - means the keyfile is read in but nothing moreis done (i.e. no defaults set,
            no santization performed etc. This is useful if keyfile_parsers() is used in reading in pre-run keyfiles where rigerous
            assessment is not needed
           

        """
        
        # expected keywords contains a list of possible keywords which can be read form the keyfile. Importantly if these
        # keywords are *not* included then they are set to their default values
        self.expected_keywords = CONFIG.EXPECTED_KEYWORDS


        # required keywords are those which MUST be included in the keyfile - i.e. PIMMS can't set default values for these
        # keywords. Note that required_keywords is a subset of expected_keywords
        self.required_keywords = CONFIG.REQUIRED_KEYWORDS 

        # list of keywords that can support multiple entries in a keyfile
        self.keywords_with_multiple_entries = ['CHAIN', 'EXTRA_CHAIN', 'ANA_RESIDUE_PAIRS']
        self.keyword_lookup = {}
        self.DEFAULTS = {}


        ## IF PARSE ONLY mode
        if parse_only:
            self.parse(filename)        
            return 
            
        IO_utils.horizontal_line(hzlen=40, linechar='*')
        IO_utils.status_message("Parsing keyfile [%s]" %(filename),'startup')
        IO_utils.status_message("Default values set are explicitly announced below:", 'startup')

        print("")
        self.parse(filename)        # parse the keyfile (i.e. read in and deal with the file)
        self.assign_default()       # assigns default values to the internal system (though not YET to this keyfile)
        self.set_defaults()         # finally any values missing from the keyfile get set to the default values
        self.set_dynamic_defaults() # AND FINALLY update any values that depend on keyfile-derived parameters
        self.run_sanity_checks()    # run some sanity checks  
        self.add_derived_keywords() # some keywords are not explicitly included in the keyfile but are derived from
                                    # the keyfile words. These are set here

        IO_utils.horizontal_line(hzlen=40, linechar='*')

        # initialize logging...
        pimmslogger.initialize()



    #-----------------------------------------------------------------
    #    
    def __repr__(self):
        """
        """
        return str(self)


    #-----------------------------------------------------------------
    #    
    def __str__(self):
        msg = '\n............................\n'
        msg = msg + 'PIMMS Keyfile object:\n'
        msg = msg + '............................\n'
        for i in self.keyword_lookup:
            msg = msg + '%s => %s\n' %(i, str(self.keyword_lookup[i]))

        msg = msg + '............................\n'
        return msg


    #-----------------------------------------------------------------
    #    
    def __check_experimental_features(self, kw):
        """
        Should ONLY be used after the full keyfile has been parsed (i.e. 
        during the sanity check section).

        Parameters
        --------------
        kw : str
            Name of the keyword to be checked (i.e. a keyword where we required
            EXPERIMENTAL_FEATURE to be set to True for this parameter to be 
            useable.

        Returns
        ----------
        None
            No return type, but if EXPERIMENTAL_FEATURES is False this raises a
            KeyfileException with an appropriate error message

        """
        if self.keyword_lookup['EXPERIMENTAL_FEATURES'] == False:
            print(kw)
            raise KeyFileException(f'\n\nExperimental or non-supported keyword [{kw}] being proposed but EXPERIMENTAL_FEATURES is False.\n')

                       
    #-----------------------------------------------------------------
    #    
    def parse(self, filename, verbose=True):
        """
        Main function reads in a keyfile and extracts out the relevant details based on the keywords. Keywoerds are assigned to 
        the self.keywords_lookup

        Keywords must be defined as 

        KEYWORD : VALUE # comment 

        This function reads through all the lines in the keyfile and in a keyword-specific manner parses each keyword and assigns it
        to the self.keyword_lookup dictionary. Any keywords that are expected but NOT defined in the keyfile are then assigned as 
        default values later.

        Parameters
        --------------
        filename : str
            Name of the keyfile to be read

        verbose : bool (default = True)
            Flag which determines the level of warning messages to print. 
        

        Returns
        ---------
        None
            No return type but assigns values to the self.keyword_lookup dictionary, which itself is a key-value pair fo
            simulations keywords
        

        """
        
        # read the keyfile lines
        with open(filename, 'r') as fh:
            contents = fh.readlines()

        # for each line
        for line in contents:

            # if it's a comment line skip the whole line
            if file_utilities.is_comment_line(line):
                continue
                
            # remove comment section (at end of line)
            un_comment = file_utilities.remove_comments(line)

            # split based on the keyword/value separator
            splitline = un_comment.split(':')

            # if we find multiple keyword separators (':') characters 
            if len(splitline) > 2:
                raise KeyFileException(latticeExceptions.message_preprocess('On keyword %s - found multiple keyword-value separators...' % splitline[0].strip()))

            # if we didn't find a keyword separator
            if len(splitline) < 2:
                raise KeyFileException(latticeExceptions.message_preprocess('On line [%s] - no keyword separator found...' % line))

            # if get here must have keyword/keyvalue                
            putative_keyword = splitline[0].strip().upper()
            putative_value   = splitline[1].strip()
            
            if putative_keyword in self.expected_keywords:
                
                ## ** CHECK TO ENSURE WE DON'T OVERWRITE KEYWORDS **
                # check if we've seen this keyword before - if we're trying to overwrite raise an exception
                if putative_keyword in list(self.keyword_lookup.keys()):

                    if putative_keyword in self.keywords_with_multiple_entries:
                        # this is OK - we can have multiple chains!
                        pass
                    else:
                        raise KeyFileException(latticeExceptions.message_preprocess('Found a second occurence of the [%s] keyword. Please correct your keyfile and retry' % putative_keyword))

                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # CHAIN KEYWORD
                #
                # Chains are of the formation 
                #
                # CHAIN : X ABC  
                #
                # Where 
                #       X is the number of occurence (and integer)
                #       ABC is the polymer sequence of the chain
                #
                if putative_keyword == "CHAIN":

                    ## TO DO - we need to support sequence files (seq.in) files like ABSINTH such that 
                    ## we can define much more complex polymers. This will come later but for now is not
                    ## included.
                    
                    chainSplit = putative_value.split()
                    number_of_chains = int(chainSplit[0])
                    chain_sequence   = chainSplit[1].strip()

                    # if we already have a chain - note can have as many chain entries as we 
                    # want, where each chain entry is the number of that chain, followed by the 
                    # chain sequence
                    if putative_keyword in list(self.keyword_lookup.keys()):                    
                        self.keyword_lookup['CHAIN'].append([number_of_chains, chain_sequence])
                    else:
                        self.keyword_lookup['CHAIN'] = [[number_of_chains, chain_sequence]]                
            
                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # EXTRA_CHAIN KEYWORD
                #
                # Extra chains use exactly the same format as chains
                #
                # EXTRA_CHAIN : X ABC  
                #
                # Where 
                #       X is the number of occurence (and integer)
                #       ABC is the polymer sequence of the chain
                #
                # The designation of a chain as an EXTRA_CHAIN means additional
                # chains can be added to a simulation that's being restrated from
                # a known starting configuration. If no RESTART file is passed,
                # EXTRA_CHAINs are simply treated like normal chains but there's
                # no reason to use them otherwise
                elif putative_keyword == "EXTRA_CHAIN":

                    chainSplit = putative_value.split()
                    number_of_chains = int(chainSplit[0])
                    chain_sequence   = chainSplit[1].strip()

                    # if we already have a chain - note can have as many chain entries as we 
                    # want, where each chain entry is the number of that chain, followed by the 
                    # chain sequence
                    if putative_keyword in list(self.keyword_lookup.keys()):                    
                        self.keyword_lookup['EXTRA_CHAIN'].append([number_of_chains, chain_sequence])
                    else:
                        self.keyword_lookup['EXTRA_CHAIN'] = [[number_of_chains, chain_sequence]]                
            
                # DIMENSIONS KEYWORD (define if simulation is 2D or 3D)
                elif putative_keyword == 'DIMENSIONS':
                    self.keyword_lookup['DIMENSIONS'] = [int(i) for i in putative_value.split()]                                        
                    if not len(self.keyword_lookup['DIMENSIONS']) == 2 and not len(self.keyword_lookup['DIMENSIONS'])  == 3:
                        raise KeyFileException(latticeExceptions.message_preprocess('Unexpected number of dimensions [%s] ' % line))

                # conversion factor for PDB file writing
                elif putative_keyword == 'LATTICE_TO_ANGSTROMS':
                    self.keyword_lookup['LATTICE_TO_ANGSTROMS'] = float(putative_value)
                    
                # Dimensions of compressed equilibration box
                elif putative_keyword == 'RESIZED_EQUILIBRATION':
                    self.keyword_lookup['RESIZED_EQUILIBRATION'] = [int(i) for i in putative_value.split()]                                        


                # CASE_INSENSITIVE_CHAINS
                elif putative_keyword == "CASE_INSENSITIVE_CHAINS":
                    if putative_value.upper() == 'FALSE':
                        self.keyword_lookup['CASE_INSENSITIVE_CHAINS'] = False

                # HARDWALL
                elif putative_keyword == "HARDWALL":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['HARDWALL'] = True

                # EXPERIMENTAL_FEATURES
                elif putative_keyword == "EXPERIMENTAL_FEATURES":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['EXPERIMENTAL_FEATURES'] = True

                # TEMPERATURE
                elif putative_keyword == "TEMPERATURE":
                    self.keyword_lookup['TEMPERATURE'] = float(putative_value)
                    
                # N_STEPS
                elif putative_keyword == "N_STEPS":
                    self.keyword_lookup['N_STEPS'] = int(putative_value)
                    
                # equilibration
                elif putative_keyword == "EQUILIBRATION":
                    self.keyword_lookup['EQUILIBRATION'] = int(putative_value)

                # energy parameter 
                elif putative_keyword == "PARAMETER_FILE":
                    self.keyword_lookup['PARAMETER_FILE'] = str(putative_value)

                # PRINT_FREQUENCY
                elif putative_keyword == "PRINT_FREQ":
                    self.keyword_lookup['PRINT_FREQ'] = int(putative_value)

                # REDUCED_PRINTING 
                elif putative_keyword == "REDUCED_PRINTING":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['REDUCED_PRINTING'] = True

                # XTC_FREQUENCY
                elif putative_keyword == "XTC_FREQ":
                    self.keyword_lookup['XTC_FREQ'] = int(putative_value)

                # EN_FREQUENCY
                elif putative_keyword == "EN_FREQ":
                    self.keyword_lookup['EN_FREQ'] = int(putative_value)

                # SEED
                elif putative_keyword == "SEED":
                    self.keyword_lookup['SEED'] = int(putative_value)

                # ENERGY_CHECK
                elif putative_keyword == "ENERGY_CHECK":
                    self.keyword_lookup['ENERGY_CHECK'] = int(putative_value)

                # RESTART_FREQ
                elif putative_keyword == "RESTART_FREQ":
                    self.keyword_lookup['RESTART_FREQ'] = int(putative_value)
                    
                # RESTART FILE
                elif putative_keyword == "RESTART_FILE":
                    self.keyword_lookup['RESTART_FILE'] = str(putative_value)

                # RESTART OVERRIDE DIMENSIONS (this keywords are ONLY
                # used when we parse a restart file)
                elif putative_keyword == "RESTART_OVERRIDE_DIMENSIONS":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup["RESTART_OVERRIDE_DIMENSIONS"] = True
                    else:
                        self.keyword_lookup["RESTART_OVERRIDE_DIMENSIONS"] = False

                # RESTART OVERRIDE HARDWALL (this keywords are ONLY
                # used when we parse a restart file)
                elif putative_keyword == "RESTART_OVERRIDE_HARDWALL":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup["RESTART_OVERRIDE_HARDWALL"] = True
                    else:
                        self.keyword_lookup["RESTART_OVERRIDE_HARDWALL"] = False

                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ## QUENCHING keywords
                # Should we run a quenching simulation?
                elif putative_keyword == "QUENCH_RUN":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['QUENCH_RUN'] = True
                    else:
                        self.keyword_lookup['QUENCH_RUN'] = False
                    
                # Should the quench phase be treated as equillibration?
                elif putative_keyword == "QUENCH_AS_EQUILIBRATION":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['QUENCH_AS_EQUILIBRATION'] = True
                    else:
                        self.keyword_lookup['QUENCH_AS_EQUILIBRATION'] = False

                # Frequency (in MC steps) at which temperature is changed
                elif putative_keyword == "QUENCH_FREQ":
                    self.keyword_lookup['QUENCH_FREQ'] = int(putative_value)

                # Step size (in AU degrees) that the temperature should be stepping
                # note this is always a positive value (can be float)
                elif putative_keyword == "QUENCH_STEPSIZE":
                    self.keyword_lookup['QUENCH_STEPSIZE'] = abs(float(putative_value))

                # temperature at which the quench will start  
                elif putative_keyword == "QUENCH_START":
                    self.keyword_lookup['QUENCH_START'] = float(putative_value)

                # temperature at which the quench will end
                elif putative_keyword == "QUENCH_END":
                    self.keyword_lookup['QUENCH_END'] = float(putative_value)



                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ## TSMMC keywords , '', 'TSMMC_INTERPOLATION_MODE', 'TSMMC_NUMBER_OF_POINTS']
                elif putative_keyword == 'TSMMC_JUMP_TEMP':                           
                    self.keyword_lookup['TSMMC_JUMP_TEMP'] = float(putative_value)

                elif putative_keyword == 'TSMMC_STEP_MULTIPLIER':                           
                    self.keyword_lookup['TSMMC_STEP_MULTIPLIER'] = int(putative_value)

                elif putative_keyword == 'TSMMC_NUMBER_OF_POINTS':                           
                    self.keyword_lookup['TSMMC_NUMBER_OF_POINTS'] = int(putative_value)

                elif putative_keyword == 'TSMMC_FIXED_OFFSET':                           
                    self.keyword_lookup['TSMMC_FIXED_OFFSET'] = float(putative_value)

                elif putative_keyword == 'TSMMC_INTERPOLATION_MODE':       
                    self.keyword_lookup['TSMMC_INTERPOLATION_MODE'] = str(putative_value).upper().strip()
                    if self.keyword_lookup['TSMMC_INTERPOLATION_MODE'] not in ['LINEAR']:
                        raise KeyFileException(latticeExceptions.message_preprocess('Tried to set TSMMC_INTERPOLATION_MODE with unexpected keyword [%s]' % (self.keyword_lookup['TSMMC_INTERPOLATION_MODE'])))
                        


                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ## CRANKSHAFT keywords
        
                # Number of substeps to perform in the crankshaft system-shake
                # move
                elif putative_keyword == "CRANKSHAFT_SUBSTEPS":
                    self.keyword_lookup['CRANKSHAFT_SUBSTEPS'] = int(putative_value)

                # How are crankshaft moves dealt with?
                # > UNIFORM      -> every chain gets the CRANKSHAFT_SUBSTEPS number of moves
                # > PROPORTIONAL -> every chain gets CRANKSHAFT_SUBSTEPS * chain length number of moves
                # > PROP-SQUARED -> every chain gets CRANKSHATF_SUBSTEPS * (chain length)^2 number of moves
                # > PROP-CUBED   -> every chain gets CRANKSHATF_SUBSTEPS * (chain length)^3 number of moves
                elif putative_keyword == "CRANKSHAFT_MODE":
                    self.keyword_lookup['CRANKSHAFT_MODE'] = str(putative_value).upper().strip()
                    if self.keyword_lookup['CRANKSHAFT_MODE'] not in ['UNIFORM','PROPORTIONAL','PROP-SQUARED','PROP-CUBED']:
                        raise KeyFileException(latticeExceptions.message_preprocess('Tried to set CRANKSHAFT_MODE mode with unexpected keyword [%s]' % (self.keyword_lookup['CRANKSHAFT_MODE'])))

                # Non-interacting flag                              
                elif putative_keyword == "NON_INTERACTING":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['NON_INTERACTING'] = True
                    else:
                        self.keyword_lookup['NON_INTERACTING'] = False

                # If angles should to be used set ANGLES_OFF to true - means we don't have to provide
                # angle information in the parameter file (useful for prototyping stuff)
                elif putative_keyword == "ANGLES_OFF":
                    if putative_value.upper() == 'TRUE':
                        self.keyword_lookup['ANGLES_OFF'] = True
                    else:
                        self.keyword_lookup['ANGLES_OFF'] = False


                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ## Analysis keywords
                
                # Frequency at which cluster analysis is performed and written
                # to disk (can be IO intensive so keep low probably)
                elif putative_keyword == "ANALYSIS_FREQ":
                    self.keyword_lookup['ANALYSIS_FREQ'] = int(putative_value)
                    
                # frequency with which polymeric analysis is done
                elif putative_keyword == 'ANA_POL':
                    self.keyword_lookup[putative_keyword] = int(putative_value)
                    
                # frequency with which internal scaling analysis is done
                elif putative_keyword == 'ANA_INTSCAL':
                    self.keyword_lookup[putative_keyword] = int(putative_value)
                    
                # frequency with which distance map analysisis done
                elif putative_keyword == 'ANA_DISTMAP':
                    self.keyword_lookup[putative_keyword] = int(putative_value)
                    
                # frequency with which acceptance rate analysis is done
                elif putative_keyword == 'ANA_ACCEPTANCE':
                    self.keyword_lookup[putative_keyword] = int(putative_value)
                    
                # end-to-end and pairwise residues (as defined by ANA_RESIDUE_PAIRS)
                elif putative_keyword == 'ANA_INTER_RESIDUE':
                    self.keyword_lookup[putative_keyword] = int(putative_value)
                    
                # frequency with which cluster analysis is done
                elif putative_keyword == 'ANA_CLUSTER':
                    self.keyword_lookup[putative_keyword] = int(putative_value)

                # minimum number of chains for something to be considered a 'cluster'
                elif putative_keyword == 'ANA_CLUSTER_THRESHOLD':
                    self.keyword_lookup[putative_keyword] = int(putative_value)

                # frequency with which custom analysis routines are run
                elif putative_keyword == 'ANA_CUSTOM':
                    self.keyword_lookup[putative_keyword] = int(putative_value)

                # custom analysis code 
                elif putative_keyword == 'ANALYSIS_MODULE':
                    self.keyword_lookup[putative_keyword] = str(putative_value)
                    
                elif putative_keyword == 'ANA_RESIDUE_PAIRS':

                    # Slightly more complicated as we have to deal with multiple residue
                    # pairs we might be interested in...
                    split_residues = putative_value.split()
                    res1    = int(split_residues[0])
                    res2    = int(split_residues[1])

                    # if res1 is bigger than res 2 flip em
                    # around so the order is consistent i.e.
                    # always small - big
                    if res1 > res2:
                        tmp = res1
                        res1 = res2
                        res2 = tmp
                        
                    # if this isn't the first pair
                    if putative_keyword in self.keyword_lookup:
                        self.keyword_lookup['ANA_RESIDUE_PAIRS'].append([res1, res2])
                    else:
                        self.keyword_lookup['ANA_RESIDUE_PAIRS'] = [[res1, res2]]

                        
                ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ## MOVESET keywords
                
                # if we're defining the moveset
                elif putative_keyword[0:4] == 'MOVE':
                    if putative_keyword == 'MOVE_CRANKSHAFT':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_CHAIN_TRANSLATE':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_CHAIN_ROTATE':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_CHAIN_PIVOT':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_HEAD_PIVOT':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_SLITHER':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_CLUSTER_TRANSLATE':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_CLUSTER_ROTATE':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_CTSMMC':
                        self.keyword_lookup[putative_keyword] = float(putative_value)
                        
                    elif putative_keyword == 'MOVE_MULTICHAIN_TSMMC':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_RATCHET_PIVOT':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_SYSTEM_TSMMC':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    elif putative_keyword == 'MOVE_JUMP_AND_RELAX':
                        self.keyword_lookup[putative_keyword] = float(putative_value)

                    
                # End of move keywords...
                # -------------------------------------------
                    
                else:
                    raise KeyFileException(latticeExceptions.message_preprocess('Fail to deal with a supported keyword - [%s] - this is a bug! ' % putative_keyword))
                                    
            else:
                raise KeyFileException(latticeExceptions.message_preprocess('Found an unsupported keyword - [%s]. Valid supported keywords are\n%s ' % (putative_keyword, self.expected_keywords)))




    #-----------------------------------------------------------------
    #    
    def run_sanity_checks(self):
        """
        Function which, once the keyword has been parsed, will run an arbitrary number of 
        sanity checks to ensure things make sense. Adding new checks is simply a case of 
        adding a new 'CHECK' section in the code.

        This also sets a number of 'internal' keywords - keywords which are derived from the manually
        defined keywords but are NOT explicitly read in. This lets the system set certain configurations
        by deducing things from the manual keyfile that aid in decision making. To keep things clear a 
        list of the 'internal' keywords are below. Internal keywords always begin with a  '__'

        __TSMMC_USED # BOOLEAN --> set to true or false depending on if any TSMMC are being used.


        Specific details of each sanity check are provided as code blocks below. Before adding additional
        checks please read through these! Note also that restart file sanity checks are done in their
        own function which gets called from in here

        """

        ## ----------------------------------------------------------------------------------------
        ## First we case chain and extra chain keywords as 
        #
        if self.keyword_lookup['CASE_INSENSITIVE_CHAINS']:
            new_chains = []
            for entry in self.keyword_lookup['CHAIN']:
                new_chains.append([entry[0], entry[1].upper()])

            self.keyword_lookup['CHAIN'] = new_chains


            new_chains = []
            for entry in self.keyword_lookup['EXTRA_CHAIN']:
                new_chains.append([entry[0], entry[1].upper()])

            self.keyword_lookup['EXTRA_CHAIN'] = new_chains


        ## ------------------------------------------------------------------
        ## check values we think must be bigger than 0
        # 
        for c in ['TEMPERATURE', 'N_STEPS',  'PRINT_FREQ', 'XTC_FREQ', 'EN_FREQ', 'SEED', 'ENERGY_CHECK', 'RESTART_FREQ', 'QUENCH_STEPSIZE', 'QUENCH_START', 'QUENCH_END',  'TSMMC_STEP_MULTIPLIER', 'TSMMC_NUMBER_OF_POINTS',  'CRANKSHAFT_SUBSTEPS', 'ANALYSIS_FREQ', 'ANA_POL', 'ANA_DISTMAP', 'ANA_ACCEPTANCE', 'ANA_INTER_RESIDUE', 'ANA_CLUSTER', 'ANA_CUSTOM']:

            try:
                if self.keyword_lookup[c] <= 0:
                    raise KeyFileException(latticeExceptions.message_preprocess(f'Numerical error when parsing keyfile. Expected {c} to be larger than 0'))
            except TypeError:
                print(c)
                print(self.expected_keywords[c])
                raise Exception
        

        ## ------------------------------------------------------------------
        ## check values with think must be bigger than or equal to zero
        # 
        for c in ['ANA_CLUSTER_THRESHOLD', 'EQUILIBRATION']:

            try:
                if self.keyword_lookup[c] < 0:
                    raise KeyFileException(latticeExceptions.message_preprocess(f'Numerical error when parsing keyfile. Expected {c} to be larger than 0'))
            except TypeError:
                print(c)
                print(self.expected_keywords[c])
                raise Exception

        
        ## ------------------------------------------------------------------
        ## CHAIN CHECKS
        # if we have more than 26 unique chains print a warning about this
        tmp = []
        for c in self.keyword_lookup['CHAIN']:
            tmp.append(c[1])
        if len(set(tmp)) > 26:
            print(f"[ WARNING ] : Found {len(set(tmp))} unique chains (more than 26). This means the chain IDs for chains after 'Z' will all be set to 'Z'")


        ## ---------------------------------------------------------
        ## MOVESET CHECK
        # check moveset frequencies add up to 1
        running_total = 0.0
        message = ''                            # message defines a summary of all moves for easy debugging

        
        for i in self.expected_keywords:
            # for each MOVE keyword [note this dynamically finds them so if new MOVE_ keywords are added this will just work :-) ]
            if i[0:4] == "MOVE":
                running_total = running_total+self.keyword_lookup[i]
                message = message+'%s : %1.8f\n' % (i, self.keyword_lookup[i])

        if abs(running_total - 1.0) > 0.0000001 :            
            raise KeyFileException(latticeExceptions.message_preprocess('Moveset keywords do not add up to 1.0 (instead = %2.5e) - see below for specific details:\n%s' % (running_total, message)))
        ## ---------------------------------------------------------


        ## ---------------------------------------------------------
        ## QUENCH-RUN CHECKS
        # IF we're running a quench all the quench keywords must be included
        if self.keyword_lookup['QUENCH_RUN']:

            # first check all the keywords are there
            for i in ['QUENCH_FREQ', 'QUENCH_STEPSIZE', 'QUENCH_START', 'QUENCH_END', 'QUENCH_AS_EQUILIBRATION']:
                if 'UNSET' == self.keyword_lookup[i]:
                    raise KeyFileException('Trying to run a simulation with a quench but no [%s] defined in the keyfile' % i)
                    
            # the stepsize but be greater than zero!
            if not self.keyword_lookup['QUENCH_STEPSIZE'] > 0:
                raise KeyFileException('Trying to use a stepsize of zero for quenching simulation...')


            # check that the temperature distance being traversed isn't smaller than the temperature 
            # step size...
            dT = abs(self.keyword_lookup['QUENCH_START'] - self.keyword_lookup['QUENCH_END'])
            if dT < self.keyword_lookup['QUENCH_STEPSIZE']:
                raise KeyFileException('A single quench step overshoots the desired temperature: suggests an in the logfile')

            steps_for_quench = (1+(int(dT) / float(self.keyword_lookup['QUENCH_STEPSIZE'])))  * self.keyword_lookup['QUENCH_FREQ']
                
            if steps_for_quench >= self.keyword_lookup['N_STEPS']:
                raise KeyFileException('This quench will not complete as the quench period [%i] is longer than the number of steps in the simulation [%i]' %(steps_for_quench, self.keyword_lookup['N_STEPS'] ))
                                
            # If we're doing a heating simulation (so start is lower than end) then stepsize must become negative
            # as the temperature update operation is CURRENT - STEPSIZE
            if self.keyword_lookup['QUENCH_START'] < self.keyword_lookup['QUENCH_END']:
                self.keyword_lookup['QUENCH_STEPSIZE'] = -(self.keyword_lookup['QUENCH_STEPSIZE'])


            # Update the equilibration period so the quench process is enveloped by the quench (i.e. all
            # output is only at the target temperature)
            if self.keyword_lookup['QUENCH_AS_EQUILIBRATION'] :
                print("UPDATING EQUILIBRATION TO [%i] (temperature quench period is equilibration)" % steps_for_quench)
                self.keyword_lookup['EQUILIBRATION'] = steps_for_quench


            if not self.keyword_lookup['TEMPERATURE'] == self.keyword_lookup['QUENCH_START']:
                print("[ WARNING ] : Resetting the starting temperature to the QUENCH_START value [%i]"%self.keyword_lookup['QUENCH_START'])
                self.keyword_lookup['TEMPERATURE'] = self.keyword_lookup['QUENCH_START'] 

        ## ---------------------------------------------------------
        ## TSMMC check (sets the __TSMMC_ON keyword
        if self.keyword_lookup['MOVE_CTSMMC'] > 0 or self.keyword_lookup['MOVE_MULTICHAIN_TSMMC'] > 0 or self.keyword_lookup['MOVE_SYSTEM_TSMMC'] > 0:

            # 
            self.keyword_lookup['__TSMMC_USED'] = True

            # check if we didn't use a fixed OFFSET that the jumps make sense...
            if (self.keyword_lookup['TSMMC_JUMP_TEMP'] <= self.keyword_lookup['TEMPERATURE']) and self.keyword_lookup['TSMMC_FIXED_OFFSET'] is False:
                raise KeyFileException(latticeExceptions.message_preprocess('\n\nThe TSMMC jump temperature [%3.2f] is less than or equal to the actual simulation temperature [%3.2f], which will mean the TSMMC moves will at best hurt performance and at worst reduce sampling. Please correct your keyfile appropriately' % (self.keyword_lookup['TSMMC_JUMP_TEMP'], self.keyword_lookup['TEMPERATURE'])))

        else:
            self.keyword_lookup['__TSMMC_USED'] = False


        ## ------------------------------------------------------------
        ## residue distance pairs check - make sure non of the pair-pair
        ## analysis distances fall outside of the chain length
        for pair in self.keyword_lookup['ANA_RESIDUE_PAIRS']:            
            for chain in self.keyword_lookup['CHAIN']:
                if pair[1] >= len(chain[1]):
                    raise KeyFileException('Residue-residue distance analysis pair (%i) is outside the chain length (%i)' % (pair[1], len(chain[1])))
                    
        ## ------------------------------------------------------------
        ## Crash if cluster rotation moves are on and non-equal vertices 
        ## 
        dims = self.keyword_lookup['DIMENSIONS']
        if len(set(dims)) != 1 and self.keyword_lookup['MOVE_CLUSTER_ROTATE'] > 0:
            raise KeyFileException('CANNOT use a non-square or non-cubic box and use cluster rotation moves (dimensions = %s'%(str(dims)))


        ## Crash if dimensions are < 7 in 
        ## 
        dims = self.keyword_lookup['DIMENSIONS']
        if len(dims) == 2:
            if dims[0] < 7 or dims[1] < 7:
                raise KeyFileException('Box size is too small to correctly support super-long range interactions, must be > 7 ')
        else:
            if dims[0] < 7 or dims[1] < 7 or dims[2] < 7:
                raise KeyFileException('Box size is too small to correctly support super-long range interactions, must be > 7 ')


        ## Check out resized equilibrium variables and fix as needed
        ##
        if self.keyword_lookup['RESIZED_EQUILIBRATION']:
            if len(self.keyword_lookup['RESIZED_EQUILIBRATION']) != len(self.keyword_lookup['DIMENSIONS']):
                raise KeyFileException('Number of dimensions for compressed equilibration and final simulation are not the same')

            # resized_equilibration 
            for (real, eq) in zip(self.keyword_lookup['DIMENSIONS'], self.keyword_lookup['RESIZED_EQUILIBRATION']):
                if eq > real:
                    raise KeyFileException('Resized equilibration dimension is larger than final dimension - not yet supported (only smaller)')

            if self.keyword_lookup['EQUILIBRATION'] == 0:
                print("[WARNING]: using RESIZED_EQUILIBRATION without an equilibration period makes no sense. Deactivating RESIZED_EQUILIBRATION") 
                self.keyword_lookup['RESIZED_EQUILIBRATION'] = False

                
        ##
        ## if analysis code is provided check it can be loaded
        # first check if the file even exists 
        if not os.path.isfile(self.keyword_lookup['PARAMETER_FILE']):
                raise KeyFileException('Unable to find parameter file at location [%s]. Please verify the file exists.' %(self.keyword_lookup['PARAMETER_FILE']))

            # check if we can load as a python module...

      
        ##
        ## if analysis code is provided check it can be loaded
        if self.keyword_lookup['ANALYSIS_MODULE']:

            if not os.path.isfile(self.keyword_lookup['ANALYSIS_MODULE']):
                raise KeyFileException('Unable to find analysis module file. Passed filename is: %s. If this is a relative path Please verify the file exists.' %(self.keyword_lookup['ANALYSIS_MODULE']))
            
            # NOTE we overwrite the actual file name with a Python function
            tmp = self.keyword_lookup['ANALYSIS_MODULE']
            self.keyword_lookup['ANALYSIS_MODULE'] = file_utilities.custom_analysis_module_import(self.keyword_lookup['ANALYSIS_MODULE'])
            
            print("Loaded analysis code [%s] into [%s]" % (tmp, self.keyword_lookup['ANALYSIS_MODULE']))
            del tmp


        ##
        ## Check restart file and then read in
        if self.keyword_lookup['RESTART_FILE']:

            ## ----------------------------------------------------------------------------------------------------
            ## This block of code here is where 100% of the restart file sanity checking is going to happen

            
            # first see if we can even find the restart file...
            if not os.path.isfile(self.keyword_lookup['RESTART_FILE']):
                raise KeyFileException('Unable to find restart file. Passed filename is: %s. If this is a relative path Please verify the file exists.' %(self.keyword_lookup['RESTART_FILE']))
                
            # if yes initialize and then construct
            restart_obj = restart.RestartObject()
            restart_obj.build_from_file(self.keyword_lookup['RESTART_FILE'])
            
            print("Loaded restart information from: %s" % (self.keyword_lookup['RESTART_FILE']))
            self.keyword_lookup['RESTART_FILE'] = restart_obj

            # finally we ask is if any EXTRA_CHAIN were provided, and if yes add these
            if len(self.keyword_lookup['EXTRA_CHAIN']) > 0:
                self.__check_experimental_features('EXTRA_CHAIN')
                
                try:
                    # recall the self.keyword_lookup['EXTRA_CHAIN'] is a list where each element has two 
                    # elements, [0]= number of chains [1]  = chain sequence
                    for extra_chain in self.keyword_lookup['EXTRA_CHAIN']:
                        self.keyword_lookup['RESTART_FILE'].add_extra_chains(extra_chain)
                        
                except RestartException as e:
                    raise KeyFileException(f'\n\nError when parsing EXTRA_CHAIN line. Full error below: {e}')


                
                    

            # finally using the restart file sanity check input WRT the current keyfile to make sure everything
            # seems OK...
            self.sanity_check_and_update_with_restart_file()
            ## ----------------------------------------------------------------------------------------------------


        ## ---------------------------------------------------------        
        # Check for missuse of experimental features...

        for kw in CONFIG.EXPERIMENTAL_KEYWORDS:

            # if a keyword included in the EXPERIMENTAL_KEYWORDS list is NOT set to the default value            
            if self.keyword_lookup[kw] != self.DEFAULTS[kw]:
                self.__check_experimental_features(kw)

                ## ask if experimental features was set`
                #if self.keyword_lookup['EXPERIMENTAL_FEATURES'] == False:
                #    raise KeyFileException('\n\nExperimental or non-supported feature ({kw}) being proposed but EXPERIMENTAL_FEATURES is not set (or set to False).\n')
                    
                
                
    #-----------------------------------------------------------------
    #                    
    def set_defaults(self):
        """
        This function assigns default values from the self.DEFAULT dictionary to the required keywords 

        """

        # first check all required keywords were set
        for KW in self.required_keywords:
            if KW not in list(self.keyword_lookup.keys()):
                raise KeyFileException('ERROR: Keyfile does not define [%s]  ' % KW)

        if ('CHAIN' not in list(self.keyword_lookup.keys())) and ('RESTART_FILE' not in list(self.keyword_lookup.keys())):
            raise KeyFileException('No CHAIN keyword nor RESTART_FILE provided - no information on system provided!')

        # now cycle through ALL the keywords in the expected keyword dictionary and
        # assign using the defaults if they haven't already been asigned 
        for KW in self.expected_keywords:
            if KW not in list(self.keyword_lookup.keys()):                

                # ONE edge case, we want to print the random seed being used at all costs 
                
                if KW == "SEED": 
                    print("Using random seed [%s]" % (self.DEFAULTS[KW]))
                else:
                    print("No %s set - using default [%s]" % (KW, self.DEFAULTS[KW]))

                self.keyword_lookup[KW] = self.DEFAULTS[KW]



    #-----------------------------------------------------------------
    #    
    def set_dynamic_defaults(self):
        """

        This function allows final default values to update such that 'defaults' can respond to 
        passed parameters.

        """

        # if no analysis module was passed set custom analysis to more steps than we're running
        if self.keyword_lookup['ANALYSIS_MODULE'] is False:
            self.keyword_lookup['ANA_CUSTOM'] = self.keyword_lookup['N_STEPS'] + 10
                

        # if any of the analysis keywords have been set to 0 or -1 or anything less than 1 take this
        # to mean this analysis should not be performed, so set the frequency to a number that is larger than the
        # total number of steps
        for tmpkw in ['ANALYSIS_FREQ','ANA_POL', 'ANA_INTSCAL', 'ANA_DISTMAP', 'ANA_ACCEPTANCE', 'ANA_INTER_RESIDUE', 'ANA_CLUSTER', 'ANA_CUSTOM', 'ENERGY_CHECK']:
            if self.keyword_lookup[tmpkw] < 1:
                self.keyword_lookup[tmpkw] = self.keyword_lookup['N_STEPS'] + 10

        if self.keyword_lookup['RESTART_FREQ'] == "Every 10th-percentile":
            self.keyword_lookup['RESTART_FREQ'] = int(self.keyword_lookup['N_STEPS'] / 10) # deliberate effective floor being used here..

                
                                  
    #-----------------------------------------------------------------
    #    
    def print_summary(self):
        """
        Function that prints a full summary of the 

        """

        def section(msg):
            IO_utils.newline()
            IO_utils.horizontal_line(hzlen=40, linechar='*')
            print("--> %s"%(msg))
            IO_utils.newline()
            


        IO_utils.status_message("KEYFILE SUMMARY",'major')
        
        print("--> System Overview")
        print("Total number of steps     : %i" % self.keyword_lookup['N_STEPS'])
        print("Num. equilibration steps  : %i" % self.keyword_lookup['EQUILIBRATION'])
        print("Start temperature         : %3.2f" % self.keyword_lookup['TEMPERATURE'])
        print("Final temperature         : %3.2f" % self.keyword_lookup['EQUILIBRIUM_TEMPERATURE'])
        print("Lattice-to-Angstroms      : %5.2f" % self.keyword_lookup['LATTICE_TO_ANGSTROMS'])

        ## logging
        pimmslogger.log_status("NUMBER OF STEPS : %i" % self.keyword_lookup['N_STEPS'], timestamp=False)
        pimmslogger.log_status("EQUIL. STEPS    : %i" % self.keyword_lookup['EQUILIBRATION'], timestamp=False)
        pimmslogger.log_status("START TEMPERATURE     : %3.2f" % self.keyword_lookup['TEMPERATURE'], timestamp=False)
        pimmslogger.log_status("FINAL TEMPERATURE     : %3.2f" % self.keyword_lookup['EQUILIBRIUM_TEMPERATURE'], timestamp=False)
        

        # count total occupied volume fraction
        total=0
        chain_count = 0
        for chain in self.keyword_lookup['CHAIN']:            
            total = total + chain[0] * len(chain[1])
            chain_count = chain_count+ chain[0]
        
        
        ##
        ## BOX DIMENSIONS SECTION
        ##
        section('Box Dimensions')

        ## If we are running an initial resized simulation, provide info on conc. for the equilibrations
        if self.keyword_lookup['RESIZED_EQUILIBRATION']:

            print("Equilibration box dimensions:")
            print("")

            if len(self.keyword_lookup['RESIZED_EQUILIBRATION']) == 2:
                print("Initial equilibration (eq.) box dimensions  : %i x %i" % (self.keyword_lookup['RESIZED_EQUILIBRATION'][0],self.keyword_lookup['RESIZED_EQUILIBRATION'][1]))
            else:
                print("Initial equilibration (eq.) box dimensions  : %i x %i x %i" % (self.keyword_lookup['RESIZED_EQUILIBRATION'][0], self.keyword_lookup['RESIZED_EQUILIBRATION'][1], self.keyword_lookup['RESIZED_EQUILIBRATION'][2]))

            if len(self.keyword_lookup['DIMENSIONS']) == 2:
                print("Total occupied volume fraction during eq. = %3.5f" % (float(total)/(self.keyword_lookup['RESIZED_EQUILIBRATION'][0]*self.keyword_lookup['RESIZED_EQUILIBRATION'][1])))
            else:
                print("Total occupied volume fraction during eq. = %3.5f" % (float(total)/(self.keyword_lookup['RESIZED_EQUILIBRATION'][0]*self.keyword_lookup['RESIZED_EQUILIBRATION'][1]*self.keyword_lookup['RESIZED_EQUILIBRATION'][2])))

                # assume a conversion of 1 lattice unit = 4 angstroms - *0.4 is *4 / 10 to get in units of nm
                
                v_in_nm_3 = self.keyword_lookup['RESIZED_EQUILIBRATION'][0]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1 * self.keyword_lookup['RESIZED_EQUILIBRATION'][1]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1 * self.keyword_lookup['RESIZED_EQUILIBRATION'][2]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1
                v_in_L  = 1e-24*v_in_nm_3
                conc = (chain_count/6.02e23)/v_in_L        
                print("Total concentration of solute during equilibration  = %10.12f (M)" % (conc))

            # space between EQ and main simulation section
            print("")
            print("Main simulation box dimensions:")
            print("")



        # once we get here we're printing box information on the full simulation
        if len(self.keyword_lookup['DIMENSIONS']) == 2:
            print("BOX DIMENSIONS                 : %i x %i" % (self.keyword_lookup['DIMENSIONS'][0],self.keyword_lookup['DIMENSIONS'][1]))
        else:
            print("BOX DIMENSIONS                 : %i x %i x %i" % (self.keyword_lookup['DIMENSIONS'][0], self.keyword_lookup['DIMENSIONS'][1], self.keyword_lookup['DIMENSIONS'][2]))

        if len(self.keyword_lookup['DIMENSIONS']) == 2:
            print("Total occupied volume fraction : %3.5f" % (float(total)/(self.keyword_lookup['DIMENSIONS'][0]*self.keyword_lookup['DIMENSIONS'][1])))
        else:
            print("Total occupied volume fraction : %3.5f" % (float(total)/(self.keyword_lookup['DIMENSIONS'][0]*self.keyword_lookup['DIMENSIONS'][1]*self.keyword_lookup['DIMENSIONS'][2])))

            # assume a conversion of 1 lattice unit = 4 angstroms - *0.4 is *4 / 10 to get in units of nm
            v_in_nm_3 = self.keyword_lookup['DIMENSIONS'][0]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1*self.keyword_lookup['DIMENSIONS'][1]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1*self.keyword_lookup['DIMENSIONS'][2]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1

            v_in_L  = 1e-24*v_in_nm_3
            conc = (chain_count/6.02e23)/v_in_L      
            print("Total conc. of solute(s)       : %10.12f (M)" % (conc))
        


        ##
        ## QUENCH SETTINGS SECTION
        ##
        section('Quench Settings')

        if self.keyword_lookup['QUENCH_RUN']:

            quenchsteps = (self.keyword_lookup['QUENCH_FREQ']/self.keyword_lookup['QUENCH_STEPSIZE'])*abs(self.keyword_lookup['QUENCH_START'] - self.keyword_lookup['QUENCH_END'])

            if self.keyword_lookup['QUENCH_AS_EQUILIBRATION']:
                print("Quench running as equilibration: TRUE")
                print("Equilibration length: %i" % self.keyword_lookup['EQUILIBRATION'])
            else:                
                print("Quench running as equilibration: FALSE")
            print("Number of steps for quenching:  %i" % quenchsteps)

            # Quick comment on this: quenchsteps vs./ the equilibration quench are actually different. quenchsteps is the number of steps 
            # we spend quenching - once we've reached the target temperature we are by definition no longer quenching. HOWEVER, for 
            # equilibration we continue equilibrating for QUENCH_STEPSIZE steps before equilibration is finished (i.e. we need to run
            # *some* simulation at the production temperature for the equilibration to actually equilibrate at that temperature. 
            # Hence these two numbers differ by QUENCH_STEPSIZE steps. This is _not_ a bug!!

            print("QUENCH FREQ  : %i"     % self.keyword_lookup['QUENCH_FREQ'])
            print("QUENCH START : %3.2f " % self.keyword_lookup['QUENCH_START'])
            print("QUENCH STEP  : %3.2f " % self.keyword_lookup['QUENCH_STEPSIZE'])
            print("QUENCH END   : %3.2f " % self.keyword_lookup['QUENCH_END'])
        else:            
            print("NO TEMPERATURE QUENCH IN EFFECT")

        
        ##
        ## QUENCH SETTINGS SECTION
        ##
        section('Simulation Components')

        # if we're running a non-interacting simulation
        if self.keyword_lookup['NON_INTERACTING']:

            print("")
            print("NOTE: This is a non-interacting simulation - all energy interactions other")
            print("      than excluded volume and angle potentials have been switched OFF")
            print("")
            

        chainGroup = 1
        total=0
        for chain in self.keyword_lookup['CHAIN']:            
            print("Chain group %i:" % chainGroup)
            print("   %i copies of the following chain:" % chain[0])
            print("   %s" % chain[1])
            chainGroup = chainGroup+1
            total = total + chain[0] * len(chain[1])
            print("")
                                                                

        ##
        ## QUENCH SETTINGS SECTION
        ##
        section('Output Parameters')
        print("Print-to-screen freq.  : %i" % self.keyword_lookup['PRINT_FREQ'])
        print("XTC out freq.          : %i" % self.keyword_lookup['XTC_FREQ'])
        print("Energy out freq.       : %i" % self.keyword_lookup['EN_FREQ'])
        print("Overal analysis freq.  : %i" % self.keyword_lookup['ANALYSIS_FREQ'])
        IO_utils.newline(2)
        print("Keyfile fully parsed! Preparing to start the simulation...")
        IO_utils.horizontal_line()
        


    #-----------------------------------------------------------------
    #    
    def assign_default(self):
        """
        This is the function which defines the default values. The absolute reference
        default values are defined in


In reality this could 
        be in a configuration file somewhere, but in the interest of keeping everything 
        within this file we define those default values here.
        
        This basically builds up a self.DEFAULTS dictionary which defines default values 
        for each keyword. Note that for a few of these the default can depend on the 
        keywords parsed.

        Note that these defaults are used in two ways:

        1. If keywords are not provided then they are used as the default. Obviously.

        2. This can be used as a test to ask if a keyword has been set, because IF
           you pass a keyword and it changes the parsed keyword from the default this
           is the only time we care about a keyword being provided, hence functionaly
           this is how we define if a keyword is provided (or not). In particular, this
           is used for evaluating if EXPERIMENTAL keywords are being used.
        
        

        """

        # set defaults from the CONFIG file
        for k in CONFIG.DEFAULTS:
            self.DEFAULTS[k] = CONFIG.DEFAULTS[k]
        

        # if we defined a standard analysis frequency...
        if 'ANALYSIS_FREQ' in self.keyword_lookup:
            anafreq = self.keyword_lookup['ANALYSIS_FREQ']
        else:
            anafreq = self.DEFAULTS['ANALYSIS_FREQ']

        # assign the analysis frequencies to the default value
        self.DEFAULTS['ANA_POL']                = anafreq
        self.DEFAULTS['ANA_INTSCAL']            = anafreq
        self.DEFAULTS['ANA_DISTMAP']            = anafreq
        self.DEFAULTS['ANA_ACCEPTANCE']         = anafreq
        self.DEFAULTS['ANA_INTER_RESIDUE']      = anafreq
        self.DEFAULTS['ANA_CLUSTER']            = anafreq

        # get a real random integer
        random.seed()            
        self.DEFAULTS['SEED']    = random.randint(1,sys.maxsize-1)


    def add_derived_keywords(self):
        """
        Final function where additional derived keywords can be added

        """

        # we consider the end-to-end distance to be a polymeric property
        self.keyword_lookup['ANA_END_TO_END'] = self.keyword_lookup['ANA_POL'] 


        # set the equilibrium temperature - if we're doing a temperature run use the final temperature
        # from the run (note 'TEMPERATURE' will have already been updated to match QUENCH_START at this
        # point), if not just used the default temperature
        if self.keyword_lookup['QUENCH_RUN']:
            self.keyword_lookup['EQUILIBRIUM_TEMPERATURE'] = self.keyword_lookup['QUENCH_END']
        else:
            self.keyword_lookup['EQUILIBRIUM_TEMPERATURE'] = self.keyword_lookup['TEMPERATURE']
            
                                           
        
    def sanity_check_and_update_with_restart_file(self):    
        """
        Run sanity checks and update the CHAIN information so the box dimensions/concentration
        info is correct. Note that we DON'T sanity check the restart file input against the
        keyfile-defined chain values. This means we can be agnostic about what's in a restart
        file, and allows restarts to be more permissive. 

        This set of sanity checks is actually the last thing done, which means that we have 
        to re-check some of the things that were already checked with updated information
        read from the restart file.

        Below is the general ruberic for how restart files are dealt with:

        - The restart file fully overwrites all chain information. Chain information in a keyfile
          is completely ignored if a restart file is provided.

        - By default we assume the keyfile provided dimensions and hardwall status are to be used
          by in the restart. However, this may actually not be compatible, in which case an exception
          is thrown. The RESTART_OVERRIDE_DIMENSIONS, RESTART_OVERRIDE_HARDWALL force the simulation 
          to use the dimensions and hardwall values passed by the restart file. The one trick here is that
          if the keyfile requires an resized_equilibration then the restart file's dimensions will be
          used for the initial equilibration, and the production part of the simulation will be run
          using information from the keyfile.
        

        """


        if not self.keyword_lookup['RESTART_FILE']:
            return 

        new_chains ={}

        # grab the restart object
        restart_object = self.keyword_lookup['RESTART_FILE']
        
        print("")        
        print("############################")
        print("#   Reading Restart File   #")
        print("############################")
        print("")
        print("Restart file dimensions:    %s" % restart_object.dimensions)
        print("Restart file hardwall flag: %s"   % restart_object.hardwall)
        print("")
        

        ## ........................................................................................................................
        ##
        ## Part 1: Chains
        ## 

        # extract the chains from the restart object into a dictionary indexed
        # by type, where the tuple is  [ count, chain_sequence ]

        # Note that we do not worry about EXTRA_CHAINS here as they are built
        # de novo and need to be kept seperate from the CHAINS here. This entire
        # block of code just ensures that we can convert a RESTART file into data
        # that matches a <COUNT> <CHAIN SEQUENCE> format
                                           
                                    
        chain_type_dictionary ={}
        for chainID in restart_object.chains:
            
            # for each chainID get the chain sequence and 
            # the chain type
            c_seq  = restart_object.chains[chainID][1]
            c_type = restart_object.chains[chainID][2]

            # if this is the first time we encounter this chain type create a new 
            # entry in the chain_type_dictionary. 
            if c_type not in chain_type_dictionary:
                chain_type_dictionary[c_type] = [0,c_seq]

            # next verify that a chain of type c_type matches the other chains of 
            # type c_type that we've already seen so far
            ## NOTE - over places in the code would allow this; probably should address this at
            # some point for consistency...
            if not chain_type_dictionary[c_type][1] == c_seq:
                raise RestartException('When reading in the restart file, found a chain of type %i that did not match sequence of another chain of type %s. Chain sequences are\n: %s\n%s\n' % (c_type, c_type, chain_type_dictionary[c_type][1], c_seq))


            # if all seems good increment
            chain_type_dictionary[c_type][0] = chain_type_dictionary[c_type][0] + 1

        # finally sort the types and reconstruct a chains list, that has
        # each type as a seperate entry (i.e. [[count_1, seq_1],[count_2, seq_2]] 
        # and so on

        chain_types = list(chain_type_dictionary.keys())
        chain_types.sort()
        chains = []
        for CT in chain_types:
            chains.append(chain_type_dictionary[CT])

        self.keyword_lookup['CHAINS'] = chains
        
        if len(chain_types) > 1:
            print("--> Read in %i different chain types from the restart file" % (len(chain_types)))
        else:
            print("--> Read in a single chain type from the restart file") 
        
            
        print("--> Chain(s) read in from restart file are as follows:")
        for tmp in chains:
            if tmp[0] == 1:                
                print(f"    1 copy of {tmp[1]}")
            else:
                print(f"    {tmp[0]} copies of {tmp[1]}")

        print('')

        # note no need to actually check stuff, but, print things here...
        if len(restart_object.extra_chains) > 0:
            print("--> Also read in additiona; chains from EXTRA_CHAIN keyword")
            print("--> Chain(s) read in from EXTRA_CHAIN keyword are as follows:")
            for tmp in self.keyword_lookup['EXTRA_CHAIN']:
                if tmp[0] == 1:                
                    print(f"    1 copy of {tmp[1]}")
                else:
                    print(f"    {tmp[0]} copies of {tmp[1]}")
            
        print('')                  
        
        ## ........................................................................................................................
        ##
        ## Part 2: Hardwall rules
        ## 

        # if the keyfile says we should default to the hardwall rules associated with the restart file...
        if self.keyword_lookup['RESTART_OVERRIDE_HARDWALL']:
            print("Setting HARDWALL keyword based on restart file")
            self.keyword_lookup['HARDWALL'] = restart_object.hardwall

            # assume that the chains are consistent with the restart file mode (maybe should explicitly check this
            # in future versions...?)

        # if the restart object was a PBC simulation (i.e. not hardwall)
        elif not restart_object.hardwall:

            # if the restart file was a PBC simulation and we are trying to run a hardwall simulation - not compatible
            if self.keyword_lookup['HARDWALL']:
                raise RestartException('\n\nRestart file describes a periodic boundary simulation, but the keyfile is trying to run as a hardwall simulation. These options are incompatible with one another. Please set HARDWALL : False (or delete the HARDWALL keyword)\n')
                
            # if the restart file was a PBC simulation and we are trying to run a RESIZED_EQUILIBRIUM situation this doesn't work
            if self.keyword_lookup['RESIZED_EQUILIBRATION']:
                raise RestartException('\n\nRestart file describes a periodic boundary simulation, but the keyfile is trying to run as simulation that has a resized equilibration step (RESIZED_EQUILIBRATION : True). These options are incompatible with one another. The initial simulation box that is then resized MUST be a hardwall simulation.\n')
                
            
        # if the restart object was a hardwall simulation we can, from there, start either a PBC or a hardwall simulation. Don't need to do anything,
        # just wanted to explicitly include this case in the code to show it is considered!
        else:   
            if not self.keyword_lookup['HARDWALL']:
                print("While restart file was generated by a hardwall simulation,\nthis simulation applies periodic boundary conditions\n")
            



        ## ........................................................................................................................
        ##
        ## Part 3: Dimension Rules
        ##
        
        # if the restart file is being used to override everything
        if self.keyword_lookup['RESTART_OVERRIDE_DIMENSIONS']:
            print("Setting DIMENSIONS keyword based on restart file")
            self.keyword_lookup['DIMENSIONS'] = restart_object.dimensions
            
            if self.keyword_lookup['RESIZED_EQUILIBRATION']:
                raise RestartException("\n\nRESTART_OVERRIDE_DIMENSIONS is set to true, but the simulation also wants a resized equilibration. These options are incompatible. If the simulation being restarted is a hardwall simulation then you can use the following approach\n1) Set the RESIZED_EQUILIBRATION to equal (or bigger) than the restart file's dimensions\n2) Set the DIMENSIONS to the production dimensions desired\n3) Set RESTART_OVERRIDE_DIMENSIONS and RESTART_OVERRIDE_HARDWALL to False\n")

        else:
            print("Setting DIMENSIONS based on ...")
            dimchange_flag=False

            n_dims = len(self.keyword_lookup['DIMENSIONS'])

            # check if number of dimensions in keyfile matches number of dimensions in restart file - if not throw an error
            if n_dims != len(restart_object.dimensions):
                raise RestartException("\n\nRestart object has %i dimensions but kefile specifies %i dimensions\n" % (len(restart_object.dimensions), n_dims) )

            # check that if the restart file was a PBC simulation, the new dimensions match (note we KNOW that if we're doing a PBC we're not running a resize equilibration
            # simulation because this is dealt with in part 2)
            if not restart_object.hardwall:
                for dim in range(0, n_dims):
                    if not restart_object.dimensions[dim] == self.keyword_lookup['DIMENSIONS'][dim]:
                        raise RestartException("\n\nThe restart file was a PBC simulation, so the box dimensions must match EXACTLY to avoid breaking PBC assumptions. The restart file dimensions are [%s] but the DIMENSIONS keyword in the keyfile is set to [%s]\n" %(restart_object.dimensions, self.keyword_lookup['DIMENSIONS']))

                        
            ### If the initial part of the simulation will be running on the RESIZED_EQUILIBRATION dimensions..
            # if we're trying to use a restart file to run a RESIZED_EQUILIBRATION (note we now know that the restart file was a hardwall simulation
            # as this is checked in part 2)
            if self.keyword_lookup['RESIZED_EQUILIBRATION']:                

                # must ensure that the resize file dimenions are equal to or smaller than the resized equilibration being used
                for dim in range(0, n_dims):
                    if self.keyword_lookup['RESIZED_EQUILIBRATION']:
                        if restart_object.dimensions[dim] > self.keyword_lookup['RESIZED_EQUILIBRATION'][dim]:
                            raise RestartException("\n\nRESIZED_EQUILIBRATION dimensions (%s) are smaller than the restart file's dimensions. This is an incompatible situation.\n"%(self.keyword_lookup['RESIZED_EQUILIBRATION'][dim], restart_object.dimensions[dim]))

                # finally center lattice if needed
                for dim in range(0, n_dims):
                    if not (restart_object.dimensions[dim] == self.keyword_lookup['RESIZED_EQUILIBRATION'][dim]):
                        dimchange_flag=True

                if dimchange_flag:
                    print("Dimensions in restart file were %s, but these have been updated to %s based on the keyfile (where RESIZED_EQUILIBRATION is set to %s)" % (restart_object.dimensions, self.keyword_lookup['RESIZED_EQUILIBRATION'],self.keyword_lookup['RESIZED_EQUILIBRATION']))
                    restart_object.update_lattice_dimensions(self.keyword_lookup['RESIZED_EQUILIBRATION'])
                else:
                    print("Dimensions in restart file matched RESIZED_EQUILIBRATION keyword: %s" % (self.keyword_lookup['RESIZED_EQUILIBRATION']))

            ### If the intial part of the simulation will be running on the normal DIMENSIONS keyword
            # finally, ensure that the dimensions are equal to or bigger than the restart file (this is a catch all but is true for
            # PBC and hardwall simulations alike)
            else:
                for dim in range(0, n_dims):
                    if restart_object.dimensions[dim] > self.keyword_lookup['DIMENSIONS'][dim]:
                        raise RestartException("\n\nDIMENSIONS dimensions [%s] are smaller than the restart file's dimensions [%s]. This is an incompatible situation.\n"%(self.keyword_lookup['DIMENSIONS'][dim], restart_object.dimensions[dim]))

                    if not (restart_object.dimensions[dim] == self.keyword_lookup['DIMENSIONS'][dim]):
                        dimchange_flag=True

                if dimchange_flag:
                    print("Dimensions in restart file were %s, but these have been updated to %s based on the keyfile" % (restart_object.dimensions, self.keyword_lookup['DIMENSIONS']))
                    restart_object.update_lattice_dimensions(self.keyword_lookup['DIMENSIONS'])
                else:
                    print("Dimensions in restart file matched keyfile: %s" % (self.keyword_lookup['DIMENSIONS']))


        print("")
        print(".... Restart file processed")
        print("############################")
        print("")

        
