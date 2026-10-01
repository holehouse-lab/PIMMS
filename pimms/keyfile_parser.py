## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
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


import math
import random
import os
import sys
import os.path
import numpy as np

from . import IO_utils
from . import latticeExceptions
from .latticeExceptions import KeyFileException, RestartException
from . import file_utilities
from . import restart
from . import pimmslogger
from . import CONFIG
from . data_structures import FreezeFile


def print_keyword_info():
    """
    Print a human-readable table of every supported keyfile keyword.

    Iterates over ``CONFIG.KEYWORDS_DESCRIPTION`` and prints, for each keyword,
    its name, its expected value type, and a short description. Output is written
    to standard output in aligned columns.

    Returns
    -------
    None
        No return value; the keyword table is printed to standard output.
    """
    maxlen = 25

    for d in CONFIG.KEYWORDS_DESCRIPTION:

        spacer = " "*(maxlen - len(d))

        # each value is a [type, description] pair
        keyword_type, keyword_description = CONFIG.KEYWORDS_DESCRIPTION[d]
        print("%s%s | %s | %s" % (d, spacer, keyword_type, keyword_description))


# ===================================================================================================
#
#

def _quench_rung_count(temperature_range, stepsize):
    """
    Number of temperature updates a quench needs to traverse ``temperature_range``
    in steps of ``stepsize``.

    The runtime (``nonequilibrium_utils.update_temperature_in_quench``) snaps to
    ``QUENCH_END`` once it is within a relative ``1e-9`` of it, so a ratio that is
    an integer up to floating-point round-off (``(1.0 - 0.7) / 0.1`` evaluates to
    ``3.0000000000000004``) counts as that integer rather than being rounded up
    to one rung too many; anything else rounds up.

    Parameters
    ----------
    temperature_range : float
        ``abs(QUENCH_START - QUENCH_END)``.

    stepsize : float
        Magnitude of ``QUENCH_STEPSIZE``.

    Returns
    -------
    int
        Number of temperature changes the quench will make.
    """
    ratio = float(temperature_range) / float(abs(stepsize))
    nearest = round(ratio)
    if abs(ratio - nearest) <= 1e-9 * max(1.0, abs(ratio)):
        return int(nearest)
    return int(np.ceil(ratio))



def write_keyword_lookup(keyword_lookup, output_filename, PADDING=10, header_lines=(),
                         expected_keywords=None, keywords_with_multiple_entries=None):
    """
    Write a keyword dictionary to disk as a keyfile the parser can read back.

    This is the writer behind ``KeyFileParser.write_keyfile``, split out so that
    anything holding a resolved keyword dictionary (the Simulation, which is
    handed only the dictionary) can write one without a parser. Values are
    serialised in the syntax the parser reads: multi-entry keywords one line per
    entry, lists space-separated, a disabled frequency as ``0``, a loaded freeze
    file by its path. Derived ``__`` keywords, ``UNSET`` sentinels, the obsolete
    ``CRANKSHAFT_MODE`` and loaded objects with no keyfile form (a restart
    object, an imported analysis callable) are skipped. After restart processing
    the ``CHAIN`` list already contains the ``EXTRA_CHAIN`` chains and the restart
    object itself cannot be written, so a restart run is written as the
    equivalent de novo composition: the same chains in the same box, started
    afresh rather than from the snapshot.

    Parameters
    ----------
    keyword_lookup : dict
        The resolved keyword dictionary.

    output_filename : str
        Where to write.

    PADDING : int, optional
        Column the keyword names are padded out to. Default 10.

    header_lines : sequence of str, optional
        Lines written at the top of the file as ``#`` comments, for provenance.
        Default none.

    expected_keywords : sequence of str, optional
        The keywords that are keyfile syntax; anything else is skipped. Default
        ``CONFIG.EXPECTED_KEYWORDS``.

    keywords_with_multiple_entries : sequence of str, optional
        The keywords written one line per entry. Default the parser's list
        (``CHAIN``, ``EXTRA_CHAIN``, ``ANA_RESIDUE_PAIRS``).

    Returns
    -------
    None
        The file is written.
    """
    if expected_keywords is None:
        expected_keywords = CONFIG.EXPECTED_KEYWORDS
    if keywords_with_multiple_entries is None:
        keywords_with_multiple_entries = ['CHAIN', 'EXTRA_CHAIN', 'ANA_RESIDUE_PAIRS']
    with open(output_filename, 'w') as fh:
        for line in header_lines:
            fh.write('# ' + line + '\n' if line else '#\n')
        if header_lines:
            fh.write('\n')
    # open the file 

        # for each pair in the keyword_lookup dictionary
        for key, value in keyword_lookup.items():

            # define the padding to be used between the keyword and the value
            padding = ' ' * max((PADDING - len(key), 1))

            # serialise values in the same syntax the parser reads, so a
            # written keyfile round-trips: multi-entry keywords (CHAIN /
            # EXTRA_CHAIN / ANA_RESIDUE_PAIRS) get one line per entry,
            # int/float lists become space-separated, and post-parse objects
            # (loaded restart/freeze files, the imported analysis callable)
            # and derived __ keywords are skipped - they are not keyfile
            # syntax.

            # skip derived keywords (those starting with '__') 
            if key.startswith('__'):
                continue


            
            # skip any keyword that is not in the expected list
            if key not in expected_keywords:
                continue

            # obsolete keyword (accepted with a warning, ignored): never
            # write it back; this is maintained for backwards compatibility 
            # but we'll probrably kill at some point
            if key == 'CRANKSHAFT_MODE':
                continue

            # non-expressible sentinels: 'UNSET' placeholders, and False for
            # keywords that are NOT genuinely boolean (an unset feature) -
            # writing them verbatim made every full-init round-trip unparseable
            if isinstance(value, str) and value == 'UNSET':
                continue


            # (QUENCH_AS_EQUILIBRATION is a genuine boolean and is NOT in this
            # list: skipping its False used to leave a QUENCH_RUN keyfile without
            # a required quench keyword, so the written file did not re-parse)
            if value is False and key in ('RESIZED_EQUILIBRATION', 'EQUILIBRATION_OFFSET',
                                          'TSMMC_FIXED_OFFSET', 'ANALYSIS_MODULE',
                                          'RESTART_FILE', 'FREEZE_FILE'):
                continue


            # after restart processing CHAIN already holds the EXTRA_CHAIN
            # chains and RESTART_FILE is a loaded object (not writable), so
            # writing EXTRA_CHAIN too would re-parse as "EXTRA_CHAIN without a
            # restart file"; the CHAIN lines describe the same composition de novo
            if key == 'EXTRA_CHAIN' and keyword_lookup.get('RESTART_FILE') and \
                    not isinstance(keyword_lookup['RESTART_FILE'], str):
                continue

            # for the same reason the flags that only mean something alongside a
            # RESTART_FILE are written as False once that file is a loaded object
            # (and so not written): the overrides have already been applied to
            # DIMENSIONS / HARDWALL, and RESTART_CONTINUE : True without a
            # RESTART_FILE is refused on re-parse
            if key in ('RESTART_CONTINUE', 'RESTART_OVERRIDE_DIMENSIONS', 'RESTART_OVERRIDE_HARDWALL') and \
                    keyword_lookup.get('RESTART_FILE') and \
                    not isinstance(keyword_lookup['RESTART_FILE'], str):
                value = False
            if key in keywords_with_multiple_entries:
                if not value:
                    continue
                for entry in value:
                    fh.write(f'  {key} : {padding} {" ".join(str(x) for x in entry)}  \n')
                continue

            # a disabled analysis/energy-check frequency is stored internally as
            # the N_STEPS+10 display sentinel; written verbatim it would re-parse
            # as an ENABLED analysis at that cadence. 0 is the keyfile spelling
            # of "disabled".
            if key in keyword_lookup.get('__DISABLED_FREQUENCIES', ()):
                value = 0
                
            if isinstance(value, (list, tuple)):
                value = ' '.join(str(x) for x in value)
            elif key == 'FREEZE_FILE' and getattr(value, 'filename', None):
                # the loaded FreezeFile remembers its path: a round-tripped
                # keyfile used to drop it and run a fully mobile simulation
                value = value.filename
            elif key == 'ANALYSIS_MODULE' and callable(value) and \
                    keyword_lookup.get('__ANALYSIS_MODULE_PATH'):
                # likewise the imported analysis function: dropping it left
                # ANA_CUSTOM set with no module, so a re-run skipped the analysis
                value = keyword_lookup['__ANALYSIS_MODULE_PATH']
            elif not isinstance(value, (str, int, float, bool)):
                # loaded objects (RestartObject, callables, ...) cannot be
                # expressed in keyfile syntax - skip them
                continue
            fh.write(f'  {key} : {padding} {value}  \n')


#-----------------------------------------------------------------
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

        parse_only : bool, optional
            Flag which, if set to True, means the keyfile is read in but nothing more is done (i.e. no defaults set, no
            sanitization performed, no logging initialized). This is useful when the KeyFileParser is used to read in
            pre-run keyfiles where rigorous assessment is not needed. Default is False.


        """

        # expected keywords contains a list of possible keywords which can be read form the keyfile. Importantly if these
        # keywords are *not* included then they are set to their default values
        self.expected_keywords = CONFIG.EXPECTED_KEYWORDS

        # required keywords are those which MUST be included in the keyfile - i.e. PIMMS can't set default values for these
        # keywords. Note that required_keywords is a subset of expected_keywords
        self.required_keywords = CONFIG.REQUIRED_KEYWORDS 

        # list of keywords that can support multiple entries in a keyfile
        self.keywords_with_multiple_entries = ['CHAIN', 'EXTRA_CHAIN', 'ANA_RESIDUE_PAIRS']
        # whether SEED was in the keyfile (as opposed to generated as a default);
        # RESTART_CONTINUE refuses an explicit seed, since the generator state
        # then comes from the restart file and a seed would be a contradiction
        self._seed_was_given = False
        self.keyword_lookup = {}
        # every keyword encountered during parsing is recorded here (independent
        # of whether its handler ends up writing to keyword_lookup) so that
        # duplicate keywords are ALWAYS detected - see parse().
        self._seen_keywords = set()
        self.DEFAULTS = {}

        ## IF PARSE ONLY mode
        if parse_only:
            self.parse(filename)
            # provenance for keyfile_used.kf: where this came from, and whether the
            # seed was the user's or generated (derived keys are never written back)
            self.keyword_lookup['__KEYFILE'] = os.path.abspath(filename)
            self.keyword_lookup['__SEED_GIVEN'] = self._seed_was_given

            return 

        # print intro for the keyfile parsing
        IO_utils.horizontal_line(hzlen=40, linechar='*')
        IO_utils.status_message("Parsing keyfile [%s]" %(filename),'startup')
        IO_utils.status_message("Default values set are explicitly announced below:", 'startup')

        print("")

        # parse the keyfile (i.e. read in and deal with the file)
        self.parse(filename)
        # provenance for keyfile_used.kf: where this came from, and whether the
        # seed was the user's or generated (derived keys are never written back)
        self.keyword_lookup['__KEYFILE'] = os.path.abspath(filename)
        self.keyword_lookup['__SEED_GIVEN'] = self._seed_was_given


        # assigns default values to the internal system (though not YET to this keyfile)
        self.assign_default()      

        # any values missing from the keyfile get set to the default values
        self.set_defaults()         

        # update any values that depend on keyfile-derived parameters
        self.set_dynamic_defaults() 

        # run some sanity checks  
        self.run_sanity_checks()    

        # some keywords are not explicitly included in the keyfile but are derived from
        # the keyfile words. These are set here
        self.add_derived_keywords() 
                                    
        IO_utils.horizontal_line(hzlen=40, linechar='*')

        # initialize logging...
        pimmslogger.initialize()



    #-----------------------------------------------------------------
    #    
    def __repr__(self):
        """
        Return the developer/`repr` representation of the parser.

        This simply defers to :meth:`__str__`, so the representation is the same
        human-readable summary of all parsed keyword/value pairs.

        Returns
        -------
        str
            A formatted string listing each keyword and its current value.
        """
        return str(self)


    #-----------------------------------------------------------------
    #
    def __str__(self):
        """
        Return a human-readable summary of all parsed keyword/value pairs.

        Builds a multi-line string in which each line shows one entry from the
        internal ``keyword_lookup`` dictionary as ``keyword => value``.

        Returns
        -------
        str
            A formatted, multi-line string describing the parsed keyfile state.
        """
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
            Name of the keyword to be checked (i.e. a keyword where we require
            EXPERIMENTAL_FEATURES to be set to True for this parameter to be
            useable). Used only for the error message.

        Returns
        ----------
        None
            No return type; the function is silent if EXPERIMENTAL_FEATURES is True.

        Raises
        ------
        KeyFileException
            If EXPERIMENTAL_FEATURES is False.

        """
        if self.keyword_lookup['EXPERIMENTAL_FEATURES'] == False:
            raise KeyFileException(f'\n\nExperimental or non-supported keyword [{kw}] being proposed but EXPERIMENTAL_FEATURES is False.\n')

                       
    #-----------------------------------------------------------------
    #
    # Typed-value helpers. These convert a raw keyfile string into the expected
    # type and raise a clear, user-facing KeyFileException on malformed input
    # (instead of an opaque ValueError traceback). Using them everywhere is part
    # of sanity-checking every keyword input at startup.
    #
    def _kw_int(self, keyword, value):
        """
        Convert a raw keyfile string into an integer.

        Parameters
        ----------
        keyword : str
            Name of the keyword being parsed (used only for error messaging).
        value : str
            The raw value associated with the keyword in the keyfile.

        Returns
        -------
        int
            The value cast to an integer.

        Raises
        ------
        KeyFileException
            If ``value`` cannot be interpreted as an integer.
        """
        try:
            return int(value)
        except (ValueError, TypeError):
            raise KeyFileException(latticeExceptions.message_preprocess(
                "Keyword [%s] expects an integer value, but got [%s]" % (keyword, value)))

    def _kw_float(self, keyword, value):
        """
        Convert a raw keyfile string into a floating-point number.

        Parameters
        ----------
        keyword : str
            Name of the keyword being parsed (used only for error messaging).
        value : str
            The raw value associated with the keyword in the keyfile.

        Returns
        -------
        float
            The value cast to a float.

        Raises
        ------
        KeyFileException
            If ``value`` cannot be interpreted as a numeric value.
        """
        try:
            parsed_value = float(value)
        except (ValueError, TypeError):
            raise KeyFileException(latticeExceptions.message_preprocess(
                "Keyword [%s] expects a float-parsible value, but got [%s]" % (keyword, value)))

        # float() deliberately accepts spellings such as ``nan`` and ``inf``.
        # Those values defeat essentially every range check below (comparisons
        # with NaN are false) and can therefore reach the Monte Carlo kernels.
        if not math.isfinite(parsed_value):
            raise KeyFileException(latticeExceptions.message_preprocess(
                "Keyword [%s] expects a finite float value, but got [%s]" %
                (keyword, value)))
        return parsed_value

    def _kw_bool(self, keyword, value):
        """
        Convert a raw keyfile string into a boolean.

        The comparison is case-insensitive and whitespace-insensitive: the
        strings ``TRUE`` and ``FALSE`` (in any case) map to ``True`` and
        ``False`` respectively.

        Parameters
        ----------
        keyword : str
            Name of the keyword being parsed (used only for error messaging).
        value : str
            The raw value associated with the keyword in the keyfile.

        Returns
        -------
        bool
            ``True`` if ``value`` is ``TRUE``, ``False`` if it is ``FALSE``.

        Raises
        ------
        KeyFileException
            If ``value`` is neither ``TRUE`` nor ``FALSE``.
        """
        v = str(value).strip().upper()
        if v == 'TRUE':
            return True
        if v == 'FALSE':
            return False
        raise KeyFileException(latticeExceptions.message_preprocess(
            "Keyword [%s] expects TRUE or FALSE, but got [%s]" % (keyword, value)))

    def _kw_path(self, keyword, value):
        """
        Convert a raw keyfile string into a (non-empty) file path.

        The consumers of the path keywords (``RESTART_FILE``, ``FREEZE_FILE``,
        ``ANALYSIS_MODULE``, ``PARAMETER_FILE``) test the stored value for truth
        to decide whether the feature is active, so an empty value used to be
        indistinguishable from the keyword being absent: ``RESTART_FILE :`` with
        nothing after the colon silently ran a fresh simulation, ``FREEZE_FILE :``
        silently froze nothing. A path keyword that is present must name a file.

        Parameters
        ----------
        keyword : str
            Name of the keyword being parsed (used only for error messaging).
        value : str
            The raw value associated with the keyword in the keyfile.

        Returns
        -------
        str
            The stripped path, with a leading '~' expanded to the user's home
            directory.

        Raises
        ------
        KeyFileException
            If the value is empty.
        """
        v = str(value).strip()
        if not v:
            raise KeyFileException(latticeExceptions.message_preprocess(
                "Keyword [%s] expects a file path, but no value was given - remove the keyword "
                "if the feature is not wanted" % keyword))
        # '~' is expanded for every path keyword (it used to be expanded for
        # ANALYSIS_MODULE only, so PARAMETER_FILE : ~/x.prm was 'not found')
        return os.path.expanduser(v)

    def _kw_int_list(self, keyword, value):
        """
        Convert a space-separated keyfile string into a list of integers.

        Parameters
        ----------
        keyword : str
            Name of the keyword being parsed (used only for error messaging).
        value : str
            The raw value associated with the keyword; expected to be a
            whitespace-separated sequence of integers (e.g. box dimensions).

        Returns
        -------
        list of int
            The parsed integers, in the order they appear in ``value``.

        Raises
        ------
        KeyFileException
            If any whitespace-separated token cannot be cast to an integer.
        """
        try:
            parsed = [int(i) for i in value.split()]
        except (ValueError, TypeError):
            raise KeyFileException(latticeExceptions.message_preprocess(
                "Keyword [%s] expects a space-separated list of integers, but got [%s]" % (keyword, value)))
        if not parsed:
            # an empty value would silently become [] - which is falsy, so every
            # downstream check AND the feature itself would be skipped without a
            # word (a truncated RESIZED_EQUILIBRATION line silently disabled it)
            raise KeyFileException(latticeExceptions.message_preprocess(
                "Keyword [%s] was given no values - expected a space-separated list of integers" % (keyword,)))
        return parsed

    #-----------------------------------------------------------------
    #
    def parse(self, filename):
        """
        Main function reads in a keyfile and extracts out the relevant details 
        based on the keywords. Keywords are assigned to the self.keywords_lookup.
        
        Keywords must be defined as: 

        KEYWORD : VALUE # comment 

        This function reads through all the lines in the keyfile and in a 
        keyword-specific manner parses each keyword and assigns it to the 
        self.keyword_lookup dictionary. Any keywords that are expected but 
        NOT defined in the keyfile are then assigned as default values later.

        Parameters
        --------------
        filename : str
            Name of the keyfile to be read

        Returns
        ---------
        None
            No return type but assigns values to the self.keyword_lookup 
            dictionary, which itself is a set of key-value pairs of simulation 
            keywords            

        Raises
        ---------
        KeyFileException
            If a line has no keyword separator, if an unsupported keyword 
            is found, if a supported keyword has no handler (a bug), if 
            a keyword that does not support multiple entries appears 
            twice, or if any keyword's value cannot be parsed into the 
            expected type.
        """
        
        # parse line-by-line to avoid loading large keyfiles into memory.
        with open(filename, 'r') as fh:

            # for each line...
            for line in fh:


                # if it's a comment line skip the whole line
                if file_utilities.is_comment_line(line):
                    continue

                # remove comment section (at end of line)
                un_comment = file_utilities.remove_comments(line)

                # split once so values can legitimately contain ':' (e.g. paths).
                splitline = un_comment.split(':', 1)

                # if we didn't find a keyword separator
                if len(splitline) < 2:
                    raise KeyFileException(latticeExceptions.message_preprocess('On line [%s] - no keyword separator found...' % line))

                # if get here must have keyword/keyvalue
                putative_keyword = splitline[0].strip().upper()
                putative_value = splitline[1].strip()

                # if we recognize ye olde keyword...
                if putative_keyword in self.expected_keywords:

                    ## ** CHECK TO ENSURE WE DON'T OVERWRITE KEYWORDS **
                    # check if we've seen this keyword before - if we're trying to overwrite raise an exception.
                    # NOTE: we track this in self._seen_keywords (NOT keyword_lookup) so that a duplicate is
                    # always caught even for keywords whose handler may not write to keyword_lookup.
                    if putative_keyword in self._seen_keywords:

                        if putative_keyword in self.keywords_with_multiple_entries:
                            # this is OK - we can have multiple chains! 
                            pass
                        else:
                            raise KeyFileException(latticeExceptions.message_preprocess('Found a second occurence of the [%s] keyword. Please correct your keyfile and retry' % putative_keyword))

                    self._seen_keywords.add(putative_keyword)

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # CHAIN KEYWORD
                    if putative_keyword == "CHAIN":
                        chainSplit = putative_value.split()
                        if len(chainSplit) != 2:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid CHAIN keyword format [{putative_value}]. Expected: CHAIN : <count> <sequence>'))
                        try:
                            number_of_chains = int(chainSplit[0])
                        except ValueError:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid CHAIN keyword format [{putative_value}]. Expected integer chain count'))
                        if number_of_chains < 1:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid CHAIN keyword [{putative_value}]. Chain count must be >= 1'))
                        chain_sequence = chainSplit[1].strip()

                        if putative_keyword in self.keyword_lookup:
                            self.keyword_lookup['CHAIN'].append([number_of_chains, chain_sequence])
                        else:
                            self.keyword_lookup['CHAIN'] = [[number_of_chains, chain_sequence]]

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    # EXTRA_CHAIN KEYWORD
                    elif putative_keyword == "EXTRA_CHAIN":
                        chainSplit = putative_value.split()
                        if len(chainSplit) != 2:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid EXTRA_CHAIN keyword format [{putative_value}]. Expected: EXTRA_CHAIN : <count> <sequence>'))
                        try:
                            number_of_chains = int(chainSplit[0])
                        except ValueError:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid EXTRA_CHAIN keyword format [{putative_value}]. Expected integer chain count'))
                        if number_of_chains < 1:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid EXTRA_CHAIN keyword [{putative_value}]. Chain count must be >= 1'))
                        chain_sequence = chainSplit[1].strip()

                        if putative_keyword in self.keyword_lookup:
                            self.keyword_lookup['EXTRA_CHAIN'].append([number_of_chains, chain_sequence])
                        else:
                            self.keyword_lookup['EXTRA_CHAIN'] = [[number_of_chains, chain_sequence]]

                    # DIMENSIONS KEYWORD (define if simulation is 2D or 3D)
                    elif putative_keyword == 'DIMENSIONS':
                        self.keyword_lookup['DIMENSIONS'] = self._kw_int_list(putative_keyword, putative_value)
                        if not len(self.keyword_lookup['DIMENSIONS']) == 2 and not len(self.keyword_lookup['DIMENSIONS']) == 3:
                            raise KeyFileException(latticeExceptions.message_preprocess('Unexpected number of dimensions [%s] ' % line))

                    # conversion factor for PDB file writing
                    elif putative_keyword == 'LATTICE_TO_ANGSTROMS':
                        self.keyword_lookup['LATTICE_TO_ANGSTROMS'] = self._kw_float(putative_keyword, putative_value)
                        if not self.keyword_lookup['LATTICE_TO_ANGSTROMS'] > 0:
                            raise KeyFileException('LATTICE_TO_ANGSTROMS must be a positive length (got %s)' % putative_value)

                    # Dimensions of compressed equilibration box
                    elif putative_keyword == 'RESIZED_EQUILIBRATION':
                        self.keyword_lookup['RESIZED_EQUILIBRATION'] = self._kw_int_list(putative_keyword, putative_value)

                    # Dimensions of manual offset
                    elif putative_keyword == 'EQUILIBRATION_OFFSET':
                        self.keyword_lookup['EQUILIBRATION_OFFSET'] = self._kw_int_list(putative_keyword, putative_value)

                    # CASE_INSENSITIVE_CHAINS
                    elif putative_keyword == "CASE_INSENSITIVE_CHAINS":
                        self.keyword_lookup['CASE_INSENSITIVE_CHAINS'] = self._kw_bool(putative_keyword, putative_value)

                    # HARDWALL
                    elif putative_keyword == "HARDWALL":
                        self.keyword_lookup['HARDWALL'] = self._kw_bool(putative_keyword, putative_value)

                    # AUTOCENTER
                    elif putative_keyword == "AUTOCENTER":
                        self.keyword_lookup['AUTOCENTER'] = self._kw_bool(putative_keyword, putative_value)

                    # EXPERIMENTAL_FEATURES
                    elif putative_keyword == "EXPERIMENTAL_FEATURES":
                        self.keyword_lookup['EXPERIMENTAL_FEATURES'] = self._kw_bool(putative_keyword, putative_value)

                    # TEMPERATURE
                    elif putative_keyword == "TEMPERATURE":
                        self.keyword_lookup['TEMPERATURE'] = self._kw_float(putative_keyword, putative_value)

                    # N_STEPS
                    elif putative_keyword == "N_STEPS":
                        self.keyword_lookup['N_STEPS'] = self._kw_int(putative_keyword, putative_value)

                    # equilibration
                    elif putative_keyword == "EQUILIBRATION":
                        self.keyword_lookup['EQUILIBRATION'] = self._kw_int(putative_keyword, putative_value)

                    # energy parameter
                    elif putative_keyword == "PARAMETER_FILE":
                        self.keyword_lookup['PARAMETER_FILE'] = self._kw_path(putative_keyword, putative_value)

                    # PRINT_FREQUENCY
                    elif putative_keyword == "PRINT_FREQ":
                        self.keyword_lookup['PRINT_FREQ'] = self._kw_int(putative_keyword, putative_value)

                    # REDUCED_PRINTING
                    elif putative_keyword == "REDUCED_PRINTING":
                        self.keyword_lookup['REDUCED_PRINTING'] = self._kw_bool(putative_keyword, putative_value)

                    # XTC_FREQUENCY
                    elif putative_keyword == "XTC_FREQ":
                        self.keyword_lookup['XTC_FREQ'] = self._kw_int(putative_keyword, putative_value)

                    # EN_FREQUENCY
                    elif putative_keyword == "EN_FREQ":
                        self.keyword_lookup['EN_FREQ'] = self._kw_int(putative_keyword, putative_value)

                    # SEED
                    elif putative_keyword == "SEED":
                        self.keyword_lookup['SEED'] = self._kw_int(putative_keyword, putative_value)
                        self._seed_was_given = True

                    # ENERGY_CHECK
                    elif putative_keyword == "ENERGY_CHECK":
                        self.keyword_lookup['ENERGY_CHECK'] = self._kw_int(putative_keyword, putative_value)

                    # RESTART_FREQ
                    elif putative_keyword == "RESTART_FREQ":
                        self.keyword_lookup['RESTART_FREQ'] = self._kw_int(putative_keyword, putative_value)

                    # RESTART FILE
                    elif putative_keyword == "RESTART_FILE":
                        self.keyword_lookup['RESTART_FILE'] = self._kw_path(putative_keyword, putative_value)

                    # RESTART OVERRIDE DIMENSIONS
                    elif putative_keyword == "RESTART_OVERRIDE_DIMENSIONS":
                        self.keyword_lookup["RESTART_OVERRIDE_DIMENSIONS"] = self._kw_bool(putative_keyword, putative_value)

                    # RESTART OVERRIDE HARDWALL
                    elif putative_keyword == "RESTART_OVERRIDE_HARDWALL":
                        self.keyword_lookup["RESTART_OVERRIDE_HARDWALL"] = self._kw_bool(putative_keyword, putative_value)

                    # RESTART_CONTINUE (resume the run the restart file was written by)
                    elif putative_keyword == "RESTART_CONTINUE":
                        self.keyword_lookup["RESTART_CONTINUE"] = self._kw_bool(putative_keyword, putative_value)

                    # FREEZE_FILE
                    elif putative_keyword == 'FREEZE_FILE':
                        self.keyword_lookup['FREEZE_FILE'] = self._kw_path(putative_keyword, putative_value)

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## QUENCHING keywords
                    elif putative_keyword == "QUENCH_RUN":
                        self.keyword_lookup['QUENCH_RUN'] = self._kw_bool(putative_keyword, putative_value)

                    elif putative_keyword == "QUENCH_AS_EQUILIBRATION":
                        self.keyword_lookup['QUENCH_AS_EQUILIBRATION'] = self._kw_bool(putative_keyword, putative_value)

                    elif putative_keyword == "QUENCH_FREQ":
                        self.keyword_lookup['QUENCH_FREQ'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == "QUENCH_STEPSIZE":
                        self.keyword_lookup['QUENCH_STEPSIZE'] = abs(self._kw_float(putative_keyword, putative_value))

                    elif putative_keyword == "QUENCH_START":
                        self.keyword_lookup['QUENCH_START'] = self._kw_float(putative_keyword, putative_value)

                    elif putative_keyword == "QUENCH_END":
                        self.keyword_lookup['QUENCH_END'] = self._kw_float(putative_keyword, putative_value)

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## TSMMC keywords
                    elif putative_keyword == 'TSMMC_JUMP_TEMP':
                        self.keyword_lookup['TSMMC_JUMP_TEMP'] = self._kw_float(putative_keyword, putative_value)

                    elif putative_keyword == 'TSMMC_STEP_MULTIPLIER':
                        self.keyword_lookup['TSMMC_STEP_MULTIPLIER'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'TSMMC_NUMBER_OF_POINTS':
                        self.keyword_lookup['TSMMC_NUMBER_OF_POINTS'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'TSMMC_FIXED_OFFSET':
                        self.keyword_lookup['TSMMC_FIXED_OFFSET'] = self._kw_float(putative_keyword, putative_value)

                    elif putative_keyword == 'TSMMC_INTERPOLATION_MODE':
                        self.keyword_lookup['TSMMC_INTERPOLATION_MODE'] = str(putative_value).upper().strip()
                        if self.keyword_lookup['TSMMC_INTERPOLATION_MODE'] not in ['LINEAR']:
                            raise KeyFileException(latticeExceptions.message_preprocess('Tried to set TSMMC_INTERPOLATION_MODE with unexpected keyword [%s]' % (self.keyword_lookup['TSMMC_INTERPOLATION_MODE'])))

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## CRANKSHAFT keywords
                    elif putative_keyword == "CRANKSHAFT_SUBSTEPS":
                        self.keyword_lookup['CRANKSHAFT_SUBSTEPS'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == "SLITHER_SUBSTEPS":
                        self.keyword_lookup['SLITHER_SUBSTEPS'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == "PULL_SUBSTEPS":
                        self.keyword_lookup['PULL_SUBSTEPS'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == "VMMC_MAX_DISPLACEMENT":
                        self.keyword_lookup['VMMC_MAX_DISPLACEMENT'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == "VMMC_MAX_CLUSTER":
                        self.keyword_lookup['VMMC_MAX_CLUSTER'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == "CRANKSHAFT_MODE":
                        # obsolete: accepted for backwards compatibility with old
                        # keyfiles, announced, and ignored (the crankshaft always
                        # performs a fixed CRANKSHAFT_SUBSTEPS sub-moves per
                        # megamove, which is the old UNIFORM mode)
                        print("[ WARNING ] : CRANKSHAFT_MODE is obsolete and ignored - the crankshaft always performs a fixed CRANKSHAFT_SUBSTEPS sub-moves per megamove (the old UNIFORM mode). Remove the keyword from the keyfile.")
                        self.keyword_lookup['CRANKSHAFT_MODE'] = CONFIG.DEFAULTS['CRANKSHAFT_MODE']

                    elif putative_keyword == "NON_INTERACTING":
                        self.keyword_lookup['NON_INTERACTING'] = self._kw_bool(putative_keyword, putative_value)

                    elif putative_keyword == "ANGLES_OFF":
                        self.keyword_lookup['ANGLES_OFF'] = self._kw_bool(putative_keyword, putative_value)

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## Analysis keywords
                    elif putative_keyword == "ANALYSIS_FREQ":
                        self.keyword_lookup['ANALYSIS_FREQ'] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_POL':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_INTSCAL':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_DISTMAP':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_ACCEPTANCE':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_INTER_RESIDUE':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_CLUSTER':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_CLUSTER_THRESHOLD':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_CUSTOM':
                        self.keyword_lookup[putative_keyword] = self._kw_int(putative_keyword, putative_value)

                    elif putative_keyword == 'ANALYSIS_MODULE':
                        self.keyword_lookup[putative_keyword] = self._kw_path(putative_keyword, putative_value)

                    elif putative_keyword == 'ANA_RESIDUE_PAIRS':
                        split_residues = putative_value.split()
                        if len(split_residues) != 2:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid ANA_RESIDUE_PAIRS format [{putative_value}]. Expected: ANA_RESIDUE_PAIRS : <res1> <res2>'))
                        try:
                            res1 = int(split_residues[0])
                            res2 = int(split_residues[1])
                        except ValueError:
                            raise KeyFileException(latticeExceptions.message_preprocess(f'Invalid ANA_RESIDUE_PAIRS format [{putative_value}]. Residue indices must be integers'))

                        if res1 > res2:
                            tmp = res1
                            res1 = res2
                            res2 = tmp

                        if putative_keyword in self.keyword_lookup:
                            self.keyword_lookup['ANA_RESIDUE_PAIRS'].append([res1, res2])
                        else:
                            self.keyword_lookup['ANA_RESIDUE_PAIRS'] = [[res1, res2]]

                    elif putative_keyword == 'WRITE_CHAIN_TO_CHAINID':
                        self.keyword_lookup[putative_keyword] = self._kw_bool(putative_keyword, putative_value)

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## MOVESET keywords
                    elif putative_keyword[0:4] == 'MOVE':
                        if putative_keyword == 'MOVE_CRANKSHAFT':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_CHAIN_TRANSLATE':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_CHAIN_ROTATE':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_CHAIN_PIVOT':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_HEAD_PIVOT':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_SLITHER':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_CLUSTER_TRANSLATE':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_CLUSTER_ROTATE':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_CTSMMC':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_MULTICHAIN_TSMMC':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_PULL':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_VMMC':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_SYSTEM_TSMMC':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                        elif putative_keyword == 'MOVE_JUMP_AND_RELAX':
                            self.keyword_lookup[putative_keyword] = self._kw_float(putative_keyword, putative_value)

                    ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    ## SAVING keywords
                    elif putative_keyword == 'SAVE_AT_END':
                        self.keyword_lookup['SAVE_AT_END'] = self._kw_bool(putative_keyword, putative_value)

                    elif putative_keyword == 'SAVE_EQ':
                        self.keyword_lookup['SAVE_EQ'] = self._kw_bool(putative_keyword, putative_value)

                    elif putative_keyword == 'TRAJECTORY_PBC_UNWRAP':
                        self.keyword_lookup['TRAJECTORY_PBC_UNWRAP'] = self._kw_bool(putative_keyword, putative_value)

                    # PARALLELIZE - run the crankshaft, slither and pull moves on the parallel kernels
                    elif putative_keyword == 'PARALLELIZE':
                        self.keyword_lookup['PARALLELIZE'] = self._kw_bool(putative_keyword, putative_value)

                    # PARALLEL_THREADS - number of OpenMP threads (0 => auto)
                    elif putative_keyword == 'PARALLEL_THREADS':
                        self.keyword_lookup['PARALLEL_THREADS'] = self._kw_int(putative_keyword, putative_value)

                    else:
                        raise KeyFileException(latticeExceptions.message_preprocess('Fail to deal with a supported keyword - [%s] - this is a bug! ' % putative_keyword))

                else:
                    raise KeyFileException(latticeExceptions.message_preprocess('Found an unsupported keyword - [%s]. Valid supported keywords are\n%s ' % (putative_keyword, self.expected_keywords)))

    def update_keyfile(self, update_dictionary):
        """
        Function that takes an update dictionary of key : value pairs
        and updates the current key:value pairs with these new keywords.

        Parameters
        -----------------------

        update_dictionary : dict
            Key-value pairs where keys are keywords and values are the values to change to. Any key-value pair is
            added, overwriting an existing entry with the same key. Note that no validation of the keywords or their
            values is performed here.

        Returns
        -----------------------
        None
            No return value, but updates the self.keyword_lookup dictionary

        """
        self.keyword_lookup.update(update_dictionary)


    def write_keyfile(self, output_filename, PADDING=10):
        """
        Write the key-value pairs from the ``keyword_lookup`` dictionary to a file.

        Each keyword and its value are written on a single line in
        ``KEY : value`` format, with padding inserted between the key and the
        value for readability.

        Parameters
        ----------
        output_filename : str
            The name (path) of the file to write the key-value pairs to.
        PADDING : int, optional
            The column width the keyword name is padded out to, so that the values
            line up. Keywords at least this long still get a single space. Default
            is 10.

        Returns
        -------
        None
            No return value; the keyword/value pairs are written to disk.

        Examples
        --------
        >>> parser = KeyFileParser('input.kf')
        >>> parser.write_keyfile('output.txt')
        """

        write_keyword_lookup(self.keyword_lookup, output_filename, PADDING=PADDING,
                             expected_keywords=self.expected_keywords,
                             keywords_with_multiple_entries=self.keywords_with_multiple_entries)


    def _check_restart_continue(self, restart_object):
        """
        Refuse a RESTART_CONTINUE keyfile that could not be an exact continuation.

        Resuming a run bit for bit needs three things the restart file may or may
        not hold (the step it was written at, the generator states at that step,
        the temperature in force) and a run that is the same system: the same
        box and boundary condition, the same chains, no seed of its own, and a
        total length beyond the checkpoint. Each is checked here and each failure
        says what to change.

        Parameters
        ----------
        restart_object : restart.RestartObject
            The loaded restart file, after the box and boundary reconciliation.

        Raises
        ------
        KeyFileException
            For any of the conditions above.
        """
        kl = self.keyword_lookup
        problems = []

        # check we can get continuation state from the restart file...
        if restart_object.step is None or restart_object.rng_python is None:
            problems.append("the restart file carries no continuation state (it was written by a "
                            "PIMMS older than 1.0.8, or by a RestartObject built outside a run); "
                            "it can seed a new run but cannot be resumed exactly")

        # check no seed was provided in the keyfile
        if self._seed_was_given:
            problems.append("SEED must not be given: the random-number generators are restored from "
                            "the restart file, so a seed would be a contradiction (remove SEED)")

        # check that the run is the same system as the restart file
        if kl['RESIZED_EQUILIBRATION']:
            problems.append("RESIZED_EQUILIBRATION cannot be combined with a continuation - the box "
                            "of a resumed run is the box of the run it resumes")

        # check no extra chains were provided
        if kl['EXTRA_CHAIN']:
            problems.append("EXTRA_CHAIN cannot be combined with a continuation - adding chains "
                            "makes it a different system")

        # check that restart file dimensions and proposed dimensions match, unless the user has explicitly overridden them
        if not kl.get('RESTART_OVERRIDE_DIMENSIONS') and \
                list(kl['DIMENSIONS']) != list(restart_object.dimensions):
            problems.append("DIMENSIONS %s differ from the restart file's %s; a continuation runs in "
                            "the same box (RESTART_OVERRIDE_DIMENSIONS : True is the easy way to "
                            "guarantee that)" % (list(kl['DIMENSIONS']), list(restart_object.dimensions)))

        # check hardwall status matches
        if bool(kl['HARDWALL']) != bool(restart_object.hardwall):
            problems.append("HARDWALL %s differs from the restart file's %s; a continuation keeps the "
                            "boundary condition (RESTART_OVERRIDE_HARDWALL : True is the easy way to "
                            "guarantee that)" % (bool(kl['HARDWALL']), bool(restart_object.hardwall)))

        # check the temperature. The continuation always runs at the temperature the
        # checkpoint recorded, while the Hamiltonian (the ANGLE_PENALTY_T_NORM scaling)
        # is built from this keyfile's EQUILIBRIUM_TEMPERATURE, so a keyfile that
        # disagrees would silently give a run that is neither the original nor the
        # one requested.
        checkpoint_T = restart_object.temperature
        if checkpoint_T is not None:
            if kl['QUENCH_RUN']:
                lo = min(kl['QUENCH_START'], kl['QUENCH_END'])
                hi = max(kl['QUENCH_START'], kl['QUENCH_END'])
                tol = 1e-9 * max(1.0, abs(hi))
                if not (lo - tol <= checkpoint_T <= hi + tol):
                    problems.append("the restart file was written at temperature %.6g, which is outside "
                                    "this keyfile's quench range (QUENCH_START %.6g -> QUENCH_END %.6g); "
                                    "a continuation must use the quench settings of the run it resumes"
                                    % (checkpoint_T, kl['QUENCH_START'], kl['QUENCH_END']))
            elif abs(kl['TEMPERATURE'] - checkpoint_T) > 1e-9 * max(1.0, abs(checkpoint_T)):
                problems.append("TEMPERATURE %.6g differs from the restart file's %.6g; a continuation "
                                "runs at the temperature it was checkpointed at (set TEMPERATURE : %.6g)"
                                % (kl['TEMPERATURE'], checkpoint_T, checkpoint_T))

        # check that the run is long enough to continue
        if restart_object.step is not None and kl['N_STEPS'] <= restart_object.step:
            problems.append("N_STEPS (%d) is the TOTAL length of the run and must exceed the step the "
                            "restart file was written at (%d)" % (kl['N_STEPS'], restart_object.step))
        if problems:
            raise KeyFileException(
                "RESTART_CONTINUE : True cannot resume this run:\n  - " + "\n  - ".join(problems))

    def _warn_if_too_many_chain_types(self):
        """
        Warn when the run has more chain types than PDB chain identifiers.

        Every chain type has its own identifier in ``START.pdb``; there are 62
        (``A-Z``, ``a-z``, ``0-9``) and every type past them shares the last one.
        Before restart processing each ``CHAIN`` / ``EXTRA_CHAIN`` line is a type.
        After it, the ``CHAIN`` list has been rebuilt one line per type from the
        restart file, with any ``EXTRA_CHAIN`` already merged (one whose sequence
        matches an existing type joins that type), so the extra chains are not
        added twice.
        """

        # get number of chain types from CHAIN keyword
        n_types = len(self.keyword_lookup['CHAIN'])

        # define restart file
        restart = self.keyword_lookup['RESTART_FILE']

        # RESTART_FILE is tri-state and each state answers "are the EXTRA_CHAIN types
        # already inside the CHAIN count above?" differently, so we have to look at the
        # type and not just at truthiness:
        #
        #   False/None       - no restart run, so nothing has been folded into CHAIN and
        #                      the extra chains (if any) still have to be added here. Note
        #                      this is the state at the pre-restart call site, and the
        #                      truthiness test alone would send it down the wrong branch
        #                      because isinstance(False, str) is False.
        #   str              - the restart file has been read from the keyfile but not yet
        #                      loaded, so again CHAIN knows nothing about the extra chains.
        #   RestartObject    - the restart file has been loaded and CHAIN was rebuilt from
        #                      the restart chains PLUS the extra chains, so adding
        #                      EXTRA_CHAIN here would count those types twice.
        if not restart or isinstance(restart, str):
            n_types += len(self.keyword_lookup['EXTRA_CHAIN'])

        # nb max 62 PDB chain identifiers comes from: A-Z, a-z, 0-9
        if n_types > 62:
            print(f"[ WARNING ] : Found {n_types} chain types (more than the 62 PDB chain identifiers available). Chain types past the 62nd all share the identifier '9' in START.pdb, so analyses that read chain types from the PDB will merge them")

    def _check_cluster_rotate_box(self):
        """
        Trigger and exception if we're trying to define a keyfile where cluster rotation 
        moves are set to >0  on a non-cubic/non-square periodic production box.

        Uses the FINAL ``HARDWALL`` and ``DIMENSIONS`` values, so it must run after
        any restart-file override of either keyword.

        Raises
        ------
        KeyFileException
            If ``MOVE_CLUSTER_ROTATE`` is enabled, ``HARDWALL`` is off and the
            production box has unequal axes.
        """
        dims = self.keyword_lookup['DIMENSIONS']
        if (not self.keyword_lookup['HARDWALL']) and self.keyword_lookup['MOVE_CLUSTER_ROTATE'] > 0:
            offending_box = None
            if len(set(dims)) != 1:
                offending_box = ('production (DIMENSIONS)', dims)
            if offending_box is not None:
                which, box = offending_box
                raise KeyFileException(
                    f'MOVE_CLUSTER_ROTATE cannot be used with a non-cubic/non-square {which} box '
                    f'= {box} under periodic boundaries (HARDWALL : False). A cluster rotation is a rigid '
                    '90/180/270 degree rotation applied as an energy-neutral move; a 90/270 degree rotation '
                    'swaps two axes, and under periodic boundaries on unequal axes that swap changes '
                    'intra-cluster minimum-image distances, so the rotation is no longer energy-preserving and '
                    'the tracked energy would drift from the true energy. The production box therefore '
                    'must be cubic/square when HARDWALL is off. A RESIZED_EQUILIBRATION box is exempt '
                    'because that phase always uses hard-wall boundaries. To use cluster '
                    'rotation, either (a) make all boxes cubic/square, (b) turn HARDWALL on (a rotation is a '
                    'valid isometry when there is no periodic wrapping), or (c) set MOVE_CLUSTER_ROTATE : 0.')


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

        Returns
        -------
        None
            No return value, but ``self.keyword_lookup`` is updated in place
            (e.g. derived ``__`` keywords are set, chains may be upper-cased,
            restart/freeze files are loaded, and quench parameters normalised).

        Raises
        ------
        KeyFileException
            If any sanity check fails (e.g. out-of-range numerical values, move
            fractions that do not sum to 1.0, missing/invalid parameter or
            restart files, or incompatible box/hardwall/experimental settings).

        RestartException
            Propagated from ``sanity_check_and_update_with_restart_file()`` if the
            restart file is incompatible with the keyfile.

        """

        ## ----------------------------------------------------------------------------------------
        ## First we case chain and extra chain keywords, and if yes cast all chains to upper case
        #  (both CHAIN and EXTRA_CHAIN)
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
        ## check values we think must be bigger than 0 (or can be unset)
        # 
        # NOTE: the ANA_* frequencies and ENERGY_CHECK are deliberately ABSENT
        # from this >0 list - set_dynamic_defaults runs first and rewrites any
        # value < 1 to "beyond the run length" (the documented "0 disables this
        # analysis" convention), so a >0 check on them could never fire anyway
        for c in ['TEMPERATURE', 'N_STEPS',  'PRINT_FREQ', 'XTC_FREQ', 'EN_FREQ', 'SEED', 'RESTART_FREQ', 'QUENCH_STEPSIZE', 'QUENCH_START', 'QUENCH_END',  'TSMMC_STEP_MULTIPLIER', 'TSMMC_NUMBER_OF_POINTS',  'CRANKSHAFT_SUBSTEPS', 'SLITHER_SUBSTEPS', 'PULL_SUBSTEPS', 'VMMC_MAX_DISPLACEMENT', 'VMMC_MAX_CLUSTER']:

            try:

                # we have to check this unset because several keywords default to 'UNSET' which triggers an error if we try and compare to an int
                if self.keyword_lookup[c] != 'UNSET':
                    if self.keyword_lookup[c] <= 0:
                        raise KeyFileException(latticeExceptions.message_preprocess(f'Numerical error when parsing keyfile. Expected {c} to be larger than 0'))
            except TypeError:
                # the comparison failed because the value is a non-numeric type
                raise KeyFileException(latticeExceptions.message_preprocess(
                    f'Keyword {c} could not be checked to be greater than 0 - its value [{self.keyword_lookup[c]}] is not numeric'))
        

        ## ------------------------------------------------------------------
        ## check values with think must be bigger than or equal to zero
        #
        for c in ['ANA_CLUSTER_THRESHOLD', 'EQUILIBRATION', 'PARALLEL_THREADS']:

            try:
                if self.keyword_lookup[c] < 0:
                    raise KeyFileException(latticeExceptions.message_preprocess(f'Numerical error when parsing keyfile. Expected {c} to be 0 or larger'))
            except TypeError:
                # the comparison failed because the value is a non-numeric type
                raise KeyFileException(latticeExceptions.message_preprocess(
                    f'Keyword {c} could not be checked to be >= 0 - its value [{self.keyword_lookup[c]}] is not numeric'))

        
        ## ------------------------------------------------------------------
        ## CHAIN CHECKS
        # every CHAIN / EXTRA_CHAIN line is a chain type with its own PDB chain
        # identifier; there are 62 labels (A-Z, a-z, 0-9) and every type past
        # them shares the last one
        # (for a restart run the CHAIN list is only known after the restart file
        # has been read, so the same check is repeated at the end of that step)
        if not self.keyword_lookup['RESTART_FILE']:
            self._warn_if_too_many_chain_types()


        ## ---------------------------------------------------------
        ## MOVESET CHECK
        # check moveset frequencies add up to 1
        running_total = 0.0
        message = ''                            # message defines a summary of all moves for easy debugging

        
        for i in self.expected_keywords:
            # for each MOVE keyword [note this dynamically finds them so if new MOVE_ keywords are added this will just work :-) ]
            if i[0:4] == "MOVE":
                # an individual move fraction must be a sensible probability (>= 0);
                # a negative fraction is meaningless even if the total still sums to 1
                if self.keyword_lookup[i] < 0:
                    raise KeyFileException(latticeExceptions.message_preprocess(
                        'Move fraction %s cannot be negative (got %s)' % (i, self.keyword_lookup[i])))
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

            # the quench frequency must be a positive number of steps: 0 makes
            # steps_for_quench zero (EQUILIBRATION silently set to 0 with
            # QUENCH_AS_EQUILIBRATION) and negative values only fail deep at runtime
            if not self.keyword_lookup['QUENCH_FREQ'] > 0:
                raise KeyFileException('QUENCH_FREQ must be a positive number of steps (got %s)' % self.keyword_lookup['QUENCH_FREQ'])


            # check that the temperature distance being traversed isn't smaller than the temperature 
            # step size...
            dT = abs(self.keyword_lookup['QUENCH_START'] - self.keyword_lookup['QUENCH_END'])
            # a relative 1e-9 slack matches the runtime's snap-to-target, so
            # 0.3 -> 0.2 in steps of 0.1 (dT = 0.09999999999999998) is one rung,
            # not an overshoot
            if dT < self.keyword_lookup['QUENCH_STEPSIZE'] * (1.0 - 1e-9):
                raise KeyFileException('A single quench step (QUENCH_STEPSIZE = %s) overshoots the full QUENCH_START -> QUENCH_END temperature range (%s); reduce QUENCH_STEPSIZE or widen the temperature range' % (self.keyword_lookup['QUENCH_STEPSIZE'], dT))

            # number of temperature changes: rounded up, but an integer ratio up to
            # floating-point round-off counts as that integer (truncating dT to an
            # integer used to under-count fractional ramps; a plain ceil on binary
            # floats over-counted 1.0 -> 0.7 by 0.1 as four rungs instead of three)
            self.quench_rungs = _quench_rung_count(dT, self.keyword_lookup['QUENCH_STEPSIZE'])
            steps_for_quench = (1 + self.quench_rungs) * self.keyword_lookup['QUENCH_FREQ']
                
            if steps_for_quench >= self.keyword_lookup['N_STEPS']:
                raise KeyFileException('This quench will not complete as the quench period [%i] is longer than the number of steps in the simulation [%i]' %(steps_for_quench, self.keyword_lookup['N_STEPS'] ))
                                
            # If we're doing a heating simulation (so start is lower than end) then stepsize must become negative
            # as the temperature update operation is CURRENT - STEPSIZE
            if self.keyword_lookup['QUENCH_START'] < self.keyword_lookup['QUENCH_END']:
                self.keyword_lookup['QUENCH_STEPSIZE'] = -(self.keyword_lookup['QUENCH_STEPSIZE'])


            # Update the equilibration period so the quench process is enveloped by the quench (i.e. all
            # output is only at the target temperature)
            if self.keyword_lookup['QUENCH_AS_EQUILIBRATION'] :
                # EQUILIBRATION must be an integer number of steps (it is used in
                # integer/range contexts downstream); round up so the equilibration
                # period fully envelopes the quench.
                equilibration_steps = int(np.ceil(steps_for_quench))
                print("UPDATING EQUILIBRATION TO [%i] (temperature quench period is equilibration)" % equilibration_steps)
                self.keyword_lookup['EQUILIBRATION'] = equilibration_steps


            if not self.keyword_lookup['TEMPERATURE'] == self.keyword_lookup['QUENCH_START']:
                print("[ WARNING ] : Resetting the starting temperature to the QUENCH_START value [%s]"%self.keyword_lookup['QUENCH_START'])
                self.keyword_lookup['TEMPERATURE'] = self.keyword_lookup['QUENCH_START'] 

        ## ---------------------------------------------------------
        ## TSMMC check (sets the __TSMMC_ON keyword
        if self.keyword_lookup['MOVE_CTSMMC'] > 0 or self.keyword_lookup['MOVE_MULTICHAIN_TSMMC'] > 0 or self.keyword_lookup['MOVE_SYSTEM_TSMMC'] > 0:

            # 
            self.keyword_lookup['__TSMMC_USED'] = True

            # check if we didn't use a fixed OFFSET that the jumps make sense...
            # a fixed offset must be a positive temperature increment - a negative
            # value silently builds a COOLING "excursion" schedule
            if self.keyword_lookup['TSMMC_FIXED_OFFSET'] is not False and not self.keyword_lookup['TSMMC_FIXED_OFFSET'] > 0:
                raise KeyFileException('TSMMC_FIXED_OFFSET must be a positive temperature increment (got %s)' % self.keyword_lookup['TSMMC_FIXED_OFFSET'])

            # The excursion must climb ABOVE the current simulation temperature.
            # In a quench run that temperature moves between QUENCH_START and
            # QUENCH_END, so check against the hottest temperature the run will
            # visit - otherwise a heating quench passes the parser and aborts
            # mid-run the first time a TSMMC move is drawn above the jump temperature.
            hottest = self.keyword_lookup['TEMPERATURE']
            if self.keyword_lookup['QUENCH_RUN']:
                hottest = max(self.keyword_lookup['QUENCH_START'], self.keyword_lookup['QUENCH_END'])
            if (self.keyword_lookup['TSMMC_JUMP_TEMP'] <= hottest) and self.keyword_lookup['TSMMC_FIXED_OFFSET'] is False:
                raise KeyFileException(latticeExceptions.message_preprocess('\n\nThe TSMMC jump temperature [%3.2f] is less than or equal to the highest simulation temperature the run will visit [%3.2f], which will mean the TSMMC moves will at best hurt performance and at worst reduce sampling (in a quench run the jump temperature must exceed both QUENCH_START and QUENCH_END). Please correct your keyfile appropriately' % (self.keyword_lookup['TSMMC_JUMP_TEMP'], hottest)))

        else:
            self.keyword_lookup['__TSMMC_USED'] = False


        ## ------------------------------------------------------------
        ## residue distance pairs check - make sure none of the pair-pair
        ## analysis indices fall outside of the chain length. Negative indices
        ## must be rejected too: they parse as valid integers, but downstream
        ## they are used directly as Python list indices, so e.g. -1 would
        ## silently measure the distance from the LAST residue rather than
        ## erroring.
        for pair in self.keyword_lookup['ANA_RESIDUE_PAIRS']:
            if pair[0] < 0:
                raise KeyFileException('Residue-residue distance analysis index (%i) is negative - residue indices count from 0' % pair[0])
            # a restart run ignores the keyfile's CHAIN lines, so its bound check
            # runs after restart processing against the chains actually loaded
            if not self.keyword_lookup['RESTART_FILE']:
                for chain in self.keyword_lookup['CHAIN']:
                    if pair[1] >= len(chain[1]):
                        raise KeyFileException('Residue-residue distance analysis pair (%i) is outside the chain length (%i)' % (pair[1], len(chain[1])))
                    
        ## ------------------------------------------------------------
        ## Non-cubic / non-square boxes are fully supported (2D and 3D, hardwall OR
        ## periodic): the engine wraps and evaluates every axis independently
        ## (per-axis minimum image / grid.shape[i]), verified move-by-move against a
        ## from-scratch energy recompute (ENERGY_CHECK). The ONE restriction is
        ## cluster ROTATION under periodic boundaries.
        ##
        ## A cluster move rotates a whole (isolated) cluster rigidly and is applied as
        ## an ENERGY-NEUTRAL move (no energy recomputation). A cardinal 90/270 degree
        ## rotation swaps two axes; under periodic boundaries on a box whose axes have
        ## different periods, that swap changes the minimum-image distance of any
        ## intra-cluster pair separated by more than half a period on a swapped axis,
        ## so the cluster's internal energy actually changes while the move assumes it
        ## did not - silently corrupting the tracked energy. It is fine on a
        ## cube/square (the rotation is an exact box symmetry) and fine under HARDWALL
        ## (no periodic wrapping) - both verified via ENERGY_CHECK.
        ##
        ## Only the PRODUCTION box (DIMENSIONS) needs to satisfy this: the
        ## resized-equilibration phase is ALWAYS run with a hardwall regardless of
        ## the keyfile (Simulation.__init__ forces it), and under hardwall a rigid
        ## rotation is a valid isometry on any box shape - so a non-cubic
        ## RESIZED_EQUILIBRATION box is fine.
        # The guard is a method so it can run again AFTER restart processing,
        # which may replace HARDWALL (RESTART_OVERRIDE_HARDWALL) and DIMENSIONS
        # (RESTART_OVERRIDE_DIMENSIONS) with the restart file's values; checking
        # only the keyfile's values here let a hardwall keyfile + periodic restart
        # + override start a periodic non-cubic run with cluster rotation, and
        # rejected the converse (periodic keyfile, hardwall restart) spuriously.
        if not self.keyword_lookup['RESTART_FILE']:
            self._check_cluster_rotate_box()


        
        ## Crash if dimensions are < 7 in 
        ## 
        dims = self.keyword_lookup['DIMENSIONS']
        if len(dims) == 2:
            if dims[0] < 7 or dims[1] < 7:
                raise KeyFileException('Box size is too small to correctly support super-long range interactions, must be >= 7 ')
        else:
            if dims[0] < 7 or dims[1] < 7 or dims[2] < 7:
                raise KeyFileException('Box size is too small to correctly support super-long range interactions, must be >= 7 ')


        ## Check out resized equilibrium variables and fix as needed
        ##
        if self.keyword_lookup['RESIZED_EQUILIBRATION']:
            if len(self.keyword_lookup['RESIZED_EQUILIBRATION']) != len(self.keyword_lookup['DIMENSIONS']):
                raise KeyFileException('Number of dimensions for compressed equilibration and final simulation are not the same')

            # resized_equilibration 
            for (real, eq) in zip(self.keyword_lookup['DIMENSIONS'], self.keyword_lookup['RESIZED_EQUILIBRATION']):
                if eq > real:
                    raise KeyFileException('Resized equilibration dimension is larger than final dimension - not yet supported (only smaller)')

                # the equilibration box IS simulated, so it must satisfy the same
                # floor as DIMENSIONS (super-long-range interactions need > 7 sites)
                if eq < 7:
                    raise KeyFileException('RESIZED_EQUILIBRATION box is too small to correctly support super-long range interactions, every axis must be >= 7 (got %s)' % (self.keyword_lookup['RESIZED_EQUILIBRATION'],))

            if self.keyword_lookup['EQUILIBRATION'] == 0:
                print("[WARNING]: using RESIZED_EQUILIBRATION without an equilibration period makes no sense. Deactivating RESIZED_EQUILIBRATION"
                      + (" and EQUILIBRATION_OFFSET" if self.keyword_lookup['EQUILIBRATION_OFFSET'] else ""))
                self.keyword_lookup['RESIZED_EQUILIBRATION'] = False
                # the offset only places the equilibration box, so it goes too;
                # leaving it set made the check below blame the user for a
                # RESIZED_EQUILIBRATION they had in fact given
                self.keyword_lookup['EQUILIBRATION_OFFSET'] = False

        if self.keyword_lookup['EQUILIBRATION_OFFSET']:
            if not self.keyword_lookup['RESIZED_EQUILIBRATION']:
                raise KeyFileException('RESIZED_EQUILIBRATION MUST be turned on if EQUILIBRATION_OFFSET is specified')

            if len(self.keyword_lookup['RESIZED_EQUILIBRATION']) != len(self.keyword_lookup['EQUILIBRATION_OFFSET']):
                raise KeyFileException('Number of dimensions for compressed equilibration and offset are not the same')

            for (real, eq, offset) in zip(self.keyword_lookup['DIMENSIONS'], self.keyword_lookup['RESIZED_EQUILIBRATION'], self.keyword_lookup['EQUILIBRATION_OFFSET']):
                if offset < 0:
                    raise KeyFileException('EQUILIBRATION_OFFSET values must be >= 0 (got %s)' % (self.keyword_lookup['EQUILIBRATION_OFFSET'],))
                if eq + offset > real:
                    raise KeyFileException('Equilibration dimension + offset is larger than final dimension')       

        ##
        ## if analysis code is provided check it can be loaded
        # first check if the file even exists 
        if not os.path.isfile(self.keyword_lookup['PARAMETER_FILE']):
                raise KeyFileException('Unable to find parameter file at location [%s]. Please verify the file exists.' %(self.keyword_lookup['PARAMETER_FILE']))

            # check if we can load as a python module...

      
        ##
        ## if custom analysis code is provided, load and validate it now (fail fast
        ## at parse time rather than part-way through a run). The loader raises a
        ## clear KeyFileException on any problem and, on success, returns the
        ## validated analysis_function callable - which we store in place of the
        ## path so the rest of PIMMS holds a ready-to-call function.
        # the converse of the ANALYSIS_MODULE/ANA_CUSTOM check below: a frequency
        # without a module silently does nothing, so say so
        if (not self.keyword_lookup['ANALYSIS_MODULE']) and not getattr(self, '_ana_custom_disabled', True):
            print("[ WARNING ] : ANA_CUSTOM is set but no ANALYSIS_MODULE was given - no custom analysis will run")

        if self.keyword_lookup['ANALYSIS_MODULE']:
            original_path = self.keyword_lookup['ANALYSIS_MODULE']
            self.keyword_lookup['ANALYSIS_MODULE'] = file_utilities.custom_analysis_module_import(original_path)
            # the keyword now holds the imported callable; keep the path so a
            # written keyfile (keyfile_used.kf) still names the module
            self.keyword_lookup['__ANALYSIS_MODULE_PATH'] = original_path
            print("Loaded custom analysis code from [%s]" % original_path)

            # a validated module that never runs is almost certainly a mistake:
            # ANA_CUSTOM defaults to 0 (disabled), so require an explicit frequency.
            # NB: test the RAW pre-dynamic-defaults state - set_dynamic_defaults has
            # already rewritten a disabled ANA_CUSTOM to N_STEPS+10, so testing the
            # rewritten value here would make this check dead code.
            if getattr(self, '_ana_custom_disabled', False):
                raise KeyFileException(
                    'ANALYSIS_MODULE was given (and loaded successfully from [%s]) but '
                    'ANA_CUSTOM is 0/unset, so the custom analysis would never run. Set '
                    'ANA_CUSTOM to the desired frequency (in steps), or remove '
                    'ANALYSIS_MODULE.' % original_path)


        ##
        ## Check restart file and then read in
        if self.keyword_lookup.get('RESTART_CONTINUE') and not self.keyword_lookup.get('RESTART_FILE'):
            raise KeyFileException('RESTART_CONTINUE : True requires a RESTART_FILE to continue from')

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

            # finally we ask is if any EXTRA_CHAIN were provided, and if yes add these to the newly generated
            # RestartObject. Note that the RestartObject sanity checks the new chains as they are added, so if 
            # there is an issue with the new chains this will be caught and raised as a RestartException 
            # which we catch and re-raise as a KeyFileException to keep things consistent for the user.
            if len(self.keyword_lookup['EXTRA_CHAIN']) > 0:
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

            # re-run the ANA_RESIDUE_PAIRS bound check now the restart has rebuilt
            # the CHAIN list: the earlier check in run_sanity_checks loops over the
            # keyfile's own CHAIN lines, which are empty for a restart run, so
            # out-of-range pairs sailed through; EXTRA_CHAIN chains were never
            # examined at all.
            _all_chains = list(self.keyword_lookup['CHAIN']) + list(self.keyword_lookup.get('EXTRA_CHAIN') or [])
            for pair in self.keyword_lookup['ANA_RESIDUE_PAIRS']:
                for chain in _all_chains:
                    if pair[1] >= len(chain[1]):
                        raise KeyFileException('Residue-residue distance analysis pair (%i) is outside the chain length (%i)' % (pair[1], len(chain[1])))
            # CHAIN now holds every chain type of the run (extra chains included)
            self._warn_if_too_many_chain_types()
            ## ----------------------------------------------------------------------------------------------------

        ## else if we DID not pass a restart file, check we didn't have 
        else:
            if len(self.keyword_lookup['EXTRA_CHAIN']) > 0:
                raise KeyFileException('\n\nEXTRA_CHAIN keyword defined but no restart file was provided. This is not supported. Please provide a restart file or remove the EXTRA_CHAIN line.')

        ##
        ## Check for and parse the freeze file
        if self.keyword_lookup['FREEZE_FILE']:

            # note all sanity checking is done in the FreezeFile class
            self.keyword_lookup['FREEZE_FILE'] = FreezeFile(self.keyword_lookup['FREEZE_FILE'])

        ## ---------------------------------------------------------        
        # Check for missuse of experimental features...

        for kw in CONFIG.EXPERIMENTAL_KEYWORDS:

            # if a keyword included in the EXPERIMENTAL_KEYWORDS list is NOT set to the default value            
            if self.keyword_lookup[kw] != self.DEFAULTS[kw]:
                self.__check_experimental_features(kw)

                ## ask if experimental features was set`
                #if self.keyword_lookup['EXPERIMENTAL_FEATURES'] == False:
                #    raise KeyFileException('\n\nExperimental or non-supported feature ({kw}) being proposed but EXPERIMENTAL_FEATURES is not set (or set to False).\n')


        ## ---------------------------------------------------------
        # check simulation is not longer than eq period
        if self.keyword_lookup['EQUILIBRATION'] >= self.keyword_lookup['N_STEPS']:
            raise KeyFileException(f"Simulation will equilibrate for longer than the total number of steps; equilibration steps = {self.keyword_lookup['EQUILIBRATION']} while total steps = {self.keyword_lookup['N_STEPS']}")
                
                
    #-----------------------------------------------------------------
    #                    
    def set_defaults(self):
        """
        Assign default values from ``self.DEFAULTS`` to any unset keywords.

        First the presence of every required keyword is validated (and at least
        one of ``CHAIN`` or ``RESTART_FILE`` must be present). Then every
        expected keyword that was not explicitly supplied in the keyfile is set
        to its default value from ``self.DEFAULTS``, announcing each default as
        it is applied.

        Returns
        -------
        None
            No return value, but ``self.keyword_lookup`` is updated in place with
            default values for any keywords that were not explicitly defined.

        Raises
        ------
        KeyFileException
            If a required keyword is missing, or if neither ``CHAIN`` nor
            ``RESTART_FILE`` was provided.
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

        Specifically, any analysis frequency (or ENERGY_CHECK) set to a value below 1 is treated as
        'disabled' and rewritten to N_STEPS + 10 so it can never fire, with the set of genuinely
        disabled keywords recorded in the derived __DISABLED_FREQUENCIES keyword. A RESTART_FREQ left
        at its 'Every 10th-percentile' default is converted into an actual number of steps.

        Returns
        -----------------------
        None
            No return value, but updates the self.keyword_lookup dictionary with new default values as
            needed, and sets self._ana_custom_disabled for use by the sanity checks.

        """

        frequency_keywords = [
            'ANALYSIS_FREQ', 'ANA_POL', 'ANA_INTSCAL', 'ANA_DISTMAP',
            'ANA_ACCEPTANCE', 'ANA_INTER_RESIDUE', 'ANA_CLUSTER',
            'ANA_CUSTOM', 'ENERGY_CHECK',
        ]

        # A frequency below one is a real disabled state, not merely a large
        # cadence.  Historically it was represented only as N_STEPS + 10, which a
        # written keyfile then round-tripped as an ENABLED analysis at that
        # cadence, and which anything that moved N_STEPS could overtake.  Retain
        # the numeric sentinel for backwards-compatible summaries, while
        # recording the state explicitly for Simulation.setup_analysis/IO and
        # write_keyword_lookup.
        disabled_frequencies = {
            keyword for keyword in frequency_keywords
            if self.keyword_lookup[keyword] < 1
        }

        # Remember the raw user/default state before any rewrite.  A supplied
        # module with a disabled ANA_CUSTOM is rejected later as likely user
        # error; conversely, without a module ANA_CUSTOM is always disabled even
        # if a stray positive frequency was supplied.
        self._ana_custom_disabled = 'ANA_CUSTOM' in disabled_frequencies
        if self.keyword_lookup['ANALYSIS_MODULE'] is False:
            disabled_frequencies.add('ANA_CUSTOM')
            self.keyword_lookup['ANA_CUSTOM'] = self.keyword_lookup['N_STEPS'] + 10

        for tmpkw in frequency_keywords:
            if self.keyword_lookup[tmpkw] < 1:
                self.keyword_lookup[tmpkw] = self.keyword_lookup['N_STEPS'] + 10

        # Tuple rather than a mutable set keeps this derived parser result stable
        # when callers retain/pass around keyword_lookup.
        self.keyword_lookup['__DISABLED_FREQUENCIES'] = tuple(sorted(disabled_frequencies))

        if self.keyword_lookup['RESTART_FREQ'] == "Every 10th-percentile":
            # deliberate effective floor being used here.. but never let it reach 0: for
            # N_STEPS < 10 that produced RESTART_FREQ = 0, which the sanity checks then
            # rejected with "Expected RESTART_FREQ to be larger than 0" - about a keyword
            # the user never set. Short runs simply write a restart every step.
            self.keyword_lookup['RESTART_FREQ'] = max(1, int(self.keyword_lookup['N_STEPS'] / 10))

                
                                  
    #-----------------------------------------------------------------
    #    
    def print_summary(self):
        """
        Print a full human-readable summary of the parsed simulation setup.

        Reports the system overview (step counts, temperatures, expected number
        of frames), box dimensions and resulting occupied volume fraction /
        solute concentration, quench settings, freeze-file settings, simulation
        components (chains), and output frequencies. A subset of these values are
        also written to the PIMMS logfile via :mod:`pimmslogger`.

        Returns
        -------
        None
            No return value; the summary is printed to standard output (and key
            values are logged).
        """

        def section(msg):
            """
            Print a formatted section header for the summary output.

            Parameters
            ----------
            msg : str
                The section title to display in the header.

            Returns
            -------
            None
                No return value; the header is printed to standard output.
            """
            IO_utils.newline()
            IO_utils.horizontal_line(hzlen=40, linechar='*')
            print("--> %s"%(msg))
            IO_utils.newline()
            


        IO_utils.status_message("KEYFILE SUMMARY",'major')

        ## caclualte expected number of steps
        if self.keyword_lookup['SAVE_EQ']:
            expected_number_of_frames = int(np.floor(self.keyword_lookup['N_STEPS']/self.keyword_lookup['XTC_FREQ']))+1
        else:
            # match the writer's real cadence: production frames are written on
            # steps that are multiples of XTC_FREQ AFTER equilibration, plus the
            # initial frame that opening the writer always records (frame 0 is the
            # starting configuration regardless of SAVE_EQ - see docs). The old
            # floor((N-EQ)/XTC)+1 formula only agreed when EQ was a multiple of XTC.
            expected_number_of_frames = (int(np.floor(self.keyword_lookup['N_STEPS']/self.keyword_lookup['XTC_FREQ']))
                                         - int(np.floor(self.keyword_lookup['EQUILIBRATION']/self.keyword_lookup['XTC_FREQ']))
                                         + 1)
        

        # a resized equilibration that saves its frames splits them across
        # eq_traj.xtc (steps 0 .. EQUILIBRATION) and traj.xtc (the post-resize
        # configuration plus the production multiples of XTC_FREQ); with SAVE_EQ
        # off no eq_ files are written and the single production count above is
        # already right
        if self.keyword_lookup.get('RESIZED_EQUILIBRATION') and self.keyword_lookup['SAVE_EQ']:
            n_eq = int(np.floor(self.keyword_lookup['EQUILIBRATION'] / self.keyword_lookup['XTC_FREQ'])) + 1
            n_prod = (int(np.floor(self.keyword_lookup['N_STEPS'] / self.keyword_lookup['XTC_FREQ']))
                      - int(np.floor(self.keyword_lookup['EQUILIBRATION'] / self.keyword_lookup['XTC_FREQ'])) + 1)
            expected_number_of_frames = f"{n_eq} (eq_traj.xtc) + {n_prod} (traj.xtc)"

        # a continuation writes the checkpoint as its frame 0 and then only the
        # XTC_FREQ multiples after the checkpoint's step
        restart_object = self.keyword_lookup.get('RESTART_FILE')
        if self.keyword_lookup.get('RESTART_CONTINUE') and getattr(restart_object, 'step', None) is not None:
            first = restart_object.step
            if not self.keyword_lookup['SAVE_EQ']:
                first = max(first, self.keyword_lookup['EQUILIBRATION'])
            expected_number_of_frames = (int(np.floor(self.keyword_lookup['N_STEPS'] / self.keyword_lookup['XTC_FREQ']))
                                         - int(np.floor(first / self.keyword_lookup['XTC_FREQ'])) + 1)

        ## print the system overview
        print("--> System Overview")
        print("Total number of steps     : %i" % self.keyword_lookup['N_STEPS'])
        print("Num. equilibration steps  : %i" % self.keyword_lookup['EQUILIBRATION'])
        print(f"Save equilibration frames : {str(self.keyword_lookup['SAVE_EQ'])}")
        print(f"Expected number of frames : {expected_number_of_frames}")
        print("Start temperature         : %3.2f" % self.keyword_lookup['TEMPERATURE'])
        print("Final temperature         : %3.2f" % self.keyword_lookup['EQUILIBRIUM_TEMPERATURE'])
        print("Lattice-to-Angstroms      : %5.2f" % self.keyword_lookup['LATTICE_TO_ANGSTROMS'])
        print(f"Autocenter                : {str(self.keyword_lookup['AUTOCENTER'])}")
        print(f"Save at end only          : {str(self.keyword_lookup['SAVE_AT_END'])}")
        
            
        # log key things
        pimmslogger.log_status("NUMBER OF STEPS     : %i" % self.keyword_lookup['N_STEPS'], timestamp=False)
        pimmslogger.log_status("EQUIL. STEPS        : %i" % self.keyword_lookup['EQUILIBRATION'], timestamp=False)
        pimmslogger.log_status(f"SAVE EQ            : {self.keyword_lookup['SAVE_EQ']}", timestamp=False)
        pimmslogger.log_status(f"XTC_OUT_FREQ       : {self.keyword_lookup['XTC_FREQ']}", timestamp=False)
        pimmslogger.log_status(f"EXPECTED # FRAMES  : {expected_number_of_frames}", timestamp=False)        
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

                # lattice units -> nm via LATTICE_TO_ANGSTROMS (x 0.1 for Angstrom -> nm)
                
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

            # lattice units -> nm via LATTICE_TO_ANGSTROMS (x 0.1 for Angstrom -> nm)
            v_in_nm_3 = self.keyword_lookup['DIMENSIONS'][0]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1*self.keyword_lookup['DIMENSIONS'][1]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1*self.keyword_lookup['DIMENSIONS'][2]*self.keyword_lookup['LATTICE_TO_ANGSTROMS']*0.1

            v_in_L  = 1e-24*v_in_nm_3
            conc = (chain_count/6.02e23)/v_in_L      
            print("Total conc. of solute(s)       : %10.12f (M)" % (conc))
        


        ##
        ## QUENCH SETTINGS SECTION
        ##
        section('Quench Settings')

        if self.keyword_lookup['QUENCH_RUN']:

            # abs(stepsize): for heating quenches QUENCH_STEPSIZE has been
            # sign-flipped to negative by the parser, which used to print a
            # negative step count; use the same ceil-based rung count as the
            # sanity check so the two agree for fractional temperature ranges.
            quenchsteps = self.keyword_lookup['QUENCH_FREQ'] * _quench_rung_count(
                abs(self.keyword_lookup['QUENCH_START'] - self.keyword_lookup['QUENCH_END']),
                self.keyword_lookup['QUENCH_STEPSIZE'])

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
            # Hence these two numbers differ by QUENCH_FREQ steps. This is _not_ a bug!!

            print("QUENCH FREQ  : %i"     % self.keyword_lookup['QUENCH_FREQ'])
            print("QUENCH START : %3.2f " % self.keyword_lookup['QUENCH_START'])
            # the stored step is negated for a heating run; report the magnitude
            # the user wrote
            print("QUENCH STEP  : %3.2f " % abs(self.keyword_lookup['QUENCH_STEPSIZE']))
            print("QUENCH END   : %3.2f " % self.keyword_lookup['QUENCH_END'])
        else:            
            print("NO TEMPERATURE QUENCH IN EFFECT")


        ##
        ## FREEZE FILE SETTINGS SECTION
        ##

        if self.keyword_lookup['FREEZE_FILE']:
            section('Freeze File Settings')

            print(f"Freeze file      : {self.keyword_lookup['FREEZE_FILE'].filename}")
            print(f"Chains to freeze : {self.keyword_lookup['FREEZE_FILE'].chains}")
            print("** NB: detailed sanity checking will happen after full initiation**")
            
        
        ##
        ## SIMULATIO COMPONENTS SECTION
        ##
        section('Simulation Components')

        # if we're running a non-interacting simulation
        if self.keyword_lookup['NON_INTERACTING']:

            print("")
            print("NOTE: This is a non-interacting simulation - all energy interactions other")
            print("      than excluded volume and angle potentials have been switched OFF")
            print("")
            

        if self.keyword_lookup['RESTART_FILE']:
            print('RESTART FILE PROVIDED: System components from restart file described above')
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
        ## OUTPUT SETTINGS SECTION
        ##
        section('Output Parameters')
        print("Print-to-screen frequency    : %i" % self.keyword_lookup['PRINT_FREQ'])
        print("XTC out frequency            : %i" % self.keyword_lookup['XTC_FREQ'])
        print("Energy out frequency         : %i" % self.keyword_lookup['EN_FREQ'])
        print("Overal analysis frequency    : %i" % self.keyword_lookup['ANALYSIS_FREQ'])
        print(f"Write chain to chainID file? : {self.keyword_lookup['WRITE_CHAIN_TO_CHAINID']}")
        IO_utils.newline(2)
        print("Keyfile fully parsed! Preparing to start the simulation...")
        IO_utils.horizontal_line()
        


    #-----------------------------------------------------------------
    #    
    def assign_default(self):
        """
        Build the ``self.DEFAULTS`` dictionary of default keyword values.

        The absolute reference default values are taken from ``CONFIG.DEFAULTS``.
        In reality these could live in a separate configuration file, but in the
        interest of keeping everything within this file we define those default
        values here.

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

        Returns
        -------
        None
            No return value, but ``self.DEFAULTS`` is populated in place (including
            analysis-frequency defaults derived from ``ANALYSIS_FREQ`` and a freshly
            generated random ``SEED``).

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
        Add keywords that are derived from already-parsed keyword values.

        This runs after parsing, default assignment and sanity checking. It sets
        the end-to-end-distance analysis frequency (``ANA_END_TO_END``) equal to
        the polymer-analysis frequency, and sets ``EQUILIBRIUM_TEMPERATURE`` to
        the final simulation temperature (the quench end temperature for a quench
        run, otherwise the standard ``TEMPERATURE``).

        Returns
        -------
        None
            No return value, but ``self.keyword_lookup`` is updated in place with
            the derived keywords ``ANA_END_TO_END`` and ``EQUILIBRIUM_TEMPERATURE``.
        """

        # we consider the end-to-end distance to be a polymeric property
        self.keyword_lookup['ANA_END_TO_END'] = self.keyword_lookup['ANA_POL'] 

        if 'ANA_POL' in self.keyword_lookup.get('__DISABLED_FREQUENCIES', ()):
            disabled = set(self.keyword_lookup['__DISABLED_FREQUENCIES'])
            disabled.add('ANA_END_TO_END')
            self.keyword_lookup['__DISABLED_FREQUENCIES'] = tuple(sorted(disabled))


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
          if the keyfile requires a resized equilibration, the restart file must be a hardwall snapshot
          no larger than RESIZED_EQUILIBRATION on any axis; its box is grown to the RESIZED_EQUILIBRATION
          box for the equilibration phase (neither override may be set), and the production part of the
          simulation is run in the keyfile's DIMENSIONS.

        Returns
        -------
        None
            No return value. ``self.keyword_lookup`` is updated in place (the
            ``CHAIN`` composition is rebuilt from the restart object, and
            ``HARDWALL`` / ``DIMENSIONS`` may be overridden), and the restart
            object's lattice dimensions may be recentred/resized. Returns early
            (without changes) if no restart file is set.

        Raises
        ------
        RestartException
            If the restart file is incompatible with the keyfile (e.g. mismatched
            dimensionality, conflicting hardwall/PBC modes, box dimensions that
            are smaller than the restart file's dimensions, or two chains of the
            same type with different sequences).

        KeyFileException
            If ``RESTART_OVERRIDE_DIMENSIONS`` pulls in a box with an axis
            shorter than 7 sites, or if the final HARDWALL/DIMENSIONS combination
            is incompatible with ``MOVE_CLUSTER_ROTATE``.

        """


        if not self.keyword_lookup['RESTART_FILE']:
            return 

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

        # Build chain composition from both restart chains and EXTRA_CHAIN entries,
        # so startup concentration/occupied-volume reporting reflects the full system.
        all_restart_chain_entries = {}
        all_restart_chain_entries.update(restart_object.chains)
        all_restart_chain_entries.update(restart_object.extra_chains)

        # this cycles through each individual chain in the restart object, where
        # chain_info[0] = chain_positions
        # chain_info[1] = chain sequence
        # chain_info[2] = chain type
        for chainID in all_restart_chain_entries:
            
            # for each chainID get the chain sequence and 
            # the chain type
            c_seq  = all_restart_chain_entries[chainID][1]
            c_type = all_restart_chain_entries[chainID][2]

            # if this is the first time we encounter this chain type create a new 
            # entry in the chain_type_dictionary. 
            if c_type not in chain_type_dictionary:
                chain_type_dictionary[c_type] = [0, c_seq]

            # next verify that a chain of type c_type matches the other chains of 
            # type c_type that we've already seen so far
            ## NOTE - other places in the code would allow this; probably should address this at
            # some point for consistency...
            if not chain_type_dictionary[c_type][1] == c_seq:
                raise RestartException('When reading in the restart file, found a chain of type %i that did not match sequence of another chain of type %s. Chain sequences are\n: %s\n%s\n' % (c_type, c_type, chain_type_dictionary[c_type][1], c_seq))


            # if all seems good increment
            chain_type_dictionary[c_type][0] = chain_type_dictionary[c_type][0] + 1

        # when we get here chain_type_dictionary is a dictionary where keys
        # are unique chain types and the values are lists of the form [ count, chain_sequence ], 

        # finally sort the types and reconstruct a chains list, that has
        # each type as a seperate entry (i.e. [[count_1, seq_1],[count_2, seq_2]] 
        # and so on

        chain_types = list(chain_type_dictionary.keys())
        chain_types.sort()
        chains = []
        for CT in chain_types:
            chains.append(chain_type_dictionary[CT])

        # at this point chains is a list of lists, where each element is of the form [ count, chain_sequence ], 
        # and the order of the list is sorted by chain type. 

        # Store restart-derived chain composition under the canonical CHAIN keyword.
        self.keyword_lookup['CHAIN'] = chains
        
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

        # A periodic snapshot may contain a chain joined across a box face.  It
        # cannot be used as the compact hardwall phase of resized equilibration,
        # irrespective of whether RESTART_OVERRIDE_HARDWALL is enabled.  This
        # check used to live only in the ``elif`` below, so the override silently
        # bypassed it.
        if not restart_object.hardwall and self.keyword_lookup['RESIZED_EQUILIBRATION']:
            raise RestartException('\n\nRestart file describes a periodic boundary simulation, but the keyfile is trying to run as simulation that has a resized equilibration step (RESIZED_EQUILIBRATION : True). These options are incompatible with one another. The initial simulation box that is then resized MUST be a hardwall simulation.\n')

        # if the restart object was a PBC simulation (i.e. not hardwall). Note this
        # also runs after RESTART_OVERRIDE_HARDWALL, which has already set HARDWALL
        # from the restart file, so the check below cannot fire in that case.
        if not restart_object.hardwall:

            # if the restart file was a PBC simulation and we are trying to run a hardwall simulation - not compatible
            if self.keyword_lookup['HARDWALL']:
                raise RestartException('\n\nRestart file describes a periodic boundary simulation, but the keyfile is trying to run as a hardwall simulation. These options are incompatible with one another. Please set HARDWALL : False (or delete the HARDWALL keyword)\n')
                
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
        # A continuation must be the same system as the run it resumes, and the
        # file must carry the state to resume from. Checked here, after the
        # HARDWALL override but BEFORE the box reconciliation below, because that
        # reconciliation grows the restart object's box to the keyfile's for a
        # hardwall restart, after which a mismatch could no longer be seen.
        if self.keyword_lookup.get('RESTART_CONTINUE'):
            self._check_restart_continue(restart_object)

        if self.keyword_lookup['RESTART_OVERRIDE_DIMENSIONS']:
            print("Setting DIMENSIONS keyword based on restart file")
            self.keyword_lookup['DIMENSIONS'] = restart_object.dimensions

            # RESTART_OVERRIDE_DIMENSIONS replaces DIMENSIONS *after* the keyfile-level
            # box-floor check ran against the keyfile's original DIMENSIONS, so re-validate
            # the overridden box here (the cluster-rotate cubic-PBC rule runs once at the
            # end of restart processing, on the final HARDWALL and DIMENSIONS).
            _od = self.keyword_lookup['DIMENSIONS']
            if any(d < 7 for d in _od):
                raise KeyFileException(
                    'Box size from the restart file (used because RESTART_OVERRIDE_DIMENSIONS '
                    'is True) is too small to correctly support super-long range interactions, '
                    'must be >= 7 in every axis: %s' % (list(_od),))
            
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
                            raise RestartException("\n\nRESIZED_EQUILIBRATION dimensions [%s] are smaller than the restart file's dimensions [%s]. This is an incompatible situation.\n"%(self.keyword_lookup['RESIZED_EQUILIBRATION'][dim], restart_object.dimensions[dim]))

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


        # now that HARDWALL and DIMENSIONS are final (either may have been taken
        # from the restart file), apply the cluster-rotate box rule to the box
        # the run will actually use
        self._check_cluster_rotate_box()

        print("")
        print(".... Restart file processed")
        print("############################")
        print("")

        
