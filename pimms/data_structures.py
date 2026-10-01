## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................



##
## Small configuration containers: AnalysisSettings (built by the Simulation
## from ANA_CLUSTER_THRESHOLD) and FreezeFile (the parsed FREEZE_FILE, built by
## the keyfile parser).
##
## The module also holds the two things the three text-input parsers (keyfile,
## parameter file, freeze file) share: read_text_input, which reads a file as
## UTF-8 and turns a file in any other encoding into a clean PIMMS exception,
## and the plain-ASCII number checks is_ascii_integer / is_ascii_float.
##

import os
import re
import unicodedata
from typing import List, Optional, Type

from pimms.latticeExceptions import UnfinishedCodeException, KeyFileException
from pimms import pimmslogger


# What a number in a PIMMS input file may look like. Python's int() and float()
# accept a good deal more than this - digit-group underscores ("1_000"), every
# Unicode decimal digit (full-width, Arabic-Indic, ...), surrounding whitespace of
# any kind - so a typo or a value pasted from a formatted document was read as
# some number rather than refused. re.ASCII is belt and braces: the classes below
# are spelled out and match ASCII only either way.
_ASCII_INTEGER = re.compile(r'[+-]?[0-9]+', re.ASCII)
_ASCII_FLOAT = re.compile(r'[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', re.ASCII)


# The whitespace that separates columns and pads values in the three text inputs:
# exactly the ASCII characters str.split() treats as whitespace. str.split() and
# str.strip() with no argument also accept every Unicode space, which is how a
# no-break space pasted from a web page or a PDF used to act as an invisible
# column separator.
ASCII_WHITESPACE = ' \t\n\r\x0b\x0c\x1c\x1d\x1e\x1f'
_ASCII_WHITESPACE_RUN = re.compile('[' + re.escape(ASCII_WHITESPACE) + ']+')

# characters that take no width at all and are not whitespace to Python, but
# separate or hide text all the same: a byte-order mark that is not at the very
# start of the file (two files concatenated, say) and the zero-width space
_INVISIBLE_CHARACTERS = ('\ufeff', '\u200b')


def split_ascii_whitespace(text: str) -> List[str]:
    """
    Split a line into columns on ASCII whitespace only.

    This is ``str.split()`` restricted to the ASCII whitespace characters
    (space, tab, newline, carriage return, vertical tab, form feed and the four
    ASCII separator controls), so a no-break space or any other non-ASCII space
    stays inside a column, where the caller can see and refuse it.

    Parameters
    ----------
    text : str
        The text to split.

    Returns
    -------
    list of str
        The non-empty columns, in order.
    """
    return [column for column in _ASCII_WHITESPACE_RUN.split(text) if column]


def find_non_ascii_space(text: str) -> Optional[str]:
    """
    Find the first non-ASCII whitespace (or invisible) character in a piece of text.

    Such a character looks like an ordinary space, or like nothing at all, so a
    message that only echoes the line is no help; the code point and its Unicode
    name say what to delete.

    Parameters
    ----------
    text : str
        The text to search (a line of an input file, with its comment removed).

    Returns
    -------
    str or None
        A description such as ``'U+00A0 NO-BREAK SPACE'``, or None if the text
        has no such character.
    """
    for character in text:
        if ord(character) > 127 and (character.isspace() or character in _INVISIBLE_CHARACTERS):
            return 'U+%04X %s' % (ord(character), unicodedata.name(character, 'unnamed character'))
    return None


def is_ascii_integer(token: str) -> bool:
    """
    Report whether a token is an integer written in plain ASCII.

    That is an optional sign followed by one or more of the digits 0-9 and
    nothing else: no digit-group underscore, no decimal point, no exponent, no
    surrounding whitespace and no non-ASCII digit.

    Parameters
    ----------
    token : str
        The token to check.

    Returns
    -------
    bool
        True if ``token`` is a string of that form.
    """
    return isinstance(token, str) and _ASCII_INTEGER.fullmatch(token) is not None


def is_ascii_float(token: str) -> bool:
    """
    Report whether a token is a decimal number written in plain ASCII.

    Accepted are an optional sign, digits with an optional decimal point (``10``,
    ``10.``, ``10.5``, ``.5``) and an optional exponent (``1e5``, ``2.5E-3``).
    The words Python's ``float`` also takes (``nan``, ``inf``, ``infinity``),
    digit-group underscores, surrounding whitespace and non-ASCII digits are not.

    Parameters
    ----------
    token : str
        The token to check.

    Returns
    -------
    bool
        True if ``token`` is a string of that form.
    """
    return isinstance(token, str) and _ASCII_FLOAT.fullmatch(token) is not None


def read_text_input(filename: str, description: str, exception_class: Type[Exception]) -> List[str]:
    """
    Read a PIMMS text input file (keyfile, parameter file or freeze file) as UTF-8.

    The file is decoded as UTF-8 whatever the platform's default encoding is, so
    the same file reads the same way everywhere, and a leading byte-order mark
    (which some Windows editors add) is dropped rather than glued onto the first
    keyword. A file that is not UTF-8 - one saved as Latin-1 / Windows-1252 with
    an accented character or an Angstrom sign in a comment, a UTF-16 file, a
    binary file - is refused with an exception that names the file and where the
    first offending byte is, in place of a bare ``UnicodeDecodeError``.

    Parameters
    ----------
    filename : str
        Path of the file to read.

    description : str
        What the file is, for the error message (``'keyfile'``,
        ``'parameter file'``, ``'freeze file'``).

    exception_class : type
        The exception to raise if the file is not UTF-8 (``KeyFileException`` or
        ``ParameterFileException``).

    Returns
    -------
    list of str
        The lines of the file, each with its line ending, as ``readlines`` gives
        them.

    Raises
    ------
    exception_class
        If the file cannot be decoded as UTF-8.
    """
    try:
        with open(filename, 'r', encoding='utf-8-sig') as fh:
            return fh.readlines()
    except UnicodeDecodeError:
        pass

    # the offset a streaming decoder reports is relative to its read buffer, so
    # decode the raw bytes in one piece to find the offset in the file
    with open(filename, 'rb') as fh:
        raw = fh.read()

    offset = 0
    try:
        raw.decode('utf-8')
    except UnicodeDecodeError as error:
        offset = error.start

    line_number = raw.count(b'\n', 0, offset) + 1
    # UTF-32 first: its little-endian mark begins with the UTF-16 one
    if raw[:4] in (b'\xff\xfe\x00\x00', b'\x00\x00\xfe\xff'):
        hint = ('The file starts with a UTF-32 byte-order mark, so it was saved as UTF-32; save it '
                'again as UTF-8.')
    elif raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        hint = ('The file starts with a UTF-16 byte-order mark, so it was probably saved as '
                '"Unicode" / UTF-16; save it again as UTF-8.')
    else:
        hint = ('It was probably saved in a legacy encoding such as Latin-1 / Windows-1252 with a '
                'non-ASCII character (an accented letter, an Angstrom or micro sign) in it; save it '
                'again as UTF-8, or remove the character.')
    raise exception_class(
        'The %s [%s] is not a UTF-8 text file: byte 0x%s at byte offset %d (line %d) cannot be '
        'decoded. %s' % (description, filename, raw[offset:offset + 1].hex(), offset, line_number, hint))

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
                C 1 2 3
                C 10 11

                # Inline comments are also supported
                C 42  # freeze chain 42

        Parameters
        ----------
        filename : str
            The name (path) of the freeze file to be read.

        Raises
        ------
        KeyFileException
            If the file does not exist, if it is not a UTF-8 text file, if a line
            has a non-ASCII space outside a comment, if a C or B directive
            contains something other than integer IDs (written in plain ASCII
            digits), if a directive carries no IDs at all, or if a line starts
            with anything other than C or B.

        UnfinishedCodeException
            If a bead-level (``B``) directive is encountered, as bead freezing is
            not yet implemented.

        """

        # check that the file exists
        if not os.path.isfile(filename):
            raise KeyFileException(f'Unable to find FREEZE FILE. Passed filename is: {filename}. Please verify the file actually exists')

                    
        # initialize the chains and beads lists
        chains = []
        beads = []

        # read the contents (as UTF-8, with any byte-order mark dropped; a file
        # in another encoding is refused with a message naming the file)
        content = read_text_input(filename, 'freeze file', KeyFileException)


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

            # a non-ASCII space outside a comment is refused by name: it used to
            # act as an invisible separator between chain IDs
            culprit = find_non_ascii_space(line.split('#')[0])
            if culprit is not None:
                raise KeyFileException(
                    f'Line {idx + 1} of the freeze file [{filename}] contains the non-ASCII whitespace '
                    f'character {culprit} outside a comment. Chain IDs are separated by plain spaces or '
                    f'tabs only - replace the character (it usually arrives by pasting from a web page, a '
                    f'PDF or a word processor): {line}')

            # discard comments at the end of the line as well, if they exist
            sline = sline.split('#')[0].strip()


            # if this line is reporting on chains (C)
            if sline[0] == 'C':

                # plain ASCII integers only: int() alone would also take "1_0"
                # or full-width digits and freeze a chain the file does not name
                tokens = split_ascii_whitespace(sline[1:])
                if not all(is_ascii_integer(i) for i in tokens):
                    raise KeyFileException(
                        f'Error parsing chains in freeze file on line {idx + 1}: {line}'
                    )
                local_chains = [int(i) for i in tokens]

                if not local_chains:
                    raise KeyFileException(
                        f'Freeze-file chain directive on line {idx + 1} does not '
                        'contain a chain ID')

                chains.extend(local_chains)

            # if this line is reporting on beads [NOT YET IMPLEMENTED]
            elif sline[0] == 'B':

                # split() with no argument (as the 'C' branch above does), NOT
                # split(' '): the latter keeps the empty field produced by the
                # space right after the 'B', so int('') blew up on well-formed
                # input and reported it as a parse error
                tokens = split_ascii_whitespace(sline[1:])
                if not all(is_ascii_integer(i) for i in tokens):
                    raise KeyFileException(
                        f'Error parsing beads in freeze file on line {idx + 1}: {line}'
                    )
                local_beads = [int(i) for i in tokens]

                if not local_beads:
                    raise KeyFileException(
                        f'Freeze-file bead directive on line {idx + 1} does not '
                        'contain a bead ID')

                beads.extend(local_beads)
                raise UnfinishedCodeException('Beads not yet implemented for freezing')

            else:
                # anything else is an error, never silently dropped: a typo'd
                # freeze file (e.g. lowercase 'c') previously ran the whole
                # simulation with NOTHING frozen and no hint anything was wrong
                raise KeyFileException(
                    f'Unrecognised directive on line {idx + 1} of freeze file: '
                    f'{sline!r} (expected a line starting with C or B)')

        # remove duplicates
        # Stable ordering makes logs, summaries and frozen-index construction
        # deterministic even when the input contains duplicate or unsorted IDs.
        self._chains = sorted(set(chains))
        self._beads = sorted(set(beads))
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
            The lattice object to be validated against. Only its ``chains``
            dictionary (keyed by chain ID) is read.

        Returns
        -------
        None
            No return variable, but an exception is raised if the freeze file is not valid

        Raises
        ------
        KeyFileException
            If the freeze file names a chain ID that is not present in the
            lattice object.

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
    
