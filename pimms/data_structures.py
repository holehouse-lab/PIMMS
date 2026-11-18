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
from typing import Any, Dict, Iterable, List, Optional, Tuple, Type

from pimms.latticeExceptions import KeyFileException
from pimms import pimmslogger


# What a number in a PIMMS input file may look like. Python's int() and float()
# accept a good deal more than this - digit-group underscores ("1_000"), every
# Unicode decimal digit (full-width, Arabic-Indic, ...), surrounding whitespace of
# any kind - so a typo or a value pasted from a formatted document was read as
# some number rather than refused. re.ASCII is belt and braces: the classes below
# are spelled out and match ASCII only either way.
_ASCII_INTEGER = re.compile(r'[+-]?[0-9]+', re.ASCII)
_ASCII_FLOAT = re.compile(r'[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', re.ASCII)

# The entries of a freeze-file bead (B) line: a chain ID, then bead indices
# counted from 0 or inclusive ranges of them ("0-36"). Only unsigned plain-ASCII
# digits: a sign has no meaning for an index, and "-3" would read as a range
# with its first number missing.
_FREEZE_CHAIN_ID = re.compile(r'[0-9]+', re.ASCII)
_FREEZE_BEAD_RANGE = re.compile(r'([0-9]+)(?:-([0-9]+))?', re.ASCII)


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

# The shortest free stretch held at both ends by frozen beads that single-bead
# crankshaft moves can trap. A lone free bead between two frozen ones is not at
# risk: the crankshaft draws its new site from every site next to both of its
# neighbours, so any position it can hold is one move away. Two or more free
# beads can sit in an arrangement none of them can leave one at a time (two free
# beads closing a 2x2 square with their frozen neighbours, for example).
MIN_TRAPPABLE_FREE_STRETCH = 2


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


def contiguous_runs(indices: Iterable[int]) -> List[Tuple[int, int]]:
    """
    Collapse bead indices into runs of consecutive indices.

    Used to write and report sets of frozen beads compactly, e.g. beads 0 to 36
    and 120 become ``[(0, 36), (120, 120)]``.

    Parameters
    ----------
    indices : iterable of int
        The bead indices, in any order; duplicates are allowed.

    Returns
    -------
    list of tuple of (int, int)
        The inclusive (first, last) runs, in ascending order. Empty if
        ``indices`` is empty.
    """
    runs: List[Tuple[int, int]] = []
    for index in sorted(set(int(i) for i in indices)):
        if runs and index == runs[-1][1] + 1:
            runs[-1] = (runs[-1][0], index)
        else:
            runs.append((index, index))
    return runs


def enclosed_free_stretches(frozen_bead_indices: Iterable[int]) -> List[Tuple[int, int]]:
    """
    The free stretches of a partly frozen chain that are held at both ends.

    A stretch of free beads with a frozen bead on each side is moved by
    single-bead crankshaft moves (and, when it is enabled, by
    ``MOVE_SEGMENT_REGROWTH`` windows), and when it is
    ``MIN_TRAPPABLE_FREE_STRETCH`` beads or longer single-bead moves cannot
    reach every arrangement of it: some arrangements cannot be left at all. A
    free stretch at a chain end has one free end and is not counted here.

    Parameters
    ----------
    frozen_bead_indices : iterable of int
        The frozen beads of one chain, counted from 0, in any order.

    Returns
    -------
    list of tuple of (int, int)
        The inclusive (first, last) free stretches of at least
        ``MIN_TRAPPABLE_FREE_STRETCH`` beads that lie between two frozen beads,
        in ascending order.
    """
    frozen = sorted(set(int(index) for index in frozen_bead_indices))
    stretches = []
    for left, right in zip(frozen, frozen[1:]):
        if right - left - 1 >= MIN_TRAPPABLE_FREE_STRETCH:
            stretches.append((left + 1, right - 1))
    return stretches


def format_bead_runs(indices: Iterable[int], separator: str = ', ') -> str:
    """
    Write bead indices as runs, e.g. ``'0-36, 120'``.

    Parameters
    ----------
    indices : iterable of int
        The bead indices.
    separator : str, optional
        What goes between runs. Default ``', '``; a freeze file uses ``' '``.

    Returns
    -------
    str
        Each run as ``first-last``, or as one number for a run of one bead.
    """
    return separator.join(str(first) if first == last else '%d-%d' % (first, last)
                          for first, last in contiguous_runs(indices))


class FreezeFile:
    """
    Reader/container for a simulation freeze file.

    Parses a freeze file (see :meth:`__init__` for the file format) and exposes
    the chains it freezes whole (``C`` lines) and the beads it freezes one by one
    (``B`` lines). :meth:`resolve` combines the two, given the chain lengths, into
    what the simulation works with: the chains that hold at least one frozen bead,
    and the frozen beads of the chains that are only partly frozen.
    """


    # ...........................................................................
    #
    def __init__(self, filename: str) -> None:
        """
        Read a freeze file, which says which chains, or which beads of a chain,
        are held fixed during the simulation.

                Expected freeze file structure
                ------------------------------
                The freeze file is plain text with one directive per line.

                Supported directives:
                - `C <chain_id> <chain_id> ...`
                    Freeze every bead of one or more chains.
                - `B <chain_id> <bead> <bead> ...`
                    Freeze some beads of one chain. A bead is its index along
                    the chain counted from 0 (as in ANA_RESIDUE_PAIRS), or an
                    inclusive range of indices written `first-last`.

                Notes:
                - Either directive may appear on any number of lines, and a
                  chain or bead named more than once is frozen once.
                - Empty lines are ignored.
                - Lines beginning with `#` are ignored.
                - Inline comments are allowed after `#`.
                - A chain whose every bead is frozen by `B` lines is treated
                  exactly like a chain on a `C` line.
                - The free beads of a partly frozen chain are moved by the
                  crankshaft move (MOVE_CRANKSHAFT) only; every move that
                  displaces a whole chain or a segment of one leaves such a
                  chain alone.

                Example:
                # Freeze two groups of chains
                C 1 2 3
                C 10 11

                # Inline comments are also supported
                C 42  # freeze chain 42

                # Freeze chain 5 except its first 20 beads (0-19)
                B 5 20-134

        Parameters
        ----------
        filename : str
            The name (path) of the freeze file to be read.

        Raises
        ------
        KeyFileException
            If the file does not exist, if it is not a UTF-8 text file, if a line
            has a non-ASCII space outside a comment, if a C directive contains
            something other than integer chain IDs (written in plain ASCII
            digits), if a B directive does not start with a chain ID or holds a
            bead that is neither an index nor an ascending range of indices, if a
            directive carries no IDs at all (or a B directive names no beads), or
            if a line starts with anything other than C or B.

        """

        # check that the file exists
        if not os.path.isfile(filename):
            raise KeyFileException(f'Unable to find FREEZE FILE. Passed filename is: {filename}. Please verify the file actually exists')


        # initialize the chains list and the bead lines. The bead lines keep
        # their line numbers and their ranges unexpanded: whether a range fits
        # its chain is only known once the chain lengths are (see resolve), and
        # expanding a mistyped 0-9999999999 before that check would exhaust
        # memory instead of reporting the typo
        chains = []
        bead_lines = []

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

            # if this line is reporting on beads (B)
            elif sline[0] == 'B':
                bead_lines.append(_parse_bead_line(sline[1:], idx + 1, line))

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
        self._bead_lines = bead_lines
        self._filename = filename

    # ...........................................................................
    #
    @property
    def chains(self) -> List[int]:
        """
        list of int : The deduplicated chain IDs on ``C`` lines, which are
        frozen entirely. A chain that ``B`` lines freeze entirely is not listed
        here; :meth:`resolve` gives the full picture.
        """
        return self._chains

    # ...........................................................................
    #
    @property
    def beads(self) -> List[Tuple[int, int, int]]:
        """
        list of tuple of (int, int, int) : The beads named on ``B`` lines, as
        ``(chainID, first, last)`` runs of consecutive bead indices (inclusive,
        counted from 0), merged across lines and sorted. Empty when the file has
        no ``B`` lines.
        """
        ranges_by_chain: Dict[int, List[Tuple[int, int]]] = {}
        for _, chainID, ranges in self._bead_lines:
            ranges_by_chain.setdefault(chainID, []).extend(ranges)

        merged = []
        for chainID in sorted(ranges_by_chain):
            # merge in range space, never by expanding: the ranges have not been
            # checked against the chain length yet
            for first, last in sorted(ranges_by_chain[chainID]):
                if merged and merged[-1][0] == chainID and first <= merged[-1][2] + 1:
                    merged[-1] = (chainID, merged[-1][1], max(last, merged[-1][2]))
                else:
                    merged.append((chainID, first, last))
        return merged

    # ...........................................................................
    #
    @property
    def filename(self) -> str:
        """
        str : The path of the freeze file that was read.
        """
        return self._filename


    # ...........................................................................
    #
    @staticmethod
    def lattice_chain_lengths(latticeObject: Any) -> Dict[int, int]:
        """
        The number of beads in every chain of a lattice, the input :meth:`resolve` needs.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice. Its ``chains`` dictionary (chainID -> Chain) is read.

        Returns
        -------
        dict of int to int
            chainID -> number of beads.
        """
        return {chainID: len(chain.positions) for chainID, chain in latticeObject.chains.items()}


    # ...........................................................................
    #
    def resolve(self, chain_lengths: Dict[int, int]) -> Tuple[List[int], Dict[int, Tuple[int, ...]]]:
        """
        Combine the C and B lines into the frozen chains and the frozen beads.

        Parameters
        ----------
        chain_lengths : dict of int to int
            The number of beads in each chain, keyed by chainID. Every chain
            named on a ``B`` line must be in it; chains on ``C`` lines are not
            looked up (:meth:`validate_freeze_file` checks those exist).

        Returns
        -------
        tuple of (list of int, dict of int to tuple of int)
            ``frozen_chains``: every chain with at least one frozen bead,
            sorted. No move that displaces a whole chain or a segment of one
            acts on these chains. ``frozen_beads``: for each chain in
            ``frozen_chains`` that is only partly frozen, its frozen bead
            indices (sorted, counted from 0); the crankshaft moves its other
            beads. A chain in ``frozen_chains`` but not in ``frozen_beads`` is
            frozen entirely - it is on a ``C`` line, or its ``B`` lines cover
            every bead.

        Raises
        ------
        KeyFileException
            If a ``B`` line names a chain that is not in ``chain_lengths``, or a
            bead past the end of its chain.
        """
        named_beads: Dict[int, set] = {}
        for line_number, chainID, ranges in self._bead_lines:
            if chainID not in chain_lengths:
                raise KeyFileException(
                    f'Line {line_number} of the freeze file [{self._filename}] freezes beads of chain '
                    f'{chainID}, which is NOT present in the system. The chains present are '
                    f'{sorted(chain_lengths)}.')
            length = int(chain_lengths[chainID])
            for first, last in ranges:
                if last >= length:
                    raise KeyFileException(
                        f'Line {line_number} of the freeze file [{self._filename}] freezes bead {last} of '
                        f'chain {chainID}, which has {length} bead(s), numbered 0 to {length - 1} (bead '
                        f'indices in a freeze file count from 0).')
                named_beads.setdefault(chainID, set()).update(range(first, last + 1))

        entirely_frozen = set(self._chains)
        frozen_beads: Dict[int, Tuple[int, ...]] = {}
        for chainID in sorted(named_beads):
            if chainID in entirely_frozen:
                continue
            if len(named_beads[chainID]) == int(chain_lengths[chainID]):
                entirely_frozen.add(chainID)
            else:
                frozen_beads[chainID] = tuple(sorted(named_beads[chainID]))

        return sorted(entirely_frozen | set(frozen_beads)), frozen_beads


    # ...........................................................................
    #
    def validate_freeze_file(self, latticeObject: Any) -> None:
        """
        Function to validate that the chains and beads specified in the freeze file
        are actually present in the lattice object.

        Parameters
        ----------
        latticeObject : Lattice
            The lattice object to be validated against. Its ``chains``
            dictionary (keyed by chain ID) is read, and the length of each
            chain named on a ``B`` line.

        Returns
        -------
        None
            No return variable, but an exception is raised if the freeze file is not valid

        Raises
        ------
        KeyFileException
            If the freeze file names a chain ID that is not present in the
            lattice object, or a bead past the end of its chain.

        """

        # for each chain specified in the freeze file
        for chainID in self.chains:

            # if the chain is not present in the lattice object, raise an exception
            if chainID not in latticeObject.chains:

                raise KeyFileException(f"\n\nFreeze file {self.filename} specifies chain {chainID}, which is NOT present in the lattice object. Lattice object chains are {list(latticeObject.chains.keys())} while freeze file chains are {self.chains}.")

        # and every bead line: its chain must exist and its beads must lie on it
        if self._bead_lines:
            self.resolve(self.lattice_chain_lengths(latticeObject))


    # ...........................................................................
    #
    def log_freeze_file(self, latticeObject: Optional[Any] = None) -> None:
        """
        Write a summary of the freeze file to the simulation log.

        Logs the freeze file path, the number of frozen chains, and the
        list of frozen chain IDs as STATUS entries. When the file has ``B``
        lines it also logs each partly frozen chain and which of its beads are
        frozen; given the lattice, chains that ``B`` lines freeze entirely are
        counted with the frozen chains.

        Parameters
        ----------
        latticeObject : Lattice, optional
            The lattice the freeze file applies to, which gives the chain
            lengths. Without it, the ``B`` lines are logged as written. Default
            None.

        Returns
        -------
        None
            No return value; entries are appended to the simulation log.

        """
        pimmslogger.log_status(f"Freeze file             : {self.filename}", timestamp=False)

        if not self._bead_lines or latticeObject is None:
            pimmslogger.log_status(f"Number of frozen chains : {len(self.chains)}", timestamp=False)
            pimmslogger.log_status(f"Frozen chainIDs         : {str(self.chains)}", timestamp=False)
            for chainID, first, last in self.beads:
                pimmslogger.log_status(f"Frozen beads of chain {chainID} : {first}-{last}", timestamp=False)
            return

        lengths = self.lattice_chain_lengths(latticeObject)
        (frozen_chains, frozen_beads) = self.resolve(lengths)
        entirely = [chainID for chainID in frozen_chains if chainID not in frozen_beads]
        pimmslogger.log_status(f"Number of frozen chains : {len(entirely)}", timestamp=False)
        pimmslogger.log_status(f"Frozen chainIDs         : {str(entirely)}", timestamp=False)
        pimmslogger.log_status(f"Partly frozen chains    : {len(frozen_beads)} (their free beads move by "
                               f"MOVE_CRANKSHAFT only)", timestamp=False)
        for chainID, beads in frozen_beads.items():
            pimmslogger.log_status(f"  chain {chainID}: {len(beads)} of {lengths[chainID]} beads frozen "
                                   f"({format_bead_runs(beads)})", timestamp=False)


    # ...........................................................................
    #
    def __str__(self) -> str:
        """
        Return a human-readable summary of the frozen chains and beads.

        Returns
        -------
        str
            Multi-line string listing the frozen chain IDs and the bead runs.

        """
        s = 'Freeze file:\n'
        s = s + 'Chains: %s\n'%self.chains
        s = s + 'Beads: %s\n'%self.beads
        return s

    # ...........................................................................
    #
    def __repr__(self) -> str:
        """
        Return the same summary string as :meth:`__str__`.

        Returns
        -------
        str
            Multi-line string listing the frozen chain IDs and the bead runs.

        """
        return self.__str__()


def _parse_bead_line(text: str, line_number: int, line: str) -> Tuple[int, int, Tuple[Tuple[int, int], ...]]:
    """
    Read the entries of one freeze-file ``B`` line.

    Parameters
    ----------
    text : str
        The line after its ``B``, with any comment removed.
    line_number : int
        The line's number in the file (from 1), for the messages.
    line : str
        The whole line as written, for the messages.

    Returns
    -------
    tuple of (int, int, tuple of tuple of (int, int))
        The line number, the chainID and the inclusive (first, last) bead ranges
        in the order written (a single bead ``n`` is the range ``(n, n)``).

    Raises
    ------
    KeyFileException
        If the line names no chain, if its chain ID is not plain ASCII digits,
        if it names no beads, or if a bead is neither an index nor an ascending
        range of indices.
    """
    tokens = split_ascii_whitespace(text)
    if not tokens:
        raise KeyFileException(
            f'Freeze-file bead directive on line {line_number} does not contain a chain ID')

    chain_token, bead_tokens = tokens[0], tokens[1:]
    if _FREEZE_CHAIN_ID.fullmatch(chain_token) is None:
        raise KeyFileException(
            f'Error parsing beads in freeze file on line {line_number}: a B line starts with the ID of '
            f'the chain whose beads it freezes, not {chain_token!r}: {line}')
    if not bead_tokens:
        raise KeyFileException(
            f'Freeze-file bead directive on line {line_number} names chain {chain_token} but no beads: '
            f'list bead indices (counted from 0) or ranges such as 0-36 after the chain ID: {line}')

    ranges = []
    for token in bead_tokens:
        match = _FREEZE_BEAD_RANGE.fullmatch(token)
        if match is None:
            raise KeyFileException(
                f'Error parsing beads in freeze file on line {line_number}: {token!r} is neither a bead '
                f'index (counted from 0) nor a range such as 0-36: {line}')
        first = int(match.group(1))
        last = first if match.group(2) is None else int(match.group(2))
        if last < first:
            raise KeyFileException(
                f'Error parsing beads in freeze file on line {line_number}: the range {token} runs '
                f'backwards (write it as {last}-{first}): {line}')
        ranges.append((first, last))

    return (line_number, int(chain_token), tuple(ranges))
