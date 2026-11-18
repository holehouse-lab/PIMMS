## ...........................................................................
##
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

##
## structure_reader
##
## Reads the atom records of a PDB or PDBx/mmCIF structure file into one
## format-independent list of atoms. It exists for pdb2pimms, which turns a
## structure into a PIMMS restart file, so it reads only what that conversion
## needs: the coordinates, the names that identify each atom and residue,
## alternate locations, occupancies, elements, model numbers, and the parent
## residue that the file itself gives for each modified residue (PDB MODRES
## records, mmCIF _pdbx_struct_mod_residue). Everything else in the file is
## skipped. Both formats may be gzip-compressed.
##

import gzip
import math
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional, TextIO

from .latticeExceptions import StructureFileException


# File name endings that name the format outright (checked after any .gz is removed)
PDB_SUFFIXES = ('.pdb', '.ent')
CIF_SUFFIXES = ('.cif', '.mmcif')
GZIP_SUFFIX = '.gz'

# The two PDB records that hold atoms (record names are the first six columns)
PDB_ATOM_RECORDS = ('ATOM  ', 'HETATM')

# A PDB atom line must reach the end of the z coordinate (column 54) to be read
PDB_MIN_ATOM_LINE_LENGTH = 54

# Unquoted mmCIF values meaning "unknown" (?) and "not applicable" (.)
CIF_MISSING_VALUES = ('?', '.')

# The mmCIF categories the reader uses
CIF_ATOM_SITE = '_atom_site.'
CIF_MOD_RESIDUE = '_pdbx_struct_mod_residue.'
CIF_ENTRY_ID = '_entry.id'

# The occupancy assumed for an atom whose file leaves the column empty
DEFAULT_OCCUPANCY = 1.0

# The model number given to every atom of a file with no MODEL records
DEFAULT_MODEL = 1


@dataclass(frozen=True, slots=True)
class AtomRecord:
    """
    One atom as read from a structure file, in the same form for both formats.

    Attributes
    ----------
    record : str
        'ATOM' or 'HETATM'.
    name : str
        Atom name, e.g. ``'CA'`` or ``"C4'"``.
    altloc : str
        Alternate location identifier, ``''`` when the atom has none.
    resname : str
        Residue (chemical component) name, e.g. ``'ALA'`` or ``'G'``.
    chain : str
        Author chain identifier (in a PDB file with a blank chain column, the
        segment identifier instead).
    resseq : str
        Author residue number as written in the file.
    icode : str
        Insertion code, ``''`` when the residue has none.
    x, y, z : float
        Cartesian coordinates in Angstroms.
    occupancy : float
        Occupancy, 1.0 when the file gives none.
    element : str
        Upper-case element symbol. When the file gives none, the first letter
        of the atom name (after any leading digits), which is enough to
        recognise a hydrogen.
    model : int
        Model number the atom belongs to.
    """
    record: str
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: str
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    element: str
    model: int


@dataclass
class StructureFile:
    """
    The atoms of a structure file, with what the file says about them.

    Attributes
    ----------
    path : Path
        The file that was read.
    file_format : str
        'PDB' or 'mmCIF'.
    atoms : list of AtomRecord
        Every ATOM and HETATM record of the first data block, in file order.
    models : list of int
        The model numbers present, in the order they first appear.
    parent_residues : dict of str to str
        Modified residue name -> the standard residue it derives from, as the
        file gives it (PDB MODRES records, mmCIF _pdbx_struct_mod_residue).
    entry_id : str
        The PDB identifier from the HEADER record or _entry.id, ``''`` if absent.
    chains_from_segment_ids : bool
        True if any atom's chain came from the PDB segment identifier because
        its chain column was blank (as some simulation packages write).
    """
    path: Path
    file_format: str
    atoms: list[AtomRecord]
    models: list[int]
    parent_residues: dict[str, str] = field(default_factory=dict)
    entry_id: str = ''
    chains_from_segment_ids: bool = False


# .............................................................................
# Format detection and file access
# .............................................................................

def _format_from_name(path: Path) -> Optional[str]:
    """
    Name the format a file's extension implies, ignoring a trailing ``.gz``.

    Parameters
    ----------
    path : Path
        The structure file.

    Returns
    -------
    str or None
        'PDB', 'mmCIF', or None when the extension names neither.
    """
    name = path.name.lower()
    if name.endswith(GZIP_SUFFIX):
        name = name[:-len(GZIP_SUFFIX)]
    if name.endswith(CIF_SUFFIXES):
        return 'mmCIF'
    if name.endswith(PDB_SUFFIXES):
        return 'PDB'
    return None


def _open_text(path: Path) -> TextIO:
    """
    Open a structure file for reading as text, decompressing it if it is gzipped.

    Atom records are plain ASCII in both formats. Other records (REMARKs, author
    names) are occasionally written in another encoding, so undecodable bytes
    are replaced rather than refused: they can never be part of a field this
    reader uses. A UTF-8 byte-order mark is removed (it would otherwise become
    part of the first record and hide it).

    Parameters
    ----------
    path : Path
        The structure file.

    Returns
    -------
    TextIO
        An open text stream; the caller closes it.
    """
    if path.name.lower().endswith(GZIP_SUFFIX):
        return gzip.open(path, 'rt', encoding='utf-8-sig', errors='replace')
    return open(path, 'r', encoding='utf-8-sig', errors='replace')


def _sniff_format(path: Path) -> str:
    """
    Name the format of a file from its content, for a file whose extension names none.

    The file is read until a decisive line: an mmCIF ``data_`` block header,
    ``loop_`` or data name, or a PDB ATOM, HETATM or MODEL record. Header
    records of either format before it are skipped.

    Parameters
    ----------
    path : Path
        The structure file.

    Returns
    -------
    str
        'PDB' or 'mmCIF'.

    Raises
    ------
    StructureFileException
        If the content looks like neither format.
    """
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            lower = stripped.lower()
            if lower.startswith('data_') or lower == 'loop_' or stripped.startswith('_'):
                return 'mmCIF'
            if line.startswith(PDB_ATOM_RECORDS + ('MODEL ',)):
                return 'PDB'
    raise StructureFileException(
        "Cannot tell whether %s is a PDB or an mmCIF file: its name does not end in .pdb, .ent, "
        ".cif or .mmcif (optionally followed by .gz) and it holds neither a PDB atom record nor an "
        "mmCIF data_ block" % path)


def read_structure(path: str | Path) -> StructureFile:
    """
    Read every atom of a PDB or mmCIF file (optionally gzip-compressed).

    The format comes from the extension (``.pdb``/``.ent`` or ``.cif``/``.mmcif``,
    each optionally followed by ``.gz``) or, failing that, from the content.

    Parameters
    ----------
    path : str or Path
        The structure file.

    Returns
    -------
    StructureFile
        The atoms and the file-level information that goes with them.

    Raises
    ------
    StructureFileException
        If the file cannot be opened, is malformed, or holds no atoms.
    """
    path = Path(path)
    if not path.is_file():
        raise StructureFileException("Cannot read structure file %s: no such file" % path)

    try:
        # sniffing reads the file too, so a corrupt gzip file named only .gz fails here
        file_format = _format_from_name(path) or _sniff_format(path)
        with _open_text(path) as handle:
            if file_format == 'PDB':
                structure = _read_pdb(handle, path)
            else:
                structure = _read_cif(handle, path)
    except (OSError, EOFError, zlib.error) as error:
        # EOFError and zlib.error are what a truncated or corrupt gzip file raises
        raise StructureFileException("Cannot read structure file %s: %s" % (path, error))

    if not structure.atoms:
        raise StructureFileException(
            "%s (read as %s) holds no ATOM or HETATM records, so there is nothing to convert"
            % (path, file_format))
    return structure


def _models_in_order(atoms: list[AtomRecord]) -> list[int]:
    """
    List the model numbers of a set of atoms in the order they first appear.

    Parameters
    ----------
    atoms : list of AtomRecord
        The atoms, in file order.

    Returns
    -------
    list of int
        Each model number once.
    """
    return list(dict.fromkeys(atom.model for atom in atoms))


# .............................................................................
# PDB format
# .............................................................................

def _infer_element(atom_name: str) -> str:
    """
    Guess an element symbol from an atom name, for an atom whose file gives no element.

    That is a PDB line with no element column, or an mmCIF atom with no
    ``_atom_site.type_symbol``. It is only used to recognise hydrogens, water's
    oxygen and a CA-only model's carbon, so the first letter after any leading
    digits (``1HB`` is a hydrogen) is enough; it cannot tell a two-letter
    element from a one-letter one (a calcium ion named ``CA`` reads as a
    carbon).

    Parameters
    ----------
    atom_name : str
        The atom name, stripped.

    Returns
    -------
    str
        The guessed element symbol (one upper-case letter), or ``''``.
    """
    letters = atom_name.lstrip('0123456789')
    return letters[:1].upper()


def _pdb_float(text: str, what: str, line_number: int, path: Path) -> float:
    """
    Read one numeric PDB field.

    Parameters
    ----------
    text : str
        The field's columns.
    what : str
        What the field is, for the error message.
    line_number : int
        One-based line number, for the error message.
    path : Path
        The file being read, for the error message.

    Returns
    -------
    float
        The value.

    Raises
    ------
    StructureFileException
        If the field is not a finite number.
    """
    try:
        value = float(text)
    except ValueError:
        raise StructureFileException(
            "Line %d of %s: the %s field %r is not a number" % (line_number, path, what, text))
    if not math.isfinite(value):
        raise StructureFileException(
            "Line %d of %s: the %s field is %r, which is not a finite number" % (line_number, path, what, text))
    return value


def _pdb_atom(line: str, model: int, line_number: int, path: Path) -> tuple[AtomRecord, bool]:
    """
    Read one ATOM or HETATM line using the PDB fixed columns.

    Residue names are read from columns 18-21 rather than 18-20: column 21 is
    blank in a standard file, and some simulation packages write four-character
    names there (``NALA``, ``TIP3``).

    Parameters
    ----------
    line : str
        The record, without its line ending.
    model : int
        The model the record belongs to.
    line_number : int
        One-based line number, for error messages.
    path : Path
        The file being read, for error messages.

    Returns
    -------
    tuple of (AtomRecord, bool)
        The atom, and whether its chain came from the segment identifier
        (columns 73-76) because the chain column was blank.

    Raises
    ------
    StructureFileException
        If the line is too short to hold coordinates or a coordinate is not a number.
    """
    if len(line) < PDB_MIN_ATOM_LINE_LENGTH:
        raise StructureFileException(
            "Line %d of %s is an atom record too short to hold its coordinates (%d columns, "
            "at least %d needed)" % (line_number, path, len(line), PDB_MIN_ATOM_LINE_LENGTH))

    name = line[12:16].strip()
    chain = line[21:22].strip()
    from_segment = False
    if not chain:
        chain = line[72:76].strip()
        from_segment = bool(chain)

    occupancy_text = line[54:60].strip()
    occupancy = (_pdb_float(occupancy_text, 'occupancy', line_number, path)
                 if occupancy_text else DEFAULT_OCCUPANCY)
    element = line[76:78].strip().upper() or _infer_element(name)

    atom = AtomRecord(
        record=line[0:6].strip(),
        name=name,
        altloc=line[16:17].strip(),
        resname=line[17:21].strip(),
        chain=chain,
        resseq=line[22:26].strip(),
        icode=line[26:27].strip(),
        x=_pdb_float(line[30:38], 'x coordinate', line_number, path),
        y=_pdb_float(line[38:46], 'y coordinate', line_number, path),
        z=_pdb_float(line[46:54], 'z coordinate', line_number, path),
        occupancy=occupancy,
        element=element,
        model=model)
    return atom, from_segment


def _pdb_model_number(line: str, line_number: int, path: Path) -> int:
    """
    Read the serial number of a MODEL record.

    Parameters
    ----------
    line : str
        The MODEL record.
    line_number : int
        One-based line number, for the error message.
    path : Path
        The file being read, for the error message.

    Returns
    -------
    int
        The model number.

    Raises
    ------
    StructureFileException
        If the record has no integer serial number.
    """
    fields = line[6:].split()
    try:
        return int(fields[0])
    except (IndexError, ValueError):
        raise StructureFileException(
            "Line %d of %s: the MODEL record has no model number: %r" % (line_number, path, line))


def _read_pdb(handle: TextIO, path: Path) -> StructureFile:
    """
    Read the atoms of a PDB-format file.

    Parameters
    ----------
    handle : TextIO
        The open file.
    path : Path
        The file being read (recorded, and used in error messages).

    Returns
    -------
    StructureFile
        The atoms, model numbers, MODRES parents and entry identifier.

    Raises
    ------
    StructureFileException
        If a record is malformed, or a model number is used for a second block
        of atoms (two MODEL records with the same number, as in two
        single-model files joined together).

    Notes
    -----
    Some programs (VMD, for one) write a trajectory as frames separated by END
    or ENDMDL with no MODEL records. Atoms that follow an END or ENDMDL with no
    MODEL record in between therefore start a new model, numbered one more than
    the last, so that the frames are not read as one doubled structure. For the
    same reason a model number cannot start a second block of atoms.
    """
    atoms: list[AtomRecord] = []
    parents: dict[str, str] = {}
    entry_id = ''
    any_segment_chain = False
    model = DEFAULT_MODEL
    # True once an END or ENDMDL has closed the current model
    model_closed = False
    # every model number that has atoms, and whether the current block of atoms
    # (since the last MODEL, END or ENDMDL) has begun
    models_with_atoms: set[int] = set()
    block_started = False

    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.rstrip('\r\n')
        record = line[0:6]

        if record in PDB_ATOM_RECORDS:
            if model_closed:
                model += 1
                model_closed = False
            if not block_started:
                if model in models_with_atoms:
                    raise StructureFileException(
                        "Line %d of %s starts a second block of atoms for model %d (two MODEL records with the "
                        "same number, as in two single-model files joined together), so the two cannot be told "
                        "apart; give each model its own number" % (line_number, path, model))
                models_with_atoms.add(model)
                block_started = True
            atom, from_segment = _pdb_atom(line, model, line_number, path)
            atoms.append(atom)
            any_segment_chain = any_segment_chain or from_segment
        elif record == 'MODEL ':
            model = _pdb_model_number(line, line_number, path)
            model_closed = False
            block_started = False
        elif record.rstrip() in ('END', 'ENDMDL') and atoms:
            model_closed = True
            block_started = False
        elif record == 'MODRES':
            modified = line[12:15].strip()
            parent = line[24:27].strip()
            if modified and parent:
                parents[modified] = parent
        elif record == 'HEADER':
            entry_id = line[62:66].strip()

    return StructureFile(path=path, file_format='PDB', atoms=atoms, models=_models_in_order(atoms),
                         parent_residues=parents, entry_id=entry_id,
                         chains_from_segment_ids=any_segment_chain)


# .............................................................................
# mmCIF format
# .............................................................................

# A token is its text plus whether it was quoted: a quoted '?' is the
# character ?, an unquoted one means "unknown", and only unquoted text can be
# a data name, loop_ or a block header.
CifToken = tuple[str, bool]


def _split_cif_line(line: str, line_number: int, path: Path) -> Iterator[CifToken]:
    """
    Split one line of an mmCIF file into tokens.

    Tokens are separated by whitespace. A token that starts with a quote runs to
    the next matching quote that is followed by whitespace or the end of the
    line, so ``"C4'"`` is the atom name C4'. A ``#`` at the start of a token
    begins a comment that runs to the end of the line.

    Parameters
    ----------
    line : str
        The line, without its line ending.
    line_number : int
        One-based line number, for the error message.
    path : Path
        The file being read, for the error message.

    Yields
    ------
    tuple of (str, bool)
        Each token's text and whether it was quoted.

    Raises
    ------
    StructureFileException
        If a quoted token is not closed on its line.
    """
    position = 0
    length = len(line)
    while position < length:
        character = line[position]
        if character.isspace():
            position += 1
            continue
        if character == '#':
            return
        if character in ('"', "'"):
            end = position + 1
            while True:
                end = line.find(character, end)
                if end == -1:
                    raise StructureFileException(
                        "Line %d of %s: a quoted value is not closed on its line: %r"
                        % (line_number, path, line))
                if end + 1 == length or line[end + 1].isspace():
                    break
                end += 1
            yield line[position + 1:end], True
            position = end + 1
            continue
        end = position
        while end < length and not line[end].isspace():
            end += 1
        yield line[position:end], False
        position = end


def _cif_tokens(handle: TextIO, path: Path) -> Iterator[CifToken]:
    """
    Turn an mmCIF file into a stream of tokens.

    Besides the tokens of ordinary lines, a line that starts with ``;`` opens a
    text field that runs, newlines included, to the next line that starts with
    ``;``; the whole field is one quoted token.

    Parameters
    ----------
    handle : TextIO
        The open file.
    path : Path
        The file being read, for error messages.

    Yields
    ------
    tuple of (str, bool)
        Each token's text and whether it was quoted.

    Raises
    ------
    StructureFileException
        If a text field is never closed, or a quoted value is not closed on its line.
    """
    text_field: Optional[list[str]] = None
    opened_at = 0
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.rstrip('\r\n')
        if text_field is not None:
            if line.startswith(';'):
                yield '\n'.join(text_field), True
                text_field = None
                yield from _split_cif_line(line[1:], line_number, path)
            else:
                text_field.append(line)
            continue
        if line.startswith(';'):
            text_field = [line[1:]]
            opened_at = line_number
            continue
        yield from _split_cif_line(line, line_number, path)

    if text_field is not None:
        raise StructureFileException(
            "%s: the text field opened by the ';' on line %d is never closed" % (path, opened_at))


class _CifTokenStream:
    """
    A token stream with one token of look-ahead, which loop parsing needs.

    Parameters
    ----------
    tokens : iterator of tuple of (str, bool)
        The tokens, from _cif_tokens().
    """

    def __init__(self, tokens: Iterator[CifToken]) -> None:
        """
        Wrap a token iterator.

        Parameters
        ----------
        tokens : iterator of tuple of (str, bool)
            The tokens to read from.

        Returns
        -------
        None
        """
        self._tokens = tokens
        self._next: Optional[CifToken] = next(tokens, None)

    def peek(self) -> Optional[CifToken]:
        """
        Look at the next token without consuming it.

        Returns
        -------
        tuple of (str, bool) or None
            The next token, or None at the end of the file.
        """
        return self._next

    def take(self) -> Optional[CifToken]:
        """
        Consume and return the next token.

        Returns
        -------
        tuple of (str, bool) or None
            The token, or None at the end of the file.
        """
        token = self._next
        self._next = next(self._tokens, None)
        return token


def _starts_new_item(token: CifToken) -> bool:
    """
    Say whether a token ends the values of a loop: a data name, loop_, or a block or frame header.

    Parameters
    ----------
    token : tuple of (str, bool)
        The token.

    Returns
    -------
    bool
        True if the token begins something other than a loop value.
    """
    text, quoted = token
    if quoted:
        return False
    lower = text.lower()
    return (lower.startswith('_') or lower == 'loop_' or lower.startswith('data_')
            or lower.startswith('save_') or lower in ('global_', 'stop_'))


def _cif_value(token: CifToken) -> str:
    """
    Turn a value token into a string, mapping the unquoted missing-value markers to ``''``.

    Parameters
    ----------
    token : tuple of (str, bool)
        The value token.

    Returns
    -------
    str
        The value, ``''`` for an unquoted ``?`` or ``.``.
    """
    text, quoted = token
    if not quoted and text in CIF_MISSING_VALUES:
        return ''
    return text


def _read_loop(stream: _CifTokenStream, path: Path) -> tuple[list[str], list[list[str]]]:
    """
    Read one loop_: its data names, then its values, row by row.

    Parameters
    ----------
    stream : _CifTokenStream
        The token stream, positioned just after ``loop_``.
    path : Path
        The file being read, for the error message.

    Returns
    -------
    tuple of (list of str, list of list of str)
        The lower-cased data names and the rows of values.

    Raises
    ------
    StructureFileException
        If the number of values is not a multiple of the number of data names.
    """
    tags: list[str] = []
    while True:
        token = stream.peek()
        if token is None or token[1] or not token[0].startswith('_'):
            break
        stream.take()
        tags.append(token[0].lower())

    values: list[str] = []
    while True:
        token = stream.peek()
        if token is None or _starts_new_item(token):
            break
        stream.take()
        values.append(_cif_value(token))

    if not tags:
        raise StructureFileException("%s: a loop_ has no data names" % path)
    if len(values) % len(tags) != 0:
        raise StructureFileException(
            "%s: the loop that starts with %s has %d values, which is not a whole number of rows "
            "of %d" % (path, tags[0], len(values), len(tags)))
    width = len(tags)
    rows = [values[start:start + width] for start in range(0, len(values), width)]
    return tags, rows


def _column(tags: list[str], *names: str) -> Optional[int]:
    """
    Find the first of several data names present in a loop.

    Parameters
    ----------
    tags : list of str
        The loop's lower-cased data names.
    *names : str
        Candidate data names, most preferred first.

    Returns
    -------
    int or None
        The column index of the first name present, or None if none is.
    """
    for name in names:
        if name.lower() in tags:
            return tags.index(name.lower())
    return None


def _first_present(row: list[str], *columns: Optional[int]) -> str:
    """
    Return the first non-missing value of a row among several columns.

    mmCIF gives most identifiers twice (``auth_*`` as the depositors wrote
    them, ``label_*`` as the archive assigned them). The author versions are the
    ones that match a PDB-format file and the literature, so they come first,
    with the ``label_*`` value used where the author value is missing.

    Parameters
    ----------
    row : list of str
        One row of the loop.
    *columns : int or None
        Column indices in order of preference; None entries are skipped.

    Returns
    -------
    str
        The first non-empty value, or ``''``.
    """
    for column in columns:
        if column is not None and row[column]:
            return row[column]
    return ''


def _cif_float(row: list[str], column: int, what: str, row_number: int, path: Path) -> float:
    """
    Read one numeric _atom_site value.

    Parameters
    ----------
    row : list of str
        The atom's row.
    column : int
        The column to read.
    what : str
        The column's data name, for the error message.
    row_number : int
        One-based row number within _atom_site, for the error message.
    path : Path
        The file being read, for the error message.

    Returns
    -------
    float
        The value.

    Raises
    ------
    StructureFileException
        If the value is missing or not a finite number.
    """
    try:
        value = float(row[column])
    except ValueError:
        raise StructureFileException(
            "%s: atom %d of _atom_site has %s = %r, which is not a number"
            % (path, row_number, what, row[column]))
    if not math.isfinite(value):
        raise StructureFileException(
            "%s: atom %d of _atom_site has %s = %r, which is not a finite number"
            % (path, row_number, what, row[column]))
    return value


def _cif_atoms(tags: list[str], rows: list[list[str]], path: Path) -> list[AtomRecord]:
    """
    Build atoms from the rows of the _atom_site loop.

    Parameters
    ----------
    tags : list of str
        The loop's lower-cased data names.
    rows : list of list of str
        The loop's rows.
    path : Path
        The file being read, for error messages.

    Returns
    -------
    list of AtomRecord
        One atom per row.

    Raises
    ------
    StructureFileException
        If a required column is absent or a coordinate is not a number.
    """
    def required(*names: str) -> int:
        """
        Find a column the conversion cannot do without.

        Parameters
        ----------
        *names : str
            Candidate data names, most preferred first.

        Returns
        -------
        int
            The column index.

        Raises
        ------
        StructureFileException
            If none of the names is present.
        """
        index = _column(tags, *names)
        if index is None:
            raise StructureFileException(
                "%s: the _atom_site loop has no %s column" % (path, ' or '.join(names)))
        return index

    x_col = required('_atom_site.cartn_x')
    y_col = required('_atom_site.cartn_y')
    z_col = required('_atom_site.cartn_z')
    required('_atom_site.auth_atom_id', '_atom_site.label_atom_id')
    required('_atom_site.auth_comp_id', '_atom_site.label_comp_id')
    required('_atom_site.auth_asym_id', '_atom_site.label_asym_id')
    required('_atom_site.auth_seq_id', '_atom_site.label_seq_id')

    group_col = _column(tags, '_atom_site.group_pdb')
    auth_atom = _column(tags, '_atom_site.auth_atom_id')
    label_atom = _column(tags, '_atom_site.label_atom_id')
    auth_comp = _column(tags, '_atom_site.auth_comp_id')
    label_comp = _column(tags, '_atom_site.label_comp_id')
    auth_asym = _column(tags, '_atom_site.auth_asym_id')
    label_asym = _column(tags, '_atom_site.label_asym_id')
    auth_seq = _column(tags, '_atom_site.auth_seq_id')
    label_seq = _column(tags, '_atom_site.label_seq_id')
    altloc_col = _column(tags, '_atom_site.label_alt_id')
    icode_col = _column(tags, '_atom_site.pdbx_pdb_ins_code')
    occupancy_col = _column(tags, '_atom_site.occupancy')
    element_col = _column(tags, '_atom_site.type_symbol')
    model_col = _column(tags, '_atom_site.pdbx_pdb_model_num')

    atoms = []
    for row_number, row in enumerate(rows, start=1):
        name = _first_present(row, auth_atom, label_atom)
        # as for a PDB line with no element column, so that both formats find the same hydrogens
        element = (row[element_col].upper() if element_col is not None else '') or _infer_element(name)
        occupancy = DEFAULT_OCCUPANCY
        if occupancy_col is not None and row[occupancy_col]:
            occupancy = _cif_float(row, occupancy_col, '_atom_site.occupancy', row_number, path)
        model = DEFAULT_MODEL
        if model_col is not None and row[model_col]:
            try:
                model = int(row[model_col])
            except ValueError:
                raise StructureFileException(
                    "%s: atom %d of _atom_site has pdbx_PDB_model_num = %r, which is not an integer"
                    % (path, row_number, row[model_col]))

        atoms.append(AtomRecord(
            record=(row[group_col] if group_col is not None and row[group_col] else 'ATOM'),
            name=name,
            altloc=row[altloc_col] if altloc_col is not None else '',
            resname=_first_present(row, auth_comp, label_comp),
            chain=_first_present(row, auth_asym, label_asym),
            resseq=_first_present(row, auth_seq, label_seq),
            icode=row[icode_col] if icode_col is not None else '',
            x=_cif_float(row, x_col, '_atom_site.Cartn_x', row_number, path),
            y=_cif_float(row, y_col, '_atom_site.Cartn_y', row_number, path),
            z=_cif_float(row, z_col, '_atom_site.Cartn_z', row_number, path),
            occupancy=occupancy,
            element=element,
            model=model))
    return atoms


def _cif_parents(tags: list[str], rows: list[list[str]]) -> dict[str, str]:
    """
    Read modified residue -> parent residue from the _pdbx_struct_mod_residue loop.

    Parameters
    ----------
    tags : list of str
        The loop's lower-cased data names.
    rows : list of list of str
        The loop's rows.

    Returns
    -------
    dict of str to str
        Modified residue name -> parent residue name (empty if the loop lacks
        either column).
    """
    modified_col = _column(tags, '_pdbx_struct_mod_residue.auth_comp_id',
                           '_pdbx_struct_mod_residue.label_comp_id')
    parent_col = _column(tags, '_pdbx_struct_mod_residue.parent_comp_id')
    if modified_col is None or parent_col is None:
        return {}
    parents = {}
    for row in rows:
        if row[modified_col] and row[parent_col]:
            parents[row[modified_col]] = row[parent_col]
    return parents


def _read_cif(handle: TextIO, path: Path) -> StructureFile:
    """
    Read the atoms of the first data block of an mmCIF file.

    Categories are read whether they are written as a loop_ or, for a
    category with a single row, as data-name/value pairs.

    Parameters
    ----------
    handle : TextIO
        The open file.
    path : Path
        The file being read (recorded, and used in error messages).

    Returns
    -------
    StructureFile
        The atoms, model numbers, modified-residue parents and entry identifier.

    Raises
    ------
    StructureFileException
        If the file is malformed or has no _atom_site category.
    """
    stream = _CifTokenStream(_cif_tokens(handle, path))
    atom_tags: Optional[list[str]] = None
    atom_rows: list[list[str]] = []
    parents: dict[str, str] = {}
    entry_id = ''
    # categories written as name/value pairs: data name -> value
    pairs: dict[str, str] = {}
    seen_block = False

    while True:
        token = stream.take()
        if token is None:
            break
        text, quoted = token
        if quoted:
            continue
        lower = text.lower()
        if lower.startswith('data_'):
            if seen_block:
                break
            seen_block = True
        elif lower == 'loop_':
            tags, rows = _read_loop(stream, path)
            if tags[0].startswith(CIF_ATOM_SITE):
                atom_tags, atom_rows = tags, rows
            elif tags[0].startswith(CIF_MOD_RESIDUE):
                parents.update(_cif_parents(tags, rows))
        elif lower.startswith('_'):
            value = stream.take()
            if value is None:
                raise StructureFileException("%s: the data name %s has no value" % (path, text))
            pairs[lower] = _cif_value(value)

    entry_id = pairs.get(CIF_ENTRY_ID, '')
    if atom_tags is None:
        # a one-atom file may give _atom_site as name/value pairs
        single = sorted(name for name in pairs if name.startswith(CIF_ATOM_SITE))
        if not single:
            raise StructureFileException("%s holds no _atom_site category, so it has no atoms" % path)
        atom_tags, atom_rows = single, [[pairs[name] for name in single]]
    single_parent = sorted(name for name in pairs if name.startswith(CIF_MOD_RESIDUE))
    if single_parent:
        parents.update(_cif_parents(single_parent, [[pairs[name] for name in single_parent]]))

    atoms = _cif_atoms(atom_tags, atom_rows, path)
    return StructureFile(path=path, file_format='mmCIF', atoms=atoms, models=_models_in_order(atoms),
                         parent_residues=parents, entry_id=entry_id)
