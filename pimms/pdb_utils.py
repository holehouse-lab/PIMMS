## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


import numpy as np

from .latticeExceptions import PDBException
from . import CONFIG 
from . import IO_utils
from string import ascii_lowercase, ascii_uppercase, digits

# PDB chain identifiers, one per chain type: A-Z first (so runs with up to 26
# types are labelled exactly as before), then a-z and 0-9
ALPHABET = ascii_uppercase + ascii_lowercase + digits


ATOM_NAME='CA'

def one_to_three(res):
    """
    Function which converts the residue code into a PDB compliant
    3 letter code. Gets the amino acids right and then will build
    new residue names on the fly for those it doesn't know...

    Parameters
    ------------
    res : str
        The residue code to be converted. Normally a one-letter amino acid code, but any
        string is accepted.

    Returns
    ---------
    str
        The 3 letter PDB compliant residue code. Codes not in CONFIG.ONE_TO_THREE are
        truncated to their first three characters, or left-padded with 'X' if shorter
        than three characters.


    """
    if res in CONFIG.ONE_TO_THREE:
        return CONFIG.ONE_TO_THREE[res]
    else:
        if len(res) > 3:
            return res[0:3]
        else:

            # note we need the +res at end so that first letter
            # of the resid is a letter (X)
            return (3-len(res))*'X'+res



def write_positions_to_file(positions, filename, spacing, dimensions=False, sequence=False):
    """
    Function which takes a list of positions on a 
    lattice and writes them to a single PDB file. Note this does
    not facilitate including bead/residue types - so may be more
    useful for graphical debugging. The dimensions of the lattice 
    can simply be inferred from the positions provided, and a 4
    site cushion is provided around the min and max values in
    each dimension. Alternativly dimenions can just be provided
    directly.

    Parameters
    ----------------
    
    positions : list of list of ints
        A set of positions on the lattice - note we assume these are non
        overlapping...

    filename : string
        Name of the file to be written (include the .pdb extension please)

    spacing : float 
        Spacing between lattice sites - i.e. 4 means each site is
        4 angstroms appart.

    dimensions : list of int or bool, optional
        Lattice dimensions as a 2- or 3-element list, which must match the
        dimensionality of the positions. If False (the default) the dimensions
        are inferred from the largest position in each dimension.

    sequence : str or bool, optional
        If provided, defines the amino acid sequence using standard amino acid
        one letter code, and must be the same length as positions. If False (the
        default) every site is written as residue 'X'.

    Returns
    -------
    None
        No return value; a finalized PDB file is written to ``filename``.

    Raises
    ------
    PDBException
        If ``positions`` is empty, not 2D, of unsupported dimensionality,
        falls outside supplied ``dimensions``, or if ``sequence`` is
        provided but is not a string of matching length.

    """

    if len(positions) == 0:
        raise PDBException("Trying to write positions to PDB but positions is empty")

    positions_np = np.array(positions)
    if positions_np.ndim != 2:
        raise PDBException("Trying to write positions to PDB but positions must be a 2D array-like")

    n_dim = positions_np.shape[1]
    if n_dim not in (2, 3):
        raise PDBException(f"Trying to write positions to PDB with unsupported dimensionality ({n_dim})")

    max_pos = positions_np.max(axis=0)


    if dimensions:
        # check the number of dimensions match
        if len(dimensions) != n_dim:
            raise PDBException('Trying to write a PDB file where the positions appear to have %i dimensions, but the supplied lattice dimensions (%s) does not match'%(n_dim, str(dimensions)))

        # check the largest position is inside the dimensions max
        for dim in range(0, n_dim):
            if max_pos[dim] >= dimensions[dim]:
                raise PDBException('Trying to write a PDB file where the dimensions provided are [%s], but there is a position that lies outside this (%s) in the %i dimension' % (str(dimensions), max_pos[dim], dim))

        local_dimensions = dimensions

    else:
        local_dimensions = []
        for dim in range(0, n_dim):
            # Dimension lengths are one-indexed relative to max coordinate.
            local_dimensions.append(int(max_pos[dim]) + 1)
    
    # sanity check passed sequence 
    if sequence is not False:
        if type(sequence) is not str:
            raise PDBException('Trying to pass a sequence which is not a string')
            
        if len(sequence) != len(positions):
            raise PDBException('Trying to pass a sequence which is not the same length as the number of positions')

    UPO={}
    UPO['dimensions'] = local_dimensions
    UPO['length'] = len(positions)
    UPO['positions'] = positions
    UPO['sequence'] = sequence
    
    # write CRYST line and initialize the PDB file
    initialize_pdb_file(local_dimensions, spacing, filename)

    # add positions
    build_pdb_file([], spacing, filename=filename, usePositionsOnly=UPO)

    # end the PDB file
    finalize_pdb_file(filename)



def build_pdb_file(latticeObject, spacing, filename='lattice.pdb', sequence=False, usePositionsOnly=None, write_connect=False, autocenter=False, unwrap=False):
    """
    Function which writes a PDB file based on lattice or postition information. The normal usage
    is to pass a latticeObject and write the whole lattice to file. However, one can also just pass
    a list of arbitrary positions using the usePositionsOnly.

    The filename MUST have already been initialized using the initialize_PDB_file() function. This
    function creates a new (empty) file and writes the CRYST line (top line at the start of a PDB
    file that defines the crystalographic symmetry group and dimensions).

    Parameters
    -----------

    latticeObject : Lattice
        The lattice object being written out. Ignored (and typically passed as an empty list) when
        usePositionsOnly is used.

    spacing : float
        How lattice spacing is converted to real-world spacing.

    filename : str, optional
        Name of PDB file being written, which must already have been initialized with
        initialize_pdb_file(). Default is lattice.pdb

    sequence : list of str or bool, optional
        Vestigial argument. It is not read anywhere in the function body, so passing it has no
        effect; residue names come from the chains in latticeObject, or from the 'sequence' entry of
        the usePositionsOnly dictionary. Default is False.

    usePositionsOnly : dict, optional
        If provided, the function will use this dictionary to construct the output rather than
        reading the lattice object. A usePositionsOnly dictionary has a specific structure and has
        exactly four key/value pairs. These are:

            dimensions : box dimensions
            length     : number of residues in the chain
            positions  : a list of lists, where each sublist has the positions of a bead.
            sequence   : one-letter-code sequence string, or False to write every bead as 'X'

        Default is None.

    write_connect : bool, optional
        Flag which, if set to True, will write one CONECT record per backbone bond. Ignored
        when usePositionsOnly is used. Once the file holds more than 99999 serials (beads
        plus one per TER record) the 5-column serials wrap modulo 100000, and every bond
        with an endpoint whose written serial is shared by two records (or that lies past
        the wrap) is omitted, so no CONECT can name the wrong atom; bonds between beads
        with unique serials are still written. Default is False.

    autocenter : bool, optional
        Flag which, if set to True and there's a single chain will center the protein in the box.
        This is useful for visualization purposes but does mean any translational diffusion will
        be lost. Default = False

    unwrap : bool, optional
        Flag which, if set to True, writes each chain as a single "whole" periodic
        image (not torn across a box face); coordinates may fall outside the box.
        Ignored where ``autocenter`` applies (autocenter already unwraps). Default False.

    Returns
    --------
    None
        Does not return anything but appends a MODEL/ATOM/TER/ENDMDL block to the PDB file on disk.
        Note the file is NOT closed off with an END line here - that is done by finalize_pdb_file().

    Raises
    ------
    PDBException
        If usePositionsOnly is not a four-entry dictionary, is missing a required key, or holds a
        sequence whose length does not match its positions, or if the box is neither 2D nor 3D.

    """
    
    ## Internal functions to avoid code duplication
    ##
    ## <><><><><><><><><><><><><><><><><><><><><><><><><>
    def segupdate(resindex_num, segment):
        """
        Roll the residue counter over to a new PDB segment past 9999.

        PDB residue IDs are limited to four columns, so when
        ``resindex_num`` exceeds 9999 it is reset to 1 and the segment
        index is incremented.

        Parameters
        ----------
        resindex_num : int
            Current residue number used for the PDB RESID column.

        segment : int
            Current PDB segment index.

        Returns
        -------
        tuple
            Updated ``(resindex_num, segment)`` pair.

        """
        if resindex_num > 9999:
            resindex_num = 1
            segment = segment+1
        return (resindex_num, segment)

    ##  .................................................
    def update_increments(i, resindex, resindex_num):
        """
        Increment the per-bead (ATOM serial) and per-residue counters by one.

        Parameters
        ----------
        i : int
            Running ATOM serial, one per bead.

        resindex : int
            Running index into the chain sequence.

        resindex_num : int
            Running residue number used for the PDB RESID column.

        Returns
        -------
        tuple
            Updated ``(i, resindex, resindex_num)`` triple.

        """
        i = i + 1
        resindex = resindex + 1
        resindex_num = resindex_num + 1
        return (i, resindex, resindex_num)



    ## <><><><><><><><><><><><><><><><><><><><><><><><><>
    if usePositionsOnly is not None:
        
        # first we validate this
        if not isinstance(usePositionsOnly, dict) or len(usePositionsOnly) != 4:
            print(usePositionsOnly)
            raise PDBException("In 'build_pdb_file' trying to generate a file using the usePositionsOnly only setting but an INVALID usePositionsOnly dictionary was passed")

        try:
            dimensions = usePositionsOnly['dimensions']

            if usePositionsOnly['sequence'] is False:
                chain_seq  = list('X'*usePositionsOnly['length'])
            else:
                chain_seq = usePositionsOnly['sequence']
                
            positions  = usePositionsOnly['positions']
            chains_list = [1]

            if len(chain_seq) != len(positions):
                raise PDBException("When extracting usePositionsOnly data found mismatch between sequence and number of positions")                

        except KeyError:
            print(usePositionsOnly)
            raise PDBException("In 'build_pdb_file' trying to generate a file using the usePositionsOnly only setting but an INVALID usePositionsOnly dictionary was passed (missing one of the keywords)")
                    
    else:
        chains_list = latticeObject.chains
        dimensions  = latticeObject.dimensions

        # one identifier per chain type; past the 62 available labels every
        # further type shares the last one, which is announced rather than done
        # silently (it used to alias onto 'Z', the label of the 26th type)
        all_pdb_chain_ids = {}
        alphabet_idx = 0
        collapsed = 0
        for chainID in chains_list:
            if chains_list[chainID].chainType not in all_pdb_chain_ids:
                if alphabet_idx < len(ALPHABET):
                    all_pdb_chain_ids[chains_list[chainID].chainType] = ALPHABET[alphabet_idx]
                    alphabet_idx = alphabet_idx + 1
                else:
                    all_pdb_chain_ids[chains_list[chainID].chainType] = ALPHABET[-1]
                    collapsed += 1
        if collapsed:
            IO_utils.status_message(
                "PDB output has only %i distinct chain identifiers; %i further chain type(s) "
                "share the identifier '%s'. Analyses that read chain types from the PDB "
                "(e.g. lemonade without a keyfile) will merge them." % (len(ALPHABET), collapsed, ALPHABET[-1]),
                'warning')
  

    CONNECT_RECORDS = []

    with open(filename,'a') as fh:
        fh.write(build_model_line(1)+"\n")
        
        i=1
        segment=1
        for chainID in chains_list:

            
            if usePositionsOnly:
                # if use positions these are set at the start
                pdb_chain_ID = 'A'
            else:
                # else define for each chain

                # autocenter is only valid for a single chain; then pick the output
                # convention (autocenter / PBC-unwrap / raw) for this chain
                use_autocenter = autocenter and len(latticeObject.chains) == 1
                positions = latticeObject.chains[chainID].get_output_positions(autocenter=use_autocenter, unwrap=unwrap)
                    
                chain_seq = latticeObject.chains[chainID].sequence
                pdb_chain_ID = all_pdb_chain_ids[latticeObject.chains[chainID].chainType]

            resindex  = 1    # used to index into the sequence
            resindex_num = 1 # used to record the RESID in the PDB file
            previous_i = None
            
            if len(dimensions) == 2:                        
             
                first_in_chain = True
                for position in positions:   

                    resindex_num, segment = segupdate(resindex_num, segment)
                    fh.write(build_atom_line(i, ATOM_NAME, one_to_three(chain_seq[resindex-1]), pdb_chain_ID,   str(resindex_num),       float(position[0])*spacing, float(position[1])*spacing, 0.0,segment))                    

                    # connect record info
                    if first_in_chain:
                        first_in_chain = False
                    else:                            
                        CONNECT_RECORDS.append([previous_i, i])
                    previous_i = i
                    
                    i, resindex, resindex_num = update_increments(i, resindex, resindex_num)


            elif len(dimensions) == 3:

                first_in_chain = True
                for position in positions:                      
                    resindex_num, segment = segupdate(resindex_num, segment)
                    fh.write(build_atom_line(i, ATOM_NAME, one_to_three(chain_seq[resindex-1]), pdb_chain_ID,   str(resindex_num),       float(position[0])*spacing, float(position[1])*spacing, float(position[2])*spacing, segment))

                    # connect record info
                    if first_in_chain:
                        first_in_chain = False
                    else:                            
                        CONNECT_RECORDS.append([previous_i, i])                        
                    previous_i = i

                    i, resindex, resindex_num = update_increments(i, resindex, resindex_num)

            else:
                raise PDBException('Unusable number of dimensions...')
        
            # resindex_num was already incremented past the last written residue by
            # update_increments, so the TER record must use resindex_num - 1: the PDB
            # spec says TER carries the SAME residue number as the terminal residue,
            # and for a 9999-residue chain the off-by-one (10000) also overflowed the
            # 4-column field and crashed the write.
            fh.write(build_ter_line(i, one_to_three(chain_seq[resindex-2]), pdb_chain_ID, resindex_num - 1))
            i=i+1


        # if we want to write the connect record...
        if write_connect:
            if usePositionsOnly is None:

                # Serials wrap at 100000 in the ATOM/TER records (5-column field), so
                # once the file needs more than 99999 of them the low serials are
                # written twice: serial s and serial s + 100000 both appear as s. A
                # CONECT can only name the written serial, and readers resolve a
                # duplicated serial to one of its carriers (mdtraj takes the LAST),
                # so a bond naming a shared serial gets drawn between the wrong atoms
                # - often across chains. Skipping only the bonds past the wrap point
                # was not enough for exactly this reason. A serial is unambiguous iff
                # it is below the wrap AND its wrapped twin (s + 100000) was never
                # used; i is one past the last serial used (the final TER), so the
                # safe window is max_serial - 100000 < s < 100000. Every bond with
                # both ends in that window is written; every other bond is omitted,
                # which drops bonds rather than ever drawing a wrong one.
                max_serial = i - 1
                for record in CONNECT_RECORDS:
                    if any(s >= 100000 or s + 100000 <= max_serial for s in record):
                        continue
                    fh.write(build_conect_line(record[0], record[1]))
                         
        fh.write("ENDMDL\n")


            




#-----------------------------------------------------------------
#
def finalize_pdb_file(filename='lattice.pdb'):

    """
    Function that finalizes a PDB file with an END line

    Parameters
    --------------
    filename : str, optional
        Filename to append the END line to. Default is lattice.pdb

    Returns
    -------------
    None
        No return type, but the PDB file is written to with an
        END line.

    """

    with open(filename,'a') as fh:
        fh.write("END\n")



#-----------------------------------------------------------------
#
def initialize_pdb_file(dimensions, spacing, filename='lattice.pdb'):
    """
    Initialize PDB with box dimensions line (CYRST line)


    Parameters
    --------------
    dimensions : list
        A list of length 2 or 3, depending on the dimensionality of the system
        being studied, that reflects the lattice dimensions.

    spacing : float
        Lattice-to-realspace spacing in angstroms.

    filename : str, optional
        Filename to write to. The file is created (or truncated if it already exists),
        so this must be called before anything is appended to it. Default is lattice.pdb

    Returns
    -------------
    None
        No return type, but a new PDB file containing just the CRYST1 line is written out

    Raises
    ------
    PDBException
        If dimensions is neither 2 nor 3 elements long.

    """
    with open(filename,'w') as fh:
        fh.write(build_cryst_line(dimensions, spacing))
        
        
    
#-----------------------------------------------------------------
#
def build_section_string(content, length, justification='L'):
    """
    Generates a string for a PDB element, where you can define
    the string length, string content, and justifictaion as
    L, C, or R.

    Parameters
    --------------
    content : str
        The content of the string

    length : int
        The number of columns the returned string must fill

    justification : {'L', 'C', 'R'}, optional
        The justification of the content within those columns, either L (left), C (centred)
        or R (right). For centred content with an odd amount of padding the extra space goes
        on the right. Default is 'L'.

    Returns
    -------------
    return_string : str
        The padded string, exactly length characters long

    Raises
    ------
    PDBException
        If content is longer than length, or if justification is not one of L, C or R.

    """

    content_len = len(content)
    filler      = length - content_len

    if content_len > length:
        raise PDBException("Trying to build a PDB section string but the content [%s] is longer than the allowed column (len=%i)"%(content, length))
    

    if justification == 'L':
        return_string = content + (filler)*" "        

    elif justification == 'R':
        return_string = (filler)*" " + content
    
    elif justification == "C":
        if filler %2 == 0:
            RHS = int(filler/2)
            LHS = RHS
        else:            
            LHS = int(filler/2)
            RHS = LHS+1
        return_string = LHS*" " + content + RHS*" " 
    else:
        raise PDBException('Invalid section justification provided [%s]'%justification)

    return return_string



#-----------------------------------------------------------------
#
def build_line(section_list, section_columns):
    """
    Place content into fixed columns to build an 80-character PDB line.

    Section columns are defined by the PDB specification and are
    corrected here for the -1 (0-based vs. 1-based) offset. Each chunk
    of content is written into its corresponding column range and the
    result is padded out to a full 80-character line.

    Parameters
    ----------
    section_list : list
        Ordered list of content strings to place into the line.

    section_columns : list of [int, int]
        Ordered list of inclusive, 1-based ``[start, end]`` column ranges
        (one per entry in ``section_list``) defining where each piece of
        content is written.

    Returns
    -------
    str
        An 80-character line with the content placed into the specified
        columns and padded with spaces.

    Raises
    ------
    PDBException
        If ``section_list`` and ``section_columns`` differ in length, or if the
        assembled content exceeds 80 characters.

    """
    if len(section_list) != len(section_columns):
        raise PDBException(
            "PDB line content and column specifications must have equal length")

    line = [""]*80
    
    for (content, region) in zip(section_list, section_columns):
        if region[0] == region[1]:
            line[region[0]-1] = content
        else:
            line[region[0]-1:region[1]-1] = content
        
    init_string = "".join(line)

    extra = 80 - len(init_string)
    if extra < 0:
        raise PDBException('Line has ended up being longer than 80? - line shown below for debugging...\n%s'%(init_string))

    init_string = init_string + extra*" "
    
    return init_string        



#-----------------------------------------------------------------
#
def build_model_line(serial):
    """
    As defined by wwpdb.org
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#MODEL

    Constructs a valid MODEL line for a PDB file. Note 'serial' here is basically the 
    frame number, which for a PDB file will always be 1 but we allow it to be passed
    in for consistency with other file formats.

    Parameters
    --------------
    serial : int
        The serial number of the model

    Returns
    -------------
    line : str
        The fully-formatted MODEL line, as defined by the PDB specification.
    
    """
    
    name_section   = build_section_string('MODEL', 6, 'L')       # 1  - 6
    BREAK_1        = "    "                                      # 7  - 10
    serial_section = build_section_string(str(serial), 4,  'R')   # 11 - 14
    
    return build_line([name_section, BREAK_1, serial_section],[[1,6],[7,10],[11,14]])



#-----------------------------------------------------------------
#
def build_atom_line(atom_index, atom_name, res_name, chain, res_id, x,y,z, segment):
    """
    As defined by wwpdb.org
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    Constructs a valid ATOM line for a PDB file.

    Parameters
    --------------

    atom_index : int
        The index of the atom in the system

    atom_name : str
        The name of the atom

    res_name : str
        The name of the residue

    chain : str
        The single-character chain identifier

    res_id : int or str
        The residue ID, written into the four-column RESID field

    x : float
        The x coordinate of the atom, in angstroms

    y : float
        The y coordinate of the atom, in angstroms

    z : float
        The z coordinate of the atom, in angstroms

    segment : int or str
        The segment identifier, written into the four-column segment field

    Returns
    -------------
    line : str
        The fully-formatted ATOM line, as defined by the PDB specification, with a
        trailing newline.

    """
    
    ATOM      = build_section_string("ATOM",          6, 'L') # 1  - 6
    # serials only have 5 columns; wrap at 100000 rather than crashing at
    # topology-write time for large systems (>=100k beads). mdtraj and most other
    # readers rebuild indices sequentially, so wrapped serials load fine; only CONECT
    # records address atoms by serial, and build_pdb_file omits the ambiguous ones.
    ATOM_IDX  = build_section_string(str(atom_index % 100000), 5, 'R') # 7  - 11
    BREAK_1   = " "                                           # 12
    ATOM_NAME = build_section_string(str(atom_name),  4, 'C') # 13 - 16
    ALTLOC    = " "                                           # 17
    RES_NAME  = build_section_string(str(res_name),   3, 'L') # 18 - 20 | note we don't actually expect the residue names to be anything other than a 3 letter code
    BREAK_2   = " "                                           # 21
    CHAIN     = build_section_string(str(chain),      1, 'L') # 22 
    RES_ID    = build_section_string(str(res_id),     4, 'R') # 23 - 26
    ICODE     = " "                                           # 27
    BREAK_3   = "   "                                         # 29 - 33
    X         = build_section_string("%8.3f" % x,     8, 'R')      # 31 - 38
    Y         = build_section_string("%8.3f" % y,     8, 'R')      # 39 - 46
    Z         = build_section_string("%8.3f" % z,     8, 'R')      # 47 - 54
    BREAK_4   = "                  "                          # 55 - 72
    SEG       = build_section_string(str(segment), 4, 'L')    # 73 - 76

    return build_line([ATOM, ATOM_IDX, BREAK_1,ATOM_NAME, ALTLOC, RES_NAME, BREAK_2, CHAIN, RES_ID, ICODE, BREAK_3, X, Y, Z, BREAK_4, SEG],  [[1,6],[7,11],[12,12],[13,16],[17,17],[18,20],[21,21],[22,22],[23,26],[27,27],[28,30],[31,38],[39,46],[47,54],[55,72],[73,76]])+"\n"


#-----------------------------------------------------------------
#
def build_conect_line(atom1, atom2):
    """
    As defined by wwpdb.org
    https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html

    Constructs a valid CONECT line.

    Parameters
    --------------

    atom1 : int
        The serial number of the first atom in the bond

    atom2 : int
        The serial number of the second atom in the bond

    Returns
    -------------
    line : str
        The fully-formatted CONECT line, as defined by the PDB specification, with a
        trailing newline.

    """

    CONECT_DEF = 'CONECT' # 1 - 6
    ATOM1_LINE  = build_section_string(str(atom1), 5, 'R') # 7  - 11
    ATOM2_LINE  = build_section_string(str(atom2), 5, 'R') # 12  - 16

    return build_line([CONECT_DEF, ATOM1_LINE, ATOM2_LINE],[[1,6], [7,11],[12,16]]) +'\n'
    

#-----------------------------------------------------------------
#
def build_ter_line(atom_index, res_name, chain, res_id):
    """
    As defined by wwpdb.org
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#TER

    Constructs a valid TER line for a PDB file    

    Parameters
    --------------

    atom_index : int
        The serial number of the terminating record (wrapped at 100000 to fit the
        five-column serial field)

    res_name : str
        The three-letter name of the terminal residue

    chain : str
        The single-character chain identifier

    res_id : int or str
        The residue ID of the terminal residue, which per the PDB spec must be the
        same number as that residue's ATOM records

    Returns
    -------------
    line : str
        The fully-formatted TER line, as defined by the PDB specification, with a
        trailing newline.

    """

    TER_SEC = "TER   "                                        # 1  - 6
    ATOM_IDX  = build_section_string(str(atom_index % 100000), 5, 'R') # 7  - 11
    BREAK_1   = "      "                                      # 12 - 17
    RES_NAME  = build_section_string(str(res_name),   3, 'L') # 18 - 20
    BREAK_2   = " "                                           # 21
    CHAIN     = build_section_string(str(chain),      1, 'L') # 22 
    RES_ID    = build_section_string(str(res_id),     4, 'R') # 23 - 26
    ICODE     = " "                                           # 27
    
    return build_line([TER_SEC, ATOM_IDX, BREAK_1, RES_NAME, BREAK_2, CHAIN, RES_ID, ICODE],  [[1,6],[7,11],[12,17],[18,20],[21,21],[22,22],[23,26],[27,27]])+"\n"


#-----------------------------------------------------------------
#
def build_cryst_line(dimensions, spacing):
    """
    As defined by wwpdb.org
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    Constructs a valid CRYST1 line for a PDB file - this creates an appropriate box size in 2 or 3
    dimensions, where spacing defines how lattice sites relate to angstroms (i.e. by default each
    lattice site is 3.65 angstroms apart --> inter-amino acid distance).

    The CRYST1 record defines the *periodic unit cell*, so an ``L``-site axis has a period of
    ``L * spacing`` angstroms, NOT ``(L - 1) * spacing``: sites ``L-1`` and ``0`` are periodic
    neighbours exactly one lattice unit apart. This previously wrote ``(L - 1) * spacing``, which
    (a) disagreed with the box vectors PIMMS writes into the XTC (see
    ``lattice_utils._lattice_frame_xyz_and_box``, which has always used ``L * spacing``), (b) made
    every PBC-aware calculation performed on START.pdb by mdtraj/VMD wrong by one lattice unit, and
    (c) made ``pimms.lemonade.load(pdb=...)`` (with no keyfile, so the box has to be inferred from
    the file) infer an ``L-1`` box and then wrap the coordinates into it, silently corrupting them.

    For a 2D system the z axis has no period, so ``c`` is set to a single lattice unit - again
    matching what the XTC writer does.

    Parameters
    ----------------
    dimensions : list
        A list of 2 or 3 elements in length that defines the x, y, and maybe z dimensions

    spacing : float
        Multiplier that converts lattice spacing to real-world spacing in the PDB. In units
        of Angstroms. i.e. 4 would mean 4 per lattice unit.

    Returns
    -------------
    str
        Returns a fully-formatted valid CRYST1 line for a PDB file, with a trailing
        newline

    Raises
    ------
    PDBException
        If dimensions is neither 2 nor 3 elements long.

    """

    if len(dimensions) not in (2, 3):
        raise PDBException(f"CRYST line only supports 2D/3D dimensions, got {len(dimensions)}")

    CRYST_SECT = "CRYST1"                                               # 1  - 6
    a          = build_section_string("%9.3f" % (dimensions[0]*spacing), 9, 'R')  # 7  - 15
    b          = build_section_string("%9.3f" % (dimensions[1]*spacing), 9, 'R')  # 16 - 24


    # set the third dimension depending on lattice type
    if len(dimensions) == 3:
        c      = build_section_string("%9.3f" % (dimensions[2]*spacing), 9, 'R')  # 25 - 33
    else:
        # 2D: no periodicity in z, so use a single lattice unit (matches the XTC box)
        c      = build_section_string("%9.3f" % spacing, 9, 'R')        # 25 - 33

    alpha      = build_section_string("%7.2f" % 90.0, 7, 'R')           # 34 - 40
    beta       = build_section_string("%7.2f" % 90.0, 7, 'R')           # 41 - 47
    gamma      = build_section_string("%7.2f" % 90.0, 7, 'R')           # 48 - 54
    BREAK_1    = " "                                                    # 55 - 55
    sGroup     = build_section_string("P 1", 10,"L")                    # 56 - 66
    ZVALUE     = "    "                                                 # 67 - 70
    
    return build_line([CRYST_SECT, a, b, c, alpha, beta, gamma, BREAK_1, sGroup, ZVALUE],  [[1,6],[7,15],[16,24],[25,33],[34,40],[41,47],[48,54],[55,55],[56,66],[67,70]])+"\n"

    
    

    
    
                                          
        
            

        
        
        
            
    
    
