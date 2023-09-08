## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2023
## ...........................................................................


import numpy as np

from .latticeExceptions import PDBException
from . import CONFIG 
from string import ascii_uppercase as ALPHABET


ATOM_NAME='CA'

def one_to_three(res):
    """
    Function which converts the residue code into a PDB compliant
    3 letter code. Gets the amino acids right and then will build
    new residue names on the fly for those it doesn't know...

    Parameters
    ------------
    res : string
        The residue code to be converted

    Returns
    ---------
    string
        The 3 letter PDB compliant residue code


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



def write_positions_to_file(positions, filename, spacing, dimensions=False):
    """
    Function which takes a list of positions on a 
    lattice and writes them to a single PDB file. Note this does
    not facilitate including atom/residue types - so may be more
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

    dimensions : 2D or 3D list of ints (DEFAULT=False)
        Lattice dimensions, if set to false dimensions are inferred
        from the set of positions


    """
    
    n_dim = len(positions[0])

    max_pos = []
    min_pos = []

    for dim in range(0, n_dim):
        max_pos.append(max(np.transpose(positions)[dim]))
        min_pos.append(min(np.transpose(positions)[dim]))


    if dimensions:
        # check the number of dimensions match
        if len(dimensions) != n_dim:
            raise PDBException('Trying to write a PDB file where the positions appear to have %i dimensions, but the supplied lattice dimensions (%s) does not match'%(n_dim, str(dimensions)))

        # check the largest position is inside the dimensions max
        for dim in range(0, n_dim):            
            if max_pos[dim] > dimensions[dim]:
                raise PDBException('Trying to write a PDB file where the dimensions provided are [%s], but there is a position that lies outside this (%s) in the %i dimension' % (str(dimensions), max_pos[dim], dim))

        local_dimensions = dimensions

    else:
        local_dimensions = []
        for dim in range(0, n_dim):
            local_dimensions.append(max_pos[dim])

    UPO={}
    UPO['dimensions'] = local_dimensions
    UPO['length'] = len(positions)
    UPO['positions'] = positions
    
    # write CRYST line and initialize the PDB file
    initialize_pdb_file(local_dimensions, filename)

    # add positions
    build_pdb_file([], spacing, filename=filename, usePositionsOnly=UPO)

    # end the PDB file
    finalize_pdb_file(filename)



def build_pdb_file(latticeObject, spacing, filename='lattice.pdb', usePositionsOnly=None, write_connect=False, autocenter=False):
    """
    Function which writes a PDB file based on lattice or postition information. The normal usage
    is to pass a latticeObject and write the whole lattice to file. However, one can also just pass
    a list of arbitrary positions using the usePositionsOnly.

    The filename MUST have already been initialized using the initialize_PDB_file() function. This
    function creates a new (empty) file and writes the CRYST line (top line at the start of a PDB
    file that defines the crystalographic symmetry group and dimensions).

    Parameters
    -----------

    latticeObject : LatticeObject
        The lattice object being written out

    spacing : float 
        How lattice spacing is converted to real-world spacing. 

    filename : string
        Name of PDB file being written. Default is lattice.pdb
    
    usePositionsOnly : dict
        If provided, the function will use this dictionary to construct the output. A usePostionsOnly
        dictionary has a specific structure and has three key/value pairs. These are:

            dimensions : box dimensions
            length     : number of residues in the chain
            positions  : a list of lists, where each sublist has the positions of a bead.
      
        Default = None

    write_connect : bool
        Flag which, if set to True, will write CONNECT records if possible

    autocenter : bool
        Flag which, if set to True and there's a single chain will center the protein in the box. 
        This is useful for visualization purposes but does mean any translational diffusion will
        be lost. Default = False
    
    Returns
    --------
    None
        Does not return anything but writes a PDB file to disk based on the input information. Note
        the PDB file written is a fully finalized PDB file that 

    """
    
    ## Internal functions to avoid code duplication
    ##
    ## <><><><><><><><><><><><><><><><><><><><><><><><><>
    def segupdate(resindex_num, segment):
        if resindex_num > 9999:                        
            resindex_num = 1
            segment = segment+1
        return (resindex_num, segment)

    ##  .................................................
    def update_increments(i, resindex, resindex_num):
        i = i + 1
        resindex = resindex + 1
        resindex_num = resindex_num + 1
        return (i, resindex, resindex_num)



    ## <><><><><><><><><><><><><><><><><><><><><><><><><>
    if usePositionsOnly is not None:
        
        # first we validate this
        if len(usePositionsOnly) != 3 or type(usePositionsOnly) != dict:
            print(usePositionsOnly)
            raise PDBException("In 'build_pdb_file' trying to generate a file using the usePositionsOnly only setting but an INVALID usePositionsOnly dictionary was passed")

        try:
            dimensions = usePositionsOnly['dimensions']
            n_chains   = 1
            chain_seq  = list('X'*usePositionsOnly['length'])
            positions = usePositionsOnly['positions']
            chains_list = [1]

            if len(chain_seq) != len(positions):
                raise PDBException("When extracting usePositionsOnly data found mismatch between sequence and number of positions")                

        except KeyError:
            print(usePositionsOnly)
            raise PDBException("In 'build_pdb_file' trying to generate a file using the usePositionsOnly only setting but an INVALID usePositionsOnly dictionary was passed (missing one of the keywords)")
                    
    else:
        chains_list = latticeObject.chains
        n_chains    = len(chains_list)
        dimensions  = latticeObject.dimensions

        # if we have more than 26 different types of chains 
        all_pdb_chain_ids = {}
        alphabet_idx = 0
        for chainID in chains_list:
            

            if chains_list[chainID].chainType not in all_pdb_chain_ids:
                try:
                    all_pdb_chain_ids[chains_list[chainID].chainType] = ALPHABET[alphabet_idx]
                    alphabet_idx = alphabet_idx + 1
                except IndexError:
                    all_pdb_chain_ids[chains_list[chainID].chainType] = 'Z'
  

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

                # note if we want to and can autocenter...                
                if autocenter and len(latticeObject.chains) == 1:
                    positions = latticeObject.chains[chainID].get_ordered_positions(center_positions=True)

                else:
                    positions = latticeObject.chains[chainID].get_ordered_positions()    
                    
                
                #positions = latticeObject.chains[chainID].get_ordered_positions()
                chain_seq = latticeObject.chains[chainID].sequence
                pdb_chain_ID = all_pdb_chain_ids[latticeObject.chains[chainID].chainType]

            resindex  = 1    # used to index into the sequence
            resindex_num = 1 # used to record the RESID in the PDB file
            
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
        
            fh.write(build_ter_line(i, one_to_three(chain_seq[resindex-2]), pdb_chain_ID, resindex_num))
            i=i+1


        # if we want to write the connect record...
        if write_connect:
            if usePositionsOnly is None:
                for record in CONNECT_RECORDS:
                    fh.write(build_conect_line(record[0], record[1]))
                         
        fh.write("ENDMDL\n")


            




#-----------------------------------------------------------------
#
def finalize_pdb_file(filename='lattice.pdb'):

    """
    Function that finalizes a PDB file with an END line

    Parameters
    --------------
    filename : str
        Filename to write to

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

    filename : str
        Filename to write to

    Returns
    ------------- 
    None
        No return type, but a filename is written out

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
        The length of the string

    justification : str
        The justification of the string, either L, C, or R

    Returns
    -------------
    return_string : str
        The padded string 

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
    Section columns as defined by PDB specification - correct for -1 offset!

    Takes a set of data that defines content (section_lists) and position (section_columns)
    and constructs a string where content is placed into the right columns

    Returns a 80-character string

    """
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
    BREAK_1        = "   "                                       # 7  - 10
    serial_section = build_section_string(str(serial), 4,  'L')   # 11 - 14
    
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
        The chain identifier

    res_id : int
        The residue ID

    x : float
        The x coordinate of the atom

    y : float
        The y coordinate of the atom

    z : float
        The z coordinate of the atom

    segment : str
        The segment identifier

    Returns
    -------------
    line : str
        The fully-formatted ATOM line, as defined by the PDB specification.
    
    """
    
    ATOM      = build_section_string("ATOM",          6, 'L') # 1  - 6
    ATOM_IDX  = build_section_string(str(atom_index), 5, 'R') # 7  - 11
    BREAK_1   = " "                                           # 12
    ATOM_NAME = build_section_string(str(atom_name),  4, 'C') # 13 - 16
    ALTLOC    = " "                                           # 17
    RES_NAME  = build_section_string(str(res_name),   3, 'L') # 18 - 20 | note we don't actually expect the residue names to be anything other than a 3 letter code
    BREAK_2   = " "                                           # 21
    CHAIN     = build_section_string(str(chain),      1, 'L') # 22 
    RES_ID    = build_section_string(str(res_id),     4, 'R') # 23 - 26
    ICODE     = " "                                           # 27
    BREAK_3   = "   "                                         # 29 - 33
    X         = build_section_string(str(np.around(x,3)),     8, 'C')      # 31 - 38
    Y         = build_section_string(str(np.around(y,3)),     8, 'C')      # 39 - 46
    Z         = build_section_string(str(np.around(z,3)),     8, 'C')      # 47 - 54
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
        The index of the first atom in the system

    atom2 : int
        The index of the second atom in the system

    Returns
    -------------
    line : str
        The fully-formatted CONECT line, as defined by the PDB specification.


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
        The index of the atom in the system

    res_name : str
        The name of the residue

    chain : str
        The chain identifier

    res_id : int
        The residue ID

    Returns
    -------------
    line : str
        The fully-formatted TER line, as defined by the PDB specification.

    """    

    TER_SEC = "TER   "                                        # 1  - 6
    ATOM_IDX  = build_section_string(str(atom_index), 5, 'R') # 7  - 11
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

    Constructs a valid TER line for a PDB file - this creates an appropriate box size in 2 or 3 dimensions,
    where spacing defines how lattice sites relate to angstroms (i.e. by default each lattice site is 4
    angstroms apart --> inter-amino acid distance.


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
        Returns a fully-formatted valid CRYST line for a PDB file
    
    """    
 

    # the 4 here relfects the 4 anstroms per lattice site spacing being used..
    CRYST_SECT = "CRYST1"                                               # 1  - 6
    a          = build_section_string("%9.3f" % ( ((dimensions[0]*spacing)-spacing)), 9, 'R')  # 7  - 15
    b          = build_section_string("%9.3f" % ( ((dimensions[1]*spacing)-spacing)), 9, 'R')  # 16 - 24


    # set the third dimension depending on lattice type
    if len(dimensions) == 3:
        c      = build_section_string("%9.3f" % ( ((dimensions[2]*spacing)-spacing)), 9, 'R')  # 25 - 33
    else:
        c      = build_section_string("%9.3f" % 1.0, 9, 'R')            # 25 - 33

    alpha      = build_section_string("%7.2f" % 90.0, 7, 'R')           # 34 - 40
    beta       = build_section_string("%7.2f" % 90.0, 7, 'R')           # 41 - 47
    gamma      = build_section_string("%7.2f" % 90.0, 7, 'R')           # 48 - 54
    BREAK_1    = " "                                                    # 55 - 55
    sGroup     = build_section_string("P 1", 10,"L")                    # 56 - 66
    ZVALUE     = "    "                                                 # 67 - 70
    
    return build_line([CRYST_SECT, a, b, c, alpha, beta, gamma, BREAK_1, sGroup, ZVALUE],  [[1,6],[7,15],[16,24],[25,33],[34,40],[41,47],[48,54],[55,55],[56,66],[67,70]])+"\n"

    
    

    
    
                                          
        
            

        
        
        
            
    
    
