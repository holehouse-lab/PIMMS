## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2017
## ...........................................................................

import time

from .latticeExceptions import ParameterFileException
from . import file_utilities
from . import IO_utils
from . import CONFIG

#-----------------------------------------------------------------
#
def parse_energy(filename):
    """
    Function which reads in an energy parameter file and returns
    the interaction matrix as a redundant dictionary of dictionaries
    and the set of non-redundant particles defined therein.

    Each line in the parameter file defines the interaction between
    two residues, in the format
    
    A B X

    or

    A B X Y

    In both cases X defines the short range interaction energy
    between A and B. This occurs when A and B are adjacent to
    one another on a lattice

    Y defines the long range interaction energy (i.e. electrostatics
    though you could use it for anything) which occurs over one 
    lattice site.
    

    """

    with open(filename, 'r') as fh:
        contents = fh.readlines()

    # these two dictionaries becomes the interaction matrix for short range
    # and long range interactions
    energy_pairs = {}
    LR_energy_pairs = {}
    SLR_energy_pairs = {}
    
    # this becomes the set of all residue types defined
    # including solvent
    non_redundant_particles    = set()
    non_redundant_LR_particles = set()
    long_range_entries         = []
    
    for line in contents:
        
        # if it's a comment line skip
        if file_utilities.is_comment_line(line):
            continue
                
        # remove comment section
        un_comment = file_utilities.remove_comments(line)

        split_line = line.split()

        # reset the long-range flag
        LR_flag = False

        # skip lines which define angle penalties - we deal with them later
        if split_line[0] == "ANGLE_PENALTY" or split_line[0] == "ANGLE_PENALTY_T_NORM":
            continue

        # ------------------------------------------------
        # From here on out we're defining the pairwise residue-residue interactions
        # If we have a line which is not either 3 or 4 separate values
        linesplitlen = len(split_line)
        if linesplitlen != 3 and linesplitlen != 4 and linesplitlen !=5:
            raise ParameterFileException('ERROR: Trying to parse line ["%s"] - after comment parsing get ["%s"] - can not be broken into the format <residue> <residue> <energy>' %(line, un_comment))
            
        # first deal with all the short-range interaction stuff
        P1     = split_line[0].strip()
        P2     = split_line[1].strip()
        ENERGY = float(split_line[2].strip())

        # non reundant particles is a set of all the particle names
        non_redundant_particles.add(P1)
        non_redundant_particles.add(P2)

        # then if this line contains a long-range component - we'll deal with this after dealing with the short range 
        # components
        if linesplitlen == 4 or linesplitlen == 5:

            # IF we passed SLR value then use, else set to zero
            if linesplitlen == 4:
                long_range_entries.append([P1,P2, float(split_line[3].strip()), 0.0])

            else:
                long_range_entries.append([P1,P2, float(split_line[3].strip()), float(split_line[4].strip())])
            
        # the following code ensures we build the fully redundant square interaction matrix
        if P1 in energy_pairs:

            if P2 in energy_pairs[P1]:
                raise ParameterFileException('ERROR: Trying to assign the [%s -- %s] energy (%s) but this value was already set - parameter files MUST be non-redundant' % (P1, P2, ENERGY))
            else:
                energy_pairs[P1][P2] = ENERGY
                                                                                        
        else:
            energy_pairs[P1]={}
            energy_pairs[P1][P2] = ENERGY

        if P2 in energy_pairs:

            if P1 in energy_pairs[P2] and not (P1 == P2):
                raise ParameterFileException('ERROR: Trying to assign the [%s -- %s] energy (%s) but this value was already set - parameter files MUST be non-redundant' % (P1, P2, ENERGY))
            else:
                energy_pairs[P2][P1] = ENERGY

        else:
            energy_pairs[P2]={}
            energy_pairs[P2][P1] = ENERGY

    # ************************************************
    # END OF SHORT RANGE INTERACTIONS


    ## Long range interactions are then parsed from the list generated when we were reading the 
    ## lines in the previous section
    for LR_pair in long_range_entries:
        P1 = LR_pair[0]
        P2 = LR_pair[1]
        LR_ENERGY = LR_pair[2]
        SLR_ENERGY = LR_pair[3]
            
        non_redundant_LR_particles.add(P1)
        non_redundant_LR_particles.add(P2)

        # the following code ensures we build the fully redundant square interaction matrix
        if P1 in LR_energy_pairs:
            if P2 in LR_energy_pairs[P1]:
                raise ParameterFileException('ERROR: Trying to assign the [%s -- %s] energy (%s) but this value was already set - parameter files MUST be non-redundant' % (P1, P2, LR_ENERGY))
            else:
                LR_energy_pairs[P1][P2] = LR_ENERGY
                SLR_energy_pairs[P1][P2] = SLR_ENERGY
                                                                                        
        else:
            LR_energy_pairs[P1]={}
            LR_energy_pairs[P1][P2] = LR_ENERGY

            SLR_energy_pairs[P1]={}
            SLR_energy_pairs[P1][P2] = SLR_ENERGY


        if P2 in LR_energy_pairs:

            if P1 in LR_energy_pairs[P2] and not (P1 == P2):
                raise ParameterFileException('ERROR: Trying to assign the [%s -- %s] energy (%s) but this value was already set - parameter files MUST be non-redundant' % (P1, P2, LR_ENERGY))
            else:                
                LR_energy_pairs[P2][P1] = LR_ENERGY
                SLR_energy_pairs[P2][P1] = SLR_ENERGY

        else:
            LR_energy_pairs[P2]={}
            LR_energy_pairs[P2][P1] = LR_ENERGY

            SLR_energy_pairs[P2]={}
            SLR_energy_pairs[P2][P1] = SLR_ENERGY


    # ----------------------------------------------------------------------------------------
    # End of file reading assignment loop
                
    ##
    ## SOLVENT SOLUTE SANITY CHECKS
    ##
    
    # check we have at least one solvation interaction - we check for the full redundancy in a second, 
    # but knowing that there's at least one solute-solvent interaction energy term defined is useful
    # for the next stage. NOTE that LR interactions cannot have long range solute-solvent interactions 
    # (this assumption is hard-coded in later)
    if '0' not in non_redundant_particles:
        raise ParameterFileException('ERROR: None of the interactions in the parameter keyfile define any solvation interactions\nPLEASE ensure each partice type defined has a solvation energy defined\nThis should look like\n\nX 0 <SOLVATION ENERGY>') 

    # check IF we defined solvent-solvent interaction it was zero - if we didn't define set it to zero!
    if '0' in list(energy_pairs['0'].keys()):
        if not energy_pairs['0']['0'] == 0:
            raise ParameterFileException('ERROR: PIMMS does not support an interaction scheme where the SOLVENT-SOLVENT interaction energy is not zero')
    else:
        energy_pairs['0']['0'] = 0.0

    if '0' in LR_energy_pairs:
        raise ParameterFileException('ERROR: PIMMS does not support long range solvent-solute interactions')
                    
    ##
    ## INTERACTION REDUNDANCY CHECKS
    ##
        
    # check we defined the full non-redundant matrix for short range interactions
    for i in non_redundant_particles:
        for j in non_redundant_particles:
            if j not in list(energy_pairs[i].keys()):
                raise ParameterFileException("ERROR: The interaction between %s and %s is not defined, suggesting the parameter file is missing a pair" %(j,i))

    # NOTE - WE DO NOT make this check for long-range interactions because we don't want to force a convention on the form of 
    # those interactions. As an example, one might define LR interactions between two charged residues and between two aromatic
    # residues. HOWEVER, this does not mean charge-to-aromatic residues experience LR interactions, but for the sake of the 
    # underlying code these interactions ARE defined but set to zero. 
    for i in non_redundant_LR_particles:
        for j in non_redundant_LR_particles:
            if j not in list(LR_energy_pairs[i].keys()):

                print("WARNING: Long-range energy check: No defined long-range interaction between [%s] and [%s]. Setting to 0.0" % (i,j))
                LR_energy_pairs[i][j] = 0.0
                SLR_energy_pairs[i][j] = 0.0

                if i in list(LR_energy_pairs[j].keys()) and not (i == j):
                    raise ParameterFileException("ERROR: The interaction between %s and %s was not defined, BUT the reverse (%s to %s) was - this is indicative of a bug in the parameter file parser" %(i,j,j,i))
                    
                LR_energy_pairs[j][i] = 0.0
                SLR_energy_pairs[j][i] = 0.0

    ## WARNING:
    #
    # this ensures that in the non-redundant list of residue names
    # the solvent residue '0' is always first
    # we later make the assumption that residue_names[0] = solvent 
    #
    # *** PLEASE DO NOT BREAK THIS!! ***
    # (yes, convention assumptions are a dangerous game)
    #
    ## 

    residue_names = []
    LR_residue_names = []
    residue_names.append('0')
    for i in non_redundant_particles:
        if i == '0':
            continue
        else:
            residue_names.append(i)

    # The LR_residue_names are much less important so we don't
    # need to worry about any kind of convention there 
    for i in non_redundant_LR_particles:
        LR_residue_names.append(i)
            
    # print to STDOUT    
    # This works but is kind of a pain, so to avoid masses of STDOUT have turned this off...
    # print_energy_parameter_summary(energy_pairs, residue_names, LR_energy_pairs, LR_residue_names, filename, SLR_energy_pairs)

    print("")
    print("Writing the complete set of parameters used in this simulation out to: %s" % (CONFIG.OUTPUT_USED_PARAMETER_FILE))
    print("")
    # Finally SAVE this parameter file to the current directory, so we can ALWAYS fully reproduce any simulation
    contents.insert(0,'## This is a copy of the parameter file used for the simulation\n')
    contents.insert(0,'## This file was generated by PIMMS on %s\n' % (time.strftime("%c")))
    contents.insert(0,'##\n')
    contents.insert(0,'##\n')
    IO_utils.write_list_to_file(contents,CONFIG.OUTPUT_USED_PARAMETER_FILE)



    return (energy_pairs, residue_names, LR_energy_pairs, LR_residue_names, SLR_energy_pairs)




#-----------------------------------------------------------------
#
def parse_angles(filename, temperature=False):
    """
    Reads in a parameter file and constructs an angle penalty dictionary. We allow the definition of two types of angle penalties

    ANGLE_PENALTY         lines in the parameter file define an ABSOLUTE value 
    ANGLE_PENALTY_T_NORM
    
    """
    with open(filename, 'r') as fh:
        contents = fh.readlines()

    angle_dict = {}
    angle_dict_multiplier = {}
    
    for line in contents:

        # if it's a comment line skip
        if file_utilities.is_comment_line(line):
            continue
                
        # remove comment section
        un_comment = file_utilities.remove_comments(line)

        # split the line up 
        split_line = line.split()

        # if we found an angle penalty line using absolute values
        if split_line[0] == 'ANGLE_PENALTY':

            # check it's formatted in a valid way
            if len(split_line) != 5:
                raise ParameterFileException('ERROR: malformatted ANGLE_PENALTY line found [%s]' % (line))

            # set the residue name and try and extract residue-specific angle penalty values
            resname = split_line[1]
            try:
                AP1     = float(split_line[2])
                AP2     = float(split_line[3])
                AP3     = float(split_line[4])
            except ValueError:
                raise ParameterFileException('Unable to convert one or more values into ANGLE_PENALTY values for line [ %s ]' % line)
                
            # if we already had an entry for this residue crash!
            if resname in angle_dict:
                raise ParameterFileException('ERROR: Multiple ANGLE_PENALTY definitions for residue %s' % resname)
            
            angle_dict[resname] = [AP1,AP2,AP3]

        # if we found an angle penalty line using T normalized values (where the penality units are now in kT, assuming
        # k = 1 which is gene
        elif split_line[0] == 'ANGLE_PENALTY_T_NORM':
            
            print("Found T_NORM angle definition")

            if temperature is False:
                raise ParameterFileException('ERROR: T_NORM angle penality was found in parameter file, yet no temperature provided for parsing. THIS IS A BUG.')

            # check it's formatted in a valid way
            if len(split_line) != 5:
                raise ParameterFileException('ERROR: malformatted ANGLE_PENALTY line found [%s]' % (line))

            # set the residue name and try and extract residue-specific angle penalty values
            resname = split_line[1]
            try:
                AP1_M  = float(split_line[2])
                AP2_M  = float(split_line[3])
                AP3_M  = float(split_line[4])

                AP1    = AP1_M*temperature
                AP2    = AP2_M*temperature
                AP3    = AP3_M*temperature
            except ValueError:
                raise ParameterFileException('Unable to convert one or more values into ANGLE_PENALTY values for line [ %s ]' % line)
                
            # if we already had an entry for this residue crash!
            if resname in angle_dict:
                raise ParameterFileException('ERROR: Multiple ANGLE_PENALTY definitions for residue %s' % resname)
            
            angle_dict[resname] = [AP1,AP2,AP3]
            
            
        
    return angle_dict


#-----------------------------------------------------------------
#
def print_energy_parameter_summary(energy_pairs, residue_names, LR_energy_pairs, LR_residue_names, filename, SLR_energy_pairs):    

    print("===================================================")
    print("||           ENERGY PARAMETER SUMMARY            ||")
    print("===================================================")
    print("Parameter file used:")
    print(filename)
    print("")
    
    print("Interaction energy table to be used printed below")
    for R1 in residue_names:
        print("| %s |" % R1)
        for R2 in residue_names:
            if R1 in LR_residue_names and R2 in LR_residue_names:
                print("%s ::: %s = %3.2f,  %3.2f,  %3.2f ****" % (R1, R2, energy_pairs[R1][R2], LR_energy_pairs[R1][R2], SLR_energy_pairs[R1][R2]))
            else:
                print("%s ::: %s = %3.2f" % (R1, R2, energy_pairs[R1][R2]))
        print("")

#-----------------------------------------------------------------
#
def write_angle_parameter_summary(angle_dict, filename):    

    with open(CONFIG.OUTPUT_FULL_ANGLE_POTENTIAL,'w') as fh:
        fh.write("===================================================\n")
        fh.write("||            ANGLE PARAMETER SUMMARY            ||\n")
        fh.write("===================================================\n")
        fh.write("Parameter file used:\n")
        
        fh.write("\n")
                 
        fh.write("Angle penalties for each residues to be used printed below\n")
        for R1 in angle_dict:
            fh.write("%s -> %3.2f, %3.2f, %3.2f\n" % (R1, angle_dict[R1][0],angle_dict[R1][1],angle_dict[R1][2]))
        fh.write("\n")
    

        


            
            
            
            

        
