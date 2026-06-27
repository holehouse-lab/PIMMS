## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2026
## ...........................................................................
# 

##
## analysis_general
##
## Stand alone functions for reading/writing analysis output to disk
##

from .latticeExceptions import AcceptanceException
from . import CONFIG
import os


def _sorted_chain_types(IDtoType):
    """
    Return a deterministic, sorted ordering of the distinct chain types.

    Used so that the per-chain-type cluster output files are always generated
    in a reproducible order regardless of the (unordered) dictionary iteration.

    Parameters
    ----------
    IDtoType : dict
        Dictionary mapping each chainID to its chainType.

    Returns
    -------
    list
        Sorted list of the unique chain types present in ``IDtoType``.
    """
    return sorted(set(IDtoType.values()))


def _cluster_type_fractions(clusters, chain_types, IDtoType):
    """
    Precompute the per-cluster fraction of each chain type.

    For every cluster this computes the fraction of its member chains that
    belong to each chain type. Precomputing the fractions in a single pass
    avoids repeated O(types * cluster_size) scans when the per-chain-type
    output files are written.

    Parameters
    ----------
    clusters : list of list of int
        List of clusters, where each cluster is a list of chainIDs.

    chain_types : list
        Ordered list of the distinct chain types to compute fractions for.

    IDtoType : dict
        Dictionary mapping each chainID to its chainType.

    Returns
    -------
    dict
        Dictionary keyed by chain type. Each value is a list (one entry per
        cluster, in input order) giving the fraction of that cluster which is
        of the given chain type. Empty clusters contribute 0.0.
    """
    fractions = {chain_type: [] for chain_type in chain_types}

    for cluster in clusters:
        cluster_len = len(cluster)
        if cluster_len == 0:
            for chain_type in chain_types:
                fractions[chain_type].append(0.0)
            continue

        type_counts = {}
        for chainID in cluster:
            chain_type = IDtoType[chainID]
            type_counts[chain_type] = type_counts.get(chain_type, 0) + 1

        denom = float(cluster_len)
        for chain_type in chain_types:
            fractions[chain_type].append(type_counts.get(chain_type, 0) / denom)

    return fractions


def _prefixed_output_name(base_name, prefix):
    """
    Prepend a prefix to a filename while preserving its directory path.

    Splits ``base_name`` into directory and basename, prepends ``prefix`` to
    the basename only, and re-joins with the directory. This ensures prefixed
    output files are written into the same configured output directory rather
    than the current working directory.

    Parameters
    ----------
    base_name : str
        The original output path (may include a directory component).

    prefix : str
        String to prepend to the basename portion of ``base_name``.

    Returns
    -------
    str
        The prefixed path, with any directory component preserved.
    """
    directory = os.path.dirname(base_name)
    basename = os.path.basename(base_name)
    local_name = f"{prefix}{basename}"
    if directory:
        return os.path.join(directory, local_name)
    return local_name
    

#-----------------------------------------------------------------
#    
def write_energy(step, energy):
    """
    Append the current step and system energy to the ENERGY output file.

    Parameters
    ----------
    step : int
        Current simulation step.

    energy : float
        Current total system energy.

    Returns
    -------
    None
        Nothing is returned; a line is appended to the file defined by
        ``CONFIG.OUTNAME_ENERGY``.
    """

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    with open(CONFIG.OUTNAME_ENERGY, 'a') as fh:
        fh.write('%i\t%10.4f\n' % (step, energy))




#-----------------------------------------------------------------
#    
def write_clusters(step, clusters, IDtoType):
    """
    Write short-range cluster size distribution and composition for a step.

    Appends (1) the number of chains in each cluster, (2) the total number of
    clusters, and - when more than one chain type is present - (3) one file per
    chain type giving the fraction of each cluster made up of that chain type.

    Parameters
    ----------
    step : int
        Current simulation step.

    clusters : list of list of int
        List where each sublist represents a cluster containing the chainIDs
        of its member chains. For example ``[[1, 2], [3], [4]]`` represents a
        four-chain system where chains 1 and 2 are clustered together while
        chains 3 and 4 are isolated.

    IDtoType : dict
        Dictionary mapping each chainID to its chainType.

    Returns
    -------
    None
        Nothing is returned; results are appended to the cluster output files
        defined in ``CONFIG`` (``OUTNAME_CLUSTERS``, ``OUTNAME_NUM_CLUSTERS``
        and per-chain-type ``CHAIN_<type>_`` prefixed files).
    """

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the number of chains in each cluster
    with open(CONFIG.OUTNAME_CLUSTERS, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in clusters:
            fh.write('%i, '%(len(i)))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the number of clusters
    with open(CONFIG.OUTNAME_NUM_CLUSTERS, 'a') as fh:                
        fh.write('%i\t%i\n'%(step,len(clusters)))

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Write the fraction of each cluster which is of a certain
    # chain type (provides info on heterogenous clusters)

    # get the number of different types
    n_types = len(set(IDtoType.values()))

    if n_types == 1:
        # if  there's only one type of chain this analysis is moot
        return

    else:
        
        # else get a list of the different chain types
        typelist = _sorted_chain_types(IDtoType)
        fractions = _cluster_type_fractions(clusters, typelist, IDtoType)

        
        # for each different chain type cycle over the clusters and identify
        # what fraction of each cluster is current iteration chainType
        # - means we generate $n_types different output files. The combinatorics
        # can be derived from these output files :-)
        for chainType in typelist:

            outname = _prefixed_output_name(CONFIG.OUTNAME_CLUSTERS, f"CHAIN_{chainType}_")
            with open(outname, 'a') as fh:
                for frac in fractions[chainType]:
                    fh.write('%2.4f, ' % frac)
                fh.write('\n')


#-----------------------------------------------------------------
#    
def write_LR_clusters(step, clusters, IDtoType):
    """
    Write long-range (LR) cluster size distribution and composition for a step.

    Identical in behaviour to :func:`write_clusters` but writes to the
    long-range cluster output files. LR clusters are defined as clusters
    connected through short-range OR long-range interactions.

    Parameters
    ----------
    step : int
        Current simulation step.

    clusters : list of list of int
        List where each sublist represents a cluster containing the chainIDs
        of its member chains. For example ``[[1, 2], [3], [4]]`` represents a
        four-chain system where chains 1 and 2 are clustered together while
        chains 3 and 4 are isolated.

    IDtoType : dict
        Dictionary mapping each chainID to its chainType.

    Returns
    -------
    None
        Nothing is returned; results are appended to the LR cluster output
        files defined in ``CONFIG`` (``OUTNAME_LR_CLUSTERS``,
        ``OUTNAME_NUM_LR_CLUSTERS`` and per-chain-type ``CHAIN_<type>_``
        prefixed files).
    """

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the number of chains in each cluster
    with open(CONFIG.OUTNAME_LR_CLUSTERS, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in clusters:
            fh.write('%i, '%(len(i)))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the number of clusters
    with open(CONFIG.OUTNAME_NUM_LR_CLUSTERS, 'a') as fh:                
        fh.write('%i\t%i\n'%(step,len(clusters)))

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Write the fraction of each cluster which is of a certain
    # chain type (provides info on heterogenous clusters)

    # get the number of different types
    n_types = len(set(IDtoType.values()))

    if n_types == 1:
        # if  there's only one type of chain this analysis is moot
        return

    else:
        
        # else get a list of the different chain types
        typelist = _sorted_chain_types(IDtoType)
        fractions = _cluster_type_fractions(clusters, typelist, IDtoType)

        
        # for each different chain type cycle over the clusters and identify
        # what fraction of each cluster is current iteration chainType
        # - means we generate $n_types different output files. The combinatorics
        # can be derived from these output files :-)
        for chainType in typelist:

            outname = _prefixed_output_name(CONFIG.OUTNAME_LR_CLUSTERS, f"CHAIN_{chainType}_")
            with open(outname, 'a') as fh:
                for frac in fractions[chainType]:
                    fh.write('%2.4f, ' % frac)
                fh.write('\n')
                        

def write_cluster_properties(step, cluster_polymeric_properties_list, cluster_size_list, cluster_radial_density):
    """
    Write the per-cluster polymeric and gross properties for a single step.

    Appends one line per output file for the current step, distributing the
    supplied per-cluster properties across the dedicated short-range cluster
    property files (Rg, asphericity, volume, surface area, density and radial
    density profile).

    Parameters
    ----------
    step : int
        Current simulation step.

    cluster_polymeric_properties_list : list of list
        One entry per cluster; each entry is ``[Rg, asphericity]`` as returned
        by the polymeric-property analysis.

    cluster_size_list : list of list
        One entry per cluster; each entry is ``[volume, surface_area, density]``
        as returned by the gross-property analysis.

    cluster_radial_density : list of list of float
        One entry per cluster; each entry is the radial density profile
        (a list of densities as a function of distance from the cluster COM).

    Returns
    -------
    None
        Nothing is returned; results are appended to the cluster property
        files defined in ``CONFIG`` (``OUTNAME_CLUSTER_RG``,
        ``OUTNAME_CLUSTER_ASPH``, ``OUTNAME_CLUSTER_VOL``,
        ``OUTNAME_CLUSTER_AREA``, ``OUTNAME_CLUSTER_DENSITY`` and
        ``OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE``).
    """

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the rg of each cluster!
    with open(CONFIG.OUTNAME_CLUSTER_RG, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in cluster_polymeric_properties_list:
            fh.write('%2.4f, '%(i[0]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the asphericity of each cluster!
    with open(CONFIG.OUTNAME_CLUSTER_ASPH, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in cluster_polymeric_properties_list:
            fh.write('%2.4f, '%(i[1]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the volume of each cluster!
    with open(CONFIG.OUTNAME_CLUSTER_VOL, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in cluster_size_list:
            fh.write('%2.4f, '%(i[0]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the surface area of each cluster!
    with open(CONFIG.OUTNAME_CLUSTER_AREA, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in cluster_size_list:
            fh.write('%2.4f, '%(i[1]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the density of each cluster!
    with open(CONFIG.OUTNAME_CLUSTER_DENSITY, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in cluster_size_list:
            fh.write('%2.4f, '%(i[2]))
        fh.write('\n')


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the radial density profile for every cluster
    with open(CONFIG.OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE, 'a') as fh:          

        idx=1        
        for cluster in cluster_radial_density:
            fh.write('%i, C%i, ' % (step, idx))
            for i in cluster:
                fh.write('%1.4f, '%(i))
            fh.write('\n')
            idx=idx+1



def write_LR_cluster_properties(step, LR_cluster_polymeric_properties_list, LR_cluster_size_list, LR_cluster_radial_density):
    """
    Write the per-cluster polymeric and gross properties for long-range clusters.

    Identical in behaviour to :func:`write_cluster_properties` but writes to the
    long-range (LR) cluster property output files.

    Parameters
    ----------
    step : int
        Current simulation step.

    LR_cluster_polymeric_properties_list : list of list
        One entry per LR cluster; each entry is ``[Rg, asphericity]``.

    LR_cluster_size_list : list of list
        One entry per LR cluster; each entry is ``[volume, surface_area, density]``.

    LR_cluster_radial_density : list of list of float
        One entry per LR cluster; each entry is the radial density profile.

    Returns
    -------
    None
        Nothing is returned; results are appended to the LR cluster property
        files defined in ``CONFIG`` (``OUTNAME_LR_CLUSTER_RG``,
        ``OUTNAME_LR_CLUSTER_ASPH``, ``OUTNAME_LR_CLUSTER_VOL``,
        ``OUTNAME_LR_CLUSTER_AREA``, ``OUTNAME_LR_CLUSTER_DENSITY`` and
        ``OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE``).
    """

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the rg of each cluster!
    with open(CONFIG.OUTNAME_LR_CLUSTER_RG, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in LR_cluster_polymeric_properties_list:
            fh.write('%2.4f, '%(i[0]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the asphericity of each cluster!
    with open(CONFIG.OUTNAME_LR_CLUSTER_ASPH, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in LR_cluster_polymeric_properties_list:
            fh.write('%2.4f, '%(i[1]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the volume of each cluster!
    with open(CONFIG.OUTNAME_LR_CLUSTER_VOL, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in LR_cluster_size_list:
            fh.write('%2.4f, '%(i[0]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the surface area of each cluster!
    with open(CONFIG.OUTNAME_LR_CLUSTER_AREA, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in LR_cluster_size_list:
            fh.write('%2.4f, '%(i[1]))
        fh.write('\n')

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the density of each cluster!
    with open(CONFIG.OUTNAME_LR_CLUSTER_DENSITY, 'a') as fh:          

        fh.write('%i, ' % step)
        for i in LR_cluster_size_list:
            fh.write('%2.4f, '%(i[2]))
        fh.write('\n')


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # write the radial density profile for every cluster
    with open(CONFIG.OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE, 'a') as fh:          

        idx=1        
        for cluster in LR_cluster_radial_density:
            fh.write('%i, C%i, ' % (step, idx))
            for i in cluster:
                fh.write('%1.4f, '%(i))
            fh.write('\n')
            idx=idx+1


#-----------------------------------------------------------------
#    
def write_internal_scaling(mean_IS, mean_IS_squared, prefix=False):
    """
    Write the mean internal scaling and mean internal scaling squared profiles.

    Each profile is written (overwriting any existing file) as a two-column
    table of ``gap`` (sequence separation, starting at 1) versus the mean value
    at that gap.

    Parameters
    ----------
    mean_IS : list of float
        Mean internal scaling value for each sequence-separation gap, ordered
        by increasing gap.

    mean_IS_squared : list of float
        Mean internal scaling squared value for each sequence-separation gap,
        ordered by increasing gap.

    prefix : str or bool, optional
        If ``False`` (default) the standard output filenames are used. If a
        string is supplied it is prepended to the output basenames (e.g. to
        separate per-chain-type output).

    Returns
    -------
    None
        Nothing is returned; results are written to the files defined by
        ``CONFIG.OUTNAME_INTERNAL_SCALING`` and
        ``CONFIG.OUTNAME_INTERNAL_SCALING_SQUARED`` (optionally prefixed).
    """

    ## First normal internal scaling
    # count refers to the IS gap
    count = 1
    if prefix is False:
        FN = CONFIG.OUTNAME_INTERNAL_SCALING
    else:
        FN = _prefixed_output_name(CONFIG.OUTNAME_INTERNAL_SCALING, prefix)

    with open(FN, 'w') as fh:
        for i in mean_IS:
            fh.write('%i\t%4.4f\n' % (count, i))
            count=count+1


    ## Next internal scaling squared
    # count refers to the IS gap
    count = 1
    if prefix is False:
        FN = CONFIG.OUTNAME_INTERNAL_SCALING_SQUARED
    else:
        FN = _prefixed_output_name(CONFIG.OUTNAME_INTERNAL_SCALING_SQUARED, prefix)

    with open(FN, 'w') as fh:
        for i in mean_IS_squared:
            fh.write('%i\t%4.4f\n' % (count, i))
            count=count+1


#-----------------------------------------------------------------
#    
def write_scaling_information(all_nu, all_R0, prefix=False):
    """
    Write paired scaling exponent and prefactor values to file.

    Writes (overwriting any existing file) a two-column table where each row is
    a ``nu`` (scaling exponent) and ``R0`` (prefactor) pair extracted from the
    internal-scaling fitting.

    Parameters
    ----------
    all_nu : list of float
        Scaling exponent values.

    all_R0 : list of float
        Prefactor values, paired element-wise with ``all_nu``.

    prefix : str or bool, optional
        If ``False`` (default) the standard output filename is used. If a
        string is supplied it is prepended to the output basename.

    Returns
    -------
    None
        Nothing is returned; results are written to the file defined by
        ``CONFIG.OUTNAME_SCALING_INFORMATION`` (optionally prefixed).

    Raises
    ------
    ValueError
        If ``all_nu`` and ``all_R0`` do not have the same length.
    """

    if len(all_nu) != len(all_R0):
        raise ValueError("all_nu and all_R0 must have the same length")

    if prefix is False:
        FN = CONFIG.OUTNAME_SCALING_INFORMATION
    else:
        FN = _prefixed_output_name(CONFIG.OUTNAME_SCALING_INFORMATION, prefix)

    with open(FN, 'w') as fh:
        for i in range(0,len(all_nu)):
            fh.write('%4.4f\t%4.4f\n' % (all_nu[i], all_R0[i]))




#-----------------------------------------------------------------
#    
def write_distance_map(dMap, prefix=False):
    """
    Write a square inter-residue distance map to file.

    Writes (overwriting any existing file) the full ``seqlen x seqlen`` matrix
    as tab-separated rows, one row per source residue.

    Parameters
    ----------
    dMap : numpy.ndarray
        Square ``(seqlen, seqlen)`` distance map matrix.

    prefix : str or bool, optional
        If ``False`` (default) the standard output filename is used. If a
        string is supplied it is prepended to the output basename.

    Returns
    -------
    None
        Nothing is returned; results are written to the file defined by
        ``CONFIG.OUTNAME_DMAP`` (optionally prefixed).
    """

    seqlen = dMap.shape[0]

    if prefix is False:
        FN = CONFIG.OUTNAME_DMAP
    else:
        FN = _prefixed_output_name(CONFIG.OUTNAME_DMAP, prefix)

    with open(FN, 'w') as fh:

        for i in range(0, seqlen):
            for j in range(0, seqlen):
                fh.write('%4.4f\t' % dMap[i][j])

            fh.write('\n')


#-----------------------------------------------------------------
#    
def write_radius_of_gyration(step, RG_list):
    """
    Append the instantaneous radius of gyration values for a single step.

    Parameters
    ----------
    step : int
        Current simulation step.

    RG_list : list of float
        Radius of gyration for each chain (or analysed object) at this step.

    Returns
    -------
    None
        Nothing is returned; a line is appended to the file defined by
        ``CONFIG.OUTNAME_RG``.
    """

    with open(CONFIG.OUTNAME_RG, 'a') as fh:
        fh.write('%i\t' %(step))
        for i in RG_list:
            fh.write('%3.3f\t' % i)

        fh.write('\n')


#-----------------------------------------------------------------
#    
def write_asphericity(step, asphericity_list):
    """
    Append the instantaneous asphericity values for a single step.

    Parameters
    ----------
    step : int
        Current simulation step.

    asphericity_list : list of float
        Asphericity for each chain (or analysed object) at this step.

    Returns
    -------
    None
        Nothing is returned; a line is appended to the file defined by
        ``CONFIG.OUTNAME_ASPH``.
    """

    with open(CONFIG.OUTNAME_ASPH, 'a') as fh:
        fh.write('%i\t' %(step))
        for i in asphericity_list:
            fh.write('%3.3f\t' % i)

        fh.write('\n')


#-----------------------------------------------------------------
#    
def write_end_to_end(step, e2e_list):
    """
    Append the instantaneous end-to-end distance values for a single step.

    Parameters
    ----------
    step : int
        Current simulation step.

    e2e_list : list of float
        End-to-end distance for each chain at this step.

    Returns
    -------
    None
        Nothing is returned; a line is appended to the file defined by
        ``CONFIG.OUTNAME_E2E``.
    """

    with open(CONFIG.OUTNAME_E2E, 'a') as fh:
        fh.write('%i\t' %(step))
        for i in e2e_list:
            fh.write('%3.3f\t' % i)

        fh.write('\n')


#-----------------------------------------------------------------
#    
def write_residue_residue_distance(step, R2R_info, all_data):
    """
    Append specific residue-residue distance measurements for a single step.

    For each monitored residue pair a single line is written containing the
    step, the two residue indices, and the per-chain distances measured for
    that pair.

    Parameters
    ----------
    step : int
        Current simulation step.

    R2R_info : list of sequence of int
        List of residue pairs; each element is an indexable pair where
        element 0 and element 1 are the two residue indices.

    all_data : list of list of float
        Distance data paired element-wise with ``R2R_info``; each element is a
        list of distances (e.g. one per chain) for the corresponding pair.

    Returns
    -------
    None
        Nothing is returned; lines are appended to the file defined by
        ``CONFIG.OUTNAME_R2R``.

    Raises
    ------
    ValueError
        If ``R2R_info`` and ``all_data`` do not have the same length.
    """

    if len(R2R_info) != len(all_data):
        raise ValueError("R2R_info and all_data must have the same length")

    with open(CONFIG.OUTNAME_R2R, 'a') as fh:
        
        # cycle through each pair writing a single line with
        # STEP | PAIR1 PAIR2 | RG1 RG2 .... RGN
        for pair, data in zip(R2R_info, all_data):

            # write the step
            fh.write('%i\t' %(step))

            fh.write('%i\t' %(pair[0]))
            fh.write('%i\t' %(pair[1]))
            

            for i in data:
                fh.write('%3.3f\t' % i)

            fh.write('\n')


#-----------------------------------------------------------------
#    
def write_acceptance_statistics(step, acceptanceObject):
    """
    Write move attempt, acceptance, and total move statistics for a step.

    Writes the per-move attempt counts and accepted counts to two separate
    files (so attempted and accepted moves can be compared directly), and the
    cumulative total number of attempted MC moves (across all sub-loops) to a
    third file.

    Parameters
    ----------
    step : int
        Current simulation step.

    acceptanceObject : AcceptanceCalculator
        Object tracking move statistics. Must expose ``move_count`` and
        ``accepted_count`` (indexable by move code) and
        ``alt_Markov_chain_moves``.

    Returns
    -------
    None
        Nothing is returned; lines are appended to the files defined by
        ``CONFIG.OUTNAME_MOVES``, ``CONFIG.OUTNAME_ACCEPTANCE`` and
        ``CONFIG.OUTNAME_TOTAL_MOVES``.

    Raises
    ------
    AcceptanceException
        If the number of tracked moves is not the expected hard-coded value
        (15), which guards the explicit move list used to compute total moves.
    """
    n_moves = len(acceptanceObject.move_count)

    # first write out the moves whcih were attempted
    with open(CONFIG.OUTNAME_MOVES, 'a') as fh:
        fh.write('%i\t' %(step))

        for move in range(1, n_moves):
            m_count = acceptanceObject.move_count[move]
            fh.write('%i\t' % m_count)
            
        fh.write('\n')

    # next write out the accepted moves
    with open(CONFIG.OUTNAME_ACCEPTANCE, 'a') as fh:
        fh.write('%i\t' %(step))

        for move in range(1, n_moves):
            a_count = acceptanceObject.accepted_count[move]
            fh.write('%i\t' % a_count)
            
        fh.write('\n')

    # finally write out the TOTAL moves so far...
    with open(CONFIG.OUTNAME_TOTAL_MOVES, 'a') as fh:
        
        # first count explicit moves (note n_moves is the true number of moves +1)
        if not n_moves == 15:
            print(n_moves)
            raise AcceptanceException('\n\nWhen trying to compute total moves found a hard-coded bug - this is probably because you tried to add a new move and not update this part of the code. You must explicitly define which moves use a sub-MC chain and which do not\n\n')

        # total attempted MC moves across all sub-loops (kept in sync with
        # AcceptanceCalculator.total_attempted_moves, which performance reporting
        # uses; the n_moves guard above protects this hard-coded move list).
        total_moves = 0
        for move in [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14]:
            total_moves = total_moves + acceptanceObject.move_count[move]
        total_moves = total_moves + acceptanceObject.alt_Markov_chain_moves

        fh.write('%i\t%i\n' %(step, total_moves))


#-----------------------------------------------------------------
#    
def write_performance(step, eq_string, steps_per_second, overall_moves_per_second, time_elapsed, time_remaining):
    """
    Function which writes out the time per step (as taken
    at some specific step in the simulation) for convenient
    comparison of simulation efficiency.

    Parameters
    ----------
    step : int
        The current simulation step number

    eq_string : str
        A string which indicates if the simulation is in
        equilibrium (E) or production (P) phase.

    steps_per_second : float
        The number of outer master-loop steps per second the
        simulation is currently running at.

    overall_moves_per_second : float
        The number of individual MC accept/reject moves per second
        across ALL sub-loops (megamove substeps + TSMMC excursion
        substeps), i.e. the true MC throughput.

    time_elapsed : str
        The time elapsed in the simulation so far.

    time_remaining : str
        The time remaining in the simulation.

    Returns
    -------
    None

        No return value, but writes to the file defined in
        CONFIG.OUTNAME_PERFORMANCE.
    
    """
    
    with open(CONFIG.OUTNAME_PERFORMANCE, 'a') as fh:

        # format the strings and pad them to a fixed length so we generate
        # nice formatted OUTNAME_PERFORMANCE files. Note the padding size
        # here is defined by the header, which in turn is defined when the file
        # is initially created at the start of the Simulation class.
        sps_string = '%.2f' % steps_per_second
        sps_pad = " "*(21-len(sps_string))

        ov_string = '%.2f' % overall_moves_per_second
        ov_pad = " "*(27-len(ov_string))

        te_pad = " "*(21-len(time_elapsed))

        fh.write(f"{step}\t{eq_string}\t{sps_string}{sps_pad}\t{ov_string}{ov_pad}\t{time_elapsed}{te_pad}\t{time_remaining}\n")
        


#-----------------------------------------------------------------
#    
def write_quench_file(step, temperature, energy):
    """
    Append the current step, temperature, and energy to the quench file.

    Parameters
    ----------
    step : int
        Current simulation step.

    temperature : float
        Current simulation temperature.

    energy : float
        Current total system energy.

    Returns
    -------
    None
        Nothing is returned; a line is appended to the file defined by
        ``CONFIG.QUENCHFILE_NAME``.
    """

    with open(CONFIG.QUENCHFILE_NAME, 'a') as fh:
        fh.write('%i\t%3.2f\t%10.4f\n' % (step, temperature, energy))                        



