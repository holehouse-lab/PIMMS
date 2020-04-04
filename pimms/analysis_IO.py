## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab 
## Copyright 2015 - 2020
## ...........................................................................
# 

##
## analysis_general
##
## Stand alone functions for reading/writing analysis output to disk
##

from .latticeExceptions import AcceptanceException
from . import CONFIG
    

#-----------------------------------------------------------------
#    
def write_energy(step, energy):
    """
    Function to write the current energy and step
    to the output ENERGY file.
    """

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    with open(CONFIG.OUTNAME_ENERGY, 'a') as fh:
        fh.write('%i\t%10.4f\n' % (step, energy))




#-----------------------------------------------------------------
#    
def write_clusters(step, clusters, IDtoType):
    """
    Function to write the number of clusters
    and the distribution of chains across
    clusters out at the current step


    Arguments:

    step [int]
    Current simulation step

    clusters [list of lists of integers]
    List, where each sublist represents a cluster containing the
    sublists chains (sublist-list items are chainIDs). e.g.
    [[1,2],[3],[4]] would represent a system with four chains were
    chains 1 and 2 are clustered together while 3 and 4 are in
    isolation

    IDtoType [dict]
    Dictionary which maps a chainID to a chainType
    
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
        typelist = list(set(IDtoType.values()))

        
        # for each different chain type cycle over the clusters and identify
        # what fraction of each cluster is current iteration chainType
        # - means we generate $n_types different output files. The combinatorics
        # can be derived from these output files :-)
        for chainType in typelist:
            
            with open("CHAIN_%i_"%(chainType) + CONFIG.OUTNAME_CLUSTERS, 'a') as fh:          

                
                # for each cluster
                for cluster in clusters:

                    # initialze the fraction of the cluster which is of chainType
                    chainsOfType=0.0

                    # for each chain in a cluster
                    for chainID in cluster:
                        if IDtoType[chainID] == chainType:                            
                            chainsOfType=chainsOfType+1

                    # write what fraction of that cluster was chainType
                    fh.write('%2.4f, '%(chainsOfType/float(len(cluster))))
                fh.write('\n')


#-----------------------------------------------------------------
#    
def write_LR_clusters(step, clusters, IDtoType):
    """
    Function to write the number of clusters
    and the distribution of chains across
    clusters out at the current step


    Arguments:

    step [int]
    Current simulation step

    clusters [list of lists of integers]
    List, where each sublist represents a cluster containing the
    sublists chains (sublist-list items are chainIDs). e.g.
    [[1,2],[3],[4]] would represent a system with four chains were
    chains 1 and 2 are clustered together while 3 and 4 are in
    isolation

    IDtoType [dict]
    Dictionary which maps a chainID to a chainType
    
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
        typelist = list(set(IDtoType.values()))

        
        # for each different chain type cycle over the clusters and identify
        # what fraction of each cluster is current iteration chainType
        # - means we generate $n_types different output files. The combinatorics
        # can be derived from these output files :-)
        for chainType in typelist:
            
            with open("CHAIN_%i_"%(chainType) + CONFIG.OUTNAME_LR_CLUSTERS, 'a') as fh:          

                
                # for each cluster
                for cluster in clusters:

                    # initialze the fraction of the cluster which is of chainType
                    chainsOfType=0.0

                    # for each chain in a cluster
                    for chainID in cluster:
                        if IDtoType[chainID] == chainType:                            
                            chainsOfType=chainsOfType+1

                    # write what fraction of that cluster was chainType
                    fh.write('%2.4f, '%(chainsOfType/float(len(cluster))))
                fh.write('\n')
                        

def write_cluster_properties(step, cluster_polymeric_properties_list, cluster_size_list, cluster_radial_density):

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

    ## First normal internal scaling
    # count refers to the IS gap
    count = 1
    if prefix is False:
        FN = CONFIG.OUTNAME_INTERNAL_SCALING
    else:
        FN = prefix+CONFIG.OUTNAME_INTERNAL_SCALING

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
        FN = prefix+CONFIG.OUTNAME_INTERNAL_SCALING_SQUARED

    with open(FN, 'w') as fh:
        for i in mean_IS_squared:
            fh.write('%i\t%4.4f\n' % (count, i))
            count=count+1


#-----------------------------------------------------------------
#    
def write_scaling_information(all_nu, all_R0, prefix=False):

    if prefix is False:
        FN = CONFIG.OUTNAME_SCALING_INFORMATION
    else:
        FN = prefix+CONFIG.OUTNAME_SCALING_INFORMATION

    with open(FN, 'w') as fh:
        for i in range(0,len(all_nu)):
            fh.write('%4.4f\t%4.4f\n' % (all_nu[i], all_R0[i]))




#-----------------------------------------------------------------
#    
def write_distance_map(dMap, prefix=False):

    seqlen = dMap.shape[0]

    if prefix is False:
        FN = CONFIG.OUTNAME_DMAP
    else:
        FN = prefix+CONFIG.OUTNAME_DMAP

    with open(FN, 'w') as fh:

        for i in range(0, seqlen):
            for j in range(0, seqlen):
                fh.write('%4.4f\t' % dMap[i][j])

            fh.write('\n')


#-----------------------------------------------------------------
#    
def write_radius_of_gyration(step, RG_list):
    """
    Writes the instantenous Rg for the values in the Rg_list

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
    Writes the instantenous asphericity values
    from the list provided

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
    Writes the instantenous end-to-end distance
    for the values in the e2e_list

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
    Writes the specific residue-residue distances out
    

    """

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
    Writes information on move attempts and acceptance out
    to two different files (defined in CONFIG). This makes
    it easy to compare attempted moves with accepted moves.
    
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
        if not n_moves == 14:
            print(n_moves)
            raise AcceptanceException('\n\nWhen trying to compute total moves found a hard-coded bug - this is probably because you tried to add a new move and not update this part of the code. You must explicitly define which moves use a sub-MC chain and which do not\n\n')
            
        total_moves = 0
        for move in [1,2,3,4,5,6,7,8,12,13]:
            total_moves = total_moves + acceptanceObject.move_count[move]

        # finally update from the alt Markov chain moves
        total_moves = total_moves + acceptanceObject.alt_Markov_chain_moves 

        fh.write('%i\t%i\n' %(step, total_moves))


#-----------------------------------------------------------------
#    
def write_time_per_step(step, step_interval, dt):
    """
    Function which writes out the time per step (as taken
    at some specific step in the simulation) for convenient
    comparison of simulation efficiency
    
    """
    
    with open(CONFIG.OUTNAME_PERFORMANCE, 'a') as fh:
        fh.write('%i\t%i\t%s\n' % (step, step_interval, str(dt)))


#-----------------------------------------------------------------
#    
def write_quench_file(step, temperature, energy):
    """
    Function which writes out the quench info 

    """

    with open(CONFIG.QUENCHFILE_NAME, 'a') as fh:
        fh.write('%i\t%3.2f\t%10.4f\n' % (step, temperature, energy))                        



