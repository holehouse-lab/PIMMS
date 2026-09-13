## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


###
### 
###

## Set of tools for analysis routines. ALL routines should 
# 1) not change any of the passed data

import math

import numpy as np
from numpy import linalg as LA
from scipy.spatial import ConvexHull # compute volume of clusters
from scipy.spatial import QhullError

from . import lattice_utils
from . import cluster_utils
from .latticeExceptions import AnalysisRoutineException


def get_inter_position_distance(P1, P2, dimensions, pbc_correction=True):
    """
    Returns the distance between two positions on the lattice (in real space)
    accounting for periodic boundary conditions.
   
    Routine optimized for a single distance (i.e. doesn't perform any of the
    setup/teardown used for vectorized implementations which are important
    when a set of positions are being compared)

    Parameters
    -----------------
    
    P1 : list
        A position list (e.g. a list of integers specifying the X/Y or X/Y/Z
        coordinates of a position)

    P2 : list
        A position list (e.g. a list of integers specifying the X/Y or X/Y/Z
        coordinates of a position)

    dimensions : list
        Defines the box size in 2 or 3 dimensions (a list of 2 or 3 ints). Its
        length sets how many coordinates of P1/P2 are used.

    pbc_correction : bool, optional
        Flag which if set to true means the minimum-image PBC correction is
        applied. Default is True.

    Returns
    -------
    float
        The (optionally PBC-corrected) Euclidean distance between ``P1`` and
        ``P2`` in real space.

    Notes
    -----
    This is deliberately written in plain Python arithmetic rather than with numpy.
    It is one of the most-called functions in PIMMS - the O(L^2) distance-map and
    internal-scaling analyses and the single-image seed search all funnel through it -
    and for a 2- or 3-element distance the numpy scalar machinery (two ``np.array``
    constructions plus ``np.power`` / ``np.sqrt`` dispatch per call) costs about 16x
    more than the arithmetic itself: ~6.5 us/call versus ~0.4 us. The results are
    bit-identical: integer coordinates stay exact through the squaring, and
    ``math.sqrt`` and ``np.sqrt`` are both the correctly-rounded IEEE-754 double
    square root.

    Prefer :func:`get_inter_position_distances` (or one of the vectorized callers in
    ``chain.py``) when you have many pairs to measure.
    """
    n_dim = len(dimensions)

    total = 0
    for idx in range(n_dim):
        d = P1[idx] - P2[idx]
        if d < 0:
            d = -d

        # minimum-image convention: a separation of more than half the box is
        # shorter the other way round
        if pbc_correction and d > dimensions[idx] * 0.5:
            d = dimensions[idx] - d

        total += d * d

    return math.sqrt(total)



def _minimum_image_lengths(delta, dims, pbc_correction=True):
    """Per-axis minimum-image separations from a raw coordinate difference array.

    Mirrors the per-axis logic of :func:`get_inter_position_distance` exactly (take the
    absolute separation, and where it exceeds half the box replace it with
    ``box - separation``), but over a whole array at once. Integer input stays integer,
    so the subsequent squaring is exact.

    Parameters
    ----------
    delta : numpy.ndarray
        Raw (unsigned or signed) coordinate differences with the per-axis
        components on the last axis, i.e. any shape ``(..., n_dim)``.

    dims : numpy.ndarray
        The box dimensions as a 1D array of length ``n_dim``, broadcast against
        the last axis of ``delta``.

    pbc_correction : bool, optional
        If True (the default) the minimum-image correction is applied. If False
        the absolute separations are returned unchanged.

    Returns
    -------
    numpy.ndarray
        Array with the same shape as ``delta`` holding the per-axis
        separations. The dtype follows ``delta`` and ``dims``, so integer input
        with integer dimensions gives an integer result.

    """
    out = np.abs(delta)
    if pbc_correction:
        out = np.where(out > 0.5 * dims, dims - out, out)
    return out


def get_distance_matrix(positions, dimensions, pbc_correction=True):
    """Full ``(L, L)`` minimum-image distance matrix for a set of positions.

    The vectorized equivalent of calling :func:`get_inter_position_distance` on every
    pair, and bit-identical to it. The per-pair Python double loop this replaces was
    the single largest cost in a PIMMS run with distance-map analysis switched on.

    Work is done in row blocks so that the ``(block, L, n_dim)`` intermediate stays
    bounded for long chains rather than scaling as ``L^2 * n_dim``.

    Parameters
    ----------
    positions : list or numpy.ndarray
        ``(L, n_dim)`` set of lattice positions.

    dimensions : list of int
        Box size in 2 or 3 dimensions.

    pbc_correction : bool, optional
        Apply the minimum-image correction (default True).

    Returns
    -------
    numpy.ndarray
        ``(L, L)`` symmetric matrix of pairwise distances.
    """
    pos = np.asarray(positions)
    dims = np.asarray(dimensions)
    n_pos = pos.shape[0]
    n_dim = pos.shape[1]

    out = np.empty((n_pos, n_pos), dtype=np.float64)

    # ~2M intermediate elements per block, and always at least one row
    block = max(1, int(2_000_000 // max(1, n_pos * n_dim)))

    for start in range(0, n_pos, block):
        stop = min(start + block, n_pos)
        delta = pos[start:stop, np.newaxis, :] - pos[np.newaxis, :, :]
        d = _minimum_image_lengths(delta, dims, pbc_correction)
        out[start:stop] = np.sqrt((d * d).sum(axis=-1))

    return out


def get_internal_scaling_profile(positions, dimensions, pbc_correction=True,
                                 return_squared=False):
    """Mean inter-bead distance as a function of sequence separation.

    Returns ``(gaps, means)`` for gaps ``1 .. L-1``. Each gap is measured with one
    vectorized pass over the ``L - gap`` pairs at that separation instead of a Python
    loop over pairs, and is bit-identical to the per-pair version. When
    ``return_squared`` is true, the mean of the *squared pair distances* is returned as
    a third array. This is not generally the square of ``means``: squaring after
    averaging loses the within-snapshot variance and does not give the RMS internal
    scaling profile.

    Parameters
    ----------
    positions : list or numpy.ndarray
        ``(L, n_dim)`` set of lattice positions, in chain order.

    dimensions : list of int
        Box size in 2 or 3 dimensions.

    pbc_correction : bool, optional
        Apply the minimum-image correction (default True).

    return_squared : bool, optional
        Also return the mean squared distance at every sequence separation. Default
        False, preserving the historical two-value return signature.

    Returns
    -------
    tuple
        ``(gaps, means)`` by default. If ``return_squared=True``, returns
        ``(gaps, means, mean_squares)``.
    """
    pos = np.asarray(positions)
    dims = np.asarray(dimensions)
    n_pos = pos.shape[0]

    gaps = []
    means = []
    mean_squares = []
    for gap in range(1, n_pos):
        d = _minimum_image_lengths(pos[gap:] - pos[:-gap], dims, pbc_correction)
        squared_distances = (d * d).sum(axis=-1)
        gaps.append(gap)
        means.append(np.mean(np.sqrt(squared_distances)))
        if return_squared:
            mean_squares.append(np.mean(squared_distances))

    if return_squared:
        return (gaps, means, mean_squares)
    return (gaps, means)


def get_inter_position_distances(P1s, P2s, dimensions, pbc_correction=True):
    """
    Returns the list of distances between lists of two positions on the lattice (in real space)
    accounting for periodic boundary conditions.

    Optimized for multiple values - vectorizes calculations.

    Parameters
    ----------
    P1s : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    P2s : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    pbc_correction : bool, optional
        If True (default), the minimum-image PBC correction is applied to each
        per-dimension separation before computing distances.

    Returns
    -------
    numpy.ndarray
        1D array of (optionally PBC-corrected) Euclidean distances, one per
        position pair.

    Raises
    ------
    AnalysisRoutineException
        If ``P1s`` and ``P2s`` do not have the same length.

    """

    # Check lists are the same length!
    if not len(P1s) == len(P2s):
        raise AnalysisRoutineException('Two lists of positions for distance analysis did not match one another in length')


    # extract box size in X/Y dimensions
    x_max = dimensions[0]
    y_max = dimensions[1]

    # convert to numpy arrays
    P1s = np.array(P1s)
    P2s = np.array(P2s)

    # get all the X/Y positions for list 1 and list 2
    P1_x = P1s.transpose()[0]
    P1_y = P1s.transpose()[1]
    P2_x = P2s.transpose()[0]
    P2_y = P2s.transpose()[1]

    # get vector of differences in X and Y dimensions
    x_dif = P1_x - P2_x
    y_dif = P1_y - P2_y

    # perform PBC correction for distances 
    if pbc_correction:
        # minimum image convention: where |d| > L/2, replace with L - |d|. The
        # selection mask must be identical on both sides (previously the RHS used
        # the always-empty mask abs(d) > L, which crashed on any over-half-box
        # separation due to a shape mismatch).
        x_mask = np.abs(x_dif) > 0.5*x_max
        x_dif[x_mask] = x_max - np.abs(x_dif[x_mask])
        y_mask = np.abs(y_dif) > 0.5*y_max
        y_dif[y_mask] = y_max - np.abs(y_dif[y_mask])
    
    # if we're in 3D do all the equivalent work for the 3D dimension (Z)
    if len(dimensions) == 3:

        z_max = dimensions[2]
        P1_z = P1s.transpose()[2]
        P2_z = P2s.transpose()[2]
        z_dif = P1_z - P2_z
        
        # PBC correction in Z
        if pbc_correction:
            z_mask = np.abs(z_dif) > 0.5*z_max
            z_dif[z_mask] = z_max - np.abs(z_dif[z_mask])
        
        distance_vector = np.sqrt(np.power(x_dif,2) + np.power(y_dif, 2) + np.power(z_dif, 2) )

    else:
        distance_vector = np.sqrt(np.power(x_dif,2) + np.power(y_dif, 2))

    return distance_vector


def get_cluster_distribution(lattice_grid, chainDict, hardwall=False):
    """
    Returns a list of lists, where each sublist contains the chainIDs associated 
    with a cluster. Cluster sublists are ordered from largest cluster to smallest.

    This is a computationally expensive algorithm that probably could be ported 
    into Cython at some point...

    Parameters
    ---------------

    lattice_grid : numpy.ndarray
        Standard lattice occupancy grid: a 2D or 3D integer array holding the
        chainID at each site (0 = solvent).

    chainDict : dict
        Standard dictionary mapping chainIDs to Chain objects.

    hardwall : bool, optional
        If True the box has hard walls, so chains are not connected through
        opposite faces of the box. Default is False (periodic boundaries).

    Returns
    -------
    list of list of int
        List of clusters, where each sublist contains the chainIDs of the
        chains in that connected component. Clusters are ordered from largest
        to smallest.

    """

    allChainIDs=[]
    for chainID in chainDict:
        allChainIDs.append(chainDict[chainID].chainID)

    # will contain lists of chains belonging to each cluster
    cluster_map = []

    # list of chains we've found so we only examine the minimum
    # number of clusters to get full coverage
    unfound_chains = set(allChainIDs)
    
    # until we've found all the chains...
    while len(unfound_chains) > 0:

        # take the first chainID from the set of unfound chains (next(iter(...)) picks
        # the same element as list(...)[0] but without materialising the whole set,
        # which made this loop O(n_chains^2))
        chainID = next(iter(unfound_chains))

        # get the set of chains in the connected component associated with chainID 

        cluster_members = lattice_utils.get_all_chains_in_connected_component(chainID, lattice_grid, chainDict, useChains=True, hardwall=hardwall)
        cluster_map.append(cluster_members)        
        
        # remove the found chains from the unfound chains set
        unfound_chains = unfound_chains.difference(cluster_members)

    # sort the cluster list - this sorts cluster map by the length of each sublist
    # and then reverses the order to get a list of sublists with the largest cluster
    # first. Finally, cycle through until we find a cluster smaller than the threshold, 
    # which point we're done
    clusters = sorted(cluster_map, key=len)[::-1]

    return clusters

def get_LR_cluster_distribution(latticeObject, hardwall=False, LR_table=None, SLR_table=None):
    """
    Returns a list of lists, where each sublist contains the chainIDs associated
    with a cluster. Cluster sublists are ordered from largest cluster to smallest.
    LR clusters are defined as clusters were interactions are through short-range
    OR long-range interactions (any short-range contact, or a Chebyshev-2/3 pair
    with nonzero LR/SLR interaction energy - see
    :func:`~pimms.lattice_utils.get_all_chains_in_long_range_cluster`).

    Parameters
    ---------------

    latticeObject : Lattice
        The PIMMS Lattice object. Its ``grid`` (the lattice grid) and ``chains``
        (mapping of chainIDs to chain objects) attributes are used, along with
        long-range interaction information, to build the long-range clusters.

    hardwall : bool, optional
        If True, do not connect chains through opposite faces of the box.
        Default is False (periodic boundaries).

    LR_table, SLR_table : numpy.ndarray or None, optional
        The Hamiltonian's LR and SLR residue interaction tables, each an
        ``(n_residues, n_residues)`` integer array indexed by residue integer
        code. When given, a Chebyshev-2/3 pair connects two chains only if its
        table entry is nonzero; when omitted (the default, None) the pair
        connects if both beads are LR-capable.

    Returns
    -------
    list of list of int
        List of clusters, where each sublist contains the chainIDs of the
        chains in that long-range connected component. Clusters are ordered
        from largest to smallest.

    """
    lattice_grid = latticeObject.grid
    chainDict = latticeObject.chains

    allChainIDs=[]    
    for chainID in chainDict:        
        allChainIDs.append(chainID)

    # will contain lists of chains belonging to each cluster
    cluster_map = []

    # list of chains we've found so we only examine the minimum
    # number of clusters to get full coverage
    unfound_chains = set(allChainIDs)
    
    # until we've found all the chains...
    while len(unfound_chains) > 0:

        # take the first chainID from the set of unfound chains (next(iter(...)) picks
        # the same element as list(...)[0] but without materialising the whole set,
        # which made this loop O(n_chains^2))
        chainID = next(iter(unfound_chains))

        # get the set of chains in the connected component associated with chainID 
        #cluster_members = lattice_utils.get_all_chains_in_connected_component(chainID, lattice_grid, chainDict, useChains=True)
        cluster_members = lattice_utils.get_all_chains_in_long_range_cluster(
            chainID, latticeObject, hardwall=hardwall,
            LR_table=LR_table, SLR_table=SLR_table)
        cluster_map.append(cluster_members)        
        
        # remove the found chains from the unfound chains set
        unfound_chains = unfound_chains.difference(cluster_members)

    # sort the cluster list - this sorts cluster map by the length of each sublist
    # and then reverses the order to get a list of sublists with the largest cluster
    # first
    clusters = sorted(cluster_map, key=len)[::-1]
    
    return clusters



def get_eigenvalues_of_the_T_matrix(positions, dimensions, pbc_correction=True):
    """
    Compute the eigenvalues and eigenvectors of the gyration (T) tensor.

    Builds the gyration tensor from the supplied positions relative to their
    (optionally PBC-corrected) center of mass, then diagonalizes it. The
    eigenvalues are the principal components used downstream to compute the
    radius of gyration and asphericity.

    Parameters
    ----------
    positions : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    pbc_correction : bool, optional
        If True (default), each position is PBC-corrected relative to the
        center of mass before contributing to the gyration tensor.

    Returns
    -------
    tuple
        ``(EIG, norm)`` where ``EIG`` is the array of eigenvalues of the
        gyration tensor and ``norm`` is the matrix of corresponding
        eigenvectors (as returned by ``numpy.linalg.eigh``).

    Notes
    -----
    The tensor is built with a single vectorized pass over the positions. This used to
    be a Python loop that called :func:`~pimms.lattice_utils.pbc_correct` and allocated
    an ``np.outer`` product for *every bead*, which made it one of the costliest parts
    of both the per-chain and the per-cluster property analyses.

    The gyration tensor is real and symmetric by construction, so ``eigh`` is used
    rather than the general ``eig``. Besides being faster, ``eig`` can return a complex
    array for a matrix that is only symmetric to within rounding, which would then
    propagate into the radius of gyration; ``eigh`` is guaranteed real. Every quantity
    derived downstream (see :func:`get_polymeric_properties`) is a symmetric function of
    the eigenvalues, so the different ordering ``eigh`` returns does not matter.
    """

    # NB: the gyration tensor here IS sensitive to whether pbc_correction is applied
    # (an incorrect wrap in a box barely larger than the cluster shifts far beads by a
    # full box length and inflates Rg). Callers that pass already-single-image
    # positions must therefore use pbc_correction=False (see
    # extract_cluster_polymeric_properties).
    pos  = np.asarray(positions, dtype=np.float64)
    dims = np.asarray(dimensions, dtype=np.float64)

    if pbc_correction:
        com = np.asarray(
            lattice_utils.center_of_mass_from_positions(
                positions, dimensions, on_lattice=False),
            dtype=np.float64)
        # vectorized pbc_correct(COM, pos): shift each position by one box length where
        # it sits more than half a box from the COM, so the whole set lies in one image
        diff = com - pos
        pos = pos + np.where(diff > dims / 2.0, dims, 0.0) - np.where(diff < -dims / 2.0, dims, 0.0)

        # The circular COM is used ONLY to select the periodic image. The gyration
        # tensor must be referenced to the ARITHMETIC mean of the reconstructed
        # single-image coordinates: the circular mean generally differs from it
        # for any non-symmetric configuration, and by the parallel-axis theorem
        # referencing the tensor to the wrong point inflated every PBC Rg by
        # exactly |mean - circular|^2 (a strictly one-sided bias, up to ~1e-2
        # relative in a box-7 system) and perturbed the asphericity. With this
        # line the PBC path satisfies the definition Rg^2 = mean |p - COM|^2 and
        # is bit-identical to the Cartesian path for non-straddling chains.
        com = pos.mean(axis=0)
    else:
        # Single-image and hardwall coordinates are ordinary Cartesian data.
        # A circular/PBC COM can lie a full box away and inflate the tensor even
        # when the caller explicitly disabled PBC correction.
        com = pos.mean(axis=0)

    delta = pos - com

    # T = <delta_i delta_j> over the beads
    T = (delta.T @ delta) / len(positions)

    # get the eigenvalues of the (symmetric) T matrix
    (EIG, norm) = LA.eigh(T)

    return (EIG, norm)

    


def get_polymeric_properties(positions, dimensions, pbc_correction=True):
    r"""
    Returns a list of polymeric properties calculated over the set of positions

    [0] - radius of gyration 
    [1] - asphericity

    Rg is defined as

    \sqrt(\dfrac{1}{N}\sum_{k=1}^N(r_k-r_{mean})^2)

    Where
    N = number of residues
    r_{mean} = mean residue position (Center of Mass)


    Parameters
    ----------
    positions : list of positions
        A list of positions, where each position is a 2-length or 3-length
        list specifying X/Y/[Z] coordinate positions on the lattice.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    pbc_correction : bool, optional
        Defines whether to perform PBC correction here (default True). For
        certain types of analysis (notably cluster analysis) the PBC correction
        is dealt with by the algorithms that construct the cluster, such that
        performing it again here is redundant (and generally not possible, as
        the snakesearch algorithm re-positions the cluster in terms of
        non-periodic space).

    Returns
    -------
    list of float
        A two-element list ``[rg, asph]`` where ``rg`` is the radius of
        gyration and ``asph`` is the asphericity (acylindricity in 2D), both
        derived from the gyration-tensor eigenvalues. Degenerate cases where
        ``rg ~ 0`` return an asphericity of 0.0.

    """

    n_dim = len(dimensions)


    # compute the eigenvalues and normal of the T matrix. NOTE - the function below USED
    # to be part of this function but we extracted it out
    (EIG, norm) = get_eigenvalues_of_the_T_matrix(positions, dimensions, pbc_correction)
        
    # Numerical tolerance for degenerate chains/clusters where Rg ~ 0.
    eps = 1e-12

    # if we're doing a 2D simulation
    if n_dim == 2:
        # radius of gyration from the gyration tensor
        rg2 = max(0.0, EIG[0] + EIG[1])
        rg = np.sqrt(rg2)

        # acylindiricity. For degenerate cases (rg == 0), define asphericity as 0.
        # NOTE on conventions: this 2D value is kappa = |l1 - l2| / (l1 + l2),
        # whereas the 3D branch below reports the relative shape anisotropy
        # kappa^2. The 2D analogue of the 3D quantity is therefore the SQUARE of
        # this value (kappa^2 = ((l1-l2)/(l1+l2))^2). The unsquared form is kept
        # for backwards compatibility with existing 2D analyses; square it to
        # compare against 3D asphericities (documented in docs/output_files.rst).
        if rg2 <= eps:
            asph = 0.0
        else:
            asph = abs(EIG[0] - EIG[1]) / rg2


    else:
        # radius of gyration from the gyration tensor
        rg_sum = max(0.0, EIG[0] + EIG[1] + EIG[2])
        rg = np.sqrt(rg_sum)

        # asphericity from the gyration tensor
        denom = np.power(rg_sum, 2)
        if denom <= eps:
            asph = 0.0
        else:
            asph = 1 - 3 * ((EIG[0] * EIG[1] + EIG[1] * EIG[2] + EIG[2] * EIG[0]) / denom)

    return [rg, asph]


def extract_positions_from_clusters(cluster_list, chainDict):
    """
    Function which takes a list of clusters (i.e. a list of lists, where 
    each sublist is a list of chainIDs in a specific cluster) and returns
    a list of lists of the same length where each sublist in the return list
    contains the positon of all residues in the cluster

    Parameters
    ----------
    cluster_list : list of list of int
        List of clusters, where each sublist is a list of chainIDs in that
        cluster.

    chainDict : dict
        Dictionary mapping each chainID to its chain object (each chain object
        must expose ``get_ordered_positions()``).

    Returns
    -------
    list of list
        List of the same length as ``cluster_list`` where each sublist contains
        the ordered positions of all residues belonging to that cluster.

    """

    return_list = []
    for cluster in cluster_list:
    
        sublist = []
        for chainID in cluster:
            sublist.extend(chainDict[chainID].get_ordered_positions())

        return_list.append(sublist)

    return return_list



def extract_cluster_polymeric_properties(cluster_position_list, dimensions=False):
    """
    Function which takes a list of cluster positions (i.e. a list of lists, where 
    each sublist is a list of positions associated with the residues in a specific cluster) 
    and returns a list of lists of the same length where each sublist in the return list
    contains the polymeric properties of the actual cluster.

    Parameters
    ----------
    cluster_position_list : list of list of positions
        List where each sublist is a list of positions; each sublist is its own
        cluster. NOTE that each cluster should exist within its own single-image
        convention, so that for each cluster the properties can be computed
        naively over those positions without any further PBC handling.

    dimensions : list of int or bool, optional
        Defines the dimensions of the lattice the positions sit on. If PBC
        correction has already been performed this can be left as ``False``
        (the default), in which case a per-cluster bounding box is computed
        dynamically (+10 beyond the largest value in each dimension).

    Returns
    -------
    list of list of float
        List of the same length as ``cluster_position_list`` where each entry
        is the ``[rg, asph]`` polymeric properties of the corresponding
        cluster (empty list if no clusters are supplied).

    """
    return_list = []

    # if no positions return empty list
    if len(cluster_position_list) == 0:
        return return_list

    local_dimensions = dimensions
    n_dim = len(cluster_position_list[0][0])

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            

        if dimensions is False:   # `is`: a numpy-array dimensions must not hit an ambiguous truth test
            # The positions are already single-image (snakesearch-corrected) per this
            # function's contract, so NO periodic correction must be applied when
            # computing their gyration tensor: a further PBC wrap in a box barely
            # larger than the cluster shifts any bead more than half a (tiny) box from
            # the circular-mean COM by a full box length, which silently inflates Rg
            # and distorts the asphericity of elongated/asymmetric clusters (up to
            # ~15% Rg error measured). We therefore pass pbc_correction=False; the
            # bounding box below is then irrelevant to the result and is kept only so
            # get_polymeric_properties has a dimensions list of the right length.
            if n_dim == 2:
                local_dimensions = [max(np.transpose(cluster)[0]+10), max(np.transpose(cluster)[1])+10]
            else:
                local_dimensions = [max(np.transpose(cluster)[0])+10, max(np.transpose(cluster)[1]+10), max(np.transpose(cluster)[2])+10]

            return_list.append(get_polymeric_properties(cluster, local_dimensions, pbc_correction=False))

        else:
            # The function's input contract is already-single-image cluster
            # coordinates regardless of whether the caller supplies the original box
            # dimensions. Re-applying minimum-image wrapping here can fold an extended
            # (but non-percolating) cluster and silently alter its tensor. Explicit
            # dimensions are metadata for dimensionality only; they do not change the
            # coordinate convention.
            return_list.append(get_polymeric_properties(
                cluster, local_dimensions, pbc_correction=False))

    return return_list



def correct_cluster_positions_to_single_image(cluster_position_list, dimensions):
    """
    Function which takes a list of cluster positions (i.e. a list of lists, where 
    each sublist is a list of positions associated with the residues in a specific cluster) 
    and for EACH CLUSTER re-configures the cluster position so the cluster is in its own single 
    periodic image.

    Parameters
    -------------
    cluster_position_list : list
        A list of lists; each sublist is the list of positions belonging to one
        cluster (where, in fact, each position is itself a list of 2 or 3
        integer coordinates).

    dimensions : list
        A list of 2 or 3 ints that defines the X/Y or X/Y/Z box dimensions.

    Returns
    ----------
    list
        List of the same length as ``cluster_position_list`` where each entry
        is the cluster's positions re-expressed in a single (non-periodic)
        image, as returned by
        ``cluster_utils.convert_positions_to_single_image_snakesearch`` with a
        ``space_threshold`` of 1.

    """

    num_clusters = len(cluster_position_list)
    
    return_list = []

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            

        # then perform single image PBC correction 
        return_list.append(cluster_utils.convert_positions_to_single_image_snakesearch(cluster, dimensions, space_threshold=1))

    return return_list


def _residue_types_for_positions(cluster, dimensions, type_grid):
    """
    Read the residue integer code sitting at each of a cluster's positions.

    Parameters
    ----------
    cluster : list or numpy.ndarray
        The cluster's positions, either wrapped or single-image (they are
        reduced mod the box here, which the single-image gather guarantees is
        congruent to the wrapped position).

    dimensions : list
        A list of 2 or 3 ints defining the X/Y or X/Y/Z box dimensions.

    type_grid : numpy.ndarray
        The lattice type grid, holding the residue integer code at each site.

    Returns
    -------
    numpy.ndarray
        A 1D int64 array of residue codes, one per position, in input order.

    """
    pos = np.asarray(cluster, dtype=np.int64)
    index = tuple(pos[:, d] % int(dimensions[d]) for d in range(len(dimensions)))
    return np.asarray(type_grid)[index].astype(np.int64)


def correct_LR_cluster_positions_to_single_image(cluster_position_list, dimensions,
                                                 type_grid=None, LR_table=None,
                                                 SLR_table=None):
    """
    Function which takes a list of cluster positions (i.e. a list of lists, where
    each sublist is a list of positions associated with the residues in a specific cluster)
    and for EACH CLUSTER re-configures the cluster position so the cluster is in its own single periodic image

    Like :func:`correct_cluster_positions_to_single_image`, but the gather walks
    the relation that DEFINES a long-range cluster rather than a plain distance
    rule: a Chebyshev-1 contact, or a Chebyshev-2 / Chebyshev-3 pair with a
    nonzero LR / SLR table entry. The type grid and both tables are needed for
    that; without them the old distance-only walk at ``space_threshold`` 3 is
    used.

    The distinction is not cosmetic. Long-range cluster MEMBERSHIP has always
    been the interaction-based relation, but the gather used to link any two
    beads within Chebyshev 3 of each other. As soon as a cluster's single-image
    extent reached within three sites of the box on some axis its two extreme
    beads became "neighbours" through the periodic face even when they carry no
    interaction at all, and part of the cluster was placed a box-length away.
    The ``LR_CLUSTER_*`` files then reported a torn, more compact object than
    the ``CLUSTER_*`` files did for the very same set of chains, and the
    percolation warning fired for a cluster that touches nothing.

    Parameters
    ----------
    cluster_position_list : list
        A list of lists; each sublist is a list of cluster positions (where
        each position is itself a 2- or 3-element list).

    dimensions : list
        A list of 2 or 3 elements defining the X/Y or X/Y/Z box dimensions.

    type_grid : numpy.ndarray or None, optional
        The lattice type grid, used to look up the residue code at each bead.
        Required (together with both tables) for the interaction-based gather.
        Default is None.

    LR_table : numpy.ndarray or None, optional
        The long-range residue interaction table. Default is None.

    SLR_table : numpy.ndarray or None, optional
        The super-long-range residue interaction table. Default is None.

    Returns
    -------
    list
        List of the same length as ``cluster_position_list`` where each entry
        is the cluster's positions re-expressed in a single (non-periodic)
        image, using a ``space_threshold`` of 3.

    """
    return_list = []

    use_tables = type_grid is not None and LR_table is not None and SLR_table is not None

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:

        # then perform single image PBC correction
        if use_tables and len(cluster) > 0:
            types = _residue_types_for_positions(cluster, dimensions, type_grid)
            return_list.append(cluster_utils.convert_positions_to_single_image_snakesearch(
                cluster, dimensions, space_threshold=3, types=types,
                LR_table=LR_table, SLR_table=SLR_table))
        else:
            return_list.append(cluster_utils.convert_positions_to_single_image_snakesearch(
                cluster, dimensions, space_threshold=3))

    return return_list


def flag_percolating_clusters(cluster_position_list, dimensions, space_threshold=1,
                              type_grid=None, LR_table=None, SLR_table=None):
    """
    Flag the gathered clusters that are connected to their own periodic image.

    Such a cluster is an unbounded object in the infinite periodic system the
    simulation represents, so it has no radius of gyration, no asphericity, no
    convex hull and no radial density profile - the gather still returns
    coordinates, but they are an arbitrary finite window cut out of an infinite
    object and the numbers depend on which bead the walk started from. Callers
    use this to blank those quantities rather than write a plausible-looking
    number.

    Parameters
    ----------
    cluster_position_list : list
        List of clusters, each an array of already-gathered (single-image)
        positions.

    dimensions : list
        A list of 2 or 3 ints defining the X/Y or X/Y/Z box dimensions.

    space_threshold : int, optional
        The per-dimension threshold the gather used (1 for contact clusters,
        3 for long-range clusters). Default is 1.

    type_grid : numpy.ndarray or None, optional
        The lattice type grid. Supplying this and both tables restricts the
        test to pairs that genuinely interact through the face, which is what
        makes it exact for long-range clusters. Default is None.

    LR_table : numpy.ndarray or None, optional
        The long-range residue interaction table. Default is None.

    SLR_table : numpy.ndarray or None, optional
        The super-long-range residue interaction table. Default is None.

    Returns
    -------
    list of bool
        One flag per cluster, True where the cluster percolates the box.

    """
    use_tables = type_grid is not None and LR_table is not None and SLR_table is not None

    flags = []
    for cluster in cluster_position_list:
        if len(cluster) == 0:
            flags.append(False)
            continue

        if use_tables:
            types = _residue_types_for_positions(cluster, dimensions, type_grid)
        else:
            types = None

        axis = cluster_utils.cluster_percolates(
            cluster, dimensions, space_threshold, types=types,
            LR_table=LR_table, SLR_table=SLR_table)
        flags.append(axis is not None)

    return flags



def compute_cluster_gross_properties(cluster_position_list):
    """
    Determines the volume of each cluster in a list of clusters. NOTE that positions
    here MUST have been corrected for PBC effects as the ConvexHull algorithm will
    calculate ASSUMES a single non-periodic image. This means that when you have clusters
    that wrap around the PBC the convex hull algorithm is calculating a single instance
    (i.e. using the boundaries as edges) so take care when extrapolating cluster volume 
    for such system spanning clusters.

    Using these positions and the Complex Hull algorithm we compute the volume, area
    and density of the cluster.
    
    Parameters
    -----------------
    cluster_position_list : list
        List of clusters, where each entry is an ``(N, n_dim)`` array (or list
        of lists) of the single-image positions of the beads in that cluster.


    Returns
    -------------
    list of lists
        Returns a list of lists, where each sublist contains three elements that reflect
        the cluster gross properties. Sublist indices match indices for cluster_position_list
        indices
    
        [0] - volume
        [1] - surface area
        [2] - density

    """
    
    return_list = []

    # for each set of positions associated with each cluster
    for cluster in cluster_position_list:            

        # run convex hull - if that throws an exception then
        # set everything to -1
        try:
            CH  = ConvexHull(cluster)
        except QhullError:
            vol = -1
            SA = -1
            den = -1
            return_list.append([vol, SA, den])
            continue
        

        # scipy's ConvexHull exposes .volume/.area on every version PIMMS supports
        # (scipy >= 1.9 is the install floor; the attributes date to 0.17), so the
        # old manual tetrahedron fallback was dead code and has been removed.
        #
        # NOTE the 2D convention: for a 2D hull scipy's .volume is the polygon
        # AREA and .area is the PERIMETER, so in 2D simulations CLUSTER_VOL.dat
        # holds areas, CLUSTER_AREA.dat holds perimeters, and den is beads per
        # unit area (see docs/output_files.rst).
        vol = CH.volume
        SA  = CH.area
        den = float(len(cluster))/vol # density in residues/VOLUME (2D: per area)

        # update lists
        return_list.append([vol, SA, den])

    return return_list
        



def compute_cluster_radial_density_profile(cluster_position_list, dimensions, minimum_cluster_size_in_beads=None, hardwall=False):
    """
    Compute the radial density profile of each cluster about its center of mass.

    For each cluster the density at "shell k" is the fraction of the lattice sites
    at Chebyshev (max-norm) distance k from the cluster centre of mass that are
    occupied by a bead - i.e. (beads at distance k) / (sites in shell k). The
    profile runs outward from the COM until every bead has been placed in a shell
    (or the box half-extent is reached), and short profiles are zero-padded to a
    common length.

    The profile STARTS AT SHELL 1 - the 26 (3D) or 8 (2D) sites immediately around
    the centre - so entry ``k`` of the returned list is the shell at Chebyshev
    distance ``k + 1``. Shell 0 is the single COM site itself and is never
    emitted; a bead sitting on it is excluded from the completion count rather
    than binned.

    Under periodic boundaries "distance k" means the MINIMUM-IMAGE Chebyshev
    distance, because the shell is a set of lattice sites and the site at
    ``COM - k`` is the site at ``COM + L - k``. Under a hardwall the box does
    not wrap, so plain Cartesian distances (and wall-clipped shells) are used.

    This is computed directly by binning each bead's Chebyshev distance from the COM
    (an O(num_beads) ``np.bincount``), rather than scanning every site of every
    concentric shell (which was O(offset_max ** n_dim) and dominated the cost for
    large clusters). It also fixes an off-by-one in the previous ring-scan, which
    additionally emitted a spurious shell at offset_max+1 (whose extent spills
    outside the box); profiles are now capped at ``offset_max`` entries as intended.

    Parameters
    ----------
    cluster_position_list : list of numpy.ndarray
        List of clusters; each entry is an array of lattice positions (each
        position being a 2- or 3-element coordinate) for that cluster. Positions
        are expected to be single-image (PBC-corrected) per cluster.

    dimensions : list of int
        Defines the box size in 2 or 3 dimensions (length 2 or 3).

    minimum_cluster_size_in_beads : int or None, optional
        If supplied, clusters with fewer beads than this threshold are skipped
        (no profile is emitted for them). Default is None (no filtering).

    hardwall : bool, optional
        If True the box has hard walls and a shell is normalised by the number
        of its sites that actually lie inside ``[0, dim)`` on every axis. Under
        periodic boundaries (the default) every shell site exists, so the full
        shell size is used. Without this a cluster touching a wall - the normal
        geometry of a wetting film or a wall-pinned droplet - had its
        wall-truncated shells divided by sites that do not exist and was
        reported too dilute (a completely full first shell read as 65 %
        occupied in the reproducer).

    Returns
    -------
    list of list of float
        One radial density profile per (non-skipped) cluster; each profile is a
        list of occupied-site fractions as a function of Chebyshev distance from the
        cluster center of mass, zero-padded to a uniform length. Entry ``k`` is the
        shell at Chebyshev distance ``k + 1``.
    """

    return_densities = []
    n_dim = len(dimensions)

    # Shells are only well-defined out to where they fit within the SMALLEST box
    # axis; beyond that a shell would wrap the short axis under periodic boundaries
    # (or run off the box under a hardwall). min(dimensions) keeps every shell
    # physical and profile lengths comparable across box shapes. For a cubic/square
    # box min == max, so this is unchanged there.
    offset_max = int((min(dimensions) / 2)) - 1

    for cluster_positions_nd in cluster_position_list:

        pts = np.asarray(cluster_positions_nd)
        num_beads = len(pts)

        # IF we've defined a smallest cluster worth computing for, skip small ones
        if minimum_cluster_size_in_beads is not None and num_beads < minimum_cluster_size_in_beads:
            continue

        # cluster COM position. The input is SINGLE-IMAGE (snakesearch) so its
        # centre is the plain arithmetic mean - NOT the periodic circular-mean COM,
        # which wraps back into [0, dim) and, for a cluster whose single image
        # extends past the box edge, lands a full box away from the beads (every
        # Chebyshev distance then exceeds the profile range and the whole profile
        # silently zeroes out).
        COM = np.rint(pts.mean(axis=0)).astype(int)

        # Chebyshev (max-norm) distance of every bead from the COM, then bin it:
        # counts[k] is the number of beads sitting in shell k.
        #
        # Under periodic boundaries shell k is the set of SITES {(COM + v) mod L,
        # |v|_inf = k}, and the denominator below already counts them that way,
        # so the numerator has to be counted the same way: a bead L-4 sites past
        # the COM in its single image occupies the shell-4 site at COM-4. Binning
        # by the raw single-image distance dropped such a bead from every shell
        # (the profile could then only read too dilute), which is only invisible
        # while every bead sits within offset_max of the COM. The shell cap
        # (2*offset_max+1 < min(dimensions)) is what guarantees the minimum-image
        # shells do not overlap, so the two conventions agree site for site.
        # A HARDWALL box does not wrap at all, so its branch keeps plain
        # Cartesian distances (and wall-clipped shells) - applying the periodic
        # metric there is a bug that has been fixed once already.
        delta = np.abs(pts - COM)
        if not hardwall:
            box = np.asarray(dimensions, dtype=delta.dtype)
            delta = delta % box
            delta = np.minimum(delta, box - delta)
        cheb = delta.max(axis=1).astype(np.int64)
        counts = np.bincount(cheb, minlength=offset_max + 1)

        # a bead sitting exactly on the (integer) COM is at shell 0 and can never be
        # found by shells >= 1, so it does not count toward completion.
        max_num_beads = num_beads - int(counts[0])

        # Walk outward shell by shell, out to offset_max - the largest shell that
        # still fits inside the smallest box axis (2*offset_max+1 <= min(dimensions)).
        # The density at shell k is (beads at Chebyshev distance k) / (number of
        # lattice sites in that shell), read off the histogram in O(1). Scanning
        # stops early once every findable bead has been placed.
        #
        # NB: the previous ring-scan loop had an off-by-one that also evaluated shell
        # offset_max+1, whose extent (2k+1 = min+1) spills outside the box; that
        # spurious out-of-bounds shell is no longer emitted, so profiles are now
        # capped at offset_max entries as intended.
        # number of lattice sites in the cube of Chebyshev radius k about the COM,
        # clipped to the box under a hardwall (every site exists under PBC); the
        # shell is the difference of consecutive cubes
        def _cube_sites(k):
            """
            Number of lattice sites in the cube of Chebyshev radius k about the COM.

            Under periodic boundaries every site of the cube exists, so this is
            just ``(2k+1) ** n_dim``. Under a hardwall the cube is clipped to
            the box on every axis, so a COM near a wall gives a smaller count.
            The shell size used for the density is the difference between two
            consecutive cubes.

            Parameters
            ----------
            k : int
                The Chebyshev radius of the cube (k=0 is the single COM site).

            Returns
            -------
            int
                The number of lattice sites in the cube, clipped to the box
                when ``hardwall`` is True.

            """
            if not hardwall:
                return (2 * k + 1) ** n_dim
            sites = 1
            for d in range(n_dim):
                lo = max(int(COM[d]) - k, 0)
                hi = min(int(COM[d]) + k, int(dimensions[d]) - 1)
                sites *= max(hi - lo + 1, 0)
            return sites

        ring_density = []
        found = 0
        for offset in range(1, offset_max + 1):
            occupied = int(counts[offset])
            total = _cube_sites(offset) - _cube_sites(offset - 1)
            ring_density.append(occupied / total if total > 0 else float('nan'))
            found += occupied

            # stop once every findable bead has been placed in a shell
            if found == max_num_beads:
                break

        # zero-pad short profiles to a common length
        if len(ring_density) < offset_max:
            ring_density.extend((offset_max - len(ring_density)) * [0])

        return_densities.append(ring_density)

    return return_densities

            


                    

                    
                
                

                    
            
            
            
        
        
