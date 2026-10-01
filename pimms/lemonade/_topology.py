## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................

"""
Trajectory topology: the fixed (time-independent) description of which beads
belong to which chain, each chain's sequence and type, and the bead-type alphabet.

Stored columnar for speed: a CSR-style ``offsets`` array delimits each chain's
contiguous block of beads, so any per-chain reduction is a single
``numpy.add.reduceat`` and any type filter is a boolean mask over one ``(n_beads,)``
array - no per-chain Python objects.
"""

import warnings

import numpy as np

from pimms import CONFIG

# PIMMS writes each bead as a PDB residue whose name is one_to_three(type_char):
# known amino acids use their standard 3-letter code, everything else becomes
# 'XX<c>'. Invert that so we can read the 1-letter bead type straight back.
_THREE_TO_ONE = {three: one for one, three in CONFIG.ONE_TO_THREE.items()}

# PIMMS writes one PDB chain identifier per chainType from A-Z, a-z, 0-9 and gives
# every type past the 62nd the last one (pdb_utils.build_pdb_file), so a PDB that
# uses all 62 may hold several real types under one label.
_N_PDB_CHAIN_IDS = 62


def pdb_chain_labels(md_topology):
    """Return the PDB chain identifier of every chain, or ``None`` if any is blank.

    Parameters
    ----------
    md_topology : mdtraj.Topology
        Topology read from the PDB.

    Returns
    -------
    list or None
        One chain identifier per chain, in chain order, or ``None`` if the
        topology has no chains or any chain has a missing or blank identifier
        (mdtraj reports a blank chain column as ``' '``), in which case the PDB
        carries no chain-type information.
    """
    labels = [getattr(chain, "chain_id", None) for chain in md_topology.chains]
    if labels and all(label is not None and str(label).strip() for label in labels):
        return labels
    return None


def three_to_one(resname):
    """Decode a PDB residue name back to its 1-letter PIMMS bead type.

    Parameters
    ----------
    resname : str
        Residue name as written in the PDB. Standard amino acids use their
        3-letter code; any other bead type is written as ``'XX<c>'``.

    Returns
    -------
    str
        The 1-letter bead type. A name that is neither a known 3-letter code nor
        an ``XX``-padded type is returned unchanged.
    """
    if resname in _THREE_TO_ONE:
        return _THREE_TO_ONE[resname]
    stripped = resname.lstrip("X")
    return stripped if stripped else resname


class Topology:
    """Columnar chain/bead topology for a trajectory.

    Attributes
    ----------
    n_chains, n_beads : int
    offsets : (n_chains + 1,) int64
        Bead-index boundaries; chain ``c`` owns beads ``offsets[c]:offsets[c+1]``.
    lengths : (n_chains,) int64
        Number of beads in each chain.
    chain_types : (n_chains,) int32
        Integer type label per chain.
    sequences : list[str]
        1-letter bead sequence per chain.
    bead_chainid : (n_beads,) int32
        Owning chain index for every bead.
    bead_codes : (n_beads,) int8
        Index into ``alphabet`` for every bead's type.
    alphabet : list[str]
        Sorted unique 1-letter bead types.
    """

    __slots__ = ("offsets", "lengths", "chain_types", "sequences",
                 "bead_chainid", "bead_codes", "alphabet")

    def __init__(self, sequences, chain_types=None):
        """Build the columnar arrays from the per-chain sequences.

        Parameters
        ----------
        sequences : iterable of str
            One 1-letter bead sequence per chain, in trajectory chain order.
            Chain ``c`` of the trajectory owns ``len(sequences[c])`` beads.
        chain_types : iterable of int, optional
            Integer type label per chain, one entry per sequence (default
            ``None``, which groups identical sequences into the same type).

        Raises
        ------
        ValueError
            If a sequence is not a string, if any sequence is empty (a
            zero-length chain breaks the ``reduceat`` reductions in
            ``_analysis``), if ``chain_types`` does not have one entry per
            chain, or if a chain type is not a non-negative int32 integer.
        """
        self.sequences = list(sequences)
        if any(not isinstance(s, str) for s in self.sequences):
            raise ValueError("Topology: every chain sequence must be a string")
        # a zero-length chain has no beads: the CSR reduceat machinery in
        # _analysis would silently return garbage (reduceat on an empty segment
        # yields the NEXT element, then /0 gives inf/nan) and the unwrap kernel
        # defends against it separately - reject it at construction instead
        if any(len(s) == 0 for s in self.sequences):
            raise ValueError("Topology: chain sequences must be non-empty")
        self.lengths = np.array([len(s) for s in self.sequences], dtype=np.int64)
        self.offsets = np.zeros(len(self.sequences) + 1, dtype=np.int64)
        np.cumsum(self.lengths, out=self.offsets[1:])

        if chain_types is None:
            # group identical sequences into the same type
            seen = {}
            chain_types = []
            for s in self.sequences:
                chain_types.append(seen.setdefault(s, len(seen)))
        chain_types = list(chain_types)
        if len(chain_types) != len(self.sequences):
            raise ValueError("Topology: chain_types must contain one entry per chain")
        normalized_types = []
        limits = np.iinfo(np.int32)
        for value in chain_types:
            if (isinstance(value, (bool, np.bool_)) or
                    not isinstance(value, (int, np.integer)) or
                    value < 0 or value > limits.max):
                raise ValueError(
                    "Topology: chain types must be non-negative int32 integers")
            normalized_types.append(int(value))
        self.chain_types = np.asarray(normalized_types, dtype=np.int32)

        n_beads = int(self.offsets[-1])
        self.bead_chainid = np.empty(n_beads, dtype=np.int32)
        for c in range(len(self.sequences)):
            self.bead_chainid[self.offsets[c]:self.offsets[c + 1]] = c

        all_chars = "".join(self.sequences)
        self.alphabet = sorted(set(all_chars))
        code = {ch: i for i, ch in enumerate(self.alphabet)}
        self.bead_codes = np.fromiter((code[ch] for ch in all_chars),
                                      dtype=np.int8, count=n_beads)

    # -- construction ------------------------------------------------------
    @classmethod
    def from_mdtraj(cls, md_topology):
        """Build from an mdtraj Topology (each bead is one residue/atom).

        PIMMS assigns the same PDB chain identifier to every chain of one
        ``chainType``.  Preserve that information instead of grouping by sequence:
        different sequences may intentionally share a type, and identical sequences
        may intentionally have different types.

        Parameters
        ----------
        md_topology : mdtraj.Topology
            Topology read from the PDB. Each of its chains becomes one PIMMS
            chain and each residue one bead; the chain identifier, where the PDB
            carries one, becomes the chain type.

        Returns
        -------
        Topology
            The new topology. Chain types come from the PDB chain identifiers
            when every chain has a non-blank one, and are grouped by sequence
            otherwise.
        """
        sequences = []
        for chain in md_topology.chains:
            sequences.append("".join(three_to_one(atom.residue.name)
                                     for atom in chain.atoms))
        labels = pdb_chain_labels(md_topology)
        if labels is not None:
            seen = {}
            chain_types = [seen.setdefault(label, len(seen)) for label in labels]
            return cls(sequences, chain_types=chain_types)
        return cls(sequences)

    def with_keyfile_types(self, chain_specs, labelled=False):
        """Return a copy whose chain types follow the keyfile CHAIN order.

        ``chain_specs`` is the parser's ``[[count, sequence], ...]`` list. If it is
        consistent with this topology (same chain count, same per-chain sequences
        and, for a ``labelled`` topology, the same partition of chains into types)
        the authoritative per-spec type index is used; otherwise the topology is
        returned unchanged.

        Parameters
        ----------
        chain_specs : list of [int, str]
            The keyfile CHAIN (plus EXTRA_CHAIN) specification in keyfile order:
            each entry is a ``[count, sequence]`` pair, and its position in the
            list is the chain type index.
        labelled : bool, optional
            Whether this topology's chain types came from the PDB chain
            identifiers (default ``False``, i.e. they are a sequence grouping and
            carry no type information of their own). PIMMS writes one identifier
            per chainType, so a labelled topology already holds the run's real
            partition of chains into types, and the keyfile types are then only
            accepted if they reproduce it: each PDB label must hold exactly one
            keyfile type and each keyfile type exactly one PDB label. The one
            relaxation is a PDB that uses all 62 identifiers, where PIMMS gives
            every type past the 62nd the last identifier; there a keyfile type
            may split a label but may still not straddle two.

        Returns
        -------
        Topology
            A new topology with keyfile-ordered chain types, or ``self`` if the
            specification does not match this topology's chains.

        Notes
        -----
        The partition check is what catches a keyfile whose CHAIN lines list the
        right sequences in an order that is not the trajectory's chain order.
        The ``keyfile_used.kf`` of a restart run is exactly that: it writes one
        CHAIN line per chain type, while the trajectory keeps the snapshot's
        chains first and appends any EXTRA_CHAIN chains at the end, so when two
        types share a sequence a sequence-only check expanded the lines onto the
        wrong chains without complaint.
        """
        expanded = []                       # (sequence, type) per chain, in order
        for type_idx, (count, seq) in enumerate(chain_specs):
            for _ in range(int(count)):
                expanded.append((seq, type_idx))
        if len(expanded) != len(self.sequences):
            return self
        if any(seq != self.sequences[i] for i, (seq, _t) in enumerate(expanded)):
            return self
        keyfile_types = [t for _s, t in expanded]
        if labelled:
            pdb_to_keyfile = {}
            keyfile_to_pdb = {}
            for pdb_type, key_type in zip(self.chain_types.tolist(), keyfile_types):
                pdb_to_keyfile.setdefault(pdb_type, set()).add(key_type)
                keyfile_to_pdb.setdefault(key_type, set()).add(pdb_type)
            # a keyfile type spread over two PDB labels is always a mismatch
            if any(len(v) != 1 for v in keyfile_to_pdb.values()):
                return self
            # a PDB label holding two keyfile types is a mismatch unless the PDB
            # ran out of identifiers, in which case the keyfile is the only way
            # to split the merged types
            if (len(pdb_to_keyfile) < _N_PDB_CHAIN_IDS and
                    any(len(v) != 1 for v in pdb_to_keyfile.values())):
                return self
        return Topology(self.sequences, chain_types=keyfile_types)

    def matches_keyfile_composition(self, chain_specs):
        """Does the keyfile describe these chain types, in any order?

        Compares the multiset of per-type ``(count, sequence)`` pairs of this
        topology against the keyfile's, ignoring the order of both chains and
        CHAIN lines. A match means the keyfile and the trajectory describe the
        same system even though the CHAIN lines do not expand onto the chains in
        order - which is what the ``keyfile_used.kf`` of a restart run looks
        like, since it lists one CHAIN line per type while the trajectory has
        the EXTRA_CHAIN chains appended after the snapshot's.

        Parameters
        ----------
        chain_specs : list of [int, str]
            The keyfile CHAIN (plus EXTRA_CHAIN) specification, each entry a
            ``[count, sequence]`` pair. Lines with a zero count describe no
            chains and are ignored.

        Returns
        -------
        bool
            ``True`` if every chain type of this topology holds one sequence
            and the ``(count, sequence)`` multisets agree, ``False`` otherwise.
        """
        mine = {}
        for c, t in enumerate(self.chain_types.tolist()):
            seqs, n = mine.get(t, (set(), 0))
            seqs.add(self.sequences[c])
            mine[t] = (seqs, n + 1)
        if any(len(seqs) != 1 for seqs, _n in mine.values()):
            return False
        ours = sorted((n, next(iter(seqs))) for seqs, n in mine.values())
        theirs = sorted((int(count), str(seq)) for count, seq in chain_specs
                        if int(count) > 0)
        return ours == theirs

    # -- convenience -------------------------------------------------------
    @property
    def n_chains(self):
        """int : Number of chains in the trajectory."""
        return len(self.sequences)

    @property
    def n_beads(self):
        """int : Total number of beads across all chains."""
        return int(self.offsets[-1])

    @property
    def n_atoms(self):
        """int : Deprecated alias of :attr:`n_beads`.

        PIMMS is a coarse-grained model and its particles are beads; the name was
        inherited from the PDB/XTC vocabulary and is kept only so scripts written
        against earlier builds keep running.
        """
        warnings.warn("n_atoms is deprecated, use n_beads (PIMMS has beads, not atoms)",
                      DeprecationWarning, stacklevel=2)
        return self.n_beads

    def type_mask(self, bead_type):
        """Boolean mask selecting the beads of one bead type.

        Parameters
        ----------
        bead_type : str
            A 1-letter bead type. A type that is not in ``alphabet`` selects
            nothing rather than raising.

        Returns
        -------
        numpy.ndarray
            ``(n_beads,)`` bool mask, ``True`` at every bead of that type.
        """
        if bead_type not in self.alphabet:
            return np.zeros(self.n_beads, dtype=bool)
        return self.bead_codes == self.alphabet.index(bead_type)
