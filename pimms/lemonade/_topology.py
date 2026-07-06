"""
Trajectory topology: the fixed (time-independent) description of which beads
belong to which chain, each chain's sequence and type, and the bead-type alphabet.

Stored columnar for speed: a CSR-style ``offsets`` array delimits each chain's
contiguous block of atoms, so any per-chain reduction is a single
``numpy.add.reduceat`` and any type filter is a boolean mask over one ``(n_atoms,)``
array - no per-chain Python objects.
"""

import numpy as np

from pimms import CONFIG

# PIMMS writes each bead as a PDB residue whose name is one_to_three(type_char):
# known amino acids use their standard 3-letter code, everything else becomes
# 'XX<c>'. Invert that so we can read the 1-letter bead type straight back.
_THREE_TO_ONE = {three: one for one, three in CONFIG.ONE_TO_THREE.items()}


def three_to_one(resname):
    """Decode a PDB residue name back to its 1-letter PIMMS bead type."""
    if resname in _THREE_TO_ONE:
        return _THREE_TO_ONE[resname]
    stripped = resname.lstrip("X")
    return stripped if stripped else resname


class Topology:
    """Columnar chain/bead topology for a trajectory.

    Attributes
    ----------
    n_chains, n_atoms : int
    offsets : (n_chains + 1,) int64
        Atom-index boundaries; chain ``c`` owns atoms ``offsets[c]:offsets[c+1]``.
    lengths : (n_chains,) int64
        Number of beads in each chain.
    chain_types : (n_chains,) int32
        Integer type label per chain.
    sequences : list[str]
        1-letter bead sequence per chain.
    atom_chainid : (n_atoms,) int32
        Owning chain index for every atom.
    bead_codes : (n_atoms,) int8
        Index into ``alphabet`` for every atom's bead type.
    alphabet : list[str]
        Sorted unique 1-letter bead types.
    """

    __slots__ = ("offsets", "lengths", "chain_types", "sequences",
                 "atom_chainid", "bead_codes", "alphabet")

    def __init__(self, sequences, chain_types=None):
        self.sequences = list(sequences)
        self.lengths = np.array([len(s) for s in self.sequences], dtype=np.int64)
        self.offsets = np.zeros(len(self.sequences) + 1, dtype=np.int64)
        np.cumsum(self.lengths, out=self.offsets[1:])

        if chain_types is None:
            # group identical sequences into the same type
            seen = {}
            chain_types = []
            for s in self.sequences:
                chain_types.append(seen.setdefault(s, len(seen)))
        self.chain_types = np.asarray(chain_types, dtype=np.int32)

        n_atoms = int(self.offsets[-1])
        self.atom_chainid = np.empty(n_atoms, dtype=np.int32)
        for c in range(len(self.sequences)):
            self.atom_chainid[self.offsets[c]:self.offsets[c + 1]] = c

        all_chars = "".join(self.sequences)
        self.alphabet = sorted(set(all_chars))
        code = {ch: i for i, ch in enumerate(self.alphabet)}
        self.bead_codes = np.fromiter((code[ch] for ch in all_chars),
                                      dtype=np.int8, count=n_atoms)

    # -- construction ------------------------------------------------------
    @classmethod
    def from_mdtraj(cls, md_topology):
        """Build from an mdtraj Topology (each bead is one residue/atom)."""
        sequences = []
        for chain in md_topology.chains:
            sequences.append("".join(three_to_one(atom.residue.name)
                                     for atom in chain.atoms))
        return cls(sequences)

    def with_keyfile_types(self, chain_specs):
        """Return a copy whose chain types follow the keyfile CHAIN order.

        ``chain_specs`` is the parser's ``[[count, sequence], ...]`` list. If it is
        consistent with this topology (same chain count, same per-chain sequences)
        the authoritative per-spec type index is used; otherwise the topology is
        returned unchanged (sequence-grouped types).
        """
        expanded = []                       # (sequence, type) per chain, in order
        for type_idx, (count, seq) in enumerate(chain_specs):
            for _ in range(int(count)):
                expanded.append((seq, type_idx))
        if len(expanded) != len(self.sequences):
            return self
        if any(seq != self.sequences[i] for i, (seq, _t) in enumerate(expanded)):
            return self
        return Topology(self.sequences, chain_types=[t for _s, t in expanded])

    # -- convenience -------------------------------------------------------
    @property
    def n_chains(self):
        return len(self.sequences)

    @property
    def n_atoms(self):
        return int(self.offsets[-1])

    def type_mask(self, bead_type):
        """Boolean ``(n_atoms,)`` mask selecting beads of a given 1-letter type."""
        if bead_type not in self.alphabet:
            return np.zeros(self.n_atoms, dtype=bool)
        return self.bead_codes == self.alphabet.index(bead_type)
