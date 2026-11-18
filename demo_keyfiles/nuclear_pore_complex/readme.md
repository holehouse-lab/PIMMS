# Nuclear pore complex demo

A schematic human nuclear pore complex (NPC) sitting in a nuclear envelope, at 2 nm per lattice site. The envelope and the NPC scaffold are frozen. FG-nucleoporins hang from the scaffold and fill the central channel. Transport receptors that stick to FG motifs, and inert cargo of the same size that doesn't, start in the cytoplasm. Over the run almost every receptor gets through to the nucleus, while most of the cargo stays out.

It shows off:

- **Freezing whole structures and single beads together.** Six frozen chains (`C` lines) make the envelope and scaffold, and each of the 88 FG-Nups is held by its first bead alone (`B <chain> 0` lines). The rest of each FG domain is free.
- **Compartments.** The membranes meet the hard walls, so the nucleus, the cytoplasm and the space between the membranes are sealed off from each other except through the pore.
- **A selective barrier from nothing but contact energies.** Receptors and cargo differ only in whether they stick to FG motifs.

## What's in it

All of it is built by `build_npc.py` from simple solids sized to published dimensions. It's a schematic, not a reconstruction of any one structure.

| part | how it's drawn | beads | frozen? |
|------|----------------|-------|---------|
| nuclear envelope | outer and inner membranes, each 2 sites (4 nm) thick, joined by a curved pore membrane (half a torus) | 32,392 (one chain) | whole |
| inner ring | fills the pore from the channel wall to the pore membrane, 8 spokes; a hollow shell 2 sites deep | 15,576 (one chain) | whole |
| cytoplasmic ring | two concentric head-to-tail rings of Y-complexes, Nup93 bridges, legs onto the membrane, connectors to the inner ring, Nup358 and Nup214 domains | 3,876 (one chain) | whole |
| nuclear ring | the same rings mirrored, with ELYS instead of Nup358 and Nup214 | 3,696 (one chain) | whole |
| nuclear basket | 8 Tpr filaments converging on a distal ring | 2,044 (one chain) | whole |
| luminal ring | the Gp210 ring in the space between the membranes, reaching through to the pore membrane | 1,492 (one chain) | whole |
| FG-Nups | 32 Nup62 (40 beads) and 32 Nup98 (80) on the inner ring; 8 Nup214 (60) and 8 Nup358 (64) on the cytoplasmic ring; 8 Nup153 (100) on the nuclear ring | 5,632 (88 chains) | first bead only |
| transport receptors | 6-bead chains that bind FG motifs | 144 (24 chains) | no |
| inert cargo | 6-bead chains that bind nothing | 144 (24 chains) | no |

That's 64,996 beads in a 93 x 93 x 111 box (186 x 186 x 222 nm).

### Where the numbers come from

| quantity | value used | source |
|----------|-----------|--------|
| NPC scaffold diameter | 120 nm | Mosalaganti et al., *Science* **376**, eabm9506 (2022): "a 120 nm-wide protein scaffold of three stacked rings" |
| spokes per ring | 8 | same: "each ring comprises eight spokes" |
| central channel (inner ring) diameter | 54 nm | same: the dilated, in-cell inner ring (42 nm constricted) |
| membrane pore diameter | 92 nm | same, Fig. 1 (dilated; 82 nm constricted) |
| Y-complexes in the outer rings | two rings, head to tail | same: "a double head-to-tail arrangement of Y-complexes" |
| spacing of inner and outer membranes | 50 nm | Hetzer, *Cold Spring Harb. Perspect. Biol.* **2**, a000539 (2010): "an even distance of ~50 nm" |
| membrane thickness | 4 nm | Mitra et al., *PNAS* **101**, 4083-4088 (2004): rough ER bilayer 37.5 A (the nuclear envelope is continuous with the ER) |
| FG-Nup copies per pore | multiples of 8 | Ori et al., *Mol. Syst. Biol.* **9**, 648 (2013) |

Illustrative choices, not fitted to anything:

- the heights of the rings (the outer rings' Y-complex layers are 70 nm apart top to bottom, 82 nm counting the domains on them);
- the basket's length (about 60 nm; published lengths vary widely);
- the FG domains: about 5 residues per bead, an FG motif every third bead, the lengths above;
- the copy numbers per spoke;
- every interaction energy (below).

The inner ring here seals the pore completely, so there are no peripheral channels beside it.

### Interactions

Energies are per contact, in units of kT (`TEMPERATURE : 100`, so -20 in `params.prm` is -0.2 kT):

| pair | energy (kT) |
|------|-------------|
| FxFG-FxFG (`F`-`F`) | -0.20 |
| FxFG-GLFG (`F`-`Q`) | -0.25 |
| GLFG-GLFG (`Q`-`Q`) | -0.30 |
| receptor-FG motif (`K`-`F`, `K`-`Q`) | -0.40 |
| everything else, including cargo and all scaffold and membrane | 0 |

The receptor binding was tuned by hand. At -0.6 kT the receptors pile into the FG meshwork and stay there; at -0.4 kT they pass through.

## Run it

```bash
cd demo_keyfiles/nuclear_pore_complex
python build_npc.py          # ~15 s; optional, the files it writes are already here
PIMMS -k KEYFILE.kf          # 60,000 steps, ~10 minutes
python analyze_transport.py  # where the receptors and cargo ended up
```

`build_npc.py` checks the envelope seals before it writes anything. The cytoplasm must reach the nucleus through the open channel and not once the channel is plugged, and the space between the membranes must reach neither.

The run writes 201 frames to `traj.xtc` (about 70 MB). `analyze_transport.py` places each receptor and cargo chain by the height of its centre: in the channel within the inner ring's height, the cytoplasm above, the nucleus below. It prints the counts over the run and writes them to `transport.csv`. With the shipped seed:

```
chains in the cytoplasm / channel / nucleus:
  frame    receptors (cyt/chan/nuc)    cargo (cyt/chan/nuc)
      0        24      0      0              24      0      0
     40        18      6      0              24      0      0
     80        16      4      4              24      0      0
    120        11      5      8              23      0      1
    160        11      3     10              22      0      2
    200         8      8      8              20      1      3
receptor: 24 of 24 entered the channel, 23 reached the nucleus
cargo   :  9 of 24 entered the channel,  4 reached the nucleus
```

A second seed gave the same picture (23 receptors and 5 cargo chains reached the nucleus).

Two things worth knowing when you read this:

- **There's no net import.** Real transport is driven by the RanGTP gradient, which this model doesn't have. Without it the receptors simply equilibrate between cytoplasm and nucleus: they end up roughly evenly split, about 8 each in the cytoplasm, the channel and the nucleus. What the model shows is facilitated diffusion: binding the FG motifs lets receptors through much faster than inert cargo.
- **Some cargo leaks.** The inert cargo is small, a 6-bead chain, and small molecules do diffuse slowly through real pores. A handful get through in a run. Make the cargo chains longer (`SOLUBLE_LENGTH`) and fewer will.

## Visualize

Load the trajectory onto the topology written at start-up (coordinates are to scale, since `LATTICE_TO_ANGSTROMS : 20` makes one site 2 nm):

- **VMD:** `vmd START.pdb traj.xtc`
- **PyMOL:** `load START.pdb`, then `load_traj traj.xtc, START`

Draw the beads as spheres of radius about 10 A (half a site) and colour by residue name. To see inside, clip the view through the pore axis or select one half, e.g. in VMD `y < 920` (the axis sits at 920 A).

| part | bead | residue name |
|------|------|--------------|
| nuclear envelope | `M` | MET |
| inner ring (alternate spokes) | `I`, `V` | ILE, VAL |
| cytoplasmic ring | `C` | CYS |
| Nup358 domains | `R` | ARG |
| Nup214 domains | `W` | TRP |
| nuclear ring | `N` | ASN |
| ELYS | `E` | GLU |
| nuclear basket | `T` | THR |
| luminal ring (Gp210) | `L` | LEU |
| FG-Nup anchor bead (frozen) | `A` | ALA |
| FxFG motif | `F` | PHE |
| GLFG motif (Nup98) | `Q` | GLN |
| FG spacer | `S` | SER |
| transport receptor | `K` | LYS |
| inert cargo | `D` | ASP |

## How it's built

**Drawing.** Each component is a set of lattice sites, drawn from cylinders, tori, rods and spheres in a frame with the pore axis along z and the envelope's mid-plane at z = 0.

**One chain per component.** Each frozen component is threaded into a single chain that visits every one of its sites once, each bonded to the next. This uses the same Hamiltonian-path search as the `star_destroyer_fleet` demo. That's why thin parts are at least 2 sites wide: a 1-site-wide spike would force the chain to end at its tip.

**FG-Nups and soluble chains.** Each FG-Nup's first bead is placed on a free site touching its anchoring structure. The rest of the chain is grown as a random self-avoiding walk into free space, steering away from dead ends. Receptors and cargo are grown the same way in the top of the cytoplasm.

**Why nothing leaks through the membranes.** Only local moves are switched on (crankshaft, slithering, pivots), and every one keeps each bead bonded to its chain neighbours. A bead can only move to a site bonded to a neighbour that stays put. The membranes are 2 sites thick, so no bond can span one. The hollow inside of the inner ring is out of reach for the same reason. A full run bears this out: in none of its 201 saved frames did a mobile bead sit in a membrane or in the space between the membranes. Moves that relocate a whole chain at once (`MOVE_CHAIN_TRANSLATE`, `MOVE_CLUSTER_TRANSLATE`, `MOVE_JUMP_AND_RELAX`, `MOVE_VMMC`) could carry a chain straight through the envelope, so leave them off.

## Changing it

Everything is set by the constants at the top of `build_npc.py`. Rebuild after changing any of them, and keep `DIMENSIONS` in `KEYFILE.kf` in step with `BOX_DIMENSIONS`:

- `CHANNEL_DIAMETER_NM`, `MEMBRANE_PORE_DIAMETER_NM`: try the constricted state (42 and 82 nm).
- `FG_NUPS`: the FG-Nups, their lengths and copies per spoke.
- `INTERACTIONS`: FG cohesion and receptor binding.
- `N_RECEPTORS`, `N_CARGO`, `SOLUBLE_LENGTH`: the soluble chains.
- `BUILD_SEED`: the FG and soluble chains' starting shapes.

The same seed always gives the same files.
