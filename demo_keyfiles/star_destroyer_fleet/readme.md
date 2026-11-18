# Star Destroyer fleet demo

A detailed Imperial Star Destroyer, a squadron of 12 TIE fighters and a flight of 8 X-wings, where **every ship is a single PIMMS chain**. It shows off three things:

- **Freezing a whole structure with one line.** The Star Destroyer is one 43,725-bead chain, so `freeze.txt` is just `C 1`.
- **Rigid bodies.** The fighters aren't frozen, but the only moves switched on are rigid-body moves (VMMC translations and whole-chain rotations), which never change a chain's shape. So each fighter flies around as a rigid body, bouncing off the others and off the Star Destroyer.
- **Arbitrary shapes as chains.** `build_fleet.py` draws each ship as a set of lattice sites and threads them into one chain that visits every site exactly once.

It's the big sibling of the `star_destroyer` demo, which builds a simpler hull from about 170 short chains and flies 2-bead "ships" around it.

## Files

| file | what it is |
|------|------------|
| `build_fleet.py` | draws the three ships, threads each into one chain, places the fleet and writes the three files below |
| `fleet.restart` | the fleet: chain 1 is the Star Destroyer, chains 2-13 are TIE fighters, chains 14-21 are X-wings |
| `freeze.txt` | `C 1`: freezes the Star Destroyer |
| `params.prm` | 21 bead types, every interaction zero (pure excluded volume); the bead types are only there for colouring |
| `KEYFILE.kf` | the simulation |

## Run it

```bash
cd demo_keyfiles/star_destroyer_fleet
python build_fleet.py      # ~30 s; optional, the files it writes are already here
PIMMS -k KEYFILE.kf        # ~5 minutes
```

The keyfile's box (410 x 420 x 496 sites) is bigger than the one `build_fleet.py` builds the fleet in (210 x 120 x 96), so PIMMS puts the fleet in the middle of it, which leaves the fighters room to fly. The run needs about 1 GB of memory and writes 101 frames to `traj.xtc` (about 28 MB). In that time each fighter moves 7 to 43 sites and turns at least once. It manages only about 7 steps a second, because every VMMC step re-evaluates the energy of all 62,449 beads. PIMMS prints one start-up warning, that the move set is rigid only, so every chain keeps the shape it was built with. Here that's the point.

## Visualize

Load the trajectory onto the topology written at start-up:

- **VMD:** `vmd START.pdb traj.xtc`
- **PyMOL:** `load START.pdb`, then `load_traj traj.xtc, START`

Draw the beads as spheres (VDW) of radius about 1.8 Å, half of the 3.65 Å lattice spacing PIMMS writes into the PDB, and colour by residue name. Each bead type is written as a residue name:

| ship | part | bead | residue name |
|------|------|------|--------------|
| Star Destroyer | hull | `H` | HIS |
| | darker panels and greebles | `G` | GLY |
| | superstructure, tower and bridge | `T` | THR |
| | lit windows | `W` | TRP |
| | shield-generator domes | `D` | ASP |
| | engine housings | `K` | LYS |
| | engine glow | `E` | GLU |
| | hangar-bay lights | `L` | LEU |
| | turbolaser turrets | `A` | ALA |
| TIE fighter | cockpit | `C` | CYS |
| | viewport | `V` | VAL |
| | wing frame and spars | `F` | PHE |
| | solar panels | `P` | PRO |
| | wing pylons | `N` | ASN |
| X-wing | fuselage and wings | `S` | SER |
| | red stripes | `R` | ARG |
| | canopy | `Q` | GLN |
| | astromech | `X` | XXX |
| | laser cannons | `I` | ILE |
| | engines | `M` | MET |
| | engine glow | `Y` | TYR |

In VMD, for example, `resname HIS GLY THR` picks out the Star Destroyer's hull and superstructure.

## How the ships are made

**Drawing.** Each ship is built from simple solids on the lattice:

- **Star Destroyer (128 sites long).** The arrowhead hull is two low pyramids with an equatorial trench between them along the sides. On top sit a stepped four-tier superstructure with rows of lit windows, the bridge tower and its two shield domes, turbolaser turrets along the flanks and the hull edges, and a few hundred raised greebles and darker panels scattered over the deck. Three main engines and four smaller ones sit on the stern, with recessed glowing nozzles, and there's a lit hangar bay in the keel.
- **TIE fighter.** A ball cockpit with a forward viewport, two pylons, and two hexagonal wings with a rim, a hub and six spars.
- **X-wing.** A tapering nose, a cockpit with a canopy and an astromech behind it, four S-foils open in the X with red stripes, an engine at each foil root, and a laser cannon along each wingtip.

The fighters are drawn far larger than true scale: next to a 128-site Star Destroyer, a real TIE fighter would be smaller than one site.

**Threading.** To be one chain, a ship's sites have to be visited in an order where each site is bonded to the next (in PIMMS two beads are bonded when they differ by at most 1 on every axis). That's a Hamiltonian path through the ship. `build_fleet.py` finds one greedily: it always steps to the free neighbour with the fewest free neighbours of its own (Warnsdorff's rule), avoiding steps that would strand a site. When it gets stuck it reverses the end of the chain to free up a new end (a Posa rotation). Any site still left over is spliced in between two chain sites bonded to it. The script checks that every ship's chain covers it exactly, with every bond valid, before writing anything.

**Shapes that thread.** Thin parts have to be at least two sites wide so the chain can run out along them and back. That's why the X-wing's cannons are 2 sites wide: a 1-site-wide spike would force the chain to end at its tip, and a chain only has two ends.

**A hollow Star Destroyer.** Solid, the Star Destroyer would be about 66,000 beads. Only a skin 3 sites deep is kept, which brings it to 43,725. Every site inside is at least 4 sites from open space, so a fighter can't get in by a VMMC step of up to 2 sites, and a fighter can't straddle the skin. But a move that relocates a chain anywhere in the box, `MOVE_CHAIN_TRANSLATE` or `MOVE_CLUSTER_TRANSLATE`, could drop a fighter inside. So leave those off, and keep `VMMC_MAX_DISPLACEMENT` at 2 or less.

## Changing it

Everything is set by the constants at the top of `build_fleet.py`. Rebuild after changing any of them, and keep `DIMENSIONS` in `KEYFILE.kf` at least as big as `BOX_DIMENSIONS` on every axis:

- `N_TIE_FIGHTERS`, `N_X_WINGS`: the number of fighters.
- `TIE_REGION`, `XWING_REGION`: where they start.
- `SD_LENGTH`, `SD_TIERS`, `SD_MAIN_ENGINES` and friends: the shape of the Star Destroyer.
- `BUILD_SEED`: the greebles, the threading and the fighter placement.

The same seed always gives the same files. A bigger Star Destroyer makes every VMMC step slower, and past 99,999 beads in total the PDB's atom numbers wrap.
