# Amphiphile bilayer demo

The **membrane** analogue of the `slab_phase_separation` demo. Instead of a slab of a
single sticky homopolymer, this builds a lipid-style **bilayer** out of amphiphilic
molecules, in the same kind of z-elongated periodic box with "solvent" (empty space)
above and below.

## The molecule and the physics

Each molecule is a 5-bead amphiphile **`HHTTT`**: two hydrophilic **head** beads
(`H`) and three hydrophobic **tail** beads (`T`). The interactions in `params.prm`
encode the hydrophobic effect:

- **Tails cohere and hide from solvent** — `T–T` is strongly attractive (`-20`) and
  every tail–solvent contact is penalised (`+4`), so the tails bury themselves into a
  hydrophobic **core**.
- **Heads like solvent** — each head–solvent contact is favourable (`-4`), so the
  heads sit on the two outer **surfaces** of the membrane.
- **Heads and tails demix** — `H–T` is mildly repulsive (`+6`), keeping the two in
  separate layers.

Together these make a **tails-in / heads-out** bilayer a stable arrangement: a
hydrophobic tail core sandwiched between two head layers, with solvent on both sides.

## Seeded, then relaxed

A clean bilayer that spans the periodic plane is hard to *nucleate* from a random
start, so — exactly like the `multiphase_core_shell` demo — the membrane is **seeded**
(`build_bilayer_seed.py` lays down two leaflets of amphiphiles) and this run *relaxes*
it. The amphiphile interactions hold the bilayer together, so it settles into a
natural, rough-edged membrane rather than dissolving — a clean demonstration that
these interactions self-organise a bilayer.

## Files

| file | what it is |
|------|------------|
| `build_bilayer_seed.py` | generates `bilayer_seed.restart` (two leaflets, 512 amphiphiles) |
| `bilayer_seed.restart`  | the seeded bilayer, ready to load |
| `params.prm`            | the amphiphile interaction matrix (see above) |
| `KEYFILE.kf`            | the simulation |

## Run it

```bash
cd demo_keyfiles/amphiphile_bilayer
python build_bilayer_seed.py     # (optional - the restart file is already here)
PIMMS -k KEYFILE.kf
```

## Visualize

Load `traj.xtc` onto `START.pdb` and colour by bead name (`H` vs `T`):

- **VMD:**  `vmd START.pdb traj.xtc`
- **PyMOL:** `load START.pdb` then `load traj.xtc, START`

Viewed **edge-on** (look down x or y) you see the classic membrane cross-section:
a band of tails (`T`) sandwiched between two layers of heads (`H`), with empty
solvent above and below. The keyfile sets `TRAJECTORY_PBC_UNWRAP : True` so lipids
that cross a box face are drawn whole rather than split.

## Tuning

- **Membrane thickness / lipid length** — change the sequence and the seed z-columns
  in `build_bilayer_seed.py` (e.g. `HHTTTT` for a thicker core).
- **How "liquid" the membrane is** — raise `TEMPERATURE` for a more disordered, fluid
  membrane; lower it for a tighter, more gel-like one. Too hot and the hydrophobic
  core loses cohesion and the membrane breaks up.
- **Head/tail balance** — more tail beads (or stronger `T–T` / larger tail-solvent
  penalty) makes a more strongly-segregated membrane; weaker values let heads and
  tails intermix and blur the leaflets.
- **Area** — the seed places one amphiphile per (x, y) column per leaflet, so the
  membrane tiles the whole periodic plane; change `DIMS` in the build script to
  resize it.
