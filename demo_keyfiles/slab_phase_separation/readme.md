# Slab phase-separation demo

This demo builds a **phase-separated slab** using the standard "slab method", and doubles as a showcase of PIMMS' **resized equilibration** (`RESIZED_EQUILIBRATION`) with **periodic boundaries** and a **non-cubic production box**.

## The idea

600 sticky 5-bead polymers (`SSSSS`, 3000 beads, with an `S–S` contact energy of -20) are packed into a small **cubic** `20 × 20 × 20` box (37.5% of lattice sites occupied) and equilibrated at `T = 120`. PIMMS always runs the equilibration box under hardwall boundaries, whatever `HARDWALL` says, so during equilibration the melt condenses into one dense droplet sitting in the middle of the cube, clear of the walls. At the end of equilibration PIMMS grows **only the Z axis** to the production `DIMENSIONS` (`20 × 20 × 60`), keeping X and Y fixed, re-centres the material in the new box, and switches to the periodic boundaries the keyfile asks for (`HARDWALL : False`). The droplet is almost as wide as the 20-site cross-section, so once the X–Y faces are periodic it touches its own images and within the first hundred or so production steps spreads into a flat **slab** that spans the X–Y cross-section, sits in the middle of the long box, and has a nearly empty vapour above and below it along Z — two flat interfaces under periodic boundaries. This slab geometry is the standard way to measure coexistence densities and interfacial properties.

In our runs of this keyfile the slab is about 8 lattice sites thick, its interior is about 90% occupied, and fewer than 100 of the 3000 beads are in the vapour at any time; it stays in place for the rest of the 8500 steps.

## Files

| file | what it is |
|------|------------|
| `params.prm` | one sticky bead type `S` with a `-20` `S–S` short-range contact, no angle penalty and no solvent interaction |
| `KEYFILE.kf` | the simulation (cubic equilibration box → Z-elongated slab box) |

## Run it

```bash
cd demo_keyfiles/slab_phase_separation
PIMMS -k KEYFILE.kf
```

The first `EQUILIBRATION` (1000) steps run in the dense `20×20×20` cube; the box then resizes to `20×20×60` and the remaining 7500 steps let the slab equilibrate. Because the box changes size, the run writes two trajectories, each with its own topology:

| files | what they hold |
|-------|----------------|
| `eq_START.pdb` + `eq_traj.xtc` | the equilibration stage in the `20×20×20` box (41 frames: the starting configuration plus one every 25 steps). Written only because `SAVE_EQ : True`. |
| `START.pdb` + `traj.xtc` | the production stage in the `20×20×60` box (301 frames). Frame 0 is the re-centred configuration just after the resize, at step 1000. |

`TRAJECTORY_PBC_UNWRAP : True` makes every chain whole in each written frame, so no chain is drawn split across a periodic face (this only affects the written coordinates, not the simulation).

## Visualize

Load each trajectory onto its own topology:

- **VMD:** `vmd eq_START.pdb eq_traj.xtc` (equilibration) and `vmd START.pdb traj.xtc` (production)
- **PyMOL:** `load START.pdb` then `load traj.xtc, START` (and likewise `eq_START.pdb` / `eq_traj.xtc`)

In the production trajectory, orient the view so Z is horizontal: you should see a dense block of beads in the middle of the long box with empty (vapour) regions on either side along Z.

## Key ingredients (and one restriction)

- **`RESIZED_EQUILIBRATION : 20 20 20`** — the (smaller, cubic) box used during equilibration; the box grows to the full `DIMENSIONS` afterwards. It must be `<=` `DIMENSIONS` on every axis. The equilibration box is always hardwall.
- **`HARDWALL : False`** — periodic boundaries in the production box are essential; a slab needs the dense phase to tile across the X–Y faces and present two interfaces along Z (hardwall walls would not give a periodic slab).
- **No cluster rotation.** The production box is non-cubic and periodic, so `MOVE_CLUSTER_ROTATE` is intentionally left off — a 90° rigid rotation is only an energy-preserving symmetry of a cube/square, and PIMMS rejects that combination at startup with an explanatory error. Crankshaft + slither rearrange the dense melt perfectly well here.

## Tuning

- **Amount of material / slab thickness** — more chains give a thicker slab and fewer a thinner one; the densities of the slab and the vapour are set by the temperature, not by how much material there is. The slab only forms if the condensate is wide enough to reach across the periodic X–Y cross-section; with much less material it can stay a droplet instead.
- **Temperature** — the keyfile runs at `T = 120`, where the `S–S = -20` contact already gives a sharp slab with an almost empty vapour. Lowering `TEMPERATURE` makes the slab denser and the vapour emptier but slows rearrangement inside the dense phase (eventually it kinetically arrests); raising it weakens the demixing.
- **Z elongation** — a longer production Z (e.g. `20 20 80`) leaves more vapour room around the slab; the slab thickness itself is set by the amount of material.
