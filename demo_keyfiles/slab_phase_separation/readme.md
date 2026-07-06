# Slab phase-separation demo

This demo builds a **phase-separated slab** using the standard "slab method",
and doubles as a showcase of PIMMS' **resized equilibration** (`RESIZED_EQUILIBRATION`)
with **periodic boundaries** and a **non-cubic production box**.

## The idea

A high concentration of sticky 5-bead polymers (`SSSSS`, strong `S–S` attraction) is
packed into a small **cubic** box (`20 × 20 × 20`) and equilibrated, so the melt
condenses into one dense blob that fills the box. At the end of equilibration PIMMS
grows **only the Z axis** to the production `DIMENSIONS` (`20 × 20 × 60`), keeping X
and Y fixed, and re-centres the material. Because the dense phase already spans the
periodic X–Y cross-section, it settles into a flat **slab** in the middle of the long
box, with a dilute vapour above and below it along Z — two flat interfaces under
periodic boundaries. This slab geometry is the standard way to measure coexistence
densities and interfacial properties.

## Files

| file | what it is |
|------|------------|
| `params.prm` | one sticky bead type `S` with a strong `S–S` short-range attraction |
| `KEYFILE.kf` | the simulation (cubic equilibration box → Z-elongated slab box) |

## Run it

```bash
cd demo_keyfiles/slab_phase_separation
PIMMS -k KEYFILE.kf
```

The first `EQUILIBRATION` steps run in the dense `20×20×20` cube; the box then
resizes to `20×20×60` and the remaining steps let the slab equilibrate. With
`SAVE_EQ : True` the trajectory captures all three stages (dense cube → resize →
slab).

## Visualize

Load the trajectory `traj.xtc` onto the topology `START.pdb`:

- **VMD:** `vmd START.pdb traj.xtc`
- **PyMOL:** `load START.pdb` then `load traj.xtc, START`

Orient the view so Z is horizontal: you should see a dense block of beads in the
middle of the long box with empty (vapour) regions on either side along Z.

## Key ingredients (and one restriction)

- **`RESIZED_EQUILIBRATION : 20 20 20`** — the (smaller, cubic) box used during
  equilibration; the box grows to the full `DIMENSIONS` afterwards. It must be
  ``<=`` `DIMENSIONS` on every axis.
- **`HARDWALL : False`** — periodic boundaries are essential; a slab needs the dense
  phase to tile across the X–Y faces and present two interfaces along Z (hardwall
  walls would not give a periodic slab).
- **No cluster rotation.** The production box is non-cubic and periodic, so
  `MOVE_CLUSTER_ROTATE` is intentionally left off — a 90° rigid rotation is only an
  energy-preserving symmetry of a cube/square, and PIMMS rejects that combination at
  startup with an explanatory error. Crankshaft + slither rearrange the dense melt
  perfectly well here.

## Tuning

- **Concentration / slab thickness** — more chains (or a smaller equilibration box)
  give a thicker, denser slab; fewer give a thinner one. Keep the equilibration box
  cubic and dense so the condensate fills the X–Y cross-section.
- **Temperature** — lower `TEMPERATURE` sharpens the interfaces and raises the dense
  density (down to where the melt kinetically arrests); higher weakens the demixing.
  Around `T = 52` with the `S–S = -20` contact gives a clean slab.
- **Z elongation** — a longer production Z (e.g. `20 20 80`) leaves more vapour room
  around the slab; the slab thickness itself is set by the amount of material.
