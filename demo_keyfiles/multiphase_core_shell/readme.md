# Multiphase core-shell condensate demo

A **four-layer "onion" condensate**: four sticky 5-bead species (`A`, `B`, `C`, `D`)
form a single droplet made of four **concentric layers** — `A` in the core, wrapped
by `B`, then `C`, then `D` on the outside. It is a minimal model of the multiphase,
multilayered condensates seen in biology (nucleolus-style core-shell organisation).

## How the layering is engineered

Concentric layering needs two ingredients, both in `params.prm`:

- **Adjacent layers attract.** Each species attracts itself (`self = -25`) and its
  immediate neighbours in the stack (`A–B`, `B–C`, `C–D` = `-20`), so all four
  condense into one droplet and consecutive layers wet each other.
- **Non-adjacent layers repel.** `A–C`, `A–D` and `B–D` are set to `+15`
  (repulsive), so a layer refuses to touch any layer other than the two it belongs
  between. That repulsion is what forces the species into a nested `A|B|C|D` onion
  rather than a scrambled blob.

With these (symmetric) interactions the concentric `A|B|C|D` order is a stable
free-energy minimum; **which** species is the core is fixed by the starting
configuration (see below).

## Two-step run

The starting droplet already has the four species *roughly* banded — the simulation
then rounds the droplet and **sharpens** the four interfaces.

```bash
cd demo_keyfiles/multiphase_core_shell
python build_onion_seed.py      # writes onion_seed.restart (optional - it's already here)
PIMMS -k KEYFILE.kf
```

`build_onion_seed.py` packs 640 chains into a sphere and assigns them to four
concentric bands by radius (fewer chains in the inner bands, more in the outer ones,
so the layers come out at similar thickness). The keyfile loads that restart and runs
crankshaft + slither Monte Carlo at `T = 28`.

## Files

| file | what it is |
|------|------------|
| `build_onion_seed.py` | generates `onion_seed.restart` (the pre-banded droplet) |
| `onion_seed.restart`  | 640 chains in 4 concentric bands (`A` core → `D` shell) |
| `params.prm`          | the interaction matrix (adjacent-attract / non-adjacent-repel) |
| `KEYFILE.kf`          | the simulation |

## Visualize

Load `traj.xtc` onto `START.pdb` and colour by bead name (`A`/`B`/`C`/`D`):

- **VMD:**  `vmd START.pdb traj.xtc`
- **PyMOL:** `load START.pdb` then `load traj.xtc, START`

A slice through the centre of the droplet shows four concentric rings; a radial
composition profile (bead count vs. distance from the droplet centre) shows four
cleanly separated peaks, core (`A`) to shell (`D`).

## Tuning

- **Which species is the core** is set purely by the seed order in
  `build_onion_seed.py` (the interactions are symmetric); swap the band assignment to
  reorder the onion.
- **Layer thickness** is set by how many chains go to each band (`COUNTS` in the
  build script). An inner shell has less volume, so it needs fewer chains.
- **Sharpness / roundness** — colder (`TEMPERATURE` down) sharpens the interfaces but
  slows rounding; a longer `N_STEPS` lets the droplet relax to a cleaner sphere. The
  `+15` non-adjacent repulsion is the main knob for how strictly the layers refuse to
  intermix — lower it and the interfaces broaden, raise it and they stay crisp.
- **Number of layers** — the same adjacent-attract / non-adjacent-repel pattern
  extends to any number of species; add an `E` bead (attract `D`, repel `A`/`B`/`C`)
  for a fifth shell.
