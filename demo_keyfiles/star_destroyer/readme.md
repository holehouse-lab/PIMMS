# Star Destroyer demo

A bit of fun that doubles as a showcase of **fixed (frozen) chains loaded from a
restart file**, **`EXTRA_CHAIN`**, and **`PARALLELIZE` working alongside frozen
chains**.

A 3D Imperial Star Destroyer is built out of ~170 PIMMS chains, written to a
restart file, and held completely **fixed** for the whole run (every hull chain is
listed in `freeze.txt`). A swarm of 150 two-bead "spaceships" is then added with
`EXTRA_CHAIN` and left to **fly around** the stationary ship (the forcefield is pure
excluded volume, so the ships just zip about and bounce off the hull and each
other).

## Files

| file | what it is |
|------|------------|
| `build_star_destroyer.py` | generates `star_destroyer.restart` + `freeze.txt` (the hull geometry) |
| `star_destroyer.restart`  | the prebuilt hull (chainIDs `1..170`), ready to load |
| `freeze.txt`              | freezes all 170 hull chains |
| `params.prm`              | bead types `H` (hull), `T` (tower), `S` (ship); all interactions zero |
| `KEYFILE.kf`              | the simulation |

## Run it

```bash
cd demo_keyfiles/star_destroyer
python build_star_destroyer.py     # (optional - the restart/freeze files are already here)
PIMMS -k KEYFILE.kf
```

## Visualize

Load the trajectory `traj.xtc` onto the topology `START.pdb` (written at startup):

- **VMD:**  `vmd START.pdb traj.xtc`
- **PyMOL:** `load START.pdb` then `load traj.xtc, START`

Colour by atom/bead name to separate the hull (`H`), the command tower (`T`) and
the spaceships (`S`). You should see the frozen wedge-shaped Star Destroyer sitting
still while the little 2-bead ships dart around it.

## How the frozen ship is made

Each horizontal row of the hull is one connected chain; rows are stacked in length
(`x`) and height (`z`) to fill a solid wedge: a flat triangular **deck** (the
top-down arrowhead), a **wedge underside** that thickens toward the stern, a raised
command **superstructure**, a **bridge** block, a central **tower** spike and two
sensor **domes**. `build_star_destroyer.py` has all the geometry and is easy to
tweak (make it bigger, add a hangar, etc.).

Because the hull chains are frozen they never move, but they still occupy the
lattice and exclude the spaceships - so the ships fly *around* the ship, never
through it. `PARALLELIZE : True` is enabled to show that parallelization now
composes with frozen chains (the frozen hull beads are kept as fixed obstacles
while the moves are threaded).
