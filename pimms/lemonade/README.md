# lemonade

A fast, hierarchical analysis backend for PIMMS lattice trajectories.

`lemonade` loads a PIMMS simulation from an XTC / PDB / keyfile and exposes it as a
navigable object hierarchy — `LatticeTrajectory → Frame → Polymer` (and
`Frame → Cluster → Polymer`) — while keeping all the heavy data in contiguous numpy
arrays so loading and analysis stay fast on large systems.

## Quickstart

```python
import pimms.lemonade as lemonade

traj = lemonade.load(xtc="traj.xtc", pdb="START.pdb", keyfile="KEYFILE.kf")

# whole-trajectory analyses (vectorised; return (n_frames, n_chains) arrays)
rg   = traj.radius_of_gyration()
com  = traj.center_of_mass()          # (n_frames, n_chains, n_dim)
asph = traj.asphericity()
ete  = traj.end_to_end_distance()

# navigate the hierarchy
frame   = traj[0]                     # a Frame (int index); traj[::2] gives a sub-trajectory
polymer = frame[3]                    # a Polymer view: chain 3 in frame 0
polymer.radius_of_gyration            # scalar
polymer.whole_positions               # chain made whole across PBC (L, 3)
polymer.distance_map()                # (L, L) inter-bead distance matrix

for cluster in frame.clusters:        # connected-component clusters, largest first
    cluster.n_chains, cluster.n_beads
    cluster.radius_of_gyration, cluster.volume, cluster.density
    cluster.bead_type_composition     # {'A': ..., 'B': ...}
```

## Phase separation & droplet physics

`pimms.lemonade.phase_separation` quantifies liquid-liquid phase separation:

```python
from pimms.lemonade import phase_separation as ps

result = ps.analyze(traj)                    # auto-detects droplet vs slab geometry
result.rho_dense, result.rho_dilute          # binodal (coexistence occupied fractions)
result.condensed_fraction                    # fraction of material in the condensate
result.binodal.interface_width               # interfacial width (lattice units)
result.shape["sphericity"], result.shape["radius_of_gyration"]
result.is_phase_separated                    # heuristic yes/no
```

Individual tools:

- `condensed_fraction(traj)`, `number_of_clusters(traj)`, `largest_cluster_size(traj)`,
  `cluster_size_distribution(traj)` — per-frame order parameters.
- `radial_density_profile(traj)` — spherically averaged occupied-fraction profile about
  the condensate COM (dense core → interface → dilute background).
- `slab_density_profile(traj, axis=...)` — 1D profile along the long axis, slabs
  re-centred per frame (the `slab_phase_separation` geometry).
- `fit_radial_profile(...)` / `fit_slab_profile(...)` — `tanh` fits returning a
  `BinodalFit` (`rho_dense`, `rho_dilute`, `interface_width`, `radius`/`half_width`).
- `frame.droplet` — the largest cluster in a frame, with `.radius_of_gyration`,
  `.volume`, `.sphericity`, `.density`, `.radial_density_profile()`.

Densities are occupied lattice-site fractions in `[0, 1]`, comparable across box sizes.

### Surface tension from interfacial undulations

`pimms.lemonade.surface_tension` estimates the condensate surface tension from
capillary-wave / shape fluctuations (k_B T is taken from the trajectory temperature;
γ is in reduced units — interaction energy per lattice area):

```python
from pimms.lemonade import surface_tension as st

st.surface_tension(traj)                 # auto: slab vs droplet by box shape
st.slab_surface_tension(traj)            # <|h(q)|^2> = kT / (gamma A q^2)  (robust)
st.droplet_surface_tension(traj)         # <|u_lm|^2> = kT / (gamma R0^2 (l-1)(l+2))
```

Each returns a `SurfaceTension` (`gamma`, `gamma_std`, `n_modes`, `spectrum`). The
**slab** method (planar capillary waves off the two flat interfaces of a
box-spanning condensate) is the robust one; the **droplet** method (spherical-harmonic
shape fluctuations) needs a single, compact, reasonably large droplet sampled over
many frames, and is noisier for small lattice droplets — always check `gamma_std` and
the returned `spectrum`.

## Loading

`load()` accepts any sensible combination of inputs:

| inputs | result |
|--------|--------|
| `xtc` + `pdb` | full trajectory (topology from the PDB, exact XTC atom order) |
| `xtc` + `pdb` + `keyfile` | as above, plus authoritative spacing / dimensions / hardwall / chain types |
| `pdb` only | a single frame (e.g. `START.pdb`) |

Without a keyfile, the lattice spacing defaults to PIMMS's `3.65 Å` and the box is
inferred from the trajectory's unit cell. Overrides (`spacing=`, `dimensions=`,
`hardwall=`) and frame selection (`start`/`stop`/`step`, `n_frames`) are available.

## Why it's fast

The original lemonade was slow because it built a Python object per chain **per
frame** and eagerly painted a full grid for every frame. lemonade instead:

- stores the whole trajectory as one contiguous `(n_frames, n_atoms, 3)` int32
  array with CSR chain offsets; `Frame` / `Polymer` / `Cluster` are thin **views**,
  so navigating allocates nothing per bead;
- converts XTC coordinates back to the integer lattice in a **single vectorised**
  `round(nm / (spacing/10))`;
- makes every chain whole across periodic boundaries in a **compiled kernel**
  (`kernels/_pbc.pyx`, batched over all frames and chains — ~125× faster than the
  per-chain Python unwrap it replaces);
- computes Rg / COM / gyration tensor / end-to-end for the **whole trajectory at
  once** with `numpy.add.reduceat` (no per-chain or per-frame Python loop);
- builds grids and clusters **lazily**, per frame, only when asked.

As a reference point, a 101-frame × 250-chain (2000-bead) trajectory loads in
~40 ms and its full per-chain Rg array computes in ~3 ms.

## Correctness

- The integer-lattice recovery is exact (float32 round-off only).
- The PBC unwrap is **bit-identical** to PIMMS's `make_chain_whole`.
- Rg matches its mathematical definition exactly and agrees with PIMMS's own
  `get_polymeric_properties` in the dilute regime. For chains larger than ~half the
  box the two intentionally differ: lemonade uses the whole (contiguous) chain,
  whereas PIMMS returns the minimum-image value that collapses under finite-size
  artefacts.
- Cluster detection, single-image reconstruction and gross properties (volume,
  surface area, density, radial profile) reuse PIMMS's own (Cython-accelerated)
  routines.

## Build

lemonade ships a compiled kernel (`kernels/_pbc.pyx`), so the package must be built
before use:

```bash
./build.sh          # from the pimms repo root
```

## Layout

```
lemonade/
  __init__.py       public API: load, LatticeTrajectory, Frame, Polymer, Cluster
  _load.py          load() - orchestration, coordinate conversion, inference
  _topology.py      Topology (CSR chain offsets, sequences, types, bead codes)
  _store.py         TrajectoryStore - columnar backing store + memoised batched results
  _analysis.py      vectorised numeric core (Rg / COM / gyration / distance maps)
  trajectory.py     LatticeTrajectory
  frame.py          Frame
  polymer.py        Polymer
  cluster.py        Cluster
  kernels/_pbc.pyx  compiled PBC unwrap + grid painting
  tests/            end-to-end tests against real PIMMS output
```
