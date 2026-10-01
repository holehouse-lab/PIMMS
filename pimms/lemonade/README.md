# lemonade

A fast, hierarchical analysis backend for PIMMS lattice trajectories.

`lemonade` loads a PIMMS simulation from an XTC / PDB / keyfile and exposes it as a navigable object hierarchy — `LatticeTrajectory → Frame → Polymer` (and `Frame → Cluster → Polymer`) — while keeping all the heavy data in contiguous numpy arrays so loading and analysis stay fast on large systems.

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
frame   = traj[-1]                    # a Frame (int index; frame 0 is the start configuration); traj[::2] gives a sub-trajectory
polymer = frame[3]                    # a Polymer view: chain 3 in that frame
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
result.shape["sphericity"], result.shape["radius_of_gyration"]   # droplet geometry only - shape is None for slab (a percolating slab has no meaningful hull)
result.is_phase_separated                    # heuristic yes/no
```

Individual tools:

- `condensed_fraction(traj)`, `number_of_clusters(traj)`, `largest_cluster_size(traj)`, `cluster_size_distribution(traj)` — per-frame order parameters.
- `spanning_fraction(traj)` — fraction of frames whose largest cluster spans the box (connected to its own periodic image, or touching both walls under a hardwall); `all_axes=True` counts only box-filling networks.
- `droplet_shape(traj)` — frame-averaged Rg, asphericity, sphericity, hull volume and density of the largest cluster, leaving out frames where it spans the box.
- `radial_density_profile(traj)` / `radial_density_profile_with_site_counts(traj)` — spherically averaged occupied-fraction profile about the condensate COM (dense core → interface → dilute background), optionally with the lattice sites per shell that the fit uses to drop thin shells.
- `slab_density_profile(traj, axis=...)` — 1D profile along the long axis, slabs re-centred per frame (the `slab_phase_separation` geometry); under a hardwall only a slab clear of both walls is moved.
- `fit_radial_profile(...)` / `fit_slab_profile(...)` — `tanh` fits returning a `BinodalFit` (`rho_dense`, `rho_dilute`, `interface_width`, `radius`/`half_width`, plus `success`/`reason` — always check `success`: a flat one-phase profile otherwise produces a converged but meaningless fit).
- `frame.droplet` — the largest cluster in a frame, with `.radius_of_gyration`, `.volume`, `.sphericity`, `.density`, `.radial_density_profile()`.

Densities are occupied lattice-site fractions in `[0, 1]`, comparable across box sizes.

### Surface tension from interfacial undulations

`pimms.lemonade.surface_tension` estimates the condensate surface tension from capillary-wave / shape fluctuations (k_B T is taken from the trajectory temperature; γ is in reduced units — interaction energy per lattice area):

```python
from pimms.lemonade import surface_tension as st

st.surface_tension(traj)                 # auto: slab vs droplet by box shape
st.slab_surface_tension(traj)            # <|h(q)|^2> = kT / (gamma A q^2)  (robust; on the lattice q^2 is the exact dispersion 2-2cos(q) per axis)
st.droplet_surface_tension(traj)         # <|u_lm|^2> = kT / (gamma R0^2 (l-1)(l+2))
```

Each returns a `SurfaceTension` (`gamma`, `gamma_std`, `n_modes`, `spectrum`). Both estimators are 3D only, and need a temperature (from the keyfile, or `temperature=`). The **slab** method (planar capillary waves off the two flat interfaces of a box-spanning condensate) is the robust one, and skips any frame whose largest cluster is not a slab (a network that also spans the normal, or a droplet or strip). The **droplet** method (spherical-harmonic shape fluctuations) needs a single, compact, reasonably large droplet sampled over many frames, skips frames whose largest cluster spans the box, and reads low: 10-17 % low on lattice droplets of known γ, and about a quarter below the slab estimate for the same chains on real PIMMS droplets, with a strong dependence on the angular grid (see the phase-separation docs). Treat it as a rough, low estimate, and always check `gamma_std` and the returned `spectrum`. For the slab method, `n_modes` counts independent Fourier wavevectors: the conjugate `+q` and `-q` coefficients of a real height field count once.

## Loading

`load()` accepts any sensible combination of inputs:

| inputs | result |
|--------|--------|
| `xtc` + `pdb` | full trajectory (topology from the PDB, exact XTC bead order) |
| `xtc` + `pdb` + `keyfile` | as above, plus authoritative spacing / dimensions / hardwall / chain types |
| `pdb` only | a single frame (e.g. `START.pdb`) |

Without a keyfile, the lattice spacing defaults to PIMMS's `3.65 Å`, the box is inferred from the trajectory's unit cell, `hardwall` is `False` and `temperature` is `None`. The run's own `keyfile_used.kf` is the best keyfile to pass, since it holds the resolved configuration (restart overrides applied, merged chain list); given the original keyfile instead, `load()` resolves `RESTART_OVERRIDE_*` from the restart file. The `eq_` files of a `RESIZED_EQUILIBRATION` run are always loaded in the compact box under hard walls (keyfile or not), and a quench keyfile gives `QUENCH_END` as the temperature. Overrides (`spacing=`, `dimensions=`, `hardwall=`, `temperature=`) and frame selection (`start`/`stop`/`step`, `n_frames`) are available.

## Why it's fast

The original lemonade was slow because it built a Python object per chain **per frame** and eagerly painted a full grid for every frame. lemonade instead:

- stores the whole trajectory as one contiguous `(n_frames, n_beads, 3)` int32 array with CSR chain offsets; `Frame` / `Polymer` / `Cluster` are thin **views**, so navigating allocates nothing per bead;
- converts XTC coordinates back to the integer lattice in a **single vectorised** `round(nm / (spacing/10))`;
- makes every chain whole across periodic boundaries in a **compiled kernel** (`kernels/_pbc.pyx`, batched over all frames and chains — under 1 ms for 101 frames × 250 chains, over a hundred times faster than calling PIMMS's per-chain `make_chain_whole` in a Python loop);
- computes Rg / COM / gyration tensor / end-to-end for the **whole trajectory at once** with `numpy.add.reduceat` (no per-chain or per-frame Python loop);
- builds grids and clusters **lazily**, per frame, only when asked.

As a reference point, a 101-frame × 250-chain (2000-bead) trajectory loads in ~40-50 ms and its full per-chain Rg array computes in ~4 ms; `phase_separation.analyze()` on it takes about a second.

## Correctness

- The integer-lattice recovery is exact (float32 round-off only).
- The PBC unwrap is **bit-identical** to PIMMS's `make_chain_whole`.
- Rg matches its mathematical definition exactly and agrees with PIMMS's own `get_polymeric_properties`, which since 1.0.8 is computed on the same whole (bond-walked) chain, for chains of any size that do not percolate the box.
- Cluster detection, single-image reconstruction and gross properties (volume, surface area, density, radial profile) reuse PIMMS's own (Cython-accelerated) routines.

## Build

lemonade ships a compiled kernel (`kernels/_pbc.pyx`), so the package must be built before use:

```bash
./build.sh uv       # from the pimms repo root; clean rebuild + editable install with uv
./build.sh pip      # ... or the same with pip
```

`build.sh` takes exactly one argument naming the installer (`uv` or `pip`) and prints its usage and exits without building if given anything else.

## Layout

```
lemonade/
  __init__.py       public API: load, LatticeTrajectory, Frame, Polymer, Cluster,
                    phase_separation, surface_tension
  _load.py          load() - orchestration, coordinate conversion, inference
  _topology.py      Topology (CSR chain offsets, sequences, types, bead codes)
  _store.py         TrajectoryStore - columnar backing store + memoised batched results
  _analysis.py      vectorised numeric core (Rg / COM / gyration / distance maps)
  trajectory.py     LatticeTrajectory
  frame.py          Frame
  polymer.py        Polymer
  cluster.py        Cluster
  phase_separation.py  analyze(), density profiles, binodal fits, order parameters
  surface_tension.py   capillary-wave / shape-fluctuation surface tension
  kernels/_pbc.pyx  compiled PBC unwrap + grid painting
  tests/            end-to-end tests against real PIMMS output
```
