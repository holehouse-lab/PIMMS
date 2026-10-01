This example keyfile runs a single-chain simulation: one 80-bead homopolymer (`B`) in a 70 x 70 x 70 box at T = 400, for 1020 steps of which the first 20 are equilibration.

Run it from this directory with `PIMMS -k KEYFILE.kf`. It takes a few seconds (about 2.6 s wall time on an M3 Max MacBook Pro) and writes a trajectory (`traj.xtc`, topology `START.pdb`) with 205 frames - the starting configuration plus one frame every 5 steps - and 200 data points in each per-step polymer analysis file (`RG.dat`, `END_TO_END_DIST.dat`, `ASPH.dat`), one every 5 production steps from step 25 to step 1020.
