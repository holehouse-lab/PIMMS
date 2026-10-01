# Minimal phase-separating system

The keyfile here is a small, fast keyfile for a phase-separating mixture. It places 800 short chains in a 30 x 30 x 30 hardwall box - 200 each of `AAAAA`, `BBBBB` and `CCCCC`, and 200 of `DDDD` - and runs 200 steps (the first 50 are equilibration) at T = 60 using only slither (reptation) moves on the parallel kernel (`PARALLELIZE : True`). It finishes in a few seconds. In `params.prm` every like pair (A-A, B-B, C-C, D-D) has a -20 contact energy and every unlike pair a weaker -5, with no angle penalties and no solvent interaction.

What you should see: within the first few frames the chains condense out of the uniform starting configuration, and because like contacts are four times stronger than unlike ones the four chain types also demix from one another, so by the end of the run about 90% of inter-chain contacts are between chains of the same type (a random mix would give about 25%). Each droplet is therefore (almost) a single type, and because unlike pairs still attract, droplets of different types stick to one another in a handful of larger clusters. 200 steps is not long enough for coarsening to finish: at the end of the run there are still roughly ten separate droplets of each type. `HARDWALL` is on so that no chain is drawn split across a periodic face, which makes the trajectory easy to look at.

This is a good system to start playing with parameters in the keyfile to explore how the simulation output changes. Suggested things to play with include...

* Increase `N_STEPS` and watch the droplets keep merging.

* Change the temperature (the keyfile uses `TEMPERATURE : 60`). PIMMS works in reduced units with k = 1, so the Boltzmann factor is exp(-dE/T), and (interaction energy / T) is the relevant interaction strength.

* Change the like and unlike energies in `params.prm` - what happens as the unlike energies approach the like ones?

* Add in new chains or bead types (e.g. where do short single beads partition?).

* Change the moveset and the number of substeps per move. The keyfile uses slither alone (`MOVE_SLITHER : 1`, `SLITHER_SUBSTEPS : 200`); try mixing in `MOVE_CRANKSHAFT`, `MOVE_PULL` or `MOVE_VMMC` (the move fractions must add up to 1).

* Change `XTC_FREQ` to write more or fewer trajectory frames.

* Change the box dimensions (super small box, super big box etc).
