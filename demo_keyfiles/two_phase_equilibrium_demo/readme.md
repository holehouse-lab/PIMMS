# Minimal system
The keyfile here is a minimal keyfile that forms a two-phase coexisting system. One type of bead "A" is defined, and at this temperature, the chains should rapidly coalesce into a single liquid droplet.

This is a good sytem to shart playing with parameters in the keyfile to explore how the simulation output changes. Suggested things to play with include...

* Change the temperature (default = 70). In general wise to think about (interaction strength / T ) as the relevant interaction strength...

* Add in new chains or bead types (i.e. where do short single beads partition?)

* Change moveset probabilities and number of crankshaft moves

* change `XTC_FREQ` to see

* change box dimensions (super small boxe, super big box etc).