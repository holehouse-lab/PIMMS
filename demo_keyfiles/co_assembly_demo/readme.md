## Demo 1

This example keyfile represents a simulation for two different chain typs which co-assemble. To run, after installing PIMMS, simply run as:


	PIMMS -k KEYFILE.kf
	
from within this directory.

The simulation should taken ~5 minutes, and will generate a trajectory with 100 frames and 100 datapoints. It will run for 5K macroscopic Monte Carlo (MC) steps (`N_STEPS : 5000`) but but in reality will perform ~40 million independent accept/reject moves.

The keyfile defines two chains 

	20 AAAA
	50 B
	
The `[AAAA]` chains and the `[B]` chains will co-assemble with one another, as defined by the interaction strengths in `params.prm` file.
