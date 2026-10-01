## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


##
## Developer sanity-check script for the cluster-translate offset draw.
##
## Draws per-axis offsets exactly the way moves.py does for a cluster translation
## with a CLUSTER_MOVE_THRESHOLD - randneg(randint(1, min(d-1, threshold))) for
## each of the three axes of an 80^3 box - tallies every draw, and prints the
## count for each offset next to the count a uniform distribution over
## +/-[1, min(d-1, threshold)] would give, plus a chi-squared statistic. Every
## offset should come up about equally often, never 0, and nothing beyond the
## threshold; chi-squared should be close to its degrees of freedom.
##


import random

from . import numpy_utils


# Developer sanity-check script, run by hand with `python -m pimms.<name>`.
# The body is guarded so that merely importing this module (it ships inside the
# package) does not execute a long side-effecting loop at import time.
if __name__ == "__main__":
    dims = [80, 80, 80]
    cluster_move_threshold = 60
    n_iterations = 5000000

    results_count = {}

    for i in range(0, n_iterations):

        # one offset per axis, as in a cluster translation; every draw is counted
        for d in dims:
            r = numpy_utils.randneg(random.randint(1, min(d-1, cluster_move_threshold)))

            if r not in results_count:
                results_count[r] = 0
            results_count[r] = results_count[r] + 1

    # all three axes are the same length here, so they share one expected support
    max_offset = min(dims[0] - 1, cluster_move_threshold)
    support = [v for v in range(-max_offset, max_offset + 1) if v != 0]
    n_draws = n_iterations * len(dims)
    expected = n_draws / len(support)

    print("offset   count   (expected %.1f per offset)" % expected)
    chi2 = 0.0
    for v in support:
        observed = results_count.get(v, 0)
        chi2 = chi2 + (observed - expected)**2 / expected
        print("%6i %8i   %+.2f%%" % (v, observed, 100 * (observed - expected) / expected))

    unexpected = sorted(v for v in results_count if v not in support)
    print("draws: %i, distinct offsets seen: %i (expected %i)" % (n_draws, len(results_count), len(support)))
    print("offsets outside +/-[1, %i]: %s" % (max_offset, unexpected if unexpected else "none"))
    print("chi-squared: %.1f on %i degrees of freedom" % (chi2, len(support) - 1))

