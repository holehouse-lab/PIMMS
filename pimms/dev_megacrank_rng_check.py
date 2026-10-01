## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2026
## ...........................................................................

##
## Developer sanity-check script for the mega_crank random-number helpers.
##
## Seeds the kernel PRNG, prints a few raw draws (crand_test) next to PRNG_MAX,
## then for two inclusive ranges prints a batch of randint_python draws and
## reports how many of THOSE draws landed on each end of the range, flagging any
## draw outside it. Each draw is taken once and then printed, counted and
## range-checked, so the summary line describes exactly the numbers shown.
##
## This is NOT a test module - it is a throwaway script with top-level side
## effects, run by hand with `python -m pimms.dev_megacrank_rng_check`. It used
## to be named test_megagrank.py, which meant pytest collected it and the whole
## session died during collection on its top-level code (it also called a
## `python_randint` entry point that mega_crank has never exported - the real
## names are `randint_python` / `randint_ext`).
##

import random

from . import mega_crank


# Guarded so importing the module never runs the checks; run with
# `python -m pimms.dev_megacrank_rng_check`.
if __name__ == "__main__":

    # seed_C_rand takes a C int, so the seed has to be an integer that fits in
    # one (it used to be a float, which was silently truncated)
    local_seed = random.randint(1, 2**31 - 1)
    print("seed: %i" % local_seed)
    mega_crank.seed_C_rand(local_seed)

    prng_max = mega_crank.RAND_MAX_test()
    print("raw draws (crand_test, expected in [0, %i]):" % prng_max)
    for i in range(0, 9):
        raw = mega_crank.crand_test()
        flag = "" if 0 <= raw <= prng_max else "   <-- OUT OF RANGE"
        print("  %i%s" % (raw, flag))

    n_draws = 10
    for (start, end) in zip([0, 1], [20, 21]):
        count_start = 0
        count_end = 0
        count_outside = 0
        print("randint_python(%i, %i) x %i:" % (start, end, n_draws))
        for i in range(0, n_draws):

            # one draw per iteration, used for everything below
            value = mega_crank.randint_python(start, end)
            print("  %i" % value)
            if value == start:
                count_start = count_start + 1
            elif value == end:
                count_end = count_end + 1
            if value < start or value > end:
                count_outside = count_outside + 1

        print("Range [%i to %i] - got %i = %i and %i = %i, %i outside the range (expected %.2f at each end)" % (start, end, count_start, start, count_end, end, count_outside, n_draws / (end - start + 1)))
