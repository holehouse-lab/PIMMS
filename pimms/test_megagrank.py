## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab 
## Copyright 2015 - 2023
## ...........................................................................

import random
import sys

from . import mega_crank

print(random.random())

C_RAND_FACTOR = 0.0000000000001
local_seed = random.randint(1,sys.maxsize-1)*C_RAND_FACTOR

print(local_seed)
mega_crank.seed_C_rand(local_seed)

for i in range(0,9):
    mega_crank.crand_test()

for (start, end) in zip([0,1], [20,21]):
    count_start = 0
    count_end   = 0 
    for i in range(0, 10):
        
        print(mega_crank.python_randint(start,end))
        if mega_crank.python_randint(start,end) == start:
            count_start=count_start+1
        elif mega_crank.python_randint(start,end) == end:
            count_end = count_end+1

    print("Range [%i to %i] - got %i = %i and %i = %i" %(start, end, count_start, start, count_end, end))
    


print("Testing limits")

print(mega_crank.randint_tester(1,20, 2147483647-10))
print(mega_crank.randint_tester(1,20, 2147483647-1))
