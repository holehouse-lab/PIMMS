



cimport numpy as cnp
import numpy as np

def randint_np(int start, int end):
    return np.random.randint(start, end+1)


def seed_randint_np(int seedval):
    np.random.seed(seedval)
