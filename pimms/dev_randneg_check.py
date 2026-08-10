## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2026
## ...........................................................................


from . import mega_crank


# Developer sanity-check script, run by hand with `python -m pimms.<name>`.
# The body is guarded so that merely importing this module (it ships inside the
# package) does not execute a long side-effecting loop at import time.
if __name__ == "__main__":
    for i in range(0,1000000):

        print(mega_crank.get_random_position_python(1))
    
