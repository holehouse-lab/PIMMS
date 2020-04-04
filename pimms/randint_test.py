## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2020
## ...........................................................................


from . import mega_crank

for i in range(0,10000000):
    print(mega_crank.randint_ext(0, 3))
