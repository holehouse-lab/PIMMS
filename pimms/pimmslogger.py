## ...........................................................................
## 
## PIMMS (Polymer Interactions in Multicomponent Mixtures)
## Alex Holehouse, Pappu Lab, Holehouse Lab
## Copyright 2015 - 2021
## ...........................................................................

from time import gmtime, strftime

from . import CONFIG
from  . import IO_utils


def initialize():
    IO_utils.wipe_file(CONFIG.OUTNAME_LOGFILE)
    with open(CONFIG.OUTNAME_LOGFILE,'a') as fh:
        fh.write("PIMMS Simulation\n")
        fh.write("Simulation Start:  %s \n" % (strftime("%Y-%m-%d %H:%M:%S")))

def log_warning(msg, timestamp=True):
    with open(CONFIG.OUTNAME_LOGFILE,'a') as fh:

        if timestamp:
            fh.write("> WARNING: [ %s ]: %s \n" % (strftime("%Y-%m-%d %H:%M:%S", gmtime()), msg))            
        else:
            spacer = " "*len("%s" % (strftime("%Y-%m-%d %H:%M:%S", gmtime())))
            fh.write(">            %s    %s \n" % (spacer, msg))
        
def log_status(msg, timestamp=True):
    with open(CONFIG.OUTNAME_LOGFILE,'a') as fh:

        if timestamp:
            fh.write("> STATUS: [ %s ]: %s \n" % (strftime("%Y-%m-%d %H:%M:%S", gmtime()), msg))
        else:
            spacer = " "*len("%s" % (strftime("%Y-%m-%d %H:%M:%S", gmtime())))
            fh.write("> STATUS:   %s    %s \n" % (spacer, msg))
        
        
