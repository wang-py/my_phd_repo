import numpy as np
import sys
from pymol import cmd

###############################################################################
# this script strips chain H from the input structure                         #
# Author: Panyue Wang                                                         #
# Email: pywang@ucdavis.edu                                                   #
###############################################################################
if __name__ == "__main__":
    # read in the pdb file from command line
    cmd.load(sys.argv[1])
    # output path
    output_path = sys.argv[1].split('.')[0] + "_chainH.pdb"
    # select chain H
    chainH = "chainH"
    cmd.select(chainH, "chain H and (bb. + sc.)")
    cmd.save(output_path, selection=chainH)