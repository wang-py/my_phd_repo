import sys
from pymol import cmd

###############################################################################
# this script aligns different chain Hs with one chain H from the input       #
# Author: Panyue Wang                                                         #
# Email: pywang@ucdavis.edu                                                   #
###############################################################################
if __name__ == "__main__":
    target = "target"
    to_align = "to_align"
    # read in the pdb file from command line
    cmd.load(sys.argv[1], target)
    # structure to align
    cmd.load(sys.argv[2], to_align)
    # output
    output_path = sys.argv[2].split('.')[0] + "_aligned.pdb"
    # select chain H
    cmd.align(to_align, target)
    cmd.save(output_path, selection=to_align)