import numpy as np
import sys
from biopandas.pdb import PandasPdb

###############################################################################
# this script strips chain H from the input structure                         #
# Author: Panyue Wang                                                         #
# Email: pywang@ucdavis.edu                                                   #
###############################################################################
if __name__ == "__main__":
    ppdb = PandasPdb()
    # read in the pdb file from command line
    ppdb.read_pdb(sys.argv[1])
    # output path
    output_path = sys.argv[1].split('.')[0] + "_chainH.pdb"
    df = ppdb.df
    df = df['ATOM'][df['ATOM']['chain_id'] == 'H']
    ppdb.to_pdb(output_path)