import numpy as np
import sys
from biopandas.pdb import PandasPdb

if __name__ == "__main__":
    ppdb = PandasPdb()
    # read in the pdb file from command line
    ppdb.read_pdb(sys.argv[1])
    df = ppdb.df
    print(df['ATOM'][df['ATOM']['chain_id'] == '6'])