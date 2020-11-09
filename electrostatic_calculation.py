import numpy as np
import sys
from biopandas.pdb import PandasPdb

if __name__ == "__main__":
    ppdb = PandasPdb()
    # read in the pdb file from command line
    ppdb.read_pdb(sys.argv[1])
    df = ppdb.df
    #print(df['ATOM'][df['ATOM']['chain_id'] == '6'])
    # Two focus atoms on chain 4
    print("TYR87:")
    print(df['ATOM'][df['ATOM']['atom_number'] == 11221])
    print("ASP139:")
    print(df['ATOM'][df['ATOM']['atom_number'] == 11618])
    # all GLU on chain 6
    print(df['ATOM'][(df['ATOM']['chain_id'] == '6') & \
         (df['ATOM']['residue_name'] == 'GLU') & \
         (df['ATOM']['atom_name'] == 'OE2')])
    # all ASP on chain 6
    print(df['ATOM'][(df['ATOM']['chain_id'] == '6') & \
         (df['ATOM']['residue_name'] == 'ASP') & \
         (df['ATOM']['atom_name'] == 'OD2')])
    # all ARG on chain 6
    print(df['ATOM'][(df['ATOM']['chain_id'] == '6') & \
         (df['ATOM']['residue_name'] == 'ARG') & \
         (df['ATOM']['atom_name'] == 'NH1')])
    # all LYS on chain 6
    print(df['ATOM'][(df['ATOM']['chain_id'] == '6') & \
         (df['ATOM']['residue_name'] == 'LYS') & \
         (df['ATOM']['atom_name'] == 'NZ')])
    # N2