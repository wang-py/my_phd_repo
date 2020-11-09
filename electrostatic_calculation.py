import numpy as np
import sys
from biopandas.pdb import PandasPdb

def get_coulomb_force(q1, q2, r1, r2):
    """
    q1, q2: signs of charges 1 and 2
    r1, r2: positions of charges 1 and 2
    f1: the force on q1
    """
    r12 = np.linalg.norm(r1 - r2)
    r12_hat = (r1 - r2) / r12
    q12 = q1 * q2
    f_norm = q12 / np.square(r12)
    f1 = f_norm * r12_hat
    f1 += -0

    return f1

def convert_to_numpy(df_entry):
    return df_entry.to_numpy()

if __name__ == "__main__":
    ppdb = PandasPdb()
    # read in the pdb file from command line
    ppdb.read_pdb(sys.argv[1])
    df = ppdb.df
    #print(df['ATOM'][df['ATOM']['chain_id'] == '6'])
    # extract coordinates
    coords=['x_coord', 'y_coord', 'z_coord']
    # Two focus atoms on chain 4
    print("TYR87:")
    TYR87 = convert_to_numpy(df['ATOM'][coords][df['ATOM']['atom_number'] == 11221])
    print(TYR87)
    print("ASP139:")
    ASP139 = convert_to_numpy(df['ATOM'][coords][df['ATOM']['atom_number'] == 11618])
    print(ASP139)
    # all GLU on chain 6
    GLU = convert_to_numpy(df['ATOM'][coords][(df['ATOM']['chain_id'] == '6') \
        & (df['ATOM']['residue_name'] == 'GLU') \
        & (df['ATOM']['atom_name'] == 'OE2')])
    print(GLU)
    # all ASP on chain 6
    ASP = convert_to_numpy(df['ATOM'][coords][(df['ATOM']['chain_id'] == '6') \
        & (df['ATOM']['residue_name'] == 'ASP') \
        & (df['ATOM']['atom_name'] == 'OD2')])
    # all ARG on chain 6
    ARG = convert_to_numpy(df['ATOM'][coords][(df['ATOM']['chain_id'] == '6') \
        & (df['ATOM']['residue_name'] == 'ARG') \
        & (df['ATOM']['atom_name'] == 'NH1')])
    # all LYS on chain 6
    LYS = convert_to_numpy(df['ATOM'][coords][(df['ATOM']['chain_id'] == '6') \
        & (df['ATOM']['residue_name'] == 'LYS') \
        & (df['ATOM']['atom_name'] == 'NZ')])
    # N2

    f1 = get_coulomb_force(-1, 1, np.array([1,0,0]), np.array([0,0,0]))
    print("F1 is:", f1)
