import numpy as np
import sys
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt

# script to find out water molecules (oxygens) that are inside of a cavity
# TODO: cavity preprocessing: split cavity in half by z coordinates, exclude
# waters that are out side the two boxes.
def separate_input_structure(input_ppdb):
    """
    separate cavity coords and water from input pdb
    """
    cavity_coords = input_ppdb.df['ATOM'][atoms.df['ATOM']['atom_name' == 'C']
    water_coords = input_ppdb.df['ATOM'][atoms.df['ATOM']['atom_name' == 'O']

    cavity_coords = cavity_coords.df['ATOM']\
        [['atom_number', 'atom_name','x_coord', 'y_coord', 'z_coord']].to_numpy()
    
    water_coords = water_coords.df['ATOM']\
        [['atom_number', 'atom_name','x_coord', 'y_coord', 'z_coord']].to_numpy()

    return cavity_coords, water_coords

def preprocess_cavity(cavity_coords):
    pass
# TODO: find enclosed water molecules

if __name__ == "__main__":
    # input pdb containing water and cavity
    input_pdb = sys.argv[1]
    ppdb = PandasPdb()
    ppdb.read_pdb(input_pdb)

    cavity_coords, water_coords = separate_input_structure(ppdb)