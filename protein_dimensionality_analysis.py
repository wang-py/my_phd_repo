# Script used to analyze the dimensionality of protein structure
# Author: Panyue Wang
# Email: pywang@ucdavis.edu

import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from calculate_interaction_energy import get_distance_vec, parse_pdb, assign_params, read_topology

# distance cutoff function
def apply_cutoff(r_cutoff, atom_index, atom_df):
    """
    apply distance cutoff for coulomb interactions
    ----------------------------------------------
    r_cutoff: float
    cutoff distance in angstroms

    atom_index: int
    index of focus atom

    atom_df: ndarray
    parameters of atoms

    Returns: 
    ----------------------------------------------
    trunc_params: ndarray
    truncated vector of distances within cutoff range
    """
    atom_dist = get_distance_vec(atom_index, atom_df[:,0:3])
    atom_df_dist = np.c_[atom_df, atom_dist]
    atoms_within_cutoff = np.array([x[-1] < r_cutoff for x in atom_df_dist])
    # include self
    #atoms_within_cutoff[atom_index] = True
    atoms_in_r = atom_df[atoms_within_cutoff]

    return atoms_in_r

def plot_dimensionality(radius, number_of_atoms):
    log_r = np.log(radius)
    log_n = np.log(number_of_atoms)
    slope, y_intercept = np.polyfit(log_r, log_n, 1)
    y_fit = slope * log_r + y_intercept

    plt.figure()
    plt.plot(log_r, log_n)
    plt.plot(log_r, y_fit, 'k--', label="linear fit k = %.3f"%slope)
    plt.legend()
    plt.title("Log(# of atoms) vs. log(r)")
    plt.xlabel("log(r)")
    plt.ylabel("log(# of atoms)")
    plt.show()

# TODO: include HETATMs in the pdb to get a more accurate count of atoms
if __name__ == "__main__":
    input_pdb = sys.argv[1]
    # ID in pymol
    selected_atom = int(sys.argv[2])
    ppdb = PandasPdb()
    # convert pdb into a data frame
    input_df = parse_pdb(ppdb, input_pdb)
    # read in topology
    topology = read_topology(open("topology.dat"))
    atom_i = input_df[input_df['atom_number'] == selected_atom].index[0]
    # building data structure with properties of the atoms
    atom_df = input_df[['x_coord', 'y_coord', 'z_coord',\
        'residue_number', 'atom_number']].to_numpy()
    #coords_params = assign_params(input_df, topology)
    # array of radii
    radius_range = np.linspace(1.2,30,41)
    # applying cutoff based on distance
    atoms_within_cutoff_arr = []
    number_of_atoms_arr = []
    for radius in radius_range:
        atoms_within_cutoff = apply_cutoff(radius, atom_i, atom_df)
        atoms_within_cutoff_arr.append(atoms_within_cutoff)
        number_of_atoms = len(atoms_within_cutoff)
        number_of_atoms_arr.append(number_of_atoms)
        print("Number of atoms within %f angstroms is %d"%(radius, number_of_atoms))

    number_of_atoms_arr = np.array(number_of_atoms_arr)
    atoms_within_cutoff_arr = np.array(number_of_atoms_arr)

    plot_dimensionality(radius_range, number_of_atoms_arr)
