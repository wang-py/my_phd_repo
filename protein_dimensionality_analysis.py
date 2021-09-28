# Script used to analyze the dimensionality of protein structure
# Author: Panyue Wang
# Email: pywang@ucdavis.edu

from os import read
import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from calculate_interaction_energy import get_distance_vec, parse_pdb, assign_params, read_topology

# distance cutoff function
def apply_cutoff(r_cutoff, atom_index, coords_params):
    """
    apply distance cutoff for coulomb interactions
    ----------------------------------------------
    r_cutoff: float
    cutoff distance in angstroms

    atom_index: int
    index of focus atom

    coords_params: ndarray
    parameters of atoms

    Returns: 
    ----------------------------------------------
    trunc_params: ndarray
    truncated vector of distances within cutoff range
    """
    atom_dist = get_distance_vec(atom_index, coords_params[:,0:3])
    coords_params_dist = np.c_[coords_params, atom_dist]
    atoms_within_cutoff = np.array([x[-1] < r_cutoff for x in coords_params_dist])
    # include self
    #atoms_within_cutoff[atom_index] = True
    atoms_in_r = coords_params[atoms_within_cutoff]

    return atoms_in_r

# TODO: add a linear fit to find out the slope aka dimensionality
def plot_dimensionality(radius, number_of_atoms):
    plt.figure()
    plt.plot(np.log(radius), np.log(number_of_atoms))
    plt.title("Log(# of atoms) vs. log(r)")
    plt.xlabel("log(r)")
    plt.ylabel("log(# of atoms)")
    plt.show()

if __name__ == "__main__":
    input_pdb = sys.argv[1]
    selected_atom = int(sys.argv[2])
    ppdb = PandasPdb()
    # convert pdb into a data frame
    input_df = parse_pdb(ppdb, input_pdb)
    # read in topology
    topology = read_topology(open("topology.dat"))
    # building data structure with properties of the atoms
    coords_params = assign_params(input_df, topology)
    # array of radii
    radius_range = np.linspace(0,10,21)
    # applying cutoff based on distance
    atoms_within_cutoff_arr = []
    number_of_atoms_arr = []
    for radius in radius_range:
        atoms_within_cutoff = apply_cutoff(radius, selected_atom, coords_params)
        atoms_within_cutoff_arr.append(atoms_within_cutoff)
        number_of_atoms = len(atoms_within_cutoff)
        number_of_atoms_arr.append(number_of_atoms)
        #print("Number of atoms within %f angstroms is %d"%(radius, number_of_atoms))

    number_of_atoms_arr = np.array(number_of_atoms_arr)
    atoms_within_cutoff_arr = np.array(number_of_atoms_arr)

    plot_dimensionality(radius_range, number_of_atoms_arr)
