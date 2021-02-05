import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import sys
from biopandas.pdb import PandasPdb

# distance cutoff function
# TODO: write a function to sort the distances and filter out atoms that are
# further away than cutoff distance.

# TODO: write a function to read topology

# TODO: verify file reading

ppdb = PandasPdb()

def parse_pdb(input_pdb):
    ppdb.read_pdb(input_pdb)
    coords = ppdb.df['ATOM'][['x_coord', 'y_coord', 'z_coord']]
    coords = coords.to_numpy()
    return coords

def get_distance_vec(xi, x):
    """
    calculate the distance vector between xi and x
    ---------------------------------------- 
    xi:
    index of coordinates of atom i 

    x:
    Nx3 vector of all coordinates in angstroms
    Returns
    ---------------------------------------- 
    r: float
    Nx1 vector of all distances in angstroms

    """
    # calculate the xyz differences between all the atoms and the focus atom
    c_diff = x - x[xi] 
    # calculate the distance
    r = np.sqrt(np.sum(np.square(c_diff), axis=-1)) 
    # prevent division by zero
    r[xi] = 1.0

    return r

# lorentz berthelot combining rule
def LB_combining(si, sj, ei, ej):
    sij = 0.5 * (si + sj)
    eij = np.sqrt(ei * ej)

    return sij, eij

def get_lennard_jones_potential(R, epsilon, sigma):
    """
    calculate the LJ potential between i and j
    ---------------------------
    R:
    distance between two charges in angstrom

    epsilon:
    epsilon in LJ

    sigma:
    sigma in LJ

    Returns
    ---------------------------
    U: float
    LJ potential energy of the system in kJ/mol
    """

    U = 4 * epsilon * (np.power(sigma / R, 12) - np.power(sigma / R, 6))

    return U

def get_coulomb_potential(qij, R):
    """
    calculate the coulomb potential between qi and qj
    ---------------------------
    qij:
    product of charges on i and j in e^2

    R:
    distance between two charges in angstrom

    Returns
    ---------------------------
    Uc: float
    Coulomb potential energy of the system in kJ/mol
    """
    #Coulomb constant in kJ/mol*A/e^2
    k0 = 1389.35456 
    Uc = k0 * qij / R

    return Uc

def get_energy(x, topology):
    """
    calculate the potential energy in the system
    ---------------------------
    x:
    positions of the atoms in nm

    topology:
    properties of the atoms

    Returns
    ---------------------------
    U: float
    potential energy of the system in kJ/mol
    """
    # unpacking topology
    charge = topology[0]
    sigma = topology[1]
    epsilon = topology[2]
    atom_count = x.shape[0]
    U_sys = np.zeros(atom_count)
    U = 0.0

    for i in range(atom_count):
        # calculate the xyz differences between all the atoms and the focus atom
        c_diff = x - x[i] 
        # calculate the distance
        r_mat = np.sqrt(np.sum(np.square(c_diff), axis=-1)) 
        # prevent division by zero
        r_mat[i] = 1.0
        # figure out sigma and epsilon
        s_mat, e_mat = LB_combining(sigma[i], sigma, epsilon[i], epsilon)
        # calculate charge
        c_mat = charge[i] * charge
        # prevent the self interaction
        c_mat[i] = charge[i]
        U_sys = force_field_potential(r_mat, c_mat, e_mat, s_mat)
        # prevent the self interaction
        U_sys[i] = 0 
        U += np.sum(U_sys)

    # remove duplicate interactions
    U /= 2

    return U

if __name__ == "__main__":
    input_pdb = sys.argv[1]
    coords = parse_pdb(input_pdb)
    atom_index = int(sys.argv[2])
    D = get_distance_vec(500, coords)
