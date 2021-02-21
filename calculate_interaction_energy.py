import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import sys
from biopandas.pdb import PandasPdb

# distance cutoff function
# TODO: write a function to sort the distances and filter out atoms that are
# further away than cutoff distance.

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
    atom_dist = get_distance_vec(atom_index, coords_params[0:3])
    coords_params_dist = np.append(coords_params, atom_dist)
    trunc_params = np.array([x[3] < r_cutoff for x in coords_params_dist])
    return trunc_params

def read_topology(top_file):
    topology = {}
    entry = []
    for line in top_file:
        line_arr = line.split()
        if line_arr[0] == '#':
            continue
        else:
            atomtype = line_arr[0]
            sigma = float(line_arr[1])
            epsilon = float(line_arr[2])
            charge = float(line_arr[3])
            entry = [sigma, epsilon, charge]
            topology[atomtype] = entry
    
    return topology

# TODO: verify file reading

def assign_params(atom_df, topology):
    atom_coords = atom_df[['x_coord', 'y_coord', 'z_coord',\
        'residue_number', 'atom_number']].to_numpy()
    atom_count = atom_coords.shape[0]
    atom_xyz = atom_coords.shape[1]
    # create empty array to store atom xyz and interaction parameters
    coords_params = np.zeros([atom_count, atom_xyz + 3])
    for i in range(atom_count):
        atom_name = atom_df.iloc[i]['atom_name']
        params = np.array(topology[atom_name])
        coords_params[i] = np.append(atom_coords[i], params)
    
    return coords_params

def parse_pdb(ppdb, input_pdb):
    ppdb.read_pdb(input_pdb)
    atom_df = ppdb.df['ATOM']\
        [['atom_number', 'atom_name','x_coord', 'y_coord', 'z_coord', 'residue_number']]
    return atom_df

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

def force_field_potential(R, qij, epsilon, sigma):
    LJ_potential = get_lennard_jones_potential(R, epsilon, sigma)
    coulomb_potential = get_coulomb_potential(qij, R)
    total_potential = LJ_potential + coulomb_potential

    return total_potential

def get_energy(coords_params, resi):
    """
    calculate the potential energy in the system
    ---------------------------
    coords_params:
    positions of the atoms in A with LJ and coulomb parameters

    resi:
    index of residue to calculate interaction energy

    Returns
    ---------------------------
    U: float
    potential energy of the focus molecule in kJ/mol
    """
    # unpacking topology
    sigma = coords_params[:,5]
    epsilon = coords_params[:,6]
    charge = coords_params[:,7]
    x = coords_params[:, 0:3]
    # all atoms in residue i
    res_x = np.array([x for x in coords_params if x[3] == resi])
    # index of the first atom in the molecule
    first_atom_i = int(res_x[0,4] - 1)
    # number of atoms in this molecule
    res_atom_count = res_x.shape[0]
    # separate topology arrays for the molecule
    res_sigma = np.zeros(res_atom_count)
    res_epsilon = np.zeros(res_atom_count)
    res_charge = np.zeros(res_atom_count)
    for j in range(res_atom_count):
        res_sigma[j] = sigma[first_atom_i + j]
        res_epsilon[j] = epsilon[first_atom_i + j]
        res_charge[j] = charge[first_atom_i + j]
        #sigma[first_atom_i + j] = 0
        epsilon[first_atom_i + j] = 0
        charge[first_atom_i + j] = 0

    U = 0
    for i in range(res_atom_count):
        # indexing through the residue
        atom_i = first_atom_i + i
        # calculate the xyz differences between all the atoms and the focus atom
        c_diff = x - res_x[i, 0:3] 
        # calculate the distance
        r_mat = np.sqrt(np.sum(np.square(c_diff), axis=-1)) 
        # prevent division by zero
        r_mat[atom_i] = 1.0
        # figure out sigma and epsilon
        s_mat, e_mat = LB_combining(res_sigma[i], sigma, res_epsilon[i], epsilon)
        # calculate charge
        c_mat = res_charge[i] * charge
        # prevent the self interaction
        #c_mat[atom_i] = charge[atom_i]
        U_sys = force_field_potential(r_mat, c_mat, e_mat, s_mat)
        # prevent the self interaction
        # NOTE: might be unnecessary since epsilon and charge are already zero
        #U_sys[atom_i] = 0
        # total interaction energy
        U += np.sum(U_sys)

    return U

if __name__ == "__main__":
    ppdb = PandasPdb()
    input_pdb = sys.argv[1]
    coords = parse_pdb(ppdb, input_pdb)
    topology_file = open(sys.argv[2], 'r')
    topol = read_topology(topology_file)
    coords_params = assign_params(coords, topol)
    residue_index = int(sys.argv[3])
    #D = get_distance_vec(500, coords)
    # 7290 is at the center of the box
    H2O_i = 2431
    E_H2O = get_energy(coords_params, H2O_i)
    #print("E_O is %f kJ/mol, E_H1 is %f kJ/mol and E_H2 is %f kJ/mol"%(E_O, E_H1, E_H2))
    print("E_H2O is %f kJ/mol"%E_H2O)
