import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import sys
from biopandas.pdb import PandasPdb

def read_topology(top_file):
    """
    read in topology from file
    ----------------------------------------------
    top_file: file obj
    file that contains information about topology

    Returns:
    ----------------------------------------------
    topology: ndarray
    array that contains topology info about atoms
    """

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
    """
    assign parameters to atoms
    ----------------------------------------------
    atom_df: Pandas Dataframe
    pdb data frame of all atoms

    topology: ndarray
    array of atom topology

    Returns:
    ----------------------------------------------
    coords_params:
    array that contains all atom information
    """
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
    """
    read pdb file using biopandas
    ----------------------------------------------
    ppdb: ppdb obj
    pandas pdb object

    input_pdb: file obj
    input pdb file

    Returns:
    ----------------------------------------------
    atom_df: pandas dataframe
    dataframe that contains all information needed for energy calculation
    """
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
    #total_potential = LJ_potential + coulomb_potential

    return LJ_potential, coulomb_potential

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
    within_cutoff = np.array([x[-1] < r_cutoff for x in coords_params_dist])
    trunc_params = coords_params[within_cutoff]
    # prevent division by zero
    atom_dist[atom_index] = 1.0
    trunc_atom_dist = atom_dist[within_cutoff]

    return trunc_params, trunc_atom_dist

def turn_off_interaction(resi, coords_params):
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
        res_atom_i = first_atom_i + j
        res_sigma[j] = sigma[res_atom_i]
        res_epsilon[j] = epsilon[res_atom_i]
        res_charge[j] = charge[res_atom_i]
        # no interaction within molecule
        coords_params[res_atom_i, 6] = 0
        coords_params[res_atom_i, 7] = 0
    
    return coords_params, res_x, res_sigma, res_epsilon, res_charge

# TODO: distance cutoff should be different for every atom
# atoms of the same molecule should have zero charge wrt each other
def get_energy(coords_params, resi, r_cutoff=None):
    """
    calculate the potential energy in the system
    ---------------------------
    coords_params:
    positions of the atoms in A with LJ and coulomb parameters

    resi:
    index of residue to calculate interaction energy
    
    r_cutoff: float
    cutoff distance in angstroms

    Returns
    ---------------------------
    U_sys_LJ: float
    Lennard-Jones energy of the focus molecule in kJ/mol

    U_sys_C: float
    Coulomb energy of the focus molecule in kJ/mol
    """
    # turn off inteaction within molecule
    coords_params, res_x, res_sigma, res_epsilon, res_charge = \
        turn_off_interaction(resi, coords_params)

    res_atom_count = res_x.shape[0]
    # index of the first atom in the molecule
    first_atom_i = int(res_x[0,4] - 1)

    U_sys_LJ = 0.0
    U_sys_C = 0.0
    total_charge = 0.0
    for i in range(res_atom_count):
        # indexing through the residue
        atom_i = first_atom_i + i
        if r_cutoff:
            coords_params_cutoff, r_mat = apply_cutoff(r_cutoff, atom_i, coords_params.copy())
        else:
            coords_params_cutoff = coords_params
            x = coords_params_cutoff[:, 0:3]
            r_mat = get_distance_vec(atom_i, x)
        # unpacking topology
        sigma = coords_params_cutoff[:,5]
        epsilon = coords_params_cutoff[:,6]
        charge = coords_params_cutoff[:,7]
        # figure out total charge within cutoff range
        total_charge += np.sum(charge)
        # figure out sigma and epsilon
        s_mat, e_mat = LB_combining(res_sigma[i], sigma, res_epsilon[i], epsilon)
        # calculate charge
        c_mat = res_charge[i] * charge
        # prevent the self interaction
        #c_mat[atom_i] = charge[atom_i]
        U_LJ, U_C = force_field_potential(r_mat, c_mat, e_mat, s_mat)
        # prevent the self interaction
        # NOTE: might be unnecessary since epsilon and charge are already zero
        #U_sys[atom_i] = 0
        # total interaction energy
        U_sys_LJ += np.sum(U_LJ)
        U_sys_C += np.sum(U_C)

    return U_sys_LJ, U_sys_C, total_charge

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
    H2O_i = residue_index
    r_cutoff = np.linspace(0.5, 50.5, 51)
    E_H2O = np.zeros([r_cutoff.shape[0]])
    E_H2O_pair_wise = get_energy(coords_params.copy(), H2O_i)[1]
    charge_H2O = np.zeros([r_cutoff.shape[0]])
    result = np.zeros([3])
    for i in range(r_cutoff.shape[0]):
        result = get_energy(coords_params.copy(), H2O_i, r_cutoff[i])
        E_H2O[i] = result[1]
        charge_H2O[i] = result[2]
    plt.figure()
    plt.plot(r_cutoff, E_H2O, 'o')
    plt.hlines(E_H2O_pair_wise, r_cutoff[0], r_cutoff[-1], linestyle='--',label="pairwise coulomb:%.2f kJ/mol"%E_H2O_pair_wise)
    plt.xlabel("cutoff distance [A]")
    plt.ylabel("nonbonded interaction energy [kJ/mol]")
    plt.legend()
    plt.figure()
    plt.plot(r_cutoff, charge_H2O, 'o')
    plt.xlabel("cutoff distance [A]")
    plt.ylabel("total charge")
    plt.legend()
    plt.show()
    #print("E_O is %f kJ/mol, E_H1 is %f kJ/mol and E_H2 is %f kJ/mol"%(E_O, E_H1, E_H2))
    print("E_H2O is %f kJ/mol"%E_H2O_pair_wise)
