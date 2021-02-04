import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import sys

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

if __name__ == "__main__":
    input_pdb = open(sys.argv[1], 'r')
    atom_index = int(sys.argv[2])
    exit()
