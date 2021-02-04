import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import sys

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
