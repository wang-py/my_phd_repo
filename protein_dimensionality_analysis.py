# Script used to analyze the dimensionality of protein structure
# Author: Panyue Wang
# Email: pywang@ucdavis.edu

import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from calculate_interaction_energy import apply_cutoff, parse_pdb

# TODO: plotting function that plots log # of atoms vs log(radius)

if __name__ == "__main__":
    input_pdb = sys.argv[1]
    selected_atom = int(sys.argv[2])
    ppdb = PandasPdb()
    input_df = parse_pdb(ppdb, input_pdb)