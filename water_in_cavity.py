import numpy as np
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt

# script to find out water molecules (oxygens) that are inside of a cavity
# TODO: cavity preprocessing: split cavity in half by z coordinates, exclude
# waters that are out side the two boxes.
# TODO: find enclosed water molecules