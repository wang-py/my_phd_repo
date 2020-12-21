import os
import sys
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    # folder that contains all result pdb frames
    pdb_folder = sys.argv[1]
    # mode
    mode_pdb = open(sys.argv[2], 'r')
    # reading the pdb folder
    results = os.fsencode(pdb_folder)
    lines = mode_pdb.readlines()
    # starting line of chain H
    chainH_start = 620
    # lines of one frame
    frame_lines = 814
    # from 38 to 76
    frame_index = int(sys.argv[3]) - 1
    mode_start = frame_index * frame_lines
    #end = start + frame_lines
    for i in range(frame_lines - 3):
        input_lines[chainH_start + i] = lines[mode_start + i + 1]

    for one_result in os.scandir(results):
        filename = os.fsdecode(one_result)
    pass