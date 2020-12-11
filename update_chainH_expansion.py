import numpy as np
import sys

if __name__ == "__main__":
    # mode file
    mode = open(sys.argv[1], 'r')
    lines = mode.readlines()
    # input pdb (last step's output)
    input_pdb = open(sys.argv[2], 'r')
    input_lines = input_pdb.readlines()
    # starting line of chain H
    chainH_start = 620
    # lines of one frame
    frame_lines = 814
    # from 37 to 75
    frame_index = int(sys.argv[3])
    # output pdb (next step's input)
    output_pdb = open(sys.argv[4], 'w')
    mode_start = frame_index * frame_lines
    #end = start + frame_lines
    for i in range(frame_lines - 3):
        input_lines[chainH_start + i] = lines[mode_start + i + 1]
    
    # swapping the output chain H with the mode chain H
    for i in range(len(input_lines)):
        output_pdb.write(input_lines[i])

    pass