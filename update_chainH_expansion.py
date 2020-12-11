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
    #output_pdb = open(sys.argv[4], 'w')
    mode_start = frame_index * frame_lines
    #end = start + frame_lines
    for i in range(frame_lines - 3):
        input_line = input_lines[chainH_start + i]
        print(input_line, end="")
        #mode_line = lines[mode_start + i + 1]
        #print(mode_line, end="")
    pass