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
    for i in range(chainH_start, chainH_start+frame_lines):
        input_line = input_lines[i]
        print(input_line, end="")
    # from 37 to 75
    frame_index = int(sys.argv[3])
    # output pdb (next step's input)
    #output_pdb = open(sys.argv[4], 'w')
    start = frame_index * frame_lines
    end = start + frame_lines
    for i in range(start, end):
        line = lines[i]
        print(line, end="")
    pass