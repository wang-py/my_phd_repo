import numpy as np
import sys

if __name__ == "__main__":
    # mode file
    mode = open(sys.argv[1], 'r')
    lines = mode.readlines()
    # input pdb
    input_pdb = open(sys.argv[2], 'r')
    # lines of one frame
    frame_lines = 814
    # from 37 to 75
    frame_index = int(sys.argv[3])
    # output pdb
    output_pdb = open(sys.argv[4], 'w')
    start = frame_index * frame_lines
    end = start + frame_lines
    for i in range(start, end):
        line = lines[i]
        print(line, end="")
    pass