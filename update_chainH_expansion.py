import numpy as np
import sys

if __name__ == "__main__":
    # mode file
    mode = open(sys.argv[1], 'r')
    # lines of one frame
    frame_lines = 814
    for i in range(frame_lines):
        one_line = mode.readline()
        print(one_line, end="")
    pass