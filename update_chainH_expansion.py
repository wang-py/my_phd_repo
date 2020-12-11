import numpy as np
import sys

if __name__ == "__main__":
    # mode file
    mode = open(sys.argv[1], 'r')
    # lines of one frame
    frame_lines = 815
    i = 0
    one_frame = mode.readlines(frame_lines * i)
    print(one_frame)
    pass