import matplotlib.pyplot as plt
import pandas as pd
import sys

# Script that reads the csv output of glances and nvidia-smi and plots the 
# results 

def read_cpu_csv(cpu_csv):
    cpu_df = pd.read_csv(cpu_csv)
    pass

def read_gpu_csv(gpu_csv):
    gpu_df = pd.read_csv(gpu_csv)
    pass

def stat_plot(time, temp, usage):
    """
    function that plots the stats of selected hardware
    ----------------------------------------------------------------------------
    time: ndarray

    time of 
    """
    fig, ax1 = plt.subplots()
    ax1.scatter(time, temp, label="temperature")
    ax1.set_xlabel("time [S]")
    ax1.set_ylabel("temperature [C]")

    ax2 = ax1.twinx()
    ax2.scatter(time, usage, label="usage")
    ax2.set_ylabel("usage [%]")

    ax1.legend()
    ax2.legend()

    pass    

if __name__ == "__main__":
    cpu_csv = sys.argv[1]
    gpu_csv = sys.argv[2]

    pass