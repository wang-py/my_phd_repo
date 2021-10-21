import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# Script that reads the csv output of glances and nvidia-smi and plots the 
# results 

def read_gpu_csv(gpu_csv):
    gpu_df = pd.read_csv(gpu_csv)
    gpu_temp = gpu_df.loc[:, " temperature.gpu"].to_numpy()
    gpu_usage = gpu_df.loc[:, " utilization.gpu [%]"].str[:-1].astype(int).to_numpy()

    
    return gpu_temp, gpu_usage

def stat_plot(time, temp):
    """
    function that plots the stats of selected hardware
    ----------------------------------------------------------------------------
    time: ndarray

    time of 
    """
    fig, ax1 = plt.subplots()
    ax1.scatter(time, temp, label="temperature", color='b')
    ax1.set_xlabel("time [S]")
    ax1.set_ylabel("temperature [C]")

    ax1.legend()

    pass    

if __name__ == "__main__":
    gpu_csv = sys.argv[1]
    gpu_temp, gpu_usage = read_gpu_csv(gpu_csv)
    time_arr = np.arange(len(gpu_temp))
    stat_plot(time_arr, gpu_temp)

    plt.show()

    pass