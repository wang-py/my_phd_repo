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

def stat_plot(time, temp, usage, title):
    """
    function that plots the stats of selected hardware
    ----------------------------------------------------------------------------
    time: ndarray    
    time of monitoring

    temp: ndarray
    temperature data

    usage: ndarray
    usage in percentage

    title: str
    title of the trial

    ----------------------------------------------------------------------------

    """
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle(title)
    ax1 = ax[0]
    ax1.scatter(time, temp, label="temperature", color='b')
    ax1.set_xlabel("time [S]")
    ax1.set_ylabel("temperature [C]")

    ax2 = ax[1]
    # window for moving average
    N = 10
    usage_mean = np.convolve(gpu_usage, np.ones((N,))/N, mode = 'same')
    ax2.scatter(time, usage, label="usage", s=1)
    ax2.plot(time[N:-N], usage_mean[N:-N], label="moving average", linestyle='-', color='r')
    ax2.set_xlabel("time [S]")
    ax2.set_ylabel("usage [%]")
    ax2.legend()

    pass    

if __name__ == "__main__":
    gpu_csv = sys.argv[1]
    fig_title = gpu_csv.split('.')[0].split('/')[-1]
    gpu_temp, gpu_usage = read_gpu_csv(gpu_csv)
    N = len(gpu_temp)
    time_arr = np.arange(N)
    stat_plot(time_arr, gpu_temp, gpu_usage, fig_title)

    plt.show()

    pass