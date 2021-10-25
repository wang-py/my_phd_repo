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

def stat_plot_single_GPU(time, temp, usage, title):
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

def stat_plot_4GPU(time, temp, usage, title):
    """
    function that plots the stats of 4 GPUs
    ----------------------------------------------------------------------------
    """
    fig, ax = plt.subplots(1, 4, figsize=(24, 6))
    fig.suptitle(title)
    fig2, ax2 = plt.subplots(1, 4, figsize=(24, 6))
    fig2.suptitle(title)
    N = 10
    for i in range(4):
        time_i = time[i::4]
        temp_i = temp[i::4]
        usage_i = usage[i::4]
        ax[i].scatter(time_i, usage_i, label="usage of GPU %d"%i, s=1)
        usage_mean_i = np.convolve(usage_i, np.ones((N,))/N, mode = 'same')
        ax[i].plot(time_i[N:-N], usage_mean_i[N:-N], label="moving average", linestyle='-', color='r')
        ax[i].legend()
        ax[i].set_xlabel("time [S]")
        ax[i].set_ylabel("usage [%]")
        ax2[i].scatter(time_i, temp_i, label="temperature of GPU %d"%i)
        ax2[i].legend()
        ax2[i].set_xlabel("time [S]")
        ax2[i].set_ylabel("temperature [C]")

    pass

if __name__ == "__main__":
    gpu_csv = sys.argv[1]
    fig_title = gpu_csv.split('.')[0].split('/')[-1]
    gpu_temp, gpu_usage = read_gpu_csv(gpu_csv)
    N = len(gpu_temp)
    time_arr = np.arange(N)
    if fig_title[0] == '1':
        stat_plot_single_GPU(time_arr, gpu_temp, gpu_usage, fig_title)
    else:
        stat_plot_4GPU(time_arr, gpu_temp, gpu_usage, fig_title)

    plt.show()

    pass