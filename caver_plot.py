import pandas as pd
import matplotlib.pyplot as plt
import sys
# Script that reads the csv output of caver and plots tunnel profiles

def plotting(x, y, label):
    """
    plotting function
    ---------------------------------------------------------------------------
    x: ndarray
    x axis of the plot

    y: ndarray
    y axis of the plot

    label: str
    label of the data
    """
    #plt.figure()
    plt.scatter(x, y, label=label)
    plt.xlabel("points")
    plt.ylabel("Radius [A]")
    plt.title("sphere size along caver points")

def select_tunnel_to_plot(tunnel_index, input_csv):
    """
    search in csv file for the radius data of a selected tunnel
    ---------------------------------------------------------------------------
    tunnel_index: int
    index of the selected tunnel

    input_csv: str
    name of the input csv file

    ---------------------------------------------------------------------------
    Returns:

    R_range: ndarray
    array of radius values of selected tunnel

    """
    profile_df = pd.read_csv(input_csv, skiprows= 7 * tunnel_index - 2, nrows=1)
    R_range = profile_df.iloc[0, 13:]

    return R_range

if __name__ == "__main__":
    # input csv file
    input_csv_Q = "test_files/tunnel_profiles_Q.csv"
    input_csv_cavity = "test_files/tunnel_profiles_cavity.csv"
    # reading csv into a pandas dataframe
    R_range_Q =  select_tunnel_to_plot(4, input_csv_Q)
    R_range_cavity =  select_tunnel_to_plot(1, input_csv_cavity)
    pts_range_Q = range(len(R_range_Q))
    pts_range_cavity = range(len(R_range_cavity))
    plotting(pts_range_Q, R_range_Q, "Q10")
    plotting(pts_range_cavity, R_range_cavity, "Q_cavity")
    plt.legend()
    plt.show()
    pass