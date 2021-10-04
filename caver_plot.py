import pandas as pd
import matplotlib.pyplot as plt
import sys
# Script that reads the csv output of caver and plots tunnel profiles

def plotting(x, y, label):
    plt.figure()
    plt.scatter(x, y, label=label)
    plt.xlabel("points")
    plt.ylabel("Radius [A]")
    plt.title("sphere size along caver points")

def select_tunnel_to_plot(tunnel_index, input_csv):
    profile_df = pd.read_csv(input_csv, skiprows= 7 * tunnel_index - 2, nrows=1)
    R_range = profile_df.iloc[0, 13:]

    return R_range

if __name__ == "__main__":
    # input csv file
    input_csv = sys.argv[1]
    # reading csv into a pandas dataframe
    R_range_Q =  select_tunnel_to_plot(4, input_csv)
    pts_range = range(len(R_range_Q))
    plotting(pts_range, R_range_Q, "Q10")
    plt.legend()
    plt.show()
    pass