import pandas as pd
import matplotlib.pyplot as plt
import sys
# Script that reads the csv output of caver and plots tunnel profiles

if __name__ == "__main__":
    # input csv file
    input_csv = sys.argv[1]
    # reading csv into a pandas dataframe
    profile_df = pd.read_csv(input_csv, skiprows=26, nrows=1)
    R_range = profile_df.iloc[0, 13:]
    pts_range = range(len(R_range))
    plt.figure()
    plt.scatter(pts_range, R_range)
    plt.xlabel("points")
    plt.ylabel("Radius [A]")
    plt.title("Q10 size along caver points")
    plt.show()
    pass