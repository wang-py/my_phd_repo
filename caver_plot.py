import pandas as pd
import matplotlib.pyplot as plt
import sys
# Script that reads the csv output of caver and plots tunnel profiles

if __name__ == "__main__":
    # input csv file
    input_csv = sys.argv[1]
    # reading csv into a pandas dataframe
    profile_df = pd.read_csv(input_csv)
    pass