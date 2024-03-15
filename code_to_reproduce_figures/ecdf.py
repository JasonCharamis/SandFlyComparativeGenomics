
import argparse
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

parser = argparse.ArgumentParser()
parser.add_argument('data')
args = parser.parse_args()

def ecdf(values):
    df = pd.read_csv(values, sep="\t")
    data = df.iloc[:, 1]  # Keep only the second column that contains protein sizes

    sorted_data = np.sort(data)
    y_values = np.arange(1, len(sorted_data) + 1) / float(len(sorted_data))

    plt.rcParams['font.family'] = "arial"
    plt.plot(sorted_data, y_values, marker='.', linestyle='none', color='black')
    plt.xlabel('Gene expression')
    plt.ylabel('Empirical Cumulative Distribution Function')
    plt.title('ECDF of 379 sand fly CCEs')
    plt.grid(False)
    # plt.show()
    plt.savefig('Figure_S10_379_CCEs.ecdf.svg', format="svg")
    
ecdf ( args.data )
