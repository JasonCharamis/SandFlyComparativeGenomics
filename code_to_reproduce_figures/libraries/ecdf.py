
import argparse
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', type=str, help='Input tab-separated file with gene_id (column 1) and protein size (column 2).')
parser.add_argument('-t', '--title', type = str, default = 'ECDF of provided data', help = 'Title of produced ECDF Figure.' )
parser.add_argument('-s', '--show', action = 'store_true', default = False, help = 'Option to show the produced ECDF Figure.' )
parser.add_argument('--show-help', action = 'help', default = False, help = 'Option to show help message.' )
args = parser.parse_args()

if not any(vars(args).values()):
    parser.print_help()
    sys.exit('Error: No arguments provided.')


def ecdf(values, title, show):
    df = pd.read_csv(values, sep='\t')
    data = df.iloc[:, 1]  # Keep only the second column that contains protein sizes

    sorted_data = np.sort(data)
    y_values = np.arange(1, len(sorted_data) + 1) / float(len(sorted_data))

    plt.rcParams['font.family'] = 'arial'
    plt.plot(sorted_data, y_values, marker='.', linestyle='none', color='black')
    plt.xlabel('Gene expression')
    plt.ylabel('Empirical Cumulative Distribution Function')
    plt.title(title)
    plt.grid(False)
    plt.savefig( args.title+'.svg', format='svg')

    if show == True:
        plt.show()


if __name__ == '__main__':    
    ecdf (values = args.data, title = args.title, show = args.show)
