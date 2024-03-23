
import argparse
import statsmodels.api as sm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', type=str, help='Input tab-separated file with gene_id (column 1) and protein size (column 2).')
parser.add_argument('-t', '--title', type = str, default = 'ECDF of provided data', help = 'Title of produced ECDF Figure.' )
parser.add_argument('-xlab', '--xlabel', type = str, default = 'Protein Size (amino-acids)', help = 'Label of x-axis of produced ECDF Figure.' )
parser.add_argument('-ylab', '--ylabel', type = str, default = 'Empirical Cumulative Distribution Function', help = 'Label of y-axis of produced ECDF Figure.' )
parser.add_argument('-s', '--show', action = 'store_true', default = False, help = 'Option to show the produced ECDF Figure.' )
parser.add_argument('-o', '--output', type = str, help = 'Name of output file.' )
parser.add_argument('-f', '--format', type = str, default = "svg", help = 'Format of output file.' )
args = parser.parse_args()

if not any(vars(args).values()):
    parser.print_help()
    sys.exit('Error: No arguments provided.')


def ecdf(values, title, output, format, xlabel = 'Protein Size (amino-acids)', ylabel = 'Empirical Cumulative Distribution Function', show = False):
    df = pd.read_csv(values, sep='\t')
    data = df.iloc[:, 1]  # Keep only the second column that contains protein sizes

    sorted_data = np.sort(data)
    y_values = np.arange(1, len(sorted_data) + 1) / float(len(sorted_data))

    plt.rcParams['font.family'] = 'arial'
    plt.plot(sorted_data, y_values, marker='.', linestyle='none', color='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(False)
    plt.savefig( re.sub(" ", "_", output+".")+format, format=format)

    if show == True:
        plt.show()


if __name__ == '__main__':    
    ecdf (values = args.data,
          title = args.title,
          xlabel = args.xlabel if args.xlabel else 'Protein Size (amino-acids)',
          ylabel = args.ylabel if args.ylabel else 'Empirical Cumulative Distribution Function',
          show = args.show,
          output = args.output if args.output else args.title,
          format = args.format if args.format else 'svg'
         )
