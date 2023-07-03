import sys
import argparse
import numpy as np
import pandas as pd

def get_p_cutoff(p_values, k):
    n = len(p_values)
    t = int(k*n)
    p_values = np.array(p_values)
    if(t>0):
      t = t-1
    p_cutoff = np.sort(p_values)[t]
    return p_cutoff

def main(argv):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-k', metavar='0.6', type=float, default=0.6, help='-K hyperparameter for p-value thresholding')
    parser.add_argument('-p', metavar='spit_test_min_p_values.txt', required=True, type=str, help='Minimum p-values from SPIT Test iterations')

    args = parser.parse_args()
    p_values = pd.read_csv(args.p, names = ['p']).p.to_list()
    p_cutoff = get_p_cutoff(p_values, args.k)
    print(p_cutoff)


if __name__ == "__main__":
   main(sys.argv[1:])
