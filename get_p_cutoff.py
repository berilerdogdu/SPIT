import sys, getopt
import numpy as np
import pandas as pd

def get_p_cutoff(p_values, k):
    n = len(p_values)
    t = int(k*n)
    p_values = np.array(p_values)
    if(t==0):
      t = 0
    else:
      t = t-1
    p_cutoff = np.sort(p_values)[t]
    return p_cutoff

def main(argv):

    p_values_file = ''
    k = ''

    try:
        opts, args = getopt.getopt(argv,"hk:p:",["spit_test_p_values_file=", "k="])
    except getopt.GetoptError:
        print("Usage: python get_p_cutoff.py -k <K> -p <spit test p_values file>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: python get_p_cutoff.py -k <K> -p <spit test p_values file>")
            sys.exit()
        elif opt in ("-k", "--K"):
            k = float(arg)
        elif opt in ("-p", "--spit_test_p_values_file"):
            p_values_file = arg
    
    p_values = pd.read_csv(p_values_file, names = ['p']).p.to_list()
    p_cutoff = get_p_cutoff(p_values, k)
    print(p_cutoff)



if __name__ == "__main__":
   main(sys.argv[1:])
