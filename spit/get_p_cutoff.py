#!/usr/bin/env python3

"""
Description:
    Retrieves a p-value threshold to be used in detecting significant DTU events. This module determines p-value threshold based on the user-defined parameter K and is obtained as the "K x N"-th smallest p-value among the N sampled by SPIT-Test.
Usage:
    ./get_p_cutoff.py -k <K> -p <sampled_p_values>
Author:
    Beril Erdogdu
Date:
    July 08, 2023
"""

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

def main(args):

    p_values = pd.read_csv(args.p, names = ['p']).p.to_list()
    p_cutoff = get_p_cutoff(p_values, args.k)
    return p_cutoff
