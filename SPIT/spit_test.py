#!/usr/bin/env python3

"""
Description:
    This module performs the SPIT-Test. The test randomly splits the control group in half, and identifies the most significant difference in isoform fractions between the two halves. The observations from this process can then be used to compare candidate DTU events in terms of their significance.
Usage:
    ./spit_test.py -i <filtered_ifs> -g <filtered_gene_counts> -l <pheno> -I <dominance_selected_ifs> -G <dominance_selected_gene_counts> -P <sampled_p_values>
Author:
    Beril Erdogdu
Date:
    July 08, 2023
"""

import os
import numpy as np
import pandas as pd
import random
import math
import scipy.stats as sts
from collections import defaultdict
from tqdm import tqdm


def split_distributions(m, split_if):
    m_left, m_right = m.astype(float), m.astype(float)
    m_right[m_right <= split_if] = np.nan
    m_left[m_left > split_if] = np.nan
    return m_left, m_right

def filter_on_sample_size(m_1, m_2, n_small):
    m_1_f, m_2_f = m_1, m_2
    boolean = ((~np.isnan(m_1)).sum(axis=0)<n_small)|((~np.isnan(m_2)).sum(axis=0)<n_small)
    m_1_f[:, boolean] = np.nan
    m_2_f[:, boolean] = np.nan
    return m_1_f, m_2_f

def run_mann_whitney_u(m_1, m_2):
    mwu_results = sts.mannwhitneyu(m_1, m_2, method = 'auto', alternative='two-sided', axis = 0, nan_policy = 'omit')
    if(np.isnan(mwu_results[1]).all()):
        min_p_val = 1
    else:
        min_p_val = np.nanmin(mwu_results[1])
    return mwu_results[1], min_p_val

def compare_tails(m_1_left, m_1_right, m_2_left, m_2_right, if_):
    min_it_p_val = 1
    p_arr = []
    if(~(np.isnan(m_1_left).all()) & ~(np.isnan(m_1_right).all())):
        mwu_results_left, p_val_left = run_mann_whitney_u(m_1_left, m_2_left)
        mwu_results_right, p_val_right = run_mann_whitney_u(m_1_right, m_2_right)
        if(p_val_left < p_val_right):
            min_it_p_val = p_val_left
            p_arr = mwu_results_left
        else:
            min_it_p_val = p_val_right
            p_arr = mwu_results_right
    elif(~(np.isnan(m_1_left).all())):
        mwu_results_left, p_val_left = run_mann_whitney_u(m_1_left, m_2_left)
        min_it_p_val = p_val_left
        p_arr = mwu_results_left
    elif(~(np.isnan(m_1_right).all())):
        mwu_results_right, p_val_right = run_mann_whitney_u(m_1_right, m_2_right)
        min_it_p_val = p_val_right
        p_arr = mwu_results_right
    return min_it_p_val, p_arr


def mannwhitneyu_permutation(IFs, ctrl_samples, num_of_it, n_small):
    min_perm_p_arr = []
    perm_p_arr = []
    sampled_txs = set()
    print("Chewing in progress:")
    for it in range(1, num_of_it+1):
        out = '\r'+ str(it) + " iterations completed."
        if (it == num_of_it):
            print(out)
        else:
            print(out, end='')
        s = int(len(ctrl_samples)/2)
        tx_matrix = np.transpose(IFs[ctrl_samples].to_numpy()).astype(float)
        np.random.shuffle(tx_matrix)
        m_1 = np.sort(tx_matrix[0:s,], axis=0)
        m_2 = np.sort(tx_matrix[s:,], axis=0)
        if(s>=n_small):
            rand_if = random.uniform(0, 1)
            m_1_left, m_1_right = split_distributions(m_1, rand_if)
            m_2_left, m_2_right = split_distributions(m_2, rand_if)
            m_1_left_f, m_2_left_f = filter_on_sample_size(m_1_left, m_2_left, n_small)
            m_1_right_f, m_2_right_f = filter_on_sample_size(m_1_right, m_2_right, n_small)
            min_it_p_val, p_arr = compare_tails(m_1_left_f, m_1_right_f, m_2_left_f, m_2_right_f, rand_if)
        else:
            p_arr, min_it_p_val = run_mann_whitney_u(m_1, m_2)
        if(min_it_p_val == 1 | len(p_arr) == 0):
            pass
        perm_p_arr.append(p_arr)
        min_tx_ind = np.nanargmin(p_arr)
        i = 0
        while((min_tx_ind in sampled_txs) and (len(sampled_txs) < tx_matrix.shape[1])):
            i+=1
            min_tx_ind = np.argsort(p_arr)[i]
            if(np.isnan(p_arr[min_tx_ind])):
                break
        if(min_tx_ind in sampled_txs):
            pass
        sampled_txs.add(min_tx_ind)
        min_it_p_val = p_arr[min_tx_ind]
        min_perm_p_arr.append(min_it_p_val)
    perm_p_medians = pd.DataFrame(np.transpose(np.array(perm_p_arr))).set_index(IFs.index).median(axis = 1)
    return min_perm_p_arr, perm_p_medians

def determine_perm_p_cutoff(min_perm_p_arr, p_values_file):
    p_cutoff = np.min(np.array(min_perm_p_arr))
    with open(p_values_file, 'w') as f:
        for i in min_perm_p_arr:
            f.write("%s\n" % i)
    return p_cutoff

def main(args):

    IFs = pd.read_csv(args.i, sep='\t', index_col=0)
    tx_names = list(IFs.index)
    gene_names = list(IFs.gene_id.unique())
    IFs = IFs.round(2)
    gene_counts = pd.read_csv(args.g, sep='\t', index_col=0)
    pheno = pd.read_csv(args.l, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    IFs.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_ifs.txt"), sep = '\t')
    gene_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_gene_counts.txt"), sep = '\t')
    min_perm_p_arr, perm_p_medians = mannwhitneyu_permutation(IFs, ctrl_samples, args.n_iter, args.n_small)
    perm_p_medians.to_csv(os.path.join(args.O, "SPIT_analysis", "perm_p_medians.txt"), sep = '\t')
    p_cutoff = determine_perm_p_cutoff(min_perm_p_arr, os.path.join(args.O, "SPIT_analysis", "spit_test_min_p_values.txt"))
