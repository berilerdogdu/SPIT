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
import warnings
import scipy.stats as sts
from collections import defaultdict
from tqdm import tqdm


def split_distributions(m, split_if):
    m_left, m_right = m.astype(np.float32), m.astype(np.float32)
    m_right[m_right <= split_if] = np.nan
    m_left[m_left > split_if] = np.nan
    return m_left, m_right

def filter_on_sample_size(m_1, m_2, n_small):
    m_1_f, m_2_f = m_1, m_2
    boolean = ((~np.isnan(m_1)).sum(axis=0)<n_small)|((~np.isnan(m_2)).sum(axis=0)<n_small)
    m_1_f[:, boolean] = np.nan
    m_2_f[:, boolean] = np.nan
    return m_1_f, m_2_f

def jitter_matrix(m):
    m = m + np.random.uniform(-0.0005, 0.0005, m.shape)
    m = np.clip(m, 0, 1)
    return m

def run_mann_whitney_u(m_1, m_2):
    m_1 = jitter_matrix(m_1)
    m_2 = jitter_matrix(m_2)
    warnings.filterwarnings("ignore", message=".*small.*sample.*")
    mwu_results = sts.mannwhitneyu(m_1, m_2, method = 'auto', alternative='two-sided', axis = 0, nan_policy = 'omit')
    if(np.isnan(mwu_results[1]).all()):
        min_p_val = 1
    else:
        min_p_val = np.nanmin(mwu_results[1])
    return mwu_results[1], round(min_p_val, 16)

def compare_tails(m_1_left, m_1_right, m_2_left, m_2_right):
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

def get_random_halves(matrix, s):
    np.random.shuffle(matrix)
    m_1 = np.sort(matrix[0:s,], axis=0)
    m_2 = np.sort(matrix[s:,], axis=0)
    return m_1, m_2

def complete_w_rand_samples(m, target_sample_size):
    num_of_rand_samples = target_sample_size - m.shape[0]
    rand_samples = np.random.uniform(0, 1, (max(0, num_of_rand_samples), m.shape[1]))
    completed_m = np.concatenate((m, rand_samples), axis=0)
    return completed_m

def mannwhitneyu_permutation(IFs, samples, case_sample_size, num_of_it, n_small, quiet=False):
    min_perm_p_arr = []
    perm_p_arr = []
    sampled_txs = set()
    if not quiet:
        print("Chewing in progress:")
    for it in range(1, num_of_it+1):
        if not quiet:
            out = '\r'+ str(it) + " iterations completed."
            if (it == num_of_it):
                print(out)
            else:
                print(out, end='')
        ctrl_sample_size = len(samples)
        s = int(ctrl_sample_size/2)
        tx_matrix = np.transpose(IFs[samples].to_numpy()).astype(np.float32)
        m_1, m_2 = get_random_halves(tx_matrix, s)
        if(s <= 16):
            m_1 = complete_w_rand_samples(m_1, ctrl_sample_size)
            m_2 = complete_w_rand_samples(m_2, case_sample_size)
            p_arr, min_it_p_val = run_mann_whitney_u(m_1, m_2)
            min_it_p_val = np.random.choice(np.sort(p_arr)[0:int(m_1.shape[1]*0.01)])
        elif(s>n_small):
            rand_if = random.uniform(0, 1)
            m_1_left, m_1_right = split_distributions(m_1, rand_if)
            m_2_left, m_2_right = split_distributions(m_2, rand_if)
            m_1_left_f, m_2_left_f = filter_on_sample_size(m_1_left, m_2_left, n_small)
            m_1_right_f, m_2_right_f = filter_on_sample_size(m_1_right, m_2_right, n_small)
            min_it_p_val, p_arr = compare_tails(m_1_left_f, m_1_right_f, m_2_left_f, m_2_right_f)
        else:
            p_arr, min_it_p_val = run_mann_whitney_u(m_1, m_2)
        if(min_it_p_val == 1 | len(p_arr) == 0):
            pass
        perm_p_arr.append(p_arr)
        if (s > 16):
            i = 0
            min_tx_ind = np.nanargmin(p_arr)
            while((min_tx_ind in sampled_txs) and (len(sampled_txs) < tx_matrix.shape[1])):
                i+=1
                min_tx_ind = np.argsort(p_arr)[i]
                if(np.isnan(p_arr[min_tx_ind])):
                    break
            sampled_txs.add(min_tx_ind)
            min_it_p_val = p_arr[min_tx_ind]
        min_perm_p_arr.append(min_it_p_val)
    perm_p_medians = pd.DataFrame(np.transpose(np.array(perm_p_arr))).set_index(IFs.index).median(axis = 1)
    return min_perm_p_arr, perm_p_medians

def write_min_p_values(min_perm_p_arr, p_values_file):
    with open(p_values_file, 'w') as f:
        for i in min_perm_p_arr:
            i = format(i, '.16f')
            f.write("%s\n" % i)
    return

def main(args):
    np.random.seed(42)
    random.seed(42)
    IFs = pd.read_csv(args.i, sep='\t', index_col=0)
    if_cols = IFs.select_dtypes(include=[np.number]).columns
    IFs[if_cols] = IFs[if_cols].astype(np.float32).round(3)
    pheno = pd.read_csv(args.l, sep='\t')
    samples = pheno[pheno.condition == 0].id.to_list()
    case_sample_size = pheno[pheno.condition == 1].shape[0]
    min_perm_p_arr, perm_p_medians = mannwhitneyu_permutation(IFs, samples, case_sample_size, args.n_iter, args.n_small, args.quiet)
    if(args.exp):
        write_dir = args.exp
    else:
        write_dir = os.path.join(args.O, "SPIT_analysis")
    perm_p_medians.to_csv(os.path.join(write_dir, "perm_p_medians.txt"), sep = '\t')
    write_min_p_values(min_perm_p_arr, os.path.join(write_dir, "spit_test_min_p_values.txt"))
