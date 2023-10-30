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

import sys
import os
import argparse
import numpy as np
import pandas as pd
import random
import math
import scipy.stats as sts
from collections import defaultdict
from tqdm import tqdm


def select_genes_w_dom_iso(IFs, ctrl_samples, gene_names, p_dom):
    IFs.insert(0, "ctrl_IF_mean", IFs[ctrl_samples].mean(axis = 1))
    ctrl_IF_max = IFs.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_no_max = IFs.drop(ctrl_IF_max.index)
    ctrl_IF_second_max = ctrl_IF_no_max.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_min = IFs.sort_values('ctrl_IF_mean', ascending=True).drop_duplicates(['gene_id'])
    genes_w_dom_iso_counter = 0
    final_dom_iso_genes = set()
    print("Selecting genes with consistently dominant isoforms in control group:")
    for g in tqdm(gene_names):
        min_if = ctrl_IF_min[ctrl_IF_min.gene_id == g].ctrl_IF_mean
        min_iso = ctrl_IF_min[ctrl_IF_min.gene_id == g].index
        max_if = ctrl_IF_max[ctrl_IF_max.gene_id == g].ctrl_IF_mean
        max_iso = ctrl_IF_max[ctrl_IF_max.gene_id == g].index
        second_max_iso = ctrl_IF_second_max[ctrl_IF_second_max.gene_id == g].index
        if(np.sum(IFs.loc[max_iso, ctrl_samples].to_numpy() > IFs.loc[second_max_iso, ctrl_samples].to_numpy()) < (len(ctrl_samples) * p_dom
        )):
            pass
        else:
            final_dom_iso_genes.add(g)
            genes_w_dom_iso_counter += 1
    return list(final_dom_iso_genes)

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
    pd.DataFrame(np.transpose(np.array(perm_p_arr))).set_index(IFs.index).median(axis = 1).to_csv('perm_p_medians.txt', sep = '\t')
    return min_perm_p_arr

def determine_perm_p_cutoff(min_perm_p_arr, p_values_file):
    p_cutoff = np.min(np.array(min_perm_p_arr))
    with open(p_values_file, 'w') as f:
        for i in min_perm_p_arr:
            f.write("%s\n" % i)
    return p_cutoff

def main(argv):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-i', metavar='filtered_ifs.txt', type=str, default=os.path.join(os.getcwd(), "SPIT_analysis", "filtered_ifs.txt"), help='Isoform fractions file (tsv)')
    parser.add_argument('-g', metavar='filtered_gene_counts.txt', type=str, default = os.path.join(os.getcwd(), "SPIT_analysis", "filtered_gene_counts.txt"), help='Gene counts file (tsv)')
    parser.add_argument('-l', metavar='labels.txt', required=True, type=str, help='Labels/metadata file (tsv)')
    parser.add_argument('--n_iter', metavar='100', type=int, default=100, help='Number of iterations')
    parser.add_argument('-d', '--p_dom', metavar='0.75', type=float, default=0.75, help='Dominance selection threshold')
    parser.add_argument('-n', '--n_small', metavar='12', type=int, default=12, help='Smallest sample size for the subgroups')
    parser.add_argument('-O', '--output_dir', type=str, default = os.getcwd(), help = "Output directory path where the SPIT output folder will be written")
    args = parser.parse_args()
    IFs = pd.read_csv(args.i, sep='\t', index_col=0)
    tx_names = list(IFs.index)
    gene_names = list(IFs.gene_id.unique())
    IFs = IFs.round(2)
    gene_counts = pd.read_csv(args.g, sep='\t', index_col=0)
    pheno = pd.read_csv(args.l, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    selected_dom_iso_genes = select_genes_w_dom_iso(IFs, ctrl_samples, gene_names, args.p_dom)
    IFs_selected = IFs[IFs.gene_id.isin(selected_dom_iso_genes)]
    IFs_selected.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_ifs.txt"), sep = '\t')
    gene_counts_selected = gene_counts[gene_counts.index.isin(selected_dom_iso_genes)]
    gene_counts_selected.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_gene_counts.txt"), sep = '\t')
    min_perm_p_arr = mannwhitneyu_permutation(IFs_selected, ctrl_samples, args.n_iter, args.n_small)
    p_cutoff = determine_perm_p_cutoff(min_perm_p_arr, os.path.join(args.O, "SPIT_analysis", "spit_test_min_p_values.txt"))


if __name__ == "__main__":
   main(sys.argv[1:])
