#!/usr/bin/env python3

"""
Description:
    This module fits KDE on Isoform Fraction (IF) distributions to search for separation between- case and control samples. P-value threshold is then used to select significant events out of all candidate DTU events.
Usage:
    ./dtu_detection.py -i <dominance_selected_ifs.txt> -g <dominance_selected_gene_counts.txt> -m <tx2gene> -l <pheno> -M <spit_cluster_matrix.txt> -O <spit_out.txt>
Author:
    Beril Erdogdu
Date:
    July 08, 2023
"""

import os
import pandas as pd
import numpy as np
import scipy.stats as sts
from collections import Counter
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
from matplotlib import pyplot as plt
from spit.transform_tx_counts_to_ifs import convert_counts_to_IF_and_gene_level as convert_counts_to_IF
from tqdm import tqdm


def convert_counts_to_cpm(counts):
    lib_sizes = counts.sum(axis = 0).to_list()
    scaled_counts = counts.divide(lib_sizes)
    cpms = (scaled_counts * (10**6)).astype(np.float32)
    return cpms

def filter_samples_on_gene_cpm(f_cpm, tx, tx2gene_dict, samples, ifs, cpms):
    gene_id = tx2gene_dict[tx]
    if f_cpm:
        cpm_vals = cpms.loc[gene_id, samples].values
        keep_mask = cpm_vals > 10
        selected_samples = [s for s, k in zip(samples, keep_mask) if k]
        filtered_samples = [s for s, k in zip(samples, keep_mask) if not k]
    else:
        selected_samples = samples
        filtered_samples = []
    selected_if_arr = ifs.loc[tx, selected_samples].to_numpy(copy=True, dtype=np.float32)
    return np.array(selected_samples), selected_if_arr, filtered_samples
    
def groupsize_limit_check(ctrl, case, size_lim):
    c = 1
    if((len(ctrl) < size_lim) or (len(case) < size_lim)):
        c = 0
    return c

def index_if_arrays(arr):
    sort_ind_arr = np.argsort(arr)
    sorted_arr = arr[sort_ind_arr]
    return sort_ind_arr, sorted_arr
    
def split_vector_on_kernel(sorted_vector, b):
    p = np.linspace(0,1).reshape(-1, 1)
    kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(sorted_vector.reshape(-1, 1))
    scores = kde.score_samples(p)
    mi = argrelextrema(scores, np.less)[0]
    split_ind = len(sorted_vector) - 1
    split_val = -1
    if(len(mi) > 0):
        split_val = p[mi[np.argmin(scores[mi])]]
    return split_val
        
def subgroupsize_limit_check(group, split_in, size_lim):
    if((split_in >= size_lim-1) and ((len(group)-split_in) > size_lim)):
        return 1
    return 0

def set_ctrl_tails(group_arr, split_in, size_lim):
    left_tail = group_arr[0:split_in]
    right_tail = group_arr[split_in:]
    if(split_in < size_lim-1):
        left_tail = group_arr[0:size_lim]
    if((len(group_arr)-split_in) <= size_lim):
        right_tail = group_arr[-size_lim:]
    return left_tail, right_tail

def jitter_matrix(m):
    m = m + np.random.uniform(-0.0005, 0.0005, m.shape)
    m = np.clip(m, 0, 1)
    return m
    
def fit_kde_and_mannwhitney(ctrl_subgroup, case_subgroup, b):
    ctrl_kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(ctrl_subgroup.reshape(-1, 1))
    likelihood = ctrl_kde.score(case_subgroup.reshape(-1, 1), y = None)
    case_subgroup = jitter_matrix(case_subgroup)
    ctrl_subgroup = jitter_matrix(ctrl_subgroup)
    r = sts.mannwhitneyu(case_subgroup, ctrl_subgroup, method = 'auto', alternative='two-sided')
    return round(r[1], 16), likelihood

def compute_double_stats(ifs, gene_counts, tx2gene_dict, ctrl, case, cluster_size_lim, perm_p_cutoff, b, f_cpm, quiet=False):
    likelihood_arr = []
    wrst_p_arr = []
    cpms = convert_counts_to_cpm(gene_counts)
    genotype_cluster_df = pd.DataFrame(0, index=ifs.index, columns=case, dtype=np.int8)
    if not quiet:
        print("Detecting potential DTU transcripts:")
    iterator = ifs.index.to_list()
    for i in (tqdm(iterator) if not quiet else iterator):
        selected_ctrls, ctrl_arr, filtered_ctrls = filter_samples_on_gene_cpm(f_cpm, i, tx2gene_dict, ctrl, ifs, cpms)
        selected_cases, case_arr, filtered_cases = filter_samples_on_gene_cpm(f_cpm, i, tx2gene_dict, case, ifs, cpms)
        genotype_cluster_df.loc[i, filtered_cases] = -1
        if((groupsize_limit_check(ctrl_arr, case_arr, cluster_size_lim)) == 0):
            wrst_p_arr.append(1)
            likelihood_arr.append(0)
            continue
        sort_ind_ctrl_arr, sorted_ctrl_arr = index_if_arrays(ctrl_arr)
        sort_ind_case_arr, sorted_case_arr = index_if_arrays(case_arr)
        case_split_val = split_vector_on_kernel(sorted_case_arr, b)
        if(case_split_val > -1):
            case_split_ind = np.nonzero(sorted_case_arr <= case_split_val)[0][-1]
            if(len(np.nonzero(sorted_ctrl_arr <= case_split_val)[0])>0):
                ctrl_split_ind = np.nonzero(sorted_ctrl_arr <= case_split_val)[0][-1]
            else:
                ctrl_split_ind = len(sorted_ctrl_arr) - 1
        else:
            case_split_ind = len(sorted_case_arr) - 1
            ctrl_split_ind = len(sorted_ctrl_arr) - 1
        case_tx_likelihood = 0        
        if((subgroupsize_limit_check(sorted_case_arr, case_split_ind, cluster_size_lim)) == 0):
            r, likelihood = fit_kde_and_mannwhitney(sorted_ctrl_arr, sorted_case_arr, b)
            if(r <= perm_p_cutoff):
                genotype_cluster_df.loc[i, selected_cases] = 1
            likelihood_arr.append(likelihood)
            wrst_p_arr.append(r)
        else:
            case_left_tail = sorted_case_arr[0:case_split_ind]
            case_left_tail_samps = selected_cases[sort_ind_case_arr[0:case_split_ind]]
            case_right_tail = sorted_case_arr[case_split_ind:]
            case_right_tail_samps = selected_cases[sort_ind_case_arr[case_split_ind:]]
            ctrl_left_tail, ctrl_right_tail = set_ctrl_tails(sorted_ctrl_arr, ctrl_split_ind, cluster_size_lim)
            r_left, likelihood_left = fit_kde_and_mannwhitney(ctrl_left_tail, case_left_tail, b)
            r_right, likelihood_right = fit_kde_and_mannwhitney(ctrl_right_tail, case_right_tail, b)
            if((r_left <= perm_p_cutoff) and (r_right <= perm_p_cutoff)):
                genotype_cluster_df.loc[i, selected_cases] = 1
                likelihood_arr.append(max(likelihood_left, likelihood_right))
                wrst_p_arr.append(max(r_left, r_right))
            elif(r_left <= perm_p_cutoff):
                genotype_cluster_df.loc[i, case_left_tail_samps] = 1
                likelihood_arr.append(likelihood_left)
                wrst_p_arr.append(r_left)
            elif(r_right <= perm_p_cutoff):
                genotype_cluster_df.loc[i, case_right_tail_samps] = 1
                likelihood_arr.append(likelihood_right)
                wrst_p_arr.append(r_right)
            else:
                r, likelihood = fit_kde_and_mannwhitney(sorted_ctrl_arr, sorted_case_arr, b)
                if(r <= perm_p_cutoff):
                    genotype_cluster_df.loc[i, selected_cases] = 1
                likelihood_arr.append(likelihood)
                wrst_p_arr.append(r)
    return np.array(likelihood_arr), np.array(wrst_p_arr), genotype_cluster_df

def MADsfromMedian(arr):
    if len(arr.shape) == 1:
        arr = arr[:,None]
    median = np.median(arr, axis=0)
    diff = np.sum((arr - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    if(med_abs_deviation == 0):
        modified_z_score = float('inf')
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score

def filter_likelihood_on_mad(likelihood_arr, mad_arr):
    lklhd_mad_dict = dict(zip(likelihood_arr, mad_arr))
    for i in sorted(lklhd_mad_dict):
        if lklhd_mad_dict[i] < 3.5:
            return i
    return float('inf')

def find_and_flag_dtu_genes(ifs, likelihood_arr, likelihood_cutoff, wrst_p_arr, perm_p_cutoff):
    wrst_pos_genes = set(ifs.iloc[np.where(wrst_p_arr <= perm_p_cutoff)[0], :].gene_id.to_list())
    likelihood_outlier_genes = set(ifs.iloc[np.where(likelihood_arr <= likelihood_cutoff)[0], :].gene_id.to_list())
    flagged_genes = likelihood_outlier_genes.intersection(wrst_pos_genes)
    return wrst_pos_genes, flagged_genes

def compute_infReps_stats(infReps, IFs, tx2gene_dict, ctrl, case, cluster_size_lim, p_cutoff, likelihood_cutoff, b, f_cpm):
    infReps_stacked = np.stack(infReps, axis=2)
    infReps_transposed = np.transpose(infReps_stacked, axes=(1, 0, 2))
    infReps_dtu_genes = []
    infReps_flagged_genes = []
    infRep_counter = 0
    for inf in infReps_transposed:
        print("--Processing inferential replicate " + str(infRep_counter+1) + " out of " + str(len(infReps_transposed)) + "--")
        col_names = [col for col in IFs.columns if col != 'gene_id']
        inf_df = pd.DataFrame(inf, columns = col_names, index = IFs.index)
        inf_df['gene_id'] = IFs.gene_id
        inf_IFs, inf_gene_counts = convert_counts_to_IF(inf_df)
        inf_likelihood_arr, inf_wrst_p_arr, inf_genotype_cluster_df = compute_double_stats(inf_IFs, inf_gene_counts, tx2gene_dict, ctrl, case, cluster_size_lim, p_cutoff, b, f_cpm)
        inf_wrst_pos_genes, inf_flagged_genes = find_and_flag_dtu_genes(inf_IFs, inf_likelihood_arr, likelihood_cutoff, inf_wrst_p_arr, p_cutoff)
        infReps_dtu_genes += inf_wrst_pos_genes
        infReps_flagged_genes += inf_flagged_genes
        infRep_counter+=1
    return infRep_counter, infReps_dtu_genes, infReps_flagged_genes

def get_stable_dtu(num_of_infReps, infReps_dtu_genes, infReps_flagged_genes, wrst_pos_genes, flagged_genes):
    infReps_dtu_genes += wrst_pos_genes
    infReps_flagged_genes += flagged_genes
    dtu_genes_counts = Counter(infReps_dtu_genes)
    dtu_genes_stable = [g for g, count in dtu_genes_counts.items() if count > int(num_of_infReps*0.5)]
    flagged_gene_counts = Counter(infReps_flagged_genes)
    flagged_genes_stable = [g for g, count in flagged_gene_counts.items() if count > int(num_of_infReps*0.5)]
    return dtu_genes_stable, flagged_genes_stable

def write_final_results(wrst_pos_genes, flagged_genes, output):
    with open(output, 'w') as f:
        f.write("gene_id\tflag\n")
        for g in wrst_pos_genes:
            if(g in flagged_genes):
                f.write(g + "\t" + "F" + "\n")
            else:
                f.write(g + "\n")
    f.close()
    return

def main(args):
    np.random.seed(42)
    IFs = pd.read_csv(args.i, sep='\t', index_col=0)
    if_cols = IFs.select_dtypes(include=[np.number]).columns
    IFs[if_cols] = IFs[if_cols].astype(np.float32).round(3)
    gene_counts = pd.read_csv(args.g, sep='\t', index_col=0)
    gene_counts = gene_counts.astype(np.int32)
    tx2gene = pd.read_csv(args.m, sep = '\t', index_col=0)
    tx2gene_dict = tx2gene.gene_id.to_dict()
    pheno = pd.read_csv(args.l, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    likelihood_arr, wrst_p_arr, genotype_cluster_df = compute_double_stats(IFs, gene_counts, tx2gene_dict, ctrl_samples, case_samples, args.n_small, args.p_cutoff, args.bandwidth, args.f_cpm, args.quiet)
    mad_scores = MADsfromMedian(likelihood_arr)
    likelihood_cutoff = filter_likelihood_on_mad(likelihood_arr, mad_scores)
    wrst_pos_genes, flagged_genes = find_and_flag_dtu_genes(IFs, likelihood_arr, likelihood_cutoff, wrst_p_arr, args.p_cutoff)
    if(args.infReps):
        all_samples = pheno.id.to_list()
        infReps = []
        for i in all_samples:
            infReps_dir = os.path.join(args.O, "SPIT_analysis", "infReps", "tximport_infReps_sample_" + i + ".txt")
            infReps_sample = pd.read_csv(infReps_dir, sep = '\t', index_col = 0)
            infReps_sample = infReps_sample.loc[IFs.index.to_list()].to_numpy()
            infReps.append(infReps_sample)
        num_of_infReps, infReps_dtu_genes, infReps_flagged_genes = compute_infReps_stats(infReps, IFs, tx2gene_dict, ctrl_samples, case_samples, args.n_small, args.p_cutoff, likelihood_cutoff, args.bandwidth, args.f_cpm)
        dtu_genes_stable, flagged_genes_stable = get_stable_dtu(num_of_infReps, infReps_dtu_genes, infReps_flagged_genes, wrst_pos_genes, flagged_genes)
        wrst_pos_genes = dtu_genes_stable
        flagged_genes = flagged_genes_stable
    if(args.exp):
        write_dir = args.exp
        write_final_results(wrst_pos_genes, flagged_genes, os.path.join(write_dir, "spit_out_k" + str(round(args.k, 1)) + ".b" + str(format(args.bandwidth, '.2f')) + ".txt"))
        genotype_cluster_df.to_csv(os.path.join(write_dir, "spit_cluster_matrix_k"+ str(round(args.k, 1)) + ".b" + str(format(args.bandwidth, '.2f')) + ".txt"), sep = '\t')
        pd.DataFrame(wrst_p_arr, index = IFs.index, columns = ['p_value']).to_csv(os.path.join(write_dir, "all_p_values.txt"), sep = '\t')
    else:
        write_dir = os.path.join(args.O, "SPIT_analysis")
        write_final_results(wrst_pos_genes, flagged_genes, os.path.join(write_dir, "spit_out.txt"))
        genotype_cluster_df.to_csv(os.path.join(write_dir, "spit_cluster_matrix.txt"), sep = '\t')
        pd.DataFrame(wrst_p_arr, index = IFs.index, columns = ['p_value']).to_csv(os.path.join(write_dir, "all_p_values.txt"), sep = '\t')

