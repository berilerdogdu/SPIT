import sys
import argparse
import random
import math
import pandas as pd
import numpy as np
import scipy.stats as sts
from collections import defaultdict
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
from matplotlib import pyplot as plt
from tqdm import tqdm


def convert_counts_to_cpm(counts):
    lib_sizes = counts.sum(axis = 0).to_list()
    scaled_counts = counts.divide(lib_sizes)
    cpms = scaled_counts * (10**6)
    return cpms

def filter_samples_on_gene_cpm(f_cpm, tx, tx2gene_dict, samples, ifs, cpms):
    gene_id = tx2gene_dict[tx]
    cpm_gene_lvl = cpms.loc[gene_id, samples].to_dict()
    selected_samples = []
    filtered_samples = []
    if(f_cpm):
        for s in cpm_gene_lvl:
            if (cpm_gene_lvl[s] > 10):
                selected_samples.append(s)
            else:
                filtered_samples.append(s)
    else:
        selected_samples = samples
    selected_if_arr = ifs.loc[tx, selected_samples].to_numpy(copy = True, dtype=np.float64)
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
     
def fit_kde_and_mannwhitney(ctrl_subgroup, case_subgroup, b):
    ctrl_kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(ctrl_subgroup.reshape(-1, 1))
    likelihood = ctrl_kde.score(case_subgroup.reshape(-1, 1), y = None)
    r = sts.mannwhitneyu(case_subgroup, ctrl_subgroup, method = 'auto', alternative='two-sided')
    return r, likelihood

def compute_double_stats(ifs, gene_counts, tx2gene_dict, ctrl, case, cluster_size_lim, perm_p_cutoff, b, f_cpm):
    likelihood_arr = []
    wrst_p_arr = []
    cpms = convert_counts_to_cpm(gene_counts)
    genotype_cluster_df = pd.DataFrame(0, index=ifs.index, columns=case)
    for i in tqdm(ifs.index.to_list()):
        selected_ctrls, ctrl_arr, filtered_ctrls = filter_samples_on_gene_cpm(f_cpm, i, tx2gene_dict, ctrl, ifs, cpms)
        selected_cases, case_arr, filtered_cases = filter_samples_on_gene_cpm(f_cpm, i, tx2gene_dict, case, ifs, cpms)
        genotype_cluster_df.loc[i, filtered_cases] = None
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
            if(r[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i, selected_cases] = 1
            likelihood_arr.append(likelihood)
            wrst_p_arr.append(r[1])
        else:
            case_left_tail = sorted_case_arr[0:case_split_ind]
            case_left_tail_samps = selected_cases[sort_ind_case_arr[0:case_split_ind]]
            case_right_tail = sorted_case_arr[case_split_ind:]
            case_right_tail_samps = selected_cases[sort_ind_case_arr[case_split_ind:]]
            ctrl_left_tail, ctrl_right_tail = set_ctrl_tails(sorted_ctrl_arr, ctrl_split_ind, cluster_size_lim)
            r_left, likelihood_left = fit_kde_and_mannwhitney(ctrl_left_tail, case_left_tail, b)
            r_right, likelihood_right = fit_kde_and_mannwhitney(ctrl_right_tail, case_right_tail, b)
            if((r_left[1] < perm_p_cutoff) and (r_right[1] < perm_p_cutoff)):
                genotype_cluster_df.loc[i, selected_cases] = 1
                likelihood_arr.append(max(likelihood_left, likelihood_right))
                wrst_p_arr.append(max(r_left[1], r_right[1]))
            elif(r_left[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i, case_left_tail_samps] = 1
                likelihood_arr.append(likelihood_left)
                wrst_p_arr.append(r_left[1])
            elif(r_right[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i, case_right_tail_samps] = 1
                likelihood_arr.append(likelihood_right)
                wrst_p_arr.append(r_right[1])
            else:
                r, likelihood = fit_kde_and_mannwhitney(sorted_ctrl_arr, sorted_case_arr, b)
                if(r[1] < perm_p_cutoff):
                    genotype_cluster_df.loc[i, selected_cases] = 1
                likelihood_arr.append(likelihood)
                wrst_p_arr.append(r[1])
    return np.array(likelihood_arr), np.array(wrst_p_arr), genotype_cluster_df

def MADsfromMedian(arr):
    if len(arr.shape) == 1:
        arr = arr[:,None]
    median = np.median(arr, axis=0)
    diff = np.sum((arr - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score

def filter_likelihood_on_dmad(likelihood_arr, dmad_arr):
    lklhd_dmad_dict = dict(zip(likelihood_arr, dmad_arr))
    for i in sorted(lklhd_dmad_dict):
        if lklhd_dmad_dict[i] < 3.5:
            return i
    return

def find_and_flag_dtu_genes(ifs, likelihood_arr, likelihood_cutoff, wrst_p_arr, perm_p_cutoff):
    pd.DataFrame(wrst_p_arr, index = ifs.index, columns = ['p_value']).to_csv("all_p_values.txt", sep = '\t')
    wrst_pos_genes = set(ifs.iloc[np.where(wrst_p_arr < perm_p_cutoff)[0], :].gene_id.to_list())
    likelihood_outlier_genes = set(ifs.iloc[np.where(likelihood_arr < likelihood_cutoff)[0], :].gene_id.to_list())
    flagged_genes = likelihood_outlier_genes.intersection(wrst_pos_genes)
    return wrst_pos_genes, flagged_genes

def write_final_results(wrst_pos_genes, flagged_genes, output):
    with open(output, 'w') as f:
        for g in wrst_pos_genes:
            if(g in flagged_genes):
                f.write(g + "\t" + "F" + "\n")
            else:
                f.write(g + "\n")
    f.close()
    return

def main(argv):

    parser = argparse.ArgumentParser()
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-i', metavar='dom_selected_ifs.txt', required=True, type=str, help='Isoform fractions file (tsv)')
    parser.add_argument('-g', metavar='dom_selected_gene_counts.txt', required=True, type=str, help='Gene counts file (tsv)')
    parser.add_argument('-m', metavar='tx2gene.txt', required=True, type=str, help='Transcript to gene mapping file (tsv)')
    parser.add_argument('-l', metavar='labels.txt', required=True, type=str, help='Labels/metadata file (tsv)')
    parser.add_argument('--p_cutoff', metavar='0.05', required=True, type=float, help='p-value cutoff from SPIT Test permutations')
    parser.add_argument('-b', '--bandwidth', metavar='0.09', type=float, default=0.09, help='p-value cutoff from SPIT Test permutations')
    parser.add_argument('-n', '--n_small', metavar='12', type=int, default=12, help='Smallest sample size for the subgroups')
    parser.add_argument('--f_cpm', action='store_true', help='Apply filter-CPM thresholding')
    parser.add_argument('-A', metavar='spit_cluster_array.txt', required=True, type=str, help='Output file path for SPIT cluster array')
    parser.add_argument('-O', metavar='spit_out.txt', required=True, type=str, help='Output file path for candidate DTU genes')
    
    args = parser.parse_args()
    IFs = pd.read_csv(args.i, sep='\t', index_col=0)
    tx_names = list(IFs.index)
    IFs = IFs.round(2)
    gene_counts = pd.read_csv(args.g, sep='\t', index_col=0)
    tx2gene = pd.read_csv(args.m, sep = '\t', index_col=0)
    tx2gene_dict = tx2gene.gene_id.to_dict()
    pheno = pd.read_csv(args.l, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    likelihood_arr, wrst_p_arr, genotype_cluster_df = compute_double_stats(IFs, gene_counts, tx2gene_dict, ctrl_samples, case_samples, args.n_small, args.p_cutoff, args.bandwidth, args.f_cpm)
    double_mad_scores = MADsfromMedian(likelihood_arr)
    likelihood_cutoff = filter_likelihood_on_dmad(likelihood_arr, double_mad_scores)
    wrst_pos_genes, flagged_genes = find_and_flag_dtu_genes(IFs, likelihood_arr, likelihood_cutoff, wrst_p_arr, args.p_cutoff)
    write_final_results(wrst_pos_genes, flagged_genes, args.O)
    genotype_cluster_df.to_csv(args.A, sep = '\t')


if __name__ == "__main__":
   main(sys.argv[1:])
