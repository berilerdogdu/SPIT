import sys, getopt
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

def filter_samples_on_gene_cpm(f_par, tx, tx2gene_dict, samples, ifs, cpms):
    gene_id = tx2gene_dict[tx]
    cpm_gene_lvl = cpms.loc[gene_id, samples].to_dict()
    selected_samples = []
    filtered_samples = []
    if(f_par):
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

def index_if_arrays(arr, b):
    sort_ind_arr = np.argsort(arr)
    sorted_arr = arr[sort_ind_arr]
    split_ind = split_vector_on_kernel(sorted_arr, b)
    return sort_ind_arr, sorted_arr, split_ind
    
def split_vector_on_kernel(sorted_vector, b):
    p = np.linspace(0,1).reshape(-1, 1)
    kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(sorted_vector.reshape(-1, 1))
    scores = kde.score_samples(p)
    mi = argrelextrema(scores, np.less)[0]
    split_ind = len(sorted_vector) - 1
    split_val = -1
    if(len(mi) > 0):
        split_val = p[mi[np.argmin(scores[mi])]]
        split_ind = np.nonzero(sorted_vector <= split_val)[0][-1]
    return split_ind
        
def subgroupsize_limit_check(group, split_in, size_lim):
    if((split_in >= size_lim-1) and ((len(group)-split_in) > size_lim)):
        return 1
    return 0
 
def fit_kde_and_mannwhitney(ctrl_subgroup, case_subgroup, b):
    ctrl_kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(ctrl_subgroup.reshape(-1, 1))
    likelihood = ctrl_kde.score(case_subgroup.reshape(-1, 1), y = None)
    r = sts.mannwhitneyu(case_subgroup, ctrl_subgroup, method = 'auto', alternative='two-sided')
    return r, likelihood

def compute_double_stats(ifs, gene_counts, tx2gene_dict, ctrl, case, cluster_size_lim, perm_p_cutoff, b, f_par):
    likelihood_arr = []
    wrst_p_arr = []
    cpms = convert_counts_to_cpm(gene_counts)
    genotype_cluster_df = pd.DataFrame(0, index=ifs.index, columns=case)
    for i in ifs.index.to_list():
        selected_ctrls, ctrl_arr, filtered_ctrls = filter_samples_on_gene_cpm(f_par, i, tx2gene_dict, ctrl, ifs, cpms)
        selected_cases, case_arr, filtered_cases = filter_samples_on_gene_cpm(f_par, i, tx2gene_dict, case, ifs, cpms)
        genotype_cluster_df.loc[i, filtered_cases] = None
        if((groupsize_limit_check(ctrl_arr, case_arr, cluster_size_lim)) == 0):
            wrst_p_arr.append(1)
            continue
        sort_ind_ctrl_arr, sorted_ctrl_arr, ctrl_split_ind = index_if_arrays(ctrl_arr, b)
        sort_ind_case_arr, sorted_case_arr, case_split_ind = index_if_arrays(case_arr, b)
        case_tx_likelihood = 0
        
        if((subgroupsize_limit_check(sorted_case_arr, case_split_ind, cluster_size_lim)) == 0):
            r, likelihood = fit_kde_and_mannwhitney(sorted_ctrl_arr, sorted_case_arr, b)
            if(r[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i] = 1
            wrst_p_arr.append(r[1])
        else:
            case_left_tail = sorted_case_arr[0:case_split_ind]
            case_left_tail_samps = selected_cases[sort_ind_case_arr[0:case_split_ind]]
            case_right_tail = sorted_case_arr[case_split_ind:]
            case_right_tail_samps = selected_cases[sort_ind_case_arr[case_split_ind:]]
            if((subgroupsize_limit_check(sorted_ctrl_arr, ctrl_split_ind, cluster_size_lim))):
                ctrl_left_tail = sorted_ctrl_arr[0:ctrl_split_ind]
                ctrl_right_tail = sorted_ctrl_arr[ctrl_split_ind:]
            elif((len(sorted_ctrl_arr) >= len(case_left_tail)) and (len(sorted_ctrl_arr) >= len(case_left_tail))):
                ctrl_left_tail = sorted_ctrl_arr[0:case_split_ind]
                ctrl_right_tail = sorted_ctrl_arr[-len(sorted_case_arr[case_split_ind:]):]
            else:
                ctrl_left_tail = sorted_ctrl_arr
                ctrl_right_tail = sorted_ctrl_arr
            r_left, likelihood_left = fit_kde_and_mannwhitney(ctrl_left_tail, case_left_tail, b)
            r_right, likelihood_right = fit_kde_and_mannwhitney(ctrl_right_tail, case_right_tail, b)
            likelihood_arr.append(min(likelihood_left, likelihood_right))
            if(r_left[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i, case_left_tail_samps] = 1
            if(r_right[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i, case_right_tail_samps] = 1
            wrst_p_arr.append(min(r_left[1], r_right[1]))
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

    IFs_file = ''
    gene_counts_file = ''
    labels_file = ''
    tx2gene_file = ''
    perm_p_cutoff = ''
    output_file = ''
    cluster_size_limit = 12
    b = ''
    f_par = True
    spit_cluster_matrix = ''   
    try:
        opts, args = getopt.getopt(argv,"hi:g:l:m:p:b:n:A:O:",["IFs_file=", "gene_counts_file=", "labels_file=", "tx2gene_file=", "perm_p_cutoff=", "b=", "n_small=", "spit_cluster_matrix=", "output_file="])
    except getopt.GetoptError:
        print("Usage: dtu_detection.py -i <input isoform fractions file> -g <input gene counts file> -l <input label file> -m <transcript to gene ids mapping file> -p <permutation p cutoff> -b <KDE bandwidth> -n <n_small> -A <cluster array file> -O <output file path>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: dtu_detection.py -i <input isoform fractions file> -g <input gene counts file> -l <input label file> -m <transcript to gene ids mapping file> -p <permutation p cutoff> -b <KDE bandwidth> -n <n_small> -A <cluster array file> -O <output file path>")
            sys.exit()
        elif opt in ("-i", "--IFs_file"):
            IFs_file = arg
        elif opt in ("-g", "--gene_counts_file"):
            gene_counts_file = arg
        elif opt in ("-l", "--labels_file"):
            labels_file = arg
        elif opt in ("-m", "--tx2gene_file"):
            tx2gene_file = arg
        elif opt in ("-p", "--perm_p_cutoff"):
            perm_p_cutoff = float(arg)
        elif opt in ("-b", "--b"):
            b = float(arg)
        elif opt in ("-n", "--n_small"):
            cluster_size_limit = int(arg)
        elif opt in ("-A", "--spit_cluster_matrix"):
            spit_cluster_matrix = arg
        elif opt in ("-O", "--output_file"):
            output_file = arg           
    IFs = pd.read_csv(IFs_file, sep='\t', index_col=0)
    tx_names = list(IFs.index)
    IFs = IFs.round(2)   
    gene_counts = pd.read_csv(gene_counts_file, sep='\t', index_col=0)
    tx2gene = pd.read_csv(tx2gene_file, sep = '\t', index_col=0)
    tx2gene_dict = tx2gene.gene_id.to_dict()   
    pheno = pd.read_csv(labels_file, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    likelihood_arr, wrst_p_arr, genotype_cluster_df = compute_double_stats(IFs, gene_counts, tx2gene_dict, ctrl_samples, case_samples, cluster_size_limit, perm_p_cutoff, b, f_par)
    mad_scores = MADsfromMedian(likelihood_arr)
    likelihood_cutoff = filter_likelihood_on_dmad(likelihood_arr, mad_scores)
    wrst_pos_genes, flagged_genes = find_and_flag_dtu_genes(IFs, likelihood_arr, likelihood_cutoff, wrst_p_arr, perm_p_cutoff)
    write_final_results(wrst_pos_genes, flagged_genes, output_file)
    genotype_cluster_df.to_csv(spit_cluster_matrix, sep = '\t')
    
    
if __name__ == "__main__":
   main(sys.argv[1:])
