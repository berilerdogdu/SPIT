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


def convert_counts_to_cpm(counts):
    
    lib_sizes = counts.sum(axis = 0).to_list()
    scaled_counts = counts.divide(lib_sizes)
    cpms = scaled_counts * (10**6)
    
    return cpms


def filter_samples_on_gene_cpm(tx, tx2gene_dict, ifs, gene_counts, ctrl, case, cluster_size_lim):
    
    cpm_gene_lvl = convert_counts_to_cpm(gene_counts)
    
    gene_id = tx2gene_dict[tx]
    cpm_filtered_ifs = ifs.loc[tx, cpm_gene_lvl.columns[(cpm_gene_lvl.loc[gene_id] >= 10)]]
    cpm_filtered_ctrl = cpm_filtered_ifs[cpm_filtered_ifs.index.isin(ctrl)]
    cpm_filtered_case = cpm_filtered_ifs[cpm_filtered_ifs.index.isin(case)]
    ctrl_arr = cpm_filtered_ctrl.to_numpy(copy = True, dtype=np.float64)
    case_arr = cpm_filtered_case.to_numpy(copy = True, dtype=np.float64)

    return ctrl_arr, case_arr


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
    
  
def compute_double_stats(ifs, gene_counts, tx2gene_dict, ctrl, case, cluster_size_lim, perm_p_cutoff, b):
    likelihood_arr = []
    cond_arr = []
    wrst_p_arr = []
    genotype_cluster_df = pd.DataFrame(0, index=ifs.index, columns=[c for c in range(len(case))])

    for i in ifs.index.to_list():
        ctrl_arr, case_arr = filter_samples_on_gene_cpm(i, tx2gene_dict, ifs, gene_counts, ctrl, case, cluster_size_lim)
        if(len(ctrl_arr) < cluster_size_lim):
            wrst_p_arr.append(1)
            cond_arr.append(-1)
            continue
        if(len(case_arr) < cluster_size_lim):
            wrst_p_arr.append(1)
            cond_arr.append(-1)
            continue
        
        sort_ind_ctrl_arr = np.argsort(ctrl_arr)
        sort_ind_case_arr = np.argsort(case_arr)
        
        sorted_ctrl_arr = ctrl_arr[sort_ind_ctrl_arr]
        sorted_case_arr = case_arr[sort_ind_case_arr]
        
        ctrl_split_ind = split_vector_on_kernel(sorted_ctrl_arr, b)
        case_split_ind = split_vector_on_kernel(sorted_case_arr, b)
        
        case_tx_likelihood = 0

        if((case_split_ind >= cluster_size_lim-1) and ((len(sorted_case_arr) - case_split_ind) > cluster_size_lim)):
            if((ctrl_split_ind >= cluster_size_lim-1) and ((len(sorted_ctrl_arr) - ctrl_split_ind) > cluster_size_lim)):
            
                ctrl_kde_1 = KernelDensity(kernel='gaussian', bandwidth=b).fit(sorted_ctrl_arr[0:ctrl_split_ind].reshape(-1, 1))
                likelihood_1 = ctrl_kde_1.score(sorted_case_arr[0:case_split_ind].reshape(-1, 1), y = None)
                r_1 = sts.mannwhitneyu(sorted_case_arr[0:case_split_ind], sorted_ctrl_arr[0:ctrl_split_ind], method = 'auto', alternative='two-sided')

                ctrl_kde_2 = KernelDensity(kernel='gaussian', bandwidth=b).fit(sorted_ctrl_arr[ctrl_split_ind:].reshape(-1, 1))
                likelihood_2 = ctrl_kde_2.score(sorted_case_arr[case_split_ind:].reshape(-1, 1), y = None)
                r_2 = sts.mannwhitneyu(sorted_case_arr[case_split_ind:], sorted_ctrl_arr[ctrl_split_ind:], method = 'auto', alternative='two-sided')
                
                case_tx_likelihood = min(likelihood_1, likelihood_2)
                likelihood_arr.append(case_tx_likelihood)
                if(r_1[1] < perm_p_cutoff):
                    genotype_cluster_df.loc[i, sort_ind_case_arr[0:case_split_ind]] = 1

                if(r_2[1] < perm_p_cutoff):
                    genotype_cluster_df.loc[i, sort_ind_case_arr[case_split_ind:]] = 1
                wrst_p = min(r_1[1], r_2[1])
                wrst_p_arr.append(wrst_p)
                cond_arr.append(1)
                
            else:
                ctrl_split_ind = case_split_ind
                if((len(sorted_ctrl_arr) - case_split_ind) <= cluster_size_lim):
                    ctrl_split_ind = len(sorted_ctrl_arr) - 1
                ctrl_kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(sorted_ctrl_arr.reshape(-1, 1))
                likelihood_1 = ctrl_kde.score(sorted_case_arr[0:case_split_ind].reshape(-1, 1), y = None)
                likelihood_2 = ctrl_kde.score(sorted_case_arr[case_split_ind:].reshape(-1, 1), y = None)
                r_1 = sts.mannwhitneyu(sorted_case_arr[0:case_split_ind], sorted_ctrl_arr[0:ctrl_split_ind], method = 'auto', alternative='two-sided')
                r_2 = sts.mannwhitneyu(sorted_case_arr[case_split_ind:], sorted_ctrl_arr[ctrl_split_ind:], method = 'auto', alternative='two-sided')

                case_tx_likelihood = min(likelihood_1, likelihood_2)
                likelihood_arr.append(case_tx_likelihood)

                if(r_1[1] < perm_p_cutoff):
                    genotype_cluster_df.loc[i, sort_ind_case_arr[0:case_split_ind]] = 1
                if(r_2[1] < perm_p_cutoff):
                    genotype_cluster_df.loc[i, sort_ind_case_arr[case_split_ind:]] = 1

                wrst_p = min(r_1[1], r_2[1])
                wrst_p_arr.append(wrst_p)
                cond_arr.append(2)
        else:
            ctrl_kde = KernelDensity(kernel='gaussian', bandwidth=b).fit(sorted_ctrl_arr.reshape(-1, 1))
            case_tx_likelihood = ctrl_kde.score(sorted_case_arr.reshape(-1, 1), y = None)
            likelihood_arr.append(case_tx_likelihood)
            r = sts.mannwhitneyu(sorted_case_arr, sorted_ctrl_arr, method = 'auto', alternative='two-sided')
            if(r[1] < perm_p_cutoff):
                genotype_cluster_df.loc[i] = 1

            wrst_p_arr.append(r[1])
            cond_arr.append(3)
    genotype_cluster_df.columns = case
    return np.array(likelihood_arr), np.array(wrst_p_arr), np.array(cond_arr), genotype_cluster_df


def doubleMADsfromMedian(arr):
    m = np.median(arr)
    abs_dev = np.abs(arr - m)
    left_mad = np.median(abs_dev[arr <= m])
    right_mad = np.median(abs_dev[arr >= m])
    arr_mad = left_mad * np.ones(len(arr))
    arr_mad[arr > m] = right_mad
    modified_z_score = abs_dev / arr_mad
    modified_z_score[arr == m] = 0
    return modified_z_score


def filter_likelihood_on_dmad(likelihood_arr, dmad_arr):
    lklhd_dmad_dict = dict(zip(likelihood_arr, dmad_arr))
    for i in sorted(lklhd_dmad_dict):
        if lklhd_dmad_dict[i] < 3.5:
            return i
    return


def find_and_flag_dtu_genes(ifs, likelihood_arr, likelihood_cutoff, wrst_p_arr, perm_p_cutoff, cond_arr):

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
    b = ''
    spit_cluster_matrix = ''
    
    
    try:
        opts, args = getopt.getopt(argv,"hi:g:l:m:p:b:A:O:",["IFs_file=", "gene_counts_file=", "labels_file=", "tx2gene_file=", "perm_p_cutoff=", "b=", "spit_cluster_matrix=", "output_file="])
    except getopt.GetoptError:
        print("Usage: dtu_detection.py -i <input isoform fractions file> -g <input gene counts file> -l <input label file> -m <transcript to gene ids mapping file> -p <permutation p cutoff> -b <KDE bandwidth> -A <cluster array file> -O <output file path>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: dtu_detection.py -i <input isoform fractions file> -g <input gene counts file> -l <input label file> -m <transcript to gene ids mapping file> -p <permutation p cutoff> -b <KDE bandwidth> -A <cluster array file> -O <output file path>")
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
        elif opt in ("-d", "--case_samples"):
            case_samples_file = arg
        elif opt in ("-c", "--ctrl_samples"):
            ctrl_samples_file = arg
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
    likelihood_arr, wrst_p_arr, cond_arr, genotype_cluster_df = compute_double_stats(IFs, gene_counts, tx2gene_dict, ctrl_samples, case_samples, 12, perm_p_cutoff, b)
    double_mad_scores = doubleMADsfromMedian(likelihood_arr)
    likelihood_cutoff = filter_likelihood_on_dmad(likelihood_arr, double_mad_scores)
    wrst_pos_genes, flagged_genes = find_and_flag_dtu_genes(IFs, likelihood_arr, likelihood_cutoff, wrst_p_arr, perm_p_cutoff, cond_arr)
    write_final_results(wrst_pos_genes, flagged_genes, output_file)
    genotype_cluster_df.to_csv(spit_cluster_matrix, sep = '\t')


if __name__ == "__main__":
   main(sys.argv[1:])
