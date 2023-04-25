import sys, getopt
import numpy as np
import pandas as pd
import random
import math
import scipy.stats as sts
from collections import defaultdict


def select_genes_w_dom_iso(IFs, ctrl_samples, gene_names, p_dom):

    IFs.insert(0, "ctrl_IF_mean", IFs[ctrl_samples].mean(axis = 1))
    ctrl_IF_max = IFs.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_no_max = IFs.drop(ctrl_IF_max.index)
    ctrl_IF_second_max = ctrl_IF_no_max.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_min = IFs.sort_values('ctrl_IF_mean', ascending=True).drop_duplicates(['gene_id'])
    genes_w_dom_iso_counter = 0
    final_dom_iso_genes = set()
    
    for g in gene_names:
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
    

def mannwhitneyu_permutation(IFs, ctrl_samples, num_of_it):

    min_perm_p_arr = []
    perm_p_arr = []
    sampled_txs = set()

    for it in range(1, num_of_it+1):
        out = '\r'+ str(it) + " iterations completed."
        if (it == num_of_it):
            print(out)
        else:
            print(out, end='')
        s = int(len(ctrl_samples)/2)
        tx_matrix = np.transpose(IFs[ctrl_samples].to_numpy())
        np.random.shuffle(tx_matrix)
        m_1 = tx_matrix[0:s,]
        m_2 = tx_matrix[s:,]
        mwu_results = sts.mannwhitneyu(m_1, m_2, method = 'auto', alternative='two-sided', axis = 0)
        perm_p_arr.append(mwu_results[1])
        min_tx_ind = np.argmin(mwu_results[1])
        i = 0
        while(min_tx_ind in sampled_txs):
            i+=1
            min_tx_ind = np.argsort(mwu_results[1])[i]
        sampled_txs.add(min_tx_ind)
        min_it_p_val = mwu_results[1][min_tx_ind]
        min_perm_p_arr.append(min_it_p_val)
    
    return min_perm_p_arr








def determine_perm_p_cutoff(min_perm_p_arr, p_values_file):

    p_cutoff = np.array(min_perm_p_arr).min()
    with open(p_values_file, 'w') as f:
        for i in min_perm_p_arr:
            f.write("%s\n" % i)
    return p_cutoff



def main(argv):

    IFs_file = ''
    gene_counts_file = ''
    labels_file = ''
    num_of_it = ''
    dominance_filtered_ifs_file = ''
    dominance_filtered_gene_counts_file = ''
    p_dom = 0.75
    p_values_file = ''
    
    try:
        opts, args = getopt.getopt(argv,"hi:g:m:l:n:I:G:P:d:",["IFs_file=", "gene_counts_file=", "labels_file=", "num_of_it=", "dominance_filtered_ifs_file=", "dominance_filtered_gene_counts_file=", "p_values_file=", "dominance_filter_threshold"])
    except getopt.GetoptError:
        print("Usage: spit_test.py -i <input isoform fractions file> -g <input gene counts file> -l <input pheno file> -n <number of iterations> -I <dominance filtered ifs file> -G <dominance filtered gene counts file> -P <p values file> -d dominance_filter_threshold")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: spit_test.py -i <input isoform fractions file> -g <input gene counts file> -l <input pheno file> -n <number of iterations> -I <dominance filtered ifs file> -G <dominance filtered gene counts file> -P <p values file> -d dominance_filter_threshold")
            sys.exit()
        elif opt in ("-i", "--IFs_file"):
            IFs_file = arg
        elif opt in ("-g", "--gene_counts_file"):
            gene_counts_file = arg
        elif opt in ("-l", "--labels_file"):
            labels_file = arg
        elif opt in ("-n", "--num_of_it"):
            num_of_it = int(arg)
        elif opt in ("-I", "--dominance_filtered_ifs_file"):
            dominance_filtered_ifs_file = arg
        elif opt in ("-G", "--dominance_filtered_gene_counts_file"):
            dominance_filtered_gene_counts_file = arg
        elif opt in ("-P", "--p_values_file"):
            p_values_file = arg
        elif opt in ("-d", "--dominance_filter_threshold"):
            p_dom = int(arg)
            
    IFs = pd.read_csv(IFs_file, sep='\t', index_col=0)
    tx_names = list(IFs.index)
    gene_names = list(IFs.gene_id.unique())
    IFs = IFs.round(2)
    
    gene_counts = pd.read_csv(gene_counts_file, sep='\t', index_col=0)
    pheno = pd.read_csv(labels_file, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    
    selected_dom_iso_genes = select_genes_w_dom_iso(IFs, ctrl_samples, gene_names, p_dom)
    IFs_selected = IFs[IFs.gene_id.isin(selected_dom_iso_genes)]
    IFs_selected.to_csv(dominance_filtered_ifs_file, sep = '\t')
    gene_counts_selected = gene_counts[gene_counts.index.isin(selected_dom_iso_genes)]
    gene_counts_selected.to_csv(dominance_filtered_gene_counts_file, sep = '\t')
    min_perm_p_arr = mannwhitneyu_permutation(IFs_selected, ctrl_samples, num_of_it)
    p_cutoff = determine_perm_p_cutoff(min_perm_p_arr, p_values_file)
    

if __name__ == "__main__":
   main(sys.argv[1:])
