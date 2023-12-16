import os
import sys
import argparse
import numpy as np
import pandas as pd
import warnings
import random
import math
from collections import defaultdict


def select_dtu_genes(IFs, ctrl_samples, gene_names):
    pot_dtu_genes = random.sample(gene_names, 1000)
    dtu_IFs = IFs[IFs.gene_id.isin(pot_dtu_genes)]
    dtu_IFs.insert(0, "ctrl_IF_mean", dtu_IFs[ctrl_samples].mean(axis = 1))
    ctrl_IF_max = dtu_IFs.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_no_max = dtu_IFs.drop(ctrl_IF_max.index)
    ctrl_IF_second_max = ctrl_IF_no_max.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_min = dtu_IFs.sort_values('ctrl_IF_mean', ascending=True).drop_duplicates(['gene_id'])
    dtu_gene_counter = 0
    final_dtu_genes = set()
    for g in pot_dtu_genes:
        if(dtu_gene_counter == 100):
            return list(final_dtu_genes)
        min_if = ctrl_IF_min[ctrl_IF_min.gene_id == g].ctrl_IF_mean
        min_iso = ctrl_IF_min[ctrl_IF_min.gene_id == g].index
        max_if = ctrl_IF_max[ctrl_IF_max.gene_id == g].ctrl_IF_mean
        max_iso = ctrl_IF_max[ctrl_IF_max.gene_id == g].index
        second_max_iso = ctrl_IF_second_max[ctrl_IF_second_max.gene_id == g].index
        if(np.sum(dtu_IFs.loc[max_iso, ctrl_samples].to_numpy() > dtu_IFs.loc[second_max_iso, ctrl_samples].to_numpy()) < (len(ctrl_samples) * 0.75
        )):
            pass
        else:
            final_dtu_genes.add(g)
            dtu_gene_counter += 1
            
    return list(final_dtu_genes)

def partition_samples_to_subgroups(case_samples, n_of_genotypes, n_small):
    c = 0
    val_sample_size = False
    while(1):
        if(val_sample_size):
            return genotype_sample_dict
        if(c > 50):
            print("Cannot partition samples into genotypes, make sure (number of genotypes * n_small) <= N.")
            sys.exit(1)
        genotype_sample_dict = defaultdict(list)
        val_sample_size = True
        for i in case_samples:
            genotype = random.randint(1, n_of_genotypes)
            genotype_sample_dict[genotype].append(i)
        for i in range(1, n_of_genotypes+1):
            genotype_samples = genotype_sample_dict[i]
            if(len(genotype_samples) < n_small):
                val_sample_size = False
        c+=1

def simulate_dtu(IFs, dtu_genes, ctrl_samples, case_samples, n_of_genotypes, n_small):
    genotype_gene_dict = defaultdict(list)
    genotype_cluster_df = pd.DataFrame(0, index=IFs.index, columns=IFs.columns)
    dtu_IFs = IFs[IFs.gene_id.isin(dtu_genes)]
    dtu_IFs.insert(0, "ctrl_IF_mean", dtu_IFs[ctrl_samples].mean(axis = 1))
    ctrl_IF_max = dtu_IFs.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_min = dtu_IFs.sort_values('ctrl_IF_mean', ascending=True).drop_duplicates(['gene_id'])
    for i in range(n_of_genotypes+1):
        genotype_gene_dict[i] = random.sample(dtu_genes, 30)
    genotype_sample_dict = partition_samples_to_subgroups(case_samples, n_of_genotypes, n_small)
    for i in range(1, n_of_genotypes+1):
        genotype_samples = genotype_sample_dict[i]
        for g in genotype_gene_dict[i]:
            min_if = ctrl_IF_min[ctrl_IF_min.gene_id == g].ctrl_IF_mean
            min_iso = ctrl_IF_min[ctrl_IF_min.gene_id == g].index
            max_if = ctrl_IF_max[ctrl_IF_max.gene_id == g].ctrl_IF_mean
            max_iso = ctrl_IF_max[ctrl_IF_max.gene_id == g].index
            for s in genotype_samples:
                simulated_noise = random.uniform(-0.05, 0.05)
                simulated_max = float(max_if.iloc[0]) + simulated_noise
                simulated_min = float(min_if.iloc[0]) - simulated_noise
                if(simulated_min < 0):
                    simulated_min = 0
                if(simulated_max > 1):
                    simulated_max = 1
                IFs.loc[min_iso, s] = simulated_max
                IFs.loc[max_iso, s] = simulated_min
                genotype_cluster_df.loc[min_iso, s] = 1
                genotype_cluster_df.loc[max_iso, s] = 1

    return IFs, genotype_cluster_df

def recreate_tx_counts_matrix(simulated_dtu_IFs, gene_level_counts, ids):
    merged_df = simulated_dtu_IFs.reset_index().merge(gene_level_counts.reset_index(), on='gene_id')
    tx_counts_df = pd.DataFrame()
    for i in ids:
        if_col = i + "_x"
        gene_col = i + "_y"
        tx_counts_df[str(i)] = (merged_df[if_col] * merged_df[gene_col]).astype(int)
    tx_counts_df["tx_id"] = merged_df["tx_id"]
    tx_counts_df["gene_id"] = merged_df["gene_id"]
    return tx_counts_df.set_index("tx_id", drop = True)

def convert_counts_to_IF_and_gene_level(counts):
    counts_w_genes_multilevel = counts.set_index([counts.index, 'gene_id'])
    gene_level_counts = counts.groupby('gene_id').sum()
    gene_level_counts = gene_level_counts + 0.00001
    IFs = counts_w_genes_multilevel.div(gene_level_counts,axis='index',level='gene_id').reset_index('gene_id')
    
    return IFs, gene_level_counts


def main(args):
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    IFs = pd.read_csv(args.i, sep='\t', index_col=0)
    gene_level_counts = pd.read_csv(args.g, sep='\t', index_col=0)
    gene_names = list(IFs.gene_id.unique())
    ctrl_file = os.path.join(args.exp, 'ctrl_samples.txt')
    case_file = os.path.join(args.exp, 'case_samples.txt')
    with open(ctrl_file) as ctrl_file:
        ctrl_samples = [line.strip() for line in ctrl_file]
    with open(case_file) as case_file:
        case_samples = [line.strip() for line in case_file]
    all_ids = ctrl_samples + case_samples
    sim_pheno_condition = [0 for i in range(len(ctrl_samples))] + [1 for i in range(len(case_samples))]
    sim_pheno = (pd.DataFrame([all_ids, sim_pheno_condition])).transpose()
    sim_pheno.columns = ['id', 'condition']
    sim_pheno.to_csv(os.path.join(args.exp, 'simulation_pheno.txt'), sep = '\t', index = False)
    dtu_genes = select_dtu_genes(IFs, ctrl_samples, gene_names)
    simulated_dtu_ifs, genotype_cluster_df = simulate_dtu(IFs, dtu_genes, ctrl_samples, case_samples, args.n_splicotypes, args.n_small)
    simulated_dtu_counts = recreate_tx_counts_matrix(simulated_dtu_ifs, gene_level_counts, all_ids)
    corrected_IFs, corrected_gene_level_counts = convert_counts_to_IF_and_gene_level(simulated_dtu_counts)
    corrected_IFs.to_csv(os.path.join(args.exp, 'simulated_ifs.txt'), sep = '\t')
    simulated_dtu_counts.to_csv(os.path.join(args.exp, 'simulated_tx_counts.txt'), sep = '\t')
    corrected_gene_level_counts.to_csv(os.path.join(args.exp, 'simulated_gene_counts.txt'), sep = '\t')
    genotype_cluster_df.drop(columns=['gene_id']).to_csv(os.path.join(args.exp, 'true_cluster_matrix.txt'), sep = '\t')
