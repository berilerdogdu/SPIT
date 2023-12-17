#!/usr/bin/env python3

"""
Description:
    This module performs pre-filtering of the input transcripts. The protocol implemented in this modules follows closely the filtering criteria defined in DRIMSeq (Nowicka M. et al, 2016) and involves filtering by Counts Per Million (CPM), Sample Count, Isoform Fraction (IF), Transcrit Count Per Gene.
Usage:
    ./filter_and_transform_tx_counts.py -i <tx_counts> -m <tx2gene> -l <pheno> -T <filtered_tx_counts> -F <filtered_ifs> -G <filtered_gene_counts> --write
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
from collections import defaultdict
from tqdm import tqdm


def convert_counts_to_cpm(counts):
    lib_sizes = counts.sum(axis = 0).to_list()
    scaled_counts = counts.divide(lib_sizes)
    cpms = scaled_counts * (10**6)
    return cpms

def filter_on_cpm(cpms, n):
    filter_count = (cpms > 1).sum(axis = 1)
    filtered_cpm_row_ind = cpms[filter_count >= n].index.to_list()
    return filtered_cpm_row_ind

def filter_on_tx_count(counts, pr_fraction, ctrl, case):
    filtered_tx_row_ind = counts.index.to_list()
    ctrl_zero_txs = counts[((counts[ctrl] == 0).astype(int).sum(axis=1) >= (len(ctrl)*pr_fraction))].index.to_list()
    case_zero_txs = counts[((counts[case] == 0).astype(int).sum(axis=1) >= (len(case)*pr_fraction))].index.to_list()
    both_zero_txs = set(ctrl_zero_txs).intersection(set(case_zero_txs))
    filtered_tx_row_ind = list(set(filtered_tx_row_ind).difference(both_zero_txs))
    return filtered_tx_row_ind

def convert_counts_to_IF_and_gene_level(counts, tx2gene, genefilter_count, genefilter_sample):
    counts_w_genes = counts.join(tx2gene.set_index('tx_id'))
    gene_level_counts = counts_w_genes.groupby('gene_id').sum()
    gene_level_counts = gene_level_counts + 0.00001
    gene_threshold_bool = gene_level_counts[gene_level_counts < genefilter_count].count(axis = 1) < genefilter_sample
    valid_genes = gene_threshold_bool[gene_threshold_bool].index
    gene_level_counts = gene_level_counts[gene_level_counts.index.isin(valid_genes)]
    counts_w_genes = counts_w_genes[counts_w_genes.gene_id.isin(valid_genes)]
    counts_w_genes_multilevel = counts_w_genes.set_index([counts_w_genes.index, 'gene_id'])
    IFs = counts_w_genes_multilevel.div(gene_level_counts,axis='index',level='gene_id').droplevel('gene_id')
    return IFs, gene_level_counts
    
def filter_on_IF(IFs, n, f):
    filter_count = (IFs > f).sum(axis = 1)
    filtered_IFs_row_ind = IFs[filter_count >= n].index.to_list()
    return filtered_IFs_row_ind

def filter_on_isoform_count(counts, tx2gene):
    counts_w_genes = counts.join(tx2gene.set_index('tx_id'))
    gene_isoform_counts = counts_w_genes.groupby('gene_id').count().iloc[:, 1]
    filtered_genes = gene_isoform_counts[gene_isoform_counts > 1].index.to_list()
    return counts_w_genes[counts_w_genes.gene_id.isin(filtered_genes)].index.to_list()
    
def select_genes_w_dom_iso(IFs, ctrl_samples, p_dom):
    IFs = IFs.round(2)
    IFs.insert(0, "ctrl_IF_mean", IFs[ctrl_samples].mean(axis = 1))
    ctrl_IF_max = IFs.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_no_max = IFs.drop(ctrl_IF_max.index)
    ctrl_IF_second_max = ctrl_IF_no_max.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_min = IFs.sort_values('ctrl_IF_mean', ascending=True).drop_duplicates(['gene_id'])
    final_dom_iso_genes = set()
    print("Selecting genes with consistently dominant isoforms in control group:")
    for g in tqdm(list(IFs.gene_id.unique())):
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
    return list(final_dom_iso_genes)


def main(args):
    tx_count_data = pd.read_csv(args.i, sep='\t', index_col=0)
    tx2gene = pd.read_csv(args.m, sep = '\t')
    txs_w_genes = tx_count_data.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("Input number of transcripts: ", tx_count_data.shape[0])
        print("Input number of genes: ", len(txs_w_genes.gene_id.unique()))
    tx_count_data = tx_count_data[tx_count_data.values.sum(axis=1) != 0]
    tx_count_data_w_genes = tx_count_data.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("After filtering out all-zeros:")
        print("\tNumber of transcripts", len(tx_count_data.index))
        print("\tNumber of genes: ", len(tx_count_data_w_genes.gene_id.unique()))
    pheno = pd.read_csv(args.l, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    cpms = convert_counts_to_cpm(tx_count_data)
    filtered_inds_on_cpm = filter_on_cpm(cpms, args.n_small)
    filtered_count_data_on_cpm = tx_count_data.loc[filtered_inds_on_cpm,:]
    on_cpm_w_genes = filtered_count_data_on_cpm.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("After CPM-filter:")
        print("\tNumber of transcripts", len(filtered_count_data_on_cpm.index))
        print("\tNumber of genes: ", len(on_cpm_w_genes.gene_id.unique()))
    filtered_inds_on_tx_and_cpm = filter_on_tx_count(filtered_count_data_on_cpm, args.pr_fraction, ctrl_samples, case_samples)
    filtered_counts_data_on_tx_and_cpm = filtered_count_data_on_cpm.loc[filtered_inds_on_tx_and_cpm,:]
    on_tx_and_cpm = filtered_counts_data_on_tx_and_cpm.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("After transcript read count filter:")
        print("\tNumber of transcripts", len(filtered_counts_data_on_tx_and_cpm.index))
        print("\tNumber of genes: ", len(on_tx_and_cpm.gene_id.unique()))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_counts_data_on_tx_and_cpm, tx2gene, args.genefilter_count, args.genefilter_sample)
    filtered_inds_on_IF = filter_on_IF(IFs, args.n_small, args.if_fraction)
    filtered_count_data_on_IFs = tx_count_data.loc[filtered_inds_on_IF,:]
    on_ifs = filtered_count_data_on_IFs.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("After filtering on isoform fraction:")
        print("\tNumber of transcripts", len(filtered_count_data_on_IFs.index))
        print("\tNumber of genes: ", len(on_ifs.gene_id.unique()))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_count_data_on_IFs, tx2gene, args.genefilter_count, args.genefilter_sample)
    filtered_count_data_on_IFs = filtered_count_data_on_IFs[filtered_count_data_on_IFs.index.isin(IFs.index.to_list())]
    filtered_ind_on_isoform_count = filter_on_isoform_count(filtered_count_data_on_IFs, tx2gene)
    filtered_count_data_on_isoform_count = tx_count_data.loc[filtered_ind_on_isoform_count,:]
    on_iso_count = filtered_count_data_on_isoform_count.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("After filtering on isoform count:")
        print("\tNumber of transcripts", len(filtered_count_data_on_isoform_count.index))
        print("\tNumber of genes: ", len(on_iso_count.gene_id.unique()))
    filtered_IF_data_on_isoform_count = IFs.loc[filtered_ind_on_isoform_count,:]
    final_filtered_tx_counts = filtered_count_data_on_isoform_count.join(tx2gene.set_index('tx_id'))
    final_filtered_tx_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_tx_counts.txt"), sep = '\t')
    final_filtered_ifs = filtered_IF_data_on_isoform_count.join(tx2gene.set_index('tx_id'))
    final_filtered_ifs.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_ifs.txt"), sep = '\t')
    final_gene_ids = filtered_IF_data_on_isoform_count.join(tx2gene.set_index('tx_id')).gene_id.to_list()
    final_filtered_gene_counts = gene_level_counts[gene_level_counts.index.isin(final_gene_ids)]
    final_filtered_gene_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_gene_counts.txt"), sep = '\t')
    
    if(args.p_dom > 0):
        final_filtered_ifs = final_filtered_ifs.round(2)
        selected_dom_iso_genes = select_genes_w_dom_iso(final_filtered_ifs, ctrl_samples, args.p_dom)
        dominance_selected_ifs = final_filtered_ifs[final_filtered_ifs.gene_id.isin(selected_dom_iso_genes)]
        dominance_selected_gene_counts = final_filtered_gene_counts[final_filtered_gene_counts.index.isin(selected_dom_iso_genes)]
        dominance_selected_tx_counts = final_filtered_tx_counts[final_filtered_tx_counts.index.isin(selected_dom_iso_genes)]
        if(args.write):
            print("After selecting for genes with a dominant transcript in control group:")
            print("\tNumber of transcripts", len(dominance_selected_ifs.index))
            print("\tNumber of genes: ", len(dominance_selected_ifs.gene_id.unique()))
        dominance_selected_ifs.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_ifs.txt"), sep = '\t')
        dominance_selected_gene_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_gene_counts.txt"), sep = '\t')
        dominance_selected_tx_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_tx_counts.txt"), sep = '\t')