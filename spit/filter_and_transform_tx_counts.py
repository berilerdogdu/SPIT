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


def filter_on_cpm(counts, n):
    lib_sizes = counts.sum(axis = 0).to_list()
    cpm_values = (counts.divide(lib_sizes) * 1e6).astype(np.float32)
    mask = cpm_values > 1.0
    filtered_cpm_row_ind = counts.index[np.flatnonzero(mask.sum(axis=1) >= n)]
    return filtered_cpm_row_ind

def filter_on_tx_count(counts, pr_fraction, ctrl, case):
    filtered_tx_row_ind = counts.index.to_list()
    ctrl_zero_txs = counts[((counts[ctrl] == 0).astype(int).sum(axis=1) >= (len(ctrl)*pr_fraction))].index.to_list()
    case_zero_txs = counts[((counts[case] == 0).astype(int).sum(axis=1) >= (len(case)*pr_fraction))].index.to_list()
    both_zero_txs = set(ctrl_zero_txs).intersection(set(case_zero_txs))
    filtered_tx_row_ind = list(set(filtered_tx_row_ind).difference(both_zero_txs))
    return filtered_tx_row_ind

def convert_counts_to_IF_and_gene_level(counts, tx2gene, genefilter_count, genefilter_sample):
    counts_w_genes = counts.join(tx2gene)
    gene_level_counts = (counts_w_genes.groupby('gene_id').sum() + 0.00001).astype(np.float32)
    gene_threshold_bool = gene_level_counts[gene_level_counts < genefilter_count].count(axis = 1) < genefilter_sample
    valid_genes = gene_threshold_bool[gene_threshold_bool].index
    gene_level_counts = gene_level_counts[gene_level_counts.index.isin(valid_genes)]
    counts_w_genes = counts_w_genes[counts_w_genes.gene_id.isin(valid_genes)]
    counts_w_genes_multilevel = counts_w_genes.set_index([counts_w_genes.index, 'gene_id'])
    IFs = counts_w_genes_multilevel.div(gene_level_counts,axis='index',level='gene_id').droplevel('gene_id')
    if_num_cols = IFs.select_dtypes(include=[np.number]).columns
    IFs[if_num_cols] = IFs[if_num_cols].astype(np.float32).round(3)
    gene_num_cols = gene_level_counts.select_dtypes(include=[np.number]).columns
    gene_level_counts[gene_num_cols] = gene_level_counts[gene_num_cols].astype(np.int32)
    return IFs, gene_level_counts
    
def filter_on_IF(IFs, n, f):
    filter_count = (IFs > f).sum(axis = 1)
    filtered_IFs_row_ind = IFs[filter_count >= n].index.to_list()
    return filtered_IFs_row_ind

def filter_on_isoform_count(counts, tx2gene):
    counts_w_genes = counts.join(tx2gene)
    gene_isoform_counts = counts_w_genes.groupby('gene_id').count().iloc[:, 1]
    filtered_genes = gene_isoform_counts[gene_isoform_counts > 1].index.to_list()
    return counts_w_genes[counts_w_genes.gene_id.isin(filtered_genes)].index.to_list()
    
def select_genes_w_dom_iso(IFs, ctrl_samples, p_dom):
    print("...Selecting genes with consistently dominant isoforms in control group...")
    IFs = IFs[[*ctrl_samples, 'gene_id']]
    IFs.insert(0, "ctrl_IF_mean", IFs[ctrl_samples].mean(axis = 1))
    ctrl_IF_max = IFs.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])
    ctrl_IF_no_max = IFs.drop(ctrl_IF_max.index)
    ctrl_IF_second_max = ctrl_IF_no_max.sort_values('ctrl_IF_mean', ascending=False).drop_duplicates(['gene_id'])    
    max_isoforms = ctrl_IF_max.drop(columns=['ctrl_IF_mean']).set_index('gene_id')
    second_max_isoforms = ctrl_IF_second_max.drop(columns=['ctrl_IF_mean']).set_index('gene_id')
    second_max_isoforms = second_max_isoforms.loc[max_isoforms.index]
    dominance_counts = (max_isoforms > second_max_isoforms).sum(axis=1)
    threshold = len(ctrl_samples) * p_dom
    dominant_genes = max_isoforms.index[dominance_counts >= threshold]
    return dominant_genes.tolist()


def main(args):
    tx_count_data = pd.read_csv(args.i, sep='\t', index_col=0)
    tx_count_data = tx_count_data.astype(np.int32)
    tx2gene = pd.read_csv(args.m, sep = '\t').set_index('tx_id')
    if(args.write):
        print("Input number of transcripts: ", tx_count_data.shape[0])
        print("Input number of genes: ", len(tx2gene.loc[tx_count_data.index].gene_id.unique()))
    tx_count_data = tx_count_data[tx_count_data.values.sum(axis=1) != 0]
    if(args.write):
        print("After filtering out all-zeros:")
        print("\tNumber of transcripts", len(tx_count_data.index))
        print("\tNumber of genes: ", len(tx2gene.loc[tx_count_data.index].gene_id.unique()))
    pheno = pd.read_csv(args.l, sep='\t')
    ctrl_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = pheno[pheno.condition == 1].id.to_list()
    filtered_inds_on_cpm = filter_on_cpm(tx_count_data, args.n_small)
    filtered_count_data_on_cpm = tx_count_data.loc[filtered_inds_on_cpm,:]
    del tx_count_data
    if(args.write):
        print("After CPM-filter:")
        print("\tNumber of transcripts", len(filtered_count_data_on_cpm.index))
        print("\tNumber of genes: ", len(filtered_count_data_on_cpm.join(tx2gene).gene_id.unique()))
    filtered_inds_on_tx_and_cpm = filter_on_tx_count(filtered_count_data_on_cpm, args.pr_fraction, ctrl_samples, case_samples)
    filtered_counts_data_on_tx_and_cpm = filtered_count_data_on_cpm.loc[filtered_inds_on_tx_and_cpm,:]
    del filtered_count_data_on_cpm
    if(args.write):
        print("After transcript read count filter:")
        print("\tNumber of transcripts", len(filtered_counts_data_on_tx_and_cpm.index))
        print("\tNumber of genes: ", len(tx2gene.loc[filtered_counts_data_on_tx_and_cpm.index].gene_id.unique()))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_counts_data_on_tx_and_cpm, tx2gene, args.genefilter_count, args.genefilter_sample)
    filtered_inds_on_IF = filter_on_IF(IFs, args.n_small, args.if_fraction)
    filtered_count_data_on_IFs = filtered_counts_data_on_tx_and_cpm.loc[filtered_inds_on_IF,:]
    del filtered_counts_data_on_tx_and_cpm
    if(args.write):
        print("After filtering on isoform fraction:")
        print("\tNumber of transcripts", len(filtered_count_data_on_IFs.index))
        print("\tNumber of genes: ", len(tx2gene.loc[filtered_count_data_on_IFs.index].gene_id.unique()))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_count_data_on_IFs, tx2gene, args.genefilter_count, args.genefilter_sample)
    filtered_count_data_on_IFs = filtered_count_data_on_IFs[filtered_count_data_on_IFs.index.isin(IFs.index.to_list())]
    filtered_ind_on_isoform_count = filter_on_isoform_count(filtered_count_data_on_IFs, tx2gene)
    filtered_count_data_on_isoform_count = filtered_count_data_on_IFs.loc[filtered_ind_on_isoform_count,:]
    del filtered_count_data_on_IFs
    if(args.write):
        print("After filtering on isoform count:")
        print("\tNumber of transcripts", len(filtered_count_data_on_isoform_count.index))
        print("\tNumber of genes: ", len(tx2gene.loc[filtered_count_data_on_isoform_count.index].gene_id.unique()))
    filtered_IF_data_on_isoform_count = IFs.loc[filtered_ind_on_isoform_count,:]
    final_filtered_tx_counts = filtered_count_data_on_isoform_count.join(tx2gene)
    # Ensure integer counts
    value_cols_tx = final_filtered_tx_counts.columns.difference(['gene_id'])
    final_filtered_tx_counts[value_cols_tx] = final_filtered_tx_counts[value_cols_tx].astype(np.int32)
    final_filtered_tx_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_tx_counts.txt"), sep = '\t')
    final_filtered_ifs = filtered_IF_data_on_isoform_count.join(tx2gene)
    # Ensure float32 IFs
    value_cols_ifs = final_filtered_ifs.columns.difference(['gene_id'])
    final_filtered_ifs[value_cols_ifs] = final_filtered_ifs[value_cols_ifs].astype(np.float32)
    final_filtered_ifs.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_ifs.txt"), sep = '\t')
    final_gene_ids = filtered_IF_data_on_isoform_count.join(tx2gene).gene_id.to_list()
    final_filtered_gene_counts = gene_level_counts[gene_level_counts.index.isin(final_gene_ids)]
    final_filtered_gene_counts = final_filtered_gene_counts.astype(np.int32)
    final_filtered_gene_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_gene_counts.txt"), sep = '\t')
    
    if(args.p_dom > 0):
        selected_dom_iso_genes = select_genes_w_dom_iso(final_filtered_ifs, ctrl_samples, args.p_dom)
        dominance_selected_ifs = final_filtered_ifs[final_filtered_ifs.gene_id.isin(selected_dom_iso_genes)].copy()
        dominance_selected_gene_counts = final_filtered_gene_counts[final_filtered_gene_counts.index.isin(selected_dom_iso_genes)].copy()
        dominance_selected_tx_counts = final_filtered_tx_counts[final_filtered_tx_counts.index.isin(dominance_selected_ifs.index)].copy()
        if(args.write):
            print("After selecting for genes with a dominant transcript in control group:")
            print("\tNumber of transcripts", len(dominance_selected_ifs.index))
            print("\tNumber of genes: ", len(dominance_selected_ifs.gene_id.unique()))
        dom_if_val_cols = dominance_selected_ifs.select_dtypes(include=[np.number]).columns
        dominance_selected_ifs.loc[:, dom_if_val_cols] = dominance_selected_ifs.loc[:, dom_if_val_cols].astype(np.float32)
        dominance_selected_ifs.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_ifs.txt"), sep = '\t')
        dominance_selected_gene_counts = dominance_selected_gene_counts.astype(np.int32)
        dominance_selected_gene_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_gene_counts.txt"), sep = '\t')
        dom_tx_val_cols = dominance_selected_tx_counts.select_dtypes(include=[np.number]).columns
        dominance_selected_tx_counts.loc[:, dom_tx_val_cols] = dominance_selected_tx_counts.loc[:, dom_tx_val_cols].astype(np.int32)
        dominance_selected_tx_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "dominance_selected_tx_counts.txt"), sep = '\t')
