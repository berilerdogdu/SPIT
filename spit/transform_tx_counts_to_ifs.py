#!/usr/bin/env python3

"""
Description:
    Directly convert transcript counts to Isoform Fractions (IF). Transcripts with no non-zero counts will still be filtered out. This module is useful when skipping pre-filtering step.
Usage:
    ./transform_tx_counts_to_ifs.py -i <tx_counts.txt> -m <tx2gene.txt> -F <filtered_ifs.txt> -G <filtered_gene_counts.txt> --write
Author:
    Beril Erdogdu
Date:
    July 08, 2023
"""

import os
import numpy as np
import pandas as pd


def filter_on_isoform_count(counts_w_genes):
    gene_isoform_counts = counts_w_genes.groupby('gene_id').count().iloc[:, 1]
    filtered_genes = gene_isoform_counts[gene_isoform_counts > 1].index.to_list()
    return counts_w_genes[counts_w_genes.gene_id.isin(filtered_genes)].index.to_list()
    
def convert_counts_to_IF_and_gene_level(counts): 
    counts_w_genes = counts
    gene_level_counts = (counts_w_genes.groupby('gene_id').sum() + 0.00001).astype(np.float32)
    counts_w_genes_multilevel = counts_w_genes.set_index([counts_w_genes.index, 'gene_id'])
    IFs = counts_w_genes_multilevel.div(gene_level_counts,axis='index',level='gene_id').reset_index(level='gene_id')
    if_num_cols = IFs.select_dtypes(include=[np.number]).columns
    IFs[if_num_cols] = IFs[if_num_cols].astype(np.float32).round(3)
    gene_num_cols = gene_level_counts.select_dtypes(include=[np.number]).columns
    gene_level_counts[gene_num_cols] = gene_level_counts[gene_num_cols].astype(np.int32)
    return IFs, gene_level_counts


def main(args):
    tx_count_data = pd.read_csv(args.i, sep='\t', index_col=0)
    tx_count_data = tx_count_data.astype(np.int32)
    tx2gene = pd.read_csv(args.m, sep = '\t').set_index('tx_id')
    txs_w_genes = tx_count_data.join(tx2gene)
    
    if(args.write):
        print("Input number of transcripts: ", txs_w_genes.shape[0])
        print("Input number of genes: ", len(txs_w_genes.gene_id.unique()))
    tx_count_data = tx_count_data[tx_count_data.values.sum(axis=1) != 0]
    txs_w_genes = tx_count_data.join(tx2gene)
    if(args.write):
        print("After filtering out all-zeros:")
        print("\tNumber of transcripts", len(txs_w_genes.index))
        print("\tNumber of genes: ", len(txs_w_genes.gene_id.unique()))
    
    filtered_ind_on_isoform_count = filter_on_isoform_count(txs_w_genes)
    filtered_count_data_on_isoform_count = txs_w_genes.loc[filtered_ind_on_isoform_count,:]
    if(args.write):
        print("After filtering on isoform count:")
        print("\tNumber of transcripts", len(filtered_count_data_on_isoform_count.index))
        print("\tNumber of genes: ", len(filtered_count_data_on_isoform_count.gene_id.unique()))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_count_data_on_isoform_count)
    if_sample_cols = [c for c in IFs.columns if c != 'gene_id']
    IFs[if_sample_cols] = IFs[if_sample_cols].astype(np.float32)
    IFs.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_ifs.txt"), sep = '\t')
    filtered_count_data_on_isoform_count = filtered_count_data_on_isoform_count.astype(np.int32)
    filtered_count_data_on_isoform_count.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_tx_counts.txt"), sep = '\t')
    gene_level_counts = gene_level_counts.astype(np.int32)
    gene_level_counts.to_csv(os.path.join(args.O, "SPIT_analysis", "filtered_gene_counts.txt"), sep = '\t')
