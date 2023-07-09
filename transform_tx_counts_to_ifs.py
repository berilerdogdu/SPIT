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

import sys
import argparse
import numpy as np
import pandas as pd
import random
import math
from collections import defaultdict


def convert_counts_to_IF_and_gene_level(counts): 
    counts_w_genes = counts
    gene_level_counts = counts_w_genes.groupby('gene_id').sum()
    gene_level_counts = gene_level_counts + 0.00001
    counts_w_genes_multilevel = counts_w_genes.set_index([counts_w_genes.index, 'gene_id'])
    IFs = counts_w_genes_multilevel.div(gene_level_counts,axis='index',level='gene_id')
    return IFs, gene_level_counts


def main(argv):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-i', metavar='tx_counts.txt', required=True, type=str, help='Transcript level counts file (tsv)')
    parser.add_argument('-m', metavar='tx2gene.txt', required=True, type=str, help='Transcript to gene mapping file (tsv)')
    parser.add_argument('-F', metavar='filtered_ifs.txt', required=True, type=str, help='Output file path for isoform fractions (IFs)')
    parser.add_argument('-G', metavar='filtered_gene_counts.txt', required=True, type=str, help='Output file path for gene counts')
    parser.add_argument('-w', '--write', action='store_true', help='Write the number of transcripts & genes left after filtering all-zeroes to stdout.')
    
    args = parser.parse_args()    
    tx_count_data = pd.read_csv(args.i, sep='\t', index_col=0)
    tx2gene = pd.read_csv(args.m, sep = '\t')
    txs_w_genes = tx_count_data.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("Input number of transcripts: ", txs_w_genes.shape[0])
        print("Input number of genes: ", len(txs_w_genes.gene_id.unique()))
    tx_count_data = tx_count_data[tx_count_data.values.sum(axis=1) != 0]
    txs_w_genes = tx_count_data.join(tx2gene.set_index('tx_id'))
    if(args.write):
        print("After filtering out all-zeros:")
        print("\tNumber of transcripts", len(txs_w_genes.index))
        print("\tNumber of genes: ", len(txs_w_genes.gene_id.unique()))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(txs_w_genes)
    IFs.to_csv(args.F, sep = '\t')
    gene_level_counts.to_csv(args.G, sep = '\t')


if __name__ == "__main__":
   main(sys.argv[1:])
