import sys, getopt
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

    tx_counts_file = ''
    tx2gene_file = ''
    iso_fracs_out = ''
    gene_counts_out = ''
    
    try:
        opts, args = getopt.getopt(argv,"hi:m:F:G:",["tfile=", "tx2gene_file=", "iso_fracs_out=", "gene_counts_out="])
    except getopt.GetoptError:
        print("Usage: filter_and_transform_tx_counts.py -i  <transcript level counts file (tsv)> -m <transcript to gene mapping file (tsv)> -F <output file path for filtered fractions> -G <output file path for filtered gene counts>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("Usage: filter_and_transform_tx_counts.py -i  <transcript level counts file (tsv)> -f <output file path for filtered fractions> -g <output file path for filtered gene counts>")
            sys.exit()
        elif opt in ("-i", "--tfile"):
            tx_counts_file = arg
        elif opt in ("-m", "--tx2gene_file"):
            tx2gene_file = arg
        elif opt in ("-f", "--iso_fracs_out"):
            iso_fracs_out = arg
        elif opt in ("-g", "--gene_counts_out"):
            gene_counts_out = arg
    
    tx_count_data = pd.read_csv(tx_counts_file, sep='\t', index_col='feature_id')
    tx2gene = pd.read_csv(tx2gene_file, sep = '\t')
    txs_w_genes = tx_count_data.join(tx2gene.set_index('tx_id'))
    print("Input number of transcripts: ", tx_count_data.shape[0])
    print("Input number of genes: ", len(tx_count_data.gene_id.unique()))
    
    tx_count_data = tx_count_data[tx_count_data.values.sum(axis=1) != 0]
    print("After filtering out all-zeros:")
    print("\tNumber of transcripts", len(tx_count_data.index))
    print("\tNumber of genes: ", len(tx_count_data_w_genes.gene_id.unique()))

    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(tx_count_data)
    IFs.to_csv(iso_fracs_out, sep = '\t')
    gene_level_counts.to_csv(gene_counts_out, sep = '\t')


if __name__ == "__main__":
   main(sys.argv[1:])
