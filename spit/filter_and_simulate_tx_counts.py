import numpy as np
import pandas as pd
import random

def convert_counts_to_cpm(counts):
    lib_sizes = counts.sum(axis = 0).to_list()
    scaled_counts = counts.divide(lib_sizes)
    cpms = (scaled_counts * (10**6)).astype(np.float32)
    return cpms

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

def write_simulated_samples(samples_file, samples):
    samples_file = open(samples_file,'w')
    for i in samples:
        samples_file.write(i+'\n')
    samples_file.close()
    return

def main(args):
    def log(msg):
        if getattr(args, 'write', False):
            with open(args.log_file, 'a') as lf:
                lf.write(str(msg) + '\n')
    tx_count_data = pd.read_csv(args.i, sep='\t', index_col=0)
    tx_count_data = tx_count_data.astype(np.int32)
    tx2gene = pd.read_csv(args.m, sep = '\t').set_index('tx_id')
    txs_w_genes = tx_count_data.join(tx2gene)
    log("Input number of transcripts: " + str(tx_count_data.shape[0]))
    log("Input number of genes: " + str(len(txs_w_genes.gene_id.unique())))
    tx_count_data = tx_count_data[tx_count_data.values.sum(axis=1) != 0]
    tx_count_data_w_genes = tx_count_data.join(tx2gene)
    if getattr(args, 'write', False):
        log("After filtering out all-zeros:")
        log("\tNumber of transcripts " + str(len(tx_count_data.index)))
        log("\tNumber of genes: " + str(len(tx_count_data_w_genes.gene_id.unique())))
    pheno = pd.read_csv(args.l, sep='\t')
    all_samples = pheno[pheno.condition == 0].id.to_list()
    case_samples = random.sample(all_samples, k=int(len(all_samples)/2))
    ctrl_samples = list(set(all_samples).difference(set(case_samples)))
    write_simulated_samples(args.case, case_samples)
    write_simulated_samples(args.ctrl, ctrl_samples)
    filtered_inds_on_cpm = filter_on_cpm(tx_count_data, args.n_small)
    filtered_count_data_on_cpm = tx_count_data.loc[filtered_inds_on_cpm,:]
    del tx_count_data
    if getattr(args, 'write', False):
        log("After CPM-filter:")
        log("\tNumber of transcripts " + str(len(filtered_count_data_on_cpm.index)))
        log("\tNumber of genes: " + str(len(filtered_count_data_on_cpm.join(tx2gene).gene_id.unique())))
    filtered_inds_on_tx_and_cpm = filter_on_tx_count(filtered_count_data_on_cpm, args.pr_fraction, ctrl_samples, case_samples)
    filtered_counts_data_on_tx_and_cpm = filtered_count_data_on_cpm.loc[filtered_inds_on_tx_and_cpm,:]
    del filtered_count_data_on_cpm
    if getattr(args, 'write', False):
        log("After transcript read count filter:")
        log("\tNumber of transcripts " + str(len(filtered_counts_data_on_tx_and_cpm.index)))
        log("\tNumber of genes: " + str(len(filtered_counts_data_on_tx_and_cpm.join(tx2gene).gene_id.unique())))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_counts_data_on_tx_and_cpm, tx2gene, args.genefilter_count, args.genefilter_sample)
    filtered_inds_on_IF = filter_on_IF(IFs, args.n_small, args.if_fraction)
    filtered_count_data_on_IFs = filtered_counts_data_on_tx_and_cpm.loc[filtered_inds_on_IF,:]
    del filtered_counts_data_on_tx_and_cpm
    if getattr(args, 'write', False):
        log("After filtering on isoform fraction:")
        log("\tNumber of transcripts " + str(len(filtered_count_data_on_IFs.index)))
        log("\tNumber of genes: " + str(len(tx2gene.loc[filtered_count_data_on_IFs.index].gene_id.unique())))
    IFs, gene_level_counts = convert_counts_to_IF_and_gene_level(filtered_count_data_on_IFs, tx2gene, args.genefilter_count, args.genefilter_sample)
    filtered_count_data_on_IFs = filtered_count_data_on_IFs[filtered_count_data_on_IFs.index.isin(IFs.index.to_list())]
    filtered_ind_on_isoform_count = filter_on_isoform_count(filtered_count_data_on_IFs, tx2gene)
    filtered_count_data_on_isoform_count = filtered_count_data_on_IFs.loc[filtered_ind_on_isoform_count,:]
    del filtered_count_data_on_IFs
    if getattr(args, 'write', False):
        log("After filtering on isoform count:")
        log("\tNumber of transcripts " + str(len(filtered_count_data_on_isoform_count.index)))
        log("\tNumber of genes: " + str(len(tx2gene.loc[filtered_count_data_on_isoform_count.index].gene_id.unique())))
    filtered_IF_data_on_isoform_count = IFs.loc[filtered_ind_on_isoform_count,:]
    final_tx = filtered_count_data_on_isoform_count.join(tx2gene)
    val_cols_tx = final_tx.columns.difference(['gene_id'])
    final_tx[val_cols_tx] = final_tx[val_cols_tx].astype(np.int32)
    final_tx.to_csv(args.T, sep = '\t')
    final_ifs = filtered_IF_data_on_isoform_count.join(tx2gene)
    val_cols_ifs = final_ifs.columns.difference(['gene_id'])
    final_ifs[val_cols_ifs] = final_ifs[val_cols_ifs].astype(np.float32)
    final_ifs.to_csv(args.F, sep = '\t')
    final_gene_ids = filtered_IF_data_on_isoform_count.join(tx2gene).gene_id.to_list()
    gene_level_counts = gene_level_counts.astype(np.int32)
    gene_level_counts[gene_level_counts.index.isin(final_gene_ids)].to_csv(args.G, sep = '\t')
