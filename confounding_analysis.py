#!/usr/bin/env python3

"""
Description:
    This module helps control for the confounding factors inherently present in the experimental design. By fitting a random forest regressor on each detected DTU transcript the module determines which covariates might be contributing into observable variance in isoform fraction (IF) values.
Usage:
    ./confounding_analysis.py <dominance_selected_ifs> -l <pheno> --cluster_matrix <spit_cluster_matrix> -o <spit_out> -M <controlled_spit_cluster_matrix> -O <controlled_spit_out> -P <importance_score_plots>
Author:
    Beril Erdogdu
Date:
    July 08, 2023
"""

import sys
import os
import argparse
import warnings
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
import matplotlib.backends.backend_pdf as pdf
from tqdm import tqdm


def make_joint_cluster_mat(spit_cluster_m, pheno_df):
    cluster_m_ctrl = pd.DataFrame(0, index=spit_cluster_m.index, columns=pheno_df[pheno_df.condition == 0].id.to_list())
    cluster_m_joint = spit_cluster_m.join(cluster_m_ctrl)
    cluster_m_joint_dtu = cluster_m_joint[(cluster_m_joint == 1).any(axis = 1)]
    cluster_m_joint_dtu = cluster_m_joint_dtu[pheno_df.id]
    return cluster_m_joint_dtu

def make_if_arr(ifs_df, cluster_m_joint_dtu):
    ifs_df = ifs_df.drop(columns = ['ctrl_IF_mean', 'gene_id'])
    ifs_dtu = ifs_df.loc[cluster_m_joint_dtu.index]
    ifs_dtu_arr = ifs_dtu.to_numpy()
    return ifs_dtu_arr

def make_cov_arr(pheno_df):
    cov_names = []
    cov_labels = pheno_df.columns.to_list()
    covs_arr = np.array([])
    for c in cov_labels:
        if(c == 'id' or c == 'condition'):
            continue
        if(c[-4:] == "_cat"):
            cov_names.append(c[:-4])
            c_arr = pd.factorize(pheno_df[c])[0]
        else:
            cov_names.append(c)
            c_arr = pheno_df[c].to_numpy()
        if(len(covs_arr) == 0):
            covs_arr = np.array([c_arr])
        else:
            covs_arr = np.vstack((covs_arr, c_arr))
    return cov_names, np.transpose(covs_arr)

def get_upper_quartile(data):
    upper_quartile = np.percentile(data, 75)
    return upper_quartile

def get_lower_quartile(data):
    lower_quartile = np.percentile(data, 25)
    return lower_quartile

def build_random_forest_regr(cluster_m_joint_dtu, cov_names, covs_arr, ifs_dtu_arr, n_small, plot_file):
    regr = RandomForestRegressor(n_estimators = 200, max_depth=1, criterion='absolute_error', bootstrap = True, min_samples_split = n_small)
    cluster_m_joint_dtu_arr = cluster_m_joint_dtu.to_numpy()
    tx_ids = cluster_m_joint_dtu.index.to_list()
    sig_tx_inds = []
    if(plot_file):
        pdf_pages = pdf.PdfPages(plot_file)
        fig, axes = plt.subplots(nrows=3, ncols=2, gridspec_kw={'hspace': 0.3, 'wspace': 0.5})
    for tx in tqdm(range(len(cluster_m_joint_dtu_arr))):
        tx_cls_arr = cluster_m_joint_dtu_arr[tx]
        covs_tx = covs_arr[~np.isnan(tx_cls_arr)]
        tx_if = ifs_dtu_arr[tx]
        tx_tif = tx_if[~np.isnan(tx_cls_arr)]
        tx_cls_arr = tx_cls_arr[~np.isnan(tx_cls_arr)]
        tx_all_attr = np.append(covs_tx,np.array(tx_cls_arr).reshape(len(tx_cls_arr), 1),1)
        regr.fit(tx_all_attr, tx_tif)
        result = permutation_importance(regr, tx_all_attr, tx_tif, n_repeats=100, random_state=42)
        importances = pd.DataFrame(
        result.importances.T,
        columns = cov_names + ['spit_v'])
        upper_quartile_arr = importances.apply(get_upper_quartile, axis = 0).to_numpy()
        lower_quartile_arr = importances.apply(get_lower_quartile, axis = 0).to_numpy()
        if(lower_quartile_arr[-1] > max(upper_quartile_arr[:-1])):
            sig_tx_inds.append(tx)
            if(plot_file):
                ax = importances.plot.box(vert=False, whis=10, ax=axes[tx % 3, tx % 2],
                                        whiskerprops = dict(linestyle='-',linewidth=4.5, color='steelblue'),
                                        boxprops = dict(linestyle='-',linewidth=4.5, color='steelblue'),
                                        medianprops = dict(linestyle='-',linewidth=4.5, color='green'))
        elif(plot_file):
            ax = importances.plot.box(vert=False, whis=10, ax=axes[tx % 3, tx % 2],
                                    whiskerprops = dict(linestyle='-',linewidth=4.5, color='red'),
                                    boxprops = dict(linestyle='-',linewidth=4.5, color='red'),
                                    medianprops = dict(linestyle='-',linewidth=4.5, color='red'))
        if(plot_file):
            ax.set_title("Permutation Importances - " + tx_ids[tx])
            ax.axvline(x=0, color="k", linestyle="--")
            ax.set_xlabel("Decrease in accuracy score")
            if (tx + 1) % 6 == 0 or tx == len(cluster_m_joint_dtu_arr) - 1:
                pdf_pages.savefig(fig)
                plt.close(fig)
                fig, axes = plt.subplots(nrows=3, ncols=2, gridspec_kw={'hspace': 0.3, 'wspace': 0.5})
    if(plot_file):
        pdf_pages.close()
    return sig_tx_inds

def write_post_filter_results(sig_tx_inds, cluster_m_joint_dtu, spit_cluster_matrix, tx2gene_dict, spit_out_file, controlled_spit_cluster_matrix, controlled_spit_out):
    sig_tx_ids = cluster_m_joint_dtu.iloc[sig_tx_inds,].index
    sig_gene_ids = [tx2gene_dict[i] for i in sig_tx_ids]
    spit_out_df = pd.read_csv(spit_out_file, sep = '\t', header = None, names = ["gene_id", "flag"])
    spit_out_post_confounding = spit_out_df[spit_out_df.gene_id.isin(sig_gene_ids)]
    spit_out_post_confounding = spit_out_post_confounding.fillna('')
    spit_out_post_confounding.to_csv(controlled_spit_out, sep = '\t', index = False)
    spit_cluster_m_dtu = spit_cluster_matrix[(spit_cluster_matrix == 1).any(axis = 1)]
    gene_ids = [tx2gene_dict[i] for i in spit_cluster_m_dtu.index]
    spit_cluster_m_dtu['gene_id'] = gene_ids
    cluster_m_post_conf = spit_cluster_m_dtu.loc[sig_tx_ids]
    cluster_m_post_conf.drop(columns=['gene_id']).to_csv(controlled_spit_cluster_matrix, sep = '\t')
    return


def main(argv):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-i', metavar='dom_selected_ifs.txt', required=True, type=str, help='Isoform fractions file (tsv)')
    parser.add_argument('-l', metavar='labels.txt', required=True, type=str, help='Labels/metadata file (tsv)')
    parser.add_argument('--cluster_matrix', metavar='spit_cluster_matrix.txt', required=True, type=str, help='SPIT cluster matrix')
    parser.add_argument('-n', '--n_small', metavar='12', type=int, default=12, help='Smallest sample size for the subgroups')
    parser.add_argument('-o', metavar='spit_out.txt', required=True, type=str, help='SPIT candidate DTU genes')
    parser.add_argument('-M', metavar='controlled_spit_cluster_matrix.txt', type=str, default = os.path.join(os.getcwd(), "SPIT_analysis", "controlled_spit_cluster_matrix.txt"), help='Output file path for updated SPIT cluster matrix after confounding analysis')
    parser.add_argument('-O', metavar='controlled_spit_out.txt', type=str, default = os.path.join(os.getcwd(), "SPIT_analysis", "controlled_spit_out.txt"), help='Output file path for updated SPIT candidate DTU genes after confounding analysis')
    parser.add_argument('-P', metavar='importance_score_plots.pdf', type=str, const=os.path.join(os.getcwd(), "SPIT_analysis", "controlled_spit_cluster_matrix.txt"), nargs='?', default=None, help='PDF file for plots')

    args = parser.parse_args()
    warnings.simplefilter(action='ignore', category=FutureWarning)
    if(args.P):
        plt.rcParams['figure.figsize'] = [60, 45]
        plt.rcParams.update({'font.size': 50})
        plt.rcParams['axes.labelsize'] = 50
        plt.rcParams['xtick.labelsize'] = 50
        plt.rcParams['ytick.labelsize'] = 50
        
    IFs_df = pd.read_csv(args.i, sep = '\t', index_col = 0)
    tx2gene_dict = IFs_df.to_dict()['gene_id']
    pheno_df = pd.read_csv(args.l, sep = '\t')
    spit_cluster_m = pd.read_csv(args.cluster_matrix, sep = '\t', index_col = 0)
    cluster_m_joint_dtu = make_joint_cluster_mat(spit_cluster_m, pheno_df)
    ifs_dtu_arr = make_if_arr(IFs_df, cluster_m_joint_dtu)
    cov_names, covs_arr = make_cov_arr(pheno_df)
    sig_tx_inds = build_random_forest_regr(cluster_m_joint_dtu, cov_names, covs_arr, ifs_dtu_arr, args.n_small, args.P)
    write_post_filter_results(sig_tx_inds, cluster_m_joint_dtu, spit_cluster_m, tx2gene_dict, args.o, args.M, args.O)


if __name__ == "__main__":
   main(sys.argv[1:])
