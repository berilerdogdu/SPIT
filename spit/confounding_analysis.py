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

import os
import warnings
import pandas as pd
import numpy as np
import matplotlib
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
    ifs_df = ifs_df.drop(columns = ['gene_id'])
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

def build_random_forest_regr(cluster_m_joint_dtu, cov_names, covs_arr, ifs_dtu_arr, n_small, plot, plot_dir):
    regr = RandomForestRegressor(n_estimators = 50, max_depth=1, criterion='absolute_error', bootstrap = True, min_samples_split = n_small)
    cluster_m_joint_dtu_arr = cluster_m_joint_dtu.to_numpy()
    tx_ids = cluster_m_joint_dtu.index.to_list()
    sig_tx_inds = []
    num_of_txs = len(cluster_m_joint_dtu_arr)
    for tx in tqdm(range(num_of_txs)):
        if(plot):
            fig, ax = plt.subplots()
        tx_cls_arr = cluster_m_joint_dtu_arr[tx]
        covs_tx = covs_arr[~np.isnan(tx_cls_arr)]
        tx_if = ifs_dtu_arr[tx]
        tx_tif = tx_if[~np.isnan(tx_cls_arr)]
        tx_cls_arr = tx_cls_arr[~np.isnan(tx_cls_arr)]
        tx_all_attr = np.append(covs_tx,np.array(tx_cls_arr).reshape(len(tx_cls_arr), 1),1)
        regr.fit(tx_all_attr, tx_tif)
        result = permutation_importance(regr, tx_all_attr, tx_tif, n_repeats=50, random_state=42)
        importances = pd.DataFrame(
        result.importances.T,
        columns = cov_names + ['spit_v'])
        upper_quartile_arr = importances.apply(get_upper_quartile, axis = 0).to_numpy()
        lower_quartile_arr = importances.apply(get_lower_quartile, axis = 0).to_numpy()
        if(lower_quartile_arr[-1] > max(upper_quartile_arr[:-1])):
            sig_tx_inds.append(tx)
            if(plot):
                ax = importances.plot.box(vert=False, ax = ax,
                                        whiskerprops = dict(linestyle='-', color='steelblue'),
                                        boxprops = dict(linestyle='-', color='steelblue'),
                                        medianprops = dict(linestyle='-', color='green'))
        elif(plot):
            ax = importances.plot.box(vert=False, ax = ax,
                                        whiskerprops = dict(linestyle='-', color='red'),
                                        boxprops = dict(linestyle='-', color='red'),
                                        medianprops = dict(linestyle='-', color='red'))
        if(plot):
            ax.set_title("Permutation Importances - " + tx_ids[tx])
            ax.axvline(x=0, color="k", linestyle="--")
            ax.set_xlabel("Decrease in accuracy score")
            plt.savefig(os.path.join(plot_dir, tx_ids[tx]+'.png'))
            plt.close(fig)
    return sig_tx_inds

def write_post_filter_results(sig_tx_inds, cluster_m_joint_dtu, spit_cluster_matrix, tx2gene_dict, spit_out_file, controlled_spit_cluster_matrix, controlled_spit_out):
    sig_tx_ids = cluster_m_joint_dtu.iloc[sig_tx_inds,].index
    sig_gene_ids = [tx2gene_dict[i] for i in sig_tx_ids]
    spit_out_df = pd.read_csv(spit_out_file, sep = '\t', header = None, names = ["gene_id", "flag"])
    spit_out_post_confounding = spit_out_df[spit_out_df.gene_id.isin(sig_gene_ids)]
    spit_out_post_confounding = spit_out_post_confounding.fillna('')
    spit_out_post_confounding.to_csv(controlled_spit_out, sep = '\t', index = False)
    spit_cluster_m_dtu = spit_cluster_matrix[(spit_cluster_matrix == 1).any(axis = 1)].copy(deep=True)
    gene_ids = [tx2gene_dict[i] for i in spit_cluster_m_dtu.index]
    spit_cluster_m_dtu['gene_id'] = gene_ids
    cluster_m_post_conf = spit_cluster_m_dtu.loc[sig_tx_ids]
    cluster_m_post_conf.drop(columns=['gene_id']).to_csv(controlled_spit_cluster_matrix, sep = '\t')
    return


def main(args):

    matplotlib.use('Agg')
    warnings.simplefilter(action='ignore', category=FutureWarning)
    IFs_df = pd.read_csv(args.i, sep = '\t', index_col = 0)
    tx2gene_dict = IFs_df.to_dict()['gene_id']
    pheno_df = pd.read_csv(args.l, sep = '\t')
    spit_cluster_m = pd.read_csv(os.path.join(args.O, "SPIT_analysis", "spit_cluster_matrix.txt"), sep = '\t', index_col = 0)
    cluster_m_joint_dtu = make_joint_cluster_mat(spit_cluster_m, pheno_df)
    ifs_dtu_arr = make_if_arr(IFs_df, cluster_m_joint_dtu)
    cov_names, covs_arr = make_cov_arr(pheno_df)
    confounding_analysis_plot_dir = os.path.join(args.O, "SPIT_analysis", "confounding_analysis_plots")
    if(args.plot and (os.path.exists(confounding_analysis_plot_dir) == False)):
        os.mkdir(confounding_analysis_plot_dir)
    sig_tx_inds = build_random_forest_regr(cluster_m_joint_dtu, cov_names, covs_arr, ifs_dtu_arr, args.n_small, args.plot, confounding_analysis_plot_dir)
    write_post_filter_results(sig_tx_inds, cluster_m_joint_dtu, spit_cluster_m, tx2gene_dict, os.path.join(args.O, "SPIT_analysis", "spit_out.txt"), os.path.join(args.O, "SPIT_analysis", "controlled_spit_cluster_matrix.txt"), os.path.join(args.O, "SPIT_analysis", "controlled_spit_out.txt"))
