import os
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from matplotlib.backends.backend_pdf import PdfPages


def make_true_dtu_sets(args, exps, tx_2_gene_dict):
    exp_dtu_genes_dict = defaultdict()
    exp_dtu_txs_dict = defaultdict()
    fit_param_directory_path = os.path.join(args.O, "SPIT_analysis", "parameter_fitting")
    for e in exps:
        true_cluster_arr = pd.read_csv(os.path.join(fit_param_directory_path, e, "true_cluster_matrix.txt"), sep = '\t', index_col = 0)
        dtu_txs = true_cluster_arr[(true_cluster_arr==1).any(axis=1)].index.to_list()
        dtu_genes = set()
        exp_dtu_txs_dict[e] = dtu_txs
        
        for t in dtu_txs:
            dtu_genes.add(tx_2_gene_dict[t])
        exp_dtu_genes_dict[e] = dtu_genes
    return exp_dtu_genes_dict
    
    
def apply_and_plot_loocv(args, exps, exp_dtu_genes_dict, tx_2_gene_dict):
    fit_param_directory_path = os.path.join(args.O, "SPIT_analysis", "parameter_fitting")
    pdf_file_path = os.path.join(fit_param_directory_path, "LOOCV_venn_diags.pdf")
    f_prime = 0
    plot_counter = 0
    pdf_pages = PdfPages(pdf_file_path)
    for e_out in exps:
        l_o = e_out
        f_bar_max = 0
        k_best = 0
        b_best = 0
        test_f = 0
        for b in np.round(np.arange(0.02, 0.21, 0.02), 2):
            b = format(b, '.2f')
            for k in np.round(np.arange(0.1, 1, 0.1), 1):
                f_bar = 0
                for e_in in exps:
                    if e_in is not l_o:
                        dtu_genes = exp_dtu_genes_dict[e_in]
                        spit_cluster_arr = pd.read_csv(os.path.join(fit_param_directory_path, e_in, "spit_cluster_matrix_k" + str(k) + ".b" + str(b) + ".txt"), sep='\t', index_col=0)
                        spit_dtu_txs = spit_cluster_arr[(spit_cluster_arr == 1).any(axis=1)].index.to_list()
                        spit_genes = set()
                        for t in spit_dtu_txs:
                            spit_genes.add(tx_2_gene_dict[t])
                        tp = len(dtu_genes.intersection(spit_genes))
                        fp = len(spit_genes.difference(dtu_genes))
                        fn = len(dtu_genes.difference(spit_genes))
                        f_1 = (2 * tp) / ((2 * tp) + fp + fn)
                        f_bar += f_1
                f_bar = f_bar / (args.n_exps - 1)
                if f_bar > f_bar_max:
                    f_bar_max = f_bar
                    k_best = k
                    b_best = b
        test_dtu_genes = exp_dtu_genes_dict[e_out]
        spit_cluster_arr = pd.read_csv(os.path.join(fit_param_directory_path, e_out, "spit_cluster_matrix_k" + str(k_best) + ".b" + str(b_best) + ".txt"), sep='\t', index_col=0)
        spit_dtu_txs = spit_cluster_arr[(spit_cluster_arr == 1).any(axis=1)].index.to_list()
        spit_genes = set()
        for t in spit_dtu_txs:
            spit_genes.add(tx_2_gene_dict[t])
        tp = len(test_dtu_genes.intersection(spit_genes))
        fp = len(spit_genes.difference(test_dtu_genes))
        fn = len(test_dtu_genes.difference(spit_genes))
        f_t = (2 * tp) / ((2 * tp) + fp + fn)
        output = f"Left-out experiment: {e_out}\n"
        output += f"Best parameters on training exps (k, band): {k_best}, {b_best}\n"
        output += f"F-score with above parameters on left-out exp: {f_bar_max}\n"
        output += f"TPR: {tp / len(test_dtu_genes)}\n"
        output += f"FDR: {fp / len(spit_genes)}\n"
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        ax1.axis('off')
        ax2.axis('off')
        venn2([set(spit_genes), set(test_dtu_genes)], ('SPIT Genes', 'DTU genes'), ax=ax1, set_colors=['red', 'lightseagreen'])
        ax2.text(0.5, 0.5, output, ha='center')
        pdf_pages.savefig(fig)
        plt.close(fig)
        f_prime += f_t
    f_prime = f_prime / len(exps)
    print("Final F-Score: ", f_prime)
    pdf_pages.close()
    
    return


def main(args):
    exps = ["exp" + str(i) for i in range(1, args.n_exps+1)]
    tx_2_gene= pd.read_csv(args.m, sep = '\t')
    tx_2_gene_dict = pd.Series(tx_2_gene.gene_id.values,index=tx_2_gene.tx_id).to_dict()
    
    exp_dtu_genes_dict = make_true_dtu_sets(args, exps, tx_2_gene_dict)
    apply_and_plot_loocv(args, exps, exp_dtu_genes_dict, tx_2_gene_dict)
