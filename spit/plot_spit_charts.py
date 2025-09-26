import os
import sys
import sys
import shutil
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
try:
    from IPython.display import display, Image as IPImage
    _HAS_IPYTHON = True
except Exception:
    _HAS_IPYTHON = False
import seaborn as sns
from tqdm import tqdm


def check_confounding_control(dir):
    controlled_matrix = os.path.join(dir, "SPIT_analysis", "controlled_spit_cluster_matrix.txt")
    orgnl_matrix = os.path.join(dir, "SPIT_analysis", "spit_cluster_matrix.txt")
    if os.path.exists(controlled_matrix):
        return controlled_matrix, orgnl_matrix
    elif os.path.exists(orgnl_matrix):
        return False, orgnl_matrix
    else:
        print("Error: SPIT DTU output not found.")
        sys.exit()

def spit_chart(args):
    plt.rcParams['figure.figsize'] = [8, 10]
    confounding_controlled, orig_matrix = check_confounding_control(args.O)
    pre_conf_cluster_m = pd.read_csv(orig_matrix, sep = '\t', index_col = 0)
    pre_conf_dtu_txs = pre_conf_cluster_m[(pre_conf_cluster_m == 1).any(axis = 1)].index.to_list()
    dtu_comp_p_values = pd.read_csv(os.path.join(args.O, "SPIT_analysis", "all_p_values.txt"), sep = '\t', index_col = 0, header = 0, names = ['p_value'])
    perm_p_values = pd.read_csv(os.path.join(args.O, "SPIT_analysis", "perm_p_medians.txt"), sep = '\t', index_col = 0, header = 0, names = ['perm_min_p_value'])
    joint_p_df = perm_p_values.join(dtu_comp_p_values)
    def neglog10(series):
        s = series.astype(float).copy()
        pos = s[s > 0]
        if len(pos) == 0:
            eps = 1e-300
            s[:] = eps
        else:
            eps = pos.min()
            s[s <= 0] = eps
        return -np.log10(s)
    joint_p_dtu_pre_conf_df = joint_p_df.loc[pre_conf_dtu_txs]
    non_dtu_txs = list(set(joint_p_df.index.to_list()).difference(set(pre_conf_dtu_txs)))
    joint_p_non_dtu_df = joint_p_df.loc[non_dtu_txs]
    plt.scatter(neglog10(joint_p_non_dtu_df['perm_min_p_value']), neglog10(joint_p_non_dtu_df['p_value']), s = 2, color = 'forestgreen', label='Insignificant transcripts')
    if(confounding_controlled):
        sig_tx_ids = pd.read_csv(confounding_controlled, sep = '\t', index_col = 0).index.to_list()
        joint_p_dtu_df = joint_p_df.loc[sig_tx_ids]
        plt.scatter(neglog10(joint_p_dtu_pre_conf_df['perm_min_p_value']), neglog10(joint_p_dtu_pre_conf_df['p_value']), s = 2, marker = 'o', facecolors='none', color = 'purple', label='Transcripts eliminated by confounding control')
        plt.scatter(neglog10(joint_p_dtu_df['perm_min_p_value']), neglog10(joint_p_dtu_df['p_value']), s = 2, marker = 'o', color = 'crimson', label='DTU transcripts')
    else:
        sig_tx_ids = pre_conf_dtu_txs
        plt.scatter(neglog10(joint_p_dtu_pre_conf_df['perm_min_p_value']), neglog10(joint_p_dtu_pre_conf_df['p_value']), s = 2, marker = 'o', facecolors='none', color = 'crimson', label='DTU transcripts')
    plt.axhline(y=-np.log10(float(args.p_cutoff)), linestyle='dashed', linewidth=1, label='SPIT-Test p-value cutoff')
    plt.ylabel('$-log_{10}$(p-value)', fontsize = 12)
    plt.xlabel('$-log_{10}$(Minimum p-value from SPIT-Test iterations)', fontsize = 12)
    plt.title("SPIT-Chart for Sample Analysis", fontsize = 14)
    plt.legend(fontsize=10)
    output_png = os.path.join(args.O, "SPIT_analysis", "SPIT_chart.png")
    plt.savefig(output_png)
    plt.close()
    # If running in a notebook/IPython, also display inline
    if _HAS_IPYTHON and hasattr(sys.modules.get('IPython'), 'get_ipython'):
        try:
            display(IPImage(filename=output_png))
        except Exception:
            pass
    return

def plot_violins(args):
    print("Generating violin plots for DTU transcripts:")
    dir_name = os.path.join(args.O, "SPIT_analysis", "violin_plots")
    if os.path.exists(dir_name) == True:
        shutil.rmtree(dir_name)
    os.mkdir(dir_name)
    plt.rcParams['figure.figsize'] = [8, 5]
    pheno_df = pd.read_csv(args.l, sep='\t')
    ctrl_ids = pheno_df[pheno_df.condition == 0].id.to_list()
    ifs = pd.read_csv(args.i, sep = '\t', index_col = 0)
    if_cols = ifs.select_dtypes(include=[np.number]).columns
    ifs[if_cols] = ifs[if_cols].astype(np.float32).round(3)
    confounding_controlled, orig_matrix = check_confounding_control(args.O)
    if(confounding_controlled):
        clm_dtu = pd.read_csv(confounding_controlled, sep = '\t', index_col = 0)
    else:
        clm_dtu = pd.read_csv(orig_matrix, sep = '\t', index_col = 0)
    sig_tx_ids = clm_dtu[clm_dtu.any(axis = 1)].index.to_list()
    for i in tqdm(sig_tx_ids):
        spit_v = clm_dtu.loc[i]
        dtu_plus = spit_v.index[(spit_v == 1)]
        df_index = ctrl_ids + dtu_plus.to_list()
        df_class = ["Ctrl" for i in range(len(ctrl_ids))] + ["DTU+" for i in range(len(dtu_plus))]
        df_if = ifs.loc[i][ctrl_ids].to_list() + ifs.loc[i][dtu_plus].to_list()
        tx_df = pd.DataFrame({'sample': df_index, 'group': df_class, 'IF':df_if})
        custom_palette = ["darkgreen", "pink"]
        sns.set_palette(custom_palette)
        sns.violinplot(data=tx_df, x="IF", y="group", hue="group")
        plt.title(i + " - Isoform Fraction Violin Plots", size = 10, pad = 10)
        plt.savefig(os.path.join(args.O, "SPIT_analysis", "violin_plots", i+".png"))
        plt.close()
    return
        
def main(args):
    matplotlib.use('Agg')
    spit_chart(args)
    plot_violins(args)
