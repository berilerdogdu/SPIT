#!/usr/bin/env python3

"""
Description:
    This module performs hierarchical clustering on DTU events detected by SPIT. 
    It creates dendrograms and heatmaps to visualize sample clustering patterns 
    based on differential transcript usage events.
Usage:
    ./cluster_samples.py -l <pheno> --include_shared_dtu --color_palette <palette> --color_covariate <covariate> -O <output_dir>
Author:
    Beril Erdogdu
Date:
    December 2024
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings


def perform_hclust_python(pheno_file, include_shared_dtu=True, col_palette="BuGn", cov=None, output_dir="."):

    controlled_file_path = os.path.join(output_dir, "controlled_spit_cluster_matrix.txt")
    non_controlled_file_path = os.path.join(output_dir, "spit_cluster_matrix.txt")
    
    if os.path.exists(controlled_file_path):
        spit_cluster_m = pd.read_csv(controlled_file_path, sep='\t', index_col=0)
        print("Using controlled cluster matrix")
    elif os.path.exists(non_controlled_file_path):
        spit_cluster_m = pd.read_csv(non_controlled_file_path, sep='\t', index_col=0)
        print("Using non-controlled cluster matrix")
    else:
        print("No SPIT cluster matrices found. Please run SPIT dtu before running the cluster module.")
        return None
    
    spit_cluster_m_processed = spit_cluster_m.fillna(0).replace(-1, 0)
    num_cols = spit_cluster_m_processed.select_dtypes(include=[np.number]).columns
    spit_cluster_m_processed[num_cols] = spit_cluster_m_processed[num_cols].astype(np.int8)
    spit_cluster_m_original = spit_cluster_m.copy()
    
    if not include_shared_dtu:
        rows_all_zeros = (spit_cluster_m_processed == 0).all(axis=1)
        rows_no_zeros = (spit_cluster_m_original != 0).all(axis=1)
        rows_to_remove = rows_all_zeros | rows_no_zeros
        spit_cluster_m_processed = spit_cluster_m_processed[~rows_to_remove]
        spit_cluster_m_original = spit_cluster_m_original[~rows_to_remove]
        
        if len(spit_cluster_m_processed) == 0:
            print("No DTU events remain after filtering.")
            return None
    
    if (spit_cluster_m_processed == 1).all().all():
        print("All DTU events are shared amongst case samples, cannot apply clustering.")
        return None
    
    mask = spit_cluster_m_original == -1
    
    palette_map = {
        "Blues": "Blues",
        "Reds": "Reds",
        "Greens": "Greens", 
        "Oranges": "Oranges",
        "Purples": "Purples",
        "Greys": "Greys"
    }
    seaborn_palette = palette_map[col_palette]
    cmap = sns.color_palette(seaborn_palette, as_cmap=True)
    cmap.set_bad(color='red')
    col_colors = None
    if cov and cov != "FALSE":
        try:
            cov += '_cat'
            metadata = pd.read_csv(pheno_file, sep='\t')
            id_order = spit_cluster_m_processed.columns.tolist()
            metadata_reordered = metadata.set_index('id').reindex(id_order)
            cov_values = metadata_reordered[cov].values
            unique_cov = metadata[cov].unique()
            n_colors = len(unique_cov)
            colors = plt.cm.managua(np.linspace(0, 1, n_colors))
            color_mapping = dict(zip(unique_cov, colors))
            col_colors = [color_mapping.get(val, 'white') for val in cov_values]
        except Exception as e:
            print(f"Warning: Could not process covariate {cov}: {e}")
            col_colors = None
            
    n_rows, n_cols = spit_cluster_m_processed.shape
    max_width = 20
    max_height = 30
    min_cell_size = 0.02
    
    fig_width = min(max_width, max(8, n_cols * min_cell_size))
    fig_height = min(max_height, max(6, n_rows * min_cell_size * 1.5))
        
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        g = sns.clustermap(
            spit_cluster_m_processed, 
            method='weighted', 
            metric='jaccard',
            cmap=cmap, 
            col_cluster=True, 
            row_cluster=True,
            mask=mask, 
            col_colors=col_colors,
            linewidths=0.02,
            figsize=(fig_width, fig_height),
            cbar_kws={'aspect': 30, 'ticks': [0, 1], 'format': '%.0f', 'orientation': 'horizontal'},
            square=False
        )

        cbar = g.cax
        cbar.set_position([0.02, 0.95, 0.15, 0.02])
        cbar.set_xlabel('DTU state', fontsize=7, labelpad=0)
        
        if col_colors is not None:
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=color_mapping[val], label=str(val)) for val in unique_cov]
            legend = g.figure.legend(handles=legend_elements, 
                                title=cov.replace('_cat', ''), 
                                loc='upper left', 
                                bbox_to_anchor=(0.02, 0.90),
                                fontsize=7,
                                title_fontsize=8)
            legend.get_frame().set_facecolor('white')
            legend.get_frame().set_alpha(0.8)
        g.ax_heatmap.set_xlabel('Samples', fontsize=12)
        g.ax_heatmap.set_ylabel('Isoforms (DTU state)', fontsize=12)
    
    output_file = os.path.join(output_dir, "spit_dendrogram.pdf")
    g.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Dendrogram saved to: {output_file}")
    return g

def call_hclust(args):
    output_dir = os.path.join(args.O, "SPIT_analysis")
    
    perform_hclust_python(
        pheno_file=args.l,
        include_shared_dtu=args.include_shared_dtu,
        col_palette=args.color_palette,
        cov=args.color_covariate,
        output_dir=output_dir
    )