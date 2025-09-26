import os
from types import SimpleNamespace
from typing import Union, Optional, Dict, Any

import pandas as pd
import numpy as np
try:
    from IPython.display import display, Image as IPImage
    _HAS_IPYTHON = True
except Exception:
    _HAS_IPYTHON = False

from spit.run_spit import handle_filter as _handle_filter
from spit.run_spit import handle_dtu as _handle_dtu
from spit.run_spit import make_output_dir as _make_output_dir
from spit.cluster_samples import perform_hclust_python as _perform_hclust


def _write_if_dataframe(
    obj: Union[str, pd.DataFrame],
    path: str,
    *,
    write_index: bool,
) -> str:
    """Return a file path for the object. If a DataFrame is provided, write it to TSV.

    - If `obj` is a path-like string, it is returned unchanged.
    - If `obj` is a DataFrame, it is saved as TSV to `path` and that path is returned.
    """
    if isinstance(obj, pd.DataFrame):
        # Ensure parent directory exists
        os.makedirs(os.path.dirname(path), exist_ok=True)
        obj.to_csv(path, sep='\t', index=write_index)
        return path
    if isinstance(obj, str):
        return obj
    raise TypeError("Expected a file path (str) or a pandas.DataFrame")


def dtu(
    labels: Union[str, pd.DataFrame],
    output_dir: str = os.getcwd(),
    *,
    # DTU options
    no_clusters: bool = False,
    n_small: Optional[int] = None,
    n_iter: int = 100,
    k: float = 0.6,
    bandwidth: float = 0.09,
    f_cpm: bool = False,
    plot: bool = False,
    quiet: bool = False,
    infReps: bool = False,
    quant_path: Optional[str] = None,
    quant_type: Optional[str] = None,
) -> Dict[str, Any]:
    """Run SPIT DTU analysis on preprocessed data.

    Requires that spit.preprocess() has been run first to generate filtered files.
    Only labels parameter is needed (for sample grouping).

    Returns a dict with keys:
    - 'spit_out_path': path to `spit_out.txt`
    - 'cluster_matrix_path': path to `spit_cluster_matrix.txt`
    - 'all_p_values_path': path to `all_p_values.txt`
    - 'spit_out': DataFrame of final results
    - 'cluster_matrix': DataFrame of cluster matrix
    - 'all_p_values': DataFrame with column 'p_value'
    """

    # Ensure output structure exists
    _make_output_dir(output_dir)
    analysis_dir = os.path.join(output_dir, "SPIT_analysis")

    # Only need labels for DTU analysis
    labels_path = _write_if_dataframe(
        labels, os.path.join(analysis_dir, "pheno.txt"), write_index=False
    )

    # Check that preprocessing has been done
    filtered_ifs_path = os.path.join(analysis_dir, "filtered_ifs.txt")
    filtered_gene_counts_path = os.path.join(analysis_dir, "filtered_gene_counts.txt")
    
    if not (os.path.exists(filtered_ifs_path) and os.path.exists(filtered_gene_counts_path)):
        raise FileNotFoundError(
            "Preprocessing files not found. Please run spit.preprocess() first to generate "
            f"filtered_ifs.txt and filtered_gene_counts.txt in {analysis_dir}"
        )

    # 2) DTU step
    # Match CLI default when not provided
    _n_small = 12 if n_small is None else n_small
    
    dtu_args = SimpleNamespace(
        # The handler will set args.i/args.g to the appropriate filtered/dominance files
        i=os.path.join(analysis_dir, "filtered_ifs.txt"),
        g=os.path.join(analysis_dir, "filtered_gene_counts.txt"),
        m=os.path.join(analysis_dir, "tx2gene.txt"),  # from preprocessing
        l=labels_path,
        O=output_dir,
        infReps=infReps,
        quant_path=quant_path,
        quant_type=quant_type,
        exp=False,
        no_clusters=no_clusters,
        n_small=_n_small,
        n_iter=n_iter,
        k=k,
        bandwidth=bandwidth,
        f_cpm=f_cpm,
        plot=plot,
        quiet=quiet,
    )
    _handle_dtu(dtu_args)

    # Collect outputs
    spit_out_path = os.path.join(analysis_dir, "spit_out.txt")
    cluster_matrix_path = os.path.join(analysis_dir, "spit_cluster_matrix.txt")
    all_p_values_path = os.path.join(analysis_dir, "all_p_values.txt")

    result = {
        "spit_out_path": spit_out_path,
        "cluster_matrix_path": cluster_matrix_path,
        "all_p_values_path": all_p_values_path,
    }

    # Read outputs if available
    if os.path.exists(spit_out_path):
        result["spit_out"] = pd.read_csv(spit_out_path, sep='\t')
    if os.path.exists(cluster_matrix_path):
        result["cluster_matrix"] = pd.read_csv(cluster_matrix_path, sep='\t', index_col=0)
    if os.path.exists(all_p_values_path):
        result["all_p_values"] = pd.read_csv(all_p_values_path, sep='\t')

    return result


def preprocess(
    tx_counts: Union[str, pd.DataFrame],
    tx2gene: Union[str, pd.DataFrame],
    labels: Union[str, pd.DataFrame],
    output_dir: str = os.getcwd(),
    *,
    no_clusters: bool = False,
    n_small: Optional[int] = None,
    keep_all_nonzeros: bool = False,
    pr_fraction: float = 0.2,
    if_fraction: float = 0.1,
    genefilter_count: int = 10,
    genefilter_sample: int = 10,
    p_dom: float = 0.75,
    write: bool = False,
    quiet: bool = False,
) -> Dict[str, Any]:
    """Run only the SPIT preprocessing (filtering + IF transform).

    Accepts file paths or DataFrames. Returns paths and DataFrames for:
    - filtered_ifs.txt
    - filtered_gene_counts.txt
    - filtered_tx_counts.txt
    - dominance_selected_*.txt (if p_dom > 0)
    """

    _make_output_dir(output_dir)
    analysis_dir = os.path.join(output_dir, "SPIT_analysis")

    tx_counts_path = _write_if_dataframe(
        tx_counts, os.path.join(analysis_dir, "input_tx_counts.txt"), write_index=True
    )
    tx2gene_path = _write_if_dataframe(
        tx2gene, os.path.join(analysis_dir, "tx2gene.txt"), write_index=False
    )
    labels_path = _write_if_dataframe(
        labels, os.path.join(analysis_dir, "pheno.txt"), write_index=False
    )

    _n_small = 12 if n_small is None else n_small

    preprocess_args = SimpleNamespace(
        i=tx_counts_path,
        m=tx2gene_path,
        l=labels_path,
        no_clusters=no_clusters,
        n_small=_n_small,
        O=output_dir,
        quiet=quiet,
        write=write,
        keep_all_nonzeros=keep_all_nonzeros,
        pr_fraction=pr_fraction,
        if_fraction=if_fraction,
        genefilter_count=genefilter_count,
        genefilter_sample=genefilter_sample,
        p_dom=p_dom,
    )
    _handle_filter(preprocess_args)

    paths = {
        "filtered_ifs_path": os.path.join(analysis_dir, "filtered_ifs.txt"),
        "filtered_gene_counts_path": os.path.join(analysis_dir, "filtered_gene_counts.txt"),
        "filtered_tx_counts_path": os.path.join(analysis_dir, "filtered_tx_counts.txt"),
        "dominance_selected_ifs_path": os.path.join(analysis_dir, "dominance_selected_ifs.txt"),
        "dominance_selected_gene_counts_path": os.path.join(analysis_dir, "dominance_selected_gene_counts.txt"),
        "dominance_selected_tx_counts_path": os.path.join(analysis_dir, "dominance_selected_tx_counts.txt"),
    }

    result: Dict[str, Any] = {**paths}
    # Attach DataFrames if available
    if os.path.exists(paths["filtered_ifs_path"]):
        df = pd.read_csv(paths["filtered_ifs_path"], sep='\t', index_col=0)
        num_cols = df.select_dtypes(include=[np.number]).columns
        df[num_cols] = df[num_cols].astype(np.float32)
        result["filtered_ifs"] = df
    if os.path.exists(paths["filtered_gene_counts_path"]):
        df = pd.read_csv(paths["filtered_gene_counts_path"], sep='\t', index_col=0)
        result["filtered_gene_counts"] = df.astype(np.int32)
    if os.path.exists(paths["filtered_tx_counts_path"]):
        df = pd.read_csv(paths["filtered_tx_counts_path"], sep='\t', index_col=0)
        counts_cols = df.select_dtypes(include=[np.number]).columns
        df[counts_cols] = df[counts_cols].astype(np.int32)
        result["filtered_tx_counts"] = df
    if os.path.exists(paths["dominance_selected_ifs_path"]):
        df = pd.read_csv(paths["dominance_selected_ifs_path"], sep='\t', index_col=0)
        num_cols = df.select_dtypes(include=[np.number]).columns
        df[num_cols] = df[num_cols].astype(np.float32)
        result["dominance_selected_ifs"] = df
    if os.path.exists(paths["dominance_selected_gene_counts_path"]):
        df = pd.read_csv(paths["dominance_selected_gene_counts_path"], sep='\t', index_col=0)
        result["dominance_selected_gene_counts"] = df.astype(np.int32)
    if os.path.exists(paths["dominance_selected_tx_counts_path"]):
        df = pd.read_csv(paths["dominance_selected_tx_counts_path"], sep='\t', index_col=0)
        counts_cols = df.select_dtypes(include=[np.number]).columns
        df[counts_cols] = df[counts_cols].astype(np.int32)
        result["dominance_selected_tx_counts"] = df

    return result


def cluster(
    labels: Union[str, pd.DataFrame],
    output_dir: str = os.getcwd(),
    *,
    include_shared_dtu: bool = False,
    color_palette: str = "Blues",
    color_covariate: Optional[str] = None,
) -> Optional[Any]:
    """Run only the clustering/heatmap step on existing SPIT outputs.

    - Requires that `spit_cluster_matrix.txt` (and optionally `controlled_spit_cluster_matrix.txt`)
      exist under `<output_dir>/SPIT_analysis` from a prior dtu() run.
    - `labels` can be a path or DataFrame; DataFrame will be written under output_dir.
    - Returns the Seaborn ClusterGrid or None if clustering is skipped.
    """

    analysis_dir = os.path.join(output_dir, "SPIT_analysis")
    # Write labels if given as DataFrame so that the plotting function can load it
    labels_path = _write_if_dataframe(labels, os.path.join(analysis_dir, "pheno.txt"), write_index=False)

    cov = color_covariate if color_covariate else "FALSE"
    return _perform_hclust(
        pheno_file=labels_path,
        include_shared_dtu=include_shared_dtu,
        col_palette=color_palette,
        cov=cov,
        output_dir=analysis_dir,
    )


# Convenience viewers for notebooks
def confounding_plot_of(tx_id: str, output_dir: str = os.getcwd()) -> Optional[str]:
    """Display the confounding analysis plot for the given transcript if available.

    Returns the file path, or None if not found.
    """
    plot_path = os.path.join(output_dir, "SPIT_analysis", "confounding_analysis_plots", f"{tx_id}.png")
    if os.path.exists(plot_path):
        if _HAS_IPYTHON:
            try:
                display(IPImage(filename=plot_path))
            except Exception:
                pass
        return plot_path
    # Fallback: try PDF variant
    pdf_path = os.path.join(output_dir, "SPIT_analysis", "confounding_analysis_plots", f"{tx_id}.pdf")
    if os.path.exists(pdf_path):
        if _HAS_IPYTHON:
            try:
                display(IPImage(filename=pdf_path))
            except Exception:
                pass
        return pdf_path
    print(f"Confounding plot not found for {tx_id} under {os.path.dirname(plot_path)}")
    return None


def violin_plot_of(tx_id: str, output_dir: str = os.getcwd()) -> Optional[str]:
    """Display the violin plot for the given transcript if available.

    Returns the file path, or None if not found.
    """
    plot_path = os.path.join(output_dir, "SPIT_analysis", "violin_plots", f"{tx_id}.png")
    if os.path.exists(plot_path):
        if _HAS_IPYTHON:
            try:
                display(IPImage(filename=plot_path))
            except Exception:
                pass
        return plot_path
    # Fallback: try PDF variant
    pdf_path = os.path.join(output_dir, "SPIT_analysis", "violin_plots", f"{tx_id}.pdf")
    if os.path.exists(pdf_path):
        if _HAS_IPYTHON:
            try:
                display(IPImage(filename=pdf_path))
            except Exception:
                pass
        return pdf_path
    print(f"Violin plot not found for {tx_id} under {os.path.dirname(plot_path)}")
    return None

__all__ = [
    "preprocess",
    "dtu",
    "cluster",
    "confounding_plot_of",
    "violin_plot_of",
]


