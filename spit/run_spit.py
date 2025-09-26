import sys
import os
import argparse
from spit.transform_tx_counts_to_ifs import main as transform_counts
from spit.filter_and_transform_tx_counts import main as filter_counts
from spit.spit_test import main as spit_test
from spit.get_p_cutoff import main as get_p
from spit.dtu_detection import main as detect_dtu
from spit.confounding_analysis import main as control_confounding
from spit.cluster_samples import call_hclust as cluster_samples
from spit.plot_spit_charts import main as spit_plots
from spit.simulate_exps import main as simulate_exps_main
from spit.run_spit_search_params import main as run_param_search_main
from spit.LOOCV import main as loocv
import pandas as pd


def handle_shared_args(parser):
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-i', metavar='tx_counts.txt', type=str, required=True, help='Transcript level counts file (tsv)')
    parser.add_argument('-m', metavar='tx2gene.txt', type=str, required=True, help='Transcript to gene mapping file (tsv)')
    parser.add_argument('-l', metavar='labels.txt', type=str, required=True, help='Labels/metadata file (tsv)')
    parser.add_argument('--no_clusters', action='store_true', help='Assume case samples do not form clusters/subgroups.')
    parser.add_argument('--n_small', metavar='12', type=int, default=12, help='Smallest sample size for the subgroups. If "--no_clusters" option is used, it will be overwritten to 0.75 times the sample size of the smaller group for filtering purposes.')
    parser.add_argument('-O', metavar='/path/', type=str, default=os.getcwd(), help='Output directory path where the SPIT output folder will be written.')
    parser.add_argument('--quiet', action='store_true', default=False, help='No verbose logging. Only show warnings and errors.')
    return

def set_no_cluster_nsmall(args, param_bool):
    pheno = pd.read_csv(args.l, sep = '\t')
    case_sample_size = (pheno.condition == 1).sum()
    ctrl_sample_size = (pheno.condition == 0).sum()
    if(param_bool):
        n_small_ = int((ctrl_sample_size/2)*0.75)
    else:
        n_small_ = min(int(case_sample_size*0.75), int(ctrl_sample_size*0.75))
    return n_small_
    
def handle_filter(args):
    if(args.no_clusters):
        args.n_small = set_no_cluster_nsmall(args, False)
    if(args.keep_all_nonzeros):
        transform_counts(args)
    else:
        filter_counts(args)
    
def handle_dtu(args):
    dominance_filtered_file = os.path.join(args.O, "SPIT_analysis", "dominance_selected_ifs.txt")
    if os.path.exists(dominance_filtered_file):
        args.i = os.path.join(args.O, "SPIT_analysis", "dominance_selected_ifs.txt")
        args.g = os.path.join(args.O, "SPIT_analysis", "dominance_selected_gene_counts.txt")
    else:
        args.i = os.path.join(args.O, "SPIT_analysis", "filtered_ifs.txt")
        args.g = os.path.join(args.O, "SPIT_analysis", "filtered_gene_counts.txt")
    if(args.infReps):
        print("Warning: Inferential replicates should be imported and processed using the independent R script provided in the documentation.")
    args.exp = False
    if(args.no_clusters):
        args.n_small = set_no_cluster_nsmall(args, False)
        args.bandwidth = 1
    spit_test(args)
    args.p = os.path.join(args.O, "SPIT_analysis", "spit_test_min_p_values.txt")
    p_cutoff = round(get_p(args), 16)
    args.p_cutoff = p_cutoff
    print("Determined p-value threshold: ", p_cutoff)
    if(p_cutoff > 0.05):
        print("WARNING: Determined p-value threshold is larger than 0.05. In this case we recommend using a simple Mann-Whitney U test instead of using SPIT.")
    detect_dtu(args)
    if(check_for_covariates(args.l)):
        control_confounding(args)
    if(args.plot):
        spit_plots(args)
        
def check_for_covariates(label_file):
    pheno_df = pd.read_csv(label_file, sep = '\t')
    labels = pheno_df.columns.to_list()
    cov_labels = [label for label in labels if label not in ["id", "condition"]]
    if cov_labels:
        print("Running confounding effect control on the following covariates: ", cov_labels)
        return True
    return False

def handle_clustering(args):
    cluster_samples(args)
    
def handle_param_fit(args):
    args.quiet = True
    fit_param_directory_path = os.path.join(args.O, "SPIT_analysis", "parameter_fitting")
    if os.path.exists(fit_param_directory_path) == False:
        os.mkdir(fit_param_directory_path)
    if(args.no_clusters):
        args.n_splicotypes = 1
        args.n_small = set_no_cluster_nsmall(args, True)  
    simulate_exps_main(args)
    run_param_search_main(args)
    loocv(args)
    
def make_output_dir(target_dir):
    dir_name = os.path.join(target_dir, "SPIT_analysis")
    if os.path.exists(dir_name) == False:
        os.mkdir(dir_name)
    return


def main(argv=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(help='Available SPIT functions')
    # Filter subcommand
    parser_filter = subparsers.add_parser('preprocess', help='Apply pre-filtering on transcript counts and generate isoform fractions. For specific parameters, run "spit preprocess -h".')
    handle_shared_args(parser_filter)
    parser_filter.add_argument('-w', '--write', action='store_true', help='Write the number of transcripts & genes left after each filtering step to stdout.')
    parser_filter.add_argument('--keep_all_nonzeros', action='store_true', help='If used, this options skips all SPIT prefiltering steps and only removes transcripts that do not have any non-zero counts in any sample. Any other filtering argument becomes irrelevant.')
    parser_filter.add_argument('-p', '--pr_fraction', metavar='0.2', type=float, default=0.2, help='Each transcript must have a positive read count in at least a fraction p_r of the samples in both the case and control groups.')
    parser_filter.add_argument('-f', '--if_fraction', metavar='0.1', type=float, default=0.1, help='Each transcript must have an IF value larger than f in at least n_small samples.')
    parser_filter.add_argument('-c', '--genefilter_count', metavar='10', type=int, default=10, help='Each gene must have a read count of at least c in at least s samples.')
    parser_filter.add_argument('-s', '--genefilter_sample', metavar='10', type=int, default=10, help='Each gene must have a read count of at least c in at least s samples.')
    parser_filter.add_argument('-d', '--p_dom', metavar='0.75', type=float, default=0.75, help='Dominance selection threshold')
    parser_filter.set_defaults(func=handle_filter)
    
    # DTU subcommand
    parser_dtu = subparsers.add_parser('dtu', help='Run DTU analysis on transcript abundances. For specific parameters, run "spit dtu -h".')
    handle_shared_args(parser_dtu)
    parser_dtu.add_argument('--n_iter', metavar='100', type=int, default=100, help='Number of iterations')
    parser_dtu.add_argument('-k', metavar='0.6', type=float, default=0.6, help='-K hyperparameter for p-value thresholding')
    parser_dtu.add_argument('--infReps', action='store_true', default=None, help='Use inferential replicates in DTU analysis. Parameters "--quant_path" and "--quant_type" must be specified.')
    parser_dtu.add_argument('--quant_path', metavar='/path/', type=str, default=None, help='Path to the parent output directory of Salmon, kallisto, or other quantification tool compatible with tximport (the directory will contain a folder for each sample).')
    parser_dtu.add_argument('--quant_type', metavar='salmon', type=str, default=None, help='The quantification tool used to generate counts. Options compatible with tximport are "salmon", "sailfish", "alevin", "kallisto", "rsem", and "stringtie".')
    parser_dtu.add_argument('-b', '--bandwidth', metavar='0.09', type=float, default=0.09, help='Choice of bandwidth for kernel density estimation. It will be overwritten and set to 1 if "--no_clusters" option is used.')
    parser_dtu.add_argument('--f_cpm', action='store_true', help='Apply filtered-CPM thresholding')
    parser_dtu.add_argument('--plot', action='store_true', help='Plot permutation importances')
    parser_dtu.set_defaults(func=handle_dtu)
    
    # Clustering subcommand
    parser_cluster = subparsers.add_parser('cluster', help='Cluster case samples based on detected DTU events. For specific parameters, run "spit cluster -h".')
    parser_cluster.add_argument('-l', metavar='labels.txt', type=str, required=True, help='Labels/metadata file (tsv)')
    parser_cluster.add_argument('-O', metavar='/path/', type=str, default=os.getcwd(), help='Output directory path where the SPIT output folder will be written')
    parser_cluster.add_argument('--include_shared_dtu', action='store_true', help='Include transcripts that are DTU in all case samples in the heatmap')
    parser_cluster.add_argument('--color_palette', metavar = 'Blues', default='Blues', help='The Seaborn palette to be used for heatmap. See list of available palettes here: https://seaborn.pydata.org/tutorial/color_palettes.html')
    parser_cluster.add_argument('--color_covariate', metavar = 'batch', default = False, help='Color samples in the heatmap based on this covariate')
    parser_cluster.set_defaults(func=handle_clustering)
    
    # Parameter-fitting subcommand
    parser_fit = subparsers.add_parser('fit_parameters', help='Apply parameter-fitting on your dataset (Optional) For specific parameters, run "spit fit_parameters -h".')
    handle_shared_args(parser_fit)
    parser_fit.add_argument('--n_exps', metavar='10', type=int, default=10, help='Number of experiments to simulate.')
    parser_fit.add_argument('--n_spliceotypes', dest='n_splicotypes', metavar='5', type=int, default=5, help='Number of subgroups to simulate.')
    parser_fit.add_argument('--n_dtu_genes', metavar='30', type=int, default=30, help='Number of DTU genes to simulate per subgroup.')
    parser_fit.add_argument('-w', '--write', action='store_true', help='Write the number of transcripts & genes left after each filtering step to stdout.')
    parser_fit.add_argument('--keep_all_nonzeros', action='store_true', help='If used, this options skips all SPIT prefiltering steps and only removes transcripts that do not have any non-zero counts in any sample. Any other filtering argument becomes irrelevant.')
    parser_fit.add_argument('-p', '--pr_fraction', metavar='0.2', type=float, default=0.2, help='Each transcript must have a positive read count in at least a fraction p_r of the samples in both the case and control groups.')
    parser_fit.add_argument('-f', '--if_fraction', metavar='0.1', type=float, default=0.1, help='Each transcript must have an IF value larger than f in at least n_small samples.')
    parser_fit.add_argument('-c', '--genefilter_count', metavar='10', type=int, default=10, help='Each gene must have a read count of at least c in at least s samples.')
    parser_fit.add_argument('-s', '--genefilter_sample', metavar='10', type=int, default=10, help='Each gene must have a read count of at least c in at least s samples.')
    parser_fit.add_argument('-d', '--p_dom', metavar='0.75', type=float, default=0.75, help='Dominance selection threshold')
    parser_fit.add_argument('--n_iter', metavar='100', type=int, default=100, help='Number of iterations')
    parser_fit.add_argument('--f_cpm', action='store_true', help='Apply filtered-CPM thresholding. This option might lower the F-scores of experiments.')
    parser_fit.add_argument('--threads', metavar='auto', type=int, default=None, help='Number of parallel workers to use (default is auto-detect based on number of available cores)')
    parser_fit.set_defaults(func=handle_param_fit)
    

    args = parser.parse_args()
    if not any(vars(args).values()):
        parser.print_help()
        exit(1)
    make_output_dir(args.O)
    args.func(args)


if __name__ == "__main__":
   main(sys.argv[1:])
