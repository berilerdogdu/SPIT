import os
import argparse
import shutil
from spit.transform_tx_counts_to_ifs import main as transform_counts
from spit.filter_and_transform_tx_counts import main as filter_counts
from spit.spit_test import main as spit_test
from spit.get_p_cutoff import main as get_p
from spit.dtu_detection import main as detect_dtu
from spit.confounding_analysis import main as control_confounding
from spit.cluster_samples import call_hclust as cluster_samples
from spit.plot_spit_charts import main as spit_plots
from spit.import_infreps import call_import_infreps as import_infreps
from spit.parameter_fitting.filter_and_simulate_tx_counts import main as fit_param_filter_and_partition
from spit.parameter_fitting.simulate_dtu_exp import main as simulate_exp
from spit.parameter_fitting.LOOCV import main as loocv
import pandas as pd
import numpy as np


def handle_shared_args(parser):
    parser._optionals.title = 'Command-line arguments:'
    parser.add_argument('-i', metavar='tx_counts.txt', type=str, required=True, help='Transcript level counts file (tsv)')
    parser.add_argument('-m', metavar='tx2gene.txt', type=str, required=True, help='Transcript to gene mapping file (tsv)')
    parser.add_argument('-l', metavar='labels.txt', type=str, required=True, help='Labels/metadata file (tsv)')
    parser.add_argument('--n_small', metavar='12', type=int, default=12, help='Smallest sample size for the subgroups')
    parser.add_argument('-O', metavar='/path/', type=str, default=os.getcwd(), help="Output directory path where the SPIT output folder will be written")
    return
    
def handle_filter(args):
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
        if((args.quant_path==None) or (args.quant_type==None)):
            print("Parameters --quant_path and --quant_type are required in order to use infReps.")
            exit(1)
        print("Importing inferential replicates using tximport:")
        import_infreps(args)
    args.exp = False
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
    fit_param_directory_path = os.path.join(args.O, "SPIT_analysis", "parameter_fitting")
    if os.path.exists(fit_param_directory_path) == False:
        os.mkdir(fit_param_directory_path)
    initial_ifs = args.i
    initial_pheno = args.l
    for i in range(1, args.n_exps+1):
        print("--Simulating Experiment No " + str(i))
        args.exp = os.path.join(fit_param_directory_path, 'exp'+str(i))
        os.mkdir(args.exp)
        args.i = initial_ifs
        args.l = initial_pheno
        fit_param_filter_and_partition(args)
        args.i = os.path.join(args.exp, "filtered_ifs.txt")
        args.g = os.path.join(args.exp, "filtered_gene_counts.txt")
        simulate_exp(args)
        print("Running SPIT-Test on Experiment No " + str(i))
        args.i = os.path.join(args.exp, "simulated_ifs.txt")
        args.g = os.path.join(args.exp, "simulated_gene_counts.txt")
        args.l = os.path.join(args.exp, "simulation_pheno.txt")
        args.infReps = False
        spit_test(args)
        for b in np.arange(0.02, 0.21, 0.02):
            b = round(b, 2)
            for k in np.arange(0.1, 1, 0.1):
                k = round(k, 1)
                print("Detecting DTU on Experiment No " + str(i) + " with bandwidth = " + str(b) + " and k= " + str(k))
                args.p = os.path.join(args.exp, "spit_test_min_p_values.txt")
                args.bandwidth = b
                args.k = k
                p_cutoff = round(get_p(args), 16)
                print("Determined p-value threshold: ", p_cutoff)
                args.p_cutoff = p_cutoff
                detect_dtu(args)
    loocv(args)
    
def make_output_dir(target_dir, pheno_path):
    dir_name = os.path.join(target_dir, "SPIT_analysis")
    if os.path.exists(dir_name) == False:
        os.mkdir(dir_name)
    return


def main(argv=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(help='Available SPIT functions')
    # Filter subcommand
    parser_filter = subparsers.add_parser('filter', help='Apply pre-filtering on transcript counts. For specific parameters, run "spit filter -h".')
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
    parser_dtu.add_argument('-b', '--bandwidth', metavar='0.09', type=float, default=0.09, help='choice of bandwidth for kernel density estimation')
    parser_dtu.add_argument('--f_cpm', action='store_true', help='Apply filter-CPM thresholding')
    parser_dtu.add_argument('--plot', action='store_true', help='Plot permutation importances')
    parser_dtu.set_defaults(func=handle_dtu)
    
    # Clustering subcommand
    parser_cluster = subparsers.add_parser('cluster', help='Cluster case samples based on detected DTU events. For specific parameters, run "spit cluster -h".')
    parser_cluster.add_argument('-l', metavar='labels.txt', type=str, required=True, help='Labels/metadata file (tsv)')
    parser_cluster.add_argument('-O', metavar='/path/', type=str, default=os.getcwd(), help="Output directory path where the SPIT output folder will be written")
    parser_cluster.add_argument('--include_shared_dtu', action='store_true', help='Include transcripts that are DTU in all case samples in the heatmap')
    parser_cluster.add_argument('--color_palette', metavar = 'BuGn', default='BuGn', help='The RColorBrewer palette to be used for heatmap. See list of available palettes here: https://r-graph-gallery.com/38-rcolorbrewers-palettes.html')
    parser_cluster.add_argument('--color_covariate', metavar = 'batch', default = False, help='Color samples in the heatmap based on this covariate')
    parser_cluster.set_defaults(func=handle_clustering)
    
    # Parameter-fitting subcommand
    parser_fit = subparsers.add_parser('fit_parameters', help='Apply parameter-fitting on your dataset (Optional) For specific parameters, run "spit fit_parameters -h".')
    handle_shared_args(parser_fit)
    parser_fit.add_argument('--n_exps', metavar='10', type=int, default=10, help='Number of experiments to simulate.')
    parser_fit.add_argument('--n_splicotypes', metavar='5', type=int, default=5, help='Number of splicotypes (subgroups) to simulate.')
    parser_fit.add_argument('-w', '--write', action='store_true', help='Write the number of transcripts & genes left after each filtering step to stdout.')
    parser_fit.add_argument('--keep_all_nonzeros', action='store_true', help='If used, this options skips all SPIT prefiltering steps and only removes transcripts that do not have any non-zero counts in any sample. Any other filtering argument becomes irrelevant.')
    parser_fit.add_argument('-p', '--pr_fraction', metavar='0.2', type=float, default=0.2, help='Each transcript must have a positive read count in at least a fraction p_r of the samples in both the case and control groups.')
    parser_fit.add_argument('-f', '--if_fraction', metavar='0.1', type=float, default=0.1, help='Each transcript must have an IF value larger than f in at least n_small samples.')
    parser_fit.add_argument('-c', '--genefilter_count', metavar='10', type=int, default=10, help='Each gene must have a read count of at least c in at least s samples.')
    parser_fit.add_argument('-s', '--genefilter_sample', metavar='10', type=int, default=10, help='Each gene must have a read count of at least c in at least s samples.')
    parser_fit.add_argument('-d', '--p_dom', metavar='0.75', type=float, default=0.75, help='Dominance selection threshold')
    parser_fit.add_argument('--n_iter', metavar='100', type=int, default=100, help='Number of iterations')
    parser_fit.add_argument('--f_cpm', action='store_true', help='Apply filter-CPM thresholding')
    parser_fit.set_defaults(func=handle_param_fit)
    
    

    args = parser.parse_args()
    make_output_dir(args.O, args.l)
    args.func(args)
    if not any(vars(args).values()):
        parser.print_help()
        exit(1)


if __name__ == "__main__":
   main(sys.argv[1:])
