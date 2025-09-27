#!/usr/bin/env python3

"""
Description:
    Parallel experiment simulation using Python multiprocessing.
Usage:
    ./simulate_exps_python.py -i <tx_counts> -m <tx2gene> -l <pheno> -n <num_experiments>
Author:
    Beril Erdogdu
Date:
    September 2025
"""

import os
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial
from spit.filter_and_simulate_tx_counts import main as fit_param_filter_and_partition
from spit.simulate_dtu import main as simulate_dtu_main


def run_single_experiment(exp_num, args):
    print(f"Simulating Experiment No {exp_num}")
    exp_dir = os.path.join(args.O, "SPIT_analysis", "parameter_fitting", f'exp{exp_num}')
    os.makedirs(exp_dir, exist_ok=True)
    filt_args = argparse.Namespace(
        i=args.i,
        m=args.m,
        l=args.l,
        T=os.path.join(exp_dir, 'filtered_tx_counts.txt'),
        F=os.path.join(exp_dir, 'filtered_ifs.txt'),
        G=os.path.join(exp_dir, 'filtered_gene_counts.txt'),
        case=os.path.join(exp_dir, 'case_samples.txt'),
        ctrl=os.path.join(exp_dir, 'ctrl_samples.txt'),
        n_small=getattr(args, 'n_small', 12),
        pr_fraction=getattr(args, 'pr_fraction', 0.2),
        if_fraction=getattr(args, 'if_fraction', 0.1),
        genefilter_count=getattr(args, 'genefilter_count', 10),
        genefilter_sample=getattr(args, 'genefilter_sample', 10),
        write=getattr(args, 'write', False),
        log_file=os.path.join(exp_dir, 'filter_log.txt'),
    )
    fit_param_filter_and_partition(filt_args)

    sim_args = argparse.Namespace(
        i=os.path.join(exp_dir, 'filtered_ifs.txt'),
        g=os.path.join(exp_dir, 'filtered_gene_counts.txt'),
        case=os.path.join(exp_dir, 'case_samples.txt'),
        ctrl=os.path.join(exp_dir, 'ctrl_samples.txt'),
        F=os.path.join(exp_dir, 'ifs.txt'),
        T=os.path.join(exp_dir, 'tx_counts.txt'),
        G=os.path.join(exp_dir, 'gene_counts.txt'),
        true_cluster=os.path.join(exp_dir, 'true_cluster_array.txt'),
        sim_pheno=os.path.join(exp_dir, 'simulation_pheno.txt'),
        n_small=getattr(args, 'n_small', None),
        n_splicotypes=getattr(args, 'n_splicotypes', getattr(args, 'n_spliceotypes', None)),
        p_dom=getattr(args, 'p_dom', 0.75),
    )
    simulate_dtu_main(sim_args)

    print(f"Experiment {exp_num} completed successfully")
    return f"exp{exp_num}"


def main(args):
    num_experiments = getattr(args, 'n_exps', getattr(args, 'num_experiments', 10))
    jobs = getattr(args, 'jobs', getattr(args, 'threads', None))
    num_jobs = jobs if jobs else min(cpu_count(), num_experiments)

    print(f"Running {num_experiments} experiments using {num_jobs} parallel jobs")

    run_exp = partial(run_single_experiment, args=args)

    with Pool(processes=num_jobs) as pool:
        results = pool.map(run_exp, range(1, num_experiments + 1))
    
    successful = [r for r in results if r is not None]
    failed = len(results) - len(successful)
    
    print(f"\nSimulation complete!")
    print(f"Successful experiments: {len(successful)}")
    print(f"Failed experiments: {failed}")
    
    if successful:
        print(f"Successful experiments: {', '.join(successful)}")
