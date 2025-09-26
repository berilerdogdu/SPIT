#!/usr/bin/env python3

"""
Description:
    Parallel SPIT parameter search using Python multiprocessing.
Usage:
    ./run_SPIT_search_params_python.py -m <tx2gene_file> -j <num_jobs> -e <num_experiments>
Author:
    Beril Erdogdu
Date:
    September 2025
"""

from multiprocessing import Pool, cpu_count
from functools import partial
import numpy as np
import os
import argparse
from spit.spit_test import main as spit_test_main
from spit.get_p_cutoff import main as get_p_cutoff_main
from spit.dtu_detection import main as dtu_detection_main


def run_spit_test(exp_num, args):
    print(f"Running SPIT test for experiment {exp_num}")
    exp_dir = os.path.join(args.O, "SPIT_analysis", "parameter_fitting", f"exp{exp_num}")
    test_args = argparse.Namespace()
    test_args.i = os.path.join(exp_dir, "ifs.txt")
    test_args.g = os.path.join(exp_dir, "filtered_gene_counts.txt")
    test_args.l = os.path.join(exp_dir, "simulation_pheno.txt")
    test_args.n_iter = 100
    test_args.n_small = getattr(args, 'n_small', 12)
    test_args.O = args.O
    test_args.exp = exp_dir
    test_args.quiet = getattr(args, 'quiet', True)
    spit_test_main(test_args)
    print(f"SPIT test completed for experiment {exp_num}")
    return True


def run_parameter_search(exp_num, args):

    print(f"Starting parameter search for experiment {exp_num}")
    if hasattr(args, 'no_clusters') and args.no_clusters:
        bandwidths = np.array([np.round(1.0, 2)])
    else:
        bandwidths = np.round(np.arange(0.02, 0.21, 0.01), 2)
    k_values = np.round(np.arange(0.1, 1.1, 0.1), 1)
    
    successful_runs = 0
    total_runs = len(bandwidths) * len(k_values)
    
    for b in bandwidths:
        for k in k_values:
            try:
                exp_dir = os.path.join(args.O, "SPIT_analysis", "parameter_fitting", f"exp{exp_num}")
                gp_args = argparse.Namespace(**vars(args))
                gp_args.k = float(k)
                gp_args.p = os.path.join(exp_dir, "spit_test_min_p_values.txt")
                p_cutoff = get_p_cutoff_main(gp_args)

                dtu_args = argparse.Namespace(**vars(args))
                dtu_args.i = os.path.join(exp_dir, "ifs.txt")
                dtu_args.g = os.path.join(exp_dir, "gene_counts.txt")
                dtu_args.m = args.m
                dtu_args.l = os.path.join(exp_dir, "simulation_pheno.txt")
                dtu_args.p_cutoff = p_cutoff
                dtu_args.bandwidth = float(b)
                dtu_args.k = float(k)
                dtu_args.n_small = getattr(args, 'n_small', 12)
                dtu_args.f_cpm = getattr(args, 'f_cpm', False)
                dtu_args.infReps = getattr(args, 'infReps', False)
                dtu_args.exp = exp_dir
                dtu_args.quiet = getattr(args, 'quiet', True)
                dtu_detection_main(dtu_args)
                successful_runs += 1
            except Exception as e:
                print(f"Failed experiment {exp_num}, k={k}, b={b}: {e}")
    
    print(f"Parameter search completed for experiment {exp_num}: {successful_runs}/{total_runs} successful")
    return successful_runs, total_runs


def run_single_experiment(exp_num, args):
    print(f"Processing experiment {exp_num}")
    
    if not run_spit_test(exp_num, args):
        return f"exp{exp_num}", False, 0, 0
    
    successful, total = run_parameter_search(exp_num, args)
    
    return f"exp{exp_num}", True, successful, total


def main(args):
    experiments = getattr(args, 'n_exps', getattr(args, 'experiments', 1))
    jobs = getattr(args, 'jobs', getattr(args, 'threads', None))
    num_jobs = jobs if jobs else min(cpu_count(), experiments)
    print(f"Running SPIT parameter search on {experiments} experiments using {num_jobs} parallel jobs")
    run_exp = partial(run_single_experiment, args=args)
    
    with Pool(processes=num_jobs) as pool:
        results = pool.map(run_exp, range(1, experiments + 1))
    
    total_successful = 0
    total_runs = 0
    failed_experiments = []
    
    for exp_name, success, successful_runs, total_runs_exp in results:
        if success:
            total_successful += successful_runs
            total_runs += total_runs_exp
        else:
            failed_experiments.append(exp_name)
    
    print(f"\nSPIT parameter search complete!")
    print(f"Successful experiments: {len(results) - len(failed_experiments)}/{len(results)}")
    print(f"Successful parameter combinations: {total_successful}/{total_runs}")
    
    if failed_experiments:
        print(f"Failed experiments: {', '.join(failed_experiments)}")
