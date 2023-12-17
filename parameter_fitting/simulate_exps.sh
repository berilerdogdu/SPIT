#!/bin/sh

while getopts "i:m:l:" opt; do
  case $opt in
    i)
      tx_counts_file=$OPTARG ;;
    m)
      tx2gene_file=$OPTARG ;;
    l)
      pheno_file=$OPTARG ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1 ;;
  esac
done

for i in {1..10}; do
  echo "Simulating Experiment No $i"

  exp_dir="./exp$i"
  mkdir "$exp_dir"
  
  python filter_and_simulate_tx_counts.py -i "$tx_counts_file" -m "$tx2gene_file" -l "$pheno_file" -T "$exp_dir/filtered_tx_counts.txt" -F "$exp_dir/filtered_ifs.txt" -G "$exp_dir/filtered_gene_counts.txt" --case "$exp_dir/case_samples.txt" --ctrl "$exp_dir/ctrl_samples.txt" --write

  python simulate_dtu_exp.py -i "$exp_dir/filtered_ifs.txt" -g "$exp_dir/filtered_gene_counts.txt" -F "$exp_dir/ifs.txt" -T "$exp_dir/tx_counts.txt" -G "$exp_dir/gene_counts.txt" --case "$exp_dir/case_samples.txt" --ctrl "$exp_dir/ctrl_samples.txt" --true_cluster "$exp_dir/true_cluster_array.txt" --sim_pheno "$exp_dir/simulation_pheno.txt"
done
