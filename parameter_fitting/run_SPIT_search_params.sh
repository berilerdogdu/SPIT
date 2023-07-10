#!/bin/bash

run_spit() {
    local i=$1
    local tx2gene_file=$2

    python ../spit_test.py -i ./exp$i/ifs.txt -g ./exp$i/filtered_gene_counts.txt -l ./exp$i/simulation_pheno.txt -n 1000 --p_dom 0 -P ./exp$i/spit_test_p_values.txt

    wait

    for b in $(seq 0.02 .01 0.20)
    do
        for k in $(seq 0.1 .1 1.0)
        do
                echo "Running SPIT on experiment" $i "with bandwidth=" $b "and k=" $k

                P_CUTOFF=$(python ../get_p_cutoff.py -k $k -p ./exp$i/spit_test_p_values.txt)

                python ../dtu_detection.py -i ./exp$i/ifs.txt -g ./exp$i/gene_counts.txt -m $tx2gene_file -l ./exp$i/simulation_pheno.txt --p_cutoff $P_CUTOFF -b $b -M ./exp$i/spit_cluster_array_$k.$b.txt -O ./exp$i/spit_out_$k.$b.txt

        done

    done
}

export -f run_spit

if [ $# -lt 3 ]; then
    echo "Usage: $0 <num_threads> -m <tx2gene_file>"
    exit 1
fi

threads=$1
shift

while getopts "m:" opt; do
  case $opt in
    m)
      tx2gene_file=$OPTARG ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1 ;;
  esac
done

experiments=(1 2 3 4 5 6 7 8 9 10)

parallel --jobs $threads run_spit ::: ${experiments[@]} ::: "$tx2gene_file"

wait
