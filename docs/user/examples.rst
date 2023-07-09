.. _examples:

.. role:: img-inline
.. |colab_logo| image:: content/images/colab.png
   :class: img-inline
   :alt: SPIT Colab
   :target: https://colab.research.google.com/github/berilerdogdu/spit/blob/master/notebooks/SPIT.ipynb
.. raw:: html

   <style>
   .img-inline {
       height: 2em;
       vertical-align: middle;
   }
   </style>

Examples
======================

In this section we will demonstrate how to use **SPIT** to identify DTU events using a test dataset.
While the data in this example is provided with the software and is for demonstration purposes only,\
the same steps can be repeated using your owndatasets!
**SPIT** is very efficient and can easily handle large datasets.
The full interactive workflow is available in the |colab_logo| `Google Colab <github.com>`__ notebook.

-------------------

In this module, we illustrate the application of pre-filtering, SPIT-Test, DTU detection,
and confounding control steps of the SPIT pipeline using an example dataset.
Our demonstration highlights the process of transforming input transcript counts to detect DTU events,
controlling for confounding variables, and exploring potential subclusters in your disease group.
Please note that if you wish to employ the optional cross-validation based parameter-fitting process
on your data, you will need to clone the GitHub project and run it locally,
as there are runtime limitations with Google Colab |colab_logo|.


Pre-filtering
--------------
The default behavior of SPIT involves the conservative pre-filtering steps listed below which build upon the DRIMSeq[2] filtering criteria. We define n as the smallest sample size presumed for the subgroups within populations which is a user-set parameter of SPIT and defaults to 12.

1. Each transcript must have a Counts per million (CPM) value of at least 10 in at least n samples.
2. Each transcript must have a positive read count in at least a r fraction of the samples in both the case and control groups, respectively. r is a user-set parameter and defaults to 0.20.
3. Each gene must have a read count of at least c in at least s samples, where c and s are user-set parameters and default to 10.
4. Each transcript must have an IF value larger than f in at least n samples, where f is a user-set parameter and defaults to 0.1.
5. After the filtering steps above, there must remain at least 2 transcripts for each gene.

The pre-filtering script is called as follows:

.. code-block:: console
   :caption: Pre-filtering transcript counts and generating isoform fractions

   $ filter_and_transform_tx_counts.py -i example_analysis/tx_counts.tsv -m example_analysis/tx2gene.txt -l example_analysis/pheno.txt -T filtered_tx_counts.txt -F filtered_ifs.txt -G filtered_gene_counts.txt --write

If you would like to skip pre-filtering, you may convert your transcript counts into IFs as follows. Please note that transcripts without any non-zero counts will still be filtered out.

.. code-block:: console
    :caption:

    $ transform_tx_counts_to_ifs.py -h



.. code-block:: console
   :caption: Running SPIT Test with 1000 iterations

   $ spit_test.py -i filtered_ifs.txt -g filtered_gene_counts.txt -l example_analysis/pheno.txt -n 1000 -I dominance_selected_ifs.txt -G dominance_selected_gene_counts.txt -P spit_test_min_p_values.txt

.. code-block:: console
   :caption: Get P-value Threshold

   $ get_p_cutoff.py -k 0.06 -p spit_test_min_p_values.txt

.. code-block:: console
   :caption: Comparing the case and control samples with default parameters (k=0.6 and bandwidth=0.09)

   $ dtu_detection.py -i dominance_selected_ifs.txt -g dominance_selected_gene_counts.txt -m example_analysis/tx2gene.txt -l example_analysis/pheno.txt --p_cutoff <p-value threshold from previous step> -b 0.09 -M spit_cluster_matrix.txt -O spit_out.txt

.. code-block:: clonsole
   :caption: Applying confounding analysis based on the provided pheno-data. The resulting permutation importance plots for all candidate DTU transcripts will be outputted into 'importance_score_plots.pdf' in the current directory.

   $ confounding_analysis.py -i dominance_selected_ifs.txt -l example_analysis/pheno.txt --cluster_matrix spit_cluster_matrix.txt -n $n_small -o spit_out.txt -M controlled_spit_cluster_matrix.txt -O controlled_spit_out.txt -P importance_score_plots.pdf

