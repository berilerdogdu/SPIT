<img src="https://raw.githubusercontent.com/berilerdogdu/SPIT/main/spitting_llama.png" alt="spitting_llama" style="width:376px; height:200px; float:left;"> 


# SPIT 

A statistical tool that quantifies the heterogeneity in transcript usage within a population and identifies predominant subgroups along with their distinctive sets of DTU events. 


### Why use SPIT?

Detecting DTU events for single-gene genetic traits is relatively uncomplicated; however, the heterogeneity of populations with complex diseases presents an intricate challenge due to the presence of diverse causal events and undetermined subtypes.
SPIT can detect DTU events exclusive to subgroups as well as DTU events shared amongst all case samples. Downstream of DTU analysis, SPIT uses detected DTU events to provide insight into potentially hierarchical subgrouping patterns present in complex disease populations using hierarchical clustering.

SPIT is equally effective on relatively homogeneous populations, and proves to be applicable for diverse scenarios, including simple genetic disorders, tissue-to-tissue comparisons and other types of DTU studies. SPIT consistently maintains notably low false discovery rates regardless of the level of dispersion in the datasets.

### How to use SPIT?

An extensive step-by-step guide that demonstrates the application of SPIT using a mock dataset is provided [here](https://colab.research.google.com/drive/1u3NpleqcAfNz_0EAgO2UHItozd9PsF1w?usp=sharing).

Users can also directly upload their datasets into this Colab environment and easily run SPIT [online](https://colab.research.google.com/drive/1u3NpleqcAfNz_0EAgO2UHItozd9PsF1w?usp=sharing).

### Parameter-fitting

If you wish to run the parameter-fitting module, please clone this project and follow these steps:
- Navigate into the "parameter_fitting" directory, and generate your 10 DTU simulations by running:
```
sh simulate_exps.sh -i [tx_counts_file] -m [tx2gene_file] -l [pheno_file]
```
Please note that your input files should follow the formatting requirements described in the Colab notebook. For detailed explanations of what these files should contain, please refer to the [tutorial](https://colab.research.google.com/drive/1u3NpleqcAfNz_0EAgO2UHItozd9PsF1w?usp=sharing).

- Run SPIT with combinations of **b** and **k** parameters to search for the optimal choice for your dataset:
```
sh run_SPIT_search_params.sh [#number of threads] -m [tx2gene_file]
```
If you have multiple threads available, the run is distributed accordingly via GNU Parallel. The maximum number of threads that can be used is equal to the number of simulated experiments (10). For example, if you would like to run with 10 threads, run:

```
sh run_SPIT_search_params.sh 10 -m [tx2gene_file]
```
Otherwise, please run as:
```
sh run_SPIT_search_params.sh 1 -m [tx2gene_file]
```
- Run the leave-one-out cross-validation (LOOCV) step to see the optimal parameters in all 10 experiments as:
```
python LOOCV.py -m tx2gene_file.txt -P venns.pdf
```
The output PDF file (venns.pdf) will include the optimal parameters at each iteratioin of the LOOCV process along with corresponding true positive and false positive rates and *F*-scores. 

##### If you use SPIT, please cite:
Erdogdu, B., Varabyou, A., Hicks, S.C., Salzberg, S.L. & Pertea, M. Detecting differential transcript usage in complex diseases with SPIT. bioRxiv, 2023.2007.2010.548289 (2023)

##### Please use [this Google group](https://groups.google.com/g/spit_dtu) to post your questions, comments, or bug reports.
