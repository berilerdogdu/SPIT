<img src="https://raw.githubusercontent.com/berilerdogdu/SPIT/main/logos/spitting_llama.png" alt="spitting_llama" style="width:376px; height:200px; float:left;"> 


# SPIT 

A statistical tool that quantifies the heterogeneity in transcript usage within a population and identifies predominant subgroups along with their distinctive sets of DTU events. 


### Why use SPIT?

Detecting DTU events for single-gene genetic traits is relatively uncomplicated; however, the heterogeneity of populations with complex diseases presents an intricate challenge due to the presence of diverse causal events and undetermined subtypes.
SPIT can detect DTU events exclusive to subgroups as well as DTU events shared amongst all case samples. Downstream of DTU analysis, SPIT uses detected DTU events to provide insight into potentially hierarchical subgrouping patterns present in complex disease populations using hierarchical clustering.

SPIT is equally effective on relatively homogeneous populations, and proves to be applicable for diverse scenarios, including simple genetic disorders, tissue-to-tissue comparisons and other types of DTU studies. SPIT consistently maintains notably low false discovery rates regardless of the level of dispersion in the datasets.

### How to use SPIT?

SPIT is available as a PyPI package and can be installed by calling:
```
pip install spit
```

An extensive step-by-step guide that demonstrates the application of SPIT using a mock dataset is provided [here](https://colab.research.google.com/drive/1u3NpleqcAfNz_0EAgO2UHItozd9PsF1w?usp=sharing).

Users can also directly upload their datasets into this Colab environment and easily run SPIT [online](https://colab.research.google.com/drive/1u3NpleqcAfNz_0EAgO2UHItozd9PsF1w?usp=sharing).


##### If you use SPIT, please cite:
Erdogdu, B., Varabyou, A., Hicks, S.C., Salzberg, S.L. & Pertea, M. Detecting differential transcript usage in complex diseases with SPIT. bioRxiv, 2023.2007.2010.548289 (2023)

##### Please use [this Google group](https://groups.google.com/g/spit_dtu) to post your questions, comments, or bug reports.
