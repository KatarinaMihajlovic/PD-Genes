# PD-Genes 
PD - Parkinson's disease

This repository has been created to present the work that includes: 
1) NetSC-NMTF - NMTF-based integration model that fuses four molecular interaction networks (protein-protein interaction (PPI), gene co-expression (COEX), metabolic interaction (MI), and genetic interaction (GI) networks) with a scRNA-seq dataset of a PD or control cell line at a specific time point (i.e., cell condition) to obtain gene embedding spaces of individual cell conditions;
2) a 3-step procedure that mines gene embedding spaces of the cell conditions obtained with NetSC-NMTF to identify PD-related gene predictions;
3) methods for analyzing the PD-relevance of the gene predictions.

### What can you find in this repository?
This repository contains all data, scripts and results related to our recent work. 
In particular, you will find:
- the jupyter notebook [PD_NMTF_uniqueEs_AUTOMATED.ipynb](https://github.com/KatarinaMihajlovic/PD-Genes/blob/main/PD_NMTF_uniqueEs_AUTOMATED.ipynb) contains the complete pipeline for: intergating molecular interaction networks with individual expression matrices, mining the resulting gene embedding spaces to get PD-related gene predictions, and computational methods for analyzing the PD-relevance of the gene predictions.
- one folder [CommonData](https://github.com/KatarinaMihajlovic/PD-Genes/tree/main/CommonData) that contains commonly used data in several steps of the pipeline
- folders for each step of the analysis (Step1 - Step11; [auxStep_JustifyClusteringDistanceApproach](https://github.com/KatarinaMihajlovic/PD-Genes/tree/main/auxStep_JustifyClusteringDistanceApproach)) from the jupyter notebook

- one folder [auxStep_Bestk1k2](https://github.com/KatarinaMihajlovic/PD-Genes/tree/main/auxStep_Bestk1k2) that contains folders:
  - [step1_ComputeDispersionCoefficient](https://github.com/KatarinaMihajlovic/PD-Genes/tree/main/auxStep_Bestk1k2/step1_ComputeDispersionCoefficient) that contains the code on how to compute the dispersion coefficients for multiple combinations of cluster numbers (k<sub>1</sub> and k<sub>2</sub>) for individual cell conditions, and the corresponding results 
  - [step2_IdentifyBestk1k2](https://github.com/KatarinaMihajlovic/PD-Genes/tree/main/auxStep_Bestk1k2/step2_IdentifyBestk1k2) that contains the code which visualizes 
the dispersion coefficients of different k<sub>1</sub> and k<sub>2</sub> values, allowing us to choose k<sub>1</sub> and k<sub>2</sub> that lead to high dispersion coefficients


### How to run the notebook
pip install -r requirements.txt

Execute the jupyter notebook [PD_NMTF_uniqueEs_AUTOMATED.ipynb](https://github.com/KatarinaMihajlovic/PD-Genes/blob/main/PD_NMTF_uniqueEs_AUTOMATED.ipynb) 

