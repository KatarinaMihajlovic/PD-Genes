# PD-Genes 
PD - Parkinson's disease

This repository has been created to present the work that includes: 
1) NetSC-NMTF - NMTF-based integration model that fuses four molecular interaction networks (protein-protein interaction (PPI), gene co-expression (COEX), metabolic interaction (MI), and genetic interaction (GI) networks) with a scRNA-seq dataset of a PD or control cell line at a specific time point (i.e., cell condition) to obtain gene embedding spaces of individual cell conditions;
2) a 3-step procedure that mines gene embedding spaces of the cell conditions obtained with NetSC-NMTF to identify PD-related gene predictions;
3) methods for analyzing the PD-relevance of the gene predictions.

The jupyter notebook [PD_NMTF_uniqueEs_AUTOMATED.ipynb](https://github.com/KatarinaMihajlovic/PD-Genes) contains the complete pipeline for: intergating molecular interaction networks with individual expression matrices, mining the resulting gene embedding spaces to get PD-related gene predictions, and computational methods for analyzing the PD-relevance of the gene predictions.


### What can you find in this repository?
This repository contains all data, scripts and results related to the LIHC tumor analysis.
In particular, you will find:
- [preprocessing.py](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/preprocessing.py), [creation_matrices.py](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/creation_matrices.py), [fusion_networks.py](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/fusion_networks.py), [extraction_classification.py](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/extraction_classification.py), python files which contain the main steps of the process,
- two folders [Datasets](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Datasets) and [Datasets_created](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Datasets_created) that contain the TCGA gene expression data retrieved through [GMQL](http://gmql.eu) and [PyGMQL](https://pygmql.readthedocs.io/en/latest/),
- one folder [Graphs](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Graphs) that contains the LIHC fused network in the format required by [Gephi](https://gephi.org),
- one folder [Extracted](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Extracted) that contains the data regarding all the genes extracted with the pipeline,
- one folder [Auc_acc_f1](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Auc_acc_f1) that contains the tables with the performances values after having executed [extraction_classification.py](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/extraction_classification.py),
- one folder [Boxplots](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Boxplots) that contains the boxplots with the comparisons of the performances of the classification,
- one file [matrix_GTypeGName.xlsx](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/matrix_GTypeGName.xls) that contains the correspondences between the gene symbols and the gene type,
- one folder [Results](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Results) containing the gene symbols of the LIHC fused network, the IC gene symbols for LIHC, KIRC and PRAD, and the PubMed evaluation results for LIHC presented in [Knowledge-based evaluation of the results.ipynb](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/Knowledge-based%20evaluation%20of%20the%20results.ipynb),
- one folder [Supplementary data](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Supplementary%20data) containing the [DB_pharmacologically_active.csv](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/Supplementary%20data/DB_pharmacologically_active.csv) file downloaded from DrugBank to be used in the DrugBank evaluation part of [Knowledge-based evaluation of the results.ipynb](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/Knowledge-based%20evaluation%20of%20the%20results.ipynb). 

### How to run the notebook
pip install -r requirements.txt

Execute the jupyter notebook [Pipeline.ipynb](https://github.com/DEIB-GECO/GeneNetFusion/blob/master/Pipeline.ipynb) until the part 'After extraction of communities with Gephi'.
Use [Gephi](https://gephi.org) in order to extract the relevant communities and save the genes in the folder [Extracted](https://github.com/DEIB-GECO/GeneNetFusion/tree/master/Extracted) with the name 'IC_'+str(tumor)+'.csv'.
