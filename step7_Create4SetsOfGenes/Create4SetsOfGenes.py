# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:57:52 2021

@author: kmihajlo
"""
import os
import pandas as pd
import pickle

EXP_dir = 'input/EXP_matrices'
Gene_sets_dir = 'output/GeneSets'
if not os.path.exists(Gene_sets_dir):
    os.mkdir(Gene_sets_dir)
EXP_files = os.listdir(EXP_dir)

All_genes = pd.read_csv('input/All_genes.csv', header=None, index_col = 0)
All_genes = All_genes.index.tolist()

    
for cell_cond1_file in EXP_files:
    EXP_cellcond1 = pd.read_csv(f'{EXP_dir}/{cell_cond1_file}', index_col=0)
    Genes_cell_cond1 = EXP_cellcond1.index.tolist()
    cell_cond1 = cell_cond1_file.split('.')[0]
    print(cell_cond1)
    Gene_sets = {}
    
    for cell_cond2_file in EXP_files:
        if cell_cond1_file != cell_cond2_file:
            EXP_cellcond2 = pd.read_csv(f'{EXP_dir}/{cell_cond2_file}', index_col=0)
            Genes_cell_cond2 = EXP_cellcond2.index.tolist()
            cell_cond2 = cell_cond2_file.split('.')[0]
            
            Shared_genes = list(set(Genes_cell_cond1) & set(Genes_cell_cond2))
            Genes_cell_cond1_ONLY = [x for x in Genes_cell_cond1 if x not in Shared_genes]
            Genes_cell_cond2_ONLY = [x for x in Genes_cell_cond2 if x not in Shared_genes]
            Gene_union = set(Genes_cell_cond1 + Genes_cell_cond2)
            Non_expressedGenes = [x for x in All_genes if x not in Gene_union]
            
            Gene_sets[cell_cond2] = {}
            Gene_sets[cell_cond2][f'{cell_cond1}_ONLY'] = Genes_cell_cond1_ONLY
            Gene_sets[cell_cond2][f'{cell_cond2}_ONLY'] = Genes_cell_cond2_ONLY
            Gene_sets[cell_cond2][f'{cell_cond1} - {cell_cond2}'] = Shared_genes
            Gene_sets[cell_cond2][f'{cell_cond1} - {cell_cond2} - NonExp'] = Non_expressedGenes
            print(len(Genes_cell_cond1_ONLY), len(Genes_cell_cond2_ONLY), len(Shared_genes), len(Non_expressedGenes))
    
    cell_cond1 = 'GS_' + cell_cond1.split('_')[1] + cell_cond1.split('_')[2]
    with open(f'{Gene_sets_dir}/{cell_cond1}.pkl', 'wb') as handle:
        pickle.dump(Gene_sets, handle)


'''
with open(f'{Gene_sets_dir}/{cell_cond1}.pkl', 'rb') as handle:
    Gene_sets_dict = pickle.load(handle)
'''