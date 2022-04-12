# -*- coding: utf-8 -*-
"""
Created on Thu Oct  14 12:48:56 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input/PD_preds'):
    os.makedirs(f'{work_dir}/input/PD_preds') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_from1 = 'step5_PDEnrichClusts_InitPreds_BestClustMeth/output/All_clust/kMeans'
copy_to = f'{work_dir}/input'     


for root, dirs, files in os.walk(copy_from1):
    for dir2copy in dirs:

    
        if os.path.exists(f'{copy_to}/PD_preds/{dir2copy}'): 
            rmtree(f'{copy_to}/PD_preds/{dir2copy}') 
        copytree(f'{root}/{dir2copy}', f'{copy_to}/PD_preds/{dir2copy}') 


copy_from2 = 'step8_SimilarityBetwGeneEmbeddings/output'
dirs2copy = os.listdir(copy_from2)
dirs2copy.pop(3)
save_dir2 = f'{copy_to}/{dirs2copy[3]}'
if not os.path.exists(save_dir2):
    os.makedirs(save_dir2)   
    
for dir2copy in dirs2copy[:-1]:
    save_dir1 = f'{copy_to}/Sim_measure/{dir2copy}'
    if not os.path.exists(save_dir1):
        os.makedirs(save_dir1)


for dir2copy in dirs2copy[:-1]:    
    for root, dirs, files in os.walk(f'{copy_from2}/{dir2copy}'):
        for file in files:
            print(file)
            cell_type1 = file.split('-')[0].split('_')[0]
            cell_stage1 = file.split('-')[0].split('_')[1]
            cell_type2 = file.split('-')[1].split('_')[0]
            cell_stage2 = file.split('-')[1].split('_')[1].split('.')[0]
            
            save_dir1 = f'{copy_to}/Sim_measure/{dir2copy}'
            if cell_type1 != cell_type2 and cell_stage1 == cell_stage2:
                copyfile(f'{copy_from2}/{dir2copy}/{file}', f'{save_dir1}/{file}')
                file2 = file.split('.')[0] + '.pkl'
                copyfile(f'{copy_from2}/{dirs2copy[3]}/{file2}', f'{save_dir2}/{file2}')   


# copyfile('CommonData/Gene4PD.pkl', f'{copy_to}/Gene4PD.pkl')
# copyfile('CommonData/GWASdb_SNP.pkl', f'{copy_to}/GWASdb_SNP.pkl') 
copyfile('CommonData/LitValid_AllGenes.pkl', f'{copy_to}/LitValid_AllGenes.pkl')
copyfile('CommonData/PD_genes_DGN_DEGs.pkl', f'{copy_to}/PD_genes_DGN_DEGs.pkl') 
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info')          

os.chdir(work_dir) 
