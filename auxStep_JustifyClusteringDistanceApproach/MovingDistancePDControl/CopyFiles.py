# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:49:56 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os


work_dir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 

copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file == 'ALL_G1_with_headers.csv':
            cell_cond = root.split('\\')[1]
            if not os.path.exists( f'{copy_to}/NMTF_G1s/{cell_cond}'):
                os.makedirs(f'{copy_to}/NMTF_G1s/{cell_cond}')
            copyfile(f'{root}/{file}', f'{copy_to}/NMTF_G1s/{cell_cond}/{file}')

                   

copy_from = 'step10_PubMedDEG_Validation/output/StageSpecPreds'
copy_to = f'{work_dir}/input/StageSpecPreds'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl') and 'Unique' not in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')

copy_from = 'step10_PubMedDEG_Validation/output/Other_genes'

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl') and 'Unique' not in file:
            if not os.path.exists(copy_to):
                os.makedirs(copy_to)
            copyfile(f'{root}/{file}', f'{copy_to}/{file}')


copy_to = f'{work_dir}/input'
copy_from = 'step7_Create4SetsOfGenes/output'
dir2copy = os.listdir(copy_from)
if os.path.exists(f'{copy_to}/{dir2copy[0]}'): 
    rmtree(f'{copy_to}/{dir2copy[0]}') 
copytree(f'{copy_from}/{dir2copy[0]}', f'{copy_to}/{dir2copy[0]}')
 
copyfile('CommonData/PD_genes_DGN_DEGs.pkl', f'{copy_to}/PD_genes_DGN_DEGs.pkl')
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info') 
copyfile('CommonData/Gene4PD.pkl', f'{copy_to}/Gene4PD.pkl') 
copyfile('CommonData/GWASdb_SNP.pkl', f'{copy_to}/GWASdb_SNP.pkl') 

copyfile('step10_PubMedDEG_Validation/output/CorePreds/All_CommonGenes_LitValid.pkl', f'{copy_to}/All_CommonGenes_LitValid.pkl') 



os.chdir(work_dir) 
