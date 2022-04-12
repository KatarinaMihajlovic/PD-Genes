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


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

copy_to = f'{work_dir}/input'  
copyfile('step9_StageSpec_CorePreds_BestSimMeas/output/CorePreds_BetwStagePairs.pkl', f'{copy_to}/CorePreds_BetwStagePairs.pkl') 


copy_dir = 'step10_PubMedDEG_Validation/input/All_genes'
destination_dir = f'{copy_to}/All_genes'
if os.path.exists(destination_dir):
    rmtree(destination_dir) 
copytree(copy_dir, destination_dir)

'''
copy_from = 'step10_PubMedDEG_Validation/output/StageSpecPreds/Unique_StageSpecPreds'
copy_to = f'{work_dir}/input/Unique_StageSpecPreds'
if not os.path.exists(copy_to):
    os.makedirs(copy_to)

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl') and 'Unique' in file:
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')  
'''

copy_from = 'CommonData'
copy_to = f'{work_dir}/input'  

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if not file.endswith('.txt'):
            copyfile(f'{copy_from}/{file}', f'{copy_to}/{file}')    
copyfile('CommonData/compound_names.txt', f'{copy_to}/compound_names.txt')
copyfile('CommonData/genes_compounds.txt', f'{copy_to}/genes_compounds.txt')
copyfile('CommonData/HSA_Reactome_Pathways_meaning.lst', f'{copy_to}/HSA_Reactome_Pathways_meaning.lst')

# copyfile('PD_genesGOenrich/output/PDgenes_KPmeaning.txt', f'{copy_to}/PDgenes_KPmeaning.txt')
# copyfile('PD_genesGOenrich/input/PD_genes_DGN_DEGs.pkl', f'{copy_to}/PD_genes_DGN_DEGs.pkl')


os.chdir(work_dir)
