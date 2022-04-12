# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 16:57:53 2021

@author: kmihajlo
"""

from shutil import copyfile, copytree, rmtree
import os

workdir = os.getcwd()
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{workdir}/input'):
    os.makedirs(f'{workdir}/input') 
if not os.path.exists(f'{workdir}/output'):
    os.makedirs(f'{workdir}/output')


copy_from = 'step9_StageSpec_CorePreds_BestSimMeas/output/StageSpecPreds/Euclidean_distance'
copy_to = f'{workdir}/input'         


for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.pkl'):
            print(root, file)
            sd = copy_from.split('/')[2]
            save_dir = f'{copy_to}/{sd}'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)             
            copyfile(f'{root}/{file}', f'{save_dir}/{file}')                       

copy_dir = 'step9_StageSpec_CorePreds_BestSimMeas/output/All_genes'
destination_dir = f'{copy_to}/All_genes'
if os.path.exists(destination_dir):
    rmtree(destination_dir) 
copytree(copy_dir, destination_dir)

copy_from = 'step2_NMTF/output'
save_dir = f'{copy_to}/All_genes'
for root, dirs, files in os.walk(copy_from):
    for file in files:
        if file.endswith('.csv') and 'Geneslist' in file:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)  
            copyfile(f'{root}/{file}', f'{save_dir}/{file}')                       

copyfile('step9_StageSpec_CorePreds_BestSimMeas/output/CorePreds_BetwStagePairs.pkl', f'{copy_to}/CorePreds_BetwStagePairs.pkl') 
copyfile('CommonData/Homo_sapiens.gene_info', f'{copy_to}/Homo_sapiens.gene_info') 
copyfile('CommonData/Gene4PD.pkl', f'{copy_to}/Gene4PD.pkl') 
copyfile('CommonData/GWASdb_SNP.pkl', f'{copy_to}/GWASdb_SNP.pkl') 
copyfile('CommonData/PD_genes_DGN_DEGs.pkl', f'{copy_to}/PD_genes_DGN_DEGs.pkl')             

copy_dir = 'CommonData/StageSpecDEGs/fold0_1'
destination_dir = f'{copy_to}/StageSpecDEGs_fold0_1'
if os.path.exists(destination_dir):
    rmtree(destination_dir) 
copytree(copy_dir, destination_dir)

os.chdir(workdir)
