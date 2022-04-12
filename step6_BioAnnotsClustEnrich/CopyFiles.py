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


copy_from = 'step3_ClusterG1NMTF/output'
copy_to = f'{work_dir}/input'         

for root, dirs, files in os.walk(copy_from):
	for file in files:
		if 'kMeans' in file:
			cell_cond = root.split('\\')[1]
			nets = root.split('\\')[2]
			if not os.path.exists(f'{copy_to}/{cell_cond}/{nets}'):
				os.makedirs(f'{copy_to}/{cell_cond}/{nets}')
			copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{nets}/{file}')

copy_from = 'step4_PDgenesClustEnirch/output'
for root, dirs, files in os.walk(copy_from):
	for file in files:
		if 'indices' in file and 'kMeans' in file:
			cell_cond = root.split('\\')[1]
			nets = root.split('\\')[2]
			if not os.path.exists(f'{copy_to}/{cell_cond}/{nets}'):
				os.makedirs(f'{copy_to}/{cell_cond}/{nets}')
			copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{nets}/{file}')
            
copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input' 

for root, dirs, files in os.walk(copy_from):
	for file in files:
		if 'Geneslist' in file:
			cell_cond = root.split('\\')[1]
			if not os.path.exists(f'{copy_to}/{cell_cond}'):
				os.makedirs(f'{copy_to}/{cell_cond}')
			copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{file}')


copyfile('CommonData/HSA_GO-BP.lst', f'{copy_to}/HSA_GO-BP.lst')
copyfile('CommonData/HSA_Kegg_Pathways.lst', f'{copy_to}/HSA_Kegg_Pathways.lst')
copyfile('CommonData/HSA_Reactome_Pathways.lst', f'{copy_to}/HSA_Reactome_Pathways.lst')
copyfile('CommonData/go-basic.obo', f'{copy_to}/go-basic.obo')

os.chdir(work_dir) 
