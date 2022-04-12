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


if not os.path.exists(f'{work_dir}/input/NMTF_G1s'):
    os.makedirs(f'{work_dir}/input/NMTF_G1s') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_from = 'step2_NMTF/output'
copy_to = f'{work_dir}/input/NMTF_G1s'         

for root, dirs, files in os.walk(copy_from):
    for file in files:
        if 'G1_with_headers' in file and 'ALL' in root:
            cell_cond = root.split('\\')[1]
            if not os.path.exists( f'{copy_to}/{cell_cond}'):
                os.makedirs(f'{copy_to}/{cell_cond}')
            copyfile(f'{root}/{file}', f'{copy_to}/{cell_cond}/{file}')

copy_from = 'step7_Create4SetsOfGenes/output'
copy_to = f'{work_dir}/input'         

dir2copy = os.listdir(copy_from)
if os.path.exists(f'{copy_to}/{dir2copy[0]}'): 
    rmtree(f'{copy_to}/{dir2copy[0]}') 
copytree(f'{copy_from}/{dir2copy[0]}', f'{copy_to}/{dir2copy[0]}')

os.chdir(work_dir) 
