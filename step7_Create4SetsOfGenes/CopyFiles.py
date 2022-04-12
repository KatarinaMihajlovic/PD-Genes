from shutil import copyfile
import os

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/input/EXP_matrices'):
    os.makedirs(f'{work_dir}/input/EXP_matrices') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')

copy_from = 'step1_CreateInputNetworksMatrices/output'
copy_to = f'{work_dir}/input'

copyfile(f'{copy_from}/All_genes.csv', f'{copy_to}/All_genes.csv')

subdirs = []
cluster_folders = set()
for root, dirs, files in os.walk(f'{copy_from}/Expression_Matrix'):
    for file in files:
        if file.endswith('.csv'):
            copyfile(f'{copy_from}/Expression_Matrix/{file}', f'{copy_to}/EXP_matrices/{file}')

os.chdir(work_dir) 
