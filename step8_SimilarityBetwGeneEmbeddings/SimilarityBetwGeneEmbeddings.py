# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:38:35 2021

@author: kmihajlo
"""

import pandas as pd
import os, pickle
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
from math import sqrt
 
# calculate manhattan distance
def manhattan_distance(a, b):
	return sum(abs(e1-e2) for e1, e2 in zip(a,b))
def euclidean_distance(a, b):
	return sqrt(sum((e1-e2)**2 for e1, e2 in zip(a,b)))

def write_txt(file_path, name, list_l):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')
in_dir = 'input'   

cell_conds = []
for root, dirs, files in os.walk(in_dir):
    for file in files:
        if 'NMTF_G1s' in root:
            cell_cond_1 = root.split('\\')[2]
            cell_conds.append(cell_cond_1)

base_dir = 'input/NMTF_G1s'

comparison_types = ['Cosine_similarity', 'Manhattan_distance', 'Euclidean_distance']
for comp_type in comparison_types:
    if not os.path.exists(f'output/{comp_type}'):
        os.makedirs(f'output/{comp_type}') 

for C1_i in range(len(cell_conds)): #C1 = cell condition 1
    cell_cond_1 = cell_conds[C1_i]
    cc_1d = cell_cond_1.split('_')[1]
    G1_1 = pd.read_csv(f'{base_dir}/{cell_cond_1}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
    for C2_j in range(C1_i, len(cell_conds)):
        cell_cond_2 = cell_conds[C2_j]
        cc_2d = cell_cond_2.split('_')[1]
        G1_2 = pd.read_csv(f'{base_dir}/{cell_cond_2}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')

        save_dir = 'output/SortedGenes'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        d1 = cell_cond_1.split('_')[1]
        d2 = cell_cond_2.split('_')[1]

        if cell_cond_1 != cell_cond_2 and d1 == d2:
            print(cell_cond_1, cell_cond_2)
            cellcond_1 = cell_cond_1.split('_')[0] + cell_cond_1.split('_')[1]
            with open(f'{in_dir}/GeneSets/GS_{cellcond_1}.pkl', 'rb') as handle:
                GeneSets = pickle.load(handle)
            CC1_genes = GeneSets[f'E_{cell_cond_2}'][f'E_{cell_cond_1}_ONLY']
            CC2_genes = GeneSets[f'E_{cell_cond_2}'][f'E_{cell_cond_2}_ONLY']
            if cc_1d == cc_2d:
                write_txt('output/RemovedGenesSameTP',f'{cell_cond_1}_RGs.txt', CC1_genes)
                write_txt('output/RemovedGenesSameTP',f'{cell_cond_2}_RGs.txt', CC2_genes)
            
            G1_1AE = G1_1.drop(CC1_genes)
            G1_2AE = G1_2.drop(CC2_genes)
            G1_1AE = G1_1AE.sort_index(axis=0)
            G1_1AE = G1_2AE.sort_index(axis=0)
            
            Geneslist = G1_1AE.index.tolist()

            with open(f'{save_dir}/{cell_cond_1}-{cell_cond_2}.pkl', 'wb') as handle:
                pickle.dump(Geneslist, handle)

        
            G1_1_vals = G1_1AE.to_numpy()
            G1_2_vals = G1_2AE.to_numpy()
            n = len(G1_1_vals)

            #cosine similarity
            cos_sim_all = cosine_similarity(G1_1_vals, G1_2_vals, dense_output=True)
            cos_sim = np.empty(n, float)
            for i in range(n):
                cos_sim[i] = cos_sim_all[i][i]
            save_file = f'output/{comparison_types[0]}/{cell_cond_1}-{cell_cond_2}'
            np.save(save_file, cos_sim) 
            
            #Manhattan distance
            Manh_dist = np.empty(n, float)
            for i in range(len(G1_1_vals)):
                gene1_embed = G1_1_vals[i]
                gene2_embed = G1_2_vals[i]
                manh_dist = manhattan_distance(gene1_embed, gene2_embed)
                Manh_dist[i] = manh_dist
            save_file = f'output/{comparison_types[1]}/{cell_cond_1}-{cell_cond_2}'
            np.save(save_file, Manh_dist)  
            
            #Euclidean distance
            Eucl_dist = np.empty(n, float)
            for i in range(len(G1_1_vals)):
                gene1_embed = G1_1_vals[i]
                gene2_embed = G1_2_vals[i]
                euc_dist = euclidean_distance(gene1_embed, gene2_embed)
                Eucl_dist[i] = euc_dist
            save_file = f'output/{comparison_types[2]}/{cell_cond_1}-{cell_cond_2}'
            np.save(save_file, Eucl_dist)          
                           
                            
            #cos_sim_df = pd.DataFrame(cos_sim, index = G1_1AE.index, columns=G1_1AE.index)
            
            '''
            n = len(G1_1_vals)
            Pearson_corr = np.empty((n, n), float)
            for i in range(len(G1_1_vals)):
                gene1_embed = G1_1_vals[i]
                if i%1000 == 0:
                    print(i)
                for j in range(len(G1_2_vals)):
                    gene2_embed = G1_2_vals[j]
                    pearson_corr, pvalue_pear = scipy.stats.pearsonr(gene1_embed, gene2_embed)
                    Pearson_corr[i][j] = pearson_corr
            
            Spearman_corr = np.empty((n, n), float)
            for i in range(len(G1_1_vals)):
                gene1_embed = G1_1_vals[i]
                for j in range(len(G1_2_vals)):
                    gene2_embed = G1_2_vals[j]
                    spearman_corr, pvalue_spear = scipy.stats.spearmanr(gene1_embed, gene2_embed)
                    Spearman_corr[i][j] = spearman_corr
            '''
            '''
            fig, ax = plt.subplots()
            im = ax.imshow(cos_sim)
            plt.title(f'{cell_cond_1} - {cell_cond_2}')
            plt.show()
            
            fig, ax = plt.subplots(figsize=(11, 9))
            hm = sb.heatmap(cos_sim,annot=False)
            plt.show()
            '''