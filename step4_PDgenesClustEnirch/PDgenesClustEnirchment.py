# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:09:01 2021

@author: kmihajlo
"""

import os, pickle
from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


legend = {'Control_D06':'C6','Control_D15':'C15', 'Control_D10':'C10', 'Control_D21':'C21', 'Control_IPSCs':'C0',
          'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21', 'PINK1_IPSCs':'PD0'}

in_dir = 'input' 
out_dir = 'output'
with open(f'{in_dir}/PD_genes_DGN_DEGs.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)

print(len(PD_genes))
Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)
   
perc_enr_cluster_all = []    
PDGenes_EnrCLusts_all_perc = []
ccs = []
for root, dirs, files in os.walk(in_dir):
    for file in files:
        if file != 'PD_genes_DGN_DEGs.pkl' and '.gene' not in file:
            cell_cond = root.split('\\')[1]
            #print(cell_cond)
            nets = root.split('\\')[2]
            clustmeth = file.split('.')[0].split('_')[0]  
            
            if clustmeth == 'kMeans' and nets == 'ALL':
                ccs.append(legend[cell_cond])

            #print(clustmeth, nets)
            with open(f'{root}/{file}', 'rb') as handle:
                G1_clust = pickle.load(handle)
            genes = [str(item) for sublist in G1_clust for item in sublist]
            genes = [Entrez2Sym[gene] for gene in genes]

            Possible_PDgenes = list(set(genes) & set(PD_genes))
            #print(PDdb)
            #Possible_PDgenes = list(set(genes) & set(PD_genes_databases['DisGeNet']))
       
            save_dir = f'{out_dir}/{cell_cond}/{nets}'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir) 
            
            Enriched_clusts = []
            Pred_PDgenes_clusts = []
            pvals = []
            Enriched_clusts_i = []
            cont = 0
            PDGenes_EnrCLusts = 0
            for i in range(len(G1_clust)): 
                genes_clust = G1_clust[i]
                genes_clust = [Entrez2Sym[str(x)] for x in genes_clust]
                PDgenes_clust = [x for x in Possible_PDgenes if x in genes_clust]
                
                
                M = len(genes)
                K = len(Possible_PDgenes)
                N = len(genes_clust)
                X = len(PDgenes_clust)
                try:
                    fold = (X/N)/(K/M)
                except ZeroDivisionError:
                    fold = 0
                if fold >= 1:
                    pval = hypergeom.sf(X-1, M, K, N)
                    if pval <= 0.05: 
                        #print(f'Enriched clust {i},  fold = {fold},  pval = {pval}')
                        #print(N, X)
                        Enriched_clusts_i.append(i)
                        Enriched_clusts.append(genes_clust)
                        Pred_PDgenes = [x for x in genes_clust if x not in Possible_PDgenes]
                        Clst_PDgenes = [x for x in genes_clust if x in Possible_PDgenes]
                        Pred_PDgenes_clusts.append(Pred_PDgenes)
                        PDGenes_EnrCLusts+= len(Clst_PDgenes)
                        pvals.append(pval)
                        cont+=1
            if clustmeth == 'kMeans'  and nets == 'ALL':
                enr_cluster = len(Enriched_clusts)
                total_cluster = len(G1_clust)
                perc_enr_cluster = 100.*enr_cluster/total_cluster
                perc_enr_cluster_all.append(perc_enr_cluster)
                PDGenes_EnrCLusts_all_perc.append(PDGenes_EnrCLusts/len(Possible_PDgenes)*100)

            # Benjamini-Hochberg p-value correction   
            pvals_adj =  p_adjust_bh(pvals)
            pvals_adj_smaller = [p for p in pvals_adj if p <= 0.05]
            indexes = []
            for i in range(len(pvals_adj)):
                if pvals_adj[i] > 0.05:
                    indexes.append(i)
            
            if len(indexes) > 1:
                for index in sorted(indexes, reverse=True):
                    del Pred_PDgenes_clusts[index]
                    del Enriched_clusts[index]
                    del Enriched_clusts_i[index]
                
            Pred_PDgenes_cc = [item for sublist in Pred_PDgenes_clusts for item in sublist]
            Pred_PDgenes_cc = set(Pred_PDgenes_cc)
            if nets == 'ALL':
                print(f'{cell_cond} {clustmeth}  {len(Pred_PDgenes_clusts)}')
            
            
            with open(f'{save_dir}/{clustmeth}_Pred_PDgenes.pkl', 'wb') as handle:
                pickle.dump(Pred_PDgenes_cc, handle)
            with open(f'{save_dir}/{clustmeth}_PredPDgenesClusts.pkl', 'wb') as handle:
                pickle.dump(Pred_PDgenes_clusts, handle)
            with open(f'{save_dir}/{clustmeth}_PDEnrichClusts.pkl', 'wb') as handle:
                pickle.dump(Enriched_clusts, handle)   
            with open(f'{save_dir}/{clustmeth}_PDEnrichClusts_indices.pkl', 'wb') as handle:
                pickle.dump(Enriched_clusts_i, handle)                  

print(PDGenes_EnrCLusts_all_perc)
            
ccs = ['C6', 'C15', 'C21', 'C0', 'PD6', 'PD15', 'PD21', 'PD0']



ccs.insert(0, ccs[3])
ccs.pop(4)
ccs.insert(4, ccs[-1])
ccs.pop(-1)
perc_enr_cluster_all.insert(0, perc_enr_cluster_all[3])
perc_enr_cluster_all.pop(4)
perc_enr_cluster_all.insert(4, perc_enr_cluster_all[-1])
perc_enr_cluster_all.pop(-1)


### Plot perc of clusters (kmeans) enriched in PD genes
N = len(ccs)
ind = np.arange(N)  # the x locations for the groups
width = 0.2    

fig, ax = plt.subplots(figsize=(16, 9))
# Remove top and right border
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
# Remove y-axis tick marks
ax.yaxis.set_ticks_position('none')

# Set plot title
# ax.set_title(f'GCV-G - {cc}', fontsize = 22)
# Add major gridlines in the y-axis
ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
ax.set_facecolor('xkcd:white')
 
for i in range(N):
	rects1 = ax.bar(ind[i-1], perc_enr_cluster_all[i-1], width, color='cornflowerblue')

ax.set_ylabel('Clusters enriched in PD genes (%)', fontsize = 26, fontweight = 'bold')
ax.set_xticks(ind)
ax.set_xticklabels(ccs,  fontsize = 26, rotation=90) 
ax.tick_params(axis='y', which='major', labelsize=22)
plt.savefig('output/kmeans_ECs.png', dpi = 600)	

plt.show()
plt.close()
               
# with open('output/PINK1_D21/ALL/kMeans_PredPDgenesClusts.pkl', 'rb') as handle:
# 	kMeans_PredPDgenesClusts = pickle.load(handle)

# lst = []
# for li in kMeans_PredPDgenesClusts:
#     for gene in li:
#         lst.append(gene)