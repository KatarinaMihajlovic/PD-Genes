# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:19:19 2022

@author: kmihajlo
"""

# Euclidean distance between PD-Control stage-spec pairs, all genes expressed in them
# order according to highest distance
# intersection and rank according to avg highest distance between all

import pandas as pd
import os, pickle, re, requests, math
import numpy as np
from math import sqrt
from scipy.stats import hypergeom, mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def euclidean_distance(a, b):
	return sqrt(sum((e1-e2)**2 for e1, e2 in zip(a,b)))

def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[int(lspt[1])] = Symbol
    return(Entrez2Sym)


def logPubmedn(cc1cc2_PDpreds_LitValid):
    cc1cc2_pubmed = []
    for gene in cc1cc2_PDpreds_LitValid.keys():
        logn = math.log(cc1cc2_PDpreds_LitValid[gene][0] + 1, 10)
        cc1cc2_pubmed.append(logn)  
    return cc1cc2_pubmed

def PlotHist_MWU(cc1cc2_pubmed, cc1cc2_AllGenes_pubmed, cc_pair, save_dir):
    ##### MWU test - shows that clustered genes have higher number of citations 
    x = cc1cc2_pubmed
    y = cc1cc2_AllGenes_pubmed
    statmwu,pvalmwu = mannwhitneyu(x,y, alternative='greater')
    # print(pvalmwu)
    
    ##computing the histograms
    num_bin = 50  
    lst = cc1cc2_pubmed + cc1cc2_AllGenes_pubmed
    minv = min(lst)
    maxv = max(lst)
        
    bin_lims = np.linspace(minv,maxv,num_bin+1)
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]
    
    hist1, _ = np.histogram(cc1cc2_AllGenes_pubmed, bins=bin_lims)
    hist2, _ = np.histogram(cc1cc2_pubmed, bins=bin_lims)
    
    ##normalizing
    hist1b = hist1/np.max(hist1)
    hist2b = hist2/np.max(hist2)
    
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_title(f'{cc_pair}: PD preds vs Other genes', fontsize=32, pad=20)    
    #ax.set_title('Core PD preds vs Other Shared genes', fontsize=22, pad=20)    
    
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.5, color='red')
    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.5, color='cornflowerblue')
    
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax.set_xlabel('log(num_PubMed + 1)', fontsize=28)
    ax.legend(handles= [Patch(facecolor='red', edgecolor='red',label=f'(1) PD preds ({len(x)})'),
                        Patch(facecolor='cornflowerblue', edgecolor='cornflowerblue',label=f'(2) Other Genes ({len(y)})')],
              loc='best', labelspacing=1, prop={'size': 24})
    pvalmwu = "{:.3e}".format(pvalmwu)
    fig.text(0.57, 0.63, f'MWU (1) > (2)\np-value = {pvalmwu}', fontsize=24)
    
    plt.savefig(f'{save_dir}/{cc_pair}_PDpredsvsOGs', dpi = 600)  
    plt.show()
    plt.close()    
 

def LitEnrich(All_litValid, gene_list, outf):  
    M = 0
    K = 0
    N = 0
    X = 0
    for gene in All_litValid:
        M+=1
        if gene in gene_list:
            N+=1
        if All_litValid[gene][0] > 0 or All_litValid[gene][1] != '' or All_litValid[gene][2] != '':
            K+=1
            if gene in gene_list:
                X+=1
    
    print(M,K,N,X)
    outf.write(f'{M}\t{K}\t{N}\t{X}\n')
    
    perc = X/N*100
    print(perc)
    outf.write(f'{perc}%\n')
    
    try:
        fold = (X/N)/(K/M)
        print(fold)
        outf.write(f'{fold}\n')
    
    except ZeroDivisionError:
        fold = 0
    
    if fold >= 1:
        pval = hypergeom.sf(X-1, M, K, N)
        print(pval)
        outf.write(f'{pval}\n')
    
        if pval <= 0.05: 
            #print(f'Enriched clust {i},  fold = {fold},  pval = {pval}')
            #print(N, X)
            print('Enriched in LitValid - PubMed,Gene4PD,GWASdb')    
            outf.write('Enriched in LitValid - PubMed,Gene4PD,GWASdb\n')
 
    outf.write('\n')
    print('\n')   
    

### MAIN CODE

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)

in_dir = 'input'   
cell_conds = []
for root, dirs, files in os.walk(in_dir):
    for file in files:
        if 'NMTF_G1s' in root:
            cell_cond_1 = root.split('\\')[2]
            cell_conds.append(cell_cond_1)
 

with open(f'{in_dir}/PD_genes_DGN_DEGs.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)
    
with open('input/All_CommonGenes_LitValid.pkl', 'rb') as handle:
    All_CommonGenes_LitValid = pickle.load(handle)    
common_genes = list(All_CommonGenes_LitValid.keys())

day_dict = {'0':'IPSCs', '6':'D06', '15':'D15', '21':'D21'}

topn = 172
base_dir = 'input/NMTF_G1s'
Euclid_dist_common_all = {key:None for key in cell_conds if 'PINK1' in key}       
Euclid_dist_allccgenes_dict = {key:None for key in cell_conds if 'PINK1' in key}      

for C1_i in range(len(cell_conds)): #C1 = cell condition 1
    cell_cond_1 = cell_conds[C1_i]
    cc_1d = cell_cond_1.split('_')[1]
    G1_1 = pd.read_csv(f'{base_dir}/{cell_cond_1}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')
    for C2_j in range(C1_i, len(cell_conds)):
        cell_cond_2 = cell_conds[C2_j]
        cc_2d = cell_cond_2.split('_')[1]
        G1_2 = pd.read_csv(f'{base_dir}/{cell_cond_2}/ALL_G1_with_headers.csv', header=0, index_col=0, delimiter='\t')

        # save_dir = 'output/SortedGenes'
        # if not os.path.exists(save_dir):
        #     os.makedirs(save_dir)
        d1 = cell_cond_1.split('_')[1]
        d2 = cell_cond_2.split('_')[1]

        if cell_cond_1 != cell_cond_2 and d1 == d2:
            print(cell_cond_1, cell_cond_2)
            cellcond_1 = cell_cond_1.split('_')[0] + cell_cond_1.split('_')[1]
            
            
            with open(f'{in_dir}/GeneSets/GS_{cellcond_1}.pkl', 'rb') as handle:
                GeneSets = pickle.load(handle)
            CC1_genes = GeneSets[f'E_{cell_cond_2}'][f'E_{cell_cond_1}_ONLY']
            CC2_genes = GeneSets[f'E_{cell_cond_2}'][f'E_{cell_cond_2}_ONLY']

            G1_1AE = G1_1.drop(CC1_genes)
            G1_2AE = G1_2.drop(CC2_genes)
            G1_1AE = G1_1AE.sort_index(axis=0)
            G1_1AE = G1_2AE.sort_index(axis=0)
            
            Geneslist = G1_1AE.index.tolist()

            # with open(f'{save_dir}/{cell_cond_1}-{cell_cond_2}.pkl', 'wb') as handle:
            #     pickle.dump(Geneslist, handle)

        
            G1_1_vals = G1_1AE.to_numpy()
            G1_2_vals = G1_2AE.to_numpy()
            n = len(G1_1_vals)

            #Euclidean distance
            Euclid_dist_common = []
            Euclid_dist_allccgenes = []
            # Eucl_dist = np.empty(n, float)
            for i in range(len(G1_1_vals)):
                gene1_embed = G1_1_vals[i]
                gene2_embed = G1_2_vals[i]
                euc_dist = euclidean_distance(gene1_embed, gene2_embed)
                # Eucl_dist[i] = euc_dist
                gene = Entrez2Sym[Geneslist[i]]
                if gene not in PD_genes:
                    if gene in common_genes:
                        Euclid_dist_common.append([gene, euc_dist])
                    Euclid_dist_allccgenes.append([gene, euc_dist])
            Euclid_dist_common.sort(key=lambda x: x[1])
            Euclid_dist_allccgenes.sort(key=lambda x: x[1])            
            Euclid_dist_common_all[cell_cond_2] = Euclid_dist_common
            Euclid_dist_allccgenes_dict[cell_cond_2] = Euclid_dist_allccgenes
     
save_dir = 'output'    
           
    
### Stage Specific predictions
file_out = 'StageSpecValids.txt'
outf = open(f'{save_dir}/{file_out}', 'w')
for root, dirs, files in os.walk('input/StageSpecPreds'):
    for file in files:
        if 'OtherGenes' not in file:
            print(file)
            ss = file.split('_')[0]
            with open(f'{root}/{file}', 'rb') as handle:
                StageSpecPreds = pickle.load(handle)   
            day = file.split('P')[0].split('C')[1]
            with open(f'{root}/C{day}PD{day}_OtherGenes.pkl', 'rb') as handle:
                cc1cc2_OtherGenes = pickle.load(handle)  
            cc1cc2_AllGenes = {**cc1cc2_OtherGenes , **StageSpecPreds}
            for gene in PD_genes:
                cc1cc2_AllGenes.pop(gene, None)
            
            cc = f'PINK1_{day_dict[day]}'
            Euclid_dist_cc = Euclid_dist_allccgenes_dict[cc] 
            topPDgenes = Euclid_dist_cc[:len(StageSpecPreds)]
            topPDgenes = [x[0] for x in topPDgenes]  
            
            outf.write(f'{file}\n')
            LitEnrich(cc1cc2_AllGenes, topPDgenes, outf)
            
            #pubmed
            for gene in PD_genes:
                cc1cc2_AllGenes.pop(gene, None)
            topPDgenes_dict = { your_key: cc1cc2_AllGenes[your_key] for your_key in topPDgenes }
            for gene in topPDgenes:
                cc1cc2_AllGenes.pop(gene, None)

            topPDgenes_pubmed = logPubmedn(topPDgenes_dict)
            OtherGenes_pubmed = logPubmedn(cc1cc2_AllGenes) #all genes a OtherGenes_LitValid
            PlotHist_MWU(topPDgenes_pubmed, OtherGenes_pubmed, ss, save_dir)
outf.close()  
    
    
### Core predictions comparisson            
Euclid_dist_common_all_pd = pd.DataFrame(index= common_genes, columns = Euclid_dist_common_all.keys())
for cell_cond in Euclid_dist_common_all:
    for i in range(len(Euclid_dist_common_all[cell_cond])):
        gene = Euclid_dist_common_all[cell_cond][i][0]
        if gene not in PD_genes:
            distance = Euclid_dist_common_all[cell_cond][i][1]
            Euclid_dist_common_all_pd.loc[gene,cell_cond] = distance            

Euclid_dist_common_all_pd['average'] = Euclid_dist_common_all_pd.mean(axis=1)
Euclid_dist_common_all_pd = Euclid_dist_common_all_pd.sort_values(by='average', ascending=False)
common_genes_ordered = Euclid_dist_common_all_pd.index.tolist()
topnCGs = common_genes_ordered[:topn]    

for gene in PD_genes:
    All_CommonGenes_LitValid.pop(gene, None)       

file_out = 'CGs172_LitValid_Enirch.txt'
outf = open(f'{save_dir}/{file_out}', 'w')
print('Core Preds\n')
LitEnrich(All_CommonGenes_LitValid, topnCGs, outf)
outf.close()   

topnCGs_dict = { your_key: All_CommonGenes_LitValid[your_key] for your_key in topnCGs }
for gene in topnCGs:
    All_CommonGenes_LitValid.pop(gene, None)

topnCGs_pubmed = logPubmedn(topnCGs_dict)
OtherGenes_pubmed = logPubmedn(All_CommonGenes_LitValid) #all genes a OtherGenes_LitValid
PlotHist_MWU(topnCGs_pubmed, OtherGenes_pubmed, 'top172genes', save_dir)

                            