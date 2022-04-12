# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 20:00:30 2021

@author: kmihajlo
"""

import pickle, os
import pandas as pd 
from scipy.stats import hypergeom
from random import sample
import matplotlib.pyplot as plt


labels_key = {'Control_IPSCs':'C0', 'Control_D06':'C6', 'Control_D10':'C10', 'Control_D15':'C15', 'Control_D21':'C21', 
              'PINK1_IPSCs':'PD0', 'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21'}


def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]
        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)
 

def Sort(sub_li):
  
    sub_li.sort(key = lambda x: x[1], reverse=True)
    return sub_li


def IsolatePDenrichCls(clst_opt, All_litValidGenes, PD_genes):
    Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
    Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)

    Best_cm_allclusts = {'kMeans':[], 'kMedoids':[], 'HardClusters':[]}
    Best_cm_PINK1_allclusts = {'kMeans':[], 'kMedoids':[], 'HardClusters':[]}

    sd = f'output/All_clust/{clst_opt}'
    if not os.path.exists(sd):
        os.makedirs(sd) 
    outfile_stats = open(f'{sd}/PDgenes_stats.txt','w')
    for root, dirs, files in os.walk('input/Clusters'):   
        for file in files:
            if file.endswith('.pkl'):
                print(root, file)
                clus_meth = file.split('_')[0]          
    
                cell_cond = root.split('\\')[1]
                
                save_dirA = f'output/All_clust/{clus_meth}/{cell_cond}'
                if not os.path.exists(save_dirA):
                    os.makedirs(save_dirA) 
    
                Geneslist = pd.read_csv(f'input/Genelists/Geneslist_{cell_cond}.csv')
                Geneslist = Geneslist['genes'].tolist()
                Geneslist = [Entrez2Sym[str(x)] for x in Geneslist]    

                with open(f'{root}/{file}', 'rb') as handle:
                    ClustsPDEnrich = pickle.load(handle)
                
                
                i = 0
                for clst in range(len(ClustsPDEnrich)):
                    print(clst)
                    PD_gene_Predictions = []
                    outfile_name = f'{clst}_{labels_key[cell_cond]}.txt'
                    outfile = open(f'{save_dirA}/{outfile_name}','w')
                    EnrClust = ClustsPDEnrich[clst]
                    Litvalid_genes = []
                    for gene in EnrClust:
                        if gene not in PD_genes:
                            PD_gene_Predictions.append(gene)
                            i+=1
                            outfile.write(f'{gene}\t{LitValid_AllGenes[gene][0]}\t{LitValid_AllGenes[gene][1]}\t{LitValid_AllGenes[gene][2]}\n')

                            # check also from pubmed
                            try:
                                if LitValid_AllGenes[gene][0] > 0:
                                    Litvalid_genes.append(gene)
                                elif LitValid_AllGenes[gene][1] != '':
                                    Litvalid_genes.append(gene)
                                elif All_litValidGenes[gene][2] != '':
                                    Litvalid_genes.append(gene)
                            except KeyError:
                                print(gene)
                    outfile.close()
    
                                        
                    N = len(PD_gene_Predictions)
                    # X = len(Gene4PDGWAS_Proof_all)
                    X = len(Litvalid_genes)
                    # print(N)
    
                    Best_cm_allclusts[clus_meth].append(X/N*100)
                    if 'PINK1' in root:
                        Best_cm_PINK1_allclusts[clus_meth].append(X/N*100)  
                print(cell_cond, i)
                if clus_meth == clst_opt:
                    outfile_stats.write(f'{cell_cond}\t{i}\t{len(ClustsPDEnrich)}\n')
                       
    outfile_stats.close()    
    return Best_cm_allclusts, Best_cm_PINK1_allclusts



### Main Code

with open('input/PD_genes_DGN_DEGs.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle)  
with open('input/LitValid_AllGenes.pkl', 'rb') as handle:
    LitValid_AllGenes = pickle.load(handle)  

Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)
All_genes = pd.read_csv('input/All_genes.csv', header = None)      
All_genes = All_genes[0].tolist()  
All_genes = [str(x) for x in All_genes]
All_genes_sym = [Entrez2Sym[str(x)] for x in All_genes]

clusts_opts = ['kMeans', 'kMedoids', 'HardClusters']
for clst_opt in clusts_opts:
    Best_cm_allclusts, Best_cm_PINK1_allclusts = IsolatePDenrichCls(clst_opt, LitValid_AllGenes, PD_genes)

                 
Best_cm_list_allclusts = []
#print(Best_cm_dist)
for cm in Best_cm_allclusts.keys():
    avg = sum(Best_cm_allclusts[cm])/len(Best_cm_allclusts[cm])
    Best_cm_list_allclusts.append([cm,avg,len(Best_cm_allclusts[cm])])  

Best_cm_list_allclusts = Sort(Best_cm_list_allclusts)
with open('output/Best_cm_allclusts.txt', 'w') as f:
    for i in range(len(Best_cm_list_allclusts)):
        f.write(f'{Best_cm_list_allclusts[i][0]}\t{Best_cm_list_allclusts[i][1]}\t{Best_cm_list_allclusts[i][2]}\n')                        
                  
Best_cm_list_PINK1_allclusts = []
for cm in Best_cm_PINK1_allclusts.keys():
    avg = sum(Best_cm_PINK1_allclusts[cm])/len(Best_cm_PINK1_allclusts[cm])
    Best_cm_list_PINK1_allclusts.append([cm,avg,len(Best_cm_PINK1_allclusts[cm])])  

Best_cm_list_PINK1_allclusts = Sort(Best_cm_list_PINK1_allclusts)
with open('output/Best_cm_PINK1_allclusts.txt', 'w') as f:
    for i in range(len(Best_cm_list_PINK1_allclusts)):
        f.write(f'{Best_cm_list_PINK1_allclusts[i][0]}\t{Best_cm_list_PINK1_allclusts[i][1]}\t{Best_cm_list_PINK1_allclusts[i][2]}\n')                                    
                        
                        
# with open('input/Clusters/Control_D21/kMeans_PDEnrichClusts.pkl', 'rb') as handle:
#     dsa = pickle.load(handle)  
                        
                        