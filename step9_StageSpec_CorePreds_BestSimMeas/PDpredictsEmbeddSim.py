# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 13:13:31 2021

@author: kmihajlo
"""

import os, pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns; sns.set_theme(color_codes=True)
#from scipy.cluster.hierarchy import dendrogram, linkage
#from sklearn.cluster import AgglomerativeClustering
# from scipy.stats import zscore, hypergeom, kstest, ttest_ind, mannwhitneyu
from scipy.stats import mannwhitneyu

from matplotlib.lines import Line2D
from random import sample
from matplotlib.patches import Patch
from itertools import repeat

# with open('input/Gene4PD.pkl', 'rb') as handle:
#     Gene4PD = pickle.load(handle)  
# with open('input/GWASdb_SNP.pkl', 'rb') as handle:
#     GWASdb_SNP = pickle.load(handle)  

legend = {'Control_D06-PINK1_D06':'C6PD6', 'Control_D15-PINK1_D15':'C15PD15', 
          'Control_D21-PINK1_D21':'C21PD21', 'Control_IPSCs-PINK1_IPSCs':'C0PD0'}
legend_r = {'C6PD6':'Control_D06-PINK1_D06', 'C15PD15':'Control_D15-PINK1_D15', 
          'C21PD21':'Control_D21-PINK1_D21', 'C0PD0':'Control_IPSCs-PINK1_IPSCs'}

with open('input/LitValid_AllGenes.pkl', 'rb') as handle:
    LitValid_AllGenes = pickle.load(handle)  
    
def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]
        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

def write_txt(file_path, name, list_l):
    if not os.path.exists(file_path):
        os.makedirs(file_path) 
    with open(f'{file_path}/{name}', 'w') as f: 
        for x in list_l:
            f.write(f'{x}\n')
            
def compare2ccs(ccs2compare, PD_preds, distance_meas = 'Euclidean_distance', n_genes = 20, save_dir = 'output/PD_Predictions'):    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)   
    sim_np = np.load(f'input/Sim_measure/{distance_meas}/{ccs2compare}.npy')

    with open(f'input/SortedGenes/{ccs2compare}.pkl', 'rb') as handel:
        gene_order = pickle.load(handel)
        
    gene_order_sym = [Entrez2Sym[str(x)] for x in gene_order]
    sim_dist_df = pd.DataFrame(sim_np, index = gene_order_sym)
    Genes2drop = [gene for gene in gene_order_sym if gene not in PD_preds] 
    simdist_PDpreds_df = sim_dist_df.drop(Genes2drop)
    
    RemPDpreds = list(sim_dist_df.index)
    RemPDpreds = [x for x in PD_preds if x not in RemPDpreds]
    write_txt('output/AllPredictions/RemovedPDpreds', f'{legend[ccs2compare]}_RemPDpreds.txt', RemPDpreds)
    
    # gene2gene = []
    # for i in range(len(simdist_PDpreds_df.index)):
    #     gene = simdist_PDpreds_df.index[i]
    #     gene2gene.append(simdist_PDpreds_df[gene][gene])
    # gene2gene_df = pd.DataFrame(data = gene2gene, columns = [0], index = simdist_PDpreds_df.index)
    gene2gene_df = simdist_PDpreds_df
    if distance_meas != 'Cosine_similarity':
        gene2gene_df = gene2gene_df.sort_values(by = 0, axis=0, ascending=False)
    else:
        gene2gene_df = gene2gene_df.sort_values(by = 0, axis=0, ascending=True)
    gene2gene_df = gene2gene_df.T
    
    c1c2_genes = gene2gene_df.columns.tolist()
    if n_genes > 0:
        topxgenes_c1c2 = c1c2_genes[:n_genes]
    else:
        topxgenes_c1c2 = c1c2_genes
    # fig, ax = plt.subplots(figsize=(25, 9))
    # hm = sns.heatmap(gene2gene_df, annot=True)#, cmap = 'bwr')
    # ax.set_title(ccs2compare, fontsize=20)
    # #plt.savefig(f'{save_dir}/{ccs2compare}')
    # plt.show()
    # plt.close()
    
    fig, ax = plt.subplots(figsize=(16, 9))
    gene2gene_df.T.plot.hist(bins=50) #fix this
    ax.set_title(ccs2compare, fontsize=20)
    ax.set_xlabel(f'{distance_meas}')
    #plt.savefig(f'{save_dir}/{ccs2compare}_dist')
    plt.show()
    plt.close()
    #print(gene2gene_df.columns.tolist())
    return gene2gene_df, c1c2_genes, topxgenes_c1c2
        

def savePDPreds(cc1cc2_genes, out_filename, save_dir = 'output/PD_Predictions', LitValid_AllGenes = LitValid_AllGenes):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)        
    with open(f'{save_dir}/{out_filename}.txt', 'w') as outfile:
        for gene in cc1cc2_genes:
            outfile.write(f'{gene}\t{LitValid_AllGenes[gene][0]}\t{LitValid_AllGenes[gene][1]}\t{LitValid_AllGenes[gene][2]}\n')


def create_ccpair_df(sim_measure, ccs2compare):
    sim_np = np.load(f'input/Sim_measure/{sim_measure}/{ccs2compare}.npy')
    with open(f'input/SortedGenes/{ccs2compare}.pkl', 'rb') as handel:
        gene_order = pickle.load(handel)           
    gene_order_sym = [Entrez2Sym[str(x)] for x in gene_order]
    sim_dist_df = pd.DataFrame(sim_np, index = gene_order_sym)
    return sim_dist_df

def PDGs2Background_dist(ccs2compare, sim_measure, cc1cc2_genes, PD_Gs = 'PD preds'):
    sim_dist_df = create_ccpair_df(sim_measure, ccs2compare)
    all_genes = sim_dist_df.index.to_list()
    sd = 'output/All_genes'
    if not os.path.exists(sd):
        os.makedirs(sd) 
    with open(f'{sd}/{legend[ccs2compare]}_AllGenes.pkl', 'wb') as fp:   #Pickling
        pickle.dump(all_genes, fp)    
    
    Dist_gene2gene = []
    cc1cc2_genes_geneDist = []
    for gene in sim_dist_df.index:
        if gene in cc1cc2_genes:
            cc1cc2_genes_geneDist.append(sim_dist_df.loc[gene,0])
        else:
            Dist_gene2gene.append(sim_dist_df.loc[gene,0])
    
    print(len(cc1cc2_genes_geneDist), len(Dist_gene2gene))
    if sim_measure != 'Cosine_similarity':
        statmwu,pvalmwu = mannwhitneyu(cc1cc2_genes_geneDist,Dist_gene2gene, alternative='greater')
    else:
        statmwu,pvalmwu = mannwhitneyu(cc1cc2_genes_geneDist,Dist_gene2gene, alternative='less')
    print(pvalmwu)
    
    ##computing the histograms
    num_bin = 50  
    lst = list(Dist_gene2gene)
    minv = min(lst)
    maxv = max(lst)
        
    bin_lims = np.linspace(minv,maxv,num_bin+1)
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]

    hist1, _ = np.histogram(cc1cc2_genes_geneDist, bins=bin_lims)
    hist2, _ = np.histogram(Dist_gene2gene, bins=bin_lims)

    ##normalizing
    hist1b = hist1/np.max(hist1)
    hist2b = hist2/np.max(hist2)

    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_facecolor('xkcd:white')
    ax.grid(color='grey', axis='x', linestyle='-', linewidth=0.25, alpha=0.5)
    ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
    title = legend[ccs2compare].split('P')
    title = title[0] + '-P' + title[1]
    ax.set_title(title, fontsize=34, pad=20) #': {PD_Gs} vs Other genes'   

    ax.bar(bin_centers, hist1b, width = bin_widths, align = 'center', alpha = 0.33, color='red')
    ax.bar(bin_centers, hist2b, width = bin_widths, align = 'center', alpha = 0.33, color='cornflowerblue')

    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax.set_xlabel(sim_measure, fontsize=28)
    
    pvalmwu = "{:.3e}".format(pvalmwu)
    if sim_measure != 'Cosine_similarity':
        fig.text(0.6, 0.6, f'MWU (1) > (2)\np-value = {pvalmwu}\n', fontsize=24)
    else:
        fig.text(0.6, 0.6, f'MWU (1) < (2)\np-value = {pvalmwu}\n', fontsize=24)
    ax.legend(handles= [Patch(facecolor='red', edgecolor='red',label=f'(1) {PD_Gs} ({len(cc1cc2_genes_geneDist)})'),
                Patch(facecolor='cornflowerblue', edgecolor='cornflowerblue',label=f'(2) Other genes ({len(Dist_gene2gene)})')],
              loc=1, labelspacing=1, prop={'size': 24}, facecolor ='white')
    
    plt.savefig(f'{save_dir}/{legend[ccs2compare]}_{PD_Gs}vsOG', dpi = 600)                        
    plt.show()
    plt.close()
        
def Plot_precorrecall(pd_DF, PR_RL, condit = 'no_cond'):
    sim_col = {'Cosine_similarity': 'green', 'Euclidean_distance' : 'red', 'Manhattan_distance': 'blue'}
    cc_pair_marker = {'Control_D06-PINK1_D06':'s', 'Control_D15-PINK1_D15':'P', 
              'Control_D21-PINK1_D21':'*', 'Control_IPSCs-PINK1_IPSCs':'o'}
   
    fig, ax = plt.subplots(figsize=(16, 9))
    if condit == 'no_cond':
        for sim_measure in pd_DF.columns:    
            for lst_ccpair in pd_DF.index:
                plt.plot(ks,pd_DF[sim_measure][lst_ccpair], marker = cc_pair_marker[lst_ccpair], color = sim_col[sim_measure], alpha = 0.4, markersize=12)
    else:
        for sim_measure in pd_DF.columns:
            plt.plot(ks,pd_DF[sim_measure][condit], marker = cc_pair_marker[condit], color = sim_col[sim_measure], alpha = 0.4, markersize=12)
        
        
    ax.set_ylabel(PR_RL, fontsize = 24)    
    ax.set_xlabel('top N-ranked genes', fontsize = 24)  
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_facecolor('xkcd:white')
    ax.grid(color='grey', axis='x', linestyle='-', linewidth=0.25, alpha=0.5)
    ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
    
    
    if 'Cosine_similarity' not in pd_DF.columns:
        legend = ax.legend(handles= [Patch(facecolor='red', edgecolor='red',label='Euclidean dist'),
                            Patch(facecolor='blue', edgecolor='blue',label='Manhattan dist'),
                            Line2D([0], [0], marker='o', color='b', label='C0PD0', markerfacecolor='w', markersize=12), 
                            Line2D([0], [0], marker='s', color='b', label='C6PD6', markerfacecolor='w', markersize=12), 
                            Line2D([0], [0], marker='P', color='b', label='C15PD15', markerfacecolor='w', markersize=12),
                            Line2D([0], [0], marker='*', color='b', label='C21PD21', markerfacecolor='w', markersize=12)],

                           loc='best', labelspacing=1, prop={'size': 18}) 
    else:
        legend = ax.legend(handles= [Patch(facecolor='green', edgecolor='green',label='Cosine sim'),
                        Patch(facecolor='red', edgecolor='red',label='Euclidean dist'),
                        Patch(facecolor='blue', edgecolor='blue',label='Manhattan dist'),
                        Line2D([0], [0], marker='o', color='b', label='C0PD0', markerfacecolor='w', markersize=12), 
                        Line2D([0], [0], marker='s', color='b', label='C6PD6', markerfacecolor='w', markersize=12), 
                        Line2D([0], [0], marker='P', color='b', label='C15PD15', markerfacecolor='w', markersize=12),
                        Line2D([0], [0], marker='*', color='b', label='C21PD21', markerfacecolor='w', markersize=12)],

                       loc='best', labelspacing=1, prop={'size': 18}) 
    legend.get_frame().set_facecolor('white')
    if not os.path.exists('output/Precision_Recall'):
        os.makedirs('output/Precision_Recall') 
    plt.savefig(f'output/Precision_Recall/{PR_RL}_{condit}', dpi = 600)
    plt.show()
    


Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)                       

in_dir = 'input/Sim_measure'
    
sim_measures = os.listdir(in_dir)
ccs2compare_all = os.listdir(f'{in_dir}/{sim_measures[0]}')
#### direct comparison of similarities between pairs of control and pink1

PDvalidations_pos = {key:[] for key in sim_measures}
for sim_measure in sim_measures:
    print(sim_measure)
    for ccs2compare in ccs2compare_all:
        ccs2compare = ccs2compare.split('.')[0]
        print(ccs2compare)
        save_dir = f'output/AllPredictions/{sim_measure}'
        
        PD_cond = ccs2compare.split('-')[1]
        in_dir = 'input/PD_preds'
        PD_preds = []
        for root, dirs, files in os.walk(in_dir):
            if PD_cond in root:
                for file in files:
                    with open(f'{root}/{file}') as f:
                        for line in f:
                            if '\t' not in line:
                                PD_preds.append(line[:-1])
                            else:
                                lspt = line.split('\t')[0]
                                PD_preds.append(lspt)

        gene2gene_df, cc1cc2_genes, top20genes_cc1cc2 = compare2ccs(ccs2compare, PD_preds, distance_meas = sim_measure, save_dir = save_dir)
        with open(f'{save_dir}/{legend[ccs2compare]}_PD_preds.pkl', 'wb') as handel:
            pickle.dump(cc1cc2_genes, handel)

        savePDPreds(cc1cc2_genes, f'{legend[ccs2compare]}_PDpreds', save_dir = save_dir)
        
        PD_pos = []
        for gene in cc1cc2_genes:
            if LitValid_AllGenes[gene][0] > 0:
                PD_pos.append(cc1cc2_genes.index(gene))      
            elif LitValid_AllGenes[gene][1] != '':
                PD_pos.append(cc1cc2_genes.index(gene))      
            elif LitValid_AllGenes[gene][2] != '':
                PD_pos.append(cc1cc2_genes.index(gene))      
        PDvalidations_pos[sim_measure].append([ccs2compare, PD_pos])

        # PD_pos = []
        # for gene in cc1cc2_genes:
        #     if gene in Gene4PD or gene in GWASdb_SNP:
        #         PD_pos.append(cc1cc2_genes.index(gene))      
        # PDvalidations_pos[sim_measure].append([ccs2compare, PD_pos])


        ### check if distances are significant compared to background
        PDGs2Background_dist(ccs2compare, sim_measure, cc1cc2_genes)
        ### check if PD known gene distances are significant compared to background
        with open('input/PD_genes_DGN_DEGs.pkl', 'rb') as handle:
            DisGeNet_PDgenes = pickle.load(handle)  
        PDGs2Background_dist(ccs2compare, sim_measure, DisGeNet_PDgenes, PD_Gs = 'PD genes')
    


### finding out the best distance measure -> cosine #make a statistic based on precision and recall, like Marinka
# precision - number of positive outcomes among top k-rated genes
# recall - fraction of all positive outcomes among the top k-rated genes
ks = [10,20, 30, 40, 50, 60, 70, 80, 90, 100,150]

Precision_pd = pd.DataFrame(index=legend.keys(), columns = PDvalidations_pos.keys())
Recall_pd = pd.DataFrame(index=legend.keys(), columns = PDvalidations_pos.keys())
for col in Precision_pd.columns:
    Precision_pd[col] = [[] for i in repeat(None, len(Precision_pd))]
for col in Recall_pd.columns:
    Recall_pd[col] = [[] for i in repeat(None, len(Recall_pd))]

for k in ks:
    for sim_measure in PDvalidations_pos.keys():
        for lst_ccpair in PDvalidations_pos[sim_measure]:  
            cc_pair = lst_ccpair[0]
            PO = 0
            for gene_pos in lst_ccpair[1]:
                if gene_pos <= k:
                    PO+=1
            precision = PO/k
            recall = PO/len(lst_ccpair[1])
            Precision_pd.loc[cc_pair, sim_measure].append(precision)
            Recall_pd.loc[cc_pair, sim_measure].append(recall)


Plot_precorrecall(Precision_pd, 'Precision')
Plot_precorrecall(Recall_pd, 'Recall')
# Plot_precorrecall(Precision_pd, 'Precision', condit = 'Control_D21-PINK1_D21')
# Plot_precorrecall(Recall_pd, 'Recall', condit = 'Control_D21-PINK1_D21')

Precision_pd_noCos = Precision_pd.drop('Cosine_similarity', axis=1)  
Recall_pd_noCos = Recall_pd.drop('Cosine_similarity', axis=1)
Plot_precorrecall(Precision_pd_noCos, 'Precision')
Plot_precorrecall(Recall_pd_noCos, 'Recall')

Plot_precorrecall(Precision_pd_noCos, 'Precision', condit = 'Control_D21-PINK1_D21')
Plot_precorrecall(Recall_pd_noCos, 'Recall', condit = 'Control_D21-PINK1_D21')
Plot_precorrecall(Precision_pd_noCos, 'Precision', condit = 'Control_D15-PINK1_D15')
Plot_precorrecall(Recall_pd_noCos, 'Recall', condit = 'Control_D15-PINK1_D15')
Plot_precorrecall(Precision_pd_noCos, 'Precision', condit = 'Control_D06-PINK1_D06')
Plot_precorrecall(Recall_pd_noCos, 'Recall', condit = 'Control_D06-PINK1_D06')
Plot_precorrecall(Precision_pd_noCos, 'Precision', condit = 'Control_IPSCs-PINK1_IPSCs')
Plot_precorrecall(Recall_pd_noCos, 'Recall', condit = 'Control_IPSCs-PINK1_IPSCs')



### Core pd preds - genes between all pairs of control and Pink1
save_dir_SG = 'output/SharedGenes'
PD_preds_ccs = []
PD_preds_ccs_dict = {}
for root, dirs, files in os.walk('output/AllPredictions/Euclidean_distance'): 
    for file in files:
        if file.endswith('.pkl'):
            pair = file.split('_')[0]
            print(pair)
            with open(f'{root}/{file}', 'rb') as handle:
                cc_genes = pickle.load(handle)  
            PD_preds_ccs.append(cc_genes)
            PD_preds_ccs_dict[pair] = cc_genes
    common_genes = set.intersection(*map(set,PD_preds_ccs))
    print(len(common_genes))

   
# is the number of Core pd preds better than random  
reps = 10000
Successes = 0 
for rep in range(reps):
    PD_preds_ccs_rand = []  
    for pair in PD_preds_ccs_dict.keys():
        with open(f'input/SortedGenes/{legend_r[pair]}.pkl', 'rb') as handle:
            cc1cc2_SGs = pickle.load(handle) 
        cc1cc2_SGs = [Entrez2Sym[str(x)] for x in cc1cc2_SGs]
        n_samp = len(PD_preds_ccs_dict[pair])
        lst = sample(cc1cc2_SGs, n_samp)
        PD_preds_ccs_rand.append(lst)
    common_genes_rand = set.intersection(*map(set,PD_preds_ccs_rand))
    if len(common_genes_rand) >= len(common_genes):
        Successes+=1
pval = (Successes + 1)/(reps+1) 
print(pval)

# sort Common_genes by the avg distance across all stage pairs - from most distant to least
cc_pairs_noIPSC = list(legend.keys())[:-1]
CGs_dist_pd = pd.DataFrame(index= common_genes, columns = legend.keys())
CGs_dist_pd_noIPSC = pd.DataFrame(index= common_genes, columns = cc_pairs_noIPSC)
for sim_measure in sim_measures:
    if sim_measure == 'Euclidean_distance':
        print(sim_measure)
        for ccs2compare in ccs2compare_all:
            ccs2compare = ccs2compare.split('.')[0]
            print(ccs2compare)
            sim_dist_df = create_ccpair_df(sim_measure, ccs2compare)
            for gene in common_genes:
                CGs_dist_pd.loc[gene,ccs2compare] = sim_dist_df.loc[gene,0]
                if ccs2compare != 'Control_IPSCs-PINK1_IPSCs':
                    CGs_dist_pd_noIPSC.loc[gene,ccs2compare] = sim_dist_df.loc[gene,0]
                    
CGs_dist_pd['average'] = CGs_dist_pd.mean(axis=1)
CGs_dist_pd = CGs_dist_pd.sort_values(by='average', ascending=False)
CGs_dist_pd_noIPSC['average'] = CGs_dist_pd_noIPSC.mean(axis=1)
CGs_dist_pd_noIPSC = CGs_dist_pd_noIPSC.sort_values(by='average', ascending=False)

common_genes_ordered = CGs_dist_pd.index.tolist()
# with open('output/CommonGenes_BetwStagePairs.txt', 'w') as f:
#     for gene in common_genes_ordered:
#         f.write(f'{gene}\n')
with open('output/CorePreds_BetwStagePairs.pkl', 'wb') as fp:   
    pickle.dump(common_genes_ordered, fp)
savePDPreds(common_genes_ordered, 'CorePreds_BetwStagePairs_Valid', save_dir = 'output')

# ### ALL genes statistics
# =============================================================================
# PD21inC21_geneDist = []
# for gene in C21PD21_genes:
#     PD21inC21_geneDist.append(PD21inC21_geneDist_dic[gene][1])
# print(len(PD21inC21_geneDist))

# #statistical tests
# statistics,pvalue = kstest(PD21inC21_geneDist,'norm')       
# plt.hist(PD21inC21_geneDist,bins = 20)
# statistics_all,pvalue_all = kstest(Dist_gene2gene,'norm')       
# plt.hist(Dist_gene2gene,bins = 50)

# statistics,pvalue = ttest_ind(PD21inC21_geneDist,Dist_gene2gene)
# stat,pval = mannwhitneyu(PD21inC21_geneDist,Dist_gene2gene)
# =============================================================================


'''
# permutation test
reps = 100
Successes = 0
N = len(C21PD21_genes)
for rep in range(reps):
    RandDists = sample(Dist_gene2gene,N)
    _,pval_tt = mannwhitneyu(PD21inC21_geneDist,RandDists)
    if pval_tt<0.05:
        Successes += 1
print(Successes)
#pval = (Successes + 1)/(reps+1) 
#if pval < 0.05:
    #print('yay')

# enrich 
M = len(Geneslist)
K = len(Thresh_genes)
N = len(C21PD21_genes)
X = len(PD_pred_thresh)

fold = (X/N)/(K/M) #how many of background genes compared to cluster is pd genes - just percentage
if fold >= 1:
    pval = hypergeom.sf(X-1, M, K, N) #background
    if pval < 0.05:   
        print('great')  

'''
    