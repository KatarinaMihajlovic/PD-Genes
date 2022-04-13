# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:12:57 2021

@author: kmihajlo
"""
import os, requests, re, math
import pickle, sys
from scipy.stats import hypergeom, mannwhitneyu
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from random import sample


def Genesym2Ensembl(filename = 'input/Homo_sapiens.gene_info'):
    Genesym2Ensembl = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            lspt = line.split('\t')
            GeneNames = lspt[5].split('|')
            for Genename in GeneNames:
                if 'Ensembl' in Genename:
                    #print(Genename)
                    Ensembl = Genename.split(':')[1]
                    Genesym2Ensembl[lspt[2]] = Ensembl
    return Genesym2Ensembl

Genesym2Ensembl = Genesym2Ensembl()
with open('input/Gene4PD.pkl', 'rb') as handle:
    Gene4PD = pickle.load(handle)  
with open('input/GWASdb_SNP.pkl', 'rb') as handle:
    GWASdb_SNP = pickle.load(handle)  
    
#Gaias code
def extract_single(drug):
    try:
        contains = requests.get("https://pubmed.ncbi.nlm.nih.gov/?term=%28{}%29".format(drug))
        m = re.search('<meta name="log_resultcount" content="[0-9]+"', str(contains.content))
        return int(m.group(0)[38:-1])
    except:
        return 1

def extract_or(*args):
    query = "https://pubmed.ncbi.nlm.nih.gov/?term=%28{}%29".format(args[0].replace(" ", "+"))
    for a in args[1:]:
        query += "+OR+%28{}%29".format(a.replace(" ", "+"))
    #print(query)
    try:
        contains = requests.get(query)
        m = re.search('<meta name="log_resultcount" content="[0-9]+"', str(contains.content))
        return int(m.group(0)[38:-1])
    except:
        return 0
    
def extract_and(drugs, prots):
    drug_sentence = "%28%28{}%29".format(drugs[0].replace(" ", "+"))
    for d in drugs[1:]:
        drug_sentence += "+OR+%28{}%29".format(d.replace(" ", "+"))
    drug_sentence += "%29"
    prot_sentence = "%28%28{}%29".format(prots[0].replace(" ", "+"))
    for d in prots[1:]:
        prot_sentence += "+OR+%28{}%29".format(d.replace(" ", "+"))
    prot_sentence += "%29"
    
    query = "https://pubmed.ncbi.nlm.nih.gov/?term={}+AND+{}".format(drug_sentence,prot_sentence)
    #print(query)
    contains = requests.get(query)
    
    m = re.search('<meta name="log_resultcount" content="[0-9]+"', str(contains.content))
    if m is not None:
        res = int(m.group(0)[38:-1])
    else:
        res = 1
    #print(res)
    return res

def Pubmed_genes(cc1cc2_PD_preds, search_term = 'Parkinson\'s disease'):
    cc1cc2_PDpreds_LitValid = {}
    for gene in cc1cc2_PD_preds:
        numPublications1 = extract_and([gene], [search_term])
        try:
            numPublications2 = extract_and([Genesym2Ensembl[gene]], [search_term])
        except KeyError:
            numPublications2 = 0            
        numPublications = numPublications1 + numPublications2
        print(gene, numPublications)
        if gene in Gene4PD:
            in_Gene4PD = 'Gene4PD'
        else:
            in_Gene4PD = ''
        if gene in GWASdb_SNP:
            in_GWASdb_SNP = 'GWASdb_SNP'
        else:
            in_GWASdb_SNP = ''
        cc1cc2_PDpreds_LitValid[gene] = [numPublications, in_Gene4PD, in_GWASdb_SNP]
    return cc1cc2_PDpreds_LitValid

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
 
def nvalid_genes(genes_LitValid):
    n_valids = 0
    for gene in genes_LitValid:
        if genes_LitValid[gene][0] > 0:
            n_valids +=1
        elif genes_LitValid[gene][1] != '':
            n_valids +=1
        elif genes_LitValid[gene][2] != '':
            n_valids +=1
    return n_valids

def BarPlotValidGenes(PDpreds_enrich, file_out, filens):
    ### plot DEGs enrichment of pairwise PD preds
    cc_pairs = [x[0] for x in PDpreds_enrich]
    PD_preds = [x[1] for x in PDpreds_enrich]
    ValidGenes = [x[2] for x in PDpreds_enrich]
    p_vals =  [x[3] for x in PDpreds_enrich]
    cc_pairs.insert(1, cc_pairs[-1])
    cc_pairs.pop(-1)
    PD_preds.insert(1, PD_preds[-1])
    PD_preds.pop(-1)
    ValidGenes.insert(1, ValidGenes[-1])
    ValidGenes.pop(-1)
    p_vals.insert(1, p_vals[-1])
    p_vals.pop(-1)
    
    width = 0.3    
    opacity = 0.7
    
    fig, ax = plt.subplots(figsize=(16, 9))
    # Remove top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # Remove y-axis tick marks
    ax.yaxis.set_ticks_position('none')
    # Set plot title
    
    if file_out == filens[0]:
        ax.set_title('Enrichment of PD preds in LitValid genes', fontsize = 28)
    elif file_out == filens[1]:
        ax.set_title('Enrichment of PD preds in DEGs', fontsize = 28)
    # Add major gridlines in the y-axis
    ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
    ax.set_facecolor('xkcd:white')
    
    bar1 = ax.bar(np.arange(len(PD_preds)), PD_preds, width, align='center', alpha=opacity, color='r', label='PD preds')
    if file_out == filens[0]:
        bar2 = ax.bar(np.arange(len(ValidGenes)) + width, ValidGenes, width, align='center', alpha=opacity,  color='limegreen', label='LitValid PD preds')
    elif file_out == filens[1]:
        bar2 = ax.bar(np.arange(len(ValidGenes)) + width, ValidGenes, width, align='center', alpha=opacity,  color='limegreen', label='DEGs PD preds')
     
    # for i in range(N):
    #     rects1 = ax.bar(ind[i-1], PD_preds[i-1], width, color='red')
    #     rects2 = ax.bar(ind[i-1]+width, ValidGenes[i-1], width, color='limegreen')
        # fig.text(ind+width/2, 0.9, p_vals[i-1], fontsize=16)
    # ax.legend( (bar1[0], bar2[0]), ('PD preds', 'DEGs PD preds'),  fontsize = 18)
    ax.legend(fontsize=22, loc = 'upper left')
    ax.set_ylabel('#genes', fontsize = 24, fontweight = 'bold')
    ax.set_xticks(np.arange(len(PD_preds))+width/2)
    ax.set_xticklabels(cc_pairs,  fontsize = 24) 
    ax.tick_params(axis='y', which='major', labelsize=24)
    plt.ylim(0,1350)
    i = 0
    for rect in bar1:
        height = rect.get_height()
        if float(p_vals[i]) < 0.05:
            style = 'bold'
        else: 
            style = 'normal'
        plt.text(rect.get_x() + rect.get_width(), height+0.1, f'{p_vals[i]}', ha='center', va='bottom', fontsize = 24, fontweight = style)
        i+=1
    
    if file_out == filens[0]:
        plt.savefig('output/StageSpecPreds/LitValid_Enrich.png', dpi = 600, bbox_inches='tight')	
    elif file_out == filens[1]:
        plt.savefig('output/StageSpecPreds/StageSpecDEGs_Enirch.png', dpi = 600, bbox_inches='tight')	
    plt.show()
    plt.close()


def Entrez2Symbols(gene_info): #dict with entrez as keys and symbol as values (key:value)
    Entrez2Sym = {}
    gene_info.readline()
    for line in gene_info.readlines():
        lspt=line.strip().split('\t')
        Symbol = lspt[2]

        Entrez2Sym[lspt[1]] = Symbol
    return(Entrez2Sym)

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
    print('\n')   
        

####### MAIN CODE
day_dict = {'0':'IPSCs', '6':'D06', '15':'D15', '21':'D21'}

Flag1 = str(sys.argv[1])
Flag2 = str(sys.argv[2])

with open('input/PD_genes_DGN_DEGs.pkl', 'rb') as handle:
    PD_genes = pickle.load(handle) 
# remove PD genes from background to make it consistent with the background from where you choose the predictions

# All_litValidGenes = {}
    
### Stage Specific validation 
filens = ['LitValid_Enirch.txt', 'StageSpecDEGs_Enirch.txt']
for file_out in filens: 
    PDpreds_enrich = []
    if not os.path.exists('output/StageSpecPreds'):
        os.makedirs('output/StageSpecPreds')  
    outf = open(f'output/StageSpecPreds/{file_out}','w')
    
    for root, dirs, files in os.walk('input/StageSpecPreds'):
        for file in files:
            if file.endswith('.pkl'):
                print(root, file)
                sd = root.split('/')[1]
                filen = file.split('.')[0]
                cc1cc2_PDpreds_LitValid = {}
                save_dir = f'output/{sd}'
                day = filen.split('C')[1].split('P')[0]
                
                with open(f'input/StageSpecDEGs_fold0_1/Control_{day_dict[day]}_PINK1_{day_dict[day]}.txt') as f:
                    DEGs = f.readlines()
                DEGs = [x[:-1] for x in DEGs]
                
                if Flag1:
                    with open(f'{root}/{file}', 'rb') as handel:
                        cc1cc2_PD_preds = pickle.load(handel) #StageSPec Preds
                    cc1cc2_PDpreds_LitValid = Pubmed_genes(cc1cc2_PD_preds)
                    if not os.path.exists(save_dir):
                        os.makedirs(save_dir)  
                    with open(f'{save_dir}/{filen}.txt', 'w') as f:
                        for gene in cc1cc2_PDpreds_LitValid.keys():
                            f.write(f'{gene}\t{cc1cc2_PDpreds_LitValid[gene][0]}\t{cc1cc2_PDpreds_LitValid[gene][1]}\t{cc1cc2_PDpreds_LitValid[gene][2]}\n')
                    with open(f'{save_dir}/{filen}.pkl', 'wb') as fp:   
                        pickle.dump(cc1cc2_PDpreds_LitValid, fp)
                    
                else:
                    with open(f'{save_dir}/{filen}.pkl', 'rb') as handle:
                        cc1cc2_PDpreds_LitValid = pickle.load(handle) 
                    cc1cc2_PD_preds = cc1cc2_PDpreds_LitValid.keys()
                # All_litValidGenes = {**All_litValidGenes , **cc1cc2_PDpreds_LitValid}
                    
                #check with background
                
                save_dir2 = 'output/Other_genes'
                cc_pair = file.split('_')[0]
                if Flag2:               
                    with open(f'input/All_genes/{cc_pair}_AllGenes.pkl', 'rb') as handel:
                        cc1cc2_AllGenes = pickle.load(handel)                
                    cc1cc2_OtherGenes = [gene for gene in cc1cc2_AllGenes if gene not in cc1cc2_PD_preds]
                    
                    cc1cc2_OtherGenes_LitValid = Pubmed_genes(cc1cc2_OtherGenes)                
                    if not os.path.exists(save_dir2):
                        os.makedirs(save_dir2)        
                    with open(f'{save_dir2}/{cc_pair}_OtherGenes.pkl', 'wb') as fp:   
                        pickle.dump(cc1cc2_OtherGenes_LitValid, fp)
                    
                else:  
                    with open(f'{save_dir2}/{cc_pair}_OtherGenes.pkl', 'rb') as handle:
                        cc1cc2_OtherGenes_LitValid = pickle.load(handle) 
                # All_litValidGenes = {**All_litValidGenes , **cc1cc2_OtherGenes_LitValid}
                
                ##### MWU test - shows if clustered genes have higher number of citations 
                for gene in PD_genes:
                    cc1cc2_OtherGenes_LitValid.pop(gene, None)
                cc1cc2_pubmed = logPubmedn(cc1cc2_PDpreds_LitValid)
                cc1cc2_OtherGenes_pubmed = logPubmedn(cc1cc2_OtherGenes_LitValid)
                PlotHist_MWU(cc1cc2_pubmed, cc1cc2_OtherGenes_pubmed, cc_pair, save_dir)
                
                ### check enrichment of of PD preds with PD proof from 1) Pubmed and Gene4PD, GWASdb; 2) DEGs
                
                outf.write(f'{filen}\n')

                # x = []
                # y = []
                M = 0
                K = 0
                N = 0
                X = 0
                for gene in cc1cc2_PDpreds_LitValid:
                    M+=1
                    N+=1
                    # x.append(N)
                    if file_out == filens[0]:
                        if cc1cc2_PDpreds_LitValid[gene][0] > 0 or cc1cc2_PDpreds_LitValid[gene][1] != '' or cc1cc2_PDpreds_LitValid[gene][2] != '':                    
                            K+=1
                            X+=1
                            # y.append(1)
                        # else:
                            # y.append(0)
                    elif file_out == filens[1]:
                        if gene in DEGs:
                            K+=1
                            X+=1
                            # y.append(1)
                        # else:
                            # y.append(0)
                for gene in cc1cc2_OtherGenes_LitValid:
                    M+=1
                    if file_out == filens[0]:
                        if cc1cc2_OtherGenes_LitValid[gene][0] > 0 or cc1cc2_OtherGenes_LitValid[gene][1] != '' or cc1cc2_OtherGenes_LitValid[gene][2] != '':
                            K+=1
                    elif file_out == filens[1]:
                        if gene in DEGs:
                            K+=1
                print(M,K,N,X)
                outf.write(f'{M}\t{K}\t{N}\t{X}\n')
                
                try:
                    fold = (X/N)/(K/M)
                    print(fold)
                except ZeroDivisionError:
                    fold = 0
                outf.write(f'{fold}\n')
    
                if fold >= 1:
                    pval = hypergeom.sf(X-1, M, K, N)
                    print(pval)
                    outf.write(f'{pval}\n')
                    if pval <= 0.05: 
                        #print(f'Enriched clust {i},  fold = {fold},  pval = {pval}')
                        #print(N, X)
                        if file_out == filens[0]:
                            print('Enriched in LitValid - PubMed,Gene4PD,GWASdb')
                            outf.write('Enriched in LitValid - PubMed,Gene4PD,GWASdb\n')
                        elif file_out == filens[1]:    
                            print('Enriched in StageSpec DEGs')
                            outf.write('Enriched in StageSpec DEGs\n')
                else:
                     pval = 1
                pval = "{:.3e}".format(pval)
                PDpreds_enrich.append([filen.split('_')[0], N, X, pval])
                print(f'percentage of enriched genes: {X/N*100}')
                print('\n')  
                # plt.bar(x,y)
    outf.close()

    
    ### plot Validated Genes and p values of enrichments of pairwise PD preds
    BarPlotValidGenes(PDpreds_enrich, file_out, filens)
    Flag1 = False
    Flag2 = False

######### Core PD predictions Validation
save_dir = 'output/CorePreds'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)  

# Create Background set of Other genes; Genes expressed across all conditinos - includes PD genes from DisGeNet

import csv


Flag = str(sys.argv[3])
if Flag:
    Entrez2Symbol_file = open('input/Homo_sapiens.gene_info', 'r')
    Entrez2Sym = Entrez2Symbols(Entrez2Symbol_file)

    all_genes = []
    for root, dirs, files in os.walk('input/All_genes'):
        for file in files:
            if file.endswith('.csv') and 'Geneslist' in file:
                cc_genes = []
                with open(f'{root}/{file}') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter='\n')
                    next(csv_reader)                
                    for row in csv_reader:
                        cc_genes.append(Entrez2Sym[row[0]])
                all_genes.append(cc_genes)
        
    common_genes = set.intersection(*map(set,all_genes)) #genes expressed across all cell conditions

    # with open('input/PD_genes_DGN_DEGs.pkl', 'rb') as handle:
    #     PD_genes = pickle.load(handle)
    # cgs_noPDgenes = [x for x in common_genes if x not in PD_genes]

    with open('input/CorePreds_BetwStagePairs.pkl', 'rb') as handel:
        CGs_PD_preds = pickle.load(handel)
    Other_CommonGenes = [x for x in common_genes if x not in CGs_PD_preds] #Genes expressed across all conditinos - includes PD genes from DisGeNet

    CGs_PD_preds_LitValid = Pubmed_genes(CGs_PD_preds)            
    with open(f'{save_dir}/CorePreds_LitValid.txt', 'w') as f:
        for gene in CGs_PD_preds_LitValid.keys():
            f.write(f'{gene}\t{CGs_PD_preds_LitValid[gene][0]}\t{CGs_PD_preds_LitValid[gene][1]}\t{CGs_PD_preds_LitValid[gene][2]}\n')
    with open(f'{save_dir}/CorePreds_LitValid.pkl', 'wb') as fp:   
        pickle.dump(CGs_PD_preds_LitValid, fp)
        
    Other_CommonGenes_LitValid = Pubmed_genes(Other_CommonGenes)            
    with open(f'{save_dir}/Other_CommonGenes_LitValid.pkl', 'wb') as fp:   
        pickle.dump(Other_CommonGenes_LitValid, fp)
        
    All_CommonGenes_LitValid = {**CGs_PD_preds_LitValid,**Other_CommonGenes_LitValid}
    with open(f'{save_dir}/All_CommonGenes_LitValid.pkl', 'wb') as fp:   
        pickle.dump(All_CommonGenes_LitValid, fp)
else:
    with open(f'{save_dir}/CorePreds_LitValid.pkl', 'rb') as handle:
        CGs_PD_preds_LitValid = pickle.load(handle)     
    with open(f'{save_dir}/Other_CommonGenes_LitValid.pkl', 'rb') as handle:
        Other_CommonGenes_LitValid = pickle.load(handle)     
    with open(f'{save_dir}/All_CommonGenes_LitValid.pkl', 'rb') as handle:
        All_CommonGenes_LitValid = pickle.load(handle)        


for gene in PD_genes:
    All_CommonGenes_LitValid.pop(gene, None)

for gene in PD_genes:
    Other_CommonGenes_LitValid.pop(gene, None)
        
CGs_PD_preds_pubmed = logPubmedn(CGs_PD_preds_LitValid)
OtherGenes_pubmed = logPubmedn(Other_CommonGenes_LitValid) #all genes a OtherGenes_LitValid
PlotHist_MWU(CGs_PD_preds_pubmed, OtherGenes_pubmed, 'CorePreds', save_dir)


### check enrichment of PD preds with PD proof from Pubmed and Gene4PD, GWASdb
print('Core Preds')
file_out = 'CorePreds_LitValid_Enirch.txt'
outf = open(f'{save_dir}/{file_out}','w')
LitEnrich(All_CommonGenes_LitValid, CGs_PD_preds_LitValid, outf)


 # check with totally random, Sampling with replacement 
reps = 10000
Successes = 0 
Othergenes = list(Other_CommonGenes_LitValid.keys()) 
n_samp = len(CGs_PD_preds_LitValid)
n_valids = nvalid_genes(CGs_PD_preds_LitValid)        
for rep in range(reps):
    lst = sample(Othergenes, n_samp) 
    n_valids_r = 0       
    for gene in lst:
        if Other_CommonGenes_LitValid[gene][0] > 0 or Other_CommonGenes_LitValid[gene][1] != '' or Other_CommonGenes_LitValid[gene][2] != '': 
            n_valids_r +=1
    #print(n_valids_r)
    if n_valids_r > n_valids:
        Successes+=1
pval = (Successes + 1)/(reps+1) 
print(pval)
outf.write(f'{pval} - Sampling with replacement {reps} reps')            
outf.close()    

'''
### save StageSPec genes unique to each set of SS preds, i.e. not in Core predictions
save_dir = 'output/StageSpecPreds/Unique_StageSpecPreds'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)  
    
file_out = 'Unique_StageSpecPreds_LitValid_Enirch.txt'
outf = open(f'{save_dir}/{file_out}','w')
    
for root, dirs, files in os.walk('output/StageSpecPreds'):
    for file in files:
        if file.endswith('.pkl') and 'Unique' not in file:
            print(file)
            outf.write(f'{file}\n')
            cc_pair = file.split('_')[0]
            with open(f'output/StageSpecPreds/{file}', 'rb') as handle:
                StageSpecPreds = pickle.load(handle)   
            Unique_StageSpecPreds = StageSpecPreds.copy()
            for gene in CGs_PD_preds_LitValid:
                Unique_StageSpecPreds.pop(gene, None)
            with open(f'{save_dir}/Unique_{file}', 'wb') as fp:   
                pickle.dump(Unique_StageSpecPreds, fp)
                
            with open(f'output/Other_genes/{cc_pair}_OtherGenes.pkl', 'rb') as handel:
                cc1cc2_OtherGenes = pickle.load(handel)                
            
            cc1cc2_AllGenes = {**cc1cc2_OtherGenes , **StageSpecPreds}#, **CGs_PD_preds_LitValid}           
            LitEnrich(cc1cc2_AllGenes, Unique_StageSpecPreds, outf)

            Unique_StageSpecPreds_pubmed = logPubmedn(Unique_StageSpecPreds)
            cc1cc2_OtherGenes = cc1cc2_AllGenes.copy()
            for gene in Unique_StageSpecPreds:
                cc1cc2_OtherGenes.pop(gene, None)
            cc1cc2_OtherGenes_pubmed = logPubmedn(cc1cc2_OtherGenes) #all genes a OtherGenes_LitValid
            PlotHist_MWU(Unique_StageSpecPreds_pubmed, cc1cc2_OtherGenes_pubmed, f'{cc_pair}_unique', save_dir)
outf.close()    
'''
# work_dir = os.getcwd()
# path_parent = os.path.dirname(os.getcwd())
# os.chdir(path_parent) 
# All_litValidGenes = {**All_litValidGenes , **All_CommonGenes_LitValid}
# with open('CommonData/All_litValidGenes.pkl', 'wb') as fp:   
#     pickle.dump(All_litValidGenes, fp)
# os.chdir(work_dir) 
    
    