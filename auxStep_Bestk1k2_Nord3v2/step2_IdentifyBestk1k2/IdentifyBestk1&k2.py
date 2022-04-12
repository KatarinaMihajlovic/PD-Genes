# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:58:07 2021

@author: kmihajlo
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from kneed import KneeLocator


in_dir = 'input'
out_dir = 'output'
cell_conds = os.listdir(in_dir)
columns = ['Net_comb', 'k1', 'k2', 'dispersionCoefG1', 'dispersionCoefG2', 'DispCoeff_avg']

labels_key = {'Control_IPSCs':'C0', 'Control_D06':'C6', 'Control_D15':'C15', 'Control_D21':'C21', 
              'PINK1_IPSCs':'PD0', 'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21'}

for cell_cond in cell_conds:
    print(cell_cond)
    DispCoeff_df = pd.DataFrame(columns=columns)
    cnt = 0
    for file in os.listdir(f'{in_dir}/{cell_cond}'):
        #print(file)
        with open(f'{in_dir}/{cell_cond}/{file}') as f:
            lines = f.readlines()
            if lines:
                lines.pop(0)
                for line in lines:
                    lspt = line.split('\t')
                    Net_comb = lspt[0]
                    k1 = lspt[1]
                    k2 = lspt[2]
                    dispersionCoefG1 = lspt[3]
                    dispersionCoefG2 = lspt[4]
                    DispCoeff_avg = (float(dispersionCoefG1) + float(dispersionCoefG2[:-1]))/2
                    DispCoeff_df.loc[cnt] = [Net_comb, k1, k2, dispersionCoefG1, dispersionCoefG2, DispCoeff_avg] 
                    cnt+=1
    
    DispCoeff_df = DispCoeff_df.sort_values(by = 'DispCoeff_avg', ascending=False)  
    # print(DispCoeff_df)      
    if not os.path.exists(f'{out_dir}/{cell_cond}'):
        os.makedirs(f'{out_dir}/{cell_cond}')
    DispCoeff_df.to_csv(f'{out_dir}/{cell_cond}/DispCoeff_{cell_cond}.csv')
    DispCoeff_df.to_csv(f'{out_dir}/AllHMs/DispCoeff_{cell_cond}.csv')
    DispCoeff_df.head(10).to_csv(f'{out_dir}/{cell_cond}/DispCoeff_{cell_cond}_top10.csv')
    
    #find elbow point
    y = list(DispCoeff_df['DispCoeff_avg'])
    x = range(1, len(y)+1)
    kn = KneeLocator(x, y, curve='concave', direction='decreasing')
    elbow = kn.knee
    print(DispCoeff_df.iloc[elbow]['DispCoeff_avg'])
    print(DispCoeff_df.iloc[elbow]['k1'])
    print(DispCoeff_df.iloc[elbow]['k2'])
    dc_elbow = DispCoeff_df.iloc[elbow]['DispCoeff_avg']
    # from kneed import KneeLocator, DataGenerator as dg
    # x, y = dg.concave_increasing()
    # kl = KneeLocator(x, y, curve="concave", direction="increasing")
    # kn.plot_knee()

    #plot DispCoeff_avg
    fig, ax = plt.subplots(figsize=(11, 9)) 
    DispCoeff_avglist = DispCoeff_df.loc[:,'DispCoeff_avg'].tolist()
    DispCoeff_avglist.reverse()
    counts = []
    for count, value in enumerate(DispCoeff_avglist):
        counts.append(count)
    plt.plot(counts, DispCoeff_avglist)
    plt.axhline(y=dc_elbow, color='r', linestyle='-')
    plt.title(labels_key[cell_cond], fontsize=20)
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeff_{cell_cond}',dpi = 600)
    plt.show()
    plt.close()
    
    # extract values for ALL comb and sort for k1 and k2
    DispCoeffALL_df = DispCoeff_df.loc[DispCoeff_df['Net_comb'] == 'ALL']
    DispCoeffALL_df = DispCoeffALL_df.sort_values(by = 'DispCoeff_avg', ascending=False)  
    
    DispCoeffALL_df_k1sorted = DispCoeffALL_df.sort_values(by = 'k1', ascending=False)  
    DispCoeffALL_df_k1sorted['k1'] = DispCoeffALL_df['k1'].astype(int)
    DispCoeffALL_df_k1sorted = DispCoeffALL_df_k1sorted.sort_values(by = 'k1', ascending=True)  
    DispCoeffALL_df_k2sorted = DispCoeffALL_df.sort_values(by = 'k2', ascending=False)  
    DispCoeffALL_df_k2sorted['k2'] = DispCoeffALL_df['k2'].astype(int)
    DispCoeffALL_df_k2sorted = DispCoeffALL_df_k2sorted.sort_values(by = 'k2', ascending=True)     
    
    #plot DispCoeff of ALL vs k1 values
    fig, ax = plt.subplots(figsize=(11, 9))    
    DispCoeff_avglistk1 = DispCoeffALL_df_k1sorted.loc[:,'DispCoeff_avg'].tolist()
    k1s = DispCoeffALL_df_k1sorted.loc[:,'k1'].tolist()   
    plt.plot(k1s, DispCoeff_avglistk1, marker='o')
    plt.title(labels_key[cell_cond], fontsize = 28)
    plt.xlabel('k1', fontsize = 22)
    plt.ylabel('Dispersion Coefficient', fontsize = 22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(set(k1s))
    plt.axhline(y=dc_elbow, color='r', linestyle='-')
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeffvsk1_{cell_cond}', dpi = 600)
    plt.show()
    plt.close()
    
    #plot DispCoeff of ALL vs max k1 values    
    k1s_set = list(set(k1s))
    k1s_set.sort()
    max_DispCoeff_avg = []
    for k1 in k1s_set:
        max_DispCoeff = DispCoeffALL_df_k1sorted[DispCoeffALL_df_k1sorted['k1']==k1]['DispCoeff_avg'].max()
        max_DispCoeff_avg.append(max_DispCoeff)
    
    fig, ax = plt.subplots(figsize=(11, 9))    
    plt.plot(k1s_set, max_DispCoeff_avg, marker='o')
    plt.title(labels_key[cell_cond], fontsize = 28)
    plt.xlabel('k1', fontsize = 22)
    plt.ylabel('Dispersion Coefficient', fontsize = 22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(set(k1s))
    plt.savefig(f'{out_dir}/{cell_cond}/Max_DispCoeffvsk1_{cell_cond}', dpi = 600)
    plt.show()
    plt.close()    
    
    #plot DispCoeff of ALL vs k2 values
    fig, ax = plt.subplots(figsize=(11, 9)) 
    DispCoeff_avglistk2 = DispCoeffALL_df_k2sorted.loc[:,'DispCoeff_avg'].tolist()
    k2s = DispCoeffALL_df_k2sorted.loc[:,'k2'].tolist()
    plt.plot(k2s, DispCoeff_avglistk2, marker='o')
    plt.title(labels_key[cell_cond], fontsize = 28)
    plt.xlabel('k2', fontsize = 22)
    plt.ylabel('Dispersion Coefficient', fontsize = 22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(set(k2s))
    plt.axhline(y=dc_elbow, color='r', linestyle='-')
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeffvsk2_{cell_cond}', dpi = 600)
    plt.show()
    plt.close()
    
    #plot DispCoeff of ALL vs max k2 values    
    k2s_set = list(set(k2s))
    k2s_set.sort()
    max_DispCoeff_avg = []
    for k2 in k2s_set:
        max_DispCoeff = DispCoeffALL_df_k2sorted[DispCoeffALL_df_k2sorted['k2']==k2]['DispCoeff_avg'].max()
        max_DispCoeff_avg.append(max_DispCoeff)
    
    fig, ax = plt.subplots(figsize=(11, 9))    
    plt.plot(k2s_set, max_DispCoeff_avg, marker='o')
    plt.title(labels_key[cell_cond], fontsize = 28)
    plt.xlabel('k2', fontsize = 22)
    plt.ylabel('Dispersion Coefficient', fontsize = 22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.grid(set(k2s))
    plt.savefig(f'{out_dir}/{cell_cond}/Max_DispCoeffvsk2_{cell_cond}', dpi = 600)
    plt.show()
    plt.close()   
    
    #heatmap of DispCoeff with k1 and k2 on axes
    k1s = set(k1s)
    k2s = set(k2s)
    DispCoeff_df_k1k2 = pd.DataFrame(columns=k2s,index=k1s)

    for k1 in k1s:
        for k2 in k2s:
            try:
                DispCoeff = DispCoeffALL_df.loc[(DispCoeffALL_df['k1'] == str(k1)) & (DispCoeffALL_df['k2'] == str(k2)), 'DispCoeff_avg'].values[0]
                DispCoeff_df_k1k2.loc[k1,k2] = DispCoeff
            except IndexError:
                DispCoeff_df_k1k2.loc[k1,k2] = None
    
    DispCoeff_df_k1k2 = DispCoeff_df_k1k2.sort_index(axis=1)
    DispCoeff_df_k1k2 = DispCoeff_df_k1k2.sort_index(axis=0)
    DispCoeff_df_k1k2.fillna(value=np.nan, inplace=True)
    fig, ax = plt.subplots(figsize=(11, 9)) 
    heatmap_DC = sb.heatmap(DispCoeff_df_k1k2, cmap='Spectral', annot=True, fmt='.3g', annot_kws={"fontsize":9})
    ax.set_title(f'Dispersion coefficient - {labels_key[cell_cond]}', fontsize = 22)
    ax.set_xlabel('k2', fontsize = 20)
    ax.set_ylabel('k1', fontsize = 20)
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeffHM_{cell_cond}', dpi = 600)
    plt.savefig(f'{out_dir}/AllHMs/DispCoeffHM_{cell_cond}', dpi = 600)