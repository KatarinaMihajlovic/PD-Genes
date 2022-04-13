# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:58:07 2021

@author: kmihajlo
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np


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
    DispCoeff_df.head(10).to_csv(f'{out_dir}/{cell_cond}/DispCoeff_{cell_cond}_top10.csv')
    
    
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
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeffvsk1_{cell_cond}', dpi = 600)
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
    plt.savefig(f'{out_dir}/{cell_cond}/DispCoeffvsk2_{cell_cond}', dpi = 600)
    plt.show()
    plt.close()
    
