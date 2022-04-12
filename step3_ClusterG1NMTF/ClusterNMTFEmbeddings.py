# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:02:34 2021

@author: kmihajlo
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 11:28:29 2021

@author: bscuser
"""
import Matrix_Factorization as mf
import pickle, os
from sklearn_extra.cluster import KMedoids
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd

def kMedoids(M, nodes, n_clusters):
    M = M.tolist()
    nodes2coordinates = dict(zip(nodes, M))  
    Cluster_belonging = {k: [] for k in range(n_clusters)}
    
    kmedoids = KMedoids(n_clusters=n_clusters).fit(M)
    kmedoids_labels = list(kmedoids.labels_)
    
    for cluster_index in range(len(kmedoids_labels)):
        cluster = kmedoids_labels[cluster_index]
        node_coords = M[cluster_index]
        for node, coordinates in nodes2coordinates.items():
            if node_coords == coordinates:
                Cluster_belonging[cluster].append(node)
    
    #print(Cluster_belonging)
    Cluster_belonging_list = []
    for _, values in Cluster_belonging.items():
        Cluster_belonging_list.append(values)
    #print(Cluster_belonging_list)
    return Cluster_belonging_list

def kMeans(M, nodes, n_clusters):
    M = M.tolist()
    nodes2coordinates = dict(zip(nodes, M))  
    Cluster_belonging = {k: [] for k in range(n_clusters)}
    
    kMeans = KMeans(n_clusters=n_clusters, random_state=0).fit(M)
    KMeans_labels = list(kMeans.labels_)
    
    for cluster_index in range(len(KMeans_labels)):
        cluster = KMeans_labels[cluster_index]
        node_coords = M[cluster_index]
        for node, coordinates in nodes2coordinates.items():
            if node_coords == coordinates:
                Cluster_belonging[cluster].append(node)
    
    #print(Cluster_belonging)
    Cluster_belonging_list = []
    for _, values in Cluster_belonging.items():
        Cluster_belonging_list.append(values)
    #print(Cluster_belonging_list)
    return Cluster_belonging_list

def Get_Hard_Clusters(M, nodes):
	n,k = np.shape(M)
	Clusters = [[] for i in range(k)]
	for i in range(n):
		idx = np.argmax(M[i])
		Clusters[idx].append(nodes[i])
	return Clusters

def CosSim_kMedoids(M, nodes, n_clusters):
    cos_sim = cosine_similarity(M, dense_output=True)
    print(cos_sim.shape)
    Clusters = kMedoids(cos_sim, nodes, n_clusters)
    return Clusters

cluster_choices = ['kMedoids', 'kMeans', 'HardClusters'] #when number of clusters is the same for all custerings
# cluster_choices = ['kMeans'] #when number of clusters is the same for all custerings

in_dir = 'input' 
out_dir = 'output'
for root, dirs, files in os.walk(in_dir):
    for file in files:
        cell_cond = root.split('\\')[1]
        nets = file.split('_')[0]
        G1_df = pd.read_csv(f'{root}/{file}', header=0, index_col=0, sep='\t')
        G1_NumClusters = len(list(G1_df.columns))
        genes = list(G1_df.index)
        G1 = G1_df.values
        print(cell_cond, nets)
        save_dir = f'{out_dir}/{cell_cond}/{nets}'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir) 
        for cluster_choice in cluster_choices:
            print(cluster_choice)
            if cluster_choice == 'kMedoids': #doesn't work                
                clustersG1_kmed = kMedoids(G1, genes, G1_NumClusters)
                #print(clustersG1_kmed)
                with open(f'{save_dir}/{cluster_choice}_G1.pkl', 'wb') as handle:
                    pickle.dump(clustersG1_kmed, handle)
            elif cluster_choice == 'kMeans':
                clustersG1_kmean = kMeans(G1, genes, G1_NumClusters)
                #print(clustersG1_kmean)
                with open(f'{save_dir}/{cluster_choice}_G1.pkl', 'wb') as handle:
                    pickle.dump(clustersG1_kmean, handle)
            elif cluster_choice == 'HardClusters':
                clustersG1_HC =  mf.Get_Hard_Clusters(G1, genes)
                #print(clustersG1_HC)
                with open(f'{save_dir}/{cluster_choice}_G1.pkl', 'wb') as handle:
                    pickle.dump(clustersG1_HC, handle)
