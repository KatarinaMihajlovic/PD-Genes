import math
from scipy.stats import hypergeom
import numpy as np
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import statistics 

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


# Read GO terms annotations
def load_GO(fname, genes):
	geneset=set(genes)
	GO_counter = {}
	gene2GO = {}
	#loading raw gene & go relationships
	fRead = open(fname, "r")
	for line in fRead.readlines():
		lspt = line.strip().split()
		if len(lspt)>1:
			geneid = lspt[0]
			term = lspt[1]
			if geneid in geneset:
				if term not in GO_counter:
					GO_counter[term] = 0
				GO_counter[term] += 1
				if geneid not in gene2GO:
					gene2GO[geneid] = set()
				gene2GO[geneid].add(term)
	fRead.close()
	#print(GO_counter)
	#filter out GO terms that annotates only one gene
	GO_counter2 = {}
	gene2GO2 = {}
	removed = set()
	for term in GO_counter:
		if GO_counter[term] > 1:
			GO_counter2[term] = GO_counter[term]
		else:
			removed.add(term)
	for gene in gene2GO:
		genego = gene2GO[gene].difference(removed)
		if len(genego)>0:
			gene2GO2[gene]=genego
	return [GO_counter2, gene2GO2]



def go_enrichment(clusters, gene2go, go_counts):
	M = len(gene2go) # total number of annotated genes in the dataset
	data = []
	i = -1
	enrichment = [[] for j in range(len(clusters))]
	cts = 0
	
	clts = []
	gos = []
	pvals = []
	
	NE_list = []
	
	for cluster in clusters:
		i+=1
		annotated_genes = []
		for gene in cluster:
			if gene in gene2go:
				annotated_genes.append(gene)
		N = len(annotated_genes) # number of annotated genes in the cluster
		#print N
		annotation_set = set()
		for gene in annotated_genes:
			for term in gene2go[gene]:
				annotation_set.add(term)
		#print len(annotation_set)
		for term in annotation_set:
			K = go_counts[term] # number of genes annotated with the given term in all the data
			X = 0	# number of gene in the clusters that are annotated with the go term
			for gene in annotated_genes:
				if term in gene2go[gene]:
					X += 1
			pval = hypergeom.sf(X-1, M, K, N) # faster and more accurate than 1-cdf
			clts.append(i)
			gos.append(term)
			pvals.append(pval)
			if pval <= 0.05:
				cts += 1
				#print "Cluster %i, term %s: X=%i, K=%i, pval = %s"%(i,term,X,K,str(pval))
				enrichment[i].append([term,pval])
			#print "%s %s"%(term,prb)

		#print(len(enrichment))    
		#Mean normalized entropy:
		d = float(len(annotation_set))	# number of different annotations in cluster c
		nec = 0.
		if d>1.:
			Cl = float(len(cluster))	# number of gene in cluster c
			sum_pi = 0.
			#counting map
			counting_map = {}
			for gene in cluster:
				if gene in gene2go:
					for go in gene2go[gene]:
						if go not in counting_map:
							counting_map[go] = 0.
						counting_map[go] += 1.
			for go in counting_map.keys():
				pi = counting_map[go]/Cl	# percentage of genes in c annotated with the considered go term
				sum_pi += pi*math.log(pi)
			nec = (-1./(math.log(d)))*sum_pi
		NE_list.append( nec )
	
	#applying BH correction (Benjamini-Hochner correction)
	BHpvals = p_adjust_bh(pvals)

	BHcts = 0
	BH_enrichment = [[] for j in range(len(clusters))]
	enriched_genes = []
	for i in range(len(clts)):
		#print pvals[i], BHpvals[i]
		if BHpvals[i]<=0.05:
			BHcts += 1
			BH_enrichment[clts[i]].append([gos[i],BHpvals[i]])
	for i in range(len(clusters)):
			cluster_set = set()
			enriched_gos = BH_enrichment[i]
			for gene in clusters[i]:
				for go, pval in enriched_gos:
					if gene in gene2go:
						if go in gene2go[gene]:
							cluster_set.add(gene)
			enriched_genes.append(list(cluster_set)) #genes that are enriched in each cluster
	
	#print cts, BHcts
	MNE = sum(NE_list)/float(len(NE_list))
	#print "Mean Normalized Entropy = %s"%(str(MNE))
	#print(BH_enrichment)
	enr_cluster = 0
	total_cluster = 0
	for i in range(len(clusters)):
		if len(clusters[i])>0:
			total_cluster += 1.
			if len(BH_enrichment[i])>0:
				enr_cluster += 1
	perc_enr_cluster = 100.*enr_cluster/total_cluster
	perc_genes = sum([len(enriched_genes[i]) for i in range(len(Cp))]) #number of genes enrihced in all clusters together (sum number of genes for each individual cluster)
	#print(perc_genes)
	perc_genes = 100.*perc_genes/float(len(gene2go))
	#print(perc_genes, float(len(gene2go)))    
	return [BH_enrichment, MNE, perc_genes, perc_enr_cluster]



def GO_meaning(clusters_GOs):
    GOterms_meaning = []
    for clust in clusters_GOs:
        GOterms_meaning_clst = []
        for GO_term in clust:
            try:
                GOterm_meaning = godag[GO_term].name
            except KeyError:
                GOterm_meaning = GO_term
            GOterms_meaning_clst.append(GOterm_meaning)
        GOterms_meaning.append(GOterms_meaning_clst)    
    return GOterms_meaning

def avg_stdv_enrich(ER_ccs):
	for i in range(len(ER_ccs)):
		mean = statistics.mean(ER_ccs[i])
		stdv = statistics.stdev(ER_ccs[i])
		ER_ccs[i].append(mean)
		ER_ccs[i].append(stdv)
	return ER_ccs
	        
"""
	Main code starts here
"""
from goatools.obo_parser import GODag
import pickle

indir = 'input'
cell_conds = os.listdir(indir)
cell_conds = [x for x in cell_conds if '.' not in x]
tissue = 'Skupin'
godag = GODag(obo_file="input/go-basic.obo") 

nets_1 = ['PPI','GI','MI','COEX']
nets_2 = ['PPI+GI','PPI+MI','PPI+COEX','COEX+GI','COEX+MI','GI+MI']
nets_3 = ['GI+COEX+MI','PPI+GI+MI','PPI+COEX+MI','PPI+GI+COEX']
nets_everything = ['EXP'] + nets_1 + nets_2 + nets_3 + ['ALL']

labels_key = {'Control_IPSCs':'C0', 'Control_D06':'C6', 'Control_D10':'C10', 'Control_D15':'C15', 'Control_D21':'C21', 
              'PINK1_IPSCs':'PD0', 'PINK1_D06':'PD6', 'PINK1_D15':'PD15', 'PINK1_D21':'PD21'}

nets_leg = {'EXP':'E', 'PPI':'E+P','GI':'E+G','MI':'E+M','COEX':'E+C',
            'PPI+GI':'E+P+G','PPI+MI':'E+P+M','PPI+COEX':'E+P+C','COEX+GI':'E+C+G','COEX+MI':'E+C+M','GI+MI':'E+G+M',
            'GI+COEX+MI':'E+G+C+M','PPI+GI+MI':'E+P+G+M','PPI+COEX+MI':'E+P+C+M','PPI+GI+COEX':'E+P+G+C', 'ALL':'E+P+G+C+M'}

Flag = False

ER_clusters_BP_ccs = x = [[] for i in range(len(nets_leg))]
ER_clusters_RP_ccs = x = [[] for i in range(len(nets_leg))]
ER_clusters_KP_ccs = x = [[] for i in range(len(nets_leg))]

ER_genes_BP_ccs = x = [[] for i in range(len(nets_leg))]
ER_genes_RP_ccs = x = [[] for i in range(len(nets_leg))]
ER_genes_KP_ccs = x = [[] for i in range(len(nets_leg))]

for cell_cond in cell_conds:
	if Flag:
		all_nets = []
		genes = pd.read_csv(f'{indir}/{cell_cond}/Geneslist_{cell_cond}.csv', header = 0)
		genes = genes.values
		genes = [str(item) for sublist in genes for item in sublist]
		PPIlist = genes

		fgoname_kp = f'{indir}/HSA_Kegg_Pathways.lst'
		go_counts_kp,gene2go_kp = load_GO(fgoname_kp, PPIlist)
		print("%s: %i KP annotated genes"%(tissue, len(gene2go_kp)))

		fgoname_rp = f'{indir}/HSA_Reactome_Pathways.lst'
		go_counts_rp,gene2go_rp = load_GO(fgoname_rp, PPIlist)
		print("%s: %i RP annotated genes"%(tissue, len(gene2go_rp)))

		fgoname_bp = f'{indir}/HSA_GO-BP.lst'
		go_counts_bp,gene2go_bp = load_GO(fgoname_bp, PPIlist)
		print("%s: %i BP annotated genes"%(tissue, len(gene2go_bp)))

		ER_clusters_KP = []
		ER_clusters_RP = []
		ER_clusters_BP = []

		ER_genes_KP = []
		ER_genes_RP = []
		ER_genes_BP = []   
		for nets in nets_everything:
			for root, dirs, files in os.walk(f'{indir}/{cell_cond}'):
				for file in files:
					if file == 'kMeans_G1.pkl' and 'indices' not in file:
						lspt = root.split('/')
						lspt1 = lspt[1].split('\\')
						lspt1 = lspt1 + [lspt[0]] 
						clust_meth = file.split('_')[0]
						if len(lspt1) > 2 and lspt1[1] == nets:
							print(root)
							nets = root.split('\\')[1]
							all_nets.append(nets_leg[nets])
	                        
							with open(f'{root}/{file}', 'rb') as handle:
							    G1clusters = pickle.load(handle)
							with open(f'{root}/{clust_meth}_PDEnrichClusts_indices.pkl', 'rb') as handle:
							    PDEnrClusts_indices = pickle.load(handle)
							Cp = [[str(y) for y in x] for x in G1clusters]           

							bh_kp, mne_kp, eg_kp, perc_cluster_kp = go_enrichment(Cp, gene2go_kp, go_counts_kp)
							bh_rp, mne_rp, eg_rp, perc_cluster_rp = go_enrichment(Cp, gene2go_rp, go_counts_rp)
							bh_bp, mne_bp, eg_bp, perc_cluster_bp = go_enrichment(Cp, gene2go_bp, go_counts_bp)
							#print(f'Function Enirchment : {len(bh_bp)}')

							KP_terms = [[i[0] for i in nested] for nested in bh_kp]
							RP_terms = [[i[0] for i in nested] for nested in bh_rp]
							BP_terms = [[i[0] for i in nested] for nested in bh_bp]

							BP_terms_meaning = GO_meaning(BP_terms)                 

							print(perc_cluster_kp, eg_kp)
							print(perc_cluster_rp, eg_rp)
							print(perc_cluster_bp, eg_bp)

							ER_clusters_KP.append(perc_cluster_kp)
							ER_clusters_RP.append(perc_cluster_rp)
							ER_clusters_BP.append(perc_cluster_bp)

							ER_genes_KP.append(eg_kp)
							ER_genes_RP.append(eg_rp)
							ER_genes_BP.append(eg_bp)


	else:
		all_nets = []
		ER_clusters_KP, ER_clusters_RP, ER_clusters_BP = [],[],[]
		ER_genes_KP, ER_genes_RP, ER_genes_BP = [],[],[]
		with open(f'output/{cell_cond}/{cell_cond}_ECs.txt') as ECf:
            # print(ECf.readlines())
			for i in range(len(nets_everything)):
				all_nets.append(nets_leg[ECf.readline()[:-1]])
				ER_clusters_KP.append(float(ECf.readline()[:-1]))
				ER_clusters_RP.append(float(ECf.readline()[:-1]))
				ER_clusters_BP.append(float(ECf.readline()[:-1]))
		with open(f'output/{cell_cond}/{cell_cond}_EGs.txt') as EGf:
			for i in range(len(nets_everything)):
				EGf.readline()
				ER_genes_KP.append(float(EGf.readline()[:-1]))
				ER_genes_RP.append(float(EGf.readline()[:-1]))
				ER_genes_BP.append(float(EGf.readline()[:-1]))
		all_nets = all_nets[1:] + [all_nets[0]]
		ER_clusters_KP = ER_clusters_KP[1:] + [ER_clusters_KP[0]]
		ER_clusters_RP = ER_clusters_RP[1:] + [ER_clusters_RP[0]]
		ER_clusters_BP = ER_clusters_BP[1:] + [ER_clusters_BP[0]]
		ER_genes_KP = ER_genes_KP[1:] + [ER_genes_KP[0]]
		ER_genes_RP = ER_genes_RP[1:] + [ER_genes_RP[0]]
		ER_genes_BP = ER_genes_BP[1:] + [ER_genes_BP[0]]
        
		for i in range(len(ER_clusters_KP)):
		    ER_clusters_BP_ccs[i].append(ER_clusters_BP[i])
		    ER_clusters_RP_ccs[i].append(ER_clusters_RP[i])
		    ER_clusters_KP_ccs[i].append(ER_clusters_KP[i])
		    ER_genes_BP_ccs[i].append(ER_genes_BP[i])
		    ER_genes_RP_ccs[i].append(ER_genes_RP[i])
		    ER_genes_KP_ccs[i].append(ER_genes_KP[i])
            

            
	    
	###########
    #
    #  plotting enrichments, individual cell conds
    #
    ##################
	    
	N = len(all_nets)
	ind = np.arange(N)  # the x locations for the groups
	width = 0.2       # the width of the bars
	save_dir = f'output/{cell_cond}'

	if not os.path.exists(save_dir):
		os.makedirs(save_dir) 

	if Flag:
		clusters_output_file1 = open(f'{save_dir}/{cell_cond}_ECs.txt','w')
		genes_output_file1 = open(f'{save_dir}/{cell_cond}_EGs.txt','w')

	fig = plt.figure(figsize=(16, 9))
	ax = fig.add_subplot(111)
	ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
    
	for i in range(N):
		rects1 = ax.bar(ind[i-1], ER_clusters_KP[i-1], width, color='r')
		rects2 = ax.bar(ind[i-1]+width, ER_clusters_RP[i-1], width, color='g')
		rects3 = ax.bar(ind[i-1]+width*2, ER_clusters_BP[i-1], width, color='b')

		ax.set_ylabel('Clusters with enriched annotations (%)', fontsize = 20, fontweight = 'bold')
		ax.set_xticks(ind+width)
		if labels_key[cell_cond] == 'C21' or labels_key[cell_cond] == 'PD21':
			ax.set_xticklabels(all_nets,  fontsize = 18, rotation=90) #nets_everythinga should be all_nets
		else:
			ax.set_xticklabels(len(all_nets)*[''])
		ax.tick_params(axis='y', which='major', labelsize=16)
		ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 18)
		plt.ylim(0,100)
		plt.title(f'{labels_key[cell_cond]}',  fontsize = 30, pad = 24)

		print(all_nets[i-1])
		print(ER_clusters_KP[i-1])
		print(ER_clusters_RP[i-1])
		print(ER_clusters_BP[i-1])
		print('\n')
		if Flag:
			clusters_output_file1.write(f'{all_nets[i-1]}\n{ER_clusters_KP[i-1]}\n{ER_clusters_RP[i-1]}\n{ER_clusters_BP[i-1]}\n')
    
	plt.tight_layout()    
	plt.savefig(f'{save_dir}/{cell_cond}_ECs.png', dpi = 600)	
	plt.show()	
	plt.close()

	fig = plt.figure(figsize=(16, 9))
	ax = fig.add_subplot(111)
	ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	for i in range(N):
		rects1 = ax.bar(ind[i-1], ER_genes_KP[i-1], width, color='r')
		rects2 = ax.bar(ind[i-1]+width, ER_genes_RP[i-1], width, color='g')
		rects3 = ax.bar(ind[i-1]+width*2, ER_genes_BP[i-1], width, color='b')

		ax.set_ylabel('Genes with enriched annotations (%)', fontsize = 20, fontweight = 'bold')
		ax.set_xticks(ind+width)
		if labels_key[cell_cond] == 'C21' or labels_key[cell_cond] == 'PD21':
			ax.set_xticklabels(all_nets,  fontsize = 18, rotation=90) #nets_everythinga should be all_nets
		else:
			ax.set_xticklabels(len(all_nets)*[''])
		ax.tick_params(axis='y', which='major', labelsize=16)
		ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 18, loc="upper left") 
		plt.ylim(0,60)
		plt.title(f'{labels_key[cell_cond]}',  fontsize = 30, pad = 20)

		print(all_nets[i-1])
		print(ER_genes_KP[i-1])
		print(ER_genes_RP[i-1])
		print(ER_genes_BP[i-1])
		print('\n')
		if Flag:
			genes_output_file1.write(f'{all_nets[i-1]}\n{ER_genes_KP[i-1]}\n{ER_genes_RP[i-1]}\n{ER_genes_BP[i-1]}\n')

	plt.tight_layout()
	plt.savefig(f'{save_dir}/{cell_cond}_EGs.png', dpi = 600)	
	plt.show()	
	plt.close()
	if Flag:
		clusters_output_file1.close()
		genes_output_file1.close()
        
  
        
  
#### merging results of all cell conds, for one plot for genes and other for clusters	
ER_clusters_BP_ccs = avg_stdv_enrich(ER_clusters_BP_ccs)
ER_clusters_RP_ccs = avg_stdv_enrich(ER_clusters_RP_ccs)
ER_clusters_KP_ccs = avg_stdv_enrich(ER_clusters_KP_ccs)
ER_genes_BP_ccs = avg_stdv_enrich(ER_genes_BP_ccs)
ER_genes_RP_ccs = avg_stdv_enrich(ER_genes_RP_ccs)
ER_genes_KP_ccs = avg_stdv_enrich(ER_genes_KP_ccs)        
   

capsize = 5 
alpha = 0.7
# clusters
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

for i in range(N):
	rects1 = ax.bar(ind[i-1], ER_clusters_KP_ccs[i-1][-2], yerr=ER_clusters_KP_ccs[i-1][-1], width=width, color='r', alpha = alpha, ecolor='black', capsize=capsize)
	rects2 = ax.bar(ind[i-1]+width, ER_clusters_RP_ccs[i-1][-2], yerr=ER_clusters_RP_ccs[i-1][-1], width=width, color='g', alpha = alpha, ecolor='black', capsize=capsize)
	rects3 = ax.bar(ind[i-1]+width*2, ER_clusters_BP_ccs[i-1][-2], yerr=ER_clusters_BP_ccs[i-1][-1], width=width, color='b', alpha = alpha, ecolor='black', capsize=capsize)

	ax.set_ylabel('Clusters with enriched annotations (%)', fontsize = 22, fontweight = 'bold')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(all_nets,  fontsize = 22, rotation=90) #nets_everythinga should be all_nets
	ax.tick_params(axis='y', which='major', labelsize=20)
	ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 22, loc = 'upper left')
	# plt.ylim(0,100)
# 	plt.title('Average across all cell conditions',  fontsize = 30, pad = 24)

# 	print(all_nets[i-1])
# 	print(ER_clusters_KP[i-1])
# 	print(ER_clusters_RP[i-1])
# 	print(ER_clusters_BP[i-1])
# 	print('\n')

plt.tight_layout()    
plt.savefig('output/AVG_ECs.png', dpi = 600)	
plt.show()	
plt.close()


# 3 scenarios
netsplot = ['E','E+C','E+P+G+C+M']
N = len(netsplot)
ind = np.arange(N)  # the x locations for the groups
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

for i in range(N):
    net = netsplot[i]
    index_r = all_nets.index(net)
    rects1 = ax.bar(ind[i], ER_clusters_KP_ccs[index_r][-2], yerr=ER_clusters_KP_ccs[index_r][-1], width=width, color='r', alpha = alpha, ecolor='black', capsize=capsize)
    rects2 = ax.bar(ind[i]+width, ER_clusters_RP_ccs[index_r][-2], yerr=ER_clusters_RP_ccs[index_r][-1], width=width, color='g', alpha = alpha, ecolor='black', capsize=capsize)
    rects3 = ax.bar(ind[i]+width*2, ER_clusters_BP_ccs[index_r][-2], yerr=ER_clusters_BP_ccs[index_r][-1], width=width, color='b', alpha = alpha, ecolor='black', capsize=capsize)
ax.set_ylabel('Clusters with enriched annotations (%)', fontsize = 24, fontweight = 'bold')
ax.set_xticks(ind+width)
ax.set_xticklabels(netsplot,  fontsize = 26) #nets_everythinga should be all_nets
ax.tick_params(axis='y', which='major', labelsize=24)
ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 24, loc = 'upper left')

plt.tight_layout()    
plt.savefig('output/AVG_ECs_3scs.png', dpi = 600)	
plt.show()	
plt.close()


#genes
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

for i in range(N):
	rects1 = ax.bar(ind[i-1], ER_genes_KP_ccs[i-1][-2], yerr=ER_genes_KP_ccs[i-1][-1], width=width, color='r', alpha = alpha, ecolor='black', capsize=capsize)
	rects2 = ax.bar(ind[i-1]+width, ER_genes_RP_ccs[i-1][-2], yerr=ER_genes_RP_ccs[i-1][-1], width=width, color='g', alpha = alpha, ecolor='black', capsize=capsize)
	rects3 = ax.bar(ind[i-1]+width*2, ER_genes_BP_ccs[i-1][-2], yerr=ER_genes_BP_ccs[i-1][-1], width=width, color='b', alpha = alpha, ecolor='black', capsize=capsize)

	ax.set_ylabel('Genes with enriched annotations (%)', fontsize = 22, fontweight = 'bold')
	ax.set_xticks(ind+width)
	ax.set_xticklabels(all_nets,  fontsize = 22, rotation=90) #nets_everythinga should be all_nets
	ax.tick_params(axis='y', which='major', labelsize=20)
	ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 22, loc = 'upper left')
	plt.ylim(0,60)
# 	plt.title('Average across all cell conditions',  fontsize = 30, pad = 24)

# 	print(all_nets[i-1])
# 	print(ER_clusters_KP[i-1])
# 	print(ER_clusters_RP[i-1])
# 	print(ER_clusters_BP[i-1])
# 	print('\n')

plt.tight_layout()    
plt.savefig('output/AVG_EGs.png', dpi = 600)	
plt.show()	
plt.close()
        
# 3 scenarios
netsplot = ['E','E+C','E+P+G+C+M']
N = len(netsplot)
ind = np.arange(N)  # the x locations for the groups
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.grid(color='grey', axis='y', linestyle='-', linewidth=0.25, alpha=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

for i in range(N):
	net = netsplot[i]
	index_r = all_nets.index(net)
	rects1 = ax.bar(ind[i], ER_genes_KP_ccs[index_r][-2], yerr=ER_genes_KP_ccs[index_r][-1], width=width, color='r', alpha = alpha, ecolor='black', capsize=capsize)
	rects2 = ax.bar(ind[i]+width, ER_genes_RP_ccs[index_r][-2], yerr=ER_genes_RP_ccs[index_r][-1], width=width, color='g', alpha = alpha, ecolor='black', capsize=capsize)
	rects3 = ax.bar(ind[i]+width*2, ER_genes_BP_ccs[index_r][-2], yerr=ER_genes_BP_ccs[index_r][-1], width=width, color='b', alpha = alpha, ecolor='black', capsize=capsize)

	print(all_nets[index_r])
	print(ER_genes_KP_ccs[index_r][-2])
	print(ER_genes_RP_ccs[index_r][-2])
	print(ER_genes_BP_ccs[index_r][-2])
ax.set_ylabel('Genes with enriched annotations (%)', fontsize = 24, fontweight = 'bold')
ax.set_xticks(ind+width)
ax.set_xticklabels(netsplot,  fontsize = 26) #nets_everythinga should be all_nets
ax.tick_params(axis='y', which='major', labelsize=24)
ax.legend( (rects1[0], rects2[0], rects3[0]), ('KP', 'RP', 'GO-BP'),  fontsize = 24, loc = 'upper left')

plt.tight_layout()    
plt.savefig('output/AVG_EGs_3scs.png', dpi = 600)	
plt.show()	
plt.close()

        
        
        
	# 						save_dir = f'output/{cell_cond}/{nets}/All_genes'
	# 						if not os.path.exists(save_dir):
	# 							os.makedirs(save_dir)
	# 						with open(f'{save_dir}/{clust_meth}_KP_terms.pkl', 'wb') as handle:
	# 							pickle.dump(KP_terms, handle)
	# 						with open(f'{save_dir}/{clust_meth}_RP_terms.pkl', 'wb') as handle:
	# 							pickle.dump(RP_terms, handle)            
	# 						with open(f'{save_dir}/{clust_meth}_BP_terms.pkl', 'wb') as handle:
	# 							pickle.dump(BP_terms, handle)
	# 						with open(f'{save_dir}/{clust_meth}_BP_terms_meaning.pkl', 'wb') as handle:
	# 							pickle.dump(BP_terms_meaning, handle)

	# 						### isolate just PD enriched clusts
	# 						save_dir = f'output/{cell_cond}/{nets}/PD_genes'
	# 						if not os.path.exists(save_dir):
	# 							os.makedirs(save_dir)
	# 						PDEnrClusts = {}
	# 						for clust in range(len(KP_terms)):	
	# 							if clust in PDEnrClusts_indices:
	# 								PDEnrClusts[clust] = KP_terms[clust]
	# 						with open(f'{save_dir}/{clust_meth}_KP_terms_PDEnrCls.pkl', 'wb') as handle:
	# 							pickle.dump(PDEnrClusts, handle)

	# 						PDEnrClusts = {}
	# 						for clust in range(len(RP_terms)):	
	# 							if clust in PDEnrClusts_indices:
	# 								PDEnrClusts[clust] = RP_terms[clust]
	# 						with open(f'{save_dir}/{clust_meth}_RP_terms_PDEnrCls.pkl', 'wb') as handle:
	# 							pickle.dump(PDEnrClusts, handle)

	# 						PDEnrClusts = {}
	# 						for clust in range(len(BP_terms)):	
	# 							if clust in PDEnrClusts_indices:
	# 								PDEnrClusts[clust] = BP_terms[clust]
	# 						with open(f'{save_dir}/{clust_meth}_BP_terms_PDEnrCls.pkl', 'wb') as handle:
	# 							pickle.dump(PDEnrClusts, handle)   
	# 						    
	# 						PDEnrClusts = {}
	# 						for clust in range(len(BP_terms_meaning)):	
	# 							if clust in PDEnrClusts_indices:
	# 								PDEnrClusts[clust] = BP_terms_meaning[clust]
	# 						with open(f'{save_dir}/{clust_meth}_BP_terms_meaning_PDEnrCls.pkl', 'wb') as handle:
	# 							pickle.dump(PDEnrClusts, handle)                        