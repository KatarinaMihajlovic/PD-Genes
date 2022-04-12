import pandas as pd
import os
import networkx as nx

def CreateNetwork(GeneralNet, genes_SC):

	fRead = open(GeneralNet, 'r')
	network = nx.Graph()
	for line in fRead.readlines():
		try:
			lspt = line.strip().split('\t')
			if (int(lspt[0]) in genes_SC) and (int(lspt[1]) in genes_SC):
				network.add_edge(lspt[0],lspt[1])

		except ValueError:
			lspt = line.strip().split(' ')
			if (int(lspt[0]) in genes_SC) and (int(lspt[1]) in genes_SC):
				network.add_edge(lspt[0],lspt[1])
	fRead.close()

	print(GeneralNet.split('/')[2].split('.')[0])
	print(f"genes: {network.number_of_nodes()} \n interactions: {network.number_of_edges()}\n density: {nx.density(network)}\n")
	return(network)

def cell_conds2singlecells(Metadata_file):
	cell_conds_2_singlecells = {}
	Metadata_file.readline()
	cell_conds = set()
	for line in Metadata_file.readlines():
		lspt = line.strip().split(',')
		key = lspt[1][1:-1]
		if key in cell_conds_2_singlecells:
			cell_conds_2_singlecells[key].append(lspt[0][1:-1])
		else:
			cell_conds_2_singlecells[key] = [lspt[0][1:-1]]
		cell_conds.add(lspt[1][1:-1])
	cell_conds = list(cell_conds)
	return(cell_conds_2_singlecells, cell_conds)

if not os.path.exists('output/MolecularEXPSpecificNetworks'):
    os.makedirs('output/MolecularEXPSpecificNetworks') 
if not os.path.exists('output/Expression_Matrix'):
    os.makedirs('output/Expression_Matrix') 
    

basefolderIN = 'input/'
basefolderOUT = 'output/'   

EXPs_Entrez = pd.read_csv(f'{basefolderIN}/ExpressionMatrix/EXPs_Skupin_Entrez_Curated.csv', index_col=0)
genes_SC = EXPs_Entrez.index.tolist()

Metadata_file = open(f'{basefolderIN}/ExpressionMatrix/Metadata.csv','r')
cell_conds_2_singlecells, cell_conds = cell_conds2singlecells(Metadata_file)

GeneralNets = [basefolderIN + 'Skupin_GeneralNetworks/PPI_General_Biogrid.edgelist', basefolderIN + 'Skupin_GeneralNetworks/MI_General_KEGG.edgelist', basefolderIN + 'Skupin_GeneralNetworks/GI_General_BIOGRID+Rohusher.edgelist', basefolderIN + 'Skupin_GeneralNetworks/COEX_General_CoexpressDB.edgelist']

All_genes = set()
All_SCs = []

for cell_cond in cell_conds:
    if cell_cond != 'Control_D10':
        print(cell_cond)

        if not os.path.exists(f'output/MolecularEXPSpecificNetworks/{cell_cond}'):
            os.makedirs(f'output/MolecularEXPSpecificNetworks/{cell_cond}') 

        savedir_nets = f'{basefolderOUT}/MolecularEXPSpecificNetworks/{cell_cond}'
        MolecularEXPSpecificNetworks_info = open(f'{savedir_nets}/MolecularEXPSpecificNetworks_Statistics_{cell_cond}.txt','w')
        
        EXP_cellcond = pd.DataFrame(index=EXPs_Entrez.index)
        for SC in EXPs_Entrez:
            if SC in cell_conds_2_singlecells[cell_cond]:
                EXP_cellcond[SC] = EXPs_Entrez[SC]      
        EXP_values = EXP_cellcond.values.tolist()

        zero_genes = []
        for i in range(len(EXP_values)):
            row = EXP_values[i]
            sum_row = sum(row)
            if sum_row == 0:
                zero_genes.append(genes_SC[i])
        expressed_genes = [x for x in genes_SC if x not in zero_genes]   
        EXP_cellcond_unique = EXP_cellcond.loc[expressed_genes, :]
        
        for i in range(len(GeneralNets)):
            SpecificNet = CreateNetwork(GeneralNets[i], expressed_genes)
            net_name = GeneralNets[i].split('/')[2].split('_')[0] + f'_Skupin_{cell_cond}.edgelist'
            nx.write_edgelist(SpecificNet, f'{savedir_nets}/{net_name}', data=False)
            net = GeneralNets[i].split('/')[2].split('.')[0]
            MolecularEXPSpecificNetworks_info.write(f'{net} \n genes: {SpecificNet.number_of_nodes()} \n interactions: {SpecificNet.number_of_edges()}\n density: {nx.density(SpecificNet)}\n\n')
        MolecularEXPSpecificNetworks_info.close()   
        
        PPI = nx.read_edgelist(f'{savedir_nets}/PPI_Skupin_{cell_cond}.edgelist')
        genes_PPI = list(PPI.nodes())
        print(len(genes_PPI))

        RemovedGenes = []
        for gene in expressed_genes:
            if str(gene) not in genes_PPI:
                EXP_cellcond_unique = EXP_cellcond_unique.drop(gene)
                RemovedGenes.append(gene)
        RemovedGenes = RemovedGenes + zero_genes #no match with PPI and no expression value for the cell cond
        RemovedGenes_file = open(f'{basefolderOUT}/Expression_Matrix/{cell_cond}_RemovedGenes.txt', 'w')
        for gene in RemovedGenes:
            RemovedGenes_file.write(str(gene) + '\n')
        print(EXP_cellcond_unique)
        EXP_cellcond_unique.to_csv(f'{basefolderOUT}/Expression_Matrix/E_{cell_cond}.csv' , header = True, index = True)
        
        EXP_cellcond_unique_genes = EXP_cellcond_unique.index.tolist()
        # if cell_cond != 'ControlD10':
        for gene in EXP_cellcond_unique_genes:
                All_genes.add(gene)
        EXP_cellcond_unique_SCs = EXP_cellcond_unique.columns.tolist()
        for SC in EXP_cellcond_unique_SCs:
            All_SCs.append(SC)
            
        genes_PPI = [int(x) for x in genes_PPI]
        print(len(list(set(genes_PPI) & set(EXP_cellcond_unique_genes))))
    
with open('output/All_genes.csv', 'w') as f:
    for gene in All_genes:
        f.write(f'{gene}\n')    
with open('output/ALL_SCs.csv', 'w') as f:
    for SC in All_SCs:
        f.write(f'{SC}\n')   