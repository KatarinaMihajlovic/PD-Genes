# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""

import networkx as nx
import numpy as np

#Load network
def Load_Network(fname):
	net = nx.read_edgelist(fname)
	nodes = [n for n in net.nodes()]
	nb_nodes = len(nodes)
	nodes2ind = {}
	for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n
	return net, nodes, nodes2ind



#Computing Adjacency matrix, A
def Make_Adj(net, net_nodes, net_n2i, MOL_or_EXP):
	if MOL_or_EXP == 'MOL':
		nb_nodes = len(net_nodes)
		A = np.zeros((nb_nodes,nb_nodes))
		
		for e in net.edges():
			if e[0] in net_n2i and e[1] in net_n2i:
				n1=net_n2i[ e[0] ]
				n2=net_n2i[ e[1] ]
				A[n1][n2] = 1.
				A[n2][n1] = 1.
		return A

	elif MOL_or_EXP == 'EXP':
		gene_nodes = net_nodes
		SC_nodes = list(net)
		i2SC = dict(list(enumerate(list(net))))
		A = np.zeros((len(gene_nodes),len(SC_nodes)))
		
		for gene in range(len(gene_nodes)):
			for SC in range(len(SC_nodes)):
				A[gene][SC] = net[SC_nodes[SC]][int(gene_nodes[gene])]
		return A, i2SC


def Save_Matrix_Factor(M, fname, rownames, colnames):
	n,m = M.shape
	ofile = open(fname, 'w')
	#column header
	for i in range(m):
		ofile.write("\t%s"%(colnames[i]))
	ofile.write("\n")
	#rows
	for i in range(n):
		ofile.write("%s"%(rownames[i]))
		for j in range(m):
			ofile.write("\t%s"%(str(M[i][j])))
		ofile.write("\n")
	ofile.close()

def Save_Matrix_Factor_no_headers(M, fname, rownames, colnames):
	n,m = M.shape
	ofile = open(fname, 'w')
	#column header
	# for i in range(m):
	# 	ofile.write("\t%s"%(colnames[i]))
	# ofile.write("\n")
	#rows
	for i in range(n):
		# ofile.write("%s"%(rownames[i]))
		for j in range(m):
			if j > 0: ofile.write("\t")
			ofile.write("%s"%(str(M[i][j])))
		ofile.write("\n")
	ofile.close()
	
	
def Load_Matrix_Factor(fname):
	ifile = open(fname, 'r')
	#column header
	line = ifile.readline()
	colnames = line.strip().split('\t')
	rownames = []
	matrix = []
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		rownames.append(lspt[0])
		matrix.append( [float(i) for i in lspt[1:]])
	ifile.close()
	return matrix, rownames, colnames
