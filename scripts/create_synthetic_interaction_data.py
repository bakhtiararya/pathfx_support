# written to create model interactome data for creating a new
# network
# written 22 Feb 2019, JLW

import pickle,sys,os
import networkx as nx

# old_net_file = '/Users/jenwilson/Documents/Stanford_CERSI/IntegratorCode/PathFX_clean/rscs/drug_interact_netx.pkl'
old_net_file = '/Users/jenwilson/Documents/Stanford_CERSI/IntegratorCode/forAngela/forAngela/rscs/merged_interact_netx.pkl'
G = pickle.load(open(old_net_file,'rb'))
outf = open('../data/scored_interactome_data_merged.txt','w')

for (a,b) in G.edges():
	wt_str = str(G[a][b]['weight'])
	outdata = '\t'.join([a,b,wt_str,'\n'])
	outf.write(outdata)

outf.close()
