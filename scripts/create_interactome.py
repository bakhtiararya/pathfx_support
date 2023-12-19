# Creates the interactome and returns the interactome size
# written 22 Feb 2019 JLW

import pickle,os,argparse
import networkx as nx

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-data_file', action='store', dest='fname', help='the name of the tab-delimited interaction network file')
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
args = parser.parse_args()

G = nx.Graph()
for l in open(args.fname,'rU').readlines():
	[a,b,s] = l.strip().split('\t')
	G.add_edge(a, b, weight=float(s))

intome_size = len(G.nodes())
unn = [x for x in G.nodes()]

od = [(intome_size,'intome_size'),(G,'interactome'),(unn,'unique_nodes')]
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)
