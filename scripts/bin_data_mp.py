# written to bin data by the number of cores
# written 12-18-19 JLW

import pickle,os,argparse
import networkx as nx
from collections import defaultdict

parser = argparse.ArgumentParser(description='Process input files.')
parser.add_argument('-data_file', action='store', dest='fname', help='the name of the pickled list of all unique network nodes')
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
parser.add_argument('-num_cores',action='store', dest='ncores',help='the number of cores available')
args = parser.parse_args()

stack = pickle.load(open(args.fname,'rb'))
n = int(args.ncores)

bassign = defaultdict(list) # bin assignments
bassign_r = {}
bin_num = 0
while stack:
	x = stack[0]
	stack.remove(x)
	bassign[bin_num].append(x)
	bassign_r[x] = bin_num
	if bin_num != (n-1):
		bin_num+=1
	else:
		bin_num=0

od = [(bassign,'genes_by_bin'),(bassign_r,'bins_by_genes'),]
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)
