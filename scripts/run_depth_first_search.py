# runs the depth-first search on all nodes in the network
# written 22 Feb 2019 JLW

import pickle,os,argparse
import networkx as nx
import find_network_dfs_noFT as dfs #updated to skip fast tracking, 6 JUL 2019

def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-intome', action='store', dest='intome', help='the pickled interactome')
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
parser.add_argument('-thr',action='store',dest='thr',help='the score threshold for stopping the dfs')
args = parser.parse_args()

G = pickle.load(open(args.intome,'rb'))
dfs.G = G # set up the network for dfs
gnodes = [(n,i) for (i,n) in enumerate(sorted(G.nodes()))]

# store hashing as uppercase (human only?)
node_to_hashID = dict([(n,"{:06d}".format(i)) for (n,i) in gnodes])
node_to_nbhd = {}

rdir = '../results/'+args.aname+'_dfs_all_nodes/'

thr_dir = os.path.join(rdir,'thr_'+str(args.thr))
check_dir(thr_dir)

for (n,i) in gnodes:
#for (n,i) in gnodes[0:100]:
	nfname = "rn"+"{:06d}".format(i)
	ndir_root = nfname[0:6]
	ndir = os.path.join(thr_dir,'raw_networks',ndir_root)
	check_dir(ndir)
	neighborhood = dfs.find_neighborhood(n,float(args.thr))
	nbhd_fpath = os.path.join(ndir,nfname+'_nbhd.pkl') 
	pickle.dump(neighborhood,open(nbhd_fpath,'wb'))
	node_to_nbhd[nfname] = nbhd_fpath 	

od = [(node_to_hashID,'node_to_hashID'),(thr_dir,'thr_dir'),(node_to_nbhd,'node_to_nbhd')] 
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)
