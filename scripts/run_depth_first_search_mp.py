# runs the depth-first search on all nodes in the network
# added multiprocessing 
# re-written 12-18-19 JLW

import pickle,os,argparse
import networkx as nx
import find_network_dfs_noFT as dfs #updated to skip fast tracking, 6 JUL 2019
import multiprocessing as mp

def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def remove_copies(flist):
	for f in flist:
		cmd = "rm %s"%(f)
		os.system(cmd)

parser = argparse.ArgumentParser(description='Process input arguments.')
parser.add_argument('-intome', action='store', dest='intome', help='the pickled interactome')
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
parser.add_argument('-thr',action='store',dest='thr',help='the score threshold for stopping the dfs')
parser.add_argument('-gbbin',action='store',dest='gbbin',help='genes by bin assignment')
args = parser.parse_args()

gbbin = pickle.load(open(args.gbbin,'rb'))

# save the hash based on the position in the original graph
G = pickle.load(open(args.intome,'rb'))
gnodes = [(n,i) for (i,n) in enumerate(sorted(G.nodes()))]
gnode_dic = dict(gnodes)

# directories for saving data
rdir = '../results/'+args.aname+'_dfs_all_nodes/'
thr_dir = os.path.join(rdir,'thr_'+str(args.thr))
check_dir(thr_dir)

# store hashing as uppercase (human only?)
node_to_hashID = dict([(n,"{:06d}".format(i)) for (n,i) in gnodes])

# first copy the network for each core
nets_per_bin = {}
for k in gbbin.keys():
	oldfname = args.intome
	newfname = oldfname.replace('.pkl','_'+str(k)+'.pkl')
	cmd = "cp %s %s"%(oldfname,newfname)
	print(cmd)
	os.system(cmd)
	nets_per_bin[k] = newfname

def f(gfile,gnodes,bnum):
	#G = pickle.load(open(args.intome,'rb'))
	#dfs.G = G # set up the network for dfs
	#gnodes = [(n,i) for (i,n) in enumerate(sorted(G.nodes()))]

	G = pickle.load(open(gfile,'rb'))
	dfs.G = G
	node_to_nbhd = {}

	for n in gnodes:
		i = gnode_dic[n]
		nfname = "rn"+"{:06d}".format(i)
		ndir_root = nfname[0:6]
		ndir = os.path.join(thr_dir,'raw_networks',ndir_root)
		check_dir(ndir)
		neighborhood = dfs.find_neighborhood(n,float(args.thr))
		nbhd_fpath = os.path.join(ndir,nfname+'_nbhd.pkl') 
		pickle.dump(neighborhood,open(nbhd_fpath,'wb'))
		node_to_nbhd[nfname] = nbhd_fpath 	
	pickle.dump(node_to_nbhd,open(os.path.join(thr_dir,'temp_n2nbhd_'+str(bnum)+'.pkl'),'wb'))

for (bnum,gnodes) in gbbin.items(): 
	gfile = nets_per_bin[bnum]
	p = mp.Process(target=f, args=(gfile,gnodes,bnum))
	p.start()
	p.join()

# merge all node-to-neighborhood files
node_to_nbhd = {}
allf = [f for f in os.listdir(thr_dir) if 'temp_n2nbhd_' in f]
print(allf)
for df in allf:
	d = pickle.load(open(os.path.join(thr_dir,df),'rb'))
	node_to_nbhd.update(d)

print('cleaning')
remove_copies([v for v in nets_per_bin.values()])

od = [(node_to_hashID,'node_to_hashID'),(thr_dir,'thr_dir'),(node_to_nbhd,'node_to_nbhd')] 
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+str(args.thr)+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)
