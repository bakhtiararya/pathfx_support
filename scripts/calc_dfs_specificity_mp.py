# written to calculate the pathway specificity
# for all nodes during dfs
# written 3-5-19 JLW
# re-written to accommodate MP and possibly writing to RAM, 12-18-19

import pickle,os,csv,argparse
import numpy as np
import multiprocessing as mp

def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def get_psl_fpath(n,thr_dir):
	hashID = n2hash[n]
	list_fname = "psl"+hashID
	listf_root = list_fname[0:7]
	listf_dir = os.path.join(thr_dir,'path_score_lists',listf_root)
	check_dir(listf_dir)
	listf_name = list_fname+'.pkl'
	listfpath = os.path.join(listf_dir,listf_name)
	
	return listfpath	

def store_score(n,scr,thr_dir):
	listfpath = get_psl_fpath(n,thr_dir)
	if os.path.exists(listfpath):
		scr_list = pickle.load(open(listfpath,'rb'))
		scr_list.append(scr)
		pickle.dump(scr_list,open(listfpath,'wb'))
	else:
		scr_list = [scr]
		pickle.dump(scr_list,open(listfpath,'wb'))

def get_avg_score(n):
	psl_fpath = get_psl_fpath(n,thr_dir) 
	if os.path.exists(psl_fpath):
		scr_list = pickle.load(open(psl_fpath,'rb'))
	else:
		scr_list = [1]
	avg_score = np.mean(scr_list)
	return avg_score
	
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-n2hash', action='store', dest='n2hash', help="The dictionary mapping a node ID to it's numerical hash ID")
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
parser.add_argument('-thr',action='store',dest='thr',help='the score threshold for stopping the dfs')
parser.add_argument('-n2nbhd',action='store',dest='n2nbhd',help='the dictionary mapped the hashed node ID to the file path for the neighborhood file')
parser.add_argument('-thr_dir',action='store',dest='thr_dir',help='the directory where all other results for the threshold are stored')
args = parser.parse_args()

### create vectors of all path scores
n2hash = pickle.load(open(args.n2hash,'rb'))
node_to_nbhd = pickle.load(open(args.n2nbhd,'rb'))
thr_dir = pickle.load(open(args.thr_dir,'rb'))

for (nfname,nfpath) in node_to_nbhd.items():
	nbhd = pickle.load(open(nfpath,'rb'))
	for (pth,scr) in nbhd.items():
		if '@' in pth:
			lastn = pth.split('@')[-1]
			store_score(lastn,scr,thr_dir)
		else:
			lastn = pth
			store_score(lastn,scr,thr_dir)

# Return to the neighborhood and calculate specific networks, store path scores for plotting
thr_values = []
spec_nbhd_hash = {}
for (nfname,nfpath) in node_to_nbhd.items():
	nbhd = pickle.load(open(nfpath,'rb'))
	spec_ints = []
	for (pth,scr) in nbhd.items():
		if '@' in pth:
			lastn = pth.split('@')[-1]
			
		else:
			lastn = pth
		avs = get_avg_score(lastn)
		if scr > avs:
			spec_ints.append((pth,scr))
		thr_values.append(scr-avs)
	spec_int_d = dict(spec_ints)

	spn = nfname.replace("rn","spn")
	spn_froot = spn[0:7]
	spn_dir = os.path.join(thr_dir,'specific_nbhds',spn_froot)
	check_dir(spn_dir)
	spn_fname = spn+'.pkl'
	spn_fpath = os.path.join(spn_dir,spn_fname)

	spec_nbhd_hash[spn] = spn_fpath
	pickle.dump(spec_int_d,open(spn_fpath,'wb'))
		
			
od = [(thr_values,'thr_values'),(spec_nbhd_hash,'spec_nbhd_hash')]
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+args.thr+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)
