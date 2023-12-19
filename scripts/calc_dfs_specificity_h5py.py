# written to calculate the pathway specificity
# for all nodes during dfs
# written 3-5-19 JLW
# written to test h5py as a data structure for saving path-scores

import pickle,os,csv,argparse, collections, datetime
import numpy as np
import multiprocessing as mp
from tables import *

def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def get_netw(gname,n2hash,node_to_nbhd):
	netw_r = 'no network'
	if gname in n2hash:
		gid = n2hash[gname]
		net_id = 'rn'+gid
		if net_id in node_to_nbhd:
			net_path = node_to_nbhd[net_id]
			netw = pickle.load(open(net_path,'rb'))
			netw_r =  netw
	return netw_r

def get_spec_net(pth_dic,asd,n2h):
	asd = dict([(k.decode(),v) for (k,v) in asd.items()])# convert to regular strings to be compatible with Python3
	temp_scores = [(pth,pth.split('@')[-1],scr) if '@' in pth else (pth,pth,scr) for (pth,scr) in pth_dic.items()]
	convert_scores = [(pth,n2h[termPth],scr) for (pth,termPth,scr) in temp_scores] # just convert terminal to hash id number
	keep_scores = [(pth,scr) for (pth,termPth,scr) in convert_scores if scr > asd[termPth]] # now asd should recognize the hashID of the termPath
	return dict(keep_scores)

def get_scores(pth_dic):
	temp_scores = [(pth.split('@')[-1],scr) if '@' in pth else (pth,scr) for (pth,scr) in pth_dic.items()]
	return temp_scores

def get_spdic_fpath(gname,n2hash,savedir2):
	gid = n2hash[gname]
	net_id = 'spn'+gid
	spn_froot = net_id[0:7]

	spn_dir = os.path.join(savedir2,spn_froot)
	check_dir(spn_dir)
	spn_fname = net_id+'.pkl'
	spn_fpath = os.path.join(spn_dir,spn_fname)
	return spn_fpath

def make_copies(fname,bnum):
	copy_dic = {}
	for k in range(bnum):
		oldfname = fname
		newfname = oldfname.replace('.pkl','_'+str(k)+'.pkl')
		cmd = "cp %s %s"%(oldfname,newfname)
		os.system(cmd)
		copy_dic[k] = newfname
	return copy_dic

def remove_copies(flist):
	for f in flist:
		cmd = "rm %s"%(f)
		os.system(cmd)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-n2hash', action='store', dest='n2hash', help="The dictionary mapping a node ID to it's numerical hash ID")
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
parser.add_argument('-thr',action='store',dest='thr',help='the score threshold for stopping the dfs')
parser.add_argument('-n2nbhd',action='store',dest='n2nbhd',help='the dictionary mapped the hashed node ID to the file path for the neighborhood file')
parser.add_argument('-thr_dir',action='store',dest='thr_dir',help='the directory where all other results for the threshold are stored')
parser.add_argument('-gbbin',action='store',dest='gbbin',help='genes by bin assignment')
args = parser.parse_args()
print('started analysis: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

global n2hash
global node_to_nbhd
global thr_dir
gbbin = pickle.load(open(args.gbbin,'rb')) # genes binned by number of processes
n2hash = pickle.load(open(args.n2hash,'rb'))
node_to_nbhd = pickle.load(open(args.n2nbhd,'rb'))
thr_dir = pickle.load(open(args.thr_dir,'rb'))
class TermPath(IsDescription):
	name = StringCol(10)
	origGene = StringCol(10)
	termScore = Float64Col()

def fnxn(ii,glist,h5fname,n2hf,n2nf):
	n2n = pickle.load(open(n2nf,'rb'))
	n2h = pickle.load(open(n2hf,'rb'))
	h5f = open_file(h5fname,mode="a")
	rsg = h5f.create_group(h5f.root, "rawScores", "Raw Scores") #raw scores group
	rstable = h5f.create_table(rsg, 'saveData', TermPath, "Save Raw Scores")
	for g in glist:
		g_id = n2h[g]
		netw = get_netw(g,n2h,n2n)
		if netw != 'no network':
			termpath = rstable.row
			for (tg,pth_scr) in get_scores(netw):
				tg_id = n2h[tg]
				termpath['name'] = tg_id
				termpath['origGene'] = g_id
				termpath['termScore'] = pth_scr
				termpath.append()

			rstable.flush()

	# print(h5f)
	h5f.close()

# create multiple h5 files to dump all the path scores
print('distributing summary')
all_h5 = []
n2h_dic = make_copies(args.n2hash,len(gbbin)) # node to Hash number
n2n_dic = make_copies(args.n2nbhd,len(gbbin)) # node to neighborhood file path (using the 'nXXXX' notation

savedir = os.path.join(thr_dir,'h5files')
check_dir(savedir)   
for (bnum,gnodes) in gbbin.items():
	n2n = n2n_dic[bnum]
	n2h = n2h_dic[bnum] 
	h5fname = os.path.join(savedir,"pthwy_spec_"+str(bnum)+".h5") 
	h5file = open_file(h5fname,mode="w",title="Pathway Specificty")
	h5file.close()

	p = mp.Process(target=fnxn, args=(bnum,gnodes,h5fname,n2h,n2n))
	p.start()
	p.join()

	all_h5.append(h5fname)

# merge all Pytables into one
h5fname = os.path.join(savedir,"pthwy_spec_merged.h5")
h5f = open_file(h5fname,mode="w",title="Merged Pathway Specificity")
h5f.create_group(h5f.root, "rawScores", "Merged Raw Scores")
h5f.close()

# initialize with the first file
h5f = open_file(h5fname,mode="a")
sh5f = open_file(all_h5[0],'r')
sx = sh5f.root.rawScores.saveData
rst = sh5f.copy_node(sh5f.root.rawScores, name = 'saveData', newparent = h5f.root.rawScores, newName = 'saveData')
sh5f.close()
# extend with remaining files
for sh5fname in all_h5[1:]:
	sh5f = open_file(sh5fname,'r')
	sx = sh5f.root.rawScores.saveData
	rst.append(sx[:])
	rst.flush()
	sh5f.close()
h5f.close()
print('merged summary data')

# create summary values, save average in dictionary, and keep differences in separate table
h5f = open_file(h5fname,mode="a")
table = h5f.root.rawScores.saveData
all_genes = list(set([x['name'] for x in table]))

class AvgScore(IsDescription):
	name = StringCol(10)
	avgScore = Float64Col()
class ScoreDiff(IsDescription): 
	diffScore = Float64Col()

ssg = h5f.create_group(h5f.root,"sumScores","Average Path Score Values")
sstable = h5f.create_table(ssg,"avgScores",AvgScore,"Saved Average Scores")
difftable = h5f.create_table(ssg,"scoreDiff",ScoreDiff,"All differences between scores and averages")

as_dic = {}
for gname in all_genes: # note: gname is now the hash id for the gene
	all_scores = [x['termScore'] for x in table.where("(name==gname)")]
	ascore = np.mean(all_scores) 

	avgscore = sstable.row
	avgscore['name'] = gname
	avgscore['avgScore'] = ascore
	avgscore.append()
	sstable.flush()

	scrdiff = difftable.row
	for rs in all_scores:
		scrdiff['diffScore'] = (rs-ascore)
		scrdiff.append()
	difftable.flush()
	as_dic[gname] = ascore # save averages by hash id

# get all score differences and store the array
thr_values = [x['diffScore'] for x in difftable]
outfname = os.path.join('../rscs/',args.aname+'_'+args.thr+'_thr_values.pkl')
pickle.dump(thr_values,open(outfname,'wb'))
print('thr_values=='+outfname)
print('saved score differences for all paths')

h5f.close()
print('saved all average scores')

asd_fname =os.path.join(savedir,'avg_score_dic.pkl') 
pickle.dump(as_dic,open(asd_fname,'wb'))

# make copies of the summary dictionary
asd_by_bin = make_copies(asd_fname,len(gbbin)) # returns a dictionary with bin# as the key, and a filepath as the value

global savedir2
savedir2 = os.path.join(thr_dir,'specific_nbhds_mp')
check_dir(savedir2)

def fxn2(asd_fname,glist,bnum,n2nf,n2hf,savedir2):
	asd = pickle.load(open(asd_fname,'rb'))
	n2n = pickle.load(open(n2nf,'rb'))
	n2h = pickle.load(open(n2hf,'rb'))
	spn_hash = {}
	for g in glist:
		g_id = n2h[g]
		netw = get_netw(g,n2h,n2n)
		spec_net = get_spec_net(netw,asd,n2h) # keys for asd are now hash id numbers
		sd_fname = get_spdic_fpath(g,n2h,savedir2)
		pickle.dump(spec_net,open(sd_fname,'wb'))
		sd_name = os.path.split(sd_fname)[-1].replace('.pkl','')
		spn_hash[sd_name] = sd_fname
	pickle.dump(spn_hash,open(os.path.join(savedir2,'spn_hash_batch_'+str(bnum)+'.pkl'),'wb'))

for (bnum,gnodes) in gbbin.items():
	savedir_temp = savedir2
	asd_fname = asd_by_bin[bnum]
	n2nf = n2n_dic[bnum]
	n2hf = n2h_dic[bnum]
	p = mp.Process(target=fxn2, args=(asd_fname,gnodes,bnum,n2nf,n2hf,savedir_temp))
	p.start()
	p.join()

# merge all specific hash dictionaries
spec_nbhd_hash = {}
batch_spn = [f for f in os.listdir(savedir2) if 'spn_hash_batch_' in f]
for spf in batch_spn:
	spd = pickle.load(open(os.path.join(savedir2,spf),'rb'))
	spec_nbhd_hash.update(spd)

od = [(spec_nbhd_hash,'spec_nbhd_hash')]
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+args.thr+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)

print('cleaning')
# remove copied files
remove_copies([v for v in n2h_dic.values()])
remove_copies([v for v in n2n_dic.values()])
remove_copies([v for v in asd_by_bin.values()])
remove_copies(all_h5)
print('finished at: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
			
#### Development code from original script ###
#def get_psl_fpath(n,thr_dir):
#	hashID = n2hash[n]
#	list_fname = "psl"+hashID
#	listf_root = list_fname[0:7]
#	listf_dir = os.path.join(thr_dir,'path_score_lists',listf_root)
#	check_dir(listf_dir)
#	listf_name = list_fname+'.pkl'
#	listfpath = os.path.join(listf_dir,listf_name)
#	
#	return listfpath	
#
#def store_score(n,scr,thr_dir):
#	listfpath = get_psl_fpath(n,thr_dir)
#	if os.path.exists(listfpath):
#		scr_list = pickle.load(open(listfpath,'rb'))
#		scr_list.append(scr)
#		pickle.dump(scr_list,open(listfpath,'wb'))
#	else:
#		scr_list = [scr]
#		pickle.dump(scr_list,open(listfpath,'wb'))
#
#def get_avg_score(n):
#	psl_fpath = get_psl_fpath(n,thr_dir) 
#	if os.path.exists(psl_fpath):
#		scr_list = pickle.load(open(psl_fpath,'rb'))
#	else:
#		scr_list = [1]
#	avg_score = np.mean(scr_list)
#	return avg_score
	
## Return to the neighborhood and calculate specific networks, store path scores for plotting
#thr_values = []
#spec_nbhd_hash = {}
#for (nfname,nfpath) in node_to_nbhd.items():
#	nbhd = pickle.load(open(nfpath,'rb'))
#	spec_ints = []
#	for (pth,scr) in nbhd.items():
#		if '@' in pth:
#			lastn = pth.split('@')[-1]
#			
#		else:
#			lastn = pth
#		avs = get_avg_score(lastn)
#		if scr > avs:
#			spec_ints.append((pth,scr))
#		thr_values.append(scr-avs)
#	spec_int_d = dict(spec_ints)
#
#	spn = nfname.replace("rn","spn")
#	spn_froot = spn[0:7]
#	spn_dir = os.path.join(thr_dir,'specific_nbhds',spn_froot)
#	check_dir(spn_dir)
#	spn_fname = spn+'.pkl'
#	spn_fpath = os.path.join(spn_dir,spn_fname)
#
#	spec_nbhd_hash[spn] = spn_fpath
#	pickle.dump(spec_int_d,open(spn_fpath,'wb'))
#	
