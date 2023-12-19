# written to calculate and store 'specific'
# networks for each node in the interactome
# written 19 Feb 2019 JLW

import pickle,os,csv
import numpy as np

def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def clean_node_name(sn):
	if '+' in sn:
		sn = sn.replace('+','-')
	if ':' in sn:
		sn = sn.replace(':','-')
	if '/' in sn:
		sn = sn.replace('/','-')
	if ' ' in sn:
		sn = sn.replace(' ','')
	if '*' in sn:
		sn = sn.replace('*','star')
	 
	return sn

rscs_dir = '../rscs/'
res_dir = '../results/'
new_hash_net = pickle.load(open(os.path.join(rscs_dir,'gene_to_hash_map.pkl'),'rb'))
new_hash_sum_files = pickle.load(open(os.path.join(rscs_dir,'gene_to_sum_hash_map.pkl'),'rb'))

sp_net_dir = '../results/node_specific_results/'
check_dir(sp_net_dir)

hash_specific_nets = {}
for (i,(t,tfile)) in [(i,(t,tfile)) for (i,(t,tfile)) in enumerate(sorted(new_hash_net.items()))]:

	net_index = "{:06d}".format(i)
	spd_fname = 'sp_nets_'+str(net_index)+'.pkl' #fname reflects the index
	spd_dir_name = os.path.join(res_dir,sp_net_dir,'sp_nets_'+net_index[0:4])

	check_dir(spd_dir_name) 

	# print((t,net_index))
	
	pth_dic = pickle.load(open(tfile,'rb'))
	spec_ints = [] # list for storing specific interactions
	for (pth,s) in pth_dic.items():
		if '@' in pth:
			pg = pth.split('@')[-1].upper()
		else:
			continue
			# spec_dic.append((pth,1)) # skip adding the target alone
		grf = None
		if pg.upper() in new_hash_sum_files:
			grf  = new_hash_sum_files[pg.upper()]
		elif clean_node_name(pg).upper() in new_hash_sum_files:
			grf  = new_hash_sum_files[clean_node_name(pg).upper()]
		else:
			print(pg)
		pg_scores = pickle.load(open(grf,'rb'))
		avg = np.mean(pg_scores)
		if (s-avg) >0:
			spec_ints.append((pth,s))
	spec_dic = dict(spec_ints)

	outfname = os.path.join(spd_dir_name,spd_fname)
	pickle.dump(spec_dic,open(outfname,'wb'))

	hash_specific_nets[t]=outfname

pickle.dump(hash_specific_nets,open(os.path.join(rscs_dir,'hash_specific_nets.pkl'),'wb'))


