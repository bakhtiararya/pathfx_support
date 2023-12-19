# Written to calculate the expected pvalue distributions
# for random input network sizes
# written 19 Feb 2019 JLW

import pickle,os,functools
from optparse import OptionParser
from collections import defaultdict
from scipy.stats import hypergeom
import numpy as np

# methods
dd_float = functools.partial(defaultdict, float)

def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def get_nodes_from_merged(m_dic):
	node_set = set()
	for pth in m_dic.keys():
		for n in pth.split('@'):
			node_set.add(n)
	return list(node_set)

def get_assoc(node_list):
	# pull assocationes
	assoc_genes = defaultdict(list)
	for n in node_list:
		if n in genes_to_cuis:
			for cui in genes_to_cuis[n]:
				assoc_genes[cui].append(n)

	return assoc_genes

def calc_hyp(node_list):
	n = len(node_list)
	assoc_genes = get_assoc(node_list)

	assoc_analy = []
	for (a,gene_list) in assoc_genes.items():
		K = len(cui_to_genes[a])
		k = len(gene_list)
		prb = 1 - hypergeom.cdf(k,intome_size,K,n)	
		assoc_analy.append([a,k,K,prb])
	# Q = 0.001
	sort_assoc = sorted(assoc_analy,key = lambda x:(x[3],x[0]))
	m = len(sort_assoc)
	mhc_assoc = {}
	for (i,[a,k,K,prb]) in enumerate(sort_assoc):
		BH = (float(i+1)/m)*Q # calculate Benjamini-Hochberg based on ranked data
		mhc_assoc[a] = BH

	# return associations after multiple hypothesis testing
	return mhc_assoc

def main():
	parser=OptionParser()

	parser.add_option('-s','--spec_dic_hash_file',dest='spdhf',help='Pickled dictionary of hashed specific network files {gene_name:full_file_path')
	parser.add_option('-t','--num_targ',dest='numt',help='Max number of targets to choose')
	parser.add_option('-r','--num_rand',dest='numr',help='The number of random networks to create for the distribution, default is 100 if no input provided')
	parser.add_option('-c','--cui_2_genes',dest='c2g',help='Pickled dict of {cui:[gene_list]}')
	parser.add_option('-g','--genes_2_cui',dest='g2c',help='Pickled dict of {gene:[cui_list]}')
	parser.add_option('-i','--intome_size',dest='ins',help='Pickled int of total nodes in the interactome')
	parser.add_option('-a','--dup_name',dest='dupname',help='Data update name')

	(options,args) = parser.parse_args()

	# global variables
	global intome_size
	global Q
	global cui_to_genes
	global genes_to_cuis
	Q = 0.001

	cui_to_genes = pickle.load(open(options.c2g,'rb'))
	genes_to_cuis = pickle.load(open(options.g2c,'rb'))
	intome_size = pickle.load(open(options.ins,'rb'))

	sp_dic_hash = pickle.load(open(options.spdhf,'rb'))
	all_nodes = [k for k in sp_dic_hash.keys()]

	num_rand = 100
	if options.numr:
		num_rand = int(options.numr)

	# save results
	rdir = '../results/'+options.dupname+'random_networks/'
	check_dir(rdir)

	input_sizes = range(int(options.numt))
	# input_sizes = [2,5,6] # for debugging
	for nt in input_sizes:
		targ_dir = os.path.join(rdir,'targ_'+str(nt))
		check_dir(targ_dir)

		targ_results = defaultdict(list)

		counter = 0
		while counter < num_rand:
		# while counter < 5: # for debugging
			# generate a random input set
			rand_targs = np.random.choice(all_nodes,nt,replace=False)

			# get individual networks, merge networks
			merged_net = {}
			for nname in rand_targs:
				spnf = sp_dic_hash[nname] # specific network file
				merged_net.update(pickle.load(open(spnf,'rb'))) # this will add the dictionary but will replace previous path_scores 
				
			# save merged network for reproducibility
			outfname = 'merged_random_'+str(nt)+'targets_repeat'+str(counter)+'.pkl'
			pickle.dump(merged_net,open(os.path.join(targ_dir,outfname),'wb'))
			
			node_list = get_nodes_from_merged(merged_net)

			# store multiple-hypothesis-corrected  p-values for phenotypes 
			mhc_assoc = calc_hyp(node_list)

			for (cui,mhc_pv) in mhc_assoc.items():
				targ_results[cui].append(mhc_pv)

			# advance counter
			counter += 1	

		# store distribution for this number of targets
		outfname = os.path.join(targ_dir,'summarize_pval_dist_'+str(nt)+'targets.pkl')
		pickle.dump(targ_results,open(outfname,'wb'))

	# After running all randomizations, compile median values and save
	exp_med_pv = defaultdict(dd_float)
	for nt in input_sizes:
		targ_dir = os.path.join(rdir,'targ_'+str(nt))
		pv_dis_file = os.path.join(targ_dir,'summarize_pval_dist_'+str(nt)+'targets.pkl')
		pv_dis = pickle.load(open(pv_dis_file,'rb'))
		for (cui,ph_list) in pv_dis.items():
			exp_pval = np.median(ph_list)
			exp_med_pv[nt][cui] = exp_pval
	outfname = os.path.join(rdir,options.dupname+'_expected_pvalue_summary.pkl')
	pickle.dump(exp_med_pv,open(outfname,'wb'))		

	print('exp_med_pv=='+outfname)


if __name__ == "__main__":
	main()
