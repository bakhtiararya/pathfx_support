# Outline for commands for data update
# 2-15-19 JLW

import os, pickle, datetime

def get_data_files(output):
	data_dic = {}
	for l in output.split('\n'):
		if '==' in l:
			[dfname,dfpath] = l.split('==')
			data_dic[dfname] = dfpath
	return data_dic

#########################################################################
# Enter the name of the update, could be the date or specific data sources
dup_name = 'pfxapp_mike_test'

#########################################################################
### 0.1 Merge interactions with scores, pickle network file
int_fname = '../data/scored_interactome_data.txt'
cmd = "python create_interactome.py -data_file %s -a_name %s"%(int_fname,dup_name)
print(cmd)
output = os.popen(cmd).read()
# Return the intome size
intome_out = get_data_files(output)
intome_size = intome_out['intome_size']
interactome = intome_out['interactome']
print('Created interactome: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

### 0.2 Run depth-first search for all interactome nodes
thr_hashes = {}
thr_values = {}
# for thr in [0.75,0.77,0.79]:
# for thr in [0.99,0.95,0.9,0.89,0.88,0.87,0.86]:
for thr in [0.99,0.95,0.9,0.89,0.88]:
# for thr in [0.89,0.88]:
# for thr in [0.87, 0.86, 0.85]:
	# run DFS at each threshold
	cmd = "python run_depth_first_search.py -intome %s -a_name %s -thr %s"%(interactome,dup_name,thr)
	print(cmd)
	output = os.popen(cmd).read()
	dfs_dic = get_data_files(output)
	n2hash = dfs_dic['node_to_hashID']
	thr_dir = dfs_dic['thr_dir']
	node_to_nbhd = dfs_dic['node_to_nbhd']
	print('finished network at: '+str(thr)+'\t'+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

	# calc average path scores and specific networks
	cmd = "python calc_dfs_specificity.py -n2hash %s -a_name %s -thr %s -n2nbhd %s -thr_dir %s"%(n2hash,dup_name,str(thr),node_to_nbhd,thr_dir)
	print(cmd)
	output = os.popen(cmd).read()
	calc_spc_dic = get_data_files(output)
	thr_hashes[str(thr)] = calc_spc_dic['spec_nbhd_hash']
	thr_values[str(thr)] = calc_spc_dic['thr_values']
	print('calculated specificity at: '+str(thr)+'\t'+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

thr_hf = '../rscs/'+dup_name+'_thresh_hash_dic.pkl'
if os.path.isfile(thr_hf): #update dictionary if previous thresholds were already completed
	thr_hf_o = pickle.load(open(thr_hf,'rb'))
	thr_hashes.update(thr_hf_o)
	pickle.dump(thr_hashes,open(thr_hf,'wb'))
else:
	pickle.dump(thr_hashes,open(thr_hf,'wb'))
thr_vf = '../rscs/'+dup_name+'_thresh_val_dic.pkl'
if os.path.isfile(thr_vf):
	thr_v_o = pickle.load(open(thr_vf,'rb'))
	thr_values.update(thr_v_o)
	pickle.dump(thr_values,open(thr_vf,'wb'))
else:
	pickle.dump(thr_values,open(thr_vf,'wb'))
# plot all values, return the optimal threshold
cmd = "python plot_int_spec.py -a_name %s -thr_values %s -thr_hash %s"%(dup_name,thr_vf,thr_hf) 
print(cmd)
output = os.popen(cmd).read()
sp_dic_hash = output.replace('\n','') 
print('finished plotting: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

# OLD - REMOVE SUBLOCK - ### 0.3 Calculate Specific Networks ### NEW: get plotting script to return the optimal threshold and then get the hashed dictionary to pass for phenotype calculations
#cmd = "python calculate_specific_nets.py"
##print(cmd)
## output = os.popen(cmd).read()
#### Ideally return the hashed directory for all specific networks
#sp_dic_hash = '../rscs/hash_specific_nets.pkl' # summary specific dictionaries, hashed

#########################################################################
### 1.1 Pull new data from source and create formatted text files
# possibly involve Nicholas?

### 1.2 Merge gene-phenotype associations ###
data_source_dir = '../data/'
data_source_list = 'ClinVar.txt,DisGeNet.txt,HumPhenOnt.txt,PheGenI.txt'
cmd = "python merge_data_sources.py -data_list %s -data_dir %s -a_name %s"%(data_source_list,data_source_dir,dup_name)
print(cmd)
output = os.popen(cmd).read()
merge_out_dic = get_data_files(output)

### 1.2 Derive expected p-value distributions ###
cui_to_genes_file = merge_out_dic['merged_unique_cuis2genes']
gene_to_cui_file = merge_out_dic['merged_genes_to_cuis']
all_assoc_to_nodes = merge_out_dic['all_assoc_to_nodes']
sourced_phens = merge_out_dic['sourced_phens']
all_phens_to_cuis = merge_out_dic['all_phens_to_cuis']

num_targ = '40'

#for input_check in (sp_dic_hash,num_targ,cui_to_genes_file,gene_to_cui_file,intome_size,dup_name):
#	print(input_check)
#	print(type(input_check))

cmd = "python calc_exp_pvale_dist.py --spec_dic_hash_dile %s --num_targ %s --cui_2_genes %s --genes_2_cui %s --intome_size %s --dup_name %s"%(sp_dic_hash,num_targ,cui_to_genes_file,gene_to_cui_file,intome_size,dup_name)
print(cmd)
output = os.popen(cmd).read()
exp_pval_out = get_data_files(output) # {'exp_med_pv':exp_med_pv_file_path}
exp_med_pv = exp_pval_out['exp_med_pv']

#########################################################################
# PathFXapp needs the hashed dictionary for specific networks 
cmd = "python push_new_data_to_pathfx.py --exp_pval %s --spec_net_has %s"%s(exp_med_pv,sp_dic_hash)

# and get_network_associations needsi the expected pvalue distributions 
