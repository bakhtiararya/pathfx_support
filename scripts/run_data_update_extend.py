# Outline for commands for data update
# 2-15-19 JLW

import os, pickle, datetime 
from functools import partial
from collections import defaultdict

#########################################################################
# method for grabbing output at each stage
def get_data_files(output):
	data_dic = {}
	for l in output.split('\n'):
		if '==' in l:
			[dfname,dfpath] = l.split('==')
			data_dic[dfname] = dfpath
	return data_dic

#########################################################################
# Enter the name of the update, could be the date or specific data sources
dup_name = 'pfxapp_geneTest'

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
unn = intome_out['unique_nodes']
print('\nCreated interactome: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

### 0.2 Import new DrugBank binding data - OPTION: wget and unzip file?
parsed_drugbank = '../data/proteins.tsv'
mapping_drugbank = '../data/drugbank.tsv'
uniprot_table = '../data/uniprot_to_gene_name_table.txt'
cmd = 'python create_drug_bank_data.py -a_name %s -proteins_file %s -drug_file %s -map_file %s'%(dup_name,parsed_drugbank,mapping_drugbank,uniprot_table) 
print(cmd)
output = os.popen(cmd).read()
dbank_dic = get_data_files(output)
dbid2name = dbank_dic['dbid2name']
drug_to_targets = dbank_dic['dint']
print('\nfinished new DrugBank data')

### 0.3 Run depth-first search for all interactome nodes
print('\nStarting to draw networks at input thresholds')
thr_hashes = {}
thr_values = {}
# for thr in [0.99,0.95,0.9]:
# for thr in [0.89,0.88]:
for thr in [0.87, 0.86, 0.85]:
# for thr in [0.75,0.77,0.79]:
	### run DFS at each threshold
	cmd = "python run_depth_first_search.py -intome %s -a_name %s -thr %s"%(interactome,dup_name,thr)
	print(cmd)
	output = os.popen(cmd).read()
	dfs_dic = get_data_files(output)
	n2hash = dfs_dic['node_to_hashID']
	thr_dir = dfs_dic['thr_dir']
	node_to_nbhd = dfs_dic['node_to_nbhd']
	print('finished network at: '+str(thr)+'\t'+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

	### 0.31 calc average path scores and specific networks
	cmd = "python calc_dfs_specificity.py -n2hash %s -a_name %s -thr %s -n2nbhd %s -thr_dir %s"%(n2hash,dup_name,str(thr),node_to_nbhd,thr_dir)
	print(cmd)
	output = os.popen(cmd).read()
	calc_spc_dic = get_data_files(output)
	thr_hashes[str(thr)] = calc_spc_dic['spec_nbhd_hash']
	thr_values[str(thr)] = calc_spc_dic['thr_values']
	print('calculated specificity at: '+str(thr)+'\t'+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

### Update results for a threshold, aware that other thresholds may have already been completed
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

### 0.4 plot all values, return the optimal threshold
cmd = "python plot_int_spec.py -a_name %s -thr_values %s -thr_hash %s"%(dup_name,thr_vf,thr_hf) 
print(cmd)
output = os.popen(cmd).read()
sp_dic_hash = output.replace('\n','') 
print('\nfinished plotting: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

#########################################################################
### 1.1 Pull new data from source and create formatted text files

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
cui_to_phens = merge_out_dic['cui_to_phens']

num_targ = '40'

print('\Calculating expected p-value distributions')
cmd = "python calc_exp_pvale_dist.py --spec_dic_hash_dile %s --num_targ %s --cui_2_genes %s --genes_2_cui %s --intome_size %s --dup_name %s"%(sp_dic_hash,num_targ,cui_to_genes_file,gene_to_cui_file,intome_size,dup_name)
print(cmd)
output = os.popen(cmd).read()
exp_pval_out = get_data_files(output) # {'exp_med_pv':exp_med_pv_file_path}
exp_med_pv = exp_pval_out['exp_med_pv']

#########################################################################
### 2.0 Incorporate updated variables and create updated versions of PathFX
cmd = "python push_new_data_to_pathfx.py -phen2cui %s -cui2phen %s -cui2gene %s -intome_size %s -gene2cui %s -source_phen %s -exp_pvals %s -drug_targets %s -drugname_to_DBid %s -gene_hash %s -sum_hash %s -unique_nodes %s -update_name %s"%(all_phens_to_cuis, cui_to_phens, cui_to_genes_file, intome_size, gene_to_cui_file, sourced_phens, exp_med_pv, drug_to_targets, dbid2name, n2hash, sp_dic_hash, unn, dup_name)

output = os.popen(cmd).read()

### TEST!? ###
print('\n\nfinished creating PathFX, example usage: \n\n')
cmd = 'python phenotype_enrichment_pathway_pfxapp_july19.py -d Metformin -a test_data_update'
print(cmd)
