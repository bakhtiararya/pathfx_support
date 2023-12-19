# Outline for commands for data update
# 2-15-19 JLW

import os,pickle,datetime

def get_data_files(output):
	data_dic = {}
	for l in output.split('\n'):
		if '==' in l:
			[dfname,dfpath] = l.split('==')
			data_dic[dfname] = dfpath
	return data_dic

#########################################################################
# Enter the name of the update, could be the date or specific data sources
dup_name = 'pfxapp_feb19'

#########################################################################
### 0.1 Merge interactions with scores, pickle network file
int_fname = '/Users/jenwilson/Documents/Stanford_CERSI/IntegratorCode/pathfx_data_update/data/scored_interactome_data_merged.txt'
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
for thr in [0.75,0.77,0.79]:
# for thr in [0.70,0.73,0.8,0.85]:
	# run DFS at each threshold
	cmd = "python run_depth_first_search.py -intome %s -a_name %s -thr %s"%(interactome,dup_name,thr)
	output = os.popen(cmd).read()
	dfs_dic = get_data_files(output)
	n2hash = dfs_dic['node_to_hashID']
	thr_dir = dfs_dic['thr_dir']
	node_to_nbhd = dfs_dic['node_to_nbhd']
	print('finished network at: '+str(thr)+'\t'+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
	
	# calc average path scores and specific networks
	cmd = "python calc_dfs_specificity.py -n2hash %s -a_name %s -thr %s -n2nbhd %s -thr_dir %s"%(n2hash,dup_name,str(thr),node_to_nbhd,thr_dir)
	output = os.popen(cmd).read()
	calc_spc_dic = get_data_files(output)
	thr_hashes[str(thr)] = calc_spc_dic['spec_nbhd_hash']
	thr_values[str(thr)] = calc_spc_dic['thr_values']
	print('calculated specificity at: '+str(thr)+'\t'+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

thr_hf = '../rscs/'+dup_name+'_thresh_hash_dic.pkl'
pickle.dump(thr_hashes,open(thr_hf,'wb'))
thr_vf = '../rscs/'+dup_name+'_thresh_val_dic.pkl'
pickle.dump(thr_values,open(thr_vf,'wb'))
# plot all values, return the optimal threshold
cmd = "python plot_int_spec.py -a_name %s -thr_values %s -thr_hash %s"%(dup_name,thr_vf,thr_hf) 
output = os.popen(cmd).read()

print('finished plotting: '+'{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
