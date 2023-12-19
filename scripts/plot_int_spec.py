# written to look at the proportion of specific
# interaction paths at increasing threshold values
# written 3-6-19 JLW

import matplotlib, pickle, os, csv, argparse, math
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal

parser = argparse.ArgumentParser(description='Handle data for plotting.')
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')
parser.add_argument('-thr_values',action='store', dest='tval',help='a dictionary where the index is the threshold and the value is a list of [pth_score-avg_path_scr]')
parser.add_argument('-thr_hash',action='store', dest='thash',help='a dictionary where the key is the threshold and the value is a dictionary mapping node:spec_nbhd_fpath')
args = parser.parse_args()

thr_values = pickle.load(open(args.tval,'rb'))
num_thr = len(thr_values)
thrs = sorted([k for k in thr_values.keys()])
thr_index = dict([(t,i) for (i,t) in enumerate(thrs)])
thr_hash = pickle.load(open(args.thash,'rb'))

spec_vs_total = []  # for plotting specificity against total number of paths
find_thr = {} # save values with threshold
label_thr = {}

fig,ax_arr = plt.subplots(num_thr,sharex=True,figsize=(5,10))
for (thr_v,vlistf) in thr_values.items():
	vlist = pickle.load(open(vlistf,'rb'))
	#print(thr_v)
	#print('max: '+str(max(vlist)))
	#print('min: '+str(min(vlist)))
	thr_ind = thr_index[thr_v]
	num_paths = len(vlist)
	post_paths = [x for x in vlist if x > 0]
	specificity = float(len(post_paths))/num_paths

	spec_vs_total.append((num_paths,specificity))
	find_thr[thr_v] = num_paths*specificity 
	label_thr[thr_v] = (num_paths,specificity)

	ax = ax_arr[thr_ind]
	leg_label = 'thr: '+str(thr_v)+'\n'+'%.2E' % Decimal(num_paths)+' paths\n'+'{0:.2f}'.format(specificity)+' spec'
	weights = np.ones_like(vlist)/float(len(vlist))
	ax.hist(vlist, weights=weights,histtype='step',fill=True,alpha=0.3,label=leg_label,color='k',bins=25)
	# ax.set_ylim([0.1,0.5])
	ax.set_xlim([-0.05,0.25])
	ax.axvline(x=0,color='k', linestyle='dashed')
	ax.legend(prop={'size': 10},loc='upper right')
	plt.setp(ax, yticks=[0.1,0.5], yticklabels = ['0.1','0.5'],)
	plt.setp(ax, xticks=[-0.05,0,0.05,0.15,0.25], xticklabels = ['-0.05','0','0.05','0.15','0.25'],)

	if thr_ind==(num_thr-1):
		ax.set_xlabel('path_score - avg_path_score')
		plt.setp(ax,xticks=[-.2, 0.0, .2], xticklabels=['-0.2', '0.0', '0.2'],)	
	else:
		ax.tick_params(labelbottom=False)

	if thr_ind == math.floor(num_thr/2.0):
		ax.set_ylabel('proportion of paths')

fig_name = os.path.join('../rscs/',args.aname+'_specificity_scross_thr.png')
plt.savefig(fig_name,format='png')

# determine optimal threshold
max_thr = max(find_thr.items(),key = lambda x: x[1])[0]

fig,ax = plt.subplots()
(x,y) = zip(*spec_vs_total)
ax.scatter(x,y,s = 150,alpha=0.5,color='indigo')
for (thrn,(npx,spy)) in label_thr.items():
	ax.annotate(thrn,(npx,spy))
ax.set_ylim([0.1,0.5])
ax.set_xscale('log')
ax.set_ylabel('Proportion of specific paths')
ax.set_xlabel('Total number of paths')
plt.setp(ax, yticks=[0.1,0.5], yticklabels = ['0.1','0.5'],xticks=[1000,10000,100000,1000000],xticklabels=['10e4','10e5','10e6','10e7'])
ax.set_title('Specificity vs. number of total paths\noptimal threshold: '+str(max_thr)) 
fig_name = os.path.join('../rscs/',args.aname+'_spec_vs_total_paths.png')
plt.savefig(fig_name,format='png')

# get the optimal hash dictionary based on threshold
opt_sp_hash = thr_hash[max_thr] # the hash dictionary from the optimal threshold, directs to specific neighborhoods
print(opt_sp_hash)
