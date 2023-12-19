# written to run all single targets
# considered for novel target discovering for two Sanofi 
# phenotypes. Written 6-2-20 JLW

import pickle,os,csv
from collections import defaultdict
from functools import partial
import multiprocessing as mp
import networkx as nx
num_cores = mp.cpu_count()
dd_dic = partial(defaultdict,dict)

def check_dir(dir_name):
        if not os.path.exists(dir_name):
                os.makedirs(dir_name)

# bin targets for multiprocessing
def get_binned_targs(tlist,nc,rdir):
	fdic = {}
	for j in range(nc):
		binv = [x for (i,x) in enumerate(tlist) if i%nc==j]
		binf = os.path.join(rdir,'all_targets_'+str(j)+'.pkl')
		pickle.dump(binv,open(binf,'wb'))
		fdic[j]=binf
	return fdic

aname = 'Leukemia_SingleTargets'
rdir = os.path.join('../results/',aname)
check_dir(rdir)
cui_set_file = '../results/leukemia/cuis_of_interest.pkl'
cui_list = list(pickle.load(open(cui_set_file,'rb')))

G = pickle.load(open('../rscs/pfx041520_interactome.pkl','rb'))
cui_to_genes_file = '../rscs/Pfx050120_merged_unique_cuis2genes.pkl'
cui2genes = pickle.load(open(cui_to_genes_file,'rb'))

drug_to_targets = '../rscs/pfxDB050620_dint.pkl'
dint = pickle.load(open(drug_to_targets,'rb'))
all_targets = list(set([dt for dtset in dint.values() for dt in dtset]))
all_targets.remove('pgk/tpi')

# generate short list of targets of interest
# all_disease_genes = set([cg for cui_term in cui_list for cg in cui2genes[cui_term]])
disease_gene_file = '../results/leukemia/genes_of_interest.pkl'
all_disease_genes = pickle.load(open(disease_gene_file,'rb'))

dis_in_net = [gn for gn in all_disease_genes if gn in G]
disease_neighbors =set([x for gname in dis_in_net for x in G.neighbors(gname) if x in all_targets])
short_list = all_disease_genes.union(disease_neighbors)
pickle.dump(short_list,open(os.path.join(rdir,'all_single_targets.pkl'),'wb'))
short_list = [[x] for x in short_list] # need to pass a list of lists to the calc_exp_pvale_allTargs_v2 script
targs_by_bin = get_binned_targs(short_list,num_cores,rdir)

# other files needed for processing
sp_dic_hash = '../rscs/pfx041520_0.82_spec_nbhd_hash.pkl'
gene_to_cui_file = '../rscs/Pfx050120_merged_genes_to_cuis.pkl'
intome_size = '../rscs/pfx041520_intome_size.pkl'
n2hash = '../rscs/pfx041520_0.8_node_to_hashID.pkl'

# add dmes to cuis of interest
rsdir = '../rscs/'
dme_c_dic = pickle.load(open(os.path.join(rsdir,'dme_2_cuis.pkl'),'rb'))
dme_cuis = list(set([v for vlist in dme_c_dic.values() for v in vlist]))
cui_to_keep = cui_list+dme_cuis
clf = os.path.join(rdir,'cuis_to_keep.pkl')

def copy_files(f1,f2):
	cmd = "cp %s %s"%(f1,f2)
	os.system(cmd)

def make_copies(fname,nc,newdir):
	(old_path,froot) = os.path.split(fname)
	new_path = os.path.join(old_path,newdir)
	check_dir(new_path)
	(fnm,fext) = os.path.splitext(froot)
	fdic = {}
	for i in range(nc):
		fnewroot = ''.join([fnm,'copy',str(i),fext])
		fnew = os.path.join(new_path,fnewroot)
		fdic[i] = fnew
		copy_files(fname,fnew)
	return fdic

# do something - I guess do the same function, and then process them for only top 5%, otherwise assign 1.
mp_dic = defaultdict(dd_dic)
for fn in [sp_dic_hash,cui_to_genes_file,gene_to_cui_file,intome_size,n2hash,clf]:
        fdic = make_copies(fn,num_cores,aname + 'Temp')
        for (i,newf) in fdic.items():
                mp_dic[i][fn]=newf

# function to use all the copies
def fxn(sf,cf,gf,intf,aname,lf,nhf,core_num,clf):
        cmd = 'python calc_exp_pvale_allTargs_v3.py -s %s -c %s -g %s -i %s -a %s -l %s -b %s -n %s -p %s'%(sf,cf,gf,intf,aname,lf,nhf,core_num,clf)
        print(cmd)
        os.system(cmd)

for cn in range(num_cores):
        [sf,cf,gf,intf,nhf,cuif] =[mp_dic[cn][sp_dic_hash],mp_dic[cn][cui_to_genes_file],mp_dic[cn][gene_to_cui_file],mp_dic[cn][intome_size],mp_dic[cn][n2hash],mp_dic[cn][clf]]
        lf = targs_by_bin[cn]
        p = mp.Process(target=fxn, args=(sf,cf,gf,intf,aname,lf,nhf,cn,cuif))
        p.start()
        p.join()

