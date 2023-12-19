# written to bundle the components
# for the release of PathFX
# written 4-28-20 JLW

from zipfile import ZipFile
import pickle,os,csv
from optparse import OptionParser

def copy_files(old_dir,f,new_dir):
	op = os.path.join(old_dir,f)
	np = os.path.join(new_dir,f)
	cmd = "cp %s %s"%(op,np)
	print(cmd)
	os.system(cmd)

def check_dir(dname):
	if not os.path.exists(dname):
		os.makedirs(dname)

def find_all_files(dupnm,pdir,thr):
	# Resources Files
	rscs_dir = '../rscs/'
	nr_dir = os.path.join(pdir,'rscs')
	check_dir(nr_dir)

	# hard code some files relevant to the patch
	nonDupFiles = ['pfx041520_unique_nodes.pkl','pfx041520_0.82_spec_nbhd_hash.pkl','pfx041520_0.8_node_to_hashID.pkl','pfx041520_interactome.pkl','pfx041520_intome_size.pkl']

	
	assimg = [f for f in os.listdir(rscs_dir) if 'pfx041520' in f and 'png' in f] # assessment images
	rsfs = [f for f in os.listdir(rscs_dir) if dupnm in f and 'pkl' in f] # get resources files
	rm_files = ['_r.pkl','thresh_hash_dic.pkl','thresh_val_dic.pkl','thr_dir.pkl','thr_values.pkl','bins_by_genes.pkl','genes_by_bin.pkl']
	cl_rsfs = [f for f in rsfs for fp in rm_files if fp in f] # these are development files not needing for packaging
	
	for f in rsfs+assimg:
		if f not in cl_rsfs:
			copy_files(rscs_dir,f,nr_dir)

	# expected phenotype distributions
	expf = '../results/DUPNAMErandom_networks/DUPNAME_expected_pvalue_summary.pkl'
	fexpf = expf.replace('DUPNAME',dupnm)
	(rdir,f) = os.path.split(fexpf)
	new_rdir = rdir.replace("..",pdir)
	check_dir(new_rdir)
	copy_files(rdir,f,new_rdir)

	# get relevant scripts
	ns_dir = os.path.join(pdir,'scripts')
	check_dir(ns_dir)
	scpfs = [f for f in os.listdir('.') if dupnm in f and '.swp' not in f]
	for f in scpfs:
		copy_files('.',f,ns_dir)

	# get specific neighborhood files
	#spnbhd_dir_root = '../results/DUPNAME_dfs_all_nodes/thr_THR/specific_nbhds_mp'
	spnbhd_dir_root = '../results/pfx041520_dfs_all_nodes/thr_0.82/specific_nbhds_mp'
	spnbhd_dir = spnbhd_dir_root.replace('DUPNAME',dupnm).replace('THR',thr)
	new_sndir = spnbhd_dir.replace('..',pdir)
	check_dir(new_sndir)
	cmd = 'cp -r %s %s'%(spnbhd_dir,new_sndir)
	print(cmd)
	os.system(cmd)


def main():
	parser=OptionParser()
	parser.add_option('-d','--dupname',dest='dupnm',help='Data update name')
	parser.add_option('-t','--threshld',dest='thr',help='Optimal threshold needed for packaging')

	(options,args) = parser.parse_args()
	dupnm = options.dupnm
	thr = options.thr

	# create home for the project
	proj_nm = dupnm + '_exec_code'
	proj_dir = os.path.join('..',proj_nm)
	check_dir(proj_dir)

	# move over relevant file paths
	find_all_files(dupnm,proj_dir,thr)

	# archiv the folder
	cmd = 'zip -r %s.zip %s'%(proj_dir,proj_dir)
	print(cmd)
	os.system(cmd)

if __name__ == "__main__": 
    main()
