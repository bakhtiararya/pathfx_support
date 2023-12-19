# creates the pickled dictionaries using DrugBank data
# written 7-22-19 JLW

import pickle,os,argparse, csv
from collections import defaultdict

parser = argparse.ArgumentParser(description='Process some file paths for data')
parser.add_argument('-proteins_file', action='store', dest='pfile', help='tab-delimited text file, with drugbank to target information')
parser.add_argument('-drug_file',action='store', dest='dfile',help='tab-delimited file with drugbank to name mapping')
parser.add_argument('-map_file',action='store', dest='mfile',help='tab-delimited file mapping entrex gene IDs to gene symbol')
parser.add_argument('-a_name',action='store', dest='aname',help='the name of the data update')

args = parser.parse_args()

# read in mapping for uniprot
u2g = {}
dR = csv.DictReader(open(args.mfile,'r'),delimiter='\t')
for r in dR:
	u2g[r['From']] = r['To']

dbid2name = {}
dR = csv.DictReader(open(args.dfile,'r'),delimiter='\t')
for r in dR:
	[dbid,drug_name] = [r['drugbank_id'],r['name']]
	dbid2name[dbid] = drug_name

drug_to_targets = defaultdict(set)
dR = csv.DictReader(open(args.pfile,'r'),delimiter='\t')
for r in dR:
	[dbid,uniid] = [r['drugbank_id'],r['uniprot_id']]
	drug_name = dbid2name[dbid]

	# convert from uniprot to gene symbol
	if uniid in u2g:
		gene_sym = u2g[uniid]
		drug_to_targets[dbid].add(gene_sym)
		drug_to_targets[drug_name].add(gene_sym)

od = [(drug_to_targets,'dint'),(dbid2name,'dbid2name')]
for (dat_obj,obj_name) in od:
	outfname = os.path.join('../rscs/',args.aname+'_'+obj_name+'.pkl')
	pickle.dump(dat_obj,open(outfname,'wb'))
	print(obj_name+'=='+outfname)
