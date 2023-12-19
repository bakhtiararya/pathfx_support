# written to run all of DrugBank through
# each threshold of PathFX
# written 1/21/20 JLW

import pickle, os, csv
from collections import defaultdict

# run all of DrugBank???
dbf = '../rscs/pathfx_mp_dbid2name.pkl'
db2name = pickle.load(open(dbf,'rb'))
n2db = dict([(dn.lower(),db_id) for (db_id,dn) in db2name.items()])

# copy/paste code to isolate drugs with interactome targets
f = os.path.join('..','data','Drugs_labeled_for_AEs.txt') #Each column header is a DME, and the rows are drugs associated with the DME
dR = csv.DictReader(open(f,'r'),delimiter='\t')

drugs_by_dme = defaultdict(list)
dmes_by_drug = defaultdict(list)
for r in dR:
        for (dme,drug) in r.items():
                if drug != '':
			if drug.lower() in n2db:
				dbid = n2db[drug.lower()]
				drugs_by_dme[dme].append(dbid)
				dmes_by_drug[dbid].append(dme)
unique_drugs = [k for k in dmes_by_drug.keys()]
unique_dmes = [k for k in drugs_by_dme.keys()]

fpdmes_by_drug = defaultdict(list)
for (dbid,dme_list) in dmes_by_drug.items():
	fp_dmes = list(set(unique_dmes).difference(set(dme_list)))# DMES not listed on the drug's label
	fpdmes_by_drug[dbid] = fp_dmes

pickle.dump(unique_drugs,open('../data/tpfp_unique_drugs.pkl','wb'))
pickle.dump(unique_dmes,open('../data/tpfp_unique_dmes.pkl','wb'))
pickle.dump(drugs_by_dme,open('../data/tpfp_drugs_by_dme.pkl','wb')) 
pickle.dump(dmes_by_drug,open('../data/tpfp_dmes_by_drug.pkl','wb')) 
pickle.dump(fpdmes_by_drug,open('../data/tpfp_fpdmes_by_drug.pkl','wb')) 
