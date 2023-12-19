# written to generate an excel file of
# drugbank to tarkget gene names 
# written 5-18-20 JLW

import pickle,os,csv
import pandas as pd
from collections import defaultdict

df1 = '../data/drugbank_050120.tsv'
df2 = '../data/proteins_050120.tsv'
mapf = '../data/uniprot_to_gene_name_table_050620.txt'

# just the DrugBank Data
outf = '../data/Drugbank050120.xlsx'
writer = pd.ExcelWriter(outf, engine='xlsxwriter')
df = pd.read_csv(df1,sep="\t")
df.to_excel(writer,sheet_name="DrugBankSummary")
db2name = dict(zip(df.drugbank_id, df.name))

# get drug->target, also map to gene symbol
mapd = dict([(r['From'],r['To']) for r in csv.DictReader(open(mapf,'r'),delimiter='\t')])
df = pd.read_csv(df2,sep="\t")
df['geneSym'] = df['uniprot_id'].map(mapd).fillna('NotMapped')
df['DrugName'] = df['drugbank_id'].map(db2name).fillna('NotMapped')
df.to_excel(writer,sheet_name="DrugsToTargets")

# now assemble all targets for a drug, or all drugs for a target
df['allTargets'] = df[['drugbank_id','geneSym']].groupby(['drugbank_id'])['geneSym'].transform(lambda x: ','.join(x))
df2 = df[['drugbank_id','allTargets']].drop_duplicates()
df2['DrugName'] = df2['drugbank_id'].map(db2name).fillna('NotMapped')
df2.to_excel(writer,sheet_name="Drugs2Targets")

df['allDrugs'] = df[['DrugName','geneSym']].groupby(['geneSym'])['DrugName'].transform(lambda x: ','.join(x))
df3 = df[['geneSym','allDrugs']].drop_duplicates()
df3.to_excel(writer,sheet_name="Targets2Drugs")
writer.save()
