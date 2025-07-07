#!/usr/bin/env python
from urllib.request import urlopen
import json,os,glob

# biochemistry repository cloned from https://github.com/ModelSEED/ModelSEEDDatabase
# Using dev branch as of 01/09/25 (commit a00786)

biochem_root = "/Users/seaver/Seaver_Lab/Git_Repos/ModelSEEDDatabase/Biochemistry/"

# reactions for sample data is a couple of steps in the glycolytic pathway
# rxn00781 rxn01100 rxn01106 rxn00459
sample_reactions = ["rxn00459","rxn00781","rxn01100","rxn01106"]
sample_reactions_list = list()
 
search_path = os.path.join(biochem_root,"reaction_*.json")
all_rxns_list = list()
for reactions_file in sorted(glob.glob(search_path)):
	with open(reactions_file) as json_file_handle:
		rxns_list = json.load(json_file_handle)
		all_rxns_list.extend(rxns_list)

		for rxn in rxns_list:
			if(rxn['id'] in sample_reactions):
				sample_reactions_list.append(rxn)

with open("MSD_Reactions.json",'w') as fh:
	fh.write(json.dumps(all_rxns_list))

with open("Sample_MSD_Reactions.json",'w') as fh:
	fh.write(json.dumps(sample_reactions_list))

# sample compounds will be taken from sample reactions
sample_compounds = list()
for rxn in sample_reactions_list:
	for cpd in rxn['compound_ids'].split(';'):
		if(cpd not in sample_compounds):
			sample_compounds.append(cpd)

sample_compounds = sorted(sample_compounds)
sample_compounds_list = list()

search_path = os.path.join(biochem_root,"compound_*.json")
all_cpds_list = list()
for compounds_file in sorted(glob.glob(search_path)):
	with open(compounds_file) as json_file_handle:
		cpds_list = json.load(json_file_handle)
		all_cpds_list.extend(cpds_list)

		# no need to keep doing this if we've found them all
		if(len(sample_compounds_list)==len(sample_compounds)):
			continue

		for cpd in cpds_list:
			if(cpd['id'] in sample_compounds):
				sample_compounds_list.append(cpd)

with open("MSD_Compounds.json",'w') as fh:
	fh.write(json.dumps(all_cpds_list))

with open("Sample_MSD_Compounds.json",'w') as fh:
	fh.write(json.dumps(sample_compounds_list))
