#!/usr/bin/env python

"""Converter for KEGG Biochemistry."""

from collections.abc import Iterable

from pyobo.struct import Obo, Reference, Term

from pyobo.struct.typedef import (
	has_participant,
	participates_in
)

import json

__all__ = [
	"KEGGBiochemistryGetter",
]

PATHWAY_PREFIX = "kegg.pathway"
REACTION_PREFIX = "kegg.reaction"
PREFIX=PATHWAY_PREFIX

# AS OF 2023-09-18
# https://www.kegg.jp/kegg/docs/relnote.html
VERSION = "107.1"

class KEGGBiochemistryGetter(Obo):
	"""An ontology representation of the KEGG Biochemistry database."""

	ontology = bioversions_key = PREFIX
	data_version = VERSION
	name = PREFIX

	idspaces = {"kegg.pathway":"https://rest.kegg.jp/find/pathway/",
			 	"kegg.reaction":"http://rest.kegg.jp/find/reaction/"}
	
	typedefs = [
		has_participant,
		participates_in
	]

	def iter_terms(self, force: bool = False) -> Iterable[Term]:
		"""Iterate over terms in the ontology."""
		return get_terms(version=self._version_or_raise, force=force)

def get_obo(force: bool = False) -> Obo:
	"""Get KEGG Biochemistry OBO."""
	return KEGGBiochemistryGetter(force=force)

def get_terms(*, version: str, use_tqdm: bool = True, force: bool = False) -> Iterable[Term]:
	"""Get KEGG Biochemistry terms."""

	header = 1
	kg_pwys = list()
	kg_rxns = dict()
	kg_pwy_parents = dict()
	with open('/Users/seaver/Seaver_Lab/Git_Repos/ModelSEEDDatabase/Scripts/Provenance/KEGG/KEGG_pathways.tsv') as fh:
		for line in fh.readlines():
			
			if(header==1):
				header=0
				continue

			line=line.strip('\r\n')
			tmp_lst=line.split('\t')

			# for i in range(len(tmp_lst)):
			#	print(i,tmp_lst[i])

			if(' ' in tmp_lst[0]):
				tmp_lst[0] = tmp_lst[0].replace(' ','_')

			pwy_dict = {'id':tmp_lst[0],'name':tmp_lst[1],'description':tmp_lst[2],'reactions':[]}

			if(tmp_lst[3] != ''):
				reactions = tmp_lst[3].split('|')

				pwy_dict['reactions']=reactions

				for rxn in pwy_dict['reactions']:
					if(rxn not in kg_rxns):
						rxn_ref = Reference(prefix=REACTION_PREFIX, identifier=rxn)
						rxn_term = Term(reference=rxn_ref)
						kg_rxns[rxn]={'ref':rxn_ref,'term':rxn_term,'pathways':[pwy_dict['id']]}
					else:
						kg_rxns[rxn]['pathways'].append(pwy_dict['id'])

			if(tmp_lst[4] != ''):
				for parent in tmp_lst[4].split('|'):
					if(' ' in parent):
						parent = parent.replace(' ','_')
					if(pwy_dict['id'] not in kg_pwy_parents):
						kg_pwy_parents[pwy_dict['id']]=[]
					if(parent not in kg_pwy_parents[pwy_dict['id']]):
						kg_pwy_parents[pwy_dict['id']].append(parent)

			kg_pwys.append(pwy_dict)

	terms=list()
	pathway_refs = dict()
	pathway_terms = dict()
	for pwy in kg_pwys:
		pwy_ref = Reference(prefix=PATHWAY_PREFIX, identifier=pwy['id'], name=pwy['name'])
		pathway_refs[pwy['id']]=pwy_ref
		pwy_term = Term(reference=pwy_ref)
		pathway_terms[pwy['id']]=pwy_term

		for rxn in pwy['reactions']:
			rxn_ref = kg_rxns[rxn]['ref']
			pwy_term.append_relationship(has_participant,rxn_ref)
		terms.append(pwy_term)

	for rxn in kg_rxns:
		rxn_term = kg_rxns[rxn]['term']
		for pwy in kg_rxns[rxn]['pathways']:
			pwy_ref = pathway_refs[pwy]
			rxn_term.append_relationship(participates_in,pwy_ref)

		terms.append(rxn_term)

	for child in kg_pwy_parents:
		child_term = pathway_terms[child]
		child_ref = pathway_refs[child]
		for parent in kg_pwy_parents[pwy]:
			parent_ref = pathway_refs[parent]
			child_term.append_relationship(participates_in,parent_ref) 

			parent_term = pathway_terms[parent]
			parent_term.append_relationship(has_participant,child_ref)

	print(len(kg_pwys),len(kg_rxns))
	return terms

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('file')
	# args = parser.parse_args()
	# print(args['file'])

	ontology = get_obo()
	with open('kegg.obo','w') as fh:
		for line in ontology.iterate_obo_lines():
			print(line, file=fh)
