#!/usr/bin/env python

"""Converter for MetaCyc Biochemistry."""

from collections.abc import Iterable

from pyobo.struct import Obo, Reference, Term

from pyobo.struct.typedef import (
	has_participant,
	participates_in
)

import json

__all__ = [
	"MetaCycBiochemistryGetter",
]

PATHWAY_PREFIX = "metacyc.pathway"
REACTION_PREFIX = "metacyc.reaction"
PREFIX=PATHWAY_PREFIX
VERSION = "26.1"

class MetaCycBiochemistryGetter(Obo):
	"""An ontology representation of the MetaCyc Biochemistry database."""

	ontology = bioversions_key = PREFIX
	name = PREFIX

	def __post_init__(self):
		self.data_version = VERSION

	typedefs = [
		has_participant,
		participates_in
	]

	def iter_terms(self, force: bool = False) -> Iterable[Term]:
		"""Iterate over terms in the ontology."""
		return get_terms(version=self._version_or_raise, force=force)

def get_obo(force: bool = False) -> Obo:
	"""Get MetaCyc Biochemistry OBO."""
	return MetaCycBiochemistryGetter(force=force)

def get_terms(*, version: str, use_tqdm: bool = True, force: bool = False) -> Iterable[Term]:
	"""Get MetaCyc Biochemistry terms."""

	header = 1
	mc_pwys = list()
	mc_rxns = dict()
	mc_pwy_parents = dict()
	with open('../Data/MetaCyc_pathways.tsv') as fh:
		for line in fh.readlines():
			
			if(header==1):
				header=0
				continue

			line=line.strip('\r\n')
			tmp_lst=line.split('\t')

			# for i in range(len(tmp_lst)):
			#	print(i,tmp_lst[i])

			pwy_dict = {'id':tmp_lst[0],'name':tmp_lst[1],'reactions':[]}
			
			if(tmp_lst[2] != ''):
				reactions = tmp_lst[2].split('|')
				
				pwy_dict['reactions']=reactions

				for rxn in pwy_dict['reactions']:
					if(rxn not in mc_rxns):
						rxn_ref = Reference(prefix=REACTION_PREFIX, identifier=rxn)
						rxn_term = Term(reference=rxn_ref)
						mc_rxns[rxn]={'ref':rxn_ref,'term':rxn_term,'pathways':[pwy_dict['id']]}
					else:
						mc_rxns[rxn]['pathways'].append(pwy_dict['id'])

			if(tmp_lst[3] != ''):
				for parent in tmp_lst[3].split('|'):
					if(pwy_dict['id'] not in mc_pwy_parents):
						mc_pwy_parents[pwy_dict['id']]=[]
					if(parent not in mc_pwy_parents[pwy_dict['id']]):
						mc_pwy_parents[pwy_dict['id']].append(parent)

			mc_pwys.append(pwy_dict)

	terms=list()
	pathway_refs = dict()
	pathway_terms = dict()
	for pwy in mc_pwys:
		pwy_ref = Reference(prefix=PATHWAY_PREFIX, identifier=pwy['id'], name=pwy['name'])
		pathway_refs[pwy['id']]=pwy_ref
		pwy_term = Term(reference=pwy_ref)
		pathway_terms[pwy['id']]=pwy_term

		for rxn in pwy['reactions']:
			rxn_ref = mc_rxns[rxn]['ref']
			pwy_term.append_relationship(has_participant,rxn_ref)
		terms.append(pwy_term)

	for rxn in mc_rxns:
		rxn_term = mc_rxns[rxn]['term']
		for pwy in mc_rxns[rxn]['pathways']:
			pwy_ref = pathway_refs[pwy]
			rxn_term.append_relationship(participates_in,pwy_ref)

		terms.append(rxn_term)

	for child in mc_pwy_parents:
		child_term = pathway_terms[child]
		child_ref = pathway_refs[child]
		for parent in mc_pwy_parents[pwy]:
			parent_ref = pathway_refs[parent]
			child_term.append_relationship(participates_in,parent_ref) 

			parent_term = pathway_terms[parent]
			parent_term.append_relationship(has_participant,child_ref)

	print(len(mc_pwys),len(mc_rxns))

	return terms

if __name__ == "__main__":
	ontology = get_obo()
	with open('metacyc.obo','w') as fh:
		for line in ontology.iterate_obo_lines():
			print(line, file=fh)