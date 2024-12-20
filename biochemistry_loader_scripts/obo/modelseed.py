#!/usr/bin/env python

"""Converter for ModelSEED Biochemistry."""

import logging
from collections.abc import Iterable

from tqdm.auto import tqdm

from pyobo.struct import Obo, Reference, Term

from pyobo.struct.typedef import (
	has_participant,
	participates_in
)

import json

__all__ = [
	"ModelSEEDBiochemistryGetter",
]

logger = logging.getLogger(__name__)

REACTION_PREFIX = "seed.reaction"
COMPOUND_PREFIX = "seed.compound"
PREFIX = REACTION_PREFIX
VERSION = "2.0"

class ModelSEEDBiochemistryGetter(Obo):
	"""An ontology representation of the ModelSEED Biochemistry database."""

	ontology = bioversions_key = PREFIX
	name = PREFIX

	def __post_init__ (self):
		self.data_version = VERSION

	typedefs = [
		has_participant,
		participates_in
	]

	def iter_terms(self, force: bool = False) -> Iterable[Term]:
		"""Iterate over terms in the ontology."""
		return get_terms(version=self._version_or_raise, force=force)

def get_obo(force: bool = False) -> Obo:
	"""Get ModelSEED Biochemistry OBO."""
	return ModelSEEDBiochemistryGetter(force=force)

def get_terms(*, version: str, use_tqdm: bool = True, force: bool = False) -> Iterable[Term]:
	"""Get ModelSEED Biochemistry terms."""

	ms_cpds = list()
	with open('../Data/MSD_Compounds.json') as fh:
		ms_cpds = json.load(fh)

	cpds_in_refs = dict()
	cpds_in_terms = dict()
	cpd_aliases = dict()
	for cpd in ms_cpds:

		cpd_ref = Reference(prefix=COMPOUND_PREFIX, identifier=cpd['id'], name=cpd['name'])
		cpds_in_refs[cpd['id']]=cpd_ref
		cpd_term = Term(reference=cpd_ref)

		if(cpd['aliases'] is not None):
			for alias in cpd['aliases']:
				xref = alias.split(': ')[-1]
				for db in ['MetaCyc','KEGG','Rhea','ChEBI']:
					if(db in alias):
						if(db not in cpd_aliases):
							cpd_aliases[db]=0
						
						db_prefix = db.lower()
						if(db=='MetaCyc' or db=='KEGG'):
							db_prefix+='.compound'

						xref_ref = Reference(prefix=db_prefix, identifier=xref)
						cpd_term.append_xref(xref_ref)

						cpd_aliases[db]+=1

		cpds_in_terms[cpd['id']]=cpd_term

	ms_rxns = list()
	with open('../Data/MSD_Reactions.json') as fh:
		ms_rxns = json.load(fh)

	rxns_in_refs=dict()
	rxns_in_terms=dict()
	rxn_aliases=dict()
	for rxn in ms_rxns:
		rxn_ref = Reference(prefix=REACTION_PREFIX, identifier=rxn['id'], name=rxn['name'])
		rxns_in_refs[rxn['id']]=rxn_ref

		rxn_term = Term(reference=rxn_ref)

		if(rxn['aliases'] is not None):
			for alias in rxn['aliases']:
				xref = alias.split(': ')[-1]
				for db in ['MetaCyc','KEGG','Rhea']:
					if(db in alias):

						if(db not in rxn_aliases):
							rxn_aliases[db]=0
						
						db_prefix = db.lower()
						if(db=='MetaCyc' or db=='KEGG'):
							db_prefix+='.reaction'

						xref_ref = Reference(prefix=db_prefix, identifier=xref)
						rxn_term.append_xref(xref_ref)

						rxn_aliases[db]+=1
					
		rxns_in_terms[rxn['id']]=rxn_term

	rxn_cpds = dict()
	for rxn in ms_rxns:
		rxn_cpds[rxn['id']]=list()
		rxn_term = rxns_in_terms[rxn['id']]
		rxn_ref = rxns_in_refs[rxn['id']]
		if(rxn['compound_ids'] is None):
			continue

		for cpd in rxn['compound_ids'].split(';'):
			rxn_cpds[rxn['id']].append(cpd)
			cpd_ref = cpds_in_refs[cpd]
			rxn_term.append_relationship(has_participant,cpd_ref)

			cpd_term = cpds_in_terms[cpd]
			cpd_term.append_relationship(participates_in,rxn_ref)

	terms = list()
	cpds_to_use = list()
	for rxn in rxns_in_terms:
		terms.append(rxns_in_terms[rxn])

		for cpd in rxn_cpds[rxn]:
			cpds_to_use.append(cpd)
	
	for cpd in cpds_in_terms:
		terms.append(cpds_in_terms[cpd])
	
	print(len(rxns_in_terms),len(cpds_in_terms))

	for db in cpd_aliases:
		print("cpd:"+db+":"+str(cpd_aliases[db]))

	for db in rxn_aliases:
		print("rxn:"+db+":"+str(rxn_aliases[db]))

	return terms

if __name__ == "__main__":
	ontology = get_obo()
	with open('modelseed.obo','w') as fh:
		for line in ontology.iterate_obo_lines():
			print(line, file=fh)

	# https://pyobo.readthedocs.io/en/stable/api/pyobo.get_sssom_df.html
	# ontology.get_sssom_df()

	with open('modelseed_sssom.tsv','w') as fh:
		print('\t'.join(['subject_id','object_id','predicate_id','mapping_justification']), file=fh)
		for entity in ontology.iterate_xref_rows():
			# default values as per get_sssom_df
			row = ['','','oboInOwl:hasDbXref','sempav:UnspecifiedMatching']
			if('cpd' in entity[0]):
				row[0] = 'seed.compound:'+entity[0]
			if('rxn' in entity[0]):
				row[0] = 'seed.reaction:'+entity[0]
			row[1] = entity[1]+':'+entity[2]
			print('\t'.join(row), file=fh)