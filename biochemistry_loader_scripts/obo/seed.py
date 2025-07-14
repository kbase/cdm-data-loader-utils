#!/usr/bin/env python

"""Converter for SEED Subsystems."""

from collections.abc import Iterable

from pyobo.struct import Obo, Reference, Term

from pyobo.struct.typedef import (
	has_participant,
	participates_in
)

import json

__all__ = [
	"SEEDSubsystemGetter",
]

SUBSYSTEM_PREFIX = "seed.subsystem"
ROLE_PREFIX = "seed.role"
PREFIX=SUBSYSTEM_PREFIX

# Date that data was retrieved
VERSION = "2025-01-06"

class SEEDSubsystemGetter(Obo):
	"""An ontology representation of the SEED Subsystems and Functional Roles."""

	ontology = bioversions_key = PREFIX
	name = PREFIX
	idspaces = {"seed.role":"https://pubseed.theseed.org/RoleEditor.cgi?page=ShowRole&Role=",
				"seed.subsystem":"https://pubseed.theseed.org/SubsysEditor.cgi?page=ShowSubsystem&subsystem=",
				"seed.reaction":"https://modelseed.org/biochem/reactions/"}

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
	"""Get SEED Subsystems OBO."""
	return SEEDSubsystemGetter(force=force)

def get_terms(*, version: str, use_tqdm: bool = True, force: bool = False) -> Iterable[Term]:
	"""Get SEED Subsystem terms."""

	terms=list()
	
	# retrieve roles reactions dbxrefs
	roles_reactions = dict()
	with open('Roles_Reactions_20250224.txt') as fh:
		for line in fh.readlines():
			line=line.strip('\r\n')
			tmp_lst=line.split('\t')
			if(tmp_lst[0] not in roles_reactions):
				roles_reactions[tmp_lst[0]]=list()
			if(tmp_lst[1] not in roles_reactions[tmp_lst[0]]):
				roles_reactions[tmp_lst[0]].append(tmp_lst[1])

	# retrieve identifiers from Chris' file
	subsystem_id_dict = dict()
	role_id_dict = dict()
	last_role_id = 0
	last_subsystem_id = 0
	with open('SSO_dictionary.json') as fh:
		tmp_dict=json.load(fh)

		for sso_id in tmp_dict['term_hash']:
			role_id = int(sso_id.split(':')[1])
			if(role_id > last_role_id):
				last_role_id = role_id
			role_id_dict[tmp_dict['term_hash'][sso_id]['name']]=str(role_id).zfill(13)

	subsystems_dict = dict()
	roles_dict = dict()
	with open('2018-0531-Subsystem-Roles.txt') as fh:
		for line in fh.readlines():
			line=line.strip('\r\n')

			(subsystem,role)=line.split('\t')

			if(subsystem != "-" and subsystem not in subsystems_dict):

				subsystem_id = None
				if(subsystem in subsystem_id_dict):
					subsystem_id = subsystem_id_dict[subsystem]
				else:
					last_subsystem_id +=1
					subsystem_id = str(last_subsystem_id).zfill(10)
					subsystem_id_dict[subsystem]=subsystem_id

				ss_ref = Reference(prefix=SUBSYSTEM_PREFIX, identifier=subsystem_id, name=subsystem)
				ss_term = Term(reference=ss_ref)
				subsystems_dict[subsystem]={'ref':ss_ref,'term':ss_term,'roles':[]}
			
			if(role not in roles_dict):
				role_id = None
				if(role in role_id_dict):
					role_id = role_id_dict[role]
				else:
					last_role_id += 1
					role_id = str(last_role_id).zfill(13)
					role_id_dict[role]=role_id

				role_ref = Reference(prefix=ROLE_PREFIX, identifier=role_id, name=role)
				role_term = Term(reference=role_ref)

				if(role in roles_reactions):
					for rxn in roles_reactions[role]:
						xref_ref = Reference(prefix='seed.reaction', identifier=rxn)
						role_term.append_xref(xref_ref)
				
				roles_dict[role]={'ref':role_ref,'term':role_term}

				if(subsystem != "-" and role not in subsystems_dict[subsystem]['roles']):
					subsystems_dict[subsystem]['roles'].append(role)
	
	for subsystem in subsystems_dict:
		ss_term = subsystems_dict[subsystem]['term']
		ss_ref = subsystems_dict[subsystem]['ref']

		for role in subsystems_dict[subsystem]['roles']:
			role_ref = roles_dict[role]['ref']
			ss_term.append_relationship(has_participant, role_ref)

			role_term = roles_dict[role]['term']
			role_term.append_relationship(participates_in, ss_ref)
		terms.append(ss_term)

	for role in roles_dict:
		role_term = roles_dict[role]['term']
		terms.append(role_term)

	print(len(list(subsystems_dict.keys())),len(list(roles_dict.keys())))
	return terms

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('file')
	# args = parser.parse_args()
	# print(args['file'])

	ontology = get_obo()
	with open('seed.obo','w') as fh:
		for line in ontology.iterate_obo_lines():
			print(line, file=fh)

	with open('seed_sssom.tsv','w') as fh:
		print('\t'.join(['subject_id','object_id','predicate_id','mapping_justification']), file=fh)
		for entity in ontology.iterate_xref_rows():
			# default values as per get_sssom_df
			row = ['','','oboInOwl:hasDbXref','sempav:UnspecifiedMatching']
			row[0] = 'seed.role:'+entity[0]
			row[1] = 'seed.reaction:'+entity[2]
			print('\t'.join(row), file=fh)
