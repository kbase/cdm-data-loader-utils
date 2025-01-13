import pyarrow as pa
from modelseedpy import MSGenome
from cdm_data_loader_utils.core.genome import CDMContig, CDMContigSet, CDMFeature, CDMProtein
from tqdm import tqdm


class EtlAssemblyAndGenome:

    def __init__(self):

        self.data_contig_set = []
        self.data_contigs = []

        self.contig_name_to_hash = {}

        self.data_proteins = {}
        self.data_features = {}
        self.data_feature_to_protein = {}

        self.names = set()

    def etl_assembly(self, g_assembly):
        contig_set_hash = g_assembly.hash_value

        cdm_contig_set = CDMContigSet(contig_set_hash)

        for feature in g_assembly.features:
            cdm_contig = CDMContig(contig_set_hash, feature.seq)
            cdm_contig.names.append(('ncbi', feature.id))
            cdm_contig_set.contigs.append(cdm_contig)
            self.data_contigs.append(cdm_contig)

            for _t, _n in cdm_contig.names:
                self.names.add((cdm_contig.hash, _t, _n))
                if _n not in self.contig_name_to_hash:
                    self.contig_name_to_hash[_n] = cdm_contig.hash
                elif self.contig_name_to_hash[_n] != cdm_contig.hash:
                    raise ValueError('dup')

        self.data_contig_set.append(cdm_contig_set)

    def etl_genome(self, genome: MSGenome):
        for f in genome.features:
            _id = f.id
            _contig_name = _id[::-1].split('_', 1)[1][::-1]  # reverse split once reverse again
            _contig_name_hash = self.contig_name_to_hash[_contig_name]

            start, end, strand_s, attr = f.description.split(' # ')
            if start.startswith('# '):
                start = int(start[1:].strip())
            end = int(end)
            strand = '?'
            if strand_s == '1':
                strand = '+'
            elif strand_s == '-1':
                strand = '-'
            else:
                raise ValueError(f'bad strand: {strand_s}')

            cdm_feature = CDMFeature(f.id, _contig_name_hash, start, end, strand,
                                     dict(x.split('=') for x in attr.split(';')))
            cdm_protein = CDMProtein(f.seq)
            self.data_feature_to_protein[cdm_feature.id] = (cdm_protein.hash, cdm_protein.stop_codon)
            if cdm_feature.id not in self.data_features:
                self.data_features[cdm_feature.id] = cdm_feature
            else:
                raise ValueError('dup feature id !~!')
            if cdm_protein.hash not in self.data_proteins:
                self.data_proteins[cdm_protein.hash] = cdm_protein

    def export_protein(self):
        names = ['hash', 'description', 'length', 'sequence', 'evidence_for_existence']
        _data = {k: [] for k in names}
        for i in tqdm(self.data_proteins):
            o = self.data_proteins[i]
            _data['hash'].append(o.hash)
            _data['length'].append(len(o.seq))
            _data['sequence'].append(o.seq)
            _data['description'].append(None)
            _data['evidence_for_existence'].append("Protein predicted (Prodigal)")

        pa_table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return pa_table

    def export_names(self):
        _label_names = ['entity_id', 'type', 'name']
        _data_names = {k: [] for k in _label_names}
        for _hash, _type, _value in tqdm(self.names):
            _data_names['entity_id'].append(_hash)
            _data_names['type'].append(_type)
            _data_names['name'].append(_value)

        table = pa.Table.from_arrays([pa.array(_data_names[k]) for k in _label_names], names=_label_names)
        return table

    def export_data_contigset(self):
        names = ['hash', 'n_contigs']
        _data = {k: [] for k in names}
        for o in tqdm(self.data_contig_set):
            _data['hash'].append(o.sha256)
            _data['n_contigs'].append(len(o.contigs))

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_contigs(self):
        names = ['contig_set_hash', 'contig_hash', 'length', 'gc_content', 'base_count']
        _data = {k: [] for k in names}

        _label_names = ['entity_id', 'type', 'name']
        _data_names = {k: [] for k in _label_names}
        for o in tqdm(self.data_contigs):
            _data['contig_set_hash'].append(o.contig_set_id)
            _data['contig_hash'].append(o.hash)
            _data['length'].append(len(o.seq))
            _data['gc_content'].append(o.gc)
            _data['base_count'].append(o.base_count)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_feature(self):
        names = ['id', 'contig_hash', 'hash',
                 'source_database', 'source_algorithm',
                 'type', 'start', 'end', 'strand', 'cds_phase', 'attributes']
        _data = {k: [] for k in names}
        for i in tqdm(self.data_features):
            o = self.data_features[i]
            _data['id'].append(o.id)
            _data['contig_hash'].append(o.contig_hash)
            _data['hash'].append(None)
            _data['source_database'].append('ke-pangenomes')
            _data['source_algorithm'].append('prodigal_ke_pangenomes')
            _data['type'].append('CDS')
            _data['start'].append(o.start)
            _data['end'].append(o.end)
            _data['strand'].append(o.strand)
            _data['cds_phase'].append(None)
            _data['attributes'].append(o.attributes)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_encoded_feature(self):
        names = ['feature_id', 'protein_hash', 'has_stop_codon']
        _data = {k: [] for k in names}
        for feature_id in tqdm(self.data_feature_to_protein):
            protein_hash, has_stop_codon = self.data_feature_to_protein[feature_id]
            _data['feature_id'].append(feature_id)
            _data['protein_hash'].append(protein_hash)
            _data['has_stop_codon'].append(has_stop_codon)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table
