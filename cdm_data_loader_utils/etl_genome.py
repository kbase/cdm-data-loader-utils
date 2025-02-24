import pyarrow as pa
import pandas as pd
import os
from modelseedpy import MSGenome
from cdm_data_loader_utils.core.genome import (
    CDMContig,
    CDMContigSet,
    CDMFeature,
    CDMProtein,
)
from tqdm import tqdm


def _get_feature_contig_kbase(feature):
    contigs = {x[0] for x in feature.location}
    if len(contigs) == 1:
        return list(contigs)[0]
    raise ValueError(f'Feature: {feature.id} - invalid contig location {contigs}')


def _get_feature_contig(feature):
    _id = feature.id
    _contig_name = _id[::-1].split("_", 1)[1][::-1]  # reverse split once reverse again
    return _contig_name


def _get_start_end_strand_s_attr(feature):
    start, end, strand_s, attr = feature.description.split(" # ")
    if start.startswith("# "):
        start = int(start[1:].strip())
    end = int(end)
    strand = "?"
    if strand_s == "1":
        strand = "+"
    elif strand_s == "-1":
        strand = "-"
    else:
        raise ValueError(f"bad strand: {strand_s}")

    if attr:
        attr = dict(x.split("=") for x in attr.split(";"))
    else:
        attr = {}

    return start, end, strand, attr


def _get_start_end_strand_s_attr_kbase(feature):
    for l in feature.location:
        # print(l)
        contig, start, strand, end = l
        return start, end, strand, {}

    # start, end, strand_s, attr = f.description.split(" # ")
    return


def _get_feature_id(feature):
    return feature.id


class ETLGbffGenome:

    def __init__(self):
        self._data_cdm_contig = {}
        self._data_cdm_contigset_x_contig = set()
        self._cdm_contigset = None
        self._data_cdm_protein = {}
        self._data_cdm_feature = {}
        self._data_cdm_feature_x_protein = set()

    def etl(self, genome_id, protocol_id, gbff_records):
        contigs = []
        h_to_gbff_record = {}
        for contig in gbff_records:
            cdm_contig = CDMContig(str(contig.seq))
            if cdm_contig.hash_contig not in self._data_cdm_contig:
                self._data_cdm_contig[cdm_contig.hash_contig] = cdm_contig
            if cdm_contig.hash_contig not in h_to_gbff_record:
                h_to_gbff_record[cdm_contig.hash_contig] = []
            h_to_gbff_record[cdm_contig.hash_contig].append(contig)
            contigs.append(cdm_contig)
        self._cdm_contigset = CDMContigSet.from_contigs(contigs)

        repeat_counter = {}
        for hash_contig in h_to_gbff_record:
            cdm_contig = self._data_cdm_contig[hash_contig]
            index_contig = 0
            contigset_x_contig_id = ETLGbffGenome.get_contigset_x_contig_id(self._cdm_contigset, cdm_contig,
                                                                            index_contig)

            self._data_cdm_contigset_x_contig.add((
                contigset_x_contig_id,
                self._cdm_contigset.hash_contigset,
                cdm_contig.hash_contig,
                index_contig
            ))

            for gbff_record in h_to_gbff_record[hash_contig]:
                for gbff_feature in list(gbff_record.features)[:]:
                    feature_id = ETLGbffGenome._get_feature_id(genome_id, gbff_feature)
                    if feature_id in self._data_cdm_feature:
                        counter = repeat_counter.get(feature_id, 1)
                        repeat_counter[feature_id] = counter + 1
                        feature_id = f'{feature_id}_{counter}'
                    cdm_feature = ETLGbffGenome.convert_to_cdm_feature(feature_id, contigset_x_contig_id, gbff_feature)
                    # set protocol
                    cdm_feature.protocol = protocol_id
                    # set source
                    cdm_feature.source = 'gbff'
                    if len(cdm_feature.id) < 62 and cdm_feature.id not in self._data_cdm_feature:
                        self._data_cdm_feature[cdm_feature.id] = cdm_feature
                    else:
                        if len(cdm_feature.id) >= 62:
                            raise ValueError(f'x feature id {cdm_feature.id}')
                        else:
                            raise ValueError(f'duplicate feature id {cdm_feature.id}')
                    if gbff_feature.type == 'CDS':
                        cdm_protein = ETLGbffGenome.get_protein(gbff_feature)
                        if cdm_protein.hash not in self._data_cdm_protein:
                            self._data_cdm_protein[cdm_protein.hash] = cdm_protein
                            self._data_cdm_feature_x_protein.add((cdm_feature.id, cdm_protein.hash))
                index_contig += 1

    @staticmethod
    def _get_feature_id(genome_id, f):
        start = int(f.location.start)
        end = int(f.location.end)
        feature_id = f'{genome_id}_{f.type}_{start}_{end}'
        if 'locus_tag' in f.qualifiers and len(f.qualifiers['locus_tag']) == 1:
            locus = f.qualifiers['locus_tag'][0]
            return f'{genome_id}_{locus}_{f.type}_{start}_{end}'
        return feature_id

    @staticmethod
    def get_contigset_x_contig_id(cdm_contigset, cdm_contig, index_contig):
        contigset_x_contig_str = f'{cdm_contigset.hash_contigset}_{cdm_contig.hash_contig}_{index_contig}'
        return ETLGbffGenome._hash_string(contigset_x_contig_str)

    @staticmethod
    def get_protein(f):
        translation = f.qualifiers.get('translation')
        if translation and type(translation) == list and len(translation) == 1:
            return CDMProtein(translation[0])
        raise ValueError('unable to get translation')

    @staticmethod
    def get_strand(f_strand):
        if f_strand == 1:
            return '+'
        elif f_strand == -1:
            return '-'
        raise ValueError('strand error')

    @staticmethod
    def convert_to_cdm_feature(feature_id, contigset_x_contig_id, f):
        start = int(f.location.start)
        end = int(f.location.end)
        strand = ETLGbffGenome.get_strand(f.location.strand)
        feature_type = f.type
        cdm_feature = CDMFeature(feature_id, contigset_x_contig_id, None,
                                 start, end, strand,
                                 feature_type)
        return cdm_feature

    @staticmethod
    def _hash_string(s):
        import hashlib
        return hashlib.sha256(s.encode("utf-8")).hexdigest()


class MapAssemblyToGff:

    def __init__(self, assembly, gff_records):
        from modelseedpy_ext.re.hash_seq import HashSeq
        self.contig_id_to_feature_contig = {}
        self.gff_index_to_contig = {}
        self.contig_to_hash_seq = {}

        from collections import Counter
        counter = Counter()
        for gff_record in gff_records:
            counter.update([gff_record.feature_type])

        for gff_record in gff_records:
            if gff_record.contig_id not in self.contig_id_to_feature_contig:
                self.contig_id_to_feature_contig[gff_record.contig_id] = None

        if 'DNA' in counter:
            for gff_record in gff_records:
                if gff_record.feature_type == 'DNA':
                    feature_contig = MapAssemblyToGff.get_feature_contig(gff_record, assembly)
                    if gff_record.contig_id in self.contig_id_to_feature_contig:
                        self.contig_id_to_feature_contig[gff_record.contig_id] = feature_contig
                    if 'Alias' in gff_record.attr and gff_record.attr.get('Alias') in self.contig_id_to_feature_contig:
                        alias = gff_record.attr.get('Alias')
                        if self.contig_id_to_feature_contig[alias] is None:
                            self.contig_id_to_feature_contig[alias] = feature_contig
                        elif self.contig_id_to_feature_contig[alias] != feature_contig:
                            raise ValueError('!')
        else:
            for f in assembly.features:
                if f.id in self.contig_id_to_feature_contig:
                    self.contig_id_to_feature_contig[f.id] = f

        for k in self.contig_id_to_feature_contig:
            feature_contig = self.contig_id_to_feature_contig[k]
            if feature_contig is None:
                raise ValueError(f'Unmapped contig: {k}')
            if feature_contig.id not in self.contig_to_hash_seq:
                hash_contig = HashSeq(feature_contig.seq).hash_value
                self.contig_to_hash_seq[feature_contig.id] = hash_contig

        for i in range(len(gff_records)):
            gff_record = gff_records[i]
            feature_contig = self.contig_id_to_feature_contig[gff_record.contig_id]
            if feature_contig is None:
                raise ValueError(f'Unable to find contig: {gff_record}')

            self.gff_index_to_contig[i] = feature_contig

    @staticmethod
    def get_feature_contig(gff_record, assembly):
        contig_id = gff_record.contig_id
        if contig_id in assembly.features:
            return assembly.features.get_by_id(contig_id)
        else:
            contig_id = gff_record.attr['Note'].split(',')[0][:-1]
            if contig_id in assembly.features:
                return assembly.features.get_by_id(contig_id)
            else:
                contig_id = gff_record.attr['Note'].split(" ")[0]
                return assembly.features.get_by_id(contig_id)


class ETLAlexeyGenome:
    _S = len('/scratch/shared/CDM/alexey_genomes_ranjan/genomes_final/genome_processing/')
    _STRIP_LEN = -1 * len('_contigs.fasta')

    def __init__(self, sql_engine, exclude=None):
        self.cdm_contigsets = {}
        self.cdm_name_contigset = set()
        self.cdm_name_contig = set()
        self.genome_scan = {}
        self.data_contigs = {}
        self.data_features = {}
        self.data_feature_to_protein = {}
        self.cdm_name_feature = set()
        self.cdm_feature_annotation = {}
        self.sql_engine = sql_engine
        self.exclude = exclude
        if self.exclude is None:
            self.exclude = set()

    @staticmethod
    def swap_dict(d):
        res = {}
        for k, v in d.items():
            if v not in res:
                res[v] = k
            else:
                raise ValueError('not one to one')
        return res

    @staticmethod
    def _a_sz(a):
        _sz = 0
        for f in a.features:
            _sz += len(f.seq)
        return _sz

    def scan_genomes(self, path_contigs):
        for path_contig in tqdm(path_contigs):
            self.scan_genome(path_contig)

        df_alexey_genome = pd.read_sql("SELECT * FROM browser_genome;", self.sql_engine)

        alexey_unique_genome_names = set(df_alexey_genome['name'])
        missing = alexey_unique_genome_names - set(self.genome_scan)
        not_in_alexey = set(self.genome_scan) - alexey_unique_genome_names
        both = set(self.genome_scan) & alexey_unique_genome_names
        print(len(alexey_unique_genome_names), len(missing), len(not_in_alexey), len(both))

    def scan_genome(self, path_contig):
        _prefix_filename = path_contig[:ETLAlexeyGenome._STRIP_LEN]
        path_faa = _prefix_filename + '_Prodigal_proteins.faa'
        path_gff = _prefix_filename + '_Prodigal.gff'

        if os.path.exists(path_faa) and os.path.exists(path_gff):
            _p = path_contig[ETLAlexeyGenome._S:].split('/')
            base_name = _p[0]
            base_name_v = _p[2]
            if base_name_v not in self.genome_scan:
                self.genome_scan[base_name_v] = {
                    'path_contigs.fasta': path_contig,
                    'path_faa': path_faa,
                    'path_gff': path_gff
                }
            else:
                print('!0', path_contig)

    def process_genomes(self):
        df_alexey_genome = pd.read_sql("SELECT * FROM browser_genome;", self.sql_engine)

        alexey_unique_genome_names = set(df_alexey_genome['name'])
        missing = alexey_unique_genome_names - set(self.genome_scan)
        not_in_alexey = set(self.genome_scan) - alexey_unique_genome_names
        both = set(self.genome_scan) & alexey_unique_genome_names
        print(len(alexey_unique_genome_names), len(missing), len(not_in_alexey), len(both))

        for row_id, d in tqdm(df_alexey_genome.iterrows()):
            genome_name = d['name']
            if genome_name in both and genome_name not in self.exclude:
                _scan_data = self.genome_scan[genome_name]
                assembly = REAssembly.from_fasta(_scan_data['path_contigs.fasta'])
                assembly_size = ETLAlexeyGenome._a_sz(assembly)
                contigs = len(assembly.features)
                contigs_ids = [x.id for x in assembly.features]

                match = assembly_size == d['size'] and contigs == d['contigs']
                if match:
                    id_genome_alexey = d['id']
                    self.process_contigset(assembly, genome_name, id_genome_alexey)

                    df_contigs = pd.read_sql(
                        f"SELECT * FROM browser_contig WHERE genome_id = {id_genome_alexey};",
                        self.sql_engine)
                    d_id_contig_alexey_to_contig_name = ETLAlexeyGenome.swap_dict(
                        df_contigs.set_index('id')['contig_id'].to_dict())
                    #genome = MSGenome.from_fasta(_scan_data['path_faa'])
                    gff_records = _read_gff_features(_scan_data['path_gff'])

                    from collections import Counter
                    counter = Counter()
                    for gff_record in gff_records:
                        counter.update([gff_record.feature_type])
                    gff_feature_types = dict(counter)

                    try:
                        if 'DNA' in gff_feature_types:
                            self.process_contigs(gff_records, assembly, d_id_contig_alexey_to_contig_name)
                        else:
                            self.process_contig_no_gff(assembly, d_id_contig_alexey_to_contig_name)
                    except Exception as ex:
                        print(f'I failed to parse this: {genome_name}')
                        raise ValueError(ex)
                else:
                    print(genome_name)

    def process_contigset(self, assembly, contigset_name, alexey_id):
        hash_contigset = assembly.hash_value
        if hash_contigset not in self.cdm_contigsets:
            self.cdm_contigsets[hash_contigset] = {
                'n_contigs': len(assembly.features),
                'length': ETLAlexeyGenome._a_sz(assembly),
                'alexey_db_id': [alexey_id]
            }
        else:
            self.cdm_contigsets[hash_contigset]['alexey_db_id'].append(alexey_id)

        self.cdm_name_contigset.add((hash_contigset, 'contigset', contigset_name, 'alexey_db'))

    def process_contig_no_gff(self, assembly, contig_alias_to_alexey_id):
        """
        process contigs without gff

        The alexey dict maps the alternative alias to the surrogate key from the alexey database
        :param gff_records:
        :param assembly:
        :param contig_alias_to_alexey_id:
        :return:
        """
        scaffold_alias_to_contig_id = {}

        hash_contigset = assembly.hash_value

        # iterate all records look for DNA type feature and extract Alias Name and the original feature_id
        #  from the Note attr
        for feature_contig in assembly.features:
            contig_id = feature_contig.id
            cdm_contig = CDMContig(hash_contigset, feature_contig.seq)
            pk = (cdm_contig.hash_contig, cdm_contig.contig_set_id)

            if pk not in self.data_contigs:
                self.data_contigs[pk] = cdm_contig
                self.data_contigs[pk].alexey_id = set()
                self.data_contigs[pk].names = set()

            self.data_contigs[pk].alexey_id.add(contig_alias_to_alexey_id[contig_id])
            self.data_contigs[pk].names.add((
                cdm_contig.hash_contig, cdm_contig.contig_set_id,
                'contig', contig_id, 'enigma', 'from _contig.fasta'))

    @staticmethod
    def get_feature_contig(gff_record, assembly):
        contig_id = gff_record.contig_id
        if contig_id in assembly.features:
            return assembly.features.get_by_id(contig_id)
        else:
            contig_id = gff_record.attr['Note'].split(',')[0][:-1]
            if contig_id in assembly.features:
                return assembly.features.get_by_id(contig_id)
            else:
                contig_id = gff_record.attr['Note'].split(" ")[0]
                return assembly.features.get_by_id(contig_id)

    def process_contigs(self, gff_records, assembly, contig_alias_to_alexey_id):
        """
        process contigs we need both the assembly object and the gff records
        the gff records contain the alternative aliases given by KBase Prokka while
        the assembly object has the contig sequence and the original feature_id

        The alexey dict maps the alternative alias to the surrogate key from the alexey database
        :param gff_records:
        :param assembly:
        :param contig_alias_to_alexey_id:
        :return:
        """
        scaffold_alias_to_contig_id = {}

        hash_contigset = assembly.hash_value

        # iterate all records look for DNA type feature and extract Alias Name and the original feature_id
        #  from the Note attr
        for gff_record in gff_records:
            if gff_record.feature_type == 'DNA':

                alias = gff_record.attr.get('Alias')
                gff_contig_id = gff_record.contig_id
                alexey_db_id = contig_alias_to_alexey_id.get(alias)
                if alexey_db_id is None:
                    alexey_db_id = contig_alias_to_alexey_id.get(gff_contig_id)

                if alias not in scaffold_alias_to_contig_id:
                    feature_contig = ETLAlexeyGenome.get_feature_contig(gff_record, assembly)

                    cdm_contig = CDMContig(hash_contigset, feature_contig.seq)
                    pk = (cdm_contig.hash_contig, cdm_contig.contig_set_id)

                    if pk not in self.data_contigs:
                        self.data_contigs[pk] = cdm_contig
                        self.data_contigs[pk].alexey_id = set()
                        self.data_contigs[pk].names = set()

                    self.data_contigs[pk].alexey_id.add(alexey_db_id)

                    self.data_contigs[pk].names.add((
                        cdm_contig.hash_contig, cdm_contig.contig_set_id,
                        'contig', alias, 'enigma', 'from _Prodigal.gff'
                    ))
                    self.data_contigs[pk].names.add((
                        cdm_contig.hash_contig, cdm_contig.contig_set_id,
                        'contig', feature_contig.id, 'enigma', 'from _contigs.fasta'
                    ))
                    if gff_record.attr['Name'] != alias:
                        self.data_contigs[pk].names.add((
                            cdm_contig.hash_contig, cdm_contig.contig_set_id,
                            'contig',
                            gff_record.attr['Name'],
                            'enigma',
                            'from _contigs.fasta'
                        ))
                else:
                    raise ValueError('!')

    @staticmethod
    def build_feature_name(genome_id, gff_record):
        """
        if 'Name' in gff_record.attr:
            attr_name = gff_record.attr.get('Name')
            return f"{genome_id}_{gff_record.feature_type}_{attr_name}"
        elif 'ID' in gff_record.attr:
            attr_id = gff_record.attr.get('ID')
            return f"{genome_id}_{gff_record.feature_type}_{attr_id}"
        else:
        """
        return f"{genome_id}_{gff_record.contig_id}_{gff_record.feature_type}_{gff_record.start}_{gff_record.end}_{gff_record.strand}"

    @staticmethod
    def find_feature(genome, gff_record):
        if 'Name' in gff_record.attr:
            _name = gff_record.attr.get('Name')
            if _name in genome.features:
                return genome.features.get_by_id(_name)
            else:
                _id = gff_record.attr.get('ID')
                if _id in genome.features:
                    return genome.features.get_by_id(_id)
                else:
                    return genome.features.get_by_id(_id.split('.')[0])

        _id = gff_record.attr.get('ID')
        if _id in genome.features:
            return genome.features.get_by_id(_id)

        search_id = '_'.join(gff_record.contig_id.split('_')[:-1]) + '_' + gff_record.attr.get('ID')
        if search_id in genome.features:
            return genome.features.get_by_id(search_id)

        search_id = gff_record.contig_id[:-1] + gff_record.attr.get('ID')
        return genome.features.get_by_id(search_id)

    def process_features(self, genome_id, gff_records, assembly, genome, protocol_id):
        from modelseedpy_ext.re.hash_seq import HashSeq
        feature_function = {}  # product

        def _add_name(cdm_feature, name):
            if name:
                self.cdm_name_feature.add((cdm_feature.id, 'feature', name, 'enigma', 'from _Prodigal.gff'))

        hash_contigset = assembly.hash_value
        assembly_to_gff = MapAssemblyToGff(assembly, gff_records)

        for i in range(len(gff_records)):
            gff_record = gff_records[i]

            feature_contig = assembly_to_gff.gff_index_to_contig[i]
            hash_contig = assembly_to_gff.contig_to_hash_seq[feature_contig.id]

            feature_id = ETLAlexeyGenome.build_feature_name(genome_id, gff_record)
            cdm_feature = CDMFeature(feature_id, hash_contig, hash_contigset,
                                     gff_record.start, gff_record.end, gff_record.strand,
                                     gff_record.feature_type, gff_record.source, gff_record.phase,
                                     protocol_id,
                                     gff_record.attr)

            _add_name(cdm_feature, gff_record.attr.get('Name'))
            _add_name(cdm_feature, gff_record.attr.get('gene_synonym'))
            annotation = gff_record.attr.get('product')
            if annotation:
                self.cdm_feature_annotation[cdm_feature.id] = annotation

            if cdm_feature.id not in self.data_features:
                self.data_features[cdm_feature.id] = cdm_feature
            else:
                raise ValueError(f'id error {cdm_feature.id}')

            if gff_record.feature_type == 'CDS':
                feature_protein = ETLAlexeyGenome.find_feature(genome, gff_record)
                _add_name(cdm_feature, feature_protein.id)

                cdm_protein = CDMProtein(feature_protein.seq)
                # map encoded protein
                self.data_feature_to_protein[(cdm_feature.id, cdm_protein.hash)] = cdm_protein.stop_codon

            # print(feature_id, name_to_contig_hash[gff_record.contig_id])
            # print(gff_record.attr)

    def export_data_contigs(self):
        names = ["hash_contig", "hash_contigset", "length", "gc_content", "base_count"]
        _data = {k: [] for k in names}

        for o in tqdm(self.data_contigs.values()):
            _data["hash_contigset"].append(o.contig_set_id)
            _data["hash_contig"].append(o.hash_contig)
            _data["length"].append(len(o.seq))
            _data["gc_content"].append(o.gc)
            _data["base_count"].append(o.base_count)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_name_contigs(self):
        _table_name_names = ["id_entity", "hash_contigset", "name_table", "name", "source", "description"]
        _table_name_data = {k: [] for k in _table_name_names}
        names = {}
        for o in tqdm(self.data_contigs.values()):
            for hash_contig, hash_contigset, name_table, name, source, description in o.names:
                names[(hash_contig, name, source)] = (
                    hash_contig, hash_contigset, name_table, name, source, description)
        for hash_contig, hash_contigset, name_table, name, source, description in names.values():
            _table_name_data['id_entity'].append(hash_contig)
            _table_name_data['hash_contigset'].append(hash_contigset)
            _table_name_data['name_table'].append(name_table)
            _table_name_data['name'].append(name)
            _table_name_data['source'].append(source)
            _table_name_data['description'].append(description)
        table = pa.Table.from_arrays([pa.array(_table_name_data[k]) for k in _table_name_names],
                                     names=_table_name_names)
        return table


class EtlAssemblyAndGenome:

    def __init__(self, fn_get_feature_contig=None, fn_get_start_end_strand_s_attr=None, fn_get_feature_id=None):

        self.data_contig_set = []
        self.data_contigs = []

        self.contig_name_to_hash = {}

        self.data_proteins = {}
        self.data_features = {}
        self.data_feature_to_protein = {}

        self.names = set()

        self._fn_get_feature_contig = fn_get_feature_contig
        self._fn_get_start_end_strand_s_attr = fn_get_start_end_strand_s_attr
        if self._fn_get_feature_contig is None:
            self._fn_get_feature_contig = _get_feature_contig
        if self._fn_get_start_end_strand_s_attr is None:
            self._fn_get_start_end_strand_s_attr = _get_start_end_strand_s_attr
        self._fn_get_feature_id = _get_feature_id
        if fn_get_feature_id:
            self._fn_get_feature_id = fn_get_feature_id

    def etl_assembly(self, g_assembly, contig_id_name_source='ncbi'):
        contig_set_hash = g_assembly.hash_value

        cdm_contig_set = CDMContigSet(contig_set_hash)

        for feature in g_assembly.features:
            cdm_contig = CDMContig(contig_set_hash, feature.seq)
            if contig_id_name_source:
                cdm_contig.names.append((contig_id_name_source, feature.id))
            cdm_contig_set.contigs.append(cdm_contig)
            self.data_contigs.append(cdm_contig)

            for _t, _n in cdm_contig.names:
                self.names.add((cdm_contig.hash_contig, _t, _n))
                if _n not in self.contig_name_to_hash:
                    self.contig_name_to_hash[_n] = cdm_contig.hash_contig
                elif self.contig_name_to_hash[_n] != cdm_contig.hash_contig:
                    raise ValueError("dup")

        self.data_contig_set.append(cdm_contig_set)

        return contig_set_hash

    def features(self, genome: MSGenome, gff):
        for f in genome.features:
            yield f

    def etl_genome(self, genome: MSGenome, hash_contigset, gff=None):

        for f in self.features():
            feature_id = self._fn_get_feature_id(f)
            _contig_name = self._fn_get_feature_contig(f)
            _contig_name_hash = self.contig_name_to_hash[_contig_name]

            start, end, strand, attr = self._fn_get_start_end_strand_s_attr(f, gff)

            cdm_feature = CDMFeature(
                feature_id,
                _contig_name_hash,
                hash_contigset,
                start,
                end,
                strand,
                attr,
            )
            cdm_protein = CDMProtein(f.seq)
            self.data_feature_to_protein[cdm_feature.id] = (
                cdm_protein.hash,
                cdm_protein.stop_codon,
            )
            if cdm_feature.id not in self.data_features:
                self.data_features[cdm_feature.id] = cdm_feature
            else:
                raise ValueError("dup feature id !~!")
            if cdm_protein.hash not in self.data_proteins:
                self.data_proteins[cdm_protein.hash] = cdm_protein

    def export_protein(self):
        names = ["hash_protein_sequence", "description", "length", "sequence", "evidence_for_existence"]
        _data = {k: [] for k in names}
        for i in tqdm(self.data_proteins):
            o = self.data_proteins[i]
            _data["hash_protein_sequence"].append(o.hash)
            _data["length"].append(len(o.seq))
            _data["sequence"].append(o.seq)
            _data["description"].append(None)
            _data["evidence_for_existence"].append(None)

        pa_table = pa.Table.from_arrays(
            [pa.array(_data[k]) for k in names], names=names
        )
        return pa_table

    def export_names(self):
        _label_names = ["entity_id", "type", "name"]
        _data_names = {k: [] for k in _label_names}
        for _hash, _type, _value in tqdm(self.names):
            _data_names["entity_id"].append(_hash)
            _data_names["type"].append(_type)
            _data_names["name"].append(_value)

        table = pa.Table.from_arrays(
            [pa.array(_data_names[k]) for k in _label_names], names=_label_names
        )
        return table

    def export_data_contigset(self):
        names = ["hash_contigset", "n_contigs"]
        _data = {k: [] for k in names}
        for o in tqdm(self.data_contig_set):
            _data["hash_contigset"].append(o.hash_contigset)
            _data["n_contigs"].append(len(o.contigs))

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_contigs(self):
        names = ["hash_contigset", "hash_contig", "length", "gc_content", "base_count"]
        _data = {k: [] for k in names}

        _label_names = ["entity_id", "type", "name"]
        _data_names = {k: [] for k in _label_names}
        for o in tqdm(self.data_contigs):
            _data["hash_contigset"].append(o.contig_set_id)
            _data["hash_contig"].append(o.hash_contig)
            _data["length"].append(len(o.seq))
            _data["gc_content"].append(o.gc)
            _data["base_count"].append(o.base_count)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_feature(self):
        names = [
            "feature_id",
            "hash_contig",
            "hash_contigset",
            "hash",
            "source",
            "type",
            "start",
            "end",
            "strand",
            "cds_phase",
            "attributes",
        ]
        _data = {k: [] for k in names}
        for i in tqdm(self.data_features):
            o = self.data_features[i]
            _data["feature_id"].append(o.id)
            _data["hash_contig"].append(o.hash_contig)
            _data["hash_contigset"].append(o.hash_contigset)
            _data["hash"].append(None)
            _data["source"].append(o.source)
            _data["type"].append("CDS")
            _data["start"].append(o.start)
            _data["end"].append(o.end)
            _data["strand"].append(o.strand)
            _data["cds_phase"].append(None)
            _data["attributes"].append(o.attributes)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

    def export_data_feature_x_protein(self):
        names = ["feature_id", "hash_protein_sequence", "has_stop_codon"]
        _data = {k: [] for k in names}
        for feature_id in tqdm(self.data_feature_to_protein):
            protein_hash, has_stop_codon = self.data_feature_to_protein[feature_id]
            _data["feature_id"].append(feature_id)
            _data["hash_protein_sequence"].append(protein_hash)
            _data["has_stop_codon"].append(has_stop_codon)

        table = pa.Table.from_arrays([pa.array(_data[k]) for k in names], names=names)
        return table

