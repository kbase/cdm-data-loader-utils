from Bio import SeqIO
from cdm_data_loader_utils.core.genome import CDMContig


def read_gbff_records_from_file(filename: str):
    if filename.endswith(".gbff"):
        with open(filename, 'r') as fh:
            return read_gbff_records(fh)
    elif filename.endswith(".gz"):
        import gzip
        from io import StringIO
        with gzip.open(filename, 'rb') as fh:
            return read_gbff_records(StringIO(fh.read().decode('utf-8')))


def read_gbff_records(handler):
    gbff_records = []
    for record in SeqIO.parse(handler, "gb"):
        gbff_records.append(record)
    return gbff_records


class ContigFeaturePairsBuilder:

    @staticmethod
    def get_feature_contig_id(f):
        ids = {x[0] for x in f.location}
        if len(ids) == 1:
            return list(ids)[0]
        else:
            raise ValueError('unable to fetch contig id')

    @staticmethod
    def from_kbase(genome_assembly, genome_features):
        contig_to_features = {}
        for contig in genome_assembly.features:
            contig_to_features[contig.id] = []
        for feature in genome_features.features:
            contig_id = ContigFeaturePairsBuilder.get_feature_contig_id(feature)
            contig_to_features[contig_id].append(feature)
        contig_feature_pairs = []
        for contig_id, features in contig_to_features.items():
            contig = genome_assembly.features.get_by_id(contig_id)
            cdm_contig = CDMContig(str(contig.seq))
            contig_feature_pairs.append((cdm_contig, features))

        return contig_feature_pairs

    @staticmethod
    def from_gbff(filename):
        gbff_records = read_gbff_records_from_file(filename)
        contig_feature_pairs = []
        for rec in gbff_records:
            cdm_contig = CDMContig(str(rec.seq))
            contig_feature_pairs.append((cdm_contig, [f for f in rec.features]))
        return contig_feature_pairs
