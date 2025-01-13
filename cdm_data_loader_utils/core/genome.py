from modelseedpy.core.msgenome import MSGenome, read_fasta2
from cdm_data_loader_utils.core.hash_seq import HashSeq, HashSeqList
from collections import Counter


class CDMContigSet:

    def __init__(self, sha256):
        self.sha256 = sha256
        self.contigs = []


class CDMContig:

    def __init__(self, contig_set_id: str, seq: str):
        self.seq = seq
        self.contig_set_id = contig_set_id
        self.hash = HashSeq(self.seq).hash_value
        self.base_count = dict(Counter(list(self.seq.upper())))
        self.length = len(self.seq)
        self.gc = (self.base_count.get('G', 0) + self.base_count.get('C', 0)) / self.length

        self.names = []

    def __repr__(self):
        return f'len: {self.length}, gc: {self.gc}, base_count: {self.base_count}, names: {self.names}'


class CDMProtein:

    def __init__(self, seq: str):
        self.stop_codon = False
        _seq = seq
        if _seq[-1] == '*':
            _seq = _seq[:-1]
            self.stop_codon = True

        self.seq = _seq
        self.hash = HashSeq(self.seq).hash_value
        self.length = len(self.seq)

        self.names = []

    def __repr__(self):
        return f'len: {self.length}, hash: {self.hash}'


class CDMFeature:

    def __init__(self, feature_id: str, contig_hash, start, end, strand, attributes=None):
        self.id = feature_id
        self.contig_hash = contig_hash
        self.start = start
        self.end = end
        self.strand = strand
        self.type = None
        self.source = None
        self.cds_phase = None
        self.attributes = {} if attributes is None else attributes

        self.names = []


class GffRecord:

    def __init__(self, contig_id: str, source: str,
                 feature_type,
                 start: int, end: int, score, strand, phase, attr):
        self.contig_id = contig_id
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attr = attr

    def get_attribute_string(self):
        attr_values = []
        for k, v in self.attr.items():
            attr_values.append(f"{k}={v}")
        return ';'.join(attr_values)

    def __str__(self):
        return '\t'.join([str(x) for x in [self.contig_id, self.source, self.feature_type,
                                           self.start, self.end, self.score, self.strand, self.phase,
                                           self.get_attribute_string()]])

    @staticmethod
    def from_str(s):
        contig_id, source, feature_type, start, end, score, strand, phase, attr_str = s.strip().split('\t')
        attr = dict([x.split('=') for x in attr_str.split(';')])
        return GffRecord(contig_id, source, feature_type, int(start), int(end), score, strand, phase, attr)


class REAssembly(MSGenome):

    def __init__(self):
        super().__init__()
        self.hash_list = HashSeqList()
        for contig in self.features:
            seq = HashSeq(contig.seq)
            self.hash_list.append(seq)

    def re(self):
        pass

    def ke(self):
        pass

    @staticmethod
    def from_fasta(filename, split=" ", h_func=None):
        genome = REAssembly()
        genome.features += read_fasta2(filename, split, h_func)
        return genome

    @property
    def hash_value(self):
        hl = HashSeqList()
        for contig in self.features:
            seq = HashSeq(contig.seq)
            hl.append(seq)
        return hl.hash_value

    @staticmethod
    def _process_contigs(contigs):
        hash_list = HashSeqList()
        contig_h_d = []
        for contig in contigs.features:
            seq = HashSeq(contig.seq)
            hash_list.append(seq)
            seq_h = seq.hash_value
            contig_h_d.append([seq_h, contig.id, contig.description])
        return {
            'genome_h': hash_list.hash_value,
            'contig_h': contig_h_d
        }
