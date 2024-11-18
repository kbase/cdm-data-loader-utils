from collections.abc import Iterable
import hashlib

DEFAULT_SPLIT = " "

def _hash_string(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()

class HashSeq(str):
    def __new__(cls, v):
        instance = super().__new__(cls, v.upper())
        return instance

    @property
    def hash_value(self):
        return _hash_string(self)

class HashSeqList(list):
    def append(self, o):
        if isinstance(o, str):
            super().append(HashSeq(o))
        elif isinstance(o, HashSeq):
            super().append(o)
        else:
            raise ValueError("bad type")

    @property
    def hash_value(self):
        h_list = [x.hash_value for x in self]
        hash_seq = "_".join(sorted(h_list))
        return _hash_string(hash_seq)

def extract_features(faa_str, split=DEFAULT_SPLIT, h_func=None):
    features = []
    active_seq = None
    seq_lines = []
    for line in faa_str.split("\n"):
        if line.startswith(">"):
            if active_seq is not None:
                active_seq.seq = "".join(seq_lines)
                features.append(active_seq)
                seq_lines = []
            seq_id = line[1:]
            desc = None
            if h_func:
                seq_id, desc = h_func(seq_id)
            elif split:
                header_data = line[1:].split(split, 1)
                seq_id = header_data[0]
                if len(header_data) > 1:
                    desc = header_data[1]
            active_seq = Feature(seq_id, "", desc)
        else:
            seq_lines.append(line.strip())
    if len(seq_lines) > 0:
        active_seq.seq = "".join(seq_lines)
        features.append(active_seq)
    return features

def read_fasta2(f, split=DEFAULT_SPLIT, h_func=None):
    if f.endswith(".gz"):
        import gzip
        with gzip.open(f, "rb") as fh:
            return extract_features(fh.read().decode("utf-8"), split, h_func)
    else:
        with open(f, "r") as fh:
            return extract_features(fh.read(), split, h_func)

class Feature:
    def __init__(self, feature_id, sequence, description=None, aliases=None):
        self.id = feature_id
        self.seq = sequence
        self.description = description
        self.ontology_terms = {}
        self.aliases = aliases

    def add_ontology_term(self, ontology_term, value):
        if ontology_term not in self.ontology_terms:
            self.ontology_terms[ontology_term] = []
        if value not in self.ontology_terms[ontology_term]:
            self.ontology_terms[ontology_term].append(value)

def contig_set_hash(features):
    hl = HashSeqList()
    for contig in features:
        seq = HashSeq(contig.seq)
        hl.append(seq)
    return hl.hash_value

