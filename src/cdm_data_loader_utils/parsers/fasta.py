"""FASTA file parser."""

import gzip

from cdm_data_loader_utils.model.feature import Feature

DEFAULT_SPLIT = " "


def extract_features(faa_str: str, split: str = DEFAULT_SPLIT, h_func=None) -> list[Feature]:
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


def read_fasta(f: str, split: str = DEFAULT_SPLIT, h_func=None) -> list[Feature]:
    if f.endswith(".gz"):
        with gzip.open(f, "rb") as fh:
            return extract_features(fh.read().decode("utf-8"), split, h_func)

    with open(f) as fh:
        return extract_features(fh.read(), split, h_func)
