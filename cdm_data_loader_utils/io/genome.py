from Bio import SeqIO


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
