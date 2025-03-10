from cdm_data_loader_utils.io.genome import read_gbff_records_from_file


def test_read_gbff_records_from_file():
    records = read_gbff_records_from_file('tests/data/genome/GCF_000005845.2_ASM584v2_genomic.gbff.gz')

    assert len(records) == 1
