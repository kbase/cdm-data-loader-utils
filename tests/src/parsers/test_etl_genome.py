from cdm_data_loader_utils.io.genome import ContigFeaturePairsBuilder
from cdm_data_loader_utils.etl_genome import ETLGbffGenome


def test_etl_genome_gbff():
    pairs = ContigFeaturePairsBuilder.from_gbff('tests/data/genome/GCF_000005845.2_ASM584v2_genomic.gbff.gz')
    etl = ETLGbffGenome()
    etl.etl('my_genome', 'some_protocol', pairs)

    assert etl._cdm_contigset is not None
    assert len(etl._data_cdm_contig) == 1
    assert len(etl._data_cdm_contigset_x_contig) == 1
    assert len(etl._data_cdm_protein) == 4257
