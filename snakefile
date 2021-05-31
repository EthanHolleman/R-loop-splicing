import pandas as pd

tNet_samples = pd.read_table(
    'samples/tNetSeq_samples.tsv'
    ).set_index('sample_name', drop=False)


include: 'rules/link_intron_exons.smk'
include: 'rules/map_intron_exons.smk'
include: 'rules/prep_intron_exon_regions.smk'
include: 'rules/prep_RNA_data.smk'
include: 'rules/map_drip_signal.smk'
include: 'rules/calculate_se.smk'



rule all:
    input:
        'data/hg19/knownGene_labeled.csv',
        expand(
            'data/hg19/knownGene_labeled.{region_type}.{strand}.formated.bedgraph',
            region_type=['intron', 'exon'], strand=['fwd', 'rev']
        ),
        'output/linked_introns/link_all_intron_files.done',
        'output/map_drip/map_drip_to_all_samples.done'