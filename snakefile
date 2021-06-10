import pandas as pd


wildcard_constraints:
    strand='\w+',
    calc_type='\w+',
    pol_type='\w+',
    intron_inclusion='\w+'


include: 'rules/map_drip_signal.smk'
include: 'rules/plots.smk'
include: 'rules/prep_precalc_data.smk'
include: 'rules/prep_hg19_genes.smk'
include: 'rules/truncate_first_last_introns.smk'


rule all:
    input:
        'output/plots/plot.pdf',