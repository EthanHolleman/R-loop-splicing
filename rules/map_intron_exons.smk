 

rule map_intron_exons_over_nascent_transcripts:
    conda:
        '../envs/bedtools.yml'
    input:
        genic_region='data/hg19/knownGene_labeled.{region_type}.{strand}.formated.bedgraph',
        nacent_read_count='data/RNAss/{tNet_sample}.bedgraph'
    output:
        'output/intron_exons_tNetseq/{region_type}.{strand}.{tNet_sample}.bedgraph'
    shell:'''
    bedtools map -a {input.genic_region} -b {input.nacent_read_count} \
    -c 4 -o sum > {output}
    '''

# equilavent for drip peaks but use mean
rule map_intron_exons_over_drip_signal:
    input:
        genic_region='data/hg19/knownGene_labeled.{region_type}.{strand}.formated.bedgraph',


