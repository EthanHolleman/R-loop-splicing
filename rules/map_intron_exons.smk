

rule sort_intron_exons:
    conda:
        '../envs/bedtools.yml'
    input:
        'data/hg19/knownGene_labeled.{region_type}.{strand}.formated.clean.bedgraph'
    output:
        'data/hg19/knownGene_labeled.{region_type}.{strand}.formated.sorted.bedgraph'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''


rule sort_tnet_seq:
    conda:
        '../envs/bedtools.yml'
    input:
        'data/RNAss/{tNet_sample}.bedgraph'
    output:
        'data/RNAss/{tNet_sample}.sorted.bedgraph'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''


rule map_intron_exons_over_nascent_transcripts:
    # determine the number of nascent transcripts measured by tNet-seq over
    # intron and exon regions in hg19
    conda:
        '../envs/bedtools.yml'
    input:
        genic_region='data/hg19/knownGene_labeled.{region_type}.{strand}.formated.sorted.bedgraph',
        nacent_read_count='data/RNAss/{tNet_sample}.sorted.bedgraph'
    output:
        'output/intron_exons_tNetseq/{region_type}.{strand}.{tNet_sample}.bedgraph'
    shell:'''
    mkdir -p output/intron_exons_tNetseq
    bedtools map -a {input.genic_region} -b {input.nacent_read_count} \
    -c 4 -o sum > {output} && [[ -s {output} ]]
    '''

# equilavent for drip peaks but use mean
rule map_intron_exons_over_drip_signal:
    input:
        genic_region='data/hg19/knownGene_labeled.{region_type}.{strand}.formated.bedgraph',


