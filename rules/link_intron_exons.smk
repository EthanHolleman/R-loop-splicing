

rule intersect_intron_exon_bounds:
    input:
        scored_introns='output/intron_exons_tNetseq/intron.{strand}.{tNet_sample}.bedgraph',
        scored_exons='output/intron_exons_tNetseq/exon.{strand}.{tNet_sample}.bedgraph'
    output:
        'output/linked_introns/linked_introns.{tNet_sample}.{strand}.bedgraph'
    shell:'''
    bedtools window -a {input.scored_introns} -b {input.scored_exons} -w 1 > {output} && [[ -s {output} ]]
    '''

rule link_all_intron_files:
    input:
        expand(
            'output/linked_introns/linked_introns.{tNet_sample}.{strand}.bedgraph',
            tNet_sample=tNet_samples.index.values.tolist(),
            strand=['fwd', 'rev']
        )
    output:
        'output/linked_introns/link_all_intron_files.done'
    shell:'''
    touch {output}
    '''



