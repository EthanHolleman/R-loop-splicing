
rule sort_drip_files:
    conda:
        '../envs/bedtools.yml'
    input:
        'data/DRIPc/DRIPc_{strand}_macs.bed'
    output:
        'data/DRIPc/DRIPc_{strand}_macs.sorted.bed'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''


rule map_drip_signal:
    conda:
        '../envs/bedtools.yml'
    input:
        SE='output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.{intron_inclusion}.sorted.bed',
        DRIP='data/DRIPc/DRIPc_{strand}_macs.sorted.bed'
        #bedgraph='output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.bedgraph'
    output:
        'output/map_drip/map_drip.{calc_type}.{pol_type}.{strand}.{intron_inclusion}.bed'
    shell:'''
    bedtools map -a {input.SE} -b {input.DRIP} -c 5 -null 0 -o mean > {output} && [[ -s {output} ]]
    '''
