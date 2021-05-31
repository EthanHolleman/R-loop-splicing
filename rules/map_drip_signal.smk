

rule map_drip_siganl:
    conda:
        '../envs/bedtools.yml'
    input:
        SE='output/SE_calcs/SE_calc.{tNet_sample}.{strand}.bed',
        DRIP='data/DRIPc/DRIPc_{strand}_macs.bed'
    output:
        'output/map_drip/map_drip.{tNet_sample}.{strand}.bed'
    shell:'''
    bedtools map -a {input.SE} -b {input.DRIP} -c 5 -o mean /
    > {output} && [[ -s {output} ]]
    '''

rule map_drip_to_all_samples:
    input:
        expand(
            'output/map_drip/map_drip.{tNet_sample}.{strand}.bed',
            tNet_sample=tNet_samples.index.values.tolist(),
            strand=['fwd', 'rev']
            )
    output:
        'output/map_drip/map_drip_to_all_samples.done'
    shell:'''
    touch {output}
    '''
