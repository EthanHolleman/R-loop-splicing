

rule map_drip_siganl:
    input:
        SE='output/SE_calcs/SE_calc.{tNet_sample}.{strand}.bed',
        DRIP='data/DRIPc/DRIPc_{strand}_macs.bed'
    output:
        ''
    shell:'''
    bedtools map -a {input.SE} -b {input.DRIP} -c 5 -o mean /
    > {output}
    '''
