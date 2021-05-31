
rule calculate_splicing_efficiency:
    input:
        'output/linked_introns/linked_introns.{tNet_sample}.{strand}.bedgraph'
    output:
        'output/SE_calcs/SE_calc.{tNet_sample}.{strand}.bed'
    script:'../scripts/calculate_splicing_efficiency.py'

