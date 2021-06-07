rule get_intron_adar_editing_WT:
    input:
        'data/RNAss_processed_data/intron_editing.csv'
    output:
        'output/proccessed_SE/ADAR_editing.WT.scores.bed'
    shell:'''
    awk  -F',' '{{print $1 "\t" $2 "\t" $3 "\t" "intron_" NR "\t" $4 "\t" $5}}' {input} > {output}
    '''


rule get_intron_adar_editing_slow_pol:
    input:
        'data/RNAss_processed_data/SE_calcs.csv'
    output:
        'output/proccessed_SE/SE_intron_coords.slow.scores.bed'
    shell:'''
    awk  -F',' '{{print $1 "\t" $2 "\t" $3 "\t" "intron_" NR "\t" $4 "\t" $6}}' {input} > {output}
    '''


rule seperate_adar_editing_strands:
    input:
        'output/proccessed_SE/SE_intron_coords.{pol_type}.scores.bed'
    output:
        fwd='output/proccessed_SE/SE_intron_coords.{pol_type}.scores.fwd.bed',
        rev='output/proccessed_SE/SE_intron_coords.{pol_type}.scores.rev.bed'
    shell:'''
    awk '$6 == "+" {{print $0}}' {input} > {output.fwd} && [[ -s {output.fwd} ]]
    awk '$6 == "+" {{print $0}}' {input} > {output.rev} && [[ -s {output.rev} ]]
    '''


    
