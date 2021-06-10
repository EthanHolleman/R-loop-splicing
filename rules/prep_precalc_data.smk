# found that SE calcs in sup data
rule get_intron_cords_WT:
    input:
        'data/RNAss_processed_data/SE_calcs.csv'
    output:
        'output/pre_calc_processed/splicingEfficiency.WT.bed'
    shell:'''
    awk  -F',' '{{print $1 "\t" $2 "\t" $3 "\t" "intron_" NR "\t" $9 "\t" $8}}' {input} > {output}
    '''


rule get_intron_adar_editing_WT:
    # seperate wildtype edits into their own file
    input:
        'data/RNAss_processed_data/intron_editing.csv'
    output:
        'output/pre_calc_processed/ADARediting.WT.bed'
    shell:'''
    awk  -F',' '{{print $1 "\t" $2 "\t" $3 "\t" "intron_" NR "\t" $6 "\t" $5}}' {input} > {output}
    '''


rule get_intron_adar_editing_slow:
    input:
        'data/RNAss_processed_data/intron_editing.csv'
    output:
        'output/pre_calc_processed/ADARediting.slow.bed'
    shell:'''
    awk  -F',' '{{print $1 "\t" $2 "\t" $3 "\t" "intron_" NR "\t" $7 "\t" $5}}' {input} > {output}
    '''

rule get_intron_adar_editing_WT_total:
    input:
        'data/RNAss_processed_data/intron_editing.csv'
    output:
        'output/pre_calc_processed/ADAReditingTotalEdits.WT.bed'
    shell:'''
    awk  -F',' '{{print $1 "\t" $2 "\t" $3 "\t" "intron_" NR "\t" $4 "\t" $5}}' {input} > {output}
    '''


rule seperate_pre_calc_strands_fwd:
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.fwd.bed'
    shell:'''
    awk '$6 == "+" {{print $0}}' {input} > {output} && [[ -s {output} ]]
    '''


rule seperate_pre_calc_strands_rev:
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.rev.bed'
    shell:'''
    awk '$6 == "+" {{print $0}}' {input} > {output} && [[ -s {output} ]]
    '''


rule sort_precalc_strands:
    # sort and add file identified indicating that all introns
    # are included to differentiate from files where introns have
    # been truncated. see remove_first_last_introns.smk
    conda:
        '../envs/bedtools.yml'
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.allIntrons.sorted.bed'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''


# rule convert_to_bedgraph:
#     input:
#         'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.sorted.bed'
#     output:
#         'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.bedgraph'
#     shell:'''
#     awk '{{print $1 "\t" $2 "\t" $3 "\t" $5}}' {input} > {output}
#     '''



    


