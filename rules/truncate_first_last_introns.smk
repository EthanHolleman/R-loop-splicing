# Remove the first and last intron from a file. This is done when considering
# correlations to DRIP data as most signal will be at the start and end
# of genes.


rule intersect_introns_and_genes:
    # intersect intron regions over genes in a strand specific manor
    conda:
        '../envs/bedtools.yml'
    input:
        introns='output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.allIntrons.sorted.bed',
        genes='data/hg19/hg19_apprisplus_gene.sorted.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.geneIntronIntersect.bed'
    shell:'''
    bedtools intersect -a {input.genes} -b {input.introns} -s -sorted -wa -wb > {output}
    '''


rule label_introns_with_gene_and_number:
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.geneIntronIntersect.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.labledAllIntrons.bed'
    script:'../scripts/label_introns.py'


rule sort_label_introns_with_gene_and_number:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.labledAllIntrons.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.labledAllIntrons.sorted.bed'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''


rule remove_first_last_introns:
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.labledAllIntrons.sorted.bed'
    output:
        no_first_last='output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.noFirstLast.bed',
    script:'../scripts/remove_first_last_introns_from_intersection.py'



rule sort_truncated_introns:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.noFirstLast.bed'
    output:
        'output/pre_calc_processed/{calc_type}.{pol_type}.{strand}.noFirstLast.sorted.bed'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''





