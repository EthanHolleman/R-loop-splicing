

rule sort_genes:
    conda:
        '../envs/bedtools.yml'
    input:
        'data/hg19/hg19_apprisplus_gene.bed'  # do not delete
    output:
        'data/hg19/hg19_apprisplus_gene.sorted.bed'
    shell:'''
    sort-bed --max-mem 8G {input} > {output} && [[ -s {output} ]]
    '''

