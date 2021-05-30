
rule download_intron_exon_regions_hg19:
    # https://www.biostars.org/p/411934/ formated for snakemake
    output:
        'data/hg19/knownGene_labeled.csv'
    shell:'''
    mkdir -p data/hg19
    curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" /
    | gunzip -c /
    | awk '{{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {{printf("%s,%s,%s,%s,%s,EXON\\n",$1,$2,$3,S[i],E[i]); if(i+1<=n) printf("%s,%s,%s,%s,%s,INTRON\\n",$1,$2,$3,int(E[i]),int(S[i+1]));  }}}}' /
    > {output}
    '''


rule seperate_introns_exons:
    input:
        'data/hg19/knownGene_labeled.csv'
    output:
        introns='data/hg19/knownGene_labeled.intron.csv',
        exons='data/hg19/knownGene_labeled.exon.csv'
    shell:'''
    awk -F, '$6 == "EXON" {{print $0}}' {input} > {output.exons} && [[ -s {output.exons} ]]
    awk -F, '$6 == "INTRON" {{print $0}}' {input} > {output.introns} && [[ -s {output.introns} ]]
    '''


rule seperate_intron_exon_strands_fwd:
    input:
        'data/hg19/knownGene_labeled.{region_type}.csv'
    output:
        'data/hg19/knownGene_labeled.{region_type}.fwd.csv'
    shell:'''
    awk -F, '$3 == "+" {{print $0}}' {input} > {output} && [[ -s {output} ]]
    '''


rule seperate_intron_exon_strands_rev:
    input:
        'data/hg19/knownGene_labeled.{region_type}.csv'
    output:
        'data/hg19/knownGene_labeled.{region_type}.rev.csv'
    shell:'''
     awk -F, '$3 == "-" {{print $0}}' {input} > {output} && [[ -s {output} ]]
    '''


rule format_intron_exon_csv_as_bed:
    input:
        'data/hg19/knownGene_labeled.{region_type}.{strand}.csv'
    output:
        'data/hg19/knownGene_labeled.{region_type}.{strand}.formated.bedgraph'
    shell:'''
    awk -F, '{{print $2 "\t" $4 "\t" $5 "\t" $1 "\t" $3 "\t" $6}}' {input} > {output} && [[ -s {output} ]]
    '''