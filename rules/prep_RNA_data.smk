# Download data from RNA structure score paper


rule download_processed_data_tar:
    output:
        compressed=temp('GSE149018_processed_data_files.tar.gz'),
        uncompressed='data/RNAss/done.txt',
        rep1='data/RNAss/WT_tNetSeq_rep1.bedgraph',
        rep2='data/RNAss/WT_tNetSeq_rep2.bedgraph',
        rep3='data/RNAss/WT_tNetSeq_rep3.bedgraph'
    shell:'''
    rm -rf data/RNAss/
    mkdir -p data/RNAss/
    cd data/RNAss/
    wget -O {output.compressed} \
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149018/suppl/GSE149018%5Fprocessed%5Fdata%5Ffiles%2Etar%2Egz"
    tar -xf {output.compressed}
    touch {output.uncompressed}
    '''
    