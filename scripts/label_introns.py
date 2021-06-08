import csv
from remove_first_last_introns_from_intersection import group_introns_by_gene, write_introns_to_bed_file


def label_introns(intron_groups):
    # change intron names to genename_intronOrderFromStart_intronOrderFromEnd
    # ordering from start and end allows finding first and last introns
    # faster and easier in R
    # reverse numbering on reverse stranded introns

    sort_func = lambda i: int(i[1])
    for gene_name in intron_groups:
        introns = intron_groups[gene_name]
        strands = [each_intron[5] for each_intron in introns]
        assert len(set(strands)) == 1
        strand = strands[0]
        if strand == '+':
            introns = sorted(introns, key=sort_func)
        else:
            assert strand == '-'
            introns = sorted(introns, key=sort_func, reversed=True)
        
        # create intron names that 
        intron_names = [f'{gene_name}_{i}_{len(introns) - i}' for i, each_intron in enumerate(introns)]
        for each_intron, intron_name in zip(introns, intron_names):
            each_intron[3] = intron_name
        
        intron_groups[gene_name] = introns

    return intron_groups


def write_labled_introns_to_bed_file(labled_introns, output_path):
    with open(output_path, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(labled_introns)


def main():
    input_file = str(snakemake.input)
    output_file = str(snakemake.output)

    intron_group_dict = group_introns_by_gene(input_file)
    intron_group_dict = label_introns(intron_group_dict)

    write_introns_to_bed_file(intron_group_dict, output_file)


if __name__ == '__main__':
    main()
