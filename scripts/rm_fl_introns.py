
# Read in an intersected bed file with first 6 columns being a gene and
# next 6 being an intersected intron. Remove introns mapping to genes with < 3
# introns and for those with more remove the first and last introns.
# Write the remaining introns to a bed file

import csv

def group_introns_by_gene(intersection_filepath):
    intron_groups = {}
    with open(intersection_filepath) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            gene_name = row[3]
            intron = row[6:]
            if gene_name in intron_groups:
                intron_groups[gene_name].append(intron)
            else:
                intron_groups[gene_name] = [intron]

    assert intron_groups
    return intron_groups


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
        
        intron_names = [f'{gene_name}_{i+1}_{len(introns) - i}' for i, each_intron in enumerate(introns)]
        for each_intron, intron_name in zip(introns, intron_names):
            each_intron[3] = intron_name
        
        intron_groups[gene_name] = introns

    return intron_groups


def cull_introns(intron_groups):
    culled_intron_groups = {}
    for gene_name in intron_groups:
        introns = intron_groups[gene_name]
        if len(introns) < 3:
            continue  # remove introns mapping to genes with less than 3 mapped introns
        else:
            introns = sorted(introns, key=lambda i: int(i[1]))
            assert int(introns[0][1]) <= int(introns[-1][1])
            introns = introns[1:-1]  # remove first and last by position
            assert len(introns) > 0
            culled_intron_groups[gene_name] = introns
    
    assert culled_intron_groups
    return culled_intron_groups


def write_introns_to_bed_file(culled_introns_dict, output_filepath):
    with open(output_filepath, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for gene_name in culled_introns_dict:
            introns = culled_introns_dict[gene_name]
            assert len(introns) > 0
            writer.writerows(introns)
    return output_filepath


def main():
    input_file = str(snakemake.input)
    output_file_no_first_last = str(snakemake.output.no_first_last)


    intron_group_dict = group_introns_by_gene(input_file)
    intron_group_dict = label_introns(intron_group_dict)
    culled_intron_dict = cull_introns(intron_group_dict)

    write_introns_to_bed_file(culled_intron_dict, output_file_no_first_last)

if __name__ == '__main__':
    main()

