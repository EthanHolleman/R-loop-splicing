
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


def write_introns_to_bed_file(introns_dict, output_filepath):
    with open(output_filepath, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for gene_name in introns_dict:
            introns = introns_dict[gene_name]
            assert len(introns) > 0
            writer.writerows(introns)
    return output_filepath


def main():
    
    input_file = str(snakemake.input)
    output_file_no_first_last = str(snakemake.output.no_first_last)

    intron_group_dict = group_introns_by_gene(input_file)
    culled_intron_dict = cull_introns(intron_group_dict)

    write_introns_to_bed_file(culled_intron_dict, output_file_no_first_last)

if __name__ == '__main__':
    main()

