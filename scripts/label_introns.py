import csv
from rm_fl_introns import group_introns_by_gene, write_introns_to_bed_file, label_introns



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
