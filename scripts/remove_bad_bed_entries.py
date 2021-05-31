# remove bed file entries that DONT FOLLOW THE RULES

import csv

input_file = str(snakemake.input)
output_file = str(snakemake.output)

with open(input_file) as input_handle:
    reader = csv.reader(input_handle, delimiter='\t')
    with open(output_file, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        for row in reader:
            if row[1] >= row[2]:
                continue
            else:
                writer.writerow(row)



